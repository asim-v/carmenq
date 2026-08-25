"""Extend an arbitrary-depth common-instrument spectral frontier.

Legacy two- and three-separator artifacts are normalised to a common schema:
each cell stores parallel ``branches`` and ``caps`` arrays, while the artifact
stores parallel separator coefficient and cube-face-grid arrays.  A new valid
real coefficient vector is crossed with every currently open cell through its
two scalar branches and a proved Bloch cover.

Distinct spectral branch tuples reuse one DPP SOCP; only cap parameters vary.
This removes model-construction growth with the number of angular cells, while
retaining the exact finite union required by the reverse-convex qubit norm.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import cvxpy as cp
import numpy as np

from fourier_behavior_cap_cover import cube_face_caps
from pairwise_inellipse_box_cover import deserialise_box
from terminal_reconstruction_enclosure import reconstruction_anchor_and_errors
from ternary_multicolumn_branch_tree import build_oracle, fixed_fourier_contractions
from ternary_shared_separator_cover import SOLVED_STATUSES, _cap


BRANCHES = ("scalar-positive", "scalar-negative", "bloch")


def normalise_frontier(
    payload: dict[str, Any], target: float
) -> tuple[tuple[np.ndarray, ...], tuple[int, ...], tuple[dict[str, Any], ...]]:
    """Return coefficient/grid arrays and canonical open spectral cells."""

    if "separator_coefficients" in payload:
        coefficients = tuple(
            np.asarray(item, dtype=float) for item in payload["separator_coefficients"]
        )
        grids = tuple(int(item) for item in payload["separator_grids"])
        cells = tuple(
            {
                "branches": tuple(str(item) for item in cell["branches"]),
                "caps": tuple(cell["caps"]),
                **{key: cell.get(key) for key in ("status", "bound", "audit", "return")},
            }
            for cell in payload.get("cells", ())
            if float(cell.get("bound", math.inf)) >= target
        )
    else:
        coefficients_list = [
            np.asarray(payload["first_separator"], dtype=float),
            np.asarray(payload["shared_second_separator"], dtype=float),
        ]
        grids_list = [int(payload["contraction_grid"])] * 2
        branch_fields = ["first_branch", "second_branch"]
        cap_fields = ["first_cap", "second_cap"]
        if "new_separator" in payload:
            coefficients_list.append(np.asarray(payload["new_separator"], dtype=float))
            grids_list.append(int(payload["new_grid"]))
            branch_fields.append("new_branch")
            cap_fields.append("new_cap")
        coefficients = tuple(coefficients_list)
        grids = tuple(grids_list)
        cells = tuple(
            {
                "branches": tuple(str(cell[field]) for field in branch_fields),
                "caps": tuple(cell.get(field) for field in cap_fields),
                **{key: cell.get(key) for key in ("status", "bound", "audit", "return")},
            }
            for cell in payload.get("cells", ())
            if float(cell.get("bound", math.inf)) >= target
        )
    if not coefficients or len(coefficients) != len(grids):
        raise ValueError("separator coefficient and grid arrays must be nonempty and parallel")
    for coefficient in coefficients:
        if coefficient.shape != (4,) or np.linalg.norm(coefficient) <= 1e-14:
            raise ValueError("every stored separator must be a nonzero four-vector")
    for cell in cells:
        if len(cell["branches"]) != len(coefficients) or len(cell["caps"]) != len(coefficients):
            raise ValueError("cell branch/cap depth does not match separator depth")
        for branch, cap in zip(cell["branches"], cell["caps"], strict=True):
            if branch not in BRANCHES:
                raise ValueError("unknown stored spectral branch")
            if (branch == "bloch") != (cap is not None):
                raise ValueError("Bloch branches, and only Bloch branches, require caps")
    return coefficients, grids, cells


def extend_frontier(
    payload: dict[str, Any],
    new_coefficients: np.ndarray,
    target: float = 0.758,
    new_grid: int = 2,
    capture_top: bool = False,
) -> dict[str, Any]:
    value = np.asarray(new_coefficients, dtype=float)
    if value.shape != (4,) or np.linalg.norm(value) <= 1e-14:
        raise ValueError("one nonzero four-component separator is required")
    value = value / np.linalg.norm(value)
    if new_grid < 2:
        raise ValueError("new_grid must be at least two")
    coefficients, grids, source_cells = normalise_frontier(payload, target)
    cap_families = tuple(cube_face_caps(grid) for grid in grids)
    new_caps = cube_face_caps(new_grid)
    box = deserialise_box(payload["box"])
    base = fixed_fourier_contractions(
        str(payload["base_code"]),
        _cap(payload.get("base_plane")),
        _cap(payload.get("base_sphere")),
    )
    anchor, errors, reconstruction_audit = reconstruction_anchor_and_errors(
        box["terminal_alpha"], box["terminal_beta"]
    )
    support_weight = float(payload["support_weight"])
    safety = float(payload["safety"])
    maximum_weight_floor = float(payload["maximum_weight_floor"])
    projective_support_upper = float(payload["projective_support_upper"])
    projective_support_lines = tuple(
        tuple(float(x) for x in line) for line in payload["projective_support_lines"]
    )
    prefix_order = tuple(int(x) for x in payload["prefix_order"])
    base_bloch_count = sum(item["branch"] == "bloch" for item in base)
    existing_branch_tuples = tuple(
        dict.fromkeys(tuple(cell["branches"]) for cell in source_cells)
    )

    oracles: dict[
        tuple[tuple[str, ...], str],
        tuple[Any, tuple[cp.Parameter | None, ...], cp.Parameter | None],
    ] = {}
    for existing_branches in existing_branch_tuples:
        for new_branch in BRANCHES:
            bloch_count = base_bloch_count
            items: list[dict[str, object]] = []
            parameters: list[cp.Parameter | None] = []
            for coefficient, branch in zip(coefficients, existing_branches, strict=True):
                parameter: cp.Parameter | None = None
                item: dict[str, object] = {
                    "coefficients": coefficient,
                    "branch": branch,
                }
                if branch == "bloch":
                    parameter = cp.Parameter(3)
                    item.update({"gauge_rank": bloch_count, "cap": parameter})
                    bloch_count += 1
                items.append(item)
                parameters.append(parameter)
            new_parameter: cp.Parameter | None = None
            new_item: dict[str, object] = {
                "coefficients": value,
                "branch": new_branch,
            }
            if new_branch == "bloch":
                new_parameter = cp.Parameter(3)
                new_item.update({"gauge_rank": bloch_count, "cap": new_parameter})
            oracle = build_oracle(
                support_weight,
                prefix_order,  # type: ignore[arg-type]
                maximum_weight_floor,
                projective_support_upper,
                projective_support_lines,  # type: ignore[arg-type]
                base + tuple(items) + (new_item,),
                (anchor, errors),
            )
            oracles[(existing_branches, new_branch)] = (
                oracle, tuple(parameters), new_parameter
            )

    rows: list[dict[str, Any]] = []
    maximum = -math.inf
    top_key: tuple[int, str, int | None] | None = None
    for source_index, cell in enumerate(source_cells):
        existing_branches = tuple(cell["branches"])
        for new_branch in BRANCHES:
            oracle, parameters, new_parameter = oracles[
                (existing_branches, new_branch)
            ]
            for parameter, cap_index, cap_family in zip(
                parameters, cell["caps"], cap_families, strict=True
            ):
                if parameter is not None and cap_index is not None:
                    normal, cosine = cap_family[int(cap_index)]
                    parameter.value = normal / cosine
            options: tuple[tuple[int | None, tuple[np.ndarray, float] | None], ...]
            options = tuple(enumerate(new_caps)) if new_parameter is not None else ((None, None),)
            for cap_index, cap in options:
                if new_parameter is not None and cap is not None:
                    new_parameter.value = cap[0] / cap[1]
                result = oracle.solve(box, safety)
                bound = float(result["bound"])
                row = {
                    "source_cell": source_index,
                    "branches": list(existing_branches + (new_branch,)),
                    "caps": list(tuple(cell["caps"]) + (cap_index,)),
                    "status": result["status"],
                    "bound": bound,
                    "audit": result.get("audit"),
                    "return": result.get("return"),
                }
                rows.append(row)
                if bound > maximum:
                    maximum = bound
                    top_key = (source_index, new_branch, cap_index)
        if (source_index + 1) % 25 == 0 or source_index + 1 == len(source_cells):
            print(
                json.dumps(
                    {
                        "source_cells_done": source_index + 1,
                        "source_cells_total": len(source_cells),
                        "solve_count": len(rows),
                        "maximum_bound": maximum,
                    }
                ),
                flush=True,
            )

    top_solution = None
    if capture_top and top_key is not None and math.isfinite(maximum):
        source_index, new_branch, cap_index = top_key
        cell = source_cells[source_index]
        existing_branches = tuple(cell["branches"])
        oracle, parameters, new_parameter = oracles[(existing_branches, new_branch)]
        for parameter, stored_index, cap_family in zip(
            parameters, cell["caps"], cap_families, strict=True
        ):
            if parameter is not None and stored_index is not None:
                normal, cosine = cap_family[int(stored_index)]
                parameter.value = normal / cosine
        if new_parameter is not None and cap_index is not None:
            normal, cosine = new_caps[cap_index]
            new_parameter.value = normal / cosine
        top_solution = oracle.solve(box, safety, capture=True)

    open_rows = [row for row in rows if float(row["bound"]) >= target]
    statuses_complete = all(row["status"] in SOLVED_STATUSES for row in rows)
    passthrough = {
        key: payload[key]
        for key in (
            "support_weight",
            "maximum_weight_floor",
            "projective_support_upper",
            "projective_support_lines",
            "prefix_order",
            "safety",
            "base_code",
            "box",
            "base_plane",
            "base_sphere",
        )
    }
    return {
        **passthrough,
        "target": target,
        "terminal_reconstruction": reconstruction_audit,
        "separator_coefficients": [item.tolist() for item in (*coefficients, value)],
        "separator_grids": list((*grids, new_grid)),
        "separator_depth": len(coefficients) + 1,
        "new_cap_count": len(new_caps),
        "new_covering_cosine": new_caps[0][1],
        "source_open_cell_count": len(source_cells),
        "crossed_solve_count": len(rows),
        "maximum_crossed_bound": maximum,
        "top_cell": None if not rows else max(rows, key=lambda row: float(row["bound"])),
        "open_crossed_cells": len(open_rows),
        "source_statuses_complete": bool(payload.get("statuses_complete", False)),
        "statuses_complete": statuses_complete,
        "complete": bool(payload.get("statuses_complete", False)) and statuses_complete and not open_rows,
        "cells": rows,
        "top_solution": top_solution,
        "scope": (
            "one terminal box and one fixed Fourier angular cell; arbitrary-depth "
            "common-instrument separator frontier extended by one exhaustive "
            "spectral/angular cover; solver-conditional"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--frontier-json", type=Path, required=True)
    parser.add_argument("--coefficients", type=float, nargs=4, required=True)
    parser.add_argument("--target", type=float, default=0.758)
    parser.add_argument("--new-grid", type=int, default=2)
    parser.add_argument("--capture-top", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = json.loads(args.frontier_json.read_text(encoding="utf-8"))
    result = extend_frontier(
        payload,
        np.asarray(args.coefficients, dtype=float),
        args.target,
        args.new_grid,
        args.capture_top,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "complete": result["complete"],
                "separator_depth": result["separator_depth"],
                "crossed_solve_count": result["crossed_solve_count"],
                "open_crossed_cells": result["open_crossed_cells"],
                "maximum_crossed_bound": result["maximum_crossed_bound"],
                "top_cell": result["top_cell"],
            }
        )
    )


if __name__ == "__main__":
    main()

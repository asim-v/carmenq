"""Cross one new separator with an accumulated ternary angular frontier.

This is the reusable continuation of ``ternary_shared_separator_cover.py``.
Its input stores the open product cells of two globally valid real
common-instrument contractions.  A third coefficient vector is supplied on
the command line and is covered by its two scalar branches plus a proved
cube-face Bloch cover.  DPP cap parameters let each distinct spectral branch
tuple reuse one canonical SOCP.

The result is a solver-conditional cover of exactly the source artifact's
fixed terminal box and Fourier angular cell.  A result marked ``complete``
also relies on the source artifact's already-closed cells.
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


def open_source_cells(
    payload: dict[str, Any], target: float
) -> tuple[dict[str, Any], ...]:
    cells = tuple(
        cell for cell in payload.get("cells", ())
        if float(cell.get("bound", math.inf)) >= target
    )
    for cell in cells:
        if cell.get("first_branch") not in BRANCHES:
            raise ValueError("unknown first source branch")
        if cell.get("second_branch") not in BRANCHES:
            raise ValueError("unknown second source branch")
        if cell["first_branch"] == "bloch" and cell.get("first_cap") is None:
            raise ValueError("a first Bloch branch is missing its cap")
        if cell["second_branch"] == "bloch" and cell.get("second_cap") is None:
            raise ValueError("a second Bloch branch is missing its cap")
    return cells


def cover_frontier(
    payload: dict[str, Any],
    coefficients: np.ndarray,
    target: float = 0.758,
    new_grid: int = 2,
    capture_top: bool = False,
) -> dict[str, Any]:
    value = np.asarray(coefficients, dtype=float)
    if value.shape != (4,) or np.linalg.norm(value) <= 1e-14:
        raise ValueError("one nonzero four-component separator is required")
    value = value / np.linalg.norm(value)
    if new_grid < 2:
        raise ValueError("new_grid must be at least two")
    source_cells = open_source_cells(payload, target)
    old_grid = int(payload["contraction_grid"])
    old_caps = cube_face_caps(old_grid)
    new_caps = cube_face_caps(new_grid)
    box = deserialise_box(payload["box"])
    base = fixed_fourier_contractions(
        str(payload["base_code"]),
        _cap(payload.get("base_plane")),
        _cap(payload.get("base_sphere")),
    )
    first = np.asarray(payload["first_separator"], dtype=float)
    second = np.asarray(payload["shared_second_separator"], dtype=float)
    reconstruction_anchor, reconstruction_errors, reconstruction_audit = (
        reconstruction_anchor_and_errors(
            box["terminal_alpha"], box["terminal_beta"]
        )
    )
    support_weight = float(payload["support_weight"])
    safety = float(payload["safety"])
    # The preceding driver inherits these constants from its adaptive tree.
    maximum_weight_floor = float(payload.get("maximum_weight_floor", 0.79))
    projective_support_upper = float(
        payload.get("projective_support_upper", 0.7573)
    )
    projective_support_lines = tuple(
        tuple(float(x) for x in line)
        for line in payload.get("projective_support_lines", ((0.6, 0.76591),))
    )
    prefix_order = tuple(int(x) for x in payload.get("prefix_order", (0, 1, 2, 3)))
    base_bloch_count = sum(item["branch"] == "bloch" for item in base)
    source_branch_pairs = tuple(
        dict.fromkeys(
            (str(cell["first_branch"]), str(cell["second_branch"]))
            for cell in source_cells
        )
    )

    oracles: dict[
        tuple[str, str, str],
        tuple[Any, cp.Parameter | None, cp.Parameter | None, cp.Parameter | None],
    ] = {}
    for first_branch, second_branch in source_branch_pairs:
        for new_branch in BRANCHES:
            bloch_count = base_bloch_count
            first_cap: cp.Parameter | None = None
            first_item: dict[str, object] = {
                "coefficients": first,
                "branch": first_branch,
            }
            if first_branch == "bloch":
                first_cap = cp.Parameter(3)
                first_item.update(
                    {"gauge_rank": bloch_count, "cap": first_cap}
                )
                bloch_count += 1
            second_cap: cp.Parameter | None = None
            second_item: dict[str, object] = {
                "coefficients": second,
                "branch": second_branch,
            }
            if second_branch == "bloch":
                second_cap = cp.Parameter(3)
                second_item.update(
                    {"gauge_rank": bloch_count, "cap": second_cap}
                )
                bloch_count += 1
            new_cap: cp.Parameter | None = None
            new_item: dict[str, object] = {
                "coefficients": value,
                "branch": new_branch,
            }
            if new_branch == "bloch":
                new_cap = cp.Parameter(3)
                new_item.update({"gauge_rank": bloch_count, "cap": new_cap})
            oracle = build_oracle(
                support_weight,
                prefix_order,  # type: ignore[arg-type]
                maximum_weight_floor,
                projective_support_upper,
                projective_support_lines,  # type: ignore[arg-type]
                base + (first_item, second_item, new_item),
                (reconstruction_anchor, reconstruction_errors),
            )
            oracles[(first_branch, second_branch, new_branch)] = (
                oracle, first_cap, second_cap, new_cap
            )

    rows: list[dict[str, Any]] = []
    maximum = -math.inf
    top_key: tuple[int, str, int | None] | None = None
    for source_index, cell in enumerate(source_cells):
        first_branch = str(cell["first_branch"])
        second_branch = str(cell["second_branch"])
        for new_branch in BRANCHES:
            oracle, first_cap, second_cap, new_cap = oracles[
                (first_branch, second_branch, new_branch)
            ]
            if first_cap is not None:
                normal, cosine = old_caps[int(cell["first_cap"])]
                first_cap.value = normal / cosine
            if second_cap is not None:
                normal, cosine = old_caps[int(cell["second_cap"])]
                second_cap.value = normal / cosine
            options: tuple[tuple[int | None, tuple[np.ndarray, float] | None], ...]
            options = tuple(enumerate(new_caps)) if new_cap is not None else ((None, None),)
            for cap_index, cap in options:
                if new_cap is not None and cap is not None:
                    new_cap.value = cap[0] / cap[1]
                result = oracle.solve(box, safety)
                bound = float(result["bound"])
                row = {
                    "source_cell": source_index,
                    "first_branch": first_branch,
                    "first_cap": cell.get("first_cap"),
                    "second_branch": second_branch,
                    "second_cap": cell.get("second_cap"),
                    "new_branch": new_branch,
                    "new_cap": cap_index,
                    "status": result["status"],
                    "bound": bound,
                    "audit": result.get("audit"),
                    "return": result.get("return"),
                }
                rows.append(row)
                if bound > maximum:
                    maximum = bound
                    top_key = (source_index, new_branch, cap_index)
        if (source_index + 1) % 10 == 0 or source_index + 1 == len(source_cells):
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
        key = (str(cell["first_branch"]), str(cell["second_branch"]), new_branch)
        oracle, first_cap, second_cap, new_cap = oracles[key]
        if first_cap is not None:
            normal, cosine = old_caps[int(cell["first_cap"])]
            first_cap.value = normal / cosine
        if second_cap is not None:
            normal, cosine = old_caps[int(cell["second_cap"])]
            second_cap.value = normal / cosine
        if new_cap is not None and cap_index is not None:
            normal, cosine = new_caps[cap_index]
            new_cap.value = normal / cosine
        top_solution = oracle.solve(box, safety, capture=True)

    open_rows = [row for row in rows if float(row["bound"]) >= target]
    statuses_complete = all(row["status"] in SOLVED_STATUSES for row in rows)
    source_statuses_complete = bool(payload.get("statuses_complete", False))
    return {
        "support_weight": support_weight,
        "maximum_weight_floor": maximum_weight_floor,
        "projective_support_upper": projective_support_upper,
        "projective_support_lines": [list(line) for line in projective_support_lines],
        "prefix_order": list(prefix_order),
        "target": target,
        "safety": safety,
        "base_code": payload["base_code"],
        "box": payload["box"],
        "base_plane": payload.get("base_plane"),
        "base_sphere": payload.get("base_sphere"),
        "terminal_reconstruction": reconstruction_audit,
        "contraction_grid": old_grid,
        "first_separator": first.tolist(),
        "shared_second_separator": second.tolist(),
        "new_separator": value.tolist(),
        "new_grid": new_grid,
        "new_cap_count": len(new_caps),
        "new_covering_cosine": new_caps[0][1],
        "source_cell_count": len(source_cells),
        "source_closed_cell_count": len(payload.get("cells", ())) - len(source_cells),
        "crossed_solve_count": len(rows),
        "maximum_crossed_bound": maximum,
        "top_cell": None if not rows else max(rows, key=lambda row: float(row["bound"])),
        "open_crossed_cells": len(open_rows),
        "source_statuses_complete": source_statuses_complete,
        "statuses_complete": statuses_complete,
        "complete": source_statuses_complete and statuses_complete and not open_rows,
        "cells": rows,
        "top_solution": top_solution,
        "scope": (
            "one terminal box and one fixed Fourier angular cell; a third "
            "global common-instrument separator crossed with every open cell "
            "of the two-separator source frontier; solver-conditional"
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
    result = cover_frontier(
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
                "crossed_solve_count": result["crossed_solve_count"],
                "open_crossed_cells": result["open_crossed_cells"],
                "maximum_crossed_bound": result["maximum_crossed_bound"],
                "top_cell": result["top_cell"],
            }
        )
    )


if __name__ == "__main__":
    main()

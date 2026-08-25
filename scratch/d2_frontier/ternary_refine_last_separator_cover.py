"""Refine only the open Bloch cells of the newest separator frontier.

Cube-face grids are nested whenever the fine grid is an integer multiple of
the coarse grid.  A coarse face-chart square is then the disjoint union of
``factor**2`` fine squares.  This driver replaces each open coarse Bloch cell
by exactly those children and re-solves open scalar cells unchanged.  It avoids
crossing every source cell with the entire fine spherical cover.
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


def cube_face_children(index: int, coarse_grid: int, fine_grid: int) -> tuple[int, ...]:
    """Return nested fine-grid indices covering one coarse cube-face cell."""

    if coarse_grid < 1 or fine_grid % coarse_grid:
        raise ValueError("fine_grid must be a positive multiple of coarse_grid")
    face, local = divmod(index, coarse_grid * coarse_grid)
    if not 0 <= face < 6:
        raise ValueError("coarse cap index is outside the cube-face cover")
    first, second = divmod(local, coarse_grid)
    factor = fine_grid // coarse_grid
    return tuple(
        face * fine_grid * fine_grid
        + (first * factor + di) * fine_grid
        + second * factor
        + dj
        for di in range(factor)
        for dj in range(factor)
    )


def refine_last(
    payload: dict[str, Any],
    target: float = 0.758,
    fine_grid: int = 4,
    capture_top: bool = False,
) -> dict[str, Any]:
    coarse_grid = int(payload["new_grid"])
    if fine_grid <= coarse_grid or fine_grid % coarse_grid:
        raise ValueError("fine_grid must be a larger multiple of the source grid")
    old_grid = int(payload["contraction_grid"])
    old_caps = cube_face_caps(old_grid)
    fine_caps = cube_face_caps(fine_grid)
    source_cells = tuple(
        cell for cell in payload.get("cells", ())
        if float(cell.get("bound", math.inf)) >= target
    )
    allowed_last = {"scalar-positive", "scalar-negative", "bloch"}
    if any(cell.get("new_branch") not in allowed_last for cell in source_cells):
        raise ValueError("source contains an unknown newest spectral branch")

    box = deserialise_box(payload["box"])
    base = fixed_fourier_contractions(
        str(payload["base_code"]),
        _cap(payload.get("base_plane")),
        _cap(payload.get("base_sphere")),
    )
    coefficients = (
        np.asarray(payload["first_separator"], dtype=float),
        np.asarray(payload["shared_second_separator"], dtype=float),
        np.asarray(payload["new_separator"], dtype=float),
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
    branch_triples = tuple(
        dict.fromkeys(
            (
                str(cell["first_branch"]),
                str(cell["second_branch"]),
                str(cell["new_branch"]),
            )
            for cell in source_cells
        )
    )
    oracles: dict[
        tuple[str, str, str],
        tuple[Any, cp.Parameter | None, cp.Parameter | None, cp.Parameter | None],
    ] = {}
    for branches in branch_triples:
        bloch_count = base_bloch_count
        items: list[dict[str, object]] = []
        parameters: list[cp.Parameter | None] = []
        for position, (coefficient, branch) in enumerate(zip(coefficients, branches, strict=True)):
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
        oracle = build_oracle(
            support_weight,
            prefix_order,  # type: ignore[arg-type]
            maximum_weight_floor,
            projective_support_upper,
            projective_support_lines,  # type: ignore[arg-type]
            base + tuple(items),
            (anchor, errors),
        )
        oracles[branches] = (oracle, *parameters)  # type: ignore[assignment]

    rows: list[dict[str, Any]] = []
    maximum = -math.inf
    top_key: tuple[int, int | None] | None = None
    for source_index, cell in enumerate(source_cells):
        branches = (
            str(cell["first_branch"]),
            str(cell["second_branch"]),
            str(cell["new_branch"]),
        )
        oracle, first_cap, second_cap, last_cap = oracles[branches]
        if first_cap is not None:
            normal, cosine = old_caps[int(cell["first_cap"])]
            first_cap.value = normal / cosine
        if second_cap is not None:
            normal, cosine = old_caps[int(cell["second_cap"])]
            second_cap.value = normal / cosine
        if last_cap is None:
            options: tuple[int | None, ...] = (None,)
        else:
            options = cube_face_children(
                int(cell["new_cap"]), coarse_grid, fine_grid
            )
        for fine_index in options:
            if last_cap is not None and fine_index is not None:
                normal, cosine = fine_caps[fine_index]
                last_cap.value = normal / cosine
            result = oracle.solve(box, safety)
            bound = float(result["bound"])
            row = {
                "source_cell": source_index,
                "first_branch": branches[0],
                "first_cap": cell.get("first_cap"),
                "second_branch": branches[1],
                "second_cap": cell.get("second_cap"),
                "new_branch": branches[2],
                "new_cap": fine_index,
                "status": result["status"],
                "bound": bound,
                "audit": result.get("audit"),
                "return": result.get("return"),
            }
            rows.append(row)
            if bound > maximum:
                maximum = bound
                top_key = (source_index, fine_index)
        if (source_index + 1) % 50 == 0 or source_index + 1 == len(source_cells):
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
        source_index, fine_index = top_key
        cell = source_cells[source_index]
        branches = (
            str(cell["first_branch"]),
            str(cell["second_branch"]),
            str(cell["new_branch"]),
        )
        oracle, first_cap, second_cap, last_cap = oracles[branches]
        if first_cap is not None:
            normal, cosine = old_caps[int(cell["first_cap"])]
            first_cap.value = normal / cosine
        if second_cap is not None:
            normal, cosine = old_caps[int(cell["second_cap"])]
            second_cap.value = normal / cosine
        if last_cap is not None and fine_index is not None:
            normal, cosine = fine_caps[fine_index]
            last_cap.value = normal / cosine
        top_solution = oracle.solve(box, safety, capture=True)

    open_rows = [row for row in rows if float(row["bound"]) >= target]
    statuses_complete = all(row["status"] in SOLVED_STATUSES for row in rows)
    return {
        **{
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
                "contraction_grid",
                "first_separator",
                "shared_second_separator",
                "new_separator",
            )
        },
        "target": target,
        "terminal_reconstruction": reconstruction_audit,
        "source_new_grid": coarse_grid,
        "new_grid": fine_grid,
        "new_cap_count": len(fine_caps),
        "new_covering_cosine": fine_caps[0][1],
        "source_open_cell_count": len(source_cells),
        "refined_solve_count": len(rows),
        "maximum_refined_bound": maximum,
        "top_cell": None if not rows else max(rows, key=lambda row: float(row["bound"])),
        "open_refined_cells": len(open_rows),
        "source_statuses_complete": bool(payload.get("statuses_complete", False)),
        "statuses_complete": statuses_complete,
        "complete": bool(payload.get("statuses_complete", False)) and statuses_complete and not open_rows,
        "cells": rows,
        "top_solution": top_solution,
        "scope": (
            "one terminal box and one fixed Fourier angular cell; exact nested "
            "cube-face refinement of only the newest open Bloch caps; "
            "solver-conditional"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--frontier-json", type=Path, required=True)
    parser.add_argument("--target", type=float, default=0.758)
    parser.add_argument("--fine-grid", type=int, default=4)
    parser.add_argument("--capture-top", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = json.loads(args.frontier_json.read_text(encoding="utf-8"))
    result = refine_last(payload, args.target, args.fine_grid, args.capture_top)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "complete": result["complete"],
                "refined_solve_count": result["refined_solve_count"],
                "open_refined_cells": result["open_refined_cells"],
                "maximum_refined_bound": result["maximum_refined_bound"],
                "top_cell": result["top_cell"],
            }
        )
    )


if __name__ == "__main__":
    main()

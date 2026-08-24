"""Exhaust one multicolumn contraction over one Fourier/pair cap cell."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import cvxpy as cp
import numpy as np

from fourier_behavior_cap_cover import cube_face_caps, plane_caps
from fourier_behavior_upper import solve_behavior_outer
from fourier_branch_upper import PRIOR_BOX


def solve(problem: cp.Problem) -> tuple[float, str]:
    try:
        problem.solve(
            solver="CLARABEL",
            tol_gap_abs=2e-8,
            tol_gap_rel=2e-8,
            tol_feas=2e-8,
            max_iter=1000,
            warm_start=True,
        )
    except cp.SolverError:
        return float("inf"), "clarabel_solver_error"
    if problem.status in {cp.INFEASIBLE, cp.INFEASIBLE_INACCURATE}:
        return float("-inf"), problem.status
    if problem.status not in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}:
        return float("inf"), problem.status
    return float(problem.value), problem.status


def run_cover(
    coefficients: np.ndarray,
    plane_cells: int,
    plane_index: int,
    face_grid: int,
    sphere_index: int,
    pair_cap_index: int,
    contraction_grid: int,
    parent_contractions: tuple[dict[str, object], ...] = (),
    target: float = 0.758,
    pair_branch: str = "bloch",
) -> dict[str, object]:
    coefficients = np.asarray(coefficients, dtype=float)
    if coefficients.shape != (4,) or np.linalg.norm(coefficients) <= 1e-14:
        raise ValueError("one nonzero four-component coefficient vector is required")
    coefficients = coefficients / np.linalg.norm(coefficients)
    plane = plane_caps(plane_cells)[plane_index]
    sphere = cube_face_caps(face_grid)[sphere_index]
    pair_axes = [
        (normal, 1.0 / np.sqrt(3.0))
        for axis in range(3)
        for sign in (-1.0, 1.0)
        for normal in [np.eye(3)[axis] * sign]
    ]
    common_caps = (
        None,
        (*plane[0], plane[1]),
        (*sphere[0], sphere[1]),
    )
    base_cut = {
        "pair": (2, 3),
        "scale": 1.0,
        "branch": pair_branch,
    }
    if pair_branch == "bloch":
        pair_cap = pair_axes[pair_cap_index]
        base_cut["cap"] = (*pair_cap[0], pair_cap[1])
    elif pair_branch != "scalar-positive":
        raise ValueError("the base pair branch must be scalar-positive or bloch")
    base_cuts = (base_cut, *parent_contractions)

    scalar_models = []
    for branch in ("scalar-positive", "scalar-negative"):
        scalar_models.append(
            (
                branch,
                solve_behavior_outer(
                    ("bloch", "bloch", "bloch"),
                    PRIOR_BOX,
                    common_caps,
                    build_only=True,
                    pairwise_contractions=(
                        *base_cuts,
                        {"coefficients": coefficients, "branch": branch},
                    ),
                ),
            )
        )
    scaled_cap = cp.Parameter(3)
    vector_model = solve_behavior_outer(
        ("bloch", "bloch", "bloch"),
        PRIOR_BOX,
        common_caps,
        build_only=True,
        pairwise_contractions=(
            *base_cuts,
            {
                "coefficients": coefficients,
                "branch": "bloch",
                "cap": scaled_cap,
            },
        ),
    )

    rows: list[dict[str, object]] = []
    for branch, model in scalar_models:
        bound, status = solve(model["problem"])
        rows.append({"branch": branch, "cap": None, "bound": bound, "status": status})
    caps = cube_face_caps(contraction_grid)
    for index, (normal, cosine) in enumerate(caps):
        scaled_cap.value = normal / cosine
        bound, status = solve(vector_model["problem"])
        rows.append(
            {
                "branch": "bloch",
                "cap": index,
                "bound": bound,
                "status": status,
            }
        )
    maximum = max(rows, key=lambda item: float(item["bound"]))
    open_rows = [item for item in rows if float(item["bound"]) >= target]
    return {
        "scope": "one fixed Fourier/pair cell with one exhaustive multicolumn contraction",
        "target": target,
        "coefficients": coefficients.tolist(),
        "parent_contractions": list(parent_contractions),
        "base_cell": {
            "plane_cells": plane_cells,
            "plane_index": plane_index,
            "face_grid": face_grid,
            "sphere_index": sphere_index,
            "pair_branch": pair_branch,
            "pair_cap_index": pair_cap_index,
        },
        "contraction_grid": contraction_grid,
        "contraction_caps": len(caps),
        "contraction_covering_cosine": caps[0][1],
        "branch_count": len(rows),
        "maximum_bound": maximum["bound"],
        "maximum_branch": maximum,
        "closed_branch_count": len(rows) - len(open_rows),
        "open_branch_count": len(open_rows),
        "complete_cell_closure": not open_rows,
        "branches": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("separator", type=Path)
    parser.add_argument("--plane-cells", type=int, default=8)
    parser.add_argument("--plane-index", type=int, default=4)
    parser.add_argument("--face-grid", type=int, default=4)
    parser.add_argument("--sphere-index", type=int, default=18)
    parser.add_argument("--pair-cap-index", type=int, default=3)
    parser.add_argument("--contraction-grid", type=int, default=4)
    parser.add_argument(
        "--parent-spec",
        type=Path,
        help="JSON list of fixed parent contraction dictionaries",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    separator = json.loads(args.separator.read_text(encoding="utf-8"))
    parent_contractions = (
        ()
        if args.parent_spec is None
        else tuple(json.loads(args.parent_spec.read_text(encoding="utf-8")))
    )
    payload = run_cover(
        np.asarray(separator["coefficients"], dtype=float),
        args.plane_cells,
        args.plane_index,
        args.face_grid,
        args.sphere_index,
        args.pair_cap_index,
        args.contraction_grid,
        parent_contractions,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "maximum_bound": payload["maximum_bound"],
                "maximum_branch": payload["maximum_branch"],
                "closed_branch_count": payload["closed_branch_count"],
                "open_branch_count": payload["open_branch_count"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

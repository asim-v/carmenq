"""Fourier cap cover augmented by one exact pairwise channel contraction."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import cvxpy as cp
import numpy as np

from fourier_behavior_cap_cover import cube_face_caps, plane_caps
from fourier_behavior_upper import solve_behavior_outer
from fourier_branch_upper import PRIOR_BOX


def solve_parameterised(problem: cp.Problem) -> tuple[float, str]:
    problem.solve(
        solver="CLARABEL",
        tol_gap_abs=2e-8,
        tol_gap_rel=2e-8,
        tol_feas=2e-8,
        max_iter=1000,
        warm_start=True,
    )
    if problem.status not in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}:
        return float("inf"), problem.status
    return float(problem.value), problem.status


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--plane-cells", type=int, default=8)
    parser.add_argument("--face-grid", type=int, default=4)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    plane_scaled = cp.Parameter(3)
    sphere_scaled = cp.Parameter(3)
    pair_scaled = cp.Parameter(3)
    common_caps = (None, plane_scaled, sphere_scaled)
    scalar_model = solve_behavior_outer(
        ("bloch", "bloch", "bloch"),
        PRIOR_BOX,
        common_caps,
        build_only=True,
        pairwise_contractions=(
            {"pair": (2, 3), "scale": 1.0, "branch": "scalar-positive"},
        ),
    )
    vector_model = solve_behavior_outer(
        ("bloch", "bloch", "bloch"),
        PRIOR_BOX,
        common_caps,
        build_only=True,
        pairwise_contractions=(
            {
                "pair": (2, 3),
                "scale": 1.0,
                "branch": "bloch",
                "cap": pair_scaled,
            },
        ),
    )
    planes = plane_caps(args.plane_cells)
    spheres = cube_face_caps(args.face_grid)
    axes = [
        (normal, 1.0 / np.sqrt(3.0))
        for axis in range(3)
        for sign in (-1.0, 1.0)
        for normal in [np.eye(3)[axis] * sign]
    ]
    maximum = -np.inf
    maximum_cell: dict[str, object] | None = None
    rows = []
    for plane_index, (plane_normal, plane_cosine) in enumerate(planes):
        plane_scaled.value = plane_normal / plane_cosine
        for sphere_index, (sphere_normal, sphere_cosine) in enumerate(spheres):
            sphere_scaled.value = sphere_normal / sphere_cosine
            bound, status = solve_parameterised(scalar_model["problem"])
            candidates = [(bound, status, "scalar-positive", None)]
            for pair_index, (pair_normal, pair_cosine) in enumerate(axes):
                pair_scaled.value = pair_normal / pair_cosine
                bound, status = solve_parameterised(vector_model["problem"])
                candidates.append((bound, status, "bloch", pair_index))
            bound, status, pair_branch, pair_index = max(
                candidates, key=lambda item: item[0]
            )
            row = {
                "plane": plane_index,
                "sphere": sphere_index,
                "bound": bound,
                "status": status,
                "pair_branch": pair_branch,
                "pair_cap": pair_index,
                "branches": [
                    {
                        "bound": candidate_bound,
                        "status": candidate_status,
                        "pair_branch": candidate_branch,
                        "pair_cap": candidate_cap,
                    }
                    for (
                        candidate_bound,
                        candidate_status,
                        candidate_branch,
                        candidate_cap,
                    ) in candidates
                ],
            }
            rows.append(row)
            if bound > maximum:
                maximum = bound
                maximum_cell = row.copy()
        print(
            json.dumps(
                {
                    "plane": plane_index,
                    "running_maximum": maximum,
                }
            ),
            flush=True,
        )
    payload = {
        "scope": "bbb Fourier cover plus the (2,3), t=1 flagged contraction",
        "target": 0.758,
        "plane_cells": args.plane_cells,
        "plane_covering_cosine": planes[0][1],
        "cube_face_grid": args.face_grid,
        "sphere_caps": len(spheres),
        "sphere_covering_cosine": spheres[0][1],
        "pair_vector_caps": len(axes),
        "pair_covering_cosine": axes[0][1],
        "cell_count": len(rows),
        "conic_solve_count": len(rows) * (1 + len(axes)),
        "maximum_bound": maximum,
        "closed": maximum < 0.758,
        "maximum_cell": maximum_cell,
        "cells": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

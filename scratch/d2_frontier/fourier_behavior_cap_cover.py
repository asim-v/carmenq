"""Dense, provable angular cover using one parameterised SOCP."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import cvxpy as cp
import numpy as np

from fourier_behavior_upper import solve_behavior_outer
from fourier_branch_upper import PRIOR_BOX


def plane_caps(count: int) -> list[tuple[np.ndarray, float]]:
    if count < 2:
        raise ValueError("at least two planar caps are required")
    cosine = float(np.cos(np.pi / (2.0 * count)))
    return [
        (
            np.asarray([np.cos(angle), 0.0, np.sin(angle)]),
            cosine,
        )
        for angle in (
            -np.pi / 2.0 + (index + 0.5) * np.pi / count
            for index in range(count)
        )
    ]


def cube_face_caps(grid: int) -> list[tuple[np.ndarray, float]]:
    """Cover S^2 by normalised cube-face cells with a proved cosine."""

    if grid < 1:
        raise ValueError("face grid must be positive")
    centers = -1.0 + (np.arange(grid) + 0.5) * 2.0 / grid
    cosine = 1.0 - 1.0 / grid**2
    result: list[tuple[np.ndarray, float]] = []
    for axis in range(3):
        other = [coordinate for coordinate in range(3) if coordinate != axis]
        for sign in (-1.0, 1.0):
            for first in centers:
                for second in centers:
                    vector = np.zeros(3)
                    vector[axis] = sign
                    # Ratios are taken relative to the signed dominant
                    # coordinate, so every face uses the same [-1,1]^2 chart.
                    vector[other[0]] = sign * first
                    vector[other[1]] = sign * second
                    vector /= np.linalg.norm(vector)
                    result.append((vector, cosine))
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--plane-cells", type=int, default=16)
    parser.add_argument("--face-grid", type=int, default=8)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    plane_scaled = cp.Parameter(3)
    sphere_scaled = cp.Parameter(3)
    model = solve_behavior_outer(
        ("bloch", "bloch", "bloch"),
        PRIOR_BOX,
        (None, plane_scaled, sphere_scaled),
        build_only=True,
    )
    problem = model["problem"]
    planes = plane_caps(args.plane_cells)
    spheres = cube_face_caps(args.face_grid)
    rows: list[dict[str, object]] = []
    maximum = -np.inf
    maximum_cell: dict[str, object] | None = None
    for plane_index, (plane_normal, plane_cosine) in enumerate(planes):
        plane_scaled.value = plane_normal / plane_cosine
        for sphere_index, (sphere_normal, sphere_cosine) in enumerate(spheres):
            sphere_scaled.value = sphere_normal / sphere_cosine
            problem.solve(
                solver="CLARABEL",
                tol_gap_abs=2e-8,
                tol_gap_rel=2e-8,
                tol_feas=2e-8,
                max_iter=1000,
                warm_start=True,
            )
            if problem.status not in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}:
                bound = float("inf")
            else:
                bound = float(problem.value)
            row = {
                "plane": plane_index,
                "sphere": sphere_index,
                "bound": bound,
                "status": problem.status,
            }
            rows.append(row)
            if bound > maximum:
                maximum = bound
                maximum_cell = {
                    **row,
                    "plane_normal": plane_normal.tolist(),
                    "plane_cosine": plane_cosine,
                    "sphere_normal": sphere_normal.tolist(),
                    "sphere_cosine": sphere_cosine,
                    "prior": [
                        float(item.value) for item in model["prior_expressions"]
                    ],
                    "audit": float(model["audit_expression"].value),
                    "return": float(model["return_expression"].value),
                    "flagged_norms": [
                        float(item.value) for item in model["flagged_expressions"]
                    ],
                    "input_fourier_vectors": [
                        np.asarray(item.value, dtype=float).tolist()
                        for item in model["input_fourier_expressions"]
                    ],
                }
        print(
            json.dumps(
                {
                    "plane": plane_index,
                    "completed": len(spheres),
                    "running_maximum": maximum,
                }
            ),
            flush=True,
        )
    payload = {
        "scope": "parameterised Fourier-behavior cap cover of the bbb branch",
        "target": 0.758,
        "plane_cells": args.plane_cells,
        "plane_covering_cosine": planes[0][1],
        "cube_face_grid": args.face_grid,
        "sphere_caps": len(spheres),
        "sphere_covering_cosine": spheres[0][1],
        "cell_count": len(rows),
        "maximum_bound": maximum,
        "closed": maximum < 0.758,
        "maximum_cell": maximum_cell,
        "cells": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

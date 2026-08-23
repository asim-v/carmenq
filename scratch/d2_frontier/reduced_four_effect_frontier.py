"""Symmetric four-effect branch of the interleaved qubit frontier.

The coarse effect has spectrum ``(p,1-p)``.  Each binary group is split by
congruence with a rank-one projector; the complementary group is its bit-flip
image.  Two symmetry-related pairs of pure prefix states suffice.  Optimising
their total prior weights is the largest eigenvalue of a 2-by-2
diagonal-plus-rank-one matrix.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import differential_evolution, minimize


def evaluate(
    p: float,
    theta: float,
    u_small: float,
    u_large: float,
    weight: float,
) -> dict[str, object]:
    root_p = math.sqrt(max(0.0, p))
    root_one_p = math.sqrt(max(0.0, 1.0 - p))
    sine = math.sin(theta)
    cosine = math.cos(theta)

    q_small = (1.0 - p) + (2.0 * p - 1.0) * u_small
    q_large = (1.0 - p) + (2.0 * p - 1.0) * u_large
    d_small = (
        root_p * sine * math.sqrt(max(0.0, u_small))
        + root_one_p * cosine * math.sqrt(max(0.0, 1.0 - u_small))
    ) ** 2
    d_large = (
        root_p * cosine * math.sqrt(max(0.0, u_large))
        + root_one_p * sine * math.sqrt(max(0.0, 1.0 - u_large))
    ) ** 2
    c_small = math.sqrt(max(0.0, q_small)) + math.sqrt(
        max(0.0, 1.0 - q_small)
    )
    c_large = math.sqrt(max(0.0, q_large)) + math.sqrt(
        max(0.0, 1.0 - q_large)
    )

    hellinger = np.asarray((c_small, c_large), dtype=float)
    matrix = weight * np.diag((d_small, d_large))
    matrix += (1.0 - weight) * np.outer(hellinger, hellinger) / 4.0
    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    vector = np.abs(eigenvectors[:, -1])
    vector /= np.linalg.norm(vector)
    small_total = float(vector[0] ** 2)
    large_total = float(vector[1] ** 2)
    audit = small_total * d_small + large_total * d_large
    returned = (vector[0] * c_small + vector[1] * c_large) ** 2 / 4.0
    score = weight * audit + (1.0 - weight) * returned
    return {
        "weight": weight,
        "p": p,
        "theta": theta,
        "u_small": u_small,
        "u_large": u_large,
        "score": float(score),
        "audit": float(audit),
        "return": float(returned),
        "priors": [
            small_total / 2.0,
            large_total / 2.0,
            large_total / 2.0,
            small_total / 2.0,
        ],
        "q_small": q_small,
        "q_large": q_large,
        "d_small": d_small,
        "d_large": d_large,
        "c_small": c_small,
        "c_large": c_large,
        "largest_eigenvalue": float(eigenvalues[-1]),
    }


def optimise(weight: float) -> dict[str, object]:
    def objective(point: np.ndarray) -> float:
        p, angle_coordinate, u_small, u_large = map(float, point)
        return -float(
            evaluate(
                p,
                angle_coordinate * math.pi / 2.0,
                u_small,
                u_large,
                weight,
            )["score"]
        )

    global_result = differential_evolution(
        objective,
        bounds=((0.5, 1.0), (0.0, 1.0), (0.0, 1.0), (0.0, 1.0)),
        seed=20260822,
        popsize=32,
        maxiter=1400,
        tol=2e-12,
        polish=False,
    )
    starts = [
        global_result.x,
        np.asarray((0.952, 0.0142, 0.477, 0.998)),
        np.asarray((0.90, 0.03, 0.50, 0.98)),
        np.asarray((0.75, 0.10, 0.50, 0.90)),
    ]
    candidates = []
    for start in starts:
        local = minimize(
            objective,
            start,
            method="Nelder-Mead",
            bounds=((0.5, 1.0), (0.0, 1.0), (0.0, 1.0), (0.0, 1.0)),
            options={"xatol": 2e-13, "fatol": 2e-15, "maxiter": 50000},
        )
        p, angle_coordinate, u_small, u_large = map(float, local.x)
        candidates.append(
            evaluate(
                p,
                angle_coordinate * math.pi / 2.0,
                u_small,
                u_large,
                weight,
            )
        )
    return max(candidates, key=lambda row: float(row["score"]))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lambda", dest="weight", type=float, default=0.6)
    parser.add_argument("--sweep", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    if args.sweep:
        payload: object = [
            optimise(weight)
            for weight in (0.45, 0.48, 0.5, 0.52, 0.55, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99)
        ]
    else:
        payload = optimise(args.weight)
    rendered = json.dumps(payload, indent=2) + "\n"
    print(rendered)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

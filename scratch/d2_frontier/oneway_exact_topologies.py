"""Direct support optimisation of every extreme one-way qubit topology.

The four correct-decision effects are grouped into two binary POVMs with
coarse effects ``E`` and ``I-E``.  For a fixed coarse effect, convexity of the
Perron secular functional reduces an optimal binary split to either an
endpoint ``(E,0)`` or a congruence rank-one split
``E^(1/2) (P,I-P) E^(1/2)``.  This file optimises the remaining scalar state
coordinates and reports the Perron-optimal prefix priors.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import differential_evolution, minimize


def state_data(
    high: float,
    low: float,
    topology: str,
    coordinates: tuple[float, ...],
) -> tuple[list[float], list[float]]:
    if topology == "endpoint":
        active_coordinate, null_coordinate = coordinates
        q_active = low + (high - low) * active_coordinate
        q_null = low + (high - low) * null_coordinate
        return [q_active, q_null], [q_active, 0.0]

    angle_coordinate, first_coordinate, second_coordinate = coordinates
    angle = (math.pi / 2.0) * angle_coordinate
    cosine = math.cos(angle)
    sine = math.sin(angle)

    def pair(coordinate: float, swapped: bool) -> tuple[float, float]:
        q_value = low + (high - low) * coordinate
        first = math.sqrt(max(0.0, high * coordinate))
        second = math.sqrt(max(0.0, low * (1.0 - coordinate)))
        if swapped:
            d_value = (sine * first + cosine * second) ** 2
        else:
            d_value = (cosine * first + sine * second) ** 2
        return q_value, d_value

    first = pair(first_coordinate, False)
    second = pair(second_coordinate, True)
    return [first[0], second[0]], [first[1], second[1]]


def evaluate(
    point: np.ndarray,
    first_topology: str,
    second_topology: str,
    weight: float,
) -> dict[str, object]:
    high = float(point[0])
    low = high * float(point[1])
    cursor = 2
    first_count = 2 if first_topology == "endpoint" else 3
    first_coordinates = tuple(map(float, point[cursor : cursor + first_count]))
    cursor += first_count
    second_count = 2 if second_topology == "endpoint" else 3
    second_coordinates = tuple(map(float, point[cursor : cursor + second_count]))

    first_q, first_d = state_data(
        high, low, first_topology, first_coordinates
    )
    second_q, second_d = state_data(
        1.0 - low, 1.0 - high, second_topology, second_coordinates
    )
    q_values = np.asarray(first_q + second_q, dtype=float)
    d_values = np.asarray(first_d + second_d, dtype=float)
    hellinger = np.sqrt(np.maximum(q_values, 0.0)) + np.sqrt(
        np.maximum(1.0 - q_values, 0.0)
    )
    matrix = weight * np.diag(d_values)
    matrix += (1.0 - weight) * np.outer(hellinger, hellinger) / 8.0
    values, vectors = np.linalg.eigh(matrix)
    vector = np.abs(vectors[:, -1])
    vector /= np.linalg.norm(vector)
    priors = vector**2
    audit = float(priors @ d_values)
    returned = float((vector @ hellinger) ** 2 / 8.0)
    return {
        "weight": weight,
        "score": float(values[-1]),
        "audit": audit,
        "return": returned,
        "high_eigenvalue": high,
        "low_eigenvalue": low,
        "first_topology": first_topology,
        "second_topology": second_topology,
        "first_coordinates": list(first_coordinates),
        "second_coordinates": list(second_coordinates),
        "q": q_values.tolist(),
        "d": d_values.tolist(),
        "priors": priors.tolist(),
        "stationarity_residual": abs(
            float(values[-1]) - (weight * audit + (1.0 - weight) * returned)
        ),
    }


def optimise_topology(
    first: str,
    second: str,
    weight: float,
    seed: int,
) -> dict[str, object]:
    dimension = 2 + (2 if first == "endpoint" else 3) + (
        2 if second == "endpoint" else 3
    )

    def objective(point: np.ndarray) -> float:
        return -float(evaluate(point, first, second, weight)["score"])

    result = differential_evolution(
        objective,
        bounds=((1e-10, 1.0 - 1e-10),) + ((0.0, 1.0),) * (dimension - 1),
        seed=seed,
        popsize=36,
        maxiter=1800,
        tol=2e-12,
        polish=False,
        updating="immediate",
    )
    starts = [result.x]
    # Known basins: endpoint/rank (3E) and symmetric rank/rank (4E).
    if first == "endpoint" and second == "rank":
        starts.append(np.asarray((0.46, 0.0, 1.0, 0.5, 0.5, 0.5, 0.5)))
    if first == "rank" and second == "rank":
        starts.append(
            np.asarray((0.952, 0.05, 0.014, 0.998, 0.477, 0.014, 0.998, 0.477))
        )
    candidates = []
    bounds = ((1e-12, 1.0 - 1e-12),) + ((0.0, 1.0),) * (dimension - 1)
    for start in starts:
        if len(start) != dimension:
            continue
        local = minimize(
            objective,
            start,
            method="Nelder-Mead",
            bounds=bounds,
            options={"xatol": 2e-13, "fatol": 2e-15, "maxiter": 100000},
        )
        candidates.append(evaluate(local.x, first, second, weight))
    return max(candidates, key=lambda row: float(row["score"]))


def optimise(weight: float) -> list[dict[str, object]]:
    rows = []
    for first_index, first in enumerate(("endpoint", "rank")):
        for second_index, second in enumerate(("endpoint", "rank")):
            rows.append(
                optimise_topology(
                    first,
                    second,
                    weight,
                    20260822 + 17 * first_index + second_index,
                )
            )
    rows.sort(key=lambda row: float(row["score"]), reverse=True)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lambda", dest="weight", type=float, default=0.6)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = {"weight": args.weight, "topologies": optimise(args.weight)}
    rendered = json.dumps(payload, indent=2) + "\n"
    print(rendered)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

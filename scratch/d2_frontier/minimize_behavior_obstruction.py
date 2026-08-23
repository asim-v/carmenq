"""Greedily minimize a nested-ellipsoid qubit-behaviour obstruction.

The full four-preparation certificate may use every outcome column.  Each
supported column becomes one branch of the resulting valid disjunction, so a
small irreducible support is materially cheaper than the raw certificate.
This script tries many deletion orders, caches every conic solve, and archives
the strongest smallest support it finds.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path
from typing import Any

import cvxpy as cp
import numpy as np

from full_behavior_psd_rank_certificate import solve_dual, solve_primal


def load_behavior(path: Path, leaf: int) -> np.ndarray:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if isinstance(payload, dict) and "behavior" in payload:
        behavior = payload["behavior"]
    elif isinstance(payload, dict) and "leaves" in payload:
        behavior = payload["leaves"][leaf]["behavior"]
    else:
        raise ValueError("input must contain behavior or leaves[*].behavior")
    value = np.asarray(behavior, dtype=float)
    if value.ndim != 2 or value.shape[0] != 4:
        raise ValueError("expected a 4 x m behavior")
    return value


def valid_dual(payload: dict[str, Any], tolerance: float) -> bool:
    return (
        payload.get("status") in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}
        and float(payload.get("certified_common_margin", -math.inf)) > tolerance
        and float(payload.get("state_dual_min_eigenvalue", -math.inf)) > -tolerance
        and min(payload.get("containment_dual_min_eigenvalues", [-math.inf]))
        > -tolerance
        and float(payload.get("stationarity_frobenius_residual", math.inf))
        < 10.0 * tolerance
    )


def normalized_quality(payload: dict[str, Any]) -> float:
    coefficients = np.asarray(payload["halfspace_linear_coefficients"], dtype=float)
    scale = max(float(np.max(np.sum(np.abs(coefficients), axis=1))), 1e-300)
    return float(payload["certified_common_margin"]) / scale


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--leaf", type=int, default=0)
    parser.add_argument("--orders", type=int, default=32)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--robust-budget", type=float, default=1000.0)
    parser.add_argument("--tolerance", type=float, default=2e-8)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    behavior = load_behavior(args.input, args.leaf)
    columns = tuple(range(behavior.shape[1]))
    cache: dict[tuple[int, ...], dict[str, Any]] = {}

    def audit(selected: tuple[int, ...]) -> dict[str, Any]:
        key = tuple(sorted(selected))
        if key not in cache:
            cache[key] = solve_dual(
                behavior,
                list(key),
                robust_budget=args.robust_budget,
            )
        return cache[key]

    root = audit(columns)
    if not valid_dual(root, args.tolerance):
        raise RuntimeError("the full column set has no validated obstruction")

    generator = np.random.default_rng(args.seed)
    candidates: list[tuple[tuple[int, ...], dict[str, Any]]] = []
    for order_index in range(args.orders):
        active = list(columns)
        order = list(columns)
        if order_index:
            generator.shuffle(order)
        changed = True
        while changed:
            changed = False
            for column in order:
                if column not in active or len(active) <= 1:
                    continue
                trial = tuple(item for item in active if item != column)
                dual = audit(trial)
                if valid_dual(dual, args.tolerance):
                    active.remove(column)
                    changed = True
        selected = tuple(sorted(active))
        candidates.append((selected, audit(selected)))
        print(
            {
                "order": order_index,
                "support": list(selected),
                "size": len(selected),
                "quality": normalized_quality(candidates[-1][1]),
                "cached_solves": len(cache),
            },
            flush=True,
        )

    validated_small_supports: list[tuple[tuple[int, ...], dict[str, Any]]] = []
    for size in (1, 2):
        for selected in itertools.combinations(columns, size):
            dual = audit(selected)
            if valid_dual(dual, args.tolerance):
                validated_small_supports.append((selected, dual))
    candidates.extend(validated_small_supports)
    candidates.sort(key=lambda item: (len(item[0]), -normalized_quality(item[1])))
    selected, dual = candidates[0]
    primal = solve_primal(behavior, list(selected))
    payload = {
        "source": str(args.input),
        "leaf": args.leaf,
        "behavior": behavior.tolist(),
        "orders": args.orders,
        "conic_solves": len(cache),
        "selected_columns": list(selected),
        "support_size": len(selected),
        "all_singletons_and_pairs_checked": True,
        "validated_smaller_support_count": len(validated_small_supports),
        "minimum_cardinality_certified": (
            len(selected) == 3 and not validated_small_supports
        ),
        "primal_status": primal["status"],
        "robust_quality": normalized_quality(dual),
        "dual": dual,
        "scope_note": (
            "Greedy irreducibility is not a proof of minimum cardinality. "
            "The archived conic witness is a numerical certificate until "
            "its PSD data are outward-rounded or interval-validated."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "selected_columns": payload["selected_columns"],
                "support_size": payload["support_size"],
                "primal_status": payload["primal_status"],
                "robust_quality": payload["robust_quality"],
                "minimum_cardinality_certified": payload[
                    "minimum_cardinality_certified"
                ],
                "conic_solves": payload["conic_solves"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

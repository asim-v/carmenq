"""Global probability upper bound for fixed terminal rank-one POVM weights.

If ``P_s = w_s |n_s><n_s|`` is a rank-one qubit POVM, then ``0 <= w_s <= 1``
and ``sum_s w_s = 2``.  For every path state,

    Tr(P_s sigma_zy) <= w_s p_zy.

Moreover, pulling the terminal POVM back to the prefix qubit gives a POVM
whose every effect has operator norm at most ``w_max``.  Its maximum success
for prefix priors ``a`` is bounded by filling total effect trace two, in
decreasing prior order, with per-effect cap ``w_max``.  The script combines
these two AUDIT bounds with the exact Hellinger return and globally solves the
resulting probability relaxation by enumerating the 24 prefix-prior orders.
"""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import cvxpy as cp
import numpy as np


PERMUTATIONS = tuple(itertools.permutations(range(4)))


def filled_effect_weights(maximum: float) -> np.ndarray:
    remaining = 2.0
    values = []
    for _ in range(4):
        value = min(maximum, remaining)
        values.append(value)
        remaining -= value
    if abs(remaining) > 1e-10:
        raise ValueError("effect-norm cap is too small to complete a qubit POVM")
    return np.asarray(values)


def hellinger_return(probabilities: cp.Expression) -> cp.Expression:
    flat = cp.reshape(probabilities, (16,), order="C")
    cross = [
        cp.geo_mean(cp.hstack((flat[i], flat[j])))
        for i in range(16)
        for j in range(i + 1, 16)
    ]
    return (cp.sum(flat) + 2.0 * cp.sum(cp.hstack(cross))) / 16.0


def solve_order(
    terminal_weights: np.ndarray,
    order: tuple[int, int, int, int],
    weight: float,
    audit_upper: float | None = None,
    projective_support_upper: float | None = None,
    projective_support_lines: tuple[tuple[float, float], ...] = (),
) -> dict[str, object]:
    probabilities = cp.Variable((4, 4), nonneg=True)
    audit = cp.Variable()
    returned = cp.Variable(nonneg=True)
    score = cp.Variable()
    prefix = cp.sum(probabilities, axis=1)
    syndrome = cp.hstack(
        [sum(probabilities[z, z ^ s] for z in range(4)) for s in range(4)]
    )
    cap_weights = filled_effect_weights(float(terminal_weights.max()))
    constraints: list[cp.Constraint] = [
        cp.sum(probabilities) == 1.0,
        returned <= hellinger_return(probabilities),
        score <= weight * audit + (1.0 - weight) * returned,
    ]
    if audit_upper is not None:
        constraints.append(audit <= audit_upper)
    constraints.extend(
        prefix[order[index]] >= prefix[order[index + 1]]
        for index in range(3)
    )
    constraints.extend(
        (
            audit
            <= sum(cap_weights[index] * prefix[order[index]] for index in range(4)),
            audit <= terminal_weights @ syndrome,
        )
    )
    if projective_support_upper is not None:
        projective_support_lines = (
            *projective_support_lines,
            (weight, projective_support_upper),
        )
    for line_weight, line_upper in projective_support_lines:
        if not 0.0 < line_weight <= 1.0:
            raise ValueError("projective support-line weights must lie in (0,1]")
        for retained in range(4):
            for complement in range(4):
                if retained == complement:
                    continue
                bonus = (1.0 - terminal_weights[retained]) * syndrome[complement]
                bonus += sum(
                    terminal_weights[deleted] * syndrome[deleted]
                    for deleted in range(4)
                    if deleted not in {retained, complement}
                )
                constraints.append(
                    line_weight * audit
                    + (1.0 - line_weight) * returned
                    <= line_upper + line_weight * bonus
                )
    problem = cp.Problem(cp.Maximize(score), constraints)
    problem.solve(
        solver="CLARABEL",
        tol_gap_abs=2e-10,
        tol_gap_rel=2e-10,
        tol_feas=2e-10,
        max_iter=1000,
    )
    if problem.status not in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}:
        raise RuntimeError(f"prefix order {order} failed: {problem.status}")
    point = np.asarray(probabilities.value, dtype=float)
    prefix_value = point.sum(axis=1)
    syndrome_value = np.asarray(
        [sum(point[z, z ^ s] for z in range(4)) for s in range(4)]
    )
    return {
        "prefix_order": list(order),
        "score": float(problem.value),
        "audit_bound": float(audit.value),
        "return": float(returned.value),
        "prefix_priors": prefix_value.tolist(),
        "syndrome_priors": syndrome_value.tolist(),
        "probabilities": point.tolist(),
        "status": problem.status,
    }


def solve_weights(
    terminal_weights: np.ndarray,
    weight: float,
    audit_upper: float | None = None,
    projective_support_upper: float | None = None,
    projective_support_lines: tuple[tuple[float, float], ...] = (),
) -> dict[str, object]:
    rows = [
        solve_order(
            terminal_weights,
            order,
            weight,
            audit_upper,
            projective_support_upper,
            projective_support_lines,
        )
        for order in PERMUTATIONS
    ]
    rows.sort(key=lambda row: float(row["score"]), reverse=True)
    return {
        "weight": weight,
        "terminal_effect_weights": terminal_weights.tolist(),
        "audit_upper": audit_upper,
        "projective_support_upper": projective_support_upper,
        "projective_support_lines": [list(line) for line in projective_support_lines],
        "prefix_effect_cap_weights": filled_effect_weights(
            float(terminal_weights.max())
        ).tolist(),
        "bound": rows[0]["score"],
        "best_cell": rows[0],
        "cells": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lambda", dest="weight", type=float, default=0.6)
    parser.add_argument("--effect-weights", type=float, nargs=4, required=True)
    parser.add_argument("--audit-upper", type=float)
    parser.add_argument("--projective-upper", type=float)
    parser.add_argument(
        "--projective-line",
        type=float,
        nargs=2,
        action="append",
        default=[],
        metavar=("WEIGHT", "UPPER"),
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    terminal_weights = np.sort(np.asarray(args.effect_weights, dtype=float))[::-1]
    if (
        np.any(terminal_weights < 0.0)
        or np.any(terminal_weights > 1.0)
        or abs(float(terminal_weights.sum()) - 2.0) > 1e-9
    ):
        raise ValueError("terminal weights must lie in [0,1] and sum to two")
    payload = solve_weights(
        terminal_weights,
        args.weight,
        args.audit_upper,
        args.projective_upper,
        tuple(tuple(map(float, line)) for line in args.projective_line),
    )
    rendered = json.dumps(payload, indent=2) + "\n"
    print(rendered)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

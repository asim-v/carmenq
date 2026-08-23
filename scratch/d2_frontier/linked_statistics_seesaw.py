"""Column-generation relaxation for common qubit path statistics.

A physical twelve-effect POVM ``G[y,t]`` and four common prefix states are
retained.  Only selected probability columns are tied to the corresponding
Born products; all other path/terminal statistics remain free.  Selected
columns may be

* ``b_J`` for ``J=3*y+t``, representing ``G[y,t]``; or
* ``d_y_t``, representing the domination residual ``w_t Q[y]-G[y,t]``.

Adding columns monotonically tightens the model.  With every ``b_J`` linked,
it becomes the full pulled-statistics outer model.  Each alternating block is
convex; the seesaw is a lower diagnostic for the relaxation, not a global
upper certificate.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import cvxpy as cp
import numpy as np

from joint_effect_dimension_seesaw import random_point
from joint_statistics_helstrom_seesaw import terminal_kkt
from terminal_weight_upper import filled_effect_weights
from two_block_choi_seesaw import (
    IDENTITY,
    OUTCOMES,
    PATHS,
    canonical_three_effect_povm,
    hellinger_hypograph,
)


def solve(problem: cp.Problem) -> None:
    try:
        problem.solve(
            solver="CLARABEL",
            tol_gap_abs=2e-9,
            tol_gap_rel=2e-9,
            tol_feas=2e-9,
            max_iter=1000,
        )
    except cp.SolverError:
        problem.solve(
            solver="SCS",
            eps=2e-6,
            max_iters=300_000,
            acceleration_lookback=20,
        )
    if problem.status not in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}:
        raise RuntimeError(f"linked-statistics block failed: {problem.status}")


def parse_column(name: str) -> tuple[str, int, int]:
    parts = name.split("_")
    if len(parts) == 2 and parts[0] == "b":
        index = int(parts[1])
        if not 0 <= index < 12:
            raise ValueError(f"invalid effect column {name!r}")
        return "b", index // 3, index % 3
    if len(parts) == 3 and parts[0] == "d":
        y, t = int(parts[1]), int(parts[2])
        if not (0 <= y < 4 and 0 <= t < 3):
            raise ValueError(f"invalid residual column {name!r}")
        return "d", y, t
    raise ValueError(f"invalid linked column {name!r}")


def selected_statistic(
    statistics: cp.Expression,
    probability: cp.Expression,
    weights: np.ndarray,
    kind: str,
    z: int,
    y: int,
    t: int,
) -> cp.Expression:
    if kind == "b":
        return statistics[z, y, t]
    return weights[t] * probability[z, y] - statistics[z, y, t]


def selected_effect(
    joint: list[list[cp.Expression]] | np.ndarray,
    weights: np.ndarray,
    kind: str,
    y: int,
    t: int,
) -> cp.Expression:
    if kind == "b":
        return joint[y][t]
    return weights[t] * sum(joint[y][u] for u in range(3)) - joint[y][t]


def common_model(
    states: list[cp.Expression] | np.ndarray,
    joint: list[list[cp.Expression]] | np.ndarray,
    terminal_effects: np.ndarray,
    support_weight: float,
    prefix_order: tuple[int, int, int, int],
    linked: tuple[tuple[str, int, int], ...],
    constraints: list[cp.Constraint],
) -> tuple[cp.Expression, dict[str, cp.Expression]]:
    weights = np.trace(terminal_effects[:3], axis1=1, axis2=2).real
    statistics = cp.Variable((4, 4, 3), nonneg=True)
    probability = cp.sum(statistics, axis=2)
    constraints.append(cp.sum(statistics) == 1.0)
    constraints.extend(
        statistics[:, :, t] <= weights[t] * probability for t in range(3)
    )
    prefix = cp.sum(probability, axis=1)
    constraints.extend(prefix[z] == cp.real(cp.trace(states[z])) for z in OUTCOMES)
    constraints.extend(
        prefix[prefix_order[index]] >= prefix[prefix_order[index + 1]]
        for index in range(3)
    )

    for kind, y, t in linked:
        effect = selected_effect(joint, weights, kind, y, t)
        for z in OUTCOMES:
            constraints.append(
                selected_statistic(statistics, probability, weights, kind, z, y, t)
                == cp.real(cp.trace(effect @ states[z]))
            )

    terminal_statistics = [
        [
            sum(statistics[z, y, t] for z, y in PATHS if (z ^ y) == syndrome)
            for t in range(3)
        ]
        for syndrome in OUTCOMES
    ]
    audit = sum(terminal_statistics[s][s] for s in range(3))
    cap = filled_effect_weights(float(weights.max()))
    constraints.append(
        audit
        <= sum(cap[index] * prefix[prefix_order[index]] for index in OUTCOMES)
    )
    terminal_kkt(terminal_statistics, terminal_effects, audit, constraints)
    returned = hellinger_hypograph(
        [probability[z, y] for z, y in PATHS], constraints
    )
    return support_weight * audit + (1.0 - support_weight) * returned, {
        "statistics": statistics,
        "probability": probability,
        "audit": audit,
        "return": returned,
    }


def optimise_joint(
    states: np.ndarray,
    terminal: np.ndarray,
    support_weight: float,
    order: tuple[int, int, int, int],
    linked: tuple[tuple[str, int, int], ...],
) -> np.ndarray:
    weights = np.trace(terminal[:3], axis1=1, axis2=2).real
    joint = [[cp.Variable((2, 2), hermitian=True) for _ in range(3)] for _ in OUTCOMES]
    constraints: list[cp.Constraint] = [joint[y][t] >> 0 for y in OUTCOMES for t in range(3)]
    constraints.append(sum(joint[y][t] for y in OUTCOMES for t in range(3)) == IDENTITY)
    for y in OUTCOMES:
        coarse = sum(joint[y])
        constraints.extend(weights[t] * coarse - joint[y][t] >> 0 for t in range(3))
    objective, _ = common_model(states, joint, terminal, support_weight, order, linked, constraints)
    problem = cp.Problem(cp.Maximize(objective), constraints)
    solve(problem)
    return np.asarray([[np.asarray(joint[y][t].value) for t in range(3)] for y in OUTCOMES])


def optimise_states(
    joint: np.ndarray,
    terminal: np.ndarray,
    support_weight: float,
    order: tuple[int, int, int, int],
    linked: tuple[tuple[str, int, int], ...],
) -> np.ndarray:
    states = [cp.Variable((2, 2), hermitian=True) for _ in OUTCOMES]
    constraints: list[cp.Constraint] = [state >> 0 for state in states]
    constraints.append(sum(cp.real(cp.trace(state)) for state in states) == 1.0)
    objective, _ = common_model(states, joint, terminal, support_weight, order, linked, constraints)
    problem = cp.Problem(cp.Maximize(objective), constraints)
    solve(problem)
    return np.asarray([np.asarray(state.value) for state in states])


def evaluate(
    states: np.ndarray,
    joint: np.ndarray,
    terminal: np.ndarray,
    support_weight: float,
    order: tuple[int, int, int, int],
    linked: tuple[tuple[str, int, int], ...],
) -> dict[str, object]:
    constraints: list[cp.Constraint] = []
    objective, expressions = common_model(
        states, joint, terminal, support_weight, order, linked, constraints
    )
    problem = cp.Problem(cp.Maximize(objective), constraints)
    solve(problem)
    probability = np.asarray(expressions["probability"].value)
    statistics = np.asarray(expressions["statistics"].value)
    return {
        "score": float(problem.value),
        "audit": float(expressions["audit"].value),
        "return": float(np.sqrt(np.maximum(probability, 0.0)).sum() ** 2 / 16.0),
        "prefix_priors": probability.sum(axis=1).tolist(),
        "path_probabilities": probability.tolist(),
        "path_terminal_statistics": statistics.tolist(),
        "terminal_effect_weights": np.trace(terminal, axis1=1, axis2=2).real.tolist(),
    }


def sorted_initial(seed: int, weights: np.ndarray, order: tuple[int, ...]) -> tuple[np.ndarray, np.ndarray]:
    padded = np.zeros(4)
    padded[:3] = weights
    states, joint = random_point(seed, padded)
    ranked = states[np.argsort(np.trace(states, axis1=1, axis2=2).real)[::-1]]
    ordered = np.empty_like(states)
    for index, label in enumerate(order):
        ordered[label] = ranked[index]
    return ordered, np.asarray(joint[:, :3])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lambda", dest="weight", type=float, default=0.6)
    parser.add_argument("--fixed-three-povm-weights", type=float, nargs=3, required=True)
    parser.add_argument("--prefix-order", type=int, nargs=4, required=True)
    parser.add_argument("--linked-column", action="append", default=[])
    parser.add_argument("--rounds", type=int, default=8)
    parser.add_argument("--restarts", type=int, default=2)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    order = tuple(args.prefix_order)
    linked = tuple(parse_column(name) for name in dict.fromkeys(args.linked_column))
    weights = np.asarray(args.fixed_three_povm_weights)
    terminal = canonical_three_effect_povm(weights)
    rows = []
    arrays = []
    for restart in range(args.restarts):
        states, joint = sorted_initial(args.seed + restart, weights, order)
        best: dict[str, object] = {"score": -np.inf}
        best_arrays = (states.copy(), joint.copy())
        history = []
        for round_index in range(args.rounds):
            joint = optimise_joint(states, terminal, args.weight, order, linked)
            states = optimise_states(joint, terminal, args.weight, order, linked)
            point = evaluate(states, joint, terminal, args.weight, order, linked)
            point["round"] = round_index + 1
            history.append(point)
            if float(point["score"]) > float(best["score"]):
                best = point.copy()
                best_arrays = (states.copy(), joint.copy())
            print(restart, round_index + 1, point["score"], point["audit"], point["return"], flush=True)
        rows.append({"restart": restart, "linked_columns": args.linked_column, **best, "history": history})
        arrays.append(best_arrays)
    indices = np.argsort([-float(row["score"]) for row in rows])
    rows = [rows[int(index)] for index in indices]
    arrays = [arrays[int(index)] for index in indices]
    rendered = json.dumps(rows, indent=2) + "\n"
    print(rendered)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
        states, joint = arrays[0]
        np.savez_compressed(args.output.with_suffix(".npz"), states=states, joint=joint, terminal=terminal)


if __name__ == "__main__":
    main()

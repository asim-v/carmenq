"""Joint-effect outer seesaw with exact terminal Helstrom optimality.

This strengthens ``joint_effect_dimension_seesaw.py`` in two ways for a
fixed rank-one terminal POVM ``P_s = w_s Pi_s``:

* every pulled-back joint effect obeys ``G[y,s] <= w_s Q[y]``; and
* the syndrome priors admit terminal states for which the fixed POVM obeys
  the exact Helstrom primal-dual and complementary-slackness conditions.

The programme still drops the requirement that the joint effects and the
terminal states come from one common output instrument.  It is therefore an
outer diagnostic.  Each state/effect block is solved globally, while the
alternation is not a global certificate.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import cvxpy as cp
import numpy as np

from joint_effect_dimension_seesaw import (
    OUTCOMES,
    PATHS,
    audit_success,
    path_probabilities,
    random_point,
)
from two_block_choi_seesaw import (
    IDENTITY,
    canonical_three_effect_povm,
    hellinger_hypograph,
    solve_problem,
)


def helstrom_constraints(
    syndrome: list[cp.Expression],
    audit: cp.Expression,
    effects: np.ndarray,
    constraints: list[cp.Constraint],
) -> None:
    dual = cp.Variable((2, 2), hermitian=True)
    terminal = [cp.Variable((2, 2), hermitian=True) for _ in OUTCOMES]
    constraints.extend(item >> 0 for item in terminal)
    constraints.extend(dual - item >> 0 for item in terminal)
    constraints.extend(
        cp.real(cp.trace(terminal[s])) == syndrome[s] for s in OUTCOMES
    )
    constraints.append(cp.real(cp.trace(dual)) == audit)
    constraints.append(
        cp.sum(
            cp.hstack(
                [
                    cp.real(cp.trace(effects[s] @ terminal[s]))
                    for s in OUTCOMES
                ]
            )
        )
        == audit
    )
    constraints.extend(
        cp.real(cp.trace(effects[s] @ (dual - terminal[s]))) == 0.0
        for s in OUTCOMES
        if np.trace(effects[s]).real > 1e-9
    )


def evaluate(
    states: np.ndarray,
    joint: np.ndarray,
    effects: np.ndarray,
    weight: float,
) -> dict[str, object]:
    probabilities = path_probabilities(states, joint)
    syndromes = np.asarray(
        [sum(probabilities[z, z ^ s] for z in OUTCOMES) for s in OUTCOMES]
    )
    audit = audit_success(states, joint)
    returned = float(np.sqrt(np.maximum(probabilities, 0.0)).sum() ** 2 / 16.0)
    return {
        "score": weight * audit + (1.0 - weight) * returned,
        "audit": audit,
        "return_upper": returned,
        "normalisation": float(probabilities.sum()),
        "syndrome_priors": syndromes.tolist(),
        "terminal_effect_weights": np.trace(
            effects, axis1=1, axis2=2
        ).real.tolist(),
    }


def optimise_effects(
    states: np.ndarray,
    effects: np.ndarray,
    weight: float,
) -> np.ndarray:
    weights = np.trace(effects, axis1=1, axis2=2).real
    variables = [
        [cp.Variable((2, 2), hermitian=True) for _ in OUTCOMES]
        for _ in OUTCOMES
    ]
    constraints: list[cp.Constraint] = [
        variables[y][s] >> 0 for y, s in PATHS
    ]
    constraints.append(sum(variables[y][s] for y, s in PATHS) == IDENTITY)
    for y in OUTCOMES:
        coarse = sum(variables[y])
        constraints.extend(
            weights[s] * coarse - variables[y][s] >> 0 for s in OUTCOMES
        )

    probabilities = {
        (z, y): cp.real(cp.trace(sum(variables[y]) @ states[z]))
        for z, y in PATHS
    }
    returned = hellinger_hypograph(
        [probabilities[z, y] for z, y in PATHS], constraints
    )
    syndrome = [
        sum(probabilities[z, z ^ s] for z in OUTCOMES) for s in OUTCOMES
    ]
    audit = sum(
        cp.real(cp.trace(variables[y][z ^ y] @ states[z]))
        for z, y in PATHS
    )
    helstrom_constraints(syndrome, audit, effects, constraints)
    problem = cp.Problem(
        cp.Maximize(weight * audit + (1.0 - weight) * returned), constraints
    )
    solve_problem(problem)
    return np.stack(
        [
            np.stack([np.asarray(variables[y][s].value) for s in OUTCOMES])
            for y in OUTCOMES
        ]
    )


def optimise_states(
    joint: np.ndarray,
    effects: np.ndarray,
    weight: float,
) -> np.ndarray:
    variables = [cp.Variable((2, 2), hermitian=True) for _ in OUTCOMES]
    constraints: list[cp.Constraint] = [item >> 0 for item in variables]
    constraints.append(
        cp.sum(cp.hstack([cp.real(cp.trace(item)) for item in variables])) == 1.0
    )
    coarse = joint.sum(axis=1)
    probabilities = {
        (z, y): cp.real(cp.trace(coarse[y] @ variables[z]))
        for z, y in PATHS
    }
    returned = hellinger_hypograph(
        [probabilities[z, y] for z, y in PATHS], constraints
    )
    syndrome = [
        sum(probabilities[z, z ^ s] for z in OUTCOMES) for s in OUTCOMES
    ]
    audit = sum(
        cp.real(cp.trace(joint[y, z ^ y] @ variables[z])) for z, y in PATHS
    )
    helstrom_constraints(syndrome, audit, effects, constraints)
    problem = cp.Problem(
        cp.Maximize(weight * audit + (1.0 - weight) * returned), constraints
    )
    solve_problem(problem)
    return np.stack([np.asarray(item.value) for item in variables])


def optimise_seed(
    seed: int,
    effects: np.ndarray,
    weight: float,
    rounds: int,
) -> tuple[dict[str, object], tuple[np.ndarray, np.ndarray]]:
    weights = np.trace(effects, axis1=1, axis2=2).real
    states, joint = random_point(seed + 1_700_003, weights)
    best: dict[str, object] = {"score": -np.inf}
    best_arrays = (states.copy(), joint.copy())
    history: list[dict[str, object]] = []
    for round_index in range(rounds):
        joint = optimise_effects(states, effects, weight)
        states = optimise_states(joint, effects, weight)
        point = evaluate(states, joint, effects, weight)
        point["round"] = round_index + 1
        history.append(point)
        if float(point["score"]) > float(best["score"]):
            best = point.copy()
            best_arrays = (states.copy(), joint.copy())
        print(seed, point, flush=True)
    return {"seed": seed, **best, "history": history}, best_arrays


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lambda", dest="weight", type=float, default=0.6)
    parser.add_argument("--fixed-three-povm-weights", type=float, nargs=3, required=True)
    parser.add_argument("--rounds", type=int, default=8)
    parser.add_argument("--restarts", type=int, default=2)
    parser.add_argument("--seed-offset", type=int, default=0)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    weights = np.asarray(args.fixed_three_povm_weights, dtype=float)
    if (
        np.any(weights <= 0.0)
        or np.any(weights > 1.0)
        or abs(float(weights.sum()) - 2.0) > 1e-9
    ):
        raise ValueError("three-effect weights must lie in (0,1] and sum to two")
    effects = canonical_three_effect_povm(weights)
    results = [
        optimise_seed(seed, effects, args.weight, args.rounds)
        for seed in range(args.seed_offset, args.seed_offset + args.restarts)
    ]
    results.sort(key=lambda item: float(item[0]["score"]), reverse=True)
    rows = [item[0] for item in results]
    rendered = json.dumps(rows, indent=2) + "\n"
    print(rendered)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
        states, joint = results[0][1]
        np.savez_compressed(
            args.output.with_suffix(".npz"),
            states=states,
            joint=joint,
            effects=effects,
        )


if __name__ == "__main__":
    main()

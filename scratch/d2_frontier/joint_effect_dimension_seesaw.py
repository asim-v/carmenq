"""Joint-effect plus terminal-dimension outer relaxation.

The effective AUDIT POVM on the prefix qubit and the coarse suffix POVM are
marginals of one joint 16-effect qubit POVM ``G[y,s]``.  This script retains
that compatibility, the exact path Hellinger term, and the independent
terminal-qubit dimension bound.  It drops the requirement that every
``G[y,s]`` arise by pulling one common terminal POVM through a common output
instrument, so its optimum is an upper relaxation of the Choi programme.

For a fixed top-two terminal-syndrome pair, each state/effect update is a
globally solved conic programme.  The outer seesaw is a diagnostic rather
than a joint global certificate.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path

import cvxpy as cp
import numpy as np

from two_block_choi_seesaw import (
    IDENTITY,
    OUTCOMES,
    hellinger_hypograph,
    pulled_effect,
    solve_problem,
)


PAIRS = tuple(itertools.combinations(OUTCOMES, 2))
PATHS = tuple((z, y) for z in OUTCOMES for y in OUTCOMES)


def inverse_square_root(matrix: np.ndarray) -> np.ndarray:
    values, vectors = np.linalg.eigh(matrix)
    return (vectors * np.maximum(values.real, 1e-14) ** -0.5) @ vectors.conj().T


def random_point(
    seed: int, effect_weights: np.ndarray | None = None
) -> tuple[np.ndarray, np.ndarray]:
    generator = np.random.default_rng(seed)
    roots = (
        generator.standard_normal((4, 2, 2))
        + 1j * generator.standard_normal((4, 2, 2))
    ) / math.sqrt(2.0)
    states = np.einsum("zai,zaj->zij", roots.conj(), roots)
    states /= np.trace(states, axis1=1, axis2=2).real.sum()

    effect_roots = (
        generator.standard_normal((4, 4, 2, 2))
        + 1j * generator.standard_normal((4, 4, 2, 2))
    ) / math.sqrt(2.0)
    effects = np.einsum("ysai,ysaj->ysij", effect_roots.conj(), effect_roots)
    normalizer = inverse_square_root(effects.sum(axis=(0, 1)))
    effects = np.einsum("ai,ysij,jb->ysab", normalizer, effects, normalizer)
    if effect_weights is not None:
        coarse = effects.sum(axis=1)
        effects = np.einsum("s,yij->ysij", effect_weights / 2.0, coarse)
    return states, effects


def checkpoint_point(path: Path) -> tuple[np.ndarray, np.ndarray]:
    arrays = np.load(path)
    states = np.asarray(arrays["states"])
    choi = np.asarray(arrays["choi"])
    terminal = np.asarray(arrays["effects"])
    joint = np.stack(
        [
            np.stack([pulled_effect(choi[y], terminal[s]) for s in OUTCOMES])
            for y in OUTCOMES
        ]
    )
    return states, joint


def path_probabilities(states: np.ndarray, joint: np.ndarray) -> np.ndarray:
    coarse = joint.sum(axis=1)
    return np.einsum("yij,zji->zy", coarse, states).real


def audit_success(states: np.ndarray, joint: np.ndarray) -> float:
    return float(
        sum(
            np.trace(joint[y, z ^ y] @ states[z]).real for z, y in PATHS
        )
    )


def syndrome_priors(probabilities: np.ndarray) -> np.ndarray:
    return np.asarray(
        [sum(probabilities[z, z ^ s] for z in OUTCOMES) for s in OUTCOMES]
    )


def evaluate(
    states: np.ndarray,
    joint: np.ndarray,
    pair: tuple[int, int],
    weight: float,
) -> dict[str, float | list[float]]:
    probabilities = path_probabilities(states, joint)
    syndromes = syndrome_priors(probabilities)
    direct_audit = audit_success(states, joint)
    dimension_audit = float(syndromes[list(pair)].sum())
    audit = min(direct_audit, dimension_audit)
    returned = float(np.sqrt(np.maximum(probabilities, 0.0)).sum() ** 2 / 16.0)
    return {
        "score": weight * audit + (1.0 - weight) * returned,
        "audit_bound": audit,
        "direct_audit": direct_audit,
        "terminal_dimension_audit": dimension_audit,
        "return_upper": returned,
        "normalisation": float(probabilities.sum()),
        "syndrome_priors": syndromes.tolist(),
    }


def top_pair_constraints(
    values: list[cp.Expression], pair: tuple[int, int]
) -> list[cp.Constraint]:
    outside = tuple(index for index in OUTCOMES if index not in pair)
    return [values[i] >= values[j] for i in pair for j in outside]


def optimise_effects(
    states: np.ndarray,
    pair: tuple[int, int],
    weight: float,
    effect_weights: np.ndarray | None,
) -> np.ndarray:
    variables = [
        [cp.Variable((2, 2), hermitian=True) for _ in OUTCOMES]
        for _ in OUTCOMES
    ]
    constraints: list[cp.Constraint] = [
        variables[y][s] >> 0 for y in OUTCOMES for s in OUTCOMES
    ]
    constraints.append(sum(variables[y][s] for y, s in PATHS) == IDENTITY)
    if effect_weights is not None:
        for y in OUTCOMES:
            coarse_y = sum(variables[y])
            constraints.extend(
                effect_weights[s] * coarse_y - variables[y][s] >> 0
                for s in OUTCOMES
            )

    probabilities: dict[tuple[int, int], cp.Expression] = {}
    for z, y in PATHS:
        probabilities[z, y] = cp.real(
            cp.trace(sum(variables[y]) @ states[z])
        )
    flat = [probabilities[z, y] for z, y in PATHS]
    returned = hellinger_hypograph(flat, constraints)
    syndromes = [
        sum(probabilities[z, z ^ s] for z in OUTCOMES) for s in OUTCOMES
    ]
    constraints += top_pair_constraints(syndromes, pair)
    direct_audit = sum(
        cp.real(cp.trace(variables[y][z ^ y] @ states[z])) for z, y in PATHS
    )
    audit = cp.Variable()
    constraints += [
        audit <= direct_audit,
        audit <= sum(syndromes[s] for s in pair),
    ]
    problem = cp.Problem(
        cp.Maximize(weight * audit + (1.0 - weight) * returned), constraints
    )
    solve_problem(problem)
    return np.stack(
        [np.stack([np.asarray(variables[y][s].value) for s in OUTCOMES]) for y in OUTCOMES]
    )


def optimise_states(
    joint: np.ndarray, pair: tuple[int, int], weight: float
) -> np.ndarray:
    variables = [cp.Variable((2, 2), hermitian=True) for _ in OUTCOMES]
    constraints: list[cp.Constraint] = [item >> 0 for item in variables]
    constraints.append(
        cp.sum(cp.hstack([cp.real(cp.trace(item)) for item in variables])) == 1.0
    )
    coarse = joint.sum(axis=1)
    probabilities: dict[tuple[int, int], cp.Expression] = {}
    for z, y in PATHS:
        probabilities[z, y] = cp.real(cp.trace(coarse[y] @ variables[z]))
    flat = [probabilities[z, y] for z, y in PATHS]
    returned = hellinger_hypograph(flat, constraints)
    syndromes = [
        sum(probabilities[z, z ^ s] for z in OUTCOMES) for s in OUTCOMES
    ]
    constraints += top_pair_constraints(syndromes, pair)
    direct_audit = sum(
        cp.real(cp.trace(joint[y, z ^ y] @ variables[z])) for z, y in PATHS
    )
    audit = cp.Variable()
    constraints += [
        audit <= direct_audit,
        audit <= sum(syndromes[s] for s in pair),
    ]
    problem = cp.Problem(
        cp.Maximize(weight * audit + (1.0 - weight) * returned), constraints
    )
    solve_problem(problem)
    return np.stack([np.asarray(item.value) for item in variables])


def optimise_seed(
    seed: int,
    pair: tuple[int, int],
    weight: float,
    rounds: int,
    checkpoint: Path | None,
    effect_weights: np.ndarray | None,
) -> dict[str, object]:
    if checkpoint is not None and seed == 0:
        states, joint = checkpoint_point(checkpoint)
    else:
        states, joint = random_point(seed + 1_400_003, effect_weights)
    best = evaluate(states, joint, pair, weight)
    history: list[dict[str, object]] = []
    for round_index in range(rounds):
        joint = optimise_effects(states, pair, weight, effect_weights)
        states = optimise_states(joint, pair, weight)
        point = evaluate(states, joint, pair, weight)
        point["round"] = round_index + 1
        history.append(point)
        if float(point["score"]) > float(best["score"]):
            best = point.copy()
        print(seed, pair, point, flush=True)
    return {
        "seed": seed,
        "top_pair": list(pair),
        "effect_weights": (
            None if effect_weights is None else effect_weights.tolist()
        ),
        **best,
        "history": history,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lambda", dest="weight", type=float, default=0.6)
    parser.add_argument("--rounds", type=int, default=8)
    parser.add_argument("--restarts", type=int, default=2)
    parser.add_argument("--seed-offset", type=int, default=0)
    parser.add_argument("--pair", type=int, nargs=2)
    parser.add_argument("--effect-weights", type=float, nargs=4)
    parser.add_argument("--checkpoint", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    pairs = (tuple(args.pair),) if args.pair is not None else PAIRS
    if any(pair not in PAIRS for pair in pairs):
        raise ValueError("pair must contain two distinct syndrome labels")
    effect_weights = (
        None
        if args.effect_weights is None
        else np.asarray(args.effect_weights, dtype=float)
    )
    if effect_weights is not None and (
        np.any(effect_weights < 0.0)
        or np.any(effect_weights > 1.0)
        or abs(float(effect_weights.sum()) - 2.0) > 1e-9
    ):
        raise ValueError("effect weights must lie in [0,1] and sum to two")
    if effect_weights is not None and args.checkpoint is not None:
        raise ValueError("a checkpoint is not guaranteed feasible for fixed weights")

    rows = [
        optimise_seed(
            seed,
            pair,
            args.weight,
            args.rounds,
            args.checkpoint,
            effect_weights,
        )
        for pair in pairs
        for seed in range(args.seed_offset, args.seed_offset + args.restarts)
    ]
    rows.sort(key=lambda row: float(row["score"]), reverse=True)
    rendered = json.dumps(rows, indent=2) + "\n"
    print(rendered)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

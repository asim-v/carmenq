"""Pulled-statistics outer model with exact terminal Helstrom optimality.

For a fixed nondegenerate three-effect terminal POVM P_t, define the pulled
effects G[y,t] = Phi_y^*(P_t).  Unlike the earlier joint-effect outer model,
this programme retains *all* three terminal-measurement probabilities on every
path.  Their syndrome sums reconstruct the in-plane Bloch coordinates of the
terminal states, so the states entering the Helstrom KKT system cannot be
chosen independently of AUDIT.

The model is still an outer relaxation: it asks only that each family of
pulled effects is a qubit POVM refinement obeying G[y,t] <= w_t Q[y], not that
one completely positive output map realizes the full operator system.  Each
alternating block is a global conic problem; the seesaw is a diagnostic rather
than a global certificate.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import cvxpy as cp
import numpy as np

from joint_effect_dimension_seesaw import OUTCOMES, PATHS, random_point
from two_block_choi_seesaw import (
    IDENTITY,
    canonical_three_effect_povm,
    hellinger_hypograph,
    solve_problem,
)


SIGMA_X = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
SIGMA_Y = np.array([[0.0, -1j], [1j, 0.0]], dtype=complex)
SIGMA_Z = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)


def reconstruction_matrix(effects: np.ndarray) -> np.ndarray:
    """Map [trace, x, y] to the three POVM probabilities."""
    rows = []
    for effect in effects[:3]:
        rows.append(
            [
                0.5 * float(np.trace(effect).real),
                0.5 * float(np.trace(effect @ SIGMA_X).real),
                0.5 * float(np.trace(effect @ SIGMA_Y).real),
            ]
        )
    matrix = np.asarray(rows, dtype=float)
    if abs(float(np.linalg.det(matrix))) < 1e-10:
        raise ValueError("the terminal three-effect POVM must be nondegenerate")
    return np.linalg.inv(matrix)


def terminal_kkt(
    terminal_statistics: list[list[cp.Expression]],
    effects: np.ndarray,
    audit: cp.Expression,
    constraints: list[cp.Constraint],
) -> list[cp.Expression]:
    inverse = reconstruction_matrix(effects)
    normal_coordinates = cp.Variable(4)
    terminal: list[cp.Expression] = []
    for s in OUTCOMES:
        reconstructed = inverse @ cp.hstack(terminal_statistics[s])
        trace = reconstructed[0]
        x_coordinate = reconstructed[1]
        y_coordinate = reconstructed[2]
        state = 0.5 * (
            trace * IDENTITY
            + x_coordinate * SIGMA_X
            + y_coordinate * SIGMA_Y
            + normal_coordinates[s] * SIGMA_Z
        )
        constraints.append(state >> 0)
        terminal.append(state)

    dual = cp.Variable((2, 2), hermitian=True)
    constraints.extend(dual - terminal[s] >> 0 for s in OUTCOMES)
    constraints.append(cp.real(cp.trace(dual)) == audit)
    return terminal


def model_expressions(
    states: list[cp.Expression] | np.ndarray,
    joint: list[list[cp.Expression]] | np.ndarray,
) -> tuple[
    dict[tuple[int, int], cp.Expression],
    dict[tuple[int, int, int], cp.Expression],
]:
    path_probabilities: dict[tuple[int, int], cp.Expression] = {}
    statistics: dict[tuple[int, int, int], cp.Expression] = {}
    for z, y in PATHS:
        path_probabilities[z, y] = cp.real(
            cp.trace(sum(joint[y][t] for t in range(3)) @ states[z])
        )
        for t in range(3):
            statistics[z, y, t] = cp.real(cp.trace(joint[y][t] @ states[z]))
    return path_probabilities, statistics


def add_common_constraints(
    states: list[cp.Expression] | np.ndarray,
    joint: list[list[cp.Expression]] | np.ndarray,
    effects: np.ndarray,
    weight: float,
    constraints: list[cp.Constraint],
) -> cp.Expression:
    probabilities, statistics = model_expressions(states, joint)
    returned = hellinger_hypograph(
        [probabilities[z, y] for z, y in PATHS], constraints
    )
    terminal_statistics = [
        [
            sum(
                statistics[z, y, t]
                for z, y in PATHS
                if (z ^ y) == s
            )
            for t in range(3)
        ]
        for s in OUTCOMES
    ]
    audit = sum(terminal_statistics[s][s] for s in range(3))
    terminal_kkt(terminal_statistics, effects, audit, constraints)
    return weight * audit + (1.0 - weight) * returned


def optimise_joint(states: np.ndarray, effects: np.ndarray, weight: float) -> np.ndarray:
    weights = np.trace(effects[:3], axis1=1, axis2=2).real
    variables = [
        [cp.Variable((2, 2), hermitian=True) for _ in range(3)]
        for _ in OUTCOMES
    ]
    constraints: list[cp.Constraint] = [
        variables[y][t] >> 0 for y in OUTCOMES for t in range(3)
    ]
    constraints.append(
        sum(variables[y][t] for y in OUTCOMES for t in range(3)) == IDENTITY
    )
    for y in OUTCOMES:
        coarse = sum(variables[y])
        constraints.extend(
            weights[t] * coarse - variables[y][t] >> 0 for t in range(3)
        )
    objective = add_common_constraints(states, variables, effects, weight, constraints)
    problem = cp.Problem(cp.Maximize(objective), constraints)
    solve_problem(problem)
    return np.asarray(
        [[np.asarray(variables[y][t].value) for t in range(3)] for y in OUTCOMES]
    )


def optimise_states(joint: np.ndarray, effects: np.ndarray, weight: float) -> np.ndarray:
    variables = [cp.Variable((2, 2), hermitian=True) for _ in OUTCOMES]
    constraints: list[cp.Constraint] = [item >> 0 for item in variables]
    constraints.append(sum(cp.real(cp.trace(item)) for item in variables) == 1.0)
    objective = add_common_constraints(variables, joint, effects, weight, constraints)
    problem = cp.Problem(cp.Maximize(objective), constraints)
    solve_problem(problem)
    return np.asarray([np.asarray(item.value) for item in variables])


def evaluate(
    states: np.ndarray,
    joint: np.ndarray,
    effects: np.ndarray,
    weight: float,
) -> dict[str, object]:
    probabilities = np.asarray(
        [
            [float(np.trace(joint[y].sum(axis=0) @ states[z]).real) for y in OUTCOMES]
            for z in OUTCOMES
        ]
    )
    statistics = np.asarray(
        [
            [
                [float(np.trace(joint[y, t] @ states[z]).real) for t in range(3)]
                for y in OUTCOMES
            ]
            for z in OUTCOMES
        ]
    )
    terminal_statistics = np.asarray(
        [
            [
                sum(statistics[z, y, t] for z, y in PATHS if (z ^ y) == s)
                for t in range(3)
            ]
            for s in OUTCOMES
        ]
    )
    audit = float(sum(terminal_statistics[s, s] for s in range(3)))
    returned = float(np.sqrt(np.maximum(probabilities, 0.0)).sum() ** 2 / 16.0)
    return {
        "score": weight * audit + (1.0 - weight) * returned,
        "audit": audit,
        "return_upper": returned,
        "normalisation": float(probabilities.sum()),
        "prefix_priors": probabilities.sum(axis=1).tolist(),
        "syndrome_priors": [
            float(sum(probabilities[z, z ^ s] for z in OUTCOMES)) for s in OUTCOMES
        ],
        "terminal_statistics": terminal_statistics.tolist(),
        "terminal_effect_weights": np.trace(
            effects, axis1=1, axis2=2
        ).real.tolist(),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lambda", dest="weight", type=float, default=0.6)
    parser.add_argument("--fixed-three-povm-weights", type=float, nargs=3, required=True)
    parser.add_argument("--rounds", type=int, default=10)
    parser.add_argument("--restarts", type=int, default=4)
    parser.add_argument("--seed-offset", type=int, default=0)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    weights = np.asarray(args.fixed_three_povm_weights, dtype=float)
    effects = canonical_three_effect_povm(weights)
    rows = []
    arrays = []
    for seed in range(args.seed_offset, args.seed_offset + args.restarts):
        states, initial_joint = random_point(seed + 4_000_003, weights)
        joint = np.asarray(initial_joint[:, :3])
        history = []
        best: dict[str, object] = {"score": -np.inf}
        best_arrays = (states.copy(), joint.copy())
        for round_index in range(args.rounds):
            joint = optimise_joint(states, effects, args.weight)
            states = optimise_states(joint, effects, args.weight)
            point = evaluate(states, joint, effects, args.weight)
            point["round"] = round_index + 1
            history.append(point)
            if float(point["score"]) > float(best["score"]):
                best = point.copy()
                best_arrays = (states.copy(), joint.copy())
            print(seed, point, flush=True)
        rows.append({"seed": seed, **best, "history": history})
        arrays.append(best_arrays)
    order = np.argsort([-float(row["score"]) for row in rows])
    rows = [rows[int(index)] for index in order]
    arrays = [arrays[int(index)] for index in order]
    rendered = json.dumps(rows, indent=2) + "\n"
    print(rendered)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
        states, joint = arrays[0]
        np.savez_compressed(
            args.output.with_suffix(".npz"), states=states, joint=joint, effects=effects
        )


if __name__ == "__main__":
    main()

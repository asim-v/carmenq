"""Global first-moment SDP bound for the joint-effect/Helstrom model.

For a fixed rank-one terminal qubit POVM, write every prefix state and every
nonzero joint effect in Bloch coordinates.  All physical probabilities are
quadratic products of those real coordinates.  This script replaces the
rank-one moment matrix ``[1;x][1;x]^T`` by one global positive-semidefinite
matrix, then adds:

* RLT products of every POVM-completeness equality with every variable;
* lifted Lorentz constraints for states, effects, and effect domination;
* variable-bound and state/effect McCormick inequalities;
* all terminal-POVM statistics and the exact terminal Helstrom dual; and
* the exact conic Hellinger hypograph.

Every joint-effect strategy is feasible in this SDP, so its optimum is a
global upper bound for the fixed terminal POVM.  This is an outer relaxation
of the common-instrument Choi programme, not an equality claim.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import cvxpy as cp
import numpy as np

from terminal_weight_upper import filled_effect_weights
from two_block_choi_seesaw import (
    IDENTITY,
    OUTCOMES,
    PATHS,
    canonical_three_effect_povm,
    hellinger_hypograph,
)


class Registry:
    def __init__(self) -> None:
        self.names: list[str] = []
        self.lower: list[float] = []
        self.upper: list[float] = []

    def add(self, name: str, lower: float, upper: float) -> int:
        index = len(self.names)
        self.names.append(name)
        self.lower.append(lower)
        self.upper.append(upper)
        return index


def moment_product(
    second: cp.Expression,
    first_terms: list[tuple[float, int]],
    second_terms: list[tuple[float, int]],
) -> cp.Expression:
    return sum(
        first_coefficient
        * second_coefficient
        * second[first_index, second_index]
        for first_coefficient, first_index in first_terms
        for second_coefficient, second_index in second_terms
    )


def add_mccormick(
    first: cp.Expression,
    second_value: cp.Expression,
    product: cp.Expression,
    first_bounds: tuple[float, float],
    second_bounds: tuple[float, float],
    constraints: list[cp.Constraint],
) -> None:
    first_lower, first_upper = first_bounds
    second_lower, second_upper = second_bounds
    constraints.extend(
        (
            product
            >= first_lower * second_value
            + second_lower * first
            - first_lower * second_lower,
            product
            >= first_upper * second_value
            + second_upper * first
            - first_upper * second_upper,
            product
            <= first_upper * second_value
            + second_lower * first
            - first_upper * second_lower,
            product
            <= first_lower * second_value
            + second_upper * first
            - first_lower * second_upper,
        )
    )


def solve_povm(
    effects: np.ndarray,
    weight: float,
    prefix_order: tuple[int, int, int, int],
) -> dict[str, object]:
    traces = np.trace(effects, axis1=1, axis2=2).real
    active = tuple(int(s) for s in OUTCOMES if traces[s] > 1e-9)
    directions = np.zeros((4, 3), dtype=float)
    paulis = (
        np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex),
        np.array([[0.0, -1j], [1j, 0.0]], dtype=complex),
        np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex),
    )
    for s in active:
        projector = effects[s] / traces[s]
        directions[s] = [
            float(np.trace(projector @ pauli).real) for pauli in paulis
        ]

    registry = Registry()
    state_index: dict[tuple[int, int], int] = {}
    for z in OUTCOMES:
        state_index[z, 0] = registry.add(f"a_{z}", 0.0, 1.0)
        for axis in range(1, 4):
            state_index[z, axis] = registry.add(f"r_{z}_{axis}", -1.0, 1.0)

    effect_index: dict[tuple[int, int, int], int] = {}
    for y in OUTCOMES:
        for s in active:
            upper = 2.0 * float(traces[s])
            effect_index[y, s, 0] = registry.add(
                f"g_{y}_{s}_0", 0.0, upper
            )
            for axis in range(1, 4):
                effect_index[y, s, axis] = registry.add(
                    f"g_{y}_{s}_{axis}", -upper, upper
                )

    size = len(registry.names)
    moment = cp.Variable((size + 1, size + 1), symmetric=True)
    first = moment[0, 1:]
    second = moment[1:, 1:]
    constraints: list[cp.Constraint] = [moment >> 0, moment[0, 0] == 1.0]
    constraints.extend(
        (
            first[index] >= registry.lower[index],
            first[index] <= registry.upper[index],
            second[index, index]
            <= (registry.lower[index] + registry.upper[index]) * first[index]
            - registry.lower[index] * registry.upper[index],
        )
        for index in range(size)
    )
    # Flatten the tuples introduced above.
    constraints = [item for group in constraints for item in (group if isinstance(group, tuple) else (group,))]

    equalities: list[tuple[list[tuple[float, int]], float]] = []
    equalities.append(
        ([(1.0, state_index[z, 0]) for z in OUTCOMES], 1.0)
    )
    equalities.append(
        (
            [
                (1.0, effect_index[y, s, 0])
                for y in OUTCOMES
                for s in active
            ],
            2.0,
        )
    )
    for axis in range(1, 4):
        equalities.append(
            (
                [
                    (1.0, effect_index[y, s, axis])
                    for y in OUTCOMES
                    for s in active
                ],
                0.0,
            )
        )
    for terms, constant in equalities:
        constraints.append(
            sum(coefficient * first[index] for coefficient, index in terms)
            == constant
        )
        for column in range(size):
            constraints.append(
                sum(
                    coefficient * second[index, column]
                    for coefficient, index in terms
                )
                == constant * first[column]
            )

    for z in OUTCOMES:
        scalar = state_index[z, 0]
        constraints.append(
            second[scalar, scalar]
            >= sum(
                second[state_index[z, axis], state_index[z, axis]]
                for axis in range(1, 4)
            )
        )
    for y in OUTCOMES:
        for s in active:
            scalar = effect_index[y, s, 0]
            constraints.append(
                second[scalar, scalar]
                >= sum(
                    second[
                        effect_index[y, s, axis], effect_index[y, s, axis]
                    ]
                    for axis in range(1, 4)
                )
            )

    def coarse_terms(y: int, axis: int) -> list[tuple[float, int]]:
        return [(1.0, effect_index[y, s, axis]) for s in active]

    # Lift positivity of w_s Q_y - G_ys.
    for y in OUTCOMES:
        for s in active:
            residual = {
                axis: [
                    *( (traces[s], index) for _, index in coarse_terms(y, axis) ),
                    (-1.0, effect_index[y, s, axis]),
                ]
                for axis in range(4)
            }
            residual_first = {
                axis: sum(coefficient * first[index] for coefficient, index in residual[axis])
                for axis in range(4)
            }
            constraints.append(residual_first[0] >= 0.0)
            constraints.append(
                moment_product(second, residual[0], residual[0])
                >= sum(
                    moment_product(second, residual[axis], residual[axis])
                    for axis in range(1, 4)
                )
            )

    # McCormick consistency for every state/effect coordinate pair.
    for z in OUTCOMES:
        for state_axis in range(4):
            first_index = state_index[z, state_axis]
            for y in OUTCOMES:
                for s in active:
                    for effect_axis in range(4):
                        second_index = effect_index[y, s, effect_axis]
                        add_mccormick(
                            first[first_index],
                            first[second_index],
                            second[first_index, second_index],
                            (
                                registry.lower[first_index],
                                registry.upper[first_index],
                            ),
                            (
                                registry.lower[second_index],
                                registry.upper[second_index],
                            ),
                            constraints,
                        )

    statistics: dict[tuple[int, int, int], cp.Expression] = {}
    probabilities: dict[tuple[int, int], cp.Expression] = {}
    correct: dict[tuple[int, int], cp.Expression] = {}
    for z, y in PATHS:
        for terminal_label in active:
            statistics[z, y, terminal_label] = 0.5 * sum(
                second[
                    state_index[z, axis], effect_index[y, terminal_label, axis]
                ]
                for axis in range(4)
            )
        probabilities[z, y] = sum(
            statistics[z, y, terminal_label] for terminal_label in active
        )
        syndrome = z ^ y
        if syndrome in active:
            correct[z, y] = statistics[z, y, syndrome]
        else:
            correct[z, y] = 0.0
        constraints.append(probabilities[z, y] >= 0.0)
        for terminal_label in active:
            constraints.extend(
                (
                    statistics[z, y, terminal_label] >= 0.0,
                    statistics[z, y, terminal_label]
                    <= traces[terminal_label] * probabilities[z, y],
                )
            )
    constraints.append(sum(probabilities.values()) == 1.0)
    audit = sum(correct.values())

    prefix = [sum(probabilities[z, y] for y in OUTCOMES) for z in OUTCOMES]
    constraints.extend(
        prefix[prefix_order[index]] >= prefix[prefix_order[index + 1]]
        for index in range(3)
    )
    cap = filled_effect_weights(float(traces.max()))
    constraints.append(
        audit
        <= sum(cap[index] * prefix[prefix_order[index]] for index in OUTCOMES)
    )
    syndrome_priors_from_paths = [
        sum(probabilities[z, z ^ s] for z in OUTCOMES) for s in OUTCOMES
    ]

    terminal_statistics = [
        [
            sum(
                statistics[z, y, terminal_label]
                for z, y in PATHS
                if (z ^ y) == syndrome
            )
            for terminal_label in active
        ]
        for syndrome in OUTCOMES
    ]
    # The canonical ternary POVM spans I, sigma_x, and sigma_y.  Hence its
    # three probabilities reconstruct the trace and both in-plane Bloch
    # coordinates of every terminal state.  The normal coordinate remains a
    # free physical variable, but the terminal states can no longer be chosen
    # independently of the pulled-effect products.
    reconstruction = np.asarray(
        [
            [
                0.5 * traces[s],
                0.5 * traces[s] * directions[s, 0],
                0.5 * traces[s] * directions[s, 1],
            ]
            for s in active
        ],
        dtype=float,
    )
    reconstruction_inverse = np.linalg.inv(reconstruction)
    reconstructed = [
        reconstruction_inverse @ cp.hstack(terminal_statistics[s])
        for s in OUTCOMES
    ]
    syndrome_priors = [reconstructed[s][0] for s in OUTCOMES]
    constraints.extend(
        syndrome_priors[s] == syndrome_priors_from_paths[s] for s in OUTCOMES
    )

    dual_trace = cp.Variable(nonneg=True)
    dual_vector = cp.Variable(3)
    terminal_normal = cp.Variable(4)
    terminal_vector = [
        cp.hstack(
            [reconstructed[s][1], reconstructed[s][2], terminal_normal[s]]
        )
        for s in OUTCOMES
    ]
    constraints.append(cp.SOC(dual_trace, dual_vector))
    for s in OUTCOMES:
        constraints.extend(
            (
                cp.SOC(syndrome_priors[s], terminal_vector[s]),
                cp.SOC(
                    dual_trace - syndrome_priors[s],
                    dual_vector - terminal_vector[s],
                ),
            )
        )
    constraints.append(audit == dual_trace)
    constraints.append(
        audit <= sum(traces[s] * syndrome_priors[s] for s in OUTCOMES)
    )

    returned = hellinger_hypograph(
        [probabilities[z, y] for z, y in PATHS], constraints
    )
    problem = cp.Problem(
        cp.Maximize(weight * audit + (1.0 - weight) * returned), constraints
    )
    try:
        problem.solve(
            solver="CLARABEL",
            tol_gap_abs=2e-8,
            tol_gap_rel=2e-8,
            tol_feas=2e-8,
            max_iter=1000,
        )
    except cp.error.SolverError:
        problem.solve(
            solver="SCS",
            eps=2e-5,
            max_iters=100_000,
            acceleration_lookback=20,
        )
    if problem.status not in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}:
        raise RuntimeError(f"moment relaxation failed: {problem.status}")

    point = np.asarray(
        [[float(probabilities[z, y].value) for y in OUTCOMES] for z in OUTCOMES]
    )
    moment_value = np.asarray(moment.value)
    return {
        "weight": weight,
        "terminal_effect_weights": traces.tolist(),
        "prefix_order": list(prefix_order),
        "moment_size": size + 1,
        "bound": float(problem.value),
        "audit": float(audit.value),
        "return": float(np.sqrt(np.maximum(point, 0.0)).sum() ** 2 / 16.0),
        "normalisation": float(point.sum()),
        "prefix_priors": point.sum(axis=1).tolist(),
        "syndrome_priors": [
            float(sum(point[z, z ^ s] for z in OUTCOMES)) for s in OUTCOMES
        ],
        "terminal_statistics": [
            [float(item.value) for item in terminal_statistics[s]]
            for s in OUTCOMES
        ],
        "terminal_bloch_vectors": [
            np.asarray(terminal_vector[s].value, dtype=float).tolist()
            for s in OUTCOMES
        ],
        "path_probabilities": point.tolist(),
        "moment_min_eigenvalue": float(np.linalg.eigvalsh(moment_value).min()),
        "moment_numerical_rank_1e-7": int(
            np.count_nonzero(np.linalg.eigvalsh(moment_value) > 1e-7)
        ),
        "status": problem.status,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lambda", dest="weight", type=float, default=0.6)
    parser.add_argument("--fixed-three-povm-weights", type=float, nargs=3, required=True)
    parser.add_argument("--prefix-order", type=int, nargs=4, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    weights = np.asarray(args.fixed_three_povm_weights, dtype=float)
    order = tuple(args.prefix_order)
    if sorted(order) != list(OUTCOMES):
        raise ValueError("prefix order must be a permutation of 0,1,2,3")
    payload = solve_povm(canonical_three_effect_povm(weights), args.weight, order)
    rendered = json.dumps(payload, indent=2) + "\n"
    print(rendered)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

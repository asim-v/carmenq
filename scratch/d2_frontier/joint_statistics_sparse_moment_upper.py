"""Sparse cross-localised moment upper bound for the ternary interior gate.

For every prefix label z and coarse outcome y, form a local moment matrix on

    {1, rho_z coordinates, G[y,*] coordinates,
     rho_z coordinate * G[y,*] coordinate}.

The sixteen matrices share moments whenever their monomials coincide.  This
is a principal sparse relaxation of the order-two commutative moment matrix.
Matrix-valued localisers impose positivity of states, pulled effects, and the
domination w_t Q_y-G[y,t] after multiplication by variables from the opposite
block.  Completeness equalities are multiplied by all degree-at-most-two
monomials available in the opposite block.

Every physical rank-one evaluation gives feasible moments, so the optimum is
a rigorous global upper bound for the pulled-statistics outer model (up to
the stated conic-solver tolerances).  The relaxation is deliberately sparse;
a high value can still be slack.
"""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import cvxpy as cp
import numpy as np

from analyze_two_block_leaf import discrimination_geometry
from terminal_weight_upper import filled_effect_weights
from two_block_choi_seesaw import (
    OUTCOMES,
    PATHS,
    canonical_three_effect_povm,
    hellinger_hypograph,
)


PAULIS = (
    np.eye(2, dtype=complex),
    np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex),
    np.array([[0.0, -1j], [1j, 0.0]], dtype=complex),
    np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex),
)
Monomial = tuple[int, ...]


def multiply(*monomials: Monomial) -> Monomial:
    return tuple(sorted(itertools.chain.from_iterable(monomials)))


def degree_two_monomials(ids: list[int]) -> list[Monomial]:
    result = [()]
    result.extend((index,) for index in ids)
    result.extend((ids[i], ids[j]) for i in range(len(ids)) for j in range(i, len(ids)))
    return result


class MomentRegistry:
    def __init__(self) -> None:
        self.keys: set[Monomial] = {()}
        self.index: dict[Monomial, int] = {}
        self.variable: cp.Variable | None = None

    def add(self, monomial: Monomial) -> None:
        self.keys.add(tuple(sorted(monomial)))

    def add_matrix(self, basis: list[Monomial]) -> None:
        for left in basis:
            for right in basis:
                self.add(multiply(left, right))

    def add_localiser(
        self,
        multipliers: list[Monomial],
        coefficient_ids: list[int],
    ) -> None:
        for left in multipliers:
            for right in multipliers:
                prefix = multiply(left, right)
                self.add(prefix)
                for index in coefficient_ids:
                    self.add(multiply(prefix, (index,)))

    def finalise(self) -> None:
        ordered = sorted(self.keys, key=lambda item: (len(item), item))
        self.index = {key: index for index, key in enumerate(ordered)}
        self.variable = cp.Variable(len(ordered))

    def moment(self, monomial: Monomial) -> cp.Expression:
        if self.variable is None:
            raise RuntimeError("moment registry has not been finalised")
        return self.variable[self.index[tuple(sorted(monomial))]]


def moment_matrix(registry: MomentRegistry, basis: list[Monomial]) -> cp.Expression:
    if registry.variable is None:
        raise RuntimeError("moment registry has not been finalised")
    # Advanced indexing produces the same affine matrix as a scalar cp.bmat,
    # but avoids hundreds of thousands of expression-tree nodes at order two.
    indices = np.asarray(
        [
            [registry.index[multiply(left, right)] for right in basis]
            for left in basis
        ],
        dtype=int,
    )
    return registry.variable[indices]


def linear_coefficient_moment(
    registry: MomentRegistry,
    prefix: Monomial,
    terms: list[tuple[float, int]],
) -> cp.Expression:
    return sum(
        coefficient * registry.moment(multiply(prefix, (index,)))
        for coefficient, index in terms
    )


def qubit_localiser(
    registry: MomentRegistry,
    multipliers: list[Monomial],
    coefficients: dict[int, list[tuple[float, int]]],
) -> cp.Expression:
    if registry.variable is None:
        raise RuntimeError("moment registry has not been finalised")
    size = len(multipliers)
    result: cp.Expression = cp.Constant(
        np.zeros((2 * size, 2 * size), dtype=complex)
    )
    for mu in range(4):
        coefficient_matrix: cp.Expression = cp.Constant(np.zeros((size, size)))
        for coefficient, variable_id in coefficients[mu]:
            indices = np.asarray(
                [
                    [
                        registry.index[
                            multiply(left, right, (variable_id,))
                        ]
                        for right in multipliers
                    ]
                    for left in multipliers
                ],
                dtype=int,
            )
            coefficient_matrix = (
                coefficient_matrix + coefficient * registry.variable[indices]
            )
        result = result + cp.kron(coefficient_matrix, PAULIS[mu]) / 2.0
    return result


def solve_povm(
    effects: np.ndarray,
    support_weight: float,
    prefix_order: tuple[int, int, int, int],
    localisers: str,
    rlt: str,
    solver: str,
    verbose: bool,
    scs_eps: float,
    scs_iters: int,
) -> dict[str, object]:
    traces = np.trace(effects, axis1=1, axis2=2).real
    active = tuple(int(s) for s in OUTCOMES if traces[s] > 1e-9)
    if active != (0, 1, 2):
        raise ValueError("this sparse model expects active labels 0,1,2")
    directions = np.asarray(
        [
            [float(np.trace((effects[s] / traces[s]) @ pauli).real) for pauli in PAULIS[1:]]
            for s in active
        ]
    )
    reconstruction = np.asarray(
        [
            [
                0.5 * traces[s],
                0.5 * traces[s] * directions[s, 0],
                0.5 * traces[s] * directions[s, 1],
            ]
            for s in active
        ]
    )
    reconstruction_inverse = np.linalg.inv(reconstruction)

    next_id = 0
    state_ids: dict[tuple[int, int], int] = {}
    for z in OUTCOMES:
        for mu in range(4):
            state_ids[z, mu] = next_id
            next_id += 1
    effect_ids: dict[tuple[int, int, int], int] = {}
    for y in OUTCOMES:
        for t in active:
            for mu in range(4):
                effect_ids[y, t, mu] = next_id
                next_id += 1

    state_group = {z: [state_ids[z, mu] for mu in range(4)] for z in OUTCOMES}
    effect_group = {
        y: [effect_ids[y, t, mu] for t in active for mu in range(4)]
        for y in OUTCOMES
    }
    bases: dict[tuple[int, int], list[Monomial]] = {}
    registry = MomentRegistry()
    for z, y in PATHS:
        cross = [
            multiply((state_id,), (effect_id,))
            for state_id in state_group[z]
            for effect_id in effect_group[y]
        ]
        basis = [(), *((index,) for index in state_group[z]), *((index,) for index in effect_group[y]), *cross]
        bases[z, y] = basis
        registry.add_matrix(basis)

    if localisers not in {"none", "basic", "full"}:
        raise ValueError("localisers must be none, basic, or full")
    if localisers != "none":
        for z, y in PATHS:
            state_multipliers = [()]
            if localisers == "basic":
                state_multipliers.extend(
                    (effect_ids[y, t, 0],) for t in active
                )
            else:
                state_multipliers.extend((index,) for index in effect_group[y])
            registry.add_localiser(state_multipliers, state_group[z])
            effect_multipliers = [(), *((index,) for index in state_group[z])]
            for t in active:
                registry.add_localiser(
                    effect_multipliers,
                    [effect_ids[y, t, mu] for mu in range(4)],
                )
                registry.add_localiser(effect_multipliers, effect_group[y])

    # Bounds and RLT equalities only require moments already present in the
    # local cross bases, but register them explicitly for clarity.
    variable_bounds: dict[int, tuple[float, float]] = {}
    for z in OUTCOMES:
        variable_bounds[state_ids[z, 0]] = (0.0, 1.0)
        for mu in range(1, 4):
            variable_bounds[state_ids[z, mu]] = (-1.0, 1.0)
    for y in OUTCOMES:
        for t in active:
            bound = 2.0 * traces[t]
            variable_bounds[effect_ids[y, t, 0]] = (0.0, bound)
            for mu in range(1, 4):
                variable_bounds[effect_ids[y, t, mu]] = (-bound, bound)
    for index in variable_bounds:
        registry.add((index,))
        registry.add((index, index))
    for z in OUTCOMES:
        for monomial in degree_two_monomials(state_group[z]):
            registry.add(monomial)
            for y in OUTCOMES:
                for t in active:
                    for mu in range(4):
                        registry.add(multiply(monomial, (effect_ids[y, t, mu],)))
    for y in OUTCOMES:
        for monomial in degree_two_monomials(effect_group[y]):
            registry.add(monomial)
            for z in OUTCOMES:
                registry.add(multiply(monomial, (state_ids[z, 0],)))
    registry.finalise()

    constraints: list[cp.Constraint] = [registry.moment(()) == 1.0]
    local_matrices = []
    for z, y in PATHS:
        matrix = moment_matrix(registry, bases[z, y])
        local_matrices.append(matrix)
        constraints.append(matrix >> 0)

    for index, (lower, upper) in variable_bounds.items():
        first = registry.moment((index,))
        second = registry.moment((index, index))
        constraints.extend((first >= lower, first <= upper))
        constraints.append(second <= (lower + upper) * first - lower * upper)

    def state_coefficients(z: int) -> dict[int, list[tuple[float, int]]]:
        return {mu: [(1.0, state_ids[z, mu])] for mu in range(4)}

    def effect_coefficients(y: int, t: int) -> dict[int, list[tuple[float, int]]]:
        return {mu: [(1.0, effect_ids[y, t, mu])] for mu in range(4)}

    def residual_coefficients(y: int, t: int) -> dict[int, list[tuple[float, int]]]:
        return {
            mu: [
                *((traces[t], effect_ids[y, other, mu]) for other in active),
                (-1.0, effect_ids[y, t, mu]),
            ]
            for mu in range(4)
        }

    if localisers != "none":
        for z, y in PATHS:
            state_multipliers = [()]
            if localisers == "basic":
                state_multipliers.extend((effect_ids[y, t, 0],) for t in active)
            else:
                state_multipliers.extend((index,) for index in effect_group[y])
            constraints.append(
                qubit_localiser(
                    registry, state_multipliers, state_coefficients(z)
                )
                >> 0
            )
            effect_multipliers = [(), *((index,) for index in state_group[z])]
            for t in active:
                constraints.append(
                    qubit_localiser(
                        registry,
                        effect_multipliers,
                        effect_coefficients(y, t),
                    )
                    >> 0
                )
                constraints.append(
                    qubit_localiser(
                        registry,
                        effect_multipliers,
                        residual_coefficients(y, t),
                    )
                    >> 0
                )
    else:
        for z in OUTCOMES:
            constraints.append(
                qubit_localiser(registry, [()], state_coefficients(z)) >> 0
            )
        for y in OUTCOMES:
            for t in active:
                constraints.extend(
                    (
                        qubit_localiser(
                            registry, [()], effect_coefficients(y, t)
                        )
                        >> 0,
                        qubit_localiser(
                            registry, [()], residual_coefficients(y, t)
                        )
                        >> 0,
                    )
                )

    if rlt not in {"none", "single", "degree2"}:
        raise ValueError("rlt must be none, single, or degree2")
    # Normalisation and completeness.  Optional RLT products use monomials
    # from the opposite block; omitting them preserves validity and often
    # restores strict feasibility for interior-point solvers.
    constraints.append(
        sum(registry.moment((state_ids[z, 0],)) for z in OUTCOMES) == 1.0
    )
    if rlt != "none":
        for y in OUTCOMES:
            effect_monomials = [(index,) for index in effect_group[y]]
            if rlt == "degree2":
                effect_monomials.extend(
                    (effect_group[y][i], effect_group[y][j])
                    for i in range(len(effect_group[y]))
                    for j in range(i, len(effect_group[y]))
                )
            for monomial in effect_monomials:
                constraints.append(
                    sum(
                        registry.moment(multiply(monomial, (state_ids[z, 0],)))
                        for z in OUTCOMES
                    )
                    == registry.moment(monomial)
                )
    for mu in range(4):
        target = 2.0 if mu == 0 else 0.0
        constraints.append(
            sum(
                registry.moment((effect_ids[y, t, mu],))
                for y in OUTCOMES
                for t in active
            )
            == target
        )
        if rlt != "none":
            for z in OUTCOMES:
                state_monomials = [(index,) for index in state_group[z]]
                if rlt == "degree2":
                    state_monomials.extend(
                        (state_group[z][i], state_group[z][j])
                        for i in range(len(state_group[z]))
                        for j in range(i, len(state_group[z]))
                    )
                for monomial in state_monomials:
                    constraints.append(
                        sum(
                            registry.moment(
                                multiply(monomial, (effect_ids[y, t, mu],))
                            )
                            for y in OUTCOMES
                            for t in active
                        )
                        == target * registry.moment(monomial)
                    )

    statistics: dict[tuple[int, int, int], cp.Expression] = {}
    probabilities: dict[tuple[int, int], cp.Expression] = {}
    for z, y in PATHS:
        for t in active:
            statistics[z, y, t] = 0.5 * sum(
                registry.moment(
                    multiply((state_ids[z, mu],), (effect_ids[y, t, mu],))
                )
                for mu in range(4)
            )
            constraints.append(statistics[z, y, t] >= 0.0)
        probabilities[z, y] = sum(statistics[z, y, t] for t in active)
        constraints.append(probabilities[z, y] >= 0.0)
        for t in active:
            constraints.append(
                statistics[z, y, t] <= traces[t] * probabilities[z, y]
            )
    constraints.append(sum(probabilities.values()) == 1.0)

    prefix = [registry.moment((state_ids[z, 0],)) for z in OUTCOMES]
    constraints.extend(
        sum(probabilities[z, y] for y in OUTCOMES) == prefix[z]
        for z in OUTCOMES
    )
    constraints.extend(
        prefix[prefix_order[index]] >= prefix[prefix_order[index + 1]]
        for index in range(3)
    )

    terminal_statistics = [
        [
            sum(
                statistics[z, y, t]
                for z, y in PATHS
                if (z ^ y) == syndrome
            )
            for t in active
        ]
        for syndrome in OUTCOMES
    ]
    reconstructed = [
        reconstruction_inverse @ cp.hstack(terminal_statistics[s])
        for s in OUTCOMES
    ]
    syndrome_priors_paths = [
        sum(probabilities[z, z ^ s] for z in OUTCOMES) for s in OUTCOMES
    ]
    syndrome_priors = [reconstructed[s][0] for s in OUTCOMES]
    constraints.extend(
        syndrome_priors[s] == syndrome_priors_paths[s] for s in OUTCOMES
    )
    terminal_normal = cp.Variable(4)
    terminal_vectors = [
        cp.hstack([reconstructed[s][1], reconstructed[s][2], terminal_normal[s]])
        for s in OUTCOMES
    ]
    constraints.extend(
        cp.SOC(syndrome_priors[s], terminal_vectors[s]) for s in OUTCOMES
    )

    audit = sum(terminal_statistics[s][s] for s in active)
    cap = filled_effect_weights(float(traces.max()))
    constraints.append(
        audit
        <= sum(cap[index] * prefix[prefix_order[index]] for index in OUTCOMES)
    )
    dual_trace = cp.Variable(nonneg=True)
    dual_vector = cp.Variable(3)
    constraints.append(cp.SOC(dual_trace, dual_vector))
    constraints.extend(
        cp.SOC(
            dual_trace - syndrome_priors[s],
            dual_vector - terminal_vectors[s],
        )
        for s in OUTCOMES
    )
    constraints.append(audit == dual_trace)

    returned = hellinger_hypograph(
        [probabilities[z, y] for z, y in PATHS], constraints
    )
    problem = cp.Problem(
        cp.Maximize(support_weight * audit + (1.0 - support_weight) * returned),
        constraints,
    )
    if solver == "clarabel":
        problem.solve(
            solver="CLARABEL",
            tol_gap_abs=5e-8,
            tol_gap_rel=5e-8,
            tol_feas=5e-8,
            max_iter=1000,
            verbose=verbose,
        )
    elif solver == "scs":
        problem.solve(
            solver="SCS",
            eps=scs_eps,
            max_iters=scs_iters,
            acceleration_lookback=20,
            verbose=verbose,
        )
    else:
        raise ValueError("solver must be clarabel or scs")
    if problem.status not in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}:
        raise RuntimeError(f"sparse moment relaxation failed: {problem.status}")

    point = np.asarray(
        [[float(probabilities[z, y].value) for y in OUTCOMES] for z in OUTCOMES]
    )
    terminal_prior_value = np.asarray([float(item.value) for item in syndrome_priors])
    terminal_vector_value = np.asarray(
        [np.asarray(item.value, dtype=float) for item in terminal_vectors]
    )
    terminal_states = np.asarray(
        [
            0.5
            * (
                terminal_prior_value[s] * PAULIS[0]
                + sum(
                    terminal_vector_value[s, axis] * PAULIS[axis + 1]
                    for axis in range(3)
                )
            )
            for s in OUTCOMES
        ]
    )
    geometry = discrimination_geometry(terminal_states)
    audit_value = float(audit.value)
    return_value = float(np.sqrt(np.maximum(point, 0.0)).sum() ** 2 / 16.0)
    local_minimum = min(
        float(np.linalg.eigvalsh(np.asarray(matrix.value, dtype=float)).min())
        for matrix in local_matrices
    )
    return {
        "weight": support_weight,
        "terminal_effect_weights": traces.tolist(),
        "prefix_order": list(prefix_order),
        "localisers": localisers,
        "rlt": rlt,
        "solver": solver,
        "scalar_moment_count": len(registry.index),
        "local_moment_size": len(next(iter(bases.values()))),
        "bound": float(problem.value),
        "objective_from_reported": support_weight * audit_value
        + (1.0 - support_weight) * return_value,
        "audit": audit_value,
        "return": return_value,
        "normalisation": float(point.sum()),
        "prefix_priors": [float(item.value) for item in prefix],
        "syndrome_priors": terminal_prior_value.tolist(),
        "terminal_statistics": [
            [float(item.value) for item in terminal_statistics[s]] for s in OUTCOMES
        ],
        "terminal_bloch_vectors": terminal_vector_value.tolist(),
        "independent_terminal_geometry": geometry,
        "helstrom_independent_residual": float(
            audit_value - float(geometry["optimal_guess_probability"])
        ),
        "path_probabilities": point.tolist(),
        "local_moment_min_eigenvalue": local_minimum,
        "status": problem.status,
        "solver_stats": {
            "solve_time": problem.solver_stats.solve_time,
            "num_iters": problem.solver_stats.num_iters,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lambda", dest="weight", type=float, default=0.6)
    parser.add_argument("--fixed-three-povm-weights", type=float, nargs=3, required=True)
    parser.add_argument("--prefix-order", type=int, nargs=4, required=True)
    parser.add_argument("--localisers", choices=("none", "basic", "full"), default="basic")
    parser.add_argument("--rlt", choices=("none", "single", "degree2"), default="single")
    parser.add_argument("--solver", choices=("clarabel", "scs"), default="clarabel")
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--scs-eps", type=float, default=1e-5)
    parser.add_argument("--scs-iters", type=int, default=300_000)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    order = tuple(args.prefix_order)
    if sorted(order) != list(OUTCOMES):
        raise ValueError("prefix order must be a permutation of 0,1,2,3")
    payload = solve_povm(
        canonical_three_effect_povm(np.asarray(args.fixed_three_povm_weights)),
        args.weight,
        order,
        args.localisers,
        args.rlt,
        args.solver,
        args.verbose,
        args.scs_eps,
        args.scs_iters,
    )
    rendered = json.dumps(payload, indent=2) + "\n"
    print(rendered)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

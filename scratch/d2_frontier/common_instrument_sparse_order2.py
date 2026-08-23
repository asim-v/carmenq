"""Sparse order-two Choi hierarchy for one common quantum instrument.

The first Choi--Helstrom relaxation stores only first and second scalar
moments of all state and Choi coordinates.  It can therefore replace a
single deterministic instrument by a correlated pseudo-distribution.  This
file lifts the products that define the conditioned outputs to local
order-two moment matrices.

There are four subnormalised input states ``rho[z]`` and four Choi matrices
``J[y]``.  Every pair ``(z,y)`` receives the cross basis

    {1, rho[z], J[y], rho[z] * J[y]},

while one bridge per outcome contains *all* four input blocks and the same
``J[y]`` block.  Choi positivity is localised against every input coordinate,
and trace preservation is multiplied by every state monomial of degree at
most two.  Thus the four conditioned families cannot silently use four
independent instruments at this relaxation level.

Every physical common instrument induces a rank-one feasible moment sequence,
so the optimum is a valid upper relaxation for the fixed terminal POVM and
prefix-prior order, up to conic-solver tolerance.  Failure to reach a target
is not a certificate until the numerical dual has also been audited.
"""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path
from typing import Iterable

import cvxpy as cp
import numpy as np

from analyze_two_block_leaf import discrimination_geometry
from joint_statistics_sparse_moment_upper import (
    PAULIS,
    MomentRegistry,
    moment_matrix,
    multiply,
    qubit_localiser,
)
from terminal_weight_upper import filled_effect_weights
from two_block_choi_seesaw import (
    OUTCOMES,
    PATHS,
    canonical_three_effect_povm,
    hellinger_hypograph,
)


Monomial = tuple[int, ...]
TRANSPOSE_SIGN = np.asarray([1.0, 1.0, -1.0, 1.0])


def monomials(ids: Iterable[int], maximum_degree: int) -> list[Monomial]:
    """Return commutative monomials of degree at most ``maximum_degree``."""

    ordered = list(ids)
    result: list[Monomial] = [()]
    for degree in range(1, maximum_degree + 1):
        result.extend(itertools.combinations_with_replacement(ordered, degree))
    return result


def unique_monomials(items: Iterable[Monomial]) -> list[Monomial]:
    """Deduplicate monomials while retaining deterministic insertion order."""

    return list(dict.fromkeys(tuple(sorted(item)) for item in items))


def choi_output_coefficients(
    state_coefficients: np.ndarray, choi_coefficients: np.ndarray
) -> np.ndarray:
    """Apply Pauli-expanded Choi data and return output Pauli coefficients."""

    state = np.asarray(state_coefficients, dtype=float)
    choi = np.asarray(choi_coefficients, dtype=float)
    if state.shape != (4,) or choi.shape != (4, 4):
        raise ValueError("expected state shape (4,) and Choi shape (4,4)")
    return 0.5 * np.einsum("m,m,mn->n", TRANSPOSE_SIGN, state, choi)


def register_linear_ideal(
    registry: MomentRegistry,
    multipliers: Iterable[Monomial],
    variable_ids: Iterable[int],
) -> None:
    for multiplier in multipliers:
        registry.add(multiplier)
        for variable_id in variable_ids:
            registry.add(multiply(multiplier, (variable_id,)))


def add_linear_ideal(
    registry: MomentRegistry,
    constraints: list[cp.Constraint],
    multipliers: Iterable[Monomial],
    terms: list[tuple[float, int]],
    constant: float,
) -> None:
    for multiplier in multipliers:
        constraints.append(
            sum(
                coefficient
                * registry.moment(multiply(multiplier, (variable_id,)))
                for coefficient, variable_id in terms
            )
            == constant * registry.moment(multiplier)
        )


def matrix_localiser(
    registry: MomentRegistry,
    multipliers: list[Monomial],
    coefficient_ids: dict[tuple[int, int], int],
) -> cp.Expression:
    """Return the localising matrix for a Pauli-expanded two-qubit matrix."""

    if registry.variable is None:
        raise RuntimeError("moment registry has not been finalised")
    size = len(multipliers)
    result: cp.Expression = cp.Constant(
        np.zeros((4 * size, 4 * size), dtype=complex)
    )
    for alpha in range(4):
        for nu in range(4):
            variable_id = coefficient_ids[alpha, nu]
            indices = np.asarray(
                [
                    [
                        registry.index[multiply(left, right, (variable_id,))]
                        for right in multipliers
                    ]
                    for left in multipliers
                ],
                dtype=int,
            )
            coefficient_matrix = registry.variable[indices]
            result = result + cp.kron(
                coefficient_matrix, np.kron(PAULIS[alpha], PAULIS[nu])
            ) / 4.0
    return cp.hermitian_wrap(result)


def register_matrix_localiser(
    registry: MomentRegistry,
    multipliers: list[Monomial],
    variable_ids: Iterable[int],
) -> None:
    registry.add_localiser(multipliers, list(variable_ids))


def solve_povm(
    terminal: np.ndarray,
    support_weight: float,
    prefix_order: tuple[int, int, int, int],
    bridge_basis_mode: str,
    localisers: str,
    ideal_degree: int,
    solver: str,
    verbose: bool,
    scs_eps: float,
    scs_iters: int,
) -> dict[str, object]:
    """Solve the sparse common-instrument upper relaxation."""

    weights = np.trace(terminal[:3], axis1=1, axis2=2).real
    if np.any(weights <= 0.0):
        raise ValueError("all three terminal effects must be active")
    directions = np.asarray(
        [
            [
                float(np.trace((terminal[t] / weights[t]) @ pauli).real)
                for pauli in PAULIS[1:]
            ]
            for t in range(3)
        ]
    )
    if bridge_basis_mode not in {"minimal", "cross", "full"}:
        raise ValueError("bridge_basis_mode must be minimal, cross, or full")
    if localisers not in {"basic", "common"}:
        raise ValueError("localisers must be basic or common")
    if ideal_degree not in {0, 1, 2}:
        raise ValueError("ideal_degree must be zero, one, or two")

    next_id = 0
    state_ids: dict[tuple[int, int], int] = {}
    for z in OUTCOMES:
        for mu in range(4):
            state_ids[z, mu] = next_id
            next_id += 1
    choi_ids: dict[tuple[int, int, int], int] = {}
    for y in OUTCOMES:
        for alpha in range(4):
            for nu in range(4):
                choi_ids[y, alpha, nu] = next_id
                next_id += 1

    state_group = {
        z: [state_ids[z, mu] for mu in range(4)] for z in OUTCOMES
    }
    all_state_ids = [state_ids[z, mu] for z in OUTCOMES for mu in range(4)]
    choi_group = {
        y: [choi_ids[y, alpha, nu] for alpha in range(4) for nu in range(4)]
        for y in OUTCOMES
    }
    trace_preservation_ids = [
        choi_ids[y, alpha, 0] for y in OUTCOMES for alpha in range(4)
    ]

    state_degree2 = monomials(all_state_ids, 2)
    choi_degree2 = {y: monomials(choi_group[y], 2) for y in OUTCOMES}
    pair_bases: dict[tuple[int, int], list[Monomial]] = {}
    for z, y in PATHS:
        pair_bases[z, y] = [
            (),
            *((variable_id,) for variable_id in state_group[z]),
            *((variable_id,) for variable_id in choi_group[y]),
            *(
                (state_id, choi_id)
                for state_id in state_group[z]
                for choi_id in choi_group[y]
            ),
        ]

    bridge_bases: dict[int, list[Monomial]] = {}
    for y in OUTCOMES:
        if bridge_basis_mode == "minimal":
            bridge = [
                *state_degree2,
                *((variable_id,) for variable_id in choi_group[y]),
            ]
        elif bridge_basis_mode == "cross":
            bridge = [
                *state_degree2,
                *((variable_id,) for variable_id in choi_group[y]),
                *(
                    (state_id, choi_id)
                    for state_id in all_state_ids
                    for choi_id in choi_group[y]
                ),
            ]
        else:
            bridge = monomials([*all_state_ids, *choi_group[y]], 2)
        bridge_bases[y] = unique_monomials(bridge)

    registry = MomentRegistry()
    registry.add_matrix(state_degree2)
    instrument_bridge_basis = unique_monomials(
        [
            *state_degree2,
            *((variable_id,) for variable_id in trace_preservation_ids),
        ]
    )
    registry.add_matrix(instrument_bridge_basis)
    for y in OUTCOMES:
        registry.add_matrix(choi_degree2[y])
        registry.add_matrix(bridge_bases[y])
    if bridge_basis_mode == "minimal":
        for basis in pair_bases.values():
            registry.add_matrix(basis)

    state_coefficients = {
        z: {mu: [(1.0, state_ids[z, mu])] for mu in range(4)}
        for z in OUTCOMES
    }
    choi_coefficients = {
        y: {(alpha, nu): choi_ids[y, alpha, nu] for alpha in range(4) for nu in range(4)}
        for y in OUTCOMES
    }

    localiser_specs: list[tuple[str, int, int | None, list[Monomial]]] = []
    if localisers == "basic":
        for z in OUTCOMES:
            localiser_specs.append(("state", z, None, [()]))
        for y in OUTCOMES:
            localiser_specs.append(("choi", y, None, [()]))
    else:
        for z, y in PATHS:
            multipliers = [(), *((variable_id,) for variable_id in choi_group[y])]
            localiser_specs.append(("state", z, y, multipliers))
            registry.add_localiser(multipliers, state_group[z])
        for y in OUTCOMES:
            multipliers = [(), *((variable_id,) for variable_id in all_state_ids)]
            localiser_specs.append(("choi", y, None, multipliers))
            register_matrix_localiser(registry, multipliers, choi_group[y])
    if localisers == "basic":
        for z in OUTCOMES:
            registry.add_localiser([()], state_group[z])
        for y in OUTCOMES:
            register_matrix_localiser(registry, [()], choi_group[y])

    variable_bounds: dict[int, tuple[float, float]] = {}
    trace_upper_by_rank = (1.0, 0.5, 1.0 / 3.0, 0.25)
    trace_upper = {
        prefix_order[rank]: trace_upper_by_rank[rank] for rank in range(4)
    }
    for rank, z in enumerate(prefix_order):
        lower = 0.25 if rank == 0 else 0.0
        upper = trace_upper[z]
        variable_bounds[state_ids[z, 0]] = (lower, upper)
        for mu in range(1, 4):
            variable_bounds[state_ids[z, mu]] = (-upper, upper)
    for y in OUTCOMES:
        variable_bounds[choi_ids[y, 0, 0]] = (0.0, 2.0)
        for alpha in range(4):
            for nu in range(4):
                if (alpha, nu) != (0, 0):
                    variable_bounds[choi_ids[y, alpha, nu]] = (-2.0, 2.0)
    for variable_id in variable_bounds:
        registry.add((variable_id,))
        registry.add((variable_id, variable_id))

    # State normalisation is multiplied by local Choi monomials, while trace
    # preservation is multiplied by the full four-input state block.  These
    # are the identities that directly prevent input-dependent instruments.
    state_normalisation_multipliers = unique_monomials(
        itertools.chain(
            monomials(all_state_ids, ideal_degree),
            *(
                monomials(choi_group[y], ideal_degree)
                for y in OUTCOMES
            ),
            *(
                (
                    (state_id, choi_id)
                    for state_id in all_state_ids
                    for choi_id in choi_group[y]
                )
                for y in OUTCOMES
                if ideal_degree >= 2
            ),
        )
    )
    trace_preservation_multipliers = unique_monomials(
        [
            *monomials(all_state_ids, ideal_degree),
            *(
                (variable_id,)
                for variable_id in trace_preservation_ids
                if ideal_degree >= 1
            ),
        ]
    )
    state_trace_ids = [state_ids[z, 0] for z in OUTCOMES]
    register_linear_ideal(
        registry, state_normalisation_multipliers, state_trace_ids
    )
    for alpha in range(4):
        register_linear_ideal(
            registry,
            trace_preservation_multipliers,
            [choi_ids[y, alpha, 0] for y in OUTCOMES],
        )

    registry.finalise()
    constraints: list[cp.Constraint] = [registry.moment(()) == 1.0]
    moment_matrices: list[cp.Expression] = []

    state_matrix = moment_matrix(registry, state_degree2)
    instrument_bridge_matrix = moment_matrix(registry, instrument_bridge_basis)
    constraints.extend((state_matrix >> 0, instrument_bridge_matrix >> 0))
    moment_matrices.extend((state_matrix, instrument_bridge_matrix))
    for y in OUTCOMES:
        choi_matrix = moment_matrix(registry, choi_degree2[y])
        bridge_matrix = moment_matrix(registry, bridge_bases[y])
        constraints.extend((choi_matrix >> 0, bridge_matrix >> 0))
        moment_matrices.extend((choi_matrix, bridge_matrix))
    if bridge_basis_mode == "minimal":
        for basis in pair_bases.values():
            matrix = moment_matrix(registry, basis)
            constraints.append(matrix >> 0)
            moment_matrices.append(matrix)

    for variable_id, (lower, upper) in variable_bounds.items():
        first = registry.moment((variable_id,))
        second = registry.moment((variable_id, variable_id))
        constraints.extend((first >= lower, first <= upper))
        constraints.append(second <= (lower + upper) * first - lower * upper)

    localiser_matrices: list[cp.Expression] = []
    for kind, first_index, _, multipliers in localiser_specs:
        if kind == "state":
            matrix = qubit_localiser(
                registry, multipliers, state_coefficients[first_index]
            )
        else:
            matrix = matrix_localiser(
                registry, multipliers, choi_coefficients[first_index]
            )
        constraints.append(matrix >> 0)
        localiser_matrices.append(matrix)

    add_linear_ideal(
        registry,
        constraints,
        state_normalisation_multipliers,
        [(1.0, variable_id) for variable_id in state_trace_ids],
        1.0,
    )
    for alpha in range(4):
        add_linear_ideal(
            registry,
            constraints,
            trace_preservation_multipliers,
            [(1.0, choi_ids[y, alpha, 0]) for y in OUTCOMES],
            2.0 if alpha == 0 else 0.0,
        )

    probabilities: dict[tuple[int, int], cp.Expression] = {}
    output_vectors: dict[tuple[int, int], cp.Expression] = {}
    for z, y in PATHS:
        probabilities[z, y] = 0.5 * sum(
            TRANSPOSE_SIGN[mu]
            * registry.moment(
                multiply((state_ids[z, mu],), (choi_ids[y, mu, 0],))
            )
            for mu in range(4)
        )
        output_vectors[z, y] = cp.hstack(
            [
                0.5
                * sum(
                    TRANSPOSE_SIGN[mu]
                    * registry.moment(
                        multiply(
                            (state_ids[z, mu],),
                            (choi_ids[y, mu, nu],),
                        )
                    )
                    for mu in range(4)
                )
                for nu in range(1, 4)
            ]
        )
        constraints.append(cp.SOC(probabilities[z, y], output_vectors[z, y]))

    prefix = [registry.moment((state_ids[z, 0],)) for z in OUTCOMES]
    constraints.extend(
        sum(probabilities[z, y] for y in OUTCOMES) == prefix[z]
        for z in OUTCOMES
    )
    constraints.extend(
        prefix[prefix_order[index]] >= prefix[prefix_order[index + 1]]
        for index in range(3)
    )
    constraints.append(sum(probabilities.values()) == 1.0)

    syndrome_priors = [
        sum(probabilities[z, z ^ syndrome] for z in OUTCOMES)
        for syndrome in OUTCOMES
    ]
    terminal_vectors = [
        sum(output_vectors[z, z ^ syndrome] for z in OUTCOMES)
        for syndrome in OUTCOMES
    ]
    constraints.extend(
        cp.SOC(syndrome_priors[syndrome], terminal_vectors[syndrome])
        for syndrome in OUTCOMES
    )

    audit = sum(
        0.5
        * weights[syndrome]
        * (
            syndrome_priors[syndrome]
            + directions[syndrome] @ terminal_vectors[syndrome]
        )
        for syndrome in range(3)
    )
    cap = filled_effect_weights(float(weights.max()))
    constraints.append(
        audit
        <= sum(cap[index] * prefix[prefix_order[index]] for index in OUTCOMES)
    )

    # The fixed nonzero effects are rank one.  This exact facial form of the
    # Helstrom conditions is both stronger numerically and less singular than
    # four bare Lorentz inequalities plus an objective equality.
    dual_scalar = cp.Variable(nonneg=True)
    dual_vector = cp.Variable(3)
    slack = cp.Variable(3, nonneg=True)
    for syndrome in range(3):
        constraints.extend(
            (
                dual_scalar - syndrome_priors[syndrome] == slack[syndrome],
                dual_vector - terminal_vectors[syndrome]
                == -slack[syndrome] * directions[syndrome],
            )
        )
    constraints.append(
        cp.SOC(
            dual_scalar - syndrome_priors[3],
            dual_vector - terminal_vectors[3],
        )
    )
    constraints.append(audit == dual_scalar)

    returned = hellinger_hypograph(
        [probabilities[z, y] for z, y in PATHS], constraints
    )
    objective = support_weight * audit + (1.0 - support_weight) * returned
    problem = cp.Problem(cp.Maximize(objective), constraints)
    if solver == "clarabel":
        problem.solve(
            solver="CLARABEL",
            tol_gap_abs=2e-7,
            tol_gap_rel=2e-7,
            tol_feas=2e-7,
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
        raise RuntimeError(f"sparse Choi hierarchy failed: {problem.status}")

    point = np.asarray(
        [[float(probabilities[z, y].value) for y in OUTCOMES] for z in OUTCOMES]
    )
    terminal_prior_value = np.asarray(
        [float(item.value) for item in syndrome_priors]
    )
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
    hypograph_value = float(returned.value)
    prefix_bloch_value = np.asarray(
        [
            [float(registry.moment((state_ids[z, mu],)).value) for mu in range(4)]
            for z in OUTCOMES
        ]
    )
    conditioned_vector_value = np.asarray(
        [
            [np.asarray(output_vectors[z, y].value, dtype=float) for y in OUTCOMES]
            for z in OUTCOMES
        ]
    )
    moment_minimum = min(
        float(np.linalg.eigvalsh(np.asarray(matrix.value, dtype=float)).min())
        for matrix in moment_matrices
    )
    localiser_minimum = min(
        float(np.linalg.eigvalsh(np.asarray(matrix.value, dtype=complex)).min())
        for matrix in localiser_matrices
    )
    return {
        "weight": support_weight,
        "terminal_effect_weights": [*weights.tolist(), 0.0],
        "prefix_order": list(prefix_order),
        "bridge_basis_mode": bridge_basis_mode,
        "localisers": localisers,
        "ideal_degree": ideal_degree,
        "solver": solver,
        "state_variable_count": len(all_state_ids),
        "choi_variable_count": sum(len(group) for group in choi_group.values()),
        "pair_moment_size": len(next(iter(pair_bases.values()))),
        "state_moment_size": len(state_degree2),
        "choi_moment_size": len(next(iter(choi_degree2.values()))),
        "bridge_moment_size": len(next(iter(bridge_bases.values()))),
        "instrument_bridge_moment_size": len(instrument_bridge_basis),
        "scalar_moment_count": len(registry.index),
        "bound": float(problem.value),
        "objective_from_hypograph": support_weight * audit_value
        + (1.0 - support_weight) * hypograph_value,
        "objective_from_reported": support_weight * audit_value
        + (1.0 - support_weight) * return_value,
        "audit": audit_value,
        "return": return_value,
        "return_hypograph_value": hypograph_value,
        "normalisation": float(point.sum()),
        "prefix_priors": [float(item.value) for item in prefix],
        "prefix_bloch_coefficients": prefix_bloch_value.tolist(),
        "path_probabilities": point.tolist(),
        "conditioned_output_bloch_vectors": conditioned_vector_value.tolist(),
        "syndrome_priors": terminal_prior_value.tolist(),
        "terminal_bloch_vectors": terminal_vector_value.tolist(),
        "independent_terminal_geometry": geometry,
        "helstrom_independent_residual": float(
            audit_value - float(geometry["optimal_guess_probability"])
        ),
        "moment_min_eigenvalue": moment_minimum,
        "localiser_min_eigenvalue": localiser_minimum,
        "status": problem.status,
        "solver_stats": {
            "solve_time": problem.solver_stats.solve_time,
            "num_iters": problem.solver_stats.num_iters,
        },
        "scope": (
            "fixed terminal POVM and prefix-prior order; sparse order-two "
            "common-instrument upper relaxation"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lambda", dest="weight", type=float, default=0.55)
    parser.add_argument(
        "--fixed-three-povm-weights", type=float, nargs=3, required=True
    )
    parser.add_argument("--prefix-order", type=int, nargs=4, required=True)
    parser.add_argument(
        "--bridge-basis", choices=("minimal", "cross", "full"), default="minimal"
    )
    parser.add_argument(
        "--localisers", choices=("basic", "common"), default="common"
    )
    parser.add_argument("--ideal-degree", type=int, choices=(0, 1, 2), default=2)
    parser.add_argument("--solver", choices=("clarabel", "scs"), default="scs")
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--scs-eps", type=float, default=1e-4)
    parser.add_argument("--scs-iters", type=int, default=200_000)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    order = tuple(args.prefix_order)
    if sorted(order) != list(OUTCOMES):
        raise ValueError("prefix order must be a permutation of 0,1,2,3")
    payload = solve_povm(
        canonical_three_effect_povm(
            np.asarray(args.fixed_three_povm_weights, dtype=float)
        ),
        args.weight,
        order,
        args.bridge_basis,
        args.localisers,
        args.ideal_degree,
        args.solver,
        args.verbose,
        args.scs_eps,
        args.scs_iters,
    )
    rendered = json.dumps(payload, indent=2) + "\n"
    print(rendered, end="")
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

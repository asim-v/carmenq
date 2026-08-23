"""Order-two clique moment bound for the nonprojective terminal sector.

The earlier pathwise sparse relaxation gave every prefix label its own local
state geometry.  This file instead builds one degree-two moment matrix for
each coarse continuation outcome ``y``.  Every clique contains *all four*
prefix states and all three pulled terminal effects ``G[y,t]``.  Pure-state
moments are shared across all four cliques, so incompatible copies of the
same four-dimensional Bloch geometry are no longer possible.

The physical variables are subnormalised qubit states

    rho_z = (a_z I + r_z.sigma) / 2

and a twelve-effect qubit POVM

    G[y,t] = (g[y,t,0] I + g[y,t].sigma) / 2.

The common SO(3) gauge ``r0_y=r0_z=r1_z=0`` is eliminated exactly, leaving
thirteen state coordinates.  The default anisotropic basis contains every
state monomial of degree at most two and every degree-one effect coordinate,
giving four 117 by 117 clique matrices.  A targeted basis additionally lifts
only the state--effect products named by selected behaviour columns.  It
retains the degree-three/four moments needed to attack a certified obstruction
without paying for the much larger 351 by 351 full order-two matrices.  Exact
multiplied normalisation/completeness ideals and matrix-valued localisers
impose qubit positivity for states, effects, and domination residuals.

Selected probability columns are linked to their Born moments.  The
remaining path/terminal statistics stay free, so this is an upper relaxation
of the selected-column model.  A value below a target is meaningful only
after solver residuals and a dual witness have been audited; a high value
merely means that this hierarchy level is still slack.
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
from joint_statistics_helstrom_seesaw import reconstruction_matrix
from joint_statistics_sparse_moment_upper import (
    PAULIS,
    MomentRegistry,
    moment_matrix,
    multiply,
    qubit_localiser,
)
from linked_statistics_seesaw import parse_column
from terminal_weight_upper import filled_effect_weights
from two_block_choi_seesaw import (
    OUTCOMES,
    PATHS,
    canonical_three_effect_povm,
    hellinger_hypograph,
)


Monomial = tuple[int, ...]
LinearCoefficients = dict[int, list[tuple[float, int]]]


def monomials(ids: Iterable[int], maximum_degree: int) -> list[Monomial]:
    """Return all commutative monomials of degree at most ``maximum_degree``."""

    ordered = list(ids)
    result: list[Monomial] = [()]
    for degree in range(1, maximum_degree + 1):
        result.extend(itertools.combinations_with_replacement(ordered, degree))
    return result


def add_linear_ideal(
    registry: MomentRegistry,
    constraints: list[cp.Constraint],
    multipliers: Iterable[Monomial],
    terms: list[tuple[float, int]],
    constant: float,
) -> None:
    """Impose ``(sum coefficient*x - constant) * multiplier == 0``."""

    for multiplier in multipliers:
        constraints.append(
            sum(
                coefficient
                * registry.moment(multiply(multiplier, (variable_id,)))
                for coefficient, variable_id in terms
            )
            == constant * registry.moment(multiplier)
        )


def register_linear_ideal(
    registry: MomentRegistry,
    multipliers: Iterable[Monomial],
    variable_ids: Iterable[int],
) -> None:
    for multiplier in multipliers:
        registry.add(multiplier)
        for variable_id in variable_ids:
            registry.add(multiply(multiplier, (variable_id,)))


def coefficient_ids(coefficients: LinearCoefficients) -> list[int]:
    return sorted(
        {
            variable_id
            for terms in coefficients.values()
            for _, variable_id in terms
        }
    )


def add_qubit_localiser_registry(
    registry: MomentRegistry,
    multipliers: list[Monomial],
    coefficients: LinearCoefficients,
) -> None:
    registry.add_localiser(multipliers, coefficient_ids(coefficients))


def solve_povm(
    terminal: np.ndarray,
    support_weight: float,
    prefix_order: tuple[int, int, int, int],
    linked: tuple[tuple[str, int, int], ...],
    basis_mode: str,
    localiser_scope: str,
    ideal_degree: int,
    helstrom_form: str,
    bridge_selected: bool,
    bridge_basis_mode: str,
    solver: str,
    verbose: bool,
    scs_eps: float,
    scs_iters: int,
) -> dict[str, object]:
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

    next_id = 0
    state_scalar: dict[int, int] = {}
    state_vector: dict[tuple[int, int], int | None] = {}
    for z in OUTCOMES:
        state_scalar[z] = next_id
        next_id += 1
    # Exact rotational gauge: r0=(x,0,0), r1=(x,y,0).
    for z in OUTCOMES:
        for axis in range(3):
            eliminated = (z == 0 and axis in (1, 2)) or (z == 1 and axis == 2)
            if eliminated:
                state_vector[z, axis] = None
            else:
                state_vector[z, axis] = next_id
                next_id += 1
    state_ids = list(range(next_id))
    if len(state_ids) != 13:
        raise AssertionError("the gauge-reduced state block must have 13 coordinates")

    effect_ids: dict[tuple[int, int, int], int] = {}
    for y in OUTCOMES:
        for t in range(3):
            for mu in range(4):
                effect_ids[y, t, mu] = next_id
                next_id += 1
    effect_groups = {
        y: [effect_ids[y, t, mu] for t in range(3) for mu in range(4)]
        for y in OUTCOMES
    }

    def state_coefficients(z: int) -> LinearCoefficients:
        result: LinearCoefficients = {0: [(1.0, state_scalar[z])]}
        for axis in range(3):
            variable_id = state_vector[z, axis]
            result[axis + 1] = [] if variable_id is None else [(1.0, variable_id)]
        return result

    def effect_coefficients(y: int, t: int) -> LinearCoefficients:
        return {mu: [(1.0, effect_ids[y, t, mu])] for mu in range(4)}

    def residual_coefficients(y: int, t: int) -> LinearCoefficients:
        return {
            mu: [
                *((weights[t], effect_ids[y, other, mu]) for other in range(3)),
                (-1.0, effect_ids[y, t, mu]),
            ]
            for mu in range(4)
        }

    def selected_coefficients(kind: str, y: int, t: int) -> LinearCoefficients:
        return effect_coefficients(y, t) if kind == "b" else residual_coefficients(y, t)

    clique_ids = {y: [*state_ids, *effect_groups[y]] for y in OUTCOMES}
    pure_state_degree2 = monomials(state_ids, 2)
    selected_effect_ids: dict[int, set[int]] = {y: set() for y in OUTCOMES}
    for kind, y, t in linked:
        selected_effect_ids[y].update(
            coefficient_ids(selected_coefficients(kind, y, t))
        )
    if basis_mode == "anisotropic":
        bases = {
            y: [
                *pure_state_degree2,
                *((variable_id,) for variable_id in effect_groups[y]),
            ]
            for y in OUTCOMES
        }
    elif basis_mode == "targeted":
        bases = {
            y: [
                *pure_state_degree2,
                *((variable_id,) for variable_id in effect_groups[y]),
                *(
                    (state_id, effect_id)
                    for state_id in state_ids
                    for effect_id in sorted(selected_effect_ids[y])
                ),
            ]
            for y in OUTCOMES
        }
    elif basis_mode == "full":
        bases = {y: monomials(clique_ids[y], 2) for y in OUTCOMES}
    else:
        raise ValueError("basis_mode must be anisotropic, targeted, or full")
    registry = MomentRegistry()
    for y in OUTCOMES:
        registry.add_matrix(bases[y])
    bridge_basis: list[Monomial] | None = None
    if bridge_selected:
        bridge_effect_ids = sorted(
            set().union(*(selected_effect_ids[y] for y in OUTCOMES))
        )
        selected_branches = sum(bool(selected_effect_ids[y]) for y in OUTCOMES)
        if selected_branches < 2:
            raise ValueError("a selected bridge requires columns from multiple y branches")
        if bridge_basis_mode == "minimal":
            bridge_basis = [
                *pure_state_degree2,
                *((variable_id,) for variable_id in bridge_effect_ids),
            ]
        elif bridge_basis_mode == "cross":
            bridge_basis = [
                *pure_state_degree2,
                *((variable_id,) for variable_id in bridge_effect_ids),
                *(
                    (state_id, effect_id)
                    for state_id in state_ids
                    for effect_id in bridge_effect_ids
                ),
            ]
        elif bridge_basis_mode == "full":
            bridge_basis = monomials([*state_ids, *bridge_effect_ids], 2)
        else:
            raise ValueError("bridge_basis_mode must be minimal, cross, or full")
        registry.add_matrix(bridge_basis)

    # Degree-one localisers are the matrix-polynomial part of the relaxation.
    # In anisotropic mode, effect localisers use all degree-one state
    # multipliers; their entries x_i*x_j*g_k occur in the clique matrix.
    # Pure-state localisers are shared and need be imposed only once.
    if localiser_scope not in {"selected", "full"}:
        raise ValueError("localiser_scope must be selected or full")
    localiser_specs: list[tuple[list[Monomial], LinearCoefficients]] = []
    if basis_mode in {"anisotropic", "targeted"}:
        state_multipliers = monomials(state_ids, 1)
        for z in OUTCOMES:
            localiser_specs.append((state_multipliers, state_coefficients(z)))
        for y in OUTCOMES:
            for t in range(3):
                localiser_specs.extend(
                    (
                        (state_multipliers, effect_coefficients(y, t)),
                        (state_multipliers, residual_coefficients(y, t)),
                    )
                )
    else:
        for y in OUTCOMES:
            if localiser_scope == "full":
                state_multipliers = monomials(clique_ids[y], 1)
                effect_multipliers = state_multipliers
            else:
                state_multipliers = monomials(effect_groups[y], 1)
                effect_multipliers = monomials(state_ids, 1)
            for z in OUTCOMES:
                localiser_specs.append((state_multipliers, state_coefficients(z)))
            for t in range(3):
                localiser_specs.extend(
                    (
                        (effect_multipliers, effect_coefficients(y, t)),
                        (effect_multipliers, residual_coefficients(y, t)),
                    )
                )
    for multipliers, coefficients in localiser_specs:
        add_qubit_localiser_registry(registry, multipliers, coefficients)

    if ideal_degree < 0:
        raise ValueError("ideal_degree must be nonnegative")
    if basis_mode in {"anisotropic", "targeted"}:
        state_ideal_multipliers = {
            y: monomials(effect_groups[y], min(ideal_degree, 1))
            for y in OUTCOMES
        }
        effect_ideal_multipliers = monomials(
            state_ids, min(ideal_degree, 2)
        )
    else:
        state_ideal_multipliers = {
            y: monomials(effect_groups[y], min(ideal_degree, 3))
            for y in OUTCOMES
        }
        effect_ideal_multipliers = monomials(
            state_ids, min(ideal_degree, 3)
        )

    state_trace_ids = [state_scalar[z] for z in OUTCOMES]
    for y in OUTCOMES:
        register_linear_ideal(
            registry, state_ideal_multipliers[y], state_trace_ids
        )
    for mu in range(4):
        register_linear_ideal(
            registry,
            effect_ideal_multipliers,
            [effect_ids[y, t, mu] for y in OUTCOMES for t in range(3)],
        )

    # First and second moments used by scalar coordinate bounds.
    variable_bounds: dict[int, tuple[float, float]] = {}
    trace_upper_by_rank = (1.0, 0.5, 1.0 / 3.0, 0.25)
    trace_upper = {
        prefix_order[rank]: trace_upper_by_rank[rank] for rank in range(4)
    }
    for z in OUTCOMES:
        variable_bounds[state_scalar[z]] = (0.0, trace_upper[z])
        for axis in range(3):
            variable_id = state_vector[z, axis]
            if variable_id is not None:
                lower = 0.0 if (z, axis) in {(0, 0), (1, 1)} else -trace_upper[z]
                variable_bounds[variable_id] = (lower, trace_upper[z])
    for y in OUTCOMES:
        for t in range(3):
            # G[y,t] <= w_t I implies coefficient bounds |g_mu| <= 2 w_t.
            variable_bounds[effect_ids[y, t, 0]] = (0.0, 2.0 * weights[t])
            for mu in range(1, 4):
                variable_bounds[effect_ids[y, t, mu]] = (
                    -2.0 * weights[t],
                    2.0 * weights[t],
                )
    for variable_id in variable_bounds:
        registry.add((variable_id,))
        registry.add((variable_id, variable_id))

    registry.finalise()
    constraints: list[cp.Constraint] = [registry.moment(()) == 1.0]
    clique_matrices: list[cp.Expression] = []
    for y in OUTCOMES:
        matrix = moment_matrix(registry, bases[y])
        clique_matrices.append(matrix)
        constraints.append(matrix >> 0)
    if bridge_basis is not None:
        bridge_matrix = moment_matrix(registry, bridge_basis)
        clique_matrices.append(bridge_matrix)
        constraints.append(bridge_matrix >> 0)

    # Archimedean scalar bounds and their first RLT secant.
    for variable_id, (lower, upper) in variable_bounds.items():
        first = registry.moment((variable_id,))
        second = registry.moment((variable_id, variable_id))
        constraints.extend((first >= lower, first <= upper))
        constraints.append(second <= (lower + upper) * first - lower * upper)

    # Exact trace/completeness ideals at every supported multiplier degree.
    for y in OUTCOMES:
        add_linear_ideal(
            registry,
            constraints,
            state_ideal_multipliers[y],
            [(1.0, variable_id) for variable_id in state_trace_ids],
            1.0,
        )
    for mu in range(4):
        add_linear_ideal(
            registry,
            constraints,
            effect_ideal_multipliers,
            [
                (1.0, effect_ids[y, t, mu])
                for y in OUTCOMES
                for t in range(3)
            ],
            2.0 if mu == 0 else 0.0,
        )

    # Matrix-valued positivity localisers.
    localiser_matrices: list[cp.Expression] = []
    for multipliers, coefficients in localiser_specs:
        localiser = qubit_localiser(registry, multipliers, coefficients)
        localiser_matrices.append(localiser)
        constraints.append(localiser >> 0)

    # Path/terminal statistics.  Only requested columns are forced to equal
    # Born moments; unselected entries make this a monotone upper relaxation.
    statistics = cp.Variable((4, 4, 3), nonneg=True)
    probability = cp.sum(statistics, axis=2)
    constraints.append(cp.sum(statistics) == 1.0)
    constraints.extend(statistics[:, :, t] <= weights[t] * probability for t in range(3))
    prefix = [registry.moment((state_scalar[z],)) for z in OUTCOMES]
    constraints.extend(cp.sum(probability[z, :]) == prefix[z] for z in OUTCOMES)
    constraints.extend(
        prefix[prefix_order[index]] >= prefix[prefix_order[index + 1]]
        for index in range(3)
    )

    def coefficient_moment_product(
        z: int, coefficients: LinearCoefficients
    ) -> cp.Expression:
        coordinate_ids: list[int | None] = [
            state_scalar[z],
            state_vector[z, 0],
            state_vector[z, 1],
            state_vector[z, 2],
        ]
        terms: list[cp.Expression] = []
        for mu, state_id in enumerate(coordinate_ids):
            if state_id is None:
                continue
            terms.extend(
                coefficient
                * registry.moment(multiply((state_id,), (effect_id,)))
                for coefficient, effect_id in coefficients[mu]
            )
        return 0.5 * sum(terms)

    for kind, y, t in linked:
        coefficients = selected_coefficients(kind, y, t)
        for z in OUTCOMES:
            target = (
                statistics[z, y, t]
                if kind == "b"
                else weights[t] * probability[z, y] - statistics[z, y, t]
            )
            constraints.append(target == coefficient_moment_product(z, coefficients))

    terminal_statistics = [
        [
            sum(statistics[z, y, t] for z, y in PATHS if (z ^ y) == syndrome)
            for t in range(3)
        ]
        for syndrome in OUTCOMES
    ]
    inverse = reconstruction_matrix(terminal)
    terminal_prior: list[cp.Expression] = []
    terminal_vector: list[cp.Expression] = []
    terminal_normal = cp.Variable(4)
    for syndrome in OUTCOMES:
        reconstructed = inverse @ cp.hstack(terminal_statistics[syndrome])
        terminal_prior.append(reconstructed[0])
        vector = cp.hstack(
            [reconstructed[1], reconstructed[2], terminal_normal[syndrome]]
        )
        terminal_vector.append(vector)
        constraints.append(cp.SOC(reconstructed[0], vector))

    audit = sum(terminal_statistics[s][s] for s in range(3))
    cap = filled_effect_weights(float(weights.max()))
    constraints.append(
        audit
        <= sum(cap[index] * prefix[prefix_order[index]] for index in OUTCOMES)
    )
    dual_scalar = cp.Variable(nonneg=True)
    dual_vector = cp.Variable(3)
    if helstrom_form == "dual":
        constraints.append(cp.SOC(dual_scalar, dual_vector))
        constraints.extend(
            cp.SOC(
                dual_scalar - terminal_prior[s],
                dual_vector - terminal_vector[s],
            )
            for s in OUTCOMES
        )
    elif helstrom_form == "facial":
        # For an active rank-one effect E_s=w_s Pi_s, equality of the
        # Helstrom primal and dual objectives forces
        #
        #   Y - tau_s = kappa_s (I - Pi_s),  kappa_s >= 0.
        #
        # Parameterising this exposed face removes three singular Lorentz
        # constraints without changing the feasible physical set.  The zero
        # fourth effect contributes only the ordinary dual inequality.
        slack = cp.Variable(3, nonneg=True)
        for s in range(3):
            constraints.extend(
                (
                    dual_scalar - terminal_prior[s] == slack[s],
                    dual_vector - terminal_vector[s]
                    == -slack[s] * directions[s],
                )
            )
        constraints.append(
            cp.SOC(
                dual_scalar - terminal_prior[3],
                dual_vector - terminal_vector[3],
            )
        )
    else:
        raise ValueError("helstrom_form must be dual or facial")
    constraints.append(audit == dual_scalar)

    returned = hellinger_hypograph(
        [probability[z, y] for z, y in PATHS], constraints
    )
    objective = support_weight * audit + (1.0 - support_weight) * returned
    problem = cp.Problem(cp.Maximize(objective), constraints)
    if solver == "clarabel":
        problem.solve(
            solver="CLARABEL",
            tol_gap_abs=2e-7,
            tol_gap_rel=2e-7,
            tol_feas=2e-7,
            max_iter=500,
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
        raise RuntimeError(f"clique moment relaxation failed: {problem.status}")

    point = np.asarray(probability.value, dtype=float)
    terminal_prior_value = np.asarray([float(item.value) for item in terminal_prior])
    terminal_vector_value = np.asarray(
        [np.asarray(item.value, dtype=float) for item in terminal_vector]
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
    audit_value = float(audit.value)
    return_value = float(np.sqrt(np.maximum(point, 0.0)).sum() ** 2 / 16.0)
    clique_minimum = min(
        float(np.linalg.eigvalsh(np.asarray(matrix.value, dtype=float)).min())
        for matrix in clique_matrices
    )
    localiser_minimum = min(
        float(np.linalg.eigvalsh(np.asarray(matrix.value, dtype=complex)).min())
        for matrix in localiser_matrices
    )
    geometry = discrimination_geometry(terminal_states)
    return {
        "weight": support_weight,
        "terminal_effect_weights": weights.tolist(),
        "prefix_order": list(prefix_order),
        "linked_columns": [
            f"b_{3 * y + t}" if kind == "b" else f"d_{y}_{t}"
            for kind, y, t in linked
        ],
        "basis_mode": basis_mode,
        "localiser_scope": localiser_scope,
        "ideal_degree": ideal_degree,
        "helstrom_form": helstrom_form,
        "bridge_selected": bridge_selected,
        "bridge_basis_mode": bridge_basis_mode if bridge_selected else None,
        "solver": solver,
        "state_coordinate_count": len(state_ids),
        "clique_variable_count": len(next(iter(clique_ids.values()))),
        "clique_moment_sizes": [len(bases[y]) for y in OUTCOMES],
        "clique_moment_size": max(len(bases[y]) for y in OUTCOMES),
        "bridge_moment_size": None if bridge_basis is None else len(bridge_basis),
        "scalar_moment_count": len(registry.index),
        "bound": float(problem.value),
        "objective_from_reported": support_weight * audit_value
        + (1.0 - support_weight) * return_value,
        "audit": audit_value,
        "return": return_value,
        "normalisation": float(point.sum()),
        "prefix_priors": [float(item.value) for item in prefix],
        "path_probabilities": point.tolist(),
        "path_terminal_statistics": np.asarray(statistics.value, dtype=float).tolist(),
        "syndrome_priors": terminal_prior_value.tolist(),
        "terminal_bloch_vectors": terminal_vector_value.tolist(),
        "independent_terminal_geometry": geometry,
        "helstrom_independent_residual": float(
            audit_value - float(geometry["optimal_guess_probability"])
        ),
        "clique_moment_min_eigenvalue": clique_minimum,
        "localiser_min_eigenvalue": localiser_minimum,
        "status": problem.status,
        "solver_stats": {
            "solve_time": problem.solver_stats.solve_time,
            "num_iters": problem.solver_stats.num_iters,
        },
        "scope": (
            "global order-two upper relaxation of the selected-column qubit "
            "model; conic tolerance and dual-witness audit remain"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lambda", dest="weight", type=float, default=0.6)
    parser.add_argument("--fixed-three-povm-weights", type=float, nargs=3, required=True)
    parser.add_argument("--prefix-order", type=int, nargs=4, required=True)
    parser.add_argument("--linked-column", action="append", required=True)
    parser.add_argument(
        "--basis-mode",
        choices=("anisotropic", "targeted", "full"),
        default="anisotropic",
    )
    parser.add_argument(
        "--localiser-scope", choices=("selected", "full"), default="selected"
    )
    parser.add_argument("--ideal-degree", type=int, default=2)
    parser.add_argument(
        "--helstrom-form", choices=("dual", "facial"), default="facial"
    )
    parser.add_argument("--bridge-selected", action="store_true")
    parser.add_argument(
        "--bridge-basis",
        choices=("minimal", "cross", "full"),
        default="minimal",
        help=(
            "moment basis used by the selected-column bridge; cross adds all "
            "state--effect rows without enlarging the four branch cliques"
        ),
    )
    parser.add_argument("--solver", choices=("clarabel", "scs"), default="clarabel")
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--scs-eps", type=float, default=2e-5)
    parser.add_argument("--scs-iters", type=int, default=200_000)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    order = tuple(args.prefix_order)
    if sorted(order) != list(OUTCOMES):
        raise ValueError("prefix order must be a permutation of 0,1,2,3")
    linked = tuple(parse_column(name) for name in dict.fromkeys(args.linked_column))
    payload = solve_povm(
        canonical_three_effect_povm(np.asarray(args.fixed_three_povm_weights)),
        args.weight,
        order,
        linked,
        args.basis_mode,
        args.localiser_scope,
        args.ideal_degree,
        args.helstrom_form,
        args.bridge_selected,
        args.bridge_basis,
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

"""Facially reduced Choi--Helstrom moment upper bound.

The earlier first-order Choi moment model imposed normalisation and
trace-preservation together with all of their RLT multiples.  Those equations
place the moment matrix on a large exposed face and make first-order conic
solvers unnecessarily fragile.  This implementation eliminates five logical
coordinates instead:

* ``Tr(rho_3) = 1 - sum_{z<3} Tr(rho_z)``; and
* the four output-trace Pauli coefficients of ``J_3`` are fixed by
  ``sum_y Tr_out(J_y) = I``.

All products involving an eliminated coordinate are expanded affinely in the
same moment matrix.  Thus the reduction loses no rank-one physical point and
every returned optimum remains a global upper bound for the fixed terminal
POVM and prefix-prior ordering.  The terminal POVM is required to be exactly
Helstrom optimal by primal feasibility, dual feasibility, and equality of the
two objective values; explicit complementary-slackness equations would be
redundant and numerically harmful.

This is still only the first Shor level.  A value above the physical lower
bound is a relaxation gap, not evidence for a better strategy.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import cvxpy as cp
import numpy as np

from analyze_two_block_leaf import discrimination_geometry
from moment_helstrom_upper import Registry
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
TRANSPOSE_SIGN = np.asarray([1.0, 1.0, -1.0, 1.0])


@dataclass(frozen=True)
class Affine:
    """A sparse affine function of the retained scalar coordinates."""

    constant: float
    coefficients: tuple[tuple[int, float], ...]


def coordinate(index: int) -> Affine:
    return Affine(0.0, ((index, 1.0),))


def affine_sum(terms: list[tuple[float, Affine]], constant: float = 0.0) -> Affine:
    accumulated: dict[int, float] = {}
    total_constant = float(constant)
    for multiplier, form in terms:
        total_constant += multiplier * form.constant
        for index, coefficient in form.coefficients:
            accumulated[index] = accumulated.get(index, 0.0) + multiplier * coefficient
    return Affine(
        total_constant,
        tuple(
            (index, coefficient)
            for index, coefficient in sorted(accumulated.items())
            if abs(coefficient) > 1e-15
        ),
    )


def value(form: Affine, first: cp.Expression) -> cp.Expression:
    return form.constant + sum(
        coefficient * first[index] for index, coefficient in form.coefficients
    )


def product(
    left: Affine,
    right: Affine,
    first: cp.Expression,
    second: cp.Expression,
) -> cp.Expression:
    expression: cp.Expression = left.constant * right.constant
    expression += left.constant * sum(
        coefficient * first[index] for index, coefficient in right.coefficients
    )
    expression += right.constant * sum(
        coefficient * first[index] for index, coefficient in left.coefficients
    )
    expression += sum(
        left_coefficient * right_coefficient * second[left_index, right_index]
        for left_index, left_coefficient in left.coefficients
        for right_index, right_coefficient in right.coefficients
    )
    return expression


def add_logical_bounds(
    form: Affine,
    lower: float,
    upper: float,
    first: cp.Expression,
    second: cp.Expression,
    constraints: list[cp.Constraint],
) -> None:
    expression = value(form, first)
    constraints.extend((expression >= lower, expression <= upper))
    constraints.append(
        product(form, form, first, second)
        <= (lower + upper) * expression - lower * upper
    )


def add_mccormick_forms(
    left: Affine,
    right: Affine,
    lifted: cp.Expression,
    left_bounds: tuple[float, float],
    right_bounds: tuple[float, float],
    first: cp.Expression,
    constraints: list[cp.Constraint],
) -> None:
    left_value = value(left, first)
    right_value = value(right, first)
    left_lower, left_upper = left_bounds
    right_lower, right_upper = right_bounds
    constraints.extend(
        (
            lifted
            >= left_lower * right_value
            + right_lower * left_value
            - left_lower * right_lower,
            lifted
            >= left_upper * right_value
            + right_upper * left_value
            - left_upper * right_upper,
            lifted
            <= left_upper * right_value
            + right_lower * left_value
            - left_upper * right_lower,
            lifted
            <= left_lower * right_value
            + right_upper * left_value
            - left_lower * right_upper,
        )
    )


def solve_povm(
    effects: np.ndarray,
    weight: float,
    prefix_order: tuple[int, int, int, int],
    mccormick: str,
    data_processing: str,
    data_processing_scales: tuple[float, ...],
    witness_cut: dict[str, object] | None,
    input_trace_radii: tuple[float, ...] | None,
    solver: str,
    verbose: bool,
) -> dict[str, object]:
    traces = np.trace(effects, axis1=1, axis2=2).real
    active = tuple(int(s) for s in OUTCOMES if traces[s] > 1e-9)
    directions = np.zeros((4, 3), dtype=float)
    for s in active:
        projector = effects[s] / traces[s]
        directions[s] = [
            float(np.trace(projector @ pauli).real) for pauli in PAULIS[1:]
        ]

    registry = Registry()
    state: dict[tuple[int, int], Affine] = {}
    for z in range(3):
        state[z, 0] = coordinate(registry.add(f"rho_{z}_0", 0.0, 1.0))
    state[3, 0] = affine_sum(
        [(-1.0, state[z, 0]) for z in range(3)], constant=1.0
    )
    for z in OUTCOMES:
        for mu in range(1, 4):
            state[z, mu] = coordinate(
                registry.add(f"rho_{z}_{mu}", -1.0, 1.0)
            )

    choi: dict[tuple[int, int, int], Affine] = {}
    for y in range(3):
        for mu in range(4):
            for nu in range(4):
                bounds = (0.0, 2.0) if (mu, nu) == (0, 0) else (-2.0, 2.0)
                choi[y, mu, nu] = coordinate(
                    registry.add(f"j_{y}_{mu}_{nu}", *bounds)
                )
    for mu in range(4):
        target = 2.0 if mu == 0 else 0.0
        choi[3, mu, 0] = affine_sum(
            [(-1.0, choi[y, mu, 0]) for y in range(3)], constant=target
        )
    for mu in range(4):
        for nu in range(1, 4):
            choi[3, mu, nu] = coordinate(
                registry.add(f"j_3_{mu}_{nu}", -2.0, 2.0)
            )

    size = len(registry.names)
    moment = cp.Variable((size + 1, size + 1), symmetric=True)
    first = moment[0, 1:]
    second = moment[1:, 1:]
    constraints: list[cp.Constraint] = [moment >> 0, moment[0, 0] == 1.0]
    for index in range(size):
        lower = registry.lower[index]
        upper = registry.upper[index]
        constraints.extend((first[index] >= lower, first[index] <= upper))
        constraints.append(
            second[index, index]
            <= (lower + upper) * first[index] - lower * upper
        )

    # Bounds and quadratic interval localisers for the five eliminated
    # logical coordinates.
    add_logical_bounds(state[3, 0], 0.0, 1.0, first, second, constraints)
    for mu in range(4):
        bounds = (0.0, 2.0) if mu == 0 else (-2.0, 2.0)
        add_logical_bounds(choi[3, mu, 0], *bounds, first, second, constraints)

    # Subnormalised-state positivity, both at first order and at degree two.
    for z in OUTCOMES:
        scalar = value(state[z, 0], first)
        vector = cp.hstack([value(state[z, mu], first) for mu in range(1, 4)])
        constraints.append(cp.SOC(scalar, vector))
        constraints.append(
            product(state[z, 0], state[z, 0], first, second)
            >= sum(
                product(state[z, mu], state[z, mu], first, second)
                for mu in range(1, 4)
            )
        )

    choi_expressions: list[cp.Expression] = []
    for y in OUTCOMES:
        expression = sum(
            value(choi[y, mu, nu], first) * np.kron(PAULIS[mu], PAULIS[nu])
            for mu in range(4)
            for nu in range(4)
        ) / 4.0
        choi_expressions.append(expression)
        constraints.append(expression >> 0)
        trace_form = choi[y, 0, 0]
        trace_value = value(trace_form, first)
        for mu in range(4):
            for nu in range(4):
                if (mu, nu) == (0, 0):
                    continue
                coefficient_value = value(choi[y, mu, nu], first)
                constraints.extend(
                    (-trace_value <= coefficient_value, coefficient_value <= trace_value)
                )
                constraints.append(
                    product(choi[y, mu, nu], choi[y, mu, nu], first, second)
                    <= product(trace_form, trace_form, first, second)
                )
        constraints.append(
            sum(
                product(choi[y, mu, nu], choi[y, mu, nu], first, second)
                for mu in range(4)
                for nu in range(4)
            )
            <= 4.0 * product(trace_form, trace_form, first, second)
        )

    if mccormick not in {"none", "used", "all"}:
        raise ValueError("mccormick must be none, used, or all")
    if mccormick != "none":
        for z in OUTCOMES:
            for y in OUTCOMES:
                for mu in range(4):
                    state_bounds = (0.0, 1.0) if mu == 0 else (-1.0, 1.0)
                    nus = range(4)
                    for nu in nus:
                        lifted = product(state[z, mu], choi[y, mu, nu], first, second)
                        choi_bounds = (0.0, 2.0) if (mu, nu) == (0, 0) else (-2.0, 2.0)
                        add_mccormick_forms(
                            state[z, mu],
                            choi[y, mu, nu],
                            lifted,
                            state_bounds,
                            choi_bounds,
                            first,
                            constraints,
                        )
                if mccormick == "all":
                    for mu in range(4):
                        state_bounds = (0.0, 1.0) if mu == 0 else (-1.0, 1.0)
                        for alpha in range(4):
                            if alpha == mu:
                                continue
                            for nu in range(4):
                                choi_bounds = (
                                    (0.0, 2.0)
                                    if (alpha, nu) == (0, 0)
                                    else (-2.0, 2.0)
                                )
                                add_mccormick_forms(
                                    state[z, mu],
                                    choi[y, alpha, nu],
                                    product(
                                        state[z, mu],
                                        choi[y, alpha, nu],
                                        first,
                                        second,
                                    ),
                                    state_bounds,
                                    choi_bounds,
                                    first,
                                    constraints,
                                )

    # Degree-two localisers: a_z J_y >= 0 and Tr(J_y) rho_z >= 0.
    for z in OUTCOMES:
        for y in OUTCOMES:
            scaled_choi = sum(
                product(state[z, 0], choi[y, mu, nu], first, second)
                * np.kron(PAULIS[mu], PAULIS[nu])
                for mu in range(4)
                for nu in range(4)
            ) / 4.0
            constraints.append(scaled_choi >> 0)
            scaled_state_scalar = product(
                state[z, 0], choi[y, 0, 0], first, second
            )
            scaled_state_vector = cp.hstack(
                [
                    product(state[z, mu], choi[y, 0, 0], first, second)
                    for mu in range(1, 4)
                ]
            )
            constraints.append(cp.SOC(scaled_state_scalar, scaled_state_vector))

    probabilities: dict[tuple[int, int], cp.Expression] = {}
    output_vectors: dict[tuple[int, int], cp.Expression] = {}
    for z, y in PATHS:
        probabilities[z, y] = 0.5 * sum(
            TRANSPOSE_SIGN[mu]
            * product(state[z, mu], choi[y, mu, 0], first, second)
            for mu in range(4)
        )
        output_vectors[z, y] = cp.hstack(
            [
                0.5
                * sum(
                    TRANSPOSE_SIGN[mu]
                    * product(state[z, mu], choi[y, mu, nu], first, second)
                    for mu in range(4)
                )
                for nu in range(1, 4)
            ]
        )
        constraints.append(cp.SOC(probabilities[z, y], output_vectors[z, y]))

    witness_expression: cp.Expression | None = None
    witness_upper: float | None = None
    if witness_cut is not None:
        if input_trace_radii is None or len(input_trace_radii) != 4:
            raise ValueError("a witness cut requires four input trace radii")
        if any(not np.isfinite(radius) or radius < 0.0 for radius in input_trace_radii):
            raise ValueError("input trace radii must be finite and nonnegative")
        reference_bloch = np.asarray(witness_cut["reference_bloch"], dtype=float)
        witness_bloch = np.asarray(witness_cut["witness_bloch"], dtype=float)
        lipschitz = np.asarray(witness_cut["lipschitz"], dtype=float)
        if reference_bloch.shape != (4, 4) or witness_bloch.shape != (4, 4, 4):
            raise ValueError("invalid common-instrument witness dimensions")
        if lipschitz.shape != (4,):
            raise ValueError("invalid common-instrument Lipschitz constants")
        for z in OUTCOMES:
            scalar_difference = value(state[z, 0], first) - reference_bloch[z, 0]
            vector_difference = cp.hstack(
                [
                    value(state[z, mu], first) - reference_bloch[z, mu]
                    for mu in range(1, 4)
                ]
            )
            radius = input_trace_radii[z]
            constraints.extend(
                (
                    scalar_difference <= radius,
                    scalar_difference >= -radius,
                    cp.SOC(radius, vector_difference),
                )
            )
        witness_expression = sum(
            0.5
            * (
                witness_bloch[z, y, 0] * probabilities[z, y]
                + witness_bloch[z, y, 1:] @ output_vectors[z, y]
            )
            for z, y in PATHS
        )
        witness_upper = float(witness_cut["reference_support"]) + float(
            np.dot(lipschitz, np.asarray(input_trace_radii))
        )
        constraints.append(witness_expression <= witness_upper)

    if data_processing not in {"none", "prior", "quadratic"}:
        raise ValueError("data_processing must be none, prior, or quadratic")
    if not data_processing_scales:
        raise ValueError("at least one data-processing scale is required")
    if any(not np.isfinite(scale) or scale < 0.0 for scale in data_processing_scales):
        raise ValueError("data-processing scales must be finite and nonnegative")
    if data_processing != "none":
        # The flagged map rho -> direct_sum_y Phi_y(rho) is trace preserving.
        # Trace norm is therefore contractive.  For a qubit Hermitian matrix
        # (d I + v.sigma)/2 its trace norm is max(|d|, ||v||_2).  Replacing
        # the unknown input distance by Tr(rho_z)+Tr(rho_z') is weaker but
        # keeps a globally valid convex constraint at this relaxation level.
        for first_z in OUTCOMES:
            for second_z in range(first_z + 1, 4):
                for scale in data_processing_scales:
                    output_distances: list[cp.Variable] = []
                    for y in OUTCOMES:
                        distance = cp.Variable(nonneg=True)
                        probability_difference = (
                            probabilities[first_z, y]
                            - scale * probabilities[second_z, y]
                        )
                        vector_difference = (
                            output_vectors[first_z, y]
                            - scale * output_vectors[second_z, y]
                        )
                        constraints.extend(
                            (
                                distance >= probability_difference,
                                distance >= -probability_difference,
                                cp.SOC(distance, vector_difference),
                            )
                        )
                        output_distances.append(distance)
                    total_output_distance = sum(output_distances)
                    if data_processing == "prior":
                        constraints.append(
                            total_output_distance
                            <= value(state[first_z, 0], first)
                            + scale * value(state[second_z, 0], first)
                        )
                    else:
                        # For X=(d I+v.sigma)/2,
                        # ||X||_1^2=max(d^2,||v||^2).  The sum d^2+||v||^2
                        # is therefore a valid polynomial upper bound.  Each
                        # scale is a necessary Alberti--Uhlmann-type cut for
                        # the one common flagged channel; no pair is treated
                        # as if it came from an independent instrument.
                        differences = [
                            affine_sum(
                                [
                                    (1.0, state[first_z, mu]),
                                    (-scale, state[second_z, mu]),
                                ]
                            )
                            for mu in range(4)
                        ]
                        input_distance_square_upper = sum(
                            product(item, item, first, second)
                            for item in differences
                        )
                        constraints.append(
                            cp.square(total_output_distance)
                            <= input_distance_square_upper
                        )

    constraints.append(sum(probabilities.values()) == 1.0)
    prefix = [value(state[z, 0], first) for z in OUTCOMES]
    for z in OUTCOMES:
        constraints.append(sum(probabilities[z, y] for y in OUTCOMES) == prefix[z])
    constraints.extend(
        prefix[prefix_order[index]] >= prefix[prefix_order[index + 1]]
        for index in range(3)
    )

    syndrome_priors = [
        sum(probabilities[z, z ^ s] for z in OUTCOMES) for s in OUTCOMES
    ]
    terminal_vectors = [
        sum(output_vectors[z, z ^ s] for z in OUTCOMES) for s in OUTCOMES
    ]
    for s in OUTCOMES:
        constraints.append(cp.SOC(syndrome_priors[s], terminal_vectors[s]))

    audit = sum(
        0.5
        * traces[s]
        * (syndrome_priors[s] + directions[s] @ terminal_vectors[s])
        for s in active
    )
    cap = filled_effect_weights(float(traces.max()))
    constraints.append(
        audit
        <= sum(cap[index] * prefix[prefix_order[index]] for index in OUTCOMES)
    )

    # Exact fixed-POVM Helstrom optimality.  Weak duality plus equality of
    # primal and dual values already implies complementary slackness.
    dual_trace = cp.Variable(nonneg=True)
    dual_vector = cp.Variable(3)
    constraints.append(cp.SOC(dual_trace, dual_vector))
    for s in OUTCOMES:
        constraints.append(
            cp.SOC(
                dual_trace - syndrome_priors[s],
                dual_vector - terminal_vectors[s],
            )
        )
    constraints.append(audit == dual_trace)

    returned = hellinger_hypograph(
        [probabilities[z, y] for z, y in PATHS], constraints
    )
    problem = cp.Problem(
        cp.Maximize(weight * audit + (1.0 - weight) * returned), constraints
    )
    if solver == "clarabel":
        problem.solve(
            solver="CLARABEL",
            tol_gap_abs=2e-8,
            tol_gap_rel=2e-8,
            tol_feas=2e-8,
            max_iter=1000,
            verbose=verbose,
        )
    elif solver == "scs":
        problem.solve(
            solver="SCS",
            eps=5e-6,
            max_iters=300_000,
            acceleration_lookback=20,
            verbose=verbose,
        )
    else:
        raise ValueError("solver must be clarabel or scs")
    if problem.status not in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}:
        raise RuntimeError(f"reduced Choi moment relaxation failed: {problem.status}")

    point = np.asarray(
        [[float(probabilities[z, y].value) for y in OUTCOMES] for z in OUTCOMES]
    )
    terminal_prior_value = np.asarray(
        [float(item.value) for item in syndrome_priors], dtype=float
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
    returned_value = float(np.sqrt(np.maximum(point, 0.0)).sum() ** 2 / 16.0)
    moment_value = np.asarray(moment.value)
    eigenvalues = np.linalg.eigvalsh(moment_value)
    choi_minimum = min(
        float(np.linalg.eigvalsh(np.asarray(item.value)).real.min())
        for item in choi_expressions
    )
    output_margins = terminal_prior_value - np.linalg.norm(
        terminal_vector_value, axis=1
    )
    prefix_bloch_value = np.asarray(
        [
            [float(value(state[z, mu], first).value) for mu in range(4)]
            for z in OUTCOMES
        ]
    )
    conditioned_vector_value = np.asarray(
        [
            [np.asarray(output_vectors[z, y].value, dtype=float) for y in OUTCOMES]
            for z in OUTCOMES
        ]
    )
    flagged_diagnostics = []
    for first_z in OUTCOMES:
        for second_z in range(first_z + 1, 4):
            for scale in data_processing_scales:
                input_difference = (
                    prefix_bloch_value[first_z]
                    - scale * prefix_bloch_value[second_z]
                )
                input_norm = max(
                    abs(float(input_difference[0])),
                    float(np.linalg.norm(input_difference[1:])),
                )
                output_norm = 0.0
                for y in OUTCOMES:
                    probability_difference = point[first_z, y] - scale * point[second_z, y]
                    vector_difference = (
                        conditioned_vector_value[first_z, y]
                        - scale * conditioned_vector_value[second_z, y]
                    )
                    output_norm += max(
                        abs(float(probability_difference)),
                        float(np.linalg.norm(vector_difference)),
                    )
                flagged_diagnostics.append(
                    {
                        "input_pair": [first_z, second_z],
                        "scale": scale,
                        "input_trace_norm_at_first_moment": input_norm,
                        "output_trace_norm": output_norm,
                        "first_moment_slack": input_norm - output_norm,
                    }
                )
    return {
        "weight": weight,
        "terminal_effect_weights": traces.tolist(),
        "prefix_order": list(prefix_order),
        "mccormick": mccormick,
        "data_processing": data_processing,
        "data_processing_scales": list(data_processing_scales),
        "common_instrument_witness_cut": (
            None
            if witness_cut is None
            else {
                "input_trace_radii": list(input_trace_radii or ()),
                "reference_support": float(witness_cut["reference_support"]),
                "lipschitz_constants": np.asarray(
                    witness_cut["lipschitz"], dtype=float
                ).tolist(),
                "robust_upper": witness_upper,
                "value": float(witness_expression.value),
                "slack": float(witness_upper - float(witness_expression.value)),
            }
        ),
        "solver": solver,
        "moment_size": size + 1,
        "bound": float(problem.value),
        "objective_from_reported": weight * audit_value + (1.0 - weight) * returned_value,
        "audit": audit_value,
        "return": returned_value,
        "normalisation": float(point.sum()),
        "prefix_priors": [float(item.value) for item in prefix],
        "syndrome_priors": terminal_prior_value.tolist(),
        "terminal_bloch_vectors": terminal_vector_value.tolist(),
        "terminal_state_min_lorentz_margin": float(output_margins.min()),
        "dual_trace": float(dual_trace.value),
        "dual_vector": np.asarray(dual_vector.value, dtype=float).tolist(),
        "independent_terminal_geometry": geometry,
        "helstrom_primal_dual_residual": float(audit_value - float(dual_trace.value)),
        "helstrom_independent_residual": float(
            audit_value - float(geometry["optimal_guess_probability"])
        ),
        "path_probabilities": point.tolist(),
        "prefix_bloch_coefficients": prefix_bloch_value.tolist(),
        "conditioned_output_bloch_vectors": conditioned_vector_value.tolist(),
        "flagged_trace_norm_diagnostics": flagged_diagnostics,
        "worst_first_moment_flagged_slack": min(
            item["first_moment_slack"] for item in flagged_diagnostics
        ),
        "moment_min_eigenvalue": float(eigenvalues.min()),
        "moment_numerical_rank_1e-7": int(np.count_nonzero(eigenvalues > 1e-7)),
        "choi_first_moment_min_eigenvalue": choi_minimum,
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
    parser.add_argument("--mccormick", choices=("none", "used", "all"), default="used")
    parser.add_argument(
        "--data-processing",
        choices=("none", "prior", "quadratic"),
        default="quadratic",
    )
    parser.add_argument(
        "--data-processing-scale",
        type=float,
        action="append",
        default=None,
        help="repeatable t in ||Phi(rho)-t Phi(sigma)||_1 <= ||rho-t sigma||_1",
    )
    parser.add_argument("--solver", choices=("clarabel", "scs"), default="clarabel")
    parser.add_argument(
        "--common-instrument-audit",
        type=Path,
        help="JSON made by audit_common_instrument_candidate.py",
    )
    parser.add_argument(
        "--input-trace-radius",
        type=float,
        nargs="+",
        help="one shared radius or four inputwise radii for the witness box",
    )
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    order = tuple(args.prefix_order)
    if sorted(order) != list(OUTCOMES):
        raise ValueError("prefix order must be a permutation of 0,1,2,3")
    weights = np.asarray(args.fixed_three_povm_weights, dtype=float)
    witness_cut = None
    radii = None
    if args.common_instrument_audit is not None:
        if args.input_trace_radius is None:
            raise ValueError("--common-instrument-audit requires --input-trace-radius")
        raw_radii = tuple(map(float, args.input_trace_radius))
        if len(raw_radii) == 1:
            radii = raw_radii * 4
        elif len(raw_radii) == 4:
            radii = raw_radii
        else:
            raise ValueError("provide one shared radius or four radii")
        audit_payload = json.loads(
            args.common_instrument_audit.read_text(encoding="utf-8")
        )
        arrays = np.load(args.common_instrument_audit.with_suffix(".npz"))
        reference_states = np.asarray(arrays["prefix_states"])
        witness = np.asarray(arrays["separating_witness"])
        reference_bloch = np.asarray(
            [
                [float(np.trace(state_matrix @ pauli).real) for pauli in PAULIS]
                for state_matrix in reference_states
            ]
        )
        witness_bloch = np.asarray(
            [
                [
                    [float(np.trace(witness[z, y] @ pauli).real) for pauli in PAULIS]
                    for y in OUTCOMES
                ]
                for z in OUTCOMES
            ]
        )
        projection_payload = audit_payload["choi_projection"]
        witness_cut = {
            "reference_bloch": reference_bloch,
            "witness_bloch": witness_bloch,
            "lipschitz": np.asarray(
                projection_payload["input_lipschitz_constants"], dtype=float
            ),
            "reference_support": float(
                projection_payload["compatible_support_value"]
            ),
        }
    elif args.input_trace_radius is not None:
        raise ValueError("--input-trace-radius requires --common-instrument-audit")
    payload = solve_povm(
        canonical_three_effect_povm(weights),
        args.weight,
        order,
        args.mccormick,
        args.data_processing,
        tuple(
            (0.25, 0.5, 1.0, 2.0, 4.0)
            if args.data_processing_scale is None
            else args.data_processing_scale
        ),
        witness_cut,
        radii,
        args.solver,
        args.verbose,
    )
    rendered = json.dumps(payload, indent=2) + "\n"
    print(rendered)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

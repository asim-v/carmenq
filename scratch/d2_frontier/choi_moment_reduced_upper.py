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
    return cp.Constant(form.constant) + sum(
        coefficient * first[index] for index, coefficient in form.coefficients
    )


def product(
    left: Affine,
    right: Affine,
    first: cp.Expression,
    second: cp.Expression,
) -> cp.Expression:
    expression: cp.Expression = cp.Constant(left.constant * right.constant)
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
    state_coordinate_bounds: np.ndarray | None,
    witness_cuts: tuple[dict[str, object], ...],
    solver: str,
    verbose: bool,
    tensor_localizers: str = "none",
    gauge_fix: bool = False,
    return_linear_weights: np.ndarray | None = None,
    fourier_trace_branches: tuple[str, str, str] | None = None,
    fourier_bloch_caps: tuple[
        tuple[float, float, float, float] | None,
        tuple[float, float, float, float] | None,
        tuple[float, float, float, float] | None,
    ]
    | None = None,
) -> dict[str, object]:
    traces = np.trace(effects, axis1=1, axis2=2).real
    active = tuple(int(s) for s in OUTCOMES if traces[s] > 1e-9)
    directions = np.zeros((4, 3), dtype=float)
    for s in active:
        projector = effects[s] / traces[s]
        directions[s] = [
            float(np.trace(projector @ pauli).real) for pauli in PAULIS[1:]
        ]

    if state_coordinate_bounds is None:
        state_bounds = np.zeros((4, 4, 2), dtype=float)
        state_bounds[:, 0] = (0.0, 1.0)
        state_bounds[:, 1:, 0] = -1.0
        state_bounds[:, 1:, 1] = 1.0
    else:
        state_bounds = np.asarray(state_coordinate_bounds, dtype=float)
        if state_bounds.shape != (4, 4, 2):
            raise ValueError("state-coordinate bounds must have shape (4,4,2)")
        if np.any(~np.isfinite(state_bounds)) or np.any(
            state_bounds[..., 0] > state_bounds[..., 1]
        ):
            raise ValueError("invalid state-coordinate bounds")
    if tensor_localizers not in {"none", "state-choi", "state-choi-ppt"}:
        raise ValueError(
            "tensor_localizers must be none, state-choi, or state-choi-ppt"
        )

    registry = Registry()
    state: dict[tuple[int, int], Affine] = {}
    for z in range(3):
        state[z, 0] = coordinate(
            registry.add(f"rho_{z}_0", *state_bounds[z, 0])
        )
    state[3, 0] = affine_sum(
        [(-1.0, state[z, 0]) for z in range(3)], constant=1.0
    )
    for z in OUTCOMES:
        for mu in range(1, 4):
            if gauge_fix and ((z, mu) == (0, 2) or (z == 3 and mu in (1, 2))):
                continue
            state[z, mu] = coordinate(
                registry.add(f"rho_{z}_{mu}", *state_bounds[z, mu])
            )
    if gauge_fix:
        # Simultaneous input conjugation and inverse pre-conjugation of every
        # instrument map leave all conditioned outputs unchanged.  Use this
        # SU(2) gauge to align the total input Bloch vector with +z and the
        # transverse part of rho_0 with +x.  Three coordinates can therefore
        # be eliminated before forming the moment matrix.
        state[0, 2] = Affine(0.0, ())
        state[3, 1] = affine_sum(
            [(-1.0, state[z, 1]) for z in range(3)]
        )
        state[3, 2] = affine_sum(
            [(-1.0, state[z, 2]) for z in range(3)]
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

    # Bounds and quadratic interval localisers for eliminated logical
    # coordinates.
    add_logical_bounds(
        state[3, 0], *state_bounds[3, 0], first, second, constraints
    )
    if gauge_fix:
        for z, mu in ((0, 2), (3, 1), (3, 2)):
            add_logical_bounds(
                state[z, mu], *state_bounds[z, mu], first, second, constraints
            )
        constraints.extend(
            (
                value(state[0, 1], first) >= 0.0,
                sum(value(state[z, 3], first) for z in OUTCOMES) >= 0.0,
            )
        )
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
                    coordinate_bounds = tuple(state_bounds[z, mu])
                    nus = range(4)
                    for nu in nus:
                        lifted = product(state[z, mu], choi[y, mu, nu], first, second)
                        choi_bounds = (0.0, 2.0) if (mu, nu) == (0, 0) else (-2.0, 2.0)
                        add_mccormick_forms(
                            state[z, mu],
                            choi[y, mu, nu],
                            lifted,
                            coordinate_bounds,
                            choi_bounds,
                            first,
                            constraints,
                        )
                if mccormick == "all":
                    for mu in range(4):
                        coordinate_bounds = tuple(state_bounds[z, mu])
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
                                    coordinate_bounds,
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

    if tensor_localizers in {"state-choi", "state-choi-ppt"}:
        # At every physical rank-one moment point this matrix is exactly
        # rho_z^T tensor J_y.  Positivity of the two factors therefore makes
        # the 8-by-8 matrix positive semidefinite.  The scalar localisers
        # above are only coarse marginals of this matrix-valued constraint;
        # retaining the full tensor product couples *all* state--Choi cross
        # moments and is a substantially stronger first-level condition.
        triple_paulis = {
            (mu, alpha, nu): np.kron(
                PAULIS[mu], np.kron(PAULIS[alpha], PAULIS[nu])
            )
            for mu in range(4)
            for alpha in range(4)
            for nu in range(4)
        }
        for z in OUTCOMES:
            for y in OUTCOMES:
                tensor_product = sum(
                    TRANSPOSE_SIGN[mu]
                    * product(
                        state[z, mu], choi[y, alpha, nu], first, second
                    )
                    * triple_paulis[mu, alpha, nu]
                    for mu in range(4)
                    for alpha in range(4)
                    for nu in range(4)
                ) / 8.0
                constraints.append(cp.hermitian_wrap(tensor_product) >> 0)
                if tensor_localizers == "state-choi-ppt":
                    # The physical matrix is a product across the auxiliary
                    # input-state / Choi split, hence it is separable across
                    # 2 x 4 and must remain positive under partial transpose.
                    # Removing TRANSPOSE_SIGN performs that transpose on the
                    # first qubit factor.  This PPT constraint rejects lifted
                    # state--instrument correlations that mere positivity
                    # cannot see.
                    partial_transpose = sum(
                        product(
                            state[z, mu],
                            choi[y, alpha, nu],
                            first,
                            second,
                        )
                        * triple_paulis[mu, alpha, nu]
                        for mu in range(4)
                        for alpha in range(4)
                        for nu in range(4)
                    ) / 8.0
                    constraints.append(cp.hermitian_wrap(partial_transpose) >> 0)

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

    # Exact finite-group data processing for the one common flagged channel.
    # With G=Z_2^2 and a convolutional readout,
    #
    #   tau_s = sum_z Phi_{z xor s}(rho_z),
    #
    # every nontrivial Fourier component is the image of the corresponding
    # input component under the same flagged CPTP map.  Hence
    #
    #   sum_y ||Phi_y(rho_hat_chi)||_1 <= ||rho_hat_chi||_1.
    #
    # For a Hermitian qubit operator the norm on the right is the maximum of
    # the absolute scalar coefficient and the Bloch-vector length.  The three
    # possible active pieces are covered by scalar-positive, scalar-negative,
    # and bloch.  Scalar pieces are SOC representable.  For all Bloch pieces
    # in a branch, Parseval and state positivity give the joint convex cut
    #
    #   sum_chi q_chi^2 + ||sum_z r_z||^2 <= 4 sum_z a_z^2,
    #
    # and each a_z^2 is majorised by its secant on the current prior cell.
    fourier_flagged: list[cp.Expression] = []
    if fourier_trace_branches is not None:
        allowed_fourier_branches = {
            "scalar-positive",
            "scalar-negative",
            "bloch",
        }
        if len(fourier_trace_branches) != 3 or any(
            branch not in allowed_fourier_branches
            for branch in fourier_trace_branches
        ):
            raise ValueError(
                "fourier_trace_branches must contain three spectral branches"
            )
        if fourier_bloch_caps is None:
            spectral_caps = (None, None, None)
        else:
            if len(fourier_bloch_caps) != 3:
                raise ValueError("fourier_bloch_caps must contain three entries")
            spectral_caps = fourier_bloch_caps
        characters = np.asarray(
            [
                [1.0, 1.0, -1.0, -1.0],
                [1.0, -1.0, 1.0, -1.0],
                [1.0, -1.0, -1.0, 1.0],
            ]
        )
        bloch_flagged: list[cp.Expression] = []
        bloch_seen = 0
        for character_index, (character, branch) in enumerate(
            zip(characters, fourier_trace_branches)
        ):
            scalar_form = affine_sum(
                [(character[z], state[z, 0]) for z in OUTCOMES]
            )
            scalar = value(scalar_form, first)
            vector_forms = [
                affine_sum(
                    [(character[z], state[z, mu]) for z in OUTCOMES]
                )
                for mu in range(1, 4)
            ]
            vector = cp.hstack(
                [value(vector_form, first) for vector_form in vector_forms]
            )
            scalar_square = product(scalar_form, scalar_form, first, second)
            vector_square = sum(
                product(vector_form, vector_form, first, second)
                for vector_form in vector_forms
            )
            block_norms: list[cp.Variable] = []
            for y in OUTCOMES:
                output_scalar = sum(
                    character[z] * probabilities[z, y] for z in OUTCOMES
                )
                output_vector = sum(
                    character[z] * output_vectors[z, y] for z in OUTCOMES
                )
                block_norm = cp.Variable(nonneg=True)
                constraints.extend(
                    (
                        block_norm >= output_scalar,
                        block_norm >= -output_scalar,
                        cp.SOC(block_norm, output_vector),
                    )
                )
                block_norms.append(block_norm)
            flagged = sum(block_norms)
            fourier_flagged.append(flagged)
            if branch == "scalar-positive":
                constraints.extend(
                    (
                        cp.SOC(scalar, vector),
                        flagged <= scalar,
                        scalar_square >= vector_square,
                    )
                )
            elif branch == "scalar-negative":
                constraints.extend(
                    (
                        cp.SOC(-scalar, vector),
                        flagged <= -scalar,
                        scalar_square >= vector_square,
                    )
                )
            else:
                constraints.extend(
                    (
                        vector_square >= scalar_square,
                        cp.square(flagged) <= vector_square,
                    )
                )
                if bloch_seen == 0:
                    # The input side has a free simultaneous SU(2) gauge:
                    # rotate every prefix state and precompose the common
                    # instrument with the inverse rotation.  Outputs and the
                    # objective are unchanged.  Align the first vector-active
                    # Fourier component with +z; its trace-norm contraction
                    # is then the exact linear inequality below.
                    constraints.extend(
                        (
                            vector[0] == 0.0,
                            vector[1] == 0.0,
                            vector[2] >= 0.0,
                            flagged <= vector[2],
                        )
                    )
                elif bloch_seen == 1:
                    # The residual rotation around z can put the transverse
                    # component of the second active vector on the +x axis.
                    constraints.extend((vector[1] == 0.0, vector[0] >= 0.0))
                cap = spectral_caps[character_index]
                if cap is not None:
                    cap_array = np.asarray(cap, dtype=float)
                    if cap_array.shape != (4,):
                        raise ValueError("a Fourier Bloch cap is (nx,ny,nz,cosine)")
                    normal = cap_array[:3]
                    cosine = float(cap_array[3])
                    if (
                        not np.all(np.isfinite(cap_array))
                        or abs(float(np.linalg.norm(normal)) - 1.0) > 1e-9
                        or not 0.0 < cosine <= 1.0
                    ):
                        raise ValueError("invalid Fourier Bloch cap")
                    projection = normal @ vector
                    # The cap condition is n.v >= cos(delta)||v||.  Inside
                    # the cap, ||v|| <= n.v/cos(delta), yielding a linear
                    # outer bound on the otherwise reverse-convex trace-norm
                    # contraction.  A finite cap cover therefore gives a
                    # finite family of rigorous conic outer problems.
                    constraints.extend(
                        (
                            cp.SOC(projection / cosine, vector),
                            flagged <= projection / cosine,
                        )
                    )
                bloch_seen += 1
                bloch_flagged.append(flagged)

        if bloch_flagged:
            total_input_vector = sum(
                cp.hstack(
                    [value(state[z, mu], first) for mu in range(1, 4)]
                )
                for z in OUTCOMES
            )
            prior_square_secants = sum(
                (state_bounds[z, 0, 0] + state_bounds[z, 0, 1])
                * value(state[z, 0], first)
                - state_bounds[z, 0, 0] * state_bounds[z, 0, 1]
                for z in OUTCOMES
            )
            constraints.append(
                cp.sum_squares(cp.hstack(bloch_flagged))
                + cp.sum_squares(total_input_vector)
                <= 4.0 * prior_square_secants
            )

    witness_expressions: list[cp.Expression] = []
    witness_uppers: list[float] = []
    for witness_cut in witness_cuts:
        input_trace_radii = np.asarray(witness_cut["input_trace_radii"], dtype=float)
        if input_trace_radii.shape != (4,):
            raise ValueError("a witness cut requires four input trace radii")
        if np.any(~np.isfinite(input_trace_radii)) or np.any(input_trace_radii < 0.0):
            raise ValueError("input trace radii must be finite and nonnegative")
        reference_bloch = np.asarray(witness_cut["reference_bloch"], dtype=float)
        witness_bloch = np.asarray(witness_cut["witness_bloch"], dtype=float)
        lipschitz = np.asarray(witness_cut["lipschitz"], dtype=float)
        if reference_bloch.shape != (4, 4) or witness_bloch.shape != (4, 4, 4):
            raise ValueError("invalid common-instrument witness dimensions")
        if lipschitz.shape != (4,):
            raise ValueError("invalid common-instrument Lipschitz constants")
        if bool(witness_cut.get("restrict_to_balls", False)):
            for z in OUTCOMES:
                scalar_difference = (
                    value(state[z, 0], first) - reference_bloch[z, 0]
                )
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
            np.dot(lipschitz, input_trace_radii)
        )
        constraints.append(witness_expression <= witness_upper)
        witness_expressions.append(witness_expression)
        witness_uppers.append(witness_upper)

    if data_processing not in {"none", "prior", "quadratic", "cell"}:
        raise ValueError(
            "data_processing must be none, prior, quadratic, or cell"
        )
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
                    elif data_processing == "cell":
                        # On an axis-aligned state cell, interval arithmetic
                        # gives a deterministic upper bound on the exact qubit
                        # trace norm.  Unlike the quadratic moment surrogate,
                        # this bound cannot borrow artificial variance and it
                        # converges to the true input distance as the spatial
                        # branch shrinks.
                        lower = (
                            state_bounds[first_z, :, 0]
                            - scale * state_bounds[second_z, :, 1]
                        )
                        upper = (
                            state_bounds[first_z, :, 1]
                            - scale * state_bounds[second_z, :, 0]
                        )
                        absolute = np.maximum(np.abs(lower), np.abs(upper))
                        input_distance_upper = max(
                            float(absolute[0]),
                            float(np.linalg.norm(absolute[1:])),
                        )
                        constraints.extend(
                            (
                                total_output_distance <= input_distance_upper,
                                total_output_distance
                                <= value(state[first_z, 0], first)
                                + scale * value(state[second_z, 0], first),
                            )
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
    # The active fixed effects are nonzero multiples of rank-one projectors.
    # Equality of the primal and dual Helstrom objectives therefore exposes a
    # one-dimensional face of every active dual slack.  Writing
    # Pi_s=(I+n_s.sigma)/2 gives the exact affine identity
    #
    #   t_s = d + (Tr(Y)-b_s) n_s.
    #
    # This is stronger numerically than leaving complementary slackness
    # implicit in a conic solver, while remaining logically redundant at
    # every physical feasible point.
    for s in active:
        constraints.append(
            terminal_vectors[s]
            == dual_vector
            + (dual_trace - syndrome_priors[s]) * directions[s]
        )

    if return_linear_weights is None:
        returned = hellinger_hypograph(
            [probabilities[z, y] for z, y in PATHS], constraints
        )
        return_mode = "hellinger"
        linear_return_weights = None
    else:
        linear_return_weights = np.asarray(return_linear_weights, dtype=float)
        if linear_return_weights.shape != (4, 4):
            raise ValueError("return_linear_weights must have shape (4,4)")
        if np.any(~np.isfinite(linear_return_weights)) or np.any(
            linear_return_weights < 0.0
        ):
            raise ValueError("return linear weights must be finite and nonnegative")
        returned = sum(
            linear_return_weights[z, y] * probabilities[z, y]
            for z, y in PATHS
        )
        return_mode = "linear-majorant"
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
    hellinger_value = float(np.sqrt(np.maximum(point, 0.0)).sum() ** 2 / 16.0)
    returned_value = float(returned.value)
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
    witness_reports = [
        {
            "input_trace_radii": np.asarray(
                witness_cut["input_trace_radii"], dtype=float
            ).tolist(),
            "reference_support": float(witness_cut["reference_support"]),
            "lipschitz_constants": np.asarray(
                witness_cut["lipschitz"], dtype=float
            ).tolist(),
            "restrict_to_balls": bool(
                witness_cut.get("restrict_to_balls", False)
            ),
            "robust_upper": witness_uppers[index],
            "value": float(witness_expressions[index].value),
            "slack": float(
                witness_uppers[index] - float(witness_expressions[index].value)
            ),
        }
        for index, witness_cut in enumerate(witness_cuts)
    ]
    return {
        "weight": weight,
        "terminal_effect_weights": traces.tolist(),
        "prefix_order": list(prefix_order),
        "mccormick": mccormick,
        "data_processing": data_processing,
        "data_processing_scales": list(data_processing_scales),
        "tensor_localizers": tensor_localizers,
        "gauge_fix": gauge_fix,
        "return_mode": return_mode,
        "return_linear_weights": (
            None if linear_return_weights is None else linear_return_weights.tolist()
        ),
        "fourier_trace_branches": (
            None
            if fourier_trace_branches is None
            else list(fourier_trace_branches)
        ),
        "fourier_flagged_norms": (
            None
            if fourier_trace_branches is None
            else [float(item.value) for item in fourier_flagged]
        ),
        "fourier_bloch_caps": (
            None
            if fourier_bloch_caps is None
            else [None if cap is None else list(cap) for cap in fourier_bloch_caps]
        ),
        "state_coordinate_bounds": state_bounds.tolist(),
        "common_instrument_witness_cuts": witness_reports,
        # Backward-compatible singular field used by the first local audit.
        "common_instrument_witness_cut": (
            witness_reports[0] if len(witness_reports) == 1 else None
        ),
        "solver": solver,
        "moment_size": size + 1,
        "bound": float(problem.value),
        "objective_from_reported": weight * audit_value + (1.0 - weight) * returned_value,
        "audit": audit_value,
        "return": returned_value,
        "hellinger_return_at_reported": hellinger_value,
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
        choices=("none", "prior", "quadratic", "cell"),
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
        "--tensor-localizers",
        choices=("none", "state-choi", "state-choi-ppt"),
        default="none",
        help="matrix-valued degree-two positivity localizers",
    )
    parser.add_argument(
        "--gauge-fix",
        action="store_true",
        help="quotient the global SU(2) input gauge before relaxation",
    )
    parser.add_argument(
        "--return-tangent-from",
        type=Path,
        help="use the global Hellinger tangent at path probabilities in this JSON",
    )
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
    witness_cuts: tuple[dict[str, object], ...] = ()
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
        witness_cuts = (
            {
                "reference_bloch": reference_bloch,
                "witness_bloch": witness_bloch,
                "lipschitz": np.asarray(
                    projection_payload["input_lipschitz_constants"], dtype=float
                ),
                "reference_support": float(
                    projection_payload["compatible_support_value"]
                ),
                "input_trace_radii": np.asarray(radii, dtype=float),
                "restrict_to_balls": True,
            },
        )
    elif args.input_trace_radius is not None:
        raise ValueError("--input-trace-radius requires --common-instrument-audit")
    return_linear_weights = None
    if args.return_tangent_from is not None:
        tangent_payload = json.loads(
            args.return_tangent_from.read_text(encoding="utf-8")
        )
        tangent_point = np.asarray(
            tangent_payload["path_probabilities"], dtype=float
        )
        if tangent_point.shape != (4, 4) or np.any(tangent_point <= 0.0):
            raise ValueError("return tangent requires sixteen positive probabilities")
        root_sum = float(np.sqrt(tangent_point).sum())
        return_linear_weights = root_sum / (16.0 * np.sqrt(tangent_point))

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
        None,
        witness_cuts,
        args.solver,
        args.verbose,
        args.tensor_localizers,
        args.gauge_fix,
        return_linear_weights,
    )
    rendered = json.dumps(payload, indent=2) + "\n"
    print(rendered)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

"""Spatial branch-and-bound for the joint-effect/Helstrom outer model.

All qubit states and effects are represented by Bloch four-vectors.  Their
positive-semidefinite constraints are Lorentz-cone quadratic inequalities;
the only nonconvex equalities are the physical state/effect products that
produce path and AUDIT probabilities.  SCIP therefore supplies a global
dual bound for a *fixed* terminal rank-one POVM, subject to its numerical
tolerances.

By default this is still an outer model: ``G[y,s] <= w_s Q[y]`` is necessary
but not sufficient for all ``G[y,s]`` to be pullbacks of the same terminal
POVM by a common output instrument.  The optional complete-positive
completion reconstructs the three planar Pauli pullbacks and requires one
literal positive Choi matrix per outcome.  With that option the pulled-effect
description is exactly equivalent to a shared qubit instrument for the score.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np
from pyscipopt import Model, quicksum

OUTCOMES = range(4)
PATHS = tuple((z, y) for z in OUTCOMES for y in OUTCOMES)
IDENTITY = np.eye(2, dtype=complex)
PAULIS = (
    np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex),
    np.array([[0.0, -1j], [1j, 0.0]], dtype=complex),
    np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex),
)
FULL_PAULIS = (IDENTITY, *PAULIS)
CHOI_BASIS = np.asarray(
    [[np.kron(left, right) for right in FULL_PAULIS] for left in FULL_PAULIS]
)


def canonical_three_effect_povm(weights: np.ndarray) -> np.ndarray:
    """Return the planar rank-one qubit POVM fixed by three traces."""
    w0, w1, w2 = (float(value) for value in weights)
    cosine = (w2 * w2 - w0 * w0 - w1 * w1) / (2.0 * w0 * w1)
    cosine = float(np.clip(cosine, -1.0, 1.0))
    sine = math.sqrt(max(0.0, 1.0 - cosine * cosine))
    directions = np.asarray(
        [
            [1.0, 0.0, 0.0],
            [cosine, sine, 0.0],
            [-(w0 + w1 * cosine) / w2, -(w1 * sine) / w2, 0.0],
        ]
    )
    effects = np.zeros((4, 2, 2), dtype=complex)
    for s in range(3):
        effects[s] = (
            0.5
            * weights[s]
            * (IDENTITY + sum(directions[s, axis] * PAULIS[axis] for axis in range(3)))
        )
    return effects


def bloch(matrix: np.ndarray) -> np.ndarray:
    return np.asarray(
        [
            float(np.trace(matrix).real),
            *(float(np.trace(matrix @ pauli).real) for pauli in PAULIS),
        ]
    )


def rotate_common_instrument_input_gauge(
    states: np.ndarray, choi: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Rotate input states and the Choi input legs into the model gauge."""

    state_coefficients = np.asarray([bloch(item) for item in states])
    first = state_coefficients[0, 1:]
    second = state_coefficients[1, 1:]
    first_norm = float(np.linalg.norm(first))
    if first_norm <= 1e-12:
        return np.asarray(states), np.asarray(choi)
    x_axis = first / first_norm
    transverse = second - float(np.dot(second, x_axis)) * x_axis
    transverse_norm = float(np.linalg.norm(transverse))
    if transverse_norm <= 1e-12:
        trial = np.asarray([0.0, 0.0, 1.0])
        if abs(float(np.dot(trial, x_axis))) > 0.9:
            trial = np.asarray([0.0, 1.0, 0.0])
        transverse = trial - float(np.dot(trial, x_axis)) * x_axis
        transverse_norm = float(np.linalg.norm(transverse))
    y_axis = transverse / transverse_norm
    z_axis = np.cross(x_axis, y_axis)
    rotation = np.vstack([x_axis, y_axis, z_axis])

    state_coefficients[:, 1:] = state_coefficients[:, 1:] @ rotation.T
    choi_coefficients = np.asarray(
        [
            [
                [
                    float(np.trace(item @ CHOI_BASIS[mu, nu]).real)
                    for nu in OUTCOMES
                ]
                for mu in OUTCOMES
            ]
            for item in choi
        ]
    )
    transpose_sign = np.diag([1.0, -1.0, 1.0])
    choi_rotation = transpose_sign @ rotation @ transpose_sign
    choi_coefficients[:, 1:, :] = np.einsum(
        "ab,ybn->yan", choi_rotation, choi_coefficients[:, 1:, :]
    )
    rotated_states = np.asarray(
        [
            sum(state_coefficients[z, mu] * FULL_PAULIS[mu] for mu in OUTCOMES)
            / 2.0
            for z in OUTCOMES
        ]
    )
    rotated_choi = np.asarray(
        [
            sum(
                choi_coefficients[y, mu, nu] * CHOI_BASIS[mu, nu]
                for mu in OUTCOMES
                for nu in OUTCOMES
            )
            / 4.0
            for y in OUTCOMES
        ]
    )
    return rotated_states, rotated_choi


def add_lorentz(
    model: Model,
    scalar: object,
    vector: tuple[object, object, object] | list[object],
) -> None:
    model.addCons(scalar >= 0.0)
    model.addCons(scalar * scalar >= quicksum(item * item for item in vector))


def add_complex_cholesky(
    model: Model,
    coefficients: dict[tuple[int, int], object],
    label: str,
) -> dict[tuple[int, int, str], object]:
    """Impose positivity of a two-qubit Pauli expansion by ``J=L L*``."""

    root_two = math.sqrt(2.0)
    factor: dict[tuple[int, int, str], object] = {}
    for row in OUTCOMES:
        for column in range(row + 1):
            factor[row, column, "real"] = model.addVar(
                lb=0.0 if row == column else -root_two,
                ub=root_two,
                name=f"{label}_L_{row}_{column}_re",
            )
            if row != column:
                factor[row, column, "imag"] = model.addVar(
                    lb=-root_two,
                    ub=root_two,
                    name=f"{label}_L_{row}_{column}_im",
                )

    def component(row: int, column: int, part: str) -> object:
        if part == "imag" and row == column:
            return 0.0
        return factor[row, column, part]

    for row in OUTCOMES:
        for column in range(row + 1):
            target_real = quicksum(
                float(CHOI_BASIS[mu, nu, row, column].real) * coefficients[mu, nu] / 4.0
                for mu in OUTCOMES
                for nu in OUTCOMES
            )
            target_imaginary = quicksum(
                float(CHOI_BASIS[mu, nu, row, column].imag) * coefficients[mu, nu] / 4.0
                for mu in OUTCOMES
                for nu in OUTCOMES
            )
            model.addCons(
                target_real
                == quicksum(
                    component(row, inner, "real") * component(column, inner, "real")
                    + component(row, inner, "imag") * component(column, inner, "imag")
                    for inner in range(column + 1)
                )
            )
            if row != column:
                model.addCons(
                    target_imaginary
                    == quicksum(
                        component(row, inner, "imag") * component(column, inner, "real")
                        - component(row, inner, "real")
                        * component(column, inner, "imag")
                        for inner in range(column + 1)
                    )
                )
    return factor


def lower_psd_factor(matrix: np.ndarray, tolerance: float = 1e-10) -> np.ndarray:
    """Return a lower-triangular ``L`` with ``matrix = L L*``.

    Unlike ``numpy.linalg.cholesky``, this construction also accepts
    rank-deficient positive-semidefinite matrices.  The unpivoted QR step
    triangularises a square-root factor without changing its Gram matrix.
    """

    hermitian = 0.5 * (np.asarray(matrix) + np.asarray(matrix).conj().T)
    eigenvalues, eigenvectors = np.linalg.eigh(hermitian)
    scale = max(1.0, float(np.linalg.norm(hermitian, ord=2)))
    if float(eigenvalues.min()) < -tolerance * scale:
        raise np.linalg.LinAlgError("matrix is not positive semidefinite")
    square_root = eigenvectors @ np.diag(np.sqrt(np.maximum(eigenvalues, 0.0)))
    _, upper = np.linalg.qr(square_root.conj().T)
    lower = upper.conj().T
    for column in range(lower.shape[1]):
        diagonal = lower[column, column]
        if abs(diagonal) > tolerance:
            lower[:, column] *= np.exp(-1j * np.angle(diagonal))
    return lower


def polynomial_determinant(matrix: list[list[object]]) -> object:
    """Leibniz determinant for a small matrix of SCIP expressions."""

    dimension = len(matrix)
    if dimension < 1 or any(len(row) != dimension for row in matrix):
        raise ValueError("a square expression matrix is required")
    terms = []
    for permutation in itertools.permutations(range(dimension)):
        inversions = sum(
            permutation[first] > permutation[second]
            for first in range(dimension)
            for second in range(first + 1, dimension)
        )
        term: object = -1.0 if inversions % 2 else 1.0
        for row, column in enumerate(permutation):
            term = term * matrix[row][column]
        terms.append(term)
    return quicksum(terms)


def polynomial_adjugate(matrix: list[list[object]]) -> list[list[object]]:
    """Adjugate by complementary minors, with no division."""

    dimension = len(matrix)
    if dimension < 2 or any(len(row) != dimension for row in matrix):
        raise ValueError("a square expression matrix of dimension at least two is required")
    result: list[list[object]] = []
    for row in range(dimension):
        result_row = []
        for column in range(dimension):
            # adj(A)[row,column] is cofactor(A)[column,row].
            minor = [
                [
                    matrix[source_row][source_column]
                    for source_column in range(dimension)
                    if source_column != row
                ]
                for source_row in range(dimension)
                if source_row != column
            ]
            sign = -1.0 if (row + column) % 2 else 1.0
            result_row.append(sign * polynomial_determinant(minor))
        result.append(result_row)
    return result


def fibonacci_sphere(count: int) -> np.ndarray:
    """Return deterministic pure-state Bloch directions for Ando cuts."""

    if count < 1:
        return np.zeros((0, 3), dtype=float)
    golden_angle = math.pi * (3.0 - math.sqrt(5.0))
    directions = []
    for index in range(count):
        z_coordinate = 1.0 - 2.0 * (index + 0.5) / count
        radius = math.sqrt(max(0.0, 1.0 - z_coordinate * z_coordinate))
        angle = golden_angle * index
        directions.append(
            [radius * math.cos(angle), radius * math.sin(angle), z_coordinate]
        )
    return np.asarray(directions)


def build(
    effects: np.ndarray,
    weight: float,
    prefix_order: tuple[int, int, int, int] | None,
    target: float | None,
    fix_rotation_gauge: bool,
    linked_columns: tuple[str, ...] | None,
    require_cp_completion: bool = False,
    ando_direction_count: int = 0,
    fourier_trace_branches: tuple[str, str, str] | None = None,
    prefix_prior_bounds: np.ndarray | None = None,
    basis_determinant_sign: int = 0,
    basis_determinant_floor: float = 0.0,
    basis_choi_witnesses: np.ndarray | None = None,
) -> tuple[Model, dict[str, object]]:
    if ando_direction_count < 0:
        raise ValueError("ando_direction_count must be nonnegative")
    if basis_determinant_sign not in {-1, 0, 1}:
        raise ValueError("basis_determinant_sign must be -1, 0, or 1")
    if not np.isfinite(basis_determinant_floor) or basis_determinant_floor < 0.0:
        raise ValueError("basis_determinant_floor must be finite and nonnegative")
    if basis_determinant_sign and not require_cp_completion:
        raise ValueError("basis determinant completion requires CP completion")
    if basis_choi_witnesses is not None:
        basis_choi_witnesses = np.asarray(basis_choi_witnesses, dtype=complex)
        if (
            basis_choi_witnesses.ndim != 3
            or basis_choi_witnesses.shape[0] != 4
            or basis_choi_witnesses.shape[2] != 4
        ):
            raise ValueError("basis Choi witnesses must have shape (4,cuts,4)")
        if not basis_determinant_sign:
            raise ValueError("basis Choi witnesses require a determinant-sign branch")
    branch_choices = {"scalar-positive", "scalar-negative", "bloch"}
    if fourier_trace_branches is not None and (
        len(fourier_trace_branches) != 3
        or any(choice not in branch_choices for choice in fourier_trace_branches)
    ):
        raise ValueError("Fourier trace branches must contain three valid cases")
    if prefix_prior_bounds is not None:
        prefix_prior_bounds = np.asarray(prefix_prior_bounds, dtype=float)
        if prefix_prior_bounds.shape != (4, 2):
            raise ValueError("prefix prior bounds must have shape (4,2)")
        if np.any(prefix_prior_bounds[:, 0] < 0.0) or np.any(
            prefix_prior_bounds[:, 1] > 1.0
        ):
            raise ValueError("prefix prior bounds must lie in [0,1]")
        if np.any(prefix_prior_bounds[:, 0] > prefix_prior_bounds[:, 1]):
            raise ValueError("prefix prior lower bounds exceed upper bounds")
    traces = np.trace(effects, axis1=1, axis2=2).real
    directions = np.zeros((4, 3), dtype=float)
    for s in OUTCOMES:
        if traces[s] > 1e-12:
            directions[s] = bloch(effects[s] / traces[s])[1:]

    model = Model("joint-effect-helstrom")
    variables: dict[str, object] = {}

    state_scalar = []
    state_vector = []
    for z in OUTCOMES:
        lower, upper = (
            (0.0, 1.0) if prefix_prior_bounds is None else tuple(prefix_prior_bounds[z])
        )
        scalar = model.addVar(lb=float(lower), ub=float(upper), name=f"a_{z}")
        vector = tuple(
            model.addVar(lb=-1.0, ub=1.0, name=f"r_{z}_{axis}") for axis in range(3)
        )
        add_lorentz(model, scalar, vector)
        state_scalar.append(scalar)
        state_vector.append(vector)
        variables[f"a_{z}"] = scalar
        for axis, item in enumerate(vector):
            variables[f"r_{z}_{axis}"] = item
    model.addCons(quicksum(state_scalar) == 1.0)
    if fix_rotation_gauge:
        # A common SO(3) rotation of all input Bloch vectors and pulled
        # effects is physically irrelevant.  Rotate r_0 to the positive
        # x-axis and r_1 into the upper xy half-plane.  The convention also
        # covers zero/collinear vectors and removes three continuous gauge
        # directions without excluding a physical orbit.
        model.addCons(state_vector[0][1] == 0.0)
        model.addCons(state_vector[0][2] == 0.0)
        model.addCons(state_vector[0][0] >= 0.0)
        model.addCons(state_vector[1][2] == 0.0)
        model.addCons(state_vector[1][1] >= 0.0)

    effect_scalar: dict[tuple[int, int], object] = {}
    effect_vector: dict[tuple[int, int], tuple[object, object, object]] = {}
    for y, s in PATHS:
        upper = 2.0 * float(traces[s])
        scalar = model.addVar(lb=0.0, ub=upper, name=f"g_{y}_{s}_0")
        vector = tuple(
            model.addVar(lb=-upper, ub=upper, name=f"g_{y}_{s}_{axis + 1}")
            for axis in range(3)
        )
        add_lorentz(model, scalar, vector)
        effect_scalar[y, s] = scalar
        effect_vector[y, s] = vector
        variables[f"g_{y}_{s}_0"] = scalar
        for axis, item in enumerate(vector, start=1):
            variables[f"g_{y}_{s}_{axis}"] = item

    model.addCons(quicksum(effect_scalar.values()) == 2.0)
    for axis in range(3):
        model.addCons(quicksum(effect_vector[y, s][axis] for y, s in PATHS) == 0.0)

    coarse_scalar = {
        y: quicksum(effect_scalar[y, s] for s in OUTCOMES) for y in OUTCOMES
    }
    coarse_vector = {
        y: tuple(
            quicksum(effect_vector[y, s][axis] for s in OUTCOMES) for axis in range(3)
        )
        for y in OUTCOMES
    }
    for y, s in PATHS:
        residual_scalar = traces[s] * coarse_scalar[y] - effect_scalar[y, s]
        residual_vector = tuple(
            traces[s] * coarse_vector[y][axis] - effect_vector[y, s][axis]
            for axis in range(3)
        )
        add_lorentz(model, residual_scalar, residual_vector)

    cp_factors: dict[int, dict[tuple[int, int, str], object]] = {}
    cp_missing_pullback: dict[tuple[int, int], object] = {}
    planar_pullback: dict[tuple[int, int, int], object] = {}
    if require_cp_completion or ando_direction_count > 0:
        active_for_completion = tuple(int(s) for s in OUTCOMES if traces[s] > 1e-9)
        if len(active_for_completion) != 3:
            raise ValueError("CP completion requires three active planar effects")
        reconstruction = np.asarray(
            [
                [
                    0.5 * traces[s],
                    0.5 * traces[s] * directions[s, 0],
                    0.5 * traces[s] * directions[s, 1],
                ]
                for s in active_for_completion
            ]
        )
        if abs(float(np.linalg.det(reconstruction))) < 1e-10:
            raise ValueError("active terminal effects do not span the plane")
        inverse = np.linalg.inv(reconstruction)

        def pulled_component(outcome: int, terminal: int, mu: int) -> object:
            if mu == 0:
                return effect_scalar[outcome, terminal]
            return effect_vector[outcome, terminal][mu - 1]

        def coarse_component(outcome: int, mu: int) -> object:
            if mu == 0:
                return coarse_scalar[outcome]
            return coarse_vector[outcome][mu - 1]

        for y in OUTCOMES:
            for mu in OUTCOMES:
                reconstructed = tuple(
                    quicksum(
                        float(inverse[row, column])
                        * pulled_component(y, active_for_completion[column], mu)
                        for column in range(3)
                    )
                    for row in range(3)
                )
                # The zeroth reconstructed operator is Phi_y^*(I).  Keeping
                # this redundant identity materially strengthens spatial
                # relaxations and audits the floating-point basis inverse.
                model.addCons(reconstructed[0] == coarse_component(y, mu))
                planar_pullback[y, 0, mu] = coarse_component(y, mu)
                planar_pullback[y, 1, mu] = reconstructed[1]
                planar_pullback[y, 2, mu] = reconstructed[2]

        for y in OUTCOMES:
            for direction in fibonacci_sphere(ando_direction_count):
                expectations = [
                    0.5
                    * (
                        planar_pullback[y, domain, 0]
                        + quicksum(
                            float(direction[axis])
                            * planar_pullback[y, domain, axis + 1]
                            for axis in range(3)
                        )
                    )
                    for domain in range(3)
                ]
                add_lorentz(model, expectations[0], expectations[1:])

    if require_cp_completion:
        for y in OUTCOMES:
            coefficients: dict[tuple[int, int], object] = {}
            for mu in OUTCOMES:
                coefficients[0, mu] = planar_pullback[y, 0, mu]
                coefficients[1, mu] = planar_pullback[y, 1, mu]
                # Choi expansion transposes the domain Pauli Y.
                coefficients[2, mu] = -planar_pullback[y, 2, mu]
                missing = model.addVar(
                    lb=-2.0,
                    ub=2.0,
                    name=f"cp_missing_{y}_{mu}",
                )
                cp_missing_pullback[y, mu] = missing
                coefficients[3, mu] = missing
            cp_factors[y] = add_complex_cholesky(model, coefficients, f"cp_choi_{y}")

    active = tuple(int(s) for s in OUTCOMES if traces[s] > 1e-9)
    path_reconstruction = np.asarray(
        [
            [
                0.5 * traces[s],
                0.5 * traces[s] * directions[s, 0],
                0.5 * traces[s] * directions[s, 1],
            ]
            for s in active
        ]
    )
    path_reconstruction_inverse = np.linalg.inv(path_reconstruction)
    linked = (
        {f"b_{3 * y + t}" for y in OUTCOMES for t in active}
        if linked_columns is None
        else set(linked_columns)
    )
    statistics: dict[tuple[int, int, int], object] = {}
    probability: dict[tuple[int, int], object] = {}
    correct: dict[tuple[int, int], object] = {}
    path_output_planar: dict[tuple[int, int, int], object] = {}
    path_output_full: dict[tuple[int, int, int], object] = {}
    for z, y in PATHS:
        p = model.addVar(lb=0.0, ub=1.0, name=f"p_{z}_{y}")
        model.addCons(
            2.0 * p
            == state_scalar[z] * coarse_scalar[y]
            + quicksum(
                state_vector[z][axis] * coarse_vector[y][axis] for axis in range(3)
            )
        )
        for terminal_label in active:
            item = model.addVar(
                lb=0.0,
                ub=float(traces[terminal_label]),
                name=f"q_{z}_{y}_{terminal_label}",
            )
            if f"b_{3 * y + terminal_label}" in linked:
                model.addCons(
                    2.0 * item
                    == state_scalar[z] * effect_scalar[y, terminal_label]
                    + quicksum(
                        state_vector[z][axis] * effect_vector[y, terminal_label][axis]
                        for axis in range(3)
                    )
                )
            model.addCons(item <= traces[terminal_label] * p)
            statistics[z, y, terminal_label] = item
            variables[f"q_{z}_{y}_{terminal_label}"] = item
        model.addCons(p == quicksum(statistics[z, y, t] for t in active))
        reconstructed_output = tuple(
            quicksum(
                float(path_reconstruction_inverse[row, column])
                * statistics[z, y, active[column]]
                for column in range(3)
            )
            for row in range(3)
        )
        model.addCons(reconstructed_output[0] == p)
        if require_cp_completion:
            missing_output = 0.5 * (
                state_scalar[z] * cp_missing_pullback[y, 0]
                + quicksum(
                    state_vector[z][axis] * cp_missing_pullback[y, axis + 1]
                    for axis in range(3)
                )
            )
            add_lorentz(
                model,
                reconstructed_output[0],
                [reconstructed_output[1], reconstructed_output[2], missing_output],
            )
            full_output = (*reconstructed_output, missing_output)
        else:
            add_lorentz(model, reconstructed_output[0], list(reconstructed_output[1:]))
            full_output = (*reconstructed_output, 0.0)
        for domain, item in enumerate(reconstructed_output):
            path_output_planar[z, y, domain] = item
        for domain, item in enumerate(full_output):
            path_output_full[z, y, domain] = item
        s = z ^ y
        d = statistics[z, y, s] if s in active else 0.0
        probability[z, y] = p
        correct[z, y] = d
        variables[f"p_{z}_{y}"] = p

    if basis_determinant_sign:
        input_basis = [
            [state_scalar[z], *state_vector[z]] for z in OUTCOMES
        ]
        determinant = polynomial_determinant(input_basis)
        signed_determinant = float(basis_determinant_sign) * determinant
        model.addCons(signed_determinant >= basis_determinant_floor)
        # L_y = Q_y^T adj(R^T) / det(R).  Twice its numerator gives the
        # Pauli coefficients expected by add_complex_cholesky, whose /4 Choi
        # convention then constructs sign(det(R))*det(R)*J_y = |det(R)|J_y.
        input_transpose = [
            [input_basis[column][row] for column in OUTCOMES]
            for row in OUTCOMES
        ]
        adjugate_transpose = polynomial_adjugate(input_transpose)
        basis_factors = {}
        for y in OUTCOMES:
            coefficients: dict[tuple[int, int], object] = {}
            for mu in OUTCOMES:
                for nu in OUTCOMES:
                    transfer_numerator = quicksum(
                        path_output_full[z, y, nu]
                        * adjugate_transpose[z][mu]
                        for z in OUTCOMES
                    )
                    coefficients[mu, nu] = (
                        2.0 * basis_determinant_sign * transfer_numerator
                    )
            # The helper transposes the domain Y internally through the sign
            # convention already used for adjoint-map CP completion.
            coefficients[2, 0] = -coefficients[2, 0]
            coefficients[2, 1] = -coefficients[2, 1]
            coefficients[2, 2] = -coefficients[2, 2]
            coefficients[2, 3] = -coefficients[2, 3]
            if basis_choi_witnesses is not None:
                for witness_index, witness in enumerate(basis_choi_witnesses[y]):
                    norm = float(np.linalg.norm(witness))
                    if norm <= 1e-14:
                        raise ValueError("basis Choi witnesses must be nonzero")
                    direction = witness / norm
                    expectation = quicksum(
                        float(
                            np.vdot(
                                direction,
                                CHOI_BASIS[mu, nu] @ direction,
                            ).real
                        )
                        * coefficients[mu, nu]
                        / 4.0
                        for mu in OUTCOMES
                        for nu in OUTCOMES
                    )
                    # The determinant numerator is naturally small on the
                    # target prior box. Scaling improves feasibility checks
                    # without changing the mathematical half-space.
                    model.addCons(
                        1000.0 * expectation >= 0.0,
                        name=f"basis_witness_{y}_{witness_index}",
                    )
            basis_factors[y] = add_complex_cholesky(
                model, coefficients, f"basis_choi_{y}"
            )
        variables["basis_determinant"] = determinant
        variables["basis_determinant_sign"] = int(basis_determinant_sign)
        variables["basis_cholesky_factors"] = basis_factors

    for name in sorted(linked):
        parts = name.split("_")
        if parts[0] != "d":
            if parts[0] != "b" or len(parts) != 2:
                raise ValueError(f"invalid linked column {name!r}")
            continue
        if len(parts) != 3:
            raise ValueError(f"invalid linked residual {name!r}")
        y, terminal_label = int(parts[1]), int(parts[2])
        if y not in OUTCOMES or terminal_label not in active:
            raise ValueError(f"invalid linked residual {name!r}")
        residual_scalar = (
            traces[terminal_label] * coarse_scalar[y] - effect_scalar[y, terminal_label]
        )
        residual_vector = tuple(
            traces[terminal_label] * coarse_vector[y][axis]
            - effect_vector[y, terminal_label][axis]
            for axis in range(3)
        )
        for z in OUTCOMES:
            model.addCons(
                2.0
                * (
                    traces[terminal_label] * probability[z, y]
                    - statistics[z, y, terminal_label]
                )
                == state_scalar[z] * residual_scalar
                + quicksum(
                    state_vector[z][axis] * residual_vector[axis] for axis in range(3)
                )
            )

    model.addCons(quicksum(probability.values()) == 1.0)
    for z in OUTCOMES:
        # This follows algebraically from POVM completeness.  Stating it
        # explicitly tightens node relaxations before all effect products
        # have been resolved.
        model.addCons(quicksum(probability[z, y] for y in OUTCOMES) == state_scalar[z])
    audit = model.addVar(lb=0.0, ub=1.0, name="audit")
    model.addCons(audit == quicksum(correct.values()))
    variables["audit"] = audit

    prefix = {z: quicksum(probability[z, y] for y in OUTCOMES) for z in OUTCOMES}
    if prefix_order is not None:
        for index in range(3):
            model.addCons(
                prefix[prefix_order[index]] >= prefix[prefix_order[index + 1]]
            )
        maximum = float(traces.max())
        cap = [maximum, maximum, 2.0 - 2.0 * maximum, 0.0]
        model.addCons(
            audit
            <= quicksum(cap[index] * prefix[prefix_order[index]] for index in OUTCOMES)
        )

    syndrome_from_paths = {
        s: quicksum(probability[z, z ^ s] for z in OUTCOMES) for s in OUTCOMES
    }
    terminal_statistics = {
        (s, t): quicksum(statistics[z, y, t] for z, y in PATHS if (z ^ y) == s)
        for s in OUTCOMES
        for t in active
    }
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
    syndrome = {}
    terminal_vector = []
    for s in OUTCOMES:
        reconstructed = tuple(
            quicksum(
                float(reconstruction_inverse[row, column])
                * terminal_statistics[s, active[column]]
                for column in range(3)
            )
            for row in range(3)
        )
        syndrome[s] = reconstructed[0]
        model.addCons(syndrome[s] == syndrome_from_paths[s])
        normal = model.addVar(lb=-1.0, ub=1.0, name=f"tau_{s}_2")
        vector = (reconstructed[1], reconstructed[2], normal)
        add_lorentz(model, syndrome[s], vector)
        if require_cp_completion:
            # The planar POVM statistics determine only I, X, and Y.  The
            # actual Z component must use the very same missing pullback that
            # completes each positive Choi matrix; otherwise the Helstrom
            # certificate could choose an unphysical terminal state.
            model.addCons(
                2.0 * normal
                == quicksum(
                    state_scalar[z] * cp_missing_pullback[z ^ s, 0]
                    + quicksum(
                        state_vector[z][axis] * cp_missing_pullback[z ^ s, axis + 1]
                        for axis in range(3)
                    )
                    for z in OUTCOMES
                )
            )
        terminal_vector.append(vector)
        variables[f"tau_{s}_2"] = normal

    dual_vector = tuple(
        model.addVar(lb=-1.0, ub=1.0, name=f"dual_{axis}") for axis in range(3)
    )
    add_lorentz(model, audit, dual_vector)
    for axis, item in enumerate(dual_vector):
        variables[f"dual_{axis}"] = item
    for s in OUTCOMES:
        add_lorentz(
            model,
            audit - syndrome[s],
            tuple(dual_vector[axis] - terminal_vector[s][axis] for axis in range(3)),
        )
        if s in active:
            # Strong duality plus optimality of the fixed rank-one effect
            # P_s = w_s Pi_s exposes a one-dimensional Helstrom face:
            # Y - tau_s = (audit - Tr(tau_s)) (I - Pi_s). These linear
            # identities follow from the existing primal/dual constraints,
            # but stating them explicitly removes a large artificial root
            # relaxation before the bilinear terminal coordinates are fixed.
            for axis in range(3):
                model.addCons(
                    terminal_vector[s][axis]
                    == dual_vector[axis] + (audit - syndrome[s]) * directions[s, axis]
                )
    model.addCons(audit <= quicksum(traces[s] * syndrome[s] for s in OUTCOMES))

    # Fourier-energy common-channel cuts. For every nontrivial character of
    # Z_2^2, convolution gives tau_hat(k) = Phi_hat(k)(rho_hat(k)). The signed
    # instrument map is trace-norm contractive on Hermitian inputs, hence
    # ||tau_hat(k)||_1 <= ||rho_hat(k)||_1. For a qubit Bloch matrix
    # X=(x0 I+x.sigma)/2, ||X||_1^2=max(x0^2,||x||^2) <= x0^2+||x||^2.
    # The variables below retain the exact terminal max and only relax the
    # input side by this last, explicit factor-of-at-most-two inequality.
    fourier_norm_square = []
    fourier_input_energy = []
    fourier_flagged_sum = []
    fourier_flagged_norm_variables: dict[tuple[int, int], object] = {}
    for character in range(1, 4):
        signs = tuple(
            -1.0 if (character & label).bit_count() % 2 else 1.0 for label in OUTCOMES
        )
        terminal_scalar = quicksum(signs[s] * syndrome[s] for s in OUTCOMES)
        terminal_bloch = tuple(
            quicksum(signs[s] * terminal_vector[s][axis] for s in OUTCOMES)
            for axis in range(3)
        )
        prefix_scalar = quicksum(signs[z] * state_scalar[z] for z in OUTCOMES)
        prefix_bloch = tuple(
            quicksum(signs[z] * state_vector[z][axis] for z in OUTCOMES)
            for axis in range(3)
        )
        norm_square = model.addVar(
            lb=0.0, ub=1.0, name=f"fourier_terminal_norm2_{character}"
        )
        input_energy = prefix_scalar * prefix_scalar + quicksum(
            item * item for item in prefix_bloch
        )
        model.addCons(norm_square >= terminal_scalar * terminal_scalar)
        model.addCons(norm_square >= quicksum(item * item for item in terminal_bloch))
        model.addCons(norm_square <= input_energy)

        # Do not allow cancellations between classical outcomes. The flagged
        # CPTP channel obeys
        #   sum_y ||Phi_y(rho_hat(k))||_1 <= ||rho_hat(k)||_1.
        # Only planar output coordinates are observed, so their trace norms
        # are certified lower bounds on the full output norms. The input
        # trace norm is again relaxed by its Bloch Euclidean energy.
        flagged_norms = []
        for y in OUTCOMES:
            output_scalar = quicksum(
                signs[z] * path_output_planar[z, y, 0] for z in OUTCOMES
            )
            output_planar = tuple(
                quicksum(signs[z] * path_output_planar[z, y, domain] for z in OUTCOMES)
                for domain in (1, 2)
            )
            flagged_norm = model.addVar(
                lb=0.0,
                ub=1.0,
                name=f"fourier_flagged_norm_{character}_{y}",
            )
            model.addCons(flagged_norm >= output_scalar)
            model.addCons(flagged_norm >= -output_scalar)
            add_lorentz(model, flagged_norm, [*output_planar, 0.0])
            flagged_norms.append(flagged_norm)
            fourier_flagged_norm_variables[character, y] = flagged_norm
        flagged_sum = quicksum(flagged_norms)
        fourier_flagged_sum.append(flagged_sum)
        model.addCons(flagged_sum * flagged_sum <= input_energy)
        if fourier_trace_branches is not None:
            branch = fourier_trace_branches[character - 1]
            prefix_bloch_square = quicksum(item * item for item in prefix_bloch)
            if branch == "scalar-positive":
                model.addCons(prefix_scalar >= 0.0)
                model.addCons(prefix_bloch_square <= prefix_scalar * prefix_scalar)
                model.addCons(flagged_sum <= prefix_scalar)
            elif branch == "scalar-negative":
                model.addCons(prefix_scalar <= 0.0)
                model.addCons(prefix_bloch_square <= prefix_scalar * prefix_scalar)
                model.addCons(flagged_sum <= -prefix_scalar)
            else:
                model.addCons(prefix_bloch_square >= prefix_scalar * prefix_scalar)
                model.addCons(flagged_sum * flagged_sum <= prefix_bloch_square)
        fourier_norm_square.append(norm_square)
        fourier_input_energy.append(input_energy)

    # Parseval states the total input energy in the original state
    # coordinates. Keeping both forms helps propagation in spatial nodes.
    total_input_bloch = tuple(
        quicksum(state_vector[z][axis] for z in OUTCOMES) for axis in range(3)
    )
    parseval_energy = (
        4.0
        * quicksum(
            state_scalar[z] * state_scalar[z]
            + quicksum(item * item for item in state_vector[z])
            for z in OUTCOMES
        )
        - 1.0
        - quicksum(item * item for item in total_input_bloch)
    )
    model.addCons(quicksum(fourier_norm_square) <= parseval_energy)
    if prefix_prior_bounds is not None and fourier_trace_branches == (
        "bloch",
        "bloch",
        "bloch",
    ):
        prior_square_secant = quicksum(
            float(prefix_prior_bounds[z, 0] + prefix_prior_bounds[z, 1])
            * state_scalar[z]
            - float(prefix_prior_bounds[z, 0] * prefix_prior_bounds[z, 1])
            for z in OUTCOMES
        )
        model.addCons(
            quicksum(item * item for item in fourier_flagged_sum)
            <= 4.0 * prior_square_secant
            - quicksum(item * item for item in total_input_bloch)
        )
    variables["fourier_terminal_norm_square"] = tuple(fourier_norm_square)
    variables["fourier_input_energy"] = tuple(fourier_input_energy)
    variables["fourier_flagged_norms"] = fourier_flagged_norm_variables

    flat = [probability[z, y] for z, y in PATHS]
    hellinger = []
    for first in range(16):
        for second in range(first + 1, 16):
            item = model.addVar(lb=0.0, ub=0.5, name=f"h_{first}_{second}")
            model.addCons(item * item <= flat[first] * flat[second])
            hellinger.append(item)
            variables[f"h_{first}_{second}"] = item
    returned = (1.0 + 2.0 * quicksum(hellinger)) / 16.0
    score = model.addVar(lb=0.0, ub=1.0, name="score")
    model.addCons(score <= weight * audit + (1.0 - weight) * returned)
    if target is not None:
        # A proof that this constrained model is infeasible certifies the
        # desired strict upper bound and is usually easier than closing the
        # full optimisation gap.
        model.addCons(score >= target)
    model.setObjective(score, "maximize")
    variables.update(
        {
            "score": score,
            "return": returned,
            "terminal_effect_weights": traces,
            "terminal_directions": directions,
            "cp_factors": cp_factors,
            "cp_missing_pullback": cp_missing_pullback,
            "planar_pullback": planar_pullback,
            "ando_direction_count": ando_direction_count,
            "rotation_gauge_fixed": bool(fix_rotation_gauge),
        }
    )
    return model, variables


def seed_solution(
    model: Model,
    variables: dict[str, object],
    checkpoint: Path,
    effects: np.ndarray,
    weight: float,
) -> None:
    arrays = np.load(checkpoint)
    states = np.asarray(arrays["states"])
    joint = np.asarray(arrays["joint"])
    coarse = joint.sum(axis=1)
    probabilities = np.einsum("yij,zji->zy", coarse, states).real
    correct = np.asarray(
        [
            [np.trace(joint[y, z ^ y] @ states[z]).real for y in OUTCOMES]
            for z in OUTCOMES
        ]
    )
    audit = float(correct.sum())
    if "terminal_vectors" not in arrays or "dual_vector" not in arrays:
        return
    terminal_vectors = np.asarray(arrays["terminal_vectors"])
    dual_vector = np.asarray(arrays["dual_vector"])

    values: dict[str, float] = {}
    for z in OUTCOMES:
        vector = bloch(states[z])
        values[f"a_{z}"] = vector[0]
        for axis in range(3):
            values[f"r_{z}_{axis}"] = vector[axis + 1]
    for y, s in PATHS:
        vector = bloch(joint[y, s])
        for axis in range(4):
            values[f"g_{y}_{s}_{axis}"] = vector[axis]
    for z, y in PATHS:
        values[f"p_{z}_{y}"] = probabilities[z, y]
        values[f"d_{z}_{y}"] = correct[z, y]
    values["audit"] = audit
    for s in OUTCOMES:
        for axis in range(3):
            values[f"tau_{s}_{axis}"] = terminal_vectors[s, axis]
    for axis in range(3):
        values[f"dual_{axis}"] = dual_vector[axis]
    flat = probabilities.reshape(16)
    for first in range(16):
        for second in range(first + 1, 16):
            values[f"h_{first}_{second}"] = math.sqrt(
                max(0.0, flat[first] * flat[second])
            )
    returned = float(np.sqrt(np.maximum(probabilities, 0.0)).sum() ** 2 / 16.0)
    values["score"] = weight * audit + (1.0 - weight) * returned

    solution = model.createSol()
    for name, value in values.items():
        if name in variables:
            model.setSolVal(solution, variables[name], float(value))
    model.addSol(solution)


def apply_choi(choi: np.ndarray, state: np.ndarray) -> np.ndarray:
    """Apply an input-major Choi matrix to one input operator."""

    blocks = np.asarray(choi).reshape(2, 2, 2, 2)
    output = np.einsum("ij,iajb->ab", state, blocks)
    return 0.5 * (output + output.conj().T)


def apply_choi_adjoint(choi: np.ndarray, effect: np.ndarray) -> np.ndarray:
    """Apply the Hilbert--Schmidt adjoint of an input-major Choi map."""

    blocks = np.asarray(choi).reshape(2, 2, 2, 2)
    pulled = np.einsum("ba,iajb->ji", effect, blocks)
    return 0.5 * (pulled + pulled.conj().T)


def seed_from_common_instrument(
    model: Model,
    variables: dict[str, object],
    checkpoint: Path,
    effects: np.ndarray,
    weight: float,
) -> bool:
    """Seed the CP-complete model from literal prefix states and Choi maps."""

    arrays = np.load(checkpoint)
    states = np.asarray(arrays["states"], dtype=complex)
    choi = np.asarray(arrays["choi"], dtype=complex)
    if states.shape != (4, 2, 2) or choi.shape != (4, 4, 4):
        raise ValueError("common-instrument seed must contain states and choi")
    if variables.get("rotation_gauge_fixed", False):
        states, choi = rotate_common_instrument_input_gauge(states, choi)

    input_basis = np.asarray([bloch(state) for state in states])
    input_determinant = float(np.linalg.det(input_basis))
    if "basis_determinant_sign" in variables:
        required_sign = int(variables["basis_determinant_sign"])
        if required_sign * input_determinant <= 0.0:
            return False

    joint = np.asarray(
        [[apply_choi_adjoint(choi[y], effects[s]) for s in OUTCOMES] for y in OUTCOMES]
    )
    coarse = joint.sum(axis=1)
    probabilities = np.einsum("yij,zji->zy", coarse, states).real
    statistics = np.asarray(
        [
            [
                [float(np.trace(joint[y, s] @ states[z]).real) for s in OUTCOMES]
                for y in OUTCOMES
            ]
            for z in OUTCOMES
        ]
    )
    audit = float(sum(statistics[z, y, z ^ y] for z, y in PATHS))
    terminal_states = np.asarray(
        [
            sum(apply_choi(choi[z ^ syndrome], states[z]) for z in OUTCOMES)
            for syndrome in OUTCOMES
        ]
    )
    dual_matrix = sum(effects[s] @ terminal_states[s] for s in OUTCOMES)
    dual_matrix = 0.5 * (dual_matrix + dual_matrix.conj().T)
    dual = bloch(dual_matrix)
    flat = probabilities.reshape(16)
    returned = float(np.sqrt(np.maximum(flat, 0.0)).sum() ** 2 / 16.0)
    conditioned = np.asarray(
        [[apply_choi(choi[y], states[z]) for y in OUTCOMES] for z in OUTCOMES]
    )

    solution = model.createSol()
    for z in OUTCOMES:
        state_vector = bloch(states[z])
        model.setSolVal(solution, variables[f"a_{z}"], float(state_vector[0]))
        for axis in range(3):
            model.setSolVal(
                solution,
                variables[f"r_{z}_{axis}"],
                float(state_vector[axis + 1]),
            )
    for y, s in PATHS:
        vector = bloch(joint[y, s])
        for mu in OUTCOMES:
            model.setSolVal(solution, variables[f"g_{y}_{s}_{mu}"], float(vector[mu]))
    for z, y in PATHS:
        model.setSolVal(solution, variables[f"p_{z}_{y}"], probabilities[z, y])
        for s in OUTCOMES:
            name = f"q_{z}_{y}_{s}"
            if name in variables:
                model.setSolVal(solution, variables[name], statistics[z, y, s])
    model.setSolVal(solution, variables["audit"], audit)
    for syndrome in OUTCOMES:
        terminal = bloch(terminal_states[syndrome])
        model.setSolVal(solution, variables[f"tau_{syndrome}_2"], float(terminal[3]))
    for axis in range(3):
        model.setSolVal(solution, variables[f"dual_{axis}"], float(dual[axis + 1]))
    for first in range(16):
        for second in range(first + 1, 16):
            model.setSolVal(
                solution,
                variables[f"h_{first}_{second}"],
                float(math.sqrt(max(0.0, flat[first] * flat[second]))),
            )

    terminal_norms = variables["fourier_terminal_norm_square"]
    flagged_norms = variables["fourier_flagged_norms"]
    for character in range(1, 4):
        signs = np.asarray(
            [
                -1.0 if (character & label).bit_count() % 2 else 1.0
                for label in OUTCOMES
            ]
        )
        terminal_fourier = bloch(np.einsum("s,sij->ij", signs, terminal_states))
        terminal_norm_square = max(
            terminal_fourier[0] ** 2,
            float(np.dot(terminal_fourier[1:], terminal_fourier[1:])),
        )
        model.setSolVal(
            solution,
            terminal_norms[character - 1],
            float(terminal_norm_square),
        )
        for y in OUTCOMES:
            output_fourier = bloch(
                np.einsum("z,zij->ij", signs, conditioned[:, y])
            )
            planar_norm = math.sqrt(
                float(output_fourier[1] ** 2 + output_fourier[2] ** 2)
            )
            model.setSolVal(
                solution,
                flagged_norms[character, y],
                max(abs(float(output_fourier[0])), planar_norm),
            )
    model.setSolVal(
        solution,
        variables["score"],
        weight * audit + (1.0 - weight) * returned,
    )

    factors = variables["cp_factors"]
    missing = variables["cp_missing_pullback"]
    for y in OUTCOMES:
        pullbacks = [apply_choi_adjoint(choi[y], FULL_PAULIS[mu]) for mu in OUTCOMES]
        missing_vector = bloch(pullbacks[3])
        for mu in OUTCOMES:
            model.setSolVal(solution, missing[y, mu], float(missing_vector[mu]))
        adjoint_choi = 0.5 * sum(
            (1.0 if domain != 2 else -1.0)
            * np.kron(FULL_PAULIS[domain], pullbacks[domain])
            for domain in OUTCOMES
        )
        adjoint_choi = 0.5 * (adjoint_choi + adjoint_choi.conj().T)
        factor = np.linalg.cholesky(adjoint_choi)
        for row in OUTCOMES:
            for column in range(row + 1):
                model.setSolVal(
                    solution,
                    factors[y][row, column, "real"],
                    float(factor[row, column].real),
                )
                if row != column:
                    model.setSolVal(
                        solution,
                        factors[y][row, column, "imag"],
                        float(factor[row, column].imag),
                    )

    if "basis_cholesky_factors" in variables:
        basis_factors = variables["basis_cholesky_factors"]
        for y in OUTCOMES:
            # The operator-basis constraint factors |det(R)| J_y.  Since this
            # seed already supplies literal Choi matrices, it can populate the
            # redundant determinant-numerator factor without an inverse.
            signed_numerator = abs(input_determinant) * choi[y]
            factor = lower_psd_factor(signed_numerator)
            for row in OUTCOMES:
                for column in range(row + 1):
                    model.setSolVal(
                        solution,
                        basis_factors[y][row, column, "real"],
                        float(factor[row, column].real),
                    )
                    if row != column:
                        model.setSolVal(
                            solution,
                            basis_factors[y][row, column, "imag"],
                            float(factor[row, column].imag),
                        )
    seed_feasible = bool(
        model.checkSol(solution, printreason=False, completely=True)
    )
    if not seed_feasible:
        return False
    return bool(model.addSol(solution, free=True))


def extract_common_instrument(
    model: Model,
    variables: dict[str, object],
    effects: np.ndarray,
    weight: float,
) -> tuple[dict[str, object], dict[str, np.ndarray]] | None:
    """Reconstruct one Schrödinger-picture instrument from the best solution."""

    solution = model.getBestSol()
    if solution is None:
        return None

    def value(variable: object) -> float:
        return float(model.getSolVal(solution, variable))

    states = np.asarray(
        [
            0.5
            * sum(
                value(variables[f"a_{z}"]) * FULL_PAULIS[0]
                if mu == 0
                else value(variables[f"r_{z}_{mu - 1}"]) * FULL_PAULIS[mu]
                for mu in OUTCOMES
            )
            for z in OUTCOMES
        ]
    )
    traces = np.trace(effects, axis1=1, axis2=2).real
    directions = np.asarray(variables["terminal_directions"], dtype=float)
    active = tuple(int(s) for s in OUTCOMES if traces[s] > 1e-9)
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
    inverse = np.linalg.inv(reconstruction)
    missing = variables["cp_missing_pullback"]
    pullback_coefficients = np.zeros((4, 4, 4), dtype=float)
    for y in OUTCOMES:
        pulled = np.asarray(
            [[value(variables[f"g_{y}_{s}_{mu}"]) for mu in OUTCOMES] for s in active]
        )
        pullback_coefficients[y, :3, :] = inverse @ pulled
        pullback_coefficients[y, 3, :] = [value(missing[y, mu]) for mu in OUTCOMES]

    transpose_sign = np.asarray([1.0, 1.0, -1.0, 1.0])
    choi_coefficients = np.zeros((4, 4, 4), dtype=float)
    for y in OUTCOMES:
        choi_coefficients[y] = transpose_sign[:, None] * pullback_coefficients[y].T
    choi = np.asarray(
        [
            sum(
                choi_coefficients[y, mu, nu] * CHOI_BASIS[mu, nu]
                for mu in OUTCOMES
                for nu in OUTCOMES
            )
            / 4.0
            for y in OUTCOMES
        ]
    )

    outputs = np.asarray(
        [[apply_choi(choi[y], states[z]) for y in OUTCOMES] for z in OUTCOMES]
    )
    probabilities = np.trace(outputs, axis1=2, axis2=3).real
    statistics = np.einsum("tij,zyji->zyt", effects, outputs).real
    audit = float(sum(statistics[z, y, z ^ y] for z, y in PATHS))
    returned = float(np.sqrt(np.maximum(probabilities, 0.0)).sum() ** 2 / 16.0)
    score = weight * audit + (1.0 - weight) * returned
    trace_preservation = choi.sum(axis=0).reshape(2, 2, 2, 2).trace(axis1=1, axis2=3)
    report: dict[str, object] = {
        "score_from_reconstructed_instrument": score,
        "audit_from_reconstructed_instrument": audit,
        "return_from_reconstructed_instrument": returned,
        "normalisation": float(probabilities.sum()),
        "minimum_state_eigenvalue": float(
            min(np.linalg.eigvalsh(item).min() for item in states)
        ),
        "minimum_choi_eigenvalue": float(
            min(np.linalg.eigvalsh(item).min() for item in choi)
        ),
        "minimum_output_eigenvalue": float(
            min(np.linalg.eigvalsh(item).min() for row in outputs for item in row)
        ),
        "trace_preservation_residual": float(
            np.linalg.norm(trace_preservation - np.eye(2))
        ),
        "reported_score_residual": abs(score - value(variables["score"])),
    }
    return report, {"states": states, "choi": choi, "effects": effects}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lambda", dest="weight", type=float, default=0.6)
    parser.add_argument(
        "--fixed-three-povm-weights", type=float, nargs=3, required=True
    )
    parser.add_argument("--prefix-order", type=int, nargs=4)
    parser.add_argument("--seed-npz", type=Path)
    parser.add_argument("--common-instrument-seed-npz", type=Path)
    parser.add_argument("--solution-npz", type=Path)
    parser.add_argument("--seconds", type=float, default=300.0)
    parser.add_argument("--gap", type=float, default=1e-5)
    parser.add_argument(
        "--feasibility-tolerance",
        type=float,
        default=1e-9,
        help="SCIP primal feasibility tolerance, set before checking warm starts",
    )
    parser.add_argument("--target", type=float)
    parser.add_argument("--no-rotation-gauge", action="store_true")
    parser.add_argument(
        "--require-cp-completion",
        action="store_true",
        help="require every pulled-effect triple to have one positive Choi completion",
    )
    parser.add_argument(
        "--ando-directions",
        type=int,
        default=0,
        help="finite pure-state cover of the exact numerical-radius inequality",
    )
    parser.add_argument(
        "--fourier-trace-branches",
        nargs=3,
        choices=("scalar-positive", "scalar-negative", "bloch"),
        help=("exact trace-norm cases for the three nontrivial Z2^2 Fourier modes"),
    )
    parser.add_argument(
        "--linked-column",
        action="append",
        default=[],
        help=(
            "link one b_J or d_y_t statistic to its Born product; "
            "when omitted, all twelve b_J columns are linked"
        ),
    )
    parser.add_argument(
        "--prefix-prior-bounds",
        type=float,
        nargs=8,
        metavar=("L0", "U0", "L1", "U1", "L2", "U2", "L3", "U3"),
        help="optional four prior intervals used for spatial secant cuts",
    )
    parser.add_argument(
        "--basis-determinant-sign",
        type=int,
        choices=(-1, 0, 1),
        default=0,
        help="add redundant operator-basis Choi positivity on one determinant-sign branch",
    )
    parser.add_argument(
        "--basis-determinant-floor",
        type=float,
        default=0.0,
        help="minimum absolute input-basis determinant on the selected sign branch",
    )
    parser.add_argument(
        "--basis-choi-witnesses-npz",
        type=Path,
        help="NPZ containing witnesses with shape (4,cuts,4)",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    weights = np.asarray(args.fixed_three_povm_weights, dtype=float)
    effects = canonical_three_effect_povm(weights)
    prefix_order = None if args.prefix_order is None else tuple(args.prefix_order)
    if prefix_order is not None and sorted(prefix_order) != list(OUTCOMES):
        raise ValueError("prefix order must be a permutation of 0,1,2,3")
    model, variables = build(
        effects,
        args.weight,
        prefix_order,
        args.target,
        not args.no_rotation_gauge,
        None if not args.linked_column else tuple(dict.fromkeys(args.linked_column)),
        args.require_cp_completion,
        args.ando_directions,
        None
        if args.fourier_trace_branches is None
        else tuple(args.fourier_trace_branches),
        None
        if args.prefix_prior_bounds is None
        else np.asarray(args.prefix_prior_bounds, dtype=float).reshape(4, 2),
        args.basis_determinant_sign,
        args.basis_determinant_floor,
        None
        if args.basis_choi_witnesses_npz is None
        else np.load(args.basis_choi_witnesses_npz)["witnesses"],
    )
    if not 1e-9 <= args.feasibility_tolerance <= 1e-3:
        raise ValueError("feasibility tolerance must lie in [1e-9,1e-3]")
    model.setRealParam("limits/time", args.seconds)
    model.setRealParam("limits/gap", args.gap)
    model.setRealParam("numerics/feastol", args.feasibility_tolerance)
    model.setRealParam("numerics/dualfeastol", args.feasibility_tolerance)
    if args.seed_npz is not None and args.common_instrument_seed_npz is not None:
        raise ValueError("choose only one seed format")
    seed_accepted = None
    if args.seed_npz is not None:
        if args.require_cp_completion:
            raise ValueError(
                "legacy joint-effect seeds do not contain CP-completion factors"
            )
        seed_solution(model, variables, args.seed_npz, effects, args.weight)
    if args.common_instrument_seed_npz is not None:
        if not args.require_cp_completion:
            raise ValueError("common-instrument seeds require CP completion")
        seed_accepted = seed_from_common_instrument(
            model,
            variables,
            args.common_instrument_seed_npz,
            effects,
            args.weight,
        )
    model.setIntParam("display/verblevel", 2)
    model.optimize()

    payload = {
        "weight": args.weight,
        "terminal_effect_weights": weights.tolist(),
        "prefix_order": None if prefix_order is None else list(prefix_order),
        "target": args.target,
        "rotation_gauge_fixed": not args.no_rotation_gauge,
        "cp_completion_required": args.require_cp_completion,
        "ando_directions": args.ando_directions,
        "fourier_trace_branches": args.fourier_trace_branches,
        "prefix_prior_bounds": (
            None
            if args.prefix_prior_bounds is None
            else np.asarray(args.prefix_prior_bounds, dtype=float)
            .reshape(4, 2)
            .tolist()
        ),
        "basis_determinant_sign": args.basis_determinant_sign,
        "basis_determinant_floor": args.basis_determinant_floor,
        "feasibility_tolerance": args.feasibility_tolerance,
        "basis_choi_witnesses": (
            None
            if args.basis_choi_witnesses_npz is None
            else str(args.basis_choi_witnesses_npz)
        ),
        "common_instrument_seed": (
            None
            if args.common_instrument_seed_npz is None
            else str(args.common_instrument_seed_npz)
        ),
        "seed_accepted": seed_accepted,
        "linked_columns": (
            "all_b"
            if not args.linked_column
            else list(dict.fromkeys(args.linked_column))
        ),
        "status": str(model.getStatus()),
        "primal_bound": float(model.getPrimalbound()),
        "dual_bound": float(model.getDualbound()),
        "absolute_gap": float(model.getDualbound() - model.getPrimalbound()),
        "relative_gap": float(model.getGap()),
        "nodes": int(model.getNNodes()),
        "solving_time": float(model.getSolvingTime()),
    }
    solution = model.getBestSol()
    if solution is not None:
        payload["incumbent_score"] = float(
            model.getSolVal(solution, variables["score"])
        )
        payload["incumbent_audit"] = float(
            model.getSolVal(solution, variables["audit"])
        )
        payload["incumbent_return"] = float(
            model.getSolVal(solution, variables["return"])
        )
    if args.require_cp_completion:
        extracted = extract_common_instrument(model, variables, effects, args.weight)
        if extracted is not None:
            report, arrays = extracted
            payload["reconstructed_common_instrument"] = report
            if args.solution_npz is not None:
                args.solution_npz.parent.mkdir(parents=True, exist_ok=True)
                np.savez(args.solution_npz, **arrays)
    rendered = json.dumps(payload, indent=2) + "\n"
    print(rendered)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

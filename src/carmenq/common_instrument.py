"""Compatibility diagnostics for a common flagged quantum instrument.

For prefix states ``rho[z]`` and conditioned outputs ``sigma[z, y]``, common
instrument compatibility means that one collection of completely positive
maps ``Phi[y]`` satisfies

    sigma[z, y] = Phi[y](rho[z])

for every input label ``z`` and outcome ``y``, while ``sum_y Phi[y]`` is trace
preserving.  This module supplies two complementary numerical checks:

* flagged trace-norm data-processing inequalities, which are inexpensive
  necessary conditions; and
* an SDP projection onto the exact common-instrument output set for fixed
  prefix states.

The SDP routines import CVXPY lazily because it is an optional ``frontier``
dependency rather than a requirement for the core CARMEN-Q package.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable

import numpy as np


Array = np.ndarray


@dataclass(frozen=True)
class FlaggedTraceNormCut:
    """One flagged-channel data-processing comparison.

    ``slack`` is the input trace norm minus the sum of the conditioned output
    trace norms.  A negative value is a common-instrument obstruction.
    """

    input_pair: tuple[int, int]
    scale: float
    input_trace_norm: float
    output_trace_norm: float

    @property
    def slack(self) -> float:
        return self.input_trace_norm - self.output_trace_norm

    @property
    def violation(self) -> float:
        return max(0.0, -self.slack)


@dataclass(frozen=True)
class CommonInstrumentProjection:
    """Numerical projection and separating witness for fixed prefix states.

    The witness is normalized in the direct-sum Frobenius norm.  Its reported
    ``separation_gap`` compares the candidate value with an independently
    optimized support value over the exact Choi-compatible set.  The result is
    solver-conditional, not an interval-arithmetic certificate.
    """

    status: str
    solver: str
    distance: float
    squared_distance: float
    choi_matrices: Array
    projected_outputs: Array
    witness: Array
    candidate_witness_value: float
    compatible_support_value: float
    separation_gap: float
    trace_preservation_residual: float
    minimum_choi_eigenvalue: float

    def is_compatible(self, tolerance: float = 1e-7) -> bool:
        return self.distance <= tolerance and self.separation_gap <= tolerance

    @property
    def input_lipschitz_constants(self) -> Array:
        """Return ``max_y ||W[z,y]||_infinity`` for each input label."""
        return common_instrument_witness_lipschitz_constants(self.witness)

    @property
    def uniform_input_radius_budget(self) -> float:
        """Equal per-input trace-radius that preserves positive separation."""
        denominator = float(self.input_lipschitz_constants.sum())
        if denominator <= 1e-15:
            return math.inf
        return max(0.0, self.separation_gap) / denominator


@dataclass(frozen=True)
class BasisInstrumentReconstruction:
    """Exact linear reconstruction from four independent qubit inputs.

    ``transfer_matrices[y, nu, mu]`` is the Pauli-transfer coefficient in
    ``Phi_y(sigma_mu) = sum_nu L[y,nu,mu] sigma_nu``.  The reconstruction is
    algebraic; only the reported floating-point residuals are numerical.
    """

    input_pauli_matrix: Array
    input_determinant: float
    input_condition_number: float
    transfer_matrices: Array
    choi_matrices: Array
    signed_choi_numerators: Array
    minimum_choi_eigenvalues: Array
    output_residual: float
    trace_preservation_residual: float

    @property
    def minimum_choi_eigenvalue(self) -> float:
        return float(np.min(self.minimum_choi_eigenvalues))

    def is_compatible(self, tolerance: float = 1e-8) -> bool:
        """Return whether the unique reconstructed maps form an instrument."""

        return (
            self.minimum_choi_eigenvalue >= -tolerance
            and self.output_residual <= tolerance
            and self.trace_preservation_residual <= tolerance
        )


@dataclass(frozen=True)
class BasisPovmReconstruction:
    """Unique effective POVM reconstructed from four qubit input states.

    ``effect_pauli_coordinates[k]`` contains ``(a0, ax, ay, az)`` for
    ``F[k] = a0 I + ax X + ay Y + az Z``.  A common instrument followed by
    any terminal POVM necessarily induces such a common effective POVM.
    """

    input_pauli_matrix: Array
    input_determinant: float
    input_condition_number: float
    effect_pauli_coordinates: Array
    effect_matrices: Array
    signed_effect_numerators: Array
    minimum_effect_eigenvalues: Array
    probability_residual: float
    completeness_residual: float

    @property
    def minimum_effect_eigenvalue(self) -> float:
        return float(np.min(self.minimum_effect_eigenvalues))

    def is_compatible(self, tolerance: float = 1e-8) -> bool:
        """Return whether the unique effects form a POVM."""

        return (
            self.minimum_effect_eigenvalue >= -tolerance
            and self.probability_residual <= tolerance
            and self.completeness_residual <= tolerance
        )


def _hermitian(matrix: Array) -> Array:
    return 0.5 * (matrix + matrix.conj().T)


_PAULIS = np.asarray(
    [
        [[1.0, 0.0], [0.0, 1.0]],
        [[0.0, 1.0], [1.0, 0.0]],
        [[0.0, -1j], [1j, 0.0]],
        [[1.0, 0.0], [0.0, -1.0]],
    ],
    dtype=complex,
)


def _pauli_coefficients(matrix: Array) -> Array:
    return np.asarray(
        [float(np.trace(matrix @ pauli).real) for pauli in _PAULIS]
    )


def _adjugate(matrix: Array) -> Array:
    """Return the adjugate using minors, without calling an inverse."""

    value = np.asarray(matrix)
    if value.shape != (4, 4):
        raise ValueError("the qubit operator-basis matrix must be 4 by 4")
    cofactors = np.empty((4, 4), dtype=value.dtype)
    for row in range(4):
        for column in range(4):
            minor = np.delete(np.delete(value, row, axis=0), column, axis=1)
            cofactors[row, column] = (-1.0) ** (row + column) * np.linalg.det(minor)
    return cofactors.T


def _choi_from_pauli_transfer(transfer: Array) -> Array:
    value = np.asarray(transfer, dtype=float)
    if value.shape != (4, 4):
        raise ValueError("a qubit Pauli-transfer matrix must be 4 by 4")
    choi = sum(
        0.5
        * value[nu, mu]
        * np.kron(_PAULIS[mu].T, _PAULIS[nu])
        for mu in range(4)
        for nu in range(4)
    )
    return _hermitian(choi)


def _validate_families(prefix_states: Array, conditioned_outputs: Array) -> tuple[Array, Array]:
    states = np.asarray(prefix_states, dtype=complex)
    outputs = np.asarray(conditioned_outputs, dtype=complex)
    if states.ndim != 3 or states.shape[1:] != (2, 2):
        raise ValueError("prefix_states must have shape (inputs, 2, 2)")
    if outputs.ndim != 4 or outputs.shape[0] != states.shape[0] or outputs.shape[2:] != (2, 2):
        raise ValueError(
            "conditioned_outputs must have shape (inputs, outcomes, 2, 2)"
        )
    if outputs.shape[1] < 1:
        raise ValueError("at least one instrument outcome is required")
    states = np.asarray([_hermitian(item) for item in states])
    outputs = np.asarray(
        [[_hermitian(item) for item in row] for row in outputs]
    )
    return states, outputs


def trace_norm_hermitian(matrix: Array) -> float:
    """Return the Schatten one-norm of a Hermitian matrix."""
    values = np.linalg.eigvalsh(_hermitian(np.asarray(matrix, dtype=complex)))
    return float(np.abs(values).sum())


def common_instrument_witness_lipschitz_constants(witness: Array) -> Array:
    r"""Return inputwise continuity constants for a flagged witness.

    If ``W[z,y]`` is a Hermitian witness, then for every trace-preserving
    instrument ``Phi`` and Hermitian perturbation ``Delta[z]``,

    .. math::

       \left|\sum_y \operatorname{Tr}
       W_{zy}\Phi_y(\Delta_z)\right|
       \leq M_z\|\Delta_z\|_1,

    where ``M[z] = max_y ||W[z,y]||_infinity``.
    """
    value = np.asarray(witness, dtype=complex)
    if value.ndim != 4 or value.shape[2:] != (2, 2):
        raise ValueError("witness must have shape (inputs, outcomes, 2, 2)")
    return np.asarray(
        [
            max(
                float(np.max(np.abs(np.linalg.eigvalsh(_hermitian(item)))))
                for item in value[z]
            )
            for z in range(value.shape[0])
        ]
    )


def robust_common_instrument_witness_bound(
    reference_support: float,
    witness: Array,
    reference_states: Array,
    perturbed_states: Array,
) -> float:
    """Extend a fixed-input witness support bound to new prefix states.

    The returned value is valid for every common instrument at
    ``perturbed_states`` whenever ``reference_support`` upper-bounds the same
    witness at ``reference_states``.  It is useful as a branch-local Benders
    cut: trace-radius bounds for a state box may replace the explicit matrix
    distances used here.
    """
    reference = np.asarray(reference_states, dtype=complex)
    perturbed = np.asarray(perturbed_states, dtype=complex)
    if reference.shape != perturbed.shape or reference.ndim != 3 or reference.shape[1:] != (2, 2):
        raise ValueError("state families must share shape (inputs, 2, 2)")
    constants = common_instrument_witness_lipschitz_constants(witness)
    if constants.shape != (reference.shape[0],):
        raise ValueError("witness and state families use different input counts")
    correction = sum(
        constants[z]
        * trace_norm_hermitian(perturbed[z] - reference[z])
        for z in range(reference.shape[0])
    )
    return float(reference_support + correction)


def choi_from_kraus(kraus_operators: Iterable[Array]) -> Array:
    """Construct an input-major qubit Choi matrix."""
    result = np.zeros((4, 4), dtype=complex)
    for operator in kraus_operators:
        item = np.asarray(operator, dtype=complex)
        if item.shape != (2, 2):
            raise ValueError("every Kraus operator must be 2 by 2")
        vector = item.T.reshape(4)
        result += np.outer(vector, vector.conj())
    return _hermitian(result)


def apply_choi(choi: Array, state: Array) -> Array:
    """Apply an input-major qubit Choi matrix to one state."""
    matrix = np.asarray(choi, dtype=complex)
    rho = np.asarray(state, dtype=complex)
    if matrix.shape != (4, 4) or rho.shape != (2, 2):
        raise ValueError("expected Choi shape (4, 4) and state shape (2, 2)")
    blocks = matrix.reshape(2, 2, 2, 2)
    return _hermitian(np.einsum("ij,iajb->ab", rho, blocks))


def conditioned_outputs(prefix_states: Array, choi_matrices: Array) -> Array:
    """Evaluate every prefix state through every outcome of one instrument."""
    states = np.asarray(prefix_states, dtype=complex)
    choi = np.asarray(choi_matrices, dtype=complex)
    if states.ndim != 3 or states.shape[1:] != (2, 2):
        raise ValueError("prefix_states must have shape (inputs, 2, 2)")
    if choi.ndim != 3 or choi.shape[1:] != (4, 4):
        raise ValueError("choi_matrices must have shape (outcomes, 4, 4)")
    return np.asarray(
        [[apply_choi(item, state) for item in choi] for state in states]
    )


def reconstruct_common_instrument_from_basis(
    prefix_states: Array,
    conditioned_output_states: Array,
    rank_tolerance: float = 1e-10,
) -> BasisInstrumentReconstruction:
    r"""Reconstruct the unique instrument maps fixed by four qubit inputs.

    Four linearly independent Hermitian qubit operators form a basis of the
    real vector space of Hermitian matrices.  Their images therefore determine
    each Hermiticity-preserving map uniquely.  The supplied family comes from
    one quantum instrument exactly when every reconstructed Choi matrix is
    positive semidefinite and their sum is trace preserving.

    The polynomial numerator reported in the result avoids division by the
    basis determinant.  If ``delta`` is that determinant and ``K_y`` the
    numerator, then

    .. math::

       |\delta| J_y = \operatorname{sign}(\delta) K_y.

    This form is useful for branch-and-bound or sum-of-squares models.  Singular
    input families require the general fixed-input SDP instead.
    """

    if not math.isfinite(rank_tolerance) or rank_tolerance <= 0.0:
        raise ValueError("rank_tolerance must be finite and positive")
    states, outputs = _validate_families(
        prefix_states, conditioned_output_states
    )
    if states.shape[0] != 4:
        raise ValueError("basis reconstruction requires exactly four inputs")
    input_pauli = np.asarray([_pauli_coefficients(state) for state in states])
    singular_values = np.linalg.svd(input_pauli, compute_uv=False)
    if singular_values[-1] <= rank_tolerance * singular_values[0]:
        raise np.linalg.LinAlgError(
            "the four prefix states do not span the Hermitian qubit operators"
        )
    determinant = float(np.linalg.det(input_pauli))
    condition_number = float(singular_values[0] / singular_values[-1])
    output_pauli = np.asarray(
        [
            [_pauli_coefficients(outputs[z, y]) for z in range(4)]
            for y in range(outputs.shape[1])
        ]
    )
    transfer = np.asarray(
        [np.linalg.solve(input_pauli, output_pauli[y]).T for y in range(outputs.shape[1])]
    )
    choi = np.asarray([_choi_from_pauli_transfer(item) for item in transfer])

    # L_y = Q_y^T adj(R^T) / det(R).  Forming the numerator by minors makes
    # the exact polynomial matrix inequality visible to symbolic optimizers.
    adjugate_transpose = _adjugate(input_pauli.T)
    transfer_numerators = np.asarray(
        [output_pauli[y].T @ adjugate_transpose for y in range(outputs.shape[1])]
    )
    choi_numerators = np.asarray(
        [_choi_from_pauli_transfer(item) for item in transfer_numerators]
    )
    signed_numerators = np.sign(determinant) * choi_numerators
    reconstructed = conditioned_outputs(states, choi)
    output_residual = float(np.linalg.norm(reconstructed - outputs))
    total_choi = choi.sum(axis=0).reshape(2, 2, 2, 2)
    partial_trace = np.einsum("iaja->ij", total_choi)
    trace_residual = float(np.linalg.norm(partial_trace - np.eye(2)))
    minimum_eigenvalues = np.asarray(
        [float(np.linalg.eigvalsh(item).min()) for item in choi]
    )
    return BasisInstrumentReconstruction(
        input_pauli_matrix=input_pauli,
        input_determinant=determinant,
        input_condition_number=condition_number,
        transfer_matrices=transfer,
        choi_matrices=choi,
        signed_choi_numerators=signed_numerators,
        minimum_choi_eigenvalues=minimum_eigenvalues,
        output_residual=output_residual,
        trace_preservation_residual=trace_residual,
    )


def reconstruct_effective_povm_from_basis(
    prefix_states: Array,
    probabilities: Array,
    rank_tolerance: float = 1e-10,
) -> BasisPovmReconstruction:
    r"""Reconstruct the unique common POVM fixed by four qubit inputs.

    Let ``R[z]`` contain the Pauli coordinates of the four input states and
    let ``Q[z,k]`` be the observed outcome probabilities.  When ``R`` is
    nonsingular, the unique effect coordinates are ``A = R^{-1} Q``.  They
    form a POVM exactly when every reconstructed effect is positive and the
    effects sum to identity.

    With ``delta = det(R)``, the returned signed numerators equal
    ``sign(delta) * adj(R) Q`` in matrix form.  Their positivity is the finite
    determinant-scaled polynomial condition useful in global optimisation.
    A common instrument followed by a terminal measurement must pass this
    test, although a common effective POVM need not admit the specified
    sequential factorisation.
    """

    if not math.isfinite(rank_tolerance) or rank_tolerance <= 0.0:
        raise ValueError("rank_tolerance must be finite and positive")
    states = np.asarray(prefix_states, dtype=complex)
    if states.shape != (4, 2, 2):
        raise ValueError("basis POVM reconstruction requires states of shape (4,2,2)")
    states = np.asarray([_hermitian(item) for item in states])
    raw_probabilities = np.asarray(probabilities)
    if raw_probabilities.ndim != 2 or raw_probabilities.shape[0] != 4:
        raise ValueError("probabilities must have shape (4, outcomes)")
    if raw_probabilities.shape[1] < 1:
        raise ValueError("at least one POVM outcome is required")
    if np.iscomplexobj(raw_probabilities) and np.max(np.abs(raw_probabilities.imag)) > 1e-12:
        raise ValueError("probabilities must be real")
    probability = np.asarray(raw_probabilities.real, dtype=float)
    if not np.all(np.isfinite(probability)):
        raise ValueError("probabilities must be finite")

    input_pauli = np.asarray([_pauli_coefficients(state) for state in states])
    singular_values = np.linalg.svd(input_pauli, compute_uv=False)
    if singular_values[-1] <= rank_tolerance * singular_values[0]:
        raise np.linalg.LinAlgError(
            "the four prefix states do not span the Hermitian qubit operators"
        )
    determinant = float(np.linalg.det(input_pauli))
    condition_number = float(singular_values[0] / singular_values[-1])
    coordinates = np.linalg.solve(input_pauli, probability).T
    effects = np.asarray(
        [
            _hermitian(
                sum(
                    coordinate[mu] * _PAULIS[mu]
                    for mu in range(4)
                )
            )
            for coordinate in coordinates
        ]
    )
    signed_numerator_coordinates = abs(determinant) * coordinates
    signed_numerators = np.asarray(
        [
            _hermitian(
                sum(
                    coordinate[mu] * _PAULIS[mu]
                    for mu in range(4)
                )
            )
            for coordinate in signed_numerator_coordinates
        ]
    )
    reconstructed = input_pauli @ coordinates.T
    probability_residual = float(np.linalg.norm(reconstructed - probability))
    completeness_residual = float(np.linalg.norm(effects.sum(axis=0) - np.eye(2)))
    minimum_eigenvalues = np.asarray(
        [float(np.linalg.eigvalsh(effect).min()) for effect in effects]
    )
    return BasisPovmReconstruction(
        input_pauli_matrix=input_pauli,
        input_determinant=determinant,
        input_condition_number=condition_number,
        effect_pauli_coordinates=coordinates,
        effect_matrices=effects,
        signed_effect_numerators=signed_numerators,
        minimum_effect_eigenvalues=minimum_eigenvalues,
        probability_residual=probability_residual,
        completeness_residual=completeness_residual,
    )


def flagged_trace_norm_cut(
    prefix_states: Array,
    conditioned_output_states: Array,
    first_input: int,
    second_input: int,
    scale: float = 1.0,
) -> FlaggedTraceNormCut:
    r"""Evaluate one necessary common-instrument inequality.

    The inequality is

    .. math::

       \sum_y \|\sigma_{zy}-t\sigma_{z'y}\|_1
       \leq \|\rho_z-t\rho_{z'}\|_1.
    """
    if not math.isfinite(scale) or scale < 0.0:
        raise ValueError("scale must be finite and nonnegative")
    states, outputs = _validate_families(
        prefix_states, conditioned_output_states
    )
    if not 0 <= first_input < states.shape[0] or not 0 <= second_input < states.shape[0]:
        raise IndexError("input label outside the supplied family")
    if first_input == second_input:
        raise ValueError("the two input labels must differ")
    input_norm = trace_norm_hermitian(
        states[first_input] - scale * states[second_input]
    )
    output_norm = sum(
        trace_norm_hermitian(
            outputs[first_input, outcome]
            - scale * outputs[second_input, outcome]
        )
        for outcome in range(outputs.shape[1])
    )
    return FlaggedTraceNormCut(
        input_pair=(first_input, second_input),
        scale=float(scale),
        input_trace_norm=input_norm,
        output_trace_norm=output_norm,
    )


def comparison_scale_grid(
    minimum: float = 1e-3,
    maximum: float = 1e3,
    logarithmic_samples: int = 25,
) -> Array:
    """Return a deterministic scale grid containing ``0`` and ``1``."""
    if minimum <= 0.0 or maximum < minimum:
        raise ValueError("require 0 < minimum <= maximum")
    if logarithmic_samples < 2:
        raise ValueError("logarithmic_samples must be at least two")
    values = np.concatenate(
        (
            np.asarray([0.0, 1.0]),
            np.geomspace(minimum, maximum, logarithmic_samples),
        )
    )
    return np.unique(values)


def scan_flagged_trace_norm_cuts(
    prefix_states: Array,
    conditioned_output_states: Array,
    scales: Iterable[float] | None = None,
) -> tuple[FlaggedTraceNormCut, ...]:
    """Evaluate all unordered input pairs on a scale grid."""
    states, outputs = _validate_families(
        prefix_states, conditioned_output_states
    )
    grid = comparison_scale_grid() if scales is None else np.asarray(tuple(scales), dtype=float)
    cuts = []
    for first in range(states.shape[0]):
        for second in range(first + 1, states.shape[0]):
            for scale in grid:
                cuts.append(
                    flagged_trace_norm_cut(
                        states, outputs, first, second, float(scale)
                    )
                )
    return tuple(cuts)


def _cvxpy_output_expression(cp: object, choi: object, state: Array) -> object:
    entries = [
        [
            sum(
                state[i, j] * choi[2 * i + a, 2 * j + b]
                for i in range(2)
                for j in range(2)
            )
            for b in range(2)
        ]
        for a in range(2)
    ]
    output = cp.bmat(entries)
    return 0.5 * (output + output.H)


def _cvxpy_instrument(cp: object, outcomes: int) -> tuple[list[object], list[object]]:
    variables = [cp.Variable((4, 4), hermitian=True) for _ in range(outcomes)]
    constraints = [item >> 0 for item in variables]
    for i in range(2):
        for j in range(2):
            partial_trace = sum(
                variables[y][2 * i, 2 * j]
                + variables[y][2 * i + 1, 2 * j + 1]
                for y in range(outcomes)
            )
            constraints.append(partial_trace == (1.0 if i == j else 0.0))
    return variables, constraints


def _solve_cvxpy(problem: object, cp: object, solver: str | None) -> str:
    chosen = solver
    if chosen is None:
        installed = set(cp.installed_solvers())
        chosen = "CLARABEL" if "CLARABEL" in installed else "SCS"
    chosen = chosen.upper()
    if chosen == "CLARABEL":
        problem.solve(
            solver=chosen,
            tol_gap_abs=1e-9,
            tol_gap_rel=1e-9,
            tol_feas=1e-9,
            max_iter=1000,
        )
    elif chosen == "SCS":
        problem.solve(solver=chosen, eps=2e-7, max_iters=300_000)
    else:
        problem.solve(solver=chosen)
    if problem.status not in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}:
        raise RuntimeError(f"common-instrument SDP failed: {problem.status}")
    return chosen


def project_to_common_instrument(
    prefix_states: Array,
    candidate_outputs: Array,
    solver: str | None = None,
) -> CommonInstrumentProjection:
    """Project outputs onto those generated by one common instrument.

    Prefix states are fixed.  The first SDP finds the closest compatible
    output family in Frobenius norm.  If the residual is nonzero, a second
    SDP maximizes the induced linear witness over the compatible set.  A
    positive ``separation_gap`` numerically certifies incompatibility.
    """
    try:
        import cvxpy as cp
    except ImportError as error:  # pragma: no cover - depends on extras
        raise ImportError(
            "project_to_common_instrument requires the 'frontier' extra"
        ) from error

    states, candidate = _validate_families(prefix_states, candidate_outputs)
    outcomes = candidate.shape[1]
    variables, constraints = _cvxpy_instrument(cp, outcomes)
    expressions: list[list[object]] = []
    residual_entries = []
    for z in range(states.shape[0]):
        row = []
        for y in range(outcomes):
            output = _cvxpy_output_expression(cp, variables[y], states[z])
            row.append(output)
            difference = output - candidate[z, y]
            for a in range(2):
                for b in range(2):
                    residual_entries.extend(
                        (cp.real(difference[a, b]), cp.imag(difference[a, b]))
                    )
        expressions.append(row)
    residual_vector = cp.hstack(residual_entries)
    projection_problem = cp.Problem(
        cp.Minimize(cp.norm(residual_vector, 2)), constraints
    )
    chosen = _solve_cvxpy(projection_problem, cp, solver)
    choi_value = np.asarray([np.asarray(item.value) for item in variables])
    projected = conditioned_outputs(states, choi_value)
    raw_witness = candidate - projected
    distance = float(np.linalg.norm(raw_witness.reshape(-1)))
    witness = (
        np.zeros_like(raw_witness)
        if distance <= 1e-14
        else raw_witness / distance
    )

    candidate_value = float(
        sum(
            np.trace(witness[z, y].conj().T @ candidate[z, y]).real
            for z in range(states.shape[0])
            for y in range(outcomes)
        )
    )
    if distance <= 1e-14:
        support_value = candidate_value
    else:
        support_variables, support_constraints = _cvxpy_instrument(cp, outcomes)
        support_terms = []
        for z in range(states.shape[0]):
            for y in range(outcomes):
                output = _cvxpy_output_expression(
                    cp, support_variables[y], states[z]
                )
                support_terms.append(
                    cp.real(cp.trace(witness[z, y].conj().T @ output))
                )
        support_problem = cp.Problem(cp.Maximize(cp.sum(cp.hstack(support_terms))), support_constraints)
        support_solver = _solve_cvxpy(support_problem, cp, chosen)
        if support_solver != chosen:  # pragma: no cover - defensive
            raise RuntimeError("projection and support solvers unexpectedly differ")
        support_value = float(support_problem.value)

    partial_trace = np.zeros((2, 2), dtype=complex)
    for matrix in choi_value:
        blocks = matrix.reshape(2, 2, 2, 2)
        partial_trace += np.einsum("iaja->ij", blocks)
    tp_residual = float(np.linalg.norm(partial_trace - np.eye(2)))
    minimum_eigenvalue = min(
        float(np.linalg.eigvalsh(_hermitian(item)).min())
        for item in choi_value
    )
    return CommonInstrumentProjection(
        status=str(projection_problem.status),
        solver=chosen,
        distance=distance,
        squared_distance=distance * distance,
        choi_matrices=choi_value,
        projected_outputs=projected,
        witness=witness,
        candidate_witness_value=candidate_value,
        compatible_support_value=support_value,
        separation_gap=candidate_value - support_value,
        trace_preservation_residual=tp_residual,
        minimum_choi_eigenvalue=minimum_eigenvalue,
    )

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


def _hermitian(matrix: Array) -> Array:
    return 0.5 * (matrix + matrix.conj().T)


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

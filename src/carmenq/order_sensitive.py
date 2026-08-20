"""Exact utilities for the canonical order-sensitive AUDIT--RETURN instances.

This module covers binary rank-two syndrome checks streamed through one
persistent qubit.  It exposes the grouped four-slot frontier and the exact
perfect-AUDIT endpoint of its interleaved column permutation.  The latter is
an endpoint theorem only.  An exact two-parameter achievable construction is
also exposed as a lower bound, but the unknown interleaved interior is never
presented as an asserted upper bound or solved frontier.

All ranks and cut profiles are over :math:`GF(2)`.  A cut index ``i`` means
that columns ``[:i]`` have arrived and columns ``[i:]`` remain.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import sqrt
from typing import Final

import numpy as np
from numpy.typing import ArrayLike
from scipy.optimize import minimize


# Immutable tuple representations keep public canonical instances safe from
# accidental in-place mutation while remaining directly accepted by NumPy.
GROUPED_CHECK_MATRIX: Final[tuple[tuple[int, ...], ...]] = (
    (1, 1, 0, 0),
    (0, 0, 1, 1),
)
"""Four-slot grouped check matrix with column order ``A, A, B, B``."""

INTERLEAVED_CHECK_MATRIX: Final[tuple[tuple[int, ...], ...]] = (
    (1, 0, 1, 0),
    (0, 1, 0, 1),
)
"""Four-slot interleaved check matrix with column order ``A, B, A, B``."""


def _binary_integer_matrix(matrix: ArrayLike) -> np.ndarray:
    """Return a two-dimensional integer array reduced modulo two."""
    values = np.asarray(matrix)
    if values.ndim != 2:
        raise ValueError("matrix must be two-dimensional")
    if not (
        np.issubdtype(values.dtype, np.integer)
        or np.issubdtype(values.dtype, np.bool_)
    ):
        raise TypeError("matrix entries must be integers or booleans")
    return np.asarray(values, dtype=np.uint8) & 1


def gf2_rank(matrix: ArrayLike) -> int:
    """Compute the matrix rank over :math:`GF(2)` by row reduction.

    Integer entries are interpreted modulo two.  Empty two-dimensional
    matrices have rank zero.
    """
    reduced = _binary_integer_matrix(matrix).copy()
    n_rows, n_columns = reduced.shape
    rank = 0
    for column in range(n_columns):
        pivot_candidates = np.flatnonzero(reduced[rank:, column])
        if pivot_candidates.size == 0:
            continue
        pivot = rank + int(pivot_candidates[0])
        if pivot != rank:
            reduced[[rank, pivot]] = reduced[[pivot, rank]]
        other_rows = np.flatnonzero(reduced[:, column])
        other_rows = other_rows[other_rows != rank]
        reduced[other_rows] ^= reduced[rank]
        rank += 1
        if rank == n_rows:
            break
    return rank


def trellis_connectivity_profile(matrix: ArrayLike) -> tuple[int, ...]:
    """Return the binary trellis-connectivity exponent at every cut.

    For a check matrix ``H`` with ``n`` columns, entry ``i`` is

    ``rank(H[:, :i]) + rank(H[:, i:]) - rank(H)``.

    The returned tuple includes the two trivial cuts ``i=0`` and ``i=n``.
    This is a structural code/order invariant; by itself it is not a formula
    for the complete AUDIT--RETURN frontier.
    """
    checks = _binary_integer_matrix(matrix)
    total_rank = gf2_rank(checks)
    return tuple(
        gf2_rank(checks[:, :cut])
        + gf2_rank(checks[:, cut:])
        - total_rank
        for cut in range(checks.shape[1] + 1)
    )


def trellis_connectivity_tau(matrix: ArrayLike) -> int:
    """Return the maximum of :func:`trellis_connectivity_profile`."""
    return max(trellis_connectivity_profile(matrix))


def full_crossing_cuts(matrix: ArrayLike) -> tuple[int, ...]:
    """Return interior cuts where both halves retain the full binary rank.

    If ``r = rank(H) > 0``, a cut ``i`` is returned exactly when both
    ``H[:, :i]`` and ``H[:, i:]`` have rank ``r``.  For a rank-two check
    matrix this is the full-crossing condition used by the exact
    perfect-AUDIT qubit endpoint theorem.  The function itself also accepts
    other positive ranks, without asserting that theorem for them.
    """
    checks = _binary_integer_matrix(matrix)
    rank = gf2_rank(checks)
    if rank == 0:
        return ()
    return tuple(
        cut
        for cut in range(1, checks.shape[1])
        if gf2_rank(checks[:, :cut]) == rank
        and gf2_rank(checks[:, cut:]) == rank
    )


def _audit_weight(value: float) -> float:
    weight = float(value)
    if not 0.0 <= weight <= 1.0:
        raise ValueError("audit_weight must lie in [0, 1]")
    return weight


def rank_two_static_qubit_support(audit_weight: float = 0.5) -> float:
    """Return the static four-label/qubit support relaxation.

    The value is

    ``(1 + sqrt(lambda**2 + (1-lambda)**2)) / 2``.

    It is a universal terminal-dimension ceiling for the rank-two syndrome
    game, not a claim that every streamed column order attains the ceiling.
    The grouped order does attain it; the interleaved order does not near the
    perfect-AUDIT endpoint.
    """
    weight = _audit_weight(audit_weight)
    return (1.0 + sqrt(weight**2 + (1.0 - weight) ** 2)) / 2.0


@dataclass(frozen=True)
class GroupedFrontierPoint:
    """One exact exposed point of the grouped one-qubit frontier."""

    audit_weight: float
    weak_measurement_strength: float
    audit_probability: float
    return_fidelity: float
    support_value: float


def grouped_frontier(audit_weight: float = 0.5) -> GroupedFrontierPoint:
    """Return the exact grouped-order support point for one coherent qubit.

    The operational scope is the four-slot complete-syndrome AUDIT and
    all-pair EPR RETURN game with sequestered earlier outputs, unrestricted
    finite classical transcript, and a persistent coherent dimension of two.
    The attainer weakly records the first grouped parity and hands the qubit
    over to the remaining parity.
    """
    weight = _audit_weight(audit_weight)
    norm = sqrt(weight**2 + (1.0 - weight) ** 2)
    strength = weight / norm
    audit_probability = (1.0 + strength) / 2.0
    return_fidelity = (1.0 + sqrt(max(0.0, 1.0 - strength**2))) / 2.0
    support_value = (
        weight * audit_probability + (1.0 - weight) * return_fidelity
    )
    return GroupedFrontierPoint(
        audit_weight=weight,
        weak_measurement_strength=strength,
        audit_probability=audit_probability,
        return_fidelity=return_fidelity,
        support_value=support_value,
    )


@dataclass(frozen=True)
class PerfectAuditEndpoint:
    """Metadata for an exact constrained perfect-AUDIT endpoint."""

    check_matrix: tuple[tuple[int, ...], ...]
    coherent_dimension: int
    audit_probability: float
    maximum_return_fidelity: float
    bound_is_attained: bool
    full_crossing_cuts: tuple[int, ...]
    interior_frontier_known: bool


@dataclass(frozen=True)
class InterleavedCandidatePoint:
    """One achievable point from the analytic interleaved candidate family.

    ``support_is_globally_optimal`` is deliberately false: the construction
    is an exact lower bound, while the arbitrary-instrument interior converse
    remains open.
    """

    audit_weight: float
    strategy: str
    q: float | None
    v: float | None
    audit_probability: float
    return_fidelity: float
    support_value: float
    support_is_globally_optimal: bool = False


INTERLEAVED_PERFECT_AUDIT_ENDPOINT: Final[PerfectAuditEndpoint] = (
    PerfectAuditEndpoint(
        check_matrix=INTERLEAVED_CHECK_MATRIX,
        coherent_dimension=2,
        audit_probability=1.0,
        maximum_return_fidelity=0.25,
        bound_is_attained=True,
        full_crossing_cuts=(2,),
        interior_frontier_known=False,
    )
)
"""Exact interleaved endpoint under the streamed one-qubit protocol.

The statement is ``P_A = 1 => F_R <= 1/4``, with equality attained.  It
covers adaptive non-QND finite-outcome instruments, a genuinely classical
transcript, sequestered earlier outputs, and arbitrary transcript-conditioned
RETURN decoding.  It does not supply the unknown interior frontier.
"""


def interleaved_candidate_scores(q: float, v: float) -> tuple[float, float]:
    """Evaluate the exact two-parameter interleaved construction.

    The returned pair is ``(P_A, F_R)`` with

    ``P_A = 1/2 + q*v*sqrt(1-v**2) - q*(1-q)*v**2``

    and the corresponding flagged polar-recovery fidelity.  The parameters
    must lie in the unit square.  This is a physical streamed one-qubit
    construction, not an asserted characterization of every possible
    interleaved strategy.
    """
    q_value = float(q)
    v_value = float(v)
    if not 0.0 <= q_value <= 1.0:
        raise ValueError("q must lie in [0, 1]")
    if not 0.0 <= v_value <= 1.0:
        raise ValueError("v must lie in [0, 1]")
    audit_probability = (
        0.5
        + q_value * v_value * sqrt(max(0.0, 1.0 - v_value**2))
        - q_value * (1.0 - q_value) * v_value**2
    )
    return_fidelity = 0.25 * (
        sqrt(max(0.0, 1.0 - (1.0 - q_value**2) * v_value**2))
        + v_value
        * (
            1.0
            - q_value
            + 2.0 * sqrt(max(0.0, q_value * (1.0 - q_value)))
        )
    ) ** 2
    return audit_probability, return_fidelity


def interleaved_candidate_lower_bound(
    audit_weight: float = 0.5,
) -> InterleavedCandidatePoint:
    """Optimize the analytic interleaved family as an achievable lower bound.

    A deterministic multistart two-variable search is compared with the exact
    no-record strategy ``(P_A,F_R)=(1/2,1)``.  Numerical optimization affects
    only which explicit achievable parameters are returned.  It does not
    convert the still-open arbitrary-instrument converse into a theorem.
    """
    weight = _audit_weight(audit_weight)
    no_record_score = 1.0 - weight / 2.0

    def negative_score(point: np.ndarray) -> float:
        audit_probability, return_fidelity = interleaved_candidate_scores(
            float(point[0]), float(point[1])
        )
        return -(
            weight * audit_probability + (1.0 - weight) * return_fidelity
        )

    starts = (
        (0.15, 0.25),
        (0.35, 0.60),
        (0.58, 0.81),
        (0.75, 0.77),
        (0.90, 0.73),
        (0.99, 1.0 / sqrt(2.0)),
    )
    best_result = min(
        (
            minimize(
                negative_score,
                np.asarray(start, dtype=float),
                method="L-BFGS-B",
                bounds=((0.0, 1.0), (0.0, 1.0)),
                options={"ftol": 1e-15, "gtol": 1e-11, "maxiter": 2000},
            )
            for start in starts
        ),
        key=lambda result: float(result.fun),
    )
    q_value, v_value = map(float, best_result.x)
    audit_probability, return_fidelity = interleaved_candidate_scores(
        q_value, v_value
    )
    candidate_score = (
        weight * audit_probability + (1.0 - weight) * return_fidelity
    )
    if no_record_score >= candidate_score - 5e-13:
        return InterleavedCandidatePoint(
            audit_weight=weight,
            strategy="no_record",
            q=None,
            v=None,
            audit_probability=0.5,
            return_fidelity=1.0,
            support_value=no_record_score,
        )
    return InterleavedCandidatePoint(
        audit_weight=weight,
        strategy="two_parameter",
        q=q_value,
        v=v_value,
        audit_probability=audit_probability,
        return_fidelity=return_fidelity,
        support_value=candidate_score,
    )


__all__ = [
    "GROUPED_CHECK_MATRIX",
    "INTERLEAVED_CHECK_MATRIX",
    "INTERLEAVED_PERFECT_AUDIT_ENDPOINT",
    "GroupedFrontierPoint",
    "InterleavedCandidatePoint",
    "PerfectAuditEndpoint",
    "full_crossing_cuts",
    "gf2_rank",
    "grouped_frontier",
    "interleaved_candidate_lower_bound",
    "interleaved_candidate_scores",
    "rank_two_static_qubit_support",
    "trellis_connectivity_profile",
    "trellis_connectivity_tau",
]

"""Exact utilities for the canonical order-sensitive AUDIT--RETURN instances.

This module covers binary rank-two syndrome checks streamed through one
persistent qubit and the finite-field temporal power law.  It exposes the
grouped four-slot frontier, the exact perfect-AUDIT endpoint of its interleaved
column permutation, and a rigorous approximate-AUDIT upper certificate based
on causal list decoding.  Two exact achievable constructions are exposed as
lower bounds.  The stronger one follows from a compact three-effect Choi-MPS
family and matches unrestricted variational searches, but the still-unknown
global interleaved converse is never presented as solved.

All ranks and cut profiles are over :math:`GF(2)`.  A cut index ``i`` means
that columns ``[:i]`` have arrived and columns ``[i:]`` remain.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import sqrt
from operator import index
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

INTERLEAVED_ORDER_GAP_WEIGHT_THRESHOLD: Final[float] = 3.0 / 7.0
"""Lowest audit weight certified by the linear-tail order-gap theorem.

The certificate is strict for weights in ``(3/7, 1)``.  It is not a claim
that the true order gap begins only at this threshold.
"""


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


def full_rank_block_packing_number(matrix: ArrayLike) -> int:
    """Return the maximum number of consecutive full-rank binary blocks.

    For a positive-rank ordered check matrix ``H``, this is the largest ``m``
    for which the full column sequence can be partitioned into ``m`` nonempty
    consecutive blocks, each having rank ``rank(H)``.  The greedy algorithm is
    optimal: close each block at the earliest column where full rank is
    reached, and append any final rank-deficient tail to the last block.

    A rank-zero matrix returns zero.  The quantity is a structural descriptor,
    not a claim of a new classical trellis invariant.
    """
    checks = _binary_integer_matrix(matrix)
    rank = gf2_rank(checks)
    if rank == 0:
        return 0
    blocks = 0
    start = 0
    for stop in range(1, checks.shape[1] + 1):
        if gf2_rank(checks[:, start:stop]) == rank:
            blocks += 1
            start = stop
    return blocks


def _positive_integer(value: int, name: str) -> int:
    """Return an integer-like value after enforcing strict positivity."""
    try:
        result = index(value)
    except TypeError as error:
        raise TypeError(f"{name} must be an integer") from error
    if result <= 0:
        raise ValueError(f"{name} must be positive")
    return result


def _is_prime_power(value: int) -> bool:
    """Return whether ``value`` is the order of a finite field."""
    divisor = 2
    while divisor * divisor <= value and value % divisor:
        divisor += 1
    if divisor * divisor > value:
        return True
    remainder = value
    while remainder % divisor == 0:
        remainder //= divisor
    return remainder == 1


def full_crossing_perfect_audit_return_bound(
    syndrome_rank: int,
    coherent_dimension: int,
    alphabet_size: int = 2,
) -> float:
    """Bound perfect-AUDIT EPR return at one full-crossing syndrome cut.

    Consider a rank-``syndrome_rank`` linear check over a finite field with
    prime-power ``alphabet_size`` elements.  At the declared cut, both the arrived and
    remaining column blocks must retain the full check rank.  If at most
    ``coherent_dimension`` coherent dimensions cross the cut and reach the
    terminal decoder, perfect full-syndrome AUDIT implies

    ``F_R <= min(1, (d/N)**2)``,

    where ``N = alphabet_size**syndrome_rank``.  Unlimited finite genuinely
    classical transcript is allowed.  The function returns an upper bound;
    it does not claim attainability for every matrix or dimension.

    For the canonical binary rank-two/qubit interleaved instance, the value
    is exactly ``1/4`` and is attained.
    """
    return full_rank_block_perfect_audit_return_bound(
        syndrome_rank=syndrome_rank,
        coherent_dimension=coherent_dimension,
        block_count=2,
        alphabet_size=alphabet_size,
    )


def full_rank_block_perfect_audit_return_bound(
    syndrome_rank: int,
    coherent_dimension: int,
    block_count: int,
    alphabet_size: int = 2,
) -> float:
    """Return the perfect-AUDIT temporal power-law upper bound.

    The ordered check matrix must admit ``block_count`` consecutive blocks,
    each retaining the full ``syndrome_rank`` over a finite field of
    prime-power size ``alphabet_size``.  If at most ``coherent_dimension``
    coherent dimensions cross every block boundary and reach terminal AUDIT,
    then

    ``F_R <= min(1, (d/N)**block_count)``,

    where ``N = alphabet_size**syndrome_rank``.  Unlimited finite genuinely
    classical transcript is allowed.  The theorem is tight on repeated
    identity blocks when ``d`` is a power of the alphabet size, but this
    function does not assert attainability for every matrix.
    """
    rank = _positive_integer(syndrome_rank, "syndrome_rank")
    dimension = _positive_integer(coherent_dimension, "coherent_dimension")
    blocks = _positive_integer(block_count, "block_count")
    alphabet = _positive_integer(alphabet_size, "alphabet_size")
    if alphabet < 2:
        raise ValueError("alphabet_size must be at least two")
    if not _is_prime_power(alphabet):
        raise ValueError("alphabet_size must be a prime power")
    syndrome_count = alphabet**rank
    if dimension >= syndrome_count:
        return 1.0
    return (dimension / syndrome_count) ** blocks


def full_rank_block_approximate_audit_return_bound(
    audit_probability: float,
    syndrome_rank: int,
    coherent_dimension: int,
    block_count: int,
    alphabet_size: int = 2,
) -> float:
    """Bound RETURN for approximate AUDIT across full-rank temporal blocks.

    The ordered check must split into ``block_count`` consecutive blocks of
    full ``syndrome_rank`` over the finite field of size ``alphabet_size``.
    With coherent dimension ``d`` at every block boundary, put

    ``alpha = min(1, (d / alphabet_size**syndrome_rank)**block_count)``

    and ``theta = min(block_count * (1 - P_A), 1 - alpha)``.  Causal list
    decoding and a Ky Fan rank-tail bound give

    ``F_R <= alpha + (1 - 2*alpha)*theta``
    ``       + 2*sqrt(alpha*(1-alpha)*theta*(1-theta))``.

    At ``P_A=1`` this reduces to the exact temporal power-law endpoint.  The
    result allows arbitrary adaptive non-QND instruments and unrestricted
    finite genuinely classical transcript under the declared sequestration
    interface.  It is an upper certificate, not generally an attainable
    frontier.
    """
    probability = float(audit_probability)
    if not 0.0 <= probability <= 1.0:
        raise ValueError("audit_probability must lie in [0, 1]")
    rank = _positive_integer(syndrome_rank, "syndrome_rank")
    dimension = _positive_integer(coherent_dimension, "coherent_dimension")
    blocks = _positive_integer(block_count, "block_count")
    alphabet = _positive_integer(alphabet_size, "alphabet_size")
    if alphabet < 2:
        raise ValueError("alphabet_size must be at least two")
    if not _is_prime_power(alphabet):
        raise ValueError("alphabet_size must be a prime power")
    syndrome_count = alphabet**rank
    alpha = min(1.0, (dimension / syndrome_count) ** blocks)
    if alpha == 1.0:
        return 1.0
    theta = min(blocks * (1.0 - probability), 1.0 - alpha)
    return (
        alpha
        + (1.0 - 2.0 * alpha) * theta
        + 2.0
        * sqrt(max(0.0, alpha * (1.0 - alpha) * theta * (1.0 - theta)))
    )


def ordered_check_perfect_audit_return_bound(
    matrix: ArrayLike,
    coherent_dimension: int,
) -> float:
    """Bound perfect-AUDIT return for an ordered binary check matrix.

    If ``r = rank(H)`` and ``mu`` is the maximum number of consecutive
    full-rank blocks returned by :func:`full_rank_block_packing_number`, the
    temporal product theorem gives

    ``F_R <= min(1, (d / 2**r)**mu)``.

    The zero-rank case returns one.  This endpoint bound does not assert
    attainability for every matrix or solve an approximate-AUDIT frontier.
    """
    checks = _binary_integer_matrix(matrix)
    rank = gf2_rank(checks)
    if rank == 0:
        _positive_integer(coherent_dimension, "coherent_dimension")
        return 1.0
    return full_rank_block_perfect_audit_return_bound(
        syndrome_rank=rank,
        coherent_dimension=coherent_dimension,
        block_count=full_rank_block_packing_number(checks),
        alphabet_size=2,
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


def interleaved_return_upper_bound(audit_probability: float) -> float:
    """Return the rigorous linear-tail bound for the interleaved benchmark.

    This is the four-slot binary specialization of
    :func:`full_rank_block_approximate_audit_return_bound`.  It is exact at
    perfect AUDIT, where it returns ``1/4``, and becomes trivial at one once
    ``P_A <= 5/8``.  The square-root approach to the endpoint has the correct
    exponent, but the complete curve is not claimed optimal.
    """
    return full_rank_block_approximate_audit_return_bound(
        audit_probability=audit_probability,
        syndrome_rank=2,
        coherent_dimension=2,
        block_count=2,
        alphabet_size=2,
    )


def interleaved_support_upper_bound(audit_weight: float = 0.5) -> float:
    """Return the explicit causal-list upper certificate on interleaved support.

    For ``lambda`` equal to ``audit_weight``, the certificate is

    ``1/2 + lambda/4 + sqrt(7*lambda**2 - 10*lambda + 4)/4``.

    It lies strictly below the grouped/static support value for
    ``3/7 < lambda < 1`` and equals ``5/8 + sqrt(3)/8`` at balanced weight.
    It remains an upper bound rather than the exact interleaved frontier.
    """
    weight = _audit_weight(audit_weight)
    return (
        0.5
        + weight / 4.0
        + sqrt(max(0.0, 7.0 * weight**2 - 10.0 * weight + 4.0)) / 4.0
    )


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
    is an exact lower bound.  A finite-outcome non-QND counterexample exceeds
    this family, while the arbitrary-instrument interior converse remains
    open.
    """

    audit_weight: float
    strategy: str
    q: float | None
    v: float | None
    audit_probability: float
    return_fidelity: float
    support_value: float
    support_is_globally_optimal: bool = False


@dataclass(frozen=True)
class StoredCounterexamplePoint:
    """Framework-neutral verified point outside the two-parameter family."""

    audit_weight: float
    local_outcome_arity: int
    audit_probability: float
    return_fidelity: float
    support_value: float
    independently_verified: bool
    support_is_globally_optimal: bool = False


@dataclass(frozen=True)
class InterleavedCompactCandidatePoint:
    """One exposed point of the compact three-effect MPS construction.

    The construction is a complete physical streamed strategy after local
    Pauli completion.  ``support_is_globally_optimal`` remains false because
    agreement with unrestricted MPS and cq-instrument searches is numerical,
    not a certified global converse.
    """

    audit_weight: float
    strategy: str
    t: float | None
    r: float | None
    priors: tuple[float, float, float, float] | None
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


INTERLEAVED_BALANCED_COUNTEREXAMPLE: Final[StoredCounterexamplePoint] = (
    StoredCounterexamplePoint(
        audit_weight=0.5,
        local_outcome_arity=3,
        audit_probability=0.6257545618203884,
        return_fidelity=0.8931433788141326,
        support_value=0.7594489703172604,
        independently_verified=True,
        support_is_globally_optimal=False,
    )
)
"""Stored complete instrument that falsifies the older restricted candidate.

The local Kraus tree and AUDIT effects are in
``data/interleaved_ternary_counterexample.npz``.  The value is an achievable
lower bound, not a claimed optimum of the unrestricted interior game.  It is
superseded by :func:`interleaved_compact_lower_bound`.
"""


def interleaved_candidate_scores(q: float, v: float) -> tuple[float, float]:
    """Evaluate the exact two-parameter interleaved construction.

    The returned pair is ``(P_A, F_R)`` with

    ``P_A = 1/2 + q*v*sqrt(1-v**2) - q*(1-q)*v**2``

    and the corresponding flagged polar-recovery fidelity.  The parameters
    must lie in the unit square.  This is a physical streamed one-qubit
    construction, not an asserted characterization of every possible
    interleaved strategy.  In particular, the stored ternary-outcome
    counterexample has a strictly larger balanced score.
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


def _compact_candidate_evaluate(
    t: float,
    r: float,
    audit_weight: float,
) -> tuple[float, float, float, tuple[float, float, float, float]]:
    """Evaluate the three-effect family after analytic prior optimisation."""

    t_value = float(t)
    r_value = float(r)
    if not 0.0 <= t_value <= 1.0:
        raise ValueError("t must lie in [0, 1]")
    if not -1.0 <= r_value <= 1.0:
        raise ValueError("r must lie in [-1, 1]")
    weight = _audit_weight(audit_weight)

    other_trace = 1.0 - t_value / 2.0
    effect_z = -t_value / (2.0 - t_value)
    effect_x = sqrt(max(0.0, 1.0 - effect_z**2))
    state_x = sqrt(max(0.0, 1.0 - r_value**2))
    other_correct = (
        other_trace
        * (1.0 + effect_z * r_value + effect_x * state_x)
        / 2.0
    )
    coarse_other = t_value * (1.0 + r_value) / 2.0
    c_reference = sqrt(t_value) + sqrt(max(0.0, 1.0 - t_value))
    c_null = sqrt(2.0) if t_value >= 0.5 else c_reference
    c_other = sqrt(max(0.0, coarse_other)) + sqrt(
        max(0.0, 1.0 - coarse_other)
    )

    audit_diagonal = np.asarray((t_value, 0.0, other_correct), dtype=float)
    hellinger = np.asarray(
        (c_reference, c_null, sqrt(2.0) * c_other), dtype=float
    )
    matrix = weight * np.diag(audit_diagonal)
    matrix += (1.0 - weight) * np.outer(hellinger, hellinger) / 8.0
    _, eigenvectors = np.linalg.eigh(matrix)
    eigenvector = np.abs(eigenvectors[:, -1])
    eigenvector /= np.linalg.norm(eigenvector)
    priors = (
        float(eigenvector[0] ** 2),
        float(eigenvector[1] ** 2),
        float(eigenvector[2] ** 2 / 2.0),
        float(eigenvector[2] ** 2 / 2.0),
    )
    audit_probability = (
        priors[0] * t_value + (priors[2] + priors[3]) * other_correct
    )
    return_fidelity = (
        c_reference * sqrt(priors[0])
        + c_null * sqrt(priors[1])
        + c_other * (sqrt(priors[2]) + sqrt(priors[3]))
    ) ** 2 / 8.0
    support_value = (
        weight * audit_probability + (1.0 - weight) * return_fidelity
    )
    return support_value, audit_probability, return_fidelity, priors


def interleaved_compact_lower_bound(
    audit_weight: float = 0.5,
) -> InterleavedCompactCandidatePoint:
    """Optimise the compact three-effect interleaved MPS construction.

    The local geometry has two real parameters.  Completeness fixes a
    three-effect extremal qubit POVM from ``t``; ``r`` is the symmetric signal
    state's Bloch coordinate.  The four prior weights are optimised exactly
    by a three-by-three Perron eigenproblem.  A deterministic multistart search
    over ``sqrt(t)`` and ``r`` is compared with the exact no-record point.

    Every returned nontrivial point is physically attainable by a bond-two
    Choi MPS and a four-outcome-per-slot Pauli completion.  Numerical
    optimisation selects the best point inside this explicit family; it does
    not prove the unrestricted MPS maximum.
    """

    weight = _audit_weight(audit_weight)
    no_record_score = 1.0 - weight / 2.0

    def negative_score(point: np.ndarray) -> float:
        u_value = min(1.0, max(0.0, float(point[0])))
        r_value = min(1.0, max(-1.0, float(point[1])))
        score, _, _, _ = _compact_candidate_evaluate(
            u_value**2, r_value, weight
        )
        return -score

    starts = (
        (sqrt(0.56), 0.03),
        (sqrt(0.46), 0.016),
        (sqrt(0.15), 0.03),
        (0.20, 0.0),
        (0.05, 0.0),
        (0.005, 0.0),
        (0.999, -0.60),
    )
    results = (
        minimize(
            negative_score,
            np.asarray(start, dtype=float),
            method="Nelder-Mead",
            bounds=((0.0, 1.0), (-1.0, 1.0)),
            options={"xatol": 1e-12, "fatol": 1e-14, "maxiter": 10000},
        )
        for start in starts
    )
    best_result = min(results, key=lambda result: float(result.fun))
    u_value, r_value = map(float, best_result.x)
    t_value = u_value**2
    score, audit_probability, return_fidelity, priors = (
        _compact_candidate_evaluate(t_value, r_value, weight)
    )
    if no_record_score >= score - 5e-12:
        return InterleavedCompactCandidatePoint(
            audit_weight=weight,
            strategy="no_record",
            t=None,
            r=None,
            priors=None,
            audit_probability=0.5,
            return_fidelity=1.0,
            support_value=no_record_score,
        )
    return InterleavedCompactCandidatePoint(
        audit_weight=weight,
        strategy="three_effect_mps",
        t=t_value,
        r=r_value,
        priors=priors,
        audit_probability=audit_probability,
        return_fidelity=return_fidelity,
        support_value=score,
    )


__all__ = [
    "GROUPED_CHECK_MATRIX",
    "INTERLEAVED_CHECK_MATRIX",
    "INTERLEAVED_ORDER_GAP_WEIGHT_THRESHOLD",
    "INTERLEAVED_BALANCED_COUNTEREXAMPLE",
    "INTERLEAVED_PERFECT_AUDIT_ENDPOINT",
    "GroupedFrontierPoint",
    "InterleavedCandidatePoint",
    "InterleavedCompactCandidatePoint",
    "PerfectAuditEndpoint",
    "StoredCounterexamplePoint",
    "full_crossing_cuts",
    "full_crossing_perfect_audit_return_bound",
    "full_rank_block_packing_number",
    "full_rank_block_approximate_audit_return_bound",
    "full_rank_block_perfect_audit_return_bound",
    "gf2_rank",
    "grouped_frontier",
    "interleaved_candidate_lower_bound",
    "interleaved_candidate_scores",
    "interleaved_compact_lower_bound",
    "interleaved_return_upper_bound",
    "interleaved_support_upper_bound",
    "ordered_check_perfect_audit_return_bound",
    "rank_two_static_qubit_support",
    "trellis_connectivity_profile",
    "trellis_connectivity_tau",
]

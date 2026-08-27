"""Analytic state elimination for binary-projective secular bounds.

For a coarse qubit effect ``E`` and a rank-one congruence split ``G``, a
projective secular summand has the form

``(1 + 2*sqrt(q*(1-q))) / (level - weight*d)``

with ``q = <psi|E|psi>`` and ``d = <psi|G|psi>``.  The AM--GM tangent

``1 + 2*sqrt(q*(1-q)) <= 1 + 1/t + (t - 1/t)*q``

turns its maximisation over ``psi`` into a two-dimensional generalized
eigenvalue.  The closed form below removes one continuous state coordinate
from every term.  Any positive tangent parameter gives a valid upper bound;
choosing ``t = sqrt((1-q*)/q*)`` makes the numerator tangent at ``q*``.
"""

from __future__ import annotations

import math


def hellinger_tangent(
    probability: float, tangent: float
) -> tuple[float, float, float]:
    """Return the exact Hellinger numerator, its tangent, and their gap."""

    if not 0.0 <= probability <= 1.0:
        raise ValueError("probability must lie in [0, 1]")
    if not math.isfinite(tangent) or tangent <= 0.0:
        raise ValueError("tangent must be finite and strictly positive")
    exact = 1.0 + 2.0 * math.sqrt(probability * (1.0 - probability))
    affine = 1.0 + 1.0 / tangent + (
        tangent - 1.0 / tangent
    ) * probability
    return exact, affine, affine - exact


def endpoint_split_state_term(
    high: float,
    low: float,
    state: float,
    label: int,
    weight: float,
    level: float,
) -> float:
    """Evaluate one endpoint-split secular term at a pure-state coordinate."""

    _validate_parameters(high, low, 0.0, label, weight, level)
    if not 0.0 <= state <= 1.0:
        raise ValueError("state must lie in [0, 1]")
    probability = low + (high - low) * state
    decision = probability if label == 0 else 0.0
    numerator = 1.0 + 2.0 * math.sqrt(probability * (1.0 - probability))
    return numerator / (level - weight * decision)


def endpoint_split_tangent_upper(
    high: float,
    low: float,
    label: int,
    tangent: float,
    weight: float,
    level: float,
) -> float:
    """Upper-bound an endpoint-split term uniformly over pure states."""

    _validate_parameters(high, low, 0.0, label, weight, level)
    if not math.isfinite(tangent) or tangent <= 0.0:
        raise ValueError("tangent must be finite and strictly positive")
    inverse = 1.0 / tangent
    intercept = 1.0 + inverse
    slope = tangent - inverse
    numerator_high = intercept + slope * high
    numerator_low = intercept + slope * low
    if label == 0:
        high_value = numerator_high / (level - weight * high)
        low_value = numerator_low / (level - weight * low)
    else:
        high_value = numerator_high / level
        low_value = numerator_low / level
    return max(high_value, low_value)


def rank_split_state_term(
    high: float,
    low: float,
    sine: float,
    state: float,
    label: int,
    weight: float,
    level: float,
) -> float:
    """Evaluate one rank-split secular term at a pure-state coordinate."""

    _validate_parameters(high, low, sine, label, weight, level)
    if not 0.0 <= state <= 1.0:
        raise ValueError("state must lie in [0, 1]")
    cosine = math.sqrt(1.0 - sine * sine)
    root_state = math.sqrt(state)
    root_complement = math.sqrt(1.0 - state)
    probability = low + (high - low) * state
    if label == 0:
        amplitude = (
            math.sqrt(high) * cosine * root_state
            + math.sqrt(low) * sine * root_complement
        )
    else:
        amplitude = (
            math.sqrt(high) * sine * root_state
            + math.sqrt(low) * cosine * root_complement
        )
    decision = amplitude * amplitude
    numerator = 1.0 + 2.0 * math.sqrt(probability * (1.0 - probability))
    return numerator / (level - weight * decision)


def rank_split_tangent_upper(
    high: float,
    low: float,
    sine: float,
    label: int,
    tangent: float,
    weight: float,
    level: float,
) -> float:
    """Upper-bound a rank-split term uniformly over all pure states.

    The result is the larger generalized eigenvalue of

    ``(A*I + B*E) v = T * (level*I - weight*G) v``,

    where ``A = 1 + 1/t`` and ``B = t - 1/t``.  The discriminant is evaluated
    as a sum of nonnegative terms, avoiding cancellation near repeated roots.
    """

    _validate_parameters(high, low, sine, label, weight, level)
    if not math.isfinite(tangent) or tangent <= 0.0:
        raise ValueError("tangent must be finite and strictly positive")

    sine_squared = sine * sine
    cosine_squared = 1.0 - sine_squared
    if label == 0:
        first_split, second_split = cosine_squared, sine_squared
    else:
        first_split, second_split = sine_squared, cosine_squared

    inverse = 1.0 / tangent
    intercept = 1.0 + inverse
    slope = tangent - inverse
    numerator_first = intercept + slope * high
    numerator_second = intercept + slope * low

    denominator_first = level - weight * high * first_split
    denominator_second = level - weight * low * second_split
    split_trace = high * first_split + low * second_split
    denominator_determinant = level * level - level * weight * split_trace
    if denominator_determinant <= 0.0:
        raise ValueError("level*I - weight*G is not positive definite")

    cross = (
        numerator_first * denominator_second
        + numerator_second * denominator_first
    )
    diagonal_gap = (
        numerator_first * denominator_second
        - numerator_second * denominator_first
    )
    off_diagonal_squared = (
        weight
        * weight
        * high
        * low
        * sine_squared
        * cosine_squared
    )
    discriminant = diagonal_gap * diagonal_gap + (
        4.0
        * numerator_first
        * numerator_second
        * off_diagonal_squared
    )
    return (
        cross + math.sqrt(max(0.0, discriminant))
    ) / (2.0 * denominator_determinant)


def tangent_for_probability(probability: float) -> float:
    """Return the positive AM--GM tangent that touches at ``probability``."""

    if not 0.0 < probability < 1.0:
        raise ValueError("probability must lie strictly between zero and one")
    return math.sqrt((1.0 - probability) / probability)


def _validate_parameters(
    high: float,
    low: float,
    sine: float,
    label: int,
    weight: float,
    level: float,
) -> None:
    if not 0.0 <= low <= high <= 1.0:
        raise ValueError("effect eigenvalues must satisfy 0 <= low <= high <= 1")
    if not 0.0 <= sine <= 1.0 / math.sqrt(2.0):
        raise ValueError("canonical split sine must lie in [0, 1/sqrt(2)]")
    if label not in (0, 1):
        raise ValueError("label must be zero or one")
    if not 0.0 < weight < 1.0:
        raise ValueError("weight must lie strictly between zero and one")
    if not math.isfinite(level) or level <= weight * high:
        raise ValueError("level must exceed weight times the high eigenvalue")


__all__ = [
    "endpoint_split_state_term",
    "endpoint_split_tangent_upper",
    "hellinger_tangent",
    "rank_split_state_term",
    "rank_split_tangent_upper",
    "tangent_for_probability",
]

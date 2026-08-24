"""Outward interval enclosure of ternary planar reconstruction matrices.

For the canonical rank-one ternary qubit POVM, Horwitz parameters ``A,B``
determine the effect weights and the angle between the first two effects.  If
``q`` is the three-outcome probability vector of a Hermitian operator, its
trace is ``sum(q)`` and its two visible Bloch coordinates are ``R(A,B) q``.

This module encloses every coefficient of ``R`` over an axis-aligned ``A,B``
box using elementary outward-rounded interval arithmetic.  The returned
column errors certify

    ||(R(A,B)-R(anchor)) q||_2 <= sum_t error[t] |q[t]|.

Boxes touching a degenerate projective edge have zero sine and are rejected;
they must be handled by the independent projective sector or subdivided away.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np


def down(value: float) -> float:
    return float(np.nextafter(value, -math.inf))


def up(value: float) -> float:
    return float(np.nextafter(value, math.inf))


@dataclass(frozen=True)
class Interval:
    lower: float
    upper: float

    def __post_init__(self) -> None:
        if not math.isfinite(self.lower) or not math.isfinite(self.upper):
            raise ValueError("interval endpoints must be finite")
        if self.lower > self.upper:
            raise ValueError("interval endpoints are reversed")

    @staticmethod
    def point(value: float) -> "Interval":
        return Interval(float(value), float(value))

    def __add__(self, other: "Interval") -> "Interval":
        return Interval(down(self.lower + other.lower), up(self.upper + other.upper))

    def __neg__(self) -> "Interval":
        return Interval(down(-self.upper), up(-self.lower))

    def __sub__(self, other: "Interval") -> "Interval":
        return self + (-other)

    def __mul__(self, other: "Interval") -> "Interval":
        values = (
            self.lower * other.lower,
            self.lower * other.upper,
            self.upper * other.lower,
            self.upper * other.upper,
        )
        return Interval(down(min(values)), up(max(values)))

    def reciprocal(self) -> "Interval":
        if self.lower <= 0.0 <= self.upper:
            raise ValueError("cannot invert an interval containing zero")
        values = (1.0 / self.lower, 1.0 / self.upper)
        return Interval(down(min(values)), up(max(values)))

    def __truediv__(self, other: "Interval") -> "Interval":
        return self * other.reciprocal()

    def square(self) -> "Interval":
        upper = max(self.lower * self.lower, self.upper * self.upper)
        lower = 0.0 if self.lower <= 0.0 <= self.upper else min(
            self.lower * self.lower, self.upper * self.upper
        )
        return Interval(max(0.0, down(lower)), up(upper))

    def sqrt(self) -> "Interval":
        if self.lower < 0.0:
            raise ValueError("cannot take the square root of a negative interval")
        return Interval(
            max(0.0, down(math.sqrt(self.lower))),
            up(math.sqrt(self.upper)),
        )


ONE = Interval.point(1.0)
TWO = Interval.point(2.0)


def planar_reconstruction(alpha: float, beta: float) -> np.ndarray:
    """Return the exact 2-by-3 visible-Bloch reconstruction matrix."""

    denominator = alpha + beta - 1.0
    w0 = alpha / denominator
    w1 = beta / denominator
    cosine = 1.0 - 2.0 / alpha - 2.0 / beta + 2.0 / (alpha * beta)
    sine = math.sqrt(max(0.0, 1.0 - cosine * cosine))
    if sine <= 0.0:
        raise ValueError("the terminal POVM is projectively degenerate")
    x = np.asarray([2.0 / w0 - 1.0, -1.0, -1.0])
    y = np.asarray(
        [
            -1.0 - cosine * x[0],
            2.0 / w1 - 1.0 + cosine,
            -1.0 + cosine,
        ]
    ) / sine
    return np.vstack([x, y])


def reconstruction_intervals(
    alpha_bounds: tuple[float, float],
    beta_bounds: tuple[float, float],
) -> tuple[tuple[Interval, ...], tuple[Interval, ...]]:
    """Enclose every coefficient of the two reconstruction rows."""

    alpha = Interval(*map(float, alpha_bounds))
    beta = Interval(*map(float, beta_bounds))
    denominator = alpha + beta - ONE
    w0 = alpha / denominator
    w1 = beta / denominator
    cosine = ONE - TWO / alpha - TWO / beta + TWO / (alpha * beta)
    sine_squared = ONE - cosine.square()
    if sine_squared.lower <= 0.0:
        raise ValueError("terminal box reaches a degenerate projective edge")
    sine = sine_squared.sqrt()
    x = (TWO / w0 - ONE, -ONE, -ONE)
    y = (
        (-ONE - cosine * x[0]) / sine,
        (TWO / w1 - ONE + cosine) / sine,
        (-ONE + cosine) / sine,
    )
    return x, y


def reconstruction_anchor_and_errors(
    alpha_bounds: tuple[float, float],
    beta_bounds: tuple[float, float],
) -> tuple[np.ndarray, np.ndarray, dict[str, object]]:
    """Return a midpoint anchor and safe per-column Euclidean errors."""

    alpha = 0.5 * (alpha_bounds[0] + alpha_bounds[1])
    beta = 0.5 * (beta_bounds[0] + beta_bounds[1])
    anchor = planar_reconstruction(alpha, beta)
    rows = reconstruction_intervals(alpha_bounds, beta_bounds)
    errors = np.empty(3, dtype=float)
    interval_payload: list[list[list[float]]] = [[], []]
    for row_index, row in enumerate(rows):
        interval_payload[row_index] = [
            [item.lower, item.upper] for item in row
        ]
    for column in range(3):
        deviations = []
        for row in range(2):
            interval = rows[row][column]
            deviations.append(
                max(
                    abs(interval.lower - anchor[row, column]),
                    abs(interval.upper - anchor[row, column]),
                )
            )
        errors[column] = up(math.hypot(*deviations))
    return anchor, errors, {
        "anchor_parameters": [alpha, beta],
        "anchor_matrix": anchor.tolist(),
        "coefficient_intervals": interval_payload,
        "column_errors": errors.tolist(),
        "rounding": "IEEE-754 nextafter outward after every interval operation",
    }

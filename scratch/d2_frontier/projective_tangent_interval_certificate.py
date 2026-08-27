"""Solver-independent interval certificate for rank/rank projective boxes.

Each secular state coordinate is eliminated by the AM--GM tangent

    1 + 2 sqrt(q (1-q)) <= 1 + 1/t + (t - 1/t) q.

The resulting state maximum is a closed 2-by-2 generalized eigenvalue.  The
code evaluates the sum with outward-expanded binary64 intervals and uses
mean-value and second-order Taylor forms to retain cancellations between the
four terms.  SciPy only
chooses positive tangent parameters and branching coordinates; neither choice
is trusted by the proof kernel.

The certificate is independent of SCIP.  Its current trust boundary is the
IEEE-754 round-to-nearest/``nextafter`` contract of Python and NumPy.  Exact
decimal inputs are first enclosed from ``Fraction`` objects.  A future dyadic
integer backend can replay the same algebra without changing the reduction.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
import heapq
import json
import math
from pathlib import Path
from typing import Iterable

from scipy.optimize import minimize_scalar

from carmenq.projective_secular import (
    endpoint_split_state_term,
    rank_split_state_term,
)


def down(value: float) -> float:
    return math.nextafter(float(value), -math.inf)


def up(value: float) -> float:
    return math.nextafter(float(value), math.inf)


def rational_down(value: Fraction) -> float:
    candidate = float(value)
    if Fraction.from_float(candidate) > value:
        candidate = down(candidate)
    if Fraction.from_float(candidate) > value:
        raise ArithmeticError("failed to round a rational downward")
    return candidate


def rational_up(value: Fraction) -> float:
    candidate = float(value)
    if Fraction.from_float(candidate) < value:
        candidate = up(candidate)
    if Fraction.from_float(candidate) < value:
        raise ArithmeticError("failed to round a rational upward")
    return candidate


@dataclass(frozen=True)
class Interval:
    lo: float
    hi: float

    def __post_init__(self) -> None:
        if not (math.isfinite(self.lo) and math.isfinite(self.hi)):
            raise ValueError("interval endpoints must be finite")
        if self.lo > self.hi:
            raise ValueError("empty interval")

    @classmethod
    def exact_float(cls, value: float) -> "Interval":
        """Treat one finite binary64 number as an exact dyadic rational."""

        if not math.isfinite(value):
            raise ValueError("exact float must be finite")
        return cls(float(value), float(value))

    @classmethod
    def rational(cls, value: Fraction) -> "Interval":
        return cls(rational_down(value), rational_up(value))

    @classmethod
    def decimal(cls, value: str | float) -> "Interval":
        return cls.rational(Fraction(str(value)))

    @property
    def midpoint(self) -> float:
        return self.lo + (self.hi - self.lo) / 2.0

    @property
    def width(self) -> float:
        return self.hi - self.lo


ZERO = Interval.exact_float(0.0)
ONE = Interval.exact_float(1.0)


def add(left: Interval, right: Interval) -> Interval:
    return Interval(down(left.lo + right.lo), up(left.hi + right.hi))


def neg(value: Interval) -> Interval:
    return Interval(down(-value.hi), up(-value.lo))


def sub(left: Interval, right: Interval) -> Interval:
    return add(left, neg(right))


def mul(left: Interval, right: Interval) -> Interval:
    candidates = (
        left.lo * right.lo,
        left.lo * right.hi,
        left.hi * right.lo,
        left.hi * right.hi,
    )
    return Interval(down(min(candidates)), up(max(candidates)))


def square(value: Interval) -> Interval:
    if value.lo <= 0.0 <= value.hi:
        lower = 0.0
    else:
        lower = min(value.lo * value.lo, value.hi * value.hi)
    upper = max(value.lo * value.lo, value.hi * value.hi)
    return Interval(0.0 if lower == 0.0 else down(lower), up(upper))


def divide(left: Interval, right: Interval) -> Interval:
    if right.lo <= 0.0 <= right.hi:
        raise ZeroDivisionError("interval denominator crosses zero")
    reciprocal = Interval(down(1.0 / right.hi), up(1.0 / right.lo))
    return mul(left, reciprocal)


def root(value: Interval) -> Interval:
    if value.hi < 0.0:
        raise ValueError("square root of a negative interval")
    lower = math.sqrt(max(0.0, value.lo))
    return Interval(
        0.0 if lower == 0.0 else down(lower),
        up(math.sqrt(value.hi)),
    )


def scale(value: Interval, scalar: float) -> Interval:
    return mul(value, Interval.exact_float(scalar))


@dataclass(frozen=True)
class Jet:
    value: Interval
    gradient: tuple[Interval, ...]


DIMENSION = 4


def constant(value: Interval | float) -> Jet:
    interval = value if isinstance(value, Interval) else Interval.exact_float(value)
    return Jet(interval, (ZERO,) * DIMENSION)


def variable(value: Interval, index: int) -> Jet:
    gradient = [ZERO] * DIMENSION
    gradient[index] = ONE
    return Jet(value, tuple(gradient))


def jadd(left: Jet, right: Jet) -> Jet:
    return Jet(
        add(left.value, right.value),
        tuple(add(a, b) for a, b in zip(left.gradient, right.gradient)),
    )


def jneg(value: Jet) -> Jet:
    return Jet(neg(value.value), tuple(neg(item) for item in value.gradient))


def jsub(left: Jet, right: Jet) -> Jet:
    return jadd(left, jneg(right))


def jmul(left: Jet, right: Jet) -> Jet:
    return Jet(
        mul(left.value, right.value),
        tuple(
            add(mul(a, right.value), mul(left.value, b))
            for a, b in zip(left.gradient, right.gradient)
        ),
    )


def jsquare(value: Jet) -> Jet:
    return Jet(
        square(value.value),
        tuple(scale(mul(value.value, item), 2.0) for item in value.gradient),
    )


def jdivide(left: Jet, right: Jet) -> Jet:
    denominator = square(right.value)
    return Jet(
        divide(left.value, right.value),
        tuple(
            divide(
                sub(mul(a, right.value), mul(left.value, b)),
                denominator,
            )
            for a, b in zip(left.gradient, right.gradient)
        ),
    )


def jroot(value: Jet) -> Jet:
    rooted = root(value.value)
    denominator = scale(rooted, 2.0)
    if denominator.lo <= 0.0:
        raise ZeroDivisionError("square-root derivative reaches zero")
    return Jet(
        rooted,
        tuple(divide(item, denominator) for item in value.gradient),
    )

def jmaximum(left: Jet, right: Jet) -> Jet:
    """Enclose the pointwise maximum of two scalar jets."""
    if left.value.lo >= right.value.hi:
        return left
    if right.value.lo >= left.value.hi:
        return right
    value = Interval(
        max(left.value.lo, right.value.lo),
        max(left.value.hi, right.value.hi),
    )
    gradient = tuple(
        Interval(min(a.lo, b.lo), max(a.hi, b.hi))
        for a, b in zip(left.gradient, right.gradient)
    )
    return Jet(value, gradient)


@dataclass(frozen=True)
class Jet2:
    """Second-order interval jet in the four exterior coordinates."""

    value: Interval
    gradient: tuple[Interval, ...]
    hessian: tuple[tuple[Interval, ...], ...]


ZERO_HESSIAN = tuple(
    tuple(ZERO for _ in range(DIMENSION)) for _ in range(DIMENSION)
)


def constant2(value: Interval | float) -> Jet2:
    interval = value if isinstance(value, Interval) else Interval.exact_float(value)
    return Jet2(interval, (ZERO,) * DIMENSION, ZERO_HESSIAN)


def variable2(value: Interval, index: int) -> Jet2:
    gradient = [ZERO] * DIMENSION
    gradient[index] = ONE
    return Jet2(value, tuple(gradient), ZERO_HESSIAN)


def j2add(left: Jet2, right: Jet2) -> Jet2:
    return Jet2(
        add(left.value, right.value),
        tuple(add(a, b) for a, b in zip(left.gradient, right.gradient)),
        tuple(
            tuple(
                add(left.hessian[i][j], right.hessian[i][j])
                for j in range(DIMENSION)
            )
            for i in range(DIMENSION)
        ),
    )


def j2neg(value: Jet2) -> Jet2:
    return Jet2(
        neg(value.value),
        tuple(neg(item) for item in value.gradient),
        tuple(tuple(neg(item) for item in row) for row in value.hessian),
    )


def j2sub(left: Jet2, right: Jet2) -> Jet2:
    return j2add(left, j2neg(right))


def j2mul(left: Jet2, right: Jet2) -> Jet2:
    gradient = tuple(
        add(
            mul(left.gradient[i], right.value),
            mul(left.value, right.gradient[i]),
        )
        for i in range(DIMENSION)
    )
    hessian = tuple(
        tuple(
            add(
                add(
                    mul(left.hessian[i][j], right.value),
                    mul(left.gradient[i], right.gradient[j]),
                ),
                add(
                    mul(left.gradient[j], right.gradient[i]),
                    mul(left.value, right.hessian[i][j]),
                ),
            )
            for j in range(DIMENSION)
        )
        for i in range(DIMENSION)
    )
    return Jet2(mul(left.value, right.value), gradient, hessian)


def j2inverse(value: Jet2) -> Jet2:
    inverse_value = divide(ONE, value.value)
    inverse_square = square(inverse_value)
    inverse_cube = mul(inverse_square, inverse_value)
    gradient = tuple(neg(mul(item, inverse_square)) for item in value.gradient)
    hessian = tuple(
        tuple(
            sub(
                scale(
                    mul(mul(value.gradient[i], value.gradient[j]), inverse_cube),
                    2.0,
                ),
                mul(value.hessian[i][j], inverse_square),
            )
            for j in range(DIMENSION)
        )
        for i in range(DIMENSION)
    )
    return Jet2(inverse_value, gradient, hessian)


def j2divide(left: Jet2, right: Jet2) -> Jet2:
    return j2mul(left, j2inverse(right))


def j2square(value: Jet2) -> Jet2:
    return j2mul(value, value)


def j2root(value: Jet2) -> Jet2:
    rooted = root(value.value)
    if rooted.lo <= 0.0:
        raise ZeroDivisionError("square-root Hessian reaches zero")
    first = divide(ONE, scale(rooted, 2.0))
    second = neg(divide(ONE, scale(mul(value.value, rooted), 4.0)))
    gradient = tuple(mul(first, item) for item in value.gradient)
    hessian = tuple(
        tuple(
            add(
                mul(first, value.hessian[i][j]),
                mul(second, mul(value.gradient[i], value.gradient[j])),
            )
            for j in range(DIMENSION)
        )
        for i in range(DIMENSION)
    )
    return Jet2(rooted, gradient, hessian)


@dataclass(frozen=True)
class RankRankBox:
    y: Interval
    residual: Interval
    first_sine: Interval
    second_sine: Interval

    @property
    def coordinates(self) -> tuple[Interval, ...]:
        return (self.y, self.residual, self.first_sine, self.second_sine)

    def serialise(self) -> dict[str, list[float]]:
        return {
            name: [getattr(self, name).lo, getattr(self, name).hi]
            for name in ("y", "residual", "first_sine", "second_sine")
        }

    @classmethod
    def deserialise(cls, payload: dict[str, list[float]]) -> "RankRankBox":
        names = ("y", "residual", "first_sine", "second_sine")
        if set(payload) != set(names):
            raise ValueError("serialized box has the wrong coordinates")
        coordinates = []
        for name in names:
            values = payload[name]
            if len(values) != 2:
                raise ValueError(f"serialized {name} is not an interval")
            coordinates.append(Interval(float(values[0]), float(values[1])))
        return cls(*coordinates)


@dataclass(frozen=True)
class BoundResult:
    upper: float
    method: str
    tangents: tuple[float, ...]
    gradient: tuple[Interval, ...] | None


class RankRankCertifier:
    """Directed-interval engine for a projective topology on a symmetric chart.

    The historical class name is retained for compatibility. ``first_kind``
    and ``second_kind`` may independently be ``"endpoint"`` or ``"rank"``.
    """

    def __init__(
        self,
        weight: Fraction,
        level: Fraction,
        x_bounds: Interval,
        first_kind: str = "rank",
        second_kind: str = "rank",
    ):
        if not Fraction(0) < weight < Fraction(1):
            raise ValueError("weight must lie strictly between zero and one")
        if first_kind not in {"endpoint", "rank"}:
            raise ValueError("unknown first projective topology")
        if second_kind not in {"endpoint", "rank"}:
            raise ValueError("unknown second projective topology")
        self.weight_fraction = weight
        self.level_fraction = level
        self.weight = Interval.rational(weight)
        self.level = Interval.rational(level)
        self.scale = divide(
            sub(ONE, self.weight), Interval.rational(Fraction(8))
        )
        self.weight_float = float(weight)
        self.level_float = float(level)
        self.x_bounds = x_bounds
        self.first_kind = first_kind
        self.second_kind = second_kind

    def contract(self, box: RankRankBox) -> RankRankBox | None:
        ylo, yhi = box.y.lo, box.y.hi
        rlo, rhi = box.residual.lo, box.residual.hi
        sum_bounds = sub(ONE, self.x_bounds)
        for _ in range(8):
            ylo = max(ylo, down(sum_bounds.lo - rhi))
            yhi = min(yhi, up(sum_bounds.hi - rlo))
            rlo = max(rlo, down(sum_bounds.lo - yhi), 0.0)
            # The canonical eigenvalue order x >= y is 2*y + residual <= 1.
            yhi = min(yhi, up((1.0 - rlo) / 2.0), 0.5)
            rhi = min(rhi, up(sum_bounds.hi - ylo), up(1.0 - 2.0 * ylo))
            if (
                ylo > yhi or rlo > rhi or ylo + rlo > sum_bounds.hi
                or 2.0 * ylo + rlo > 1.0
            ):
                return None
        return RankRankBox(
            Interval(ylo, yhi),
            Interval(rlo, rhi),
            box.first_sine,
            box.second_sine,
        )

    def _select_tangents(self, box: RankRankBox) -> tuple[float, ...]:
        y = box.y.midpoint
        residual = box.residual.midpoint
        first_sine = box.first_sine.midpoint
        second_sine = box.second_sine.midpoint
        x_value = 1.0 - y - residual
        # This point only chooses tangents. Clip roundoff at the canonical
        # boundary so the untrusted numerical helper always sees x >= y.
        x_value = min(1.0, max(y, x_value))
        tangents: list[float] = []
        for kind, high, low, sine in (
            (self.first_kind, x_value, y, first_sine),
            (self.second_kind, 1.0 - y, 1.0 - x_value, second_sine),
        ):
            for label in (0, 1):
                def objective(state: float) -> float:
                    if kind == "rank":
                        return -rank_split_state_term(
                        high, low, sine, float(state), label,
                        self.weight_float, self.level_float,
                        )
                    return -endpoint_split_state_term(
                        high, low, float(state), label,
                        self.weight_float, self.level_float,
                    )

                optimum = minimize_scalar(
                    objective,
                    bounds=(0.0, 1.0),
                    method="bounded",
                    options={"xatol": 2e-12, "maxiter": 80},
                )
                state = float(optimum.x)
                probability = low + (high - low) * state
                probability = min(1.0 - 1e-14, max(1e-14, probability))
                tangents.append(math.sqrt((1.0 - probability) / probability))
        return tuple(tangents)

    def _endpoint_trial_residual_upper(
        self,
        high: Interval,
        low: Interval,
        label: int,
        tangent: float,
        trial: float,
    ) -> float:
        tangent_i = Interval.exact_float(tangent)
        inverse = divide(ONE, tangent_i)
        intercept = add(ONE, inverse)
        slope = sub(tangent_i, inverse)
        trial_i = Interval.exact_float(trial)
        coefficient = mul(trial_i, self.weight) if label == 0 else ZERO
        linear = add(slope, coefficient)
        constant_part = sub(intercept, mul(trial_i, self.level))
        at_high = add(constant_part, mul(linear, high))
        at_low = add(constant_part, mul(linear, low))
        return max(at_high.hi, at_low.hi)

    def _trial_residual_upper(
        self,
        high: Interval,
        low: Interval,
        sine: Interval,
        label: int,
        tangent: float,
        trial: float,
    ) -> float:
        tangent_i = Interval.exact_float(tangent)
        inverse = divide(ONE, tangent_i)
        intercept = add(ONE, inverse)
        slope = sub(tangent_i, inverse)
        trial_i = Interval.exact_float(trial)
        coefficient = mul(trial_i, self.weight)
        z_value = square(sine)
        one_minus_z = sub(ONE, z_value)
        first_split, second_split = (
            (one_minus_z, z_value) if label == 0 else (z_value, one_minus_z)
        )
        m11 = mul(high, add(slope, mul(coefficient, first_split)))
        m22 = mul(low, add(slope, mul(coefficient, second_split)))
        off_squared = mul(
            square(coefficient),
            mul(mul(high, low), mul(z_value, one_minus_z)),
        )
        trace = add(m11, m22)
        discriminant = add(
            square(sub(m11, m22)), scale(off_squared, 4.0)
        )
        eigenvalue = scale(add(trace, root(discriminant)), 0.5)
        residual = sub(
            add(intercept, eigenvalue), mul(trial_i, self.level)
        )
        return residual.hi

    def _trial_term_upper(
        self,
        high: Interval,
        low: Interval,
        sine: Interval,
        label: int,
        kind: str,
        tangent: float,
    ) -> float:
        def residual_upper(trial: float) -> float:
            if kind == "rank":
                return self._trial_residual_upper(
                    high, low, sine, label, tangent, trial
                )
            return self._endpoint_trial_residual_upper(
                high, low, label, tangent, trial
            )

        # Bracketing by doubling needs only logarithmically many rigorous
        # residual evaluations.  The former factor 1.05 spent roughly thirty
        # evaluations before bisection on the typical range of this proof.
        upper = 1.0
        while residual_upper(upper) > 0.0:
            upper *= 2.0
            if upper > 1e6:
                raise RuntimeError("failed to bracket a secular term")
        lower = 0.0
        # Thirty-six steps leave an absolute bracket below 1e-10 throughout
        # the certified regime.  Soundness does not depend on that estimate:
        # the outward interval residual is checked again below.
        for _ in range(36):
            trial = lower + (upper - lower) / 2.0
            if residual_upper(trial) <= 0.0:
                upper = trial
            else:
                lower = trial
        if residual_upper(upper) > 0.0:
            upper *= 1.0 + 1e-12
        if residual_upper(upper) > 0.0:
            raise ArithmeticError("could not preserve the trial bracket")
        return upper

    def _trial_sum_upper(
        self, box: RankRankBox, tangents: tuple[float, ...]
    ) -> float:
        y, residual, first_sine, second_sine = box.coordinates
        x_value = sub(sub(ONE, y), residual)
        groups = (
            (self.first_kind, x_value, y, first_sine),
            (self.second_kind, sub(ONE, y), add(y, residual), second_sine),
        )
        values: list[float] = []
        cursor = 0
        for kind, high, low, sine in groups:
            for label in (0, 1):
                values.append(
                    self._trial_term_upper(
                        high, low, sine, label, kind, tangents[cursor]
                    )
                )
                cursor += 1
        total = ZERO
        for value in values:
            total = add(total, Interval.exact_float(value))
        return mul(self.scale, total).hi

    def _endpoint_generalized_value(
        self,
        high: Interval,
        low: Interval,
        label: int,
        tangent: float,
    ) -> Interval:
        tangent_i = Interval.exact_float(tangent)
        inverse = divide(ONE, tangent_i)
        intercept = add(ONE, inverse)
        slope = sub(tangent_i, inverse)
        n_high = add(intercept, mul(slope, high))
        n_low = add(intercept, mul(slope, low))
        if label == 0:
            d_high = sub(self.level, mul(self.weight, high))
            d_low = sub(self.level, mul(self.weight, low))
        else:
            d_high = self.level
            d_low = self.level
        high_value = divide(n_high, d_high)
        low_value = divide(n_low, d_low)
        return Interval(
            max(high_value.lo, low_value.lo),
            max(high_value.hi, low_value.hi),
        )

    def _generalized_value(
        self,
        high: Interval,
        low: Interval,
        sine: Interval,
        label: int,
        tangent: float,
        kind: str = "rank",
    ) -> Interval:
        if kind == "endpoint":
            return self._endpoint_generalized_value(high, low, label, tangent)
        tangent_i = Interval.exact_float(tangent)
        inverse = divide(ONE, tangent_i)
        intercept = add(ONE, inverse)
        slope = sub(tangent_i, inverse)
        z_value = square(sine)
        one_minus_z = sub(ONE, z_value)
        first_split, second_split = (
            (one_minus_z, z_value) if label == 0 else (z_value, one_minus_z)
        )
        n1 = add(intercept, mul(slope, high))
        n2 = add(intercept, mul(slope, low))
        d11 = sub(self.level, mul(self.weight, mul(high, first_split)))
        d22 = sub(self.level, mul(self.weight, mul(low, second_split)))
        trace_g = add(mul(high, first_split), mul(low, second_split))
        determinant_d = sub(
            square(self.level), mul(mul(self.level, self.weight), trace_g)
        )
        cross = add(mul(n1, d22), mul(n2, d11))
        diagonal_gap = sub(mul(n1, d22), mul(n2, d11))
        off_squared = mul(
            square(self.weight),
            mul(mul(high, low), mul(z_value, one_minus_z)),
        )
        discriminant = add(
            square(diagonal_gap), scale(mul(mul(n1, n2), off_squared), 4.0)
        )
        return divide(add(cross, root(discriminant)), scale(determinant_d, 2.0))

    def _rank_diagonal_dominance_value(
        self,
        high: Interval,
        low: Interval,
        sine: Interval,
        label: int,
        tangent: float,
    ) -> Interval:
        """Upper-bound one rank term without a discriminant subtraction.

        Write the positive denominator as ``D=[[a,c],[c,b]]``.  Since
        ``D >= diag(a-|c|, b-|c|)``, the generalized Rayleigh quotient with
        diagonal positive numerator is at most the larger of the two
        diagonal ratios.  This is especially effective near either
        sine-zero face and remains regular at a diagonal-branch crossing.
        """

        tangent_i = Interval.exact_float(tangent)
        inverse = divide(ONE, tangent_i)
        intercept = add(ONE, inverse)
        slope = sub(tangent_i, inverse)
        z_value = square(sine)
        one_minus_z = sub(ONE, z_value)
        first_split, second_split = (
            (one_minus_z, z_value) if label == 0 else (z_value, one_minus_z)
        )
        n1 = add(intercept, mul(slope, high))
        n2 = add(intercept, mul(slope, low))
        d11 = sub(self.level, mul(self.weight, mul(high, first_split)))
        d22 = sub(self.level, mul(self.weight, mul(low, second_split)))
        off_squared = mul(
            square(self.weight),
            mul(mul(high, low), mul(z_value, one_minus_z)),
        )
        off = root(off_squared)
        first = divide(n1, sub(d11, off))
        second = divide(n2, sub(d22, off))
        return Interval(
            max(first.lo, second.lo), max(first.hi, second.hi)
        )

    def _endpoint_generalized_jet(
        self,
        high: Jet,
        low: Jet,
        label: int,
        tangent: float,
    ) -> Jet:
        tangent_j = constant(Interval.exact_float(tangent))
        inverse = jdivide(constant(ONE), tangent_j)
        intercept = jadd(constant(ONE), inverse)
        slope = jsub(tangent_j, inverse)
        n_high = jadd(intercept, jmul(slope, high))
        n_low = jadd(intercept, jmul(slope, low))
        if label == 0:
            d_high = jsub(
                constant(self.level), jmul(constant(self.weight), high)
            )
            d_low = jsub(
                constant(self.level), jmul(constant(self.weight), low)
            )
        else:
            d_high = constant(self.level)
            d_low = constant(self.level)
        return jmaximum(
            jdivide(n_high, d_high),
            jdivide(n_low, d_low),
        )

    def _rank_zero_generalized_jet(
        self,
        high: Jet,
        low: Jet,
        label: int,
        tangent: float,
    ) -> Jet:
        """Exact diagonal rank-split jet on the sine-zero face."""

        tangent_j = constant(Interval.exact_float(tangent))
        inverse = jdivide(constant(ONE), tangent_j)
        intercept = jadd(constant(ONE), inverse)
        slope = jsub(tangent_j, inverse)
        n_high = jadd(intercept, jmul(slope, high))
        n_low = jadd(intercept, jmul(slope, low))
        if label == 0:
            d_high = jsub(
                constant(self.level), jmul(constant(self.weight), high)
            )
            d_low = constant(self.level)
        else:
            d_high = constant(self.level)
            d_low = jsub(
                constant(self.level), jmul(constant(self.weight), low)
            )
        return jmaximum(
            jdivide(n_high, d_high),
            jdivide(n_low, d_low),
        )

    def _generalized_jet(
        self,
        high: Jet,
        low: Jet,
        sine: Jet,
        label: int,
        tangent: float,
        kind: str = "rank",
    ) -> Jet:
        if kind == "endpoint":
            return self._endpoint_generalized_jet(high, low, label, tangent)
        tangent_j = constant(Interval.exact_float(tangent))
        inverse = jdivide(constant(ONE), tangent_j)
        intercept = jadd(constant(ONE), inverse)
        slope = jsub(tangent_j, inverse)
        z_value = jsquare(sine)
        one_minus_z = jsub(constant(ONE), z_value)
        first_split, second_split = (
            (one_minus_z, z_value) if label == 0 else (z_value, one_minus_z)
        )
        n1 = jadd(intercept, jmul(slope, high))
        n2 = jadd(intercept, jmul(slope, low))
        d11 = jsub(
            constant(self.level), jmul(constant(self.weight), jmul(high, first_split))
        )
        d22 = jsub(
            constant(self.level), jmul(constant(self.weight), jmul(low, second_split))
        )
        trace_g = jadd(jmul(high, first_split), jmul(low, second_split))
        determinant_d = jsub(
            jsquare(constant(self.level)),
            jmul(jmul(constant(self.level), constant(self.weight)), trace_g),
        )
        cross = jadd(jmul(n1, d22), jmul(n2, d11))
        diagonal_gap = jsub(jmul(n1, d22), jmul(n2, d11))
        off_squared = jmul(
            jsquare(constant(self.weight)),
            jmul(jmul(high, low), jmul(z_value, one_minus_z)),
        )
        discriminant = jadd(
            jsquare(diagonal_gap),
            jmul(constant(4.0), jmul(jmul(n1, n2), off_squared)),
        )
        return jdivide(
            jadd(cross, jroot(discriminant)), jmul(constant(2.0), determinant_d)
        )

    def _generalized_jet2(
        self,
        high: Jet2,
        low: Jet2,
        sine: Jet2,
        label: int,
        tangent: float,
    ) -> Jet2:
        """Second-order jet of one rank-split generalized eigenvalue."""

        tangent_j = constant2(Interval.exact_float(tangent))
        inverse = j2divide(constant2(ONE), tangent_j)
        intercept = j2add(constant2(ONE), inverse)
        slope = j2sub(tangent_j, inverse)
        z_value = j2square(sine)
        one_minus_z = j2sub(constant2(ONE), z_value)
        first_split, second_split = (
            (one_minus_z, z_value) if label == 0 else (z_value, one_minus_z)
        )
        n1 = j2add(intercept, j2mul(slope, high))
        n2 = j2add(intercept, j2mul(slope, low))
        d11 = j2sub(
            constant2(self.level),
            j2mul(constant2(self.weight), j2mul(high, first_split)),
        )
        d22 = j2sub(
            constant2(self.level),
            j2mul(constant2(self.weight), j2mul(low, second_split)),
        )
        trace_g = j2add(j2mul(high, first_split), j2mul(low, second_split))
        determinant_d = j2sub(
            j2square(constant2(self.level)),
            j2mul(
                j2mul(constant2(self.level), constant2(self.weight)), trace_g
            ),
        )
        cross = j2add(j2mul(n1, d22), j2mul(n2, d11))
        diagonal_gap = j2sub(j2mul(n1, d22), j2mul(n2, d11))
        off_squared = j2mul(
            j2square(constant2(self.weight)),
            j2mul(j2mul(high, low), j2mul(z_value, one_minus_z)),
        )
        discriminant = j2add(
            j2square(diagonal_gap),
            j2mul(
                constant2(4.0), j2mul(j2mul(n1, n2), off_squared)
            ),
        )
        return j2divide(
            j2add(cross, j2root(discriminant)),
            j2mul(constant2(2.0), determinant_d),
        )

    def _value_expression(
        self, box: RankRankBox, tangents: tuple[float, ...]
    ) -> tuple[Interval, tuple[Interval, ...]]:
        y, residual, first_sine, second_sine = box.coordinates
        x_value = sub(sub(ONE, y), residual)
        groups = (
            (self.first_kind, x_value, y, first_sine),
            (self.second_kind, sub(ONE, y), add(y, residual), second_sine),
        )
        terms: list[Interval] = []
        cursor = 0
        for kind, high, low, sine in groups:
            for label in (0, 1):
                terms.append(
                    self._generalized_value(
                        high, low, sine, label, tangents[cursor], kind
                    )
                )
                cursor += 1
        total = ZERO
        for term in terms:
            total = add(total, term)
        return mul(self.scale, total), tuple(terms)

    def _jet_expression(
        self, box: RankRankBox, tangents: tuple[float, ...]
    ) -> tuple[Jet, tuple[Jet, ...]]:
        y = variable(box.y, 0)
        residual = variable(box.residual, 1)
        first_sine = variable(box.first_sine, 2)
        second_sine = variable(box.second_sine, 3)
        x_value = jsub(jsub(constant(ONE), y), residual)
        groups = (
            (self.first_kind, x_value, y, first_sine),
            (self.second_kind, jsub(constant(ONE), y), jadd(y, residual), second_sine),
        )
        terms: list[Jet] = []
        cursor = 0
        for kind, high, low, sine in groups:
            for label in (0, 1):
                terms.append(
                    self._generalized_jet(
                        high, low, sine, label, tangents[cursor], kind
                    )
                )
                cursor += 1
        total = constant(ZERO)
        for term in terms:
            total = jadd(total, term)
        return jmul(constant(self.scale), total), tuple(terms)

    def _diagonal_dominance_sum_upper(
        self, box: RankRankBox, tangents: tuple[float, ...]
    ) -> float:
        y, residual, first_sine, second_sine = box.coordinates
        x_value = sub(sub(ONE, y), residual)
        groups = (
            (self.first_kind, x_value, y, first_sine),
            (self.second_kind, sub(ONE, y), add(y, residual), second_sine),
        )
        total = ZERO
        cursor = 0
        for kind, high, low, sine in groups:
            for label in (0, 1):
                if kind == "rank":
                    term = self._rank_diagonal_dominance_value(
                        high, low, sine, label, tangents[cursor]
                    )
                else:
                    term = self._endpoint_generalized_value(
                        high, low, label, tangents[cursor]
                    )
                total = add(total, term)
                cursor += 1
        return mul(self.scale, total).hi

    def _rank_center_lipschitz_value(
        self,
        high: Interval,
        low: Interval,
        sine: Interval,
        label: int,
        tangent: float,
    ) -> Interval:
        """Bound one rank term by a regular matrix perturbation estimate.

        For ``lambda(N,D) = max_v (v* N v)/(v* D v)``, anchor the two
        matrices at the center of the box.  If ``D,D0 >= d_* I``, then

            lambda(N,D)
              <= lambda(N0,D0) (1 + ||D-D0||/d_*)
                 + ||N-N0||/d_*.

        This follows directly by comparing the two Rayleigh quotients.  It
        avoids differentiating the closed-form square root and therefore
        stays finite when the two generalized eigenvalues cross.  Every
        quantity below is enclosed with directed binary64 intervals.
        """

        high_center = Interval.exact_float(high.midpoint)
        low_center = Interval.exact_float(low.midpoint)
        sine_center = Interval.exact_float(sine.midpoint)
        anchor = self._generalized_value(
            high_center, low_center, sine_center, label, tangent
        )

        tangent_i = Interval.exact_float(tangent)
        inverse = divide(ONE, tangent_i)
        intercept = add(ONE, inverse)
        slope = sub(tangent_i, inverse)

        def numerator_components(
            high_value: Interval, low_value: Interval
        ) -> tuple[Interval, Interval]:
            return (
                add(intercept, mul(slope, high_value)),
                add(intercept, mul(slope, low_value)),
            )

        def denominator_components(
            high_value: Interval,
            low_value: Interval,
            sine_value: Interval,
        ) -> tuple[Interval, Interval, Interval]:
            z_value = square(sine_value)
            one_minus_z = sub(ONE, z_value)
            first_split, second_split = (
                (one_minus_z, z_value)
                if label == 0
                else (z_value, one_minus_z)
            )
            d11 = sub(
                self.level,
                mul(self.weight, mul(high_value, first_split)),
            )
            d22 = sub(
                self.level,
                mul(self.weight, mul(low_value, second_split)),
            )
            off = root(
                mul(
                    square(self.weight),
                    mul(
                        mul(high_value, low_value),
                        mul(z_value, one_minus_z),
                    ),
                )
            )
            return d11, d22, off

        n1, n2 = numerator_components(high, low)
        n1_center, n2_center = numerator_components(high_center, low_center)
        delta_n1 = sub(n1, n1_center)
        delta_n2 = sub(n2, n2_center)
        delta_n = Interval(
            0.0,
            up(
                max(
                    abs(delta_n1.lo),
                    abs(delta_n1.hi),
                    abs(delta_n2.lo),
                    abs(delta_n2.hi),
                )
            ),
        )

        d11, d22, off = denominator_components(high, low, sine)
        d11_center, d22_center, off_center = denominator_components(
            high_center, low_center, sine_center
        )
        delta_d11 = sub(d11, d11_center)
        delta_d22 = sub(d22, d22_center)
        delta_off = sub(off, off_center)
        # The spectral norm is bounded by the Frobenius norm.  Both
        # off-diagonal entries contribute the same squared difference.
        delta_d = root(
            add(
                add(square(delta_d11), square(delta_d22)),
                scale(square(delta_off), 2.0),
            )
        )

        # The rank-one split matrix is positive with its only nonzero
        # eigenvalue no larger than max(high, low).  Hence every denominator
        # in the box is at least this positive scalar multiple of identity.
        spectral_cap = Interval(0.0, up(max(high.hi, low.hi)))
        denominator_floor = sub(
            self.level, mul(self.weight, spectral_cap)
        )
        if denominator_floor.lo <= 0.0:
            raise ZeroDivisionError("matrix denominator floor is not positive")

        anchor_upper = Interval(0.0, up(max(0.0, anchor.hi)))
        relative_d = divide(delta_d, denominator_floor)
        perturbation = add(
            mul(anchor_upper, relative_d),
            divide(delta_n, denominator_floor),
        )
        return add(anchor, perturbation)

    def _center_matrix_lipschitz_sum_upper(
        self, box: RankRankBox, tangents: tuple[float, ...]
    ) -> float:
        """Regular center-anchored upper bound for the complete topology."""

        y, residual, first_sine, second_sine = box.coordinates
        x_value = sub(sub(ONE, y), residual)
        groups = (
            (self.first_kind, x_value, y, first_sine),
            (self.second_kind, sub(ONE, y), add(y, residual), second_sine),
        )
        total = ZERO
        cursor = 0
        for kind, high, low, sine in groups:
            for label in (0, 1):
                if kind == "rank":
                    term = self._rank_center_lipschitz_value(
                        high, low, sine, label, tangents[cursor]
                    )
                else:
                    term = self._endpoint_generalized_value(
                        high, low, label, tangents[cursor]
                    )
                total = add(total, term)
                cursor += 1
        return mul(self.scale, total).hi

    def _zero_face_jet_expression(
        self,
        box: RankRankBox,
        tangents: tuple[float, ...],
        face: int = 1,
    ) -> Jet:
        """Jet of the full sum with either rank split diagonalised."""

        if self.first_kind != "rank" or self.second_kind != "rank":
            raise ValueError("zero-face form currently covers rank/rank only")
        if face not in (0, 1):
            raise ValueError("zero face must be the first or second split")
        y = variable(box.y, 0)
        residual = variable(box.residual, 1)
        sines = (variable(box.first_sine, 2), variable(box.second_sine, 3))
        groups = (
            (jsub(jsub(constant(ONE), y), residual), y),
            (jsub(constant(ONE), y), jadd(y, residual)),
        )
        terms: list[Jet] = []
        for group, ((high, low), sine) in enumerate(zip(groups, sines)):
            for label in (0, 1):
                tangent = tangents[2 * group + label]
                if group == face:
                    terms.append(
                        self._rank_zero_generalized_jet(
                            high, low, label, tangent
                        )
                    )
                else:
                    terms.append(
                        self._generalized_jet(
                            high, low, sine, label, tangent
                        )
                    )
        total = constant(ZERO)
        for term in terms:
            total = jadd(total, term)
        return jmul(constant(self.scale), total)

    def _zero_face_lipschitz_upper(
        self, box: RankRankBox, face: int = 1
    ) -> float:
        """Bound a box near either sine-zero face without ``sqrt(Delta)``.

        On the face the second generalized eigenproblems are diagonal.  Away
        from it, the Rayleigh quotient changes by at most

            n_max * weight * ||G(s)-G(0)|| / d_min**2.

        The Frobenius norm gives an explicit outward interval for the matrix
        difference.  This remains finite when the two diagonal branches cross.
        """

        if face not in (0, 1):
            raise ValueError("zero face must be the first or second split")
        coordinates = list(box.coordinates)
        coordinates[2 + face] = ZERO
        zero_box = RankRankBox(*coordinates)
        selected = self._select_tangents(zero_box)
        enclosure = self._zero_face_jet_expression(zero_box, selected, face)
        anchors: list[float] = []
        deltas: list[Interval] = []
        for coordinate, derivative in zip(
            zero_box.coordinates, enclosure.gradient
        ):
            if derivative.hi < 0.0:
                anchors.append(coordinate.lo)
                deltas.append(Interval(0.0, up(coordinate.width)))
            elif derivative.lo > 0.0:
                anchors.append(coordinate.hi)
                deltas.append(Interval(down(-coordinate.width), 0.0))
            else:
                anchor = coordinate.midpoint
                anchors.append(anchor)
                deltas.append(
                    Interval(
                        down(coordinate.lo - anchor),
                        up(coordinate.hi - anchor),
                    )
                )
        anchor_box = RankRankBox(
            *(Interval.exact_float(value) for value in anchors)
        )
        mean_value = self._zero_face_jet_expression(
            anchor_box, selected, face
        ).value
        for derivative, delta in zip(enclosure.gradient, deltas):
            mean_value = add(mean_value, mul(derivative, delta))

        y, residual, first_sine, second_sine = box.coordinates
        if face == 0:
            high = sub(sub(ONE, y), residual)
            low = y
            sine = first_sine
        else:
            high = sub(ONE, y)
            low = add(y, residual)
            sine = second_sine
        z_value = square(sine)
        z_squared = square(z_value)
        delta_norm_squared = add(
            mul(add(square(high), square(low)), z_squared),
            scale(
                mul(
                    mul(high, low),
                    mul(z_value, sub(ONE, z_value)),
                ),
                2.0,
            ),
        )
        delta_norm = root(delta_norm_squared)
        denominator_floor = sub(self.level, mul(self.weight, high))
        numerator_sum = ZERO
        for tangent in selected[2 * face : 2 * face + 2]:
            tangent_i = Interval.exact_float(tangent)
            inverse = divide(ONE, tangent_i)
            intercept = add(ONE, inverse)
            slope = sub(tangent_i, inverse)
            n_high = add(intercept, mul(slope, high))
            n_low = add(intercept, mul(slope, low))
            numerator_sum = add(
                numerator_sum,
                Interval(0.0, up(max(n_high.hi, n_low.hi))),
            )
        perturbation = divide(
            mul(mul(self.weight, delta_norm), numerator_sum),
            square(denominator_floor),
        )
        return add(mean_value, mul(self.scale, perturbation)).hi

    def _jet2_expression(
        self, box: RankRankBox, tangents: tuple[float, ...]
    ) -> tuple[Jet2, tuple[Jet2, ...]]:
        if self.first_kind != "rank" or self.second_kind != "rank":
            raise ValueError("second-order form currently covers rank/rank only")
        y = variable2(box.y, 0)
        residual = variable2(box.residual, 1)
        first_sine = variable2(box.first_sine, 2)
        second_sine = variable2(box.second_sine, 3)
        x_value = j2sub(j2sub(constant2(ONE), y), residual)
        groups = (
            (x_value, y, first_sine),
            (j2sub(constant2(ONE), y), j2add(y, residual), second_sine),
        )
        terms: list[Jet2] = []
        cursor = 0
        for high, low, sine in groups:
            for label in (0, 1):
                terms.append(
                    self._generalized_jet2(
                        high, low, sine, label, tangents[cursor]
                    )
                )
                cursor += 1
        total = constant2(ZERO)
        for term in terms:
            total = j2add(total, term)
        return j2mul(constant2(self.scale), total), tuple(terms)

    def bound(
        self, box: RankRankBox, tangents: tuple[float, ...] | None = None
    ) -> BoundResult:
        contracted = self.contract(box)
        if contracted is None:
            return BoundResult(-math.inf, "domain-empty", (), None)
        selected = tangents or self._select_tangents(contracted)
        try:
            natural, _ = self._value_expression(contracted, selected)
            natural_upper = natural.hi
        except (ValueError, ZeroDivisionError):
            natural_upper = math.inf
        try:
            trial_upper = self._trial_sum_upper(contracted, selected)
        except (ArithmeticError, RuntimeError):
            trial_upper = math.inf
        upper = min(natural_upper, trial_upper)
        base_method = "trial" if trial_upper < natural_upper else "natural"
        if upper <= 1.0:
            return BoundResult(upper, base_method, selected, None)
        try:
            diagonal_upper = self._diagonal_dominance_sum_upper(
                contracted, selected
            )
            if diagonal_upper < upper:
                upper = diagonal_upper
                base_method = "diagonal-dominance"
        except (ArithmeticError, ValueError, ZeroDivisionError):
            pass
        if upper <= 1.0:
            return BoundResult(upper, base_method, selected, None)
        try:
            matrix_lipschitz_upper = self._center_matrix_lipschitz_sum_upper(
                contracted, selected
            )
            if matrix_lipschitz_upper < upper:
                upper = matrix_lipschitz_upper
                base_method = "center-matrix-lipschitz"
        except (ArithmeticError, ValueError, ZeroDivisionError):
            pass
        if upper <= 1.0:
            return BoundResult(upper, base_method, selected, None)
        for face in (0, 1):
            try:
                zero_face_upper = self._zero_face_lipschitz_upper(contracted, face)
                if zero_face_upper < upper:
                    upper = zero_face_upper
                    base_method = "zero-face-lipschitz"
            except (ArithmeticError, ValueError, ZeroDivisionError):
                pass
        if upper <= 1.0:
            return BoundResult(upper, base_method, selected, None)
        try:
            enclosure, _ = self._jet_expression(contracted, selected)
            anchors: list[float] = []
            deltas: list[Interval] = []
            for coordinate, derivative in zip(
                contracted.coordinates, enclosure.gradient
            ):
                if derivative.hi < 0.0:
                    anchors.append(coordinate.lo)
                    deltas.append(Interval(0.0, up(coordinate.width)))
                elif derivative.lo > 0.0:
                    anchors.append(coordinate.hi)
                    deltas.append(Interval(down(-coordinate.width), 0.0))
                else:
                    anchor = coordinate.midpoint
                    anchors.append(anchor)
                    deltas.append(
                        Interval(
                            down(coordinate.lo - anchor),
                            up(coordinate.hi - anchor),
                        )
                    )
            anchor_box = RankRankBox(
                *(Interval.exact_float(value) for value in anchors)
            )
            anchor_value, _ = self._value_expression(anchor_box, selected)
            mean_value = anchor_value
            for derivative, delta in zip(enclosure.gradient, deltas):
                mean_value = add(mean_value, mul(derivative, delta))
            best_upper = upper
            best_method = base_method
            if mean_value.hi < best_upper:
                best_upper = mean_value.hi
                best_method = "mean-value"
            if best_upper <= 1.0:
                return BoundResult(
                    best_upper, best_method, selected, enclosure.gradient
                )

            # Near the narrow stationary frontier, evaluate value and gradient
            # at one exact dyadic center and enclose the Hessian on the full
            # box.  The complete i,j sum followed by 1/2 is the standard
            # multivariate Taylor remainder, including all mixed derivatives.
            try:
                centered = [
                    coordinate.midpoint for coordinate in contracted.coordinates
                ]
                center_box = RankRankBox(
                    *(Interval.exact_float(value) for value in centered)
                )
                # Only the value and first derivatives are needed at the
                # exact dyadic centre.  A second-order jet here duplicated
                # the expensive Hessian algebra that is required only on the
                # full box below.
                center_jet, _ = self._jet_expression(center_box, selected)
                box_jet, _ = self._jet2_expression(contracted, selected)
                center_deltas = [
                    Interval(
                        down(coordinate.lo - center),
                        up(coordinate.hi - center),
                    )
                    for coordinate, center in zip(
                        contracted.coordinates, centered
                    )
                ]
                taylor = center_jet.value
                for derivative, delta in zip(
                    center_jet.gradient, center_deltas
                ):
                    taylor = add(taylor, mul(derivative, delta))
                remainder = ZERO
                for i in range(DIMENSION):
                    for j in range(DIMENSION):
                        displacement = (
                            square(center_deltas[i])
                            if i == j
                            else mul(center_deltas[i], center_deltas[j])
                        )
                        remainder = add(
                            remainder,
                            mul(box_jet.hessian[i][j], displacement),
                        )
                taylor = add(taylor, scale(remainder, 0.5))
                if taylor.hi < best_upper:
                    best_upper = taylor.hi
                    best_method = "second-order"
            except (ValueError, ZeroDivisionError):
                pass
            return BoundResult(
                best_upper, best_method, selected, enclosure.gradient
            )
        except (ValueError, ZeroDivisionError):
            return BoundResult(upper, base_method, selected, None)

    def split(
        self,
        box: RankRankBox,
        scales: tuple[float, ...],
        tangents: tuple[float, ...] | None = None,
    ) -> tuple[str, RankRankBox, RankRankBox]:
        contracted = self.contract(box)
        if contracted is None:
            raise ValueError("cannot split an empty box")
        names = ("y", "residual", "first_sine", "second_sine")
        coordinates = contracted.coordinates
        bisectable = tuple(
            coordinate.lo < coordinate.midpoint < coordinate.hi
            for coordinate in coordinates
        )
        if not any(bisectable):
            raise FloatingPointError("no representable interior midpoint remains")
        try:
            selected = tangents or self._select_tangents(contracted)
            enclosure, _ = self._jet_expression(contracted, selected)
            # A sign-definite derivative has already moved the mean-value
            # anchor to its maximizing endpoint, so bisecting that coordinate
            # cannot sharpen the upper contribution.  Concentrate the tree on
            # genuinely stationary coordinates whose derivative interval
            # crosses zero.
            scores = [
                0.0
                if not bisectable[index]
                or derivative.hi < 0.0
                or derivative.lo > 0.0
                else coordinate.width
                * max(abs(derivative.lo), abs(derivative.hi))
                for index, (coordinate, derivative) in enumerate(
                    zip(coordinates, enclosure.gradient)
                )
            ]
            if max(scores) == 0.0:
                index = max(
                    (item for item in range(DIMENSION) if bisectable[item]),
                    key=lambda item: coordinates[item].width / scales[item],
                )
            else:
                index = max(range(DIMENSION), key=lambda item: scores[item])
        except (ValueError, ZeroDivisionError):
            # At an eigenvalue crossing the raw square-root jet is singular,
            # although the matrix problem is Lipschitz.  Use the diagonal
            # zero-face jet for the regular sensitivities and compare the
            # regularised bound with its face value for the singular angular
            # score.  Both rank groups are tried symmetrically.
            # This prevents endless bisection toward a representable zero.
            try:
                if self.first_kind != "rank" or self.second_kind != "rank":
                    raise ValueError("zero-face split applies to rank/rank")
                candidates: list[tuple[float, float, int, list[float]]] = []
                for face in (0, 1):
                    try:
                        zero_coordinates = list(coordinates)
                        zero_coordinates[2 + face] = ZERO
                        zero_box = RankRankBox(*zero_coordinates)
                        selected = self._select_tangents(zero_box)
                        enclosure = self._zero_face_jet_expression(
                            zero_box, selected, face
                        )
                        scores = [
                            (
                                coordinate.width
                                * max(abs(derivative.lo), abs(derivative.hi))
                                if bisectable[index] else 0.0
                            )
                            for index, (coordinate, derivative) in enumerate(
                                zip(coordinates, enclosure.gradient)
                            )
                        ]
                        face_upper = self._zero_face_lipschitz_upper(
                            zero_box, face
                        )
                        near_upper = self._zero_face_lipschitz_upper(
                            contracted, face
                        )
                        scores[2 + face] = max(0.0, near_upper - face_upper)
                        candidates.append(
                            (near_upper, face_upper, face, scores)
                        )
                    except (ArithmeticError, ValueError, ZeroDivisionError):
                        continue
                if not candidates:
                    raise ValueError("no regular zero-face split candidate")
                critical = [
                    item for item in candidates
                    if bisectable[2 + item[2]]
                    and coordinates[2 + item[2]].lo == 0.0
                    and item[1] < 1.0 <= item[0]
                ]
                if critical:
                    _, _, face, _ = max(
                        critical, key=lambda item: item[0] - item[1]
                    )
                    index = 2 + face
                else:
                    _, _, _, scores = min(candidates, key=lambda item: item[0])
                    if max(scores) == 0.0:
                        raise ValueError("zero-face split scores vanished")
                    index = max(range(DIMENSION), key=lambda item: scores[item])
            except (ArithmeticError, ValueError, ZeroDivisionError):
                index = max(
                    (item for item in range(DIMENSION) if bisectable[item]),
                    key=lambda item: coordinates[item].width / scales[item],
                )
        # A square-root eigenvalue jet can have an arbitrarily large interval
        # derivative on a crossing surface.  Do not let that heuristic refine
        # one coordinate many orders of magnitude beyond the rest: the trial
        # and matrix-Lipschitz bounds converge only when the whole box does.
        normalized_widths = tuple(
            coordinates[item].width / scales[item] for item in range(DIMENSION)
        )
        widest = max(
            (item for item in range(DIMENSION) if bisectable[item]),
            key=lambda item: normalized_widths[item],
        )
        if normalized_widths[index] < normalized_widths[widest] / 16.0:
            index = widest
        selected_coordinate = coordinates[index]
        middle = selected_coordinate.midpoint
        if not selected_coordinate.lo < middle < selected_coordinate.hi:
            index = max(
                (item for item in range(DIMENSION) if bisectable[item]),
                key=lambda item: coordinates[item].width / scales[item],
            )
            selected_coordinate = coordinates[index]
            middle = selected_coordinate.midpoint
        left = list(coordinates)
        right = list(coordinates)
        left[index] = Interval(selected_coordinate.lo, middle)
        right[index] = Interval(middle, selected_coordinate.hi)
        children = (RankRankBox(*left), RankRankBox(*right))
        contracted_children = tuple(self.contract(child) or child for child in children)
        return names[index], contracted_children[0], contracted_children[1]

    def certify(
        self,
        root: RankRankBox,
        max_boxes: int,
        *,
        resume: dict[str, object] | None = None,
    ) -> dict[str, object]:
        root = self.contract(root) or root
        scales = tuple(max(item.width, 1e-15) for item in root.coordinates)
        method_names = (
            "mean-value",
            "sqrt-y-mean-value",
            "second-order",
            "zero-face-lipschitz",
            "center-matrix-lipschitz",
            "trial",
            "natural",
            "diagonal-dominance",
            "domain-empty",
        )
        saved_frontier = None if resume is None else resume.get("open_frontier")
        if saved_frontier is None:
            first = self.bound(root)
            initial_upper = first.upper
            queue: list[tuple[float, int, RankRankBox]] = []
            if first.upper > 1.0:
                queue.append((-first.upper, 0, root))
            counter = 0
            split_count = 0
            closed = 1 if first.upper <= 1.0 else 0
            methods = {name: 0 for name in method_names}
            if first.upper <= 1.0:
                methods[first.method] += 1
        else:
            if resume is None:
                raise RuntimeError("unreachable resume state")
            if resume.get("weight") != str(self.weight_fraction):
                raise ValueError("resume certificate has the wrong weight")
            if resume.get("level") != str(self.level_fraction):
                raise ValueError("resume certificate has the wrong level")
            if resume.get("root_box") != root.serialise():
                raise ValueError("resume certificate has the wrong root box")
            if bool(resume.get("complete")):
                raise ValueError("complete certificates should be reused directly")
            if not isinstance(saved_frontier, list) or not saved_frontier:
                raise ValueError("incomplete resume certificate has no frontier")
            initial_upper = float(resume["initial_upper"])
            split_count = int(resume["boxes_split"])
            closed = int(resume["boxes_closed"])
            recorded_methods = resume["closed_methods"]
            if not isinstance(recorded_methods, dict):
                raise ValueError("resume certificate has malformed method counts")
            methods = {
                name: int(recorded_methods.get(name, 0))
                for name in method_names
            }
            queue = []
            for counter, item in enumerate(saved_frontier):
                if not isinstance(item, dict):
                    raise ValueError("resume frontier entry is malformed")
                upper = float(item["upper"])
                box = RankRankBox.deserialise(item["box"])
                if not math.isfinite(upper) or upper <= 1.0:
                    raise ValueError("resume frontier contains a closed bound")
                heapq.heappush(queue, (-upper, counter, box))
            counter = len(queue)
        while queue and split_count < max_boxes:
            _, _, box = heapq.heappop(queue)
            selected = self._select_tangents(box)
            _, left, right = self.split(box, scales, selected)
            for child in (left, right):
                result = self.bound(child, selected)
                if result.upper <= 1.0:
                    closed += 1
                    methods[result.method] += 1
                elif math.isfinite(result.upper):
                    counter += 1
                    heapq.heappush(queue, (-result.upper, counter, child))
            split_count += 1
            if split_count % 1000 == 0:
                print(
                    json.dumps(
                        {
                            "boxes_split": split_count,
                            "open": len(queue),
                            "maximum_open_upper": -queue[0][0] if queue else -math.inf,
                            "leading_box": queue[0][2].serialise() if queue else None,
                        }
                    ),
                    flush=True,
                )
        return {
            "weight": str(self.weight_fraction),
            "level": str(self.level_fraction),
            "topology": f"{self.first_kind}/{self.second_kind}",
            "state_elimination": "AM-GM tangent and 2-by-2 generalized eigenvalue",
            "root_box": root.serialise(),
            "initial_upper": initial_upper,
            "boxes_split": split_count,
            "boxes_closed": closed,
            "closed_methods": methods,
            "boxes_remaining": len(queue),
            "maximum_open_upper": -queue[0][0] if queue else -math.inf,
            "leading_box": queue[0][2].serialise() if queue else None,
            "open_frontier": [
                {"upper": -priority, "box": box.serialise()}
                for priority, _, box in sorted(
                    queue, key=lambda item: (item[0], item[1])
                )
            ],
            "complete": not queue,
            "proof_kernel": "outward-expanded IEEE-754 binary64 intervals",
            "trusted_optimizers": [],
            "untrusted_search_helpers": ["SciPy tangent selection", "branch heuristic"],
        }


def decimal_box(lower: str | float, upper: str | float) -> Interval:
    low = Interval.decimal(lower).lo
    high = Interval.decimal(upper).hi
    return Interval(low, high)


def root_from_bounds(
    x_bounds: Interval,
    y_bounds: Interval,
    first_sine: Interval,
    second_sine: Interval,
) -> RankRankBox:
    residual_lower = max(0.0, sub(sub(ONE, x_bounds), y_bounds).lo)
    residual_upper = max(0.0, sub(sub(ONE, x_bounds), y_bounds).hi)
    return RankRankBox(
        y_bounds,
        Interval(residual_lower, residual_upper),
        first_sine,
        second_sine,
    )


def parse_pair(values: Iterable[str]) -> Interval:
    lower, upper = tuple(values)
    return decimal_box(lower, upper)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lambda", dest="weight", default="0.6")
    parser.add_argument("--level", default="0.76662")
    parser.add_argument("--x", nargs=2, default=("0.94", "0.96"))
    parser.add_argument("--y", nargs=2, default=("0.04", "0.06"))
    parser.add_argument("--first-sine", nargs=2, default=("0.03", "0.06"))
    parser.add_argument("--second-sine", nargs=2, default=("0", "0.03"))
    parser.add_argument("--max-boxes", type=int, default=3000)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    x_bounds = parse_pair(args.x)
    y_bounds = parse_pair(args.y)
    first_sine = parse_pair(args.first_sine)
    second_sine = parse_pair(args.second_sine)
    certifier = RankRankCertifier(
        Fraction(args.weight), Fraction(args.level), x_bounds
    )
    payload = certifier.certify(
        root_from_bounds(x_bounds, y_bounds, first_sine, second_sine),
        args.max_boxes,
    )
    rendered = json.dumps(payload, indent=2) + "\n"
    print(rendered, end="")
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

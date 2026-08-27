"""Exact-dual McCormick--SOCP cover of the four-active readout sector.

This strengthens :mod:`four_active_socp_exact_cover` by retaining one shared
weight vector inside every conic relaxation.  The only nonlinear products are

``w_i p_i``, ``w_0 g`` and ``k_i(A-p_i)``,

where ``g`` is the ordered-prefix gap and
``k_i=(1-w_i)/(2-w_i)``.  A fourth product ``k_i w_i`` encodes the latter
linear-fractional identity as ``2 k_i-k_i w_i=1-w_i``.  Every product receives
all four McCormick facets on an exact rational weight box.  Thus the model is
an SOCP, converges to the original convex consequences under bisection, and
admits the same solver-free exact-residual dual checker as the ternary cover.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from functools import cache
from fractions import Fraction
import heapq
import itertools
import json
import math
from pathlib import Path
from typing import Any

import cvxpy as cp
from cvxpy.reductions.solvers.conic_solvers.clarabel_conif import CLARABEL
import numpy as np

from four_active_socp_exact_cover import (
    MAXIMUM_WEIGHT_FLOOR,
    MINIMUM_ACTIVE_WEIGHT,
    NONZERO_PERMUTATIONS,
    PROJECTIVE_LINES,
    SUPPORT_WEIGHT,
    WeightBox,
    binary64_down,
    binary64_up,
    choose_split,
    decode_dual,
    encode_dual,
    fraction_pair,
    initial_box,
    polytope_vertices,
    split_box,
    validate_leaf_tree,
    weight_hull,
)
from ternary_socp_exact_dual_probe import (
    canonical_hash,
    exact_dot,
    exact_sparse_stationarity,
    fraction_decimal,
    q,
    repair_dual_cones,
)


ROOT = Path(__file__).resolve().parent
SCHEMA = "carmenq.four-active-mccormick-socp-exact-dual-cover.v26"
FOUR_ACTIVE_TARGET = Fraction(76670, 100000)
PERTURBATIONS = (Fraction(1, 20), Fraction(3, 10))
PERTURBATION = PERTURBATIONS[0]
ORDERED_EFFECT_PAIRS = tuple(itertools.permutations(range(4), 2))
PERTURBATION_PAIRS = ORDERED_EFFECT_PAIRS
# Three small rational support directions found by common-bias separation.
# Their complete orbits are label neutral.  Unlike subset-reserve cuts, every
# member couples all four priors and detects incompatible simultaneous uses
# of the single Helstrom spectral bias.
COMMON_BIAS_COEFFICIENT_REPRESENTATIVES = (
    (Fraction(1), Fraction(1, 3), Fraction(2, 5), Fraction(2, 5)),
    (Fraction(1), Fraction(1, 2), Fraction(1, 3), Fraction(1, 3)),
    (Fraction(1), Fraction(10, 29), Fraction(1, 2), Fraction(16, 27)),
)
COMMON_BIAS_COEFFICIENTS = tuple(
    sorted(
        {
            coefficients
            for representative in COMMON_BIAS_COEFFICIENT_REPRESENTATIVES
            for coefficients in itertools.permutations(representative)
        }
    )
)



def loss_factor(weight: Fraction) -> Fraction:
    return (1 - weight) / (2 - weight)



def sqrt_fraction_down(value: Fraction) -> Fraction:
    """Return a dyadic lower bound on a nonnegative rational square root."""

    if value < 0:
        raise ValueError("cannot take the square root of a negative rational")
    candidate = math.sqrt(float(value))
    while q(candidate) * q(candidate) > value:
        candidate = float(np.nextafter(candidate, -math.inf))
    return q(candidate)


def reserve_audit_cap_upper(
    maximum_weight_upper: Fraction,
    second_weight_upper: Fraction,
) -> tuple[Fraction, Fraction]:
    """Closed prior-reserve audit cap, rounded in the weak direction.

    The exact minimizing projection vertex is ``(1,2r-1,-1,-1)`` with
    ``r=(1-w0)/w1``.  The minimizing spectral bias is
    ``1/(1+sqrt(2r))``.  Lower-bounding that square root lower-bounds the
    total reserve and therefore upper-bounds the admissible audit.
    """

    ratio_lower = 2 * (1 - maximum_weight_upper) / second_weight_upper
    root_lower = sqrt_fraction_down(max(Fraction(0), ratio_lower))
    unit = 1 + root_lower
    reserve_lower = (5 - 2 / unit - 1 / (unit * unit)) / 2
    if reserve_lower <= 0:
        raise ArithmeticError("prior-reserve lower bound must be positive")
    return 1 / reserve_lower, root_lower


def complement_reserve_lower(
    excluded_weight_upper: Fraction,
    largest_other_weight_upper: Fraction,
) -> tuple[Fraction, Fraction, Fraction]:
    """Lower-bound the prior mass outside one readout effect.

    With ``u_i=(1+x_i)/2``, longitudinal closure is
    ``sum_i w_i u_i=1``.  The reserve is increasing and concave in ``u``.
    Giving the excluded effect zero cost fills that coordinate first, then
    concentrates the residual in the largest other effect.  Write the
    resulting fractional demand as ``r``.  Minimising the remaining scalar
    bias gives, for ``0 <= r < 1/2``,

        sum_{j != i} p_j / A >=
          2 [r(4r-3) + sqrt(r(4r^2-7r+3))] / (1-2r)^2.

    The increasing envelope equals one for ``r >= 1/2``.  A box lower bound
    on ``r`` and a downward dyadic square root therefore give a rigorous
    lower coefficient even if the true fractional demand crosses the branch.
    """

    ratio_lower = (1 - excluded_weight_upper) / largest_other_weight_upper
    if ratio_lower <= 0:
        return Fraction(0), Fraction(0), max(Fraction(0), ratio_lower)
    if ratio_lower >= Fraction(1, 2):
        return Fraction(1), Fraction(0), ratio_lower
    radicand = ratio_lower * (
        4 * ratio_lower * ratio_lower - 7 * ratio_lower + 3
    )
    root_lower = sqrt_fraction_down(radicand)
    denominator = (1 - 2 * ratio_lower) ** 2
    reserve_lower = 2 * (
        ratio_lower * (4 * ratio_lower - 3) + root_lower
    ) / denominator
    if reserve_lower < 0:
        raise ArithmeticError("complement reserve lower bound became negative")
    return reserve_lower, root_lower, ratio_lower


def pair_reserve_lower(
    hull: tuple[tuple[Fraction, Fraction], ...],
    first: int,
    second: int,
) -> tuple[Fraction, Fraction, Fraction]:
    """Lower-bound the normalised prior reserve on two sorted effects.

    The two zero-cost coordinates are filled first.  Their missing weight is
    ``(w_i+w_j-1)_+`` and is concentrated in the larger cost-bearing effect.
    If ``r`` is that fractional fill, exact minimisation over the common
    Helstrom bias gives

        h(r) = [2 sqrt(2r)(1-r) + r(6r-5)] / (1-2r)^2,

    with the removable value ``h(1/2)=7/8``.  The function is increasing on
    ``[0,1]``; lower box endpoints in the numerator, an upper endpoint in the
    denominator, and a downward square root therefore give a safe cut.
    """

    if not 0 <= first < second < 4:
        raise ValueError("pair indices must satisfy 0 <= first < second < 4")
    residual_lower = max(
        Fraction(0), hull[first][0] + hull[second][0] - 1
    )
    ratio_lower = residual_lower / hull[first][1]
    if not Fraction(0) <= ratio_lower <= Fraction(1):
        raise ArithmeticError("pair fractional demand escaped [0,1]")
    if ratio_lower == Fraction(1, 2):
        return Fraction(7, 8), Fraction(1), ratio_lower
    root_lower = sqrt_fraction_down(2 * ratio_lower)
    denominator = (1 - 2 * ratio_lower) ** 2
    reserve_lower = (
        2 * root_lower * (1 - ratio_lower)
        + ratio_lower * (6 * ratio_lower - 5)
    ) / denominator
    # Downward square-root rounding can be amplified near the removable
    # singularity.  Zero remains a valid lower bound in that extreme case.
    reserve_lower = max(Fraction(0), reserve_lower)
    return reserve_lower, root_lower, ratio_lower


def reserve_value_exact(t_value: Fraction, x_value: Fraction) -> Fraction:
    """Evaluate ``f_t(x)`` exactly, including its removable corner."""

    denominator = 1 + t_value * x_value
    if denominator == 0:
        if t_value == 1 and x_value == -1:
            return Fraction(0)
        raise ArithmeticError("unexpected zero reserve denominator")
    return (
        1 + 2 * t_value * x_value + t_value * t_value
    ) / (2 * denominator)


def affine_term_lower(
    coefficient: Fraction, lower: Fraction, upper: Fraction
) -> Fraction:
    return coefficient * (lower if coefficient >= 0 else upper)


def scalar_reserve_minimum_lower(
    x_value: Fraction,
    constant_term: Fraction,
    linear_term: Fraction,
    free_coefficient: Fraction,
) -> Fraction:
    """Exactly enclose the minimum of ``a+b*t+c*f_t(x)`` on ``[0,1]``.

    Its derivative has a quadratic numerator.  Exact rational monotonicity
    intervals isolate every stationary point; 96 dyadic bisections then make
    direct rational interval evaluation effectively sharp while remaining a
    one-sided lower bound.
    """

    if not -1 <= x_value <= 1 or free_coefficient <= 0:
        raise ValueError("invalid scalar reserve parameters")

    def objective(t_value: Fraction) -> Fraction:
        return (
            constant_term
            + linear_term * t_value
            + free_coefficient * reserve_value_exact(t_value, x_value)
        )

    quadratic = 2 * linear_term * x_value * x_value + free_coefficient * x_value
    linear = 4 * linear_term * x_value + 2 * free_coefficient
    constant = 2 * linear_term + free_coefficient * x_value

    def derivative_polynomial(t_value: Fraction) -> Fraction:
        return quadratic * t_value * t_value + linear * t_value + constant

    knots = [Fraction(0), Fraction(1)]
    if quadratic != 0:
        turning = -linear / (2 * quadratic)
        if 0 < turning < 1:
            knots.append(turning)
    knots = sorted(set(knots))
    root_intervals: list[tuple[Fraction, Fraction]] = []
    for lower, upper in zip(knots, knots[1:]):
        lower_value = derivative_polynomial(lower)
        upper_value = derivative_polynomial(upper)
        if lower_value == 0:
            root_intervals.append((lower, lower))
        if upper_value == 0:
            root_intervals.append((upper, upper))
        if lower_value * upper_value >= 0:
            continue
        left, right = lower, upper
        left_value = lower_value
        for _ in range(96):
            middle = (left + right) / 2
            middle_value = derivative_polynomial(middle)
            if middle_value == 0:
                left = right = middle
                break
            if left_value * middle_value <= 0:
                right = middle
            else:
                left = middle
                left_value = middle_value
        root_intervals.append((left, right))

    candidates = [objective(Fraction(0)), objective(Fraction(1))]
    q0 = 2 * constant_term + free_coefficient
    q1 = 2 * (
        constant_term * x_value
        + linear_term
        + free_coefficient * x_value
    )
    q2 = 2 * linear_term * x_value + free_coefficient
    for lower, upper in root_intervals:
        if lower == upper:
            candidates.append(objective(lower))
            continue
        numerator_lower = (
            q0
            + affine_term_lower(q1, lower, upper)
            + affine_term_lower(q2, lower * lower, upper * upper)
        )
        denominator_upper = 1 + x_value * (
            upper if x_value >= 0 else lower
        )
        if denominator_upper <= 0:
            candidates.append(Fraction(0))
        else:
            candidates.append(
                max(Fraction(0), numerator_lower)
                / (2 * denominator_upper)
            )
    return min(candidates)


@cache
def weighted_reserve_cut_lower(
    hull: tuple[tuple[Fraction, Fraction], ...],
    coefficients: tuple[Fraction, ...],
    weight_vertices: tuple[tuple[Fraction, ...], ...] | None = None,
) -> Fraction:
    """Lower-bound ``min sum_i c_i f_t(x_i)`` on a weight box.

    For fixed ``t`` the objective is concave in the four longitudinal
    coordinates, so a minimum occurs at a vertex of the closure slice: three
    coordinates are signs and one is free.  The free reserve is increasing
    in its coordinate.  Independent exact hull arithmetic therefore supplies
    a safe lower coordinate, after which the scalar minimiser above is exact.
    """

    if len(hull) != 4 or len(coefficients) != 4:
        raise ValueError("expected four weights and four coefficients")
    if weight_vertices is not None and not weight_vertices:
        raise ValueError("an exact weight-vertex list cannot be empty")
    if weight_vertices is not None and any(len(row) != 4 for row in weight_vertices):
        raise ValueError("each exact weight vertex must have four coordinates")
    candidates: list[Fraction] = []
    for free in range(4):
        fixed = [index for index in range(4) if index != free]
        for signs in itertools.product((-1, 1), repeat=3):
            if weight_vertices is None:
                sum_lower = Fraction(0)
                sum_upper = Fraction(0)
                for index, sign in zip(fixed, signs):
                    lower, upper = hull[index]
                    if sign > 0:
                        sum_lower += lower
                        sum_upper += upper
                    else:
                        sum_lower -= upper
                        sum_upper -= lower
                numerator = (-sum_upper, -sum_lower)
                denominator = hull[free]
                quotients = [
                    item / divisor
                    for item in numerator
                    for divisor in denominator
                ]
            else:
                # A linear-fractional function with positive denominator
                # reaches both extrema on vertices of the exact weight
                # polytope.  This retains sum(w)=2 and all sorting/box
                # correlations instead of multiplying independent intervals.
                quotients = [
                    -sum(
                        weights[index] * sign
                        for index, sign in zip(fixed, signs)
                    )
                    / weights[free]
                    for weights in weight_vertices
                ]
            x_lower = min(quotients)
            x_upper = max(quotients)
            if x_upper < -1 or x_lower > 1:
                continue
            x_lower = max(Fraction(-1), x_lower)
            plus = sum(
                coefficients[index]
                for index, sign in zip(fixed, signs)
                if sign > 0
            )
            minus = sum(
                coefficients[index]
                for index, sign in zip(fixed, signs)
                if sign < 0
            )
            candidates.append(
                scalar_reserve_minimum_lower(
                    x_lower,
                    (plus + minus) / 2,
                    (plus - minus) / 2,
                    coefficients[free],
                )
            )
    if not candidates:
        raise ArithmeticError("weight hull produced no closure vertex")
    return min(candidates)


def pair_loss_factor_upper(
    hull: tuple[tuple[Fraction, Fraction], ...],
    first: int,
    second: int,
) -> Fraction:
    """Upper-bound the direct projective loss for any ordered effect pair.

    Full Bloch closure bounds the resultant of a selected pair by the total
    weight of its complementary pair.  The overlap factor decreases in both
    selected weights, so exact lower hull endpoints give a box upper.  The
    universal overlap cap one covers pairs whose lower weight sum has not yet
    crossed one.
    """

    if first == second or first not in range(4) or second not in range(4):
        raise ValueError("effect indices must be distinct elements of range(4)")
    first_weight = hull[first][0]
    second_weight = hull[second][0]
    factor = (1 - first_weight) * (1 - second_weight) / (
        first_weight * second_weight
    )
    return min(Fraction(1), factor)


def dominant_pair_loss_factor_upper(
    hull: tuple[tuple[Fraction, Fraction], ...],
) -> Fraction:
    """Upper-bound the direct projective loss for the two largest effects.

    Full Bloch-vector closure gives

        ||w0*n0 + w1*n1|| <= w2+w3 = 2-w0-w1.

    Hence ``(1+n0.n1)/2`` is at most
    ``(1-w0)(1-w1)/(w0*w1)``.  This expression decreases in both weights,
    so the two exact lower hull endpoints give a rigorous box upper bound.
    """

    return pair_loss_factor_upper(hull, 0, 1)


@cache
def physical_weight_vertices(
    box: WeightBox,
) -> tuple[tuple[Fraction, Fraction, Fraction, Fraction], ...]:
    """Return every exact physical four-weight vertex of ``box``."""

    return tuple(
        (point[0], Fraction(2) - sum(point), point[1], point[2])
        for point in polytope_vertices(box)
    )


def safe_linear_row(
    size: int,
    coefficients: dict[int, Fraction],
    right: Fraction,
) -> tuple[np.ndarray, float, Fraction]:
    """Encode ``a.x <= b`` weakly for variables known to lie in [0,1]."""

    row = np.zeros(size, dtype=float)
    correction = Fraction(0)
    for index, exact in coefficients.items():
        encoded = float(exact)
        row[index] = encoded
        correction += max(Fraction(0), q(encoded) - exact)
    safe_right = binary64_up(right + correction)
    if q(safe_right) < right + correction:
        raise ArithmeticError("failed to relax a linear-row right side")
    return row, safe_right, correction


def mccormick_rows(
    size: int,
    x: int,
    y: int,
    product: int,
    x_bounds: tuple[Fraction, Fraction],
    y_bounds: tuple[Fraction, Fraction],
) -> list[tuple[np.ndarray, float, Fraction]]:
    lx, ux = x_bounds
    ly, uy = y_bounds
    return [
        safe_linear_row(
            size,
            {product: -1, y: lx, x: ly},
            lx * ly,
        ),
        safe_linear_row(
            size,
            {product: -1, y: ux, x: uy},
            ux * uy,
        ),
        safe_linear_row(
            size,
            {product: 1, y: -ux, x: -ly},
            -ux * ly,
        ),
        safe_linear_row(
            size,
            {product: 1, y: -lx, x: -uy},
            -lx * uy,
        ),
    ]


class McCormickOracle:
    """Persistent conic model with a parameterised exact outer polytope."""

    LINEAR_SIZE = 32
    PARAMETER_ROW_COUNT = (
        114
        + 12 * len(PERTURBATIONS)
        + len(COMMON_BIAS_COEFFICIENTS)
    )

    def __init__(self) -> None:
        self.path = cp.Variable((4, 4), nonneg=True, name="path")
        self.audit = cp.Variable(nonneg=True, name="audit")
        self.returned = cp.Variable(nonneg=True, name="return")
        self.cross = cp.Variable(120, nonneg=True, name="hellinger_cross")
        self.weight = cp.Variable(4, nonneg=True, name="weight")
        self.loss_factor = cp.Variable(4, nonneg=True, name="loss_factor")
        self.weight_loss = cp.Variable(4, nonneg=True, name="weight_loss")
        self.audit_gap = cp.Variable(4, nonneg=True, name="audit_gap")
        self.projective_loss = cp.Variable(4, nonneg=True, name="projective_loss")
        self.syndrome_product = cp.Variable(4, nonneg=True, name="syndrome_product")
        self.prefix_gap = cp.Variable(nonneg=True, name="prefix_gap")
        self.prefix_product = cp.Variable(nonneg=True, name="prefix_product")
        self.box_matrix = cp.Parameter(
            (self.PARAMETER_ROW_COUNT, self.LINEAR_SIZE)
        )
        self.box_right = cp.Parameter(self.PARAMETER_ROW_COUNT)

        flat = cp.reshape(self.path, (16,), order="C")
        prefix = cp.sum(self.path, axis=1)
        self.syndrome = cp.hstack(
            [sum(self.path[z, z ^ s] for z in range(4)) for s in range(4)]
        )
        # Coordinate inventory used by every safe parameter row.
        pieces = (
            self.weight,
            self.loss_factor,
            self.weight_loss,
            self.audit_gap,
            self.projective_loss,
            self.syndrome_product,
            cp.hstack(
                [
                    self.prefix_gap,
                    self.prefix_product,
                    self.audit,
                    self.returned,
                ]
            ),
            self.syndrome,
        )
        self.linear = cp.hstack(pieces)
        if self.linear.size != self.LINEAR_SIZE:
            raise AssertionError("wrong safe-linear coordinate inventory")

        constraints: list[cp.Constraint] = [
            cp.sum(flat) == 1,
            cp.sum(self.weight) == 2,
            self.audit <= 1,
            self.returned <= 1,
            self.syndrome <= self.audit,
            self.audit_gap == self.audit - self.syndrome,
            2 * self.loss_factor - self.weight_loss == 1 - self.weight,
            self.audit <= cp.sum(self.syndrome_product),
            *(prefix[index] >= prefix[index + 1] for index in range(3)),
            self.prefix_gap == prefix[0] + prefix[1] - 2 * prefix[2],
            self.audit <= 2 * prefix[2] + self.prefix_product,
            self.box_matrix @ self.linear <= self.box_right,
        ]
        cursor = 0
        for first in range(16):
            for second in range(first + 1, 16):
                constraints.append(
                    cp.SOC(
                        flat[first] + flat[second],
                        cp.hstack(
                            [
                                2 * self.cross[cursor],
                                flat[first] - flat[second],
                            ]
                        ),
                    )
                )
                cursor += 1
        constraints.append(
            16 * self.returned <= cp.sum(flat) + 2 * cp.sum(self.cross)
        )
        objective = cp.Maximize(
            binary64_up(SUPPORT_WEIGHT) * self.audit
            + binary64_up(1 - SUPPORT_WEIGHT) * self.returned
        )
        base = cp.Problem(objective, constraints)
        variables = base.variables()
        if not all(bool(variable.attributes.get("nonneg")) for variable in variables):
            raise RuntimeError("every residual-controlled variable must be nonnegative")
        self.problem = cp.Problem(
            objective,
            [*constraints, *(variable <= 1 for variable in variables)],
        )
        if not self.problem.is_dpp():
            raise RuntimeError("McCormick four-active SOCP must be DPP")

    @staticmethod
    def indices() -> dict[str, tuple[int, ...] | int]:
        return {
            "weight": tuple(range(0, 4)),
            "factor": tuple(range(4, 8)),
            "weight_factor": tuple(range(8, 12)),
            "gap": tuple(range(12, 16)),
            "loss": tuple(range(16, 20)),
            "syndrome_product": tuple(range(20, 24)),
            "prefix_gap": 24,
            "prefix_product": 25,
            "audit": 26,
            "return": 27,
            "syndrome": tuple(range(28, 32)),
        }

    def assign(
        self,
        hull: tuple[tuple[Fraction, Fraction], ...],
        permutation: tuple[int, int, int],
        weight_vertices: tuple[tuple[Fraction, ...], ...] | None = None,
    ) -> dict[str, Any]:
        order = (0, *permutation)
        physical = tuple(hull[index] for index in order)
        # Enlarge each exact hull to dyadic binary64 endpoints.  Every later
        # McCormick facet is exact for this enlarged rectangle.
        weight_bounds = tuple(
            (
                q(binary64_down(lower)),
                q(binary64_up(upper)),
            )
            for lower, upper in physical
        )
        factor_bounds = tuple(
            (
                q(binary64_down(loss_factor(upper))),
                q(binary64_up(loss_factor(lower))),
            )
            for lower, upper in weight_bounds
        )
        idx = self.indices()
        rows: list[tuple[np.ndarray, float, Fraction]] = []
        position = {base: syndrome for syndrome, base in enumerate(order)}
        for larger, smaller in zip(range(3), range(1, 4)):
            rows.append(
                safe_linear_row(
                    self.LINEAR_SIZE,
                    {
                        idx["weight"][position[smaller]]: 1,  # type: ignore[index]
                        idx["weight"][position[larger]]: -1,  # type: ignore[index]
                    },
                    Fraction(0),
                )
            )

        audit_cap, reserve_root_lower = reserve_audit_cap_upper(
            hull[0][1], hull[1][1]
        )
        rows.append(
            safe_linear_row(
                self.LINEAR_SIZE, {int(idx["audit"]): 1}, audit_cap
            )
        )
        complement_reserves = []
        for excluded in range(4):
            largest_other = 1 if excluded == 0 else 0
            reserve_lower, root_lower, ratio_lower = complement_reserve_lower(
                hull[excluded][1], hull[largest_other][1]
            )
            excluded_syndrome = position[excluded]
            coefficients = {int(idx["audit"]): reserve_lower}
            coefficients.update(
                {
                    idx["syndrome"][syndrome]: -1  # type: ignore[index]
                    for syndrome in range(4)
                    if syndrome != excluded_syndrome
                }
            )
            rows.append(
                safe_linear_row(
                    self.LINEAR_SIZE, coefficients, Fraction(0)
                )
            )
            complement_reserves.append(
                {
                    "excluded_effect": excluded,
                    "excluded_syndrome": excluded_syndrome,
                    "reserve_lower": fraction_pair(reserve_lower),
                    "ratio_lower": fraction_pair(ratio_lower),
                    "sqrt_lower": fraction_pair(root_lower),
                }
            )
        pair_reserves = []
        for first in range(4):
            for second in range(first + 1, 4):
                reserve_lower, root_lower, ratio_lower = pair_reserve_lower(
                    hull, first, second
                )
                first_syndrome = position[first]
                second_syndrome = position[second]
                rows.append(
                    safe_linear_row(
                        self.LINEAR_SIZE,
                        {
                            int(idx["audit"]): reserve_lower,
                            idx["syndrome"][first_syndrome]: -1,  # type: ignore[index]
                            idx["syndrome"][second_syndrome]: -1,  # type: ignore[index]
                        },
                        Fraction(0),
                    )
                )
                pair_reserves.append(
                    {
                        "effects": [first, second],
                        "syndromes": [first_syndrome, second_syndrome],
                        "reserve_lower": fraction_pair(reserve_lower),
                        "ratio_lower": fraction_pair(ratio_lower),
                        "sqrt_lower": fraction_pair(root_lower),
                    }
                )
        weighted_reserves = []
        for perturbation in PERTURBATIONS:
            for elevated, discounted in PERTURBATION_PAIRS:
                coefficients = tuple(
                    1 + perturbation
                    if effect == elevated
                    else 1 - perturbation
                    if effect == discounted
                    else Fraction(1)
                    for effect in range(4)
                )
                reserve_lower = weighted_reserve_cut_lower(
                    hull, coefficients, weight_vertices
                )
                row_coefficients = {int(idx["audit"]): reserve_lower}
                row_coefficients.update(
                    {
                        idx["syndrome"][position[effect]]: -coefficient  # type: ignore[index]
                        for effect, coefficient in enumerate(coefficients)
                    }
                )
                rows.append(
                    safe_linear_row(
                        self.LINEAR_SIZE, row_coefficients, Fraction(0)
                    )
                )
                weighted_reserves.append(
                    {
                        "perturbation": fraction_pair(perturbation),
                        "elevated_effect": elevated,
                        "discounted_effect": discounted,
                        "coefficients": list(map(fraction_pair, coefficients)),
                        "reserve_lower": fraction_pair(reserve_lower),
                        "scalar_minimisation": (
                            "exact derivative isolation and rational intervals"
                        ),
                    }
                )
        common_bias_reserves = []
        for coefficients in COMMON_BIAS_COEFFICIENTS:
            reserve_lower = weighted_reserve_cut_lower(
                hull, coefficients, weight_vertices
            )
            row_coefficients = {int(idx["audit"]): reserve_lower}
            row_coefficients.update(
                {
                    idx["syndrome"][position[effect]]: -coefficient  # type: ignore[index]
                    for effect, coefficient in enumerate(coefficients)
                }
            )
            rows.append(
                safe_linear_row(
                    self.LINEAR_SIZE, row_coefficients, Fraction(0)
                )
            )
            common_bias_reserves.append(
                {
                    "coefficients": list(map(fraction_pair, coefficients)),
                    "reserve_lower": fraction_pair(reserve_lower),
                    "derivation": (
                        "single common Helstrom bias, longitudinal POVM "
                        "closure, and exact scalar support minimisation"
                    ),
                }
            )
        def bounds(variable: int, lower: Fraction, upper: Fraction) -> None:
            rows.append(
                safe_linear_row(self.LINEAR_SIZE, {variable: -1}, -lower)
            )
            rows.append(
                safe_linear_row(self.LINEAR_SIZE, {variable: 1}, upper)
            )

        for syndrome in range(4):
            weight_index = idx["weight"][syndrome]  # type: ignore[index]
            factor_index = idx["factor"][syndrome]  # type: ignore[index]
            weight_factor = idx["weight_factor"][syndrome]  # type: ignore[index]
            gap = idx["gap"][syndrome]  # type: ignore[index]
            loss = idx["loss"][syndrome]  # type: ignore[index]
            syndrome_product = idx["syndrome_product"][syndrome]  # type: ignore[index]
            syndrome_prior = idx["syndrome"][syndrome]  # type: ignore[index]
            bounds(weight_index, *weight_bounds[syndrome])
            bounds(factor_index, *factor_bounds[syndrome])
            rows.extend(
                mccormick_rows(
                    self.LINEAR_SIZE,
                    weight_index,
                    factor_index,
                    weight_factor,
                    weight_bounds[syndrome],
                    factor_bounds[syndrome],
                )
            )
            rows.extend(
                mccormick_rows(
                    self.LINEAR_SIZE,
                    factor_index,
                    gap,
                    loss,
                    factor_bounds[syndrome],
                    (Fraction(0), Fraction(1)),
                )
            )
            rows.extend(
                mccormick_rows(
                    self.LINEAR_SIZE,
                    weight_index,
                    syndrome_prior,
                    syndrome_product,
                    weight_bounds[syndrome],
                    (Fraction(0), Fraction(1)),
                )
            )
        rows.extend(
            mccormick_rows(
                self.LINEAR_SIZE,
                idx["weight"][0],  # type: ignore[index]
                int(idx["prefix_gap"]),
                int(idx["prefix_product"]),
                weight_bounds[0],
                (Fraction(0), Fraction(1)),
            )
        )
        for line_weight, line_upper in PROJECTIVE_LINES:
            for syndrome in range(4):
                rows.append(
                    safe_linear_row(
                        self.LINEAR_SIZE,
                        {
                            int(idx["audit"]): line_weight,
                            idx["loss"][syndrome]: -line_weight,  # type: ignore[index]
                            int(idx["return"]): 1 - line_weight,
                        },
                        line_upper,
                    )
                )
        pair_factors = {
            pair: pair_loss_factor_upper(hull, *pair)
            for pair in ORDERED_EFFECT_PAIRS
        }
        full_closure_pair_cuts = []
        for line_weight, line_upper in PROJECTIVE_LINES:
            for retained, other in ORDERED_EFFECT_PAIRS:
                pair_factor = pair_factors[(retained, other)]
                rows.append(
                    safe_linear_row(
                        self.LINEAR_SIZE,
                        {
                            int(idx["audit"]): line_weight * (1 - pair_factor),
                            idx["syndrome"][position[other]]: line_weight * pair_factor,  # type: ignore[index]
                            int(idx["return"]): 1 - line_weight,
                        },
                        line_upper,
                    )
                )
                full_closure_pair_cuts.append(
                    {
                        "line_weight": fraction_pair(line_weight),
                        "line_upper": fraction_pair(line_upper),
                        "retained_effect": retained,
                        "other_effect": other,
                        "other_syndrome": position[other],
                        "loss_factor_upper": fraction_pair(pair_factor),
                        "derivation": (
                            "full Bloch closure and the complementary-pair triangle bound"
                        ),
                    }
                )
        if len(rows) != self.PARAMETER_ROW_COUNT:
            raise AssertionError(
                f"expected {self.PARAMETER_ROW_COUNT} parameter rows, got {len(rows)}"
            )
        self.box_matrix.value = np.vstack([row for row, _, _ in rows])
        self.box_right.value = np.asarray([right for _, right, _ in rows])
        corrections = [correction for _, _, correction in rows]
        return {
            "syndrome_order": list(order),
            "weight_bounds": [
                [fraction_pair(lower), fraction_pair(upper)]
                for lower, upper in weight_bounds
            ],
            "loss_factor_bounds": [
                [fraction_pair(lower), fraction_pair(upper)]
                for lower, upper in factor_bounds
            ],
            "maximum_linear_rounding_correction": fraction_pair(max(corrections)),
            "prior_reserve_audit_cap": fraction_pair(audit_cap),
            "prior_reserve_sqrt_lower": fraction_pair(reserve_root_lower),
            "complement_reserve_cuts": complement_reserves,
            "pair_reserve_cuts": pair_reserves,
            "weighted_reserve_cuts": weighted_reserves,
            "common_bias_support_cuts": common_bias_reserves,
            "exact_weight_vertex_count": (
                len(weight_vertices) if weight_vertices is not None else None
            ),
            "full_closure_projective_pair_cuts": full_closure_pair_cuts,
        }

    def canonical_data(self) -> dict[str, Any]:
        data, _, _ = self.problem.get_problem_data(cp.CLARABEL)
        variable_count = sum(variable.size for variable in self.problem.variables())
        if data["A"].shape[1] != variable_count:
            raise RuntimeError("canonicalisation introduced an unbounded variable")
        return data


def exact_upper(
    data: dict[str, Any], dual: np.ndarray
) -> tuple[Fraction, Fraction, Fraction]:
    residuals, correction = exact_sparse_stationarity(
        data["A"], dual, data["c"]
    )
    upper = exact_dot(np.asarray(data["b"]), dual) + correction
    return upper, correction, max(map(abs, residuals), default=Fraction(0))


class ExactCertifier:
    def __init__(self, target: Fraction) -> None:
        self.oracle = McCormickOracle()
        self.solver = CLARABEL()
        self.target = target

    def certify(
        self,
        box: WeightBox,
        hull: tuple[tuple[Fraction, Fraction], ...],
        permutation: tuple[int, int, int],
    ) -> dict[str, Any]:
        enclosure = self.oracle.assign(hull, permutation, physical_weight_vertices(box))
        data = self.oracle.canonical_data()
        result = self.solver.solve_via_data(
            data,
            warm_start=False,
            verbose=False,
            solver_opts={
                "tol_gap_abs": 1e-11,
                "tol_gap_rel": 1e-11,
                "tol_feas": 1e-11,
                "max_iter": 500,
            },
            solver_cache=None,
        )
        repaired, repaired_blocks = repair_dual_cones(
            np.asarray(result.z), data["dims"]
        )
        raw_upper, _, _ = exact_upper(data, repaired)
        storage_dtype = "f32"
        encoded_dual = encode_dual(repaired, storage_dtype)
        stored = decode_dual(encoded_dual, storage_dtype)
        stored, _ = repair_dual_cones(stored, data["dims"])
        upper, correction, residual = exact_upper(data, stored)
        if upper > self.target and raw_upper <= self.target:
            storage_dtype = "f64"
            encoded_dual = encode_dual(repaired, storage_dtype)
            stored = decode_dual(encoded_dual, storage_dtype)
            stored, _ = repair_dual_cones(stored, data["dims"])
            upper, correction, residual = exact_upper(data, stored)
        closed = upper <= self.target
        report: dict[str, Any] = {
            "syndrome_permutation": list(permutation),
            "coefficient_enclosure": enclosure,
            "canonical_shape": [int(data["A"].shape[0]), int(data["A"].shape[1])],
            "canonical_nonzeros": int(data["A"].nnz),
            "canonical_sha256": canonical_hash(data),
            "cone_dimensions": {
                "zero": int(data["dims"].zero),
                "nonnegative": int(data["dims"].nonneg),
                "soc": list(map(int, data["dims"].soc)),
            },
            "untrusted_solver_status": str(result.status),
            "untrusted_primal_objective": float(result.obj_val),
            "untrusted_dual_objective": float(result.obj_val_dual),
            "soc_heads_repaired": repaired_blocks,
            "certified_upper_fraction": fraction_pair(upper),
            "certified_upper_decimal": fraction_decimal(upper),
            "exact_residual_correction": fraction_pair(correction),
            "maximum_stationarity_residual_decimal": fraction_decimal(residual),
            "closed": closed,
            "trusted_optimizers": [],
            "untrusted_search_helpers": ["Clarabel dual-vector proposal"],
        }
        if closed:
            report.update(
                {
                    "dual_storage_dtype": storage_dtype,
                    "dual_zlib_base64": encoded_dual,
                }
            )
        return report


_WORKER_CERTIFIER: ExactCertifier | None = None


def initialise_certifier_worker(target: Fraction) -> None:
    global _WORKER_CERTIFIER
    _WORKER_CERTIFIER = ExactCertifier(target)


def certify_order_worker(
    payload: tuple[
        WeightBox,
        tuple[tuple[Fraction, Fraction], ...],
        tuple[int, int, int],
    ],
) -> dict[str, Any]:
    if _WORKER_CERTIFIER is None:
        raise RuntimeError("four-active worker was not initialised")
    box, hull, permutation = payload
    return _WORKER_CERTIFIER.certify(box, hull, permutation)


def assess_box(
    certifier: ExactCertifier,
    box: WeightBox,
    executor: ProcessPoolExecutor | None = None,
) -> dict[str, Any]:
    hull = weight_hull(box)
    if hull is None:
        return {"kind": "domain-empty", "box": box.serialise()}
    if executor is None:
        reports = [
            certifier.certify(box, hull, permutation)
            for permutation in NONZERO_PERMUTATIONS
        ]
    else:
        reports = list(
            executor.map(
                certify_order_worker,
                [(box, hull, permutation) for permutation in NONZERO_PERMUTATIONS],
            )
        )
    maximum = max(
        Fraction(*report["certified_upper_fraction"]) for report in reports
    )
    return {
        "kind": "closed" if all(report["closed"] for report in reports) else "open",
        "box": box.serialise(),
        "exact_weight_hull": [
            [fraction_pair(lower), fraction_pair(upper)]
            for lower, upper in hull
        ],
        "maximum_certified_upper_fraction": fraction_pair(maximum),
        "maximum_certified_upper_decimal": fraction_decimal(maximum),
        "order_certificates": reports,
    }


def write_payload(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", default=str(FOUR_ACTIVE_TARGET))
    parser.add_argument("--max-splits", type=int, default=5000)
    parser.add_argument("--checkpoint-every", type=int, default=10)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument(
        "--resume", action="store_true"
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "four_active_mccormick_socp_exact_cover_l060.json",
    )
    args = parser.parse_args()
    target = Fraction(args.target)
    certifier = ExactCertifier(target)
    if args.workers < 1:
        parser.error("--workers must be positive")
    executor = (
        ProcessPoolExecutor(
            max_workers=args.workers,
            initializer=initialise_certifier_worker,
            initargs=(target,),
        )
        if args.workers > 1 else None
    )
    root = initial_box()
    leaves: list[dict[str, Any]] = []
    queue: list[tuple[float, int, WeightBox]] = []
    counter = 0
    splits = 0
    if args.resume and args.output.exists():
        prior = json.loads(args.output.read_text(encoding="utf-8"))
        if prior.get("schema") != SCHEMA:
            parser.error("resume artifact has an incompatible schema")
        if prior.get("target") != fraction_pair(target):
            parser.error("resume artifact has a different target")
        if prior.get("initial_box") != root.serialise():
            parser.error("resume artifact has a different root box")
        leaves = list(prior.get("leaves", []))
        splits = int(prior.get("boxes_split", 0))
        open_payloads = prior.get("open_boxes")
        if open_payloads is None:
            parser.error("resume artifact does not store its open frontier")
        for serialised in open_payloads:
            box = WeightBox.deserialise(serialised)
            report = assess_box(certifier, box, executor)
            if report["kind"] == "open":
                maximum = Fraction(
                    *report["maximum_certified_upper_fraction"]
                )
                heapq.heappush(
                    queue, (-float(maximum), counter, box)
                )
                counter += 1
            else:
                leaves.append(report)
    else:
        first = assess_box(certifier, root, executor)
        if first["kind"] == "open":
            maximum = Fraction(*first["maximum_certified_upper_fraction"])
            heapq.heappush(queue, (-float(maximum), counter, root))
            counter += 1
        else:
            leaves.append(first)

    def payload() -> dict[str, Any]:
        complete = not queue
        closed = [leaf for leaf in leaves if leaf["kind"] == "closed"]
        maximum = max(
            (
                Fraction(*leaf["maximum_certified_upper_fraction"])
                for leaf in closed
            ),
            default=Fraction(0),
        )
        result = {
            "schema": SCHEMA,
            "support_weight": fraction_pair(SUPPORT_WEIGHT),
            "target": fraction_pair(target),
            "maximum_weight_floor": fraction_pair(MAXIMUM_WEIGHT_FLOOR),
            "minimum_active_weight": fraction_pair(MINIMUM_ACTIVE_WEIGHT),
            "projective_lines": [
                [fraction_pair(weight), fraction_pair(upper)]
                for weight, upper in PROJECTIVE_LINES
            ],
            "reserve_perturbations": list(map(fraction_pair, PERTURBATIONS)),
            "common_bias_coefficient_representatives": [
                list(map(fraction_pair, representative))
                for representative in COMMON_BIAS_COEFFICIENT_REPRESENTATIVES
            ],
            "common_bias_coefficient_orbit_size": len(COMMON_BIAS_COEFFICIENTS),
            "weighted_reserve_geometry": "exact correlated weight-polytope vertices",
            "projective_pair_geometry": (
                "all ordered effect-pair projective comparisons from full Bloch closure"
            ),
            "prefix_order_reduction": (
                "AGL(2,2)=S4: translations canonicalise the prefix order; "
                "the six GL(2,2) parts permute the nonzero syndrome labels"
            ),
            "relaxation": (
                "shared sorted weights; four-facet McCormick envelopes for "
                "w*p, w0*prefix_gap, k*w, and k*(A-p), plus exact total, "
                "four complementary, six pair, and twenty-four multiscale "
                "common-bias Helstrom prior-reserve cuts, forty-eight rational "
                "all-prior common-bias support cuts in three label orbits, plus "
                "twenty-four full-closure projective pair comparisons"
            ),
            "initial_box": root.serialise(),
            "boxes_split": splits,
            "leaf_count": len(leaves),
            "closed_leaf_count": len(closed),
            "domain_empty_leaf_count": sum(
                leaf["kind"] == "domain-empty" for leaf in leaves
            ),
            "boxes_remaining": len(queue),
            "maximum_open_upper": -queue[0][0] if queue else None,
            "leading_open_box": queue[0][2].serialise() if queue else None,
            "open_boxes": [
                item[2].serialise()
                for item in sorted(queue, key=lambda item: (item[0], item[1]))
            ],
            "maximum_certified_upper_decimal": fraction_decimal(maximum),
            "complete": complete,
            "all_cells_closed": complete,
            "leaves": sorted(leaves, key=lambda leaf: leaf["box"]["path"]),
            "trusted_optimizers": [],
            "untrusted_search_helpers": [
                "Clarabel dual-vector proposals",
                "best-first box ordering",
            ],
        }
        if complete:
            validate_leaf_tree(result["leaves"])
        return result

    while queue and splits < args.max_splits:
        _, _, box = heapq.heappop(queue)
        coordinate = choose_split(box)
        for child in split_box(box, coordinate):
            report = assess_box(certifier, child, executor)
            if report["kind"] == "open":
                counter += 1
                maximum = Fraction(*report["maximum_certified_upper_fraction"])
                heapq.heappush(queue, (-float(maximum), counter, child))
            else:
                leaves.append(report)
        splits += 1
        if splits % args.checkpoint_every == 0:
            print(
                json.dumps(
                    {
                        "boxes_split": splits,
                        "closed_or_empty": len(leaves),
                        "open": len(queue),
                        "maximum_open_upper": -queue[0][0] if queue else None,
                    }
                ),
                flush=True,
            )
            write_payload(args.output, payload())
    result = payload()
    write_payload(args.output, result)
    if executor is not None:
        executor.shutdown()
    print(json.dumps({key: value for key, value in result.items() if key != "leaves"}, indent=2))


if __name__ == "__main__":
    main()

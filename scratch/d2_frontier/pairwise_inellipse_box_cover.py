"""SOCP cover for the clean-POVM inellipse model.

For fixed reciprocal inellipse parameters ``alpha=1/w`` and ``beta=1/t``,
the homogenised Horwitz inequality is one Lorentz-cone constraint in each
conditional probability point.  A parameter box is relaxed by enlarging the
midpoint ellipse with a uniform coefficient-error allowance.  Therefore every
node is an SOCP upper relaxation and the allowance converges to zero with the
box widths.

The necessity of an inellipse follows from the preprocessing order on qubit
POVMs.  Every ternary qubit POVM is below a clean ternary qubit POVM; the
qubit clean-POVM classification makes the latter rank one (or a degenerate
rank-one limit).  Its probability range is an inellipse containing the range
of the original POVM.  Relabel the outcomes so that the smallest barycentric
coordinate of the inellipse centre is residual.  Horwitz's parameters then
obey ``w,t >= 1/2``, so the three residual choices cover the closed model with
``1 <= alpha,beta <= 2``.

Each invocation handles one residual choice per selected pair.  A complete
union certificate must cover ``xy``, ``xr``, and ``yr``.  Solver bounds remain
numerical until their conic duals are outward rounded.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import heapq
import json
import math
from pathlib import Path
from typing import Any

import cvxpy as cp
import numpy as np

from pairwise_qubit_helstrom_scip import Column, parse_pair, render_column
from terminal_weight_upper import filled_effect_weights


OUTCOMES = range(4)
PATHS = tuple((z, y) for z in OUTCOMES for y in OUTCOMES)
IDENTITY = np.eye(2, dtype=complex)
SIGMA_X = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
SIGMA_Y = np.array([[0.0, -1j], [1j, 0.0]], dtype=complex)


def canonical_three_effect_povm(weights: np.ndarray) -> np.ndarray:
    """Return the planar rank-one qubit POVM with the requested traces."""

    w0, w1, w2 = (float(value) for value in weights)
    cosine = (w2 * w2 - w0 * w0 - w1 * w1) / (2.0 * w0 * w1)
    cosine = float(np.clip(cosine, -1.0, 1.0))
    sine = math.sqrt(max(0.0, 1.0 - cosine * cosine))
    vectors = np.asarray(
        [
            [1.0, 0.0],
            [cosine, sine],
            [-(w0 + w1 * cosine) / w2, -(w1 * sine) / w2],
        ]
    )
    effects = np.zeros((4, 2, 2), dtype=complex)
    for index in range(3):
        effects[index] = 0.5 * weights[index] * (
            IDENTITY + vectors[index, 0] * SIGMA_X + vectors[index, 1] * SIGMA_Y
        )
    if np.linalg.norm(effects.sum(axis=0) - IDENTITY) > 2e-10:
        raise RuntimeError("canonical three-effect POVM failed completeness")
    return effects


def reconstruction_matrix(effects: np.ndarray) -> np.ndarray:
    """Map ``[trace, x, y]`` to the three terminal probabilities."""

    matrix = np.asarray(
        [
            [
                0.5 * float(np.trace(effect).real),
                0.5 * float(np.trace(effect @ SIGMA_X).real),
                0.5 * float(np.trace(effect @ SIGMA_Y).real),
            ]
            for effect in effects[:3]
        ],
        dtype=float,
    )
    if abs(float(np.linalg.det(matrix))) < 1e-10:
        raise ValueError("the terminal three-effect POVM must be nondegenerate")
    return np.linalg.inv(matrix)


def hellinger_hypograph(
    probabilities: list[cp.Expression], constraints: list[cp.Constraint]
) -> cp.Expression:
    """Conic hypograph of ``(sum_i sqrt(p_i))**2 / 16``."""

    cross_terms: list[cp.Expression] = []
    for first in range(len(probabilities)):
        constraints.append(probabilities[first] >= 0.0)
        for second in range(first + 1, len(probabilities)):
            geometric = cp.Variable(nonneg=True)
            constraints.append(
                cp.SOC(
                    probabilities[first] + probabilities[second],
                    cp.hstack(
                        (
                            2.0 * geometric,
                            probabilities[first] - probabilities[second],
                        )
                    ),
                )
            )
            cross_terms.append(geometric)
    return (
        cp.sum(cp.hstack(probabilities))
        + 2.0 * cp.sum(cp.hstack(cross_terms))
    ) / 16.0


Box = dict[str, tuple[float, float]]


def initial_box(pair_count: int, reciprocal_cap: float) -> Box:
    return {
        f"k{pair_index}_{coordinate}": (1.0, reciprocal_cap)
        for pair_index in range(pair_count)
        for coordinate in ("alpha", "beta")
    }


def selected_statistic(
    statistics: cp.Expression,
    probability: cp.Expression,
    weights: np.ndarray,
    column: Column,
    z: int,
) -> cp.Expression:
    kind, y, t = column
    if kind == "b":
        return statistics[z, y, t]
    return weights[t] * probability[z, y] - statistics[z, y, t]


def cross_coefficient(alpha: float, beta: float) -> float:
    return -2.0 * (alpha * beta - 2.0 * alpha - 2.0 * beta + 2.0)


def quadratic_maximum_on_unit_triangle(
    x_squared: float,
    y_squared: float,
    xy: float,
    x_linear: float,
    y_linear: float,
    constant: float = 0.0,
) -> float:
    """Exact maximum of a bivariate quadratic on ``x,y>=0, x+y<=1``."""

    def value(x_value: float, y_value: float) -> float:
        return (
            x_squared * x_value * x_value
            + y_squared * y_value * y_value
            + xy * x_value * y_value
            + x_linear * x_value
            + y_linear * y_value
            + constant
        )

    candidates = [value(0.0, 0.0), value(1.0, 0.0), value(0.0, 1.0)]

    def add_univariate_stationary(
        quadratic: float,
        linear: float,
        evaluator: Any,
    ) -> None:
        if abs(quadratic) <= 1e-15:
            return
        point = -linear / (2.0 * quadratic)
        if 0.0 < point < 1.0:
            candidates.append(float(evaluator(point)))

    add_univariate_stationary(
        x_squared, x_linear, lambda point: value(point, 0.0)
    )
    add_univariate_stationary(
        y_squared, y_linear, lambda point: value(0.0, point)
    )
    diagonal_quadratic = x_squared + y_squared - xy
    diagonal_linear = -2.0 * y_squared + xy + x_linear - y_linear
    add_univariate_stationary(
        diagonal_quadratic,
        diagonal_linear,
        lambda point: value(point, 1.0 - point),
    )

    hessian = np.asarray(
        [[2.0 * x_squared, xy], [xy, 2.0 * y_squared]], dtype=float
    )
    linear = np.asarray([x_linear, y_linear], dtype=float)
    if abs(float(np.linalg.det(hessian))) > 1e-14:
        stationary = np.linalg.solve(hessian, -linear)
        if (
            stationary[0] > 0.0
            and stationary[1] > 0.0
            and stationary.sum() < 1.0
        ):
            candidates.append(value(float(stationary[0]), float(stationary[1])))
    return max(candidates)


def tangent_coefficient_error_at(
    alpha: float,
    beta: float,
    alpha_bounds: tuple[float, float],
    beta_bounds: tuple[float, float],
) -> float:
    """Convex-tangent upper bound on ``q_anchor-q_true``.

    The Horwitz polynomial is jointly convex in ``(alpha,beta)``.  Hence
    ``q_true >= q_anchor + grad(q_anchor) dot delta``.  Maximising the
    negative tangent over the parameter box needs only its four corners;
    for each corner the remaining maximum is one quadratic on the unit
    probability triangle and is evaluated exactly.
    """

    errors = []
    for endpoint_alpha in alpha_bounds:
        delta_alpha = endpoint_alpha - alpha
        for endpoint_beta in beta_bounds:
            delta_beta = endpoint_beta - beta
            errors.append(
                quadratic_maximum_on_unit_triangle(
                    x_squared=-2.0 * delta_beta * beta,
                    y_squared=-2.0 * delta_alpha * alpha,
                    xy=(
                        -2.0 * delta_alpha * (2.0 - beta)
                        - 2.0 * delta_beta * (2.0 - alpha)
                    ),
                    x_linear=2.0 * delta_beta,
                    y_linear=2.0 * delta_alpha,
                )
            )
    return max(0.0, max(errors))


def coefficient_error_at(
    alpha: float,
    beta: float,
    alpha_bounds: tuple[float, float],
    beta_bounds: tuple[float, float],
) -> float:
    """Return a one-sided uniform error around one anchor ellipse.

    Write a normalised probability point as ``(X,Y)`` with
    ``X,Y >= 0`` and ``X+Y <= 1``.  The one-sided perturbation
    ``q_anchor-q_true`` is bounded above by

    ``A X^2 + B Y^2 + C XY + D X + E Y``.

    All five coefficients are nonnegative, so its maximum lies on
    ``X+Y=1``.  Maximising the resulting univariate quadratic is both
    rigorous and substantially sharper than separately replacing every
    monomial by one.
    """

    al, au = alpha_bounds
    bl, bu = beta_bounds
    if not (al <= alpha <= au and bl <= beta <= bu):
        raise ValueError("the anchor must lie inside its parameter box")
    # Only the one-sided difference q_midpoint - q_true is needed: from
    # q_true <= 0 and q_midpoint - q_true <= error we obtain
    # q_midpoint <= error.  Negative coefficient differences can therefore
    # be discarded because every point monomial below is nonnegative.
    alpha2_error = max(0.0, alpha * alpha - al * al)
    beta2_error = max(0.0, beta * beta - bl * bl)
    centre_cross = cross_coefficient(alpha, beta)
    cross_error = max(
        0.0,
        centre_cross
        - min(
            cross_coefficient(a, b)
            for a in (al, au)
            for b in (bl, bu)
        ),
    )
    beta_linear_error = 2.0 * max(0.0, bu - beta)
    alpha_linear_error = 2.0 * max(0.0, au - alpha)

    # On Y=1-X the positive error polynomial is
    #     q2 X^2 + q1 X + q0,  0 <= X <= 1.
    q2 = beta2_error + alpha2_error - cross_error
    q1 = -2.0 * alpha2_error + cross_error + beta_linear_error - alpha_linear_error
    q0 = alpha2_error + alpha_linear_error
    candidates = [q0, q2 + q1 + q0]
    if abs(q2) > 1e-15:
        stationary = -q1 / (2.0 * q2)
        if 0.0 < stationary < 1.0:
            candidates.append(q2 * stationary * stationary + q1 * stationary + q0)
    positive_coefficient_error = max(candidates)
    tangent_error = tangent_coefficient_error_at(
        alpha, beta, alpha_bounds, beta_bounds
    )
    return min(positive_coefficient_error, tangent_error)


def coefficient_error(
    alpha_bounds: tuple[float, float],
    beta_bounds: tuple[float, float],
) -> tuple[float, float, float]:
    """Return midpoint parameters and their one-sided uniform error."""

    alpha = 0.5 * sum(alpha_bounds)
    beta = 0.5 * sum(beta_bounds)
    return (
        alpha,
        beta,
        coefficient_error_at(alpha, beta, alpha_bounds, beta_bounds),
    )


# Five interior anchors in each direction are supplemented by the four box
# corners.  Corner anchors are especially effective after an adaptive parent
# cell has already moved away from the degenerate alpha=1 or beta=1 boundary.
# A boundary corner whose quadratic is singular is simply made vacuous below.
ANCHOR_LOCATIONS = tuple(
    (alpha_fraction, beta_fraction)
    for alpha_fraction in (index / 6.0 for index in range(1, 6))
    for beta_fraction in (index / 6.0 for index in range(1, 6))
) + ((0.0, 0.0), (0.0, 1.0), (1.0, 0.0), (1.0, 1.0))


def box_anchor_relaxations(
    alpha_bounds: tuple[float, float],
    beta_bounds: tuple[float, float],
) -> list[tuple[float, float, float]]:
    """Valid outer ellipses whose intersection contains a box union.

    A point lying in any true ellipse from the parameter box obeys every
    anchor inequality after the corresponding one-sided error enlargement.
    Intersecting several such inequalities is therefore still an outer
    relaxation, and is much tighter than a lone midpoint ellipse on broad
    boxes.  Interior anchors avoid the degenerate ``alpha=1`` or ``beta=1``
    boundary of the closed inellipse family.
    """

    al, au = alpha_bounds
    bl, bu = beta_bounds
    result = []
    for alpha_fraction, beta_fraction in ANCHOR_LOCATIONS:
        alpha = al + alpha_fraction * (au - al)
        beta = bl + beta_fraction * (bu - bl)
        error = coefficient_error_at(
            alpha, beta, alpha_bounds, beta_bounds
        )
        result.append((alpha, beta, error))
    return result


def inellipse_soc_data(
    alpha: float, beta: float, error: float
) -> tuple[np.ndarray, np.ndarray, float]:
    cross = cross_coefficient(alpha, beta)
    quadratic = np.asarray(
        [[beta * beta, 0.5 * cross], [0.5 * cross, alpha * alpha]],
        dtype=float,
    )
    linear = np.asarray([-beta, -alpha], dtype=float)
    return quadratic_soc_data(quadratic, linear, 1.0 - error)


def quadratic_soc_data(
    quadratic: np.ndarray, linear: np.ndarray, constant: float
) -> tuple[np.ndarray, np.ndarray, float]:
    values, vectors = np.linalg.eigh(quadratic)
    if values.min() <= 0.0:
        raise ValueError(
            f"ellipse quadratic is singular: eigenvalues={values}"
        )
    square_root = (vectors * np.sqrt(values)) @ vectors.T
    shift = np.linalg.solve(quadratic, linear)
    radius_squared = float(linear @ shift - constant)
    if radius_squared < -1e-10:
        raise ValueError("invalid inellipse completion radius")
    radius = math.sqrt(max(0.0, radius_squared))
    # ``square_root`` and ``shift`` are computed in binary64. Treating their
    # entries as exact dyadic coefficients can otherwise make the completed
    # SOC microscopically smaller than the intended quadratic (source cell
    # 15818 exposed deficits up to 7.60e-16). Inflating only the scalar radius
    # is monotone: it weakens the SOC without changing its centre or axes.
    # The exact-rational enclosure audit checks the combined coefficient,
    # anchor-error, and completion error for every proof-producing source
    # anchor; 32 ULPs leave a visible exact margin while remaining negligible
    # on the optimization scale.
    for _ in range(32):
        radius = math.nextafter(radius, math.inf)
    return square_root, shift, radius


def center_inellipse_coefficients(
    center_x: float, center_y: float
) -> tuple[float, float, float, float, float, float]:
    """Polynomial inellipse coefficients in one compact barycentric chart.

    Put ``center_r=1-center_x-center_y``.  On the medial triangle
    ``0<=center_i<=1/2``, every inellipse is described by

    ``cxx*x^2+cyy*y^2+cxy*x*y+cx*x*p+cy*y*p+c0*p^2 <= 0``.

    The formula is Horwitz's equation after clearing its reciprocal
    denominators.  Unlike the reciprocal representation, it remains compact
    over all three former residual-coordinate charts.
    """

    a = float(center_x)
    b = float(center_y)
    u = 2.0 * a + 2.0 * b - 1.0
    return (
        4.0 * b * b,
        4.0 * a * a,
        -8.0 * a * b + 8.0 * a + 8.0 * b - 4.0,
        -4.0 * b * u,
        -4.0 * a * u,
        u * u,
    )


def center_box_coefficient_minima(
    center_x_bounds: tuple[float, float],
    center_y_bounds: tuple[float, float],
) -> tuple[float, float, float, float, float, float]:
    """Safe coefficientwise minima over a medial-triangle box intersection."""

    al, au = center_x_bounds
    bl, bu = center_y_bounds
    lower_u = max(0.0, 2.0 * al + 2.0 * bl - 1.0)
    upper_u = max(0.0, 2.0 * au + 2.0 * bu - 1.0)
    return (
        4.0 * bl * bl,
        4.0 * al * al,
        -8.0 * al * bl + 8.0 * al + 8.0 * bl - 4.0,
        -4.0 * bu * upper_u,
        -4.0 * au * upper_u,
        lower_u * lower_u,
    )


def center_coefficient_error_at(
    center_x: float,
    center_y: float,
    center_x_bounds: tuple[float, float],
    center_y_bounds: tuple[float, float],
) -> float:
    """One-sided uniform bound on ``q_anchor-q_true`` over a box."""

    anchor = center_inellipse_coefficients(center_x, center_y)
    lower = center_box_coefficient_minima(center_x_bounds, center_y_bounds)
    differences = [max(0.0, a - b) for a, b in zip(anchor, lower, strict=True)]
    return max(
        0.0,
        quadratic_maximum_on_unit_triangle(
            x_squared=differences[0],
            y_squared=differences[1],
            xy=differences[2],
            x_linear=differences[3],
            y_linear=differences[4],
            constant=differences[5],
        ),
    )


def center_box_anchor_relaxations(
    center_x_bounds: tuple[float, float],
    center_y_bounds: tuple[float, float],
) -> list[tuple[float, float, float]]:
    """Return a fixed-size family of valid single-chart outer ellipses."""

    al, au = center_x_bounds
    bl, bu = center_y_bounds
    result = []
    for x_fraction, y_fraction in ANCHOR_LOCATIONS:
        center_x = al + x_fraction * (au - al)
        center_y = bl + y_fraction * (bu - bl)
        if (
            center_x <= 1e-10
            or center_y <= 1e-10
            or center_x + center_y <= 0.5 + 1e-10
        ):
            # The anchor need not lie in the parameter box.  The Steiner
            # centre is a fixed positive-definite fallback, and the one-sided
            # coefficient error keeps its relaxation valid for the box.
            center_x = center_y = 1.0 / 3.0
        error = center_coefficient_error_at(
            center_x,
            center_y,
            center_x_bounds,
            center_y_bounds,
        )
        result.append((center_x, center_y, error))
    return result


def center_inellipse_soc_data(
    center_x: float, center_y: float, error: float
) -> tuple[np.ndarray, np.ndarray, float]:
    cxx, cyy, cxy, cx, cy, c0 = center_inellipse_coefficients(
        center_x, center_y
    )
    quadratic = np.asarray(
        [[cxx, 0.5 * cxy], [0.5 * cxy, cyy]], dtype=float
    )
    linear = np.asarray([0.5 * cx, 0.5 * cy], dtype=float)
    return quadratic_soc_data(quadratic, linear, c0 - error)


def center_coefficientwise_box_soc_data(
    center_x_bounds: tuple[float, float],
    center_y_bounds: tuple[float, float],
) -> tuple[np.ndarray, np.ndarray, float] | None:
    coefficients = center_box_coefficient_minima(
        center_x_bounds, center_y_bounds
    )
    cxx, cyy, cxy, cx, cy, c0 = coefficients
    quadratic = np.asarray(
        [[cxx, 0.5 * cxy], [0.5 * cxy, cyy]], dtype=float
    )
    if np.linalg.eigvalsh(quadratic).min() <= 1e-10:
        return None
    linear = np.asarray([0.5 * cx, 0.5 * cy], dtype=float)
    try:
        return quadratic_soc_data(quadratic, linear, c0)
    except ValueError:
        return None


def coefficientwise_box_soc_data(
    alpha_bounds: tuple[float, float],
    beta_bounds: tuple[float, float],
) -> tuple[np.ndarray, np.ndarray, float] | None:
    """Convex quadratic lower envelope of every ellipse in a box.

    All homogenised point monomials are nonnegative.  Taking the minimum of
    each coefficient therefore gives a necessary inequality.  It is an SOCP
    whenever the resulting two-by-two quadratic block is positive definite.
    """

    al, au = alpha_bounds
    bl, bu = beta_bounds
    cross_minimum = min(
        cross_coefficient(alpha, beta)
        for alpha in (al, au)
        for beta in (bl, bu)
    )
    quadratic = np.asarray(
        [[bl * bl, 0.5 * cross_minimum], [0.5 * cross_minimum, al * al]],
        dtype=float,
    )
    if np.linalg.eigvalsh(quadratic).min() <= 1e-10:
        return None
    linear = np.asarray([-bu, -au], dtype=float)
    return quadratic_soc_data(quadratic, linear, 1.0)


class ReusableBoxOracle:
    """DPP SOCP reused across every parameter box of one chart cell."""

    def __init__(
        self,
        terminal: np.ndarray,
        support_weight: float,
        prefix_order: tuple[int, int, int, int],
        pairs: tuple[tuple[Column, Column], ...],
        coordinate_cases: tuple[str, ...],
    ) -> None:
        self.weights = np.trace(terminal[:3], axis1=1, axis2=2).real
        self.pairs = pairs
        self.coordinate_cases = coordinate_cases
        constraints: list[cp.Constraint] = []
        self.statistics = cp.Variable((4, 4, 3), nonneg=True)
        self.probability = cp.sum(self.statistics, axis=2)
        constraints.append(cp.sum(self.statistics) == 1.0)
        constraints.extend(
            self.statistics[:, :, t] <= self.weights[t] * self.probability
            for t in range(3)
        )
        self.prefix = [cp.sum(self.probability[z, :]) for z in OUTCOMES]
        constraints.extend(
            self.prefix[prefix_order[index]]
            >= self.prefix[prefix_order[index + 1]]
            for index in range(3)
        )
        constraints.append(self.prefix[prefix_order[0]] >= 0.25)
        for rank in range(1, 4):
            constraints.append(
                self.prefix[prefix_order[rank]] <= 1.0 / (rank + 1.0)
            )

        # For each selected pair, a fixed anchor grid describes outer
        # ellipses and the final triple describes the coefficientwise lower
        # envelope.
        self.soc_parameters: list[
            tuple[
                list[tuple[cp.Parameter, cp.Parameter, cp.Parameter]],
                tuple[cp.Parameter, cp.Parameter, cp.Parameter],
            ]
        ] = []
        for (first, second), coordinate_case in zip(
            pairs, coordinate_cases, strict=True
        ):
            anchors = [
                (
                    cp.Parameter((2, 2)),
                    cp.Parameter(2),
                    cp.Parameter(nonneg=True),
                )
                for _ in range(len(ANCHOR_LOCATIONS))
            ]
            lower = (
                cp.Parameter((2, 2)),
                cp.Parameter(2),
                cp.Parameter(nonneg=True),
            )
            self.soc_parameters.append((anchors, lower))
            first_cap = self.weights[first[2]]
            second_cap = self.weights[second[2]]
            for z in OUTCOMES:
                first_value = selected_statistic(
                    self.statistics, self.probability, self.weights, first, z
                )
                second_value = selected_statistic(
                    self.statistics, self.probability, self.weights, second, z
                )
                x_value = first_value / first_cap
                y_value = second_value / second_cap
                residual = self.prefix[z] - x_value - y_value
                constraints.append(residual >= 0.0)
                if coordinate_case == "xy":
                    point = cp.hstack([x_value, y_value])
                elif coordinate_case == "xr":
                    point = cp.hstack([x_value, residual])
                elif coordinate_case == "yr":
                    point = cp.hstack([y_value, residual])
                else:
                    raise ValueError(f"unknown coordinate case {coordinate_case!r}")
                for root, offset, radius in (*anchors, lower):
                    transformed = root @ point + offset * self.prefix[z]
                    constraints.append(
                        cp.SOC(radius * self.prefix[z], transformed)
                    )

        terminal_statistics = [
            [
                sum(
                    self.statistics[z, y, t]
                    for z, y in PATHS
                    if (z ^ y) == syndrome
                )
                for t in range(3)
            ]
            for syndrome in OUTCOMES
        ]
        inverse = reconstruction_matrix(terminal)
        terminal_prior: list[cp.Expression] = []
        terminal_vector: list[cp.Expression] = []
        normal = cp.Variable(4)
        for syndrome in OUTCOMES:
            reconstructed = inverse @ cp.hstack(terminal_statistics[syndrome])
            vector = cp.hstack(
                [reconstructed[1], reconstructed[2], normal[syndrome]]
            )
            constraints.extend(
                (
                    reconstructed[0]
                    == sum(
                        self.probability[z, z ^ syndrome]
                        for z in OUTCOMES
                    ),
                    cp.SOC(reconstructed[0], vector),
                )
            )
            terminal_prior.append(reconstructed[0])
            terminal_vector.append(vector)

        self.audit = sum(terminal_statistics[s][s] for s in range(3))
        cap = filled_effect_weights(float(self.weights.max()))
        constraints.append(
            self.audit
            <= sum(
                cap[index] * self.prefix[prefix_order[index]]
                for index in OUTCOMES
            )
        )
        dual_scalar = cp.Variable(nonneg=True)
        dual_vector = cp.Variable(3)
        constraints.append(cp.SOC(dual_scalar, dual_vector))
        constraints.extend(
            cp.SOC(
                dual_scalar - terminal_prior[s],
                dual_vector - terminal_vector[s],
            )
            for s in OUTCOMES
        )
        constraints.append(self.audit == dual_scalar)
        constraints.append(
            self.audit
            <= sum(self.weights[t] * terminal_prior[t] for t in range(3))
        )
        self.returned = hellinger_hypograph(
            [self.probability[z, y] for z, y in PATHS], constraints
        )
        self.score = (
            support_weight * self.audit
            + (1.0 - support_weight) * self.returned
        )
        self.problem = cp.Problem(cp.Maximize(self.score), constraints)
        if not self.problem.is_dpp():
            raise RuntimeError("reusable box oracle is not DPP")

    @staticmethod
    def _assign_soc(
        targets: tuple[cp.Parameter, cp.Parameter, cp.Parameter],
        data: tuple[np.ndarray, np.ndarray, float] | None,
    ) -> None:
        root, offset, radius = targets
        if data is None:
            # ||0|| <= a is vacuous because every prefix prior is nonnegative.
            root.value = np.zeros((2, 2))
            offset.value = np.zeros(2)
            radius.value = 1.0
            return
        root_value, shift_value, radius_value = data
        root.value = root_value
        offset.value = root_value @ shift_value
        radius.value = radius_value

    def solve(self, box: Box, safety: float, capture: bool) -> dict[str, Any]:
        ellipse_report = []
        for pair_index, (anchor_targets, lower) in enumerate(self.soc_parameters):
            alpha_bounds = box[f"k{pair_index}_alpha"]
            beta_bounds = box[f"k{pair_index}_beta"]
            anchor_values = box_anchor_relaxations(alpha_bounds, beta_bounds)
            anchor_report = []
            for (alpha, beta, error), target in zip(
                anchor_values, anchor_targets, strict=True
            ):
                try:
                    anchor_data = inellipse_soc_data(alpha, beta, error)
                except ValueError:
                    # Clean ternary POVMs include degenerate rank-one limits.
                    # The other anchors still give a valid outer cover; making
                    # this singular anchor vacuous cannot exclude a physical
                    # point.
                    anchor_data = None
                self._assign_soc(target, anchor_data)
                anchor_report.append(
                    {
                        "alpha": alpha,
                        "beta": beta,
                        "coefficient_error": error,
                        "active": anchor_data is not None,
                        "radius": None if anchor_data is None else anchor_data[2],
                    }
                )
            lower_data = coefficientwise_box_soc_data(
                alpha_bounds, beta_bounds
            )
            self._assign_soc(lower, lower_data)
            ellipse_report.append(
                {
                    "anchors": anchor_report,
                    "coefficientwise_soc_active": lower_data is not None,
                    "coordinate_case": self.coordinate_cases[pair_index],
                }
            )
        try:
            self.problem.solve(
                solver="CLARABEL",
                tol_gap_abs=2e-8,
                tol_gap_rel=2e-8,
                tol_feas=2e-8,
                max_iter=1000,
                warm_start=True,
                ignore_dpp=False,
            )
        except cp.SolverError as error_value:
            return {
                "status": "solver_error",
                "error": str(error_value),
                "bound": math.inf,
            }
        if self.problem.status in {cp.INFEASIBLE, cp.INFEASIBLE_INACCURATE}:
            return {"status": self.problem.status, "bound": -math.inf}
        if self.problem.status not in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}:
            return {"status": self.problem.status, "bound": math.inf}
        result: dict[str, Any] = {
            "status": self.problem.status,
            "raw_value": float(self.problem.value),
            "bound": float(self.problem.value) + safety,
            "audit": float(self.audit.value),
            "return": float(self.returned.value),
            "ellipse_relaxations": ellipse_report,
            "iterations": self.problem.solver_stats.num_iters,
            "solve_time": self.problem.solver_stats.solve_time,
        }
        if capture:
            result.update(
                {
                    "prefix": [float(value.value) for value in self.prefix],
                    "probability": np.asarray(self.probability.value).tolist(),
                    "statistics": np.asarray(self.statistics.value).tolist(),
                }
            )
        return result


_ORACLE_CACHE: dict[tuple[Any, ...], ReusableBoxOracle] = {}


def solve_node_reusable(
    box: Box,
    terminal: np.ndarray,
    support_weight: float,
    prefix_order: tuple[int, int, int, int],
    pairs: tuple[tuple[Column, Column], ...],
    coordinate_cases: tuple[str, ...],
    safety: float,
    capture: bool = False,
) -> dict[str, Any]:
    key = (
        terminal.shape,
        terminal.tobytes(),
        float(support_weight),
        prefix_order,
        pairs,
        coordinate_cases,
    )
    oracle = _ORACLE_CACHE.get(key)
    if oracle is None:
        oracle = ReusableBoxOracle(
            terminal,
            support_weight,
            prefix_order,
            pairs,
            coordinate_cases,
        )
        _ORACLE_CACHE[key] = oracle
    return oracle.solve(box, safety, capture)


def solve_node(
    box: Box,
    terminal: np.ndarray,
    support_weight: float,
    prefix_order: tuple[int, int, int, int],
    pairs: tuple[tuple[Column, Column], ...],
    coordinate_cases: tuple[str, ...],
    safety: float,
    capture: bool = False,
) -> dict[str, Any]:
    weights = np.trace(terminal[:3], axis1=1, axis2=2).real
    constraints: list[cp.Constraint] = []
    statistics = cp.Variable((4, 4, 3), nonneg=True)
    probability = cp.sum(statistics, axis=2)
    constraints.append(cp.sum(statistics) == 1.0)
    constraints.extend(
        statistics[:, :, t] <= weights[t] * probability for t in range(3)
    )
    prefix = [cp.sum(probability[z, :]) for z in OUTCOMES]
    constraints.extend(
        prefix[prefix_order[index]] >= prefix[prefix_order[index + 1]]
        for index in range(3)
    )
    constraints.append(prefix[prefix_order[0]] >= 0.25)
    for rank in range(1, 4):
        constraints.append(prefix[prefix_order[rank]] <= 1.0 / (rank + 1.0))

    ellipse_report = []
    for pair_index, ((first, second), coordinate_case) in enumerate(
        zip(pairs, coordinate_cases, strict=True)
    ):
        alpha, beta, error = coefficient_error(
            box[f"k{pair_index}_alpha"], box[f"k{pair_index}_beta"]
        )
        square_root, shift, radius = inellipse_soc_data(alpha, beta, error)
        lower_envelope = coefficientwise_box_soc_data(
            box[f"k{pair_index}_alpha"], box[f"k{pair_index}_beta"]
        )
        first_cap, second_cap = weights[first[2]], weights[second[2]]
        for z in OUTCOMES:
            first_value = selected_statistic(
                statistics, probability, weights, first, z
            )
            second_value = selected_statistic(
                statistics, probability, weights, second, z
            )
            x_value = first_value / first_cap
            y_value = second_value / second_cap
            residual = prefix[z] - x_value - y_value
            constraints.append(residual >= 0.0)
            if coordinate_case == "xy":
                point = cp.hstack([x_value, y_value])
            elif coordinate_case == "xr":
                point = cp.hstack([x_value, residual])
            elif coordinate_case == "yr":
                point = cp.hstack([y_value, residual])
            else:
                raise ValueError(f"unknown coordinate case {coordinate_case!r}")
            transformed = square_root @ (point + shift * prefix[z])
            constraints.append(cp.SOC(radius * prefix[z], transformed))
            if lower_envelope is not None:
                lower_root, lower_shift, lower_radius = lower_envelope
                lower_transformed = lower_root @ (
                    point + lower_shift * prefix[z]
                )
                constraints.append(
                    cp.SOC(lower_radius * prefix[z], lower_transformed)
                )
        ellipse_report.append(
            {
                "alpha_midpoint": alpha,
                "beta_midpoint": beta,
                "coefficient_error": error,
                "radius": radius,
                "coefficientwise_soc_active": lower_envelope is not None,
                "coordinate_case": coordinate_case,
            }
        )

    terminal_statistics = [
        [
            sum(statistics[z, y, t] for z, y in PATHS if (z ^ y) == syndrome)
            for t in range(3)
        ]
        for syndrome in OUTCOMES
    ]
    inverse = reconstruction_matrix(terminal)
    terminal_prior: list[cp.Expression] = []
    terminal_vector: list[cp.Expression] = []
    normal = cp.Variable(4)
    for syndrome in OUTCOMES:
        reconstructed = inverse @ cp.hstack(terminal_statistics[syndrome])
        vector = cp.hstack([reconstructed[1], reconstructed[2], normal[syndrome]])
        constraints.extend(
            (
                reconstructed[0]
                == sum(probability[z, z ^ syndrome] for z in OUTCOMES),
                cp.SOC(reconstructed[0], vector),
            )
        )
        terminal_prior.append(reconstructed[0])
        terminal_vector.append(vector)
    audit = sum(terminal_statistics[s][s] for s in range(3))
    cap = filled_effect_weights(float(weights.max()))
    constraints.append(
        audit
        <= sum(cap[index] * prefix[prefix_order[index]] for index in OUTCOMES)
    )
    dual_scalar = cp.Variable(nonneg=True)
    dual_vector = cp.Variable(3)
    constraints.append(cp.SOC(dual_scalar, dual_vector))
    constraints.extend(
        cp.SOC(
            dual_scalar - terminal_prior[s],
            dual_vector - terminal_vector[s],
        )
        for s in OUTCOMES
    )
    constraints.append(audit == dual_scalar)
    constraints.append(
        audit
        <= sum(weights[t] * terminal_prior[t] for t in range(3))
    )
    returned = hellinger_hypograph(
        [probability[z, y] for z, y in PATHS], constraints
    )
    objective = support_weight * audit + (1.0 - support_weight) * returned
    problem = cp.Problem(cp.Maximize(objective), constraints)
    try:
        problem.solve(
            solver="CLARABEL",
            tol_gap_abs=2e-8,
            tol_gap_rel=2e-8,
            tol_feas=2e-8,
            max_iter=1000,
        )
    except cp.SolverError as error_value:
        return {
            "status": "solver_error",
            "error": str(error_value),
            "bound": math.inf,
        }
    if problem.status in {cp.INFEASIBLE, cp.INFEASIBLE_INACCURATE}:
        return {"status": problem.status, "bound": -math.inf}
    if problem.status not in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}:
        return {"status": problem.status, "bound": math.inf}
    result: dict[str, Any] = {
        "status": problem.status,
        "raw_value": float(problem.value),
        "bound": float(problem.value) + safety,
        "audit": float(audit.value),
        "return": float(returned.value),
        "ellipse_relaxations": ellipse_report,
        "iterations": problem.solver_stats.num_iters,
        "solve_time": problem.solver_stats.solve_time,
    }
    if capture:
        result.update(
            {
                "prefix": [float(value.value) for value in prefix],
                "probability": np.asarray(probability.value, dtype=float).tolist(),
                "statistics": np.asarray(statistics.value, dtype=float).tolist(),
            }
        )
    return result


def branch_variable(
    box: Box,
    root: Box,
    active_variables: tuple[str, ...] | None = None,
) -> str:
    names = active_variables or tuple(box)
    return max(
        names,
        key=lambda name: (box[name][1] - box[name][0])
        / max(root[name][1] - root[name][0], 1e-15),
    )


def deserialise_box(payload: dict[str, Any]) -> Box:
    return {
        name: (float(bounds[0]), float(bounds[1]))
        for name, bounds in payload.items()
    }


def split_box(box: Box, name: str) -> tuple[Box, Box]:
    lower, upper = box[name]
    middle = 0.5 * (lower + upper)
    first, second = box.copy(), box.copy()
    first[name] = (lower, middle)
    second[name] = (middle, upper)
    return first, second


def serialise_box(box: Box) -> dict[str, list[float]]:
    return {name: [bounds[0], bounds[1]] for name, bounds in box.items()}


def write_payload(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lambda", dest="weight", type=float, default=0.6)
    parser.add_argument("--fixed-three-povm-weights", type=float, nargs=3, required=True)
    parser.add_argument("--prefix-order", type=int, nargs=4, required=True)
    parser.add_argument("--pair", action="append", default=[])
    parser.add_argument(
        "--coordinate-case",
        action="append",
        choices=("xy", "xr", "yr"),
        default=[],
        help="supply one displayed-coordinate chart per selected pair",
    )
    parser.add_argument("--reciprocal-cap", type=float, default=2.0)
    parser.add_argument("--target", type=float, default=0.76591)
    parser.add_argument("--max-nodes", type=int, default=101)
    parser.add_argument("--safety", type=float, default=2e-6)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument(
        "--branch-rule",
        choices=("width", "strong"),
        default="width",
        help="strong branching tests every parameter split before choosing one",
    )
    parser.add_argument("--resume", type=Path)
    parser.add_argument(
        "--root-box",
        type=Path,
        help=(
            "JSON object of parameter intervals; useful for conditioning an "
            "adaptive child cover on an already certified parent cell"
        ),
    )
    parser.add_argument(
        "--branch-variable",
        action="append",
        default=[],
        help=(
            "parameter name eligible for subdivision; repeat as needed. "
            "Unspecified means all parameters"
        ),
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    order = tuple(args.prefix_order)
    if sorted(order) != list(OUTCOMES):
        raise ValueError("prefix order must be a permutation")
    if args.reciprocal_cap < 2.0:
        raise ValueError("the clean-POVM relabelling proof requires reciprocal-cap >= 2")
    pair_text = args.pair or [
        "b_0,d_1_0",
        "b_1,b_9",
        "b_3,b_6",
        "b_4,d_0_2",
    ]
    pairs = tuple(parse_pair(text) for text in pair_text)
    coordinate_cases = tuple(args.coordinate_case or ["xy"] * len(pairs))
    if len(coordinate_cases) != len(pairs):
        raise ValueError("supply exactly one --coordinate-case per pair")
    terminal = canonical_three_effect_povm(
        np.asarray(args.fixed_three_povm_weights)
    )
    expected_names = set(initial_box(len(pairs), args.reciprocal_cap))
    if args.root_box is None:
        root = initial_box(len(pairs), args.reciprocal_cap)
    else:
        root_payload = json.loads(args.root_box.read_text(encoding="utf-8"))
        if "box" in root_payload:
            root_payload = root_payload["box"]
        root = deserialise_box(root_payload)
        if set(root) != expected_names:
            raise ValueError(
                "root-box keys must be exactly " + ", ".join(sorted(expected_names))
            )
        if any(lower < 1.0 or upper > args.reciprocal_cap or lower > upper
               for lower, upper in root.values()):
            raise ValueError("root-box intervals lie outside the clean-POVM domain")
    active_variables = tuple(args.branch_variable) or tuple(root)
    if not active_variables or not set(active_variables) <= set(root):
        raise ValueError("branch variables must be names from the root box")
    solved = 0
    leaves: list[dict[str, Any]] = []
    queue: list[tuple[float, int, Box, dict[str, Any]]] = []
    counter = 0
    if args.resume is None:
        result = solve_node_reusable(
            root,
            terminal,
            args.weight,
            order,
            pairs,
            coordinate_cases,
            args.safety,
        )
        solved = 1
        heapq.heappush(queue, (-float(result["bound"]), counter, root, result))
        print("root", result, flush=True)
    else:
        previous = json.loads(args.resume.read_text(encoding="utf-8"))
        solved = int(previous["solved_nodes"])
        leaves = list(previous["leaves"])
        for item in previous["open_nodes"]:
            counter += 1
            box = {
                name: (float(bounds[0]), float(bounds[1]))
                for name, bounds in item["box"].items()
            }
            result = {key: value for key, value in item.items() if key != "box"}
            heapq.heappush(queue, (-float(result["bound"]), counter, box, result))
        print("resumed", solved, len(queue), len(leaves), flush=True)

    executor = ProcessPoolExecutor(max_workers=args.workers) if args.workers > 1 else None
    oracle_solves = solved
    try:
        while queue and solved < args.max_nodes:
            neg_bound, _, box, result = heapq.heappop(queue)
            if -neg_bound <= args.target:
                leaves.append({"box": serialise_box(box), **result, "reason": "target"})
                continue
            candidate_names = (
                [branch_variable(box, root, active_variables)]
                if args.branch_rule == "width"
                else [
                    name
                    for name in active_variables
                    for bounds in (box[name],)
                    if bounds[1] - bounds[0] > 1e-14
                ]
            )
            candidate_children = [
                child
                for name in candidate_names
                for child in split_box(box, name)
            ]
            if executor is None:
                candidate_results = [
                    solve_node_reusable(
                        child,
                        terminal,
                        args.weight,
                        order,
                        pairs,
                        coordinate_cases,
                        args.safety,
                    )
                    for child in candidate_children
                ]
            else:
                candidate_results = list(
                    executor.map(
                        solve_node_reusable,
                        candidate_children,
                        [terminal] * len(candidate_children),
                        [args.weight] * len(candidate_children),
                        [order] * len(candidate_children),
                        [pairs] * len(candidate_children),
                        [coordinate_cases] * len(candidate_children),
                        [args.safety] * len(candidate_children),
                    )
                )
            oracle_solves += len(candidate_results)
            grouped = [
                (
                    name,
                    tuple(candidate_children[2 * index : 2 * index + 2]),
                    candidate_results[2 * index : 2 * index + 2],
                )
                for index, name in enumerate(candidate_names)
            ]
            name, children, results = min(
                grouped,
                key=lambda item: (
                    max(float(result["bound"]) for result in item[2]),
                    sum(float(result["bound"]) for result in item[2]),
                ),
            )
            for child, child_result in zip(children, results):
                solved += 1
                counter += 1
                child_bound = float(child_result["bound"])
                if child_bound <= args.target:
                    leaves.append(
                        {"box": serialise_box(child), **child_result, "reason": "target"}
                    )
                else:
                    heapq.heappush(queue, (-child_bound, counter, child, child_result))
            current = -queue[0][0] if queue else max(
                [float(leaf["bound"]) for leaf in leaves], default=-math.inf
            )
            print(solved, name, current, len(queue), len(leaves), flush=True)
    finally:
        if executor is not None:
            executor.shutdown()

    top_open_solution = None
    if queue:
        _, _, top_box, _ = queue[0]
        top_open_solution = {
            "box": serialise_box(top_box),
            **solve_node_reusable(
                top_box,
                terminal,
                args.weight,
                order,
                pairs,
                coordinate_cases,
                args.safety,
                capture=True,
            ),
        }
    open_nodes = [
        {"box": serialise_box(box), **result}
        for _, _, box, result in sorted(queue)
    ]
    payload = {
        "weight": args.weight,
        "terminal_weights": args.fixed_three_povm_weights,
        "prefix_order": list(order),
        "pairs": [
            [render_column(first), render_column(second)] for first, second in pairs
        ],
        "coordinate_cases": list(coordinate_cases),
        "reciprocal_cap": args.reciprocal_cap,
        "target": args.target,
        "safety": args.safety,
        "branch_rule": args.branch_rule,
        "root_box": serialise_box(root),
        "branch_variables": list(active_variables),
        "oracle_solves": oracle_solves,
        "complete": not open_nodes,
        "solved_nodes": solved,
        "maximum_open_bound": max(
            [float(item["bound"]) for item in open_nodes], default=-math.inf
        ),
        "maximum_leaf_bound": max(
            [float(item["bound"]) for item in leaves], default=-math.inf
        ),
        "open_nodes": open_nodes,
        "leaves": leaves,
        "top_open_solution": top_open_solution,
        "scope": (
            "one coordinate-chart cell for every selected pair; a complete "
            "fixed-pair certificate must cover the Cartesian product of the "
            "xy/xr/yr chart unions"
        ),
    }
    write_payload(args.output, payload)
    print(
        json.dumps(
            {key: value for key, value in payload.items() if key not in {"open_nodes", "leaves"}},
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

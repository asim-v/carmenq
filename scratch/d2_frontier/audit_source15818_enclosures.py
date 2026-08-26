"""Exact-rational audit of the inellipse SOCs used by source cell 15818.

The source oracle constructs ellipse coefficients and a completed-square SOC
with binary64 arithmetic.  Reinterpreting those numbers as exact dyadic
rationals is not automatically an outer approximation.  This audit bounds
both errors together on the full normalized probability triangle.

No optimizer is used.  Every maximum is the exact maximum of a rational
quadratic on ``x >= 0, y >= 0, x + y <= 1``.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path
from typing import Callable

import numpy as np

from pairwise_inellipse_box_cover import (
    box_anchor_relaxations,
    coefficientwise_box_soc_data,
    cross_coefficient,
    inellipse_soc_data,
)


SOURCE_INDEX = 15818
DEFAULT_FRONTIER = Path(
    "scratch/d2_frontier/"
    "ternary_reconstructed_depth4_g2_top_leaf_bbb_p1_s92_l055.json"
)


def q(value: object) -> Fraction:
    """Exact rational denoted by one finite binary64 value."""

    return Fraction(*float(value).as_integer_ratio())


def cross(alpha: Fraction, beta: Fraction) -> Fraction:
    return -2 * (alpha * beta - 2 * alpha - 2 * beta + 2)


def quadratic_value(
    coefficients: tuple[Fraction, Fraction, Fraction, Fraction, Fraction, Fraction],
    x: Fraction,
    y: Fraction,
) -> Fraction:
    xx, yy, xy, x_linear, y_linear, constant = coefficients
    return (
        xx * x * x
        + yy * y * y
        + xy * x * y
        + x_linear * x
        + y_linear * y
        + constant
    )


def rational_quadratic_maximum_on_unit_triangle(
    coefficients: tuple[Fraction, Fraction, Fraction, Fraction, Fraction, Fraction],
) -> Fraction:
    """Return the exact maximum of a rational quadratic on the unit triangle."""

    xx, yy, xy, x_linear, y_linear, _ = coefficients
    value: Callable[[Fraction, Fraction], Fraction] = lambda x, y: quadratic_value(
        coefficients, x, y
    )
    candidates = [value(Fraction(0), Fraction(0)), value(Fraction(1), Fraction(0)),
                  value(Fraction(0), Fraction(1))]

    def add_edge_stationary(
        quadratic: Fraction,
        linear: Fraction,
        evaluator: Callable[[Fraction], Fraction],
    ) -> None:
        if quadratic == 0:
            return
        point = -linear / (2 * quadratic)
        if 0 < point < 1:
            candidates.append(evaluator(point))

    add_edge_stationary(xx, x_linear, lambda point: value(point, Fraction(0)))
    add_edge_stationary(yy, y_linear, lambda point: value(Fraction(0), point))
    diagonal_quadratic = xx + yy - xy
    diagonal_linear = -2 * yy + xy + x_linear - y_linear
    add_edge_stationary(
        diagonal_quadratic,
        diagonal_linear,
        lambda point: value(point, 1 - point),
    )

    determinant = 4 * xx * yy - xy * xy
    if determinant != 0:
        stationary_x = (-2 * yy * x_linear + xy * y_linear) / determinant
        stationary_y = (xy * x_linear - 2 * xx * y_linear) / determinant
        if stationary_x > 0 and stationary_y > 0 and stationary_x + stationary_y < 1:
            candidates.append(value(stationary_x, stationary_y))
    return max(candidates)


def positive_coefficient_error(
    alpha: Fraction,
    beta: Fraction,
    alpha_bounds: tuple[Fraction, Fraction],
    beta_bounds: tuple[Fraction, Fraction],
) -> Fraction:
    al, au = alpha_bounds
    bl, bu = beta_bounds
    alpha2_error = max(Fraction(0), alpha * alpha - al * al)
    beta2_error = max(Fraction(0), beta * beta - bl * bl)
    cross_error = max(
        Fraction(0),
        cross(alpha, beta)
        - min(cross(a, b) for a in (al, au) for b in (bl, bu)),
    )
    beta_linear_error = 2 * max(Fraction(0), bu - beta)
    alpha_linear_error = 2 * max(Fraction(0), au - alpha)
    return rational_quadratic_maximum_on_unit_triangle(
        (
            beta2_error,
            alpha2_error,
            cross_error,
            beta_linear_error,
            alpha_linear_error,
            Fraction(0),
        )
    )


def tangent_error(
    alpha: Fraction,
    beta: Fraction,
    alpha_bounds: tuple[Fraction, Fraction],
    beta_bounds: tuple[Fraction, Fraction],
) -> Fraction:
    errors: list[Fraction] = []
    for endpoint_alpha in alpha_bounds:
        delta_alpha = endpoint_alpha - alpha
        for endpoint_beta in beta_bounds:
            delta_beta = endpoint_beta - beta
            errors.append(
                rational_quadratic_maximum_on_unit_triangle(
                    (
                        -2 * delta_beta * beta,
                        -2 * delta_alpha * alpha,
                        -2 * delta_alpha * (2 - beta)
                        - 2 * delta_beta * (2 - alpha),
                        2 * delta_beta,
                        2 * delta_alpha,
                        Fraction(0),
                    )
                )
            )
    return max(Fraction(0), max(errors))


def soc_polynomial_coefficients(
    root: np.ndarray, offset: np.ndarray, radius: float
) -> tuple[Fraction, Fraction, Fraction, Fraction, Fraction, Fraction]:
    r00, r01 = q(root[0, 0]), q(root[0, 1])
    r10, r11 = q(root[1, 0]), q(root[1, 1])
    o0, o1 = q(offset[0]), q(offset[1])
    r = q(radius)
    return (
        r00 * r00 + r10 * r10,
        r01 * r01 + r11 * r11,
        2 * (r00 * r01 + r10 * r11),
        2 * (r00 * o0 + r10 * o1),
        2 * (r01 * o0 + r11 * o1),
        o0 * o0 + o1 * o1 - r * r,
    )


def subtract_coefficients(
    first: tuple[Fraction, ...], second: tuple[Fraction, ...]
) -> tuple[Fraction, Fraction, Fraction, Fraction, Fraction, Fraction]:
    return tuple(a - b for a, b in zip(first, second, strict=True))  # type: ignore[return-value]


def encode(value: Fraction) -> list[int]:
    return [value.numerator, value.denominator]


def audit_box(
    alpha_bounds_float: tuple[float, float],
    beta_bounds_float: tuple[float, float],
) -> dict[str, object]:
    alpha_bounds = tuple(q(value) for value in alpha_bounds_float)
    beta_bounds = tuple(q(value) for value in beta_bounds_float)
    anchors: list[dict[str, object]] = []
    minimum_margin: Fraction | None = None

    for index, (alpha_float, beta_float, error_float) in enumerate(
        box_anchor_relaxations(alpha_bounds_float, beta_bounds_float)
    ):
        alpha, beta, error = q(alpha_float), q(beta_float), q(error_float)
        positive = positive_coefficient_error(alpha, beta, alpha_bounds, beta_bounds)
        tangent = tangent_error(alpha, beta, alpha_bounds, beta_bounds)
        analytic = min(positive, tangent)
        root, shift, radius = inellipse_soc_data(
            alpha_float, beta_float, error_float
        )
        offset = root @ shift
        actual = soc_polynomial_coefficients(root, offset, radius)
        intended = (
            beta * beta,
            alpha * alpha,
            cross(alpha, beta),
            -2 * beta,
            -2 * alpha,
            1 - error,
        )
        completion_excess = rational_quadratic_maximum_on_unit_triangle(
            subtract_coefficients(actual, intended)
        )
        margin = error - analytic - completion_excess
        minimum_margin = margin if minimum_margin is None else min(minimum_margin, margin)
        anchors.append(
            {
                "index": index,
                "alpha": encode(alpha),
                "beta": encode(beta),
                "stored_error": encode(error),
                "positive_bound": encode(positive),
                "tangent_bound": encode(tangent),
                "analytic_bound": encode(analytic),
                "completion_excess": encode(completion_excess),
                "margin": encode(margin),
                "certified_outer": margin >= 0,
            }
        )

    lower_data = coefficientwise_box_soc_data(
        alpha_bounds_float, beta_bounds_float
    )
    if lower_data is None:
        lower_report: dict[str, object] = {"present": False, "certified_outer": True}
    else:
        root, shift, radius = lower_data
        offset = root @ shift
        actual = soc_polynomial_coefficients(root, offset, radius)
        al, _au = alpha_bounds
        bl, _bu = beta_bounds
        cross_minimum = min(
            cross(a, b) for a in alpha_bounds for b in beta_bounds
        )
        ideal = (
            bl * bl,
            al * al,
            cross_minimum,
            -2 * beta_bounds[1],
            -2 * alpha_bounds[1],
            Fraction(1),
        )
        excess = rational_quadratic_maximum_on_unit_triangle(
            subtract_coefficients(actual, ideal)
        )
        lower_report = {
            "present": True,
            "completion_excess": encode(excess),
            "certified_outer": excess <= 0,
        }

    return {
        "alpha_bounds": [encode(value) for value in alpha_bounds],
        "beta_bounds": [encode(value) for value in beta_bounds],
        "anchor_count": len(anchors),
        "all_anchors_certified": all(bool(item["certified_outer"]) for item in anchors),
        "minimum_anchor_margin": encode(minimum_margin or Fraction(0)),
        "anchors": anchors,
        "coefficientwise_lower": lower_report,
        "all_certified": all(bool(item["certified_outer"]) for item in anchors)
        and bool(lower_report["certified_outer"]),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--frontier", type=Path, default=DEFAULT_FRONTIER)
    parser.add_argument("--source-index", type=int, default=SOURCE_INDEX)
    args = parser.parse_args()
    source = json.loads(args.frontier.read_text(encoding="utf-8"))
    if args.source_index >= len(source["cells"]):
        raise SystemExit("source index is outside the frontier cell array")
    box = source["box"]
    report = audit_box(tuple(box["terminal_alpha"]), tuple(box["terminal_beta"]))
    payload = {
        "schema": "carmenq.source15818-inellipse-exact-audit.v1",
        "frontier": str(args.frontier),
        "source_index": args.source_index,
        **report,
    }
    print(json.dumps(payload, indent=2))
    if not payload["all_certified"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

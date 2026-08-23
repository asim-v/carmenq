from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import numpy as np
import pytest


RESEARCH = Path(__file__).resolve().parents[1] / "scratch" / "d2_frontier"
sys.path.insert(0, str(RESEARCH))
SPEC = importlib.util.spec_from_file_location(
    "pairwise_inellipse_box_cover",
    RESEARCH / "pairwise_inellipse_box_cover.py",
)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def horwitz(alpha: float, beta: float, x: float, y: float) -> float:
    return (
        beta * beta * x * x
        + alpha * alpha * y * y
        + MODULE.cross_coefficient(alpha, beta) * x * y
        - 2.0 * beta * x
        - 2.0 * alpha * y
        + 1.0
    )


def test_quadratic_maximum_detects_vertices_edge_and_interior() -> None:
    assert MODULE.quadratic_maximum_on_unit_triangle(1, 1, 0, 0, 0) == 1

    # -(x-.2)^2-(y-.3)^2 has its unique maximum at an interior point.
    interior = MODULE.quadratic_maximum_on_unit_triangle(
        -1.0, -1.0, 0.0, 0.4, 0.6, -0.13
    )
    assert interior == pytest.approx(0.0, abs=1e-14)

    # On x+y=1 this equals -(x-.4)^2, while the other edges are smaller.
    edge = MODULE.quadratic_maximum_on_unit_triangle(
        -2.0, -1.0, -2.0, 2.8, 2.0, -1.16
    )
    assert edge == pytest.approx(0.0, abs=1e-14)


def test_tangent_error_contains_random_parameter_box_union() -> None:
    generator = np.random.default_rng(20260823)
    for _ in range(80):
        lower_alpha = generator.uniform(1.0, 1.8)
        lower_beta = generator.uniform(1.0, 1.8)
        upper_alpha = generator.uniform(lower_alpha, 2.0)
        upper_beta = generator.uniform(lower_beta, 2.0)
        anchor_alpha = generator.uniform(lower_alpha, upper_alpha)
        anchor_beta = generator.uniform(lower_beta, upper_beta)
        error = MODULE.tangent_coefficient_error_at(
            anchor_alpha,
            anchor_beta,
            (lower_alpha, upper_alpha),
            (lower_beta, upper_beta),
        )
        for _ in range(100):
            alpha = generator.uniform(lower_alpha, upper_alpha)
            beta = generator.uniform(lower_beta, upper_beta)
            x = generator.uniform()
            y = generator.uniform(0.0, 1.0 - x)
            difference = horwitz(anchor_alpha, anchor_beta, x, y) - horwitz(
                alpha, beta, x, y
            )
            assert difference <= error + 2e-13


def test_anchor_grid_includes_interior_grid_and_corners() -> None:
    anchors = MODULE.box_anchor_relaxations((1.1, 1.3), (1.2, 1.6))
    assert len(anchors) == 29
    locations = {(round(alpha, 12), round(beta, 12)) for alpha, beta, _ in anchors}
    assert (1.1, 1.2) in locations
    assert (1.3, 1.6) in locations
    assert all(error >= 0.0 for _, _, error in anchors)


def test_compact_center_polynomial_matches_reciprocal_chart() -> None:
    center_x, center_y = 0.31, 0.41
    center_r = 1.0 - center_x - center_y
    numerator = 1.0 - 2.0 * center_r
    alpha = 2.0 * center_x / numerator
    beta = 2.0 * center_y / numerator
    coefficients = MODULE.center_inellipse_coefficients(center_x, center_y)
    for x_value, y_value in ([0.12, 0.23], [0.4, 0.1], [0.05, 0.7]):
        centered = (
            coefficients[0] * x_value * x_value
            + coefficients[1] * y_value * y_value
            + coefficients[2] * x_value * y_value
            + coefficients[3] * x_value
            + coefficients[4] * y_value
            + coefficients[5]
        )
        reciprocal = horwitz(alpha, beta, x_value, y_value)
        assert centered == pytest.approx(
            numerator * numerator * reciprocal, abs=2e-14
        )


def test_center_box_error_is_one_sided() -> None:
    x_bounds = (0.22, 0.38)
    y_bounds = (0.28, 0.46)
    anchor_x, anchor_y = 1.0 / 3.0, 1.0 / 3.0
    error = MODULE.center_coefficient_error_at(
        anchor_x, anchor_y, x_bounds, y_bounds
    )
    anchor = np.asarray(
        MODULE.center_inellipse_coefficients(anchor_x, anchor_y)
    )
    generator = np.random.default_rng(91)
    for _ in range(500):
        true_x = generator.uniform(*x_bounds)
        true_y = generator.uniform(*y_bounds)
        if true_x + true_y < 0.5:
            continue
        point = generator.dirichlet(np.ones(3))[:2]
        true = np.asarray(MODULE.center_inellipse_coefficients(true_x, true_y))
        monomials = np.asarray(
            [
                point[0] ** 2,
                point[1] ** 2,
                point[0] * point[1],
                point[0],
                point[1],
                1.0,
            ]
        )
        assert np.dot(anchor - true, monomials) <= error + 2e-12

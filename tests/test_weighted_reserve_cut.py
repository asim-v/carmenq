from __future__ import annotations

import itertools
from fractions import Fraction
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.optimize import minimize_scalar


RESEARCH = Path(__file__).resolve().parents[1] / "scratch" / "d2_frontier"
sys.path.insert(0, str(RESEARCH))

from active_readout_audit_cap import projection_vertices, reserve  # noqa: E402
from four_active_socp_exact_cover import WeightBox, weight_hull  # noqa: E402
from four_active_mccormick_socp_exact_cover import (  # noqa: E402
    COMMON_BIAS_COEFFICIENT_REPRESENTATIVES,
    COMMON_BIAS_COEFFICIENTS,
    McCormickOracle,
    PERTURBATION,
    dominant_pair_loss_factor_upper,
    pair_loss_factor_upper,
    physical_weight_vertices,
    scalar_reserve_minimum_lower,
    weighted_reserve_cut_lower,
)


def numerical_scalar_minimum(x_value: float, a: float, b: float, c: float) -> float:
    def objective(bias: float) -> float:
        return float(a + b * bias + c * reserve(bias, x_value))

    interior = minimize_scalar(
        objective,
        bounds=(0.0, 1.0 - 1e-12),
        method="bounded",
        options={"xatol": 1e-14},
    )
    return min(objective(0.0), objective(1.0), float(interior.fun))


@pytest.mark.parametrize(
    "x_value,a,b,c",
    [
        (-1.0, 1.2, -0.4, 0.95),
        (-0.9872575864, 1.5, -0.5, 0.95),
        (-0.3, 1.1, 0.2, 1.05),
        (0.0, 0.8, -0.1, 1.0),
        (0.72, 1.4, 0.35, 1.05),
        (1.0, 0.9, -0.3, 0.95),
    ],
)
def test_scalar_rational_minimum_is_sharp_and_one_sided(
    x_value: float, a: float, b: float, c: float
) -> None:
    exact = scalar_reserve_minimum_lower(
        Fraction(str(x_value)),
        Fraction(str(a)),
        Fraction(str(b)),
        Fraction(str(c)),
    )
    numerical = numerical_scalar_minimum(x_value, a, b, c)
    assert float(exact) <= numerical + 2e-11
    assert float(exact) == pytest.approx(numerical, abs=3e-8)


def numerical_weighted_reserve(weights: np.ndarray, coefficients: np.ndarray) -> float:
    best = np.inf
    for point in projection_vertices(weights):
        def objective(bias: float) -> float:
            return float(np.dot(coefficients, reserve(bias, point)))

        interior = minimize_scalar(
            objective,
            bounds=(0.0, 1.0 - 1e-12),
            method="bounded",
            options={"xatol": 1e-14},
        )
        best = min(
            best,
            objective(0.0),
            objective(1.0),
            float(interior.fun),
        )
    return float(best)


def test_weighted_cut_matches_full_vertex_enumeration_at_fixed_weights() -> None:
    random = np.random.default_rng(20260826)
    checked = 0
    while checked < 12:
        weights = np.sort(2.0 * random.dirichlet(np.ones(4)))[::-1]
        if weights[0] > 1.0 or weights[0] <= 0.88325:
            continue
        exact_weights = tuple(Fraction(str(item)) for item in weights)
        hull = tuple((item, item) for item in exact_weights)
        for elevated, discounted in ((0, 1), (1, 0), (0, 3), (3, 1)):
            coefficients = tuple(
                1 + PERTURBATION
                if index == elevated
                else 1 - PERTURBATION
                if index == discounted
                else Fraction(1)
                for index in range(4)
            )
            certified = weighted_reserve_cut_lower(hull, coefficients)
            numerical = numerical_weighted_reserve(
                weights, np.asarray(list(map(float, coefficients)))
            )
            assert float(certified) <= numerical + 3e-10
            assert float(certified) == pytest.approx(numerical, abs=7e-8)
        checked += 1


def test_weight_box_enclosure_is_monotone_under_widening() -> None:
    weights = tuple(
        map(
            Fraction,
            ("0.994", "0.948", "0.041", "0.017"),
        )
    )
    narrow = tuple((item, item) for item in weights)
    radius = Fraction(1, 100000)
    wide = tuple((item - radius, item + radius) for item in weights)
    coefficients = (
        1 + PERTURBATION,
        1 - PERTURBATION,
        Fraction(1),
        Fraction(1),
    )
    assert weighted_reserve_cut_lower(wide, coefficients) <= (
        weighted_reserve_cut_lower(narrow, coefficients)
    )

def test_all_prior_common_bias_orbit_is_complete_and_one_sided() -> None:
    expected = {
        coefficients
        for representative in COMMON_BIAS_COEFFICIENT_REPRESENTATIVES
        for coefficients in itertools.permutations(representative)
    }
    assert len(COMMON_BIAS_COEFFICIENTS) == 48
    assert set(COMMON_BIAS_COEFFICIENTS) == expected
    weights = np.asarray([0.9983, 0.62354, 0.29649, 0.08167])
    exact_weights = tuple(Fraction(str(item)) for item in weights)
    hull = tuple((item, item) for item in exact_weights)
    for coefficients in COMMON_BIAS_COEFFICIENTS[::5]:
        certified = weighted_reserve_cut_lower(hull, coefficients)
        numerical = numerical_weighted_reserve(
            weights, np.asarray(list(map(float, coefficients)))
        )
        assert float(certified) <= numerical + 3e-10
        assert float(certified) == pytest.approx(numerical, abs=7e-8)


def test_common_bias_cut_encloses_interior_weights_of_a_correlated_box() -> None:
    box = WeightBox(
        (
            (Fraction("0.91"), Fraction("0.99")),
            (Fraction("0.42"), Fraction("0.67")),
            (Fraction("0.12"), Fraction("0.34")),
        )
    )
    hull = weight_hull(box)
    assert hull is not None
    vertices = physical_weight_vertices(box)
    coefficients = COMMON_BIAS_COEFFICIENT_REPRESENTATIVES[2]
    certified = weighted_reserve_cut_lower(hull, coefficients, vertices)
    random = np.random.default_rng(76138331)
    for _ in range(24):
        mixture = random.dirichlet(np.ones(len(vertices)))
        weights = sum(
            (amount * np.asarray(list(map(float, vertex)))
             for amount, vertex in zip(mixture, vertices)),
            start=np.zeros(4),
        )
        assert weights.sum() == pytest.approx(2.0, abs=2e-14)
        assert np.all(weights[:-1] >= weights[1:] - 2e-14)
        numerical = numerical_weighted_reserve(
            weights, np.asarray(list(map(float, coefficients)))
        )
        assert float(certified) <= numerical + 3e-10


def test_oracle_emits_the_full_common_bias_support_orbit() -> None:
    box = WeightBox(
        ((Fraction("0.99"), Fraction("0.995")),
         (Fraction("0.008"), Fraction("0.012")),
         (Fraction("0.002"), Fraction("0.004")))
    )
    hull = weight_hull(box)
    assert hull is not None
    enclosure = McCormickOracle().assign(hull, (1, 2, 3), physical_weight_vertices(box))
    assert len(enclosure["common_bias_support_cuts"]) == 48

def test_dominant_pair_factor_is_the_exact_triangle_expression() -> None:
    weights = tuple(map(Fraction, ("0.99", "0.98", "0.02", "0.01")))
    hull = tuple((item, item) for item in weights)
    factor = dominant_pair_loss_factor_upper(hull)
    transverse_triangle = (
        (weights[2] + weights[3]) ** 2 - (weights[0] - weights[1]) ** 2
    ) / (4 * weights[0] * weights[1])
    assert factor == transverse_triangle == Fraction(1, 4851)


def test_dominant_pair_factor_uses_a_safe_box_upper() -> None:
    box = WeightBox(
        (
            (Fraction("0.98"), Fraction("0.995")),
            (Fraction("0.008"), Fraction("0.018")),
            (Fraction("0.002"), Fraction("0.007")),
        )
    )
    hull = weight_hull(box)
    assert hull is not None
    upper = dominant_pair_loss_factor_upper(hull)
    for weights in physical_weight_vertices(box):
        exact = (1 - weights[0]) * (1 - weights[1]) / (
            weights[0] * weights[1]
        )
        assert exact <= upper


def test_pair_factor_caps_noninformative_pairs_at_one() -> None:
    weights = tuple(map(Fraction, ("0.99", "0.51", "0.30", "0.20")))
    hull = tuple((item, item) for item in weights)
    assert pair_loss_factor_upper(hull, 0, 1) == Fraction(49, 5049)
    assert pair_loss_factor_upper(hull, 2, 3) == 1


def test_every_pair_factor_is_safe_on_exact_weight_vertices() -> None:
    box = WeightBox(
        ((Fraction("0.94"), Fraction("0.99")),
         (Fraction("0.25"), Fraction("0.45")),
         (Fraction("0.05"), Fraction("0.25")))
    )
    hull = weight_hull(box)
    assert hull is not None
    for first, second in ((i, j) for i in range(4) for j in range(4) if i != j):
        upper = pair_loss_factor_upper(hull, first, second)
        for weights in physical_weight_vertices(box):
            exact = min(Fraction(1), (1 - weights[first]) * (1 - weights[second]) / (weights[first] * weights[second]))
            assert exact <= upper


def test_projective_pair_rows_cover_all_orientations_and_lines() -> None:
    box = WeightBox(
        (
            (Fraction("0.99"), Fraction("0.995")),
            (Fraction("0.008"), Fraction("0.012")),
            (Fraction("0.002"), Fraction("0.004")),
        )
    )
    hull = weight_hull(box)
    assert hull is not None
    permutation = (2, 1, 3)
    enclosure = McCormickOracle().assign(
        hull, permutation, physical_weight_vertices(box)
    )
    cuts = enclosure["full_closure_projective_pair_cuts"]
    assert len(cuts) == 24
    assert {(cut["retained_effect"], cut["other_effect"]) for cut in cuts} == {
        (i, j) for i in range(4) for j in range(4) if i != j
    }
    assert {tuple(cut["line_weight"]) for cut in cuts} == {(11, 20), (3, 5)}
    position = {effect: syndrome for syndrome, effect in enumerate((0, *permutation))}
    assert all(cut["other_syndrome"] == position[cut["other_effect"]] for cut in cuts)

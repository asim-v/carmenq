"""Proof-kernel checks for the rank/rank tangent interval certificate."""

from __future__ import annotations

from fractions import Fraction
import math
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
FRONTIER = ROOT / "scratch" / "d2_frontier"
if str(FRONTIER) not in sys.path:
    sys.path.insert(0, str(FRONTIER))

from projective_tangent_interval_certificate import (  # noqa: E402
    Interval,
    RankRankCertifier,
    RankRankBox,
    decimal_box,
    root_from_bounds,
)
from projective_low_eigenvalue_face import LowEigenvalueFaceCertifier  # noqa: E402
from carmenq.projective_secular import (  # noqa: E402
    endpoint_split_state_term,
    rank_split_state_term,
)



def test_exact_decimal_constants_are_outward_enclosed() -> None:
    for text in ("0", "0.03", "0.6", "0.76662", "1"):
        interval = Interval.decimal(text)
        exact = Fraction(text)
        assert Fraction.from_float(interval.lo) <= exact
        assert exact <= Fraction.from_float(interval.hi)


def test_tangent_choice_is_not_part_of_the_trusted_bound() -> None:
    p_value = 0.9523404308856647
    certifier = RankRankCertifier(
        Fraction(3, 5),
        Fraction(38331, 50000),
        Interval.exact_float(p_value),
    )
    box = root_from_bounds(
        Interval.exact_float(p_value),
        Interval.exact_float(1.0 - p_value),
        Interval.exact_float(0.022341589656722696),
        Interval.exact_float(0.022341589656722696),
    )
    selected = certifier.bound(box)
    deliberately_bad = certifier.bound(box, tangents=(0.4, 1.4, 0.4, 1.4))

    assert selected.upper < 1.0
    assert deliberately_bad.upper >= selected.upper
    assert deliberately_bad.upper < float("inf")


def test_small_critical_neighborhood_closes_without_branching() -> None:
    p_value = 0.9523404308856647
    radius = 2e-5
    x_bounds = decimal_box(str(p_value - radius), str(p_value + radius))
    y_bounds = decimal_box(
        str(1.0 - p_value - radius), str(1.0 - p_value + radius)
    )
    angle = 0.022341589656722696
    angle_bounds = decimal_box(str(angle - radius), str(angle + radius))
    certifier = RankRankCertifier(
        Fraction(3, 5), Fraction(38331, 50000), x_bounds
    )
    payload = certifier.certify(
        root_from_bounds(x_bounds, y_bounds, angle_bounds, angle_bounds),
        max_boxes=0,
    )

    assert payload["complete"] is True
    assert payload["boxes_split"] == 0
    assert payload["maximum_open_upper"] == float("-inf")


def test_rank_rank_box_round_trips_through_json_shape() -> None:
    box = RankRankBox(
        decimal_box("0", "0.01"),
        decimal_box("0.02", "0.03"),
        decimal_box("0.04", "0.05"),
        decimal_box("0.06", "0.07"),
    )
    assert RankRankBox.deserialise(box.serialise()) == box


def test_incomplete_cell_resumes_from_its_saved_open_frontier() -> None:
    certifier = LowEigenvalueFaceCertifier(
        Fraction(11, 20),
        Fraction(7573, 10000),
        decimal_box("0.92", "0.94"),
    )
    root = root_from_bounds(
        decimal_box("0.92", "0.94"),
        decimal_box("0", "0.02"),
        decimal_box("0", "0.1"),
        decimal_box("0", "0.1"),
    )
    first = certifier.certify(root, max_boxes=1)
    assert first["complete"] is False
    assert first["boxes_split"] == 1
    assert first["open_frontier"]
    resumed = certifier.certify(root, max_boxes=2, resume=first)
    assert resumed["boxes_split"] == 2
    assert (
        resumed["boxes_closed"] + resumed["boxes_remaining"]
        == resumed["boxes_split"] + 1
    )


def test_all_topologies_bound_a_dense_state_grid() -> None:
    weight = 0.6
    level = 0.76662
    x_value = 0.9
    y_value = 0.05
    first_sine = 0.08
    second_sine = 0.14
    box = root_from_bounds(
        Interval.exact_float(x_value),
        Interval.exact_float(y_value),
        Interval.exact_float(first_sine),
        Interval.exact_float(second_sine),
    )
    states = np.linspace(0.0, 1.0, 1025)
    for first_kind in ("endpoint", "rank"):
        for second_kind in ("endpoint", "rank"):
            certifier = RankRankCertifier(
                Fraction(3, 5),
                Fraction(38331, 50000),
                Interval.exact_float(x_value),
                first_kind,
                second_kind,
            )
            upper = certifier.bound(box).upper
            exact_terms = []
            for kind, high, low, sine in (
                (first_kind, x_value, y_value, first_sine),
                (second_kind, 1.0 - y_value, 1.0 - x_value, second_sine),
            ):
                for label in (0, 1):
                    if kind == "rank":
                        values = (
                            rank_split_state_term(
                                high, low, sine, float(state), label,
                                weight, level,
                            )
                            for state in states
                        )
                    else:
                        values = (
                            endpoint_split_state_term(
                                high, low, float(state), label, weight, level
                            )
                            for state in states
                        )
                    exact_terms.append(max(values))
            sampled = (1.0 - weight) * math.fsum(exact_terms) / 8.0
            assert sampled <= upper + 3e-11


def test_zero_face_lipschitz_bound_contains_nearby_exact_points() -> None:
    certifier = RankRankCertifier(
        Fraction(3, 5),
        Fraction(76591, 100000),
        decimal_box("0.945", "0.955"),
    )
    box = RankRankBox(
        decimal_box("0.045", "0.055"),
        decimal_box("0.0025", "0.005"),
        decimal_box("0.03", "0.045"),
        decimal_box("0.0001", "0.0004"),
    )
    upper = certifier._zero_face_lipschitz_upper(box)
    assert math.isfinite(upper)
    zero_box = RankRankBox(box.y, box.residual, box.first_sine, Interval(0.0, 0.0))
    tangents = certifier._select_tangents(zero_box)
    random = np.random.default_rng(20260826)
    for _ in range(128):
        point = RankRankBox(
            *(
                Interval.exact_float(random.uniform(item.lo, item.hi))
                for item in box.coordinates
            )
        )
        exact_interval, _ = certifier._value_expression(point, tangents)
        assert exact_interval.hi <= upper + 2e-12


def test_near_zero_face_is_used_even_when_interval_excludes_zero() -> None:
    certifier = RankRankCertifier(
        Fraction(3, 5),
        Fraction(76591, 100000),
        decimal_box("0.945", "0.955"),
    )
    touching = RankRankBox(
        decimal_box("0.0525", "0.055"),
        decimal_box("0.0025", "0.00375"),
        decimal_box("0.03", "0.03375"),
        decimal_box("0", "0.000234375"),
    )
    separated = RankRankBox(
        touching.y,
        touching.residual,
        touching.first_sine,
        decimal_box("0.000000000001", "0.000000000002"),
    )
    assert certifier._zero_face_lipschitz_upper(touching) < 1.0
    assert math.isfinite(certifier._zero_face_lipschitz_upper(separated))
    assert certifier.bound(separated).upper < 1.0


def test_zero_face_split_does_not_chase_an_infinitesimal_angle() -> None:
    certifier = RankRankCertifier(
        Fraction(3, 5),
        Fraction(76591, 100000),
        decimal_box("0.94", "0.96"),
    )
    root = RankRankBox(
        decimal_box("0.04", "0.06"),
        decimal_box("0", "0.02"),
        decimal_box("0.03", "0.06"),
        decimal_box("0", "0.03"),
    )
    infinitesimal = RankRankBox(
        root.y,
        root.residual,
        root.first_sine,
        Interval(0.0, math.ldexp(1.0, -500)),
    )
    scales = tuple(max(item.width, 1e-15) for item in root.coordinates)
    coordinate, _, _ = certifier.split(infinitesimal, scales)
    assert coordinate in {"y", "residual", "first_sine"}


def test_first_zero_face_lipschitz_bound_contains_nearby_exact_points() -> None:
    certifier = RankRankCertifier(
        Fraction(3, 5),
        Fraction(76591, 100000),
        decimal_box("0.945", "0.955"),
    )
    box = RankRankBox(
        decimal_box("0.045", "0.055"),
        decimal_box("0.0025", "0.005"),
        decimal_box("0.0001", "0.0004"),
        decimal_box("0.03", "0.045"),
    )
    upper = certifier._zero_face_lipschitz_upper(box, 0)
    assert math.isfinite(upper)
    zero_box = RankRankBox(
        box.y, box.residual, Interval(0.0, 0.0), box.second_sine
    )
    tangents = certifier._select_tangents(zero_box)
    random = np.random.default_rng(20260827)
    for _ in range(128):
        point = RankRankBox(
            *(
                Interval.exact_float(random.uniform(item.lo, item.hi))
                for item in box.coordinates
            )
        )
        exact_interval, _ = certifier._value_expression(point, tangents)
        assert exact_interval.hi <= upper + 2e-12


def test_first_zero_face_split_does_not_chase_an_infinitesimal_angle() -> None:
    certifier = RankRankCertifier(
        Fraction(3, 5),
        Fraction(76591, 100000),
        decimal_box("0.94", "0.96"),
    )
    root = RankRankBox(
        decimal_box("0.04", "0.06"),
        decimal_box("0", "0.02"),
        decimal_box("0", "0.03"),
        decimal_box("0.03", "0.06"),
    )
    infinitesimal = RankRankBox(
        root.y,
        root.residual,
        Interval(0.0, math.ldexp(1.0, -500)),
        root.second_sine,
    )
    scales = tuple(max(item.width, 1e-15) for item in root.coordinates)
    coordinate, _, _ = certifier.split(infinitesimal, scales)
    assert coordinate in {"y", "residual", "second_sine"}


def test_low_eigenvalue_hessian_split_targets_the_taylor_remainder() -> None:
    certifier = LowEigenvalueFaceCertifier(
        Fraction(11, 20),
        Fraction(7573, 10000),
        decimal_box("0.92", "0.94"),
    )
    leader = RankRankBox(
        decimal_box("0.005", "0.01"),
        decimal_box("0.05", "0.065"),
        decimal_box("0.075", "0.0765625"),
        decimal_box("0.05", "0.1"),
    )
    coordinate, _, _ = certifier.split(
        leader, (0.02, 0.04, 0.1, 0.1)
    )
    assert coordinate == "residual"


def test_center_matrix_lipschitz_bound_is_regular_at_crossings() -> None:
    certifier = RankRankCertifier(
        Fraction(3, 5),
        Fraction(76591, 100000),
        decimal_box("0.98", "1"),
    )
    box = root_from_bounds(
        decimal_box("0.98", "1"),
        decimal_box("0", "0.02"),
        decimal_box("0", "0.1"),
        decimal_box("0", "0.1"),
    )
    contracted = certifier.contract(box)
    assert contracted is not None
    tangents = certifier._select_tangents(contracted)
    upper = certifier._center_matrix_lipschitz_sum_upper(contracted, tangents)
    assert math.isfinite(upper)


def test_center_matrix_lipschitz_bound_contains_random_points() -> None:
    certifier = RankRankCertifier(
        Fraction(3, 5),
        Fraction(76591, 100000),
        decimal_box("0.98", "1"),
    )
    box = root_from_bounds(
        decimal_box("0.98", "1"),
        decimal_box("0", "0.02"),
        decimal_box("0", "0.1"),
        decimal_box("0", "0.1"),
    )
    contracted = certifier.contract(box)
    assert contracted is not None
    tangents = certifier._select_tangents(contracted)
    upper = certifier._center_matrix_lipschitz_sum_upper(contracted, tangents)
    random = np.random.default_rng(20260828)
    for _ in range(256):
        point = RankRankBox(
            *(
                Interval.exact_float(random.uniform(item.lo, item.hi))
                for item in contracted.coordinates
            )
        )
        exact_interval, _ = certifier._value_expression(point, tangents)
        assert exact_interval.hi <= upper + 2e-12


def test_low_eigenvalue_face_bound_contains_random_points() -> None:
    certifier = LowEigenvalueFaceCertifier(
        Fraction(11, 20),
        Fraction(761, 1000),
        decimal_box("0.85", "0.9"),
    )
    random = np.random.default_rng(20260829)
    for y_lower, y_upper in (("0", "0.00625"), ("0.0125", "0.025")):
        box = RankRankBox(
            decimal_box(y_lower, y_upper),
            decimal_box("0.125", "0.15"),
            decimal_box("0.08125", "0.084375"),
            decimal_box("0.1", "0.15"),
        )
        upper = certifier._low_eigenvalue_face_lipschitz_upper(box)
        assert math.isfinite(upper)
        face_box = RankRankBox(
            Interval.exact_float(0.0),
            box.residual,
            box.first_sine,
            box.second_sine,
        )
        tangents = certifier._select_tangents(face_box)
        for _ in range(128):
            point = RankRankBox(
                *(
                    Interval.exact_float(random.uniform(item.lo, item.hi))
                    for item in box.coordinates
                )
            )
            exact_interval, _ = certifier._value_expression(point, tangents)
            assert exact_interval.hi <= upper + 2e-12


def test_sqrt_y_mean_value_bound_contains_random_points() -> None:
    certifier = LowEigenvalueFaceCertifier(
        Fraction(11, 20),
        Fraction(7573, 10000),
        decimal_box("0.8", "0.9"),
    )
    box = RankRankBox(
        decimal_box("0.001", "0.00625"),
        decimal_box("0.09375", "0.125"),
        decimal_box("0.053125", "0.0546875"),
        decimal_box("0.05", "0.0625"),
    )
    face_box = RankRankBox(
        Interval.exact_float(0.0),
        box.residual,
        box.first_sine,
        box.second_sine,
    )
    tangents = certifier._select_tangents(face_box)
    upper = certifier._sqrt_y_mean_value_upper(box, tangents)
    assert math.isfinite(upper)
    random = np.random.default_rng(20260830)
    for _ in range(256):
        point = RankRankBox(
            *(
                Interval.exact_float(random.uniform(item.lo, item.hi))
                for item in box.coordinates
            )
        )
        exact_interval, _ = certifier._value_expression(point, tangents)
        assert exact_interval.hi <= upper + 2e-12

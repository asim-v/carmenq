from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest


RESEARCH = Path(__file__).resolve().parents[1] / "scratch" / "d2_frontier"
SPEC = importlib.util.spec_from_file_location(
    "active_readout_audit_cap", RESEARCH / "active_readout_audit_cap.py"
)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_reserve_endpoint_and_unbiased_dual() -> None:
    assert MODULE.reserve(0.0, -0.37) == pytest.approx(0.5)
    assert MODULE.reserve(1.0, -1.0) == pytest.approx(0.0)
    assert MODULE.reserve(1.0, 0.25) == pytest.approx(1.0)


def test_projection_vertices_obey_weighted_closure() -> None:
    weights = np.asarray([0.9, 0.5, 0.4, 0.2])
    vertices = MODULE.projection_vertices(weights)
    assert vertices
    for point in vertices:
        assert np.max(np.abs(point)) <= 1.0 + 1e-12
        assert np.dot(weights, point) == pytest.approx(0.0, abs=2e-12)


def test_balanced_four_active_readout_has_half_cap() -> None:
    result = MODULE.active_audit_cap(np.full(4, 0.5))
    assert result["minimum_total_prior_reserve"] == pytest.approx(2.0, abs=2e-10)
    assert result["audit_upper"] == pytest.approx(0.5, abs=5e-11)


def test_trine_cap_recovers_two_thirds() -> None:
    result = MODULE.active_audit_cap(np.full(3, 2.0 / 3.0))
    assert result["audit_upper"] == pytest.approx(2.0 / 3.0, abs=5e-10)


def test_asymmetric_four_active_cap_is_nontrivial() -> None:
    result = MODULE.active_audit_cap(np.asarray([0.9, 0.5, 0.4, 0.2]))
    assert result["audit_upper"] < 0.59
    assert result["projection_residual"] == pytest.approx(0.0, abs=2e-12)


@pytest.mark.parametrize(
    "weights",
    [
        [0.5, 0.5, 0.5, 0.5],
        [0.9, 0.5, 0.4, 0.2],
        [0.95, 0.5, 0.4, 0.15],
        [0.99, 0.5, 0.49, 0.02],
        [0.9994, 0.50015, 0.50015, 0.0003],
    ],
)
def test_closed_four_active_cap_matches_vertex_enumeration(weights) -> None:
    value = np.asarray(weights)
    closed = MODULE.closed_four_active_audit_cap(value)
    enumerated = MODULE.active_audit_cap(value)
    assert closed["minimum_total_prior_reserve"] == pytest.approx(
        enumerated["minimum_total_prior_reserve"], abs=2e-9
    )
    assert closed["audit_upper"] == pytest.approx(
        enumerated["audit_upper"], abs=2e-9
    )
    assert closed["projection_residual"] == pytest.approx(0.0, abs=2e-12)


def test_closed_four_active_cap_matches_random_sorted_weights() -> None:
    random = np.random.default_rng(20260826)
    checked = 0
    while checked < 100:
        weights = np.sort(2.0 * random.dirichlet(np.ones(4)))[::-1]
        if weights[0] > 1.0:
            continue
        closed = MODULE.closed_four_active_audit_cap(weights)
        enumerated = MODULE.active_audit_cap(weights)
        assert closed["audit_upper"] == pytest.approx(
            enumerated["audit_upper"], abs=3e-8
        )
        checked += 1


def enumerated_complement_reserve(
    weights: np.ndarray, excluded: int = 0
) -> float:
    coefficients = np.ones(4)
    coefficients[excluded] = 0.0
    best = np.inf
    for point in MODULE.projection_vertices(weights):
        def objective(bias: float) -> float:
            return float(np.dot(coefficients, MODULE.reserve(bias, point)))

        interior = MODULE.minimize_scalar(
            objective,
            bounds=(0.0, 1.0 - 1e-12),
            method="bounded",
            options={"xatol": 2e-14},
        )
        best = min(
            best,
            objective(0.0),
            objective(1.0),
            float(interior.fun),
        )
    return float(best)


@pytest.mark.parametrize(
    "weights",
    [
        [0.9, 0.7, 0.3, 0.1],
        [0.95, 0.5, 0.4, 0.15],
        [0.99, 0.5, 0.49, 0.02],
        [0.994, 0.948, 0.041, 0.017],
        [0.9994, 0.50015, 0.50015, 0.0003],
    ],
)
def test_closed_complement_reserve_matches_vertex_enumeration(weights) -> None:
    value = np.asarray(weights)
    closed = MODULE.closed_four_active_complement_reserve(value)
    assert closed["minimum_complement_prior_reserve"] == pytest.approx(
        enumerated_complement_reserve(value), abs=3e-8
    )
    assert closed["projection_residual"] == pytest.approx(0.0, abs=2e-12)


def test_closed_complement_reserve_matches_random_near_projective_weights() -> None:
    random = np.random.default_rng(20260826)
    checked = 0
    while checked < 40:
        weights = np.sort(2.0 * random.dirichlet(np.ones(4)))[::-1]
        if weights[0] > 1.0 or weights[0] <= 0.88325:
            continue
        closed = MODULE.closed_four_active_complement_reserve(weights)
        assert closed["minimum_complement_prior_reserve"] == pytest.approx(
            enumerated_complement_reserve(weights), abs=5e-8
        )
        checked += 1


@pytest.mark.parametrize(
    "weights",
    [
        [0.9, 0.7, 0.3, 0.1],
        [0.95, 0.5, 0.4, 0.15],
        [0.99, 0.5, 0.49, 0.02],
        [0.994, 0.948, 0.041, 0.017],
        [0.9994, 0.50015, 0.50015, 0.0003],
    ],
)
def test_all_four_complement_reserves_match_enumeration(weights) -> None:
    value = np.sort(np.asarray(weights))[::-1]
    for excluded in range(4):
        closed = MODULE.closed_four_active_complement_reserve(
            value, excluded=excluded
        )
        enumerated = enumerated_complement_reserve(value, excluded)
        certified = closed["minimum_complement_prior_reserve"]
        assert certified <= enumerated + 5e-8
        if closed["fractional_fill"] < 0.5:
            assert certified == pytest.approx(enumerated, abs=5e-8)


@pytest.mark.parametrize(
    "weights",
    [
        [0.9, 0.7, 0.3, 0.1],
        [0.95, 0.5, 0.4, 0.15],
        [0.99, 0.5, 0.49, 0.02],
        [0.994, 0.948, 0.041, 0.017],
        [0.9994, 0.50015, 0.50015, 0.0003],
    ],
)
def test_all_six_pair_reserves_match_vertex_enumeration(weights) -> None:
    value = np.sort(np.asarray(weights))[::-1]
    for first in range(4):
        for second in range(first + 1, 4):
            pair = (first, second)
            coefficients = np.zeros(4)
            coefficients[list(pair)] = 1.0
            best = np.inf
            for point in MODULE.projection_vertices(value):

                def objective(bias: float) -> float:
                    return float(
                        np.dot(coefficients, MODULE.reserve(bias, point))
                    )

                interior = MODULE.minimize_scalar(
                    objective, bounds=(0.0, 1.0 - 1e-12), method="bounded"
                )
                best = min(best, objective(0.0), objective(1.0), interior.fun)
            closed = MODULE.closed_four_active_pair_reserve(value, pair)
            assert closed["minimum_pair_prior_reserve"] == pytest.approx(
                best, abs=5e-8
            )
            assert closed["projection_residual"] == pytest.approx(0.0, abs=3e-12)

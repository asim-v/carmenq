"""Checks for the analytic projective-state elimination."""

from __future__ import annotations

import math

import numpy as np
from scipy.optimize import minimize_scalar

from carmenq.projective_secular import (
    endpoint_split_state_term,
    endpoint_split_tangent_upper,
    hellinger_tangent,
    rank_split_state_term,
    rank_split_tangent_upper,
    tangent_for_probability,
)


def test_hellinger_tangent_is_global_and_touches() -> None:
    for probability in np.linspace(0.01, 0.99, 33):
        tangent = tangent_for_probability(float(probability))
        exact, affine, gap = hellinger_tangent(float(probability), tangent)
        assert math.isclose(exact, affine, rel_tol=0.0, abs_tol=2e-14)
        assert gap >= -2e-14
        for trial in np.linspace(0.0, 1.0, 101):
            _, _, trial_gap = hellinger_tangent(float(trial), tangent)
            assert trial_gap >= -2e-14


def test_generalized_eigenvalue_bounds_every_state() -> None:
    rng = np.random.default_rng(20260826)
    for _ in range(40):
        low = float(rng.uniform(0.01, 0.35))
        high = float(rng.uniform(max(0.55, low), 0.99))
        sine = float(rng.uniform(0.0, 0.4))
        weight = float(rng.uniform(0.3, 0.7))
        level = float(rng.uniform(weight * high + 0.03, 0.99))
        tangent = float(math.exp(rng.uniform(-1.5, 1.5)))
        for label in (0, 1):
            upper = rank_split_tangent_upper(
                high, low, sine, label, tangent, weight, level
            )
            sampled = max(
                rank_split_state_term(
                    high, low, sine, float(state), label, weight, level
                )
                for state in np.linspace(0.0, 1.0, 513)
            )
            assert sampled <= upper + 3e-12


def test_endpoint_tangent_bounds_every_state() -> None:
    rng = np.random.default_rng(20260827)
    for _ in range(40):
        low = float(rng.uniform(0.0, 0.4))
        high = float(rng.uniform(max(0.5, low), 0.99))
        weight = float(rng.uniform(0.3, 0.7))
        level = float(rng.uniform(weight * high + 0.03, 0.99))
        tangent = float(math.exp(rng.uniform(-1.5, 1.5)))
        for label in (0, 1):
            upper = endpoint_split_tangent_upper(
                high, low, label, tangent, weight, level
            )
            sampled = max(
                endpoint_split_state_term(
                    high, low, float(state), label, weight, level
                )
                for state in np.linspace(0.0, 1.0, 513)
            )
            assert sampled <= upper + 3e-12


def test_optimal_probability_makes_the_tangent_tight() -> None:
    cases = (
        (0.9523404308856647, 1.0 - 0.9523404308856647, 0.022341589656722696),
        (0.94, 0.04, 0.06),
        (0.90, 0.08, 0.15),
    )
    weight = 0.6
    level = 0.76662
    for high, low, sine in cases:
        for label in (0, 1):
            optimum = minimize_scalar(
                lambda state: -rank_split_state_term(
                    high, low, sine, float(state), label, weight, level
                ),
                bounds=(0.0, 1.0),
                method="bounded",
                options={"xatol": 1e-14},
            )
            state = float(optimum.x)
            probability = low + (high - low) * state
            tangent = tangent_for_probability(probability)
            upper = rank_split_tangent_upper(
                high, low, sine, label, tangent, weight, level
            )
            assert upper >= -float(optimum.fun) - 2e-12
            assert upper <= -float(optimum.fun) + 2e-9

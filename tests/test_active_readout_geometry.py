from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pytest


RESEARCH = Path(__file__).resolve().parents[1] / "scratch" / "d2_frontier"
sys.path.insert(0, str(RESEARCH))

from active_readout_geometry_probe import (  # noqa: E402
    aligned_projective_audit,
    averaged_projective_audit,
    bloch_vectors,
    diagonal_interval,
    pair_components,
    projections,
)
from four_active_geometry_support_probe import (  # noqa: E402
    support_for_reserve,
    support_for_reserve_fast,
    water_filled_priors,
)
from four_active_spectral_upper import named_upper_matrices  # noqa: E402
from four_active_exact_projective_probe import weight_chart  # noqa: E402


def test_closed_quadrilateral_projection_parameterisation() -> None:
    weights = np.asarray([0.9, 0.5, 0.4, 0.2])
    assert diagonal_interval(weights) == pytest.approx((0.4, 0.6))
    variables = np.asarray([0.51, 0.27, 0.41, 1.19, 0.63])
    point = projections(weights, variables)
    assert np.max(np.abs(point)) <= 1.0 + 2e-12
    assert np.dot(weights, point) == pytest.approx(0.0, abs=2e-12)
    vectors = bloch_vectors(weights, variables)
    assert np.linalg.norm(vectors, axis=1) == pytest.approx(np.ones(4))
    assert np.linalg.norm((weights[:, None] * vectors).sum(axis=0)) < 2e-12


def test_pair_components_reconstruct_side_lengths() -> None:
    first, second, diagonal = 0.9, 0.5, 0.55
    a0, a1, transverse = pair_components(first, second, diagonal)
    assert a0 + a1 == pytest.approx(diagonal)
    assert a0 * a0 + transverse * transverse == pytest.approx(first * first)
    assert a1 * a1 + transverse * transverse == pytest.approx(second * second)


def test_exact_aligned_projective_audit_identity() -> None:
    rng = np.random.default_rng(20260823)
    for _ in range(100):
        audit = rng.uniform(0.2, 0.9)
        prior = rng.uniform(0.05, audit)
        first = rng.normal(size=3)
        first /= np.linalg.norm(first)
        second = rng.normal(size=3)
        second /= np.linalg.norm(second)
        dual_axis = rng.normal(size=3)
        dual_axis /= np.linalg.norm(dual_axis)
        bias = rng.uniform(0.0, 1.0)
        identity = np.eye(2, dtype=complex)
        pauli = np.asarray(
            [
                [[0, 1], [1, 0]],
                [[0, -1j], [1j, 0]],
                [[1, 0], [0, -1]],
            ],
            dtype=complex,
        )
        def projector(vector: np.ndarray) -> np.ndarray:
            return (identity + np.tensordot(vector, pauli, axes=1)) / 2.0

        y_dual = audit * (
            identity + bias * np.tensordot(dual_axis, pauli, axes=1)
        ) / 2.0
        rho_second = y_dual - (audit - prior) * (identity - projector(second))
        direct = (
            np.trace(projector(first) @ y_dual).real
            + prior
            - np.trace(projector(first) @ rho_second).real
        )
        overlap = float(first @ second)
        assert direct == pytest.approx(
            aligned_projective_audit(audit, prior, overlap), abs=2e-12
        )


def test_weighted_projective_average_eliminates_pair_geometry() -> None:
    weights = np.asarray([0.9, 0.5, 0.4, 0.2])
    variables = np.asarray([0.51, 0.27, 0.41, 1.19, 0.63])
    vectors = bloch_vectors(weights, variables)
    audit = 0.72
    priors = np.asarray([0.31, 0.27, 0.23, 0.19])
    for complement in range(4):
        retained = [index for index in range(4) if index != complement]
        direct = sum(
            weights[index]
            * aligned_projective_audit(
                audit,
                priors[complement],
                float(vectors[index] @ vectors[complement]),
            )
            for index in retained
        ) / (2.0 - weights[complement])
        assert direct == pytest.approx(
            averaged_projective_audit(
                audit, priors[complement], weights[complement]
            ),
            abs=2e-12,
        )


def test_water_filling_is_normalised_and_dominates_floors() -> None:
    lower = np.asarray([0.4, 0.2, 0.1, 0.0])
    priors = water_filled_priors(lower)
    assert priors == pytest.approx([0.4, 0.2, 0.2, 0.2])
    assert priors.sum() == pytest.approx(1.0)
    assert np.all(priors >= lower - 1e-14)


def test_uniform_reserve_joint_support() -> None:
    score, audit, priors, returned = support_for_reserve(np.full(4, 0.5), 0.6)
    assert audit == pytest.approx(0.5, abs=5e-8)
    assert priors == pytest.approx(np.full(4, 0.25), abs=5e-8)
    assert returned == pytest.approx(1.0, abs=5e-8)
    assert score == pytest.approx(0.7, abs=5e-8)


def test_closed_form_water_filling_matches_scalar_solver() -> None:
    rng = np.random.default_rng(424242)
    for _ in range(100):
        reserve_vector = rng.uniform(0.01, 1.0, size=4)
        slow = support_for_reserve(reserve_vector, 0.6)
        fast = support_for_reserve_fast(reserve_vector, 0.6)
        assert fast[0] == pytest.approx(slow[0], abs=2e-8)
        assert fast[1] == pytest.approx(slow[1], abs=2e-7)
        assert fast[3] == pytest.approx(slow[3], abs=2e-7)


def test_four_effect_matrix_family_contains_all_comparisons() -> None:
    weights = np.asarray([0.9, 0.5, 0.4, 0.2])
    named = named_upper_matrices(weights, 0.59, 0.6, 0.76591)
    assert len(named) == 15
    assert [name for name, _ in named[:3]] == [
        "prefix",
        "syndrome",
        "active_prior_reserve",
    ]
    assert all(matrix.shape == (16, 16) for _, matrix in named)


def test_four_active_weight_chart_is_sorted_and_complete() -> None:
    for maximum in (0.88325, 0.94, 1.0):
        for smallest_fraction in (0.0, 0.37, 1.0):
            for middle_fraction in (0.0, 0.61, 1.0):
                weights = weight_chart(
                    maximum, smallest_fraction, middle_fraction, 1e-6
                )
                assert weights.sum() == pytest.approx(2.0)
                assert weights[0] == pytest.approx(maximum)
                assert np.all(weights[:-1] >= weights[1:] - 2e-14)
                assert np.all(weights > 0.0)
                assert np.all(weights <= 1.0 + 2e-14)

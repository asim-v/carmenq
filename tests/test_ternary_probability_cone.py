from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scratch" / "d2_frontier"))

from pairwise_inellipse_box_cover import cross_coefficient  # noqa: E402
from fourier_behavior_cap_cover import cube_face_caps  # noqa: E402
from ternary_probability_cone_cover import (  # noqa: E402
    TernaryConeOracle,
    projective_comparison_bonus,
    terminal_weight_intervals,
    terminal_weights,
)
from terminal_reconstruction_enclosure import (  # noqa: E402
    planar_reconstruction,
    reconstruction_anchor_and_errors,
)
from validate_terminal_reconstructed_frontier import validate  # noqa: E402


def canonical_three_effect_povm(weights: np.ndarray) -> np.ndarray:
    """Construct the planar rank-one POVM without importing the Torch probes."""

    w0, w1, w2 = map(float, weights)
    cosine = np.clip(
        (w2 * w2 - w0 * w0 - w1 * w1) / (2.0 * w0 * w1),
        -1.0,
        1.0,
    )
    sine = math.sqrt(max(0.0, 1.0 - float(cosine) ** 2))
    vectors = np.asarray(
        [
            [1.0, 0.0, 0.0],
            [cosine, sine, 0.0],
            [-(w0 + w1 * cosine) / w2, -(w1 * sine) / w2, 0.0],
        ]
    )
    identity = np.eye(2, dtype=complex)
    paulis = np.asarray(
        [
            [[0.0, 1.0], [1.0, 0.0]],
            [[0.0, -1j], [1j, 0.0]],
            [[1.0, 0.0], [0.0, -1.0]],
        ],
        dtype=complex,
    )
    effects = np.asarray(
        [
            0.5 * weights[index]
            * (identity + np.tensordot(vectors[index], paulis, axes=1))
            for index in range(3)
        ]
    )
    assert np.linalg.norm(effects.sum(axis=0) - identity) < 2e-10
    return effects


def random_positive_operator(rng: np.random.Generator) -> np.ndarray:
    matrix = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
    return matrix @ matrix.conj().T


def test_horwitz_weights_and_interval_enclosure() -> None:
    rng = np.random.default_rng(20260823)
    for _ in range(200):
        beta = rng.uniform(1.0, 1.25)
        alpha = rng.uniform(beta, 2.0)
        weights = np.asarray(terminal_weights(alpha, beta))
        assert np.isclose(weights.sum(), 2.0, atol=2e-15)
        assert weights[0] >= weights[1] >= weights[2] >= 0.0

        half_alpha = rng.uniform(0.0, min(0.05, alpha - 1.0))
        half_beta = rng.uniform(0.0, min(0.03, beta - 1.0))
        box = {
            "terminal_alpha": (alpha - half_alpha, alpha + half_alpha),
            "terminal_beta": (beta - half_beta, beta + half_beta),
        }
        intervals = terminal_weight_intervals(box)
        for weight, (lower, upper) in zip(weights, intervals, strict=True):
            assert lower - 2e-15 <= weight <= upper + 2e-15


def test_probability_range_obeys_homogenized_inellipse() -> None:
    rng = np.random.default_rng(314159)
    for _ in range(100):
        beta = rng.uniform(1.01, 1.35)
        alpha = rng.uniform(max(beta, 1.01), 1.95)
        weights = np.asarray(terminal_weights(alpha, beta))
        effects = canonical_three_effect_povm(weights)
        operator = random_positive_operator(rng)
        q = np.asarray(
            [np.trace(operator @ effects[t]).real for t in range(3)]
        )
        scale = float(q.sum())
        x, y = float(q[0]), float(q[1])
        polynomial = (
            beta**2 * x**2
            + alpha**2 * y**2
            + cross_coefficient(alpha, beta) * x * y
            - 2.0 * beta * x * scale
            - 2.0 * alpha * y * scale
            + scale**2
        )
        assert polynomial <= 2e-10 * max(1.0, scale**2)


def test_aligned_projective_comparison_for_random_states() -> None:
    rng = np.random.default_rng(271828)
    for _ in range(100):
        beta = rng.uniform(1.01, 1.25)
        alpha = rng.uniform(beta, 1.8)
        weights = np.asarray(terminal_weights(alpha, beta))
        effects = canonical_three_effect_povm(weights)
        states = [random_positive_operator(rng) for _ in range(3)]
        priors = np.asarray([np.trace(state).real for state in states])
        ternary = sum(
            np.trace(effects[label] @ states[label]).real for label in range(3)
        )
        for retained in range(3):
            projector = effects[retained] / weights[retained]
            for complement in range(3):
                if complement == retained:
                    continue
                deleted = 3 - retained - complement
                projective = (
                    np.trace(projector @ states[retained]).real
                    + np.trace((np.eye(2) - projector) @ states[complement]).real
                )
                bonus = projective_comparison_bonus(
                    weights, priors, retained, complement
                )
                assert ternary - projective <= bonus + 2e-12
                assert deleted not in {retained, complement}


def test_planar_reconstruction_recovers_visible_bloch_coordinates() -> None:
    rng = np.random.default_rng(161803)
    paulis = np.asarray(
        [
            [[0.0, 1.0], [1.0, 0.0]],
            [[0.0, -1j], [1j, 0.0]],
            [[1.0, 0.0], [0.0, -1.0]],
        ],
        dtype=complex,
    )
    identity = np.eye(2, dtype=complex)
    for _ in range(100):
        beta = rng.uniform(1.02, 1.30)
        alpha = rng.uniform(max(beta, 1.05), 1.95)
        effects = canonical_three_effect_povm(
            np.asarray(terminal_weights(alpha, beta))
        )
        coefficients = rng.normal(size=4)
        operator = 0.5 * (
            coefficients[0] * identity
            + np.tensordot(coefficients[1:], paulis, axes=1)
        )
        measured = np.asarray(
            [np.trace(effect @ operator).real for effect in effects]
        )
        assert np.allclose(measured.sum(), coefficients[0], atol=2e-12)
        assert np.allclose(
            planar_reconstruction(alpha, beta) @ measured,
            coefficients[1:3],
            atol=3e-12,
        )


def test_outward_reconstruction_column_errors_cover_random_box_points() -> None:
    rng = np.random.default_rng(141421)
    for _ in range(30):
        beta_center = rng.uniform(1.04, 1.24)
        alpha_center = rng.uniform(beta_center + 0.03, 1.90)
        alpha_half = rng.uniform(1e-5, 5e-3)
        beta_half = rng.uniform(1e-5, 3e-3)
        alpha_bounds = (alpha_center - alpha_half, alpha_center + alpha_half)
        beta_bounds = (beta_center - beta_half, beta_center + beta_half)
        anchor, errors, _ = reconstruction_anchor_and_errors(
            alpha_bounds, beta_bounds
        )
        for _ in range(30):
            alpha = rng.uniform(*alpha_bounds)
            beta = rng.uniform(*beta_bounds)
            difference = planar_reconstruction(alpha, beta) - anchor
            assert np.all(
                np.linalg.norm(difference, axis=0) <= errors + 2e-15
            )
            signed = rng.normal(size=3)
            assert np.linalg.norm(difference @ signed) <= (
                errors @ np.abs(signed) + 2e-14
            )


def test_terminal_reconstructed_frontier_artifacts_are_self_consistent() -> None:
    summary = validate()
    assert summary["logical_status"] == (
        "strict local strengthening; target remains open"
    )
    assert summary["global_status"] == "continuous terminal strip still open"


def test_spectral_cover_builds_one_mixed_integer_selector_family() -> None:
    caps = tuple(
        np.append(normal, cosine) for normal, cosine in cube_face_caps(2)
    )
    anchor, errors, _ = reconstruction_anchor_and_errors(
        (1.92, 1.93), (1.14, 1.15)
    )
    oracle = TernaryConeOracle(
        0.55,
        (0, 1, 2, 3),
        (),
        (),
        0.79,
        0.7573,
        projective_support_lines=((0.6, 0.76591),),
        common_contractions=(
            {
                "coefficients": [1.0, 1.0, -1.0, -1.0],
                "branch": "spectral-cover",
                "scalar_signs": (1,),
                "caps": caps,
            },
        ),
        terminal_reconstruction=(anchor, errors),
    )
    assert oracle.problem.is_mixed_integer()
    assert len(oracle.common_contraction_selectors) == 1
    assert oracle.common_contraction_selectors[0].shape == (25,)

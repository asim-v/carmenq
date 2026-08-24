from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pytest


RESEARCH = Path(__file__).resolve().parents[1] / "scratch" / "d2_frontier"
sys.path.insert(0, str(RESEARCH))

from planar_cp_completion import (  # noqa: E402
    IDENTITY,
    PAULIS,
    cp_completion_radius,
    numerical_radius,
    numerical_radius_witness,
    reconstruct_planar_pullbacks,
)


def canonical_effects() -> np.ndarray:
    weights = np.asarray([0.92, 0.64, 0.44])
    cosine = (weights[2] ** 2 - weights[0] ** 2 - weights[1] ** 2) / (
        2.0 * weights[0] * weights[1]
    )
    sine = np.sqrt(1.0 - cosine * cosine)
    directions = np.asarray(
        [
            [1.0, 0.0, 0.0],
            [cosine, sine, 0.0],
            [
                -(weights[0] + weights[1] * cosine) / weights[2],
                -weights[1] * sine / weights[2],
                0.0,
            ],
        ]
    )
    effects = np.zeros((4, 2, 2), dtype=complex)
    for index in range(3):
        effects[index] = (
            0.5
            * weights[index]
            * (
                IDENTITY
                + sum(directions[index, axis] * PAULIS[axis] for axis in range(3))
            )
        )
    return effects


def test_identity_channel_is_on_ando_boundary() -> None:
    effects = canonical_effects()
    report = cp_completion_radius(effects, effects)
    assert report["compatible"]
    assert report["radius"] == pytest.approx(1.0, abs=2e-12)


def test_random_cp_map_is_completable_and_reconstructed() -> None:
    rng = np.random.default_rng(2026082304)
    effects = canonical_effects()
    kraus = rng.normal(size=(3, 2, 2)) + 1j * rng.normal(size=(3, 2, 2))
    scale = np.linalg.eigvalsh(np.einsum("kji,kjl->il", kraus.conj(), kraus))[-1]
    kraus /= np.sqrt(1.3 * scale)

    def adjoint(operator: np.ndarray) -> np.ndarray:
        return sum(item.conj().T @ operator @ item for item in kraus)

    pulled = np.asarray([adjoint(effect) for effect in effects])
    reconstructed = reconstruct_planar_pullbacks(pulled, effects)
    expected = (adjoint(IDENTITY), adjoint(PAULIS[0]), adjoint(PAULIS[1]))
    for actual, target in zip(reconstructed, expected, strict=True):
        assert np.linalg.norm(actual - target) < 2e-14
    report = cp_completion_radius(pulled, effects)
    assert report["compatible"]
    assert report["radius"] <= 1.0 + 2e-12


def test_singular_total_uses_support_completion() -> None:
    effects = canonical_effects()
    total = 0.5 * (IDENTITY + PAULIS[2])
    weights = np.trace(effects, axis1=1, axis2=2).real
    directions_x = np.asarray(
        [
            np.trace(effect @ PAULIS[0]).real / weights[index]
            for index, effect in enumerate(effects[:3])
        ]
    )
    pulled = np.zeros_like(effects)
    for index in range(3):
        pulled[index] = 0.5 * weights[index] * (1.0 + directions_x[index]) * total
    report = cp_completion_radius(pulled, effects)
    assert report["compatible"]
    assert report["support_rank"] == 1
    assert report["support_residual"] < 2e-14
    assert report["radius"] == pytest.approx(1.0, abs=2e-12)


def test_individually_valid_effects_can_fail_common_cp_completion() -> None:
    effects = canonical_effects()
    second = np.asarray(
        [
            [0.24175927315093138, -0.30528594092156486 + 0.05548766814312024j],
            [-0.30528594092156475 - 0.05548766814312024j, 0.39824073004173],
        ]
    )
    third = np.asarray(
        [
            [0.24654148066124223, 0.11892952988987542 - 0.18317018377261665j],
            [0.11892952988987542 + 0.18317018377261671j, 0.19345852058459123],
        ]
    )
    first = IDENTITY - second - third
    pulled = np.asarray([first, second, third, np.zeros((2, 2))])
    weights = np.trace(effects, axis1=1, axis2=2).real
    for index in range(3):
        assert np.linalg.eigvalsh(pulled[index]).min() >= -2e-15
        assert (
            np.linalg.eigvalsh(weights[index] * IDENTITY - pulled[index]).min() >= -2e-9
        )
    assert np.linalg.norm(pulled.sum(axis=0) - IDENTITY) < 2e-15
    report = cp_completion_radius(pulled, effects)
    assert not report["compatible"]
    assert report["radius"] == pytest.approx(1.72444066, abs=2e-8)
    assert report["witness_ratio"] == pytest.approx(report["radius"], abs=2e-12)
    witness_state = 0.5 * (
        IDENTITY
        + sum(report["witness_input_bloch"][axis] * PAULIS[axis] for axis in range(3))
    )
    total, x_pullback, y_pullback = reconstruct_planar_pullbacks(pulled, effects)
    denominator = np.trace(witness_state @ total).real
    numerator = np.hypot(
        np.trace(witness_state @ x_pullback).real,
        np.trace(witness_state @ y_pullback).real,
    )
    assert numerator > denominator
    assert numerator / denominator == pytest.approx(report["radius"], abs=2e-12)


def test_numerical_radius_known_cases() -> None:
    nilpotent = np.asarray([[0.0, 2.0], [0.0, 0.0]])
    assert numerical_radius(nilpotent) == pytest.approx(1.0, abs=2e-12)
    diagonal = np.diag([1.0 + 2.0j, -0.25])
    assert numerical_radius(diagonal) == pytest.approx(np.sqrt(5.0), abs=2e-12)
    radius, _, vector = numerical_radius_witness(diagonal)
    assert abs(np.vdot(vector, diagonal @ vector)) == pytest.approx(radius, abs=2e-12)

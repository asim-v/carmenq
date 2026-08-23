from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pytest

from carmenq.common_instrument import apply_choi, choi_from_kraus


RESEARCH = Path(__file__).resolve().parents[1] / "scratch" / "d2_frontier"
sys.path.insert(0, str(RESEARCH))

from common_instrument_exact_scip import (  # noqa: E402
    CHOI_BASIS,
    PAULIS,
    choi_pauli_coefficients,
    output_coefficients,
    pauli_coefficients,
    rotate_input_gauge,
)


def reconstruct_state(coefficients: np.ndarray) -> np.ndarray:
    return sum(coefficients[mu] * PAULIS[mu] for mu in range(4)) / 2.0


def reconstruct_choi(coefficients: np.ndarray) -> np.ndarray:
    return sum(
        coefficients[mu, nu] * CHOI_BASIS[mu, nu]
        for mu in range(4)
        for nu in range(4)
    ) / 4.0


def test_exact_pauli_application_matches_direct_choi() -> None:
    state = np.asarray([[0.31, 0.04 + 0.03j], [0.04 - 0.03j, 0.19]])
    operator = np.asarray([[0.7, 0.2j], [0.1, 0.5]], dtype=complex)
    choi = choi_from_kraus((operator,))
    result = output_coefficients(
        pauli_coefficients(state), choi_pauli_coefficients(choi)
    )
    expected = pauli_coefficients(apply_choi(choi, state))
    assert result == pytest.approx(expected)


def test_input_gauge_rotation_preserves_every_conditioned_output() -> None:
    generator = np.random.default_rng(2718)
    raw_states = []
    for trace in (0.31, 0.27, 0.23, 0.19):
        vector = generator.normal(size=3)
        vector *= 0.8 * trace / np.linalg.norm(vector)
        raw_states.append(np.r_[trace, vector])
    state = np.asarray(raw_states)

    operators = []
    for _ in range(4):
        item = generator.normal(size=(2, 2)) + 1j * generator.normal(size=(2, 2))
        operators.append(0.12 * item)
    choi = np.asarray(
        [choi_pauli_coefficients(choi_from_kraus((item,))) for item in operators]
    )
    rotated_state, rotated_choi = rotate_input_gauge(state, choi)

    before = np.asarray(
        [
            [output_coefficients(state[z], choi[y]) for y in range(4)]
            for z in range(4)
        ]
    )
    after = np.asarray(
        [
            [
                output_coefficients(rotated_state[z], rotated_choi[y])
                for y in range(4)
            ]
            for z in range(4)
        ]
    )
    assert after == pytest.approx(before, abs=1e-12)
    assert rotated_state[0, 2:] == pytest.approx([0.0, 0.0], abs=1e-12)
    assert rotated_state[0, 1] >= 0.0
    assert rotated_state[1, 3] == pytest.approx(0.0, abs=1e-12)
    assert rotated_state[1, 2] >= 0.0

    for original, rotated in zip(choi, rotated_choi, strict=True):
        assert np.linalg.eigvalsh(reconstruct_choi(original)) == pytest.approx(
            np.linalg.eigvalsh(reconstruct_choi(rotated)), abs=1e-12
        )
    for original, rotated in zip(state, rotated_state, strict=True):
        assert np.linalg.eigvalsh(reconstruct_state(original)) == pytest.approx(
            np.linalg.eigvalsh(reconstruct_state(rotated)), abs=1e-12
        )

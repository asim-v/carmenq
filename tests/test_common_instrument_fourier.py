from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
FRONTIER = ROOT / "scratch" / "d2_frontier"
if str(FRONTIER) not in sys.path:
    sys.path.insert(0, str(FRONTIER))

from fourier_behavior_cap_cover import cube_face_caps, plane_caps  # noqa: E402
from fourier_behavior_upper import CHARACTERS  # noqa: E402


def trace_norm(matrix: np.ndarray) -> float:
    return float(np.linalg.svd(matrix, compute_uv=False).sum())


def random_instrument(rng: np.random.Generator) -> list[np.ndarray]:
    raw = [
        rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
        for _ in range(4)
    ]
    total = sum(kraus.conj().T @ kraus for kraus in raw)
    eigenvalues, eigenvectors = np.linalg.eigh(total)
    inverse_root = (eigenvectors * eigenvalues**-0.5) @ eigenvectors.conj().T
    return [kraus @ inverse_root for kraus in raw]


def test_flagged_fourier_contraction_for_one_common_instrument() -> None:
    rng = np.random.default_rng(20260823)
    kraus = random_instrument(rng)
    priors = rng.dirichlet(np.ones(4))
    states = []
    for prior in priors:
        vector = rng.normal(size=3)
        vector *= rng.uniform() / np.linalg.norm(vector)
        matrix = np.array(
            [
                [1.0 + vector[2], vector[0] - 1j * vector[1]],
                [vector[0] + 1j * vector[1], 1.0 - vector[2]],
            ]
        ) / 2.0
        states.append(prior * matrix)

    outputs = np.asarray(
        [[operator @ state @ operator.conj().T for operator in kraus] for state in states]
    )
    terminal = np.asarray(
        [sum(outputs[z, z ^ s] for z in range(4)) for s in range(4)]
    )
    for character in CHARACTERS:
        input_fourier = sum(character[z] * states[z] for z in range(4))
        blocks = [sum(character[z] * outputs[z, y] for z in range(4)) for y in range(4)]
        terminal_fourier = sum(character[s] * terminal[s] for s in range(4))
        signed_blocks = sum(character[y] * blocks[y] for y in range(4))
        assert np.allclose(terminal_fourier, signed_blocks, atol=2e-12)
        flagged = sum(trace_norm(block) for block in blocks)
        assert flagged <= trace_norm(input_fourier) + 2e-12
        assert trace_norm(terminal_fourier) <= flagged + 2e-12


def test_qubit_hermitian_trace_norm_is_scalar_vector_maximum() -> None:
    rng = np.random.default_rng(7)
    paulis = (
        np.eye(2),
        np.array([[0.0, 1.0], [1.0, 0.0]]),
        np.array([[0.0, -1j], [1j, 0.0]]),
        np.array([[1.0, 0.0], [0.0, -1.0]]),
    )
    for _ in range(30):
        coefficients = rng.normal(size=4)
        matrix = sum(coefficients[index] * paulis[index] for index in range(4)) / 2.0
        expected = max(abs(coefficients[0]), np.linalg.norm(coefficients[1:]))
        assert np.isclose(trace_norm(matrix), expected, atol=2e-12)


def test_order_fixes_two_fourier_scalar_signs() -> None:
    rng = np.random.default_rng(11)
    for _ in range(100):
        priors = np.sort(rng.dirichlet(np.ones(4)))[::-1]
        scalars = CHARACTERS @ priors
        assert scalars[0] >= -1e-15
        assert scalars[1] >= -1e-15


def test_proved_plane_and_cube_face_caps_cover_sampled_directions() -> None:
    rng = np.random.default_rng(13)
    planar = plane_caps(16)
    spherical = cube_face_caps(8)
    for _ in range(1000):
        angle = rng.uniform(-np.pi / 2.0, np.pi / 2.0)
        vector = np.asarray([np.cos(angle), 0.0, np.sin(angle)])
        assert max(normal @ vector / cosine for normal, cosine in planar) >= 1.0 - 1e-12

        vector = rng.normal(size=3)
        vector /= np.linalg.norm(vector)
        assert max(normal @ vector / cosine for normal, cosine in spherical) >= 1.0 - 1e-12

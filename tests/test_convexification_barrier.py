from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pytest


RESEARCH = Path(__file__).resolve().parents[1] / "scratch" / "d2_frontier"
sys.path.insert(0, str(RESEARCH))

from convexification_barrier import (  # noqa: E402
    block_trace_coherence_return,
    convexified_leaf_barrier,
    pure_pinched_return,
)


def test_l055_fixed_readout_barrier_is_reproduced() -> None:
    result = convexified_leaf_barrier(np.asarray([0.92, 0.64, 0.44, 0.0]), 0.55)
    assert result["bound"] == pytest.approx(0.790265609741724, abs=1e-15)
    assert result["score_from_distribution"] == pytest.approx(
        result["bound"], abs=2e-15
    )
    assert sum(result["path_probabilities"]) == pytest.approx(1.0)
    assert result["eigen_residual"] < 1e-14


def test_rayleigh_bound_dominates_random_path_distributions() -> None:
    rng = np.random.default_rng(20260823)
    result = convexified_leaf_barrier(np.asarray([0.92, 0.64, 0.44, 0.0]), 0.55)
    expanded = np.repeat(np.asarray(result["effect_norms"]), 4)
    for _ in range(100):
        probability = rng.dirichlet(np.ones(16))
        audit = float(expanded @ probability)
        returned = float(np.square(np.sqrt(probability).sum()) / 16.0)
        assert 0.55 * audit + 0.45 * returned <= result["bound"] + 1e-14


def test_input_validation() -> None:
    with pytest.raises(ValueError, match=r"\[0,1\]"):
        convexified_leaf_barrier(np.asarray([1.1]), 0.5)
    with pytest.raises(ValueError, match="positive"):
        convexified_leaf_barrier(np.asarray([0.5]), 0.5, multiplicity=0)


def test_block_trace_coherence_equals_return_on_pure_states() -> None:
    rng = np.random.default_rng(1729)
    block_sizes = (2, 3, 1, 4)
    vector = rng.normal(size=10) + 1j * rng.normal(size=10)
    vector /= np.linalg.norm(vector)
    density = np.outer(vector, vector.conj())
    assert block_trace_coherence_return(density, block_sizes) == pytest.approx(
        pure_pinched_return(vector, block_sizes), abs=2e-15
    )


def test_block_trace_coherence_is_convex() -> None:
    rng = np.random.default_rng(314159)
    block_sizes = (2, 2, 2)
    vectors = []
    for _ in range(2):
        vector = rng.normal(size=6) + 1j * rng.normal(size=6)
        vectors.append(vector / np.linalg.norm(vector))
    states = [np.outer(vector, vector.conj()) for vector in vectors]
    mixture = 0.37 * states[0] + 0.63 * states[1]
    mixed_value = block_trace_coherence_return(mixture, block_sizes)
    average = 0.37 * block_trace_coherence_return(states[0], block_sizes)
    average += 0.63 * block_trace_coherence_return(states[1], block_sizes)
    assert mixed_value <= average + 2e-15

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pytest


RESEARCH = Path(__file__).resolve().parents[1] / "scratch" / "d2_frontier"
sys.path.insert(0, str(RESEARCH))

from coherence_polar_program import (  # noqa: E402
    coherence_from_witness,
    contraction_residuals,
    polar_coherence_witness,
    trace_norm_polar,
)
from convexification_barrier import block_trace_coherence_return  # noqa: E402


def random_density(rng: np.random.Generator, dimension: int) -> np.ndarray:
    factor = rng.normal(size=(dimension, dimension))
    factor = factor + 1j * rng.normal(size=(dimension, dimension))
    density = factor @ factor.conj().T
    return density / np.trace(density)


def test_rectangular_polar_factor_attains_trace_norm() -> None:
    rng = np.random.default_rng(2026082301)
    block = rng.normal(size=(3, 5)) + 1j * rng.normal(size=(3, 5))
    contraction, trace_norm = trace_norm_polar(block)
    attained = np.trace(contraction.conj().T @ block).real
    assert attained == pytest.approx(trace_norm, abs=2e-14)
    assert np.linalg.svd(contraction, compute_uv=False)[0] <= 1.0 + 2e-15


def test_polar_program_equals_block_trace_coherence() -> None:
    rng = np.random.default_rng(2026082302)
    block_sizes = (2, 3, 1, 4)
    density = random_density(rng, sum(block_sizes))
    witness, _ = polar_coherence_witness(density, block_sizes)
    affine_value = coherence_from_witness(density, witness) / len(block_sizes)
    direct_value = block_trace_coherence_return(density, block_sizes)
    assert affine_value == pytest.approx(direct_value, abs=3e-15)


def test_polar_witness_is_feasible() -> None:
    rng = np.random.default_rng(2026082303)
    block_sizes = (2, 2, 3)
    density = random_density(rng, sum(block_sizes))
    witness, contractions = polar_coherence_witness(density, block_sizes)
    audit = contraction_residuals(witness, block_sizes)
    assert len(contractions) == 3
    assert audit["hermiticity"] < 1e-14
    assert audit["diagonal"] == 0.0
    assert audit["maximum_pair_operator_norm"] <= 1.0 + 2e-15


def test_path_diagonal_mixture_gets_no_hellinger_bonus() -> None:
    probability = np.asarray([0.55, 0.25, 0.15, 0.05])
    density = np.diag(probability)
    witness, _ = polar_coherence_witness(density, (1, 1, 1, 1))
    value = coherence_from_witness(density, witness) / 4.0
    assert value == pytest.approx(0.25)
    assert np.linalg.norm(witness) == pytest.approx(0.0)

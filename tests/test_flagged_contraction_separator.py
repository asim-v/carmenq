from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pytest


RESEARCH = Path(__file__).resolve().parents[1] / "scratch" / "d2_frontier"
sys.path.insert(0, str(RESEARCH))

from flagged_contraction_separator import (  # noqa: E402
    contraction_report,
    qubit_trace_norm,
)


def test_qubit_trace_norm_uses_the_larger_scalar_or_bloch_radius() -> None:
    assert qubit_trace_norm(np.asarray([0.7, 0.1, 0.2, 0.3])) == pytest.approx(0.7)
    assert qubit_trace_norm(np.asarray([0.1, 0.3, 0.4, 0.0])) == pytest.approx(0.5)


def test_identity_flagged_channel_saturates_every_contraction() -> None:
    prefix = np.asarray(
        [
            [0.25, 0.10, 0.00, 0.00],
            [0.25, 0.00, 0.10, 0.00],
            [0.25, 0.00, 0.00, 0.10],
            [0.25, -0.05, -0.05, -0.05],
        ]
    )
    conditioned = np.zeros((4, 4, 4))
    conditioned[:, 0, :] = prefix
    report = contraction_report(
        prefix, conditioned, np.asarray([0.3, -0.7, 0.2, 0.4])
    )
    assert report["violation"] == pytest.approx(0.0, abs=1e-12)


def test_nonphysical_output_amplification_is_separated() -> None:
    prefix = np.repeat(np.asarray([[0.25, 0.0, 0.0, 0.0]]), 4, axis=0)
    conditioned = np.zeros((4, 4, 4))
    conditioned[0, 0] = [0.25, 0.20, 0.0, 0.0]
    conditioned[1, 0] = [0.25, -0.20, 0.0, 0.0]
    report = contraction_report(prefix, conditioned, np.asarray([1.0, -1.0, 0.0, 0.0]))
    assert report["input_trace_norm"] == pytest.approx(0.0)
    assert report["violation"] == pytest.approx(0.4 / np.sqrt(2.0))


def test_all_contractions_deliberately_retain_the_positive_non_cp_transpose() -> None:
    tetrahedron = np.asarray(
        [
            [1.0, 1.0, 1.0],
            [1.0, -1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [-1.0, -1.0, 1.0],
        ]
    ) / np.sqrt(3.0)
    prefix = np.column_stack((np.full(4, 0.25), 0.25 * tetrahedron))
    conditioned = np.zeros((4, 4, 4))
    conditioned[:, 0, :] = prefix * np.asarray([1.0, 1.0, -1.0, 1.0])
    rng = np.random.default_rng(19)
    for _ in range(30):
        report = contraction_report(prefix, conditioned, rng.normal(size=4))
        assert report["violation"] == pytest.approx(0.0, abs=2e-15)

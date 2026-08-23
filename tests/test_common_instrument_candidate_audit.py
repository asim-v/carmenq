from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pytest


RESEARCH = Path(__file__).resolve().parents[1] / "scratch" / "d2_frontier"
sys.path.insert(0, str(RESEARCH))

from audit_common_instrument_candidate import (  # noqa: E402
    load_reported_family,
    matrix_from_bloch,
)


def test_matrix_from_bloch_uses_trace_first_convention() -> None:
    matrix = matrix_from_bloch(np.asarray([0.4, 0.1, -0.2, 0.3]))
    assert np.trace(matrix).real == pytest.approx(0.4)
    assert np.linalg.eigvalsh(matrix) == pytest.approx(
        np.asarray([0.4 - np.sqrt(0.14), 0.4 + np.sqrt(0.14)]) / 2.0
    )


def test_reported_family_loader_reconstructs_all_conditioned_states() -> None:
    payload = {
        "prefix_bloch_coefficients": [[0.25, 0.0, 0.0, 0.0]] * 4,
        "path_probabilities": [[0.0625] * 4] * 4,
        "conditioned_output_bloch_vectors": [[[[0.0, 0.0, 0.0]][0] for _ in range(4)]] * 4,
    }
    states, outputs = load_reported_family(payload)
    assert states.shape == (4, 2, 2)
    assert outputs.shape == (4, 4, 2, 2)
    assert np.trace(states, axis1=1, axis2=2).real == pytest.approx([0.25] * 4)
    assert np.trace(outputs, axis1=2, axis2=3).real == pytest.approx(
        np.full((4, 4), 0.0625)
    )

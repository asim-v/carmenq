from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pytest


pytest.importorskip("cvxpy")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scratch" / "d2_frontier"))

from qubit_discrimination_geometry import discrimination_geometry  # noqa: E402


def test_weighted_ball_recovers_perfect_binary_discrimination() -> None:
    states = np.zeros((4, 2, 2), dtype=complex)
    states[0, 0, 0] = 0.5
    states[1, 1, 1] = 0.5
    result = discrimination_geometry(states)
    assert result["optimal_guess_probability"] == pytest.approx(1.0, abs=2e-10)
    assert set(result["active_indices"]) == {0, 1}


def test_weighted_ball_rejects_wrong_shape() -> None:
    with pytest.raises(ValueError, match="shape"):
        discrimination_geometry(np.zeros((3, 2, 2), dtype=complex))

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scratch" / "d2_frontier"))

from validate_projective_line_l055 import (  # noqa: E402
    WEIGHT,
    secular_rescaling_factor,
)


def test_secular_rescaling_is_pointwise_valid() -> None:
    old_level = 0.7568534
    new_level = 0.75730
    factor = secular_rescaling_factor(old_level, new_level)
    for decision in np.linspace(0.0, 1.0, 1001):
        new_term = 1.0 / (new_level - WEIGHT * decision)
        rescaled_old_term = factor / (old_level - WEIGHT * decision)
        assert new_term <= rescaled_old_term + 2e-15


def test_secular_rescaling_rejects_invalid_levels() -> None:
    with pytest.raises(ValueError):
        secular_rescaling_factor(WEIGHT, 0.8)
    with pytest.raises(ValueError):
        secular_rescaling_factor(0.8, 0.7)

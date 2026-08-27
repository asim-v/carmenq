from __future__ import annotations

import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from check_data_reproducibility import _compare_json_value  # noqa: E402


def test_json_comparison_handles_exact_integers_beyond_float_range() -> None:
    enormous = 10**1000 + 123

    _compare_json_value(enormous, enormous, "certificate.numerator")
    with pytest.raises(AssertionError, match="certificate.numerator changed"):
        _compare_json_value(enormous, enormous + 1, "certificate.numerator")


def test_json_comparison_uses_tolerance_only_for_floats() -> None:
    _compare_json_value(1.0, 1.0 + 5e-13, "estimate")
    with pytest.raises(AssertionError, match="estimate changed"):
        _compare_json_value(1.0, 1.0 + 5e-9, "estimate")


def test_json_comparison_does_not_tolerate_integer_drift() -> None:
    with pytest.raises(AssertionError, match="count changed"):
        _compare_json_value(10**18, 10**18 + 1, "count")

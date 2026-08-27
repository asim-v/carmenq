"""Checks for exact RHS transfer of ternary projective premises."""

from __future__ import annotations

from fractions import Fraction
from pathlib import Path
from types import SimpleNamespace
import sys

import numpy as np
import pytest
from scipy import sparse


FRONTIER = Path(__file__).resolve().parents[1] / "scratch" / "d2_frontier"
if str(FRONTIER) not in sys.path:
    sys.path.insert(0, str(FRONTIER))

from ternary_projective_line_sensitivity import (  # noqa: E402
    exact_delta_objective,
    require_line_only_change,
)


def canonical_data(
    b: list[float], *, c: list[float] | None = None
) -> dict[str, object]:
    return {
        "A": sparse.csc_matrix(
            np.asarray([[1.0, 0.0], [0.0, -2.0], [3.0, 4.0]])
        ),
        "b": np.asarray(b),
        "c": np.asarray(c or [5.0, 6.0]),
        "dims": SimpleNamespace(zero=1, nonneg=2, soc=[3]),
    }


def test_line_only_change_extracts_exact_dyadic_rhs_delta() -> None:
    baseline = canonical_data([0.25, 0.5, -0.75])
    modified = canonical_data([0.5, 0.5, -0.5])
    delta = require_line_only_change(baseline, modified)
    assert delta == [Fraction(1, 4), Fraction(0), Fraction(1, 4)]
    support = [(index, value) for index, value in enumerate(delta) if value]
    dual = np.asarray([0.5, 99.0, -0.25])
    assert exact_delta_objective(len(delta), support, dual) == Fraction(1, 16)


def test_line_only_change_rejects_an_objective_change() -> None:
    baseline = canonical_data([0.0, 0.0, 0.0])
    modified = canonical_data([0.0, 0.0, 0.0], c=[5.0, 7.0])
    with pytest.raises(RuntimeError, match="altered A or c"):
        require_line_only_change(baseline, modified)


def test_sparse_delta_rejects_a_wrong_dual_dimension() -> None:
    with pytest.raises(RuntimeError, match="incompatible lengths"):
        exact_delta_objective(3, [(0, Fraction(1))], np.zeros(2))

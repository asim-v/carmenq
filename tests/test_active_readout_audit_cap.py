from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest


RESEARCH = Path(__file__).resolve().parents[1] / "scratch" / "d2_frontier"
SPEC = importlib.util.spec_from_file_location(
    "active_readout_audit_cap", RESEARCH / "active_readout_audit_cap.py"
)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_reserve_endpoint_and_unbiased_dual() -> None:
    assert MODULE.reserve(0.0, -0.37) == pytest.approx(0.5)
    assert MODULE.reserve(1.0, -1.0) == pytest.approx(0.0)
    assert MODULE.reserve(1.0, 0.25) == pytest.approx(1.0)


def test_projection_vertices_obey_weighted_closure() -> None:
    weights = np.asarray([0.9, 0.5, 0.4, 0.2])
    vertices = MODULE.projection_vertices(weights)
    assert vertices
    for point in vertices:
        assert np.max(np.abs(point)) <= 1.0 + 1e-12
        assert np.dot(weights, point) == pytest.approx(0.0, abs=2e-12)


def test_balanced_four_active_readout_has_half_cap() -> None:
    result = MODULE.active_audit_cap(np.full(4, 0.5))
    assert result["minimum_total_prior_reserve"] == pytest.approx(2.0, abs=2e-10)
    assert result["audit_upper"] == pytest.approx(0.5, abs=5e-11)


def test_trine_cap_recovers_two_thirds() -> None:
    result = MODULE.active_audit_cap(np.full(3, 2.0 / 3.0))
    assert result["audit_upper"] == pytest.approx(2.0 / 3.0, abs=5e-10)


def test_asymmetric_four_active_cap_is_nontrivial() -> None:
    result = MODULE.active_audit_cap(np.asarray([0.9, 0.5, 0.4, 0.2]))
    assert result["audit_upper"] < 0.59
    assert result["projection_residual"] == pytest.approx(0.0, abs=2e-12)

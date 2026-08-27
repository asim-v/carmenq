"""Regression checks for the exact global-sector arithmetic."""

from __future__ import annotations

from fractions import Fraction
import sys
from pathlib import Path

import pytest


FRONTIER = Path(__file__).resolve().parents[1] / "scratch" / "d2_frontier"
if str(FRONTIER) not in sys.path:
    sys.path.insert(0, str(FRONTIER))

from four_active_socp_exact_cover import PROJECTIVE_LINES  # noqa: E402
from validate_projective_tangent_global import KERNEL, validate_payload  # noqa: E402
from verify_global_frontier_l060 import (  # noqa: E402
    DELETION_UPPER,
    FOUR_ACTIVE_UPPER,
    REPOSITORY_ROOT,
    portable_path,
    LOW_WEIGHT_TARGET,
    EXPECTED_PROJECTIVE_PREMISES,
    PROJECTIVE_055_UPPER,
    PROJECTIVE_060_UPPER,
    REPORTED_UPPER,
    TERNARY_UPPER,
    validate_named_projective_premises,
)


def test_four_active_uses_only_completed_projective_targets() -> None:
    assert PROJECTIVE_LINES == (
        (Fraction(11, 20), PROJECTIVE_055_UPPER),
        (Fraction(3, 5), PROJECTIVE_060_UPPER),
    )

def test_global_manifest_paths_are_repository_relative() -> None:
    artifact = REPOSITORY_ROOT / "data" / "certificate.json"
    assert portable_path(artifact) == "data/certificate.json"



def test_global_endpoint_is_exactly_the_deletion_sector() -> None:
    assert LOW_WEIGHT_TARGET == Fraction(76591, 100000)
    assert TERNARY_UPPER == Fraction(76652, 100000)
    assert FOUR_ACTIVE_UPPER == Fraction(76670, 100000)
    assert DELETION_UPPER == Fraction(76670, 100000)
    assert max(PROJECTIVE_060_UPPER, DELETION_UPPER) == REPORTED_UPPER
    assert REPORTED_UPPER == Fraction(76670, 100000)


def test_global_assembly_requires_the_exact_named_projective_premises() -> None:
    validate_named_projective_premises(
        {"named_projective_premises": EXPECTED_PROJECTIVE_PREMISES}, "component"
    )


def test_global_assembly_rejects_a_dependency_mismatch() -> None:
    with pytest.raises(RuntimeError, match="incompatible projective premises"):
        validate_named_projective_premises(
            {
                "named_projective_premises": {
                    "11/20": "761/1000",
                    "3/5": "38331/50000",
                }
            },
            "component",
        )


def _single_cell_payload(level: str) -> tuple[dict[str, object], dict[str, object]]:
    expected = {"cell": {"bounds": {"x": ("0", "1")}}}
    certificate = {
        "complete": True,
        "boxes_remaining": 0,
        "weight": "3/5",
        "level": level,
        "proof_kernel": KERNEL,
        "trusted_optimizers": False,
        "boxes_split": 0,
        "closed_methods": {"trace": 1},
    }
    payload = {
        "weight": "3/5",
        "level": level,
        "run_complete": True,
        "all_cells_complete": True,
        "total_boxes_split": 0,
        "cells": [
            {
                "certificate_name": "cell",
                "bounds": {"x": ["0", "1"]},
                "certificate": certificate,
            }
        ],
    }
    return payload, expected


def test_projective_assembly_accepts_a_stronger_component_certificate() -> None:
    payload, expected = _single_cell_payload("3/4")
    summary = validate_payload(
        payload, expected, "certificate_name", "3/5", Fraction(4, 5)
    )
    assert summary["certified_level"] == "3/4"


def test_projective_assembly_rejects_a_weaker_component_certificate() -> None:
    payload, expected = _single_cell_payload("4/5")
    with pytest.raises(RuntimeError, match="weaker level"):
        validate_payload(
            payload, expected, "certificate_name", "3/5", Fraction(3, 4)
        )

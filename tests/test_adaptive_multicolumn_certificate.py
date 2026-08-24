from __future__ import annotations

import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
RESEARCH = ROOT / "scratch" / "d2_frontier"
sys.path.insert(0, str(RESEARCH))

from audit_adaptive_multicolumn_certificate import audit_certificate  # noqa: E402
from adaptive_multicolumn_branch_tree import json_bound  # noqa: E402


def test_recorded_worst_cell_tree_is_exhaustive_and_closed() -> None:
    certificate = RESEARCH / "adaptive_multicolumn_worstcell_l055_auditable.json"
    payload = json.loads(certificate.read_text(encoding="utf-8"))
    summary = audit_certificate(payload)
    assert summary == {
        "certificate_complete": True,
        "expanded_nodes": 48,
        "source_closed_nodes": 0,
        "infeasible_source_nodes": 0,
        "closed_leaves": 4657,
        "open_leaves": 0,
        "maximum_depth": 3,
        "branches_per_expansion": 98,
    }


def test_solver_failure_bound_stays_open_in_standard_json() -> None:
    assert json_bound(float("inf")) == "+inf"
    assert json_bound(float("-inf")) is None


def test_recorded_source_closed_tree_is_auditable() -> None:
    certificate = (
        RESEARCH
        / "regular_multicolumn_forest_l055"
        / "p004_s037_bloch_c3.json"
    )
    payload = json.loads(certificate.read_text(encoding="utf-8"))
    summary = audit_certificate(payload)
    assert summary["certificate_complete"] is True
    assert summary["source_closed_nodes"] == 1
    assert summary["open_leaves"] == 0

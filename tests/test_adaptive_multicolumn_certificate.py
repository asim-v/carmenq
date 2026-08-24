from __future__ import annotations

import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
RESEARCH = ROOT / "scratch" / "d2_frontier"
sys.path.insert(0, str(RESEARCH))

from audit_adaptive_multicolumn_certificate import audit_certificate  # noqa: E402


def test_recorded_worst_cell_tree_is_exhaustive_and_closed() -> None:
    certificate = RESEARCH / "adaptive_multicolumn_worstcell_l055_auditable.json"
    payload = json.loads(certificate.read_text(encoding="utf-8"))
    summary = audit_certificate(payload)
    assert summary == {
        "certificate_complete": True,
        "expanded_nodes": 48,
        "closed_leaves": 4657,
        "open_leaves": 0,
        "maximum_depth": 3,
        "branches_per_expansion": 98,
    }

from __future__ import annotations

import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
RESEARCH = ROOT / "scratch" / "d2_frontier"
sys.path.insert(0, str(RESEARCH))

from adaptive_multicolumn_regular_forest import open_branch_orbits  # noqa: E402
from aggregate_multicolumn_regular_forest import aggregate  # noqa: E402


def test_open_base_branches_form_exact_conjugation_pairs() -> None:
    path = RESEARCH / "fourier_behavior_pair_cover_p8_g4_l055_auditable.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    orbits = open_branch_orbits(payload, 0.758)
    assert len(orbits) == 353
    assert sum(len(item["orbit"]) for item in orbits) == 706
    assert all(len(item["orbit"]) == 2 for item in orbits)
    assert all(item["representative"][2] == "bloch" for item in orbits)


def test_recorded_global_regular_forest_is_complete_and_portable() -> None:
    base = RESEARCH / "fourier_behavior_pair_cover_p8_g4_l055_auditable.json"
    manifests = [
        RESEARCH / "regular_multicolumn_forest_top3_l055.json",
        RESEARCH / "regular_multicolumn_forest_batch0_l055.json",
        RESEARCH / "regular_multicolumn_forest_batch1_l055.json",
        RESEARCH / "regular_multicolumn_forest_batch2_l055.json",
    ]
    result = aggregate(base, manifests, 0.758, require_complete=True)
    assert result["certificate_complete"] is True
    assert result["audited_complete_orbits"] == 353
    assert result["missing_or_open_orbits"] == 0
    assert result["expanded_nodes"] == 2698
    assert result["closed_leaves"] == 262059
    assert result["source_closed_leaves"] == 236
    assert result["infeasible_leaves"] == 151733
    assert result["maximum_depth"] == 7
    assert result["maximum_finite_terminal_bound"] == 0.7579983961104495
    assert result["leaf_identity_verified"] is True
    assert all("\\" not in row["certificate"] for row in result["orbits"])

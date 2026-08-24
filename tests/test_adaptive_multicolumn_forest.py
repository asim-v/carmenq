from __future__ import annotations

import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
RESEARCH = ROOT / "scratch" / "d2_frontier"
sys.path.insert(0, str(RESEARCH))

from adaptive_multicolumn_regular_forest import open_branch_orbits  # noqa: E402


def test_open_base_branches_form_exact_conjugation_pairs() -> None:
    path = RESEARCH / "fourier_behavior_pair_cover_p8_g4_l055_auditable.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    orbits = open_branch_orbits(payload, 0.758)
    assert len(orbits) == 353
    assert sum(len(item["orbit"]) for item in orbits) == 706
    assert all(len(item["orbit"]) == 2 for item in orbits)
    assert all(item["representative"][2] == "bloch" for item in orbits)

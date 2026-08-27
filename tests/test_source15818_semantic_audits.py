from __future__ import annotations

import importlib
import json
import subprocess
import sys
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FRONTIER = ROOT / "scratch" / "d2_frontier"
if str(FRONTIER) not in sys.path:
    sys.path.insert(0, str(FRONTIER))


def test_source15818_completed_ellipse_socs_are_exact_outer_relaxations() -> None:
    enclosure = importlib.import_module("audit_source15818_enclosures")
    source = enclosure.json.loads(
        (ROOT / enclosure.DEFAULT_FRONTIER).read_text(encoding="utf-8")
    )
    report = enclosure.audit_box(
        tuple(source["box"]["terminal_alpha"]),
        tuple(source["box"]["terminal_beta"]),
    )

    assert report["anchor_count"] == 29
    assert report["all_anchors_certified"] is True
    assert report["all_certified"] is True
    assert Fraction(*report["minimum_anchor_margin"]) > 0
    lower = report["coefficientwise_lower"]
    assert lower["present"] is True
    assert lower["certified_outer"] is True
    assert Fraction(*lower["completion_excess"]) < 0


def test_source15818_every_source_constraint_has_canonical_provenance() -> None:
    # CVXPY allocates process-global expression identifiers.  Run this
    # bitwise-reproducibility audit in a fresh interpreter so earlier tests
    # cannot perturb the canonical column order.
    completed = subprocess.run(
        [
            sys.executable,
            str(FRONTIER / "audit_source15818_canonicalization.py"),
            "--basis-columns",
            "8",
        ],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    report = json.loads(completed.stdout)

    assert report["source_constraints"] == 1142
    assert report["canonical_constraints"] == 1306
    assert report["canonical_rows"] == 3173
    assert report["canonical_variables"] == 274
    assert report["canonical_nonzeros"] == 20962
    assert report["canonical_data_sha256"] == (
        "0861e28c987a2fdf03864ec8f753f70698e8cd3e8ba3b241ba715379acf0f1cf"
    )
    assert report["origin_counts"] == {
        "source": 1142,
        "implicit_variable_domain": 164,
    }
    assert report["all_source_constraints_mapped"] is True
    assert report["row_spans_complete"] is True
    assert report["cone_order_ok"] is True
    assert report["manifest_sha256"] == (
        "6f1ec1704c952520b3f677d8ad2388bd5730d3080e9ad98c8f0243e390b14064"
    )
    direct = report["direct_evaluation"]
    assert direct["columns_checked"] == 8
    assert direct["right_exact"] is True
    assert direct["matrix_failures"] == 0
    assert direct["objective_failures"] == 0
    assert report["passed"] is True

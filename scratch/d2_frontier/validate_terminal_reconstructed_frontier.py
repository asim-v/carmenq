"""Audit the local terminal-reconstructed common-instrument artifacts."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
TARGET = 0.758
EXPECTED_BOX = {
    "terminal_alpha": [1.923828125, 1.92578125],
    "terminal_beta": [1.1453718354430378, 1.149525316455696],
}


def load(name: str) -> dict[str, Any]:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def validate() -> dict[str, object]:
    measured_all = load("ternary_common_instrument_top_leaf_coarse_l055.json")
    measured_pbb = load("ternary_common_instrument_top_leaf_pbb_p16_l055.json")
    reconstructed_pbb = load(
        "ternary_reconstructed_fourier_top_leaf_pbb_p16_l055.json"
    )
    reconstructed_bbb = load(
        "ternary_reconstructed_fourier_top_leaf_bbb_p4_g4_l055.json"
    )
    adaptive_bbb = load(
        "ternary_reconstructed_multicolumn_top_leaf_bbb_p1_s92_l055.json"
    )
    payloads = (
        measured_all,
        measured_pbb,
        reconstructed_pbb,
        reconstructed_bbb,
        adaptive_bbb,
    )
    for payload in payloads:
        if payload["box"] != EXPECTED_BOX:
            raise RuntimeError("frontier artifacts do not use the same terminal box")
        if not math.isclose(float(payload["support_weight"]), 0.55, abs_tol=1e-15):
            raise RuntimeError("wrong support weight")
        if not payload["complete"]:
            raise RuntimeError("an angular cover is incomplete")
        if "solver-conditional" not in payload["scope"]:
            raise RuntimeError("artifact scope lost its numerical qualifier")

    measured_pbb_max = float(measured_pbb["maximum_bound"])
    measured_bbb_max = max(
        float(row["bound"])
        for row in measured_all["cells"]
        if row["branch"] == "bbb"
    )
    if measured_pbb_max < TARGET or measured_bbb_max < TARGET:
        raise RuntimeError("measured-only controls unexpectedly close the target")

    if reconstructed_pbb["branches"] != ["pbb"]:
        raise RuntimeError("wrong reconstructed pbb branch")
    if reconstructed_bbb["branches"] != ["bbb"]:
        raise RuntimeError("wrong reconstructed bbb branch")
    if len(reconstructed_pbb["cells"]) != 16:
        raise RuntimeError("pbb plane cover is incomplete")
    if len(reconstructed_bbb["cells"]) != 384:
        raise RuntimeError("bbb plane-sphere cover is incomplete")
    for payload in (reconstructed_pbb, reconstructed_bbb):
        audit = payload.get("terminal_reconstruction")
        if audit is None or len(audit.get("column_errors", ())) != 3:
            raise RuntimeError("missing terminal reconstruction enclosure")
        if any(float(value) < 0.0 for value in audit["column_errors"]):
            raise RuntimeError("negative reconstruction error")
    reconstructed_pbb_max = float(reconstructed_pbb["maximum_bound"])
    reconstructed_bbb_max = float(reconstructed_bbb["maximum_bound"])
    adaptive_bbb_max = float(adaptive_bbb["maximum_bound"])
    if not reconstructed_pbb_max < measured_pbb_max - 1e-4:
        raise RuntimeError("planar reconstruction did not strengthen pbb")
    if reconstructed_bbb_max > measured_bbb_max + 1e-7:
        raise RuntimeError("the stronger bbb relaxation increased its bound")
    if not adaptive_bbb_max < reconstructed_bbb_max - 1e-5:
        raise RuntimeError("the adaptive contraction did not strengthen bbb")
    if min(reconstructed_pbb_max, reconstructed_bbb_max, adaptive_bbb_max) < TARGET:
        raise RuntimeError("artifact unexpectedly claims a local target closure")

    return {
        "target": TARGET,
        "terminal_box": EXPECTED_BOX,
        "measured_pbb_maximum": measured_pbb_max,
        "measured_bbb_maximum": measured_bbb_max,
        "reconstructed_pbb_maximum": reconstructed_pbb_max,
        "reconstructed_bbb_maximum": reconstructed_bbb_max,
        "adaptive_bbb_maximum": adaptive_bbb_max,
        "logical_status": "strict local strengthening; target remains open",
        "global_status": "continuous terminal strip still open",
    }


if __name__ == "__main__":
    print(json.dumps(validate(), indent=2))

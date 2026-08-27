"""Validate the committed compact spectral-cap cluster audit."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any


HERE = Path(__file__).resolve().parent
SUMMARY = HERE / "spectral_cap_cluster_cover_l055_summary.json"
SOURCE = HERE / "ternary_reconstructed_depth4_g2_top_leaf_bbb_p1_s92_l055.json"


def validate() -> dict[str, Any]:
    summary = json.loads(SUMMARY.read_bytes())
    source_raw = SOURCE.read_bytes()
    if summary["schema"] != "carmenq.spectral-cap-cluster-cover-summary.v1":
        raise ValueError("unexpected compact summary schema")
    canonical_source = source_raw.replace(b"\r\n", b"\n")
    if summary["source"]["sha256"] != hashlib.sha256(canonical_source).hexdigest():
        raise ValueError("compact summary source hash mismatch")
    target = float(summary["target"])
    bound = float(summary["aggregate_upper_bound"])
    if not summary["complete"] or not summary["selected_base_angular_cell_closed"]:
        raise ValueError("compact summary is not complete")
    if not bound < target:
        raise ValueError("compact summary does not close the target")
    cover = summary["cluster_cover"]
    if cover["unresolved_nodes"] != 0:
        raise ValueError("compact summary records unresolved nodes")
    if cover["closed_source_open_cells"] != cover["source_open_cells"]:
        raise ValueError("compact summary does not cover every source-open cell")
    patterns = cover["patterns"]
    if sum(item["closed_clusters"] for item in patterns.values()) != cover[
        "closed_clusters"
    ]:
        raise ValueError("pattern cluster counts do not sum")
    if sum(item["angular_splits"] for item in patterns.values()) != cover[
        "angular_split_nodes"
    ]:
        raise ValueError("pattern split counts do not sum")
    if sum(item["source_open_cells_closed"] for item in patterns.values()) != cover[
        "source_open_cells"
    ]:
        raise ValueError("pattern source-cell counts do not sum")
    if summary["numerical_audit"]["requires_precision_recheck"] is not True:
        raise ValueError("compact summary must retain its precision caveat")
    return summary


def main() -> None:
    summary = validate()
    print(
        json.dumps(
            {
                "schema": summary["schema"],
                "aggregate_upper_bound": summary["aggregate_upper_bound"],
                "margin_below_target": summary["margin_below_target"],
                "closed_source_open_cells": summary["cluster_cover"][
                    "closed_source_open_cells"
                ],
                "closed_clusters": summary["cluster_cover"]["closed_clusters"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

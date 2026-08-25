"""Validate and compact an Ando-strengthened spatial-cover checkpoint.

The large solver checkpoint is regenerable and intentionally untracked. This
script rechecks tree and witness accounting, audits the leading determinant
cells, and records a direct comparison with the committed common-POVM
baseline summary.
"""

from __future__ import annotations

import argparse
import collections
import hashlib
import json
import math
import platform
from pathlib import Path
from typing import Any

import cvxpy as cp
import numpy as np

from summarize_determinant_povm_cover import (
    audit_leading_cells,
    validate_accounting,
)


SCHEMA = "carmenq.ando-instrument-cover-summary.v1"


def validate_ando_accounting(payload: dict[str, Any]) -> dict[str, Any]:
    """Recompute all Ando-specific counters from node records."""

    records = payload["records"]
    count = sum(int(record.get("new_ando_witnesses", 0)) for record in records)
    if int(payload["planar_ando_witness_count"]) != count:
        raise ValueError(
            "planar_ando_witness_count: "
            f"recorded {payload['planar_ando_witness_count']!r}, recomputed {count!r}"
        )
    witness_nodes = [
        record for record in records if int(record.get("new_ando_witnesses", 0)) > 0
    ]
    audits = [
        record["ando_audit"]
        for record in records
        if record.get("ando_audit") is not None
    ]
    violated = [
        audit for audit in audits if int(audit.get("violated_direction_count", 0)) > 0
    ]
    branching = collections.Counter(
        str(record.get("branching_rule", "none")) for record in records
    )
    return {
        "planar_ando_witness_count": count,
        "records_with_new_ando_witnesses": len(witness_nodes),
        "ando_audits": len(audits),
        "audits_with_exact_ando_violation": len(violated),
        "near_ando_margin_nodes": sum(
            bool(record.get("near_ando_margin")) for record in records
        ),
        "positive_ando_split_score_nodes": sum(
            float(record.get("maximum_ando_split_score", 0.0)) > 0.0
            for record in records
        ),
        "maximum_ando_split_score": max(
            (
                float(record.get("maximum_ando_split_score", 0.0))
                for record in records
            ),
            default=0.0,
        ),
        "branching_rules": dict(sorted(branching.items())),
        "witness_nodes": [
            {
                key: record[key]
                for key in (
                    "identifier",
                    "parent",
                    "depth",
                    "parent_bound",
                    "bound",
                    "new_witnesses",
                    "new_ando_witnesses",
                )
            }
            for record in witness_nodes
        ],
    }


def summarize(
    checkpoint: Path,
    baseline_summary: Path,
    leading_cells: int = 20,
) -> dict[str, Any]:
    """Return a compact validated comparison for one strengthened run."""

    raw = checkpoint.read_bytes()
    payload = json.loads(raw)
    baseline = json.loads(baseline_summary.read_text(encoding="utf-8"))
    accounting = validate_accounting(payload)
    ando_accounting = validate_ando_accounting(payload)
    leading = audit_leading_cells(payload, leading_cells)
    target = float(payload["target"])
    base_bound = float(payload["localisation"]["base_result"]["bound"])
    maximum_pending = float(payload["maximum_pending_bound"])
    baseline_pending = float(baseline["run"]["maximum_pending_bound"])
    if int(baseline["run"]["solved_nodes"]) != int(payload["solved_nodes"]):
        raise ValueError("baseline and Ando runs must have equal node budgets")
    if not math.isclose(
        float(baseline["problem"]["target"]),
        target,
        rel_tol=0.0,
        abs_tol=1e-15,
    ):
        raise ValueError("baseline and Ando runs must use the same target")
    improvement = baseline_pending - maximum_pending
    if not bool(payload.get("planar_ando_witnesses")):
        raise ValueError("checkpoint did not enable planar Ando witnesses")
    target_closed = bool(payload["complete"] and maximum_pending < target)
    if target_closed:
        logical_status = "the strengthened cover closes the recorded target"
    elif improvement > 0.0:
        logical_status = (
            "common-instrument positivity gives a strict bound improvement "
            "but remains insufficient at the recorded node budget"
        )
    else:
        logical_status = "no strict equal-budget bound improvement was measured"
    return {
        "schema": SCHEMA,
        "checkpoint": {
            "filename": checkpoint.name,
            "sha256": hashlib.sha256(raw).hexdigest(),
            "bytes": len(raw),
        },
        "baseline": {
            "filename": baseline_summary.name,
            "schema": baseline["schema"],
            "checkpoint_sha256": baseline["checkpoint"]["sha256"],
        },
        "runtime": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "cvxpy": cp.__version__,
            "installed_solvers": cp.installed_solvers(),
        },
        "problem": {
            "support_weight": float(payload["support_weight"]),
            "target": target,
            "base_code": payload["base_code"],
            "top_spectral_cell": bool(payload["top_spectral_cell"]),
            "scope": payload["scope"],
        },
        "configuration": {
            key: payload[key]
            for key in (
                "bound_safety",
                "minimum_width",
                "max_witnesses",
                "max_new_witnesses_per_node",
                "witness_tolerance",
                "determinant_bounds_method",
                "determinant_near_relative_gap",
                "ando_near_relative_gap",
                "ando_split_shortlist",
                "maximum_determinant_branch_streak",
                "max_nodes",
            )
        },
        "run": {
            key: payload[key]
            for key in (
                "solved_nodes",
                "closed_nodes",
                "split_nodes",
                "pending_nodes",
                "unresolved_nodes",
                "maximum_depth",
                "determinant_witness_count",
                "planar_ando_witness_count",
                "maximum_pending_bound",
                "statuses_complete",
                "complete",
            )
        },
        "accounting_audit": accounting,
        "ando_accounting": ando_accounting,
        "leading_pending_determinants": leading,
        "comparison_to_baseline": {
            "baseline_maximum_pending_bound": baseline_pending,
            "ando_maximum_pending_bound": maximum_pending,
            "absolute_bound_improvement": improvement,
            "baseline_closed_nodes": int(baseline["run"]["closed_nodes"]),
            "ando_closed_nodes": int(payload["closed_nodes"]),
            "baseline_pending_nodes": int(baseline["run"]["pending_nodes"]),
            "ando_pending_nodes": int(payload["pending_nodes"]),
            "strict_bound_improvement": bool(improvement > 0.0),
            "tree_dominance": bool(
                improvement > 0.0
                and int(payload["closed_nodes"]) >= int(baseline["run"]["closed_nodes"])
            ),
        },
        "comparison_to_target": {
            "base_localisation_bound": base_bound,
            "base_gap": base_bound - target,
            "remaining_gap": maximum_pending - target,
            "bound_reduction": base_bound - maximum_pending,
            "fraction_of_initial_gap_removed": (
                (base_bound - maximum_pending) / (base_bound - target)
            ),
        },
        "conclusion": {
            "target_closed": target_closed,
            "logical_status": logical_status,
            "next_bottleneck": (
                "the convergent input-cell relaxation and McCormick products, "
                "not missing separate-effect or planar common-instrument positivity"
            ),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("--baseline-summary", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--leading-cells", type=int, default=20)
    args = parser.parse_args()
    result = summarize(
        args.checkpoint,
        args.baseline_summary,
        args.leading_cells,
    )
    encoded = json.dumps(result, indent=2) + "\n"
    if args.output is None:
        print(encoded, end="")
    else:
        args.output.write_text(encoded, encoding="utf-8")


if __name__ == "__main__":
    main()

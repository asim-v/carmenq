"""Validate and compact a common-instrument product-localizer cover.

The full checkpoint contains solver vectors and extended-real JSON values.
This script independently checks the binary tree, recomputes the leaf cover
bound, records the raw-file digest, and emits strict, reviewable JSON.
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


SCHEMA = "carmenq.common-instrument-product-localizer-cover-summary.v1"
ALLOWED_STATUSES = {
    cp.OPTIMAL,
    cp.OPTIMAL_INACCURATE,
    cp.INFEASIBLE,
    cp.INFEASIBLE_INACCURATE,
}


def require_equal(label: str, actual: Any, expected: Any) -> None:
    if actual != expected:
        raise ValueError(f"{label}: recorded {actual!r}, recomputed {expected!r}")


def require_close(label: str, actual: float, expected: float) -> None:
    if actual == expected:
        return
    if not math.isclose(actual, expected, rel_tol=2e-13, abs_tol=2e-13):
        raise ValueError(f"{label}: recorded {actual!r}, recomputed {expected!r}")


def strict_extended_real(value: float) -> tuple[float | None, str]:
    """Return a strict-JSON value and its extended-real class."""

    if math.isfinite(value):
        return value, "finite"
    if value == math.inf:
        return None, "positive-infinity"
    if value == -math.inf:
        return None, "negative-infinity"
    return None, "not-a-number"


def validate_tree(payload: dict[str, Any]) -> dict[str, Any]:
    """Recompute finite-tree accounting and exact parent/child incidence."""

    records = payload["records"]
    pending = payload["pending"]
    unresolved = payload["unresolved"]
    all_nodes = [*records, *pending, *unresolved]
    identifiers = [int(node["identifier"]) for node in all_nodes]
    require_equal("unique identifiers", len(set(identifiers)), len(identifiers))

    records_by_identifier = {
        int(record["identifier"]): record for record in records
    }
    roots = [node for node in all_nodes if node.get("parent") is None]
    require_equal("root count", len(roots), 1)
    require_equal("root identifier", int(roots[0]["identifier"]), 0)

    children: dict[int, list[int]] = collections.defaultdict(list)
    for node in all_nodes:
        parent = node.get("parent")
        if parent is not None:
            if int(parent) not in records_by_identifier:
                raise ValueError(f"node {node['identifier']} has missing parent {parent}")
            children[int(parent)].append(int(node["identifier"]))

    dispositions = collections.Counter(
        str(record.get("disposition", "missing")) for record in records
    )
    for record in records:
        identifier = int(record["identifier"])
        child_count = len(children.get(identifier, []))
        expected = 2 if record.get("disposition") == "split" else 0
        require_equal(f"child count of node {identifier}", child_count, expected)

    statuses = collections.Counter(str(record.get("status")) for record in records)
    statuses_complete = all(status in ALLOWED_STATUSES for status in statuses)
    maximum_pending = max(
        (float(node["parent_bound"]) for node in pending), default=-math.inf
    )
    closed_bounds = [
        float(record["bound"])
        for record in records
        if record.get("disposition") == "closed"
    ]
    maximum_closed = max(closed_bounds, default=-math.inf)
    cover_upper_bound = max(maximum_pending, maximum_closed)

    require_equal("solved_nodes", payload["solved_nodes"], len(records))
    require_equal("closed_nodes", payload["closed_nodes"], dispositions["closed"])
    require_equal("split_nodes", payload["split_nodes"], dispositions["split"])
    require_equal("pending_nodes", payload["pending_nodes"], len(pending))
    require_equal("unresolved_nodes", payload["unresolved_nodes"], len(unresolved))
    require_equal(
        "binary-tree node count",
        len(all_nodes),
        1 + 2 * dispositions["split"],
    )
    require_close(
        "maximum_pending_bound",
        float(payload["maximum_pending_bound"]),
        maximum_pending,
    )
    require_close(
        "maximum_closed_bound",
        float(payload["maximum_closed_bound"]),
        maximum_closed,
    )
    require_close(
        "cover_upper_bound",
        float(payload["cover_upper_bound"]),
        cover_upper_bound,
    )
    require_equal("statuses_complete", payload["statuses_complete"], statuses_complete)
    recomputed_complete = not pending and not unresolved and statuses_complete
    require_equal("complete", payload["complete"], recomputed_complete)

    branch_rules = collections.Counter(
        str(record.get("branching_rule", "none")) for record in records
    )
    maximum_closed_json, maximum_closed_class = strict_extended_real(
        maximum_closed
    )
    return {
        "identifiers_unique": True,
        "root_identifier": 0,
        "all_parent_child_links_valid": True,
        "binary_tree_identity_valid": True,
        "dispositions": dict(sorted(dispositions.items())),
        "statuses": dict(sorted(statuses.items())),
        "branching_rules": dict(sorted(branch_rules.items())),
        "recomputed_maximum_closed_bound": maximum_closed_json,
        "recomputed_maximum_closed_bound_class": maximum_closed_class,
        "recomputed_cover_upper_bound": cover_upper_bound,
    }


def compact_leaves(payload: dict[str, Any]) -> list[dict[str, Any]]:
    """Return the finite closed leaves in descending upper-bound order."""

    leaves = []
    for record in payload["records"]:
        if record.get("disposition") != "closed":
            continue
        bound = float(record["bound"])
        if not math.isfinite(bound):
            continue
        leaves.append(
            {
                "identifier": int(record["identifier"]),
                "parent": (
                    int(record["parent"])
                    if record.get("parent") is not None
                    else None
                ),
                "depth": int(record["depth"]),
                "status": str(record["status"]),
                "bound": bound,
                "raw_value": (
                    float(record["raw_value"])
                    if record.get("raw_value") is not None
                    and math.isfinite(float(record["raw_value"]))
                    else None
                ),
                "maximum_coordinate_width": float(
                    record["maximum_coordinate_width"]
                ),
                "parent_bound": float(record["parent_bound"]),
            }
        )
    return sorted(leaves, key=lambda leaf: (-leaf["bound"], leaf["identifier"]))


def summarize(
    checkpoint: Path,
    baseline_summary: Path,
    sandwich_checkpoint: Path,
) -> dict[str, Any]:
    """Validate one complete localizer checkpoint and compare it to Ando."""

    raw = checkpoint.read_bytes()
    payload = json.loads(raw)
    baseline_raw = baseline_summary.read_bytes()
    baseline = json.loads(baseline_raw)
    sandwich_raw = sandwich_checkpoint.read_bytes()
    sandwich = json.loads(sandwich_raw)
    audit = validate_tree(payload)
    sandwich_audit = validate_tree(sandwich)

    if not payload.get("common_instrument_product_trace_rules"):
        raise ValueError("checkpoint did not enable product trace rules")
    if not payload.get("common_instrument_product_psd_sandwiches"):
        raise ValueError("checkpoint did not enable product PSD sandwiches")
    if sandwich.get("common_instrument_product_trace_rules"):
        raise ValueError("sandwich-only checkpoint unexpectedly enabled trace rules")
    if not sandwich.get("common_instrument_product_psd_sandwiches"):
        raise ValueError("sandwich-only checkpoint did not enable PSD sandwiches")
    require_equal("baseline schema", baseline["schema"], "carmenq.ando-instrument-cover-summary.v1")
    require_close(
        "target agreement",
        float(payload["target"]),
        float(baseline["problem"]["target"]),
    )
    require_close(
        "sandwich target agreement",
        float(payload["target"]),
        float(sandwich["target"]),
    )

    target = float(payload["target"])
    cover_upper = float(audit["recomputed_cover_upper_bound"])
    baseline_upper = float(baseline["run"]["maximum_pending_bound"])
    leaves = compact_leaves(payload)
    target_closed = bool(payload["complete"] and cover_upper < target)
    require_equal("target closure", target_closed, True)

    return {
        "schema": SCHEMA,
        "checkpoint": {
            "filename": checkpoint.name,
            "sha256": hashlib.sha256(raw).hexdigest(),
            "bytes": len(raw),
        },
        "baseline": {
            "filename": baseline_summary.name,
            "sha256": hashlib.sha256(baseline_raw).hexdigest(),
            "schema": baseline["schema"],
            "checkpoint_sha256": baseline["checkpoint"]["sha256"],
        },
        "sandwich_only_checkpoint": {
            "filename": sandwich_checkpoint.name,
            "sha256": hashlib.sha256(sandwich_raw).hexdigest(),
            "bytes": len(sandwich_raw),
        },
        "runtime": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "cvxpy": cp.__version__,
            "installed_solvers": cp.installed_solvers(),
            "recorded_wall_seconds": float(payload["runtime_seconds"]),
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
                "planar_ando_witnesses",
                "common_povm_product_sum_rules",
                "common_instrument_product_trace_rules",
                "common_instrument_product_psd_sandwiches",
                "max_nodes",
            )
        },
        "run": {
            "solved_nodes": int(payload["solved_nodes"]),
            "closed_nodes": int(payload["closed_nodes"]),
            "split_nodes": int(payload["split_nodes"]),
            "pending_nodes": int(payload["pending_nodes"]),
            "unresolved_nodes": int(payload["unresolved_nodes"]),
            "maximum_depth": int(payload["maximum_depth"]),
            "maximum_closed_bound": float(payload["maximum_closed_bound"]),
            "cover_upper_bound": cover_upper,
            "statuses_complete": bool(payload["statuses_complete"]),
            "complete": bool(payload["complete"]),
        },
        "accounting_audit": audit,
        "sandwich_only_ablation": {
            "run": {
                "solved_nodes": int(sandwich["solved_nodes"]),
                "closed_nodes": int(sandwich["closed_nodes"]),
                "split_nodes": int(sandwich["split_nodes"]),
                "pending_nodes": int(sandwich["pending_nodes"]),
                "unresolved_nodes": int(sandwich["unresolved_nodes"]),
                "maximum_pending_bound": float(
                    sandwich["maximum_pending_bound"]
                ),
                "cover_upper_bound": float(sandwich["cover_upper_bound"]),
                "runtime_seconds": float(sandwich["runtime_seconds"]),
                "complete": bool(sandwich["complete"]),
            },
            "accounting_audit": sandwich_audit,
        },
        "finite_closed_leaves": leaves,
        "comparison": {
            "ando_reference_nodes": int(baseline["run"]["solved_nodes"]),
            "ando_maximum_pending_bound": baseline_upper,
            "sandwich_only_nodes": int(sandwich["solved_nodes"]),
            "sandwich_only_maximum_pending_bound": float(
                sandwich["maximum_pending_bound"]
            ),
            "sandwich_only_target_gap": (
                float(sandwich["maximum_pending_bound"]) - target
            ),
            "sandwich_only_complete": bool(sandwich["complete"]),
            "localizer_nodes": int(payload["solved_nodes"]),
            "localizer_cover_upper_bound": cover_upper,
            "absolute_bound_improvement": baseline_upper - cover_upper,
            "target_margin": target - cover_upper,
            "node_reduction_factor": (
                float(baseline["run"]["solved_nodes"])
                / float(payload["solved_nodes"])
            ),
        },
        "conclusion": {
            "target_closed": target_closed,
            "logical_status": (
                "the selected continuous-terminal/Fourier cell is covered below "
                "target by the recorded solver-conditional SDP tree"
            ),
            "global_claim": False,
            "remaining_scope": (
                "extend the same sound localizers across every terminal and "
                "spectral cell required by the global frontier partition"
            ),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("--baseline-summary", type=Path, required=True)
    parser.add_argument("--sandwich-checkpoint", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = summarize(
        args.checkpoint,
        args.baseline_summary,
        args.sandwich_checkpoint,
    )
    encoded = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if args.output is None:
        print(encoded, end="")
    else:
        args.output.write_text(encoded, encoding="utf-8")


if __name__ == "__main__":
    main()

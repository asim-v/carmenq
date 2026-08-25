"""Validate and compact a determinant-witness spatial-cover checkpoint.

The full branch-and-bound checkpoints are intentionally not committed: a
1,000-node run contains solver vectors and can exceed twenty megabytes.  This
script checks the combinatorial accounting, recomputes the determinant bounds
of the leading open cells, and emits a compact, reviewable JSON result.
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

from ternary_bilinear_instrument_input_cover import (
    determinant_interval,
    determinant_vertex_bounds,
)


SCHEMA = "carmenq.determinant-povm-cover-summary.v1"


def _require_equal(label: str, actual: Any, expected: Any) -> None:
    if actual != expected:
        raise ValueError(f"{label}: recorded {actual!r}, recomputed {expected!r}")


def _require_close(label: str, actual: float, expected: float) -> None:
    if not math.isclose(actual, expected, rel_tol=2e-13, abs_tol=2e-13):
        raise ValueError(f"{label}: recorded {actual!r}, recomputed {expected!r}")


def validate_accounting(payload: dict[str, Any]) -> dict[str, Any]:
    """Recompute the checkpoint's finite tree and witness accounting."""

    records = payload["records"]
    pending = payload["pending"]
    unresolved = payload["unresolved"]
    dispositions = collections.Counter(
        str(record.get("disposition", "missing")) for record in records
    )
    branch_rules = collections.Counter(
        str(record.get("branching_rule", "none")) for record in records
    )
    witness_count = sum(int(record.get("new_witnesses", 0)) for record in records)
    maximum_pending_bound = max(
        (float(node["parent_bound"]) for node in pending), default=-math.inf
    )

    _require_equal("solved_nodes", payload["solved_nodes"], len(records))
    _require_equal("closed_nodes", payload["closed_nodes"], dispositions["closed"])
    _require_equal("split_nodes", payload["split_nodes"], dispositions["split"])
    _require_equal("pending_nodes", payload["pending_nodes"], len(pending))
    _require_equal("unresolved_nodes", payload["unresolved_nodes"], len(unresolved))
    _require_equal(
        "determinant_witness_count",
        payload["determinant_witness_count"],
        witness_count,
    )
    _require_close(
        "maximum_pending_bound",
        float(payload["maximum_pending_bound"]),
        maximum_pending_bound,
    )

    # Starting from one root, every split consumes one pending node and adds
    # two children; every closed or unresolved node only consumes one.
    expected_pending = 1 + 2 * dispositions["split"] - len(records)
    _require_equal("binary-tree pending count", len(pending), expected_pending)
    return {
        "dispositions": dict(sorted(dispositions.items())),
        "branching_rules": dict(sorted(branch_rules.items())),
        "records_with_new_witnesses": sum(
            int(record.get("new_witnesses", 0)) > 0 for record in records
        ),
        "sign_definite_audits": sum(
            bool((record.get("determinant_audit") or {}).get("sign_definite"))
            for record in records
        ),
        "robust_witness_nodes": sum(
            int(
                (record.get("determinant_audit") or {}).get(
                    "robust_witness_count", 0
                )
            )
            > 0
            for record in records
        ),
        "near_margin_nodes": sum(
            bool(record.get("near_determinant_margin")) for record in records
        ),
    }


def audit_leading_cells(
    payload: dict[str, Any], count: int = 20
) -> dict[str, Any]:
    """Recompute dependency-prone and exhaustive bounds on leading cells."""

    leading = sorted(
        payload["pending"],
        key=lambda node: float(node["parent_bound"]),
        reverse=True,
    )[:count]
    cells: list[dict[str, Any]] = []
    for node in leading:
        lower = np.asarray(node["lower"], dtype=float)
        upper = np.asarray(node["upper"], dtype=float)
        dependency = determinant_interval(lower, upper)
        vertices = determinant_vertex_bounds(lower, upper)
        if vertices.lower > 0.0:
            sign = "positive"
        elif vertices.upper < 0.0:
            sign = "negative"
        else:
            sign = "indefinite"
        cells.append(
            {
                "identifier": int(node["identifier"]),
                "parent": int(node["parent"]),
                "depth": int(node["depth"]),
                "parent_bound": float(node["parent_bound"]),
                "ordinary_interval": [dependency.lower, dependency.upper],
                "exhaustive_vertex_enclosure": [vertices.lower, vertices.upper],
                "vertex_sign": sign,
                "inherited_witnesses": len(node.get("determinant_witnesses", [])),
            }
        )
    return {
        "cell_count": len(cells),
        "ordinary_interval_sign_indefinite": sum(
            cell["ordinary_interval"][0] <= 0.0 <= cell["ordinary_interval"][1]
            for cell in cells
        ),
        "vertex_sign_counts": dict(
            sorted(collections.Counter(cell["vertex_sign"] for cell in cells).items())
        ),
        "minimum_certified_absolute_determinant": min(
            (
                min(abs(bound) for bound in cell["exhaustive_vertex_enclosure"])
                for cell in cells
                if cell["vertex_sign"] != "indefinite"
            ),
            default=0.0,
        ),
        "cells": cells,
    }


def summarize(checkpoint: Path, leading_cells: int = 20) -> dict[str, Any]:
    """Return a compact validated summary of one full checkpoint."""

    raw = checkpoint.read_bytes()
    payload = json.loads(raw)
    accounting = validate_accounting(payload)
    leading = audit_leading_cells(payload, leading_cells)
    target = float(payload["target"])
    maximum_pending = float(payload["maximum_pending_bound"])
    base_bound = float(payload["localisation"]["base_result"]["bound"])
    closed_target = bool(payload["complete"] and maximum_pending < target)
    return {
        "schema": SCHEMA,
        "checkpoint": {
            "filename": checkpoint.name,
            "sha256": hashlib.sha256(raw).hexdigest(),
            "bytes": len(raw),
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
                "maximum_pending_bound",
                "statuses_complete",
                "complete",
            )
        },
        "accounting_audit": accounting,
        "leading_pending_determinants": leading,
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
            "target_closed": closed_target,
            "logical_status": (
                "selected cell closed at target"
                if closed_target
                else "determinant-POVM cuts valid but insufficient at node budget"
            ),
            "next_relaxation": (
                "coupled matrix-valued common-instrument positivity, not additional "
                "scalar common-POVM cuts"
            ),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--leading-cells", type=int, default=20)
    args = parser.parse_args()
    result = summarize(args.checkpoint, args.leading_cells)
    encoded = json.dumps(result, indent=2) + "\n"
    if args.output is None:
        print(encoded, end="")
    else:
        args.output.write_text(encoded, encoding="utf-8")


if __name__ == "__main__":
    main()

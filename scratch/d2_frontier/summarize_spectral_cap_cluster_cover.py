"""Validate and compact the adaptive spectral-cap cluster checkpoint."""

from __future__ import annotations

import argparse
import collections
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import cvxpy as cp

from spectral_product_localizer_batch import enclosing_scaled_cap, pattern_code


SCHEMA = "carmenq.spectral-cap-cluster-cover-summary.v1"
ACCEPTED = {
    cp.OPTIMAL,
    cp.OPTIMAL_INACCURATE,
    cp.INFEASIBLE,
    cp.INFEASIBLE_INACCURATE,
}


def closure_result(node: dict[str, Any]) -> tuple[float, str, str]:
    method = str(node["closure_method"])
    if method == "base-relaxation":
        result = node["localisation"]["base_result"]
        bound = -math.inf if result["bound"] is None else float(result["bound"])
        return bound, str(result["status"]), method
    cover = node.get("leaf_cover", node.get("root_cover"))
    if cover is None or not cover.get("complete") or not cover.get("target_closed"):
        raise ValueError(f"closed node {node['identifier']} lacks a complete cover")
    if cover["cover_upper_bound"] is None:
        if cover.get("cover_upper_bound_class") == "negative-infinity":
            return -math.inf, "infeasible", method
        raise ValueError(f"closed cover {node['identifier']} lacks a finite bound")
    statuses = list(cover["statuses"])
    status = statuses[0] if len(statuses) == 1 else "mixed"
    return float(cover["cover_upper_bound"]), status, method


def solve_statuses(node: dict[str, Any]) -> collections.Counter[str]:
    result: collections.Counter[str] = collections.Counter()
    localisation = node["localisation"]
    result[str(localisation["base_result"]["status"])] += 1
    for support in localisation.get("supports", []):
        result[str(support["positive"]["status"])] += 1
        result[str(support["negative"]["status"])] += 1
    for key in ("root_cover", "leaf_cover"):
        cover = node.get(key)
        if cover is not None:
            result.update({str(status): int(count) for status, count in cover["statuses"].items()})
    return result


def summarize(
    source_path: Path,
    checkpoint_path: Path,
) -> dict[str, Any]:
    source_raw = source_path.read_bytes()
    checkpoint_raw = checkpoint_path.read_bytes()
    source = json.loads(source_raw)
    checkpoint = json.loads(checkpoint_raw)
    if checkpoint.get("schema") != "carmenq.spectral-cap-cluster-cover.v1":
        raise ValueError("unexpected cluster checkpoint schema")
    if checkpoint["source"]["sha256"] != hashlib.sha256(source_raw).hexdigest():
        raise ValueError("checkpoint source hash mismatch")
    if not source.get("statuses_complete"):
        raise ValueError("source spectral cover has unresolved statuses")
    target = float(source["target"])
    grids = tuple(int(value) for value in source["separator_grids"])
    nodes = {int(key): value for key, value in checkpoint["nodes"].items()}
    if checkpoint["pending"]:
        raise ValueError("cluster checkpoint still has pending nodes")

    source_open = {
        index
        for index, cell in enumerate(source["cells"])
        if str(cell["status"]) in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}
        and float(cell["bound"]) >= target
    }
    root_nodes = [node for node in nodes.values() if node["parent"] is None]
    root_counts: collections.Counter[int] = collections.Counter(
        int(index) for node in root_nodes for index in node["source_indices"]
    )
    if set(root_counts) != source_open or any(count != 1 for count in root_counts.values()):
        raise ValueError("root clusters do not partition the source-open cells")

    cap_audits = 0
    for node in nodes.values():
        pattern = tuple(node["pattern"])
        for position, (branch, raw_caps) in enumerate(
            zip(pattern, node["caps"], strict=True)
        ):
            if branch != "bloch":
                if raw_caps is not None:
                    raise ValueError("scalar cluster branch carries a cap")
                continue
            indices = (
                tuple(int(value) for value in raw_caps)
                if isinstance(raw_caps, list)
                else (int(raw_caps),)
            )
            _, audit = enclosing_scaled_cap(grids[position], indices)
            if float(audit["cosine"]) <= 0.0:
                raise ValueError("cluster cap is not in an open hemisphere")
            cap_audits += 1
        if node.get("disposition") == "angular-split":
            children = [nodes[int(identifier)] for identifier in node["children"]]
            if any(int(child["parent"]) != int(node["identifier"]) for child in children):
                raise ValueError("child cluster has the wrong parent")
            child_counts: collections.Counter[int] = collections.Counter(
                int(index) for child in children for index in child["source_indices"]
            )
            if (
                set(child_counts) != {int(index) for index in node["source_indices"]}
                or any(count != 1 for count in child_counts.values())
            ):
                raise ValueError("angular children do not partition their parent")

    closed_nodes = [node for node in nodes.values() if node.get("disposition") == "closed"]
    if any(
        node.get("disposition") not in {"closed", "angular-split"}
        for node in nodes.values()
    ):
        raise ValueError("checkpoint contains an unresolved disposition")
    leaf_counts: collections.Counter[int] = collections.Counter(
        int(index) for node in closed_nodes for index in node["source_indices"]
    )
    if set(leaf_counts) != source_open or any(count != 1 for count in leaf_counts.values()):
        raise ValueError("closed clusters do not partition the source-open cells")

    status_counts: collections.Counter[str] = collections.Counter()
    cluster_bounds: list[tuple[float, dict[str, Any], str]] = []
    methods: collections.Counter[str] = collections.Counter()
    patterns: dict[str, dict[str, int]] = collections.defaultdict(
        lambda: {
            "closed_clusters": 0,
            "angular_splits": 0,
            "source_open_cells_closed": 0,
            "maximum_closed_cluster_size": 0,
        }
    )
    for node in nodes.values():
        status_counts.update(solve_statuses(node))
        code = str(node["pattern_code"])
        if code != pattern_code(tuple(node["pattern"])):
            raise ValueError("stored pattern code is inconsistent")
        if node["disposition"] == "angular-split":
            patterns[code]["angular_splits"] += 1
            continue
        bound, _, method = closure_result(node)
        if not bound < target:
            raise ValueError(f"closed cluster {node['identifier']} reaches the target")
        methods[method] += 1
        patterns[code]["closed_clusters"] += 1
        size = int(node["source_open_cell_count"])
        patterns[code]["source_open_cells_closed"] += size
        patterns[code]["maximum_closed_cluster_size"] = max(
            patterns[code]["maximum_closed_cluster_size"], size
        )
        cluster_bounds.append((bound, node, method))
    if any(status not in ACCEPTED for status in status_counts):
        raise ValueError("checkpoint contains an unaccepted solver status")

    source_closed: list[tuple[float, int, dict[str, Any]]] = []
    for index, cell in enumerate(source["cells"]):
        status = str(cell["status"])
        if status in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}:
            bound = float(cell["bound"])
            if index not in source_open:
                if not bound < target:
                    raise ValueError("source cell is neither open nor strictly closed")
                source_closed.append((bound, index, cell))
        elif status not in {cp.INFEASIBLE, cp.INFEASIBLE_INACCURATE}:
            raise ValueError("source cover contains an unaccepted status")

    cluster_maximum = max(cluster_bounds, key=lambda item: item[0])
    source_maximum = max(source_closed, key=lambda item: item[0])
    aggregate = max(cluster_maximum[0], source_maximum[0])
    if not aggregate < target:
        raise ValueError("aggregate base-angular-cell bound reaches the target")

    monotonicity: list[dict[str, Any]] = []
    for node in nodes.values():
        if node["parent"] is None or "root_cover" not in node:
            continue
        parent = nodes[int(node["parent"])]
        if "root_cover" not in parent:
            continue
        child_bound = node["root_cover"]["cover_upper_bound"]
        parent_bound = parent["root_cover"]["cover_upper_bound"]
        if child_bound is None or parent_bound is None:
            continue
        excess = float(child_bound) - float(parent_bound)
        if excess > 1e-8:
            monotonicity.append(
                {
                    "parent": int(parent["identifier"]),
                    "child": int(node["identifier"]),
                    "excess": excess,
                }
            )

    worst_clusters = sorted(cluster_bounds, key=lambda item: item[0], reverse=True)[:10]
    return {
        "schema": SCHEMA,
        "scope": (
            "one fixed terminal box and one fixed base Fourier angular cell; "
            "complete adaptive cover of its reconstructed spectral subcells"
        ),
        "epistemic_status": (
            "solver-conditional numerical certificate; not a solver-independent "
            "interval or rational-dual certificate"
        ),
        "source": {
            "path": str(source_path),
            "sha256": hashlib.sha256(source_raw).hexdigest(),
            "cell_count": len(source["cells"]),
            "statuses_complete": bool(source["statuses_complete"]),
        },
        "checkpoint": {
            "path": str(checkpoint_path),
            "sha256": hashlib.sha256(checkpoint_raw).hexdigest(),
            "schema": checkpoint["schema"],
        },
        "target": target,
        "complete": True,
        "selected_base_angular_cell_closed": True,
        "aggregate_upper_bound": aggregate,
        "margin_below_target": target - aggregate,
        "limiting_component": "source-preclosed-spectral-cell",
        "limiting_source_cell": {
            "source_index": int(source_maximum[1]),
            "bound": float(source_maximum[0]),
            "source_cell": int(source_maximum[2]["source_cell"]),
            "branches": source_maximum[2]["branches"],
            "caps": source_maximum[2]["caps"],
            "status": source_maximum[2]["status"],
        },
        "cluster_cover": {
            "source_open_cells": len(source_open),
            "cluster_nodes": len(nodes),
            "root_face_clusters": len(root_nodes),
            "angular_split_nodes": sum(
                node["disposition"] == "angular-split" for node in nodes.values()
            ),
            "closed_clusters": len(closed_nodes),
            "closed_source_open_cells": sum(
                int(node["source_open_cell_count"]) for node in closed_nodes
            ),
            "unresolved_nodes": 0,
            "maximum_cluster_upper_bound": cluster_maximum[0],
            "maximum_cluster_margin": target - cluster_maximum[0],
            "cap_containment_audits": cap_audits,
            "closure_methods": dict(sorted(methods.items())),
            "patterns": dict(sorted(patterns.items())),
            "solver_statuses": dict(sorted(status_counts.items())),
        },
        "numerical_audit": {
            "monotonicity_inversions_above_1e-8": len(monotonicity),
            "maximum_monotonicity_inversion": max(
                (item["excess"] for item in monotonicity), default=0.0
            ),
            "largest_inversions": sorted(
                monotonicity, key=lambda item: item["excess"], reverse=True
            )[:10],
            "requires_precision_recheck": True,
        },
        "worst_closed_clusters": [
            {
                "bound": bound,
                "identifier": int(node["identifier"]),
                "pattern": node["pattern_code"],
                "source_open_cell_count": int(node["source_open_cell_count"]),
                "cartesian_child_cell_count": int(node["cartesian_child_cell_count"]),
                "closure_method": method,
            }
            for bound, node, method in worst_clusters
        ],
        "interpretation": {
            "proved_at_this_numerical_level": (
                "the selected base angular cell is strictly below 0.758 in the "
                "stated common-instrument relaxation"
            ),
            "not_proved": (
                "the complete terminal leaf, the complete terminal strip, or a "
                "solver-independent theorem"
            ),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = summarize(args.source, args.checkpoint)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()

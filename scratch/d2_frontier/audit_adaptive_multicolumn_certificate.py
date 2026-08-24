"""Audit or numerically replay an adaptive multicolumn certificate."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from adaptive_multicolumn_branch_tree import child_contraction
from fourier_behavior_cap_cover import cube_face_caps
from multicolumn_contraction_cell_cover import run_cover


def _close(left: float, right: float, tolerance: float) -> bool:
    return abs(left - right) <= tolerance * max(1.0, abs(left), abs(right))


def audit_certificate(payload: dict[str, object], tolerance: float = 5e-6) -> dict[str, object]:
    """Check topology, exhaustive branch labels, bounds, and closure decisions."""
    target = float(payload["target"])
    grid = int(payload["contraction_grid"])
    cap_count = len(cube_face_caps(grid))
    expected_keys = {
        ("scalar-positive", None),
        ("scalar-negative", None),
        *(("bloch", index) for index in range(cap_count)),
    }
    records = list(payload["expanded_nodes"])
    by_identifier = {int(record["id"]): record for record in records}
    if len(by_identifier) != len(records) or 0 not in by_identifier:
        raise ValueError("expanded node identifiers are not a unique rooted set")
    if by_identifier[0]["parent"] is not None:
        raise ValueError("node zero is not the root")

    referenced_children: dict[int, tuple[int, float]] = {}
    closed_leaves = 0
    maximum_depth = 0
    for record in records:
        identifier = int(record["id"])
        maximum_depth = max(maximum_depth, int(record["depth"]))
        branches = list(record["branches"])
        keys = {(str(row["branch"]), row["cap"]) for row in branches}
        if keys != expected_keys or len(branches) != len(expected_keys):
            raise ValueError(f"node {identifier} does not contain the exhaustive branch cover")

        finite_bounds = []
        open_count = 0
        for row in branches:
            bound = row["bound"]
            is_closed = bound is None or float(bound) < target
            if bool(row["closed"]) != is_closed:
                raise ValueError(f"node {identifier} has an inconsistent closure flag")
            if bound is not None:
                finite_bounds.append(float(bound))
            child_id = row["child_id"]
            if is_closed:
                closed_leaves += 1
                if child_id is not None:
                    raise ValueError(f"closed branch of node {identifier} has a child")
            else:
                open_count += 1
                if child_id is None:
                    raise ValueError(f"open branch of node {identifier} has no child")
                child = int(child_id)
                if child in referenced_children:
                    raise ValueError(f"child {child} has more than one incoming branch")
                referenced_children[child] = (identifier, float(bound))

        if int(record["closed_children"]) != len(branches) - open_count:
            raise ValueError(f"node {identifier} has a wrong closed-child count")
        if int(record["open_children"]) != open_count:
            raise ValueError(f"node {identifier} has a wrong open-child count")
        maximum = max(finite_bounds, default=-math.inf)
        if not _close(maximum, float(record["maximum_child_bound"]), tolerance):
            raise ValueError(f"node {identifier} has a wrong maximum child bound")

    pending = {int(node["id"]): node for node in payload["pending_nodes"]}
    all_nonroot = (set(by_identifier) | set(pending)) - {0}
    if set(referenced_children) != all_nonroot:
        raise ValueError("the recorded edges do not equal the expanded and pending nodes")
    for child, (parent, incoming_bound) in referenced_children.items():
        node = by_identifier.get(child, pending.get(child))
        if int(node["parent"]) != parent:
            raise ValueError(f"child {child} has an inconsistent parent")
        if child in by_identifier:
            if not _close(float(node["source_bound"]), incoming_bound, tolerance):
                raise ValueError(f"child {child} has an inconsistent source bound")
        elif not _close(float(node["bound"]), incoming_bound, tolerance):
            raise ValueError(f"pending child {child} has an inconsistent bound")

    complete = not pending
    if bool(payload["certificate_complete"]) != complete:
        raise ValueError("certificate completeness does not match the pending-node set")
    if int(payload["closed_leaf_count"]) != closed_leaves:
        raise ValueError("the global closed-leaf count is inconsistent")
    if int(payload["open_leaf_count"]) != len(pending):
        raise ValueError("the global open-leaf count is inconsistent")
    return {
        "certificate_complete": complete,
        "expanded_nodes": len(records),
        "closed_leaves": closed_leaves,
        "open_leaves": len(pending),
        "maximum_depth": maximum_depth,
        "branches_per_expansion": len(expected_keys),
    }


def _paths(payload: dict[str, object]) -> dict[int, tuple[dict[str, object], ...]]:
    records = {int(record["id"]): record for record in payload["expanded_nodes"]}
    cache: dict[int, tuple[dict[str, object], ...]] = {0: ()}

    def visit(identifier: int) -> tuple[dict[str, object], ...]:
        if identifier in cache:
            return cache[identifier]
        record = records[identifier]
        parent_id = int(record["parent"])
        parent = records[parent_id]
        incoming = next(
            row for row in parent["branches"] if row["child_id"] == identifier
        )
        contraction = child_contraction(
            list(parent["separator_coefficients"]),
            incoming,
            int(payload["contraction_grid"]),
        )
        cache[identifier] = (*visit(parent_id), contraction)
        return cache[identifier]

    for identifier in records:
        visit(identifier)
    return cache


def replay_certificate(payload: dict[str, object], tolerance: float = 5e-6) -> None:
    """Re-solve every recorded conic branch and compare its upper bound."""
    paths = _paths(payload)
    for record in payload["expanded_nodes"]:
        identifier = int(record["id"])
        replay = run_cover(
            np.asarray(record["separator_coefficients"], dtype=float),
            int(payload["plane_cells"]),
            int(payload["plane_index"]),
            int(payload["face_grid"]),
            int(payload["sphere_index"]),
            int(payload["pair_cap_index"]),
            int(payload["contraction_grid"]),
            paths[identifier],
            target=float(payload["target"]),
            pair_branch=str(payload.get("pair_branch", "bloch")),
        )
        stored = {
            (str(row["branch"]), row["cap"]): row for row in record["branches"]
        }
        for row in replay["branches"]:
            expected = stored[(str(row["branch"]), row["cap"])]["bound"]
            observed = float(row["bound"])
            if expected is None:
                if math.isfinite(observed):
                    raise ValueError(f"node {identifier} changed infeasibility status")
            elif not _close(observed, float(expected), tolerance):
                raise ValueError(f"node {identifier} changed a stored conic bound")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("certificate", type=Path)
    parser.add_argument("--recompute", action="store_true")
    parser.add_argument("--tolerance", type=float, default=5e-6)
    args = parser.parse_args()
    payload = json.loads(args.certificate.read_text(encoding="utf-8"))
    summary = audit_certificate(payload, args.tolerance)
    if args.recompute:
        replay_certificate(payload, args.tolerance)
        summary["numerically_replayed"] = True
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

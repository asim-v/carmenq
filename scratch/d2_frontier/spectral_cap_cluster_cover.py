"""Adaptive angular-cluster cover of one reconstructed spectral frontier."""

from __future__ import annotations

import argparse
import collections
import gc
import hashlib
import json
from pathlib import Path
import platform
from typing import Any

import cvxpy as cp
import numpy as np

from fourier_behavior_cap_cover import cube_face_caps
from spectral_product_localizer_batch import (
    build_localisation_oracle,
    build_product_oracle,
    cover_localised_cell,
    enclosing_scaled_cap,
    localise_cell,
    pattern_code,
    source_open_cells,
)


SCHEMA = "carmenq.spectral-cap-cluster-cover.v1"


def cap_sets(
    cells: list[dict[str, Any]],
    pattern: tuple[str, ...],
) -> tuple[int | tuple[int, ...] | None, ...]:
    result: list[int | tuple[int, ...] | None] = []
    for position, branch in enumerate(pattern):
        if branch != "bloch":
            result.append(None)
            continue
        indices = tuple(sorted({int(cell["caps"][position]) for cell in cells}))
        result.append(indices[0] if len(indices) == 1 else indices)
    return tuple(result)


def cartesian_size(caps: tuple[int | tuple[int, ...] | None, ...]) -> int:
    size = 1
    for indices in caps:
        if isinstance(indices, tuple):
            size *= len(indices)
    return size


def make_node(
    identifier: int,
    parent: int | None,
    depth: int,
    cells: list[dict[str, Any]],
    pattern: tuple[str, ...],
) -> dict[str, Any]:
    caps = cap_sets(cells, pattern)
    return {
        "identifier": int(identifier),
        "parent": parent,
        "depth": int(depth),
        "pattern": list(pattern),
        "pattern_code": pattern_code(pattern),
        "source_indices": [int(cell["source_index"]) for cell in cells],
        "source_open_cell_count": len(cells),
        "maximum_source_bound": max(float(cell["source_bound"]) for cell in cells),
        "caps": [list(item) if isinstance(item, tuple) else item for item in caps],
        "cartesian_child_cell_count": cartesian_size(caps),
    }


def node_cells(
    node: dict[str, Any],
    open_by_index: dict[int, dict[str, Any]],
) -> list[dict[str, Any]]:
    return [open_by_index[int(index)] for index in node["source_indices"]]


def initial_face_nodes(
    open_cells: list[dict[str, Any]],
    grids: tuple[int, ...],
) -> list[dict[str, Any]]:
    groups: dict[tuple[tuple[str, ...], tuple[int | None, ...]], list[dict[str, Any]]]
    groups = collections.defaultdict(list)
    for cell in open_cells:
        pattern = tuple(cell["branches"])
        faces = tuple(
            (
                int(cell["caps"][position]) // grids[position] ** 2
                if branch == "bloch"
                else None
            )
            for position, branch in enumerate(pattern)
        )
        groups[(pattern, faces)].append(cell)
    ordered = sorted(
        groups.items(),
        key=lambda item: (
            -max(float(cell["source_bound"]) for cell in item[1]),
            pattern_code(item[0][0]),
            str(item[0][1]),
        ),
    )
    return [
        make_node(identifier, None, 0, cells, pattern)
        for identifier, ((pattern, _), cells) in enumerate(ordered)
    ]


def split_dimension(
    source: dict[str, Any],
    node: dict[str, Any],
) -> tuple[int, set[int], set[int], dict[str, Any]]:
    pattern = tuple(node["pattern"])
    grids = tuple(int(value) for value in source["separator_grids"])
    choices: list[tuple[float, int, tuple[int, ...], dict[str, Any]]] = []
    for position, (branch, raw_indices) in enumerate(
        zip(pattern, node["caps"], strict=True)
    ):
        if branch != "bloch" or not isinstance(raw_indices, list) or len(raw_indices) < 2:
            continue
        indices = tuple(int(index) for index in raw_indices)
        _, audit = enclosing_scaled_cap(grids[position], indices)
        choices.append((float(audit["angular_radius"]), position, indices, audit))
    if not choices:
        raise ValueError("cluster has no splittable Bloch cap set")
    _, position, indices, audit = max(choices, key=lambda item: (item[0], len(item[2])))
    child_caps = cube_face_caps(grids[position])
    normals = {index: np.asarray(child_caps[index][0], dtype=float) for index in indices}
    first, second = min(
        (
            (first, second)
            for offset, first in enumerate(indices)
            for second in indices[offset + 1 :]
        ),
        key=lambda pair: float(normals[pair[0]] @ normals[pair[1]]),
    )
    left: set[int] = set()
    right: set[int] = set()
    for index in indices:
        if float(normals[index] @ normals[first]) >= float(
            normals[index] @ normals[second]
        ):
            left.add(index)
        else:
            right.add(index)
    if not left or not right:
        raise RuntimeError("farthest-centre cap split produced an empty child")
    return position, left, right, {
        "position": position,
        "parent_cap": audit,
        "seed_indices": [first, second],
        "left_indices": sorted(left),
        "right_indices": sorted(right),
    }


def synthetic_cell(node: dict[str, Any]) -> dict[str, Any]:
    return {
        "source_index": int(node["identifier"]),
        "source_cell": -1,
        "branches": tuple(node["pattern"]),
        "caps": tuple(
            tuple(int(value) for value in item) if isinstance(item, list) else item
            for item in node["caps"]
        ),
        "source_status": "cluster",
        "source_bound": float(node["maximum_source_bound"]),
        "source_audit": 0.0,
        "source_return": 0.0,
    }


def write_payload(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


def update_summary(payload: dict[str, Any], source_open_count: int) -> None:
    nodes = payload["nodes"]
    dispositions = collections.Counter(
        node.get("disposition", "pending") for node in nodes.values()
    )
    closed_indices: set[int] = set()
    unresolved_indices: set[int] = set()
    for node in nodes.values():
        if node.get("disposition") == "closed":
            closed_indices.update(int(index) for index in node["source_indices"])
        elif node.get("disposition") in {"leaf-unresolved", "solver-unresolved"}:
            unresolved_indices.update(int(index) for index in node["source_indices"])
    pending_indices = {
        int(index)
        for identifier in payload["pending"]
        for index in nodes[str(identifier)]["source_indices"]
    }
    covered = closed_indices | unresolved_indices | pending_indices
    payload["summary"] = {
        "source_open_cells": int(source_open_count),
        "cluster_nodes": len(nodes),
        "solved_cluster_nodes": sum(
            node.get("disposition", "pending") != "pending" for node in nodes.values()
        ),
        "closed_cluster_nodes": int(dispositions["closed"]),
        "angular_split_nodes": int(dispositions["angular-split"]),
        "pending_cluster_nodes": len(payload["pending"]),
        "unresolved_cluster_nodes": int(
            dispositions["leaf-unresolved"] + dispositions["solver-unresolved"]
        ),
        "closed_source_open_cells": len(closed_indices),
        "pending_source_open_cells": len(pending_indices - closed_indices),
        "unresolved_source_open_cells": len(unresolved_indices - closed_indices),
        "source_partition_accounted": len(covered) == source_open_count,
        "selected_base_angular_cell_closed": bool(
            len(closed_indices) == source_open_count
            and not pending_indices
            and not unresolved_indices
        ),
    }


def initial_payload(
    source_path: Path,
    source_raw: bytes,
    source: dict[str, Any],
    open_cells: list[dict[str, Any]],
    coordinate_safety: float,
    bound_safety: float,
    minimum_width: float,
    leaf_max_nodes: int,
) -> dict[str, Any]:
    grids = tuple(int(value) for value in source["separator_grids"])
    roots = initial_face_nodes(open_cells, grids)
    payload = {
        "schema": SCHEMA,
        "source": {
            "path": str(source_path),
            "sha256": hashlib.sha256(source_raw).hexdigest(),
            "target": float(source["target"]),
            "source_open_cells": len(open_cells),
        },
        "configuration": {
            "initial_partition": "branch-pattern and cube-face tuple",
            "angular_split": "largest parent radius, farthest-centre bipartition",
            "coordinate_safety": coordinate_safety,
            "bound_safety": bound_safety,
            "minimum_width": minimum_width,
            "leaf_max_nodes": leaf_max_nodes,
            "state_choi_psd": True,
            "state_choi_ppt": True,
        },
        "runtime": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "cvxpy": cp.__version__,
            "installed_solvers": cp.installed_solvers(),
        },
        "group_builds": {},
        "next_identifier": len(roots),
        "nodes": {str(node["identifier"]): node for node in roots},
        "pending": [int(node["identifier"]) for node in roots],
        "summary": {},
    }
    update_summary(payload, len(open_cells))
    return payload


def process_pattern(
    source: dict[str, Any],
    payload: dict[str, Any],
    pattern: tuple[str, ...],
    open_by_index: dict[int, dict[str, Any]],
    node_budget: int,
    coordinate_safety: float,
    bound_safety: float,
    minimum_width: float,
    leaf_max_nodes: int,
    output: Path,
) -> int:
    localisation_oracle, localisation_caps, box, localisation_build = (
        build_localisation_oracle(source, pattern)
    )
    (
        product_oracle,
        product_caps,
        product_box,
        lower,
        upper,
        purity,
        product_build,
    ) = build_product_oracle(source, pattern)
    code = pattern_code(pattern)
    payload["group_builds"][code] = {
        "localisation": localisation_build,
        "product": product_build,
    }
    solved = 0
    while solved < node_budget:
        candidates = [
            payload["nodes"][str(identifier)]
            for identifier in payload["pending"]
            if tuple(payload["nodes"][str(identifier)]["pattern"]) == pattern
        ]
        if not candidates:
            break
        node = max(
            candidates,
            key=lambda item: (
                float(item["maximum_source_bound"]),
                int(item["source_open_cell_count"]),
                -int(item["identifier"]),
            ),
        )
        identifier = int(node["identifier"])
        payload["pending"].remove(identifier)
        cell = synthetic_cell(node)
        localisation = localise_cell(
            source,
            cell,
            localisation_oracle,
            localisation_caps,
            box,
            coordinate_safety,
        )
        node["localisation"] = localisation
        cover: dict[str, Any] | None = None
        if localisation["status"] == "base-closed":
            node["disposition"] = "closed"
            node["closure_method"] = "base-relaxation"
        elif localisation["status"] != "localized":
            node["disposition"] = "solver-unresolved"
        else:
            cover = cover_localised_cell(
                source,
                cell,
                localisation,
                product_oracle,
                product_caps,
                product_box,
                lower,
                upper,
                purity,
                1,
                bound_safety,
                minimum_width,
            )
            node["root_cover"] = cover
            if cover["target_closed"]:
                node["disposition"] = "closed"
                node["closure_method"] = "cluster-root-state-choi-ppt"
            elif int(node["source_open_cell_count"]) > 1:
                position, left, right, split = split_dimension(source, node)
                cells = node_cells(node, open_by_index)
                child_groups = [
                    [cell for cell in cells if int(cell["caps"][position]) in side]
                    for side in (left, right)
                ]
                if any(not child for child in child_groups):
                    raise RuntimeError("angular split lost one child cluster")
                children = []
                for child_cells in child_groups:
                    child_identifier = int(payload["next_identifier"])
                    payload["next_identifier"] += 1
                    child = make_node(
                        child_identifier,
                        identifier,
                        int(node["depth"]) + 1,
                        child_cells,
                        pattern,
                    )
                    payload["nodes"][str(child_identifier)] = child
                    payload["pending"].append(child_identifier)
                    children.append(child_identifier)
                node["disposition"] = "angular-split"
                node["angular_split"] = split
                node["children"] = children
            else:
                leaf_cover = cover_localised_cell(
                    source,
                    cell,
                    localisation,
                    product_oracle,
                    product_caps,
                    product_box,
                    lower,
                    upper,
                    purity,
                    leaf_max_nodes,
                    bound_safety,
                    minimum_width,
                    cover,
                )
                node["leaf_cover"] = leaf_cover
                if leaf_cover["target_closed"]:
                    node["disposition"] = "closed"
                    node["closure_method"] = "singleton-spatial-state-choi-ppt"
                else:
                    node["disposition"] = "leaf-unresolved"
        solved += 1
        update_summary(payload, len(open_by_index))
        write_payload(output, payload)
        active_cover = node.get("leaf_cover", cover)
        print(
            json.dumps(
                {
                    "identifier": identifier,
                    "pattern": code,
                    "depth": node["depth"],
                    "source_open_cells": node["source_open_cell_count"],
                    "cartesian_cells": node["cartesian_child_cell_count"],
                    "disposition": node["disposition"],
                    "root_bound": (
                        active_cover.get("cover_upper_bound")
                        if active_cover is not None
                        else localisation["base_result"].get("bound")
                    ),
                    "pending_cluster_nodes": len(payload["pending"]),
                }
            ),
            flush=True,
        )
    del (
        localisation_oracle,
        localisation_caps,
        product_oracle,
        product_caps,
        lower,
        upper,
        purity,
    )
    gc.collect()
    return solved


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--frontier-json", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-cluster-nodes", type=int, default=10)
    parser.add_argument("--leaf-max-nodes", type=int, default=100)
    parser.add_argument("--coordinate-safety", type=float, default=2e-6)
    parser.add_argument("--bound-safety", type=float, default=2e-6)
    parser.add_argument("--minimum-width", type=float, default=1e-6)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    source_raw = args.frontier_json.read_bytes()
    source = json.loads(source_raw)
    if not source.get("statuses_complete"):
        raise ValueError("source spectral frontier has unresolved statuses")
    open_cells = source_open_cells(source)
    open_by_index = {int(cell["source_index"]): cell for cell in open_cells}
    if args.resume:
        payload = json.loads(args.output.read_bytes())
        if payload.get("schema") != SCHEMA:
            raise ValueError("resume payload has the wrong schema")
        if payload["source"]["sha256"] != hashlib.sha256(source_raw).hexdigest():
            raise ValueError("resume source hash does not match")
    else:
        payload = initial_payload(
            args.frontier_json,
            source_raw,
            source,
            open_cells,
            args.coordinate_safety,
            args.bound_safety,
            args.minimum_width,
            args.leaf_max_nodes,
        )
        write_payload(args.output, payload)

    remaining = int(args.max_cluster_nodes)
    patterns = sorted(
        {
            tuple(payload["nodes"][str(identifier)]["pattern"])
            for identifier in payload["pending"]
        },
        key=lambda pattern: -max(
            float(payload["nodes"][str(identifier)]["maximum_source_bound"])
            for identifier in payload["pending"]
            if tuple(payload["nodes"][str(identifier)]["pattern"]) == pattern
        ),
    )
    for pattern in patterns:
        if remaining <= 0:
            break
        solved = process_pattern(
            source,
            payload,
            pattern,
            open_by_index,
            remaining,
            args.coordinate_safety,
            args.bound_safety,
            args.minimum_width,
            args.leaf_max_nodes,
            args.output,
        )
        remaining -= solved
    update_summary(payload, len(open_cells))
    write_payload(args.output, payload)
    print(json.dumps(payload["summary"], indent=2))


if __name__ == "__main__":
    main()

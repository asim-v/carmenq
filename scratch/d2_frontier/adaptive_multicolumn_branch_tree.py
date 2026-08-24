"""Adaptive flagged-contraction tree for one open behavior-cover cell.

Each expansion selects a violated real coefficient vector, then partitions the
input trace norm into two scalar branches and a proved cube-face cap cover of
the Bloch branch.  Therefore every expansion is an exhaustive disjunction,
even though the separating coefficient itself is chosen numerically.
"""

from __future__ import annotations

import argparse
import heapq
import json
import math
from pathlib import Path
import numpy as np

from flagged_contraction_separator import family_from_payload, find_worst_contraction
from fourier_behavior_cap_cover import cube_face_caps, plane_caps
from fourier_behavior_upper import solve_behavior_outer
from fourier_branch_upper import PRIOR_BOX
from multicolumn_contraction_cell_cover import run_cover


def fixed_base_data(
    plane_cells: int,
    plane_index: int,
    face_grid: int,
    sphere_index: int,
    pair_cap_index: int,
    pair_branch: str = "bloch",
) -> tuple[tuple[object | None, ...], dict[str, object]]:
    plane = plane_caps(plane_cells)[plane_index]
    sphere = cube_face_caps(face_grid)[sphere_index]
    pair_axes = [
        (normal, 1.0 / np.sqrt(3.0))
        for axis in range(3)
        for sign in (-1.0, 1.0)
        for normal in [np.eye(3)[axis] * sign]
    ]
    caps = (None, (*plane[0], plane[1]), (*sphere[0], sphere[1]))
    base_cut = {
        "pair": (2, 3),
        "scale": 1.0,
        "branch": pair_branch,
    }
    if pair_branch == "bloch":
        pair = pair_axes[pair_cap_index]
        base_cut["cap"] = (*pair[0], pair[1])
    elif pair_branch != "scalar-positive":
        raise ValueError("the base pair branch must be scalar-positive or bloch")
    return caps, base_cut


def solve_node(
    caps: tuple[object | None, ...],
    base_cut: dict[str, object],
    contractions: tuple[dict[str, object], ...],
) -> dict[str, object]:
    return solve_behavior_outer(
        ("bloch", "bloch", "bloch"),
        PRIOR_BOX,
        caps,
        pairwise_contractions=(base_cut, *contractions),
    )


def child_contraction(
    coefficients: list[float],
    row: dict[str, object],
    contraction_grid: int,
) -> dict[str, object]:
    result: dict[str, object] = {
        "coefficients": list(coefficients),
        "branch": str(row["branch"]),
    }
    if row["branch"] == "bloch":
        cap = cube_face_caps(contraction_grid)[int(row["cap"])]
        result["cap"] = [*cap[0], cap[1]]
    return result


def run_tree(
    initial_coefficients: np.ndarray | None,
    target: float,
    max_expansions: int,
    separator_samples: int,
    separator_starts: int,
    seed: int,
    plane_cells: int,
    plane_index: int,
    face_grid: int,
    sphere_index: int,
    pair_cap_index: int,
    contraction_grid: int,
    checkpoint: Path | None,
    resume: bool,
    pair_branch: str = "bloch",
) -> dict[str, object]:
    settings = {
        "target": target,
        "separator_samples": separator_samples,
        "separator_starts": separator_starts,
        "seed": seed,
        "plane_cells": plane_cells,
        "plane_index": plane_index,
        "face_grid": face_grid,
        "sphere_index": sphere_index,
        "pair_cap_index": pair_cap_index,
        "pair_branch": pair_branch,
        "contraction_grid": contraction_grid,
        "initial_coefficients": (
            None
            if initial_coefficients is None
            else (
                np.asarray(initial_coefficients, dtype=float)
                / np.linalg.norm(initial_coefficients)
            ).tolist()
        ),
    }
    if resume:
        if checkpoint is None or not checkpoint.exists():
            raise ValueError("resume requires an existing checkpoint")
        prior = json.loads(checkpoint.read_text(encoding="utf-8"))
        for key, expected in settings.items():
            if prior.get(key) != expected:
                raise ValueError(f"checkpoint mismatch for {key}")
        records = list(prior["expanded_nodes"])
        pending = list(prior["pending_nodes"])
        next_identifier = int(prior["next_identifier"])
    else:
        records = []
        pending = [
            {
                "id": 0,
                "parent": None,
                "depth": 0,
                "bound": math.inf,
                "contractions": [],
                "separator_coefficients": settings["initial_coefficients"],
            }
        ]
        next_identifier = 1

    caps, base_cut = fixed_base_data(
        plane_cells,
        plane_index,
        face_grid,
        sphere_index,
        pair_cap_index,
        pair_branch,
    )
    queue = [(-float(item["bound"]), int(item["id"]), item) for item in pending]
    heapq.heapify(queue)

    def snapshot() -> dict[str, object]:
        leaves = [item[-1] for item in sorted(queue)]
        return {
            **settings,
            "max_expansions": max_expansions,
            "solver_conditional": True,
            "expanded_nodes": records,
            "pending_nodes": leaves,
            "next_identifier": next_identifier,
            "expansion_count": len(records),
            "closed_leaf_count": sum(int(item["closed_children"]) for item in records),
            "open_leaf_count": len(leaves),
            "maximum_open_bound": max(
                (float(item["bound"]) for item in leaves), default=None
            ),
            "certificate_complete": not leaves,
        }

    while queue and len(records) < max_expansions:
        node = heapq.heappop(queue)[-1]
        contractions = tuple(node["contractions"])
        if node.get("separator_coefficients") is None:
            family = solve_node(caps, base_cut, contractions)
            prefix, conditioned = family_from_payload(family)
            separation = find_worst_contraction(
                prefix,
                conditioned,
                samples=separator_samples,
                starts=separator_starts,
                seed=seed + int(node["id"]),
            )
            coefficients = np.asarray(separation["coefficients"], dtype=float)
            source_bound = float(family["bound"])
            violation = float(separation["violation"])
        else:
            coefficients = np.asarray(node["separator_coefficients"], dtype=float)
            source_bound = float(node["bound"])
            if not math.isfinite(source_bound):
                source_bound = float(
                    solve_node(caps, base_cut, contractions)["bound"]
                )
            violation = None

        cover = run_cover(
            coefficients,
            plane_cells,
            plane_index,
            face_grid,
            sphere_index,
            pair_cap_index,
            contraction_grid,
            contractions,
            target=target,
            pair_branch=pair_branch,
        )
        open_rows = [
            row for row in cover["branches"] if float(row["bound"]) >= target
        ]
        child_identifiers: dict[tuple[str, int | None], int] = {}
        children: list[dict[str, object]] = []
        for row in open_rows:
            child = {
                "id": next_identifier,
                "parent": int(node["id"]),
                "depth": int(node["depth"]) + 1,
                "bound": float(row["bound"]),
                "contractions": [
                    *contractions,
                    child_contraction(coefficients.tolist(), row, contraction_grid),
                ],
                "separator_coefficients": None,
            }
            child_identifiers[(str(row["branch"]), row["cap"])] = next_identifier
            children.append(child)
            next_identifier += 1

        audited_branches = []
        for row in cover["branches"]:
            bound = float(row["bound"])
            branch_key = (str(row["branch"]), row["cap"])
            audited_branches.append(
                {
                    "branch": row["branch"],
                    "cap": row["cap"],
                    "bound": bound if math.isfinite(bound) else None,
                    "status": row["status"],
                    "closed": bound < target,
                    "child_id": child_identifiers.get(branch_key),
                }
            )
        record = {
            "id": int(node["id"]),
            "parent": node["parent"],
            "depth": int(node["depth"]),
            "source_bound": source_bound,
            "separator_coefficients": coefficients.tolist(),
            "separator_violation": violation,
            "maximum_child_bound": float(cover["maximum_bound"]),
            "closed_children": len(cover["branches"]) - len(open_rows),
            "open_children": len(open_rows),
            "branches": audited_branches,
        }
        records.append(record)
        print(json.dumps({key: value for key, value in record.items() if key != "branches"}), flush=True)
        for child in children:
            heapq.heappush(queue, (-float(child["bound"]), int(child["id"]), child))
        if checkpoint is not None:
            checkpoint.parent.mkdir(parents=True, exist_ok=True)
            checkpoint.write_text(
                json.dumps(snapshot(), indent=2) + "\n", encoding="utf-8"
            )

    result = snapshot()
    if checkpoint is not None:
        checkpoint.parent.mkdir(parents=True, exist_ok=True)
        checkpoint.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("initial_separator", type=Path, nargs="?")
    parser.add_argument("--target", type=float, default=0.758)
    parser.add_argument("--max-expansions", type=int, default=20)
    parser.add_argument("--separator-samples", type=int, default=20_000)
    parser.add_argument("--separator-starts", type=int, default=12)
    parser.add_argument("--seed", type=int, default=260830)
    parser.add_argument("--plane-cells", type=int, default=8)
    parser.add_argument("--plane-index", type=int, default=4)
    parser.add_argument("--face-grid", type=int, default=4)
    parser.add_argument("--sphere-index", type=int, default=18)
    parser.add_argument("--pair-cap-index", type=int, default=3)
    parser.add_argument(
        "--pair-branch", choices=("scalar-positive", "bloch"), default="bloch"
    )
    parser.add_argument("--contraction-grid", type=int, default=4)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    initial_coefficients = None
    if args.initial_separator is not None:
        separator = json.loads(args.initial_separator.read_text(encoding="utf-8"))
        initial_coefficients = np.asarray(separator["coefficients"], dtype=float)
    result = run_tree(
        initial_coefficients,
        args.target,
        args.max_expansions,
        args.separator_samples,
        args.separator_starts,
        args.seed,
        args.plane_cells,
        args.plane_index,
        args.face_grid,
        args.sphere_index,
        args.pair_cap_index,
        args.contraction_grid,
        args.output,
        args.resume,
        args.pair_branch,
    )
    print(
        json.dumps(
            {
                "certificate_complete": result["certificate_complete"],
                "expansion_count": result["expansion_count"],
                "open_leaf_count": result["open_leaf_count"],
                "maximum_open_bound": result["maximum_open_bound"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

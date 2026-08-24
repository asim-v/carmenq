"""Adaptive measured-contraction tree inside one ternary Fourier cap cell.

Every node accumulates necessary trace-norm contractions for a common quantum
instrument.  At the highest relaxed node, a deterministic numerical search
selects a real coefficient vector whose measured flagged output most exceeds
the qubit input trace norm.  The coefficient is then imposed on an exhaustive
scalar-positive/scalar-negative/Bloch-active angular cover.

The separator search need not be globally certified: every nonzero coefficient
defines a valid contraction, and the child spectral cover is exhaustive.  The
SOCP bounds remain solver-conditional, as do the parent terminal and angular
covers.  This script handles one fixed terminal-parameter box and one fixed
Fourier angular cell; a forest wrapper is required for global closure.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
import math
from pathlib import Path
from typing import Any

import cvxpy as cp
import numpy as np
from scipy.optimize import differential_evolution

from fourier_behavior_cap_cover import cube_face_caps, plane_caps
from fourier_behavior_upper import CHARACTERS
from pairwise_inellipse_box_cover import Box, serialise_box
from ternary_common_instrument_cover import BRANCH_NAME
from ternary_multicolumn_cell_cover import load_ranked_leaf
from ternary_probability_cone_cover import TernaryConeOracle
from terminal_reconstruction_enclosure import reconstruction_anchor_and_errors


@dataclass
class PendingNode:
    identifier: int
    parent: int | None
    depth: int
    contractions: tuple[dict[str, object], ...]
    result: dict[str, Any]

    @property
    def bound(self) -> float:
        return float(self.result["bound"])


def fixed_fourier_contractions(
    code: str,
    plane: tuple[np.ndarray, float] | None,
    sphere: tuple[np.ndarray, float] | None,
) -> tuple[dict[str, object], ...]:
    if len(code) != 3 or any(letter not in BRANCH_NAME for letter in code):
        raise ValueError("invalid Fourier branch code")
    result: list[dict[str, object]] = []
    bloch_seen = 0
    for coefficients, letter in zip(CHARACTERS, code, strict=True):
        contraction: dict[str, object] = {
            "coefficients": np.asarray(coefficients, dtype=float),
            "branch": BRANCH_NAME[letter],
        }
        if letter == "b":
            contraction["gauge_rank"] = bloch_seen
            if bloch_seen == 1:
                if plane is None:
                    raise ValueError("the second Bloch Fourier branch needs a plane cap")
                contraction["cap"] = np.append(plane[0], plane[1])
            elif bloch_seen == 2:
                if sphere is None:
                    raise ValueError("the third Bloch Fourier branch needs a sphere cap")
                contraction["cap"] = np.append(sphere[0], sphere[1])
            bloch_seen += 1
        result.append(contraction)
    return tuple(result)


def render_contraction(item: dict[str, object]) -> dict[str, object]:
    rendered: dict[str, object] = {
        "coefficients": np.asarray(item["coefficients"], dtype=float).tolist(),
        "branch": str(item["branch"]),
    }
    if "gauge_rank" in item:
        rendered["gauge_rank"] = int(item["gauge_rank"])
    if "cap" in item:
        cap = item["cap"]
        if isinstance(cap, cp.Parameter):
            rendered["cap"] = np.asarray(cap.value, dtype=float).tolist()
            rendered["cap_is_scaled"] = True
        else:
            rendered["cap"] = np.asarray(cap, dtype=float).tolist()
    return rendered


def compact_result(result: dict[str, Any]) -> dict[str, object]:
    return {
        key: result[key]
        for key in (
            "status",
            "raw_value",
            "bound",
            "audit",
            "return",
            "terminal_weight_intervals",
            "prior_intervals",
            "iterations",
        )
        if key in result
    }


def build_oracle(
    support_weight: float,
    prefix_order: tuple[int, int, int, int],
    maximum_weight_floor: float,
    projective_support_upper: float,
    projective_support_lines: tuple[tuple[float, float], ...],
    contractions: tuple[dict[str, object], ...],
    terminal_reconstruction: tuple[np.ndarray, np.ndarray] | None = None,
) -> TernaryConeOracle:
    return TernaryConeOracle(
        support_weight,
        prefix_order,
        (),
        (),
        maximum_weight_floor,
        projective_support_upper,
        projective_support_lines=projective_support_lines,
        common_contractions=contractions,
        terminal_reconstruction=terminal_reconstruction,
    )


def measured_gap(
    coefficients: np.ndarray,
    statistics: np.ndarray,
    priors: np.ndarray,
    input_vectors: np.ndarray,
) -> float:
    value = np.asarray(coefficients, dtype=float)
    norm = float(np.linalg.norm(value))
    if norm <= 1e-14:
        return -math.inf
    value = value / norm
    measured = float(
        np.abs(np.tensordot(value, statistics, axes=(0, 0))).sum()
    )
    scalar = abs(float(value @ priors))
    bloch = float(np.linalg.norm(value @ input_vectors))
    return measured - max(scalar, bloch)


def find_separator(
    result: dict[str, Any],
    seed: int,
    terminal_reconstruction: tuple[np.ndarray, np.ndarray] | None = None,
) -> dict[str, object]:
    statistics = np.asarray(result["statistics"], dtype=float)
    priors = np.asarray(result["prefix"], dtype=float)
    input_vectors = np.asarray(result["input_bloch_vectors"], dtype=float)

    def contraction_gap(raw: np.ndarray) -> float:
        if terminal_reconstruction is None:
            return measured_gap(raw, statistics, priors, input_vectors)
        norm = float(np.linalg.norm(raw))
        if norm <= 1e-14:
            return -math.inf
        coefficients = raw / norm
        anchor, errors = terminal_reconstruction
        flagged = 0.0
        for y in range(4):
            signed = np.tensordot(coefficients, statistics[:, y, :], axes=(0, 0))
            trace = abs(float(signed.sum()))
            planar = float(np.linalg.norm(anchor @ signed))
            error = float(
                np.sum(
                    np.abs(coefficients)[:, None]
                    * statistics[:, y, :]
                    * errors[None, :]
                )
            )
            flagged += max(trace, planar - error, 0.0)
        scalar = abs(float(coefficients @ priors))
        bloch = float(np.linalg.norm(coefficients @ input_vectors))
        return flagged - max(scalar, bloch)

    def objective(raw: np.ndarray) -> float:
        return -contraction_gap(raw)

    optimum = differential_evolution(
        objective,
        [(-1.0, 1.0)] * 4,
        seed=seed,
        popsize=20,
        maxiter=800,
        tol=1e-9,
        polish=True,
        workers=1,
    )
    coefficients = optimum.x / np.linalg.norm(optimum.x)
    measured = float(
        np.abs(np.tensordot(coefficients, statistics, axes=(0, 0))).sum()
    )
    scalar = float(coefficients @ priors)
    bloch = float(np.linalg.norm(coefficients @ input_vectors))
    return {
        "coefficients": coefficients,
        "gap": contraction_gap(coefficients),
        "measured_l1": measured,
        "input_scalar": scalar,
        "input_bloch_norm": bloch,
        "optimizer_success": bool(optimum.success),
        "optimizer_iterations": int(optimum.nit),
    }


def run_tree(
    box: Box,
    base_code: str,
    base_plane: tuple[np.ndarray, float] | None,
    base_sphere: tuple[np.ndarray, float] | None,
    support_weight: float = 0.55,
    maximum_weight_floor: float = 0.79,
    projective_support_upper: float = 0.7573,
    projective_support_lines: tuple[tuple[float, float], ...] = ((0.6, 0.76591),),
    prefix_order: tuple[int, int, int, int] = (0, 1, 2, 3),
    contraction_grid: int = 4,
    target: float = 0.758,
    safety: float = 2e-6,
    max_expansions: int = 3,
    seed: int = 20260824,
    checkpoint: Path | None = None,
    use_reconstruction: bool = False,
) -> dict[str, Any]:
    if contraction_grid < 2:
        raise ValueError("contraction_grid must be at least two")
    base = fixed_fourier_contractions(base_code, base_plane, base_sphere)
    terminal_reconstruction = None
    reconstruction_audit = None
    if use_reconstruction:
        anchor, errors, reconstruction_audit = reconstruction_anchor_and_errors(
            box["terminal_alpha"], box["terminal_beta"]
        )
        terminal_reconstruction = (anchor, errors)
    root_oracle = build_oracle(
        support_weight,
        prefix_order,
        maximum_weight_floor,
        projective_support_upper,
        projective_support_lines,
        base,
        terminal_reconstruction,
    )
    root_result = root_oracle.solve(box, safety, capture=True)
    pending = [PendingNode(0, None, 0, (), root_result)]
    records: list[dict[str, object]] = []
    next_identifier = 1
    solved_nodes = 1
    caps = cube_face_caps(contraction_grid)

    def snapshot() -> dict[str, Any]:
        return {
            "support_weight": support_weight,
            "maximum_weight_floor": maximum_weight_floor,
            "projective_support_upper": projective_support_upper,
            "projective_support_lines": [list(line) for line in projective_support_lines],
            "prefix_order": list(prefix_order),
            "box": serialise_box(box),
            "base_code": base_code,
            "base_plane": None if base_plane is None else {
                "normal": base_plane[0].tolist(), "cosine": base_plane[1]
            },
            "base_sphere": None if base_sphere is None else {
                "normal": base_sphere[0].tolist(), "cosine": base_sphere[1]
            },
            "contraction_grid": contraction_grid,
            "contraction_cap_count": len(caps),
            "contraction_covering_cosine": caps[0][1],
            "terminal_reconstruction": reconstruction_audit,
            "target": target,
            "safety": safety,
            "max_expansions": max_expansions,
            "expanded_nodes": len(records),
            "solved_nodes": solved_nodes,
            "open_nodes": len(pending),
            "maximum_open_bound": max(
                (node.bound for node in pending), default=-math.inf
            ),
            "complete": not pending,
            "nodes": records,
            "pending": [
                {
                    "id": node.identifier,
                    "parent": node.parent,
                    "depth": node.depth,
                    "contractions": [
                        render_contraction(item) for item in node.contractions
                    ],
                    **compact_result(node.result),
                }
                for node in sorted(pending, key=lambda item: item.bound, reverse=True)
            ],
            "scope": (
                "one terminal box and one Fourier angular cell; adaptive "
                "measured common-instrument contractions; solver-conditional"
            ),
        }

    def save() -> None:
        if checkpoint is not None:
            checkpoint.parent.mkdir(parents=True, exist_ok=True)
            checkpoint.write_text(
                json.dumps(snapshot(), indent=2) + "\n", encoding="utf-8"
            )

    while pending and len(records) < max_expansions:
        node = max(pending, key=lambda item: item.bound)
        pending.remove(node)
        if node.bound < target:
            records.append(
                {
                    "id": node.identifier,
                    "parent": node.parent,
                    "depth": node.depth,
                    "status": "closed_by_upper_bound",
                    **compact_result(node.result),
                }
            )
            save()
            continue
        separator = find_separator(
            node.result,
            seed + node.identifier,
            terminal_reconstruction,
        )
        coefficients = np.asarray(separator["coefficients"], dtype=float)
        inherited = base + node.contractions
        bloch_rank = sum(
            str(item["branch"]) == "bloch" for item in inherited
        )
        children: list[PendingNode] = []
        child_summaries: list[dict[str, object]] = []
        for branch in ("scalar-positive", "scalar-negative"):
            new = {
                "coefficients": coefficients,
                "branch": branch,
            }
            oracle = build_oracle(
                support_weight,
                prefix_order,
                maximum_weight_floor,
                projective_support_upper,
                projective_support_lines,
                inherited + (new,),
                terminal_reconstruction,
            )
            result = oracle.solve(box, safety, capture=True)
            child = PendingNode(
                next_identifier,
                node.identifier,
                node.depth + 1,
                node.contractions + (new,),
                result,
            )
            next_identifier += 1
            solved_nodes += 1
            children.append(child)
            child_summaries.append(
                {
                    "id": child.identifier,
                    "branch": branch,
                    "cap": None,
                    **compact_result(result),
                }
            )
        cap_parameter = cp.Parameter(3)
        new_bloch: dict[str, object] = {
            "coefficients": coefficients,
            "branch": "bloch",
            "gauge_rank": bloch_rank,
            "cap": cap_parameter,
        }
        oracle = build_oracle(
            support_weight,
            prefix_order,
            maximum_weight_floor,
            projective_support_upper,
            projective_support_lines,
            inherited + (new_bloch,),
            terminal_reconstruction,
        )
        for cap_index, (normal, cosine) in enumerate(caps):
            cap_parameter.value = normal / cosine
            result = oracle.solve(box, safety, capture=True)
            stored_bloch = {
                "coefficients": coefficients,
                "branch": "bloch",
                "gauge_rank": bloch_rank,
                "cap": np.append(normal, cosine),
            }
            child = PendingNode(
                next_identifier,
                node.identifier,
                node.depth + 1,
                node.contractions + (stored_bloch,),
                result,
            )
            next_identifier += 1
            solved_nodes += 1
            children.append(child)
            child_summaries.append(
                {
                    "id": child.identifier,
                    "branch": "bloch",
                    "cap": cap_index,
                    **compact_result(result),
                }
            )
        closed = [child for child in children if child.bound < target]
        pending.extend(child for child in children if child.bound >= target)
        records.append(
            {
                "id": node.identifier,
                "parent": node.parent,
                "depth": node.depth,
                "status": "expanded",
                **compact_result(node.result),
                "separator": {
                    key: (
                        value.tolist() if isinstance(value, np.ndarray) else value
                    )
                    for key, value in separator.items()
                },
                "branches_per_expansion": len(children),
                "closed_children": len(closed),
                "open_children": len(children) - len(closed),
                "maximum_child_bound": max(child.bound for child in children),
                "children": child_summaries,
            }
        )
        save()
        print(
            json.dumps(
                {
                    "expanded": len(records),
                    "separator_gap": separator["gap"],
                    "closed_children": len(closed),
                    "open_nodes": len(pending),
                    "maximum_open_bound": max(
                        (item.bound for item in pending), default=-math.inf
                    ),
                }
            ),
            flush=True,
        )
    return snapshot()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--box-json", type=Path, required=True)
    parser.add_argument("--leaf-rank", type=int, default=0)
    parser.add_argument("--base-code", default="pbb")
    parser.add_argument("--base-plane-cells", type=int, default=16)
    parser.add_argument("--base-plane-index", type=int)
    parser.add_argument("--base-face-grid", type=int, default=2)
    parser.add_argument("--base-sphere-index", type=int)
    parser.add_argument("--contraction-grid", type=int, default=4)
    parser.add_argument("--target", type=float, default=0.758)
    parser.add_argument("--max-expansions", type=int, default=3)
    parser.add_argument("--seed", type=int, default=20260824)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--reconstruction", action="store_true")
    args = parser.parse_args()
    planes = plane_caps(args.base_plane_cells)
    spheres = cube_face_caps(args.base_face_grid)
    base_plane = (
        None if args.base_plane_index is None else planes[args.base_plane_index]
    )
    base_sphere = (
        None if args.base_sphere_index is None else spheres[args.base_sphere_index]
    )
    payload = run_tree(
        load_ranked_leaf(args.box_json, args.leaf_rank),
        args.base_code,
        base_plane,
        base_sphere,
        contraction_grid=args.contraction_grid,
        target=args.target,
        max_expansions=args.max_expansions,
        seed=args.seed,
        checkpoint=args.output,
        use_reconstruction=args.reconstruction,
    )
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "complete": payload["complete"],
                "expanded_nodes": payload["expanded_nodes"],
                "solved_nodes": payload["solved_nodes"],
                "open_nodes": payload["open_nodes"],
                "maximum_open_bound": payload["maximum_open_bound"],
            }
        )
    )


if __name__ == "__main__":
    main()

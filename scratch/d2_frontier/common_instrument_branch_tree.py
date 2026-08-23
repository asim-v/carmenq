"""Adaptive common-instrument branch-and-cut for one fixed terminal POVM.

This is the first spatial wrapper around the facially reduced Choi moment
upper model.  A node is an axis-aligned cell in the four prefix-state Bloch
coordinates.  At a high relaxed maximizer the script:

1. projects the reported conditioned outputs onto the exact fixed-input
   common-instrument set;
2. extracts a separating witness;
3. partitions the node into a central cell on which the robust witness is
   strong and a disjoint family of side cells; and
4. propagates every witness to every descendant using a cellwise trace-radius
   correction.

The tree is a solver-conditional upper certificate for one fixed terminal POVM
and prefix-prior order.  It is not a cover over terminal POVM geometry.
"""

from __future__ import annotations

import argparse
from collections import deque
from dataclasses import dataclass
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from audit_common_instrument_candidate import load_reported_family
from carmenq.common_instrument import project_to_common_instrument
from choi_moment_reduced_upper import PAULIS, solve_povm
from common_instrument_cells import (
    StateCell,
    cover_volume_residual,
    partition_one_trace_ball_coordinate,
)
from two_block_choi_seesaw import canonical_three_effect_povm


@dataclass(frozen=True)
class WitnessTemplate:
    reference_bloch: np.ndarray
    witness_bloch: np.ndarray
    lipschitz: np.ndarray
    reference_support: float
    source_gap: float

    def __post_init__(self) -> None:
        reference = np.asarray(self.reference_bloch, dtype=float)
        witness = np.asarray(self.witness_bloch, dtype=float)
        lipschitz = np.asarray(self.lipschitz, dtype=float)
        if reference.shape != (4, 4) or witness.shape != (4, 4, 4):
            raise ValueError("invalid witness-template Bloch dimensions")
        if lipschitz.shape != (4,):
            raise ValueError("a witness template requires four Lipschitz constants")
        object.__setattr__(self, "reference_bloch", reference.copy())
        object.__setattr__(self, "witness_bloch", witness.copy())
        object.__setattr__(self, "lipschitz", lipschitz.copy())

    def cut_for_cell(self, cell: StateCell) -> dict[str, object]:
        return {
            "reference_bloch": self.reference_bloch,
            "witness_bloch": self.witness_bloch,
            "lipschitz": self.lipschitz,
            "reference_support": self.reference_support,
            "input_trace_radii": cell.maximum_trace_radii(
                self.reference_bloch
            ),
            "restrict_to_balls": False,
        }

    @property
    def uniform_radius_budget(self) -> float:
        denominator = float(self.lipschitz.sum())
        return math.inf if denominator <= 1e-15 else self.source_gap / denominator

    def to_json(self) -> dict[str, object]:
        return {
            "reference_bloch": self.reference_bloch.tolist(),
            "witness_bloch": self.witness_bloch.tolist(),
            "lipschitz": self.lipschitz.tolist(),
            "reference_support": self.reference_support,
            "source_gap": self.source_gap,
            "uniform_radius_budget": self.uniform_radius_budget,
        }

    @staticmethod
    def from_json(payload: dict[str, Any]) -> "WitnessTemplate":
        return WitnessTemplate(
            reference_bloch=np.asarray(payload["reference_bloch"], dtype=float),
            witness_bloch=np.asarray(payload["witness_bloch"], dtype=float),
            lipschitz=np.asarray(payload["lipschitz"], dtype=float),
            reference_support=float(payload["reference_support"]),
            source_gap=float(payload["source_gap"]),
        )


@dataclass(frozen=True)
class PendingNode:
    identifier: int
    parent: int | None
    depth: int
    cell: StateCell
    witness_indices: tuple[int, ...]
    inherited_upper: float


def witness_template_from_projection(
    prefix_bloch: np.ndarray,
    witness: np.ndarray,
    lipschitz: np.ndarray,
    reference_support: float,
    gap: float,
) -> WitnessTemplate:
    witness_bloch = np.asarray(
        [
            [
                [float(np.trace(witness[z, y] @ pauli).real) for pauli in PAULIS]
                for y in range(4)
            ]
            for z in range(4)
        ]
    )
    return WitnessTemplate(
        reference_bloch=np.asarray(prefix_bloch, dtype=float),
        witness_bloch=witness_bloch,
        lipschitz=np.asarray(lipschitz, dtype=float),
        reference_support=float(reference_support),
        source_gap=float(gap),
    )


def solve_cell(
    terminal: np.ndarray,
    support_weight: float,
    prefix_order: tuple[int, int, int, int],
    cell: StateCell,
    templates: list[WitnessTemplate],
    witness_indices: tuple[int, ...],
    data_processing_scales: tuple[float, ...],
    data_processing_mode: str,
    solver: str,
) -> dict[str, object]:
    cuts = tuple(templates[index].cut_for_cell(cell) for index in witness_indices)
    return solve_povm(
        terminal,
        support_weight,
        prefix_order,
        "used",
        data_processing_mode,
        data_processing_scales,
        cell.bounds,
        cuts,
        solver,
        False,
    )


def run_tree(
    support_weight: float,
    terminal_weights: np.ndarray,
    prefix_order: tuple[int, int, int, int],
    target: float,
    safety: float,
    max_nodes: int,
    radius_fraction: float,
    data_processing_scales: tuple[float, ...],
    data_processing_mode: str,
    solver: str,
    checkpoint: Path | None,
    resume: bool,
) -> dict[str, object]:
    if max_nodes < 1:
        raise ValueError("max_nodes must be positive")
    if not 0.0 < radius_fraction < 1.0:
        raise ValueError("radius_fraction must lie strictly between zero and one")
    terminal = canonical_three_effect_povm(terminal_weights)
    if resume:
        if checkpoint is None or not checkpoint.exists():
            raise ValueError("resume requires an existing output checkpoint")
        previous = json.loads(checkpoint.read_text(encoding="utf-8"))
        expected = {
            "support_weight": support_weight,
            "terminal_effect_weights": np.append(terminal_weights, 0.0).tolist(),
            "prefix_order": list(prefix_order),
            "target": target,
            "safety": safety,
            "radius_fraction": radius_fraction,
            "data_processing_scales": list(data_processing_scales),
            "data_processing_mode": data_processing_mode,
            "solver": solver,
        }
        for key, value_expected in expected.items():
            if previous.get(key) != value_expected:
                raise ValueError(f"resume checkpoint mismatch for {key}")
        templates = [
            WitnessTemplate.from_json(item)
            for item in previous["witness_templates"]
        ]
        records = list(previous["nodes"])
        record_by_id = {int(item["id"]): item for item in records}
        pending_items = []
        for item in previous["pending"]:
            parent = None if item["parent"] is None else int(item["parent"])
            inherited_upper = float(item["inherited_upper"])
            ancestor = parent
            while ancestor is not None:
                ancestor_record = record_by_id[ancestor]
                if ancestor_record.get("bound") is not None:
                    inherited_upper = min(
                        inherited_upper, float(ancestor_record["bound"])
                    )
                ancestor = (
                    None
                    if ancestor_record["parent"] is None
                    else int(ancestor_record["parent"])
                )
            pending_items.append(
                PendingNode(
                    int(item["id"]),
                    parent,
                    int(item["depth"]),
                    StateCell.from_json(item["cell"]),
                    tuple(map(int, item["witness_indices"])),
                    inherited_upper,
                )
            )
        pending = deque(pending_items)
        identifiers = [int(item["id"]) for item in records] + [
            item.identifier for item in pending
        ]
        next_identifier = max(identifiers, default=-1) + 1
        closed_nodes = int(previous["closed_nodes"])
    else:
        root = StateCell.root(prefix_order)
        pending = deque([PendingNode(0, None, 0, root, (), math.inf)])
        next_identifier = 1
        templates = []
        records = []
        closed_nodes = 0

    def snapshot() -> dict[str, object]:
        finite_open = [item.inherited_upper for item in pending if math.isfinite(item.inherited_upper)]
        maximum_open = math.inf if any(
            not math.isfinite(item.inherited_upper) for item in pending
        ) else (max(finite_open) if finite_open else -math.inf)
        return {
            "support_weight": support_weight,
            "terminal_effect_weights": np.append(terminal_weights, 0.0).tolist(),
            "prefix_order": list(prefix_order),
            "target": target,
            "safety": safety,
            "max_nodes": max_nodes,
            "radius_fraction": radius_fraction,
            "data_processing_scales": list(data_processing_scales),
            "data_processing_mode": data_processing_mode,
            "solver": solver,
            "solved_nodes": len(records),
            "closed_nodes": closed_nodes,
            "open_nodes": len(pending),
            "maximum_open_inherited_bound": maximum_open,
            "complete": not pending,
            "witness_templates": [item.to_json() for item in templates],
            "nodes": records,
            "pending": [
                {
                    "id": item.identifier,
                    "parent": item.parent,
                    "depth": item.depth,
                    "cell": item.cell.to_json(),
                    "witness_indices": list(item.witness_indices),
                    "inherited_upper": item.inherited_upper,
                }
                for item in pending
            ],
            "scope": (
                "fixed ternary terminal POVM and fixed prefix-prior order; "
                "solver-conditional state-cell branch-and-cut"
            ),
        }

    def save() -> None:
        if checkpoint is not None:
            checkpoint.parent.mkdir(parents=True, exist_ok=True)
            checkpoint.write_text(
                json.dumps(snapshot(), indent=2) + "\n", encoding="utf-8"
            )

    while pending and len(records) < max_nodes:
        node = max(
            pending,
            key=lambda item: (
                item.inherited_upper,
                item.cell.volume,
                -item.depth,
                -item.identifier,
            ),
        )
        pending.remove(node)
        active_witness_indices = tuple(range(len(templates)))
        try:
            result = solve_cell(
                terminal,
                support_weight,
                prefix_order,
                node.cell,
                templates,
                active_witness_indices,
                data_processing_scales,
                data_processing_mode,
                solver,
            )
        except RuntimeError as error:
            if "infeasible" in str(error).lower():
                records.append(
                    {
                        "id": node.identifier,
                        "parent": node.parent,
                        "depth": node.depth,
                        "status": "infeasible",
                        "cell": node.cell.to_json(),
                        "witness_indices": list(active_witness_indices),
                    }
                )
                closed_nodes += 1
                save()
                continue
            raise

        raw_bound = float(result["bound"])
        certified_bound = raw_bound + safety
        record: dict[str, Any] = {
            "id": node.identifier,
            "parent": node.parent,
            "depth": node.depth,
            "cell": node.cell.to_json(),
            "cell_volume": node.cell.volume,
            "witness_indices": list(active_witness_indices),
            "status": "solved",
            "raw_bound": raw_bound,
            "bound": certified_bound,
            "audit": float(result["audit"]),
            "return": float(result["return"]),
            "prefix_bloch_coefficients": result["prefix_bloch_coefficients"],
            "moment_rank": int(result["moment_numerical_rank_1e-7"]),
            "solver_status": result["status"],
            "solver_stats": result["solver_stats"],
            "applied_witness_cuts": result["common_instrument_witness_cuts"],
        }
        if certified_bound < target:
            record["status"] = "closed_by_upper_bound"
            records.append(record)
            closed_nodes += 1
            save()
            continue

        states, outputs = load_reported_family(result)
        projection = project_to_common_instrument(states, outputs, solver=solver)
        record["choi_audit"] = {
            "distance": projection.distance,
            "separation_gap": projection.separation_gap,
            "reference_support": projection.compatible_support_value,
            "input_lipschitz_constants": projection.input_lipschitz_constants.tolist(),
            "uniform_radius_budget": projection.uniform_input_radius_budget,
            "trace_preservation_residual": projection.trace_preservation_residual,
            "minimum_choi_eigenvalue": projection.minimum_choi_eigenvalue,
        }
        if projection.separation_gap <= 1e-7:
            left, right = node.cell.split_widest()
            children = (left, right)
            record["status"] = "split_unseparated"
            record["cover_volume_residual"] = (
                left.volume + right.volume - node.cell.volume
            )
            record["children"] = [next_identifier, next_identifier + 1]
            for child in reversed(children):
                pending.appendleft(
                    PendingNode(
                        next_identifier,
                        node.identifier,
                        node.depth + 1,
                        child,
                        active_witness_indices,
                        min(node.inherited_upper, certified_bound),
                    )
                )
                next_identifier += 1
            records.append(record)
            save()
            continue

        template = witness_template_from_projection(
            np.asarray(result["prefix_bloch_coefficients"], dtype=float),
            projection.witness,
            projection.input_lipschitz_constants,
            projection.compatible_support_value,
            projection.separation_gap,
        )
        template_index = len(templates)
        templates.append(template)
        central_radius = radius_fraction * template.uniform_radius_budget
        radii = np.full(4, central_radius)
        sides, center, split_coordinate = partition_one_trace_ball_coordinate(
            node.cell, template.reference_bloch, radii
        )
        witness_indices = tuple(range(len(templates)))
        children: list[StateCell] = list(sides)
        if center is not None and center.volume > 0.0:
            # Solve the witness-covered central child first.
            children.insert(0, center)
        if not children:
            record["status"] = "empty_partition"
            closed_nodes += 1
        else:
            record["status"] = "split_by_choi_witness"
            record["new_witness_index"] = template_index
            record["central_radius"] = central_radius
            record["split_coordinate"] = (
                None if split_coordinate is None else list(split_coordinate)
            )
            record["side_child_count"] = len(sides)
            record["cover_volume_residual"] = cover_volume_residual(
                node.cell, sides, center
            )
            identifiers = list(range(next_identifier, next_identifier + len(children)))
            record["children"] = identifiers
            queued = [
                PendingNode(
                    identifier,
                    node.identifier,
                    node.depth + 1,
                    child,
                    witness_indices,
                    min(node.inherited_upper, certified_bound),
                )
                for identifier, child in zip(identifiers, children, strict=True)
            ]
            for child_node in reversed(queued):
                pending.appendleft(child_node)
            next_identifier += len(children)
        records.append(record)
        save()

    return snapshot()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lambda", dest="support_weight", type=float, default=0.55)
    parser.add_argument(
        "--terminal-weights", type=float, nargs=3, default=(0.92, 0.64, 0.44)
    )
    parser.add_argument("--prefix-order", type=int, nargs=4, default=(0, 1, 2, 3))
    parser.add_argument("--target", type=float, default=0.758)
    parser.add_argument("--safety", type=float, default=2e-6)
    parser.add_argument("--max-nodes", type=int, default=3)
    parser.add_argument("--radius-fraction", type=float, default=0.9)
    parser.add_argument(
        "--data-processing-scale", type=float, action="append", default=None
    )
    parser.add_argument(
        "--data-processing",
        choices=("quadratic", "cell"),
        default="cell",
    )
    parser.add_argument("--solver", choices=("clarabel", "scs"), default="clarabel")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    weights = np.asarray(args.terminal_weights, dtype=float)
    if weights.shape != (3,) or np.any(weights <= 0.0) or abs(float(weights.sum()) - 2.0) > 1e-9:
        raise ValueError("three terminal effect weights must be positive and sum to two")
    order = tuple(args.prefix_order)
    if sorted(order) != [0, 1, 2, 3]:
        raise ValueError("prefix order must be a permutation")
    scales = tuple((1.0,) if args.data_processing_scale is None else args.data_processing_scale)
    payload = run_tree(
        args.support_weight,
        weights,
        order,
        args.target,
        args.safety,
        args.max_nodes,
        args.radius_fraction,
        scales,
        args.data_processing,
        args.solver,
        args.output,
        args.resume,
    )
    rendered = json.dumps(payload, indent=2) + "\n"
    print(rendered, end="")
    if args.output is not None:
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

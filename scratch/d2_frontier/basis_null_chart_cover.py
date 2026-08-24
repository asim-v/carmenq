"""Adaptive exact cover of the singular common-instrument stratum.

Let ``R`` be the 4 by 4 matrix whose rows are the Pauli coordinates of the
four subnormalised input states.  If ``det(R)=0``, a nonzero vector ``c``
satisfies ``c.T @ R = 0``.  Choose an index ``k`` on which ``abs(c[k])`` is
maximal and normalise ``c[k]=1``.  The other three coefficients then belong
to ``[-1,1]``.  The four choices of ``k`` therefore cover the singular
stratum exactly.

This script partitions each coefficient cube into disjoint axis-aligned
boxes.  On every box it asks SCIP whether the exact common-instrument model
can attain a target score.  A node reported infeasible is closed; an
unresolved node is bisected.  The resulting JSON contains the complete leaf
cover, including unresolved leaves, so that numerical solver statements are
never silently promoted to analytic proofs.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import heapq
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from joint_effect_helstrom_scip import build, canonical_three_effect_povm


DEFAULT_PRIOR_BOUNDS = np.asarray(
    [
        [0.296875, 0.42596435546875],
        [0.224609375, 0.34832000732421875],
        [0.15234375, 0.258392333984375],
        [0.1083984375, 0.201324462890625],
    ],
    dtype=float,
)


@dataclass(frozen=True)
class Cell:
    """One coefficient box in one normalised left-null-vector chart."""

    identifier: int
    pivot: int
    depth: int
    bounds: np.ndarray
    parent: int | None = None
    inherited_upper: float = math.inf

    def __post_init__(self) -> None:
        bounds = np.asarray(self.bounds, dtype=float)
        if bounds.shape != (3, 2):
            raise ValueError("cell bounds must have shape (3,2)")
        if np.any(bounds[:, 0] > bounds[:, 1]):
            raise ValueError("cell lower bounds exceed upper bounds")
        object.__setattr__(self, "bounds", bounds.copy())

    @property
    def volume(self) -> float:
        return float(np.prod(self.bounds[:, 1] - self.bounds[:, 0]))

    def to_json(self) -> dict[str, Any]:
        return {
            "id": self.identifier,
            "parent": self.parent,
            "pivot": self.pivot,
            "depth": self.depth,
            "bounds": self.bounds.tolist(),
            "volume": self.volume,
            "inherited_upper": self.inherited_upper,
        }

    @staticmethod
    def from_json(payload: dict[str, Any]) -> "Cell":
        return Cell(
            identifier=int(payload["id"]),
            parent=None if payload["parent"] is None else int(payload["parent"]),
            pivot=int(payload["pivot"]),
            depth=int(payload["depth"]),
            bounds=np.asarray(payload["bounds"], dtype=float),
            inherited_upper=float(payload["inherited_upper"]),
        )


def bisect_cell(cell: Cell, first_id: int) -> tuple[Cell, Cell]:
    """Return a disjoint binary partition along the widest coefficient."""

    widths = cell.bounds[:, 1] - cell.bounds[:, 0]
    axis = int(np.argmax(widths))
    midpoint = float(cell.bounds[axis].mean())
    lower = cell.bounds.copy()
    upper = cell.bounds.copy()
    lower[axis, 1] = midpoint
    upper[axis, 0] = midpoint
    common = {
        "pivot": cell.pivot,
        "depth": cell.depth + 1,
        "parent": cell.identifier,
        "inherited_upper": cell.inherited_upper,
    }
    return (
        Cell(identifier=first_id, bounds=lower, **common),
        Cell(identifier=first_id + 1, bounds=upper, **common),
    )


def solve_cell(
    cell: Cell,
    effects: np.ndarray,
    weight: float,
    target: float,
    prior_bounds: np.ndarray,
    seconds: float,
    gap: float,
    feasibility_tolerance: float,
) -> dict[str, Any]:
    """Solve one exact target-feasibility subproblem."""

    model, variables = build(
        effects=effects,
        weight=weight,
        prefix_order=(0, 1, 2, 3),
        target=target,
        fix_rotation_gauge=True,
        linked_columns=None,
        require_cp_completion=True,
        prefix_prior_bounds=prior_bounds,
        basis_null_pivot=cell.pivot,
        basis_null_bounds=cell.bounds,
    )
    model.setRealParam("limits/time", seconds)
    model.setRealParam("limits/gap", gap)
    model.setRealParam("numerics/feastol", feasibility_tolerance)
    model.setRealParam("numerics/dualfeastol", feasibility_tolerance)
    model.setIntParam("display/verblevel", 0)
    model.optimize()

    status = str(model.getStatus())
    primal = float(model.getPrimalbound())
    dual = float(model.getDualbound())
    result: dict[str, Any] = {
        **cell.to_json(),
        "status": status,
        "primal_bound": primal,
        "dual_bound": dual,
        "nodes": int(model.getNNodes()),
        "solving_time": float(model.getSolvingTime()),
        "closed": status == "infeasible",
    }
    solution = model.getBestSol()
    if solution is not None:
        result["incumbent_score"] = float(
            model.getSolVal(solution, variables["score"])
        )
        result["incumbent_audit"] = float(
            model.getSolVal(solution, variables["audit"])
        )
        result["incumbent_return"] = float(
            model.getSolVal(solution, variables["return"])
        )
    return result


def run_cover(
    weight: float,
    terminal_weights: np.ndarray,
    target: float,
    prior_bounds: np.ndarray,
    seconds_per_node: float,
    gap: float,
    feasibility_tolerance: float,
    max_nodes: int,
    max_depth: int,
    pivots: tuple[int, ...],
    checkpoint: Path | None,
    resume: bool,
) -> dict[str, Any]:
    """Run or resume the adaptive four-chart cover."""

    if max_nodes < 1:
        raise ValueError("max_nodes must be positive")
    if max_depth < 0:
        raise ValueError("max_depth must be nonnegative")
    if any(pivot not in range(4) for pivot in pivots) or len(set(pivots)) != len(pivots):
        raise ValueError("pivots must be distinct members of 0,1,2,3")

    settings = {
        "weight": weight,
        "terminal_effect_weights": terminal_weights.tolist(),
        "target": target,
        "prefix_order": [0, 1, 2, 3],
        "prefix_prior_bounds": prior_bounds.tolist(),
        "seconds_per_node": seconds_per_node,
        "gap": gap,
        "feasibility_tolerance": feasibility_tolerance,
        "max_depth": max_depth,
        "pivots": list(pivots),
    }
    if resume:
        if checkpoint is None or not checkpoint.exists():
            raise ValueError("resume requires an existing checkpoint")
        previous = json.loads(checkpoint.read_text(encoding="utf-8"))
        for key, expected in settings.items():
            if previous.get(key) != expected:
                raise ValueError(f"resume checkpoint mismatch for {key}")
        records = list(previous["nodes"])
        pending = [Cell.from_json(item) for item in previous["pending"]]
        identifiers = [int(item["id"]) for item in records] + [
            item.identifier for item in pending
        ]
        next_identifier = max(identifiers, default=-1) + 1
        target_feasible = bool(previous.get("target_feasible", False))
    else:
        records = []
        pending = [
            Cell(
                identifier=index,
                pivot=pivot,
                depth=0,
                bounds=np.asarray([[-1.0, 1.0]] * 3),
            )
            for index, pivot in enumerate(pivots)
        ]
        next_identifier = len(pending)
        target_feasible = False

    effects = canonical_three_effect_povm(terminal_weights)

    def priority(cell: Cell) -> tuple[float, int, int, Cell]:
        # Resolve the largest inherited upper bound first, then the shallowest
        # cell.  The identifier makes the ordering deterministic.
        return (-cell.inherited_upper, cell.depth, cell.identifier, cell)

    queue = [priority(cell) for cell in pending]
    heapq.heapify(queue)

    def snapshot() -> dict[str, Any]:
        leaves = [item[-1] for item in queue]
        by_pivot: dict[str, dict[str, Any]] = {}
        for pivot in pivots:
            open_pivot = [cell for cell in leaves if cell.pivot == pivot]
            solved_pivot = [item for item in records if int(item["pivot"]) == pivot]
            by_pivot[str(pivot)] = {
                "solved_nodes": len(solved_pivot),
                "closed_nodes": sum(bool(item["closed"]) for item in solved_pivot),
                "open_leaves": len(open_pivot),
                "open_volume": sum(cell.volume for cell in open_pivot),
                "maximum_inherited_upper": (
                    max((cell.inherited_upper for cell in open_pivot), default=None)
                ),
            }
        return {
            **settings,
            "max_nodes": max_nodes,
            "solver_conditional": True,
            "certificate_complete": not queue and not target_feasible,
            "target_feasible": target_feasible,
            "nodes": records,
            "pending": [item[-1].to_json() for item in sorted(queue)],
            "summary_by_pivot": by_pivot,
        }

    while queue and len(records) < max_nodes and not target_feasible:
        cell = heapq.heappop(queue)[-1]
        result = solve_cell(
            cell,
            effects,
            weight,
            target,
            prior_bounds,
            seconds_per_node,
            gap,
            feasibility_tolerance,
        )
        records.append(result)
        print(
            f"node={cell.identifier} pivot={cell.pivot} depth={cell.depth} "
            f"status={result['status']} dual={result['dual_bound']:.12g}",
            flush=True,
        )
        if not result["closed"]:
            if result.get("incumbent_score", -math.inf) >= target - feasibility_tolerance:
                # A feasible point at the target is a counterexample to the
                # proposed upper bound.  Preserve it and stop subdividing.
                result["target_feasible"] = True
                target_feasible = True
                queue = []
            elif cell.depth >= max_depth:
                unresolved = Cell(
                    identifier=cell.identifier,
                    pivot=cell.pivot,
                    depth=cell.depth,
                    bounds=cell.bounds,
                    parent=cell.parent,
                    inherited_upper=min(cell.inherited_upper, result["dual_bound"]),
                )
                heapq.heappush(queue, priority(unresolved))
                break
            else:
                inherited = min(cell.inherited_upper, result["dual_bound"])
                parent_with_bound = Cell(
                    identifier=cell.identifier,
                    pivot=cell.pivot,
                    depth=cell.depth,
                    bounds=cell.bounds,
                    parent=cell.parent,
                    inherited_upper=inherited,
                )
                children = bisect_cell(parent_with_bound, next_identifier)
                next_identifier += 2
                for child in children:
                    heapq.heappush(queue, priority(child))
        if checkpoint is not None:
            checkpoint.parent.mkdir(parents=True, exist_ok=True)
            checkpoint.write_text(json.dumps(snapshot(), indent=2) + "\n", encoding="utf-8")

    result = snapshot()
    if checkpoint is not None:
        checkpoint.parent.mkdir(parents=True, exist_ok=True)
        checkpoint.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lambda", dest="weight", type=float, default=0.55)
    parser.add_argument(
        "--fixed-three-povm-weights",
        type=float,
        nargs=3,
        default=(0.92, 0.64, 0.44),
    )
    parser.add_argument("--target", type=float, default=0.758)
    parser.add_argument(
        "--prefix-prior-bounds", type=float, nargs=8, default=DEFAULT_PRIOR_BOUNDS.ravel()
    )
    parser.add_argument("--seconds-per-node", type=float, default=10.0)
    parser.add_argument("--gap", type=float, default=1e-6)
    parser.add_argument("--feasibility-tolerance", type=float, default=1e-9)
    parser.add_argument("--max-nodes", type=int, default=256)
    parser.add_argument("--max-depth", type=int, default=12)
    parser.add_argument("--pivot", type=int, action="append", choices=(0, 1, 2, 3))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    result = run_cover(
        weight=args.weight,
        terminal_weights=np.asarray(args.fixed_three_povm_weights, dtype=float),
        target=args.target,
        prior_bounds=np.asarray(args.prefix_prior_bounds, dtype=float).reshape(4, 2),
        seconds_per_node=args.seconds_per_node,
        gap=args.gap,
        feasibility_tolerance=args.feasibility_tolerance,
        max_nodes=args.max_nodes,
        max_depth=args.max_depth,
        pivots=tuple(args.pivot if args.pivot is not None else range(4)),
        checkpoint=args.output,
        resume=args.resume,
    )
    print(json.dumps({
        "certificate_complete": result["certificate_complete"],
        "solved_nodes": len(result["nodes"]),
        "summary_by_pivot": result["summary_by_pivot"],
    }, indent=2))


if __name__ == "__main__":
    main()

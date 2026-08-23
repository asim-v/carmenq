"""Two-dimensional cover using exact aligned-projective comparisons.

For a ternary rank-one terminal POVM, reciprocal Horwitz parameters
``(alpha, beta)`` determine all three effect traces.  The probability-cone
oracle also carries normalized projector statistics.  Replacing the terminal
readout by ``{Pi_i, I-Pi_i}`` therefore has the exact audit value

    u[i, i] + p[j] - u[j, i].

Every such replacement has the same RETURN score and must obey every known
supporting line of the independently covered projective sector.  This driver
covers only the two terminal parameters; no auxiliary pair-inellipse chart is
needed.  The resulting artifact is numerical until every projective line and
every conic node is validated with outward rounding.
"""

from __future__ import annotations

import argparse
import heapq
import json
import math
from pathlib import Path
from typing import Any

from pairwise_inellipse_box_cover import (
    Box,
    deserialise_box,
    serialise_box,
    split_box,
    write_payload,
)
from ternary_probability_cone_cover import (
    TERMINAL_ALPHA,
    TERMINAL_BETA,
    TernaryConeOracle,
    initial_box,
    terminal_domain_intersects,
)


def compact_result(result: dict[str, Any]) -> dict[str, Any]:
    return {
        key: result[key]
        for key in (
            "status",
            "raw_value",
            "bound",
            "audit",
            "return",
            "terminal_weight_intervals",
            "iterations",
        )
        if key in result
    }


def branch_coordinate(box: Box, root: Box, alpha_branch_weight: float) -> str:
    """A fixed sensitivity-weighted bisection rule for the two coordinates."""

    alpha_relative = (box[TERMINAL_ALPHA][1] - box[TERMINAL_ALPHA][0]) / (
        root[TERMINAL_ALPHA][1] - root[TERMINAL_ALPHA][0]
    )
    beta_relative = (box[TERMINAL_BETA][1] - box[TERMINAL_BETA][0]) / (
        root[TERMINAL_BETA][1] - root[TERMINAL_BETA][0]
    )
    return (
        TERMINAL_ALPHA
        if alpha_branch_weight * alpha_relative >= beta_relative
        else TERMINAL_BETA
    )


def cover(
    support_weight: float,
    maximum_weight_floor: float,
    projective_support_upper: float,
    projective_support_lines: tuple[tuple[float, float], ...],
    prefix_order: tuple[int, int, int, int],
    target: float,
    safety: float,
    max_nodes: int,
    alpha_branch_weight: float,
    resume_payload: dict[str, Any] | None = None,
) -> dict[str, Any]:
    if alpha_branch_weight <= 0.0:
        raise ValueError("alpha_branch_weight must be positive")
    root = initial_box(0, maximum_weight_floor, include_priors=False)
    oracle = TernaryConeOracle(
        support_weight,
        prefix_order,
        (),
        (),
        maximum_weight_floor,
        projective_support_upper,
        projective_support_lines=projective_support_lines,
    )
    if resume_payload is None:
        first = oracle.solve(root, safety)
        solved = 1
        counter = 0
        queue: list[tuple[float, int, Box, dict[str, Any]]] = [
            (-float(first["bound"]), counter, root, first)
        ]
        leaves: list[dict[str, Any]] = []
        resumed_from = None
    else:
        expected = {
            "support_weight": support_weight,
            "maximum_weight_floor": maximum_weight_floor,
            "projective_support_upper": projective_support_upper,
            "safety": safety,
            "alpha_branch_weight": alpha_branch_weight,
        }
        for key, value in expected.items():
            recorded = resume_payload.get(key, 4.0 if key == "alpha_branch_weight" else None)
            if recorded is None or not math.isclose(
                float(recorded), float(value), rel_tol=0.0, abs_tol=1e-15
            ):
                raise ValueError(f"resume payload has incompatible {key}")
        recorded_target = float(resume_payload["target"])
        if target + 1e-15 < recorded_target:
            raise ValueError(
                "resume target may only be relaxed upward; a tighter target "
                "would invalidate previously closed leaves"
            )
        recorded_lines = tuple(
            tuple(map(float, line))
            for line in resume_payload.get("projective_support_lines", [])
        )
        if recorded_lines != projective_support_lines:
            raise ValueError("resume payload has incompatible projective lines")
        if tuple(resume_payload.get("prefix_order", ())) != prefix_order:
            raise ValueError("resume payload has incompatible prefix order")
        solved = int(resume_payload["solved_nodes"])
        resumed_from = solved
        leaves = list(resume_payload.get("leaves", []))
        queue = []
        counter = 0
        for item in resume_payload.get("open_nodes", []):
            result = {
                key: item[key]
                for key in (
                    "status", "raw_value", "bound", "audit", "return",
                    "terminal_weight_intervals", "iterations",
                )
                if key in item
            }
            box = deserialise_box(item["box"])
            heapq.heappush(
                queue, (-float(result["bound"]), counter, box, result)
            )
            counter += 1
        if not queue:
            raise ValueError("resume payload has no open nodes")
    while queue:
        negative_bound, node_id, box, result = heapq.heappop(queue)
        if -negative_bound <= target:
            leaves.append({"box": serialise_box(box), **compact_result(result)})
            continue
        if solved >= max_nodes:
            heapq.heappush(queue, (negative_bound, node_id, box, result))
            break
        coordinate = branch_coordinate(box, root, alpha_branch_weight)
        for child in split_box(box, coordinate):
            if not terminal_domain_intersects(child, maximum_weight_floor):
                leaves.append(
                    {
                        "box": serialise_box(child),
                        "status": "domain_empty",
                        "bound": -math.inf,
                    }
                )
                continue
            child_result = oracle.solve(child, safety)
            solved += 1
            counter += 1
            record = {
                "box": serialise_box(child),
                **compact_result(child_result),
            }
            if float(child_result["bound"]) <= target:
                leaves.append(record)
            else:
                heapq.heappush(
                    queue,
                    (-float(child_result["bound"]), counter, child, child_result),
                )
        if solved % 100 <= 1:
            print(
                json.dumps(
                    {
                        "nodes": solved,
                        "open": len(queue),
                        "maximum_open_bound": (
                            -queue[0][0] if queue else -math.inf
                        ),
                    }
                ),
                flush=True,
            )
    open_nodes = [
        {"box": serialise_box(box), **compact_result(result)}
        for _, _, box, result in sorted(queue)
    ]
    return {
        "support_weight": support_weight,
        "maximum_weight_floor": maximum_weight_floor,
        "projective_support_upper": projective_support_upper,
        "projective_support_lines": [list(line) for line in projective_support_lines],
        "prefix_order": list(prefix_order),
        "target": target,
        "safety": safety,
        "max_nodes": max_nodes,
        "alpha_branch_weight": alpha_branch_weight,
        "resumed_from_solved_nodes": resumed_from,
        "solved_nodes": solved,
        "complete": not open_nodes,
        "maximum_open_bound": max(
            (float(item["bound"]) for item in open_nodes), default=-math.inf
        ),
        "maximum_leaf_bound": max(
            (float(item["bound"]) for item in leaves), default=-math.inf
        ),
        "open_nodes": open_nodes,
        "leaves": leaves,
        "logical_scope": (
            "full sorted ternary terminal-weight strip with maximum effect "
            "trace above the supplied floor"
        ),
        "numerical_status": (
            "finite SOCP outer cover at solver tolerances; every auxiliary "
            "projective line and conic dual still requires outward validation"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lambda", dest="support_weight", type=float, default=0.6)
    parser.add_argument("--maximum-weight-floor", type=float, default=0.88325)
    parser.add_argument("--projective-support-upper", type=float, default=0.76591)
    parser.add_argument(
        "--projective-line",
        type=float,
        nargs=2,
        action="append",
        default=[],
        metavar=("WEIGHT", "UPPER"),
    )
    parser.add_argument(
        "--resume",
        type=Path,
        help="continue an incomplete cover artifact with matching parameters",
    )
    parser.add_argument("--prefix-order", type=int, nargs=4, default=(0, 1, 2, 3))
    parser.add_argument("--target", type=float, default=0.76593)
    parser.add_argument("--safety", type=float, default=2e-6)
    parser.add_argument("--max-nodes", type=int, default=2001)
    parser.add_argument(
        "--alpha-branch-weight",
        type=float,
        default=4.0,
        help="relative bisection priority for the terminal-alpha coordinate",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    resume_payload = (
        None
        if args.resume is None
        else json.loads(args.resume.read_text(encoding="utf-8"))
    )
    payload = cover(
        args.support_weight,
        args.maximum_weight_floor,
        args.projective_support_upper,
        tuple(tuple(map(float, line)) for line in args.projective_line),
        tuple(args.prefix_order),
        args.target,
        args.safety,
        args.max_nodes,
        args.alpha_branch_weight,
        resume_payload,
    )
    write_payload(args.output, payload)
    print(
        json.dumps(
            {
                "complete": payload["complete"],
                "nodes": payload["solved_nodes"],
                "open": len(payload["open_nodes"]),
                "maximum_open_bound": payload["maximum_open_bound"],
            }
        )
    )


if __name__ == "__main__":
    main()

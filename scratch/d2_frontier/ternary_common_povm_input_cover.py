"""Branch-and-bound cover of a high-score input-Pauli region.

The common-instrument relaxation represents the four subnormalised qubit
inputs by the rows of a Pauli-coordinate matrix ``R``.  A common instrument
and terminal readout necessarily induce one twelve-outcome input POVM.  On an
input box ``R in [L,U]``, positivity of every effect ``a`` gives the robust
outer relation

    |q[z,k] - C[z].a[k]| <= sum_mu D[z,mu] a0[k],

where ``C=(L+U)/2`` and ``D=(U-L)/2``.  This script first bounds every input
coordinate among relaxed points whose score can reach the target.  It then
partitions that candidate box and applies the robust common-POVM SOCP on each
node.  Closing every node certifies the selected terminal/Fourier cell,
conditional on the numerical conic solves used by the surrounding frontier.
"""

from __future__ import annotations

import argparse
import heapq
import json
import math
from pathlib import Path
from typing import Any

import cvxpy as cp
import numpy as np

from fourier_behavior_cap_cover import cube_face_caps
from pairwise_inellipse_box_cover import Box, deserialise_box
from ternary_multicolumn_branch_tree import fixed_fourier_contractions
from ternary_probability_cone_cover import TernaryConeOracle
from terminal_reconstruction_enclosure import reconstruction_anchor_and_errors


def _plane(data: dict[str, Any] | None) -> tuple[np.ndarray, float] | None:
    if data is None:
        return None
    return np.asarray(data["normal"], dtype=float), float(data["cosine"])


def _configuration(
    source: dict[str, Any],
    use_top_spectral_cell: bool,
) -> tuple[
    Box,
    tuple[dict[str, object], ...],
    tuple[np.ndarray, np.ndarray],
]:
    box = deserialise_box(source["box"])
    contractions = fixed_fourier_contractions(
        str(source["base_code"]),
        _plane(source.get("base_plane")),
        _plane(source.get("base_sphere")),
    )
    if use_top_spectral_cell:
        cell = source.get("top_cell")
        if not cell:
            raise ValueError("the source frontier has no top spectral cell")
        coefficients = tuple(
            np.asarray(item, dtype=float)
            for item in source["separator_coefficients"]
        )
        grids = tuple(int(value) for value in source["separator_grids"])
        branches = tuple(str(value) for value in cell["branches"])
        cap_indices = tuple(cell["caps"])
        if not (
            len(coefficients) == len(grids) == len(branches) == len(cap_indices)
        ):
            raise ValueError("top spectral-cell arrays are not parallel")
        additions: list[dict[str, object]] = []
        bloch_count = sum(item["branch"] == "bloch" for item in contractions)
        for coefficient, grid, branch, cap_index in zip(
            coefficients, grids, branches, cap_indices, strict=True
        ):
            item: dict[str, object] = {
                "coefficients": coefficient,
                "branch": branch,
            }
            if branch == "bloch":
                if cap_index is None:
                    raise ValueError("a Bloch spectral branch needs a cap")
                normal, cosine = cube_face_caps(grid)[int(cap_index)]
                item.update(
                    {
                        "gauge_rank": bloch_count,
                        "cap": np.append(normal, cosine),
                    }
                )
                bloch_count += 1
            elif cap_index is not None:
                raise ValueError("a scalar spectral branch cannot have a cap")
            additions.append(item)
        contractions = (*contractions, *additions)
    anchor, errors, _ = reconstruction_anchor_and_errors(
        box["terminal_alpha"], box["terminal_beta"]
    )
    return box, contractions, (anchor, errors)


def _oracle_keywords(source: dict[str, Any]) -> dict[str, Any]:
    return {
        "support_weight": float(source["support_weight"]),
        "prefix_order": tuple(int(value) for value in source["prefix_order"]),
        "pairs": (),
        "coordinate_cases": (),
        "maximum_weight_floor": float(source["maximum_weight_floor"]),
        "projective_support_upper": float(source["projective_support_upper"]),
        "projective_support_lines": tuple(
            tuple(float(value) for value in line)
            for line in source["projective_support_lines"]
        ),
    }


def _compact_result(result: dict[str, Any]) -> dict[str, Any]:
    return {
        key: result[key]
        for key in (
            "status",
            "raw_value",
            "bound",
            "audit",
            "return",
            "score",
            "iterations",
        )
        if key in result
    }


def localise_candidate_region(
    source: dict[str, Any],
    coordinate_safety: float,
    use_top_spectral_cell: bool = False,
) -> dict[str, Any]:
    """Enclose all base-relaxation inputs whose score reaches the target."""

    box, contractions, reconstruction = _configuration(
        source, use_top_spectral_cell
    )
    oracle = TernaryConeOracle(
        **_oracle_keywords(source),
        common_contractions=contractions,
        terminal_reconstruction=reconstruction,
        build_input_region_problem=True,
    )
    base = oracle.solve(box, float(source["safety"]), capture=True)
    if base["status"] not in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}:
        raise RuntimeError(f"base localisation solve failed: {base['status']}")
    lower = np.empty((4, 4), dtype=float)
    upper = np.empty((4, 4), dtype=float)
    supports: list[dict[str, Any]] = []
    for z in range(4):
        for mu in range(4):
            direction = np.zeros((4, 4), dtype=float)
            direction[z, mu] = 1.0
            positive = oracle.solve_input_support(
                float(source["target"]), direction, coordinate_safety
            )
            negative = oracle.solve_input_support(
                float(source["target"]), -direction, coordinate_safety
            )
            if positive["status"] not in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}:
                raise RuntimeError(
                    f"positive input support ({z},{mu}) failed: {positive['status']}"
                )
            if negative["status"] not in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}:
                raise RuntimeError(
                    f"negative input support ({z},{mu}) failed: {negative['status']}"
                )
            upper[z, mu] = float(positive["bound"])
            lower[z, mu] = -float(negative["bound"])
            supports.append(
                {
                    "coordinate": [z, mu],
                    "lower": float(lower[z, mu]),
                    "upper": float(upper[z, mu]),
                    "negative": _compact_result(negative),
                    "positive": _compact_result(positive),
                }
            )
            print(
                json.dumps(
                    {
                        "phase": "localise",
                        "coordinate": [z, mu],
                        "lower": lower[z, mu],
                        "upper": upper[z, mu],
                    }
                ),
                flush=True,
            )
    if np.any(lower > upper):
        raise RuntimeError("localised input box is empty")
    return {
        "base_result": _compact_result(base),
        "lower": lower.tolist(),
        "upper": upper.tolist(),
        "row_l1_widths": np.sum(upper - lower, axis=1).tolist(),
        "maximum_coordinate_width": float(np.max(upper - lower)),
        "coordinate_safety": float(coordinate_safety),
        "top_spectral_cell": bool(use_top_spectral_cell),
        "supports": supports,
    }


def _set_common_povm_box(
    lower_parameter: cp.Parameter,
    upper_parameter: cp.Parameter,
    anchor_parameter: cp.Parameter,
    radii_parameter: cp.Parameter,
    trace_radii_parameter: cp.Parameter,
    lower: np.ndarray,
    upper: np.ndarray,
) -> None:
    lower_parameter.value = lower
    upper_parameter.value = upper
    anchor_parameter.value = 0.5 * (lower + upper)
    radii = 0.5 * (upper - lower)
    radii_parameter.value = radii
    trace_radii_parameter.value = np.maximum(
        radii[:, 0], np.linalg.norm(radii[:, 1:], axis=1)
    )


def _split_coordinate(lower: np.ndarray, upper: np.ndarray) -> tuple[int, int]:
    widths = upper - lower
    row = int(np.argmax(np.sum(widths, axis=1)))
    coordinate = int(np.argmax(widths[row]))
    return row, coordinate


def _node_payload(
    identifier: int,
    parent: int | None,
    depth: int,
    lower: np.ndarray,
    upper: np.ndarray,
) -> dict[str, Any]:
    return {
        "identifier": int(identifier),
        "parent": parent,
        "depth": int(depth),
        "lower": lower.tolist(),
        "upper": upper.tolist(),
    }


def _write_checkpoint(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def cover_candidate_region(
    source: dict[str, Any],
    localisation: dict[str, Any],
    output: Path,
    max_nodes: int,
    bound_safety: float,
    minimum_width: float,
    checkpoint_every: int,
    use_top_spectral_cell: bool,
    resume: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Cover a localised input region by robust common-effective-POVM boxes."""

    box, contractions, reconstruction = _configuration(
        source, use_top_spectral_cell
    )
    anchor_parameter = cp.Parameter((4, 4))
    radii_parameter = cp.Parameter((4, 4), nonneg=True)
    trace_radii_parameter = cp.Parameter(4, nonneg=True)
    lower_parameter = cp.Parameter((4, 4))
    upper_parameter = cp.Parameter((4, 4))
    oracle = TernaryConeOracle(
        **_oracle_keywords(source),
        common_contractions=contractions,
        terminal_reconstruction=reconstruction,
        common_povm_input_anchor=anchor_parameter,
        common_povm_input_radii=radii_parameter,
        common_povm_trace_radii=trace_radii_parameter,
        input_pauli_lower=lower_parameter,
        input_pauli_upper=upper_parameter,
    )

    target = float(source["target"])
    if resume is None:
        lower = np.asarray(localisation["lower"], dtype=float)
        upper = np.asarray(localisation["upper"], dtype=float)
        pending: list[tuple[float, int, dict[str, Any]]] = [
            (-math.inf, 0, _node_payload(0, None, 0, lower, upper))
        ]
        records: list[dict[str, Any]] = []
        next_identifier = 1
        top_solution: dict[str, Any] | None = None
        maximum_bound = -math.inf
    else:
        pending = []
        for node in resume["pending"]:
            priority = -float(node.get("parent_bound", math.inf))
            heapq.heappush(pending, (priority, int(node["identifier"]), node))
        records = list(resume["records"])
        next_identifier = int(resume["next_identifier"])
        top_solution = resume.get("top_solution")
        maximum_bound = float(resume.get("maximum_bound", -math.inf))

    unresolved: list[dict[str, Any]] = []
    solved_nodes = len(records)
    while pending and solved_nodes < max_nodes:
        _, _, node = heapq.heappop(pending)
        lower = np.asarray(node["lower"], dtype=float)
        upper = np.asarray(node["upper"], dtype=float)
        _set_common_povm_box(
            lower_parameter,
            upper_parameter,
            anchor_parameter,
            radii_parameter,
            trace_radii_parameter,
            lower,
            upper,
        )
        result = oracle.solve(box, bound_safety, capture=True)
        solved_nodes += 1
        record = {
            **node,
            **_compact_result(result),
            "row_l1_radii": (0.5 * np.sum(upper - lower, axis=1)).tolist(),
        }
        records.append(record)
        bound = float(result["bound"])
        if math.isfinite(bound) and bound > maximum_bound:
            maximum_bound = bound
            top_solution = {
                **record,
                "prefix": result.get("prefix"),
                "input_bloch_vectors": result.get("input_bloch_vectors"),
                "statistics": result.get("statistics"),
                "effective_povm": result.get("effective_povm"),
            }
        if result["status"] in {cp.INFEASIBLE, cp.INFEASIBLE_INACCURATE} or bound < target:
            record["disposition"] = "closed"
        elif result["status"] not in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}:
            record["disposition"] = "solver-unresolved"
            unresolved.append(record)
        else:
            row, coordinate = _split_coordinate(lower, upper)
            width = float(upper[row, coordinate] - lower[row, coordinate])
            if width <= minimum_width:
                record["disposition"] = "resolution-limit"
                unresolved.append(record)
            else:
                midpoint = 0.5 * (lower[row, coordinate] + upper[row, coordinate])
                record["disposition"] = "split"
                record["split_coordinate"] = [row, coordinate]
                record["split_value"] = float(midpoint)
                for side in range(2):
                    child_lower = lower.copy()
                    child_upper = upper.copy()
                    if side == 0:
                        child_upper[row, coordinate] = midpoint
                    else:
                        child_lower[row, coordinate] = midpoint
                    child = _node_payload(
                        next_identifier,
                        int(node["identifier"]),
                        int(node["depth"]) + 1,
                        child_lower,
                        child_upper,
                    )
                    child["parent_bound"] = bound
                    heapq.heappush(pending, (-bound, next_identifier, child))
                    next_identifier += 1
        print(
            json.dumps(
                {
                    "phase": "cover",
                    "solved": solved_nodes,
                    "pending": len(pending),
                    "identifier": node["identifier"],
                    "depth": node["depth"],
                    "bound": bound,
                    "disposition": record["disposition"],
                }
            ),
            flush=True,
        )
        if checkpoint_every > 0 and solved_nodes % checkpoint_every == 0:
            checkpoint = _assemble_payload(
                source,
                localisation,
                records,
                pending,
                unresolved,
                next_identifier,
                maximum_bound,
                top_solution,
                max_nodes,
                bound_safety,
                minimum_width,
                use_top_spectral_cell,
            )
            _write_checkpoint(output, checkpoint)

    payload = _assemble_payload(
        source,
        localisation,
        records,
        pending,
        unresolved,
        next_identifier,
        maximum_bound,
        top_solution,
        max_nodes,
        bound_safety,
        minimum_width,
        use_top_spectral_cell,
    )
    _write_checkpoint(output, payload)
    return payload


def _assemble_payload(
    source: dict[str, Any],
    localisation: dict[str, Any],
    records: list[dict[str, Any]],
    pending_heap: list[tuple[float, int, dict[str, Any]]],
    unresolved: list[dict[str, Any]],
    next_identifier: int,
    maximum_bound: float,
    top_solution: dict[str, Any] | None,
    max_nodes: int,
    bound_safety: float,
    minimum_width: float,
    use_top_spectral_cell: bool,
) -> dict[str, Any]:
    pending = [item[2] for item in sorted(pending_heap)]
    closed = sum(record.get("disposition") == "closed" for record in records)
    split = sum(record.get("disposition") == "split" for record in records)
    statuses_complete = all(
        record.get("status")
        in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE, cp.INFEASIBLE, cp.INFEASIBLE_INACCURATE}
        for record in records
    )
    complete = not pending and not unresolved and statuses_complete
    return {
        "support_weight": source["support_weight"],
        "target": source["target"],
        "source": "ternary common-effective-POVM input-box cover",
        "source_box": source["box"],
        "base_code": source["base_code"],
        "base_plane": source.get("base_plane"),
        "base_sphere": source.get("base_sphere"),
        "top_spectral_cell": bool(use_top_spectral_cell),
        "localisation": localisation,
        "bound_safety": float(bound_safety),
        "minimum_width": float(minimum_width),
        "max_nodes": int(max_nodes),
        "solved_nodes": len(records),
        "closed_nodes": int(closed),
        "split_nodes": int(split),
        "pending_nodes": len(pending),
        "unresolved_nodes": len(unresolved),
        "maximum_bound": float(maximum_bound),
        "statuses_complete": bool(statuses_complete),
        "complete": bool(complete),
        "next_identifier": int(next_identifier),
        "pending": pending,
        "unresolved": unresolved,
        "records": records,
        "top_solution": top_solution,
        "scope": (
            "continuous cover of every high-score input Pauli matrix in the "
            "selected terminal/Fourier spectral cell, using one robust common "
            "twelve-outcome effective POVM; numerical SOCP bounds remain "
            "solver-conditional"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--frontier-json", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--localisation-json", type=Path)
    parser.add_argument("--localise-only", action="store_true")
    parser.add_argument("--max-nodes", type=int, default=1000)
    parser.add_argument("--coordinate-safety", type=float, default=2e-6)
    parser.add_argument("--bound-safety", type=float, default=2e-6)
    parser.add_argument("--minimum-width", type=float, default=1e-6)
    parser.add_argument("--checkpoint-every", type=int, default=10)
    parser.add_argument("--top-spectral-cell", action="store_true")
    args = parser.parse_args()
    source = json.loads(args.frontier_json.read_text(encoding="utf-8"))
    resumed: dict[str, Any] | None = None
    if args.resume:
        resumed = json.loads(args.output.read_text(encoding="utf-8"))
        localisation = resumed["localisation"]
    elif args.localisation_json is not None:
        localised_payload = json.loads(
            args.localisation_json.read_text(encoding="utf-8")
        )
        localisation = localised_payload.get("localisation", localised_payload)
    else:
        localisation = localise_candidate_region(
            source, args.coordinate_safety, args.top_spectral_cell
        )
        if args.localise_only:
            payload = {
                "support_weight": source["support_weight"],
                "target": source["target"],
                "source_box": source["box"],
                "base_code": source["base_code"],
                "base_plane": source.get("base_plane"),
                "base_sphere": source.get("base_sphere"),
                "localisation": localisation,
                "scope": "high-score input-coordinate localisation only",
            }
            _write_checkpoint(args.output, payload)
            print(json.dumps({"localised": True, **localisation}, indent=2))
            return
    payload = cover_candidate_region(
        source,
        localisation,
        args.output,
        args.max_nodes,
        args.bound_safety,
        args.minimum_width,
        args.checkpoint_every,
        args.top_spectral_cell,
        resumed,
    )
    print(
        json.dumps(
            {
                key: payload[key]
                for key in (
                    "solved_nodes",
                    "closed_nodes",
                    "split_nodes",
                    "pending_nodes",
                    "unresolved_nodes",
                    "maximum_bound",
                    "statuses_complete",
                    "complete",
                )
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

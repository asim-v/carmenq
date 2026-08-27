"""Replay a rank/rank projective cover with tangent intervals.

The archived SCIP hierarchy is used only as a finite geometric partition.
Every selected cell is re-bounded from its explicit ``x``, ``y``, and angle
ranges by :mod:`projective_tangent_interval_certificate`; archived primal and
dual objective values are not used in the new inequalities.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import json
import math
from pathlib import Path
from typing import Any

from projective_tangent_interval_certificate import (
    decimal_box,
    root_from_bounds,
)
from projective_low_eigenvalue_face import LowEigenvalueFaceCertifier as RankRankCertifier
import validate_projective_cover as archived


ROOT = Path(__file__).resolve().parent
SCHEMA = "carmenq.projective-tangent-interval-cover.rank-rank.v1"
ANGLE_EDGES = (0.0, 0.1, 0.3, 0.5, 1.0 / math.sqrt(2.0))
REFINED_ANGLE_EDGES = (0.0, 0.03, 0.06, 0.1)
DIRECTORY = {
    "coarse-trace": "projective_trace_boxes",
    "fine-trace": "projective_trace_fine",
    "fine-exact": "projective_full_fine",
    "angle-trace": "projective_trace_angle",
    "angle-exact": "projective_full_angle",
    "refined-trace": "projective_trace_angle2",
    "refined-exact": "projective_full_angle2",
    "refined-long": "projective_full_angle2_long",
    "refined-extended": "projective_full_angle2_extended",
}


def explicit_bounds(row: dict[str, Any]) -> dict[str, tuple[float, float]]:
    result = {
        "x": (0.0, 1.0),
        "y": (0.0, 1.0),
        "first": (0.0, 1.0 / math.sqrt(2.0)),
        "second": (0.0, 1.0 / math.sqrt(2.0)),
    }
    raw = row["bounds"]
    if isinstance(raw, dict):
        for name in ("x", "y"):
            if name in raw:
                result[name] = tuple(map(float, raw[name]))
        return result
    for item in raw:
        if item["kind"] == "squared":
            result[str(item["name"])] = (
                float(item["lower"]),
                float(item["upper"]),
            )
        elif item["name"] == "first_angle_0":
            result["first"] = (float(item["lower"]), float(item["upper"]))
        elif item["name"] == "second_angle_0":
            result["second"] = (float(item["lower"]), float(item["upper"]))
    return result


def indexed_angles(name: str, edges: tuple[float, ...]) -> dict[str, tuple[float, float]]:
    suffix = name.rsplit("_a", 1)[1]
    if len(suffix) != 2 or not suffix.isdigit():
        raise ValueError(f"could not parse angle suffix from {name!r}")
    first, second = map(int, suffix)
    return {
        "first": (edges[first], edges[first + 1]),
        "second": (edges[second], edges[second + 1]),
    }


def geometric_leaves() -> list[dict[str, Any]]:
    leaves = []
    for leaf in archived.rank_rank_leaves():
        directory = DIRECTORY[str(leaf["stage"])]
        source = ROOT / directory / f"{leaf['name']}.json"
        row = json.loads(source.read_text(encoding="utf-8"))
        bounds = explicit_bounds(row)
        if leaf["stage"] == "angle-trace":
            bounds.update(indexed_angles(str(leaf["name"]), ANGLE_EDGES))
        elif leaf["stage"] == "refined-trace":
            bounds.update(indexed_angles(str(leaf["name"]), REFINED_ANGLE_EDGES))
        leaves.append(
            {
                "stage": leaf["stage"],
                "name": leaf["name"],
                "certificate_name": f"{leaf['stage']}__{leaf['name']}",
                "source": str(source.relative_to(ROOT)),
                "archived_dual_scaled": leaf["dual_scaled"],
                "bounds": bounds,
            }
        )
    return leaves


def expand_full_angles(leaf: dict[str, Any]) -> list[dict[str, Any]]:
    bounds = leaf["bounds"]
    full = (0.0, 1.0 / math.sqrt(2.0))
    if bounds["first"] != full or bounds["second"] != full:
        return [leaf]
    children = []
    for first in range(len(ANGLE_EDGES) - 1):
        for second in range(len(ANGLE_EDGES) - 1):
            child = {**leaf, "bounds": dict(bounds)}
            child["name"] = f"{leaf['name']}_tangent_a{first}{second}"
            child["certificate_name"] = (
                f"{leaf['certificate_name']}__tangent_a{first}{second}"
            )
            child["parent_name"] = leaf["name"]
            child["bounds"]["first"] = (
                ANGLE_EDGES[first], ANGLE_EDGES[first + 1]
            )
            child["bounds"]["second"] = (
                ANGLE_EDGES[second], ANGLE_EDGES[second + 1]
            )
            children.append(child)
    return children


def select_named_cells(
    cells: list[dict[str, Any]], names: list[str]
) -> list[dict[str, Any]]:
    """Select exact certificate identifiers and reject typos."""

    requested = set(names)
    selected = [
        cell for cell in cells if str(cell["certificate_name"]) in requested
    ]
    found = {str(cell["certificate_name"]) for cell in selected}
    missing = sorted(requested - found)
    if missing:
        raise ValueError(f"unknown certificate names: {missing}")
    return selected


def certify_cell(
    leaf: dict[str, Any],
    level: Fraction,
    max_boxes: int,
    *,
    weight: Fraction = Fraction(3, 5),
    resume_certificate: dict[str, object] | None = None,
) -> dict[str, Any]:
    bounds = leaf["bounds"]
    x_bounds = decimal_box(*map(str, bounds["x"]))
    y_bounds = decimal_box(*map(str, bounds["y"]))
    first = decimal_box(*map(str, bounds["first"]))
    second = decimal_box(*map(str, bounds["second"]))
    certifier = RankRankCertifier(weight, level, x_bounds)
    result = certifier.certify(
        root_from_bounds(x_bounds, y_bounds, first, second),
        max_boxes,
        resume=resume_certificate,
    )
    return {
        **{key: value for key, value in leaf.items() if key != "bounds"},
        "bounds": {key: list(value) for key, value in bounds.items()},
        "certificate": result,
    }


def cover_payload(
    level: Fraction,
    total_cells: int,
    results: list[dict[str, Any]],
    *,
    weight: Fraction = Fraction(3, 5),
) -> dict[str, Any]:
    complete = [item for item in results if item["certificate"]["complete"]]
    return {
        "schema": SCHEMA,
        "weight": str(weight),
        "level": str(level),
        "source_geometric_leaf_count": len(geometric_leaves()),
        "expanded_cell_count": total_cells,
        "processed_cell_count": len(results),
        "complete_cell_count": len(complete),
        "run_complete": len(results) == total_cells,
        "all_cells_complete": (
            len(results) == total_cells and len(complete) == total_cells
        ),
        "total_boxes_split": sum(
            int(item["certificate"]["boxes_split"]) for item in results
        ),
        "cells": results,
        "scope": f"rank/rank projective topology at lambda={weight}",
        "archived_solver_role": (
            "geometric partition only; no archived primal or dual value enters "
            "a tangent-interval cell certificate"
        ),
    }


def write_payload(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--weight", default="0.6")
    parser.add_argument("--level", default="0.76662")
    parser.add_argument("--max-boxes-per-cell", type=int, default=3000)
    parser.add_argument("--minimum-archived-dual", type=float, default=-math.inf)
    parser.add_argument("--limit", type=int)
    parser.add_argument(
        "--cell",
        action="append",
        default=[],
        help="exact certificate_name to run; repeat for multiple cells",
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    weight = Fraction(args.weight)
    level = Fraction(args.level)
    if not Fraction(0) < weight < level < Fraction(1):
        parser.error("require 0 < weight < level < 1")
    cells = [
        child
        for leaf in geometric_leaves()
        if float(leaf["archived_dual_scaled"]) >= args.minimum_archived_dual
        for child in expand_full_angles(leaf)
    ]
    cells.sort(key=lambda item: float(item["archived_dual_scaled"]), reverse=True)
    if args.cell:
        try:
            cells = select_named_cells(cells, args.cell)
        except ValueError as exc:
            parser.error(str(exc))
    if args.limit is not None:
        cells = cells[: args.limit]

    results: list[dict[str, Any]] = []
    prior_by_name: dict[str, dict[str, Any]] = {}
    legacy_by_name: dict[str, dict[str, Any]] = {}
    name_counts = Counter(str(cell["name"]) for cell in cells)
    if args.resume:
        if args.output is None:
            parser.error("--resume requires --output")
        if args.output.exists():
            prior = json.loads(args.output.read_text(encoding="utf-8"))
            if prior.get("weight") != str(weight) or prior.get("level") != str(level):
                parser.error("resume file has incompatible weight or level")
            prior_by_name = {
                str(item["certificate_name"]): item
                for item in prior.get("cells", [])
                if "certificate_name" in item
            }
            legacy_by_name = {
                str(item["name"]): item
                for item in prior.get("cells", [])
                if "certificate_name" not in item
                and name_counts[str(item["name"])] == 1
            }

    for index, cell in enumerate(cells):
        previous = prior_by_name.get(str(cell["certificate_name"]))
        if previous is None:
            previous = legacy_by_name.get(str(cell["name"]))
        if previous is not None and previous["certificate"]["complete"]:
            result = {**previous, "certificate_name": cell["certificate_name"]}
            reused = True
        else:
            result = certify_cell(
                cell,
                level,
                args.max_boxes_per_cell,
                weight=weight,
                resume_certificate=(
                    previous["certificate"] if previous is not None else None
                ),
            )
            reused = False
        results.append(result)
        certificate = result["certificate"]
        print(
            json.dumps(
                {
                    "cell": index + 1,
                    "count": len(cells),
                    "name": result["certificate_name"],
                    "complete": certificate["complete"],
                    "split": certificate["boxes_split"],
                    "remaining": certificate["boxes_remaining"],
                    "maximum_open_upper": certificate["maximum_open_upper"],
                    "reused": reused,
                }
            ),
            flush=True,
        )
        if args.output is not None:
            write_payload(
                args.output,
                cover_payload(level, len(cells), results, weight=weight),
            )

    payload = cover_payload(level, len(cells), results, weight=weight)
    if args.output is not None:
        write_payload(args.output, payload)
    print(
        json.dumps({key: value for key, value in payload.items() if key != "cells"}, indent=2)
    )


if __name__ == "__main__":
    main()

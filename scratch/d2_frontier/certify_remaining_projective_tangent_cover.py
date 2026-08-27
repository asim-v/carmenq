"""Certify the three non-rank/rank projective topologies.

The mixed topologies are exchanged by complementing the coarse effect:
``E -> I-E`` swaps endpoint/rank with rank/endpoint.  It is therefore enough
to certify both ordered mixed topologies on the symmetric chart ``x+y <= 1``.
Endpoint/endpoint has the same complement symmetry.  The ``(y,residual)``
coordinates used by the interval kernel impose that chart exactly through
``x = 1-y-residual``.

Archived SCIP files supply only the finite geometric x/y partition.  No
archived objective value enters a cell inequality.
"""

from __future__ import annotations

import argparse
from functools import cache
from fractions import Fraction
import json
from pathlib import Path
from typing import Any

from certify_rank_rank_tangent_cover import (
    ANGLE_EDGES,
    explicit_bounds,
)
from projective_tangent_interval_certificate import (
    RankRankCertifier,
    decimal_box,
    root_from_bounds,
)
import validate_projective_cover as archived


ROOT = Path(__file__).resolve().parent
SCHEMA = "carmenq.projective-tangent-interval-cover.remaining.v1"
DIRECTORY = {
    "coarse-trace": "projective_er_trace_boxes",
    "fine-trace": "projective_er_trace_fine",
    "fine-exact": "projective_er_full_fine",
}
TOPOLOGIES = (
    ("endpoint", "endpoint"),
    ("endpoint", "rank"),
    ("rank", "endpoint"),
)


@cache
def geometric_leaves() -> list[dict[str, Any]]:
    """Return the archived endpoint/rank hierarchy as a geometric partition."""

    leaves = []
    for leaf in archived.endpoint_rank_leaves():
        source = ROOT / DIRECTORY[str(leaf["stage"])] / f"{leaf['name']}.json"
        row = json.loads(source.read_text(encoding="utf-8"))
        bounds = explicit_bounds(row)
        leaves.append(
            {
                "stage": leaf["stage"],
                "name": leaf["name"],
                "source": str(source.relative_to(ROOT)),
                "bounds": bounds,
            }
        )
    return leaves


def topology_cells() -> list[dict[str, Any]]:
    cells: list[dict[str, Any]] = []
    for first_kind, second_kind in TOPOLOGIES:
        rank_position = (
            "first" if first_kind == "rank"
            else "second" if second_kind == "rank"
            else None
        )
        for leaf in geometric_leaves():
            if rank_position is None:
                child = {**leaf, "bounds": dict(leaf["bounds"])}
                child["first_kind"] = first_kind
                child["second_kind"] = second_kind
                child["certificate_name"] = (
                    f"{first_kind}_{second_kind}__{leaf['name']}"
                )
                child["bounds"]["first"] = (0.0, 0.0)
                child["bounds"]["second"] = (0.0, 0.0)
                cells.append(child)
                continue
            for angle_index in range(len(ANGLE_EDGES) - 1):
                child = {**leaf, "bounds": dict(leaf["bounds"])}
                child["first_kind"] = first_kind
                child["second_kind"] = second_kind
                child["certificate_name"] = (
                    f"{first_kind}_{second_kind}__{leaf['name']}__a{angle_index}"
                )
                child["bounds"]["first"] = (0.0, 0.0)
                child["bounds"]["second"] = (0.0, 0.0)
                child["bounds"][rank_position] = (
                    ANGLE_EDGES[angle_index],
                    ANGLE_EDGES[angle_index + 1],
                )
                cells.append(child)
    return cells


def certify_cell(
    cell: dict[str, Any], level: Fraction, max_boxes: int,
    *, weight: Fraction = Fraction(3, 5),
) -> dict[str, Any]:
    bounds = cell["bounds"]
    x_bounds = decimal_box(*map(str, bounds["x"]))
    y_bounds = decimal_box(*map(str, bounds["y"]))
    first_sine = decimal_box(*map(str, bounds["first"]))
    second_sine = decimal_box(*map(str, bounds["second"]))
    certifier = RankRankCertifier(
        weight,
        level,
        x_bounds,
        str(cell["first_kind"]),
        str(cell["second_kind"]),
    )
    certificate = certifier.certify(
        root_from_bounds(x_bounds, y_bounds, first_sine, second_sine),
        max_boxes,
    )
    return {
        **{key: value for key, value in cell.items() if key != "bounds"},
        "bounds": {key: list(value) for key, value in bounds.items()},
        "certificate": certificate,
    }


def cover_payload(
    level: Fraction,
    total_cells: int,
    results: list[dict[str, Any]],
    *,
    weight: Fraction = Fraction(3, 5),
) -> dict[str, Any]:
    complete = [item for item in results if item["certificate"]["complete"]]
    by_topology: dict[str, dict[str, int | bool]] = {}
    for first_kind, second_kind in TOPOLOGIES:
        topology = f"{first_kind}/{second_kind}"
        selected = [
            item for item in results
            if item["certificate"]["topology"] == topology
        ]
        expected = sum(
            1 for item in topology_cells()
            if item["first_kind"] == first_kind
            and item["second_kind"] == second_kind
        )
        closed = sum(bool(item["certificate"]["complete"]) for item in selected)
        by_topology[topology] = {
            "expected_cells": expected,
            "processed_cells": len(selected),
            "complete_cells": closed,
            "complete": len(selected) == expected and closed == expected,
        }
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
        "topologies": by_topology,
        "cells": results,
        "scope": (
            "endpoint/endpoint and both mixed projective topologies on the "
            "complement-symmetric chart x+y<=1"
        ),
        "symmetry": (
            "E -> I-E exchanges endpoint/rank with rank/endpoint and maps "
            "the complementary chart into x+y<=1"
        ),
        "archived_solver_role": (
            "geometric partition only; archived primal and dual values are unused"
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
    parser.add_argument("--limit", type=int)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    weight = Fraction(args.weight)
    level = Fraction(args.level)
    if not Fraction(0) < weight < level < Fraction(1):
        parser.error("require 0 < weight < level < 1")
    cells = topology_cells()
    if args.limit is not None:
        cells = cells[: args.limit]

    prior_by_name: dict[str, dict[str, Any]] = {}
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
            }

    results: list[dict[str, Any]] = []
    for index, cell in enumerate(cells):
        previous = prior_by_name.get(str(cell["certificate_name"]))
        if previous is not None and previous["certificate"]["complete"]:
            result = previous
            reused = True
        else:
            result = certify_cell(
                cell, level, args.max_boxes_per_cell, weight=weight
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
        # A reused certificate is already durable in the resume artifact.
        # Checkpoint only new work; the complete payload is written once more
        # after the loop, which also performs cheap schema migrations.
        if args.output is not None and not reused:
            write_payload(
                args.output,
                cover_payload(level, len(cells), results, weight=weight),
            )

    payload = cover_payload(level, len(cells), results, weight=weight)
    if args.output is not None:
        write_payload(args.output, payload)
    print(
        json.dumps(
            {key: value for key, value in payload.items() if key != "cells"},
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

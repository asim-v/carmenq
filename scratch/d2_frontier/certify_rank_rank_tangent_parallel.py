"""Parallel, resumable replay of the rank/rank tangent-interval cover.

Each worker receives one complete geometric cell and returns its independent
certificate.  The parent process alone writes checkpoints, so parallelism
does not enlarge the mathematical trust boundary or introduce shared state.
"""

from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from fractions import Fraction
import json
import math
from pathlib import Path
from typing import Any

import certify_rank_rank_tangent_cover as serial


def certify_task(
    cell: dict[str, Any],
    level: Fraction,
    max_boxes: int,
    weight: Fraction,
    resume_certificate: dict[str, object] | None,
) -> dict[str, Any]:
    """Pickle-friendly entry point for an independent cell proof."""

    return serial.certify_cell(
        cell,
        level,
        max_boxes,
        weight=weight,
        resume_certificate=resume_certificate,
    )


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
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be positive")

    weight = Fraction(args.weight)
    level = Fraction(args.level)
    if not Fraction(0) < weight < level < Fraction(1):
        parser.error("require 0 < weight < level < 1")
    cells = [
        child
        for leaf in serial.geometric_leaves()
        if float(leaf["archived_dual_scaled"]) >= args.minimum_archived_dual
        for child in serial.expand_full_angles(leaf)
    ]
    cells.sort(key=lambda item: float(item["archived_dual_scaled"]), reverse=True)
    if args.cell:
        try:
            cells = serial.select_named_cells(cells, args.cell)
        except ValueError as exc:
            parser.error(str(exc))
    if args.limit is not None:
        cells = cells[: args.limit]

    prior_by_id: dict[str, dict[str, Any]] = {}
    legacy_by_name: dict[str, dict[str, Any]] = {}
    name_counts = Counter(str(cell["name"]) for cell in cells)
    if args.resume and args.output.exists():
        prior = json.loads(args.output.read_text(encoding="utf-8"))
        if prior.get("weight") != str(weight) or prior.get("level") != str(level):
            parser.error("resume file has incompatible weight or level")
        prior_by_id = {
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

    results_by_id: dict[str, dict[str, Any]] = {}
    checkpoint_by_id: dict[str, dict[str, Any]] = {}
    pending: list[dict[str, Any]] = []
    resume_by_id: dict[str, dict[str, object]] = {}
    for cell in cells:
        certificate_id = str(cell["certificate_name"])
        previous = prior_by_id.get(certificate_id)
        if previous is None:
            previous = legacy_by_name.get(str(cell["name"]))
        if previous is not None and previous["certificate"]["complete"]:
            results_by_id[certificate_id] = {
                **previous,
                "certificate_name": certificate_id,
            }
        else:
            pending.append(cell)
            if previous is not None:
                checkpoint_by_id[certificate_id] = {
                    **previous,
                    "certificate_name": certificate_id,
                }
            if (
                previous is not None
                and previous["certificate"].get("open_frontier")
            ):
                resume_by_id[certificate_id] = previous["certificate"]
        if certificate_id in results_by_id:
            checkpoint_by_id[certificate_id] = results_by_id[certificate_id]

    def ordered_results() -> list[dict[str, Any]]:
        return [
            checkpoint_by_id[str(cell["certificate_name"])]
            for cell in cells
            if str(cell["certificate_name"]) in checkpoint_by_id
        ]

    serial.write_payload(
        args.output,
        serial.cover_payload(
            level, len(cells), ordered_results(), weight=weight
        ),
    )
    print(
        json.dumps(
            {
                "cells": len(cells),
                "reused_complete": len(results_by_id),
                "pending": len(pending),
                "workers": args.workers,
                "resumable_incomplete": len(resume_by_id),
            }
        ),
        flush=True,
    )

    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(
                certify_task,
                cell,
                level,
                args.max_boxes_per_cell,
                weight,
                resume_by_id.get(str(cell["certificate_name"])),
            ): cell
            for cell in pending
        }
        for completed_now, future in enumerate(as_completed(futures), start=1):
            cell = futures[future]
            result = future.result()
            certificate_id = str(cell["certificate_name"])
            results_by_id[certificate_id] = result
            checkpoint_by_id[certificate_id] = result
            certificate = result["certificate"]
            print(
                json.dumps(
                    {
                        "completed_now": completed_now,
                        "pending_count": len(pending),
                        "processed_total": len(results_by_id),
                        "count": len(cells),
                        "name": certificate_id,
                        "complete": certificate["complete"],
                        "split": certificate["boxes_split"],
                        "remaining": certificate["boxes_remaining"],
                        "maximum_open_upper": certificate["maximum_open_upper"],
                    }
                ),
                flush=True,
            )
            serial.write_payload(
                args.output,
                serial.cover_payload(
                    level, len(cells), ordered_results(), weight=weight
                ),
            )

    payload = serial.cover_payload(
        level, len(cells), ordered_results(), weight=weight
    )
    serial.write_payload(args.output, payload)
    print(
        json.dumps(
            {key: value for key, value in payload.items() if key != "cells"},
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

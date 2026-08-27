"""Merge compatible partial rank/rank tangent-cover checkpoints.

Complete cell certificates dominate incomplete ones.  Between two incomplete
records, the checkpoint with the larger split count is retained, but only if
it carries the resumable open frontier.  The merged payload is ordered by the
canonical 944-cell catalogue and is re-summarised from its cell records.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
from typing import Any

import certify_rank_rank_tangent_cover as cover


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def preferred(
    first: dict[str, Any], second: dict[str, Any]
) -> dict[str, Any]:
    first_certificate = first["certificate"]
    second_certificate = second["certificate"]
    first_complete = bool(first_certificate["complete"])
    second_complete = bool(second_certificate["complete"])
    if first_complete != second_complete:
        return first if first_complete else second
    if first_complete:
        return min(
            (first, second),
            key=lambda row: int(row["certificate"]["boxes_split"]),
        )
    candidates = [
        row
        for row in (first, second)
        if row["certificate"].get("open_frontier")
    ]
    if not candidates:
        return max(
            (first, second),
            key=lambda row: int(row["certificate"]["boxes_split"]),
        )
    return max(
        candidates,
        key=lambda row: int(row["certificate"]["boxes_split"]),
    )


def merge_payloads(payloads: list[dict[str, Any]]) -> dict[str, Any]:
    require(bool(payloads), "no cover payloads supplied")
    weight = Fraction(payloads[0]["weight"])
    level = Fraction(payloads[0]["level"])
    expected_cells = [
        child
        for leaf in cover.geometric_leaves()
        for child in cover.expand_full_angles(leaf)
    ]
    expected_by_id = {
        str(cell["certificate_name"]): cell for cell in expected_cells
    }
    records: dict[str, dict[str, Any]] = {}
    for payload in payloads:
        require(Fraction(payload["weight"]) == weight, "weight mismatch")
        require(Fraction(payload["level"]) == level, "level mismatch")
        for record in payload.get("cells", []):
            certificate_id = str(record["certificate_name"])
            require(certificate_id in expected_by_id, "unknown certificate cell")
            require(
                record["bounds"]
                == {
                    key: list(value)
                    for key, value in expected_by_id[certificate_id]["bounds"].items()
                },
                f"bounds changed for {certificate_id}",
            )
            if certificate_id in records:
                record = preferred(records[certificate_id], record)
            records[certificate_id] = record
    ordered = [
        records[str(cell["certificate_name"])]
        for cell in expected_cells
        if str(cell["certificate_name"]) in records
    ]
    return cover.cover_payload(level, len(expected_cells), ordered, weight=weight)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payloads = [
        json.loads(path.read_text(encoding="utf-8")) for path in args.inputs
    ]
    merged = merge_payloads(payloads)
    cover.write_payload(args.output, merged)
    print(
        json.dumps(
            {key: value for key, value in merged.items() if key != "cells"},
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

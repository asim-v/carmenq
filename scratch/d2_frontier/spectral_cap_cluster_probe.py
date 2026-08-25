"""Certify a Cartesian cluster of adjacent spectral cap cells at once."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

from spectral_product_localizer_batch import (
    build_localisation_oracle,
    build_product_oracle,
    cover_localised_cell,
    enclosing_scaled_cap,
    localise_cell,
    pattern_code,
)


SCHEMA = "carmenq.spectral-cap-cluster-probe.v1"


def parse_pattern(value: str) -> tuple[str, ...]:
    symbols = {"b": "bloch", "+": "scalar-positive", "-": "scalar-negative"}
    try:
        pattern = tuple(symbols[symbol] for symbol in value)
    except KeyError as error:
        raise argparse.ArgumentTypeError(
            f"unknown branch symbol {error.args[0]!r}"
        ) from error
    if len(pattern) != 4:
        raise argparse.ArgumentTypeError("the branch pattern must have length four")
    return pattern


def parse_caps(
    value: str,
    pattern: tuple[str, ...],
) -> tuple[int | tuple[int, ...] | None, ...]:
    fields = value.split(";")
    if len(fields) != len(pattern):
        raise argparse.ArgumentTypeError("caps must have one semicolon field per branch")
    result: list[int | tuple[int, ...] | None] = []
    for field, branch in zip(fields, pattern, strict=True):
        if branch != "bloch":
            if field.strip() not in {"", "none"}:
                raise argparse.ArgumentTypeError("scalar branches use an empty cap field")
            result.append(None)
            continue
        indices = tuple(int(item) for item in field.split(",") if item.strip())
        if not indices:
            raise argparse.ArgumentTypeError("Bloch branches need at least one cap")
        result.append(indices[0] if len(indices) == 1 else indices)
    return tuple(result)


def cap_audit(
    grids: tuple[int, ...],
    pattern: tuple[str, ...],
    caps: tuple[int | tuple[int, ...] | None, ...],
) -> list[dict[str, Any] | None]:
    result: list[dict[str, Any] | None] = []
    for grid, branch, indices in zip(grids, pattern, caps, strict=True):
        if branch != "bloch":
            result.append(None)
            continue
        children = indices if isinstance(indices, tuple) else (int(indices),)
        _, audit = enclosing_scaled_cap(grid, children)
        result.append(audit)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--frontier-json", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--pattern", type=parse_pattern, required=True)
    parser.add_argument(
        "--caps",
        required=True,
        help="semicolon-separated cap sets, with comma-separated child indices",
    )
    parser.add_argument("--coordinate-safety", type=float, default=2e-6)
    parser.add_argument("--bound-safety", type=float, default=2e-6)
    parser.add_argument("--minimum-width", type=float, default=1e-6)
    parser.add_argument("--max-nodes", type=int, default=100)
    args = parser.parse_args()

    source_raw = args.frontier_json.read_bytes()
    source = json.loads(source_raw)
    pattern = tuple(args.pattern)
    caps = parse_caps(args.caps, pattern)
    cluster_size = 1
    for branch, indices in zip(pattern, caps, strict=True):
        if branch == "bloch" and isinstance(indices, tuple):
            cluster_size *= len(indices)
    cell = {
        "source_index": -1,
        "source_cell": -1,
        "branches": pattern,
        "caps": caps,
        "source_status": "cluster",
        "source_bound": max(
            float(item["bound"])
            for item in source["cells"]
            if str(item["status"]) in {"optimal", "optimal_inaccurate"}
        ),
        "source_audit": 0.0,
        "source_return": 0.0,
    }

    localisation_oracle, localisation_caps, box, localisation_build = (
        build_localisation_oracle(source, pattern)
    )
    localisation = localise_cell(
        source,
        cell,
        localisation_oracle,
        localisation_caps,
        box,
        args.coordinate_safety,
    )
    cover: dict[str, Any] | None = None
    product_build: dict[str, Any] | None = None
    if localisation["status"] == "localized":
        (
            product_oracle,
            product_caps,
            product_box,
            lower,
            upper,
            purity,
            product_build,
        ) = build_product_oracle(source, pattern)
        cover = cover_localised_cell(
            source,
            cell,
            localisation,
            product_oracle,
            product_caps,
            product_box,
            lower,
            upper,
            purity,
            args.max_nodes,
            args.bound_safety,
            args.minimum_width,
        )

    payload = {
        "schema": SCHEMA,
        "source": {
            "path": str(args.frontier_json),
            "sha256": hashlib.sha256(source_raw).hexdigest(),
            "target": float(source["target"]),
        },
        "cluster": {
            "pattern": pattern_code(pattern),
            "caps": caps,
            "cartesian_child_cell_count": cluster_size,
            "cap_audit": cap_audit(
                tuple(int(value) for value in source["separator_grids"]),
                pattern,
                caps,
            ),
        },
        "configuration": {
            "coordinate_safety": args.coordinate_safety,
            "bound_safety": args.bound_safety,
            "minimum_width": args.minimum_width,
            "max_nodes": args.max_nodes,
            "state_choi_psd": True,
            "state_choi_ppt": True,
        },
        "builds": {
            "localisation": localisation_build,
            "product": product_build,
        },
        "localisation": localisation,
        "cover": cover,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "pattern": pattern_code(pattern),
                "cartesian_child_cell_count": cluster_size,
                "localisation_status": localisation["status"],
                "solved_nodes": cover.get("solved_nodes") if cover else None,
                "cover_upper_bound": cover.get("cover_upper_bound") if cover else None,
                "target_closed": cover.get("target_closed") if cover else None,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

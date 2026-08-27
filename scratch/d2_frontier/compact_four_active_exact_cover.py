"""Compact a complete four-active exact-dual cover without weakening replay.

The production artifact records every derived coefficient enclosure and
canonical-matrix statistic beside every dual.  Those fields are useful while
developing the cover but are redundant for distribution: the solver-free
verifier reconstructs them from the exact box and current source.  The compact
format retains the complete leaf tree and, for each of the six affine orders,
only the permutation, storage dtype, and compressed dual vector.
"""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path
from typing import Any


COMPACT_SCHEMA = "carmenq.four-active-mccormick-exact-dual-cover.compact.v1"


def read_json(path: Path) -> dict[str, Any]:
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8") as stream:
            return json.load(stream)
    return json.loads(path.read_text(encoding="utf-8"))


def compact_leaf(leaf: dict[str, Any]) -> dict[str, Any]:
    result: dict[str, Any] = {
        "kind": leaf["kind"],
        "box": leaf["box"],
    }
    if leaf["kind"] == "domain-empty":
        return result
    if leaf["kind"] != "closed":
        raise RuntimeError("compact format accepts only closed or empty leaves")
    result["order_duals"] = [
        {
            "syndrome_permutation": report["syndrome_permutation"],
            "dual_storage_dtype": report["dual_storage_dtype"],
            "dual_zlib_base64": report["dual_zlib_base64"],
        }
        for report in leaf["order_certificates"]
    ]
    return result


def compact_artifact(source: dict[str, Any]) -> dict[str, Any]:
    if not source.get("complete") or not source.get("all_cells_closed"):
        raise RuntimeError("refusing to compact an incomplete cover")
    if source.get("open_boxes"):
        raise RuntimeError("complete cover unexpectedly retains open boxes")
    retained = (
        "support_weight",
        "target",
        "maximum_weight_floor",
        "minimum_active_weight",
        "projective_lines",
        "reserve_perturbations",
        "common_bias_coefficient_representatives",
        "common_bias_coefficient_orbit_size",
        "weighted_reserve_geometry",
        "projective_pair_geometry",
        "prefix_order_reduction",
        "relaxation",
        "initial_box",
        "boxes_split",
        "leaf_count",
        "closed_leaf_count",
        "domain_empty_leaf_count",
    )
    compact = {key: source[key] for key in retained}
    compact.update(
        {
            "schema": COMPACT_SCHEMA,
            "source_schema": source["schema"],
            "complete": True,
            "all_cells_closed": True,
            "boxes_remaining": 0,
            "open_boxes": [],
            "leaves": [compact_leaf(leaf) for leaf in source["leaves"]],
            "trusted_optimizers": [],
        }
    )
    return compact


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".gz":
        with gzip.open(path, "wt", encoding="utf-8", compresslevel=9) as stream:
            json.dump(payload, stream, separators=(",", ":"))
        return
    path.write_text(
        json.dumps(payload, separators=(",", ":")) + "\n", encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    compact = compact_artifact(read_json(args.source))
    write_json(args.output, compact)
    print(json.dumps({
        "leaf_count": compact["leaf_count"],
        "output": str(args.output),
        "bytes": args.output.stat().st_size,
    }))


if __name__ == "__main__":
    main()

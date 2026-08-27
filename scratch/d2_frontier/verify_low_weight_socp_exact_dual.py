"""Solver-free verifier for the low-weight SOCP dual certificate."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
from typing import Any

import numpy as np

from four_active_socp_exact_cover import decode_dual, fraction_pair
from low_weight_socp_exact_dual import (
    CAP_VECTOR,
    MAXIMUM_EFFECT_WEIGHT,
    ORDERS,
    SUPPORT_WEIGHT,
    canonical_data,
    exact_upper,
)
from ternary_socp_exact_dual_probe import canonical_hash, repair_dual_cones


ROOT = Path(__file__).resolve().parent
REPOSITORY_ROOT = ROOT.parents[1]
DEFAULT_CERTIFICATE = ROOT / "low_weight_socp_exact_dual_l060.json"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def portable_path(path: Path) -> str:
    """Return a repository-relative POSIX path or reject external input."""
    try:
        return path.resolve().relative_to(REPOSITORY_ROOT).as_posix()
    except ValueError as error:
        raise RuntimeError(f"artifact lies outside the repository: {path}") from error


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def dimensions(data: dict[str, Any]) -> dict[str, Any]:
    return {
        "zero": int(data["dims"].zero),
        "nonnegative": int(data["dims"].nonneg),
        "soc": list(map(int, data["dims"].soc)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--certificate", type=Path, default=DEFAULT_CERTIFICATE)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    artifact = json.loads(args.certificate.read_text(encoding="utf-8"))
    require(
        artifact.get("schema") == "carmenq.low-weight-socp-exact-dual.v1",
        "wrong certificate schema",
    )
    require(
        artifact.get("support_weight") == fraction_pair(SUPPORT_WEIGHT),
        "wrong support weight",
    )
    require(
        artifact.get("maximum_effect_weight")
        == fraction_pair(MAXIMUM_EFFECT_WEIGHT),
        "wrong effect-weight cap",
    )
    require(
        artifact.get("capped_simplex_extreme_vector")
        == list(map(fraction_pair, CAP_VECTOR)),
        "wrong capped-simplex support vector",
    )
    target = Fraction(*artifact["target_fraction"])
    expected = set(itertools.product(ORDERS, repeat=2))
    cells = artifact["cells"]
    require(len(cells) == len(expected), "certificate has the wrong cell count")
    seen: set[
        tuple[tuple[int, int, int, int], tuple[int, int, int, int]]
    ] = set()
    maximum = Fraction(0)
    repaired_blocks = 0
    storage_counts: dict[str, int] = {}
    for index, cell in enumerate(cells):
        prefix_order = tuple(map(int, cell["prefix_order"]))
        syndrome_order = tuple(map(int, cell["syndrome_order"]))
        key = (prefix_order, syndrome_order)
        require(key in expected, f"invalid order pair at cell {index}")
        require(key not in seen, f"duplicate order pair at cell {index}")
        seen.add(key)
        data = canonical_data(prefix_order, syndrome_order)
        require(
            cell["canonical_shape"]
            == [int(data["A"].shape[0]), int(data["A"].shape[1])],
            f"canonical shape mismatch at cell {index}",
        )
        require(
            int(cell["canonical_nonzeros"]) == int(data["A"].nnz),
            f"canonical nonzero mismatch at cell {index}",
        )
        require(
            cell["canonical_sha256"] == canonical_hash(data),
            f"canonical hash mismatch at cell {index}",
        )
        require(
            cell["cone_dimensions"] == dimensions(data),
            f"cone dimensions mismatch at cell {index}",
        )
        dtype = str(cell["dual_storage_dtype"])
        require(dtype in {"f32", "f64"}, f"invalid dual dtype at cell {index}")
        dual = decode_dual(cell["dual_zlib_base64"], dtype)
        require(
            dual.shape == (data["A"].shape[0],) and np.all(np.isfinite(dual)),
            f"invalid dual vector at cell {index}",
        )
        dual, repaired = repair_dual_cones(dual, data["dims"])
        repaired_blocks += repaired
        upper, correction, _ = exact_upper(data, dual)
        require(
            fraction_pair(upper) == cell["certified_upper_fraction"],
            f"exact upper mismatch at cell {index}",
        )
        require(
            fraction_pair(correction) == cell["exact_residual_correction"],
            f"residual correction mismatch at cell {index}",
        )
        require(upper <= target, f"cell {index} exceeds target")
        require(bool(cell["closed"]), f"cell {index} is not marked closed")
        maximum = max(maximum, upper)
        storage_counts[dtype] = storage_counts.get(dtype, 0) + 1
    require(seen == expected, "order enumeration is incomplete")
    require(bool(artifact["complete"]), "artifact is not marked complete")
    require(
        fraction_pair(maximum) == artifact["maximum_certified_upper_fraction"],
        "global maximum mismatch",
    )
    summary = {
        "schema": "carmenq.low-weight-socp-exact-dual-verification.v1",
        "certificate": portable_path(args.certificate),
        "verified_cells": len(cells),
        "certificate_sha256": sha256(args.certificate),
        "expected_cells": len(expected),
        "maximum_certified_upper_fraction": fraction_pair(maximum),
        "maximum_certified_upper_decimal": artifact[
            "maximum_certified_upper_decimal"
        ],
        "target_fraction": fraction_pair(target),
        "dual_storage_counts": storage_counts,
        "additional_soc_heads_repaired": repaired_blocks,
        "optimiser_called": False,
        "complete": True,
    }
    rendered = json.dumps(summary, indent=2) + "\n"
    print(rendered, end="")
    if args.output is not None:
        output = args.output
        if not output.is_absolute():
            output = ROOT / output
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

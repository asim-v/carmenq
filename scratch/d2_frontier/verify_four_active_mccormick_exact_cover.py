"""Solver-free verifier for the four-active weighted-reserve cover."""

from __future__ import annotations

import argparse
import gzip
import hashlib
from fractions import Fraction
import json
from pathlib import Path
from typing import Any

import numpy as np

from compact_four_active_exact_cover import COMPACT_SCHEMA
from four_active_mccormick_socp_exact_cover import (
    COMMON_BIAS_COEFFICIENT_REPRESENTATIVES,
    COMMON_BIAS_COEFFICIENTS,
    FOUR_ACTIVE_TARGET,
    MAXIMUM_WEIGHT_FLOOR,
    MINIMUM_ACTIVE_WEIGHT,
    PERTURBATIONS,
    NONZERO_PERMUTATIONS,
    PROJECTIVE_LINES,
    ROOT,
    SUPPORT_WEIGHT,
    SCHEMA,
    McCormickOracle,
    WeightBox,
    canonical_hash,
    decode_dual,
    exact_upper,
    fraction_decimal,
    fraction_pair,
    initial_box,
    physical_weight_vertices,
    repair_dual_cones,
    validate_leaf_tree,
    weight_hull,
)
REPOSITORY_ROOT = ROOT.parents[1]


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


def verify_order(
    oracle: McCormickOracle,
    box: WeightBox,
    hull: tuple[tuple[Fraction, Fraction], ...],
    report: dict[str, Any],
    target: Fraction,
    *,
    compact: bool,
) -> Fraction:
    compact_fields = {
        "syndrome_permutation",
        "dual_storage_dtype",
        "dual_zlib_base64",
    }
    full_fields = compact_fields | {
        "coefficient_enclosure",
        "canonical_shape",
        "canonical_nonzeros",
        "canonical_sha256",
        "cone_dimensions",
        "certified_upper_fraction",
        "certified_upper_decimal",
        "exact_residual_correction",
        "maximum_stationarity_residual_decimal",
        "closed",
        "trusted_optimizers",
    }
    required = compact_fields if compact else full_fields
    missing = sorted(required - report.keys())
    require(not missing, f"order report is missing required fields: {missing}")
    permutation = tuple(map(int, report["syndrome_permutation"]))
    require(permutation in NONZERO_PERMUTATIONS, "invalid syndrome permutation")
    enclosure = oracle.assign(hull, permutation, physical_weight_vertices(box))
    if not compact:
        require(
            enclosure == report["coefficient_enclosure"],
            "coefficient enclosure does not replay",
        )
    data = oracle.canonical_data()
    if not compact:
        require(
            list(map(int, data["A"].shape)) == report["canonical_shape"],
            "canonical shape mismatch",
        )
        require(int(data["A"].nnz) == report["canonical_nonzeros"], "nnz mismatch")
        require(canonical_hash(data) == report["canonical_sha256"], "matrix hash mismatch")
    dimensions = {
        "zero": int(data["dims"].zero),
        "nonnegative": int(data["dims"].nonneg),
        "soc": list(map(int, data["dims"].soc)),
    }
    if not compact:
        require(dimensions == report["cone_dimensions"], "cone dimensions changed")
    dtype = str(report["dual_storage_dtype"])
    dual = decode_dual(str(report["dual_zlib_base64"]), dtype)
    dual, _ = repair_dual_cones(np.asarray(dual), data["dims"])
    upper, correction, residual = exact_upper(data, dual)
    if not compact:
        require(
            upper == Fraction(*report["certified_upper_fraction"]),
            "exact objective does not replay",
        )
        require(
            fraction_decimal(upper) == report["certified_upper_decimal"],
            "upper decimal mismatch",
        )
        require(
            correction == Fraction(*report["exact_residual_correction"]),
            "stationarity correction mismatch",
        )
        require(
            fraction_decimal(residual)
            == report["maximum_stationarity_residual_decimal"],
            "maximum stationarity residual mismatch",
        )
        require(bool(report["closed"]), "order report is not marked closed")
    require(upper <= target, "order report exceeds target")
    if not compact:
        require(not report["trusted_optimizers"], "report trusts an optimizer")
    return upper


def load_artifact(path: Path) -> dict[str, Any]:
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8") as stream:
            return json.load(stream)
    return json.loads(path.read_text(encoding="utf-8"))


def verify_artifact(path: Path, target: Fraction) -> dict[str, Any]:
    artifact = load_artifact(path)
    compact = artifact.get("schema") == COMPACT_SCHEMA
    require(compact or artifact.get("schema") == SCHEMA, "wrong four-active schema")
    if compact:
        require(artifact.get("source_schema") == SCHEMA, "wrong compact source schema")
    require(artifact.get("support_weight") == fraction_pair(SUPPORT_WEIGHT), "wrong weight")
    require(artifact.get("target") == fraction_pair(target), "wrong target")
    require(
        artifact.get("maximum_weight_floor") == fraction_pair(MAXIMUM_WEIGHT_FLOOR),
        "wrong maximum-weight floor",
    )
    require(
        artifact.get("reserve_perturbations") == list(map(fraction_pair, PERTURBATIONS)),
        "wrong weighted-reserve perturbations",
    )
    require(
        artifact.get("common_bias_coefficient_representatives")
        == [
            list(map(fraction_pair, representative))
            for representative in COMMON_BIAS_COEFFICIENT_REPRESENTATIVES
        ],
        "wrong all-prior common-bias support directions",
    )
    require(
        int(artifact.get("common_bias_coefficient_orbit_size", -1))
        == len(COMMON_BIAS_COEFFICIENTS),
        "wrong all-prior common-bias orbit size",
    )
    require(
        artifact.get("projective_pair_geometry")
        == "all ordered effect-pair projective comparisons from full Bloch closure",
        "wrong full-closure pair geometry",
    )
    require(
        artifact.get("weighted_reserve_geometry")
        == "exact correlated weight-polytope vertices",
        "wrong common-bias weight geometry",
    )
    require(
        artifact.get("minimum_active_weight") == fraction_pair(MINIMUM_ACTIVE_WEIGHT),
        "wrong minimum active weight",
    )
    require(
        artifact.get("projective_lines")
        == [
            [fraction_pair(weight), fraction_pair(upper)]
            for weight, upper in PROJECTIVE_LINES
        ],
        "projective support lines changed",
    )
    require(artifact.get("initial_box") == initial_box().serialise(), "wrong root")
    require(bool(artifact.get("complete")), "cover is incomplete")
    require(bool(artifact.get("all_cells_closed")), "cover has open cells")
    require(int(artifact.get("boxes_remaining", -1)) == 0, "boxes remain")
    require(artifact.get("open_boxes") == [], "open frontier is not empty")
    leaves = artifact["leaves"]
    require(int(artifact["leaf_count"]) == len(leaves), "leaf count mismatch")
    validate_leaf_tree(leaves)

    oracle = McCormickOracle()
    maximum: Fraction | None = None
    closed_count = 0
    empty_count = 0
    for index, leaf in enumerate(leaves, start=1):
        box = WeightBox.deserialise(leaf["box"])
        hull = weight_hull(box)
        if leaf["kind"] == "domain-empty":
            require(hull is None, "invalid domain-empty leaf")
            if compact:
                require(
                    set(leaf) == {"kind", "box"},
                    "compact empty leaf has unexpected fields",
                )
            empty_count += 1
            continue
        require(leaf["kind"] == "closed", "unexpected leaf kind")
        require(hull is not None, "closed leaf has empty domain")
        if compact:
            require(
                set(leaf) == {"kind", "box", "order_duals"},
                "compact closed leaf has unexpected fields",
            )
            reports = leaf["order_duals"]
        else:
            expected_hull = [
                [fraction_pair(lower), fraction_pair(upper)]
                for lower, upper in hull
            ]
            require(leaf["exact_weight_hull"] == expected_hull, "weight hull changed")
            reports = leaf["order_certificates"]
        require(len(reports) == len(NONZERO_PERMUTATIONS), "missing order reports")
        require(
            {tuple(row["syndrome_permutation"]) for row in reports}
            == set(NONZERO_PERMUTATIONS),
            "order reports are not exhaustive",
        )
        bounds = [
            verify_order(oracle, box, hull, report, target, compact=compact)
            for report in reports
        ]
        leaf_maximum = max(bounds)
        if not compact:
            require(
                leaf_maximum == Fraction(*leaf["maximum_certified_upper_fraction"]),
                "leaf maximum mismatch",
            )
            require(
                fraction_decimal(leaf_maximum)
                == leaf["maximum_certified_upper_decimal"],
                "leaf maximum decimal mismatch",
            )
        maximum = leaf_maximum if maximum is None else max(maximum, leaf_maximum)
        closed_count += 1
        if index % 25 == 0:
            print(
                json.dumps({"verified": index, "count": len(leaves)}),
                flush=True,
            )

    require(int(artifact["closed_leaf_count"]) == closed_count, "closed count mismatch")
    require(
        int(artifact["domain_empty_leaf_count"]) == empty_count,
        "domain-empty count mismatch",
    )
    require(maximum is not None, "cover contains no certified leaf")
    if not compact:
        require(
            fraction_decimal(maximum) == artifact["maximum_certified_upper_decimal"],
            "artifact maximum mismatch",
        )
    return {
        "schema": "carmenq.four-active-mccormick-exact-verification.v1",
        "source": portable_path(path),
        "target": str(target),
        "source_sha256": sha256(path),
        "leaf_count": len(leaves),
        "closed_leaf_count": closed_count,
        "domain_empty_leaf_count": empty_count,
        "maximum_certified_upper": fraction_decimal(maximum),
        "compact_source": compact,
        "named_projective_premises": {
            str(weight): str(upper)
            for weight, upper in PROJECTIVE_LINES
        },
        "tree_complete": True,
        "all_six_orders_verified": True,
        "verified_without_solver": True,
        "trusted_optimizers": [],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "artifact",
        nargs="?",
        type=Path,
        default=ROOT / "four_active_weighted_reserve_exact_cover_l060.json",
    )
    parser.add_argument("--target", default=str(FOUR_ACTIVE_TARGET))
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = verify_artifact(args.artifact, Fraction(args.target))
    rendered = json.dumps(result, indent=2) + "\n"
    print(rendered, end="")
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

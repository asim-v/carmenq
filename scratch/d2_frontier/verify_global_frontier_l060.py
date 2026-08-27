"""Assemble the solver-independent lambda=0.6 global-sector enclosure.

This checker performs the final exact arithmetic and validates the summaries
emitted by the four independent proof kernels:

* outward interval projective covers at lambda 0.55 and 0.60;
* exact-residual dual replay of the capped low-weight SOCP;
* exact-residual dual replay of all ternary cells; and
* exact-residual dual replay of the four-active McCormick cover.

The component verifiers remain the proof replay entry points.  This script
does not call an optimiser and refuses incomplete or parameter-inconsistent
summaries.  The deletion sector fixes the final reported endpoint after the
two sharp projective premises close the ternary and four-active replays.
"""

from __future__ import annotations

import argparse
from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
REPOSITORY_ROOT = ROOT.parent.parent
SUPPORT_WEIGHT = Fraction(3, 5)
PROJECTIVE_055_WEIGHT = Fraction(11, 20)
PROJECTIVE_055_UPPER = Fraction(7573, 10000)
PROJECTIVE_060_UPPER = Fraction(766, 1000)
FOUR_ACTIVE_UPPER = Fraction(76670, 100000)
LOW_WEIGHT_TARGET = Fraction(76591, 100000)
TERNARY_UPPER = Fraction(76652, 100000)
MAXIMUM_WEIGHT_FLOOR = Fraction(3533, 4000)
MINIMUM_ACTIVE_WEIGHT = Fraction(3, 10000)
DELETION_UPPER = TERNARY_UPPER + SUPPORT_WEIGHT * MINIMUM_ACTIVE_WEIGHT
REPORTED_UPPER = Fraction(76670, 100000)
DECLARED_PHYSICAL_LOWER = Fraction(7_658_988_152, 10_000_000_000)
EXPECTED_PROJECTIVE_PREMISES = {
    str(PROJECTIVE_055_WEIGHT): str(PROJECTIVE_055_UPPER),
    str(SUPPORT_WEIGHT): str(PROJECTIVE_060_UPPER),
}


def load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"), parse_float=Decimal)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def fraction_pair(value: Fraction) -> list[int]:
    return [value.numerator, value.denominator]


def decimal(value: Fraction, digits: int = 40) -> str:
    with localcontext() as context:
        context.prec = digits
        return format(Decimal(value.numerator) / Decimal(value.denominator), "f")


def as_fraction(value: Any) -> Fraction:
    if isinstance(value, list) and len(value) == 2:
        return Fraction(int(value[0]), int(value[1]))
    return Fraction(str(value))


def validate_named_projective_premises(
    payload: dict[str, Any], component: str
) -> None:
    require(
        payload.get("named_projective_premises") == EXPECTED_PROJECTIVE_PREMISES,
        f"{component} uses incompatible projective premises",
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def portable_path(path: Path) -> str:
    """Return a repository-relative POSIX path or reject the artifact."""

    resolved = path.resolve()
    try:
        return resolved.relative_to(REPOSITORY_ROOT).as_posix()
    except ValueError as error:
        raise RuntimeError(
            f"certificate artifact lies outside repository root: {resolved.name}"
        ) from error

def validate_projective(
    payload: dict[str, Any], weight: Fraction, level: Fraction
) -> dict[str, Any]:
    require(
        payload.get("schema") == "carmenq.projective-tangent-global-summary.v1",
        "wrong projective summary schema",
    )
    require(Fraction(payload["weight"]) == weight, "wrong projective weight")
    require(Fraction(payload["level"]) == level, "wrong projective level")
    require(
        bool(payload["all_four_topologies_complete"]),
        "projective topology cover is incomplete",
    )
    require(
        int(payload["topology_cell_count"]) == 1448,
        "projective cover does not contain all 1,448 topology cells",
    )
    require(
        payload["proof_kernel"]
        == "outward-expanded IEEE-754 binary64 intervals",
        "wrong projective proof kernel",
    )
    return {
        "weight": str(weight),
        "upper": str(level),
        "topology_cells": int(payload["topology_cell_count"]),
        "boxes_split": int(payload["rank_rank"]["boxes_split"])
        + int(payload["remaining_topologies"]["boxes_split"]),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--physical-lower",
        type=Path,
        default=ROOT.parent.parent / "data" / "four_effect_rational_lower_l060.json",
    )
    parser.add_argument(
        "--projective-060",
        type=Path,
        default=ROOT / "projective_tangent_global_l060_L0766_summary.json",
    )
    parser.add_argument(
        "--projective-055",
        type=Path,
        default=ROOT / "projective_tangent_global_l055_L07573_summary.json",
    )
    parser.add_argument(
        "--low-weight-verification",
        type=Path,
        default=ROOT / "low_weight_socp_exact_dual_l060_verified.json",
    )
    parser.add_argument(
        "--ternary-verification",
        type=Path,
        default=ROOT / "ternary_transferred_exact_dual_full_l060_verified.json",
    )
    parser.add_argument(
        "--four-active-verification",
        type=Path,
        default=ROOT / "four_active_common_bias_fallback_exact_cover_l060_verified.json",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    paths = {
        "physical_lower": args.physical_lower,
        "projective_060": args.projective_060,
        "projective_055": args.projective_055,
        "low_weight_verification": args.low_weight_verification,
        "ternary_verification": args.ternary_verification,
        "four_active_verification": args.four_active_verification,
    }
    payloads = {name: load(path) for name, path in paths.items()}

    projective_060 = validate_projective(
        payloads["projective_060"], SUPPORT_WEIGHT, PROJECTIVE_060_UPPER
    )
    projective_055 = validate_projective(
        payloads["projective_055"],
        PROJECTIVE_055_WEIGHT,
        PROJECTIVE_055_UPPER,
    )

    low = payloads["low_weight_verification"]
    require(
        low.get("schema")
        == "carmenq.low-weight-socp-exact-dual-verification.v1",
        "wrong low-weight verification schema",
    )
    require(bool(low["complete"]), "low-weight replay is incomplete")
    require(
        int(low["verified_cells"]) == int(low["expected_cells"]) == 576,
        "low-weight order enumeration is incomplete",
    )
    require(not bool(low["optimiser_called"]), "low-weight replay called a solver")
    require(
        as_fraction(low["target_fraction"]) == LOW_WEIGHT_TARGET,
        "wrong low-weight target",
    )
    low_upper = as_fraction(low["maximum_certified_upper_fraction"])
    require(low_upper <= LOW_WEIGHT_TARGET, "low-weight cap exceeds target")

    ternary = payloads["ternary_verification"]
    require(
        ternary.get("schema")
        == "carmenq.ternary-socp-transferred-exact-dual-verification.v1",
        "wrong ternary verification schema",
    )
    require(Fraction(ternary["target"]) == TERNARY_UPPER, "wrong ternary target")
    require(
        bool(ternary["full_source_cover_verified"]),
        "ternary source cover is incomplete",
    )
    require(bool(ternary["all_cells_at_most_target"]), "ternary target failed")
    require(
        int(ternary["verified_cell_count"])
        == int(ternary["source_leaf_count"])
        == 12008,
        "ternary leaf enumeration is incomplete",
    )
    require(not ternary["trusted_optimizers"], "ternary replay trusts a solver")
    validate_named_projective_premises(ternary, "ternary replay")

    four = payloads["four_active_verification"]
    require(
        four.get("schema")
        == "carmenq.four-active-mccormick-exact-verification.v1",
        "wrong four-active verification schema",
    )
    require(Fraction(four["target"]) == FOUR_ACTIVE_UPPER, "wrong four-active target")
    require(bool(four["tree_complete"]), "four-active tree is incomplete")
    require(
        bool(four["all_six_orders_verified"]),
        "four-active affine orders are incomplete",
    )
    require(
        bool(four["verified_without_solver"]),
        "four-active replay used an optimizer",
    )
    require(not four["trusted_optimizers"], "four-active replay trusts a solver")
    validate_named_projective_premises(four, "four-active replay")
    four_upper = Fraction(four["maximum_certified_upper"])
    require(four_upper <= FOUR_ACTIVE_UPPER, "four-active cap exceeds target")

    physical = payloads["physical_lower"]
    require(
        physical.get("schema") == "carmenq.four-effect-rational-lower.v1",
        "wrong physical-lower schema",
    )
    require(
        as_fraction(physical["support_weight"]) == SUPPORT_WEIGHT,
        "wrong physical-lower support weight",
    )
    require(not bool(physical["optimiser_called"]), "lower replay called an optimizer")
    require(
        not bool(physical["floating_point_used"]),
        "lower replay used floating-point arithmetic",
    )
    certified_lower = as_fraction(physical["support_lower_fraction"])
    lower = as_fraction(physical["declared_lower_fraction"])
    require(
        certified_lower >= lower == DECLARED_PHYSICAL_LOWER,
        "physical lower witness misses the declared endpoint",
    )
    sector_bounds = {
        "binary_projective": PROJECTIVE_060_UPPER,
        "maximum_effect_at_most_0p88325": low_upper,
        "ternary": TERNARY_UPPER,
        "four_active_minimum_effect_at_least_0p0003": four_upper,
        "four_active_minimum_effect_below_0p0003_by_deletion": DELETION_UPPER,
    }
    assembled = max(sector_bounds.values())
    require(assembled == DELETION_UPPER, "unexpected dominant global sector")
    require(DELETION_UPPER == Fraction(76670, 100000), "deletion arithmetic changed")
    require(assembled == REPORTED_UPPER, "reported endpoint changed")
    require(lower <= assembled, "physical lower exceeds the certified upper")

    manifest = {
        "schema": "carmenq.global-frontier-l060-exact-assembly.v1",
        "support_weight": fraction_pair(SUPPORT_WEIGHT),
        "named_projective_premises": EXPECTED_PROJECTIVE_PREMISES,
        "explicit_physical_lower_fraction": fraction_pair(lower),
        "explicit_physical_lower_decimal": decimal(lower),
        "witness_support_lower_fraction": fraction_pair(certified_lower),
        "witness_support_lower_decimal": decimal(certified_lower),
        "assembled_upper_fraction": fraction_pair(assembled),
        "assembled_upper_decimal": decimal(assembled),
        "reported_outward_upper_fraction": fraction_pair(REPORTED_UPPER),
        "reported_outward_upper_decimal": decimal(REPORTED_UPPER),
        "reported_interval_width_decimal": decimal(REPORTED_UPPER - lower),
        "sector_bounds": {
            name: {
                "fraction": fraction_pair(value),
                "decimal": decimal(value),
            }
            for name, value in sector_bounds.items()
        },
        "auxiliary_projective_line": projective_055,
        "projective_support_line": projective_060,
        "coverage": [
            "binary-projective extreme terminal POVMs",
            "three- and four-active POVMs with maximum effect at most 0.88325",
            "three-active POVMs above the maximum-effect threshold",
            "four-active POVMs above the threshold with minimum effect at least 0.0003",
            "four-active POVMs with a smaller effect, deleted into the ternary sector",
        ],
        "component_replays": {
            "low_weight_cells": int(low["verified_cells"]),
            "ternary_cells": int(ternary["verified_cell_count"]),
            "four_active_leaves": int(four["leaf_count"]),
            "four_active_maximum": str(four_upper),
            "optimizers_called_by_replays": False,
        },
        "artifacts": {
            name: {"path": portable_path(path), "sha256": sha256(path)}
            for name, path in paths.items()
        },
        "certificate_class": (
            "solver-independent finite enclosure: outward interval projective "
            "kernels, exact-dyadic residual dual replay, and the full-closure "
            "dominant-pair projective inequality; component matrix "
            "canonicalisation and the analytic sector reductions remain in the "
            "documented Python/Lean trust boundary"
        ),
        "not_claimed": (
            "exact equality with the explicit four-effect lower strategy, a "
            "closed form for the support curve, or a fully kernel-formalised "
            "end-to-end physical theorem"
        ),
        "complete": True,
    }
    rendered = json.dumps(manifest, indent=2) + "\n"
    print(rendered, end="")
    if args.output is not None:
        output = args.output
        if not output.is_absolute():
            output = ROOT / output
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

"""Assemble and validate the lambda=0.6 interior-frontier certificate.

The unrestricted terminal qubit readout is split by arity and effect traces.

* maximum effect trace at most 0.88325: universal probability cap;
* projective readout: independent finite SCIP cover;
* three active effects above that floor: finite SOCP box cover;
* four active effects with minimum trace at least 0.0003: spatial SCIP
  Helstrom-projection relaxation;
* four active effects below that minimum: merge the smallest outcome into a
  remaining outcome.  AUDIT falls by at most its effect norm, while RETURN is
  unchanged, so the support rises by at most lambda times the threshold over
  the global ternary bound.

This validator checks the complete dyadic ternary partition, the recorded
scope parameters, and the arithmetic of the assembled bound.  The result is
a solver-conditional numerical certificate, not an interval-arithmetic proof
independent of CLARABEL and SCIP.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import tarfile
from typing import Any

from continuous_terminal_projective_envelope_cover import branch_coordinate
from pairwise_inellipse_box_cover import Box, deserialise_box, split_box
from ternary_probability_cone_cover import initial_box, terminal_domain_intersects


ROOT = Path(__file__).resolve().parent
SUPPORT_WEIGHT = 0.6
MAXIMUM_WEIGHT_FLOOR = 0.88325
MINIMUM_ACTIVE_WEIGHT = 0.0003
TERNARY_REPORTED_UPPER = 0.76643
FINAL_REPORTED_UPPER = 0.76662
TOL = 2e-12
ARTIFACTS = {
    "physical_lower": "reduced_four_effect_l060.json",
    "projective_lambda_060": "projective_cover_l060_summary.json",
    "projective_lambda_055": "projective_l055_cover_certificate.json",
    "low_maximum_weight_cap": "terminal_weight_cap_0p88325_l060.json",
    "ternary_cover": "continuous_terminal_projective_l055cert_complete.json.tar.gz",
    "four_active_spatial_bound": "four_active_projection_geometry_scip_60s.json",
}


def load(filename: str) -> dict[str, Any]:
    path = ROOT / filename
    if filename.endswith(".tar.gz"):
        with tarfile.open(path, "r:gz") as archive:
            members = [member for member in archive.getmembers() if member.isfile()]
            if len(members) != 1:
                raise RuntimeError(f"expected one JSON member in {path}")
            stream = archive.extractfile(members[0])
            if stream is None:
                raise RuntimeError(f"could not read {members[0].name} from {path}")
            return json.loads(stream.read().decode("utf-8"))
    return json.loads(path.read_text(encoding="utf-8"))


def canonical_json_sha256(filename: str) -> str:
    """Hash parsed JSON independently of indentation and line endings."""

    payload = load(filename)
    canonical = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(canonical).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def close(first: float, second: float, tolerance: float = TOL) -> bool:
    return math.isclose(first, second, rel_tol=0.0, abs_tol=tolerance)


def box_key(box: Box) -> tuple[tuple[str, float, float], ...]:
    return tuple(
        sorted((name, float(bounds[0]), float(bounds[1])) for name, bounds in box.items())
    )


def validate_ternary_partition(payload: dict[str, Any]) -> dict[str, Any]:
    require(bool(payload["complete"]), "ternary cover is not complete")
    require(not payload["open_nodes"], "complete ternary cover has open nodes")
    require(close(float(payload["support_weight"]), SUPPORT_WEIGHT), "wrong lambda")
    require(
        close(float(payload["maximum_weight_floor"]), MAXIMUM_WEIGHT_FLOOR),
        "wrong ternary maximum-weight floor",
    )
    require(
        close(float(payload["projective_support_upper"]), 0.76591),
        "wrong lambda=0.6 projective line",
    )
    require(
        tuple(tuple(map(float, row)) for row in payload["projective_support_lines"])
        == ((0.55, 0.7573),),
        "wrong auxiliary projective line",
    )
    require(
        close(float(payload["target"]), TERNARY_REPORTED_UPPER),
        "wrong ternary target",
    )
    root = initial_box(0, MAXIMUM_WEIGHT_FLOOR, include_priors=False)
    alpha_weight = float(payload["alpha_branch_weight"])
    leaves = payload["leaves"]
    leaf_map: dict[tuple[tuple[str, float, float], ...], dict[str, Any]] = {}
    for leaf in leaves:
        key = box_key(deserialise_box(leaf["box"]))
        require(key not in leaf_map, f"duplicate ternary leaf {key}")
        leaf_map[key] = leaf

    used: set[tuple[tuple[str, float, float], ...]] = set()

    def visit(box: Box, depth: int = 0) -> None:
        require(depth <= 64, "ternary cover tree is unexpectedly deep")
        key = box_key(box)
        if key in leaf_map:
            leaf = leaf_map[key]
            used.add(key)
            if leaf["status"] == "domain_empty":
                require(
                    not terminal_domain_intersects(box, MAXIMUM_WEIGHT_FLOOR),
                    "a domain-empty leaf intersects the physical strip",
                )
            else:
                require(
                    float(leaf["bound"]) <= TERNARY_REPORTED_UPPER + TOL,
                    "ternary leaf exceeds its target",
                )
            return
        coordinate = branch_coordinate(box, root, alpha_weight)
        for child in split_box(box, coordinate):
            visit(child, depth + 1)

    visit(root)
    require(len(used) == len(leaf_map), "ternary artifact contains extra leaves")
    computed_maximum = max(
        float(leaf["bound"])
        for leaf in leaves
        if leaf["status"] != "domain_empty"
    )
    require(
        close(computed_maximum, float(payload["maximum_leaf_bound"]), 5e-13),
        "ternary maximum leaf bound is inconsistent",
    )
    return {
        "leaf_count": len(leaves),
        "maximum_leaf_bound": computed_maximum,
        "target": float(payload["target"]),
        "solved_nodes": int(payload["solved_nodes"]),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    physical = load("reduced_four_effect_l060.json")
    lower = float(physical["score"])

    projective_060 = load("projective_cover_l060_summary.json")
    require(
        bool(projective_060["all_leaf_duals_at_most_one"]),
        "lambda=0.6 projective cover failed",
    )
    require(
        close(float(projective_060["certified_projective_upper"]), 0.76591),
        "wrong projective upper at lambda=0.6",
    )
    projective_055 = load("projective_l055_cover_certificate.json")
    require(
        bool(projective_055["all_leaf_duals_at_most_one_after_rescaling"]),
        "lambda=0.55 projective cover failed",
    )
    require(
        close(float(projective_055["certified_projective_upper"]), 0.7573),
        "wrong projective upper at lambda=0.55",
    )

    low_weight = load("terminal_weight_cap_0p88325_l060.json")
    low_weight_upper = float(low_weight["bound"])
    require(
        low_weight_upper <= TERNARY_REPORTED_UPPER,
        "low-maximum-weight cap exceeds ternary target",
    )

    ternary = load("continuous_terminal_projective_l055cert_complete.json.tar.gz")
    ternary_summary = validate_ternary_partition(ternary)
    global_ternary_upper = max(low_weight_upper, TERNARY_REPORTED_UPPER)

    active = load("four_active_projection_geometry_scip_60s.json")
    require(
        close(float(active["support_weight"]), SUPPORT_WEIGHT),
        "wrong four-active support weight",
    )
    require(
        close(float(active["maximum_weight_floor"]), MAXIMUM_WEIGHT_FLOOR),
        "wrong four-active maximum-weight floor",
    )
    require(
        close(float(active["minimum_active_weight"]), MINIMUM_ACTIVE_WEIGHT),
        "wrong four-active minimum effect trace",
    )
    four_active_upper = float(active["dual_bound"])

    deletion_upper = global_ternary_upper + SUPPORT_WEIGHT * MINIMUM_ACTIVE_WEIGHT
    assembled_upper = max(
        low_weight_upper,
        global_ternary_upper,
        four_active_upper,
        deletion_upper,
    )
    require(
        assembled_upper <= FINAL_REPORTED_UPPER,
        "assembled upper exceeds reported outward decimal",
    )
    require(lower <= assembled_upper, "physical lower exceeds assembled upper")

    manifest = {
        "support_weight": SUPPORT_WEIGHT,
        "explicit_physical_lower": lower,
        "assembled_solver_conditional_upper": assembled_upper,
        "reported_outward_decimal_upper": FINAL_REPORTED_UPPER,
        "reported_interval_width": FINAL_REPORTED_UPPER - lower,
        "relative_reported_interval_width": (
            FINAL_REPORTED_UPPER - lower
        ) / lower,
        "sector_bounds": {
            "projective": float(projective_060["certified_projective_upper"]),
            "maximum_effect_at_most_0p88325": low_weight_upper,
            "ternary": global_ternary_upper,
            "four_active_minimum_effect_at_least_0p0003": four_active_upper,
            "four_active_minimum_effect_below_0p0003_by_deletion": deletion_upper,
        },
        "ternary_cover": ternary_summary,
        "auxiliary_projective_line": {
            "weight": 0.55,
            "upper": float(projective_055["certified_projective_upper"]),
            "leaf_count": int(projective_055["leaf_count"]),
        },
        "artifacts": {
            role: {
                "path": filename,
                "canonical_json_sha256": canonical_json_sha256(filename),
            }
            for role, filename in ARTIFACTS.items()
        },
        "solver_environment": {
            "numpy": "2.2.6",
            "scipy": "1.15.3",
            "cvxpy": "1.7.5",
            "clarabel": "0.11.1",
            "pyscipopt": "6.2.1",
            "scip": "10.0.2",
        },
        "logical_exhaustion": [
            "projective terminal readout",
            "three-active terminal readout",
            "four-active readout with small effect deleted",
            "four-active readout with all effects above threshold",
        ],
        "certificate_class": (
            "complete finite solver-conditional numerical enclosure at "
            "lambda=0.6; analytic arity reduction and Helstrom identities, "
            "CLARABEL box cover, and spatial SCIP duals"
        ),
        "not_claimed": (
            "formal interval arithmetic, an exact closed-form optimum, or "
            "equality of the explicit 4E protocol with the unrestricted maximum"
        ),
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

"""Validate the first common-instrument strengthening experiment at lambda=.55."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


FILENAMES = {
    "single_scale": "common_instrument_l055_w092_064_044_t1.json",
    "scale_grid": "common_instrument_l055_w092_064_044_tgrid.json",
    "audit": "common_instrument_l055_w092_064_044_audit.json",
    "local_witness": "common_instrument_l055_w092_064_044_local_r010.json",
}


def load(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def validate(directory: Path) -> dict[str, object]:
    payloads = {
        name: load(directory / filename) for name, filename in FILENAMES.items()
    }
    single = payloads["single_scale"]
    grid = payloads["scale_grid"]
    audit = payloads["audit"]
    local = payloads["local_witness"]

    if float(single["weight"]) != 0.55 or float(grid["weight"]) != 0.55:
        raise RuntimeError("common-instrument benchmark uses the wrong support weight")
    expected_weights = [0.92, 0.64, 0.44, 0.0]
    if single["terminal_effect_weights"] != expected_weights:
        raise RuntimeError("single-scale benchmark uses the wrong fixed POVM")
    if grid["terminal_effect_weights"] != expected_weights:
        raise RuntimeError("scale-grid benchmark uses the wrong fixed POVM")
    single_bound = float(single["bound"])
    grid_bound = float(grid["bound"])
    if grid_bound > single_bound + 2e-6:
        raise RuntimeError("the nested scale-grid relaxation unexpectedly increased")
    if abs(float(audit["source_bound"]) - grid_bound) > 1e-10:
        raise RuntimeError("the Choi audit does not reference the scale-grid point")
    if int(audit["negative_flagged_cut_count"]) != int(
        audit["total_flagged_cut_count"]
    ):
        raise RuntimeError("not every recorded flagged comparison is violated")
    projection = audit["choi_projection"]
    if float(projection["distance"]) <= 0.2:
        raise RuntimeError("the exact Choi projection obstruction is too small")
    if float(projection["separation_gap"]) <= 0.2:
        raise RuntimeError("the exact Choi witness failed to separate the point")
    if bool(projection["compatible_at_1e-7"]):
        raise RuntimeError("the audited first-moment point was marked compatible")
    radius_budget = float(projection["uniform_input_trace_radius_budget"])
    if radius_budget <= 0.16:
        raise RuntimeError("the robust witness radius is unexpectedly small")

    witness = local["common_instrument_witness_cut"]
    witness_slack = float(witness["slack"])
    if witness_slack < -2e-7:
        raise RuntimeError("the local witness constraint is numerically violated")
    local_bound = float(local["bound"])
    if local_bound >= 0.735:
        raise RuntimeError("the branch-local witness did not remove the fake basin")
    return {
        "support_weight": 0.55,
        "fixed_terminal_effect_weights": expected_weights,
        "single_scale_bound": single_bound,
        "five_scale_bound": grid_bound,
        "five_scale_improvement": single_bound - grid_bound,
        "worst_first_moment_flagged_violation": float(
            audit["worst_flagged_cut"]["violation"]
        ),
        "exact_choi_projection_distance": float(projection["distance"]),
        "exact_choi_separation_gap": float(projection["separation_gap"]),
        "uniform_input_trace_radius_budget": radius_budget,
        "branch_radius": 0.1,
        "branch_local_bound": local_bound,
        "branch_local_drop": grid_bound - local_bound,
        "witness_constraint_slack": witness_slack,
        "status": (
            "validated solver-conditional local strengthening; the complete "
            "lambda=.55 frontier still requires a covering branch tree"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--directory", type=Path, default=Path(__file__).resolve().parent
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = validate(args.directory)
    rendered = json.dumps(payload, indent=2) + "\n"
    print(rendered, end="")
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

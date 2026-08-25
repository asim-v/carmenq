"""Cross-check the adaptive-frontier kill criterion and exact POVM audit."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from common_effective_povm_audit import audit_common_effective_povm


ROOT = Path(__file__).resolve().parent


def _load(name: str) -> dict[str, Any]:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def validate() -> dict[str, Any]:
    shared = _load(
        "ternary_reconstructed_shared_separator_top_leaf_bbb_p1_s92_l055.json"
    )
    third = _load(
        "ternary_reconstructed_third_separator_g2_top_leaf_bbb_p1_s92_l055.json"
    )
    refined = _load(
        "ternary_reconstructed_third_separator_refined_g4_top_leaf_bbb_p1_s92_l055.json"
    )
    depth4 = _load(
        "ternary_reconstructed_depth4_g2_top_leaf_bbb_p1_s92_l055.json"
    )
    archived_audit = _load("common_effective_povm_audit_depth4_top_l055.json")
    fixed_slice = _load("ternary_fixed_common_povm_depth4_top_l055.json")
    neighbourhood = _load("ternary_common_povm_neighborhood_bisection_l055.json")

    expected = (
        (shared, "crossed_solve_count", 2156, "open_crossed_cells", 246),
        (third, "crossed_solve_count", 6396, "open_crossed_cells", 826),
        (refined, "refined_solve_count", 2719, "open_refined_cells", 815),
        (depth4, "crossed_solve_count", 21190, "open_crossed_cells", 2216),
    )
    for payload, solve_key, solve_count, open_key, open_count in expected:
        if int(payload[solve_key]) != solve_count or int(payload[open_key]) != open_count:
            raise AssertionError("frontier coverage count mismatch")
        if not payload["statuses_complete"]:
            raise AssertionError("a frontier contains an unresolved solver status")
        if payload["complete"]:
            raise AssertionError("an open frontier is incorrectly marked complete")

    if not np.isclose(
        float(refined["maximum_refined_bound"]),
        0.763514590302948,
        atol=2e-10,
    ):
        raise AssertionError("unexpected refined depth-three maximum")
    if not np.isclose(
        float(depth4["maximum_crossed_bound"]),
        float(refined["maximum_refined_bound"]),
        atol=3e-10,
    ):
        raise AssertionError("depth four unexpectedly changed the frontier maximum")

    solution = depth4["top_solution"]
    recomputed = audit_common_effective_povm(
        np.asarray(solution["prefix"], dtype=float),
        np.asarray(solution["input_bloch_vectors"], dtype=float),
        np.asarray(solution["statistics"], dtype=float),
    )
    for key in (
        "determinant",
        "condition_number",
        "minimum_margin",
        "completeness_residual",
    ):
        if not np.isclose(
            float(recomputed[key]), float(archived_audit[key]), rtol=2e-12, atol=2e-14
        ):
            raise AssertionError(f"effective-POVM audit mismatch for {key}")
    if int(recomputed["negative_effect_count"]) != 10:
        raise AssertionError("unexpected negative effective-effect count")
    if recomputed["common_effective_povm"]:
        raise AssertionError("the depth-four maximiser passed the common POVM audit")

    if not fixed_slice["closed_at_target"] or not np.isclose(
        float(fixed_slice["common_povm_bound"]),
        0.7202822292251839,
        atol=2e-9,
    ):
        raise AssertionError("fixed-input common-POVM slice did not reproduce")
    rows = {float(row["row_l1_radius"]): row for row in neighbourhood["rows"]}
    if not rows[0.0871]["closed_at_target"]:
        raise AssertionError("the retained common-POVM neighbourhood did not close")
    if rows[0.0872]["closed_at_target"]:
        raise AssertionError("the first open neighbourhood radius unexpectedly closed")
    if not np.isclose(float(rows[0.0871]["bound"]), 0.7579750191382169, atol=2e-9):
        raise AssertionError("unexpected retained neighbourhood bound")

    return {
        "logical_status": "adaptive separator frontier rejected by kill criterion",
        "source_open_cells": int(refined["open_refined_cells"]),
        "depth4_open_cells": int(depth4["open_crossed_cells"]),
        "depth3_bound": float(refined["maximum_refined_bound"]),
        "depth4_bound": float(depth4["maximum_crossed_bound"]),
        "effective_povm_minimum_margin": float(recomputed["minimum_margin"]),
        "negative_effect_count": int(recomputed["negative_effect_count"]),
        "fixed_input_common_povm_bound": float(fixed_slice["common_povm_bound"]),
        "certified_row_l1_radius": 0.0871,
        "certified_neighbourhood_bound": float(rows[0.0871]["bound"]),
        "next_method": "cover the relevant input-basis region by robust common-POVM boxes",
    }


if __name__ == "__main__":
    print(json.dumps(validate(), indent=2))

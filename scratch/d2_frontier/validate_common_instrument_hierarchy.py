"""Validate the common-instrument hierarchy experiment at lambda = 0.55."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


FILENAMES = {
    "first_level": "choi_moment_ppt_l055.json",
    "behaviour_core": "behavior_core_092_064_044_l055.json",
    "selected_bridge": "clique_bridge_cross_core036_l055_eps1e4.json",
    "sparse_order_two": (
        "common_instrument_sparse_o2_instrument_bridge_eps1e5_l055.json"
    ),
    "sparse_audit": (
        "common_instrument_sparse_o2_instrument_bridge_eps1e5_l055_audit.json"
    ),
    "cell_tree": "common_instrument_cell_tree_l055_pilot1.json",
    "behaviour_disjunction": (
        "behavior_disjunction_092_064_044_l055_4w_300s.json"
    ),
    "exact_shared_instrument": (
        "common_instrument_exact_scip_0123_validation.json"
    ),
}


def load(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def validate(directory: Path) -> dict[str, object]:
    payload = {
        name: load(directory / filename) for name, filename in FILENAMES.items()
    }
    first = payload["first_level"]
    core = payload["behaviour_core"]
    bridge = payload["selected_bridge"]
    sparse = payload["sparse_order_two"]
    audit = payload["sparse_audit"]
    tree = payload["cell_tree"]
    disjunction = payload["behaviour_disjunction"]
    exact = payload["exact_shared_instrument"]

    expected_weights = [0.92, 0.64, 0.44, 0.0]
    for name, item in (("first level", first), ("sparse order two", sparse)):
        if float(item["weight"]) != 0.55:
            raise RuntimeError(f"{name} uses the wrong support weight")
        if item["terminal_effect_weights"] != expected_weights:
            raise RuntimeError(f"{name} uses the wrong terminal POVM")
        if item["prefix_order"] != [0, 1, 2, 3]:
            raise RuntimeError(f"{name} uses the wrong prefix order")

    if not bool(core["minimum_cardinality_certified"]):
        raise RuntimeError("the archived behaviour obstruction is not minimal")
    if int(core["support_size"]) != 3:
        raise RuntimeError("the minimum behaviour obstruction must use three columns")
    if bridge["linked_columns"] != ["b_0", "b_3", "b_6"]:
        raise RuntimeError("the selected-column bridge uses the wrong core")
    if bridge["bridge_basis_mode"] != "cross":
        raise RuntimeError("the selected-column bridge is not cross lifted")

    if sparse["bridge_basis_mode"] != "minimal":
        raise RuntimeError("unexpected sparse Choi bridge basis")
    if sparse["localisers"] != "common" or int(sparse["ideal_degree"]) != 2:
        raise RuntimeError("the sparse Choi model is missing common localisers")
    if int(sparse["instrument_bridge_moment_size"]) != 169:
        raise RuntimeError("the global trace-preservation bridge is missing")
    if float(sparse["moment_min_eigenvalue"]) < -5e-8:
        raise RuntimeError("sparse moment PSD residual is too large")
    if float(sparse["localiser_min_eigenvalue"]) < -2e-7:
        raise RuntimeError("sparse localiser PSD residual is too large")
    if abs(float(sparse["helstrom_independent_residual"])) > 1e-8:
        raise RuntimeError("facial Helstrom audit failed")

    projection = audit["choi_projection"]
    if float(audit["worst_flagged_cut"]["violation"]) <= 0.54:
        raise RuntimeError("the sparse pseudo-moment obstruction disappeared")
    if float(projection["distance"]) <= 0.19:
        raise RuntimeError("the exact fixed-input Choi distance is too small")
    if float(projection["separation_gap"]) <= 0.19:
        raise RuntimeError("the exact Choi witness did not separate the point")
    if bool(projection["compatible_at_1e-7"]):
        raise RuntimeError("the sparse first-moment family was marked compatible")
    if tree["data_processing_mode"] != "cell":
        raise RuntimeError("the branch pilot did not use deterministic cell cuts")
    if int(tree["solved_nodes"]) != 12 or int(tree["closed_nodes"]) != 0:
        raise RuntimeError("unexpected state-cell pilot population")
    if int(tree["open_nodes"]) != 25:
        raise RuntimeError("the state-cell pilot cover has changed")
    if float(disjunction["weight"]) != 0.55:
        raise RuntimeError("behaviour disjunction uses the wrong support weight")
    if disjunction["prefix_order"] != [0, 1, 2, 3]:
        raise RuntimeError("behaviour disjunction uses the wrong prefix order")
    disjunction_round = disjunction["rounds"][0]
    if int(disjunction_round["witness_count"]) != 4:
        raise RuntimeError("unexpected number of active behaviour witnesses")
    if disjunction_round["status"] != "timelimit":
        raise RuntimeError("canonical mixed-integer run changed status")
    if int(exact["scip_nodes"]) < 1:
        raise RuntimeError("exact shared-instrument run explored no nodes")
    repaired = exact["repaired_physical_strategy"]
    if float(repaired["minimum_state_eigenvalue"]) < 0.0:
        raise RuntimeError("repaired exact state is not positive")
    if float(repaired["minimum_choi_eigenvalue"]) < 0.0:
        raise RuntimeError("repaired exact Choi matrix is not positive")
    if float(repaired["trace_preservation_residual"]) > 1e-10:
        raise RuntimeError("repaired exact instrument is not trace preserving")

    first_bound = float(first["bound"])
    sparse_bound = float(sparse["bound"])
    return {
        "support_weight": 0.55,
        "fixed_terminal_effect_weights": expected_weights,
        "first_level_ppt_bound": first_bound,
        "minimum_behaviour_obstruction_size": int(core["support_size"]),
        "selected_cross_bridge_bound": float(bridge["bound"]),
        "sparse_order_two_bound": sparse_bound,
        "sparse_minus_first_level": sparse_bound - first_bound,
        "sparse_scalar_moments": int(sparse["scalar_moment_count"]),
        "instrument_bridge_size": int(sparse["instrument_bridge_moment_size"]),
        "worst_first_moment_flagged_violation": float(
            audit["worst_flagged_cut"]["violation"]
        ),
        "exact_choi_projection_distance": float(projection["distance"]),
        "exact_choi_separation_gap": float(projection["separation_gap"]),
        "cell_tree_solved_nodes": int(tree["solved_nodes"]),
        "cell_tree_closed_nodes": int(tree["closed_nodes"]),
        "cell_tree_open_nodes": int(tree["open_nodes"]),
        "cell_tree_maximum_open_bound": float(
            tree["maximum_open_inherited_bound"]
        ),
        "behaviour_disjunction_primal_bound": float(
            disjunction_round["primal_bound"]
        ),
        "behaviour_disjunction_dual_bound": float(
            disjunction_round["dual_bound"]
        ),
        "behaviour_disjunction_active_witnesses": int(
            disjunction_round["witness_count"]
        ),
        "exact_shared_instrument_physical_score": float(repaired["score"]),
        "exact_shared_instrument_scip_dual_bound": float(
            exact["scip_dual_bound"]
        ),
        "exact_minus_disjunctive_dual": float(exact["scip_dual_bound"])
        - float(disjunction_round["dual_bound"]),
        "physical_improvement_over_previous_checkpoint": float(
            exact["improvement_over_previous_checkpoint"]
        ),
        "status": (
            "validated one literal shared instrument and two independent "
            "solver-conditional upper formulations; neither closes the fixed "
            "interior benchmark, so stronger spatial envelopes remain necessary"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--directory", type=Path, default=Path(__file__).resolve().parent
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = validate(args.directory)
    rendered = json.dumps(result, indent=2) + "\n"
    print(rendered, end="")
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

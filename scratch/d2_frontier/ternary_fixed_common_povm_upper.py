"""Exact common-effective-POVM upper bound on one fixed input basis slice.

The input Pauli matrix is taken from a captured frontier maximiser.  With that
matrix fixed, the existence of one common twelve-outcome qubit POVM is an
ordinary SOCP: its effects are positive Lorentz-cone vectors, sum to identity,
and reproduce all path statistics.  This is a diagnostic slice, not a
neighbourhood or global certificate.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

from pairwise_inellipse_box_cover import deserialise_box
from ternary_probability_cone_cover import TernaryConeOracle


def solve_slice(
    source: dict[str, Any],
    safety: float = 2e-6,
) -> dict[str, Any]:
    solution = source.get("top_solution")
    if not solution:
        raise ValueError("source frontier lacks a captured top solution")
    fixed_input = np.column_stack(
        [
            np.asarray(solution["prefix"], dtype=float),
            np.asarray(solution["input_bloch_vectors"], dtype=float),
        ]
    )
    box = deserialise_box(source["box"])
    oracle = TernaryConeOracle(
        float(source["support_weight"]),
        tuple(int(value) for value in source["prefix_order"]),
        (),
        (),
        float(source["maximum_weight_floor"]),
        float(source["projective_support_upper"]),
        projective_support_lines=tuple(
            tuple(float(value) for value in line)
            for line in source["projective_support_lines"]
        ),
        fixed_common_povm_input=fixed_input,
    )
    result = oracle.solve(box, safety, capture=True)
    return {
        "support_weight": source["support_weight"],
        "target": source["target"],
        "source_bound": solution["bound"],
        "box": source["box"],
        "fixed_input_pauli_matrix": fixed_input.tolist(),
        "fixed_input_determinant": float(np.linalg.det(fixed_input)),
        "fixed_input_condition_number": float(np.linalg.cond(fixed_input)),
        "solver_status": result["status"],
        "common_povm_bound": result["bound"],
        "closed_at_target": float(result["bound"]) < float(source["target"]),
        "result": result,
        "scope": (
            "exact common effective POVM on the single fixed input basis "
            "taken from the source maximiser; solver-conditional; not a "
            "neighbourhood or global certificate"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--frontier-json", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source = json.loads(args.frontier_json.read_text(encoding="utf-8"))
    payload = solve_slice(source)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                key: payload[key]
                for key in (
                    "solver_status",
                    "source_bound",
                    "common_povm_bound",
                    "closed_at_target",
                    "fixed_input_determinant",
                    "fixed_input_condition_number",
                )
            }
        )
    )


if __name__ == "__main__":
    main()

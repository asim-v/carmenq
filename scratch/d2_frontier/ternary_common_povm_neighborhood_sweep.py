"""Continuous input-neighbourhood pilot for the joint common POVM bound.

For a POVM effect ``a0 I + a.v sigma``, positivity implies ``|a_mu| <= a0``.
If every Pauli coordinate of input row ``R[z]`` lies in a box of radii
``d[z,mu]`` around an anchor, then

    |(R[z]-R0[z]).a_k| <= sum_mu d[z,mu] * a0_k.

The resulting affine uncertainty envelope is an SOCP outer relaxation for all
input bases in the box.  This driver sweeps equal row-wise L1 radii, allocated
uniformly among the four Pauli coordinates.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

from pairwise_inellipse_box_cover import deserialise_box
from ternary_probability_cone_cover import TernaryConeOracle


def sweep(
    source: dict[str, Any],
    row_l1_radii: tuple[float, ...],
    safety: float = 2e-6,
) -> dict[str, Any]:
    solution = source.get("top_solution")
    if not solution:
        raise ValueError("source frontier lacks a captured top solution")
    anchor = np.column_stack(
        [
            np.asarray(solution["prefix"], dtype=float),
            np.asarray(solution["input_bloch_vectors"], dtype=float),
        ]
    )
    box = deserialise_box(source["box"])
    rows: list[dict[str, Any]] = []
    for radius in row_l1_radii:
        if not np.isfinite(radius) or radius < 0.0:
            raise ValueError("row L1 radii must be finite and nonnegative")
        coordinate_radii = np.full((4, 4), float(radius) / 4.0)
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
            common_povm_input_anchor=anchor,
            common_povm_input_radii=coordinate_radii,
        )
        result = oracle.solve(box, safety, capture=True)
        rows.append(
            {
                "row_l1_radius": float(radius),
                "coordinate_radii": coordinate_radii.tolist(),
                "status": result["status"],
                "bound": result["bound"],
                "closed_at_target": float(result["bound"]) < float(source["target"]),
                "audit": result.get("audit"),
                "return": result.get("return"),
                "prefix": result.get("prefix"),
                "effective_povm": result.get("effective_povm"),
            }
        )
        print(json.dumps(rows[-1]), flush=True)
    return {
        "support_weight": source["support_weight"],
        "target": source["target"],
        "source_bound": solution["bound"],
        "box": source["box"],
        "input_anchor": anchor.tolist(),
        "input_anchor_determinant": float(np.linalg.det(anchor)),
        "input_anchor_condition_number": float(np.linalg.cond(anchor)),
        "rows": rows,
        "largest_tested_closed_radius": max(
            (row["row_l1_radius"] for row in rows if row["closed_at_target"]),
            default=None,
        ),
        "scope": (
            "continuous coordinate-box neighbourhood around one input Pauli "
            "basis; exact common POVM positivity with a triangle-inequality "
            "probability envelope; solver-conditional; not a global input cover"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--frontier-json", type=Path, required=True)
    parser.add_argument(
        "--row-l1-radius",
        type=float,
        action="append",
        required=True,
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source = json.loads(args.frontier_json.read_text(encoding="utf-8"))
    payload = sweep(source, tuple(args.row_l1_radius))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "largest_tested_closed_radius": payload[
                    "largest_tested_closed_radius"
                ],
                "row_count": len(payload["rows"]),
            }
        )
    )


if __name__ == "__main__":
    main()

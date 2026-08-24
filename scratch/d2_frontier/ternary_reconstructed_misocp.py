"""MISOCP compression of reconstructed common-instrument contractions.

Instead of materialising the Cartesian product of scalar and angular spectral
branches, this diagnostic uses one-hot binary selectors for each contraction.
Every selector family contains the allowed scalar signs and a proved cube-face
cover of the Bloch sphere.  The model is an outer relaxation because cap
projections upper-bound the active Bloch norm.

Only an optimal SCIP result is an upper certificate.  Time-limited primal
values are diagnostic and must never be reported as upper bounds.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

from fourier_behavior_cap_cover import cube_face_caps
from fourier_behavior_upper import CHARACTERS
from terminal_reconstruction_enclosure import reconstruction_anchor_and_errors
from ternary_multicolumn_cell_cover import load_ranked_leaf
from ternary_probability_cone_cover import TernaryConeOracle


def solve(
    box_source: Path,
    leaf_rank: int,
    face_grid: int,
    extra_coefficients: tuple[tuple[float, float, float, float], ...],
    time_limit: float,
    safety: float = 2e-6,
) -> dict[str, Any]:
    if face_grid < 2:
        raise ValueError("face_grid must be at least two")
    box = load_ranked_leaf(box_source, leaf_rank)
    anchor, errors, reconstruction_audit = reconstruction_anchor_and_errors(
        box["terminal_alpha"], box["terminal_beta"]
    )
    caps = tuple(
        np.append(normal, cosine) for normal, cosine in cube_face_caps(face_grid)
    )
    contractions: list[dict[str, object]] = []
    for index, coefficients in enumerate(CHARACTERS):
        contractions.append(
            {
                "coefficients": coefficients,
                "branch": "spectral-cover",
                "scalar_signs": (1,) if index < 2 else (1, -1),
                "caps": caps,
            }
        )
    for raw in extra_coefficients:
        coefficients = np.asarray(raw, dtype=float)
        if np.linalg.norm(coefficients) <= 1e-14:
            raise ValueError("extra contraction coefficients must be nonzero")
        contractions.append(
            {
                "coefficients": coefficients / np.linalg.norm(coefficients),
                "branch": "spectral-cover",
                "scalar_signs": (1, -1),
                "caps": caps,
            }
        )
    oracle = TernaryConeOracle(
        0.55,
        (0, 1, 2, 3),
        (),
        (),
        0.79,
        0.7573,
        mip_time_limit=time_limit,
        projective_support_lines=((0.6, 0.76591),),
        common_contractions=tuple(contractions),
        terminal_reconstruction=(anchor, errors),
    )
    result = oracle.solve(box, safety, capture=True)
    optimal = result["status"] == "optimal"
    return {
        "support_weight": 0.55,
        "target": 0.758,
        "box": {key: list(value) for key, value in box.items()},
        "face_grid": face_grid,
        "cap_count": len(caps),
        "covering_cosine": float(caps[0][3]),
        "extra_coefficients": [list(item) for item in extra_coefficients],
        "contraction_count": len(contractions),
        "binary_selector_count": sum(
            len(item.get("scalar_signs", ())) + len(caps)
            for item in contractions
        ),
        "terminal_reconstruction": reconstruction_audit,
        "solver_status": result["status"],
        "solver_optimal": optimal,
        "reported_upper": result["bound"] if optimal else None,
        "raw_solver_value": result.get("raw_value"),
        "audit": result.get("audit"),
        "return": result.get("return"),
        "prefix": result.get("prefix"),
        "weights": result.get("weights"),
        "contraction_values": result.get("common_contraction_values"),
        "selectors": result.get("common_contraction_selectors"),
        "scope": (
            "one terminal box; reconstructed Fourier and optional adaptive "
            "contractions; finite spectral MISOCP cover; solver-conditional"
        ),
        "warning": (
            "a non-optimal SCIP primal value is not an upper bound"
            if not optimal
            else None
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--box-json", type=Path, required=True)
    parser.add_argument("--leaf-rank", type=int, default=0)
    parser.add_argument("--face-grid", type=int, default=2)
    parser.add_argument(
        "--extra-coefficients", type=float, nargs=4, action="append", default=[]
    )
    parser.add_argument("--time-limit", type=float, default=300.0)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = solve(
        args.box_json,
        args.leaf_rank,
        args.face_grid,
        tuple(tuple(map(float, item)) for item in args.extra_coefficients),
        args.time_limit,
    )
    rendered = json.dumps(payload, indent=2) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered, encoding="utf-8")
    print(
        json.dumps(
            {
                "status": payload["solver_status"],
                "reported_upper": payload["reported_upper"],
                "raw_solver_value": payload["raw_solver_value"],
                "binary_selector_count": payload["binary_selector_count"],
            }
        )
    )


if __name__ == "__main__":
    main()

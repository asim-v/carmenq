"""Add one adaptive measured contraction to a ternary Fourier cap cell.

This is the local branch-and-cut step used after the three fixed Fourier
contractions leave a high relaxed cell.  The supplied real coefficient vector
is normalised and its scalar-positive, scalar-negative, and Bloch-active
spectral branches are covered exhaustively.  The first three input directions
retain the rotational gauge fixed by the parent Fourier cell.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from fourier_behavior_cap_cover import cube_face_caps, plane_caps
from pairwise_inellipse_box_cover import Box, deserialise_box, serialise_box
from ternary_common_instrument_cover import BranchOracle
from terminal_reconstruction_enclosure import reconstruction_anchor_and_errors


def load_ranked_leaf(path: Path, rank: int) -> Box:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if "box" in payload:
        raw = payload["box"]
    elif "leaves" in payload:
        finite = sorted(
            (
                item
                for item in payload["leaves"]
                if item.get("bound") is not None
                and math.isfinite(float(item["bound"]))
            ),
            key=lambda item: float(item["bound"]),
            reverse=True,
        )
        if not 0 <= rank < len(finite):
            raise ValueError("leaf rank is outside the finite cover leaves")
        raw = finite[rank]["box"]
    else:
        raw = payload
    return deserialise_box(raw)


def cover_cell(
    box: Box,
    coefficients: np.ndarray,
    base_branch: str,
    base_plane: tuple[np.ndarray, float] | None,
    base_sphere: tuple[np.ndarray, float] | None,
    support_weight: float = 0.55,
    maximum_weight_floor: float = 0.79,
    projective_support_upper: float = 0.7573,
    projective_support_lines: tuple[tuple[float, float], ...] = ((0.6, 0.76591),),
    prefix_order: tuple[int, int, int, int] = (0, 1, 2, 3),
    extra_plane_cells: int = 16,
    extra_face_grid: int = 4,
    safety: float = 2e-6,
    capture_top: bool = False,
    use_reconstruction: bool = False,
) -> dict[str, Any]:
    value = np.asarray(coefficients, dtype=float)
    if value.shape != (4,) or np.linalg.norm(value) <= 1e-14:
        raise ValueError("one nonzero four-component coefficient vector is required")
    value = value / np.linalg.norm(value)
    reconstruction = None
    reconstruction_audit = None
    if use_reconstruction:
        anchor, errors, reconstruction_audit = reconstruction_anchor_and_errors(
            box["terminal_alpha"], box["terminal_beta"]
        )
        reconstruction = (anchor, errors)
    planes = plane_caps(extra_plane_cells)
    spheres = cube_face_caps(extra_face_grid)
    rows: list[dict[str, Any]] = []
    top: tuple[
        BranchOracle,
        tuple[np.ndarray, float] | None,
    ] | None = None
    maximum = -math.inf
    for extra_branch in ("scalar-positive", "scalar-negative", "bloch"):
        oracle = BranchOracle(
            base_branch,
            support_weight,
            prefix_order,
            maximum_weight_floor,
            projective_support_upper,
            projective_support_lines,
            extra_contraction={
                "coefficients": value,
                "branch": extra_branch,
            },
            terminal_reconstruction=reconstruction,
        )
        if extra_branch == "bloch" and oracle.extra_cap_kind == "plane":
            options: tuple[tuple[np.ndarray, float] | None, ...] = tuple(planes)
        elif extra_branch == "bloch" and oracle.extra_cap_kind == "sphere":
            options = tuple(spheres)
        else:
            options = (None,)
        for cap_index, cap in enumerate(options):
            result = oracle.solve(box, safety, base_plane, base_sphere, cap)
            bound = float(result["bound"])
            row = {
                "extra_branch": extra_branch,
                "extra_cap": None if cap is None else cap_index,
                "status": result["status"],
                "bound": bound,
                "audit": result.get("audit"),
                "return": result.get("return"),
            }
            rows.append(row)
            if bound > maximum:
                maximum = bound
                top = (oracle, cap)
    top_solution = None
    if capture_top and top is not None and math.isfinite(maximum):
        oracle, cap = top
        top_solution = oracle.solve(
            box, safety, base_plane, base_sphere, cap, capture=True
        )
    return {
        "support_weight": support_weight,
        "maximum_weight_floor": maximum_weight_floor,
        "projective_support_upper": projective_support_upper,
        "projective_support_lines": [list(line) for line in projective_support_lines],
        "prefix_order": list(prefix_order),
        "box": serialise_box(box),
        "base_branch": base_branch,
        "base_plane": None if base_plane is None else {
            "normal": base_plane[0].tolist(), "cosine": base_plane[1]
        },
        "base_sphere": None if base_sphere is None else {
            "normal": base_sphere[0].tolist(), "cosine": base_sphere[1]
        },
        "coefficients": value.tolist(),
        "terminal_reconstruction": reconstruction_audit,
        "extra_plane_cells": extra_plane_cells,
        "extra_face_grid": extra_face_grid,
        "extra_sphere_covering_cosine": spheres[0][1],
        "solve_count": len(rows),
        "maximum_bound": maximum,
        "top_cell": max(rows, key=lambda row: float(row["bound"])),
        "cells": rows,
        "top_solution": top_solution,
        "complete": all(
            row["status"]
            in {"optimal", "optimal_inaccurate", "infeasible", "infeasible_inaccurate"}
            for row in rows
        ),
        "scope": (
            "one terminal box and one parent Fourier angular cell; exhaustive "
            "spectral/angular cover for one additional measured contraction; "
            "solver-conditional"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--box-json", type=Path, required=True)
    parser.add_argument("--leaf-rank", type=int, default=0)
    parser.add_argument("--base-branch", default="pbb")
    parser.add_argument("--base-plane-cells", type=int, default=16)
    parser.add_argument("--base-plane-index", type=int)
    parser.add_argument("--base-face-grid", type=int, default=2)
    parser.add_argument("--base-sphere-index", type=int)
    parser.add_argument("--coefficients", type=float, nargs=4, required=True)
    parser.add_argument("--extra-plane-cells", type=int, default=16)
    parser.add_argument("--extra-face-grid", type=int, default=4)
    parser.add_argument("--capture-top", action="store_true")
    parser.add_argument("--reconstruction", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    planes = plane_caps(args.base_plane_cells)
    spheres = cube_face_caps(args.base_face_grid)
    base_plane = (
        None if args.base_plane_index is None else planes[args.base_plane_index]
    )
    base_sphere = (
        None if args.base_sphere_index is None else spheres[args.base_sphere_index]
    )
    result = cover_cell(
        load_ranked_leaf(args.box_json, args.leaf_rank),
        np.asarray(args.coefficients, dtype=float),
        args.base_branch,
        base_plane,
        base_sphere,
        extra_plane_cells=args.extra_plane_cells,
        extra_face_grid=args.extra_face_grid,
        capture_top=args.capture_top,
        use_reconstruction=args.reconstruction,
    )
    rendered = json.dumps(result, indent=2) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered, encoding="utf-8")
    print(
        json.dumps(
            {
                "complete": result["complete"],
                "solve_count": result["solve_count"],
                "maximum_bound": result["maximum_bound"],
                "top_cell": result["top_cell"],
            }
        )
    )


if __name__ == "__main__":
    main()

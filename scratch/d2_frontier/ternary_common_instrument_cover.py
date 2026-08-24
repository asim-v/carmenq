"""Fourier common-instrument cover for a continuous ternary terminal box.

The terminal probability-cone model knows the complete measured statistics
``q[z,y,t]`` but not a shared output map.  A common quantum instrument forces,
for every real coefficient vector ``c``,

    sum_{y,t} |sum_z c[z] q[z,y,t]|
        <= ||sum_z c[z] rho[z]||_1.

This driver imposes the three nontrivial Z2^2 Fourier contractions.  It covers
the reverse-convex qubit norm exactly by its scalar-positive, scalar-negative,
and Bloch-active branches.  Global rotational symmetry fixes the first active
Bloch vector to +z and the second to the xz half-plane; finite angular caps
cover the remaining directions.

The result is a solver-conditional outer bound for one terminal-parameter box.
It is a diagnostic building block, not by itself a cover of the entire
terminal strip.
"""

from __future__ import annotations

import argparse
from itertools import product
import json
import math
from pathlib import Path
from typing import Any

import cvxpy as cp
import numpy as np

from fourier_behavior_cap_cover import cube_face_caps, plane_caps
from fourier_behavior_upper import CHARACTERS
from pairwise_inellipse_box_cover import Box, deserialise_box, serialise_box
from ternary_probability_cone_cover import TernaryConeOracle, initial_box
from terminal_reconstruction_enclosure import reconstruction_anchor_and_errors


BRANCH_NAME = {
    "p": "scalar-positive",
    "n": "scalar-negative",
    "b": "bloch",
}


def branch_codes() -> tuple[str, ...]:
    """Exhaustive Fourier spectral branches under sorted prefix priors."""

    # The first two Fourier scalars are nonnegative whenever
    # p0 >= p1 >= p2 >= p3.  The third scalar can have either sign.
    return tuple(
        "".join(code)
        for code in product(("p", "b"), ("p", "b"), ("p", "n", "b"))
    )


class BranchOracle:
    """Reusable DPP oracle for one spectral branch code."""

    def __init__(
        self,
        code: str,
        support_weight: float,
        prefix_order: tuple[int, int, int, int],
        maximum_weight_floor: float,
        projective_support_upper: float,
        projective_support_lines: tuple[tuple[float, float], ...],
        extra_contraction: dict[str, object] | None = None,
        terminal_reconstruction: tuple[np.ndarray, np.ndarray] | None = None,
        parameterized_reconstruction: bool = False,
    ) -> None:
        if code not in branch_codes():
            raise ValueError(f"invalid Fourier branch code {code!r}")
        self.code = code
        self.plane_cap: cp.Parameter | None = None
        self.sphere_cap: cp.Parameter | None = None
        self.extra_cap: cp.Parameter | None = None
        self.extra_cap_kind: str | None = None
        contractions: list[dict[str, object]] = []
        bloch_seen = 0
        for coefficients, letter in zip(CHARACTERS, code, strict=True):
            contraction: dict[str, object] = {
                "coefficients": coefficients,
                "branch": BRANCH_NAME[letter],
            }
            if letter == "b":
                contraction["gauge_rank"] = bloch_seen
                if bloch_seen == 1:
                    self.plane_cap = cp.Parameter(3)
                    contraction["cap"] = self.plane_cap
                elif bloch_seen == 2:
                    self.sphere_cap = cp.Parameter(3)
                    contraction["cap"] = self.sphere_cap
                bloch_seen += 1
            contractions.append(contraction)
        if extra_contraction is not None:
            coefficients = np.asarray(
                extra_contraction["coefficients"], dtype=float
            )
            if coefficients.shape != (4,) or np.linalg.norm(coefficients) <= 1e-14:
                raise ValueError("invalid extra contraction coefficients")
            branch = str(extra_contraction["branch"])
            contraction = {
                "coefficients": coefficients / np.linalg.norm(coefficients),
                "branch": branch,
            }
            if branch == "bloch":
                contraction["gauge_rank"] = bloch_seen
                if bloch_seen >= 1:
                    self.extra_cap = cp.Parameter(3)
                    self.extra_cap_kind = "plane" if bloch_seen == 1 else "sphere"
                    contraction["cap"] = self.extra_cap
            elif branch not in {"scalar-positive", "scalar-negative"}:
                raise ValueError("invalid extra contraction branch")
            contractions.append(contraction)
        self.reconstruction_anchor: cp.Parameter | None = None
        self.reconstruction_errors: cp.Parameter | None = None
        if parameterized_reconstruction:
            if terminal_reconstruction is not None:
                raise ValueError("choose fixed or parameterized reconstruction, not both")
            self.reconstruction_anchor = cp.Parameter((2, 3))
            self.reconstruction_errors = cp.Parameter(3, nonneg=True)
            terminal_reconstruction = (
                self.reconstruction_anchor,
                self.reconstruction_errors,
            )
        self.oracle = TernaryConeOracle(
            support_weight,
            prefix_order,
            (),
            (),
            maximum_weight_floor,
            projective_support_upper,
            projective_support_lines=projective_support_lines,
            common_contractions=tuple(contractions),
            terminal_reconstruction=terminal_reconstruction,
        )

    def solve(
        self,
        box: Box,
        safety: float,
        plane: tuple[np.ndarray, float] | None,
        sphere: tuple[np.ndarray, float] | None,
        extra_cap: tuple[np.ndarray, float] | None = None,
        capture: bool = False,
    ) -> dict[str, Any]:
        if (self.plane_cap is None) != (plane is None):
            raise ValueError("plane cap assignment does not match branch")
        if (self.sphere_cap is None) != (sphere is None):
            raise ValueError("sphere cap assignment does not match branch")
        if self.plane_cap is not None and plane is not None:
            self.plane_cap.value = plane[0] / plane[1]
        if self.sphere_cap is not None and sphere is not None:
            self.sphere_cap.value = sphere[0] / sphere[1]
        if (self.extra_cap is None) != (extra_cap is None):
            raise ValueError("extra cap assignment does not match branch")
        if self.extra_cap is not None and extra_cap is not None:
            self.extra_cap.value = extra_cap[0] / extra_cap[1]
        if self.reconstruction_anchor is not None:
            anchor, errors, _ = reconstruction_anchor_and_errors(
                box["terminal_alpha"], box["terminal_beta"]
            )
            self.reconstruction_anchor.value = anchor
            self.reconstruction_errors.value = errors
        return self.oracle.solve(box, safety, capture=capture)


def cover_box(
    box: Box,
    support_weight: float = 0.55,
    maximum_weight_floor: float = 0.79,
    projective_support_upper: float = 0.7573,
    projective_support_lines: tuple[tuple[float, float], ...] = ((0.6, 0.76591),),
    prefix_order: tuple[int, int, int, int] = (0, 1, 2, 3),
    plane_cells: int = 4,
    face_grid: int = 2,
    safety: float = 2e-6,
    capture_top: bool = False,
    selected_branches: tuple[str, ...] | None = None,
    use_reconstruction: bool = False,
) -> dict[str, Any]:
    if face_grid < 2:
        raise ValueError("face_grid must be at least two for a finite cap envelope")
    planes = plane_caps(plane_cells)
    spheres = cube_face_caps(face_grid)
    rows: list[dict[str, Any]] = []
    maximum = -math.inf
    top: tuple[BranchOracle, tuple[np.ndarray, float] | None, tuple[np.ndarray, float] | None] | None = None
    solve_count = 0
    codes = branch_codes() if selected_branches is None else selected_branches
    if not codes or any(code not in branch_codes() for code in codes):
        raise ValueError("selected branches must be nonempty valid branch codes")
    reconstruction = None
    reconstruction_audit = None
    if use_reconstruction:
        anchor, errors, reconstruction_audit = reconstruction_anchor_and_errors(
            box["terminal_alpha"], box["terminal_beta"]
        )
        reconstruction = (anchor, errors)
    for code in codes:
        branch = BranchOracle(
            code,
            support_weight,
            prefix_order,
            maximum_weight_floor,
            projective_support_upper,
            projective_support_lines,
            terminal_reconstruction=reconstruction,
        )
        plane_options: tuple[tuple[np.ndarray, float] | None, ...] = (
            tuple(planes) if branch.plane_cap is not None else (None,)
        )
        sphere_options: tuple[tuple[np.ndarray, float] | None, ...] = (
            tuple(spheres) if branch.sphere_cap is not None else (None,)
        )
        for plane_index, plane in enumerate(plane_options):
            for sphere_index, sphere in enumerate(sphere_options):
                result = branch.solve(box, safety, plane, sphere)
                solve_count += 1
                bound = float(result["bound"])
                row = {
                    "branch": code,
                    "plane_cap": None if plane is None else plane_index,
                    "sphere_cap": None if sphere is None else sphere_index,
                    "status": result["status"],
                    "bound": bound,
                    "audit": result.get("audit"),
                    "return": result.get("return"),
                }
                rows.append(row)
                if bound > maximum:
                    maximum = bound
                    top = (branch, plane, sphere)
    top_solution = None
    if capture_top and top is not None and math.isfinite(maximum):
        branch, plane, sphere = top
        top_solution = branch.solve(box, safety, plane, sphere, capture=True)
    return {
        "support_weight": support_weight,
        "maximum_weight_floor": maximum_weight_floor,
        "projective_support_upper": projective_support_upper,
        "projective_support_lines": [list(line) for line in projective_support_lines],
        "prefix_order": list(prefix_order),
        "box": serialise_box(box),
        "plane_cells": plane_cells,
        "face_grid": face_grid,
        "plane_covering_cosine": planes[0][1],
        "sphere_covering_cosine": spheres[0][1],
        "branches": list(codes),
        "terminal_reconstruction": reconstruction_audit,
        "branch_count": len(codes),
        "solve_count": solve_count,
        "maximum_bound": maximum,
        "top_cell": max(rows, key=lambda row: float(row["bound"])),
        "cells": rows,
        "top_solution": top_solution,
        "complete": all(math.isfinite(float(row["bound"])) for row in rows),
        "scope": (
            "one continuous ternary terminal-parameter box; exhaustive Fourier "
            "trace-norm branches and finite angular cover; solver-conditional"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lambda", dest="support_weight", type=float, default=0.55)
    parser.add_argument("--maximum-weight-floor", type=float, default=0.79)
    parser.add_argument("--projective-support-upper", type=float, default=0.7573)
    parser.add_argument(
        "--projective-line",
        type=float,
        nargs=2,
        action="append",
        default=[(0.6, 0.76591)],
    )
    parser.add_argument("--plane-cells", type=int, default=4)
    parser.add_argument("--face-grid", type=int, default=2)
    parser.add_argument("--safety", type=float, default=2e-6)
    parser.add_argument("--box-json", type=Path)
    parser.add_argument(
        "--leaf-rank",
        type=int,
        default=0,
        help="when --box-json is a cover, select this finite leaf by descending bound",
    )
    parser.add_argument("--capture-top", action="store_true")
    parser.add_argument("--branch", action="append", default=[])
    parser.add_argument("--reconstruction", action="store_true")
    parser.add_argument(
        "--summary-only",
        action="store_true",
        help="print only the maximum while retaining the complete output artifact",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.box_json is None:
        box = initial_box(0, args.maximum_weight_floor, include_priors=False)
    else:
        payload = json.loads(args.box_json.read_text(encoding="utf-8"))
        if "box" in payload:
            raw_box = payload["box"]
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
            if not 0 <= args.leaf_rank < len(finite):
                raise ValueError("leaf rank is outside the finite cover leaves")
            raw_box = finite[args.leaf_rank]["box"]
        else:
            raw_box = payload
        box = deserialise_box(raw_box)
    result = cover_box(
        box,
        args.support_weight,
        args.maximum_weight_floor,
        args.projective_support_upper,
        tuple(tuple(map(float, line)) for line in args.projective_line),
        plane_cells=args.plane_cells,
        face_grid=args.face_grid,
        safety=args.safety,
        capture_top=args.capture_top,
        selected_branches=(tuple(args.branch) if args.branch else None),
        use_reconstruction=args.reconstruction,
    )
    rendered = json.dumps(result, indent=2) + "\n"
    if args.summary_only:
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
    else:
        print(rendered, end="")
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

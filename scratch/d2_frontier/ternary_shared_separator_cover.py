"""Share a second common-instrument separator across an open angular frontier.

The adaptive tree selects a valid real coefficient vector at its worst open
node.  Validity of the associated trace-norm contraction is global: the same
coefficient may therefore be imposed on every sibling orientation left open
by the preceding separator.  This driver performs that product cover without
re-canonicalising one SOCP per angular pair.

The input checkpoint must contain two expansions.  Expansion zero supplies
the first separator and its complete scalar/Bloch cover; expansion one
supplies the second separator.  Every still-open first-generation Bloch cap is
then crossed with the scalar-positive, scalar-negative, and exhaustive Bloch
cover of the second separator.  Together with the first-generation children
already below target, this certifies the original fixed terminal/Fourier cell
whenever every crossed child is also below target.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import cvxpy as cp
import numpy as np

from fourier_behavior_cap_cover import cube_face_caps
from pairwise_inellipse_box_cover import deserialise_box
from terminal_reconstruction_enclosure import reconstruction_anchor_and_errors
from ternary_multicolumn_branch_tree import (
    build_oracle,
    fixed_fourier_contractions,
)


SOLVED_STATUSES = {
    "optimal",
    "optimal_inaccurate",
    "infeasible",
    "infeasible_inaccurate",
}


def _cap(raw: dict[str, Any] | None) -> tuple[np.ndarray, float] | None:
    if raw is None:
        return None
    normal = np.asarray(raw["normal"], dtype=float)
    cosine = float(raw["cosine"])
    if normal.shape != (3,) or cosine <= 0.0:
        raise ValueError("invalid stored angular cap")
    return normal, cosine


def extract_shared_frontier(
    payload: dict[str, Any], target: float
) -> tuple[np.ndarray, np.ndarray, tuple[dict[str, Any], ...], dict[str, int]]:
    """Extract two separators and all open children of the first expansion."""

    nodes = payload.get("nodes", [])
    if len(nodes) < 2:
        raise ValueError("the checkpoint needs at least two adaptive expansions")
    first, second = nodes[0], nodes[1]
    # Older checkpoints let the solver status overwrite the descriptive
    # ``expanded`` label.  Presence of the separator/children payload is the
    # stable structural test.
    if "separator" not in first or "children" not in first:
        raise ValueError("the first checkpoint node must be an expansion")
    if "separator" not in second or "children" not in second:
        raise ValueError("the second checkpoint node must be an expansion")
    first_coefficients = np.asarray(
        first["separator"]["coefficients"], dtype=float
    )
    second_coefficients = np.asarray(
        second["separator"]["coefficients"], dtype=float
    )
    if first_coefficients.shape != (4,) or second_coefficients.shape != (4,):
        raise ValueError("stored separators must have four coefficients")

    children = first.get("children", [])
    open_children = [
        child for child in children if float(child.get("bound", math.inf)) >= target
    ]
    if any(
        child.get("branch")
        not in {"scalar-positive", "scalar-negative", "bloch"}
        for child in open_children
    ):
        raise ValueError("the first expansion contains an unknown spectral branch")
    bloch_indices = tuple(
        int(child["cap"])
        for child in open_children
        if child["branch"] == "bloch"
    )
    if len(set(bloch_indices)) != len(bloch_indices):
        raise ValueError("the first-generation cap list contains duplicates")
    parents = tuple(
        {
            "branch": str(child["branch"]),
            "cap": None if child["cap"] is None else int(child["cap"]),
        }
        for child in open_children
    )
    counts = {
        "first_generation_total": len(children),
        "first_generation_closed": len(children) - len(open_children),
        "first_generation_open": len(open_children),
    }
    return first_coefficients, second_coefficients, parents, counts


def cover_shared_separator(
    payload: dict[str, Any],
    target: float = 0.758,
    safety: float | None = None,
    capture_top: bool = False,
) -> dict[str, Any]:
    first, second, parents, counts = extract_shared_frontier(payload, target)
    grid = int(payload["contraction_grid"])
    caps = cube_face_caps(grid)
    if any(
        parent["cap"] is not None
        and not 0 <= int(parent["cap"]) < len(caps)
        for parent in parents
    ):
        raise ValueError("stored first-generation cap index is outside its cover")
    safety_value = float(payload["safety"] if safety is None else safety)
    box = deserialise_box(payload["box"])
    base_plane = _cap(payload.get("base_plane"))
    base_sphere = _cap(payload.get("base_sphere"))
    base = fixed_fourier_contractions(
        str(payload["base_code"]), base_plane, base_sphere
    )
    reconstruction_anchor, reconstruction_errors, reconstruction_audit = (
        reconstruction_anchor_and_errors(
            box["terminal_alpha"], box["terminal_beta"]
        )
    )
    reconstruction = (reconstruction_anchor, reconstruction_errors)
    support_weight = float(payload["support_weight"])
    prefix_order = tuple(int(value) for value in payload["prefix_order"])
    maximum_weight_floor = float(payload["maximum_weight_floor"])
    projective_support_upper = float(payload["projective_support_upper"])
    projective_support_lines = tuple(
        tuple(float(value) for value in line)
        for line in payload["projective_support_lines"]
    )
    base_bloch_count = sum(
        str(item["branch"]) == "bloch" for item in base
    )

    # Only three DPP models per surviving parent branch type are canonicalised.
    # Their cap parameters are then updated across the angular product cover.
    oracles: dict[
        tuple[str, str], tuple[Any, cp.Parameter | None, cp.Parameter | None]
    ] = {}
    parent_branches = tuple(dict.fromkeys(str(parent["branch"]) for parent in parents))
    for parent_branch in parent_branches:
        for branch in ("scalar-positive", "scalar-negative", "bloch"):
            first_cap: cp.Parameter | None = None
            inherited: dict[str, object] = {
                "coefficients": first,
                "branch": parent_branch,
            }
            inherited_bloch_count = base_bloch_count
            if parent_branch == "bloch":
                first_cap = cp.Parameter(3)
                inherited.update(
                    {
                        "gauge_rank": base_bloch_count,
                        "cap": first_cap,
                    }
                )
                inherited_bloch_count += 1
            second_cap: cp.Parameter | None = None
            added: dict[str, object] = {
                "coefficients": second,
                "branch": branch,
            }
            if branch == "bloch":
                second_cap = cp.Parameter(3)
                added.update(
                    {
                        "gauge_rank": inherited_bloch_count,
                        "cap": second_cap,
                    }
                )
            oracle = build_oracle(
                support_weight,
                prefix_order,  # type: ignore[arg-type]
                maximum_weight_floor,
                projective_support_upper,
                projective_support_lines,  # type: ignore[arg-type]
                base + (inherited, added),
                reconstruction,
            )
            oracles[(parent_branch, branch)] = (oracle, first_cap, second_cap)

    rows: list[dict[str, Any]] = []
    maximum = -math.inf
    top_key: tuple[str, int | None, str, int | None] | None = None
    for parent_number, parent in enumerate(parents, start=1):
        parent_branch = str(parent["branch"])
        parent_index = parent["cap"]
        for branch in ("scalar-positive", "scalar-negative", "bloch"):
            _, first_cap, second_cap = oracles[(parent_branch, branch)]
            if first_cap is not None and parent_index is not None:
                parent_normal, parent_cosine = caps[int(parent_index)]
                first_cap.value = parent_normal / parent_cosine
            child_options: tuple[tuple[int | None, tuple[np.ndarray, float] | None], ...]
            if branch == "bloch":
                child_options = tuple(enumerate(caps))
            else:
                child_options = ((None, None),)
            for child_index, child_cap in child_options:
                if second_cap is not None and child_cap is not None:
                    second_cap.value = child_cap[0] / child_cap[1]
                result = oracles[(parent_branch, branch)][0].solve(box, safety_value)
                bound = float(result["bound"])
                row = {
                    "first_branch": parent_branch,
                    "first_cap": parent_index,
                    "second_branch": branch,
                    "second_cap": child_index,
                    "status": result["status"],
                    "bound": bound,
                    "audit": result.get("audit"),
                    "return": result.get("return"),
                }
                rows.append(row)
                if bound > maximum:
                    maximum = bound
                    top_key = (parent_branch, parent_index, branch, child_index)
        print(
            json.dumps(
                {
                    "parents_done": parent_number,
                    "parents_total": len(parents),
                    "solve_count": len(rows),
                    "maximum_bound": maximum,
                }
            ),
            flush=True,
        )

    top_solution = None
    if capture_top and top_key is not None and math.isfinite(maximum):
        parent_branch, parent_index, branch, child_index = top_key
        oracle, first_cap, second_cap = oracles[(parent_branch, branch)]
        if first_cap is not None and parent_index is not None:
            parent_normal, parent_cosine = caps[int(parent_index)]
            first_cap.value = parent_normal / parent_cosine
        if second_cap is not None and child_index is not None:
            normal, cosine = caps[child_index]
            second_cap.value = normal / cosine
        top_solution = oracle.solve(box, safety_value, capture=True)

    open_rows = [row for row in rows if float(row["bound"]) >= target]
    statuses_complete = all(row["status"] in SOLVED_STATUSES for row in rows)
    return {
        "support_weight": support_weight,
        "maximum_weight_floor": maximum_weight_floor,
        "projective_support_upper": projective_support_upper,
        "projective_support_lines": [list(line) for line in projective_support_lines],
        "prefix_order": list(prefix_order),
        "target": target,
        "safety": safety_value,
        "base_code": payload["base_code"],
        "box": payload["box"],
        "base_plane": payload.get("base_plane"),
        "base_sphere": payload.get("base_sphere"),
        "terminal_reconstruction": reconstruction_audit,
        "contraction_grid": grid,
        "cap_count": len(caps),
        "covering_cosine": caps[0][1],
        "first_separator": first.tolist(),
        "shared_second_separator": second.tolist(),
        **counts,
        "crossed_parents": list(parents),
        "crossed_solve_count": len(rows),
        "maximum_crossed_bound": maximum,
        "top_cell": None
        if not rows
        else max(rows, key=lambda row: float(row["bound"])),
        "open_crossed_cells": len(open_rows),
        "statuses_complete": statuses_complete,
        "complete": statuses_complete and not open_rows,
        "cells": rows,
        "top_solution": top_solution,
        "scope": (
            "one terminal box and one fixed Fourier angular cell; the second "
            "adaptive common-instrument contraction is crossed with every "
            "open first-generation angular branch; solver-conditional"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree-json", type=Path, required=True)
    parser.add_argument("--target", type=float, default=0.758)
    parser.add_argument("--capture-top", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = json.loads(args.tree_json.read_text(encoding="utf-8"))
    result = cover_shared_separator(payload, args.target, capture_top=args.capture_top)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "complete": result["complete"],
                "crossed_solve_count": result["crossed_solve_count"],
                "open_crossed_cells": result["open_crossed_cells"],
                "maximum_crossed_bound": result["maximum_crossed_bound"],
                "top_cell": result["top_cell"],
            }
        )
    )


if __name__ == "__main__":
    main()

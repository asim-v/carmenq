"""Batch common-instrument product-localizer covers over spectral cells.

The depth-four reconstructed frontier contains many cells that differ only in
the branch type and angular cap of four fixed contractions.  Rebuilding an SDP
for every cell would make the next frontier step needlessly expensive.  This
driver groups cells by branch pattern, represents their scaled cap normals by
CVXPY parameters, and reuses two DPP problems per pattern:

1. a support-function problem enclosing every input basis that can reach the
   target in that spectral cell;
2. a spatial input--Choi product relaxation with PSD cone-RLT sandwiches and
   the lifted common-instrument trace identities.

The raw output retains every spatial tree.  It is a solver-conditional research
checkpoint, not a solver-independent interval certificate.
"""

from __future__ import annotations

import argparse
import collections
import gc
import hashlib
import heapq
import json
import math
from pathlib import Path
import platform
import time
from typing import Any

import cvxpy as cp
import numpy as np

from fourier_behavior_cap_cover import cube_face_caps
from pairwise_inellipse_box_cover import Box, deserialise_box
from terminal_reconstruction_enclosure import (
    reconstruction_anchor_and_errors,
    terminal_effect_anchor_and_errors,
)
from ternary_bilinear_instrument_input_cover import (
    box_purity_caps,
    product_residual_scores,
)
from ternary_common_povm_input_cover import (
    _oracle_keywords,
    _plane,
    _split_coordinate,
)
from ternary_multicolumn_branch_tree import fixed_fourier_contractions
from ternary_probability_cone_cover import TernaryConeOracle


SCHEMA = "carmenq.spectral-product-localizer-batch.v3"
ACCEPTED_STATUSES = {
    cp.OPTIMAL,
    cp.OPTIMAL_INACCURATE,
    cp.INFEASIBLE,
    cp.INFEASIBLE_INACCURATE,
}


def strict_extended_real(value: float) -> tuple[float | None, str]:
    """Encode an extended-real value in standards-compliant JSON."""

    if math.isfinite(value):
        return float(value), "finite"
    if value == math.inf:
        return None, "positive-infinity"
    if value == -math.inf:
        return None, "negative-infinity"
    return None, "not-a-number"


def pattern_code(pattern: tuple[str, ...]) -> str:
    symbols = {
        "bloch": "b",
        "scalar-positive": "+",
        "scalar-negative": "-",
    }
    try:
        return "".join(symbols[item] for item in pattern)
    except KeyError as error:
        raise ValueError(f"unknown spectral branch {error.args[0]!r}") from error


def source_open_cells(source: dict[str, Any]) -> list[dict[str, Any]]:
    """Return source cells not already below target, sorted by source bound."""

    target = float(source["target"])
    cells = []
    for index, cell in enumerate(source["cells"]):
        status = str(cell["status"])
        bound = float(cell["bound"])
        if status not in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE} or bound < target:
            continue
        cells.append(
            {
                "source_index": int(index),
                "source_cell": int(cell["source_cell"]),
                "branches": tuple(str(item) for item in cell["branches"]),
                "caps": tuple(cell["caps"]),
                "source_status": status,
                "source_bound": bound,
                "source_audit": float(cell["audit"]),
                "source_return": float(cell["return"]),
            }
        )
    return sorted(cells, key=lambda cell: (-cell["source_bound"], cell["source_index"]))


def parameterized_contractions(
    source: dict[str, Any],
    pattern: tuple[str, ...],
) -> tuple[tuple[dict[str, object], ...], dict[int, cp.Parameter]]:
    """Build one fixed-pattern contraction family with parameterized caps."""

    coefficients = tuple(
        np.asarray(item, dtype=float) for item in source["separator_coefficients"]
    )
    grids = tuple(int(value) for value in source["separator_grids"])
    if not (len(coefficients) == len(grids) == len(pattern)):
        raise ValueError("separator coefficients, grids, and pattern are not parallel")
    contractions = list(
        fixed_fourier_contractions(
            str(source["base_code"]),
            _plane(source.get("base_plane")),
            _plane(source.get("base_sphere")),
        )
    )
    cap_parameters: dict[int, cp.Parameter] = {}
    bloch_count = sum(item["branch"] == "bloch" for item in contractions)
    for position, (coefficient, branch) in enumerate(
        zip(coefficients, pattern, strict=True)
    ):
        item: dict[str, object] = {
            "coefficients": coefficient,
            "branch": branch,
        }
        if branch == "bloch":
            cap = cp.Parameter(3, name=f"spectral_cap_{position}")
            item.update(
                {
                    "gauge_rank": bloch_count,
                    "cap": cap,
                }
            )
            cap_parameters[position] = cap
            bloch_count += 1
        elif branch not in {"scalar-positive", "scalar-negative"}:
            raise ValueError(f"unsupported spectral branch {branch!r}")
        contractions.append(item)
    return tuple(contractions), cap_parameters


def enclosing_scaled_cap(
    grid: int,
    cap_indices: tuple[int, ...],
) -> tuple[np.ndarray, dict[str, Any]]:
    """Return one rigorously containing spherical cap for child caps.

    If a child cap has centre ``n_i`` and angular radius ``alpha_i``, every
    direction in it lies within ``angle(n,n_i)+alpha_i`` of any chosen parent
    centre ``n``.  The maximum of those radii therefore proves containment by
    the spherical triangle inequality.  The returned vector is ``n/cos(beta)``
    in the convention used by the contraction SOC.
    """

    indices = tuple(sorted(set(int(index) for index in cap_indices)))
    if not indices:
        raise ValueError("an enclosing cap needs at least one child")
    caps = cube_face_caps(int(grid))
    if indices[0] < 0 or indices[-1] >= len(caps):
        raise ValueError("cap index outside the cube-face cover")
    normals = np.asarray([caps[index][0] for index in indices], dtype=float)
    center = np.sum(normals, axis=0)
    norm = float(np.linalg.norm(center))
    if norm <= 1e-14:
        raise ValueError("child caps do not admit the mean-direction parent")
    center /= norm
    child_radii = [math.acos(float(caps[index][1])) for index in indices]
    parent_radius = max(
        math.acos(float(np.clip(center @ caps[index][0], -1.0, 1.0)))
        + child_radius
        for index, child_radius in zip(indices, child_radii, strict=True)
    )
    cosine = float(np.nextafter(math.cos(parent_radius), -math.inf))
    if cosine <= 0.0:
        raise ValueError("the enclosing cap is not contained in an open hemisphere")
    return center / cosine, {
        "child_indices": list(indices),
        "normal": center.tolist(),
        "cosine": cosine,
        "angular_radius": float(parent_radius),
    }


def set_cell_caps(
    source: dict[str, Any],
    pattern: tuple[str, ...],
    cap_indices: tuple[object, ...],
    parameters: dict[int, cp.Parameter],
) -> None:
    """Assign scaled angular-cap normals for one cell."""

    grids = tuple(int(value) for value in source["separator_grids"])
    if not (len(pattern) == len(cap_indices) == len(grids)):
        raise ValueError("cell branch, cap, and grid arrays are not parallel")
    for position, (branch, cap_index, grid) in enumerate(
        zip(pattern, cap_indices, grids, strict=True)
    ):
        if branch == "bloch":
            if cap_index is None or position not in parameters:
                raise ValueError("a Bloch branch requires a cap parameter")
            if isinstance(cap_index, (list, tuple)):
                scaled, _ = enclosing_scaled_cap(
                    grid, tuple(int(index) for index in cap_index)
                )
                parameters[position].value = scaled
            else:
                normal, cosine = cube_face_caps(grid)[int(cap_index)]
                if cosine <= 0.0:
                    raise ValueError("spectral cap cosine must be positive")
                parameters[position].value = np.asarray(normal, dtype=float) / float(
                    cosine
                )
        elif cap_index is not None:
            raise ValueError("a scalar branch cannot carry a cap index")


def group_box_and_reconstruction(
    source: dict[str, Any],
) -> tuple[Box, tuple[np.ndarray, np.ndarray]]:
    box = deserialise_box(source["box"])
    anchor, errors, _ = reconstruction_anchor_and_errors(
        box["terminal_alpha"], box["terminal_beta"]
    )
    return box, (anchor, errors)


def build_localisation_oracle(
    source: dict[str, Any],
    pattern: tuple[str, ...],
) -> tuple[TernaryConeOracle, dict[int, cp.Parameter], Box, dict[str, Any]]:
    """Build a reusable target-level support oracle for one branch pattern."""

    contractions, cap_parameters = parameterized_contractions(source, pattern)
    box, reconstruction = group_box_and_reconstruction(source)
    started = time.perf_counter()
    oracle = TernaryConeOracle(
        **_oracle_keywords(source),
        common_contractions=contractions,
        terminal_reconstruction=reconstruction,
        build_input_region_problem=True,
    )
    audit = {
        "build_seconds": time.perf_counter() - started,
        "problem_is_dpp": bool(oracle.problem.is_dpp()),
        "support_problem_is_dpp": bool(oracle.input_region_problem.is_dpp()),
        "mixed_integer": bool(oracle.problem.is_mixed_integer()),
    }
    if not audit["problem_is_dpp"] or not audit["support_problem_is_dpp"]:
        raise RuntimeError("parameterized localisation problems must be DPP")
    if audit["mixed_integer"]:
        raise RuntimeError("explicit branch patterns must remain continuous")
    return oracle, cap_parameters, box, audit


def localise_cell(
    source: dict[str, Any],
    cell: dict[str, Any],
    oracle: TernaryConeOracle,
    cap_parameters: dict[int, cp.Parameter],
    box: Box,
    coordinate_safety: float,
) -> dict[str, Any]:
    """Enclose every relaxed input in one cell that can reach the target."""

    pattern = tuple(cell["branches"])
    set_cell_caps(source, pattern, tuple(cell["caps"]), cap_parameters)
    started = time.perf_counter()
    base = oracle.solve(box, float(source["safety"]), capture=False)
    target = float(source["target"])
    if base["status"] in {cp.INFEASIBLE, cp.INFEASIBLE_INACCURATE} or float(
        base["bound"]
    ) < target:
        return {
            "status": "base-closed",
            "base_result": _strict_result(base),
            "runtime_seconds": time.perf_counter() - started,
        }
    if base["status"] not in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}:
        return {
            "status": "solver-unresolved",
            "base_result": _strict_result(base),
            "runtime_seconds": time.perf_counter() - started,
        }

    lower = np.empty((4, 4), dtype=float)
    upper = np.empty((4, 4), dtype=float)
    supports: list[dict[str, Any]] = []
    for z in range(4):
        for mu in range(4):
            direction = np.zeros((4, 4), dtype=float)
            direction[z, mu] = 1.0
            positive = oracle.solve_input_support(
                target, direction, coordinate_safety
            )
            negative = oracle.solve_input_support(
                target, -direction, coordinate_safety
            )
            if positive["status"] not in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE} or negative[
                "status"
            ] not in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}:
                return {
                    "status": "support-unresolved",
                    "base_result": _strict_result(base),
                    "failed_coordinate": [z, mu],
                    "positive": _strict_result(positive),
                    "negative": _strict_result(negative),
                    "runtime_seconds": time.perf_counter() - started,
                }
            upper[z, mu] = float(positive["bound"])
            lower[z, mu] = -float(negative["bound"])
            supports.append(
                {
                    "coordinate": [z, mu],
                    "lower": float(lower[z, mu]),
                    "upper": float(upper[z, mu]),
                    "positive": _strict_result(positive),
                    "negative": _strict_result(negative),
                }
            )
    if np.any(lower > upper):
        return {
            "status": "empty-box-unresolved",
            "base_result": _strict_result(base),
            "runtime_seconds": time.perf_counter() - started,
        }
    return {
        "status": "localized",
        "base_result": _strict_result(base),
        "lower": lower.tolist(),
        "upper": upper.tolist(),
        "row_l1_widths": np.sum(upper - lower, axis=1).tolist(),
        "maximum_coordinate_width": float(np.max(upper - lower)),
        "coordinate_safety": float(coordinate_safety),
        "supports": supports,
        "runtime_seconds": time.perf_counter() - started,
    }


def build_product_oracle(
    source: dict[str, Any],
    pattern: tuple[str, ...],
) -> tuple[
    TernaryConeOracle,
    dict[int, cp.Parameter],
    Box,
    cp.Parameter,
    cp.Parameter,
    cp.Parameter,
    dict[str, Any],
]:
    """Build a reusable conic-RLT product oracle for one branch pattern."""

    contractions, cap_parameters = parameterized_contractions(source, pattern)
    box, reconstruction = group_box_and_reconstruction(source)
    terminal_effects, terminal_errors, _, _ = terminal_effect_anchor_and_errors(
        box["terminal_alpha"], box["terminal_beta"]
    )
    lower = cp.Parameter((4, 4), name="input_lower")
    upper = cp.Parameter((4, 4), name="input_upper")
    purity = cp.Parameter((4, 4), name="input_purity_caps")
    started = time.perf_counter()
    oracle = TernaryConeOracle(
        **_oracle_keywords(source),
        common_contractions=contractions,
        terminal_reconstruction=reconstruction,
        input_pauli_lower=lower,
        input_pauli_upper=upper,
        input_purity_caps=purity,
        common_povm_bilinear=True,
        common_instrument_terminal_effect_anchor=terminal_effects,
        common_instrument_terminal_effect_errors=terminal_errors,
        common_instrument_product_trace_rules=True,
        common_instrument_product_psd_sandwiches=True,
        common_instrument_product_state_choi_psd=True,
        common_instrument_product_state_choi_ppt=True,
        max_common_instrument_witnesses=0,
    )
    audit = {
        "build_seconds": time.perf_counter() - started,
        "problem_is_dpp": bool(oracle.problem.is_dpp()),
        "mixed_integer": bool(oracle.problem.is_mixed_integer()),
        "constraint_count": len(oracle.problem.constraints),
    }
    if not audit["problem_is_dpp"]:
        raise RuntimeError("parameterized product problem must be DPP")
    if audit["mixed_integer"]:
        raise RuntimeError("explicit branch patterns must remain continuous")
    return oracle, cap_parameters, box, lower, upper, purity, audit


def _strict_result(result: dict[str, Any]) -> dict[str, Any]:
    bound, bound_class = strict_extended_real(float(result["bound"]))
    raw, raw_class = (
        strict_extended_real(float(result["raw_value"]))
        if result.get("raw_value") is not None
        else (None, "missing")
    )
    return {
        "status": str(result["status"]),
        "bound": bound,
        "bound_class": bound_class,
        "raw_value": raw,
        "raw_value_class": raw_class,
        "iterations": (
            int(result["iterations"]) if result.get("iterations") is not None else None
        ),
        **(
            {"error": str(result["error"])} if result.get("error") is not None else {}
        ),
    }


def _spatial_node(
    identifier: int,
    parent: int | None,
    depth: int,
    lower: np.ndarray,
    upper: np.ndarray,
    parent_bound: float,
) -> dict[str, Any]:
    return {
        "identifier": int(identifier),
        "parent": parent,
        "depth": int(depth),
        "lower": lower.tolist(),
        "upper": upper.tolist(),
        "parent_bound": float(parent_bound),
    }


def cover_localised_cell(
    source: dict[str, Any],
    cell: dict[str, Any],
    localisation: dict[str, Any],
    oracle: TernaryConeOracle,
    cap_parameters: dict[int, cp.Parameter],
    box: Box,
    lower_parameter: cp.Parameter,
    upper_parameter: cp.Parameter,
    purity_parameter: cp.Parameter,
    max_nodes: int,
    bound_safety: float,
    minimum_width: float,
    previous: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Spatially cover one localized input box with the reusable SDP."""

    if max_nodes <= 0:
        raise ValueError("max_nodes must be positive")
    pattern = tuple(cell["branches"])
    set_cell_caps(source, pattern, tuple(cell["caps"]), cap_parameters)
    previous_runtime = 0.0
    if previous is None:
        root_lower = np.asarray(localisation["lower"], dtype=float)
        root_upper = np.asarray(localisation["upper"], dtype=float)
        root_parent_bound = float(localisation["base_result"]["bound"])
        root = _spatial_node(0, None, 0, root_lower, root_upper, root_parent_bound)
        pending: list[tuple[float, int, dict[str, Any]]] = [
            (-root_parent_bound, 0, root)
        ]
        records: list[dict[str, Any]] = []
        unresolved: list[dict[str, Any]] = []
        next_identifier = 1
    else:
        records = list(previous["records"])
        unresolved = list(previous["unresolved"])
        pending = [
            (-float(node["parent_bound"]), int(node["identifier"]), node)
            for node in previous["pending"]
        ]
        heapq.heapify(pending)
        identifiers = [
            *(int(record["identifier"]) for record in records),
            *(int(node["identifier"]) for _, _, node in pending),
        ]
        next_identifier = max(identifiers, default=-1) + 1
        previous_runtime = float(previous.get("runtime_seconds", 0.0))
    target = float(source["target"])
    started = time.perf_counter()

    while pending and len(records) < max_nodes:
        _, _, node = heapq.heappop(pending)
        lower = np.asarray(node["lower"], dtype=float)
        upper = np.asarray(node["upper"], dtype=float)
        lower_parameter.value = lower
        upper_parameter.value = upper
        purity_parameter.value = box_purity_caps(lower, upper)
        result = oracle.solve(box, bound_safety, capture=False)
        bound = float(result["bound"])
        strict_result = _strict_result(result)
        widths = upper - lower
        residual_scores = (
            product_residual_scores(oracle)
            if result["status"] in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}
            and math.isfinite(bound)
            else np.zeros((4, 4), dtype=float)
        )
        record = {
            **node,
            **strict_result,
            "maximum_coordinate_width": float(np.max(widths)),
            "maximum_row_l1_width": float(np.max(np.sum(widths, axis=1))),
            "maximum_product_residual": float(np.max(residual_scores)),
        }
        records.append(record)
        if result["status"] in {cp.INFEASIBLE, cp.INFEASIBLE_INACCURATE} or (
            math.isfinite(bound) and bound < target
        ):
            record["disposition"] = "closed"
            record["branching_rule"] = "none"
            continue
        if result["status"] not in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE} or not math.isfinite(
            bound
        ):
            record["disposition"] = "solver-unresolved"
            record["branching_rule"] = "none"
            unresolved.append(record)
            continue

        eligible = residual_scores.copy()
        eligible[widths <= minimum_width] = -math.inf
        if np.any(np.isfinite(eligible)) and float(np.max(eligible)) > 1e-10:
            row, coordinate = np.unravel_index(int(np.argmax(eligible)), (4, 4))
            row, coordinate = int(row), int(coordinate)
            record["branching_rule"] = "product-residual"
        else:
            row, coordinate = _split_coordinate(lower, upper)
            record["branching_rule"] = "widest-coordinate"
        width = float(widths[row, coordinate])
        if width <= minimum_width:
            record["disposition"] = "resolution-limit"
            unresolved.append(record)
            continue
        record["disposition"] = "split"
        record["split_coordinate"] = [row, coordinate]
        midpoint = 0.5 * (lower[row, coordinate] + upper[row, coordinate])
        record["split_value"] = float(midpoint)
        for side in range(2):
            child_lower = lower.copy()
            child_upper = upper.copy()
            if side == 0:
                child_upper[row, coordinate] = midpoint
            else:
                child_lower[row, coordinate] = midpoint
            child = _spatial_node(
                next_identifier,
                int(node["identifier"]),
                int(node["depth"]) + 1,
                child_lower,
                child_upper,
                bound,
            )
            heapq.heappush(pending, (-bound, next_identifier, child))
            next_identifier += 1

    pending_nodes = [item[2] for item in sorted(pending)]
    dispositions = collections.Counter(record["disposition"] for record in records)
    statuses = collections.Counter(record["status"] for record in records)
    maximum_pending = max(
        (float(node["parent_bound"]) for node in pending_nodes), default=-math.inf
    )
    maximum_closed = max(
        (
            float(record["bound"])
            for record in records
            if record["disposition"] == "closed" and record["bound"] is not None
        ),
        default=-math.inf,
    )
    cover_upper = max(maximum_pending, maximum_closed)
    statuses_complete = all(status in ACCEPTED_STATUSES for status in statuses)
    complete = not pending_nodes and not unresolved and statuses_complete
    cover_json, cover_class = strict_extended_real(cover_upper)
    pending_json, pending_class = strict_extended_real(maximum_pending)
    closed_json, closed_class = strict_extended_real(maximum_closed)
    return {
        "solved_nodes": len(records),
        "closed_nodes": int(dispositions["closed"]),
        "split_nodes": int(dispositions["split"]),
        "pending_nodes": len(pending_nodes),
        "unresolved_nodes": len(unresolved),
        "maximum_depth": max((int(record["depth"]) for record in records), default=0),
        "maximum_pending_bound": pending_json,
        "maximum_pending_bound_class": pending_class,
        "maximum_closed_bound": closed_json,
        "maximum_closed_bound_class": closed_class,
        "cover_upper_bound": cover_json,
        "cover_upper_bound_class": cover_class,
        "statuses_complete": bool(statuses_complete),
        "complete": bool(complete),
        "target_closed": bool(complete and cover_upper < target),
        "dispositions": dict(sorted(dispositions.items())),
        "statuses": dict(sorted(statuses.items())),
        "runtime_seconds": previous_runtime + time.perf_counter() - started,
        "resumed": previous is not None,
        "records": records,
        "pending": pending_nodes,
        "unresolved": unresolved,
    }


def cell_descriptor(cell: dict[str, Any]) -> dict[str, Any]:
    return {
        "source_index": int(cell["source_index"]),
        "source_cell": int(cell["source_cell"]),
        "pattern": pattern_code(tuple(cell["branches"])),
        "branches": list(cell["branches"]),
        "caps": list(cell["caps"]),
        "source_status": str(cell["source_status"]),
        "source_bound": float(cell["source_bound"]),
        "source_audit": float(cell["source_audit"]),
        "source_return": float(cell["source_return"]),
    }


def write_payload(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def update_summary(
    payload: dict[str, Any],
    selected: list[dict[str, Any]],
    all_open: list[dict[str, Any]],
) -> None:
    localisations = payload["localisations"]
    results = payload["results"]
    selected_indices = {str(cell["source_index"]) for cell in selected}
    unprocessed = [
        cell for cell in all_open if str(cell["source_index"]) not in results
    ]
    processed_bounds = [
        float(result["cover_upper_bound"])
        for key, result in results.items()
        if key in selected_indices and result.get("cover_upper_bound") is not None
    ]
    unprocessed_bounds = [float(cell["source_bound"]) for cell in unprocessed]
    aggregate = max([*processed_bounds, *unprocessed_bounds], default=-math.inf)
    aggregate_json, aggregate_class = strict_extended_real(aggregate)
    payload["summary"] = {
        "source_open_cells": len(all_open),
        "selected_cells": len(selected),
        "localized_cells": sum(
            localisations.get(str(cell["source_index"]), {}).get("status")
            == "localized"
            for cell in selected
        ),
        "processed_cells": sum(
            str(cell["source_index"]) in results for cell in selected
        ),
        "target_closed_cells": sum(
            bool(results.get(str(cell["source_index"]), {}).get("target_closed"))
            for cell in selected
        ),
        "open_or_unresolved_selected_cells": sum(
            str(cell["source_index"]) in results
            and not bool(results[str(cell["source_index"])].get("target_closed"))
            for cell in selected
        ),
        "unprocessed_source_open_cells": len(unprocessed),
        "aggregate_source_cell_upper_bound": aggregate_json,
        "aggregate_source_cell_upper_bound_class": aggregate_class,
        "selected_base_angular_cell_closed": bool(
            len(selected) == len(all_open)
            and all(
                bool(results.get(str(cell["source_index"]), {}).get("target_closed"))
                for cell in all_open
            )
        ),
    }


def initial_payload(
    source_path: Path,
    source_raw: bytes,
    source: dict[str, Any],
    coordinate_safety: float,
    bound_safety: float,
    minimum_width: float,
    max_nodes: int,
) -> dict[str, Any]:
    return {
        "schema": SCHEMA,
        "source": {
            "filename": source_path.name,
            "sha256": hashlib.sha256(source_raw).hexdigest(),
            "bytes": len(source_raw),
            "scope": source["scope"],
            "support_weight": float(source["support_weight"]),
            "target": float(source["target"]),
            "base_code": str(source["base_code"]),
            "source_cell_count": len(source["cells"]),
            "source_statuses_complete": bool(source["statuses_complete"]),
        },
        "configuration": {
            "coordinate_safety": float(coordinate_safety),
            "bound_safety": float(bound_safety),
            "minimum_width": float(minimum_width),
            "max_nodes_per_cell": int(max_nodes),
            "common_instrument_product_trace_rules": True,
            "common_instrument_product_psd_sandwiches": True,
            "common_instrument_product_state_choi_psd": True,
            "common_instrument_product_state_choi_ppt": True,
        },
        "runtime": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "cvxpy": cp.__version__,
            "installed_solvers": cp.installed_solvers(),
        },
        "group_builds": {},
        "localisations": {},
        "results": {},
        "summary": {},
    }


def run_batch(
    source_path: Path,
    output: Path,
    limit: int,
    coordinate_safety: float,
    bound_safety: float,
    minimum_width: float,
    max_nodes: int,
    resume: bool,
) -> dict[str, Any]:
    source_raw = source_path.read_bytes()
    source = json.loads(source_raw)
    if not source.get("statuses_complete"):
        raise ValueError("source spectral frontier has unresolved solver statuses")
    all_open = source_open_cells(source)
    selected = all_open if limit <= 0 else all_open[:limit]
    if resume:
        payload = json.loads(output.read_bytes())
        if payload.get("schema") != SCHEMA:
            raise ValueError("resume payload has the wrong schema")
        if payload["source"]["sha256"] != hashlib.sha256(source_raw).hexdigest():
            raise ValueError("resume source hash does not match")
        previous_configuration = payload["configuration"]
        for key, value in (
            ("coordinate_safety", coordinate_safety),
            ("bound_safety", bound_safety),
            ("minimum_width", minimum_width),
        ):
            if float(previous_configuration[key]) != float(value):
                raise ValueError(f"resume cannot change {key}")
        payload["configuration"]["max_nodes_per_cell"] = int(max_nodes)
    else:
        payload = initial_payload(
            source_path,
            source_raw,
            source,
            coordinate_safety,
            bound_safety,
            minimum_width,
            max_nodes,
        )

    groups: dict[tuple[str, ...], list[dict[str, Any]]] = collections.defaultdict(list)
    for cell in selected:
        groups[tuple(cell["branches"])].append(cell)
    ordered_groups = sorted(
        groups.items(),
        key=lambda item: -max(cell["source_bound"] for cell in item[1]),
    )

    for pattern, cells in ordered_groups:
        code = pattern_code(pattern)
        missing_localisations = [
            cell
            for cell in cells
            if str(cell["source_index"]) not in payload["localisations"]
        ]
        if missing_localisations:
            oracle, cap_parameters, box, audit = build_localisation_oracle(
                source, pattern
            )
            payload["group_builds"].setdefault(code, {})["localisation"] = audit
            for cell in missing_localisations:
                localisation = localise_cell(
                    source,
                    cell,
                    oracle,
                    cap_parameters,
                    box,
                    coordinate_safety,
                )
                key = str(cell["source_index"])
                payload["localisations"][key] = {
                    "cell": cell_descriptor(cell),
                    **localisation,
                }
                if localisation["status"] == "base-closed":
                    payload["results"][key] = {
                        "cell": cell_descriptor(cell),
                        "closure_method": "base-relaxation",
                        "cover_upper_bound": localisation["base_result"]["bound"],
                        "cover_upper_bound_class": localisation["base_result"][
                            "bound_class"
                        ],
                        "target_closed": True,
                        "complete": True,
                    }
                update_summary(payload, selected, all_open)
                write_payload(output, payload)
                print(
                    json.dumps(
                        {
                            "phase": "localise",
                            "source_index": cell["source_index"],
                            "pattern": code,
                            "source_bound": cell["source_bound"],
                            "status": localisation["status"],
                            "maximum_coordinate_width": localisation.get(
                                "maximum_coordinate_width"
                            ),
                        }
                    ),
                    flush=True,
                )
            del oracle, cap_parameters
            gc.collect()

        pending_covers = [
            cell
            for cell in cells
            if payload["localisations"].get(str(cell["source_index"]), {}).get(
                "status"
            )
            == "localized"
            and (
                str(cell["source_index"]) not in payload["results"]
                or (
                    not payload["results"][str(cell["source_index"])].get(
                        "target_closed", False
                    )
                    and payload["results"][str(cell["source_index"])].get(
                        "pending_nodes", 0
                    )
                    > 0
                    and payload["results"][str(cell["source_index"])].get(
                        "solved_nodes", 0
                    )
                    < max_nodes
                )
            )
        ]
        if pending_covers:
            (
                oracle,
                cap_parameters,
                box,
                lower,
                upper,
                purity,
                audit,
            ) = build_product_oracle(source, pattern)
            payload["group_builds"].setdefault(code, {})["product"] = audit
            for cell in pending_covers:
                key = str(cell["source_index"])
                previous = payload["results"].get(key)
                cover = cover_localised_cell(
                    source,
                    cell,
                    payload["localisations"][key],
                    oracle,
                    cap_parameters,
                    box,
                    lower,
                    upper,
                    purity,
                    max_nodes,
                    bound_safety,
                    minimum_width,
                    previous,
                )
                payload["results"][key] = {
                    "cell": cell_descriptor(cell),
                    "closure_method": "common-instrument-product-localizer",
                    **cover,
                }
                update_summary(payload, selected, all_open)
                write_payload(output, payload)
                print(
                    json.dumps(
                        {
                            "phase": "cover",
                            "source_index": cell["source_index"],
                            "pattern": code,
                            "solved_nodes": cover["solved_nodes"],
                            "pending_nodes": cover["pending_nodes"],
                            "cover_upper_bound": cover["cover_upper_bound"],
                            "target_closed": cover["target_closed"],
                        }
                    ),
                    flush=True,
                )
            del oracle, cap_parameters, lower, upper, purity
            gc.collect()

    update_summary(payload, selected, all_open)
    write_payload(output, payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--frontier-json", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--limit", type=int, default=10)
    parser.add_argument("--coordinate-safety", type=float, default=2e-6)
    parser.add_argument("--bound-safety", type=float, default=2e-6)
    parser.add_argument("--minimum-width", type=float, default=1e-6)
    parser.add_argument("--max-nodes-per-cell", type=int, default=100)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    result = run_batch(
        args.frontier_json,
        args.output,
        args.limit,
        args.coordinate_safety,
        args.bound_safety,
        args.minimum_width,
        args.max_nodes_per_cell,
        args.resume,
    )
    print(json.dumps(result["summary"], indent=2, allow_nan=False))


if __name__ == "__main__":
    main()

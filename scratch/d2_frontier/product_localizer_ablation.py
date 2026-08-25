"""Compare common-product localizers on identical frontier cells.

The production cover uses entrywise McCormick envelopes for products between
one input Pauli coordinate and either an effective-POVM coordinate or a Choi
matrix.  This ablation keeps every box, inherited witness, objective, and
solver tolerance fixed while enabling two exact strengthenings:

* coordinatewise product-sum identities from POVM completeness;
* matrix-order sandwiches and product trace-preservation for Choi products.

The output is a solver-conditional comparison, not an interval certificate.
"""

from __future__ import annotations

import argparse
import collections
import gc
import hashlib
import json
import math
from pathlib import Path
import platform
import time
from typing import Any

import cvxpy as cp
import numpy as np

from terminal_reconstruction_enclosure import terminal_effect_anchor_and_errors
from ternary_bilinear_instrument_input_cover import (
    _configuration,
    _oracle_keywords,
    box_purity_caps,
)
from ternary_probability_cone_cover import TernaryConeOracle


SCHEMA = "carmenq.product-localizer-ablation.v1"
MODES = {
    "baseline": (False, False, False),
    "povm-sums": (True, False, False),
    "choi-trace": (False, True, False),
    "choi-sandwich": (False, False, True),
    "choi-both": (False, True, True),
}


def finite_number(value: object) -> tuple[float | None, str]:
    """Encode one extended-real solver value as strict JSON."""

    number = float(value)
    if math.isfinite(number):
        return number, "finite"
    if number == math.inf:
        return None, "positive-infinity"
    if number == -math.inf:
        return None, "negative-infinity"
    return None, "not-a-number"


def build_oracle(
    source: dict[str, Any],
    top_spectral_cell: bool,
    povm_sum_rules: bool,
    instrument_trace_rules: bool,
    instrument_psd_sandwiches: bool,
) -> tuple[TernaryConeOracle, cp.Parameter, cp.Parameter, cp.Parameter, dict[str, Any]]:
    """Build one reusable oracle for a fixed localizer ablation."""

    box, contractions, reconstruction = _configuration(source, top_spectral_cell)
    terminal_effects, terminal_errors, _, _ = terminal_effect_anchor_and_errors(
        box["terminal_alpha"], box["terminal_beta"]
    )
    lower = cp.Parameter((4, 4))
    upper = cp.Parameter((4, 4))
    purity = cp.Parameter((4, 4))
    started = time.perf_counter()
    oracle = TernaryConeOracle(
        **_oracle_keywords(source),
        common_contractions=contractions,
        terminal_reconstruction=reconstruction,
        input_pauli_lower=lower,
        input_pauli_upper=upper,
        input_purity_caps=purity,
        common_povm_bilinear=True,
        common_povm_product_sum_rules=povm_sum_rules,
        common_instrument_terminal_effect_anchor=terminal_effects,
        common_instrument_terminal_effect_errors=terminal_errors,
        common_instrument_product_trace_rules=instrument_trace_rules,
        common_instrument_product_psd_sandwiches=instrument_psd_sandwiches,
        max_common_instrument_witnesses=24,
    )
    return oracle, lower, upper, purity, {
        "build_seconds": time.perf_counter() - started,
        "is_dpp": bool(oracle.problem.is_dpp()),
        "is_mixed_integer": bool(oracle.problem.is_mixed_integer()),
    }


def solve_cell(
    oracle: TernaryConeOracle,
    lower_parameter: cp.Parameter,
    upper_parameter: cp.Parameter,
    purity_parameter: cp.Parameter,
    box: dict[str, Any],
    node: dict[str, Any],
    safety: float,
) -> dict[str, Any]:
    """Solve one cell with its complete inherited cut list."""

    lower = np.asarray(node["lower"], dtype=float)
    upper = np.asarray(node["upper"], dtype=float)
    lower_parameter.value = lower
    upper_parameter.value = upper
    purity_parameter.value = box_purity_caps(lower, upper)
    witnesses = tuple(
        (
            np.asarray(item["coefficients"], dtype=float),
            float(item["bound"]),
        )
        for item in node.get("determinant_witnesses", [])
    )
    started = time.perf_counter()
    result = oracle.solve(
        box,
        safety,
        capture=False,
        common_instrument_witnesses=witnesses,
    )
    bound, bound_class = finite_number(result["bound"])
    raw_value, raw_value_class = (
        finite_number(result["raw_value"])
        if result.get("raw_value") is not None
        else (None, "missing")
    )
    return {
        "status": result["status"],
        "bound": bound,
        "bound_class": bound_class,
        "raw_value": raw_value,
        "raw_value_class": raw_value_class,
        "iterations": result.get("iterations"),
        "solve_seconds": time.perf_counter() - started,
        "inherited_witnesses": len(witnesses),
    }


def ablate(
    source: dict[str, Any],
    checkpoint: dict[str, Any],
    limit: int,
) -> dict[str, Any]:
    """Run every localizer mode on the same leading pending cells."""

    if limit <= 0:
        raise ValueError("cell limit must be positive")
    top_spectral_cell = bool(checkpoint["top_spectral_cell"])
    box, _, _ = _configuration(source, top_spectral_cell)
    safety = float(checkpoint["bound_safety"])
    pending = sorted(
        checkpoint["pending"],
        key=lambda node: (-float(node["parent_bound"]), int(node["identifier"])),
    )[:limit]
    if not pending:
        raise ValueError("checkpoint has no pending cells")

    cells = [
        {
            "identifier": int(node["identifier"]),
            "parent": node.get("parent"),
            "depth": int(node["depth"]),
            "parent_bound": float(node["parent_bound"]),
            "maximum_coordinate_width": float(
                np.max(
                    np.asarray(node["upper"], dtype=float)
                    - np.asarray(node["lower"], dtype=float)
                )
            ),
            "results": {},
        }
        for node in pending
    ]
    construction: dict[str, Any] = {}
    for name, (povm_sums, choi_trace, choi_sandwich) in MODES.items():
        oracle, lower, upper, purity, audit = build_oracle(
            source,
            top_spectral_cell,
            povm_sums,
            choi_trace,
            choi_sandwich,
        )
        construction[name] = audit
        for cell, node in zip(cells, pending, strict=True):
            cell["results"][name] = solve_cell(
                oracle,
                lower,
                upper,
                purity,
                box,
                node,
                safety,
            )
        del oracle, lower, upper, purity
        gc.collect()

    for cell in cells:
        results = cell["results"]
        baseline = results["baseline"]["bound"]
        improvements = {
            name: (
                float(baseline) - float(results[name]["bound"])
                if baseline is not None and results[name]["bound"] is not None
                else None
            )
            for name in MODES
            if name != "baseline"
        }
        cell["improvements_from_baseline"] = improvements

    finite_improvements = {
        name: [
            float(cell["improvements_from_baseline"][name])
            for cell in cells
            if cell["improvements_from_baseline"][name] is not None
        ]
        for name in MODES
        if name != "baseline"
    }
    target = float(checkpoint["target"])
    return {
        "schema": SCHEMA,
        "runtime": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "cvxpy": cp.__version__,
            "installed_solvers": cp.installed_solvers(),
        },
        "source": {
            "support_weight": float(checkpoint["support_weight"]),
            "target": float(checkpoint["target"]),
            "top_spectral_cell": top_spectral_cell,
            "bound_safety": safety,
            "checkpoint_solved_nodes": int(checkpoint["solved_nodes"]),
            "checkpoint_maximum_pending_bound": float(
                checkpoint["maximum_pending_bound"]
            ),
        },
        "construction": construction,
        "cell_count": len(cells),
        "cells": cells,
        "summary": {
            name: {
                "comparable_cells": len(values),
                "minimum_improvement": min(values) if values else None,
                "maximum_improvement": max(values) if values else None,
                "mean_improvement": float(np.mean(values)) if values else None,
                "strict_improvement_cells_at_1e-8": sum(
                    value > 1e-8 for value in values
                ),
                "finite_bound_cells": sum(
                    cell["results"][name]["bound"] is not None for cell in cells
                ),
                "finite_bounds_below_target": sum(
                    cell["results"][name]["bound"] is not None
                    and float(cell["results"][name]["bound"]) < target
                    for cell in cells
                ),
                "status_counts": dict(
                    sorted(
                        collections.Counter(
                            str(cell["results"][name]["status"]) for cell in cells
                        ).items()
                    )
                ),
            }
            for name, values in {
                "baseline": [],
                **finite_improvements,
            }.items()
        },
        "interpretation": (
            "positive improvements tighten the solver-conditional upper bound; "
            "small negative values can arise from conic-solver tolerances"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--frontier-json", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--limit", type=int, default=8)
    args = parser.parse_args()
    frontier_raw = args.frontier_json.read_bytes()
    checkpoint_raw = args.checkpoint.read_bytes()
    source = json.loads(frontier_raw)
    checkpoint = json.loads(checkpoint_raw)
    result = ablate(source, checkpoint, args.limit)
    result["inputs"] = {
        "frontier": {
            "filename": args.frontier_json.name,
            "sha256": hashlib.sha256(frontier_raw).hexdigest(),
            "bytes": len(frontier_raw),
        },
        "checkpoint": {
            "filename": args.checkpoint.name,
            "sha256": hashlib.sha256(checkpoint_raw).hexdigest(),
            "bytes": len(checkpoint_raw),
        },
    }
    encoded = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if args.output is None:
        print(encoded, end="")
    else:
        args.output.write_text(encoded, encoding="utf-8")


if __name__ == "__main__":
    main()

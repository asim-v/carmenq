"""Rigorous local sensitivity of ternary duals to a projective support line.

This diagnostic replays stored cone-feasible dual vectors after changing the
named lambda=3/5 projective premise.  Exact stationarity residuals make every
reported number a valid upper bound for the selected cell.  It is deliberately
not a full-cover certificate: only the cells with the largest recorded bounds
are replayed, and the old dual vectors need not remain sharp.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
from typing import Any

import numpy as np

import ternary_socp_exact_dual_cover as cover


ROOT = Path(__file__).resolve().parent


def fraction_pair(value: Fraction) -> list[int]:
    return [value.numerator, value.denominator]


def selected_cells(paths: list[Path], count: int) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for path in paths:
        artifact = json.loads(path.read_text(encoding="utf-8"))
        for cell in artifact["cells"]:
            if "upper_fraction" not in cell:
                continue
            rows.append(
                {
                    "recorded_upper": Fraction(*cell["upper_fraction"]),
                    "path": path,
                    "artifact": artifact,
                    "cell": cell,
                }
            )
    rows.sort(key=lambda row: row["recorded_upper"], reverse=True)
    return rows[:count]


def replay(
    rows: list[dict[str, Any]], line_055: Fraction, line_060: Fraction
) -> dict[str, Any]:
    cover.LINE_055_UPPER = line_055
    cover.LINE_060_UPPER = line_060
    certifier = cover.ReusableExactCertifier(Fraction(1))
    reports = []
    for row in rows:
        cell = row["cell"]
        artifact = row["artifact"]
        data, _ = certifier.data(cell["box"])
        vector = cover.decode_candidate(artifact["candidates"][int(cell["candidate"])])
        dual, _ = cover.repair_dual_cones(vector, data["dims"])
        upper, correction, residual = certifier.exact_upper(data, dual)
        reports.append(
            {
                "artifact": row["path"].name,
                "candidate": int(cell["candidate"]),
                "box": cell["box"],
                "upper_fraction": fraction_pair(upper),
                "upper_decimal": cover.fraction_decimal(upper),
                "residual_correction_fraction": fraction_pair(correction),
                "maximum_stationarity_residual_fraction": fraction_pair(residual),
            }
        )
    maximum = max(Fraction(*row["upper_fraction"]) for row in reports)
    return {
        "line_055_upper": str(line_055),
        "line_060_upper": str(line_060),
        "maximum_selected_upper_fraction": fraction_pair(maximum),
        "maximum_selected_upper_decimal": cover.fraction_decimal(maximum),
        "cells": reports,
    }


def dimensions_key(data: dict[str, Any]) -> tuple[int, int, tuple[int, ...]]:
    dims = data["dims"]
    return int(dims.zero), int(dims.nonneg), tuple(map(int, dims.soc))


def require_line_only_change(
    baseline: dict[str, Any], modified: dict[str, Any]
) -> list[Fraction]:
    first = baseline["A"]
    second = modified["A"]
    same_sparse_matrix = (
        first.shape == second.shape
        and np.array_equal(first.indptr, second.indptr)
        and np.array_equal(first.indices, second.indices)
        and np.array_equal(first.data, second.data)
    )
    if not same_sparse_matrix or not np.array_equal(
        baseline["c"], modified["c"]
    ):
        raise RuntimeError("projective-line change altered A or c")
    if dimensions_key(baseline) != dimensions_key(modified):
        raise RuntimeError("projective-line change altered cone dimensions")
    return [
        cover.q(float(after)) - cover.q(float(before))
        for before, after in zip(
            np.asarray(baseline["b"]), np.asarray(modified["b"]), strict=True
        )
    ]


def exact_delta_objective(
    delta_size: int,
    delta_support: list[tuple[int, Fraction]],
    dual: np.ndarray,
) -> Fraction:
    if delta_size != len(dual):
        raise RuntimeError("dual and RHS delta have incompatible lengths")
    return sum(
        (
            coefficient * cover.q(float(dual[index]))
            for index, coefficient in delta_support
        ),
        Fraction(0),
    )


def fast_full_replay(
    paths: list[Path], line_055: Fraction, line_060: Fraction
) -> dict[str, Any]:
    loaded = [
        (path, json.loads(path.read_text(encoding="utf-8")))
        for path in paths
    ]
    expected_premises = None
    finite_rows: list[tuple[Path, dict[str, Any], dict[str, Any]]] = []
    total_cells = 0
    for path, artifact in loaded:
        premises = artifact.get("named_projective_premises")
        if expected_premises is None:
            expected_premises = premises
        if premises != expected_premises:
            raise RuntimeError("baseline artifacts use different premises")
        total_cells += len(artifact["cells"])
        finite_rows.extend(
            (path, artifact, cell)
            for cell in artifact["cells"]
            if "upper_fraction" in cell
        )
    if expected_premises is None or not finite_rows:
        raise RuntimeError("no finite ternary cells found")
    baseline_055 = Fraction(expected_premises["11/20"])
    baseline_060 = Fraction(expected_premises["3/5"])

    cover.LINE_055_UPPER = baseline_055
    cover.LINE_060_UPPER = baseline_060
    baseline = cover.ReusableExactCertifier(Fraction(1))
    cover.LINE_055_UPPER = line_055
    cover.LINE_060_UPPER = line_060
    modified = cover.ReusableExactCertifier(Fraction(1))
    if baseline.objective_correction != modified.objective_correction:
        raise RuntimeError("projective-line change altered objective correction")

    representative_indices = sorted(
        {0, len(finite_rows) // 2, len(finite_rows) - 1}
    )
    delta_b: list[Fraction] | None = None
    for index in representative_indices:
        box = finite_rows[index][2]["box"]
        baseline_data, _ = baseline.data(box)
        modified_data, _ = modified.data(box)
        current = require_line_only_change(baseline_data, modified_data)
        if delta_b is None:
            delta_b = current
        elif current != delta_b:
            raise RuntimeError("projective-line RHS delta depends on the cell")
    if delta_b is None:
        raise RuntimeError("failed to construct a projective-line RHS delta")

    delta_support = [
        (index, coefficient)
        for index, coefficient in enumerate(delta_b)
        if coefficient
    ]
    candidate_deltas: dict[tuple[str, int], Fraction] = {}
    for path, artifact in loaded:
        for candidate_id, candidate in enumerate(artifact["candidates"]):
            vector = cover.decode_candidate(candidate)
            dual, _ = cover.repair_dual_cones(
                vector, baseline_data["dims"]
            )
            candidate_deltas[(str(path), candidate_id)] = (
                exact_delta_objective(len(delta_b), delta_support, dual)
            )

    maximum: Fraction | None = None
    maximum_record: dict[str, Any] | None = None
    for path, _, cell in finite_rows:
        candidate_id = int(cell["candidate"])
        upper = (
            Fraction(*cell["upper_fraction"])
            + candidate_deltas[(str(path), candidate_id)]
        )
        if maximum is None or upper > maximum:
            maximum = upper
            maximum_record = {
                "artifact": path.name,
                "candidate": candidate_id,
                "box": cell["box"],
                "upper_fraction": fraction_pair(upper),
                "upper_decimal": cover.fraction_decimal(upper),
            }
    if maximum is None or maximum_record is None:
        raise RuntimeError("full replay produced no finite maximum")
    return {
        "baseline_projective_premises": expected_premises,
        "line_055_upper": str(line_055),
        "line_060_upper": str(line_060),
        "source_cell_count": total_cells,
        "finite_cell_count": len(finite_rows),
        "candidate_count": sum(
            len(artifact["candidates"]) for _, artifact in loaded
        ),
        "line_only_canonical_checks": len(representative_indices),
        "rhs_delta_nonzeros": len(delta_support),
        "maximum_full_upper_fraction": fraction_pair(maximum),
        "maximum_full_upper_decimal": cover.fraction_decimal(maximum),
        "maximum_cell": maximum_record,
        "optimiser_called": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("artifacts", nargs="*", type=Path)
    parser.add_argument(
        "--line-055-levels", nargs="+", default=["0.7573"]
    )
    parser.add_argument(
        "--levels",
        nargs="+",
        default=["0.76591", "0.76595", "0.766", "0.7661", "0.7662"],
    )
    parser.add_argument("--top", type=int, default=8)
    parser.add_argument(
        "--fast-full",
        action="store_true",
        help="replay every stored upper through an exact RHS sensitivity",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.top < 1:
        parser.error("--top must be positive")
    paths = args.artifacts or sorted(
        ROOT.glob("ternary_socp_exact_dual_full_l060_shard*of08.json")
    )
    if not paths:
        parser.error("no ternary artifacts found")
    original_055 = cover.LINE_055_UPPER
    original_060 = cover.LINE_060_UPPER
    try:
        pairs = [
            (Fraction(line_055), Fraction(line_060))
            for line_055 in args.line_055_levels
            for line_060 in args.levels
        ]
        if args.fast_full:
            result = {
                "schema": (
                    "carmenq.ternary-projective-line-full-sensitivity.v1"
                ),
                "source_artifacts": [str(path) for path in paths],
                "replays": [
                    fast_full_replay(paths, line_055, line_060)
                    for line_055, line_060 in pairs
                ],
                "proof_status": (
                    "exact arithmetic replay of every stored baseline upper "
                    "after a line-only canonical RHS change; final component "
                    "verification must still replay the selected new premise"
                ),
                "trusted_optimizers": [],
            }
        else:
            rows = selected_cells(paths, args.top)
            result = {
                "schema": (
                    "carmenq.ternary-projective-line-local-sensitivity.v1"
                ),
                "selected_cell_count": len(rows),
                "source_artifacts": [str(path) for path in paths],
                "replays": [
                    replay(rows, line_055, line_060)
                    for line_055, line_060 in pairs
                ],
                "proof_status": (
                    "rigorous exact-residual upper bounds for the selected "
                    "cells; not a full-cover certificate"
                ),
                "trusted_optimizers": [],
            }
    finally:
        cover.LINE_055_UPPER = original_055
        cover.LINE_060_UPPER = original_060
    rendered = json.dumps(result, indent=2) + "\n"
    print(rendered, end="")
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

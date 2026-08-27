"""Clustered exact-residual dual replay of the complete ternary cover.

The numerical cover supplies only its dyadic terminal boxes.  A persistent
CVXPY/Clarabel model proposes dual vectors, while every accepted cell is
checked with the exact-rational kernel from
``ternary_socp_exact_dual_probe``.  One cone-feasible dual may certify many
neighbouring cells; the output stores each candidate once and refers to it by
index, avoiding a multi-gigabyte vector-per-leaf artifact.

The resulting theorem is conditional only on the two named binary-projective
support lines.  Clarabel is a search helper and is absent from the verifier's
logical trust boundary.
"""

from __future__ import annotations

import argparse
import base64
from fractions import Fraction
import json
import math
from pathlib import Path
from typing import Any
import zlib

import cvxpy as cp
import numpy as np

from cvxpy.reductions.solvers.conic_solvers.clarabel_conif import CLARABEL
from ternary_probability_cone_cover import (
    OUTCOMES,
    TERMINAL_ALPHA,
    TERMINAL_BETA,
    TernaryConeOracle,
    binary64_product_down,
    binary64_product_up,
    rational_to_binary64_up,
    terminal_weight_intervals,
)
from ternary_socp_exact_dual_probe import (
    canonical_hash,
    certified_terminal_family,
    exact_dot,
    exact_sparse_stationarity,
    fraction_decimal,
    objective_correction,
    q,
    repair_dual_cones,
    safe_line_upper,
)


ROOT = Path(__file__).resolve().parent
DEFAULT_COVER = ROOT / "continuous_terminal_projective_l055cert_complete.json"
TARGET = Fraction(76643, 100000)
SUPPORT_WEIGHT = Fraction(3, 5)
LINE_055_WEIGHT = Fraction(11, 20)
LINE_055_UPPER = Fraction(7573, 10000)
LINE_060_UPPER = Fraction(76591, 100000)
MAXIMUM_WEIGHT_FLOOR = Fraction(3533, 4000)


def box_key(box: dict[str, Any]) -> str:
    return json.dumps(box, sort_keys=True, separators=(",", ":"))


def box_center(box: dict[str, Any]) -> tuple[float, float]:
    return tuple(
        sum(map(float, box[name])) / 2.0
        for name in (TERMINAL_ALPHA, TERMINAL_BETA)
    )  # type: ignore[return-value]


def exact_domain_empty(box: dict[str, Any]) -> bool:
    """Certify the only pruning rule used by the archived cover."""

    au = q(float(box[TERMINAL_ALPHA][1]))
    bl = q(float(box[TERMINAL_BETA][0]))
    # Sorted Horwitz parameters require alpha >= beta.  This was the first
    # pruning branch in the archived cover and must be replayed exactly too.
    if au < bl:
        return True
    maximum_w0 = au / (au + bl - 1)
    return maximum_w0 < MAXIMUM_WEIGHT_FLOOR


def encode_dual_f32(vector: np.ndarray) -> str:
    raw = np.asarray(vector, dtype="<f4").tobytes()
    return base64.b64encode(zlib.compress(raw, level=9)).decode("ascii")


def decode_candidate(record: dict[str, Any]) -> np.ndarray:
    if "dual_f32_zlib_base64" in record:
        raw = zlib.decompress(base64.b64decode(record["dual_f32_zlib_base64"]))
        vector = np.frombuffer(raw, dtype="<f4").astype(float)
    elif "dual_hex" in record:
        vector = np.asarray([float.fromhex(value) for value in record["dual_hex"]])
    else:
        raise ValueError("candidate has no supported dual encoding")
    if vector.shape != (3187,) or not np.all(np.isfinite(vector)):
        raise ValueError("candidate dual has the wrong shape or non-finite entries")
    return vector


def sparsify_dual_f32(
    vector: np.ndarray, dims: Any, threshold: float
) -> np.ndarray:
    result = np.asarray(vector, dtype=np.float32).astype(float)
    cursor = int(dims.zero)
    end = cursor + int(dims.nonneg)
    nonnegative = result[cursor:end]
    nonnegative[nonnegative < threshold] = 0.0
    result[cursor:end] = np.maximum(nonnegative, 0.0)
    cursor = end
    for dimension in dims.soc:
        end = cursor + int(dimension)
        if np.max(np.abs(result[cursor:end])) < threshold:
            result[cursor:end] = 0.0
        cursor = end
    return np.asarray(result, dtype=np.float32).astype(float)


SPARSITY_SCHEDULE = (
    1e-5,
    3e-6,
    1e-6,
    3e-7,
    1e-7,
    3e-8,
    1e-8,
    3e-9,
    1e-9,
    3e-10,
    1e-10,
    0.0,
)


def compact_enclosure(report: dict[str, Any]) -> dict[str, Any]:
    margins = [
        Fraction(*item["margin"])
        for item in report["anchors"]
        if item.get("margin") is not None
    ]
    lower = report["coefficientwise_lower"]
    if lower.get("present"):
        margins.append(Fraction(*lower["margin"]))
    return {
        "all_certified": bool(report["all_certified"]),
        "inflated_anchor_count": sum(
            int(item["post_sqrt_nextafter_steps"] > 0)
            for item in report["anchors"]
        ),
        "minimum_margin_decimal": fraction_decimal(min(margins)),
        "coefficientwise_lower_present": bool(lower.get("present")),
    }


class ReusableExactCertifier:
    def __init__(self, target: Fraction) -> None:
        self.target = target
        encoded_weight = float(SUPPORT_WEIGHT)
        encoded_line_weight = float(LINE_055_WEIGHT)
        line_055, correction_055 = safe_line_upper(
            LINE_055_WEIGHT, LINE_055_UPPER, encoded_line_weight
        )
        line_060, correction_060 = safe_line_upper(
            SUPPORT_WEIGHT, LINE_060_UPPER, encoded_weight
        )
        self.objective_correction = objective_correction(
            SUPPORT_WEIGHT, encoded_weight
        )
        self.line_corrections = {
            "0.55": correction_055,
            "0.60": correction_060,
        }
        self.oracle = TernaryConeOracle(
            encoded_weight,
            (0, 1, 2, 3),
            (),
            (),
            float(MAXIMUM_WEIGHT_FLOOR),
            line_060,
            projective_support_lines=((encoded_line_weight, line_055),),
        )
        variables = self.oracle.problem.variables()
        if sum(variable.size for variable in variables) != 234:
            raise RuntimeError("unexpected canonical variable inventory")
        if not all(bool(variable.attributes.get("nonneg")) for variable in variables):
            raise RuntimeError("residual-controlled variables must be nonnegative")
        self.problem = cp.Problem(
            self.oracle.problem.objective,
            [
                *self.oracle.problem.constraints,
                *(variable <= 1.0 for variable in variables),
            ],
        )

    def assign(self, box: dict[str, Any]) -> dict[str, Any]:
        intervals = terminal_weight_intervals(box)
        weight_lower = np.asarray([item[0] for item in intervals])
        weight_upper = np.asarray([item[1] for item in intervals])
        self.oracle.weight_lower.value = weight_lower
        self.oracle.weight_upper.value = weight_upper
        self.oracle.maximum_weight_floor.value = float(MAXIMUM_WEIGHT_FLOOR)
        al, au = map(float, box[TERMINAL_ALPHA])
        bl, bu = map(float, box[TERMINAL_BETA])
        self.oracle.terminal_alpha_lower.value = al
        self.oracle.terminal_alpha_upper.value = au
        self.oracle.terminal_beta_lower.value = bl
        self.oracle.terminal_beta_upper.value = bu

        rank_upper = (
            1.0,
            0.5,
            rational_to_binary64_up(Fraction(1, 3)),
            0.25,
        )
        prior_lower = np.zeros(4, dtype=float)
        prior_upper = np.ones(4, dtype=float)
        for rank, z in enumerate(self.oracle.prefix_order):
            prior_upper[z] = rank_upper[rank]
        prior_lower[self.oracle.prefix_order[0]] = 0.25
        self.oracle.prior_lower.value = prior_lower
        self.oracle.prior_upper.value = prior_upper
        self.oracle.mccormick_lower_cross.value = np.asarray(
            [
                [
                    binary64_product_down(weight_lower[t], prior_upper[z])
                    for z in OUTCOMES
                ]
                for t in range(3)
            ]
        )
        self.oracle.mccormick_upper_cross.value = np.asarray(
            [
                [
                    binary64_product_up(weight_upper[t], prior_upper[z])
                    for z in OUTCOMES
                ]
                for t in range(3)
            ]
        )
        maximum = float(intervals[0][1])
        third_cap = rational_to_binary64_up(
            max(Fraction(0), 2 - 2 * Fraction.from_float(maximum))
        )
        self.oracle.cap_weights.value = np.asarray(
            [maximum, maximum, third_cap, 0.0]
        )

        anchor_data, lower_data, enclosure = certified_terminal_family(
            (al, au), (bl, bu)
        )
        anchor_targets, lower_target = self.oracle.soc_parameters[0]
        for data, targets in zip(anchor_data, anchor_targets, strict=True):
            self.oracle.assign_soc(targets, data)
        self.oracle.assign_soc(lower_target, lower_data)
        return enclosure

    def data(self, box: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
        enclosure = self.assign(box)
        data, _, _ = self.problem.get_problem_data(cp.CLARABEL)
        if data["A"].shape[1] != 234:
            raise RuntimeError("unexpected canonical column count")
        return data, enclosure

    def solve_candidate(self, data: dict[str, Any]) -> tuple[np.ndarray, str]:
        result = CLARABEL().solve_via_data(
            data,
            warm_start=False,
            verbose=False,
            solver_opts={
                "tol_gap_abs": 1e-11,
                "tol_gap_rel": 1e-11,
                "tol_feas": 1e-11,
                "max_iter": 500,
            },
            solver_cache=None,
        )
        dual, _ = repair_dual_cones(np.asarray(result.z), data["dims"])
        return dual, str(result.status)

    def numerical_screen(self, data: dict[str, Any], dual: np.ndarray) -> float:
        residual = np.asarray(data["A"].T @ dual + data["c"])
        return float(
            np.asarray(data["b"]) @ dual
            + np.maximum(-residual, 0.0).sum()
            + float(self.objective_correction)
        )

    def exact_upper(
        self, data: dict[str, Any], dual: np.ndarray
    ) -> tuple[Fraction, Fraction, Fraction]:
        residuals, residual_correction = exact_sparse_stationarity(
            data["A"], dual, data["c"]
        )
        dual_objective = exact_dot(np.asarray(data["b"]), dual)
        upper = (
            dual_objective
            + residual_correction
            + self.objective_correction
        )
        return upper, residual_correction, max(
            map(abs, residuals), default=Fraction(0)
        )


def compact_candidate(
    certifier: ReusableExactCertifier,
    data: dict[str, Any],
    raw_dual: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, Fraction, Fraction, Fraction, float]:
    """Choose the sparsest float32 candidate that still closes exactly."""

    for threshold in SPARSITY_SCHEDULE:
        stored = sparsify_dual_f32(raw_dual, data["dims"], threshold)
        checked, _ = repair_dual_cones(stored, data["dims"])
        upper, correction, residual = certifier.exact_upper(data, checked)
        if upper <= certifier.target:
            return stored, checked, upper, correction, residual, threshold
    raise ArithmeticError("no float32 sparsification level preserved closure")


def candidate_distance(
    center: tuple[float, float], candidate: dict[str, Any]
) -> float:
    other = tuple(map(float, candidate["center"]))
    return (center[0] - other[0]) ** 2 + (center[1] - other[1]) ** 2


def payload(
    cover: dict[str, Any],
    target: Fraction,
    total_cells: int,
    results: list[dict[str, Any]],
    candidates: list[dict[str, Any]],
) -> dict[str, Any]:
    closed = [item for item in results if item["closed"]]
    finite = [
        Fraction(item["upper_fraction"][0], item["upper_fraction"][1])
        for item in closed
        if "upper_fraction" in item
    ]
    return {
        "schema": "carmenq.ternary-socp-clustered-exact-dual-cover.v1",
        "support_weight": "3/5",
        "target": str(target),
        "named_projective_premises": {
            "11/20": str(LINE_055_UPPER),
            "3/5": str(LINE_060_UPPER),
        },
        "source_cover": str(DEFAULT_COVER.name),
        "source_leaf_count": len(cover["leaves"]),
        "selected_cell_count": total_cells,
        "processed_cell_count": len(results),
        "closed_cell_count": len(closed),
        "run_complete": len(results) == total_cells,
        "all_cells_closed": len(results) == total_cells and len(closed) == total_cells,
        "candidate_count": len(candidates),
        "maximum_certified_upper": (
            fraction_decimal(max(finite)) if finite else None
        ),
        "cells": results,
        "candidates": candidates,
        "trusted_optimizers": [],
        "untrusted_search_helpers": ["Clarabel dual-vector proposals"],
    }


def write_payload(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cover", type=Path, default=DEFAULT_COVER)
    parser.add_argument("--target", default="0.76643")
    parser.add_argument("--limit", type=int)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--candidate-trials", type=int, default=2)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    if args.candidate_trials < 0:
        parser.error("--candidate-trials must be nonnegative")
    if args.shard_count < 1:
        parser.error("--shard-count must be positive")
    if not 0 <= args.shard_index < args.shard_count:
        parser.error("--shard-index must lie in [0, shard-count)")
    target = Fraction(args.target)
    source = json.loads(args.cover.read_text(encoding="utf-8"))
    leaves = list(source["leaves"])
    leaves.sort(
        key=lambda item: (
            -float(item.get("raw_value", -math.inf)),
            *box_center(item["box"]),
        )
    )
    leaves = [
        leaf
        for position, leaf in enumerate(leaves)
        if position % args.shard_count == args.shard_index
    ]
    if args.limit is not None:
        leaves = leaves[: args.limit]
    if len({box_key(item["box"]) for item in leaves}) != len(leaves):
        raise RuntimeError("duplicate terminal boxes in source cover")

    results_by_box: dict[str, dict[str, Any]] = {}
    candidates: list[dict[str, Any]] = []
    if args.resume and args.output.exists():
        prior = json.loads(args.output.read_text(encoding="utf-8"))
        if prior.get("support_weight") != "3/5" or prior.get("target") != str(target):
            parser.error("resume artifact has incompatible support weight or target")
        results_by_box = {
            box_key(item["box"]): item
            for item in prior.get("cells", [])
            if item.get("closed")
        }
        candidates = list(prior.get("candidates", []))

    candidate_vectors = [decode_candidate(item) for item in candidates]
    certifier = ReusableExactCertifier(target)

    def ordered_results() -> list[dict[str, Any]]:
        return [
            results_by_box[box_key(leaf["box"])]
            for leaf in leaves
            if box_key(leaf["box"]) in results_by_box
        ]

    pending = [leaf for leaf in leaves if box_key(leaf["box"]) not in results_by_box]
    for index, leaf in enumerate(pending, start=1):
        box = leaf["box"]
        key = box_key(box)
        if leaf.get("status") == "domain_empty":
            closed = exact_domain_empty(box)
            result = {
                "box": box,
                "source_status": "domain_empty",
                "domain_empty_exact": closed,
                "closed": closed,
            }
            if not closed:
                raise RuntimeError("numerically pruned leaf is not exactly empty")
            results_by_box[key] = result
            continue

        data, enclosure = certifier.data(box)
        center = box_center(box)
        nearest = sorted(
            range(len(candidates)),
            key=lambda item: candidate_distance(center, candidates[item]),
        )[: args.candidate_trials]
        selected: int | None = None
        upper: Fraction | None = None
        residual_correction = Fraction(0)
        maximum_residual = Fraction(0)
        for candidate_id in nearest:
            dual, _ = repair_dual_cones(
                candidate_vectors[candidate_id], data["dims"]
            )
            if certifier.numerical_screen(data, dual) > float(target) - 1e-10:
                continue
            checked, correction, residual = certifier.exact_upper(data, dual)
            if checked <= target:
                selected = candidate_id
                upper = checked
                residual_correction = correction
                maximum_residual = residual
                break
        solver_status = None
        if selected is None:
            raw_dual, solver_status = certifier.solve_candidate(data)
            (
                stored_dual,
                dual,
                checked,
                correction,
                residual,
                sparsity_threshold,
            ) = compact_candidate(certifier, data, raw_dual)
            selected = len(candidates)
            upper = checked
            residual_correction = correction
            maximum_residual = residual
            candidates.append(
                {
                    "center": list(center),
                    "source_box": box,
                    "untrusted_solver_status": solver_status,
                    "encoding": "zlib-base64-little-endian-float32",
                    "sparsity_threshold": sparsity_threshold,
                    "dual_f32_zlib_base64": encode_dual_f32(stored_dual),
                }
            )
            candidate_vectors.append(stored_dual)
        assert upper is not None and selected is not None
        closed = upper <= target
        result = {
            "box": box,
            "source_raw_value": leaf.get("raw_value"),
            "candidate": selected,
            "candidate_created_here": solver_status is not None,
            "canonical_sha256": canonical_hash(data),
            "upper_fraction": [upper.numerator, upper.denominator],
            "upper_decimal": fraction_decimal(upper),
            "residual_correction_decimal": fraction_decimal(residual_correction),
            "maximum_stationarity_residual_decimal": fraction_decimal(maximum_residual),
            "inellipse": compact_enclosure(enclosure),
            "closed": closed,
        }
        results_by_box[key] = result
        print(
            json.dumps(
                {
                    "cell": index,
                    "pending": len(pending),
                    "processed_total": len(results_by_box),
                    "candidate": selected,
                    "new_candidate": solver_status is not None,
                    "candidate_count": len(candidates),
                    "upper": result["upper_decimal"],
                    "closed": closed,
                }
            ),
            flush=True,
        )
        if not closed:
            write_payload(
                args.output,
                payload(source, target, len(leaves), ordered_results(), candidates),
            )
            raise RuntimeError("exact residual dual failed to close a source leaf")
        if index % 25 == 0:
            write_payload(
                args.output,
                payload(source, target, len(leaves), ordered_results(), candidates),
            )

    final = payload(source, target, len(leaves), ordered_results(), candidates)
    write_payload(args.output, final)
    print(
        json.dumps({key: value for key, value in final.items() if key not in {"cells", "candidates"}}, indent=2)
    )
    if not final["all_cells_closed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

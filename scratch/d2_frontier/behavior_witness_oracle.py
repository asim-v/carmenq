"""Separate-process nested-ellipsoid oracle for behaviour SCIP runs.

Keeping this CVXPY oracle in a separate process avoids loading the OpenMP
runtimes used by CVXPY and SCIP into one Windows process.
"""

from __future__ import annotations

import json
import sys
from typing import Any

import cvxpy as cp
import numpy as np

from full_behavior_psd_rank_certificate import solve_primal
from full_behavior_witness_tree import find_small_witness


def normalise_cut(cut: dict[str, Any]) -> dict[str, Any]:
    coefficients = np.asarray(cut["coefficients"], dtype=float)
    scale = float(np.max(np.abs(coefficients)))
    if not np.isfinite(scale) or scale <= 0.0:
        raise ValueError("a cut coefficient vector is zero or non-finite")
    return {
        "column": int(cut["column"]),
        "coefficients": (coefficients / scale).tolist(),
    }


def witness_from_dual(dual: dict[str, Any]) -> dict[str, Any]:
    return {
        "cuts": [
            normalise_cut({"column": column, "coefficients": coefficients})
            for column, coefficients in zip(
                dual["active_columns"],
                dual["halfspace_linear_coefficients"],
                strict=True,
            )
        ],
        "common_margin": float(dual["certified_common_margin"]),
        "stationarity_residual": float(dual["stationarity_frobenius_residual"]),
        "state_dual_min_eigenvalue": float(dual["state_dual_min_eigenvalue"]),
        "containment_dual_min_eigenvalues": [
            float(value) for value in dual["containment_dual_min_eigenvalues"]
        ],
    }


def find_witness(
    behavior: np.ndarray,
    robust_budget: float,
    tolerance: float,
) -> tuple[dict[str, Any] | None, dict[str, Any]]:
    primal = solve_primal(behavior)
    report: dict[str, Any] = {"primal_status": primal["status"]}
    if primal["status"] in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}:
        return None, report

    dual = find_small_witness(
        behavior,
        robust_budget,
        tolerance,
        exhaustive_pairs=True,
        exhaustive_triples=True,
    )
    if dual is None:
        return None, report
    report["robust_quality"] = float(dual["robust_quality"])
    report["active_columns"] = dual["active_columns"]
    return witness_from_dual(dual), report


def main() -> None:
    request = json.loads(sys.stdin.read())
    witness, report = find_witness(
        np.asarray(request["behavior"], dtype=float),
        float(request.get("robust_budget", 1000.0)),
        float(request.get("tolerance", 2e-8)),
    )
    sys.stdout.write(json.dumps({"witness": witness, "report": report}))


if __name__ == "__main__":
    main()

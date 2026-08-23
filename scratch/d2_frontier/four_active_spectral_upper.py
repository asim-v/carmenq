"""Fast outer bound for a fully active four-effect terminal qubit readout.

The bound combines four facts about a sorted rank-one qubit POVM

    E_s = w_s |n_s><n_s|,  sum_s E_s = I_2.

First, every pulled-back prefix effect has norm at most ``max(w)``.  Second,
the syndrome success is bounded by ``sum_s w_s p_s``.  Third, Helstrom
complementarity gives the prior-reserve cap implemented in
``active_readout_audit_cap.py``.  Finally, replacing one effect by its support
projector and assigning the complementary projector to a second label gives

    A_4 - A_2 <= (1-w_i) p_j + sum_{k not in {i,j}} w_k p_k.

The same RETURN value is retained, while the binary-projective support is
bounded independently.  Writing the sixteen path probabilities as squares
of a unit vector turns every one of these inequalities into a quadratic
upper bound.  Any convex combination is still valid.  A tiny SDP chooses the
combination with the smallest largest eigenvalue.

This is a solver-level numerical upper certificate.  A continuum proof needs
an outward weight-box version and outward validation of the scalar/SDP duals.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import cvxpy as cp
import numpy as np

from active_readout_audit_cap import active_audit_cap
from terminal_weight_upper import filled_effect_weights


OUTCOMES = range(4)
PATHS = tuple((z, y) for z in OUTCOMES for y in OUTCOMES)


def score_matrix(audit_diagonal: np.ndarray, support_weight: float) -> np.ndarray:
    """Return ``lambda*diag(a)+(1-lambda)J/16``."""

    diagonal = np.asarray(audit_diagonal, dtype=float)
    if diagonal.shape != (16,):
        raise ValueError("an audit diagonal must have sixteen entries")
    return (
        support_weight * np.diag(diagonal)
        + (1.0 - support_weight) * np.ones((16, 16), dtype=float) / 16.0
    )


def named_upper_matrices(
    weights: np.ndarray,
    audit_upper: float,
    support_weight: float,
    projective_support_upper: float,
) -> list[tuple[str, np.ndarray]]:
    """Build every quadratic score upper used by the certificate."""

    value = np.asarray(weights, dtype=float)
    prefix_cap = filled_effect_weights(float(value.max()))
    prefix_diagonal = np.asarray([prefix_cap[z] for z, _ in PATHS])
    syndrome_diagonal = np.asarray([value[z ^ y] for z, y in PATHS])
    named = [
        ("prefix", score_matrix(prefix_diagonal, support_weight)),
        ("syndrome", score_matrix(syndrome_diagonal, support_weight)),
        (
            "active_prior_reserve",
            score_matrix(np.full(16, audit_upper), support_weight),
        ),
    ]
    for retained in OUTCOMES:
        for complement in OUTCOMES:
            if retained == complement:
                continue
            bonus = np.zeros(4, dtype=float)
            bonus[complement] = 1.0 - value[retained]
            for deleted in OUTCOMES:
                if deleted not in {retained, complement}:
                    bonus[deleted] = value[deleted]
            diagonal = np.asarray([bonus[z ^ y] for z, y in PATHS])
            matrix = (
                projective_support_upper * np.eye(16)
                + support_weight * np.diag(diagonal)
            )
            named.append(
                (
                    f"retain_{retained}_complement_{complement}",
                    matrix,
                )
            )
    return named


def spectral_upper(
    weights: np.ndarray,
    support_weight: float = 0.6,
    projective_support_upper: float = 0.76591,
) -> dict[str, Any]:
    """Return the combined fully-active four-effect spectral upper bound."""

    value = np.sort(np.asarray(weights, dtype=float))[::-1]
    if value.shape != (4,):
        raise ValueError("expected four terminal effect traces")
    if np.any(value <= 0.0) or np.any(value > 1.0 + 1e-12):
        raise ValueError("fully active effect traces must lie in (0,1]")
    if abs(float(value.sum()) - 2.0) > 1e-9:
        raise ValueError("qubit rank-one effect traces must sum to two")

    active = active_audit_cap(value)
    named = named_upper_matrices(
        value,
        float(active["audit_upper"]),
        support_weight,
        projective_support_upper,
    )
    multipliers = cp.Variable(len(named), nonneg=True)
    level = cp.Variable()
    combined = sum(
        multipliers[index] * matrix
        for index, (_, matrix) in enumerate(named)
    )
    problem = cp.Problem(
        cp.Minimize(level),
        [cp.sum(multipliers) == 1.0, level * np.eye(16) - combined >> 0],
    )
    problem.solve(
        solver="CLARABEL",
        tol_gap_abs=1e-10,
        tol_gap_rel=1e-10,
        tol_feas=1e-10,
        max_iter=1000,
    )
    if problem.status not in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}:
        raise RuntimeError(f"spectral SDP failed with status {problem.status}")
    coefficients = np.asarray(multipliers.value, dtype=float)
    direct = sum(
        coefficients[index] * matrix
        for index, (_, matrix) in enumerate(named)
    )
    return {
        "weights": value.tolist(),
        "support_weight": float(support_weight),
        "projective_support_upper": float(projective_support_upper),
        "active_audit_cap": active,
        "bound": float(problem.value),
        "direct_maximum_eigenvalue": float(np.linalg.eigvalsh(direct)[-1]),
        "multipliers": {
            name: float(coefficients[index])
            for index, (name, _) in enumerate(named)
            if coefficients[index] > 1e-9
        },
        "solver_status": problem.status,
        "numerical_status": (
            "fixed-weight SDP upper at solver tolerances; scalar and SDP "
            "outward validation pending"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--weights", type=float, nargs=4, required=True)
    parser.add_argument("--lambda", dest="support_weight", type=float, default=0.6)
    parser.add_argument("--projective-upper", type=float, default=0.76591)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = spectral_upper(
        np.asarray(args.weights, dtype=float),
        args.support_weight,
        args.projective_upper,
    )
    rendered = json.dumps(payload, indent=2) + "\n"
    print(rendered, end="")
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

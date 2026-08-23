"""Full four-preparation qubit-behaviour certificate.

Let ``B`` be a row-stochastic ``4 x m`` behaviour.  In barycentric
coordinates the convex hull of its four rows is the fixed tetrahedron

    Delta = conv(0, e1, e2, e3),

whereas the intersection of the row affine hull with the probability
simplex is

    Q_B = {x : b_3j + sum_i (b_ij-b_3j) x_i >= 0 for every j}.

If ``B`` has a complex-qubit prepare-and-measure factorisation, an ellipsoid
(an affine image of the Bloch ball) is nested between ``Delta`` and ``Q_B``.
Conversely, a nondegenerate nested ellipsoid reconstructs such a
factorisation.  The fixed-behaviour feasibility problem is an SDP by the
S-lemma.  Its conic alternative gives an explicit obstruction.

Unlike a certificate based on selected pairs of outcomes, the tetrahedron
coordinates keep the inner vertices fixed.  For an archived dual witness,
the only behaviour-dependent checks are the *linear* quantities
``<S_j, R_j(B)>``.  This makes the witness suitable for a later rigorous
superlevel-set cover.

This script is a research certificate generator.  Solver output is numerical
until the archived matrices are rounded outwards or checked with interval
arithmetic.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import cvxpy as cp
import numpy as np


DIMENSION = 3
TRACE_SELECTOR = np.diag([1.0, 1.0, 1.0, 0.0])
TETRAHEDRON_POINTS = (
    np.asarray([0.0, 0.0, 0.0]),
    np.asarray([1.0, 0.0, 0.0]),
    np.asarray([0.0, 1.0, 0.0]),
    np.asarray([0.0, 0.0, 1.0]),
)
TETRAHEDRON_VERTICES = tuple(
    np.outer(np.r_[point, 1.0], np.r_[point, 1.0])
    for point in TETRAHEDRON_POINTS
)


def outer_halfspace_matrix(column: np.ndarray) -> np.ndarray:
    """Return R with [x;1]^T R [x;1] = -b(x).

    The fourth preparation is the barycentric origin.  Thus

        b(x) = b_3 + sum_i (b_i-b_3) x_i,

    and the probability-simplex halfspace is ``-b(x) <= 0``.
    """

    value = np.asarray(column, dtype=float)
    if value.shape != (4,):
        raise ValueError("a behaviour column must have four entries")
    normal = value[3] - value[:3]
    result = np.zeros((4, 4), dtype=float)
    result[:3, 3] = 0.5 * normal
    result[3, :3] = 0.5 * normal
    result[3, 3] = -value[3]
    return result


def validate_behavior(behavior: np.ndarray) -> np.ndarray:
    value = np.asarray(behavior, dtype=float)
    if value.ndim != 2 or value.shape[0] != 4:
        raise ValueError("expected a 4 x m behaviour")
    if np.min(value) < -1e-8:
        raise ValueError("behaviour has a negative entry")
    if np.max(np.abs(value.sum(axis=1) - 1.0)) > 2e-7:
        raise ValueError("behaviour rows are not normalized")
    return value


def active_columns(behavior: np.ndarray, tolerance: float = 1e-13) -> list[int]:
    """Omit identically zero outcomes, which define no outer halfspace."""

    return [
        index
        for index in range(behavior.shape[1])
        if float(np.max(np.abs(behavior[:, index]))) > tolerance
    ]


def solve_primal(
    behavior: np.ndarray, selected_columns: list[int] | None = None
) -> dict[str, Any]:
    """Find an ellipsoid nested between the row tetrahedron and simplex."""

    behavior = validate_behavior(behavior)
    columns = (
        active_columns(behavior)
        if selected_columns is None
        else [int(index) for index in selected_columns]
    )
    ellipse = cp.Variable((4, 4), symmetric=True)
    multipliers = cp.Variable(len(columns), nonneg=True)
    constraints: list[cp.Constraint] = [
        ellipse[:3, :3] >> 0,
        cp.trace(ellipse[:3, :3]) == 1.0,
    ]
    constraints.extend(
        cp.sum(cp.multiply(vertex, ellipse)) <= 0.0
        for vertex in TETRAHEDRON_VERTICES
    )
    halfspaces = [outer_halfspace_matrix(behavior[:, j]) for j in columns]
    constraints.extend(
        ellipse - multipliers[k] * halfspaces[k] >> 0
        for k in range(len(columns))
    )
    problem = cp.Problem(cp.Minimize(cp.norm(ellipse, "fro")), constraints)
    problem.solve(
        solver="CLARABEL",
        tol_gap_abs=1e-10,
        tol_gap_rel=1e-10,
        tol_feas=1e-10,
        max_iter=1000,
    )
    result: dict[str, Any] = {
        "status": problem.status,
        "objective": None if problem.value is None else float(problem.value),
        "active_columns": columns,
    }
    if problem.status not in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}:
        return result
    q = np.asarray(ellipse.value, dtype=float)
    lambdas = np.asarray(multipliers.value, dtype=float)
    result.update(
        {
            "ellipse": q.tolist(),
            "multipliers": lambdas.tolist(),
            "top_block_min_eigenvalue": float(np.linalg.eigvalsh(q[:3, :3]).min()),
            "vertex_quadratic_values": [
                float(np.sum(vertex * q)) for vertex in TETRAHEDRON_VERTICES
            ],
            "containment_min_eigenvalues": [
                float(np.linalg.eigvalsh(q - lambdas[k] * halfspaces[k]).min())
                for k in range(len(columns))
            ],
        }
    )
    return result


def solve_dual(
    behavior: np.ndarray,
    selected_columns: list[int] | None = None,
    robust_budget: float | None = None,
) -> dict[str, Any]:
    """Find a conic Farkas obstruction to a nested Bloch ellipsoid."""

    behavior = validate_behavior(behavior)
    columns = (
        active_columns(behavior)
        if selected_columns is None
        else [int(index) for index in selected_columns]
    )
    halfspaces = [outer_halfspace_matrix(behavior[:, j]) for j in columns]

    state_dual = cp.Variable((3, 3), symmetric=True)
    containment_duals = [cp.Variable((4, 4), symmetric=True) for _ in columns]
    vertex_dual = cp.Variable(4, nonneg=True)
    embedded = cp.bmat(
        [
            [state_dual, np.zeros((3, 1))],
            [np.zeros((1, 3)), np.zeros((1, 1))],
        ]
    )
    stationarity = (
        embedded
        + sum(containment_duals)
        - sum(vertex_dual[z] * TETRAHEDRON_VERTICES[z] for z in range(4))
    )
    constraints: list[cp.Constraint] = [state_dual >> 0]
    constraints.extend(matrix >> 0 for matrix in containment_duals)
    halfspace_expressions = [
        cp.sum(cp.multiply(containment_duals[k], halfspaces[k]))
        for k in range(len(columns))
    ]
    robust_margin = None
    if robust_budget is None:
        constraints.extend(expression >= 0.0 for expression in halfspace_expressions)
    else:
        if robust_budget <= 0.0:
            raise ValueError("robust dual budget must be positive")
        robust_margin = cp.Variable()
        constraints.extend(
            expression >= robust_margin for expression in halfspace_expressions
        )
    constraints.append(stationarity == -TRACE_SELECTOR)
    size_expression = (
        cp.trace(state_dual)
        + sum(cp.trace(matrix) for matrix in containment_duals)
        + cp.sum(vertex_dual)
    )
    if robust_budget is None:
        objective: cp.Minimize | cp.Maximize = cp.Minimize(size_expression)
    else:
        constraints.append(size_expression <= robust_budget)
        objective = cp.Maximize(robust_margin)
    problem = cp.Problem(objective, constraints)
    problem.solve(
        solver="CLARABEL",
        tol_gap_abs=1e-10,
        tol_gap_rel=1e-10,
        tol_feas=1e-10,
        max_iter=1000,
    )
    result: dict[str, Any] = {
        "status": problem.status,
        "objective": None if problem.value is None else float(problem.value),
        "active_columns": columns,
        "robust_budget": robust_budget,
    }
    if problem.status not in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}:
        return result

    state = np.asarray(state_dual.value, dtype=float)
    matrices = [np.asarray(matrix.value, dtype=float) for matrix in containment_duals]
    vertices = np.asarray(vertex_dual.value, dtype=float)
    embedded_value = np.zeros((4, 4), dtype=float)
    embedded_value[:3, :3] = state
    stationarity_value = (
        embedded_value
        + sum(matrices)
        - sum(vertices[z] * TETRAHEDRON_VERTICES[z] for z in range(4))
    )
    margins = [
        float(np.sum(matrices[k] * halfspaces[k]))
        for k in range(len(columns))
    ]
    # <S,R(b)> is linear in the four entries of b.  Store its coefficient
    # vector so later cover code need not reconstruct symbolic matrices.
    linear_coefficients = []
    for matrix in matrices:
        coefficients = []
        for row in range(4):
            basis = np.zeros(4, dtype=float)
            basis[row] = 1.0
            coefficients.append(
                float(np.sum(matrix * outer_halfspace_matrix(basis)))
            )
        linear_coefficients.append(coefficients)

    result.update(
        {
            "state_dual": state.tolist(),
            "containment_duals": [matrix.tolist() for matrix in matrices],
            "vertex_multipliers": vertices.tolist(),
            "halfspace_margins": margins,
            "certified_common_margin": min(margins),
            "optimized_common_margin": (
                None if robust_margin is None else float(robust_margin.value)
            ),
            "halfspace_linear_coefficients": linear_coefficients,
            "state_dual_min_eigenvalue": float(np.linalg.eigvalsh(state).min()),
            "containment_dual_min_eigenvalues": [
                float(np.linalg.eigvalsh(matrix).min()) for matrix in matrices
            ],
            "stationarity_frobenius_residual": float(
                np.linalg.norm(stationarity_value + TRACE_SELECTOR)
            ),
        }
    )
    return result


def physical_audit(trials: int, seed: int, weights: np.ndarray) -> dict[str, Any]:
    from qubit_behavior_lorentz_metric import physical_behavior

    rows = []
    primal_failures = 0
    dual_false_positives = 0
    for index in range(trials):
        behavior, _, _, _ = physical_behavior(seed + index, weights)
        primal = solve_primal(behavior)
        dual = solve_dual(behavior)
        primal_ok = primal["status"] in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}
        dual_ok = dual["status"] in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}
        primal_failures += int(not primal_ok)
        dual_false_positives += int(dual_ok)
        rows.append(
            {
                "trial": index,
                "rank": int(np.linalg.matrix_rank(behavior, tol=1e-10)),
                "primal_status": primal["status"],
                "dual_status": dual["status"],
            }
        )
    return {
        "family": "four random qubit states and a physical joint POVM",
        "trials": trials,
        "seed": seed,
        "weights": weights.tolist(),
        "primal_failures": primal_failures,
        "dual_false_positives": dual_false_positives,
        "rows": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", type=Path)
    parser.add_argument("--audit-physical", type=int, default=0)
    parser.add_argument("--seed", type=int, default=271828)
    parser.add_argument("--weights", type=float, nargs="+", default=[0.99, 0.8, 0.21])
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    weights = np.asarray(args.weights, dtype=float)

    if args.audit_physical:
        payload = physical_audit(args.audit_physical, args.seed, weights)
    else:
        if args.target is None:
            parser.error("--target is required unless --audit-physical is used")
        from qubit_behavior_lorentz_metric import load_moment_behavior

        behavior, loaded_weights = load_moment_behavior(args.target)
        if loaded_weights.size:
            weights = loaded_weights
        payload = {
            "source": str(args.target),
            "shape": list(behavior.shape),
            "rank": int(np.linalg.matrix_rank(behavior, tol=1e-10)),
            "weights": weights.tolist(),
            "primal": solve_primal(behavior),
            "dual": solve_dual(behavior),
        }

    rendered = json.dumps(payload, indent=2) + "\n"
    print(rendered)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

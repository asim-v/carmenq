"""Positive exact RREF recovery for the limiting source-cell certificate."""

from __future__ import annotations

from fractions import Fraction
from typing import Any

from flint import fmpq_mat
import numpy as np
import scipy.linalg
import scipy.sparse as sp
from scipy.optimize import linprog

import exact_socp_dual_certificate as base


EXTRA_PREFIXES = ("soc:762:",)

DISCOVERY_RAY_OBJECTIVE = np.empty(0, dtype=float)
original_base_linprog = base.linprog


def capture_discovery_objective(c: np.ndarray, *args: object, **kwargs: object):
    global DISCOVERY_RAY_OBJECTIVE
    DISCOVERY_RAY_OBJECTIVE = np.asarray(c, dtype=float).copy()
    return original_base_linprog(c, *args, **kwargs)


def positive_exact_recovery(
    matrix: sp.csc_matrix,
    objective: np.ndarray,
    rays: list[base.Ray],
    floating_coefficients: np.ndarray,
    threshold: float,
) -> tuple[list[int], list[Fraction], dict[str, Any]]:
    """Use exact RREF plus a small LP over its rational nullspace."""

    active = [
        index
        for index, (ray, coefficient) in enumerate(
            zip(rays, floating_coefficients, strict=True)
        )
        if ray.free
        or float(coefficient) > threshold
        or ray.label.startswith(EXTRA_PREFIXES)
        or ray.cone == "nonnegative"
        or ray.label.endswith(":clarabel")
    ]
    generators = base.floating_generators(matrix.shape[0], rays)
    active_float = np.asarray((matrix.T @ generators[:, active]).toarray())
    _, triangular, _ = scipy.linalg.qr(active_float, mode="economic", pivoting=True)
    diagonal = np.abs(np.diag(triangular))
    tolerance = (
        max(active_float.shape)
        * np.finfo(float).eps
        * max(float(diagonal[0]) if diagonal.size else 0.0, 1.0)
    )
    numerical_rank = int(np.count_nonzero(diagonal > tolerance))
    _, row_triangular, row_pivots = scipy.linalg.qr(
        active_float.T, mode="economic", pivoting=True
    )
    row_rank = int(np.count_nonzero(np.abs(np.diag(row_triangular)) > tolerance))
    if row_rank != numerical_rank:
        raise RuntimeError("row and column numerical ranks differ")
    pivot_rows = [int(value) for value in row_pivots[:numerical_rank]]

    matrix_csr = matrix.tocsr()
    exact_columns = [
        base.exact_ray_column(matrix_csr, rays[index], matrix.shape[1])
        for index in active
    ]
    augmented_flat: list[int] = []
    for row in pivot_rows:
        coefficients = [column[row] for column in exact_columns]
        target = -base.exact_float(float(objective[row]))
        denominator = max(
            [value.denominator for value in coefficients] + [target.denominator]
        )
        augmented_flat.extend(
            int(value.numerator * (denominator // value.denominator))
            for value in coefficients
        )
        augmented_flat.append(
            int(target.numerator * (denominator // target.denominator))
        )
    reduced, exact_rank = fmpq_mat(
        numerical_rank, len(active) + 1, augmented_flat
    ).rref()
    if exact_rank != numerical_rank:
        raise RuntimeError("exact and numerical ranks differ")

    pivot_columns: list[int] = []
    for row in range(exact_rank):
        pivot = next(
            (column for column in range(len(active)) if reduced[row, column] != 0),
            None,
        )
        if pivot is None:
            raise RuntimeError("RREF has no coefficient pivot")
        pivot_columns.append(pivot)
    pivot_set = set(pivot_columns)
    free_columns = [column for column in range(len(active)) if column not in pivot_set]

    # In x_p = rhs - R_pf f, require every conic x_p and conic free f to be
    # nonnegative while minimizing the same dual upper used in discovery.
    if DISCOVERY_RAY_OBJECTIVE.size != len(rays):
        raise RuntimeError("the discovery LP objective was not captured")
    ray_objectives = [float(DISCOVERY_RAY_OBJECTIVE[index]) for index in active]
    constant_objective = 0.0
    free_objective = np.asarray(
        [ray_objectives[column] for column in free_columns], dtype=float
    )
    inequalities: list[list[float]] = []
    bounds: list[float] = []
    for row, pivot in enumerate(pivot_columns):
        right = float(reduced[row, len(active)])
        coefficients = [float(reduced[row, column]) for column in free_columns]
        constant_objective += ray_objectives[pivot] * right
        free_objective -= ray_objectives[pivot] * np.asarray(coefficients)
        if not rays[active[pivot]].free:
            inequalities.append(coefficients)
            bounds.append(right - 1e-8)
    free_bounds = [
        (None, None) if rays[active[column]].free else (0.0, None)
        for column in free_columns
    ]
    nullspace_lp = linprog(
        free_objective,
        A_ub=np.asarray(inequalities) if inequalities else None,
        b_ub=np.asarray(bounds) if bounds else None,
        bounds=free_bounds,
        method="highs",
        options={
            "dual_feasibility_tolerance": 1e-10,
            "primal_feasibility_tolerance": 1e-10,
        },
    )
    if not nullspace_lp.success:
        raise RuntimeError(f"positive nullspace LP failed: {nullspace_lp.message}")
    denominator = 1 << 60
    values = [Fraction(0) for _ in active]
    for column, approximate in zip(free_columns, nullspace_lp.x, strict=True):
        values[column] = Fraction(
            int(round(float(approximate) * denominator)), denominator
        )
    for row, pivot in reversed(list(enumerate(pivot_columns))):
        right = reduced[row, len(active)]
        value = Fraction(int(right.p), int(right.q))
        for column in free_columns:
            coefficient = reduced[row, column]
            if coefficient:
                value -= (
                    Fraction(int(coefficient.p), int(coefficient.q)) * values[column]
                )
        values[pivot] = value
    negatives = [
        (values[position], rays[index].label)
        for position, index in enumerate(active)
        if not rays[index].free and values[position] < 0
    ]
    if negatives:
        worst = min(negatives)
        raise RuntimeError(
            f"exact nullspace rounding left {len(negatives)} negatives; worst {float(worst[0])} {worst[1]}"
        )
    return (
        active,
        values,
        {
            "candidate_columns": len(active),
            "active_columns": len(active),
            "rank": int(exact_rank),
            "free_parameters": len(free_columns),
            "basis_recovery": "FLINT RREF plus positive nullspace LP",
            "nullspace_lp_objective": float(constant_objective + nullspace_lp.fun),
        },
    )


# The base main assigns canonical right-hand-side data immediately before it
# calls recovery.  Wrap that call to make the data available to the exact
# ray-objective computation above without changing the stable base exporter.
DATA_RIGHT: np.ndarray
original_audit = base.exact_certificate_audit


def capture_right(
    matrix: sp.csc_matrix,
    right: np.ndarray,
    objective: np.ndarray,
    rays: list[base.Ray],
    active: list[int],
    coefficients: list[Fraction],
    target: Fraction,
):
    global DATA_RIGHT
    DATA_RIGHT = right
    return original_audit(matrix, right, objective, rays, active, coefficients, target)


# Recovery runs before the audit, so install a lightweight wrapper around the
# main recovery that receives the right-hand side through a temporary hook in
# ``exact_ray_objective``.
original_exact_ray_objective = base.exact_ray_objective


def remember_right(right: np.ndarray, ray: base.Ray) -> Fraction:
    global DATA_RIGHT
    DATA_RIGHT = right
    return original_exact_ray_objective(right, ray)


base.linprog = capture_discovery_objective
base.exact_ray_objective = remember_right
base.rationalize_active_solution = positive_exact_recovery
base.exact_certificate_audit = capture_right

if __name__ == "__main__":
    base.main()

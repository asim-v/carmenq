"""Dependency-light Helstrom geometry for subnormalised qubit ensembles."""

from __future__ import annotations

import cvxpy as cp
import numpy as np


PAULIS = (
    np.array([[0, 1], [1, 0]], dtype=complex),
    np.array([[0, -1j], [1j, 0]], dtype=complex),
    np.diag([1, -1]).astype(complex),
)


def discrimination_geometry(states: np.ndarray) -> dict[str, object]:
    """Solve the weighted smallest-enclosing-ball dual for qubit states."""

    array = np.asarray(states, dtype=complex)
    if array.shape != (4, 2, 2):
        raise ValueError("states must have shape (4, 2, 2)")
    priors = np.trace(array, axis1=1, axis2=2).real
    vectors = np.asarray(
        [[np.trace(state @ pauli).real for pauli in PAULIS] for state in array]
    )
    center = cp.Variable(3)
    radius = cp.Variable()
    constraints = [
        cp.norm(center - vectors[index], 2) <= radius - priors[index]
        for index in range(4)
    ]
    problem = cp.Problem(cp.Minimize(radius), constraints)
    problem.solve(solver="CLARABEL", tol_gap_abs=1e-11, tol_feas=1e-11)
    if problem.status not in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}:
        raise RuntimeError(f"terminal discrimination SOCP failed: {problem.status}")

    center_value = np.asarray(center.value, dtype=float)
    radius_value = float(radius.value)
    slacks = np.asarray(
        [
            radius_value
            - priors[index]
            - np.linalg.norm(center_value - vectors[index])
            for index in range(4)
        ]
    )
    weights = np.asarray(
        [float(constraint.dual_value) for constraint in constraints]
    )
    directions = []
    for index in range(4):
        displacement = vectors[index] - center_value
        norm = float(np.linalg.norm(displacement))
        directions.append(
            (displacement / norm).tolist()
            if norm > 1e-14
            else [0.0, 0.0, 0.0]
        )
    active_tolerance = 2e-8
    return {
        "optimal_guess_probability": radius_value,
        "dual_center": center_value.tolist(),
        "constraint_slacks": slacks.tolist(),
        "active_indices": np.flatnonzero(slacks <= active_tolerance).tolist(),
        "kkt_weights": weights.tolist(),
        "optimal_effect_traces": (2.0 * weights).tolist(),
        "optimal_effect_bloch": directions,
    }


__all__ = ["PAULIS", "discrimination_geometry"]

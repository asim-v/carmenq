"""Certify full-rank qubit prepare-and-measure behaviours with a Lorentz metric.

Let ``B`` be a row-stochastic matrix with four rows and full row rank.  A
qubit factorisation has the form

    B = X E,  X[z] = (1, n_z),  E[:, j] = (c_j, t_j),

where ``||n_z|| <= 1``, ``c_j >= ||t_j||``, and the POVM columns sum to
``(1, 0, 0, 0)``.  If ``A = X^{-1}`` and

    H = A.T @ diag(1, -1, -1, -1) @ A,

then effect positivity is linear in ``H`` after evaluating the fixed
quadratic forms ``b_j.T @ H @ b_j``.  State positivity is equivalent to the
four principal 3-by-3 submatrices of ``H`` being negative semidefinite,
provided ``H`` is nonsingular with Lorentz inertia.  Strict feasibility of
the SDP below therefore gives an explicit qubit factorisation; infeasibility
is a rigorous obstruction up to the chosen solver tolerances.

For the fixed terminal weights used in the frontier calculation, the same
cone test is also applied to every residual effect
``w_t Q_y - G[y,t]``.

This is a research/audit script, not yet part of the public package API.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import cvxpy as cp
import numpy as np

from joint_effect_dimension_seesaw import random_point


ETA = np.diag([1.0, -1.0, -1.0, -1.0])
PAULI = np.asarray(
    [
        [[0.0, 1.0], [1.0, 0.0]],
        [[0.0, -1.0j], [1.0j, 0.0]],
        [[1.0, 0.0], [0.0, -1.0]],
    ],
    dtype=complex,
)


def _linear_quadratic_form(h: cp.Expression, vector: np.ndarray) -> cp.Expression:
    """Return v.T H v as an affine expression when v is numerical."""

    return cp.sum(cp.multiply(h, np.outer(vector, vector)))


def load_moment_behavior(path: Path) -> tuple[np.ndarray, np.ndarray]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if isinstance(payload, list):
        if not payload:
            raise ValueError("behaviour result list is empty")
        payload = payload[0]
    joint = np.asarray(payload["path_terminal_statistics"], dtype=float)
    priors = joint.sum(axis=(1, 2))
    behavior = joint.reshape(4, -1) / priors[:, None]
    weights = np.asarray(payload["terminal_effect_weights"][: joint.shape[2]])
    return behavior, weights


def physical_behavior(
    seed: int, weights: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Generate an exactly physical behaviour and its Bloch factorisation."""

    padded = np.zeros(4, dtype=float)
    padded[: len(weights)] = weights
    states, joint = random_point(seed, padded)
    priors = np.trace(states, axis1=1, axis2=2).real
    normalized_states = states / priors[:, None, None]
    effects = joint[:, : len(weights)].reshape(-1, 2, 2)

    bloch_states = np.asarray(
        [
            [np.trace(rho @ sigma).real for sigma in PAULI]
            for rho in normalized_states
        ]
    )
    effect_scalar = np.trace(effects, axis1=1, axis2=2).real / 2.0
    effect_vector = np.asarray(
        [[np.trace(effect @ sigma).real / 2.0 for sigma in PAULI] for effect in effects]
    )
    x = np.column_stack([np.ones(4), bloch_states])
    e = np.row_stack([effect_scalar, effect_vector.T])
    behavior = x @ e
    return behavior, x, e, joint


def metric_from_factorisation(x: np.ndarray) -> np.ndarray:
    inverse = np.linalg.inv(x)
    return inverse.T @ ETA @ inverse


def domination_columns(
    behavior: np.ndarray, weights: np.ndarray
) -> list[tuple[str, np.ndarray]]:
    arity = len(weights)
    columns: list[tuple[str, np.ndarray]] = []
    for y in range(4):
        indices = [arity * y + t for t in range(arity)]
        coarse = behavior[:, indices].sum(axis=1)
        for t, index in enumerate(indices):
            columns.append((f"d_{y}_{t}", weights[t] * coarse - behavior[:, index]))
    return columns


def metric_diagnostics(
    behavior: np.ndarray, weights: np.ndarray, h: np.ndarray
) -> dict[str, Any]:
    one = np.ones(4)
    columns = [(f"b_{j}", behavior[:, j]) for j in range(behavior.shape[1])]
    columns += domination_columns(behavior, weights)
    cone = [float(vector @ h @ vector) for _, vector in columns]
    orientation = [float(one @ h @ vector) for _, vector in columns]
    principal_max = []
    for z in range(4):
        keep = [index for index in range(4) if index != z]
        principal_max.append(float(np.linalg.eigvalsh(h[np.ix_(keep, keep)]).max()))
    eigenvalues = np.linalg.eigvalsh(h)
    inverse_diagonal = (
        np.diag(np.linalg.inv(h)).tolist()
        if abs(float(np.linalg.det(h))) > 1e-12
        else None
    )
    return {
        "normalisation": float(one @ h @ one),
        "eigenvalues": eigenvalues.tolist(),
        "determinant": float(np.linalg.det(h)),
        "inverse_diagonal": inverse_diagonal,
        "largest_principal_eigenvalues": principal_max,
        "minimum_effect_or_residual_quadratic": min(cone),
        "minimum_effect_or_residual_orientation": min(orientation),
        "maximum_row_sum_error": float(np.max(np.abs(behavior.sum(axis=1) - 1.0))),
        "behavior_rank": int(np.linalg.matrix_rank(behavior, tol=1e-10)),
    }


def solve_metric_sdp(
    behavior: np.ndarray,
    weights: np.ndarray,
    cone_tolerance: float = 0.0,
    principal_margin: float = 0.0,
) -> tuple[str, np.ndarray | None, dict[str, Any]]:
    """Solve the fixed-behaviour Lorentz-metric feasibility problem."""

    h = cp.Variable((4, 4), symmetric=True)
    one = np.ones(4)
    constraints: list[cp.Constraint] = [
        _linear_quadratic_form(h, one) == 1.0
    ]
    all_columns = [behavior[:, j] for j in range(behavior.shape[1])]
    all_columns += [vector for _, vector in domination_columns(behavior, weights)]
    for vector in all_columns:
        constraints.append(_linear_quadratic_form(h, vector) >= -cone_tolerance)
        constraints.append(one @ h @ vector >= -cone_tolerance)
    for z in range(4):
        keep = [index for index in range(4) if index != z]
        constraints.append(
            h[np.ix_(keep, keep)] << -principal_margin * np.eye(3)
        )

    # The trace objective selects a bounded representative on degenerate
    # faces and makes solver comparisons reproducible.  It does not change
    # the feasible set.
    problem = cp.Problem(cp.Minimize(cp.norm(h, "fro")), constraints)
    problem.solve(
        solver="CLARABEL",
        tol_gap_abs=1e-10,
        tol_gap_rel=1e-10,
        tol_feas=1e-10,
        max_iter=500,
    )
    solution = None if h.value is None else np.asarray(h.value)
    report: dict[str, Any] = {
        "status": problem.status,
        "objective": None if problem.value is None else float(problem.value),
        "cone_tolerance": cone_tolerance,
        "principal_margin": principal_margin,
    }
    if solution is not None:
        report["diagnostics"] = metric_diagnostics(behavior, weights, solution)
    return problem.status, solution, report


def reconstruct_factorisation(
    behavior: np.ndarray, h: np.ndarray
) -> tuple[np.ndarray, np.ndarray, float]:
    """Construct X,E from a nonsingular Lorentz metric H.

    An eigendecomposition gives one congruence factor A0 with
    ``H=A0.T ETA A0``.  A Lorentz boost is then obtained numerically by a
    second congruence normalisation constrained by ``A @ 1 = e0``.  The
    simpler closed construction below uses the H-orthogonal decomposition:
    its first row is ``1.T H`` and the remaining rows form a Euclidean
    square root of ``(H 1)(H 1).T - H``.
    """

    one = np.ones(4)
    temporal = h @ one
    spatial_gram = np.outer(temporal, temporal) - h
    values, vectors = np.linalg.eigh((spatial_gram + spatial_gram.T) / 2.0)
    positive = np.where(values > 1e-9)[0]
    if len(positive) != 3:
        raise ValueError(f"expected spatial rank three, got eigenvalues {values}")
    spatial = (np.sqrt(values[positive])[:, None] * vectors[:, positive].T)
    a = np.row_stack([temporal, spatial])
    if np.linalg.det(a) < 0.0:
        a[-1] *= -1.0
    x = np.linalg.inv(a)
    e = a @ behavior
    residual = float(np.linalg.norm(x @ e - behavior))
    return x, e, residual


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", type=Path)
    parser.add_argument("--physical-seed", type=int, default=123)
    parser.add_argument("--weights", type=float, nargs="+", default=[0.99, 0.8, 0.21])
    parser.add_argument("--cone-tolerance", type=float, default=0.0)
    parser.add_argument("--principal-margin", type=float, default=0.0)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    weights = np.asarray(args.weights, dtype=float)

    if args.target is None:
        behavior, known_x, known_e, _ = physical_behavior(args.physical_seed, weights)
        known_h = metric_from_factorisation(known_x)
        source: dict[str, Any] = {
            "kind": "exact_physical_random",
            "seed": args.physical_seed,
            "factorisation_residual": float(np.linalg.norm(known_x @ known_e - behavior)),
            "known_metric": metric_diagnostics(behavior, weights, known_h),
        }
    else:
        behavior, loaded_weights = load_moment_behavior(args.target)
        if not np.allclose(weights, loaded_weights):
            weights = loaded_weights
        source = {"kind": "json", "path": str(args.target)}

    status, h, result = solve_metric_sdp(
        behavior,
        weights,
        cone_tolerance=args.cone_tolerance,
        principal_margin=args.principal_margin,
    )
    payload: dict[str, Any] = {
        "source": source,
        "weights": weights.tolist(),
        "sdp": result,
    }
    if status in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE} and h is not None:
        try:
            x, e, residual = reconstruct_factorisation(behavior, h)
            payload["reconstruction"] = {
                "residual": residual,
                "state_lorentz_values": np.diag(x @ ETA @ x.T).tolist(),
                "effect_lorentz_minimum": float(
                    np.min(np.diag(e.T @ ETA @ e))
                ),
                "effect_time_minimum": float(np.min(e[0])),
                "completeness_residual": float(
                    np.linalg.norm(e.sum(axis=1) - np.asarray([1.0, 0.0, 0.0, 0.0]))
                ),
            }
        except (ValueError, np.linalg.LinAlgError) as error:
            payload["reconstruction_error"] = str(error)

    rendered = json.dumps(payload, indent=2) + "\n"
    print(rendered)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

"""Global-search probe of the exact four-effect prior-reserve geometry.

This script strengthens :mod:`active_readout_audit_cap` diagnostically by
enforcing the full Bloch-vector closure rather than only its projection on the
Helstrom-dual axis.  Four weighted unit Bloch vectors form a closed spatial
quadrilateral.  Pairing sides ``(0,1)`` and ``(2,3)`` gives an exact
five-parameter representation: the common diagonal length, the direction of
the dual axis relative to that diagonal, two pair-plane azimuths, and the
dual spectral bias.

Differential evolution plus local polishing searches this compact domain.
The returned minimum is evidence, not an upper certificate: proving a cap
requires a lower bound on the global minimum, which will be supplied by a
subsequent interval cover.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import differential_evolution

from active_readout_audit_cap import reserve


def diagonal_interval(weights: np.ndarray) -> tuple[float, float]:
    """Return the possible length interval of ``w0*n0+w1*n1``."""

    value = np.asarray(weights, dtype=float)
    lower = max(abs(value[0] - value[1]), abs(value[2] - value[3]))
    upper = min(value[0] + value[1], value[2] + value[3])
    if lower > upper + 1e-12:
        raise ValueError("the four side lengths cannot form a closed polygon")
    return max(0.0, lower), max(0.0, upper)


def pair_components(first: float, second: float, diagonal: float) -> tuple[float, float, float]:
    """Return the two axial components and common transverse magnitude."""

    if diagonal <= 1e-12:
        if abs(first - second) > 1e-9:
            raise ValueError("a zero diagonal requires equal paired weights")
        return 0.0, 0.0, first
    axial_first = (diagonal * diagonal + first * first - second * second) / (
        2.0 * diagonal
    )
    axial_second = diagonal - axial_first
    transverse = np.sqrt(max(0.0, first * first - axial_first * axial_first))
    return axial_first, axial_second, transverse


def projections(weights: np.ndarray, variables: np.ndarray) -> np.ndarray:
    """Map exact quadrilateral parameters to four dual-axis projections."""

    diagonal, dual_cosine, first_azimuth, second_azimuth, _ = variables
    dual_sine = np.sqrt(max(0.0, 1.0 - dual_cosine * dual_cosine))
    a0, a1, h01 = pair_components(weights[0], weights[1], diagonal)
    a2, a3, h23 = pair_components(weights[2], weights[3], diagonal)
    transverse_dual = dual_sine
    return np.asarray(
        [
            (a0 * dual_cosine + h01 * transverse_dual * np.cos(first_azimuth))
            / weights[0],
            (a1 * dual_cosine - h01 * transverse_dual * np.cos(first_azimuth))
            / weights[1],
            (-a2 * dual_cosine + h23 * transverse_dual * np.cos(second_azimuth))
            / weights[2],
            (-a3 * dual_cosine - h23 * transverse_dual * np.cos(second_azimuth))
            / weights[3],
        ]
    )


def bloch_vectors(weights: np.ndarray, variables: np.ndarray) -> np.ndarray:
    """Return the four unit vectors in the closed-quadrilateral chart.

    The common diagonal is the first Cartesian axis.  The dual axis can be
    chosen in the first/second coordinate plane, so the two azimuths retain
    the complete relative three-dimensional geometry.
    """

    value = np.asarray(weights, dtype=float)
    diagonal, _, first_azimuth, second_azimuth, _ = variables
    a0, a1, h01 = pair_components(value[0], value[1], diagonal)
    a2, a3, h23 = pair_components(value[2], value[3], diagonal)
    weighted = np.asarray(
        [
            [a0, h01 * np.cos(first_azimuth), h01 * np.sin(first_azimuth)],
            [a1, -h01 * np.cos(first_azimuth), -h01 * np.sin(first_azimuth)],
            [-a2, h23 * np.cos(second_azimuth), h23 * np.sin(second_azimuth)],
            [-a3, -h23 * np.cos(second_azimuth), -h23 * np.sin(second_azimuth)],
        ]
    )
    return weighted / value[:, None]


def aligned_projective_audit(
    audit: float,
    complement_prior: float,
    bloch_overlap: float,
) -> float:
    """Exact audit after retaining one active projector and its complement.

    If every Helstrom constraint is active, complementarity gives

        rho_j = Y - (A-p_j) Pi_j^perp.

    Retaining ``Pi_i`` for answer ``i`` and ``I-Pi_i`` for answer ``j`` then
    has the audit below, where ``bloch_overlap = n_i dot n_j``.  The formula
    is independent of the spectral bias and orientation of ``Y``.
    """

    return 0.5 * (
        (1.0 - bloch_overlap) * audit
        + (1.0 + bloch_overlap) * complement_prior
    )


def averaged_projective_audit(
    audit: float,
    complement_prior: float,
    complement_weight: float,
) -> float:
    """Weighted average over every retained label other than the complement."""

    if not 0.0 <= complement_weight <= 1.0:
        raise ValueError("a qubit rank-one effect weight lies in [0,1]")
    loss_factor = (1.0 - complement_weight) / (2.0 - complement_weight)
    return audit - loss_factor * (audit - complement_prior)


def objective(weights: np.ndarray, variables: np.ndarray) -> float:
    point = projections(weights, variables)
    return float(np.sum(reserve(float(variables[4]), point)))


def geometry_probe(
    weights: np.ndarray,
    seed: int = 0,
    maxiter: int = 1200,
    popsize: int = 24,
) -> dict[str, Any]:
    value = np.sort(np.asarray(weights, dtype=float))[::-1]
    if value.shape != (4,) or np.any(value <= 0.0) or np.any(value > 1.0 + 1e-12):
        raise ValueError("expected four weights in (0,1]")
    if abs(float(value.sum()) - 2.0) > 1e-9:
        raise ValueError("qubit effect weights must sum to two")
    lower, upper = diagonal_interval(value)
    epsilon = 1e-10 if lower == 0.0 else 0.0
    result = differential_evolution(
        lambda variables: objective(value, variables),
        bounds=(
            (lower + epsilon, upper),
            (-1.0, 1.0),
            (0.0, 2.0 * np.pi),
            (0.0, 2.0 * np.pi),
            (0.0, 1.0 - 1e-10),
        ),
        seed=seed,
        maxiter=maxiter,
        popsize=popsize,
        tol=1e-11,
        polish=True,
        updating="immediate",
        workers=1,
    )
    variables = np.asarray(result.x, dtype=float)
    point = projections(value, variables)
    bloch_closure = float(np.dot(value, point))
    return {
        "weights": value.tolist(),
        "minimum_prior_reserve_found": float(result.fun),
        "candidate_audit_cap": float(1.0 / result.fun),
        "diagonal_interval": [lower, upper],
        "diagonal": float(variables[0]),
        "dual_axis_cosine": float(variables[1]),
        "pair_azimuths": variables[2:4].tolist(),
        "dual_spectral_bias": float(variables[4]),
        "bloch_projections": point.tolist(),
        "projected_closure_residual": bloch_closure,
        "optimizer_success": bool(result.success),
        "optimizer_message": str(result.message),
        "function_evaluations": int(result.nfev),
        "status": "global-search diagnostic; not a certified upper bound",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--weights", type=float, nargs=4, required=True)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--maxiter", type=int, default=1200)
    parser.add_argument("--popsize", type=int, default=24)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = geometry_probe(
        np.asarray(args.weights, dtype=float),
        args.seed,
        args.maxiter,
        args.popsize,
    )
    rendered = json.dumps(payload, indent=2) + "\n"
    print(rendered, end="")
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

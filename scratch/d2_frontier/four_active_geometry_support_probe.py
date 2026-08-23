"""Probe the joint AUDIT--RETURN cap of a fully active four-outcome readout.

For a fixed four-effect Bloch geometry and Helstrom-dual bias, positivity of
the four syndrome states gives componentwise lower bounds

    p_s >= A f_t(x_s),

where ``A`` is the AUDIT success.  After forgetting every other causal and
common-instrument constraint, the Hellinger RETURN obeys

    R <= (sum_s sqrt(p_s))**2 / 4.

For fixed ``A`` and reserve vector ``f``, the right-hand side is maximised by
water filling the residual syndrome mass as uniformly as the lower bounds
allow.  This script combines that exact scalar solution with the closed-
quadrilateral parameterisation from :mod:`active_readout_geometry_probe` and
performs a global numerical search.

The search is a diagnostic lower estimate of the maximum of a valid outer
relaxation.  It is not a proof of an upper bound until a branch-and-bound or
interval cover certifies the whole compact parameter domain.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import differential_evolution, minimize_scalar

from active_readout_audit_cap import reserve
from active_readout_geometry_probe import diagonal_interval, projections


def water_filled_priors(lower: np.ndarray) -> np.ndarray:
    """Maximise Hellinger symmetry subject to ``p>=lower, sum(p)=1``."""

    floor = np.asarray(lower, dtype=float)
    if floor.sum() > 1.0 + 1e-10:
        raise ValueError("infeasible syndrome lower bounds")
    ordered = np.sort(floor)[::-1]
    level = 0.0
    for fixed in range(len(ordered) + 1):
        remaining = 1.0 - float(ordered[:fixed].sum())
        free = len(ordered) - fixed
        if free == 0:
            level = 0.0
            break
        candidate = remaining / free
        if fixed == len(ordered) or candidate >= ordered[fixed] - 1e-14:
            level = candidate
            break
    return np.maximum(floor, level)


def support_for_reserve(
    reserve_vector: np.ndarray,
    support_weight: float,
) -> tuple[float, float, np.ndarray, float]:
    """Optimise AUDIT ``A`` and syndrome priors for one reserve vector."""

    value = np.asarray(reserve_vector, dtype=float)
    maximum_audit = min(1.0, 1.0 / float(value.sum()))

    def negative_score(audit: float) -> float:
        priors = water_filled_priors(audit * value)
        returned = float(np.square(np.sqrt(priors).sum()) / 4.0)
        return -(support_weight * audit + (1.0 - support_weight) * returned)

    candidates = [(negative_score(0.0), 0.0), (negative_score(maximum_audit), maximum_audit)]
    result = minimize_scalar(
        negative_score,
        bounds=(0.0, maximum_audit),
        method="bounded",
        options={"xatol": 1e-13},
    )
    candidates.append((float(result.fun), float(result.x)))
    negative, audit = min(candidates)
    priors = water_filled_priors(audit * value)
    returned = float(np.square(np.sqrt(priors).sum()) / 4.0)
    return -negative, audit, priors, returned


def support_for_reserve_fast(
    reserve_vector: np.ndarray,
    support_weight: float,
) -> tuple[float, float, np.ndarray, float]:
    """Exact finite water-filling solution without scalar optimisation.

    On a region where the ``k`` largest prior floors are active, put
    ``F=sum(f[:k])``, ``C=sum(sqrt(f[:k]))``, and ``m=4-k``.  With
    ``y=sqrt(F*A)`` the support is a two-dimensional Rayleigh quotient in
    ``(y,sqrt(1-y**2))``.  Its Perron vector gives the interior maximiser;
    the two water-filling breakpoints cover the constrained endpoints.
    """

    values = np.sort(np.asarray(reserve_vector, dtype=float))[::-1]
    maximum_audit = min(1.0, 1.0 / float(values.sum()))

    def evaluate(audit: float) -> tuple[float, float, np.ndarray, float]:
        priors = water_filled_priors(audit * values)
        returned = float(np.square(np.sqrt(priors).sum()) / 4.0)
        score = support_weight * audit + (1.0 - support_weight) * returned
        return score, audit, priors, returned

    candidates = [evaluate(0.0), evaluate(maximum_audit)]
    if values[0] > 0.0:
        candidates.append(evaluate(min(maximum_audit, 1.0 / (4.0 * values[0]))))
    for active in range(1, 4):
        fixed_sum = float(values[:active].sum())
        root_sum = float(np.sqrt(values[:active]).sum())
        free = 4 - active
        lower = 1.0 / (fixed_sum + free * values[active - 1])
        upper = 1.0 / (fixed_sum + free * values[active])
        lower = max(0.0, lower)
        upper = min(maximum_audit, upper)
        if lower > upper + 2e-15:
            continue
        candidates.extend((evaluate(lower), evaluate(upper)))

        matrix = np.asarray(
            [
                [
                    support_weight / fixed_sum
                    + (1.0 - support_weight)
                    * root_sum**2
                    / (4.0 * fixed_sum),
                    (1.0 - support_weight)
                    * root_sum
                    * np.sqrt(free)
                    / (4.0 * np.sqrt(fixed_sum)),
                ],
                [
                    (1.0 - support_weight)
                    * root_sum
                    * np.sqrt(free)
                    / (4.0 * np.sqrt(fixed_sum)),
                    (1.0 - support_weight) * free / 4.0,
                ],
            ]
        )
        _, vectors = np.linalg.eigh(matrix)
        perron = np.abs(vectors[:, -1])
        perron /= np.linalg.norm(perron)
        stationary_audit = float(perron[0] ** 2 / fixed_sum)
        if lower - 2e-15 <= stationary_audit <= upper + 2e-15:
            candidates.append(evaluate(np.clip(stationary_audit, lower, upper)))
    return max(candidates, key=lambda item: item[0])


def geometry_support_probe(
    weights: np.ndarray,
    support_weight: float = 0.6,
    seed: int = 0,
    maxiter: int = 1500,
    popsize: int = 28,
) -> dict[str, Any]:
    value = np.sort(np.asarray(weights, dtype=float))[::-1]
    if value.shape != (4,) or np.any(value <= 0.0) or np.any(value > 1.0 + 1e-12):
        raise ValueError("expected four weights in (0,1]")
    if abs(float(value.sum()) - 2.0) > 1e-9:
        raise ValueError("qubit effect weights must sum to two")
    lower, upper = diagonal_interval(value)
    epsilon = 1e-10 if lower == 0.0 else 0.0

    def negative_outer_support(variables: np.ndarray) -> float:
        point = projections(value, variables)
        reserve_vector = np.asarray(reserve(float(variables[4]), point), dtype=float)
        score, _, _, _ = support_for_reserve(reserve_vector, support_weight)
        return -score

    result = differential_evolution(
        negative_outer_support,
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
        tol=1e-10,
        polish=True,
        updating="immediate",
        workers=1,
    )
    variables = np.asarray(result.x, dtype=float)
    point = projections(value, variables)
    reserve_vector = np.asarray(reserve(float(variables[4]), point), dtype=float)
    score, audit, priors, returned = support_for_reserve(
        reserve_vector, support_weight
    )
    return {
        "weights": value.tolist(),
        "support_weight": float(support_weight),
        "outer_support_found": float(score),
        "audit": float(audit),
        "return": float(returned),
        "syndrome_priors": priors.tolist(),
        "prior_reserves": reserve_vector.tolist(),
        "diagonal_interval": [lower, upper],
        "diagonal": float(variables[0]),
        "dual_axis_cosine": float(variables[1]),
        "pair_azimuths": variables[2:4].tolist(),
        "dual_spectral_bias": float(variables[4]),
        "bloch_projections": point.tolist(),
        "projected_closure_residual": float(np.dot(value, point)),
        "optimizer_success": bool(result.success),
        "optimizer_message": str(result.message),
        "function_evaluations": int(result.nfev),
        "status": "global-search diagnostic of an outer model; not an upper certificate",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--weights", type=float, nargs=4, required=True)
    parser.add_argument("--lambda", dest="support_weight", type=float, default=0.6)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--maxiter", type=int, default=1500)
    parser.add_argument("--popsize", type=int, default=28)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = geometry_support_probe(
        np.asarray(args.weights, dtype=float),
        args.support_weight,
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

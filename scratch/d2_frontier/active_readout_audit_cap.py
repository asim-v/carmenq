"""Helstrom-complementarity audit cap for a fully active qubit POVM.

Suppose an optimal minimum-error readout has ``k`` nonzero rank-one effects

    E_s = w_s |n_s><n_s|,   sum_s E_s = I_2.

Let ``Y`` be the Helstrom dual, ``A = Tr(Y)``, and write the spectral bias of
``Y`` as ``t in [0,1]``.  Complementary slackness and positivity of every
state imply

    p_s >= A f_t(x_s),
    f_t(x) = (1 + 2 t x + t^2) / (2 (1 + t x)),

where ``x_s`` is the projection of the effect Bloch vector onto the Bloch
axis of ``Y``.  POVM completeness gives ``sum_s w_s x_s = 0``.  Therefore

    A <= 1 / min sum_s f_t(x_s).

For fixed ``t``, ``f_t`` is increasing and concave in ``x``.  Its minimum on
the box cut by one linear equality is attained at a vertex, so all but one
of the ``x_s`` are signs.  The remaining scalar minimisation over ``t`` is
performed separately on every finite vertex branch.

The projection constraints are necessary but not sufficient for a set of
Bloch vectors.  Consequently the returned value is a valid (occasionally
loose) audit upper bound for every fully active rank-one qubit readout.
"""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import minimize_scalar


def reserve(t: float, x: np.ndarray | float) -> np.ndarray | float:
    """Return the dimensionless prior reserve ``f_t(x)`` stably."""

    t_value = float(t)
    point = np.asarray(x, dtype=float)
    denominator = 1.0 + t_value * point
    numerator = 1.0 + 2.0 * t_value * point + t_value * t_value
    result = np.empty_like(point)
    regular = denominator > 1e-13
    result[regular] = numerator[regular] / (2.0 * denominator[regular])
    # The only singular endpoint is t=1,x=-1.  There the support expectation
    # of Y vanishes and positivity requires no prior reserve.
    result[~regular] = 0.0
    return float(result) if result.ndim == 0 else result


def projection_vertices(weights: np.ndarray) -> list[np.ndarray]:
    """Enumerate vertices of ``[-1,1]^k intersect {w.x=0}``."""

    value = np.asarray(weights, dtype=float)
    vertices: list[np.ndarray] = []
    seen: set[tuple[int, ...]] = set()
    for free in range(len(value)):
        fixed = [index for index in range(len(value)) if index != free]
        for signs in itertools.product((-1.0, 1.0), repeat=len(fixed)):
            point = np.zeros(len(value), dtype=float)
            point[fixed] = signs
            point[free] = -float(np.dot(value[fixed], point[fixed])) / value[free]
            if point[free] < -1.0 - 1e-12 or point[free] > 1.0 + 1e-12:
                continue
            point[free] = float(np.clip(point[free], -1.0, 1.0))
            key = tuple(int(round(item * 1e12)) for item in point)
            if key not in seen:
                seen.add(key)
                vertices.append(point)
    if not vertices:
        raise ValueError("the projection polytope has no vertices")
    return vertices


def active_audit_cap(weights: np.ndarray) -> dict[str, Any]:
    """Return the fully-active audit cap and its relaxed minimizer."""

    value = np.asarray(weights, dtype=float)
    if value.ndim != 1 or len(value) not in {3, 4}:
        raise ValueError("expected three or four active effect weights")
    if np.any(value <= 0.0) or np.any(value > 1.0 + 1e-12):
        raise ValueError("active rank-one effect weights must lie in (0,1]")
    if abs(float(value.sum()) - 2.0) > 1e-9:
        raise ValueError("qubit rank-one effect weights must sum to two")

    best = (np.inf, 0.0, np.zeros(len(value)))
    for point in projection_vertices(value):
        objective = lambda t: float(np.sum(reserve(t, point)))
        candidates = [(objective(0.0), 0.0), (objective(1.0), 1.0)]
        result = minimize_scalar(
            objective,
            bounds=(0.0, 1.0 - 1e-12),
            method="bounded",
            options={"xatol": 2e-14},
        )
        candidates.append((float(result.fun), float(result.x)))
        minimum, t_value = min(candidates)
        if minimum < best[0]:
            best = (minimum, t_value, point.copy())
    return {
        "weights": value.tolist(),
        "minimum_total_prior_reserve": float(best[0]),
        "audit_upper": float(min(1.0, 1.0 / best[0])),
        "dual_spectral_bias": float(best[1]),
        "relaxed_bloch_projections": best[2].tolist(),
        "projection_residual": float(np.dot(value, best[2])),
        "vertex_count": len(projection_vertices(value)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--weights", type=float, nargs="+", required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = active_audit_cap(np.asarray(args.weights, dtype=float))
    rendered = json.dumps(payload, indent=2) + "\n"
    print(rendered, end="")
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

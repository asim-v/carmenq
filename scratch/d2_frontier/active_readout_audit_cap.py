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


def closed_four_active_audit_cap(weights: np.ndarray) -> dict[str, Any]:
    """Return the exact closed-form cap for four sorted active weights.

    Put ``u_i=(x_i+1)/2``.  The projection constraint becomes
    ``sum_i w_i u_i=1``.  For every fixed dual bias, the reserve as a
    function of ``u`` is increasing and concave.  A minimum is therefore a
    fractional-knapsack vertex.  Since ``w0 >= w1 >= w2 >= w3``, the greedy
    vertex fills ``u0`` and then ``u1``:

        u = (1, (1-w0)/w1, 0, 0).

    Every vertex with at least two filled coordinates costs at least two
    endpoint increments; among one-fill vertices the displayed ratio is
    smallest.  Minimising the remaining scalar reserve gives

        s = sqrt(2(1-w0)/w1),  t = 1/(1+s),
        m = (5 - 2/(1+s) - 1/(1+s)^2)/2,  A <= 1/m.
    """

    value = np.sort(np.asarray(weights, dtype=float))[::-1]
    if value.shape != (4,):
        raise ValueError("expected four active effect weights")
    if np.any(value <= 0.0) or np.any(value > 1.0 + 1e-12):
        raise ValueError("active rank-one effect weights must lie in (0,1]")
    if abs(float(value.sum()) - 2.0) > 1e-9:
        raise ValueError("qubit rank-one effect weights must sum to two")
    fraction = (1.0 - value[0]) / value[1]
    if fraction < -1e-12 or fraction > 1.0 + 1e-12:
        raise ArithmeticError("greedy projection fraction escaped [0,1]")
    fraction = float(np.clip(fraction, 0.0, 1.0))
    root = float(np.sqrt(2.0 * fraction))
    unit = 1.0 + root
    minimum = 0.5 * (5.0 - 2.0 / unit - 1.0 / (unit * unit))
    bias = 1.0 / unit
    point = np.asarray([1.0, 2.0 * fraction - 1.0, -1.0, -1.0])
    return {
        "weights": value.tolist(),
        "minimum_total_prior_reserve": minimum,
        "audit_upper": min(1.0, 1.0 / minimum),
        "dual_spectral_bias": bias,
        "relaxed_bloch_projections": point.tolist(),
        "projection_residual": float(np.dot(value, point)),
        "proof": "concave fractional-knapsack vertex plus scalar closed form",
    }


def closed_four_active_complement_reserve(
    weights: np.ndarray, excluded: int = 0
) -> dict[str, Any]:
    """Return a closed-form reserve lower bound outside one sorted effect.

    Concavity fills the zero-cost excluded coordinate first and then the
    largest remaining coordinate.  Exact bias minimisation is a radical for
    fractional demand below one half.  Above that threshold, monotonicity
    gives the simpler (possibly non-tight) lower bound one.
    """

    value = np.sort(np.asarray(weights, dtype=float))[::-1]
    if value.shape != (4,):
        raise ValueError("expected four active effect weights")
    if np.any(value <= 0.0) or np.any(value > 1.0 + 1e-12):
        raise ValueError("active rank-one effect weights must lie in (0,1]")
    if abs(float(value.sum()) - 2.0) > 1e-9:
        raise ValueError("qubit rank-one effect weights must sum to two")
    if excluded not in range(4):
        raise ValueError("excluded effect index must be in range(4)")
    largest_other = 1 if excluded == 0 else 0
    fraction = float((1.0 - value[excluded]) / value[largest_other])
    if fraction < -1e-12:
        raise ArithmeticError("fractional demand became negative")
    fraction = max(0.0, fraction)
    if fraction >= 0.5:
        minimum = 1.0
        point = None
        residual = None
    else:
        radical = fraction * (
            4.0 * fraction * fraction - 7.0 * fraction + 3.0
        )
        minimum = (
            2.0
            * (
                fraction * (4.0 * fraction - 3.0)
                + float(np.sqrt(radical))
            )
            / ((1.0 - 2.0 * fraction) ** 2)
        )
        fill = np.zeros(4)
        fill[excluded] = 1.0
        fill[largest_other] = fraction
        point = 2.0 * fill - 1.0
        residual = float(np.dot(value, point))
    return {
        "weights": value.tolist(),
        "fractional_fill": fraction,
        "minimum_complement_prior_reserve": minimum,
        "excluded_effect": excluded,
        "relaxed_bloch_projections": None if point is None else point.tolist(),
        "projection_residual": residual,
        "proof": "concave fractional-knapsack vertex plus scalar radical",
    }


def closed_four_active_pair_reserve(
    weights: np.ndarray, pair: tuple[int, int]
) -> dict[str, Any]:
    """Return the exact reserve forced on a pair of sorted effects."""

    value = np.sort(np.asarray(weights, dtype=float))[::-1]
    if value.shape != (4,):
        raise ValueError("expected four active effect weights")
    if np.any(value <= 0.0) or np.any(value > 1.0 + 1e-12):
        raise ValueError("active rank-one effect weights must lie in (0,1]")
    if abs(float(value.sum()) - 2.0) > 1e-9:
        raise ValueError("qubit rank-one effect weights must sum to two")
    first, second = pair
    if not 0 <= first < second < 4:
        raise ValueError("pair must satisfy 0 <= first < second < 4")
    residual = max(0.0, float(value[first] + value[second] - 1.0))
    fraction = residual / float(value[first])
    if fraction < 0.0 or fraction > 1.0 + 1e-12:
        raise ArithmeticError("pair fractional fill escaped [0,1]")
    fraction = float(np.clip(fraction, 0.0, 1.0))
    if abs(fraction - 0.5) <= 1e-12:
        minimum = 7.0 / 8.0
    else:
        minimum = (
            2.0 * np.sqrt(2.0 * fraction) * (1.0 - fraction)
            + fraction * (6.0 * fraction - 5.0)
        ) / ((1.0 - 2.0 * fraction) ** 2)
    fill = np.zeros(4)
    remaining = 1.0
    outside = [index for index in range(4) if index not in pair]
    for index in sorted(outside, key=lambda item: value[item], reverse=True):
        amount = min(1.0, remaining / value[index])
        fill[index] = amount
        remaining -= amount * value[index]
        if remaining <= 1e-14:
            break
    if remaining > 1e-14:
        fill[first] = remaining / value[first]
        remaining = 0.0
    point = 2.0 * fill - 1.0
    return {
        "weights": value.tolist(),
        "pair": list(pair),
        "fractional_fill": fraction,
        "minimum_pair_prior_reserve": float(max(0.0, minimum)),
        "relaxed_bloch_projections": point.tolist(),
        "projection_residual": float(np.dot(value, point)),
        "proof": "concave subset knapsack plus scalar radical",
    }


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
        def objective(t: float) -> float:
            return float(np.sum(reserve(t, point)))

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

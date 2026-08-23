"""Diagnostic outer search for a fully active four-outcome qubit readout.

Helstrom complementarity makes the aligned-projective comparisons exact in
the fully active sector.  Write the dual as

    Y = A (I + t z.sigma) / 2

and the active POVM supports as ``Pi_i=(I+n_i.sigma)/2``.  From
``(Y-rho_i)Pi_i=0`` one obtains

    rho_i = Y - (A-p_i) Pi_i^perp,

so positivity is precisely ``p_i >= A f_t(z.n_i)`` and the projective
replacement retaining ``Pi_i`` and complementing it to label ``j`` has

    A2(i,j) = ((1-n_i.n_j) A + (1+n_i.n_j) p_j) / 2.

For fixed quadrilateral geometry the remaining maximisation over ``A``, the
four priors, and the Hellinger RETURN is an SOCP.  Differential evolution
over the five exact geometry parameters is a falsification diagnostic only;
it is not a global upper certificate until the compact domain is covered by
outward interval bounds.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import cvxpy as cp
import numpy as np
from scipy.optimize import differential_evolution

from active_readout_audit_cap import reserve
from active_readout_geometry_probe import (
    bloch_vectors,
    diagonal_interval,
    projections,
)


class ExactProjectiveOuterOracle:
    """DPP SOCP for the prior and audit optimisation at fixed geometry."""

    def __init__(
        self,
        support_weight: float,
        projective_lines: tuple[tuple[float, float], ...],
    ) -> None:
        self.support_weight = float(support_weight)
        self.projective_lines = projective_lines
        self.reserve = cp.Parameter(4, nonneg=True)
        self.one_minus_overlap = cp.Parameter((4, 4), nonneg=True)
        self.one_plus_overlap = cp.Parameter((4, 4), nonneg=True)
        self.prior = cp.Variable(4, nonneg=True)
        self.audit = cp.Variable(nonneg=True)
        self.returned = cp.Variable(nonneg=True)
        constraints: list[cp.Constraint] = [
            cp.sum(self.prior) == 1.0,
            self.prior >= self.audit * self.reserve,
            # Y-rho_i is positive for every Helstrom constraint, hence its
            # trace A-p_i is nonnegative.
            self.prior <= self.audit,
        ]
        cross = []
        for first in range(4):
            for second in range(first + 1, 4):
                geometric = cp.Variable(nonneg=True)
                constraints.append(
                    cp.SOC(
                        self.prior[first] + self.prior[second],
                        cp.hstack(
                            [
                                2.0 * geometric,
                                self.prior[first] - self.prior[second],
                            ]
                        ),
                    )
                )
                cross.append(geometric)
        constraints.append(
            self.returned
            <= (cp.sum(self.prior) + 2.0 * cp.sum(cp.hstack(cross))) / 4.0
        )
        self.projective_constraints: list[
            tuple[float, float, int, int, cp.Constraint]
        ] = []
        for retained in range(4):
            for complement in range(4):
                if retained == complement:
                    continue
                projective_audit = 0.5 * (
                    self.one_minus_overlap[retained, complement] * self.audit
                    + self.one_plus_overlap[retained, complement]
                    * self.prior[complement]
                )
                for line_weight, line_upper in projective_lines:
                    constraint = (
                        line_weight * projective_audit
                        + (1.0 - line_weight) * self.returned
                        <= line_upper
                    )
                    constraints.append(constraint)
                    self.projective_constraints.append(
                        (
                            line_weight,
                            line_upper,
                            retained,
                            complement,
                            constraint,
                        )
                    )
        self.score = (
            self.support_weight * self.audit
            + (1.0 - self.support_weight) * self.returned
        )
        self.problem = cp.Problem(cp.Maximize(self.score), constraints)
        if not self.problem.is_dpp():
            raise RuntimeError("the fixed-geometry SOCP must be DPP")

    def solve(
        self,
        reserve_vector: np.ndarray,
        overlaps: np.ndarray,
        *,
        safety: float = 0.0,
        capture: bool = False,
    ) -> dict[str, Any]:
        clipped = np.clip(np.asarray(overlaps, dtype=float), -1.0, 1.0)
        self.reserve.value = np.asarray(reserve_vector, dtype=float)
        self.one_minus_overlap.value = 1.0 - clipped
        self.one_plus_overlap.value = 1.0 + clipped
        self.problem.solve(
            solver="CLARABEL",
            tol_gap_abs=2e-9,
            tol_gap_rel=2e-9,
            tol_feas=2e-9,
            max_iter=500,
            warm_start=True,
            ignore_dpp=False,
        )
        if self.problem.status not in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}:
            return {"status": self.problem.status, "bound": float("inf")}
        result: dict[str, Any] = {
            "status": self.problem.status,
            "raw_value": float(self.problem.value),
            "bound": float(self.problem.value) + safety,
        }
        if capture:
            result.update(
                {
                    "audit": float(self.audit.value),
                    "return": float(self.returned.value),
                    "syndrome_priors": np.asarray(self.prior.value).tolist(),
                    "projective_line_duals": [
                        {
                            "weight": line_weight,
                            "upper": line_upper,
                            "retained": retained,
                            "complement": complement,
                            "dual": (
                                None
                                if constraint.dual_value is None
                                else float(constraint.dual_value)
                            ),
                        }
                        for (
                            line_weight,
                            line_upper,
                            retained,
                            complement,
                            constraint,
                        ) in self.projective_constraints
                    ],
                }
            )
        return result


def fixed_weight_probe(
    weights: np.ndarray,
    support_weight: float,
    projective_lines: tuple[tuple[float, float], ...],
    seed: int,
    maxiter: int,
    popsize: int,
) -> dict[str, Any]:
    value = np.sort(np.asarray(weights, dtype=float))[::-1]
    if value.shape != (4,) or np.any(value <= 0.0) or np.any(value > 1.0):
        raise ValueError("expected four active weights in (0,1]")
    if abs(float(value.sum()) - 2.0) > 1e-9:
        raise ValueError("qubit effect weights must sum to two")
    lower, upper = diagonal_interval(value)
    epsilon = 1e-10 if lower == 0.0 else 0.0
    oracle = ExactProjectiveOuterOracle(support_weight, projective_lines)

    def evaluate(variables: np.ndarray, capture: bool = False) -> dict[str, Any]:
        vectors = bloch_vectors(value, variables)
        point = projections(value, variables)
        reserve_vector = np.asarray(
            reserve(float(variables[4]), point), dtype=float
        )
        result = oracle.solve(
            reserve_vector, vectors @ vectors.T, capture=capture
        )
        if capture:
            result.update(
                {
                    "prior_reserves": reserve_vector.tolist(),
                    "bloch_projections": point.tolist(),
                    "bloch_overlaps": (vectors @ vectors.T).tolist(),
                    "closure_residual": float(
                        np.linalg.norm((value[:, None] * vectors).sum(axis=0))
                    ),
                }
            )
        return result

    result = differential_evolution(
        lambda variables: -float(evaluate(variables)["bound"]),
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
        tol=1e-9,
        polish=True,
        updating="immediate",
        workers=1,
    )
    variables = np.asarray(result.x, dtype=float)
    best = evaluate(variables, capture=True)
    return {
        "weights": value.tolist(),
        "support_weight": float(support_weight),
        "projective_lines": [list(line) for line in projective_lines],
        "outer_support_found": best["raw_value"],
        "geometry": {
            "diagonal": float(variables[0]),
            "dual_axis_cosine": float(variables[1]),
            "pair_azimuths": variables[2:4].tolist(),
            "dual_spectral_bias": float(variables[4]),
        },
        **{key: value for key, value in best.items() if key != "raw_value"},
        "function_evaluations": int(result.nfev),
        "optimizer_success": bool(result.success),
        "optimizer_message": str(result.message),
        "status_note": "global-search diagnostic; not an upper certificate",
    }


def weight_chart(
    maximum_weight: float,
    smallest_fraction: float,
    middle_fraction: float,
    minimum_active_weight: float,
) -> np.ndarray:
    """Map a unit square to sorted qubit weights with prescribed maximum."""

    residual = 2.0 - maximum_weight
    smallest_upper = residual / 3.0
    smallest = minimum_active_weight + smallest_fraction * (
        smallest_upper - minimum_active_weight
    )
    middle_lower = max(smallest, residual - smallest - maximum_weight)
    middle_upper = (residual - smallest) / 2.0
    middle = middle_lower + middle_fraction * (middle_upper - middle_lower)
    second = residual - middle - smallest
    return np.asarray([maximum_weight, second, middle, smallest])


def global_weight_probe(
    maximum_weight_floor: float,
    support_weight: float,
    projective_lines: tuple[tuple[float, float], ...],
    seed: int,
    maxiter: int,
    popsize: int,
    minimum_active_weight: float = 1e-6,
) -> dict[str, Any]:
    """Search all sorted active weights and exact quadrilateral geometries."""

    oracle = ExactProjectiveOuterOracle(support_weight, projective_lines)

    def decode(variables: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        weights = weight_chart(
            float(variables[0]),
            float(variables[1]),
            float(variables[2]),
            minimum_active_weight,
        )
        lower, upper = diagonal_interval(weights)
        diagonal = lower + float(variables[3]) * (upper - lower)
        geometry = np.asarray(
            [diagonal, *variables[4:]], dtype=float
        )
        return weights, geometry

    def evaluate(variables: np.ndarray, capture: bool = False) -> dict[str, Any]:
        weights, geometry = decode(variables)
        vectors = bloch_vectors(weights, geometry)
        point = projections(weights, geometry)
        reserve_vector = np.asarray(
            reserve(float(geometry[4]), point), dtype=float
        )
        result = oracle.solve(
            reserve_vector, vectors @ vectors.T, capture=capture
        )
        if capture:
            result.update(
                {
                    "weights": weights.tolist(),
                    "geometry_variables": geometry.tolist(),
                    "prior_reserves": reserve_vector.tolist(),
                    "bloch_projections": point.tolist(),
                    "bloch_overlaps": (vectors @ vectors.T).tolist(),
                    "closure_residual": float(
                        np.linalg.norm((weights[:, None] * vectors).sum(axis=0))
                    ),
                }
            )
        return result

    result = differential_evolution(
        lambda variables: -float(evaluate(variables)["bound"]),
        bounds=(
            (maximum_weight_floor, 1.0),
            (0.0, 1.0),
            (0.0, 1.0),
            (0.0, 1.0),
            (-1.0, 1.0),
            (0.0, 2.0 * np.pi),
            (0.0, 2.0 * np.pi),
            (0.0, 1.0 - 1e-10),
        ),
        seed=seed,
        maxiter=maxiter,
        popsize=popsize,
        tol=2e-8,
        polish=True,
        updating="immediate",
        workers=1,
    )
    variables = np.asarray(result.x, dtype=float)
    best = evaluate(variables, capture=True)
    geometry = np.asarray(best.pop("geometry_variables"), dtype=float)
    return {
        "maximum_weight_floor": float(maximum_weight_floor),
        "minimum_active_weight": float(minimum_active_weight),
        "support_weight": float(support_weight),
        "projective_lines": [list(line) for line in projective_lines],
        "outer_support_found": best["raw_value"],
        "weights": best.pop("weights"),
        "geometry": {
            "diagonal": float(geometry[0]),
            "dual_axis_cosine": float(geometry[1]),
            "pair_azimuths": geometry[2:4].tolist(),
            "dual_spectral_bias": float(geometry[4]),
        },
        **{key: value for key, value in best.items() if key != "raw_value"},
        "function_evaluations": int(result.nfev),
        "optimizer_success": bool(result.success),
        "optimizer_message": str(result.message),
        "status_note": "global-search diagnostic; not an upper certificate",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--weights", type=float, nargs=4)
    parser.add_argument(
        "--global-weights",
        action="store_true",
        help="also optimise the sorted four-effect weights",
    )
    parser.add_argument("--maximum-weight-floor", type=float, default=0.88325)
    parser.add_argument("--lambda", dest="support_weight", type=float, default=0.6)
    parser.add_argument(
        "--projective-line",
        type=float,
        nargs=2,
        action="append",
        default=[(0.55, 0.7573), (0.6, 0.76591)],
        metavar=("WEIGHT", "UPPER"),
    )
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--maxiter", type=int, default=300)
    parser.add_argument("--popsize", type=int, default=16)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    lines = tuple(tuple(map(float, line)) for line in args.projective_line)
    if args.global_weights:
        payload = global_weight_probe(
            args.maximum_weight_floor,
            args.support_weight,
            lines,
            args.seed,
            args.maxiter,
            args.popsize,
        )
    else:
        if args.weights is None:
            parser.error("--weights is required unless --global-weights is used")
        payload = fixed_weight_probe(
            np.asarray(args.weights, dtype=float),
            args.support_weight,
            lines,
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

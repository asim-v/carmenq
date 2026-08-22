"""Two-parameter active branch of the interleaved qubit frontier.

After relaxing the first two slots to an arbitrary four-state qubit ensemble,
the numerically optimal suffix is represented by a three-effect extremal
qubit POVM.  One effect in the nominal four-outcome POVM vanishes.  POVM
completeness fixes all effect directions and weights from a single trace
parameter ``t``.  A second Bloch coordinate ``r`` fixes the two symmetric
signal states.  Optimisation over the four prior weights is the largest
eigenvalue of a 3-by-3 arrow-plus-rank-one matrix.

This script evaluates and continues that reduced active branch.  Numerical
agreement with the unrestricted Choi-MPS and cq-instrument searches is not,
by itself, a proof that the reduction is globally exhaustive.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import differential_evolution, minimize, root


def reduced_matrix(t: float, r: float, weight: float) -> tuple[np.ndarray, dict[str, float]]:
    if not 0.0 <= t <= 1.0:
        raise ValueError("t must lie in [0,1]")
    if not -1.0 <= r <= 1.0:
        raise ValueError("r must lie in [-1,1]")

    other_trace = 1.0 - t / 2.0
    effect_z = -t / (2.0 - t)
    effect_x = np.sqrt(max(0.0, 1.0 - effect_z**2))
    state_x = np.sqrt(max(0.0, 1.0 - r**2))
    correct_other = other_trace * (
        1.0 + effect_z * r + effect_x * state_x
    ) / 2.0

    coarse_reference = t
    coarse_other = t * (1.0 + r) / 2.0
    hellinger_reference = np.sqrt(coarse_reference) + np.sqrt(1.0 - coarse_reference)
    hellinger_other = np.sqrt(coarse_other) + np.sqrt(1.0 - coarse_other)
    # The null AUDIT label carries no correctness reward.  Its state therefore
    # maximises only the coarse-outcome Hellinger factor.  When t >= 1/2 it can
    # choose Tr(E_0 rho)=1/2 exactly; below that threshold it aligns with E_0
    # and reaches only probability t.  This is the small active-phase kink
    # visible just above coexistence with the no-record strategy.
    hellinger_null = np.sqrt(2.0) if t >= 0.5 else hellinger_reference

    # Symmetry combines the last two labels into one normalised coordinate.
    audit_diagonal = np.asarray((t, 0.0, correct_other), dtype=float)
    hellinger_vector = np.asarray(
        (hellinger_reference, hellinger_null, np.sqrt(2.0) * hellinger_other),
        dtype=float,
    )
    matrix = weight * np.diag(audit_diagonal)
    matrix += (1.0 - weight) * np.outer(hellinger_vector, hellinger_vector) / 8.0
    return matrix, {
        "other_effect_trace": other_trace,
        "other_effect_z": effect_z,
        "other_correct_probability": correct_other,
        "coarse_reference_probability": coarse_reference,
        "coarse_other_probability": coarse_other,
        "hellinger_reference": hellinger_reference,
        "hellinger_null": hellinger_null,
        "hellinger_other": hellinger_other,
    }


def evaluate(t: float, r: float, weight: float) -> dict[str, float | list[float]]:
    matrix, invariants = reduced_matrix(t, r, weight)
    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    vector = np.abs(eigenvectors[:, -1])
    vector /= np.linalg.norm(vector)
    priors = np.asarray((vector[0] ** 2, vector[1] ** 2, vector[2] ** 2 / 2, vector[2] ** 2 / 2))

    t_audit = t
    other_audit = invariants["other_correct_probability"]
    audit = priors[0] * t_audit + (priors[2] + priors[3]) * other_audit
    c0 = invariants["hellinger_reference"]
    c_null = invariants["hellinger_null"]
    c1 = invariants["hellinger_other"]
    returned = (
        c0 * np.sqrt(priors[0])
        + c_null * np.sqrt(priors[1])
        + c1 * (np.sqrt(priors[2]) + np.sqrt(priors[3]))
    ) ** 2 / 8.0
    score = weight * audit + (1.0 - weight) * returned
    return {
        "weight": weight,
        "t": t,
        "r": r,
        "score": float(score),
        "audit": float(audit),
        "return": float(returned),
        "priors": [float(value) for value in priors],
        "largest_eigenvalue": float(eigenvalues[-1]),
        **invariants,
    }


def optimise(weight: float, initial: tuple[float, float] | None = None) -> dict[str, float | list[float]]:
    # Optimise in u=sqrt(t).  The high-AUDIT branch approaches t=0 through a
    # very narrow basin, which a linear t coordinate can miss numerically.
    def objective(point: np.ndarray) -> float:
        return -float(evaluate(float(point[0]) ** 2, float(point[1]), weight)["score"])

    starts: list[tuple[float, float]] = []
    if initial is None:
        global_result = differential_evolution(
            objective,
            bounds=((0.0, 1.0), (-1.0 + 1e-10, 1.0 - 1e-10)),
            seed=20260822,
            tol=1e-12,
            polish=False,
        )
        starts.append(tuple(float(value) for value in global_result.x))
    else:
        starts.append((np.sqrt(max(0.0, float(initial[0]))), float(initial[1])))
    starts.extend(((0.7, 0.0), (0.2, 0.0), (0.05, 0.0), (0.005, 0.0), (0.999, -0.6)))

    candidates = []
    for start in starts:
        local = minimize(
            objective,
            np.asarray(start),
            method="Nelder-Mead",
            bounds=((0.0, 1.0), (-1.0 + 1e-12, 1.0 - 1e-12)),
            options={"xatol": 1e-13, "fatol": 1e-15, "maxiter": 20000},
        )
        u, r = map(float, local.x)
        candidates.append(evaluate(u**2, r, weight))
    return max(candidates, key=lambda row: float(row["score"]))


def transition() -> dict[str, float | list[float]]:
    """Solve stationarity and coexistence with the no-record line."""

    def equations(point: np.ndarray) -> np.ndarray:
        t, r, weight = map(float, point)
        step = 2e-6
        base = float(evaluate(t, r, weight)["score"])
        dt = (
            float(evaluate(t + step, r, weight)["score"])
            - float(evaluate(t - step, r, weight)["score"])
        ) / (2 * step)
        dr = (
            float(evaluate(t, r + step, weight)["score"])
            - float(evaluate(t, r - step, weight)["score"])
        ) / (2 * step)
        return np.asarray((dt, dr, base - (1.0 - weight / 2.0)))

    solution = root(equations, np.asarray((0.52, -0.10, 0.44144)), tol=1e-11)
    if not solution.success:
        raise RuntimeError(solution.message)
    t, r, weight = map(float, solution.x)
    return evaluate(t, r, weight)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lambda", dest="weight", type=float)
    parser.add_argument("--sweep", action="store_true")
    parser.add_argument("--transition", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    payload: object
    if args.transition:
        payload = transition()
    elif args.sweep:
        rows = []
        initial = None
        for weight in (0.42, 0.44, 0.45, 0.48, 0.5, 0.52, 0.55, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99):
            row = optimise(weight, initial)
            initial = (float(row["t"]), float(row["r"]))
            row["no_record"] = 1.0 - weight / 2.0
            row["active_advantage"] = float(row["score"]) - float(row["no_record"])
            rows.append(row)
        payload = rows
    else:
        payload = optimise(0.5 if args.weight is None else args.weight)

    rendered = json.dumps(payload, indent=2) + "\n"
    print(rendered)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

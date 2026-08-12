"""Fast independent checks for the streaming parity audit--return theorem.

This validator is deliberately separate from ``streaming_parity.py``.  It
checks the closed forms, the phase-transition equations, stored adversarial
searches, the collective classical-record bound, and the coherent accumulator
without regenerating the expensive optimization artifact.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, minimize_scalar


HERE = Path(__file__).resolve().parent
RNG = np.random.default_rng(20260812)


def f_return(t: np.ndarray | float) -> np.ndarray | float:
    return (1.0 + np.sqrt(np.maximum(0.0, 1.0 - np.square(t)))) / 2.0


def online_score(n: int, weight: float, t: float) -> float:
    return float(
        weight * (1.0 + t**n) / 2.0
        + (1.0 - weight) * f_return(t) ** n
    )


def online_support(n: int, weight: float) -> float:
    grid = np.unique(
        np.concatenate(
            (
                np.linspace(0.0, 0.995, 20_001),
                1.0 - np.geomspace(1e-14, 5e-3, 10_001),
                np.array([1.0]),
            )
        )
    )
    values = np.asarray([online_score(n, weight, t) for t in grid])
    candidates = [float(values[0]), float(values[-1])]
    for index in np.argpartition(values, -8)[-8:]:
        lo = float(grid[max(index - 2, 0)])
        hi = float(grid[min(index + 2, grid.size - 1)])
        result = minimize_scalar(
            lambda t: -online_score(n, weight, t),
            bounds=(lo, hi),
            method="bounded",
            options={"xatol": 1e-14},
        )
        candidates.append(float(-result.fun))
    return max(candidates)


def collective_support(weight: float) -> float:
    return float(
        (1.0 + np.sqrt(weight**2 + (1.0 - weight) ** 2)) / 2.0
    )


def check_closed_forms() -> float:
    worst = 0.0
    for weight in np.linspace(0.0, 1.0, 41):
        one_slot = collective_support(float(weight))
        worst = max(worst, abs(online_support(1, float(weight)) - one_slot))
        two_slot = (
            1.0 - weight / 2.0
            if weight <= 0.5
            else weight / 2.0 + weight**2 / (3.0 * weight - 1.0)
        )
        worst = max(worst, abs(online_support(2, float(weight)) - two_slot))
    assert worst < 2e-11, worst
    return worst


def check_transition_equations() -> tuple[float, dict[str, dict[str, float]]]:
    worst = 0.0
    rows: dict[str, dict[str, float]] = {}
    for n in range(3, 13):
        lower = (n - 2.0) / n
        root = brentq(
            lambda z: (1.0 - z) * (1.0 + z) ** (n - 1) - 1.0,
            lower + 1e-13,
            1.0 - 1e-13,
        )
        weight = 1.0 / (
            1.0
            + 2.0 ** (n - 2)
            * root ** ((n - 2.0) / 2.0)
            * (1.0 - root)
        )
        t = 2.0 * np.sqrt(root) / (1.0 + root)
        baseline = 1.0 - weight / 2.0
        finite = online_score(n, weight, float(t))
        stationary = (
            weight
            * 2.0 ** (n - 2)
            * root ** ((n - 2.0) / 2.0)
            * (1.0 - root)
            - (1.0 - weight)
        )
        error = max(abs(finite - baseline), abs(stationary))
        worst = max(worst, error)
        rows[str(n)] = {
            "lambda_c": float(weight),
            "jump_strength": float(t),
        }
    assert worst < 1e-9, worst
    return worst, rows


def check_product_recovery_inequality() -> float:
    max_excess = -np.inf
    for n in range(2, 11):
        samples = RNG.random((20_000, n))
        lhs = np.prod(f_return(samples), axis=1)
        rhs = f_return(np.prod(samples, axis=1))
        max_excess = max(max_excess, float(np.max(lhs - rhs)))
    assert max_excess <= 2e-15, max_excess
    return max_excess


def check_collective_classical_record_bound() -> float:
    max_excess = -np.inf
    for n in (2, 3, 4):
        dimension = 2**n
        parity = np.arange(dimension).astype(np.uint64)
        parity = np.asarray([int(value).bit_count() & 1 for value in parity])
        for _ in range(2_000):
            outcomes = int(RNG.integers(2, 9))
            likelihood = RNG.dirichlet(np.ones(outcomes), size=dimension).T
            even = likelihood[:, parity == 0].sum(axis=1)
            odd = likelihood[:, parity == 1].sum(axis=1)
            audit = 0.5 + np.abs(even - odd).sum() / (2.0 * dimension)
            recovery_upper = np.square(np.sqrt(likelihood).sum(axis=1)).sum()
            recovery_upper /= dimension**2
            total_variation = np.abs(even - odd).sum() / dimension
            max_excess = max(
                max_excess,
                float(recovery_upper - f_return(total_variation)),
            )
            for weight in (0.2, 0.5, 0.8):
                score = weight * audit + (1.0 - weight) * recovery_upper
                max_excess = max(
                    max_excess,
                    float(score - collective_support(weight)),
                )
    assert max_excess <= 2e-14, max_excess
    return max_excess


def check_coherent_accumulator() -> float:
    worst = 0.0
    for n in range(1, 8):
        dimension = 2**n
        state = np.zeros((dimension, dimension, 2), dtype=np.complex128)
        for x in range(dimension):
            state[x, x, x.bit_count() & 1] = 1.0 / np.sqrt(dimension)
        audit = sum(
            abs(state[x, x, x.bit_count() & 1]) ** 2
            for x in range(dimension)
        )
        returned = np.zeros_like(state)
        for reference in range(dimension):
            for carrier in range(dimension):
                parity = carrier.bit_count() & 1
                returned[reference, carrier, 0] = state[
                    reference, carrier, parity
                ]
                returned[reference, carrier, 1] = state[
                    reference, carrier, parity ^ 1
                ]
        target = np.zeros_like(state)
        for x in range(dimension):
            target[x, x, 0] = 1.0 / np.sqrt(dimension)
        fidelity = abs(np.vdot(target.ravel(), returned.ravel())) ** 2
        worst = max(worst, abs(audit - 1.0), abs(fidelity - 1.0))
    assert worst < 2e-15, worst
    return worst


def check_stored_adversarial_searches() -> dict[str, float]:
    payload = json.loads((HERE / "streaming_parity.json").read_text("utf-8"))
    adaptive_excess = max(
        row["adaptive_minus_symmetric"] for row in payload["rows"]
    )
    asymmetric_excess = max(
        row["unrestricted_product"]["score"]
        - row["symmetric_product"]["score"]
        for row in payload["rows"]
    )
    minimum_collective_gap = min(
        row["global_minus_online"] for row in payload["rows"]
    )
    assert adaptive_excess < 1e-10, adaptive_excess
    assert asymmetric_excess < 1e-10, asymmetric_excess
    assert minimum_collective_gap > 0.0, minimum_collective_gap
    return {
        "maximum_adaptive_excess": float(adaptive_excess),
        "maximum_asymmetric_excess": float(asymmetric_excess),
        "minimum_tested_collective_gap": float(minimum_collective_gap),
    }


def main() -> None:
    transition_error, transitions = check_transition_equations()
    report = {
        "closed_form_max_abs_error": check_closed_forms(),
        "transition_equation_max_abs_error": transition_error,
        "transitions_n3_to_n12": transitions,
        "product_recovery_max_excess": check_product_recovery_inequality(),
        "collective_bound_max_excess": check_collective_classical_record_bound(),
        "coherent_accumulator_max_abs_error": check_coherent_accumulator(),
        "stored_adversarial_searches": check_stored_adversarial_searches(),
        "status": "PASS",
    }
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

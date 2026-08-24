"""Separate a reported qubit family by flagged trace-norm contraction.

For a common quantum instrument and every real coefficient vector ``c``,

    sum_y || sum_z c[z] sigma[z,y] ||_1
        <= || sum_z c[z] rho[z] ||_1.

The inequality is homogeneous in ``c``.  This module searches the unit
three-sphere for the largest violation of a reported outer-relaxation family.
It is a diagnostic separation oracle, not by itself a global certificate.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize


def qubit_trace_norm(coefficients: np.ndarray) -> float:
    """Trace norm of a Hermitian qubit operator in Pauli coordinates."""

    value = np.asarray(coefficients, dtype=float)
    if value.shape != (4,):
        raise ValueError("qubit Pauli coordinates must have length four")
    return max(abs(float(value[0])), float(np.linalg.norm(value[1:])))


def contraction_report(
    prefix: np.ndarray,
    conditioned: np.ndarray,
    coefficients: np.ndarray,
) -> dict[str, object]:
    """Evaluate one homogeneous flagged contraction."""

    prefix = np.asarray(prefix, dtype=float)
    conditioned = np.asarray(conditioned, dtype=float)
    coefficients = np.asarray(coefficients, dtype=float)
    if prefix.shape != (4, 4) or conditioned.shape != (4, 4, 4):
        raise ValueError("expected four inputs and sixteen qubit outputs")
    norm = float(np.linalg.norm(coefficients))
    if norm <= 1e-14:
        raise ValueError("contraction coefficients must be nonzero")
    coefficients = coefficients / norm
    input_combination = coefficients @ prefix
    output_combinations = np.einsum("z,zyd->yd", coefficients, conditioned)
    input_norm = qubit_trace_norm(input_combination)
    output_norms = np.asarray(
        [qubit_trace_norm(item) for item in output_combinations]
    )
    return {
        "coefficients": coefficients.tolist(),
        "input_combination": input_combination.tolist(),
        "output_combinations": output_combinations.tolist(),
        "input_trace_norm": input_norm,
        "output_trace_norms": output_norms.tolist(),
        "flagged_output_trace_norm": float(output_norms.sum()),
        "violation": float(output_norms.sum() - input_norm),
        "input_branch": (
            "scalar-positive"
            if input_combination[0] >= np.linalg.norm(input_combination[1:])
            else "scalar-negative"
            if -input_combination[0] >= np.linalg.norm(input_combination[1:])
            else "bloch"
        ),
        "input_bloch_direction": (
            None
            if np.linalg.norm(input_combination[1:]) <= 1e-14
            else (
                input_combination[1:] / np.linalg.norm(input_combination[1:])
            ).tolist()
        ),
    }


def find_worst_contraction(
    prefix: np.ndarray,
    conditioned: np.ndarray,
    samples: int = 50_000,
    starts: int = 32,
    seed: int = 260824,
) -> dict[str, object]:
    """Use deterministic sphere sampling followed by local refinement."""

    if samples < 1 or starts < 1:
        raise ValueError("samples and starts must be positive")
    rng = np.random.default_rng(seed)
    candidates = rng.normal(size=(samples, 4))
    candidates /= np.linalg.norm(candidates, axis=1, keepdims=True)
    # Always include simple pair differences and the Z2^2 characters.
    structured = []
    for first in range(4):
        for second in range(first + 1, 4):
            item = np.zeros(4)
            item[first], item[second] = 1.0, -1.0
            structured.append(item / np.linalg.norm(item))
    structured.extend(
        np.asarray(
            [
                [1.0, 1.0, -1.0, -1.0],
                [1.0, -1.0, 1.0, -1.0],
                [1.0, -1.0, -1.0, 1.0],
            ]
        )
        / 2.0
    )
    candidates = np.vstack([np.asarray(structured), candidates])

    def objective(raw: np.ndarray) -> float:
        norm = float(np.linalg.norm(raw))
        if norm <= 1e-14:
            return 0.0
        return -float(
            contraction_report(prefix, conditioned, raw / norm)["violation"]
        )

    sampled = np.asarray([-objective(item) for item in candidates])
    best_indices = np.argsort(sampled)[-min(starts, len(sampled)) :]
    best_vector = candidates[int(best_indices[-1])]
    best_violation = float(sampled[int(best_indices[-1])])
    refinements = []
    for index in best_indices:
        result = minimize(
            objective,
            candidates[int(index)],
            method="BFGS",
            options={"gtol": 1e-10, "maxiter": 1000},
        )
        report = contraction_report(prefix, conditioned, result.x)
        refinements.append(
            {
                "success": bool(result.success),
                "iterations": int(result.nit),
                "violation": report["violation"],
            }
        )
        if float(report["violation"]) > best_violation:
            best_violation = float(report["violation"])
            best_vector = np.asarray(report["coefficients"], dtype=float)
    report = contraction_report(prefix, conditioned, best_vector)
    report.update(
        {
            "samples": samples,
            "local_starts": starts,
            "seed": seed,
            "largest_sampled_violation": float(sampled.max()),
            "refinements": refinements,
            "epistemic_status": (
                "nonconvex numerical separation over coefficient space; "
                "a positive value is a directly checkable incompatibility witness"
            ),
        }
    )
    return report


def family_from_payload(payload: dict[str, object]) -> tuple[np.ndarray, np.ndarray]:
    """Read Pauli-coordinate inputs and outputs from a behavior payload."""

    prefix = np.asarray(payload["prefix_bloch_coefficients"], dtype=float)
    probabilities = np.asarray(payload["path_probabilities"], dtype=float)
    vectors = np.asarray(payload["conditioned_output_bloch_vectors"], dtype=float)
    conditioned = np.concatenate((probabilities[:, :, None], vectors), axis=2)
    return prefix, conditioned


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("--samples", type=int, default=50_000)
    parser.add_argument("--starts", type=int, default=32)
    parser.add_argument("--seed", type=int, default=260824)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = json.loads(args.source.read_text(encoding="utf-8"))
    prefix, conditioned = family_from_payload(payload)
    report = find_worst_contraction(
        prefix, conditioned, args.samples, args.starts, args.seed
    )
    rendered = json.dumps(report, indent=2) + "\n"
    print(rendered, end="")
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

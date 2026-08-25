"""Exact operator-basis audit for a common effective qubit POVM.

Any common instrument ``{Phi_y}`` followed by one terminal POVM ``{E_t}``
induces the twelve input effects ``F_yt = Phi_y^*(E_t)``.  Hence measured
statistics must factor as

    q[z,y,t] = Tr(rho_z F_yt)

for one POVM on the input qubit.  If four subnormalised input states span the
Hermitian operator space, the effects are unique.  In Pauli coordinates, put
the state rows ``R[z] = (p_z, r_z)`` and flatten ``q`` to a 4-by-12 matrix.
The effect coordinates are simply ``A = R^{-1} q``.  Each column is positive
iff ``A[0,k] >= ||A[1:,k]||_2``; row normalisation of q then makes the effects
sum to the identity automatically.

This is an exact necessary condition for the full common instrument and an
exact necessary-and-sufficient condition for a common effective POVM.  It is
not sufficient for the finer sequential factorisation through the specified
terminal POVM.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np


def audit_common_effective_povm(
    priors: np.ndarray,
    input_bloch_vectors: np.ndarray,
    statistics: np.ndarray,
    singular_tolerance: float = 1e-10,
) -> dict[str, Any]:
    priors_array = np.asarray(priors, dtype=float)
    bloch = np.asarray(input_bloch_vectors, dtype=float)
    measured = np.asarray(statistics, dtype=float)
    if priors_array.shape != (4,) or bloch.shape != (4, 3):
        raise ValueError("expected four qubit input Pauli vectors")
    if measured.shape != (4, 4, 3):
        raise ValueError("expected statistics with shape (4,4,3)")
    state_matrix = np.column_stack([priors_array, bloch])
    determinant = float(np.linalg.det(state_matrix))
    condition_number = float(np.linalg.cond(state_matrix))
    if abs(determinant) <= singular_tolerance:
        return {
            "status": "singular_input_basis",
            "determinant": determinant,
            "condition_number": condition_number,
            "singular_tolerance": singular_tolerance,
        }
    probability_matrix = measured.reshape(4, 12)
    effects = np.linalg.solve(state_matrix, probability_matrix)
    radii = np.linalg.norm(effects[1:, :], axis=0)
    margins = effects[0, :] - radii
    maximum_eigenvalues = effects[0, :] + radii
    interpolation_residual = float(
        np.linalg.norm(state_matrix @ effects - probability_matrix)
    )
    completeness_vector = effects.sum(axis=1)
    completeness_residual = float(
        np.linalg.norm(completeness_vector - np.asarray([1.0, 0.0, 0.0, 0.0]))
    )

    worst = int(np.argmin(margins))
    y, t = divmod(worst, 3)
    effect_vector = effects[1:, worst]
    radius = float(radii[worst])
    witness_bloch = (
        np.zeros(3) if radius <= 1e-15 else -effect_vector / radius
    )
    witness_operator = np.append(1.0, witness_bloch)
    interpolation_coefficients = np.linalg.solve(
        state_matrix.T, witness_operator
    )
    witness_expectation = float(
        interpolation_coefficients @ probability_matrix[:, worst]
    )

    # adj(R) q = det(R) R^{-1} q.  On a fixed determinant-sign branch,
    # positivity is the polynomial SOC sign(det)*n0 >= ||nvec||.
    numerators = determinant * effects
    signed_numerators = np.sign(determinant) * numerators
    numerator_margins = signed_numerators[0, :] - np.linalg.norm(
        signed_numerators[1:, :], axis=0
    )
    return {
        "status": "nonsingular",
        "determinant": determinant,
        "determinant_sign": int(np.sign(determinant)),
        "condition_number": condition_number,
        "singular_tolerance": singular_tolerance,
        "effect_pauli_coordinates": effects.T.reshape(4, 3, 4).tolist(),
        "minimum_eigenvalue_margins": margins.reshape(4, 3).tolist(),
        "maximum_eigenvalues": maximum_eigenvalues.reshape(4, 3).tolist(),
        "minimum_margin": float(margins[worst]),
        "negative_effect_count": int(np.count_nonzero(margins < -1e-9)),
        "interpolation_residual": interpolation_residual,
        "completeness_vector": completeness_vector.tolist(),
        "completeness_residual": completeness_residual,
        "worst_effect": {"y": y, "t": t, "flat_index": worst},
        "worst_effect_witness_bloch": witness_bloch.tolist(),
        "worst_effect_interpolation_coefficients": interpolation_coefficients.tolist(),
        "worst_effect_witness_expectation": witness_expectation,
        "minimum_numerator_soc_margin": float(np.min(numerator_margins)),
        "povm_positive": bool(np.min(margins) >= -1e-9),
        "povm_complete": bool(completeness_residual <= 1e-9),
        "common_effective_povm": bool(
            np.min(margins) >= -1e-9 and completeness_residual <= 1e-9
        ),
        "scope": (
            "exact common effective qubit POVM audit at one nonsingular "
            "numerical candidate; necessary for a common instrument"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--frontier-json", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    frontier = json.loads(args.frontier_json.read_text(encoding="utf-8"))
    solution = frontier.get("top_solution")
    if not solution:
        raise ValueError("frontier artifact does not contain a captured top solution")
    audit = audit_common_effective_povm(
        np.asarray(solution["prefix"], dtype=float),
        np.asarray(solution["input_bloch_vectors"], dtype=float),
        np.asarray(solution["statistics"], dtype=float),
    )
    payload = {
        "source": str(args.frontier_json),
        "source_bound": solution.get("bound"),
        **audit,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                key: payload.get(key)
                for key in (
                    "status",
                    "determinant",
                    "condition_number",
                    "minimum_margin",
                    "negative_effect_count",
                    "completeness_residual",
                    "common_effective_povm",
                )
            }
        )
    )


if __name__ == "__main__":
    main()

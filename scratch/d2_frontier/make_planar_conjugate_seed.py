"""Generate the determinant-reflected partner of a planar-POVM strategy.

For effects invariant under ``E -> X E* X``, the transformation

    rho_z -> rho_z*
    J_y   -> (I tensor X) J_y* (I tensor X)

preserves every path and terminal Born probability while reversing the
orientation of a linearly independent qubit input basis.  It therefore maps
the positive and negative determinant branches of the fixed-planar frontier
onto one another exactly.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


IDENTITY = np.eye(2, dtype=complex)
PAULI_X = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
PAULI_Y = np.array([[0.0, -1j], [1j, 0.0]], dtype=complex)
PAULI_Z = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)
PAULIS = (IDENTITY, PAULI_X, PAULI_Y, PAULI_Z)


def input_pauli_matrix(states: np.ndarray) -> np.ndarray:
    return np.asarray(
        [
            [float(np.trace(state @ pauli).real) for pauli in PAULIS]
            for state in states
        ]
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()

    arrays = np.load(args.input)
    states = np.asarray(arrays["states"], dtype=complex)
    choi = np.asarray(arrays["choi"], dtype=complex)
    effects = np.asarray(arrays["effects"], dtype=complex)
    if states.shape != (4, 2, 2) or choi.shape != (4, 4, 4):
        raise ValueError("expected states (4,2,2) and Choi matrices (4,4,4)")

    output_unitary = np.kron(IDENTITY, PAULI_X)
    reflected_states = states.conj()
    reflected_choi = np.asarray(
        [output_unitary @ item.conj() @ output_unitary for item in choi]
    )
    effect_residual = max(
        float(np.linalg.norm(PAULI_X @ item.conj() @ PAULI_X - item))
        for item in effects
    )
    if effect_residual > 1e-10:
        raise ValueError("the supplied terminal POVM is not reflection invariant")

    determinant_before = float(np.linalg.det(input_pauli_matrix(states)))
    determinant_after = float(np.linalg.det(input_pauli_matrix(reflected_states)))
    payload = {
        "source": str(args.input),
        "output": str(args.output),
        "input_determinant": determinant_before,
        "reflected_determinant": determinant_after,
        "determinant_sum_residual": abs(determinant_before + determinant_after),
        "terminal_effect_symmetry_residual": effect_residual,
        "minimum_reflected_choi_eigenvalue": min(
            float(np.linalg.eigvalsh(item).min()) for item in reflected_choi
        ),
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output,
        states=reflected_states,
        choi=reflected_choi,
        effects=effects,
    )
    rendered = json.dumps(payload, indent=2) + "\n"
    print(rendered, end="")
    if args.report is not None:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

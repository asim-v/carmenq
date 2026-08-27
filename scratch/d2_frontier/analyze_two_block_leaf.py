"""Gauge-invariant diagnostics for a general two-block leaf checkpoint."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from qubit_discrimination_geometry import PAULIS, discrimination_geometry


def matrix_data(matrix: np.ndarray) -> dict[str, object]:
    trace = float(np.trace(matrix).real)
    return {
        "trace": trace,
        "eigenvalues": np.linalg.eigvalsh(matrix).real.tolist(),
        "bloch": [
            float(np.trace(matrix @ pauli).real / trace) for pauli in PAULIS
        ]
        if trace > 1e-14
        else [0.0, 0.0, 0.0],
    }


def main() -> None:
    # Checkpoint inspection is the only path that requires the optional PyTorch
    # search stack.  Keep it lazy so the NumPy discrimination helpers remain
    # importable in the reproducibility environment.
    import torch

    from general_two_block_leaf import GeneralTwoBlockLeaf

    parser = argparse.ArgumentParser()
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("--lambda", dest="weight", type=float, default=0.6)
    args = parser.parse_args()

    model = GeneralTwoBlockLeaf(0)
    model.load_state_dict(
        torch.load(args.checkpoint, map_location="cpu", weights_only=True)
    )
    with torch.no_grad():
        columns = model.columns().numpy()
        effects = model.povm().numpy()
        score, audit, returned = model.quotient(args.weight)

    tensor = columns.reshape(4, 4, 2, 4, 4)
    paired = tensor.transpose(0, 3, 1, 2, 4).reshape(16, 32)
    left, singular, right_h = np.linalg.svd(paired, full_matrices=False)
    active = int(np.count_nonzero(singular > 1e-10))

    # Right-canonical one-way representation.  All Schmidt coefficients are
    # absorbed into the prefix tensors; the right rows remain orthonormal and
    # therefore define a trace-preserving four-outcome qubit instrument.
    prefix_amplitudes = (
        left[:, :2] * singular[:2][None, :]
    ).reshape(4, 4, 2)
    prefix_states = np.einsum(
        "bzi,bzj->zij", prefix_amplitudes, prefix_amplitudes.conj()
    )
    right = right_h[:2].reshape(2, 4, 2, 4).transpose(1, 2, 3, 0)
    # Recontracting density matrices is clearer than manipulating amplitudes.
    reconstructed_terminal_states = np.zeros((4, 2, 2), dtype=complex)
    probabilities = np.zeros((4, 4), dtype=float)
    for z in range(4):
        for y in range(4):
            sigma = np.einsum(
                "rmi,ij,rnj->mn",
                right[:, :, y, :],
                prefix_states[z],
                right[:, :, y, :].conj(),
            )
            probabilities[z, y] = np.trace(sigma).real
            reconstructed_terminal_states[z ^ y] += sigma

    direct_terminal = columns.T.reshape(16, 16, 2)
    terminal_states = np.zeros((4, 2, 2), dtype=complex)
    for z in range(4):
        for y in range(4):
            word = 4 * z + y
            terminal_states[z ^ y] += (
                direct_terminal[word].T @ direct_terminal[word].conj()
            )

    gram = columns.conj().T @ columns
    payload = {
        "weight": args.weight,
        "score": float(score),
        "audit": float(audit),
        "return": float(returned),
        "paired_schmidt_rank": active,
        "paired_singular_values": singular[:8].tolist(),
        "gram_off_diagonal_frobenius": float(
            np.linalg.norm(gram - np.diag(np.diag(gram)))
        ),
        "prefix_states": [matrix_data(state) for state in prefix_states],
        "path_probabilities": probabilities.tolist(),
        "terminal_states": [matrix_data(state) for state in terminal_states],
        "terminal_discrimination": discrimination_geometry(terminal_states),
        "one_way_reconstruction_residual": float(
            np.linalg.norm(reconstructed_terminal_states - terminal_states)
        ),
        "terminal_effects": [matrix_data(effect) for effect in effects],
    }
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()

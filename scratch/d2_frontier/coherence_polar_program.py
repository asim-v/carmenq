"""Polar-witness form of the coherence-preserving RETURN functional.

For orthogonal path blocks ``Pi_i`` and a density operator ``rho``, define

    R_coh(rho) = (Tr(rho) + sum_{i != j} ||Pi_i rho Pi_j||_1) / N.

Every off-diagonal trace norm has an independent contraction witness.  This
module constructs the optimal polar witnesses and the resulting Hermitian
operator ``W`` for which

    R_coh(rho) = Tr((I + W) rho) / N.

The identity turns the remaining pure Schmidt-rank frontier into an exact
max-max program: cover contraction witnesses, then bound one linear
Schmidt-number support function in every cell.
"""

from __future__ import annotations

import numpy as np


def _block_slices(block_sizes: tuple[int, ...]) -> tuple[slice, ...]:
    if not block_sizes or any(size < 1 for size in block_sizes):
        raise ValueError("block_sizes must be positive")
    boundaries = np.cumsum((0, *block_sizes))
    return tuple(
        slice(int(boundaries[index]), int(boundaries[index + 1]))
        for index in range(len(block_sizes))
    )


def trace_norm_polar(block: np.ndarray) -> tuple[np.ndarray, float]:
    """Return a contraction attaining ``Re Tr(U^* block)=||block||_1``."""

    matrix = np.asarray(block, dtype=complex)
    if matrix.ndim != 2:
        raise ValueError("block must be a matrix")
    left, singular_values, right_adjoint = np.linalg.svd(
        matrix, full_matrices=False
    )
    if singular_values.size == 0:
        contraction = np.zeros_like(matrix)
    else:
        tolerance = (
            max(matrix.shape)
            * np.finfo(float).eps
            * float(singular_values[0])
        )
        support = singular_values > tolerance
        contraction = left[:, support] @ right_adjoint[support, :]
    return contraction, float(singular_values.sum())


def polar_coherence_witness(
    density: np.ndarray, block_sizes: tuple[int, ...]
) -> tuple[np.ndarray, dict[tuple[int, int], np.ndarray]]:
    """Construct all optimal path-pair contractions and their Hermitian sum."""

    matrix = np.asarray(density, dtype=complex)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("density must be a square matrix")
    blocks = _block_slices(block_sizes)
    if sum(block_sizes) != matrix.shape[0]:
        raise ValueError("block sizes do not match the density dimension")

    witness = np.zeros_like(matrix)
    contractions: dict[tuple[int, int], np.ndarray] = {}
    for first in range(len(blocks)):
        for second in range(first + 1, len(blocks)):
            rows, columns = blocks[first], blocks[second]
            contraction, _ = trace_norm_polar(matrix[rows, columns])
            witness[rows, columns] = contraction
            witness[columns, rows] = contraction.conj().T
            contractions[first, second] = contraction
    return witness, contractions


def coherence_from_witness(density: np.ndarray, witness: np.ndarray) -> float:
    """Evaluate the normalized affine polar-witness functional."""

    matrix = np.asarray(density, dtype=complex)
    polar = np.asarray(witness, dtype=complex)
    if matrix.shape != polar.shape or matrix.ndim != 2:
        raise ValueError("density and witness must be same-size matrices")
    value = np.trace((np.eye(matrix.shape[0]) + polar) @ matrix)
    if abs(value.imag) > 1e-9:
        raise ValueError("witness evaluation is not real")
    return float(value.real)


def contraction_residuals(
    witness: np.ndarray, block_sizes: tuple[int, ...]
) -> dict[str, float]:
    """Audit Hermiticity, zero diagonal blocks, and pair contraction norms."""

    polar = np.asarray(witness, dtype=complex)
    blocks = _block_slices(block_sizes)
    if polar.shape != (sum(block_sizes), sum(block_sizes)):
        raise ValueError("block sizes do not match the witness dimension")
    diagonal_residual = max(
        float(np.linalg.norm(polar[block, block])) for block in blocks
    )
    pair_norms = [
        float(np.linalg.svd(polar[blocks[i], blocks[j]], compute_uv=False)[0])
        for i in range(len(blocks))
        for j in range(i + 1, len(blocks))
    ]
    return {
        "hermiticity": float(np.linalg.norm(polar - polar.conj().T)),
        "diagonal": diagonal_residual,
        "maximum_pair_operator_norm": max(pair_norms, default=0.0),
    }

"""Exact barrier for convexifying one Schmidt-rank-constrained leaf.

For a fixed terminal POVM whose largest eigenvalue on syndrome ``s`` is
``w[s]``, replace the pure two-block leaf by an arbitrary mixed state of
Schmidt number at most two after tracing the emitted registers.  The relaxed
set contains every classical state diagonal in the path labels and is already
large enough to attain the unconstrained probability optimum

    lambda * sum_i w_i p_i
      + (1-lambda) / N * (sum_i sqrt(p_i))**2.

Writing ``x_i=sqrt(p_i)`` turns this expression into a Rayleigh quotient.
Perron--Frobenius makes its leading eigenvector nonnegative, so the optimum is
the largest eigenvalue of one diagonal-plus-rank-one matrix.  Consequently no
convex Schmidt-number relaxation that retains the concave pinched RETURN term
can improve this number; the missing constraint is purity/determinism.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np


def block_trace_coherence_return(
    density: np.ndarray, block_sizes: tuple[int, ...]
) -> float:
    """Evaluate ``(1 + sum_{i!=j} ||rho_ij||_1) / block_count``.

    The normalization matches the pinched RETURN convention when the number
    of path blocks is ``block_count``.  On a pure state this is exactly the
    squared sum of block amplitudes divided by ``block_count``.
    """

    matrix = np.asarray(density, dtype=complex)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("density must be a square matrix")
    if not block_sizes or any(size < 1 for size in block_sizes):
        raise ValueError("block_sizes must be positive")
    if sum(block_sizes) != matrix.shape[0]:
        raise ValueError("block sizes do not match the density dimension")
    boundaries = np.cumsum((0, *block_sizes))
    coherence = 0.0
    for first in range(len(block_sizes)):
        rows = slice(boundaries[first], boundaries[first + 1])
        for second in range(len(block_sizes)):
            if first == second:
                continue
            columns = slice(boundaries[second], boundaries[second + 1])
            coherence += float(
                np.linalg.svd(matrix[rows, columns], compute_uv=False).sum()
            )
    return float((np.trace(matrix).real + coherence) / len(block_sizes))


def pure_pinched_return(vector: np.ndarray, block_sizes: tuple[int, ...]) -> float:
    """Evaluate pinched RETURN directly from one normalized pure vector."""

    state = np.asarray(vector, dtype=complex).reshape(-1)
    if sum(block_sizes) != state.size:
        raise ValueError("block sizes do not match the vector dimension")
    boundaries = np.cumsum((0, *block_sizes))
    amplitude_sum = sum(
        np.linalg.norm(state[boundaries[index] : boundaries[index + 1]])
        for index in range(len(block_sizes))
    )
    return float(amplitude_sum * amplitude_sum / len(block_sizes))


def convexified_leaf_barrier(
    effect_norms: np.ndarray,
    support_weight: float,
    multiplicity: int = 4,
) -> dict[str, Any]:
    """Return the exact mixed-state convexification value and maximizer."""

    weights = np.asarray(effect_norms, dtype=float)
    if weights.ndim != 1 or weights.size < 1:
        raise ValueError("effect_norms must be a nonempty vector")
    if np.any(~np.isfinite(weights)) or np.any(weights < 0.0) or np.any(weights > 1.0):
        raise ValueError("effect norms must lie in [0,1]")
    if not 0.0 <= support_weight <= 1.0:
        raise ValueError("support_weight must lie in [0,1]")
    if multiplicity < 1:
        raise ValueError("multiplicity must be positive")

    expanded = np.repeat(weights, multiplicity)
    path_count = expanded.size
    matrix = (
        support_weight * np.diag(expanded)
        + (1.0 - support_weight)
        * np.ones((path_count, path_count), dtype=float)
        / path_count
    )
    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    vector = np.asarray(eigenvectors[:, -1], dtype=float)
    if vector.sum() < 0.0:
        vector = -vector
    # The matrix is entrywise nonnegative.  Numerical diagonalisation can
    # choose arbitrary signs only on exactly degenerate endpoint faces.
    vector = np.maximum(vector, 0.0)
    vector /= np.linalg.norm(vector)
    probabilities = vector * vector
    audit = float(expanded @ probabilities)
    returned = float(np.square(np.sqrt(probabilities).sum()) / path_count)
    score = float(
        support_weight * audit + (1.0 - support_weight) * returned
    )
    residual = float(np.linalg.norm(matrix @ vector - eigenvalues[-1] * vector))
    grouped = probabilities.reshape(weights.size, multiplicity).sum(axis=1)
    return {
        "support_weight": float(support_weight),
        "effect_norms": weights.tolist(),
        "multiplicity": int(multiplicity),
        "path_count": int(path_count),
        "bound": float(eigenvalues[-1]),
        "score_from_distribution": score,
        "audit": audit,
        "return": returned,
        "path_probabilities": probabilities.tolist(),
        "syndrome_probabilities": grouped.tolist(),
        "eigen_residual": residual,
        "minimum_perron_coordinate": float(vector.min()),
        "interpretation": (
            "attained by a path-diagonal separable mixed state; any convex "
            "Schmidt-number relaxation with concave pinched RETURN is no tighter"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lambda", dest="weight", type=float, default=0.55)
    parser.add_argument(
        "--effect-norms", type=float, nargs="+", default=(0.92, 0.64, 0.44, 0.0)
    )
    parser.add_argument("--multiplicity", type=int, default=4)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = convexified_leaf_barrier(
        np.asarray(args.effect_norms, dtype=float), args.weight, args.multiplicity
    )
    rendered = json.dumps(payload, indent=2) + "\n"
    print(rendered, end="")
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

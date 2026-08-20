"""Verify the stored non-QND counterexample to the analytic QND candidate.

The artifact contains complete local Kraus instruments and terminal AUDIT
effects, not optimizer parameters.  This verifier uses only NumPy and contracts
the four-slot process independently.  RETURN is evaluated with the exact
flag-conditioned trace-norm formula and no spectral floor.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


N = 4
D = 16
ARITY = 3
NODES = 40
TRANSCRIPTS = ARITY**N
CANDIDATE_SCORE = 0.755437446228747


def syndrome(word: int) -> int:
    x1, x2, x3, x4 = ((word >> (3 - index)) & 1 for index in range(N))
    return 2 * (x1 ^ x3) + (x2 ^ x4)


def transcript_digits(index: int) -> tuple[int, ...]:
    digits = []
    for depth in range(N):
        power = ARITY ** (N - 1 - depth)
        digits.append(index // power)
        index %= power
    return tuple(digits)


def node_offset(depth: int) -> int:
    return (ARITY**depth - 1) // (ARITY - 1)


def contract(tree: np.ndarray) -> np.ndarray:
    columns = np.zeros((TRANSCRIPTS, 32, D), dtype=np.complex128)
    for transcript in range(TRANSCRIPTS):
        outcomes = transcript_digits(transcript)
        for word in range(D):
            bits = tuple((word >> (N - 1 - depth)) & 1 for depth in range(N))
            state = np.zeros((1, 2), dtype=np.complex128)
            state[0, 0] = 1.0
            prefix = 0
            for depth in range(N):
                node = node_offset(depth) + prefix
                local = tree[node, outcomes[depth], :, :, bits[depth], :]
                state = np.einsum("bno,po->pbn", local, state).reshape(-1, 2)
                prefix = ARITY * prefix + outcomes[depth]
            columns[transcript, :, word] = state.reshape(32)
    return columns


def evaluate(path: Path) -> dict[str, float | int]:
    with np.load(path, allow_pickle=False) as artifact:
        tree = artifact["kraus_tree"]
        effects = artifact["audit_effects"]
        arity = int(artifact["arity"])
        slots = int(artifact["slots"])
        audit_weight = float(artifact["audit_weight"])
    if tree.shape != (NODES, ARITY, 2, 2, 2, 2):
        raise ValueError(f"unexpected Kraus-tree shape: {tree.shape}")
    if effects.shape != (TRANSCRIPTS, 4, 2, 2):
        raise ValueError(f"unexpected AUDIT-effect shape: {effects.shape}")
    if arity != ARITY or slots != N or audit_weight != 0.5:
        raise ValueError("artifact metadata does not describe the frozen game")

    local = tree.reshape(NODES, 4 * ARITY, 4)
    identity4 = np.eye(4)
    completeness = max(
        np.linalg.norm(matrix.conj().T @ matrix - identity4) for matrix in local
    )
    identity2 = np.eye(2)
    povm_completeness = max(
        np.linalg.norm(branch.sum(axis=0) - identity2) for branch in effects
    )
    povm_minimum_eigenvalue = min(
        np.linalg.eigvalsh(effect).min().real
        for branch in effects
        for effect in branch
    )

    columns = contract(tree)
    returned = sum(
        np.linalg.svd(branch, compute_uv=False).sum() ** 2 for branch in columns
    ) / D**2
    audit = 0.0
    probabilities = np.zeros(TRANSCRIPTS)
    for transcript in range(TRANSCRIPTS):
        probabilities[transcript] = np.linalg.norm(columns[transcript]) ** 2 / D
        for word in range(D):
            vector = columns[transcript, :, word].reshape(16, 2)
            rho_m = vector.conj().T @ vector / D
            audit += np.trace(effects[transcript, syndrome(word)] @ rho_m).real
    score = audit_weight * audit + (1.0 - audit_weight) * returned
    return {
        "audit_probability": float(audit),
        "return_fidelity": float(returned),
        "support_value": float(score),
        "analytic_candidate_value": CANDIDATE_SCORE,
        "strict_excess": float(score - CANDIDATE_SCORE),
        "maximum_local_completeness_residual": float(completeness),
        "maximum_povm_completeness_residual": float(povm_completeness),
        "minimum_povm_eigenvalue": float(povm_minimum_eigenvalue),
        "total_transcript_probability": float(probabilities.sum()),
        "active_transcript_count": int(np.count_nonzero(probabilities > 1e-8)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "artifact",
        nargs="?",
        type=Path,
        default=Path("data/interleaved_ternary_counterexample.npz"),
    )
    args = parser.parse_args()
    result = evaluate(args.artifact)
    print(json.dumps(result, indent=2))
    if result["strict_excess"] <= 0.004:
        raise SystemExit("stored instrument no longer falsifies the candidate")
    if result["maximum_local_completeness_residual"] > 1e-12:
        raise SystemExit("local instrument completeness check failed")
    if result["maximum_povm_completeness_residual"] > 1e-12:
        raise SystemExit("AUDIT POVM completeness check failed")
    if result["minimum_povm_eigenvalue"] < -1e-12:
        raise SystemExit("AUDIT POVM positivity check failed")


if __name__ == "__main__":
    main()

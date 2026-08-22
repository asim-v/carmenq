"""Complete one normalized TT-rank-two Choi leaf by local Pauli twirling.

This script tests the constructive converse behind the homogeneous-leaf
reduction.  A normalized Choi MPS is put in sequential row-isometric gauge.
At slot ``i`` its seed tensor ``T_i`` obeys

    Tr_A(T_i^* T_i) = I_M.

The four local Kraus operators ``T_i P / sqrt(2)``, with ``P`` a one-qubit
Pauli acting on the fresh input, therefore form a complete instrument.  Every
global transcript is a right-Pauli translate of the seed leaf and has the
same homogeneous AUDIT--RETURN score.
"""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np
import torch

from general_single_leaf_bound import GeneralLeaf


D = 16
PAULIS = np.asarray(
    (
        ((1, 0), (0, 1)),
        ((0, 1), (1, 0)),
        ((0, -1j), (1j, 0)),
        ((1, 0), (0, -1)),
    ),
    dtype=complex,
)
X_SHIFT = np.asarray((0, 1, 1, 0), dtype=int)


def syndrome(word: int) -> int:
    bits = tuple((word >> (3 - index)) & 1 for index in range(4))
    return 2 * (bits[0] ^ bits[2]) + (bits[1] ^ bits[3])


def row_isometric_cores(columns: np.ndarray) -> tuple[list[np.ndarray], float]:
    """Return sequential seed tensors in row-isometric gauge."""

    tensor = columns.reshape(2, 2, 2, 2, 2, 2, 2, 2, 2)
    # (b1,b2,b3,b4,m,x1,x2,x3,x4) -> ((b1,x1),...,(b4,x4,m)).
    state = tensor.transpose(0, 5, 1, 6, 2, 7, 3, 8, 4).reshape(1, 4, 4, 4, 8)
    norm = float(np.linalg.norm(state))
    state = state / norm
    physical = (4, 4, 4)
    cores: list[np.ndarray] = []
    old = 1
    remainder = state
    for dimension in physical:
        matrix = remainder.reshape(old * dimension, -1)
        left, singular, right = np.linalg.svd(matrix, full_matrices=False)
        rank = int(np.sum(singular > 1e-11))
        if rank > 2:
            raise ValueError(f"TT rank {rank} exceeds the charged memory")
        core = (left[:, :rank] * singular[:rank]).reshape(old, dimension, rank)
        cores.append(core)
        remainder = right[:rank].reshape(rank, *remainder.shape[2:])
        old = rank
    cores.append(remainder.reshape(old, 2, 2, 2))  # old,b4,x4,m4

    local: list[np.ndarray] = []
    for core in cores[:3]:
        shaped = core.reshape(core.shape[0], 2, 2, core.shape[2])
        local.append(shaped.transpose(1, 3, 2, 0))  # b,new,x,old
    local.append(cores[3].transpose(1, 3, 2, 0))  # b,m4,x,old
    return local, norm


def pauli_instruments(seed: list[np.ndarray]) -> tuple[list[np.ndarray], float]:
    instruments: list[np.ndarray] = []
    residual = 0.0
    for tensor in seed:
        local = np.einsum("bnxo,vxy->v b n y o", tensor, PAULIS) / np.sqrt(2)
        local = local.reshape(4, 2, tensor.shape[1], 2, tensor.shape[3])
        matrix = local.reshape(4 * 2 * tensor.shape[1], 2 * tensor.shape[3])
        identity = np.eye(2 * tensor.shape[3])
        residual = max(residual, float(np.linalg.norm(matrix.conj().T @ matrix - identity)))
        instruments.append(local)
    return instruments, residual


def contract(instruments: list[np.ndarray]) -> np.ndarray:
    transcripts = tuple(itertools.product(range(4), repeat=4))
    columns = np.zeros((len(transcripts), 32, D), dtype=complex)
    for transcript_index, transcript in enumerate(transcripts):
        for word in range(D):
            bits = tuple((word >> (3 - index)) & 1 for index in range(4))
            state = np.ones((1, 1), dtype=complex)
            for depth, (outcome, bit) in enumerate(zip(transcript, bits, strict=True)):
                operator = instruments[depth][outcome, :, :, bit, :]
                state = np.einsum("bno,po->pbn", operator, state).reshape(-1, operator.shape[1])
            columns[transcript_index, :, word] = state.reshape(32)
    return columns


def shifted_syndrome(transcript: tuple[int, ...]) -> int:
    shift_word = sum(
        int(X_SHIFT[outcome]) << (3 - index)
        for index, outcome in enumerate(transcript)
    )
    return syndrome(shift_word)


def scores(columns: np.ndarray, seed_effects: np.ndarray) -> tuple[float, float, float]:
    transcripts = tuple(itertools.product(range(4), repeat=4))
    audit = 0.0
    for transcript_index, transcript in enumerate(transcripts):
        shift = shifted_syndrome(transcript)
        terminal = columns[transcript_index].T.reshape(D, 16, 2)
        for word in range(D):
            state = terminal[word].conj().T @ terminal[word]
            audit += np.trace(seed_effects[syndrome(word) ^ shift] @ state).real / D
    returned = sum(
        np.linalg.svd(branch, compute_uv=False).sum() ** 2 for branch in columns
    ) / D**2
    total = float(np.linalg.norm(columns) ** 2 / D)
    return float(audit), float(returned), total


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("--lambda", dest="weight", type=float, default=0.5)
    args = parser.parse_args()

    model = GeneralLeaf(0)
    model.load_state_dict(torch.load(args.checkpoint, map_location="cpu", weights_only=True))
    with torch.no_grad():
        original = model.columns().numpy()
        effects = model.povm().numpy()
        quotient, audit_leaf, return_leaf = model.quotient(args.weight)

    seed, original_norm = row_isometric_cores(original)
    instruments, completeness = pauli_instruments(seed)
    columns = contract(instruments)
    audit, returned, total = scores(columns, effects)
    score = args.weight * audit + (1 - args.weight) * returned
    payload = {
        "weight": args.weight,
        "original_hilbert_schmidt_norm": original_norm,
        "leaf_score": float(quotient),
        "leaf_audit": float(audit_leaf),
        "leaf_return": float(return_leaf),
        "completed_score": score,
        "completed_audit": audit,
        "completed_return": returned,
        "total_probability": total,
        "maximum_local_completeness_residual": completeness,
        "score_residual": score - float(quotient),
    }
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()

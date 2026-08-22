"""Independent NumPy/SciPy verifier for a fixed-arity n=4 checkpoint."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
import torch


N = 4
D = 16
PAULIS = (
    np.array([[0, 1], [1, 0]], dtype=complex),
    np.array([[0, -1j], [1j, 0]], dtype=complex),
    np.diag([1, -1]).astype(complex),
)


def syndrome(word: int) -> int:
    x1, x2, x3, x4 = ((word >> (3 - index)) & 1 for index in range(N))
    return 2 * (x1 ^ x3) + (x2 ^ x4)


def stiefel(raw: np.ndarray) -> np.ndarray:
    output = np.empty_like(raw)
    for node, matrix in enumerate(raw):
        q, r = np.linalg.qr(matrix, mode="reduced")
        diagonal = np.diag(r)
        phase = diagonal / np.maximum(np.abs(diagonal), 1e-14)
        output[node] = q * phase.conj()[None, :]
    return output


def invsqrt(matrix: np.ndarray) -> np.ndarray:
    values, vectors = np.linalg.eigh(matrix)
    return (vectors * np.maximum(values.real, 1e-13) ** -0.5) @ vectors.conj().T


def transcript_digits(index: int, arity: int) -> tuple[int, ...]:
    digits = []
    for depth in range(N):
        power = arity ** (N - 1 - depth)
        digits.append(index // power)
        index %= power
    return tuple(digits)


def node_offset(arity: int, depth: int) -> int:
    return (arity**depth - 1) // (arity - 1)


def contract(raw_iso: np.ndarray, arity: int) -> tuple[np.ndarray, float]:
    nodes = (arity**N - 1) // (arity - 1)
    transcripts = arity**N
    tree = stiefel(raw_iso).reshape(nodes, arity, 2, 2, 2, 2)
    columns = np.zeros((transcripts, 32, D), dtype=complex)
    for transcript in range(transcripts):
        outcomes = transcript_digits(transcript, arity)
        for word in range(D):
            bits = tuple((word >> (N - 1 - depth)) & 1 for depth in range(N))
            state = np.zeros((1, 2), dtype=complex)
            state[0, 0] = 1.0
            prefix = 0
            for depth in range(N):
                node = node_offset(arity, depth) + prefix
                local = tree[node, outcomes[depth], :, :, bits[depth], :]
                state = np.einsum("bno,po->pbn", local, state).reshape(-1, 2)
                prefix = arity * prefix + outcomes[depth]
            columns[transcript, :, word] = state.reshape(32)
    identity = np.eye(4)
    residual = max(
        np.linalg.norm(matrix.conj().T @ matrix - identity)
        for matrix in tree.reshape(nodes, 4 * arity, 4)
    )
    return columns, float(residual)


def tau_states(columns: np.ndarray) -> np.ndarray:
    tau = np.zeros((columns.shape[0], 4, 2, 2), dtype=complex)
    for transcript in range(columns.shape[0]):
        for word in range(D):
            matrix = columns[transcript, :, word].reshape(16, 2)
            tau[transcript, syndrome(word)] += matrix.conj().T @ matrix / D
    return tau


def qubit_dual(states: np.ndarray) -> float:
    traces = np.trace(states, axis1=1, axis2=2).real
    mass = float(traces.sum())
    if mass < 1e-14:
        return 0.0
    normalized = states / mass
    priors = traces / mass
    vectors = np.array(
        [[np.trace(normalized[s] @ pauli).real for pauli in PAULIS] for s in range(4)]
    )

    # For Y=(t I+y.sigma)/2, dual feasibility is exactly
    # t >= p_s+|y-r_s| for every label.  SLSQP can stop at the centroid of
    # this nonsmooth convex problem and return a visibly loose upper bound.
    # Multi-start Nelder--Mead reliably resolves the two-active-label cases
    # present in the counterexample.  Every trial point is dual feasible once
    # t is evaluated by the maximum below, so the result is still an upper
    # bound even before convergence.
    def objective(point: np.ndarray) -> float:
        return max(
            priors[s] + np.linalg.norm(point - vectors[s]) for s in range(4)
        )
    starts = [vectors.mean(axis=0), *vectors]
    starts.extend(
        (vectors[first] + vectors[second]) / 2
        for first in range(4)
        for second in range(first + 1, 4)
    )
    results = [
        minimize(
            objective,
            start,
            method="Nelder-Mead",
            options={"xatol": 1e-13, "fatol": 1e-14, "maxiter": 20000},
        )
        for start in starts
    ]
    successful = [result for result in results if result.success]
    if not successful:
        raise RuntimeError("qubit dual minimisation failed from every start")
    return mass * min(float(result.fun) for result in successful)


def checkpoint_povms(raw: np.ndarray) -> np.ndarray:
    gram = np.einsum("cskj,cski->csji", raw.conj(), raw)
    effects = np.empty_like(gram)
    for transcript in range(raw.shape[0]):
        normalizer = invsqrt(gram[transcript].sum(axis=0))
        effects[transcript] = np.einsum(
            "ij,sjk,kl->sil", normalizer, gram[transcript], normalizer
        )
    return effects


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("--arity", type=int, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    state = torch.load(args.checkpoint, map_location="cpu", weights_only=True)
    columns, completeness = contract(state["raw_iso"].numpy(), args.arity)
    tau = tau_states(columns)
    effects = checkpoint_povms(state["raw_povm"].numpy())
    returned = sum(
        np.linalg.svd(branch, compute_uv=False).sum() ** 2 for branch in columns
    ) / D**2
    audit_checkpoint = sum(
        np.trace(effects[c, s] @ tau[c, s]).real
        for c in range(columns.shape[0])
        for s in range(4)
    )
    audit_exact = sum(qubit_dual(states) for states in tau)
    probabilities = np.linalg.norm(columns, axis=(1, 2)) ** 2 / D
    active = np.flatnonzero(probabilities > 1e-8)
    per_leaf_audit = np.array([qubit_dual(states) for states in tau])
    payload = {
        "arity": args.arity,
        "audit_with_checkpoint_povm": float(audit_checkpoint),
        "audit_exact_qubit_dual": float(audit_exact),
        "return_fidelity_no_floor": float(returned),
        "score_with_checkpoint_povm": float((audit_checkpoint + returned) / 2),
        "score_with_exact_audit": float((audit_exact + returned) / 2),
        "candidate_score": 0.755437446228747,
        "excess_over_candidate": float((audit_exact + returned) / 2 - 0.755437446228747),
        "maximum_completeness_residual": completeness,
        "active_transcript_count": int(active.size),
        "active_transcripts": [
            {
                "transcript": "".join(map(str, transcript_digits(int(c), args.arity))),
                "probability": float(probabilities[c]),
                "syndrome_likelihood": [
                    float(4 * np.trace(tau[c, s]).real) for s in range(4)
                ],
                "return_contribution": float(
                    np.linalg.svd(columns[c], compute_uv=False).sum() ** 2 / D**2
                ),
                "audit_contribution": float(per_leaf_audit[c]),
                "homogeneous_score": float(
                    (
                        per_leaf_audit[c]
                        + np.linalg.svd(columns[c], compute_uv=False).sum() ** 2
                        / D**2
                    )
                    / (2 * probabilities[c])
                ),
            }
            for c in active
        ],
    }
    rendered = json.dumps(payload, indent=2) + "\n"
    print(rendered)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

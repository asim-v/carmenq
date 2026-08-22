"""Contract the explicit interleaved two-parameter instrument from scratch.

This verifier does not call the closed score formula until after constructing
all 64 terminal memory vectors (four classical leaves for each of 16 inputs).
It therefore checks the physical construction independently of the public
frontier helper.
"""

from __future__ import annotations

import argparse
import json
import math

import numpy as np

from carmenq import interleaved_candidate_scores


N_WORDS = 16


def terminal_vectors(q: float, v: float) -> np.ndarray:
    """Return vectors indexed by ``(c1,c3,x,m)`` for the QND instrument."""
    if not (0.0 <= q <= 1.0 and 0.0 <= v <= 1.0):
        raise ValueError("q and v must lie in [0, 1]")
    p_correct = 1.0 - (1.0 - q) * v**2
    u = q * v**2 / p_correct if p_correct > 0.0 else 0.0
    pauli_z = np.diag([1.0, -1.0])
    correct_filter = np.diag([math.sqrt(q), 1.0])
    wrong_filter = np.diag([math.sqrt(1.0 - q), 0.0])
    output = np.zeros((2, 2, N_WORDS, 2), dtype=complex)

    for c1 in range(2):
        for c3 in range(2):
            for word in range(N_WORDS):
                x1, x2, x3, x4 = (
                    (word >> (3 - index)) & 1 for index in range(4)
                )
                if c1 == x1:
                    state = math.sqrt(p_correct) * np.array(
                        [math.sqrt(u), math.sqrt(1.0 - u)], dtype=complex
                    )
                else:
                    state = math.sqrt(1.0 - p_correct) * np.array(
                        [1.0, 0.0], dtype=complex
                    )
                if x2:
                    state = pauli_z @ state
                state = (
                    correct_filter if c3 == x3 else wrong_filter
                ) @ state
                if x4:
                    state = pauli_z @ state
                output[c1, c3, word] = state
    return output


def direct_scores(q: float, v: float) -> tuple[float, float, float]:
    """Compute AUDIT, RETURN, and the maximum completeness residual."""
    vectors = terminal_vectors(q, v)
    audit_probability = 0.0
    for c1 in range(2):
        for c3 in range(2):
            guessed_first_syndrome = c1 ^ c3
            second_syndrome_states = [
                np.zeros((2, 2), dtype=complex) for _ in range(2)
            ]
            for word in range(N_WORDS):
                x1, x2, x3, x4 = (
                    (word >> (3 - index)) & 1 for index in range(4)
                )
                first_syndrome = x1 ^ x3
                second_syndrome = x2 ^ x4
                if first_syndrome == guessed_first_syndrome:
                    state = vectors[c1, c3, word]
                    second_syndrome_states[second_syndrome] += (
                        np.outer(state, state.conj()) / N_WORDS
                    )
            difference = second_syndrome_states[0] - second_syndrome_states[1]
            audit_probability += 0.5 * (
                np.trace(
                    second_syndrome_states[0] + second_syndrome_states[1]
                ).real
                + np.abs(np.linalg.eigvalsh(difference)).sum()
            )

    norms = np.linalg.norm(vectors, axis=-1)
    return_fidelity = float(np.square(norms.sum(axis=-1)).sum() / N_WORDS**2)
    completeness = np.square(norms).sum(axis=(0, 1))
    residual = float(np.max(np.abs(completeness - 1.0)))
    return float(audit_probability), return_fidelity, residual


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--q", type=float, default=0.6168956030718684)
    parser.add_argument("--v", type=float, default=0.8003177036431812)
    args = parser.parse_args()

    direct_audit, direct_return, completeness_residual = direct_scores(
        args.q, args.v
    )
    closed_audit, closed_return = interleaved_candidate_scores(args.q, args.v)
    payload = {
        "q": args.q,
        "v": args.v,
        "direct": {
            "audit_probability": direct_audit,
            "return_fidelity": direct_return,
        },
        "closed": {
            "audit_probability": closed_audit,
            "return_fidelity": closed_return,
        },
        "maximum_completeness_residual": completeness_residual,
        "maximum_formula_residual": max(
            abs(direct_audit - closed_audit),
            abs(direct_return - closed_return),
        ),
    }
    print(json.dumps(payload, indent=2))
    if completeness_residual > 1e-12 or payload["maximum_formula_residual"] > 1e-12:
        raise SystemExit("candidate verification failed")


if __name__ == "__main__":
    main()

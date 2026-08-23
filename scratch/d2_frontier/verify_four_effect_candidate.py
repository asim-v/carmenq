"""Independent NumPy verifier for the symmetric four-effect MPS branch."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from reduced_four_effect_frontier import optimise  # noqa: E402
from verify_compact_interleaved_candidate import (  # noqa: E402
    PAULIS,
    columns,
    prefix_tensors,
    syndrome,
)


def construction(weight: float) -> tuple[list[np.ndarray], np.ndarray, dict[str, object]]:
    point = optimise(weight)
    p = float(point["p"])
    theta = float(point["theta"])
    u_small = float(point["u_small"])
    u_large = float(point["u_large"])
    priors = tuple(float(value) for value in point["priors"])

    state_small = np.asarray(
        (np.sqrt(u_small), np.sqrt(1.0 - u_small)), dtype=complex
    )
    state_large = np.asarray(
        (np.sqrt(u_large), -np.sqrt(1.0 - u_large)), dtype=complex
    )
    pauli_x = np.asarray(((0, 1), (1, 0)), dtype=complex)
    states = [
        state_small,
        state_large,
        pauli_x @ state_large,
        pauli_x @ state_small,
    ]
    tensors = prefix_tensors(states, priors)

    root_effect = np.diag((np.sqrt(p), np.sqrt(1.0 - p)))
    projector_small = np.asarray((np.sin(theta), np.cos(theta)), dtype=complex)
    projector_large = np.asarray((np.cos(theta), -np.sin(theta)), dtype=complex)
    effect_vectors = [
        root_effect @ projector_small,
        root_effect @ projector_large,
    ]
    effect_vectors.extend((pauli_x @ effect_vectors[1], pauli_x @ effect_vectors[0]))

    third = np.zeros((2, 2, 2, 2), dtype=complex)  # b,new,x,old
    for label, vector in enumerate(effect_vectors):
        x3, decoded_x2 = divmod(label, 2)
        # The emitted carrier records the fresh input x3, while the persistent
        # qubit records the terminal decision x2.  The four (x3,x2) labels are
        # therefore orthogonal without duplicating either label, which makes
        # the Choi columns orthogonal and saturates the pinching bound.
        third[x3, decoded_x2, x3, :] = vector.conj()
    tensors.append(third)

    fourth = np.zeros((2, 2, 2, 2), dtype=complex)
    fourth[0, :, 0, :] = np.eye(2) / np.sqrt(2.0)
    fourth[1, :, 1, :] = pauli_x / np.sqrt(2.0)
    tensors.append(fourth)

    effects = np.zeros((4, 2, 2), dtype=complex)
    effects[0, 0, 0] = 1.0
    effects[1, 1, 1] = 1.0
    return tensors, effects, point


def evaluate(weight: float) -> dict[str, object]:
    tensors, effects, point = construction(weight)
    operator = columns(tensors)
    gram = operator.conj().T @ operator
    terminal = operator.T.reshape(16, 16, 2)
    audit = sum(
        np.trace(
            effects[syndrome(word)]
            @ (terminal[word].conj().T @ terminal[word])
        ).real
        for word in range(16)
    )
    returned = np.linalg.svd(operator, compute_uv=False).sum() ** 2 / 16.0
    score = weight * audit + (1.0 - weight) * returned

    completeness = []
    pauli_completion = []
    for tensor in tensors:
        seed = np.einsum("bnxo,bnxp->op", tensor.conj(), tensor)
        completeness.append(float(np.linalg.norm(seed - np.eye(tensor.shape[-1]))))
        local = np.einsum("bnxo,vxy->vbnyo", tensor, PAULIS) / np.sqrt(2.0)
        local = local.reshape(4, 2, tensor.shape[1], 2, tensor.shape[-1])
        matrix = local.reshape(4 * 2 * tensor.shape[1], 2 * tensor.shape[-1])
        pauli_completion.append(
            float(np.linalg.norm(matrix.conj().T @ matrix - np.eye(2 * tensor.shape[-1])))
        )

    full_tensor = operator.reshape(2, 2, 2, 2, 2, 2, 2, 2, 2).transpose(
        0, 5, 1, 6, 2, 7, 3, 8, 4
    )
    temporal_ranks = [
        int(np.linalg.matrix_rank(full_tensor.reshape(4**cut, -1), tol=1e-11))
        for cut in (1, 2, 3)
    ]
    return {
        "weight": weight,
        "parameters": point,
        "direct_audit": float(audit),
        "direct_return": float(returned),
        "direct_score": float(score),
        "normalisation": float(np.trace(gram).real),
        "gram_off_diagonal_frobenius": float(
            np.linalg.norm(gram - np.diag(np.diag(gram)))
        ),
        "temporal_ranks": temporal_ranks,
        "maximum_row_isometry_residual": max(completeness),
        "maximum_pauli_completion_residual": max(pauli_completion),
        "audit_residual": float(audit - float(point["audit"])),
        "return_residual": float(returned - float(point["return"])),
        "score_residual": float(score - float(point["score"])),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lambda", dest="weight", type=float, default=0.6)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = evaluate(args.weight)
    rendered = json.dumps(payload, indent=2) + "\n"
    print(rendered)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

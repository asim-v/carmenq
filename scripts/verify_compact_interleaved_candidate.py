"""Independently construct and verify the compact interleaved MPS candidate.

The verifier uses NumPy only for the tensor construction.  It builds four
row-isometric local tensors from the analytic three-effect geometry, contracts
the sixteen input columns, checks temporal bond ranks and local completeness,
and evaluates AUDIT and optimal polar RETURN directly.  It also verifies the
four-outcome local Pauli completion identity.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from carmenq.order_sensitive import interleaved_compact_lower_bound  # noqa: E402


PAULIS = np.asarray(
    (
        ((1, 0), (0, 1)),
        ((0, 1), (1, 0)),
        ((0, -1j), (1j, 0)),
        ((1, 0), (0, -1)),
    ),
    dtype=complex,
)


def pure_state(z_coordinate: float, transverse_sign: float = 1.0) -> np.ndarray:
    """Return a real qubit ket with Bloch vector (sign*sqrt(1-z^2),0,z)."""

    z_value = min(1.0, max(-1.0, float(z_coordinate)))
    return np.asarray(
        (
            np.sqrt((1.0 + z_value) / 2.0),
            transverse_sign * np.sqrt((1.0 - z_value) / 2.0),
        ),
        dtype=complex,
    )


def prefix_tensors(states: list[np.ndarray], priors: tuple[float, ...]) -> list[np.ndarray]:
    """Factor the two-site prefix into row-isometric bond-two tensors."""

    tensor = np.zeros((2, 2, 2, 2, 2), dtype=complex)  # b1,b2,m,x1,x2
    for x1 in (0, 1):
        for x2 in (0, 1):
            word = 2 * x1 + x2
            tensor[x1, x2, :, x1, x2] = np.sqrt(priors[word]) * states[word]
    matrix = tensor.transpose(0, 3, 1, 4, 2).reshape(4, 8)
    left, singular, right = np.linalg.svd(matrix, full_matrices=False)
    rank = int(np.sum(singular > 1e-12))
    if rank > 2:
        raise AssertionError(f"prefix Schmidt rank {rank} exceeds two")
    first = (left[:, :rank] * singular[:rank]).reshape(1, 2, 2, rank)
    second = right[:rank].reshape(rank, 2, 2, 2)
    # old, b, x, new -> b,new,x,old
    return [
        first.transpose(1, 3, 2, 0),
        second.transpose(1, 3, 2, 0),
    ]


def construction(weight: float) -> tuple[list[np.ndarray], np.ndarray, dict[str, object]]:
    point = interleaved_compact_lower_bound(weight)
    if point.strategy != "three_effect_mps":
        raise ValueError("the selected support direction uses the no-record branch")
    assert point.t is not None and point.r is not None and point.priors is not None
    t_value = point.t
    r_value = point.r
    other_trace = 1.0 - t_value / 2.0
    effect_z = -t_value / (2.0 - t_value)

    state_reference = pure_state(1.0)
    if t_value >= 0.5:
        null_z = 1.0 / t_value - 1.0
        state_null = pure_state(null_z)
    else:
        state_null = state_reference.copy()
    state_plus = pure_state(r_value, 1.0)
    state_minus = pure_state(r_value, -1.0)
    states = [state_reference, state_null, state_plus, state_minus]
    tensors = prefix_tensors(states, point.priors)

    effect_plus = pure_state(effect_z, 1.0)
    effect_minus = pure_state(effect_z, -1.0)
    third = np.zeros((2, 2, 2, 2), dtype=complex)  # b,new,x,old
    third[0, 0, 0, :] = np.sqrt(t_value) * state_reference.conj()
    third[1, 0, 1, :] = np.sqrt(other_trace) * effect_plus.conj()
    third[0, 1, 1, :] = np.sqrt(other_trace) * effect_minus.conj()
    tensors.append(third)

    pauli_x = np.asarray(((0, 1), (1, 0)), dtype=complex)
    fourth = np.zeros((2, 2, 2, 2), dtype=complex)
    fourth[0, :, 0, :] = np.eye(2) / np.sqrt(2.0)
    fourth[1, :, 1, :] = pauli_x / np.sqrt(2.0)
    tensors.append(fourth)

    terminal_effects = np.zeros((4, 2, 2), dtype=complex)
    terminal_effects[0, 0, 0] = 1.0
    terminal_effects[1, 1, 1] = 1.0
    metadata = {
        "weight": weight,
        "t": t_value,
        "r": r_value,
        "priors": list(point.priors),
        "reported_audit": point.audit_probability,
        "reported_return": point.return_fidelity,
        "reported_score": point.support_value,
    }
    return tensors, terminal_effects, metadata


def columns(tensors: list[np.ndarray]) -> np.ndarray:
    result = np.empty((32, 16), dtype=complex)
    for word in range(16):
        bits = tuple((word >> (3 - index)) & 1 for index in range(4))
        state = tensors[0][:, :, bits[0], 0]
        state = np.einsum("ao,bno->abn", state, tensors[1][:, :, bits[1], :])
        state = np.einsum("abo,cno->abcn", state, tensors[2][:, :, bits[2], :])
        state = np.einsum("abco,dno->abcdn", state, tensors[3][:, :, bits[3], :])
        result[:, word] = state.reshape(32)
    return result


def syndrome(word: int) -> int:
    bits = tuple((word >> (3 - index)) & 1 for index in range(4))
    return 2 * (bits[0] ^ bits[2]) + (bits[1] ^ bits[3])


def evaluate(weight: float = 0.5) -> dict[str, object]:
    tensors, effects, metadata = construction(weight)
    operator = columns(tensors)
    gram = operator.conj().T @ operator
    probabilities = np.diag(gram).real
    terminal = operator.T.reshape(16, 16, 2)
    audit = sum(
        np.trace(effects[syndrome(word)] @ (terminal[word].conj().T @ terminal[word])).real
        for word in range(16)
    )
    returned = np.linalg.svd(operator, compute_uv=False).sum() ** 2 / 16.0
    score = weight * audit + (1.0 - weight) * returned

    completeness_residuals = []
    pauli_residuals = []
    for tensor in tensors:
        seed = np.einsum("bnxo,bnxp->op", tensor.conj(), tensor)
        completeness_residuals.append(float(np.linalg.norm(seed - np.eye(tensor.shape[-1]))))
        local = np.einsum("bnxo,vxy->vbnyo", tensor, PAULIS) / np.sqrt(2.0)
        local = local.reshape(4, 2, tensor.shape[1], 2, tensor.shape[-1])
        matrix = local.reshape(4 * 2 * tensor.shape[1], 2 * tensor.shape[-1])
        pauli_residuals.append(
            float(np.linalg.norm(matrix.conj().T @ matrix - np.eye(2 * tensor.shape[-1])))
        )

    full_tensor = operator.reshape(2, 2, 2, 2, 2, 2, 2, 2, 2).transpose(
        0, 5, 1, 6, 2, 7, 3, 8, 4
    )
    temporal_ranks = [
        int(np.linalg.matrix_rank(full_tensor.reshape(4**cut, -1), tol=1e-11))
        for cut in (1, 2, 3)
    ]
    payload = {
        **metadata,
        "direct_audit": float(audit),
        "direct_return": float(returned),
        "direct_score": float(score),
        "normalisation": float(probabilities.sum()),
        "gram_off_diagonal_frobenius": float(
            np.linalg.norm(gram - np.diag(np.diag(gram)))
        ),
        "temporal_ranks": temporal_ranks,
        "maximum_row_isometry_residual": max(completeness_residuals),
        "maximum_pauli_completion_residual": max(pauli_residuals),
        "audit_residual": float(audit - metadata["reported_audit"]),
        "return_residual": float(returned - metadata["reported_return"]),
        "score_residual": float(score - metadata["reported_score"]),
    }
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lambda", dest="weight", type=float, default=0.5)
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

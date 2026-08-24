"""Exact complete-positive completion test for planar qubit pullbacks.

A nondegenerate three-effect planar qubit POVM spans ``{I, X, Y}``.  Its
three pullbacks through a putative subchannel therefore determine

    F = Phi*(I),  Bx = Phi*(X),  By = Phi*(Y).

There exists a Hermitian ``Bz=Phi*(Z)`` completing these data to a completely
positive map exactly when

    w(F^{-1/2} (Bx + i By) F^{-1/2}) <= 1

on the support of ``F``.  This is Ando's numerical-radius characterization
applied to the 2-by-2 block Choi matrix.  The functions below reconstruct the
operator-system data and evaluate the missing compatibility condition.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import minimize_scalar


IDENTITY = np.eye(2, dtype=complex)
PAULIS = (
    np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex),
    np.array([[0.0, -1j], [1j, 0.0]], dtype=complex),
    np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex),
)
FULL_PAULIS = (IDENTITY, *PAULIS)


def bloch(matrix: np.ndarray) -> np.ndarray:
    """Return the four real Pauli coefficients of a Hermitian qubit matrix."""

    operator = np.asarray(matrix, dtype=complex)
    if operator.shape != (2, 2):
        raise ValueError("operator must be 2 by 2")
    return np.asarray([float(np.trace(operator @ pauli).real) for pauli in FULL_PAULIS])


def from_bloch(coefficients: np.ndarray) -> np.ndarray:
    vector = np.asarray(coefficients, dtype=float)
    if vector.shape != (4,):
        raise ValueError("Bloch coefficients must have length four")
    return sum(vector[mu] * FULL_PAULIS[mu] for mu in range(4)) / 2.0


def planar_reconstruction_matrix(effects: np.ndarray) -> np.ndarray:
    """Return the inverse map from three statistics to ``I,X,Y`` data."""

    terminal = np.asarray(effects, dtype=complex)
    if terminal.ndim != 3 or terminal.shape[1:] != (2, 2):
        raise ValueError("effects must have shape (outcomes,2,2)")
    active = [
        index for index, effect in enumerate(terminal) if np.trace(effect).real > 1e-10
    ]
    if len(active) != 3:
        raise ValueError("exactly three terminal effects must be active")
    matrix = np.asarray(
        [
            [0.5 * np.trace(terminal[index] @ basis).real for basis in FULL_PAULIS[:3]]
            for index in active
        ]
    )
    if abs(float(np.linalg.det(matrix))) < 1e-12:
        raise ValueError("active effects do not span the planar operator system")
    return np.linalg.inv(matrix)


def reconstruct_planar_pullbacks(
    pulled_effects: np.ndarray, effects: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Recover ``Phi*(I)``, ``Phi*(X)``, and ``Phi*(Y)`` from pullbacks."""

    pulled = np.asarray(pulled_effects, dtype=complex)
    terminal = np.asarray(effects, dtype=complex)
    active = [
        index for index, effect in enumerate(terminal) if np.trace(effect).real > 1e-10
    ]
    if pulled.shape == (len(terminal), 2, 2):
        data = pulled[active]
    elif pulled.shape == (3, 2, 2):
        data = pulled
    else:
        raise ValueError("pulled effects have incompatible dimensions")
    inverse = planar_reconstruction_matrix(terminal)
    coefficients = inverse @ np.asarray([bloch(item) for item in data])
    return tuple(from_bloch(coefficients[index]) for index in range(3))


def numerical_radius_witness(
    matrix: np.ndarray, grid_size: int = 4096
) -> tuple[float, float, np.ndarray]:
    """Return numerical radius, a maximizing phase, and a unit vector.

    The returned vector ``u`` satisfies ``abs(u* T u) == w(T)`` up to the
    scalar optimization tolerance. It therefore supplies a directly
    checkable separating witness whenever the radius exceeds one.
    """

    operator = np.asarray(matrix, dtype=complex)
    if operator.ndim != 2 or operator.shape[0] != operator.shape[1]:
        raise ValueError("matrix must be square")
    if operator.shape[0] == 0:
        return 0.0, 0.0, np.zeros(0, dtype=complex)
    if operator.shape == (1, 1):
        phase = float(np.angle(operator[0, 0]))
        return float(abs(operator[0, 0])), phase, np.ones(1, dtype=complex)
    if grid_size < 32:
        raise ValueError("grid_size must be at least 32")

    def support(theta: float) -> float:
        phase = np.exp(-1j * theta)
        hermitian = 0.5 * (phase * operator + phase.conjugate() * operator.conj().T)
        return float(np.linalg.eigvalsh(hermitian)[-1])

    angles = np.linspace(0.0, 2.0 * math.pi, grid_size, endpoint=False)
    values = np.asarray([support(theta) for theta in angles])
    step = 2.0 * math.pi / grid_size
    candidates = np.flatnonzero(
        (values >= np.roll(values, 1)) & (values >= np.roll(values, -1))
    )
    best_index = int(np.argmax(values))
    best = float(values[best_index])
    best_angle = float(angles[best_index])
    for index in candidates:
        center = float(angles[index])
        optimum = minimize_scalar(
            lambda theta: -support(theta),
            bounds=(center - step, center + step),
            method="bounded",
            options={"xatol": 1e-14},
        )
        candidate = float(-optimum.fun)
        if candidate > best:
            best = candidate
            best_angle = float(optimum.x % (2.0 * math.pi))

    phase = np.exp(-1j * best_angle)
    hermitian = 0.5 * (phase * operator + phase.conjugate() * operator.conj().T)
    _, eigenvectors = np.linalg.eigh(hermitian)
    vector = eigenvectors[:, -1]
    witnessed = float(abs(np.vdot(vector, operator @ vector)))
    return witnessed, best_angle, vector


def numerical_radius(matrix: np.ndarray, grid_size: int = 4096) -> float:
    """Evaluate ``max_theta lambda_max(Re(exp(-i theta) T))``."""

    return numerical_radius_witness(matrix, grid_size)[0]


def cp_completion_radius(
    pulled_effects: np.ndarray,
    effects: np.ndarray,
    support_tolerance: float = 1e-10,
) -> dict[str, Any]:
    """Evaluate the exact Ando compatibility radius on the support of ``F``."""

    total, x_pullback, y_pullback = reconstruct_planar_pullbacks(
        pulled_effects, effects
    )
    eigenvalues, eigenvectors = np.linalg.eigh(total)
    support = eigenvalues > support_tolerance
    if not np.any(support):
        residual = max(
            float(np.linalg.norm(operator)) for operator in (x_pullback, y_pullback)
        )
        return {
            "radius": 0.0 if residual <= support_tolerance else math.inf,
            "support_residual": residual,
            "support_rank": 0,
            "compatible": residual <= support_tolerance,
        }

    basis = eigenvectors[:, support]
    projector = basis @ basis.conj().T
    complement = IDENTITY - projector
    support_residual = max(
        float(np.linalg.norm(complement @ operator))
        for operator in (x_pullback, y_pullback)
    )
    inverse_root = np.diag(1.0 / np.sqrt(eigenvalues[support]))
    normalized = inverse_root @ basis.conj().T
    normalized = normalized @ (x_pullback + 1j * y_pullback) @ basis @ inverse_root
    radius, phase, support_vector = numerical_radius_witness(normalized)
    functional_vector = basis @ inverse_root @ support_vector
    functional_vector /= np.linalg.norm(functional_vector)
    functional_state = np.outer(functional_vector, functional_vector.conj())
    total_expectation = float(np.trace(functional_state @ total).real)
    x_expectation = float(np.trace(functional_state @ x_pullback).real)
    y_expectation = float(np.trace(functional_state @ y_pullback).real)
    witness_ratio = math.hypot(x_expectation, y_expectation) / total_expectation
    compatible = support_residual <= support_tolerance and radius <= 1.0 + 1e-10
    return {
        "radius": radius,
        "radius_excess": radius - 1.0,
        "support_residual": support_residual,
        "support_rank": int(np.count_nonzero(support)),
        "minimum_total_eigenvalue": float(eigenvalues.min()),
        "witness_phase": phase,
        "witness_input_bloch": bloch(functional_state)[1:].tolist(),
        "witness_total_expectation": total_expectation,
        "witness_planar_norm": math.hypot(x_expectation, y_expectation),
        "witness_ratio": witness_ratio,
        "compatible": compatible,
    }


def audit_joint_pullbacks(joint: np.ndarray, effects: np.ndarray) -> dict[str, Any]:
    """Audit every outcome family in a pulled-effect array."""

    array = np.asarray(joint, dtype=complex)
    if array.ndim != 4 or array.shape[2:] != (2, 2):
        raise ValueError("joint must have shape (outcomes,terminals,2,2)")
    reports = [cp_completion_radius(array[y], effects) for y in range(array.shape[0])]
    return {
        "outcomes": reports,
        "maximum_radius": float(max(float(item["radius"]) for item in reports)),
        "all_cp_completable": all(bool(item["compatible"]) for item in reports),
        "interpretation": (
            "radius <= 1 for every outcome is equivalent to a shared CP "
            "completion of the observed planar pullbacks"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    arrays = np.load(args.checkpoint)
    payload = audit_joint_pullbacks(arrays["joint"], arrays["effects"])
    rendered = json.dumps(payload, indent=2) + "\n"
    print(rendered, end="")
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

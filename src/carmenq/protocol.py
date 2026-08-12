"""Explicit five-register simulation of reversible quantum histories.

Register order and dimensions are

``B(4) x W(4) x M(4) x G(2) x A(2)``.

``B`` is a four-valued history label, ``W`` a two-bit internal world, ``M`` a
two-bit memory, ``G`` a reversible decision/work bit, and ``A`` the phase-
kickback ancilla.  Four labels are essential: the parity predicate leaves two
histories inside each predicate class, making the residual-record diagnostic
``chi(H:R|P)`` nontrivial.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
from typing import Literal

import numpy as np

from .linalg import (
    amplitude_damping_to_ground_kraus,
    apply_imperfect_permutation,
    apply_local_kraus,
    apply_local_unitary,
    apply_permutation,
    conditional_holevo_information,
    density,
    dephasing_kraus,
    depolarizing_kraus,
    make_permutation,
    partial_trace,
    uniform_overlap_channel,
)

Array = np.ndarray

B, W, M, G, A = range(5)
REGISTER_NAMES = ("B", "W", "M", "G", "A")
DIMENSIONS = (4, 4, 4, 2, 2)
TOTAL_DIMENSION = int(np.prod(DIMENSIONS))


def bit_parity(value: int) -> int:
    """Return the parity of the two-bit history label."""
    return (int(value).bit_count()) & 1


@dataclass(frozen=True)
class NoiseModel:
    """Local Markovian channels applied after each addressed logical gate.

    Probabilities are per logical operation, not calibrated hardware error
    rates.  Amplitude damping is a qudit relaxation-to-zero channel.  The
    environment overlap is a separate, exactly interpretable leakage model.
    Inversion error is a stochastic skipped-inverse model.
    """

    dephasing: float = 0.0
    depolarizing: float = 0.0
    amplitude_damping: float = 0.0
    environment_overlap: float = 1.0
    inversion_error: float = 0.0

    def __post_init__(self) -> None:
        for name in ("dephasing", "depolarizing", "amplitude_damping", "inversion_error"):
            value = getattr(self, name)
            if not 0.0 <= value <= 1.0:
                raise ValueError(f"{name} must lie in [0, 1]")
        if not 0.0 <= self.environment_overlap <= 1.0:
            raise ValueError("environment_overlap must lie in [0, 1]")


@dataclass(frozen=True)
class ProtocolConfig:
    """Configuration for one reference-protocol execution."""

    coherent_input: bool = True
    challenges: tuple[int, ...] = (1, 0)
    uncompute: Literal["full", "leave_memory", "none"] = "full"
    direct_phase: bool = False
    enable_actions: bool = True
    noise: NoiseModel = field(default_factory=NoiseModel)

    def __post_init__(self) -> None:
        if not self.challenges:
            raise ValueError("At least one late challenge is required.")
        if any(challenge not in (0, 1) for challenge in self.challenges):
            raise ValueError("Challenges must be bits.")


@dataclass
class ProtocolResult:
    """Exact final states and scalar diagnostics for one run."""

    state_before_readout: Array
    state_after_readout: Array
    branch_state_before_readout: Array
    branch_state_after_readout: Array
    metrics: dict[str, float]


def _initial_state(coherent: bool, initial_branch: int | None = None) -> Array:
    zero_w = np.eye(4, dtype=complex)[:, 0]
    zero_m = np.eye(4, dtype=complex)[:, 0]
    zero_g = np.eye(2, dtype=complex)[:, 0]
    minus_a = np.array([1.0, -1.0], dtype=complex) / np.sqrt(2.0)

    if initial_branch is not None:
        branch = np.eye(4, dtype=complex)[:, initial_branch]
        return density(np.kron(np.kron(np.kron(np.kron(branch, zero_w), zero_m), zero_g), minus_a))

    if coherent:
        branch = np.ones(4, dtype=complex) / 2.0
        return density(np.kron(np.kron(np.kron(np.kron(branch, zero_w), zero_m), zero_g), minus_a))

    mixture = np.zeros((TOTAL_DIMENSION, TOTAL_DIMENSION), dtype=complex)
    for branch_index in range(4):
        branch = np.eye(4, dtype=complex)[:, branch_index]
        ket = np.kron(np.kron(np.kron(np.kron(branch, zero_w), zero_m), zero_g), minus_a)
        mixture += density(ket) / 4.0
    return mixture


@lru_cache(maxsize=None)
def _world_gate(sign: int) -> Array:
    return make_permutation(
        DIMENSIONS,
        lambda values: (
            values[B],
            (values[W] + sign * values[B]) % 4,
            values[M],
            values[G],
            values[A],
        ),
    )


@lru_cache(maxsize=None)
def _memory_gate(sign: int) -> Array:
    return make_permutation(
        DIMENSIONS,
        lambda values: (
            values[B],
            values[W],
            (values[M] + sign * values[W]) % 4,
            values[G],
            values[A],
        ),
    )


@lru_cache(maxsize=None)
def _decision_gate(round_index: int, challenge: int) -> Array:
    memory_bit = round_index % 2
    return make_permutation(
        DIMENSIONS,
        lambda values: (
            values[B],
            values[W],
            values[M],
            values[G] ^ (((values[M] >> memory_bit) & 1) ^ challenge),
            values[A],
        ),
    )


@lru_cache(maxsize=None)
def _action_gate(round_index: int) -> Array:
    # Each decision conditionally flips one world bit.  Every such action flips
    # world parity, which the later oracle corrects using the decision rule.
    mask = 1 << ((round_index + 1) % 2)
    return make_permutation(
        DIMENSIONS,
        lambda values: (
            values[B],
            values[W] ^ (mask if values[G] else 0),
            values[M],
            values[G],
            values[A],
        ),
    )


def _decision_parity(memory: int, challenges: tuple[int, ...]) -> int:
    answer = 0
    for round_index, challenge in enumerate(challenges):
        answer ^= ((memory >> (round_index % 2)) & 1) ^ challenge
    return answer


@lru_cache(maxsize=None)
def _history_predicate_oracle(challenges: tuple[int, ...]) -> Array:
    # On the intended causal path, parity(W) differs from parity(B) by the XOR
    # of all challenge-responsive decisions.  The oracle reads final internal
    # state, not B, and therefore restores p(B)=two-bit parity.
    return make_permutation(
        DIMENSIONS,
        lambda values: (
            values[B],
            values[W],
            values[M],
            values[G],
            values[A]
            ^ (bit_parity(values[W]) ^ _decision_parity(values[M], challenges)),
        ),
    )


@lru_cache(maxsize=None)
def _direct_parity_oracle() -> Array:
    return make_permutation(
        DIMENSIONS,
        lambda values: (
            values[B],
            values[W],
            values[M],
            values[G],
            values[A] ^ bit_parity(values[B]),
        ),
    )


def _apply_noise(rho: Array, target: int, noise: NoiseModel) -> Array:
    dimension = DIMENSIONS[target]
    if noise.dephasing:
        rho = apply_local_kraus(
            rho, DIMENSIONS, target, dephasing_kraus(dimension, noise.dephasing)
        )
    if noise.depolarizing:
        rho = apply_local_kraus(
            rho, DIMENSIONS, target, depolarizing_kraus(dimension, noise.depolarizing)
        )
    if noise.amplitude_damping:
        rho = apply_local_kraus(
            rho,
            DIMENSIONS,
            target,
            amplitude_damping_to_ground_kraus(dimension, noise.amplitude_damping),
        )
    return rho


def _forward_gate(rho: Array, gate: Array, target: int, noise: NoiseModel) -> Array:
    return _apply_noise(apply_permutation(rho, gate), target, noise)


def _inverse_gate(rho: Array, gate: Array, target: int, noise: NoiseModel) -> Array:
    # The decision and action gates are self-inverse; modular-add gates are
    # passed with the appropriate negative sign by the caller.
    rho = apply_imperfect_permutation(rho, gate, noise.inversion_error)
    return _apply_noise(rho, target, noise)


def _hadamard_four() -> Array:
    hadamard = np.array([[1.0, 1.0], [1.0, -1.0]], dtype=complex) / np.sqrt(2.0)
    return np.kron(hadamard, hadamard)


def _reset_probability(rho: Array) -> float:
    diagonal = np.real(np.diag(rho)).reshape(DIMENSIONS)
    # W=M=G=0; B and A remain unrestricted.
    return float(np.clip(diagonal[:, 0, 0, 0, :].sum(), 0.0, 1.0))


def _coherence_visibility(branch_state: Array) -> float:
    off_diagonal = np.sum(np.abs(branch_state)) - np.sum(np.abs(np.diag(branch_state)))
    return float(np.clip(off_diagonal / 3.0, 0.0, 1.0))


def _run_density(config: ProtocolConfig, initial_branch: int | None = None) -> ProtocolResult:
    rho = _initial_state(config.coherent_input, initial_branch)
    noise = config.noise

    if config.direct_phase:
        rho = _forward_gate(rho, _direct_parity_oracle(), A, noise)
        if noise.environment_overlap != 1.0:
            rho = uniform_overlap_channel(rho, DIMENSIONS, B, noise.environment_overlap)
    else:
        rho = _forward_gate(rho, _world_gate(+1), W, noise)
        rho = _forward_gate(rho, _memory_gate(+1), M, noise)

        for round_index, challenge in enumerate(config.challenges):
            decision = _decision_gate(round_index, challenge)
            rho = _forward_gate(rho, decision, G, noise)
            if config.enable_actions:
                rho = _forward_gate(rho, _action_gate(round_index), W, noise)
            rho = _forward_gate(rho, decision, G, noise)

        if noise.environment_overlap != 1.0:
            rho = uniform_overlap_channel(rho, DIMENSIONS, B, noise.environment_overlap)

        rho = _forward_gate(rho, _history_predicate_oracle(config.challenges), A, noise)

        if config.uncompute != "none":
            for round_index in reversed(range(len(config.challenges))):
                challenge = config.challenges[round_index]
                decision = _decision_gate(round_index, challenge)
                rho = _inverse_gate(rho, decision, G, noise)
                if config.enable_actions:
                    rho = _inverse_gate(rho, _action_gate(round_index), W, noise)
                rho = _inverse_gate(rho, decision, G, noise)

            if config.uncompute == "full":
                rho = _inverse_gate(rho, _memory_gate(-1), M, noise)
            rho = _inverse_gate(rho, _world_gate(-1), W, noise)

    rho_before = (rho + rho.conj().T) / 2.0
    branch_before = partial_trace(rho_before, DIMENSIONS, [B])
    reset = _reset_probability(rho_before)
    visibility = _coherence_visibility(branch_before)

    rho_after = apply_local_unitary(rho_before, DIMENSIONS, B, _hadamard_four())
    rho_after = (rho_after + rho_after.conj().T) / 2.0
    branch_after = partial_trace(rho_after, DIMENSIONS, [B])
    predicate_fidelity = float(np.clip(branch_after[3, 3].real, 0.0, 1.0))
    target_contrast = float((predicate_fidelity - 0.25) / 0.75)

    metrics = {
        "visibility": visibility,
        "reset_fidelity": reset,
        "predicate_fidelity": predicate_fidelity,
        "target_contrast": target_contrast,
        "trace": float(np.trace(rho_after).real),
    }
    return ProtocolResult(rho_before, rho_after, branch_before, branch_after, metrics)


def run_protocol(config: ProtocolConfig | None = None) -> ProtocolResult:
    """Execute the exact density-matrix reference protocol."""
    return _run_density(config or ProtocolConfig())


def conditional_record_information(config: ProtocolConfig | None = None) -> float:
    r"""Return :math:`\chi(H:WMG\mid P)` for branch-conditioned record states.

    It measures history information retained *beyond* the permitted parity
    predicate.  It is not a claim about an optimal hardware measurement; the
    conditional Holevo quantity is an information-theoretic upper bound on
    accessible classical transcript information for this record ensemble.
    """
    config = config or ProtocolConfig()
    record_states = []
    for branch in range(4):
        result = _run_density(config, initial_branch=branch)
        record_states.append(partial_trace(result.state_before_readout, DIMENSIONS, [W, M, G]))
    predicates = [bit_parity(branch) for branch in range(4)]
    return conditional_holevo_information(record_states, predicates)


def environment_conditional_information(overlap: float) -> float:
    r"""Return :math:`\chi(H:E\mid P)` for uniform-overlap pure records.

    The four environment states have Gram matrix with unit diagonal and common
    off-diagonal ``overlap``.  For each two-history parity class this reduces
    to the binary entropy of eigenvalues ``(1 +/- overlap)/2``.
    """
    overlap = float(overlap)
    if not 0.0 <= overlap <= 1.0:
        raise ValueError("overlap must lie in [0, 1]")
    gram = np.full((4, 4), overlap, dtype=float)
    np.fill_diagonal(gram, 1.0)
    eigenvalues, eigenvectors = np.linalg.eigh(gram)
    square_root = (eigenvectors * np.sqrt(np.clip(eigenvalues, 0.0, None))) @ eigenvectors.T
    states = [density(square_root[row]) for row in range(4)]
    predicates = [bit_parity(branch) for branch in range(4)]
    return conditional_holevo_information(states, predicates)

"""Unit tests for trace-preserving dense linear-algebra primitives."""

from __future__ import annotations

import numpy as np

from carmenq.linalg import (
    amplitude_damping_to_ground_kraus,
    apply_local_kraus,
    conditional_holevo_information,
    density,
    dephasing_kraus,
    depolarizing_kraus,
    partial_trace,
)


def test_partial_trace_of_bell_state_is_maximally_mixed() -> None:
    bell = np.array([1.0, 0.0, 0.0, 1.0], dtype=complex) / np.sqrt(2.0)
    reduced = partial_trace(density(bell), (2, 2), [0])
    np.testing.assert_allclose(reduced, np.eye(2) / 2.0, atol=1e-12)


def test_local_channels_are_trace_preserving_and_positive() -> None:
    ket = np.array([0.0, 1.0, 1.0j, 0.0], dtype=complex) / np.sqrt(2.0)
    initial = density(ket)
    channels = (
        dephasing_kraus(2, 0.31),
        depolarizing_kraus(2, 0.31),
        amplitude_damping_to_ground_kraus(2, 0.31),
    )
    for operators in channels:
        output = apply_local_kraus(initial, (2, 2), 0, operators)
        assert np.isclose(np.trace(output).real, 1.0, atol=1e-12)
        assert np.min(np.linalg.eigvalsh((output + output.conj().T) / 2.0)) > -1e-12


def test_conditional_holevo_resolves_information_within_coarse_classes() -> None:
    states = [density(np.eye(4)[:, index]) for index in range(4)]
    parity = [0, 1, 1, 0]
    assert np.isclose(conditional_holevo_information(states, parity), 1.0)

    class_only_states = [states[0], states[1], states[1], states[0]]
    assert np.isclose(conditional_holevo_information(class_only_states, parity), 0.0)

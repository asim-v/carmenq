"""Algebraic and control tests for the four-history reference protocol."""

from __future__ import annotations

import itertools

import numpy as np

from carmenq import (
    NoiseModel,
    ProtocolConfig,
    conditional_record_information,
    environment_conditional_information,
    run_protocol,
)
from carmenq.linalg import partial_trace
from carmenq.protocol import A, DIMENSIONS


def test_ideal_phase_kickback_uncomputation_and_parity_readout() -> None:
    result = run_protocol()
    assert np.isclose(result.metrics["trace"], 1.0)
    assert np.isclose(result.metrics["visibility"], 1.0)
    assert np.isclose(result.metrics["reset_fidelity"], 1.0)
    assert np.isclose(result.metrics["predicate_fidelity"], 1.0)
    assert np.isclose(conditional_record_information(), 0.0, atol=1e-12)

    expected_before = np.outer(
        np.array([1.0, -1.0, -1.0, 1.0]) / 2.0,
        np.array([1.0, -1.0, -1.0, 1.0]) / 2.0,
    )
    np.testing.assert_allclose(result.branch_state_before_readout, expected_before, atol=1e-12)
    np.testing.assert_allclose(np.diag(result.branch_state_after_readout), [0, 0, 0, 1], atol=1e-12)
    ancilla = partial_trace(result.state_before_readout, DIMENSIONS, [A])
    minus = np.array([1.0, -1.0]) / np.sqrt(2.0)
    np.testing.assert_allclose(ancilla, np.outer(minus, minus), atol=1e-12)


def test_late_challenge_choices_preserve_the_intended_predicate() -> None:
    for length in (1, 2, 3, 4):
        for challenges in itertools.product((0, 1), repeat=length):
            result = run_protocol(ProtocolConfig(challenges=challenges))
            assert np.isclose(result.metrics["predicate_fidelity"], 1.0, atol=1e-12)
            assert np.isclose(result.metrics["reset_fidelity"], 1.0, atol=1e-12)


def test_action_bypass_breaks_the_history_dependent_oracle_relation() -> None:
    result = run_protocol(ProtocolConfig(enable_actions=False))
    assert result.metrics["predicate_fidelity"] < 0.26
    assert np.isclose(result.metrics["reset_fidelity"], 1.0)


def test_direct_phase_is_endpoint_equivalent_to_the_ideal_history() -> None:
    ideal = run_protocol()
    direct = run_protocol(ProtocolConfig(direct_phase=True))
    np.testing.assert_allclose(
        direct.branch_state_before_readout, ideal.branch_state_before_readout, atol=1e-12
    )
    np.testing.assert_allclose(
        direct.branch_state_after_readout, ideal.branch_state_after_readout, atol=1e-12
    )
    assert np.isclose(direct.metrics["reset_fidelity"], 1.0)


def test_classical_mixture_has_no_interference_signal() -> None:
    result = run_protocol(ProtocolConfig(coherent_input=False))
    assert np.isclose(result.metrics["visibility"], 0.0, atol=1e-12)
    assert np.isclose(result.metrics["predicate_fidelity"], 0.25, atol=1e-12)
    assert np.isclose(result.metrics["target_contrast"], 0.0, atol=1e-12)


def test_retained_memory_removes_interference_and_keeps_one_extra_bit() -> None:
    config = ProtocolConfig(uncompute="leave_memory")
    result = run_protocol(config)
    assert np.isclose(result.metrics["visibility"], 0.0, atol=1e-12)
    assert np.isclose(result.metrics["predicate_fidelity"], 0.25, atol=1e-12)
    assert np.isclose(result.metrics["reset_fidelity"], 0.25, atol=1e-12)
    assert np.isclose(conditional_record_information(config), 1.0, atol=1e-12)


def test_environment_overlap_sets_visibility_and_conditional_record_tradeoff() -> None:
    overlaps = np.linspace(0.0, 1.0, 6)
    visibilities = []
    information_values = []
    for overlap in overlaps:
        config = ProtocolConfig(noise=NoiseModel(environment_overlap=float(overlap)))
        result = run_protocol(config)
        visibilities.append(result.metrics["visibility"])
        information_values.append(environment_conditional_information(float(overlap)))
        assert np.isclose(result.metrics["visibility"], overlap, atol=1e-12)
        assert np.isclose(
            result.metrics["predicate_fidelity"], 0.25 + 0.75 * overlap, atol=1e-12
        )
        assert np.isclose(result.metrics["reset_fidelity"], 1.0, atol=1e-12)
    assert np.all(np.diff(visibilities) >= -1e-12)
    assert np.all(np.diff(information_values) <= 1e-12)
    assert np.isclose(information_values[0], 1.0)
    assert np.isclose(information_values[-1], 0.0, atol=1e-12)


def test_dephasing_monotonically_reduces_visibility() -> None:
    levels = (0.0, 0.01, 0.03, 0.06)
    visibility = [
        run_protocol(ProtocolConfig(noise=NoiseModel(dephasing=level))).metrics["visibility"]
        for level in levels
    ]
    assert np.all(np.diff(visibility) <= 1e-12)
    assert visibility[-1] < visibility[0]


def test_each_reference_channel_monotonically_reduces_readout_quality() -> None:
    levels = (0.0, 0.02, 0.05)
    for channel in ("dephasing", "depolarizing", "amplitude_damping"):
        results = [
            run_protocol(
                ProtocolConfig(noise=NoiseModel(**{channel: level}))
            ).metrics
            for level in levels
        ]
        for metric in ("visibility", "predicate_fidelity"):
            values = [result[metric] for result in results]
            assert np.all(np.diff(values) <= 1e-12), (channel, metric, values)


def test_imperfect_inverse_lowers_reset_and_readout_fidelity() -> None:
    ideal = run_protocol().metrics
    imperfect = run_protocol(
        ProtocolConfig(noise=NoiseModel(inversion_error=0.08))
    ).metrics
    assert imperfect["reset_fidelity"] < ideal["reset_fidelity"]
    assert imperfect["predicate_fidelity"] < ideal["predicate_fidelity"]

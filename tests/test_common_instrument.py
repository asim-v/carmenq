from __future__ import annotations

import numpy as np
import pytest

from carmenq.common_instrument import (
    apply_choi,
    choi_from_kraus,
    comparison_scale_grid,
    conditioned_outputs,
    flagged_trace_norm_cut,
    project_to_common_instrument,
    reconstruct_common_instrument_from_basis,
    reconstruct_effective_povm_from_basis,
    robust_common_instrument_witness_bound,
    scan_flagged_trace_norm_cuts,
)


IDENTITY = np.eye(2, dtype=complex)
X = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
Y = np.array([[0.0, -1j], [1j, 0.0]], dtype=complex)
Z = np.diag([1.0, -1.0]).astype(complex)


def random_prefix_states(seed: int = 20260823) -> np.ndarray:
    generator = np.random.default_rng(seed)
    roots = generator.normal(size=(4, 2, 2)) + 1j * generator.normal(
        size=(4, 2, 2)
    )
    states = np.einsum("zai,zaj->zij", roots.conj(), roots)
    states /= np.trace(states, axis1=1, axis2=2).real.sum()
    return states


def basis_prefix_states() -> np.ndarray:
    return np.asarray(
        [
            0.125 * (IDENTITY + 0.6 * Z),
            0.125 * (IDENTITY - 0.6 * Z),
            0.125 * (IDENTITY + 0.6 * X),
            0.125 * (IDENTITY + 0.6 * Y),
        ]
    )


def test_identity_choi_convention() -> None:
    state = np.asarray([[0.3, 0.1j], [-0.1j, 0.7]], dtype=complex)
    identity_choi = choi_from_kraus((IDENTITY,))
    assert apply_choi(identity_choi, state) == pytest.approx(state)


def test_flagged_data_processing_accepts_a_physical_instrument() -> None:
    states = random_prefix_states()
    probabilities = np.asarray([0.1, 0.2, 0.3, 0.4])
    choi = np.asarray(
        [choi_from_kraus((np.sqrt(probability) * IDENTITY,)) for probability in probabilities]
    )
    outputs = conditioned_outputs(states, choi)
    cuts = scan_flagged_trace_norm_cuts(states, outputs, scales=(0.0, 0.3, 1.0, 2.0))
    assert min(cut.slack for cut in cuts) >= -2e-12
    assert comparison_scale_grid()[0] == 0.0
    assert 1.0 in comparison_scale_grid()


def test_flagged_data_processing_rejects_distinct_outputs_for_equal_inputs() -> None:
    states = np.repeat((IDENTITY / 8.0)[None, :, :], 4, axis=0)
    outputs = np.zeros((4, 4, 2, 2), dtype=complex)
    outputs[:, 0] = np.diag([0.125, 0.125])
    outputs[0, 0] = np.diag([0.25, 0.0])
    outputs[1, 0] = np.diag([0.0, 0.25])
    cut = flagged_trace_norm_cut(states, outputs, 0, 1, scale=1.0)
    assert cut.input_trace_norm == pytest.approx(0.0)
    assert cut.output_trace_norm == pytest.approx(0.5)
    assert cut.violation == pytest.approx(0.5)


def test_exact_choi_projection_separates_an_incompatible_family() -> None:
    pytest.importorskip("cvxpy")
    states = np.repeat((IDENTITY / 8.0)[None, :, :], 4, axis=0)
    outputs = np.zeros((4, 4, 2, 2), dtype=complex)
    outputs[:, 0] = np.diag([0.125, 0.125])
    outputs[0, 0] = np.diag([0.25, 0.0])
    outputs[1, 0] = np.diag([0.0, 0.25])
    projection = project_to_common_instrument(states, outputs)
    assert projection.distance > 0.1
    assert projection.separation_gap > 0.1
    assert not projection.is_compatible()
    assert projection.trace_preservation_residual < 2e-7
    assert projection.minimum_choi_eigenvalue > -2e-7
    assert projection.uniform_input_radius_budget > 0.05


def test_exact_choi_projection_recovers_a_physical_family() -> None:
    pytest.importorskip("cvxpy")
    states = random_prefix_states(8)
    probabilities = np.asarray([0.15, 0.2, 0.25, 0.4])
    choi = np.asarray(
        [choi_from_kraus((np.sqrt(probability) * IDENTITY,)) for probability in probabilities]
    )
    outputs = conditioned_outputs(states, choi)
    projection = project_to_common_instrument(states, outputs)
    assert projection.distance < 2e-7
    assert projection.is_compatible(2e-7)
    assert projection.trace_preservation_residual < 2e-7


def test_projection_witness_extends_to_perturbed_prefix_states() -> None:
    pytest.importorskip("cvxpy")
    reference = np.repeat((IDENTITY / 8.0)[None, :, :], 4, axis=0)
    incompatible = np.zeros((4, 4, 2, 2), dtype=complex)
    incompatible[:, 0] = np.diag([0.125, 0.125])
    incompatible[0, 0] = np.diag([0.25, 0.0])
    incompatible[1, 0] = np.diag([0.0, 0.25])
    projection = project_to_common_instrument(reference, incompatible)

    perturbed = reference.copy()
    perturbed[0] += np.diag([0.01, -0.01])
    probabilities = np.asarray([0.1, 0.2, 0.3, 0.4])
    choi = np.asarray(
        [choi_from_kraus((np.sqrt(probability) * IDENTITY,)) for probability in probabilities]
    )
    physical_outputs = conditioned_outputs(perturbed, choi)
    physical_witness_value = float(
        sum(
            np.trace(
                projection.witness[z, y].conj().T @ physical_outputs[z, y]
            ).real
            for z in range(4)
            for y in range(4)
        )
    )
    upper = robust_common_instrument_witness_bound(
        projection.compatible_support_value,
        projection.witness,
        reference,
        perturbed,
    )
    assert physical_witness_value <= upper + 2e-7


def test_basis_reconstruction_recovers_unique_physical_instrument() -> None:
    states = basis_prefix_states()
    probabilities = np.asarray([0.1, 0.2, 0.3, 0.4])
    choi = np.asarray(
        [choi_from_kraus((np.sqrt(probability) * IDENTITY,)) for probability in probabilities]
    )
    outputs = conditioned_outputs(states, choi)
    reconstruction = reconstruct_common_instrument_from_basis(states, outputs)
    assert reconstruction.is_compatible(2e-12)
    assert reconstruction.choi_matrices == pytest.approx(choi, abs=2e-12)
    assert reconstruction.output_residual <= 2e-15
    assert reconstruction.trace_preservation_residual <= 2e-15
    assert reconstruction.input_condition_number < 10.0
    assert reconstruction.signed_choi_numerators == pytest.approx(
        abs(reconstruction.input_determinant) * reconstruction.choi_matrices,
        abs=2e-14,
    )


def test_basis_reconstruction_detects_positive_but_incompatible_outputs() -> None:
    states = basis_prefix_states()
    probabilities = np.full(4, 0.25)
    choi = np.asarray(
        [choi_from_kraus((0.5 * IDENTITY,)) for _ in range(4)]
    )
    outputs = conditioned_outputs(states, choi)
    outputs[3, 0] += 0.004 * Z
    outputs[3, 1] -= 0.004 * Z
    assert min(np.linalg.eigvalsh(item).min() for row in outputs for item in row) >= -1e-12
    reconstruction = reconstruct_common_instrument_from_basis(states, outputs)
    assert reconstruction.output_residual <= 2e-15
    assert reconstruction.trace_preservation_residual <= 2e-15
    assert reconstruction.minimum_choi_eigenvalue < -1e-3
    assert not reconstruction.is_compatible()


def test_basis_reconstruction_rejects_singular_input_family() -> None:
    states = np.repeat((IDENTITY / 8.0)[None, :, :], 4, axis=0)
    outputs = np.zeros((4, 2, 2, 2), dtype=complex)
    with pytest.raises(np.linalg.LinAlgError, match="do not span"):
        reconstruct_common_instrument_from_basis(states, outputs)


def test_effective_povm_basis_reconstruction_recovers_a_physical_povm() -> None:
    states = basis_prefix_states()
    effects = np.repeat((IDENTITY / 12.0)[None, :, :], 12, axis=0)
    probabilities = np.asarray(
        [[np.trace(state @ effect).real for effect in effects] for state in states]
    )
    reconstruction = reconstruct_effective_povm_from_basis(states, probabilities)
    assert reconstruction.is_compatible(2e-12)
    assert reconstruction.effect_matrices == pytest.approx(effects, abs=2e-14)
    assert reconstruction.probability_residual < 2e-16
    assert reconstruction.completeness_residual < 2e-15
    assert reconstruction.signed_effect_numerators == pytest.approx(
        abs(reconstruction.input_determinant) * effects,
        abs=2e-15,
    )


def test_effective_povm_basis_reconstruction_detects_negative_unique_effects() -> None:
    states = basis_prefix_states()
    effects = np.repeat((0.09 * IDENTITY)[None, :, :], 12, axis=0)
    effects[0] = 0.05 * IDENTITY + 0.08 * X
    effects[1] = 0.05 * IDENTITY - 0.08 * X
    probabilities = np.asarray(
        [[np.trace(state @ effect).real for effect in effects] for state in states]
    )
    assert np.min(probabilities) >= 0.0
    reconstruction = reconstruct_effective_povm_from_basis(states, probabilities)
    assert reconstruction.probability_residual < 2e-16
    assert reconstruction.completeness_residual < 2e-15
    assert reconstruction.minimum_effect_eigenvalue == pytest.approx(-0.03)
    assert np.count_nonzero(reconstruction.minimum_effect_eigenvalues < -1e-12) == 2
    assert not reconstruction.is_compatible()


def test_effective_povm_basis_reconstruction_rejects_singular_inputs() -> None:
    states = np.repeat((IDENTITY / 8.0)[None, :, :], 4, axis=0)
    with pytest.raises(np.linalg.LinAlgError, match="do not span"):
        reconstruct_effective_povm_from_basis(states, np.full((4, 3), 1.0 / 24.0))

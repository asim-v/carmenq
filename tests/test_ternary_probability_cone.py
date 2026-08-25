from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scratch" / "d2_frontier"))

from pairwise_inellipse_box_cover import cross_coefficient  # noqa: E402
from common_effective_povm_audit import audit_common_effective_povm  # noqa: E402
from fourier_behavior_cap_cover import cube_face_caps  # noqa: E402
from ternary_probability_cone_cover import (  # noqa: E402
    TernaryConeOracle,
    choi_probability_coefficients,
    pauli_effect_matrix,
    pauli_state_matrix,
    projective_comparison_bonus,
    terminal_weight_intervals,
    terminal_weights,
)
from terminal_reconstruction_enclosure import (  # noqa: E402
    planar_effect_pauli,
    planar_reconstruction,
    reconstruction_anchor_and_errors,
    terminal_effect_anchor_and_errors,
)
from ternary_common_instrument_input_cover import (  # noqa: E402
    MeasuredInstrumentProjectionOracle,
    _coefficient_parameters,
    instrument_tube_data,
    robust_witness_error,
)
from ternary_bilinear_instrument_input_cover import (  # noqa: E402
    ando_witness_coefficient_bounds,
    ando_input_split_scores,
    determinant_ando_witnesses,
    planar_ando_direction,
    box_purity_caps,
    determinant_interval,
    determinant_povm_witnesses,
    determinant_split_scores,
    determinant_vertex_bounds,
    replacement_determinant_bounds,
)
from carmenq.common_instrument import apply_choi, choi_from_kraus  # noqa: E402
from validate_terminal_reconstructed_frontier import validate  # noqa: E402
from validate_spectral_cap_cluster_summary import (  # noqa: E402
    validate as validate_spectral_cap_cluster_summary,
)
from validate_common_effective_povm_frontier import (  # noqa: E402
    validate as validate_common_povm_frontier,
)
from ternary_frontier_separator_cover import open_source_cells  # noqa: E402
from spectral_product_localizer_batch import enclosing_scaled_cap  # noqa: E402
from ternary_extend_separator_frontier import normalise_frontier  # noqa: E402
from ternary_refine_last_separator_cover import cube_face_children  # noqa: E402
from ternary_shared_separator_cover import extract_shared_frontier  # noqa: E402
from summarize_ando_instrument_cover import (  # noqa: E402
    SCHEMA as ANDO_SUMMARY_SCHEMA,
)

from summarize_determinant_povm_cover import (  # noqa: E402
    SCHEMA as DETERMINANT_SUMMARY_SCHEMA,
    validate_accounting as validate_determinant_accounting,
)


def canonical_three_effect_povm(weights: np.ndarray) -> np.ndarray:
    """Construct the planar rank-one POVM without importing the Torch probes."""

    w0, w1, w2 = map(float, weights)
    cosine = np.clip(
        (w2 * w2 - w0 * w0 - w1 * w1) / (2.0 * w0 * w1),
        -1.0,
        1.0,
    )
    sine = math.sqrt(max(0.0, 1.0 - float(cosine) ** 2))
    vectors = np.asarray(
        [
            [1.0, 0.0, 0.0],
            [cosine, sine, 0.0],
            [-(w0 + w1 * cosine) / w2, -(w1 * sine) / w2, 0.0],
        ]
    )
    identity = np.eye(2, dtype=complex)
    paulis = np.asarray(
        [
            [[0.0, 1.0], [1.0, 0.0]],
            [[0.0, -1j], [1j, 0.0]],
            [[1.0, 0.0], [0.0, -1.0]],
        ],
        dtype=complex,
    )
    effects = np.asarray(
        [
            0.5
            * weights[index]
            * (identity + np.tensordot(vectors[index], paulis, axes=1))
            for index in range(3)
        ]
    )
    assert np.linalg.norm(effects.sum(axis=0) - identity) < 2e-10
    return effects


def random_positive_operator(rng: np.random.Generator) -> np.ndarray:
    matrix = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
    return matrix @ matrix.conj().T


def test_horwitz_weights_and_interval_enclosure() -> None:
    rng = np.random.default_rng(20260823)
    for _ in range(200):
        beta = rng.uniform(1.0, 1.25)
        alpha = rng.uniform(beta, 2.0)
        weights = np.asarray(terminal_weights(alpha, beta))
        assert np.isclose(weights.sum(), 2.0, atol=2e-15)
        assert weights[0] >= weights[1] >= weights[2] >= 0.0

        half_alpha = rng.uniform(0.0, min(0.05, alpha - 1.0))
        half_beta = rng.uniform(0.0, min(0.03, beta - 1.0))
        box = {
            "terminal_alpha": (alpha - half_alpha, alpha + half_alpha),
            "terminal_beta": (beta - half_beta, beta + half_beta),
        }
        intervals = terminal_weight_intervals(box)
        for weight, (lower, upper) in zip(weights, intervals, strict=True):
            assert lower - 2e-15 <= weight <= upper + 2e-15


def test_probability_range_obeys_homogenized_inellipse() -> None:
    rng = np.random.default_rng(314159)
    for _ in range(100):
        beta = rng.uniform(1.01, 1.35)
        alpha = rng.uniform(max(beta, 1.01), 1.95)
        weights = np.asarray(terminal_weights(alpha, beta))
        effects = canonical_three_effect_povm(weights)
        operator = random_positive_operator(rng)
        q = np.asarray([np.trace(operator @ effects[t]).real for t in range(3)])
        scale = float(q.sum())
        x, y = float(q[0]), float(q[1])
        polynomial = (
            beta**2 * x**2
            + alpha**2 * y**2
            + cross_coefficient(alpha, beta) * x * y
            - 2.0 * beta * x * scale
            - 2.0 * alpha * y * scale
            + scale**2
        )
        assert polynomial <= 2e-10 * max(1.0, scale**2)


def test_aligned_projective_comparison_for_random_states() -> None:
    rng = np.random.default_rng(271828)
    for _ in range(100):
        beta = rng.uniform(1.01, 1.25)
        alpha = rng.uniform(beta, 1.8)
        weights = np.asarray(terminal_weights(alpha, beta))
        effects = canonical_three_effect_povm(weights)
        states = [random_positive_operator(rng) for _ in range(3)]
        priors = np.asarray([np.trace(state).real for state in states])
        ternary = sum(
            np.trace(effects[label] @ states[label]).real for label in range(3)
        )
        for retained in range(3):
            projector = effects[retained] / weights[retained]
            for complement in range(3):
                if complement == retained:
                    continue
                deleted = 3 - retained - complement
                projective = (
                    np.trace(projector @ states[retained]).real
                    + np.trace((np.eye(2) - projector) @ states[complement]).real
                )
                bonus = projective_comparison_bonus(
                    weights, priors, retained, complement
                )
                assert ternary - projective <= bonus + 2e-12
                assert deleted not in {retained, complement}


def test_planar_reconstruction_recovers_visible_bloch_coordinates() -> None:
    rng = np.random.default_rng(161803)
    paulis = np.asarray(
        [
            [[0.0, 1.0], [1.0, 0.0]],
            [[0.0, -1j], [1j, 0.0]],
            [[1.0, 0.0], [0.0, -1.0]],
        ],
        dtype=complex,
    )
    identity = np.eye(2, dtype=complex)
    for _ in range(100):
        beta = rng.uniform(1.02, 1.30)
        alpha = rng.uniform(max(beta, 1.05), 1.95)
        effects = canonical_three_effect_povm(np.asarray(terminal_weights(alpha, beta)))
        coefficients = rng.normal(size=4)
        operator = 0.5 * (
            coefficients[0] * identity + np.tensordot(coefficients[1:], paulis, axes=1)
        )
        measured = np.asarray([np.trace(effect @ operator).real for effect in effects])
        assert np.allclose(measured.sum(), coefficients[0], atol=2e-12)
        assert np.allclose(
            planar_reconstruction(alpha, beta) @ measured,
            coefficients[1:3],
            atol=3e-12,
        )


def test_outward_reconstruction_column_errors_cover_random_box_points() -> None:
    rng = np.random.default_rng(141421)
    for _ in range(30):
        beta_center = rng.uniform(1.04, 1.24)
        alpha_center = rng.uniform(beta_center + 0.03, 1.90)
        alpha_half = rng.uniform(1e-5, 5e-3)
        beta_half = rng.uniform(1e-5, 3e-3)
        alpha_bounds = (alpha_center - alpha_half, alpha_center + alpha_half)
        beta_bounds = (beta_center - beta_half, beta_center + beta_half)
        anchor, errors, _ = reconstruction_anchor_and_errors(alpha_bounds, beta_bounds)
        for _ in range(30):
            alpha = rng.uniform(*alpha_bounds)
            beta = rng.uniform(*beta_bounds)
            difference = planar_reconstruction(alpha, beta) - anchor
            assert np.all(np.linalg.norm(difference, axis=0) <= errors + 2e-15)
            signed = rng.normal(size=3)
            assert np.linalg.norm(difference @ signed) <= (
                errors @ np.abs(signed) + 2e-14
            )


def test_terminal_effect_operator_errors_cover_random_box_points() -> None:
    rng = np.random.default_rng(173205)
    for _ in range(20):
        beta_center = rng.uniform(1.05, 1.22)
        alpha_center = rng.uniform(beta_center + 0.05, 1.88)
        alpha_half = rng.uniform(1e-5, 3e-3)
        beta_half = rng.uniform(1e-5, 2e-3)
        alpha_bounds = (alpha_center - alpha_half, alpha_center + alpha_half)
        beta_bounds = (beta_center - beta_half, beta_center + beta_half)
        anchor, errors, norm_upper, _ = terminal_effect_anchor_and_errors(
            alpha_bounds, beta_bounds
        )
        anchor_matrices = [pauli_effect_matrix(row) for row in anchor]
        assert np.linalg.norm(sum(anchor_matrices) - np.eye(2)) < 2e-13
        for _ in range(20):
            effects = planar_effect_pauli(
                rng.uniform(*alpha_bounds), rng.uniform(*beta_bounds)
            )
            for t in range(3):
                matrix = pauli_effect_matrix(effects[t])
                assert (
                    np.linalg.norm(matrix - anchor_matrices[t], ord=2)
                    <= errors[t] + 3e-14
                )
                assert np.linalg.norm(matrix, ord=2) <= norm_upper[t] + 3e-14


def test_shared_choi_probability_coefficients_match_direct_evaluation() -> None:
    rng = np.random.default_rng(223607)
    input_pauli = np.asarray(
        [
            [0.30, 0.06, 0.01, -0.02],
            [0.27, -0.03, 0.05, 0.01],
            [0.23, 0.02, -0.04, 0.03],
            [0.20, -0.01, -0.02, -0.04],
        ]
    )
    effects = planar_effect_pauli(1.75, 1.16)
    coefficients = choi_probability_coefficients(input_pauli, effects)
    operator = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
    choi = choi_from_kraus([operator])
    for z in range(4):
        state = pauli_state_matrix(input_pauli[z])
        output = apply_choi(choi, state)
        for t in range(3):
            direct = float(np.trace(pauli_effect_matrix(effects[t]) @ output).real)
            linear = float(np.sum(coefficients[z, t] * choi).real)
            assert np.isclose(direct, linear, atol=3e-13)


def test_robust_shared_instrument_tube_contains_physical_perturbations() -> None:
    rng = np.random.default_rng(244949)
    anchor = np.asarray(
        [
            [0.30, 0.06, 0.01, -0.02],
            [0.27, -0.03, 0.05, 0.01],
            [0.23, 0.02, -0.04, 0.03],
            [0.20, -0.01, -0.02, -0.04],
        ]
    )
    coordinate_radii = np.full((4, 4), 0.003)
    lower = anchor - coordinate_radii
    upper = anchor + coordinate_radii
    alpha_bounds = (1.748, 1.752)
    beta_bounds = (1.158, 1.162)
    terminal_anchor, terminal_errors, terminal_norm_upper, _ = (
        terminal_effect_anchor_and_errors(alpha_bounds, beta_bounds)
    )
    tube = instrument_tube_data(
        lower,
        upper,
        terminal_anchor,
        terminal_errors,
        terminal_norm_upper,
    )

    raw_kraus = [
        rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)) for _ in range(4)
    ]
    normalizer = sum(item.conj().T @ item for item in raw_kraus)
    eigenvalues, eigenvectors = np.linalg.eigh(normalizer)
    inverse_root = (eigenvectors / np.sqrt(eigenvalues)) @ eigenvectors.conj().T
    kraus = [item @ inverse_root for item in raw_kraus]
    choi = [choi_from_kraus([item]) for item in kraus]
    assert (
        np.linalg.norm(sum(item.conj().T @ item for item in kraus) - np.eye(2)) < 2e-13
    )

    actual_pauli = rng.uniform(lower, upper)
    actual_effects = planar_effect_pauli(
        rng.uniform(*alpha_bounds), rng.uniform(*beta_bounds)
    )
    q0 = np.empty((4, 4, 3), dtype=float)
    q = np.empty_like(q0)
    for z in range(4):
        state0 = pauli_state_matrix(tube["anchor"][z])
        state = pauli_state_matrix(actual_pauli[z])
        for y in range(4):
            output0 = apply_choi(choi[y], state0)
            output = apply_choi(choi[y], state)
            for t in range(3):
                q0[z, y, t] = np.trace(
                    pauli_effect_matrix(terminal_anchor[t]) @ output0
                ).real
                q[z, y, t] = np.trace(
                    pauli_effect_matrix(actual_effects[t]) @ output
                ).real
                assert abs(q[z, y, t] - q0[z, y, t]) <= (
                    tube["probability_radii"][z, t] * np.trace(choi[y]).real + 3e-13
                )
        assert np.sum(np.abs(q[z] - q0[z])) <= tube["row_radii"][z] + 3e-13


def test_common_instrument_tube_builds_one_realified_shared_choi_family() -> None:
    inputs = np.asarray(
        [
            [0.30, 0.06, 0.01, -0.02],
            [0.27, -0.03, 0.05, 0.01],
            [0.23, 0.02, -0.04, 0.03],
            [0.20, -0.01, -0.02, -0.04],
        ]
    )
    coefficients = choi_probability_coefficients(
        inputs, planar_effect_pauli(1.75, 1.16)
    )
    oracle = TernaryConeOracle(
        0.55,
        (0, 1, 2, 3),
        (),
        (),
        0.79,
        0.7573,
        common_instrument_probability_coefficients=coefficients,
        common_instrument_probability_radii=np.full((4, 3), 1e-3),
        common_instrument_row_radii=np.full(4, 2e-3),
    )
    assert len(oracle.common_instrument_choi) == 4
    assert all(
        real.shape == (4, 4) and imaginary.shape == (4, 4)
        for real, imaginary in oracle.common_instrument_choi
    )
    assert oracle.problem.is_dpp()
    assert not oracle.problem.is_mixed_integer()


def test_common_product_localizers_hold_for_exact_physical_products() -> None:
    rng = np.random.default_rng(224745)
    raw_kraus = [
        rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)) for _ in range(4)
    ]
    normalizer = sum(item.conj().T @ item for item in raw_kraus)
    eigenvalues, eigenvectors = np.linalg.eigh(normalizer)
    inverse_root = (eigenvectors / np.sqrt(eigenvalues)) @ eigenvectors.conj().T
    choi = [choi_from_kraus([item @ inverse_root]) for item in raw_kraus]

    lower, coordinate, upper = -0.17, 0.06, 0.21
    products = [coordinate * matrix for matrix in choi]
    for matrix, product in zip(choi, products, strict=True):
        for residual in (
            product - lower * matrix,
            upper * matrix - product,
        ):
            realification = np.block(
                [[residual.real, -residual.imag], [residual.imag, residual.real]]
            )
            assert np.linalg.eigvalsh(realification).min() >= -2e-14

    partial_trace = np.empty((2, 2), dtype=complex)
    for i in range(2):
        for j in range(2):
            partial_trace[i, j] = sum(
                product[2 * i, 2 * j] + product[2 * i + 1, 2 * j + 1]
                for product in products
            )
    np.testing.assert_allclose(
        partial_trace,
        coordinate * np.eye(2),
        atol=3e-15,
    )


def test_state_choi_tensor_localizer_is_exact_and_positive() -> None:
    rng = np.random.default_rng(20260825)
    raw_kraus = [
        rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)) for _ in range(4)
    ]
    normalizer = sum(item.conj().T @ item for item in raw_kraus)
    eigenvalues, eigenvectors = np.linalg.eigh(normalizer)
    inverse_root = (eigenvectors / np.sqrt(eigenvalues)) @ eigenvectors.conj().T
    choi = [choi_from_kraus([item @ inverse_root]) for item in raw_kraus]
    input_pauli = np.asarray([0.31, 0.08, -0.04, 0.05])
    state = pauli_state_matrix(input_pauli)

    for matrix in choi:
        products = [coordinate * matrix for coordinate in input_pauli]
        tensor_localizer = np.block(
            [
                [products[0] + products[3], products[1] - 1j * products[2]],
                [products[1] + 1j * products[2], products[0] - products[3]],
            ]
        )
        np.testing.assert_allclose(
            tensor_localizer,
            2.0 * np.kron(state, matrix),
            atol=3e-15,
        )
        realification = np.block(
            [
                [tensor_localizer.real, -tensor_localizer.imag],
                [tensor_localizer.imag, tensor_localizer.real],
            ]
        )
        assert np.linalg.eigvalsh(realification).min() >= -2e-14
        transposed_localizer = np.block(
            [
                [products[0] + products[3], products[1] + 1j * products[2]],
                [products[1] - 1j * products[2], products[0] - products[3]],
            ]
        )
        np.testing.assert_allclose(
            transposed_localizer,
            2.0 * np.kron(state.T, matrix),
            atol=3e-15,
        )
        transposed_realification = np.block(
            [
                [transposed_localizer.real, -transposed_localizer.imag],
                [transposed_localizer.imag, transposed_localizer.real],
            ]
        )
        assert np.linalg.eigvalsh(transposed_realification).min() >= -2e-14


def test_enclosing_spectral_cap_proves_child_cap_containment() -> None:
    indices = (89, 90, 93, 94)
    scaled, audit = enclosing_scaled_cap(4, indices)
    center = np.asarray(audit["normal"], dtype=float)
    cosine = float(audit["cosine"])
    radius = float(audit["angular_radius"])
    np.testing.assert_allclose(scaled, center / cosine, atol=2e-15)
    assert 0.0 < cosine < cube_face_caps(4)[0][1]
    for index in indices:
        normal, child_cosine = cube_face_caps(4)[index]
        center_distance = math.acos(float(np.clip(center @ normal, -1.0, 1.0)))
        child_radius = math.acos(float(child_cosine))
        assert center_distance + child_radius <= radius + 2e-15


def test_bilinear_oracle_builds_dpp_common_product_localizers() -> None:
    center = np.asarray(
        [
            [0.30, 0.06, 0.01, -0.02],
            [0.27, -0.03, 0.05, 0.01],
            [0.23, 0.02, -0.04, 0.03],
            [0.20, -0.01, -0.02, -0.04],
        ]
    )
    lower = center - 0.01
    upper = center + 0.01
    terminal, errors, _, _ = terminal_effect_anchor_and_errors(
        (1.74, 1.76),
        (1.15, 1.17),
    )
    common = dict(
        support_weight=0.55,
        prefix_order=(0, 1, 2, 3),
        pairs=(),
        coordinate_cases=(),
        maximum_weight_floor=0.79,
        projective_support_upper=0.7573,
        input_pauli_lower=lower,
        input_pauli_upper=upper,
        common_povm_bilinear=True,
        common_instrument_terminal_effect_anchor=terminal,
        common_instrument_terminal_effect_errors=errors,
    )
    baseline = TernaryConeOracle(**common)
    strengthened = TernaryConeOracle(
        **common,
        common_povm_product_sum_rules=True,
        common_instrument_product_trace_rules=True,
        common_instrument_product_psd_sandwiches=True,
        common_instrument_product_state_choi_psd=True,
        common_instrument_product_state_choi_ppt=True,
    )
    assert strengthened.problem.is_dpp()
    assert not strengthened.problem.is_mixed_integer()
    assert len(strengthened.common_instrument_state_choi_localizers) == 32
    assert len(strengthened.problem.constraints) >= len(baseline.problem.constraints) + 250


def test_measured_instrument_projection_accepts_a_physical_table() -> None:
    rng = np.random.default_rng(264575)
    inputs = np.asarray(
        [
            [0.30, 0.06, 0.01, -0.02],
            [0.27, -0.03, 0.05, 0.01],
            [0.23, 0.02, -0.04, 0.03],
            [0.20, -0.01, -0.02, -0.04],
        ]
    )
    effects = planar_effect_pauli(1.75, 1.16)
    raw_kraus = [
        rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)) for _ in range(4)
    ]
    normalizer = sum(item.conj().T @ item for item in raw_kraus)
    eigenvalues, eigenvectors = np.linalg.eigh(normalizer)
    inverse_root = (eigenvectors / np.sqrt(eigenvalues)) @ eigenvectors.conj().T
    choi = [choi_from_kraus([item @ inverse_root]) for item in raw_kraus]
    coefficients = choi_probability_coefficients(inputs, effects)
    statistics = np.asarray(
        [
            [
                [np.sum(coefficients[z, t] * choi[y]).real for t in range(3)]
                for y in range(4)
            ]
            for z in range(4)
        ]
    )
    parameters = _coefficient_parameters()
    for z in range(4):
        for t in range(3):
            parameters[z][t][0].value = coefficients[z, t].real
            parameters[z][t][1].value = coefficients[z, t].imag
    projection = MeasuredInstrumentProjectionOracle(parameters).project(
        statistics, 1e-9, 2e-7
    )
    assert projection["status"] in {"optimal", "optimal_inaccurate"}
    assert projection["compatible"]
    assert projection["distance"] < 2e-7


def test_witness_specific_robust_error_contains_physical_motion() -> None:
    rng = np.random.default_rng(282843)
    anchor = np.asarray(
        [
            [0.30, 0.06, 0.01, -0.02],
            [0.27, -0.03, 0.05, 0.01],
            [0.23, 0.02, -0.04, 0.03],
            [0.20, -0.01, -0.02, -0.04],
        ]
    )
    coordinate_radii = np.full((4, 4), 0.002)
    lower = anchor - coordinate_radii
    upper = anchor + coordinate_radii
    alpha_bounds = (1.749, 1.751)
    beta_bounds = (1.159, 1.161)
    terminal_anchor, terminal_errors, terminal_norm_upper, _ = (
        terminal_effect_anchor_and_errors(alpha_bounds, beta_bounds)
    )
    tube = instrument_tube_data(
        lower,
        upper,
        terminal_anchor,
        terminal_errors,
        terminal_norm_upper,
    )
    witness = rng.normal(size=(4, 4, 3))
    error, _ = robust_witness_error(witness, tube, terminal_anchor, terminal_errors)
    raw_kraus = [
        rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)) for _ in range(4)
    ]
    normalizer = sum(item.conj().T @ item for item in raw_kraus)
    eigenvalues, eigenvectors = np.linalg.eigh(normalizer)
    inverse_root = (eigenvectors / np.sqrt(eigenvalues)) @ eigenvectors.conj().T
    choi = [choi_from_kraus([item @ inverse_root]) for item in raw_kraus]
    actual_inputs = rng.uniform(lower, upper)
    actual_effects = planar_effect_pauli(
        rng.uniform(*alpha_bounds), rng.uniform(*beta_bounds)
    )
    anchor_coefficients = choi_probability_coefficients(tube["anchor"], terminal_anchor)
    actual_coefficients = choi_probability_coefficients(actual_inputs, actual_effects)
    anchor_statistics = np.asarray(
        [
            [
                [np.sum(anchor_coefficients[z, t] * choi[y]).real for t in range(3)]
                for y in range(4)
            ]
            for z in range(4)
        ]
    )
    actual_statistics = np.asarray(
        [
            [
                [np.sum(actual_coefficients[z, t] * choi[y]).real for t in range(3)]
                for y in range(4)
            ]
            for z in range(4)
        ]
    )
    motion = float(np.sum(witness * (actual_statistics - anchor_statistics)))
    assert abs(motion) <= error + 5e-13


def test_terminal_reconstructed_frontier_artifacts_are_self_consistent() -> None:
    summary = validate()
    assert summary["logical_status"] == (
        "strict local strengthening; target remains open"
    )
    assert summary["global_status"] == "continuous terminal strip still open"


def test_spectral_cap_cluster_summary_is_self_consistent() -> None:
    summary = validate_spectral_cap_cluster_summary()
    assert summary["selected_base_angular_cell_closed"] is True
    assert summary["aggregate_upper_bound"] < summary["target"]
    assert summary["cluster_cover"]["closed_source_open_cells"] == 2216
    assert summary["cluster_cover"]["unresolved_nodes"] == 0


def test_spectral_cover_builds_one_mixed_integer_selector_family() -> None:
    caps = tuple(np.append(normal, cosine) for normal, cosine in cube_face_caps(2))
    anchor, errors, _ = reconstruction_anchor_and_errors((1.92, 1.93), (1.14, 1.15))
    oracle = TernaryConeOracle(
        0.55,
        (0, 1, 2, 3),
        (),
        (),
        0.79,
        0.7573,
        projective_support_lines=((0.6, 0.76591),),
        common_contractions=(
            {
                "coefficients": [1.0, 1.0, -1.0, -1.0],
                "branch": "spectral-cover",
                "scalar_signs": (1,),
                "caps": caps,
            },
        ),
        terminal_reconstruction=(anchor, errors),
    )
    assert oracle.problem.is_mixed_integer()
    assert len(oracle.common_contraction_selectors) == 1
    assert oracle.common_contraction_selectors[0].shape == (25,)


def test_fixed_common_effective_povm_builds_one_joint_effect_family() -> None:
    fixed_input = np.asarray(
        [
            [0.25, 0.10, 0.00, 0.00],
            [0.25, 0.00, 0.10, 0.00],
            [0.25, 0.00, 0.00, 0.10],
            [0.25, 0.00, 0.00, 0.00],
        ]
    )
    oracle = TernaryConeOracle(
        0.55,
        (0, 1, 2, 3),
        (),
        (),
        0.79,
        0.7573,
        fixed_common_povm_input=fixed_input,
    )
    assert oracle.effective_povm is not None
    assert oracle.effective_povm.shape == (12, 4)
    assert len(oracle.input_vectors) == 4
    assert not oracle.problem.is_mixed_integer()


def test_common_effective_povm_neighbourhood_uses_valid_probability_envelope() -> None:
    rng = np.random.default_rng(20260824)
    radii = rng.uniform(0.0, 0.02, size=(4, 4))
    for z in range(4):
        for _ in range(50):
            a0 = rng.uniform(0.0, 1.0)
            direction = rng.normal(size=3)
            direction /= np.linalg.norm(direction)
            effect = np.append(a0, rng.uniform(0.0, a0) * direction)
            delta = rng.uniform(-radii[z], radii[z])
            assert abs(float(delta @ effect)) <= np.sum(radii[z]) * a0 + 1e-15

    anchor = np.asarray(
        [
            [0.25, 0.10, 0.00, 0.00],
            [0.25, 0.00, 0.10, 0.00],
            [0.25, 0.00, 0.00, 0.10],
            [0.25, 0.00, 0.00, 0.00],
        ]
    )
    oracle = TernaryConeOracle(
        0.55,
        (0, 1, 2, 3),
        (),
        (),
        0.79,
        0.7573,
        common_povm_input_anchor=anchor,
        common_povm_input_radii=np.full((4, 4), 1e-3),
    )
    assert oracle.effective_povm is not None
    assert not oracle.input_vectors


def test_shared_separator_extracts_all_open_first_generation_branches() -> None:
    payload = {
        "nodes": [
            {
                "status": "optimal",
                "separator": {"coefficients": [1.0, 0.0, -1.0, 0.0]},
                "children": [
                    {"branch": "scalar-positive", "cap": None, "bound": 0.74},
                    {"branch": "scalar-negative", "cap": None, "bound": 0.76},
                    {"branch": "bloch", "cap": 7, "bound": 0.759},
                    {"branch": "bloch", "cap": 8, "bound": 0.75},
                ],
            },
            {
                "status": "optimal",
                "separator": {"coefficients": [0.0, 1.0, 0.0, -1.0]},
                "children": [],
            },
        ]
    }
    first, second, parents, counts = extract_shared_frontier(payload, 0.758)
    assert np.array_equal(first, [1.0, 0.0, -1.0, 0.0])
    assert np.array_equal(second, [0.0, 1.0, 0.0, -1.0])
    assert parents == (
        {"branch": "scalar-negative", "cap": None},
        {"branch": "bloch", "cap": 7},
    )
    assert counts == {
        "first_generation_total": 4,
        "first_generation_closed": 2,
        "first_generation_open": 2,
    }


def test_frontier_filter_preserves_scalar_and_bloch_product_cells() -> None:
    payload = {
        "cells": [
            {
                "first_branch": "scalar-negative",
                "first_cap": None,
                "second_branch": "bloch",
                "second_cap": 3,
                "bound": 0.759,
            },
            {
                "first_branch": "bloch",
                "first_cap": 1,
                "second_branch": "scalar-positive",
                "second_cap": None,
                "bound": 0.757,
            },
        ]
    }
    assert open_source_cells(payload, 0.758) == (payload["cells"][0],)


def test_nested_cube_face_children_partition_each_coarse_chart_cell() -> None:
    all_children = [
        child for parent in range(24) for child in cube_face_children(parent, 2, 4)
    ]
    assert len(all_children) == 96
    assert sorted(all_children) == list(range(96))
    assert cube_face_children(0, 2, 4) == (0, 1, 4, 5)
    assert cube_face_children(23, 2, 4) == (90, 91, 94, 95)


def test_arbitrary_depth_frontier_normalises_legacy_three_separator_cells() -> None:
    payload = {
        "contraction_grid": 4,
        "new_grid": 2,
        "first_separator": [1.0, 0.0, -1.0, 0.0],
        "shared_second_separator": [0.0, 1.0, 0.0, -1.0],
        "new_separator": [1.0, -1.0, 0.0, 0.0],
        "cells": [
            {
                "first_branch": "bloch",
                "first_cap": 5,
                "second_branch": "scalar-positive",
                "second_cap": None,
                "new_branch": "bloch",
                "new_cap": 8,
                "status": "optimal",
                "bound": 0.759,
            },
            {
                "first_branch": "bloch",
                "first_cap": 6,
                "second_branch": "bloch",
                "second_cap": 7,
                "new_branch": "scalar-negative",
                "new_cap": None,
                "status": "optimal",
                "bound": 0.75,
            },
        ],
    }
    coefficients, grids, cells = normalise_frontier(payload, 0.758)
    assert len(coefficients) == 3
    assert grids == (4, 4, 2)
    assert len(cells) == 1
    assert cells[0]["branches"] == ("bloch", "scalar-positive", "bloch")
    assert cells[0]["caps"] == (5, None, 8)


def test_common_effective_povm_basis_audit_accepts_a_physical_povm() -> None:
    priors = np.full(4, 0.25)
    bloch = np.asarray(
        [[0.1, 0.0, 0.0], [0.0, 0.1, 0.0], [0.0, 0.0, 0.1], [0.0, 0.0, 0.0]]
    )
    state_matrix = np.column_stack([priors, bloch])
    effects = np.zeros((4, 12))
    effects[0, :] = 1.0 / 12.0
    statistics = (state_matrix @ effects).reshape(4, 4, 3)
    audit = audit_common_effective_povm(priors, bloch, statistics)
    assert audit["status"] == "nonsingular"
    assert audit["common_effective_povm"]
    assert audit["negative_effect_count"] == 0
    assert abs(audit["minimum_margin"] - 1.0 / 12.0) < 2e-14
    assert audit["interpolation_residual"] < 1e-14
    assert audit["completeness_residual"] < 1e-14


def test_common_effective_povm_basis_audit_exposes_negative_unique_effects() -> None:
    priors = np.full(4, 0.25)
    bloch = np.asarray(
        [[0.1, 0.0, 0.0], [0.0, 0.1, 0.0], [0.0, 0.0, 0.1], [0.0, 0.0, 0.0]]
    )
    state_matrix = np.column_stack([priors, bloch])
    effects = np.zeros((4, 12))
    effects[:, 0] = [0.05, 0.08, 0.0, 0.0]
    effects[:, 1] = [0.05, -0.08, 0.0, 0.0]
    effects[0, 2:] = 0.09
    statistics = (state_matrix @ effects).reshape(4, 4, 3)
    assert np.min(statistics) >= 0.0
    audit = audit_common_effective_povm(priors, bloch, statistics)
    assert not audit["common_effective_povm"]
    assert audit["negative_effect_count"] == 2
    assert abs(audit["minimum_margin"] + 0.03) < 2e-14
    assert abs(audit["worst_effect_witness_expectation"] + 0.03) < 2e-14
    assert audit["completeness_residual"] < 1e-14


def test_common_effective_povm_frontier_artifacts_recompute_exactly() -> None:
    summary = validate_common_povm_frontier()
    assert summary["logical_status"] == (
        "adaptive separator frontier rejected by kill criterion"
    )
    assert summary["source_open_cells"] == 815
    assert summary["depth4_open_cells"] == 2216
    assert summary["negative_effect_count"] == 10
    assert summary["certified_row_l1_radius"] == 0.0871
    assert summary["certified_neighbourhood_bound"] < 0.758


def test_pure_prefix_caps_enclose_every_box_corner() -> None:
    lower = np.asarray(
        [
            [0.20, 0.08, -0.02, 0.03],
            [0.22, -0.03, 0.10, 0.02],
            [0.24, 0.01, -0.04, 0.12],
            [0.18, -0.09, -0.03, 0.01],
        ]
    )
    upper = lower + np.asarray([0.02, 0.01, 0.015, 0.012])
    caps = box_purity_caps(lower, upper)
    for row in range(4):
        normal = caps[row, :3]
        cosine = caps[row, 3]
        assert np.linalg.norm(normal) > 0.999999999
        for bits in range(8):
            corner = np.asarray(
                [
                    upper[row, coordinate + 1]
                    if bits & (1 << coordinate)
                    else lower[row, coordinate + 1]
                    for coordinate in range(3)
                ]
            )
            assert normal @ corner + 2e-16 >= cosine * np.linalg.norm(corner)


def test_row_replacement_determinant_bounds_are_vertex_exact() -> None:
    center = np.asarray(
        [
            [0.25, 0.10, 0.00, 0.00],
            [0.25, 0.00, 0.10, 0.00],
            [0.25, 0.00, 0.00, 0.10],
            [0.25, 0.00, 0.00, 0.00],
        ]
    )
    lower = center - 0.003
    upper = center + 0.003
    replacement = np.asarray([1.0, -1.0, 0.0, 0.0])
    certified_lower, certified_upper = replacement_determinant_bounds(
        lower, upper, 2, replacement, -1
    )
    values = []
    free_rows = (0, 1, 3)
    for choices in np.ndindex(16, 16, 16):
        matrix = center.copy()
        matrix[2] = replacement
        for position, row in enumerate(free_rows):
            bits = choices[position]
            matrix[row] = [
                upper[row, coordinate]
                if bits & (1 << coordinate)
                else lower[row, coordinate]
                for coordinate in range(4)
            ]
        values.append(-float(np.linalg.det(matrix)))
    assert certified_lower <= min(values)
    assert certified_upper >= max(values)
    assert certified_lower > min(values) - 2e-14
    assert certified_upper < max(values) + 2e-14


def test_full_vertex_determinant_bound_removes_interval_dependency() -> None:
    center = np.asarray(
        [
            [0.25, 0.10, 0.00, 0.00],
            [0.25, 0.00, 0.10, 0.00],
            [0.25, 0.00, 0.00, 0.10],
            [0.25, 0.00, 0.00, 0.00],
        ]
    )
    lower = center - 0.015
    upper = center + 0.015
    dependency_bound = determinant_interval(lower, upper)
    vertex_bound = determinant_vertex_bounds(lower, upper)
    assert dependency_bound.lower < 0.0 < dependency_bound.upper
    assert vertex_bound.upper < -9.9e-6
    assert vertex_bound.lower < vertex_bound.upper


def test_determinant_witness_rejects_a_nonpositive_common_povm_robustly() -> None:
    inputs = np.asarray(
        [
            [0.25, 0.10, 0.00, 0.00],
            [0.25, 0.00, 0.10, 0.00],
            [0.25, 0.00, 0.00, 0.10],
            [0.25, 0.00, 0.00, 0.00],
        ]
    )
    effects = np.zeros((4, 12))
    effects[:, 0] = [0.05, 0.08, 0.0, 0.0]
    effects[:, 1] = [0.05, -0.08, 0.0, 0.0]
    effects[0, 2:] = 0.09
    statistics = (inputs @ effects).reshape(4, 4, 3)
    witnesses, audit = determinant_povm_witnesses(
        inputs - 1e-6,
        inputs + 1e-6,
        inputs,
        statistics,
        2e-12,
    )
    assert audit["sign_definite"]
    assert audit["negative_effect_count"] == 2
    assert audit["robust_witness_count"] == 2
    assert {item["effect_index"] for item in witnesses} == {0, 1}
    assert min(item["violation"] for item in witnesses) > 7e-6


def test_determinant_margin_branching_targets_a_useful_bisection() -> None:
    inputs = np.asarray(
        [
            [0.25, 0.10, 0.00, 0.00],
            [0.25, 0.00, 0.10, 0.00],
            [0.25, 0.00, 0.00, 0.10],
            [0.25, 0.00, 0.00, 0.00],
        ]
    )
    effects = np.zeros((4, 12))
    effects[:, 0] = [0.05, 0.08, 0.0, 0.0]
    effects[:, 1] = [0.05, -0.08, 0.0, 0.0]
    effects[0, 2:] = 0.09
    statistics = (inputs @ effects).reshape(4, 4, 3)
    lower = inputs - 0.003
    upper = inputs + 0.003
    witnesses, audit = determinant_povm_witnesses(
        lower, upper, inputs, statistics, 2e-12
    )
    assert not witnesses
    assert audit["minimum_robust_lhs"] > 0.0
    scores = determinant_split_scores(
        lower, upper, statistics, audit, maximum_effects=1
    )
    assert np.max(scores) > 0.1


def test_determinant_cover_tree_accounting_recomputes_from_records() -> None:
    payload = {
        "records": [
            {
                "identifier": 0,
                "disposition": "split",
                "branching_rule": "product-residual",
                "new_witnesses": 1,
                "determinant_audit": {"sign_definite": True},
            },
            {
                "identifier": 1,
                "disposition": "closed",
                "new_witnesses": 0,
            },
        ],
        "pending": [{"identifier": 2, "parent_bound": 0.76}],
        "unresolved": [],
        "solved_nodes": 2,
        "closed_nodes": 1,
        "split_nodes": 1,
        "pending_nodes": 1,
        "unresolved_nodes": 0,
        "determinant_witness_count": 1,
        "maximum_pending_bound": 0.76,
    }
    audit = validate_determinant_accounting(payload)
    assert audit["dispositions"] == {"closed": 1, "split": 1}
    assert audit["branching_rules"] == {"none": 1, "product-residual": 1}
    assert audit["sign_definite_audits"] == 1


def test_compact_determinant_cover_summary_records_open_status() -> None:
    path = (
        ROOT
        / "scratch"
        / "d2_frontier"
        / "determinant_povm_cover_l055_summary.json"
    )
    summary = json.loads(path.read_text(encoding="utf-8"))
    assert summary["schema"] == DETERMINANT_SUMMARY_SCHEMA
    assert summary["run"]["solved_nodes"] == 1000
    assert summary["run"]["unresolved_nodes"] == 0
    assert summary["run"]["maximum_pending_bound"] > summary["problem"]["target"]
    assert summary["leading_pending_determinants"][
        "ordinary_interval_sign_indefinite"
    ] == 20
    assert summary["leading_pending_determinants"]["vertex_sign_counts"] == {
        "positive": 20
    }
    assert not summary["conclusion"]["target_closed"]


def test_compact_ando_cover_summary_records_strict_nondominant_gain() -> None:
    path = (
        ROOT
        / "scratch"
        / "d2_frontier"
        / "ando_instrument_cover_l055_summary.json"
    )
    summary = json.loads(path.read_text(encoding="utf-8"))
    assert summary["schema"] == ANDO_SUMMARY_SCHEMA
    assert summary["run"]["solved_nodes"] == 1000
    assert summary["run"]["unresolved_nodes"] == 0
    assert summary["run"]["planar_ando_witness_count"] == 8
    assert summary["ando_accounting"]["audits_with_exact_ando_violation"] == 925
    comparison = summary["comparison_to_baseline"]
    assert comparison["strict_bound_improvement"]
    assert not comparison["tree_dominance"]
    assert (
        comparison["ando_maximum_pending_bound"]
        < comparison["baseline_maximum_pending_bound"]
    )
    assert summary["run"]["maximum_pending_bound"] > summary["problem"]["target"]
    assert not summary["conclusion"]["target_closed"]


def known_non_cp_planar_pullbacks() -> tuple[np.ndarray, float, float]:
    """Return positive pulled effects with no common planar CP completion."""

    second = np.asarray(
        [
            [0.24175927315093138, -0.30528594092156486 + 0.05548766814312024j],
            [-0.30528594092156475 - 0.05548766814312024j, 0.39824073004173],
        ]
    )
    third = np.asarray(
        [
            [0.24654148066124223, 0.11892952988987542 - 0.18317018377261665j],
            [0.11892952988987542 + 0.18317018377261671j, 0.19345852058459123],
        ]
    )
    matrices = np.asarray([np.eye(2) - second - third, second, third])
    paulis = np.asarray(
        [
            np.eye(2),
            [[0.0, 1.0], [1.0, 0.0]],
            [[0.0, -1j], [1j, 0.0]],
            [[1.0, 0.0], [0.0, -1.0]],
        ],
        dtype=complex,
    )
    coefficients = np.asarray(
        [
            [
                0.5 * np.trace(matrices[outcome] @ paulis[coordinate]).real
                for outcome in range(3)
            ]
            for coordinate in range(4)
        ]
    )
    denominator = 1.0 / (0.92 + 0.64 - 1.0)
    return coefficients, 0.92 * denominator, 0.64 * denominator


def test_planar_ando_direction_detects_common_cp_failure() -> None:
    pulled, alpha, beta = known_non_cp_planar_pullbacks()
    report = planar_ando_direction(pulled, planar_reconstruction(alpha, beta))
    assert report["violation"] > 0.1
    assert report["exact_margin"] == -report["violation"]
    np.testing.assert_allclose(
        np.linalg.norm(report["test_direction"][1:]),
        1.0,
        atol=2e-12,
    )


def test_zero_width_ando_enclosure_matches_cramers_rule() -> None:
    inputs = np.asarray(
        [
            [0.25, 0.10, 0.00, 0.00],
            [0.25, 0.00, 0.10, 0.00],
            [0.25, 0.00, 0.00, 0.10],
            [0.25, 0.00, 0.00, 0.00],
        ]
    )
    pulled, alpha, beta = known_non_cp_planar_pullbacks()
    report = planar_ando_direction(pulled, planar_reconstruction(alpha, beta))
    sign = 1 if np.linalg.det(inputs) > 0.0 else -1
    _, coefficient_upper, _ = ando_witness_coefficient_bounds(
        inputs,
        inputs,
        (alpha, alpha),
        (beta, beta),
        np.asarray(report["test_direction"]),
        float(report["phase"]),
        sign,
    )
    probabilities = inputs @ pulled
    enclosed_lhs = float(np.sum(coefficient_upper * probabilities))
    exact_lhs = abs(float(np.linalg.det(inputs))) * float(report["exact_margin"])
    np.testing.assert_allclose(
        enclosed_lhs,
        exact_lhs,
        atol=2e-13,
    )


def test_ando_witness_rejects_positive_effects_without_common_cp_map() -> None:
    inputs = np.asarray(
        [
            [0.25, 0.10, 0.00, 0.00],
            [0.25, 0.00, 0.10, 0.00],
            [0.25, 0.00, 0.00, 0.10],
            [0.25, 0.00, 0.00, 0.00],
        ]
    )
    pulled, alpha, beta = known_non_cp_planar_pullbacks()
    statistics = np.zeros((4, 4, 3), dtype=float)
    statistics[:, 0, :] = inputs @ pulled
    povm_witnesses, povm_audit = determinant_povm_witnesses(
        inputs - 1e-8,
        inputs + 1e-8,
        inputs,
        statistics,
        2e-12,
    )
    ando_witnesses, ando_audit = determinant_ando_witnesses(
        inputs - 1e-8,
        inputs + 1e-8,
        (alpha - 1e-8, alpha + 1e-8),
        (beta - 1e-8, beta + 1e-8),
        inputs,
        statistics,
        2e-12,
    )
    assert povm_audit["negative_effect_count"] == 0
    assert not povm_witnesses
    assert ando_audit["violated_direction_count"] == 1
    assert ando_audit["robust_witness_count"] == 1
    assert ando_witnesses[0]["kind"] == "planar-ando"
    assert ando_witnesses[0]["violation"] > 1e-6


def test_ando_enclosure_does_not_reject_a_physical_cp_map() -> None:
    rng = np.random.default_rng(20260825)
    inputs = np.asarray(
        [
            [0.25, 0.10, 0.00, 0.00],
            [0.25, 0.00, 0.10, 0.00],
            [0.25, 0.00, 0.00, 0.10],
            [0.25, 0.00, 0.00, 0.00],
        ]
    )
    alpha, beta = 1.61, 1.13
    terminal = np.asarray(
        [pauli_effect_matrix(row) for row in planar_effect_pauli(alpha, beta)]
    )
    kraus = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
    kraus /= 1.4 * np.linalg.svd(kraus, compute_uv=False)[0]
    pulled_matrices = np.asarray(
        [kraus.conj().T @ effect @ kraus for effect in terminal]
    )
    paulis = np.asarray(
        [
            np.eye(2),
            [[0.0, 1.0], [1.0, 0.0]],
            [[0.0, -1j], [1j, 0.0]],
            [[1.0, 0.0], [0.0, -1.0]],
        ],
        dtype=complex,
    )
    pulled = np.asarray(
        [
            [
                0.5
                * np.trace(pulled_matrices[outcome] @ paulis[coordinate]).real
                for outcome in range(3)
            ]
            for coordinate in range(4)
        ]
    )
    probabilities = inputs @ pulled
    statistics = np.zeros((4, 4, 3), dtype=float)
    statistics[:, 0, :] = probabilities
    witnesses, audit = determinant_ando_witnesses(
        inputs - 2e-6,
        inputs + 2e-6,
        (alpha - 2e-6, alpha + 2e-6),
        (beta - 2e-6, beta + 2e-6),
        inputs,
        statistics,
        2e-12,
    )
    assert not witnesses
    assert audit["robust_witness_count"] == 0
    sign = 1 if np.linalg.det(inputs) > 0.0 else -1
    for _ in range(20):
        bloch = rng.normal(size=3)
        bloch /= np.linalg.norm(bloch)
        test = np.concatenate(([1.0], bloch))
        phase = rng.uniform(0.0, 2.0 * math.pi)
        _, coefficient_upper, _ = ando_witness_coefficient_bounds(
            inputs - 2e-6,
            inputs + 2e-6,
            (alpha - 2e-6, alpha + 2e-6),
            (beta - 2e-6, beta + 2e-6),
            test,
            phase,
            sign,
        )
        assert float(np.sum(coefficient_upper * probabilities)) >= -2e-13


def test_ando_margin_branching_targets_enclosure_gap() -> None:
    inputs = np.asarray(
        [
            [0.25, 0.10, 0.00, 0.00],
            [0.25, 0.00, 0.10, 0.00],
            [0.25, 0.00, 0.00, 0.10],
            [0.25, 0.00, 0.00, 0.00],
        ]
    )
    pulled, alpha, beta = known_non_cp_planar_pullbacks()
    statistics = np.zeros((4, 4, 3), dtype=float)
    statistics[:, 0, :] = inputs @ pulled
    lower = inputs - 0.003
    upper = inputs + 0.003
    witnesses, audit = determinant_ando_witnesses(
        lower,
        upper,
        (alpha - 0.001, alpha + 0.001),
        (beta - 0.001, beta + 0.001),
        inputs,
        statistics,
        2e-12,
    )
    assert not witnesses
    assert audit["violated_direction_count"] == 1
    scores = ando_input_split_scores(
        lower,
        upper,
        (alpha - 0.001, alpha + 0.001),
        (beta - 0.001, beta + 0.001),
        statistics,
        audit,
    )
    assert np.max(scores) > 0.01

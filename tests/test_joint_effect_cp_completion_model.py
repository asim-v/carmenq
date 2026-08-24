from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pytest


RESEARCH = Path(__file__).resolve().parents[1] / "scratch" / "d2_frontier"
sys.path.insert(0, str(RESEARCH))

from joint_effect_helstrom_scip import (  # noqa: E402
    IDENTITY,
    build,
    canonical_three_effect_povm,
    extract_common_instrument,
    normalized_left_null_chart,
    seed_from_common_instrument,
)
from basis_null_chart_cover import Cell, bisect_cell  # noqa: E402


def test_exact_completion_accepts_and_reconstructs_literal_instrument(
    tmp_path: Path,
) -> None:
    effects = canonical_three_effect_povm(np.asarray([0.92, 0.64, 0.44]))
    states = np.repeat((IDENTITY / 8.0)[None, :, :], 4, axis=0)
    # Four equiprobable completely depolarizing subchannels. Their Choi
    # matrices sum to the Choi matrix of a trace-preserving channel.
    choi = np.repeat((np.eye(4) / 8.0)[None, :, :], 4, axis=0)
    checkpoint = tmp_path / "literal_instrument.npz"
    np.savez(checkpoint, states=states, choi=choi, effects=effects)

    model, variables = build(
        effects,
        weight=0.55,
        prefix_order=(0, 1, 2, 3),
        target=None,
        fix_rotation_gauge=True,
        linked_columns=None,
        require_cp_completion=True,
        ando_direction_count=16,
    )
    model.hideOutput()
    assert seed_from_common_instrument(
        model, variables, checkpoint, effects, weight=0.55
    )
    extracted = extract_common_instrument(model, variables, effects, weight=0.55)
    assert extracted is not None
    report, arrays = extracted
    assert report["score_from_reconstructed_instrument"] == pytest.approx(0.5875)
    assert report["audit_from_reconstructed_instrument"] == pytest.approx(0.25)
    assert report["return_from_reconstructed_instrument"] == pytest.approx(1.0)
    assert report["minimum_choi_eigenvalue"] >= -1e-12
    assert report["trace_preservation_residual"] <= 1e-12
    assert np.linalg.norm(arrays["choi"] - choi) <= 2e-12


def test_normalized_null_charts_cover_every_singular_operator_basis() -> None:
    matrix = np.asarray(
        [
            [0.35, 0.20, 0.00, 0.00],
            [0.27, -0.10, 0.16, 0.00],
            [0.21, 0.03, -0.08, 0.12],
            [0.35 + 0.4 * 0.27 - 0.2 * 0.21,
             0.20 + 0.4 * -0.10 - 0.2 * 0.03,
             0.4 * 0.16 - 0.2 * -0.08,
             -0.2 * 0.12],
        ]
    )
    pivot, coefficients = normalized_left_null_chart(matrix)
    assert coefficients[pivot] == pytest.approx(1.0)
    assert np.max(np.abs(coefficients)) <= 1.0 + 1e-12
    assert coefficients @ matrix == pytest.approx(np.zeros(4), abs=2e-12)


def test_normalized_null_chart_rejects_a_nonsingular_basis() -> None:
    with pytest.raises(np.linalg.LinAlgError, match="nonsingular"):
        normalized_left_null_chart(np.eye(4))


def test_null_chart_cell_bisection_is_an_exact_partition() -> None:
    parent = Cell(
        identifier=7,
        pivot=2,
        depth=3,
        bounds=np.asarray([[-0.25, 0.0], [-0.5, -0.25], [-0.5, -0.25]]),
        inherited_upper=0.758053,
    )
    left, right = bisect_cell(parent, 8)
    assert left.parent == right.parent == parent.identifier
    assert left.volume + right.volume == pytest.approx(parent.volume)
    assert left.bounds[0, 1] == pytest.approx(right.bounds[0, 0])
    assert left.bounds[1:] == pytest.approx(parent.bounds[1:])
    assert right.bounds[1:] == pytest.approx(parent.bounds[1:])


def test_inverse_basis_chart_requires_a_compact_nondegenerate_branch() -> None:
    effects = canonical_three_effect_povm(np.asarray([0.92, 0.64, 0.44]))
    with pytest.raises(ValueError, match="nonzero determinant branch"):
        build(
            effects,
            weight=0.55,
            prefix_order=(0, 1, 2, 3),
            target=None,
            fix_rotation_gauge=True,
            linked_columns=None,
            require_cp_completion=True,
            basis_inverse_bound=10.0,
        )

    model, variables = build(
        effects,
        weight=0.55,
        prefix_order=(0, 1, 2, 3),
        target=None,
        fix_rotation_gauge=True,
        linked_columns=None,
        require_cp_completion=True,
        basis_determinant_sign=1,
        basis_determinant_floor=1e-3,
        basis_inverse_bound=100.0,
    )
    assert variables["basis_inverse_bound"] == 100.0
    assert len(variables["basis_inverse"]) == 4
    assert len(variables["basis_transfer"]) == 64
    model.freeProb()


def test_lifted_adjugate_exposes_bounded_intermediate_variables() -> None:
    effects = canonical_three_effect_povm(np.asarray([0.92, 0.64, 0.44]))
    model, variables = build(
        effects,
        weight=0.55,
        prefix_order=(0, 1, 2, 3),
        target=None,
        fix_rotation_gauge=True,
        linked_columns=None,
        require_cp_completion=True,
        prefix_prior_bounds=np.asarray(
            [
                [0.296875, 0.42596435546875],
                [0.224609375, 0.34832000732421875],
                [0.15234375, 0.258392333984375],
                [0.1083984375, 0.201324462890625],
            ]
        ),
        basis_determinant_sign=1,
        basis_lifted_adjugate=True,
    )
    assert len(variables["basis_lifted_adjugate"]) == 4
    assert len(variables["basis_lifted_numerators"]) == 64
    assert variables["basis_adjugate_bound"] < 0.109
    model.freeProb()


def test_nondegenerate_basis_completion_does_not_need_a_second_channel_factor() -> None:
    effects = canonical_three_effect_povm(np.asarray([0.92, 0.64, 0.44]))
    model, variables = build(
        effects,
        weight=0.55,
        prefix_order=(0, 1, 2, 3),
        target=None,
        fix_rotation_gauge=True,
        linked_columns=None,
        require_cp_completion=False,
        basis_determinant_sign=1,
        basis_determinant_floor=1e-5,
        basis_lifted_adjugate=True,
    )
    assert not variables["cp_missing_pullback"]
    assert sum(name.startswith("basis_output_z_") for name in variables) == 16
    assert len(variables["basis_cholesky_factors"]) == 4
    model.freeProb()


@pytest.mark.parametrize("branch", ["scalar-positive", "bloch"])
def test_exact_flagged_pair_contraction_has_four_output_norms(branch: str) -> None:
    effects = canonical_three_effect_povm(np.asarray([0.92, 0.64, 0.44]))
    model, variables = build(
        effects,
        weight=0.55,
        prefix_order=(0, 1, 2, 3),
        target=None,
        fix_rotation_gauge=True,
        linked_columns=None,
        require_cp_completion=True,
        flagged_contraction_coefficients=np.asarray([0.0, 0.0, 1.0, -1.0]),
        flagged_contraction_branch=branch,
    )
    assert variables["flagged_contraction_branch"] == branch
    assert len(variables["flagged_contraction_output_norms"]) == 4
    model.freeProb()


def test_flagged_bloch_cap_replaces_the_nonconvex_input_norm_variable() -> None:
    effects = canonical_three_effect_povm(np.asarray([0.92, 0.64, 0.44]))
    model, variables = build(
        effects,
        weight=0.55,
        prefix_order=(0, 1, 2, 3),
        target=None,
        fix_rotation_gauge=True,
        linked_columns=None,
        require_cp_completion=True,
        flagged_contraction_coefficients=np.asarray([0.0, 0.0, 1.0, -1.0]),
        flagged_contraction_branch="bloch",
        flagged_bloch_cap=np.asarray([0.0, 1.0, 0.0, 0.98]),
    )
    assert variables["flagged_bloch_cap"] == pytest.approx([0.0, 1.0, 0.0, 0.98])
    assert "flagged_contraction_input_norm" not in variables
    model.freeProb()


def test_flagged_l1_orthant_uses_a_finite_linear_input_upper_bound() -> None:
    effects = canonical_three_effect_povm(np.asarray([0.92, 0.64, 0.44]))
    model, variables = build(
        effects,
        weight=0.55,
        prefix_order=(0, 1, 2, 3),
        target=None,
        fix_rotation_gauge=True,
        linked_columns=None,
        require_cp_completion=False,
        flagged_contraction_coefficients=np.asarray([0.0, 0.0, 1.0, -1.0]),
        flagged_contraction_branch="l1-upper",
        flagged_l1_signs=np.asarray([1.0, 1.0, 1.0, -1.0]),
    )
    assert variables["flagged_l1_signs"] == pytest.approx([1.0, 1.0, 1.0, -1.0])
    assert len(variables["flagged_contraction_output_norms"]) == 4
    model.freeProb()

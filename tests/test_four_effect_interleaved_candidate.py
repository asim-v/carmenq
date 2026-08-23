"""Independent tensor checks for the four-effect interleaved MPS branch."""

from math import isclose

from scripts.verify_four_effect_interleaved_candidate import evaluate


def test_four_effect_candidate_is_a_legal_bond_two_leaf() -> None:
    result = evaluate(0.6)
    assert result["temporal_ranks"] == [2, 2, 2]
    assert result["maximum_row_isometry_residual"] < 1e-13
    assert result["maximum_pauli_completion_residual"] < 1e-13
    assert result["gram_off_diagonal_frobenius"] < 1e-13
    assert isclose(result["normalisation"], 1.0, abs_tol=1e-13)
    assert isclose(result["direct_audit"], 0.8699300211, abs_tol=2e-8)
    assert isclose(result["direct_return"], 0.6098520066, abs_tol=2e-8)
    assert isclose(result["direct_score"], 0.765898815264694, abs_tol=3e-12)
    assert abs(result["audit_residual"]) < 1e-13
    assert abs(result["return_residual"]) < 1e-13
    assert abs(result["score_residual"]) < 1e-13

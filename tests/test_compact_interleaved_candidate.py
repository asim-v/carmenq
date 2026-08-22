"""Independent tensor checks for the compact interleaved MPS construction."""

from math import isclose

from scripts.verify_compact_interleaved_candidate import evaluate


def test_balanced_compact_candidate_is_a_legal_bond_two_leaf() -> None:
    result = evaluate(0.5)
    assert result["temporal_ranks"] == [2, 2, 2]
    assert result["maximum_row_isometry_residual"] < 1e-13
    assert result["maximum_pauli_completion_residual"] < 1e-13
    assert result["gram_off_diagonal_frobenius"] < 1e-13
    assert isclose(result["normalisation"], 1.0, abs_tol=1e-13)
    # The support is sharply reproducible, while its two coordinates move in
    # compensating directions along a numerically flat tangent of the optimum.
    assert isclose(result["direct_audit"], 0.620085075586, abs_tol=5e-9)
    assert isclose(result["direct_return"], 0.899520492117, abs_tol=5e-9)
    assert isclose(result["direct_score"], 0.759802783851444, abs_tol=3e-12)
    assert abs(result["audit_residual"]) < 1e-13
    assert abs(result["return_residual"]) < 1e-13
    assert abs(result["score_residual"]) < 1e-13

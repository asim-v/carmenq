"""Regression test for the stored finite-outcome counterexample."""

from pathlib import Path

from scripts.verify_interleaved_counterexample import evaluate


ROOT = Path(__file__).resolve().parents[1]


def test_ternary_instrument_strictly_beats_two_parameter_candidate() -> None:
    result = evaluate(ROOT / "data" / "interleaved_ternary_counterexample.npz")
    assert result["strict_excess"] > 0.004
    assert result["maximum_local_completeness_residual"] < 1e-12
    assert result["maximum_povm_completeness_residual"] < 1e-12
    assert result["minimum_povm_eigenvalue"] > -1e-12
    assert result["active_transcript_count"] == 16
    assert result["linear_tail_slack"] > 0.29
    assert result["margin_below_linear_tail_upper_bound"] > 0.08
    assert abs(
        result["rank_four_spectral_tail"]
        - (1.0 - result["best_four_word_list_mass"])
    ) < 1e-12
    assert result["causal_four_word_list_mass"] > 0.5477

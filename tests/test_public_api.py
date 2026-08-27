"""Tests for the stable CARMEN-Q entry points."""

from carmenq import (
    GROUPED_CHECK_MATRIX,
    INTERLEAVED_BALANCED_COUNTEREXAMPLE,
    INTERLEAVED_ORDER_GAP_WEIGHT_THRESHOLD,
    INTERLEAVED_L060_CERTIFIED_INTERVAL,
    INTERLEAVED_PERFECT_AUDIT_ENDPOINT,
    __version__,
    certify,
    certify_classical_memory,
    collective_bound,
    collective_classical_record_bound,
    full_crossing_perfect_audit_return_bound,
    full_rank_block_approximate_audit_return_bound,
    full_rank_block_packing_number,
    full_rank_block_perfect_audit_return_bound,
    grouped_frontier,
    interleaved_best_known_lower_bound,
    interleaved_candidate_lower_bound,
    interleaved_compact_lower_bound,
    interleaved_four_effect_lower_bound,
    interleaved_return_upper_bound,
    interleaved_support_upper_bound,
    ordered_check_perfect_audit_return_bound,
    plan,
    plan_experiment,
    reconstruct_common_instrument_from_basis,
    streaming_bound,
)


def test_concise_aliases_match_scientific_api() -> None:
    assert streaming_bound(8, 0.5) == 0.75
    assert collective_bound(0.5) == collective_classical_record_bound(0.5)
    assert certify is certify_classical_memory
    assert plan is plan_experiment
    assert __version__ == "2.3.1"
    assert callable(reconstruct_common_instrument_from_basis)


def test_order_sensitive_exact_results_are_public() -> None:
    assert grouped_frontier(1.0).return_fidelity == 0.5
    assert GROUPED_CHECK_MATRIX == ((1, 1, 0, 0), (0, 0, 1, 1))
    assert INTERLEAVED_PERFECT_AUDIT_ENDPOINT.maximum_return_fidelity == 0.25
    assert full_crossing_perfect_audit_return_bound(2, 2) == 0.25
    assert full_rank_block_perfect_audit_return_bound(2, 2, 3) == 0.125
    assert full_rank_block_packing_number(GROUPED_CHECK_MATRIX) == 1
    assert ordered_check_perfect_audit_return_bound(GROUPED_CHECK_MATRIX, 2) == 0.5
    assert INTERLEAVED_ORDER_GAP_WEIGHT_THRESHOLD == 3 / 7
    interval = INTERLEAVED_L060_CERTIFIED_INTERVAL
    assert interval.audit_weight == 0.6
    assert interval.lower_bound == 0.7658988152
    assert interval.upper_bound == 0.76670
    assert interval.exact_optimum_known is False
    assert interval.lower_fraction == (957_373_519, 1_250_000_000)
    assert interval.upper_fraction == (7_667, 10_000)
    assert interval.width_fraction.numerator == 1_001_481
    assert full_rank_block_approximate_audit_return_bound(1.0, 2, 2, 2) == 0.25
    assert interleaved_return_upper_bound(1.0) == 0.25
    assert interleaved_support_upper_bound(0.5) < grouped_frontier(0.5).support_value
    candidate = interleaved_candidate_lower_bound(0.5)
    assert candidate.support_value > 0.7554
    assert candidate.support_is_globally_optimal is False
    counterexample = INTERLEAVED_BALANCED_COUNTEREXAMPLE
    assert counterexample.support_value > candidate.support_value + 0.004
    assert counterexample.independently_verified is True
    assert counterexample.support_is_globally_optimal is False
    compact = interleaved_compact_lower_bound(0.5)
    assert compact.support_value > counterexample.support_value + 0.0003
    assert compact.support_is_globally_optimal is False
    four_effect = interleaved_four_effect_lower_bound(0.6)
    assert four_effect.support_value > interleaved_compact_lower_bound(0.6).support_value
    assert four_effect.support_is_globally_optimal is False
    assert interleaved_best_known_lower_bound(0.6).strategy == "four_effect_mps"

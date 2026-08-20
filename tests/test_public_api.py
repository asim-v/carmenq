"""Tests for the stable CARMEN-Q entry points."""

from carmenq import (
    GROUPED_CHECK_MATRIX,
    INTERLEAVED_BALANCED_COUNTEREXAMPLE,
    INTERLEAVED_PERFECT_AUDIT_ENDPOINT,
    __version__,
    certify,
    certify_classical_memory,
    collective_bound,
    collective_classical_record_bound,
    full_crossing_perfect_audit_return_bound,
    full_rank_block_packing_number,
    full_rank_block_perfect_audit_return_bound,
    grouped_frontier,
    interleaved_candidate_lower_bound,
    ordered_check_perfect_audit_return_bound,
    plan,
    plan_experiment,
    streaming_bound,
)


def test_concise_aliases_match_scientific_api() -> None:
    assert streaming_bound(8, 0.5) == 0.75
    assert collective_bound(0.5) == collective_classical_record_bound(0.5)
    assert certify is certify_classical_memory
    assert plan is plan_experiment
    assert __version__ == "2.1.0"


def test_order_sensitive_exact_results_are_public() -> None:
    assert grouped_frontier(1.0).return_fidelity == 0.5
    assert GROUPED_CHECK_MATRIX == ((1, 1, 0, 0), (0, 0, 1, 1))
    assert INTERLEAVED_PERFECT_AUDIT_ENDPOINT.maximum_return_fidelity == 0.25
    assert full_crossing_perfect_audit_return_bound(2, 2) == 0.25
    assert full_rank_block_perfect_audit_return_bound(2, 2, 3) == 0.125
    assert full_rank_block_packing_number(GROUPED_CHECK_MATRIX) == 1
    assert ordered_check_perfect_audit_return_bound(GROUPED_CHECK_MATRIX, 2) == 0.5
    candidate = interleaved_candidate_lower_bound(0.5)
    assert candidate.support_value > 0.7554
    assert candidate.support_is_globally_optimal is False
    counterexample = INTERLEAVED_BALANCED_COUNTEREXAMPLE
    assert counterexample.support_value > candidate.support_value + 0.004
    assert counterexample.independently_verified is True
    assert counterexample.support_is_globally_optimal is False

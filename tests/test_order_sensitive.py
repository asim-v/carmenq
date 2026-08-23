"""Tests for exact order-sensitive AUDIT--RETURN utilities."""

import itertools
from math import isclose, sqrt

import numpy as np
import pytest

from carmenq.order_sensitive import (
    GROUPED_CHECK_MATRIX,
    INTERLEAVED_CHECK_MATRIX,
    INTERLEAVED_ORDER_GAP_WEIGHT_THRESHOLD,
    INTERLEAVED_PERFECT_AUDIT_ENDPOINT,
    full_crossing_cuts,
    full_crossing_perfect_audit_return_bound,
    full_rank_block_packing_number,
    full_rank_block_approximate_audit_return_bound,
    full_rank_block_perfect_audit_return_bound,
    gf2_rank,
    grouped_frontier,
    interleaved_best_known_lower_bound,
    interleaved_candidate_lower_bound,
    interleaved_candidate_scores,
    interleaved_compact_lower_bound,
    interleaved_four_effect_lower_bound,
    interleaved_return_upper_bound,
    interleaved_support_upper_bound,
    ordered_check_perfect_audit_return_bound,
    rank_two_static_qubit_support,
    trellis_connectivity_profile,
    trellis_connectivity_tau,
)


def test_gf2_rank_reduces_integer_entries_modulo_two() -> None:
    assert gf2_rank([[1, 0], [0, 1]]) == 2
    assert gf2_rank([[1, 1, 0], [1, 1, 0]]) == 1
    assert gf2_rank([[2, 3], [4, 5]]) == 1
    assert gf2_rank(np.empty((0, 4), dtype=int)) == 0
    assert gf2_rank(np.empty((3, 0), dtype=bool)) == 0


@pytest.mark.parametrize("matrix", [[], [0, 1, 0], np.zeros((2, 2, 1), dtype=int)])
def test_gf2_rank_requires_a_matrix(matrix: object) -> None:
    with pytest.raises(ValueError, match="two-dimensional"):
        gf2_rank(matrix)


def test_gf2_rank_rejects_noninteger_entries() -> None:
    with pytest.raises(TypeError, match="integers or booleans"):
        gf2_rank([[0.0, 1.0]])


def test_order_permutation_changes_cut_profile_but_not_rank() -> None:
    assert gf2_rank(GROUPED_CHECK_MATRIX) == 2
    assert gf2_rank(INTERLEAVED_CHECK_MATRIX) == 2
    assert trellis_connectivity_profile(GROUPED_CHECK_MATRIX) == (0, 1, 0, 1, 0)
    assert trellis_connectivity_profile(INTERLEAVED_CHECK_MATRIX) == (
        0,
        1,
        2,
        1,
        0,
    )
    assert trellis_connectivity_tau(GROUPED_CHECK_MATRIX) == 1
    assert trellis_connectivity_tau(INTERLEAVED_CHECK_MATRIX) == 2


def test_full_crossing_cut_detection() -> None:
    assert full_crossing_cuts(GROUPED_CHECK_MATRIX) == ()
    assert full_crossing_cuts(INTERLEAVED_CHECK_MATRIX) == (2,)
    assert full_crossing_cuts([[0, 0, 0]]) == ()
    assert full_crossing_cuts([[1, 1]]) == (1,)


def test_full_rank_block_packing_captures_temporal_repetition() -> None:
    assert full_rank_block_packing_number(GROUPED_CHECK_MATRIX) == 1
    assert full_rank_block_packing_number(INTERLEAVED_CHECK_MATRIX) == 2
    assert full_rank_block_packing_number([[0, 0, 0]]) == 0
    assert full_rank_block_packing_number([[1, 0, 1, 0, 0]]) == 2
    for rank in range(1, 5):
        identity = np.eye(rank, dtype=int)
        for block_count in range(1, 5):
            repeated = np.tile(identity, (1, block_count))
            assert full_rank_block_packing_number(repeated) == block_count


def test_ordered_matrix_bound_unifies_grouped_and_interleaved_endpoints() -> None:
    assert ordered_check_perfect_audit_return_bound(GROUPED_CHECK_MATRIX, 2) == 0.5
    assert (
        ordered_check_perfect_audit_return_bound(INTERLEAVED_CHECK_MATRIX, 2)
        == 0.25
    )
    repeated = np.tile(np.eye(2, dtype=int), (1, 3))
    assert ordered_check_perfect_audit_return_bound(repeated, 2) == 0.125
    assert ordered_check_perfect_audit_return_bound([[0, 0]], 1) == 1.0
    with pytest.raises(ValueError, match="coherent_dimension must be positive"):
        ordered_check_perfect_audit_return_bound([[0, 0]], 0)


def test_same_columns_can_have_exponentially_different_order_bounds() -> None:
    for rank in range(2, 5):
        identity = np.eye(rank, dtype=int)
        for repetitions in range(2, 5):
            batched = np.repeat(identity, repetitions, axis=1)
            cycled = np.tile(identity, (1, repetitions))
            assert sorted(map(tuple, batched.T)) == sorted(map(tuple, cycled.T))
            assert full_rank_block_packing_number(batched) == 1
            assert full_rank_block_packing_number(cycled) == repetitions
            for retained_coordinates in range(1, rank):
                dimension = 2**retained_coordinates
                base = dimension / (2**rank)
                assert ordered_check_perfect_audit_return_bound(
                    batched, dimension
                ) == base
                assert ordered_check_perfect_audit_return_bound(
                    cycled, dimension
                ) == base**repetitions


@pytest.mark.parametrize(
    ("rank", "dimension", "alphabet", "expected"),
    [
        (1, 1, 2, 1 / 4),
        (1, 2, 2, 1.0),
        (2, 1, 2, 1 / 16),
        (2, 2, 2, 1 / 4),
        (2, 3, 2, 9 / 16),
        (3, 2, 2, 1 / 16),
        (2, 2, 3, 4 / 81),
    ],
)
def test_full_crossing_dimension_bound(
    rank: int, dimension: int, alphabet: int, expected: float
) -> None:
    assert full_crossing_perfect_audit_return_bound(
        rank, dimension, alphabet
    ) == expected


def test_full_crossing_square_law_matches_canonical_construction() -> None:
    for alphabet in (2, 3, 5):
        for rank in range(1, 5):
            for retained_coordinates in range(rank + 1):
                dimension = alphabet**retained_coordinates
                expected = alphabet ** (-2 * (rank - retained_coordinates))
                assert isclose(
                    full_crossing_perfect_audit_return_bound(
                        rank, dimension, alphabet
                    ),
                    expected,
                    abs_tol=1e-16,
                )


def test_temporal_power_law_matches_repeated_identity_construction() -> None:
    for alphabet in (2, 3, 4, 5):
        for rank in range(1, 5):
            for block_count in range(1, 5):
                for retained_coordinates in range(rank + 1):
                    dimension = alphabet**retained_coordinates
                    expected = alphabet ** (
                        -block_count * (rank - retained_coordinates)
                    )
                    assert isclose(
                        full_rank_block_perfect_audit_return_bound(
                            rank,
                            dimension,
                            block_count,
                            alphabet,
                        ),
                        expected,
                        abs_tol=1e-16,
                    )


def test_approximate_audit_bound_reduces_to_temporal_power_endpoint() -> None:
    for alphabet in (2, 3, 4, 5):
        for rank in range(1, 4):
            for block_count in range(1, 4):
                for retained_coordinates in range(rank + 1):
                    dimension = alphabet**retained_coordinates
                    assert isclose(
                        full_rank_block_approximate_audit_return_bound(
                            1.0,
                            rank,
                            dimension,
                            block_count,
                            alphabet,
                        ),
                        full_rank_block_perfect_audit_return_bound(
                            rank,
                            dimension,
                            block_count,
                            alphabet,
                        ),
                        abs_tol=1e-15,
                    )


def test_interleaved_linear_tail_return_bound() -> None:
    assert isclose(interleaved_return_upper_bound(1.0), 0.25, abs_tol=1e-15)
    assert isclose(
        interleaved_return_upper_bound(0.99),
        0.25 + 0.01 + sqrt(1.5 * 0.01 * 0.98),
        abs_tol=2e-15,
    )
    assert interleaved_return_upper_bound(0.625) == 1.0
    assert interleaved_return_upper_bound(0.5) == 1.0


@pytest.mark.parametrize("probability", [-0.01, 1.01, float("nan")])
def test_approximate_audit_bound_validates_probability(probability: float) -> None:
    with pytest.raises(ValueError, match=r"\[0, 1\]"):
        full_rank_block_approximate_audit_return_bound(
            probability, 2, 2, 2
        )


def test_interleaved_support_certificate_is_strict_above_three_sevenths() -> None:
    threshold = INTERLEAVED_ORDER_GAP_WEIGHT_THRESHOLD
    assert threshold == 3.0 / 7.0
    assert isclose(
        interleaved_support_upper_bound(threshold),
        rank_two_static_qubit_support(threshold),
        abs_tol=2e-15,
    )
    assert isclose(
        interleaved_support_upper_bound(0.5),
        5.0 / 8.0 + sqrt(3.0) / 8.0,
        abs_tol=2e-15,
    )
    for weight in (0.43, 0.5, 0.75, 0.99):
        assert (
            interleaved_support_upper_bound(weight)
            < rank_two_static_qubit_support(weight)
        )
    assert interleaved_support_upper_bound(1.0) == 1.0


def test_full_crossing_wrapper_is_the_two_block_power_law() -> None:
    for alphabet in (2, 3, 4, 5):
        for rank in range(1, 5):
            for dimension in (1, alphabet, alphabet**rank, alphabet ** (rank + 1)):
                assert full_crossing_perfect_audit_return_bound(
                    rank, dimension, alphabet
                ) == full_rank_block_perfect_audit_return_bound(
                    rank, dimension, 2, alphabet
                )


@pytest.mark.parametrize(
    ("arguments", "error", "message"),
    [
        ((0, 1, 2), ValueError, "syndrome_rank must be positive"),
        ((1, 0, 2), ValueError, "coherent_dimension must be positive"),
        ((1, 1, 1), ValueError, "alphabet_size must be at least two"),
        ((1, 1, 6), ValueError, "alphabet_size must be a prime power"),
        ((1.5, 1, 2), TypeError, "syndrome_rank must be an integer"),
    ],
)
def test_full_crossing_dimension_bound_validates_inputs(
    arguments: tuple[object, object, object],
    error: type[Exception],
    message: str,
) -> None:
    with pytest.raises(error, match=message):
        full_crossing_perfect_audit_return_bound(*arguments)  # type: ignore[arg-type]


@pytest.mark.parametrize(
    ("block_count", "error", "message"),
    [
        (0, ValueError, "block_count must be positive"),
        (1.5, TypeError, "block_count must be an integer"),
    ],
)
def test_temporal_power_law_validates_block_count(
    block_count: object,
    error: type[Exception],
    message: str,
) -> None:
    with pytest.raises(error, match=message):
        full_rank_block_perfect_audit_return_bound(
            2, 2, block_count  # type: ignore[arg-type]
        )


@pytest.mark.parametrize("weight", [0.0, 0.2, 0.5, 0.8, 1.0])
def test_grouped_frontier_attains_static_support(weight: float) -> None:
    point = grouped_frontier(weight)
    expected = rank_two_static_qubit_support(weight)
    direct_score = (
        weight * point.audit_probability
        + (1.0 - weight) * point.return_fidelity
    )
    assert isclose(point.support_value, expected, abs_tol=1e-15)
    assert isclose(point.support_value, direct_score, abs_tol=1e-15)
    assert isclose(
        (2.0 * point.audit_probability - 1.0) ** 2
        + (2.0 * point.return_fidelity - 1.0) ** 2,
        1.0,
        abs_tol=2e-15,
    )


def test_balanced_grouped_frontier() -> None:
    point = grouped_frontier(0.5)
    expected = (1.0 + 1.0 / sqrt(2.0)) / 2.0
    assert isclose(point.weak_measurement_strength, 1.0 / sqrt(2.0))
    assert isclose(point.audit_probability, expected)
    assert isclose(point.return_fidelity, expected)
    assert isclose(point.support_value, expected)


def test_grouped_frontier_endpoints_are_explicit() -> None:
    return_only = grouped_frontier(0.0)
    assert return_only.audit_probability == 0.5
    assert return_only.return_fidelity == 1.0
    assert return_only.support_value == 1.0

    audit_only = grouped_frontier(1.0)
    assert audit_only.audit_probability == 1.0
    assert audit_only.return_fidelity == 0.5
    assert audit_only.support_value == 1.0


@pytest.mark.parametrize("weight", [-0.01, 1.01, float("nan")])
def test_support_functions_validate_audit_weight(weight: float) -> None:
    with pytest.raises(ValueError, match=r"\[0, 1\]"):
        rank_two_static_qubit_support(weight)
    with pytest.raises(ValueError, match=r"\[0, 1\]"):
        grouped_frontier(weight)
    with pytest.raises(ValueError, match=r"\[0, 1\]"):
        interleaved_support_upper_bound(weight)


def test_exact_interleaved_endpoint_metadata() -> None:
    endpoint = INTERLEAVED_PERFECT_AUDIT_ENDPOINT
    assert endpoint.check_matrix == INTERLEAVED_CHECK_MATRIX
    assert endpoint.coherent_dimension == 2
    assert endpoint.audit_probability == 1.0
    assert endpoint.maximum_return_fidelity == 0.25
    assert endpoint.bound_is_attained is True
    assert endpoint.full_crossing_cuts == full_crossing_cuts(endpoint.check_matrix)
    assert endpoint.interior_frontier_known is False


def test_endpoint_exhibits_exact_order_gap_at_perfect_audit() -> None:
    grouped = grouped_frontier(1.0)
    interleaved = INTERLEAVED_PERFECT_AUDIT_ENDPOINT
    assert grouped.audit_probability == interleaved.audit_probability == 1.0
    assert grouped.return_fidelity == 0.5
    assert interleaved.maximum_return_fidelity == 0.25


def test_interleaved_candidate_closed_scores() -> None:
    audit, returned = interleaved_candidate_scores(
        0.6168956030718684, 0.8003177036431812
    )
    assert isclose(audit, 0.6446434022644623, abs_tol=2e-15)
    assert isclose(returned, 0.8662314901930318, abs_tol=2e-15)
    endpoint = interleaved_candidate_scores(1.0, 1.0 / sqrt(2.0))
    assert isclose(endpoint[0], 1.0, abs_tol=2e-15)
    assert isclose(endpoint[1], 0.25, abs_tol=2e-15)


@pytest.mark.parametrize(
    ("q", "v"), [(-0.1, 0.5), (1.1, 0.5), (0.5, -0.1), (0.5, 1.1)]
)
def test_interleaved_candidate_validates_parameters(q: float, v: float) -> None:
    with pytest.raises(ValueError, match=r"\[0, 1\]"):
        interleaved_candidate_scores(q, v)


def test_interleaved_candidate_balanced_lower_bound() -> None:
    point = interleaved_candidate_lower_bound(0.5)
    assert point.strategy == "two_parameter"
    assert point.q is not None and isclose(point.q, 0.6168956031, abs_tol=2e-7)
    assert point.v is not None and isclose(point.v, 0.8003177036, abs_tol=2e-7)
    assert isclose(point.support_value, 0.755437446228747, abs_tol=2e-13)
    assert point.support_is_globally_optimal is False


def test_interleaved_candidate_selects_no_record_below_transition() -> None:
    point = interleaved_candidate_lower_bound(0.47)
    assert point.strategy == "no_record"
    assert point.q is None and point.v is None
    assert point.audit_probability == 0.5
    assert point.return_fidelity == 1.0
    assert point.support_value == 0.765


def test_compact_interleaved_candidate_reproduces_balanced_mps_point() -> None:
    point = interleaved_compact_lower_bound(0.5)
    assert point.strategy == "three_effect_mps"
    assert point.t is not None and isclose(point.t, 0.45807398, abs_tol=3e-7)
    assert point.r is not None and isclose(point.r, 0.01637352, abs_tol=3e-7)
    assert point.priors is not None
    assert np.allclose(
        point.priors,
        (0.19921398, 0.09721290, 0.35178656, 0.35178656),
        atol=3e-7,
        rtol=0.0,
    )
    assert isclose(point.audit_probability, 0.620085075586, abs_tol=5e-9)
    assert isclose(point.return_fidelity, 0.899520492117, abs_tol=5e-9)
    assert isclose(point.support_value, 0.759802783851444, abs_tol=2e-11)
    assert point.support_is_globally_optimal is False


def test_compact_interleaved_candidate_has_first_order_coexistence() -> None:
    below = interleaved_compact_lower_bound(0.44)
    above = interleaved_compact_lower_bound(0.442)
    assert below.strategy == "no_record"
    assert below.support_value == 0.78
    assert above.strategy == "three_effect_mps"
    assert above.support_value > 1.0 - 0.442 / 2.0


def test_compact_candidate_reaches_exact_interleaved_endpoint() -> None:
    point = interleaved_compact_lower_bound(1.0)
    assert point.strategy == "three_effect_mps"
    assert isclose(point.audit_probability, 1.0, abs_tol=2e-12)
    assert isclose(point.return_fidelity, 0.25, abs_tol=2e-10)
    assert isclose(point.support_value, 1.0, abs_tol=2e-12)


def test_four_effect_phase_improves_the_compact_branch() -> None:
    compact = interleaved_compact_lower_bound(0.6)
    point = interleaved_four_effect_lower_bound(0.6)
    assert point.strategy == "four_effect_mps"
    assert point.p is not None and isclose(point.p, 0.95234043, abs_tol=4e-7)
    assert point.theta is not None and isclose(point.theta, 0.02234159, abs_tol=4e-7)
    assert point.priors is not None
    assert isclose(point.support_value, 0.765898815264694, abs_tol=3e-12)
    assert point.support_value > compact.support_value + 0.01
    assert point.support_is_globally_optimal is False


def test_best_known_envelope_switches_physical_families() -> None:
    assert interleaved_best_known_lower_bound(0.5).strategy == "three_effect_mps"
    assert interleaved_best_known_lower_bound(0.6).strategy == "four_effect_mps"


def test_four_slot_order_classification_has_both_connectivity_types() -> None:
    vectors = (1, 2, 3)

    def orbit(sequence: tuple[int, ...]) -> set[tuple[int, ...]]:
        images: set[tuple[int, ...]] = set()
        for permutation in itertools.permutations(vectors):
            relabel = dict(zip(vectors, permutation))
            image = tuple(relabel[column] for column in sequence)
            images.update((image, image[::-1]))
        return images

    sequences = tuple(
        sequence
        for sequence in itertools.product(vectors, repeat=4)
        if len(set(sequence)) >= 2
    )
    seen: set[tuple[int, ...]] = set()
    taus: list[int] = []
    for sequence in sequences:
        if sequence in seen:
            continue
        seen.update(orbit(sequence))
        matrix = (
            tuple((column >> 1) & 1 for column in sequence),
            tuple(column & 1 for column in sequence),
        )
        taus.append(trellis_connectivity_tau(matrix))

    assert len(sequences) == 78
    assert len(taus) == 9
    assert taus.count(1) == 4
    assert taus.count(2) == 5

"""Tests for exact order-sensitive AUDIT--RETURN utilities."""

import itertools
from math import isclose, sqrt

import numpy as np
import pytest

from carmenq.order_sensitive import (
    GROUPED_CHECK_MATRIX,
    INTERLEAVED_CHECK_MATRIX,
    INTERLEAVED_PERFECT_AUDIT_ENDPOINT,
    full_crossing_cuts,
    gf2_rank,
    grouped_frontier,
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

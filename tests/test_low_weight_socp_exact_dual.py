"""Checks for the capped-simplex reduction in the low-weight sector."""

from __future__ import annotations

from fractions import Fraction
import itertools
from pathlib import Path
import random
import sys


ROOT = Path(__file__).resolve().parents[1]
FRONTIER = ROOT / "scratch" / "d2_frontier"
if str(FRONTIER) not in sys.path:
    sys.path.insert(0, str(FRONTIER))

from four_active_socp_exact_cover import decode_dual  # noqa: E402
from low_weight_socp_exact_dual import (  # noqa: E402
    CAP_VECTOR,
    MAXIMUM_EFFECT_WEIGHT,
    TARGET,
    canonical_data,
    certify_cell,
    exact_upper,
)
from ternary_socp_exact_dual_probe import repair_dual_cones  # noqa: E402


def dot(left: tuple[Fraction, ...], right: tuple[Fraction, ...]) -> Fraction:
    return sum((a * b for a, b in zip(left, right, strict=True)), Fraction(0))


def test_capped_simplex_support_function_is_exact() -> None:
    vertices = sorted(set(itertools.permutations(CAP_VECTOR)))
    assert len(vertices) == 12
    random_source = random.Random(20260826)
    for _ in range(100):
        query = tuple(
            Fraction(random_source.randrange(1001), 1000) for _ in range(4)
        )
        sorted_query = tuple(sorted(query, reverse=True))
        expected = dot(CAP_VECTOR, sorted_query)
        assert max(dot(vertex, query) for vertex in vertices) == expected

        coefficients = [Fraction(random_source.randrange(1001), 1000) for _ in vertices]
        normalisation = sum(coefficients, Fraction(0))
        if normalisation == 0:
            continue
        feasible = tuple(
            sum(
                coefficient * vertex[index]
                for coefficient, vertex in zip(coefficients, vertices, strict=True)
            )
            / normalisation
            for index in range(4)
        )
        assert sum(feasible, Fraction(0)) == 2
        assert max(feasible) <= MAXIMUM_EFFECT_WEIGHT
        assert dot(feasible, query) <= expected


def test_one_low_weight_cell_has_an_exact_residual_certificate() -> None:
    order = (0, 1, 2, 3)
    report = certify_cell(order, order, TARGET)
    assert report["closed"] is True
    assert Fraction(*report["certified_upper_fraction"]) <= TARGET
    assert report["trusted_optimizers"] == []
    data = canonical_data(order, order)
    decoded = decode_dual(
        report["dual_zlib_base64"], report["dual_storage_dtype"]
    )
    repaired, _ = repair_dual_cones(decoded, data["dims"])
    replayed, correction, _ = exact_upper(data, repaired)
    assert replayed == Fraction(*report["certified_upper_fraction"])
    assert correction == Fraction(*report["exact_residual_correction"])

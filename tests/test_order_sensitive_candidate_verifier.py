"""Independent full-instrument checks for the interleaved candidate."""

from math import isclose

from scripts.verify_interleaved_candidate import direct_scores


def test_full_candidate_contraction_matches_balanced_point() -> None:
    audit, returned, completeness = direct_scores(
        0.6168956030718684, 0.8003177036431812
    )
    assert isclose(audit, 0.6446434022644625, abs_tol=2e-15)
    assert isclose(returned, 0.8662314901930314, abs_tol=2e-15)
    assert completeness < 3e-15


def test_full_candidate_contraction_reaches_exact_endpoint() -> None:
    audit, returned, completeness = direct_scores(1.0, 2**-0.5)
    assert isclose(audit, 1.0, abs_tol=2e-15)
    assert isclose(returned, 0.25, abs_tol=2e-15)
    assert completeness < 3e-15

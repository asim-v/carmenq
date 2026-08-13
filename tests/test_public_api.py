"""Tests for the stable CARMEN-Q entry points."""

from carmenq import (
    __version__,
    certify,
    certify_classical_memory,
    collective_bound,
    collective_classical_record_bound,
    plan,
    plan_experiment,
    streaming_bound,
)


def test_concise_aliases_match_scientific_api() -> None:
    assert streaming_bound(8, 0.5) == 0.75
    assert collective_bound(0.5) == collective_classical_record_bound(0.5)
    assert certify is certify_classical_memory
    assert plan is plan_experiment
    assert __version__ == "2.0.2"

"""Regression tests for the exact rational four-effect lower witness."""

from fractions import Fraction
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "verify_four_effect_rational_lower.py"
SPEC = importlib.util.spec_from_file_location("rational_lower", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
rational_lower = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(rational_lower)


def test_rational_lower_reaches_declared_endpoint_without_floats() -> None:
    payload = rational_lower.certificate()
    support = Fraction(*payload["support_lower_fraction"])
    declared = Fraction(*payload["declared_lower_fraction"])
    assert support >= declared == Fraction(7_658_988_152, 10_000_000_000)
    assert payload["optimiser_called"] is False
    assert payload["floating_point_used"] is False


def test_dyadic_square_root_is_outward_and_tight_to_one_step() -> None:
    value = Fraction(2, 7)
    lower = rational_lower.sqrt_down(value, bits=80)
    step = Fraction(1, 1 << 80)
    assert lower * lower <= value < (lower + step) ** 2


def test_half_angle_parameterisation_is_exactly_normalised() -> None:
    for coordinate in (
        rational_lower.ROOT_EFFECT_HALF_TAN,
        rational_lower.READOUT_HALF_TAN,
        rational_lower.SMALL_STATE_HALF_TAN,
        rational_lower.LARGE_STATE_HALF_TAN,
        rational_lower.PRIOR_HALF_TAN,
    ):
        sine, cosine = rational_lower.unit_pair(coordinate)
        assert sine * sine + cosine * cosine == 1

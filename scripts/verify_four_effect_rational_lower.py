"""Verify a rational four-effect lower witness at audit weight 3/5.

The optimized four-effect construction is normally reported with floating-
point angles.  This checker fixes nearby rational half-angle coordinates, so
every state, effect, and prior amplitude is rational.  The AUDIT contribution
is then rational exactly.  The only remaining radicals occur in the polar
RETURN value; dyadic square-root floors give a rigorous rational lower bound.

No optimizer or floating-point arithmetic enters the certificate.
"""

from __future__ import annotations

import argparse
from decimal import Decimal, localcontext
from fractions import Fraction
from math import isqrt
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "carmenq.four-effect-rational-lower.v1"
SUPPORT_WEIGHT = Fraction(3, 5)
DECLARED_LOWER = Fraction(7_658_988_152, 10_000_000_000)
SQRT_BITS = 192

# Rational half-angle coordinates selected near the numerical optimum.  For
# t in [0,1], (sin(phi), cos(phi)) = (2t, 1-t^2)/(1+t^2).
ROOT_EFFECT_HALF_TAN = Fraction("0.11048788031544")
READOUT_HALF_TAN = Fraction("0.011171259506937")
SMALL_STATE_HALF_TAN = Fraction("0.400875453787973")
LARGE_STATE_HALF_TAN = Fraction("0.0238892548513174")
PRIOR_HALF_TAN = Fraction("0.15098773005366")


def fraction_pair(value: Fraction) -> list[int]:
    return [value.numerator, value.denominator]


def decimal(value: Fraction, digits: int = 42) -> str:
    with localcontext() as context:
        context.prec = digits
        return format(Decimal(value.numerator) / Decimal(value.denominator), "f")


def unit_pair(half_tangent: Fraction) -> tuple[Fraction, Fraction]:
    """Return exact sine and cosine from a rational half-angle tangent."""

    if not Fraction(0) <= half_tangent <= Fraction(1):
        raise ValueError("half-angle coordinate must lie in [0,1]")
    denominator = 1 + half_tangent * half_tangent
    sine = 2 * half_tangent / denominator
    cosine = (1 - half_tangent * half_tangent) / denominator
    if sine * sine + cosine * cosine != 1:
        raise ArithmeticError("rational unit-pair identity failed")
    return sine, cosine


def sqrt_down(value: Fraction, bits: int = SQRT_BITS) -> Fraction:
    """Return a dyadic lower bound on ``sqrt(value)`` using integer arithmetic."""

    if value < 0:
        raise ValueError("cannot take the square root of a negative rational")
    scaled_square = (value.numerator << (2 * bits)) // value.denominator
    numerator = isqrt(scaled_square)
    result = Fraction(numerator, 1 << bits)
    step = Fraction(1, 1 << bits)
    if result * result > value or (result + step) ** 2 <= value:
        raise ArithmeticError("dyadic square-root enclosure failed")
    return result


def certificate() -> dict[str, Any]:
    root_one_p, root_p = unit_pair(ROOT_EFFECT_HALF_TAN)
    sine, cosine = unit_pair(READOUT_HALF_TAN)
    small_root, small_complement_root = unit_pair(SMALL_STATE_HALF_TAN)
    large_complement_root, large_root = unit_pair(LARGE_STATE_HALF_TAN)
    small_prior_root, large_prior_root = unit_pair(PRIOR_HALF_TAN)

    p_value = root_p * root_p
    u_small = small_root * small_root
    u_large = large_root * large_root
    small_total = small_prior_root * small_prior_root
    large_total = large_prior_root * large_prior_root

    q_small = (1 - p_value) + (2 * p_value - 1) * u_small
    q_large = (1 - p_value) + (2 * p_value - 1) * u_large
    if not 0 <= q_small <= 1 or not 0 <= q_large <= 1:
        raise ArithmeticError("coarse probabilities escaped [0,1]")

    d_small = (
        root_p * sine * small_root
        + root_one_p * cosine * small_complement_root
    ) ** 2
    d_large = (
        root_p * cosine * large_root
        + root_one_p * sine * large_complement_root
    ) ** 2
    audit = small_total * d_small + large_total * d_large

    q_small_root = sqrt_down(q_small)
    q_small_complement_root = sqrt_down(1 - q_small)
    q_large_root = sqrt_down(q_large)
    q_large_complement_root = sqrt_down(1 - q_large)
    c_small_lower = q_small_root + q_small_complement_root
    c_large_lower = q_large_root + q_large_complement_root
    return_lower = (
        small_prior_root * c_small_lower
        + large_prior_root * c_large_lower
    ) ** 2 / 4
    support_lower = (
        SUPPORT_WEIGHT * audit + (1 - SUPPORT_WEIGHT) * return_lower
    )
    if support_lower < DECLARED_LOWER:
        raise RuntimeError("rational witness does not reach the declared lower bound")

    coordinates = {
        "root_effect": ROOT_EFFECT_HALF_TAN,
        "readout": READOUT_HALF_TAN,
        "small_state": SMALL_STATE_HALF_TAN,
        "large_state": LARGE_STATE_HALF_TAN,
        "prior": PRIOR_HALF_TAN,
    }
    return {
        "schema": SCHEMA,
        "support_weight": fraction_pair(SUPPORT_WEIGHT),
        "half_angle_coordinates": {
            name: fraction_pair(value) for name, value in coordinates.items()
        },
        "sqrt_floor_bits": SQRT_BITS,
        "audit_exact_fraction": fraction_pair(audit),
        "audit_exact_decimal": decimal(audit),
        "return_lower_fraction": fraction_pair(return_lower),
        "return_lower_decimal": decimal(return_lower),
        "support_lower_fraction": fraction_pair(support_lower),
        "support_lower_decimal": decimal(support_lower),
        "declared_lower_fraction": fraction_pair(DECLARED_LOWER),
        "declared_lower_decimal": decimal(DECLARED_LOWER),
        "certified_at_least_declared_lower": True,
        "optimiser_called": False,
        "floating_point_used": False,
        "construction": (
            "symmetric four-effect bond-two Choi-MPS with exact rational "
            "unit amplitudes and local Pauli completion"
        ),
        "trust_boundary": (
            "Python arbitrary-precision integer and Fraction arithmetic, plus "
            "the analytic physical-construction lemma stated in the manuscript"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "data" / "four_effect_rational_lower_l060.json",
    )
    args = parser.parse_args()
    payload = certificate()
    rendered = json.dumps(payload, indent=2) + "\n"
    print(rendered, end="")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

"""Exact-residual dual certificate for one ternary SOCP cover cell.

Clarabel is used only to propose a canonical dual vector.  The checker then
repairs every second-order-cone block inward, interprets all binary64 data as
exact dyadic rationals, and evaluates both dual objective and stationarity
residual exactly.  Adding the physically redundant box ``0 <= x <= 1`` to
all 234 scalar variables converts the residual into the rigorous correction

    sum_i max(0, -(c + A^T z)_i).

For the canonical minimisation ``min c.x`` with ``A x + s = b``, this proves
an upper bound on the original maximisation of

    b.z + sum_i max(0, -(c + A^T z)_i).

The inellipse SOCs are separately checked by exact rational polynomial
arithmetic.  Solver status and reported objective values are diagnostic only.
"""

from __future__ import annotations

import argparse
from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import cvxpy as cp
import numpy as np
from scipy import sparse

from audit_source15818_enclosures import (
    positive_coefficient_error,
    rational_quadratic_maximum_on_unit_triangle,
    soc_polynomial_coefficients,
    subtract_coefficients,
    tangent_error,
)
from cvxpy.reductions.solvers.conic_solvers.clarabel_conif import CLARABEL
from pairwise_inellipse_box_cover import (
    box_anchor_relaxations,
    coefficientwise_box_soc_data,
    deserialise_box,
    inellipse_soc_data,
)
from ternary_probability_cone_cover import TernaryConeOracle


ROOT = Path(__file__).resolve().parent
DEFAULT_COVER = ROOT / "continuous_terminal_projective_l055cert_complete.json"


def q(value: float) -> Fraction:
    return Fraction.from_float(float(value))


def binary64_up(value: Fraction) -> float:
    candidate = float(value)
    while q(candidate) < value:
        candidate = float(np.nextafter(candidate, math.inf))
    return candidate


def safe_line_upper(
    exact_weight: Fraction,
    exact_upper: Fraction,
    encoded_weight: float,
) -> tuple[float, Fraction]:
    """Relax a support line enough to cover binary64 coefficient errors."""

    encoded = q(encoded_weight)
    encoded_complement = q(1.0 - encoded_weight)
    correction = max(Fraction(0), encoded - exact_weight) + max(
        Fraction(0), encoded_complement - (1 - exact_weight)
    )
    return binary64_up(exact_upper + correction), correction


def objective_correction(exact_weight: Fraction, encoded_weight: float) -> Fraction:
    encoded = q(encoded_weight)
    encoded_complement = q(1.0 - encoded_weight)
    return max(Fraction(0), exact_weight - encoded) + max(
        Fraction(0), (1 - exact_weight) - encoded_complement
    )


def inflate_radius_squared(radius: float, deficit: Fraction) -> tuple[float, int]:
    """Increase a dyadic radius until its exact square gains the deficit."""

    if deficit <= 0:
        return radius, 0
    old_squared = q(radius) * q(radius)
    candidate = max(radius, math.sqrt(float(old_squared + deficit)))
    steps = 0
    while q(candidate) * q(candidate) < old_squared + deficit:
        candidate = float(np.nextafter(candidate, math.inf))
        steps += 1
    return candidate, steps


def certified_terminal_family(
    alpha_bounds_float: tuple[float, float],
    beta_bounds_float: tuple[float, float],
) -> tuple[
    list[tuple[np.ndarray, np.ndarray, float] | None],
    tuple[np.ndarray, np.ndarray, float] | None,
    dict[str, Any],
]:
    """Return SOC data whose outer-enclosure property is rationally checked."""

    alpha_bounds = tuple(q(value) for value in alpha_bounds_float)
    beta_bounds = tuple(q(value) for value in beta_bounds_float)
    anchors: list[tuple[np.ndarray, np.ndarray, float] | None] = []
    reports: list[dict[str, Any]] = []
    for alpha_float, beta_float, error_float in box_anchor_relaxations(
        alpha_bounds_float, beta_bounds_float
    ):
        alpha, beta, error = q(alpha_float), q(beta_float), q(error_float)
        analytic = min(
            positive_coefficient_error(alpha, beta, alpha_bounds, beta_bounds),
            tangent_error(alpha, beta, alpha_bounds, beta_bounds),
        )
        try:
            root, shift, radius = inellipse_soc_data(
                alpha_float, beta_float, error_float
            )
        except ValueError:
            # At alpha=1 or beta=1 the completed ellipse can be singular.
            # Dropping that one redundant anchor is a valid outer relaxation.
            anchors.append(None)
            reports.append(
                {
                    "vacuous": True,
                    "post_sqrt_nextafter_steps": 0,
                    "margin": None,
                }
            )
            continue
        intended = (
            beta * beta,
            alpha * alpha,
            -2 * (alpha * beta - 2 * alpha - 2 * beta + 2),
            -2 * beta,
            -2 * alpha,
            1 - error,
        )
        actual = soc_polynomial_coefficients(root, root @ shift, radius)
        excess = rational_quadratic_maximum_on_unit_triangle(
            subtract_coefficients(actual, intended)
        )
        deficit = max(Fraction(0), excess - (error - analytic))
        radius, steps = inflate_radius_squared(radius, deficit)
        repaired = soc_polynomial_coefficients(root, root @ shift, radius)
        repaired_excess = rational_quadratic_maximum_on_unit_triangle(
            subtract_coefficients(repaired, intended)
        )
        margin = error - analytic - repaired_excess
        if margin < 0:
            raise ArithmeticError("failed to make an anchor SOC outward")
        anchors.append((root, shift, radius))
        reports.append(
            {
                "post_sqrt_nextafter_steps": steps,
                "margin": [margin.numerator, margin.denominator],
            }
        )

    lower = coefficientwise_box_soc_data(alpha_bounds_float, beta_bounds_float)
    lower_report: dict[str, Any]
    if lower is None:
        certified_lower = None
        lower_report = {"present": False, "certified_outer": True}
    else:
        root, shift, radius = lower
        actual = soc_polynomial_coefficients(root, root @ shift, radius)
        al, _ = alpha_bounds
        bl, _ = beta_bounds
        cross_minimum = min(
            -2 * (a * b - 2 * a - 2 * b + 2)
            for a in alpha_bounds
            for b in beta_bounds
        )
        ideal = (
            bl * bl,
            al * al,
            cross_minimum,
            -2 * beta_bounds[1],
            -2 * alpha_bounds[1],
            Fraction(1),
        )
        excess = rational_quadratic_maximum_on_unit_triangle(
            subtract_coefficients(actual, ideal)
        )
        radius, steps = inflate_radius_squared(radius, max(Fraction(0), excess))
        repaired = soc_polynomial_coefficients(root, root @ shift, radius)
        repaired_excess = rational_quadratic_maximum_on_unit_triangle(
            subtract_coefficients(repaired, ideal)
        )
        if repaired_excess > 0:
            raise ArithmeticError("failed to make the coefficientwise SOC outward")
        certified_lower = (root, shift, radius)
        lower_report = {
            "present": True,
            "certified_outer": True,
            "post_sqrt_nextafter_steps": steps,
            "margin": [-repaired_excess.numerator, repaired_excess.denominator],
        }
    return anchors, certified_lower, {
        "all_certified": True,
        "anchors": reports,
        "coefficientwise_lower": lower_report,
    }


def repair_dual_cones(z_raw: np.ndarray, dims: Any) -> tuple[np.ndarray, int]:
    """Return a binary64 vector lying exactly in the product dual cone."""

    if dims.exp or dims.psd or dims.p3d:
        raise ValueError("this checker currently accepts only zero, nonnegative, and SOC cones")
    z = np.asarray(z_raw, dtype=float).copy()
    if z.ndim != 1 or not np.all(np.isfinite(z)):
        raise ValueError("dual vector must be finite and one-dimensional")
    cursor = int(dims.zero)
    nonnegative_end = cursor + int(dims.nonneg)
    z[cursor:nonnegative_end] = np.maximum(z[cursor:nonnegative_end], 0.0)
    cursor = nonnegative_end
    repaired = 0
    for dimension in dims.soc:
        end = cursor + int(dimension)
        block = z[cursor:end]
        norm_squared = sum((q(value) * q(value) for value in block[1:]), Fraction(0))
        head = max(0.0, float(block[0]), math.sqrt(float(norm_squared)))
        while q(head) * q(head) < norm_squared:
            head = float(np.nextafter(head, math.inf))
        if head != block[0]:
            repaired += 1
            z[cursor] = head
        cursor = end
    if cursor != len(z):
        raise ValueError("cone dimensions do not consume the dual vector")
    return z, repaired


def exact_sparse_stationarity(
    matrix: sparse.spmatrix,
    dual: np.ndarray,
    objective: np.ndarray,
) -> tuple[list[Fraction], Fraction]:
    csc = sparse.csc_matrix(matrix)
    dual_q = [q(value) for value in dual]
    residuals: list[Fraction] = []
    correction = Fraction(0)
    for column in range(csc.shape[1]):
        value = q(float(objective[column]))
        for position in range(csc.indptr[column], csc.indptr[column + 1]):
            row = int(csc.indices[position])
            value += q(float(csc.data[position])) * dual_q[row]
        residuals.append(value)
        if value < 0:
            correction -= value
    return residuals, correction


def exact_dot(left: np.ndarray, right: np.ndarray) -> Fraction:
    return sum(
        (q(a) * q(b) for a, b in zip(left, right, strict=True)),
        Fraction(0),
    )


def fraction_decimal(value: Fraction, digits: int = 40) -> str:
    with localcontext() as context:
        context.prec = digits
        return format(Decimal(value.numerator) / Decimal(value.denominator), "f")


def canonical_hash(data: dict[str, Any]) -> str:
    matrix = sparse.csc_matrix(data["A"])
    digest = hashlib.sha256()
    for array in (
        matrix.indptr.astype("<i8", copy=False),
        matrix.indices.astype("<i8", copy=False),
        matrix.data.astype("<f8", copy=False),
        np.asarray(data["b"], dtype="<f8"),
        np.asarray(data["c"], dtype="<f8"),
    ):
        digest.update(array.tobytes())
    return digest.hexdigest()


def build_problem(box: dict[str, Any]) -> tuple[cp.Problem, dict[str, Any]]:
    exact_weight = Fraction(3, 5)
    exact_line_weight = Fraction(11, 20)
    encoded_weight = float(exact_weight)
    encoded_line_weight = float(exact_line_weight)
    line_055, line_055_correction = safe_line_upper(
        exact_line_weight, Fraction(7573, 10000), encoded_line_weight
    )
    line_060, line_060_correction = safe_line_upper(
        exact_weight, Fraction(76591, 100000), encoded_weight
    )
    oracle = TernaryConeOracle(
        encoded_weight,
        (0, 1, 2, 3),
        (),
        (),
        float(Fraction(3533, 4000)),
        line_060,
        projective_support_lines=((encoded_line_weight, line_055),),
    )
    # Assign every DPP parameter.  This numerical solve is dispensable for
    # validity; only the subsequently checked canonical dual is trusted.
    oracle.solve(deserialise_box(box), 0.0)
    anchor_data, lower_data, enclosure = certified_terminal_family(
        tuple(map(float, box["terminal_alpha"])),
        tuple(map(float, box["terminal_beta"])),
    )
    anchor_targets, lower_target = oracle.soc_parameters[0]
    for data, targets in zip(anchor_data, anchor_targets, strict=True):
        oracle.assign_soc(targets, data)
    oracle.assign_soc(lower_target, lower_data)
    variables = oracle.problem.variables()
    if sum(variable.size for variable in variables) != 234:
        raise RuntimeError("unexpected canonical variable inventory")
    if not all(bool(variable.attributes.get("nonneg")) for variable in variables):
        raise RuntimeError("every residual-controlled variable must be nonnegative")
    bounded = cp.Problem(
        oracle.problem.objective,
        [
            *oracle.problem.constraints,
            *(variable <= 1.0 for variable in variables),
        ],
    )
    return bounded, {
        "line_055_correction": line_055_correction,
        "line_060_correction": line_060_correction,
        "objective_correction": objective_correction(exact_weight, encoded_weight),
        "line_055_encoded_upper": q(line_055),
        "line_060_encoded_upper": q(line_060),
        "inellipse_audit": enclosure,
    }


def certify(box: dict[str, Any], target: Fraction, include_vector: bool) -> dict[str, Any]:
    problem, corrections = build_problem(box)
    enclosure = corrections["inellipse_audit"]
    data, _, _ = problem.get_problem_data(cp.CLARABEL)
    solver = CLARABEL()
    result = solver.solve_via_data(
        data,
        warm_start=False,
        verbose=False,
        solver_opts={
            "tol_gap_abs": 1e-11,
            "tol_gap_rel": 1e-11,
            "tol_feas": 1e-11,
            "max_iter": 500,
        },
        solver_cache=None,
    )
    dual, repaired_blocks = repair_dual_cones(np.asarray(result.z), data["dims"])
    residuals, residual_correction = exact_sparse_stationarity(
        data["A"], dual, data["c"]
    )
    dual_objective = exact_dot(np.asarray(data["b"]), dual)
    certified = (
        dual_objective
        + residual_correction
        + corrections["objective_correction"]
    )
    payload: dict[str, Any] = {
        "schema": "carmenq.ternary-socp-exact-residual-dual.v1",
        "box": box,
        "canonical_shape": [int(data["A"].shape[0]), int(data["A"].shape[1])],
        "canonical_nonzeros": int(data["A"].nnz),
        "canonical_sha256": canonical_hash(data),
        "variable_box": "all 234 canonical scalars satisfy 0 <= x_i <= 1",
        "cone_dimensions": {
            "zero": int(data["dims"].zero),
            "nonnegative": int(data["dims"].nonneg),
            "soc": list(map(int, data["dims"].soc)),
        },
        "untrusted_solver_status": str(result.status),
        "untrusted_primal_objective": float(result.obj_val),
        "untrusted_dual_objective": float(result.obj_val_dual),
        "soc_heads_repaired": repaired_blocks,
        "exact_dual_objective": [dual_objective.numerator, dual_objective.denominator],
        "exact_residual_correction": [
            residual_correction.numerator,
            residual_correction.denominator,
        ],
        "maximum_absolute_stationarity_residual": fraction_decimal(
            max(map(abs, residuals), default=Fraction(0))
        ),
        "objective_rounding_correction": [
            corrections["objective_correction"].numerator,
            corrections["objective_correction"].denominator,
        ],
        "safe_projective_lines": {
            "0.55": fraction_decimal(corrections["line_055_encoded_upper"]),
            "0.60": fraction_decimal(corrections["line_060_encoded_upper"]),
        },
        "certified_upper_fraction": [certified.numerator, certified.denominator],
        "certified_upper_decimal": fraction_decimal(certified),
        "target_fraction": [target.numerator, target.denominator],
        "target_decimal": fraction_decimal(target),
        "closed": certified <= target,
        "inellipse_exact_audit": enclosure,
        "trusted_optimizers": [],
        "untrusted_search_helpers": ["Clarabel dual-vector proposal"],
    }
    if include_vector:
        payload["repaired_dual_hex"] = [float(value).hex() for value in dual]
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cover", type=Path, default=DEFAULT_COVER)
    parser.add_argument("--target", default="0.76643")
    parser.add_argument("--leaf-index", type=int)
    parser.add_argument("--include-vector", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    cover = json.loads(args.cover.read_text(encoding="utf-8"))
    candidates = [leaf for leaf in cover["leaves"] if "raw_value" in leaf]
    if args.leaf_index is None:
        leaf = max(candidates, key=lambda item: float(item["raw_value"]))
    else:
        leaf = candidates[args.leaf_index]
    payload = certify(leaf["box"], Fraction(args.target), args.include_vector)
    rendered = json.dumps(payload, indent=2) + "\n"
    print(rendered, end="")
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    if not payload["closed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

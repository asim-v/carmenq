"""Construct an exact rational dual certificate for one stored SOCP cell.

Clarabel is used only to discover useful Lorentz-cone directions.  Each
direction is replaced by a dyadic ray that is exactly inside the cone.  A
linear program selects conic combinations of those rays, after which SymPy
solves the active linear system over the rationals.  The final checks use
``fractions.Fraction`` only: stationarity, coefficient signs, cone-ray
membership, and the strict upper bound.

The certificate concerns the canonical conic program whose floating-point
coefficients are interpreted as their exact IEEE-754 dyadic rationals.  It is
not yet a proof that all upstream geometric enclosures contain their intended
physical sets; that separate semantic bridge belongs in Lean.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import clarabel
import cvxpy as cp
import numpy as np
import scipy.linalg
import scipy.sparse as sp
import sympy as sy
from cvxpy.reductions.solvers.conic_solvers.clarabel_conif import (
    dims_to_solver_cones,
)
from scipy.optimize import linprog

from spectral_product_localizer_batch import (
    build_localisation_oracle,
    set_cell_caps,
)


SCHEMA = "carmenq.exact-socp-dual-certificate.v1"


def exact_float(value: float) -> Fraction:
    """Interpret one finite binary64 value as an exact rational."""

    if not math.isfinite(value):
        raise ValueError("certificate data must be finite")
    numerator, denominator = float(value).as_integer_ratio()
    return Fraction(numerator, denominator)


def encode_fraction(value: Fraction) -> list[int]:
    return [int(value.numerator), int(value.denominator)]


@dataclass(frozen=True)
class Ray:
    label: str
    entries: tuple[tuple[int, Fraction], ...]
    free: bool = False
    cone: str = "free"


def dyadic_soc_ray(
    block: np.ndarray,
    cursor: int,
    label: str,
    denominator_bits: int,
) -> Ray | None:
    """Round a numerical SOC direction to an exact, inward dyadic ray."""

    scale = float(np.max(np.abs(block)))
    if scale <= 1e-9:
        return None
    normalized = block / scale
    denominator = 1 << denominator_bits
    spatial = [int(round(float(value) * denominator)) for value in normalized[1:]]
    time = math.isqrt(sum(value * value for value in spatial))
    if time * time < sum(value * value for value in spatial):
        time += 1
    time += 1
    entries = [(cursor, Fraction(time, denominator))]
    entries.extend(
        (cursor + index + 1, Fraction(value, denominator))
        for index, value in enumerate(spatial)
        if value
    )
    return Ray(label, tuple(entries), cone="soc")


def build_rays(
    dual: np.ndarray,
    dimensions: Any,
    denominator_bits: int,
) -> tuple[list[Ray], list[tuple[int, int]]]:
    """Build an exact inner-polyhedral generating family for the dual cone."""

    rays: list[Ray] = []
    soc_blocks: list[tuple[int, int]] = []
    cursor = 0
    for local in range(int(dimensions.zero)):
        rays.append(Ray(f"zero:{local}", ((cursor + local, Fraction(1)),), True))
    cursor += int(dimensions.zero)
    for local in range(int(dimensions.nonneg)):
        rays.append(
            Ray(
                f"nonnegative:{local}",
                ((cursor + local, Fraction(1)),),
                cone="nonnegative",
            )
        )
    cursor += int(dimensions.nonneg)
    for block_index, raw_size in enumerate(dimensions.soc):
        size = int(raw_size)
        soc_blocks.append((cursor, size))
        main = dyadic_soc_ray(
            dual[cursor : cursor + size],
            cursor,
            f"soc:{block_index}:clarabel",
            denominator_bits,
        )
        if main is not None:
            rays.append(main)
        rays.append(Ray(f"soc:{block_index}:center", ((cursor, Fraction(1)),), cone="soc"))
        for axis in range(1, size):
            rays.append(
                Ray(
                    f"soc:{block_index}:axis+:{axis}",
                    ((cursor, Fraction(1)), (cursor + axis, Fraction(1))),
                    cone="soc",
                )
            )
            rays.append(
                Ray(
                    f"soc:{block_index}:axis-:{axis}",
                    ((cursor, Fraction(1)), (cursor + axis, Fraction(-1))),
                    cone="soc",
                )
            )
        cursor += size
    if dimensions.psd or dimensions.exp or dimensions.p3d or cursor != dual.size:
        raise ValueError("this exporter currently accepts zero/nonnegative/SOC cones only")
    return rays, soc_blocks


def floating_generators(row_count: int, rays: list[Ray]) -> sp.csc_matrix:
    rows: list[int] = []
    columns: list[int] = []
    values: list[float] = []
    for column, ray in enumerate(rays):
        for row, value in ray.entries:
            rows.append(row)
            columns.append(column)
            values.append(float(value))
    return sp.csc_matrix((values, (rows, columns)), shape=(row_count, len(rays)))


def exact_ray_column(
    matrix_csr: sp.csr_matrix,
    ray: Ray,
    variable_count: int,
) -> list[Fraction]:
    """Return the exact column ``A.T * ray`` from binary64 matrix entries."""

    result = [Fraction(0) for _ in range(variable_count)]
    for row, ray_value in ray.entries:
        start, stop = matrix_csr.indptr[row], matrix_csr.indptr[row + 1]
        for pointer in range(start, stop):
            column = int(matrix_csr.indices[pointer])
            result[column] += exact_float(float(matrix_csr.data[pointer])) * ray_value
    return result


def exact_ray_objective(right: np.ndarray, ray: Ray) -> Fraction:
    return sum(
        (exact_float(float(right[row])) * value for row, value in ray.entries),
        Fraction(0),
    )


def row_scaled_integer_system(
    columns: list[list[Fraction]],
    objective: np.ndarray,
) -> tuple[sy.Matrix, sy.Matrix]:
    """Convert a dyadic rational stationarity system to row-scaled integers."""

    variable_count = len(objective)
    matrix_rows: list[list[int]] = []
    right_rows: list[list[int]] = []
    for variable in range(variable_count):
        values = [column[variable] for column in columns]
        target = -exact_float(float(objective[variable]))
        denominator = max(
            [value.denominator for value in values] + [target.denominator]
        )
        if denominator & (denominator - 1):
            raise ValueError("canonical SOCP coefficients unexpectedly ceased to be dyadic")
        matrix_rows.append(
            [int(value.numerator * (denominator // value.denominator)) for value in values]
        )
        right_rows.append([int(target.numerator * (denominator // target.denominator))])
    return sy.Matrix(matrix_rows), sy.Matrix(right_rows)


def rationalize_active_solution(
    matrix: sp.csc_matrix,
    objective: np.ndarray,
    rays: list[Ray],
    floating_coefficients: np.ndarray,
    threshold: float,
) -> tuple[list[int], list[Fraction], dict[str, Any]]:
    """Solve the numerical LP's active stationarity system exactly."""

    active = [
        index
        for index, (ray, coefficient) in enumerate(zip(rays, floating_coefficients, strict=True))
        if ray.free or float(coefficient) > threshold
    ]
    matrix_csr = matrix.tocsr()
    columns = [
        exact_ray_column(matrix_csr, rays[index], matrix.shape[1]) for index in active
    ]
    integer_matrix, integer_right = row_scaled_integer_system(columns, objective)
    rank = int(integer_matrix.rank())
    augmented_rank = int(integer_matrix.row_join(integer_right).rank())
    if rank != augmented_rank:
        raise RuntimeError(
            f"active exact system is inconsistent: rank {rank}, augmented {augmented_rank}"
        )
    solution_set = sy.linsolve((integer_matrix, integer_right))
    raw_solution = next(iter(solution_set))
    symbols = sorted(
        set().union(*(value.free_symbols for value in raw_solution)), key=str
    )
    substitutions: dict[sy.Symbol, sy.Rational] = {}
    if symbols:
        # Fit the exact affine family to the floating basic solution, then
        # freeze the remaining degrees of freedom to nearby dyadic rationals.
        base = np.asarray(
            [float(value.subs({symbol: 0 for symbol in symbols})) for value in raw_solution]
        )
        directions = np.column_stack(
            [
                np.asarray(
                    [
                        float(
                            value.subs(
                                {
                                    **{other: 0 for other in symbols},
                                    symbol: 1,
                                }
                            )
                        )
                        for value in raw_solution
                    ]
                )
                - base
                for symbol in symbols
            ]
        )
        target = floating_coefficients[active]
        fitted, *_ = np.linalg.lstsq(directions, target - base, rcond=None)
        denominator = 1 << 40
        substitutions = {
            symbol: sy.Rational(int(round(float(value) * denominator)), denominator)
            for symbol, value in zip(symbols, fitted, strict=True)
        }
    values = [value.subs(substitutions) for value in raw_solution]
    if any(value.free_symbols for value in values):
        raise RuntimeError("failed to instantiate every exact nullspace parameter")
    fractions = [Fraction(int(value.p), int(value.q)) for value in values]
    return active, fractions, {
        "active_columns": len(active),
        "rank": rank,
        "free_parameters": len(symbols),
    }


def exact_certificate_audit(
    matrix: sp.csc_matrix,
    right: np.ndarray,
    objective: np.ndarray,
    rays: list[Ray],
    active: list[int],
    coefficients: list[Fraction],
    target: Fraction,
) -> tuple[Fraction, dict[str, Any]]:
    """Check the complete dual certificate using exact rational arithmetic."""

    if len(active) != len(coefficients):
        raise ValueError("active indices and coefficients are not parallel")
    for index, coefficient in zip(active, coefficients, strict=True):
        if not rays[index].free and coefficient < 0:
            raise ValueError(f"negative conic coefficient for {rays[index].label}")
    for ray in rays:
        if ray.cone != "soc":
            continue
        entries = dict(ray.entries)
        first_row = min(entries)
        time = entries.get(first_row, Fraction(0))
        spatial_square = sum(
            (value * value for row, value in entries.items() if row != first_row),
            Fraction(0),
        )
        if time < 0 or time * time < spatial_square:
            raise ValueError(f"invalid exact SOC ray {ray.label}")

    matrix_csr = matrix.tocsr()
    stationarity = [exact_float(float(value)) for value in objective]
    upper = Fraction(0)
    for index, coefficient in zip(active, coefficients, strict=True):
        if coefficient == 0:
            continue
        column = exact_ray_column(matrix_csr, rays[index], matrix.shape[1])
        for variable, value in enumerate(column):
            stationarity[variable] += coefficient * value
        upper += coefficient * exact_ray_objective(right, rays[index])
    nonzero = [index for index, value in enumerate(stationarity) if value]
    if nonzero:
        raise ValueError(f"exact stationarity failed in {len(nonzero)} coordinates")
    if not upper < target:
        raise ValueError(f"exact upper {upper} does not beat target {target}")
    return upper, {
        "stationarity_exact": True,
        "conic_coefficients_nonnegative": True,
        "soc_rays_exact": True,
        "strict_target": True,
    }


def write_payload(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--frontier-json", type=Path, required=True)
    parser.add_argument("--source-index", type=int, default=15818)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ray-bits", type=int, default=42)
    parser.add_argument("--active-threshold", type=float, default=1e-9)
    args = parser.parse_args()

    source_raw = args.frontier_json.read_bytes()
    source = json.loads(source_raw)
    stored = source["cells"][args.source_index]
    pattern = tuple(stored["branches"])
    oracle, cap_parameters, box, _ = build_localisation_oracle(source, pattern)
    set_cell_caps(source, pattern, tuple(stored["caps"]), cap_parameters)
    oracle.solve(box, safety=0.0, capture=False)
    data, _, _ = oracle.problem.get_problem_data(cp.CLARABEL)
    matrix = sp.csc_matrix(data[cp.settings.A])
    right = np.asarray(data[cp.settings.B], dtype=float)
    objective = np.asarray(data[cp.settings.C], dtype=float)
    dimensions = data["dims"]

    settings = clarabel.DefaultSettings()
    settings.verbose = False
    settings.tol_gap_abs = 1e-11
    settings.tol_gap_rel = 1e-11
    settings.tol_feas = 1e-11
    settings.max_iter = 2000
    solution = clarabel.DefaultSolver(
        sp.csc_matrix((objective.size, objective.size)),
        objective,
        matrix,
        right,
        dims_to_solver_cones(dimensions),
        settings,
    ).solve()
    dual = np.asarray(solution.z, dtype=float)
    rays, _ = build_rays(dual, dimensions, args.ray_bits)
    generators = floating_generators(matrix.shape[0], rays)
    stationarity = (matrix.T @ generators).tocsc()
    ray_objective = np.asarray(right @ generators).reshape(-1)
    lp = linprog(
        ray_objective,
        A_eq=stationarity,
        b_eq=-objective,
        bounds=[(None, None) if ray.free else (0.0, None) for ray in rays],
        method="highs",
        options={
            "dual_feasibility_tolerance": 1e-9,
            "primal_feasibility_tolerance": 1e-9,
        },
    )
    if not lp.success:
        raise RuntimeError(f"inner-ray LP failed: {lp.message}")
    active, coefficients, recovery = rationalize_active_solution(
        matrix,
        objective,
        rays,
        np.asarray(lp.x, dtype=float),
        args.active_threshold,
    )
    target = Fraction(379, 500)
    upper, audit = exact_certificate_audit(
        matrix, right, objective, rays, active, coefficients, target
    )
    selected = [
        {
            "ray_index": index,
            "label": rays[index].label,
            "free": rays[index].free,
            "cone": rays[index].cone,
            "entries": [
                {"row": row, "value": encode_fraction(value)}
                for row, value in rays[index].entries
            ],
            "coefficient": encode_fraction(coefficient),
        }
        for index, coefficient in zip(active, coefficients, strict=True)
        if coefficient
    ]
    payload = {
        "schema": SCHEMA,
        "scope": "canonical SOCP for stored source spectral cell",
        "source": {
            "path": str(args.frontier_json),
            "sha256": hashlib.sha256(source_raw).hexdigest(),
            "source_index": int(args.source_index),
            "source_cell": int(stored["source_cell"]),
            "branches": list(pattern),
            "caps": list(stored["caps"]),
        },
        "canonical_program": {
            "rows": int(matrix.shape[0]),
            "variables": int(matrix.shape[1]),
            "nonzeros": int(matrix.nnz),
            "zero": int(dimensions.zero),
            "nonnegative": int(dimensions.nonneg),
            "soc": [int(size) for size in dimensions.soc],
            "coefficient_semantics": "exact rational values of IEEE-754 binary64 canonical data",
        },
        "untrusted_discovery": {
            "clarabel_status": str(solution.status),
            "clarabel_upper": float(right @ dual),
            "inner_ray_lp_upper": float(lp.fun),
            "inner_ray_lp_residual_inf": float(
                np.linalg.norm(stationarity @ lp.x + objective, ord=np.inf)
            ),
        },
        "exact_certificate": {
            "target": encode_fraction(target),
            "upper": encode_fraction(upper),
            "upper_decimal": float(upper),
            "margin": encode_fraction(target - upper),
            "margin_decimal": float(target - upper),
            "selected_rays": selected,
            "recovery": recovery,
            "audit": audit,
        },
        "epistemic_status": (
            "solver-independent exact dual certificate for the serialized canonical "
            "SOCP; upstream physical enclosure semantics not yet formalized"
        ),
    }
    write_payload(args.output, payload)
    print(
        json.dumps(
            {
                "upper": float(upper),
                "margin": float(target - upper),
                "selected_rays": len(selected),
                **recovery,
                **audit,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

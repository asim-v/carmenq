"""Solver-free-checkable SOCP certificate for the low-weight readout sector.

For a four-outcome rank-one qubit POVM let ``w_i`` be the effect traces.  On
the sector ``0 <= w_i <= u`` and ``sum_i w_i = 2``, the support function of
the capped simplex is obtained by filling the largest coordinates first:

    c = (u, u, 2 - 2 u, 0).

Consequently ``w.q <= c.sort(q)`` for every probability vector ``q``.  The
same capped-simplex support function bounds the pulled-back prefix POVM.  We
enumerate all 24 orders of the prefix marginal and all 24 orders of the
syndrome marginal, so no unrecorded sorting or symmetry assumption remains.

Clarabel is used only to propose a dual vector for each SOCP.  Every stored
dual is repaired into the exact product cone and checked with rational
arithmetic.  All canonical variables are physically bounded by one; hence a
stationarity residual ``r`` costs at most ``sum_i max(0, -r_i)``.  The
separate verifier reconstructs every canonical problem and repeats this
calculation without calling an optimiser.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import itertools
import json
from pathlib import Path
from typing import Any

import cvxpy as cp
from cvxpy.reductions.solvers.conic_solvers.clarabel_conif import CLARABEL
import numpy as np

from four_active_socp_exact_cover import (
    binary64_up,
    decode_dual,
    encode_dual,
    fraction_pair,
)
from ternary_socp_exact_dual_probe import (
    canonical_hash,
    exact_dot,
    exact_sparse_stationarity,
    fraction_decimal,
    repair_dual_cones,
)


ROOT = Path(__file__).resolve().parent
SUPPORT_WEIGHT = Fraction(3, 5)
MAXIMUM_EFFECT_WEIGHT = Fraction(3533, 4000)
TARGET = Fraction(76591, 100000)
CAP_VECTOR = (
    MAXIMUM_EFFECT_WEIGHT,
    MAXIMUM_EFFECT_WEIGHT,
    2 - 2 * MAXIMUM_EFFECT_WEIGHT,
    Fraction(0),
)
ORDERS = tuple(itertools.permutations(range(4)))


def encoded_upper(value: Fraction) -> float:
    """Return a binary64 coefficient no smaller than ``value``."""

    return binary64_up(value)


def build_problem(
    prefix_order: tuple[int, int, int, int],
    syndrome_order: tuple[int, int, int, int],
) -> cp.Problem:
    """Build one exactly specified order cell of the universal relaxation."""

    if tuple(sorted(prefix_order)) != tuple(range(4)):
        raise ValueError("prefix_order is not a permutation")
    if tuple(sorted(syndrome_order)) != tuple(range(4)):
        raise ValueError("syndrome_order is not a permutation")

    path = cp.Variable((4, 4), nonneg=True, name="path")
    audit = cp.Variable(nonneg=True, name="audit")
    returned = cp.Variable(nonneg=True, name="return")
    cross = cp.Variable(120, nonneg=True, name="hellinger_cross")
    flat = cp.reshape(path, (16,), order="C")
    prefix = cp.sum(path, axis=1)
    syndrome = cp.hstack(
        [sum(path[z, z ^ s] for z in range(4)) for s in range(4)]
    )
    cap = tuple(encoded_upper(value) for value in CAP_VECTOR)

    constraints: list[cp.Constraint] = [
        cp.sum(flat) == 1,
        audit <= 1,
        returned <= 1,
        *(
            prefix[prefix_order[index]]
            >= prefix[prefix_order[index + 1]]
            for index in range(3)
        ),
        *(
            syndrome[syndrome_order[index]]
            >= syndrome[syndrome_order[index + 1]]
            for index in range(3)
        ),
        audit
        <= sum(
            cap[index] * prefix[prefix_order[index]] for index in range(4)
        ),
        audit
        <= sum(
            cap[index] * syndrome[syndrome_order[index]] for index in range(4)
        ),
    ]
    cursor = 0
    for first in range(16):
        for second in range(first + 1, 16):
            constraints.append(
                cp.SOC(
                    flat[first] + flat[second],
                    cp.hstack(
                        [2 * cross[cursor], flat[first] - flat[second]]
                    ),
                )
            )
            cursor += 1
    constraints.append(16 * returned <= cp.sum(flat) + 2 * cp.sum(cross))

    objective = cp.Maximize(
        encoded_upper(SUPPORT_WEIGHT) * audit
        + encoded_upper(1 - SUPPORT_WEIGHT) * returned
    )
    base = cp.Problem(objective, constraints)
    variables = base.variables()
    if not all(bool(variable.attributes.get("nonneg")) for variable in variables):
        raise RuntimeError("every residual-controlled variable must be nonnegative")
    # path, audit, return, and every geometric-mean variable are all at most
    # one on the physical feasible set.  Adding the redundant upper box makes
    # the exact residual correction immediate.
    return cp.Problem(
        objective,
        [*constraints, *(variable <= 1 for variable in variables)],
    )


def canonical_data(
    prefix_order: tuple[int, int, int, int],
    syndrome_order: tuple[int, int, int, int],
) -> dict[str, Any]:
    problem = build_problem(prefix_order, syndrome_order)
    data, _, _ = problem.get_problem_data(cp.CLARABEL)
    variable_count = sum(variable.size for variable in problem.variables())
    if data["A"].shape[1] != variable_count:
        raise RuntimeError("canonicalisation introduced an unbounded variable")
    return data


def exact_upper(
    data: dict[str, Any], dual: np.ndarray
) -> tuple[Fraction, Fraction, Fraction]:
    residuals, correction = exact_sparse_stationarity(
        data["A"], dual, data["c"]
    )
    upper = exact_dot(np.asarray(data["b"]), dual) + correction
    maximum_residual = max(map(abs, residuals), default=Fraction(0))
    return upper, correction, maximum_residual


def certify_cell(
    prefix_order: tuple[int, int, int, int],
    syndrome_order: tuple[int, int, int, int],
    target: Fraction,
) -> dict[str, Any]:
    data = canonical_data(prefix_order, syndrome_order)
    result = CLARABEL().solve_via_data(
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
    repaired, repaired_blocks = repair_dual_cones(
        np.asarray(result.z), data["dims"]
    )
    raw_upper, _, _ = exact_upper(data, repaired)
    storage_dtype = "f32"
    encoded_dual = encode_dual(repaired, storage_dtype)
    stored = decode_dual(encoded_dual, storage_dtype)
    stored, _ = repair_dual_cones(stored, data["dims"])
    upper, correction, residual = exact_upper(data, stored)
    if upper > target and raw_upper <= target:
        storage_dtype = "f64"
        encoded_dual = encode_dual(repaired, storage_dtype)
        stored = decode_dual(encoded_dual, storage_dtype)
        stored, _ = repair_dual_cones(stored, data["dims"])
        upper, correction, residual = exact_upper(data, stored)
    return {
        "prefix_order": list(prefix_order),
        "syndrome_order": list(syndrome_order),
        "canonical_shape": [int(data["A"].shape[0]), int(data["A"].shape[1])],
        "canonical_nonzeros": int(data["A"].nnz),
        "canonical_sha256": canonical_hash(data),
        "cone_dimensions": {
            "zero": int(data["dims"].zero),
            "nonnegative": int(data["dims"].nonneg),
            "soc": list(map(int, data["dims"].soc)),
        },
        "untrusted_solver_status": str(result.status),
        "untrusted_primal_objective": float(result.obj_val),
        "untrusted_dual_objective": float(result.obj_val_dual),
        "soc_heads_repaired": repaired_blocks,
        "dual_storage_dtype": storage_dtype,
        "dual_zlib_base64": encoded_dual,
        "certified_upper_fraction": fraction_pair(upper),
        "certified_upper_decimal": fraction_decimal(upper),
        "exact_residual_correction": fraction_pair(correction),
        "maximum_stationarity_residual_decimal": fraction_decimal(residual),
        "closed": upper <= target,
        "trusted_optimizers": [],
        "untrusted_search_helpers": ["Clarabel dual-vector proposal"],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", default=str(float(TARGET)))
    parser.add_argument("--max-cells", type=int)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    target = Fraction(args.target)
    pairs = list(itertools.product(ORDERS, repeat=2))
    if args.max_cells is not None:
        pairs = pairs[: args.max_cells]
    cells = []
    for index, (prefix_order, syndrome_order) in enumerate(pairs, start=1):
        cell = certify_cell(prefix_order, syndrome_order, target)
        cells.append(cell)
        if index % 24 == 0 or index == len(pairs):
            maximum = max(
                Fraction(*entry["certified_upper_fraction"]) for entry in cells
            )
            print(
                json.dumps(
                    {
                        "cells": index,
                        "closed": sum(bool(entry["closed"]) for entry in cells),
                        "maximum_upper": float(maximum),
                    }
                ),
                flush=True,
            )
    complete = len(cells) == len(ORDERS) ** 2 and all(
        bool(cell["closed"]) for cell in cells
    )
    maximum = max(
        (Fraction(*cell["certified_upper_fraction"]) for cell in cells),
        default=Fraction(0),
    )
    payload = {
        "schema": "carmenq.low-weight-socp-exact-dual.v1",
        "support_weight": fraction_pair(SUPPORT_WEIGHT),
        "maximum_effect_weight": fraction_pair(MAXIMUM_EFFECT_WEIGHT),
        "capped_simplex_extreme_vector": list(map(fraction_pair, CAP_VECTOR)),
        "target_fraction": fraction_pair(target),
        "target_decimal": fraction_decimal(target),
        "prefix_orders": len(ORDERS),
        "syndrome_orders": len(ORDERS),
        "expected_cells": len(ORDERS) ** 2,
        "cell_count": len(cells),
        "complete": complete,
        "maximum_certified_upper_fraction": fraction_pair(maximum),
        "maximum_certified_upper_decimal": fraction_decimal(maximum),
        "coverage": (
            "all order cells for both four-component marginals; the capped-"
            "simplex support function eliminates the unknown POVM weights"
        ),
        "proof_kernel": (
            "exact dyadic rational cone membership, stationarity, objective, "
            "and exhaustive permutation enumeration"
        ),
        "trusted_optimizers": [],
        "untrusted_search_helpers": ["Clarabel dual-vector proposals"],
        "cells": cells,
    }
    rendered = json.dumps(payload, indent=2) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        output = args.output
        if not output.is_absolute():
            output = ROOT / output
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(rendered, encoding="utf-8")
    if len(cells) == len(ORDERS) ** 2 and not complete:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

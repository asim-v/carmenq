"""Trace source-cell 15818 constraints into the canonical SOCP.

This audit complements the Lean theorem for the literal canonical data.  It
uses CVXPY's inverse maps to account for every source constraint and every
implicit variable-domain row, then independently evaluates the affine cone
slack on the zero vector and all standard-basis vectors.  The latter is a
second implementation path (expression evaluation rather than cone matrix
stuffing), so it detects sign, ordering, dropped-row, and coefficient errors.

The basis check is a reproducibility audit, not a replacement for a formal
proof of CVXPY itself.  Its status is reported separately from the exact Lean
certificate and from the physical enclosure premises.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import cvxpy as cp
import numpy as np
import scipy.sparse as sp

from export_exact_socp_certificate_lean import (
    DEFAULT_CERTIFICATE,
    SOURCE_INDEX,
    canonical_data_sha256,
    require_deterministic_hash_seed,
)
from spectral_product_localizer_batch import (
    build_localisation_oracle,
    set_cell_caps,
)


EXPECTED_REDUCTIONS = (
    "FlipObjective",
    "Dcp2Cone",
    "CvxAttr2Constr",
    "ConeMatrixStuffing",
    "CLARABEL",
)


def project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def canonical_json_sha256(value: object) -> str:
    raw = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def load_source(root: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    certificate = json.loads((root / DEFAULT_CERTIFICATE).read_text(encoding="utf-8"))
    source_path = root / Path(certificate["source"]["path"].replace("\\", "/"))
    source = json.loads(source_path.read_text(encoding="utf-8"))
    return source, source["cells"][SOURCE_INDEX]


def cone_argument_arrays(constraint: cp.Constraint) -> list[np.ndarray]:
    kind = type(constraint).__name__
    if kind in {"Zero", "NonNeg"}:
        return [np.asarray(constraint.args[0].value, dtype=float).reshape(-1, order="F")]
    if kind == "SOC":
        return [
            np.asarray(constraint.args[0].value, dtype=float).reshape(-1, order="F"),
            np.asarray(constraint.args[1].value, dtype=float).reshape(-1, order="F"),
        ]
    raise TypeError(f"unexpected canonical constraint type {kind!r}")


def evaluated_slack(constraints: list[cp.Constraint]) -> np.ndarray:
    pieces: list[np.ndarray] = []
    for constraint in constraints:
        sign = -1.0 if type(constraint).__name__ == "Zero" else 1.0
        pieces.extend(sign * array for array in cone_argument_arrays(constraint))
    return np.concatenate(pieces)


def zero_value(variable: cp.Variable) -> float | np.ndarray:
    return np.zeros(variable.shape, dtype=float, order="F") if variable.shape else 0.0


def assign_value_to_copies(
    variable_id: int,
    value: float | np.ndarray,
    copies: dict[int, list[cp.Variable]],
) -> None:
    for variable in copies.get(variable_id, []):
        variable.value = value


def provenance_manifest(
    source_constraints: list[cp.Constraint],
    canonical_constraints: list[cp.Constraint],
    source_id_map: dict[int, int],
) -> tuple[list[dict[str, object]], dict[str, int]]:
    source_position = {
        int(source_id_map[constraint.id]): index
        for index, constraint in enumerate(source_constraints)
    }
    manifest: list[dict[str, object]] = []
    cursor = 0
    origin_counts = {"source": 0, "implicit_variable_domain": 0}
    for canonical_position, constraint in enumerate(canonical_constraints):
        size = int(constraint.size)
        source_index = source_position.get(constraint.id)
        origin = "source" if source_index is not None else "implicit_variable_domain"
        origin_counts[origin] += 1
        manifest.append(
            {
                "canonical_position": canonical_position,
                "row_start": cursor,
                "row_stop": cursor + size,
                "cone": type(constraint).__name__,
                "size": size,
                "origin": origin,
                "source_constraint": source_index,
            }
        )
        cursor += size
    return manifest, origin_counts


def audit(basis_columns: int | None = None) -> dict[str, object]:
    root = project_root()
    source, stored = load_source(root)
    pattern = tuple(str(item) for item in stored["branches"])
    caps = tuple(stored["caps"])
    oracle, cap_parameters, box, _ = build_localisation_oracle(source, pattern)
    set_cell_caps(source, pattern, caps, cap_parameters)
    result = oracle.solve(box, safety=0.0, capture=False)
    if result.get("status") not in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}:
        raise RuntimeError(f"source reconstruction failed: {result!r}")

    data, chain, inverse = oracle.problem.get_problem_data(cp.CLARABEL)
    reductions = tuple(type(item).__name__ for item in chain.reductions)
    if reductions != EXPECTED_REDUCTIONS:
        raise RuntimeError(f"unexpected reduction chain: {reductions!r}")
    matrix = sp.csc_matrix(data[cp.settings.A], dtype=float)
    matrix.sum_duplicates()
    matrix.sort_indices()
    right = np.asarray(data[cp.settings.B], dtype=float).reshape(-1)
    objective = np.asarray(data[cp.settings.C], dtype=float).reshape(-1)
    dimensions = data["dims"]

    source_constraints = list(oracle.problem.constraints)
    inverse_dcp = inverse[1]
    inverse_stuffing = inverse[3]
    canonical_constraints = list(inverse_stuffing.constraints)
    manifest, origin_counts = provenance_manifest(
        source_constraints, canonical_constraints, inverse_dcp.cons_id_map
    )
    mapped_source = {
        item["source_constraint"]
        for item in manifest
        if item["source_constraint"] is not None
    }
    if mapped_source != set(range(len(source_constraints))):
        raise RuntimeError("not every source constraint has exactly one canonical origin")
    if manifest[-1]["row_stop"] != matrix.shape[0]:
        raise RuntimeError("provenance row spans do not cover the canonical matrix")

    reduced_variables = inverse_stuffing.id2var
    all_copies: dict[int, list[cp.Variable]] = {}
    for variable in [*reduced_variables.values(), *oracle.problem.variables()]:
        bucket = all_copies.setdefault(int(variable.id), [])
        if all(variable is not present for present in bucket):
            bucket.append(variable)
    values: dict[int, float | np.ndarray] = {
        int(variable_id): zero_value(variable)
        for variable_id, variable in reduced_variables.items()
    }
    for variable_id, value in values.items():
        assign_value_to_copies(variable_id, value, all_copies)

    direct_right = evaluated_slack(canonical_constraints)
    if not np.array_equal(direct_right, right):
        raise RuntimeError("direct zero-point slack disagrees with canonical b")

    requested_columns = matrix.shape[1] if basis_columns is None else basis_columns
    checked_columns = min(max(int(requested_columns), 0), matrix.shape[1])
    coordinate_owner: list[tuple[int, int]] = [(-1, -1)] * matrix.shape[1]
    for variable_id, variable in reduced_variables.items():
        offset = int(inverse_stuffing.var_offsets[variable_id])
        for local in range(int(variable.size)):
            coordinate_owner[offset + local] = (int(variable_id), local)
    if any(variable_id < 0 for variable_id, _ in coordinate_owner):
        raise RuntimeError("canonical variable offsets do not cover every coordinate")

    max_matrix_error = 0.0
    matrix_failures = 0
    max_objective_error = 0.0
    objective_failures = 0
    direct_objective_zero = -float(oracle.score.value)
    objective_zero_error = abs(direct_objective_zero)
    tolerance = 128 * np.finfo(float).eps

    previous: tuple[int, int] | None = None
    for column in range(checked_columns):
        if previous is not None:
            previous_id, previous_local = previous
            previous_value = values[previous_id]
            if isinstance(previous_value, np.ndarray):
                previous_value.reshape(-1, order="F")[previous_local] = 0.0
            else:
                previous_value = 0.0
                values[previous_id] = previous_value
            assign_value_to_copies(previous_id, previous_value, all_copies)

        variable_id, local = coordinate_owner[column]
        value = values[variable_id]
        if isinstance(value, np.ndarray):
            value.reshape(-1, order="F")[local] = 1.0
        else:
            value = 1.0
            values[variable_id] = value
        assign_value_to_copies(variable_id, value, all_copies)
        previous = (variable_id, local)

        direct_column = direct_right - evaluated_slack(canonical_constraints)
        stuffed_column = matrix.getcol(column).toarray().reshape(-1)
        error = float(np.max(np.abs(direct_column - stuffed_column)))
        scale = max(1.0, float(np.max(np.abs(stuffed_column))))
        max_matrix_error = max(max_matrix_error, error)
        if error > tolerance * scale:
            matrix_failures += 1

        direct_coefficient = -float(oracle.score.value) - direct_objective_zero
        objective_error = abs(direct_coefficient - objective[column])
        max_objective_error = max(max_objective_error, objective_error)
        if objective_error > tolerance * max(1.0, abs(float(objective[column]))):
            objective_failures += 1

    cone_order = [item["cone"] for item in manifest]
    order_ok = cone_order == sorted(
        cone_order, key={"Zero": 0, "NonNeg": 1, "SOC": 2}.__getitem__
    )
    payload: dict[str, object] = {
        "schema": "carmenq.source15818-canonicalization-audit.v1",
        "source_index": SOURCE_INDEX,
        "reductions": list(reductions),
        "source_constraints": len(source_constraints),
        "canonical_constraints": len(canonical_constraints),
        "canonical_rows": int(matrix.shape[0]),
        "canonical_variables": int(matrix.shape[1]),
        "canonical_nonzeros": int(matrix.nnz),
        "canonical_data_sha256": canonical_data_sha256(
            matrix, right, objective, dimensions
        ),
        "origin_counts": origin_counts,
        "all_source_constraints_mapped": True,
        "row_spans_complete": True,
        "cone_order_ok": order_ok,
        "manifest_sha256": canonical_json_sha256(manifest),
        "direct_evaluation": {
            "columns_checked": checked_columns,
            "right_exact": True,
            "objective_constant_abs": objective_zero_error,
            "max_matrix_abs_error": max_matrix_error,
            "matrix_failures": matrix_failures,
            "max_objective_abs_error": max_objective_error,
            "objective_failures": objective_failures,
            "relative_tolerance_factor": 128,
        },
    }
    payload["passed"] = bool(
        order_ok
        and matrix_failures == 0
        and objective_failures == 0
        and objective_zero_error <= tolerance
    )
    return payload


def main() -> None:
    require_deterministic_hash_seed()
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--basis-columns",
        type=int,
        default=None,
        help="check only the first N canonical columns (default: all)",
    )
    args = parser.parse_args()
    payload = audit(args.basis_columns)
    print(json.dumps(payload, indent=2))
    if not payload["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

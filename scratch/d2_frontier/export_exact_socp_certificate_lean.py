"""Export the source-cell 15818 exact SOCP certificate to Lean.

This program deliberately separates discovery from checking.  It rebuilds the
CVXPY canonical SOCP from the stored frontier cell, interprets every binary64
coefficient as its exact dyadic rational, and reads only the rational witness
from ``source_15818_exact_socp_certificate.json``.  It then independently
checks, using :class:`fractions.Fraction` arithmetic,

* ``A.T @ y + c = 0``;
* membership of ``y`` in the product dual cone;
* ``b.T @ y < 379 / 500``.

Only after those checks pass does it emit a self-contained Lean module.  Lean
repeats the checks from the serialized rational ``A``, ``b``, ``c`` and
witness, so neither this exporter nor a numerical solver is part of the final
proof's trusted base.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import subprocess
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable

import cvxpy as cp
import numpy as np
import scipy.sparse as sp

from spectral_product_localizer_batch import (
    build_localisation_oracle,
    set_cell_caps,
)


SOURCE_INDEX = 15818
TARGET = Fraction(379, 500)
CERTIFICATE_SCHEMA = "carmenq.exact-socp-dual-certificate.v1"
DEFAULT_CERTIFICATE = Path("scratch/d2_frontier/source_15818_exact_socp_certificate.json")
DATA_OUTPUT = Path("formal/CarmenQExact/Source15818Data.lean")
DUAL_DATA_OUTPUT = Path("formal/CarmenQExact/Source15818DualData.lean")
PROOF_OUTPUT = Path("formal/CarmenQExact/Source15818.lean")
DEFAULT_STATIONARITY_SHARD_SIZE = 1
DEFAULT_STATIONARITY_BRIDGE_SHARD_SIZE = 4
STATIONARITY_SINGLETON_COLUMNS = (247, 248, 249)
DUAL_CHUNK_SIZE = 32
NONNEGATIVE_SHARD_SIZE = 64
SOC_SHARD_SIZE = 32


def project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def exact_float(value: object) -> Fraction:
    """Return the exact rational denoted by one finite binary64 value."""

    binary64 = float(value)
    if not math.isfinite(binary64):
        raise ValueError("canonical SOCP data contain a non-finite coefficient")
    return Fraction(*binary64.as_integer_ratio())


def exact_json_fraction(raw: object, context: str) -> Fraction:
    """Parse one canonical ``[numerator, denominator]`` JSON rational."""

    if (
        not isinstance(raw, list)
        or len(raw) != 2
        or isinstance(raw[0], bool)
        or isinstance(raw[1], bool)
        or not isinstance(raw[0], int)
        or not isinstance(raw[1], int)
    ):
        raise ValueError(f"{context} is not an integer rational pair")
    numerator, denominator = raw
    if denominator <= 0:
        raise ValueError(f"{context} has a non-positive denominator")
    value = Fraction(numerator, denominator)
    if value.numerator != numerator or value.denominator != denominator:
        raise ValueError(f"{context} is not in canonical lowest terms")
    return value


def require_mapping(value: object, context: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ValueError(f"{context} must be a JSON object")
    return value


def require_list(value: object, context: str) -> list[Any]:
    if not isinstance(value, list):
        raise ValueError(f"{context} must be a JSON array")
    return value


def require_int(value: object, context: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{context} must be an integer")
    return value


def resolve_inside_root(root: Path, raw_path: object) -> Path:
    if not isinstance(raw_path, str):
        raise ValueError("certificate source.path must be a string")
    candidate = (root / Path(raw_path)).resolve()
    try:
        candidate.relative_to(root.resolve())
    except ValueError as error:
        raise ValueError("certificate source.path escapes the project root") from error
    return candidate


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical_text_sha256(data: bytes) -> str:
    """Hash tracked text independently of Git's CRLF checkout conversion."""

    return sha256_bytes(data.replace(b"\r\n", b"\n"))


def require_deterministic_hash_seed() -> None:
    """Re-execute a CLI entry point with deterministic Python hash ordering."""

    if os.environ.get("PYTHONHASHSEED") == "0":
        return
    environment = os.environ.copy()
    environment["PYTHONHASHSEED"] = "0"
    completed = subprocess.run([sys.executable, *sys.argv], env=environment)
    raise SystemExit(completed.returncode)


def reconstruct_canonical_socp(
    source: dict[str, Any], stored: dict[str, Any]
) -> tuple[sp.csc_matrix, np.ndarray, np.ndarray, Any]:
    """Rebuild cell 15818 and return its CLARABEL canonical conic data."""

    pattern = tuple(str(item) for item in require_list(stored.get("branches"), "cell branches"))
    caps = tuple(require_list(stored.get("caps"), "cell caps"))
    oracle, cap_parameters, box, _ = build_localisation_oracle(source, pattern)
    set_cell_caps(source, pattern, caps, cap_parameters)

    # ``solve`` assigns every remaining CVXPY Parameter from the stored box.
    # The numerical solution is discarded; only subsequently canonicalized
    # coefficient arrays are used.  Lean later checks the rational certificate
    # without trusting this call or its result.
    result = oracle.solve(box, safety=0.0, capture=False)
    if result.get("status") not in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}:
        raise RuntimeError(f"source cell reconstruction failed: {result!r}")
    data, _, _ = oracle.problem.get_problem_data(cp.CLARABEL)
    matrix = sp.csc_matrix(data[cp.settings.A], dtype=float)
    matrix.sum_duplicates()
    matrix.sort_indices()
    right = np.asarray(data[cp.settings.B], dtype=float).reshape(-1)
    objective = np.asarray(data[cp.settings.C], dtype=float).reshape(-1)
    dimensions = data["dims"]
    if dimensions.psd or dimensions.exp or dimensions.p3d:
        raise ValueError("source cell canonicalization is not a pure zero/nonnegative/SOC program")
    if matrix.shape != (right.size, objective.size):
        raise ValueError("canonical A, b, and c dimensions are inconsistent")
    return matrix, right, objective, dimensions


def canonical_matrix_entries(matrix: sp.csc_matrix) -> Iterable[tuple[int, int, Fraction]]:
    """Yield sparse entries in deterministic column-major order."""

    for column in range(matrix.shape[1]):
        start = int(matrix.indptr[column])
        stop = int(matrix.indptr[column + 1])
        previous_row = -1
        for pointer in range(start, stop):
            row = int(matrix.indices[pointer])
            if row <= previous_row:
                raise ValueError("canonical sparse matrix is not strictly row-sorted by column")
            previous_row = row
            yield row, column, exact_float(matrix.data[pointer])


def validate_source_and_metadata(
    root: Path,
    certificate: dict[str, Any],
    matrix: sp.csc_matrix,
    right: np.ndarray,
    objective: np.ndarray,
    dimensions: Any,
) -> tuple[Path, bytes, dict[str, Any]]:
    if certificate.get("schema") != CERTIFICATE_SCHEMA:
        raise ValueError("unsupported exact certificate schema")
    source_meta = require_mapping(certificate.get("source"), "certificate source")
    if require_int(source_meta.get("source_index"), "source index") != SOURCE_INDEX:
        raise ValueError(f"certificate is not for source index {SOURCE_INDEX}")
    source_path = resolve_inside_root(root, source_meta.get("path"))
    source_raw = source_path.read_bytes()
    claimed_hash = source_meta.get("sha256")
    if source_meta.get("sha256_semantics") != (
        "raw UTF-8 bytes with CRLF normalized to LF"
    ):
        raise ValueError("unsupported frontier source SHA-256 semantics")
    if not isinstance(claimed_hash, str) or canonical_text_sha256(source_raw) != claimed_hash:
        raise ValueError("frontier source SHA-256 does not match the certificate")
    source = require_mapping(json.loads(source_raw), "frontier source")
    cells = require_list(source.get("cells"), "frontier cells")
    if SOURCE_INDEX >= len(cells):
        raise ValueError("source index is outside the frontier cell array")
    stored = require_mapping(cells[SOURCE_INDEX], "stored source cell")
    expected_source_fields = {
        "source_cell": stored.get("source_cell"),
        "branches": stored.get("branches"),
        "caps": stored.get("caps"),
    }
    actual_source_fields = {
        key: source_meta.get(key) for key in expected_source_fields
    }
    if actual_source_fields != expected_source_fields:
        raise ValueError("certificate source-cell identity does not match the frontier")

    program = require_mapping(certificate.get("canonical_program"), "canonical program")
    expected_program = {
        "rows": int(matrix.shape[0]),
        "variables": int(matrix.shape[1]),
        "nonzeros": int(matrix.nnz),
        "zero": int(dimensions.zero),
        "nonnegative": int(dimensions.nonneg),
        "soc": [int(size) for size in dimensions.soc],
    }
    for key, expected in expected_program.items():
        if program.get(key) != expected:
            raise ValueError(
                f"canonical metadata {key!r} disagrees with reconstruction: "
                f"{program.get(key)!r} != {expected!r}"
            )
    if int(dimensions.zero) + int(dimensions.nonneg) + sum(
        int(size) for size in dimensions.soc
    ) != right.size:
        raise ValueError("canonical cone dimensions do not sum to the row count")
    if objective.size != matrix.shape[1]:
        raise ValueError("canonical objective length does not match A")
    return source_path, source_raw, stored


def parse_witness(
    certificate: dict[str, Any], row_count: int
) -> tuple[Fraction, list[tuple[Fraction, list[tuple[int, Fraction]]]], dict[str, Any]]:
    exact = require_mapping(certificate.get("exact_certificate"), "exact certificate")
    target = exact_json_fraction(exact.get("target"), "certificate target")
    if target != TARGET:
        raise ValueError(f"unexpected target {target}; expected {TARGET}")
    selected = require_list(exact.get("selected_rays"), "selected rays")
    witness: list[tuple[Fraction, list[tuple[int, Fraction]]]] = []
    ray_indices: set[int] = set()
    labels: set[str] = set()
    for position, raw_ray in enumerate(selected):
        ray = require_mapping(raw_ray, f"selected ray {position}")
        ray_index = require_int(ray.get("ray_index"), f"selected ray {position} index")
        if ray_index < 0 or ray_index in ray_indices:
            raise ValueError(f"selected ray {position} has an invalid or duplicate index")
        ray_indices.add(ray_index)
        label = ray.get("label")
        if not isinstance(label, str) or not label or label in labels:
            raise ValueError(f"selected ray {position} has an invalid or duplicate label")
        labels.add(label)
        coefficient = exact_json_fraction(
            ray.get("coefficient"), f"selected ray {position} coefficient"
        )
        if coefficient == 0:
            raise ValueError(f"selected ray {position} has a zero coefficient")
        entries_raw = require_list(ray.get("entries"), f"selected ray {position} entries")
        if not entries_raw:
            raise ValueError(f"selected ray {position} is empty")
        entries: list[tuple[int, Fraction]] = []
        for entry_position, raw_entry in enumerate(entries_raw):
            entry = require_mapping(
                raw_entry, f"selected ray {position} entry {entry_position}"
            )
            row = require_int(
                entry.get("row"), f"selected ray {position} entry {entry_position} row"
            )
            if not 0 <= row < row_count:
                raise ValueError(f"selected ray {position} has an out-of-range row")
            value = exact_json_fraction(
                entry.get("value"),
                f"selected ray {position} entry {entry_position} value",
            )
            if value == 0:
                raise ValueError(f"selected ray {position} serializes a zero entry")
            entries.append((row, value))
        witness.append((coefficient, entries))
    if not witness:
        raise ValueError("certificate has no selected rays")
    return target, witness, exact


def assemble_dual(
    row_count: int, witness: list[tuple[Fraction, list[tuple[int, Fraction]]]]
) -> list[Fraction]:
    dual = [Fraction(0) for _ in range(row_count)]
    for coefficient, entries in witness:
        for row, value in entries:
            dual[row] += coefficient * value
    return dual


def audit_dual_cone(dual: list[Fraction], dimensions: Any) -> dict[str, int]:
    zero_dim = int(dimensions.zero)
    nonnegative_dim = int(dimensions.nonneg)
    cursor = zero_dim
    for offset in range(nonnegative_dim):
        if dual[cursor + offset] < 0:
            raise ValueError(f"dual nonnegative cone fails at row {cursor + offset}")
    cursor += nonnegative_dim
    strict_soc = 0
    boundary_soc = 0
    for block_index, raw_size in enumerate(dimensions.soc):
        size = int(raw_size)
        if size <= 0:
            raise ValueError(f"SOC block {block_index} is empty")
        time = dual[cursor]
        spatial_square = sum(
            (dual[cursor + offset] ** 2 for offset in range(1, size)),
            Fraction(0),
        )
        gap = time * time - spatial_square
        if time < 0 or gap < 0:
            raise ValueError(f"dual SOC membership fails in block {block_index}")
        if time > 0 and gap > 0:
            strict_soc += 1
        else:
            boundary_soc += 1
        cursor += size
    if cursor != len(dual):
        raise ValueError("dual cone block traversal did not consume every row")
    return {
        "zero_rows": zero_dim,
        "nonnegative_rows": nonnegative_dim,
        "soc_blocks": len(dimensions.soc),
        "strict_soc_blocks": strict_soc,
        "boundary_soc_blocks": boundary_soc,
    }


def audit_stationarity(
    matrix: sp.csc_matrix, objective: np.ndarray, dual: list[Fraction]
) -> None:
    for column in range(matrix.shape[1]):
        residual = exact_float(objective[column])
        start = int(matrix.indptr[column])
        stop = int(matrix.indptr[column + 1])
        for pointer in range(start, stop):
            row = int(matrix.indices[pointer])
            residual += exact_float(matrix.data[pointer]) * dual[row]
        if residual:
            raise ValueError(
                f"exact stationarity fails in column {column}: {residual}"
            )


def sparse_right_entries(right: np.ndarray) -> list[tuple[int, Fraction]]:
    """Serialize dense canonical ``b`` as unique sorted nonzero entries."""

    entries = [
        (row, value)
        for row, raw in enumerate(right)
        if (value := exact_float(raw)) != 0
    ]
    rows = [row for row, _ in entries]
    if rows != sorted(set(rows)):
        raise ValueError("sparse right-hand-side rows are not unique and sorted")
    if any(value == 0 for _, value in entries):
        raise ValueError("sparse right-hand side contains an explicit zero")
    reconstructed = [Fraction(0) for _ in range(right.size)]
    for row, value in entries:
        if not 0 <= row < right.size:
            raise ValueError("sparse right-hand-side row is out of range")
        reconstructed[row] += value
    dense = [exact_float(raw) for raw in right]
    if reconstructed != dense:
        raise ValueError("sparse right-hand side does not reconstruct canonical b")
    return entries


def audit_bound(
    right: np.ndarray,
    right_entries: list[tuple[int, Fraction]],
    dual: list[Fraction],
    target: Fraction,
    exact: dict[str, Any],
) -> tuple[Fraction, Fraction]:
    dense_upper = sum(
        (exact_float(value) * dual[row] for row, value in enumerate(right)),
        Fraction(0),
    )
    sparse_upper = sum(
        (value * dual[row] for row, value in right_entries),
        Fraction(0),
    )
    if sparse_upper != dense_upper:
        raise ValueError("sparse and dense exact dual objectives disagree")
    upper = sparse_upper
    margin = target - upper
    if margin <= 0:
        raise ValueError(f"exact upper bound {upper} does not beat target {target}")
    serialized_upper = exact_json_fraction(exact.get("upper"), "serialized exact upper")
    serialized_margin = exact_json_fraction(exact.get("margin"), "serialized exact margin")
    if upper != serialized_upper or margin != serialized_margin:
        raise ValueError("recomputed exact upper bound disagrees with the certificate")
    return upper, margin

def update_fraction_hash(digest: Any, tag: str, indices: Iterable[int], value: Fraction) -> None:
    record = ":".join(
        [tag, *(str(index) for index in indices), str(value.numerator), str(value.denominator)]
    )
    digest.update(record.encode("ascii") + b"\n")


def canonical_data_sha256(
    matrix: sp.csc_matrix,
    right: np.ndarray,
    objective: np.ndarray,
    dimensions: Any,
) -> str:
    digest = hashlib.sha256()
    digest.update(
        (
            f"dims:{matrix.shape[0]}:{matrix.shape[1]}:{int(dimensions.zero)}:"
            f"{int(dimensions.nonneg)}:"
            + ",".join(str(int(size)) for size in dimensions.soc)
            + "\n"
        ).encode("ascii")
    )
    for row, column, value in canonical_matrix_entries(matrix):
        update_fraction_hash(digest, "A", (row, column), value)
    for row, raw in enumerate(right):
        update_fraction_hash(digest, "b", (row,), exact_float(raw))
    for column, raw in enumerate(objective):
        update_fraction_hash(digest, "c", (column,), exact_float(raw))
    return digest.hexdigest()


def lean_int(value: int) -> str:
    return str(value) if value >= 0 else f"({value})"


def lean_q(value: Fraction) -> str:
    return f"q {lean_int(value.numerator)} {value.denominator}"


def lean_array(items: Iterable[str], indent: str = "    ") -> list[str]:
    values = list(items)
    if not values:
        return ["#[]"]
    result = ["#["]
    result.extend(f"{indent}{value}," for value in values)
    result.append("]")
    return result


def emit_dual_function(lines: list[str], dual: list[Fraction], rows: int) -> int:
    """Emit a two-level, definitionally reducing ``Fin rows → ℚ`` witness."""

    if rows <= 0 or len(dual) != rows:
        raise ValueError("functional dual length does not match positive row count")
    chunks = [
        dual[start : min(start + DUAL_CHUNK_SIZE, rows)]
        for start in range(0, rows, DUAL_CHUNK_SIZE)
    ]
    if [value for chunk in chunks for value in chunk] != dual:
        raise ValueError("functional dual chunking changed witness order")
    for chunk_index, chunk in enumerate(chunks):
        lines.append(
            f"def source15818DualChunk{chunk_index:03d} : Nat → ℚ"
        )
        lines.extend(
            f"  | {offset} => {lean_q(value)}"
            for offset, value in enumerate(chunk)
        )
        lines.extend(["  | _ => 0", ""])
    lines.append("def source15818DualChunk : Nat → Nat → ℚ")
    lines.extend(
        f"  | {chunk_index} => source15818DualChunk{chunk_index:03d}"
        for chunk_index in range(len(chunks))
    )
    lines.extend(["  | _ => fun _ => 0", ""])
    lines.extend(
        [
            f"def source15818Dual (index : Fin {rows}) : ℚ :=",
            "  source15818DualChunk",
            f"    (index.val / {DUAL_CHUNK_SIZE})",
            f"    (index.val % {DUAL_CHUNK_SIZE})",
            "",
        ]
    )
    return len(chunks)


def stationarity_arithmetic_module_names(
    variable_count: int, shard_size: int
) -> list[str]:
    """Return the deterministic module names that own column-level atoms."""

    return [
        f"CarmenQExact.Source15818Stationarity{shard:02d}"
        for shard, _ in enumerate(stationarity_chunks(variable_count, shard_size))
    ]


def emit_dual_data_module(
    output: Path,
    dual: list[Fraction],
    rows: int,
) -> Path:
    """Emit the shared functional dual witness used by every proof shard."""

    body = [
        "private def q (numerator : Int) (denominator : Nat) : ℚ :=",
        "  (numerator : ℚ) / (denominator : ℚ)",
        "",
    ]
    emit_dual_function(body, dual, rows)
    return write_lean_module(
        output,
        ["CarmenQExact.Checker"],
        "Shared dual witness for source cell 15818",
        body,
        large_reduction=True,
    )


def emit_data_module(
    output: Path,
    source_sha256: str,
    data_sha256: str,
    matrix: sp.csc_matrix,
    right: np.ndarray,
    objective: np.ndarray,
    dimensions: Any,
    dual: list[Fraction],
    target: Fraction,
    upper: Fraction,
    margin: Fraction,
    stationarity_shard_size: int,
) -> None:
    stationarity_imports = stationarity_arithmetic_module_names(
        int(matrix.shape[1]), stationarity_shard_size
    )
    lines: list[str] = [
        "import CarmenQExact.Source15818DualData",
        *(f"import {module}" for module in stationarity_imports),
        "",
        "/-!",
        "# Exact certificate data for source spectral cell 15818",
        "",
        "This file is generated by `export_exact_socp_certificate_lean.py`.",
        "It assembles the complete rational canonical SOCP from shared generated atoms.",
        "The importing proof module checks these data with Lean kernel reduction.",
        "",
        f"Frontier source SHA-256: `{source_sha256}`",
        f"Canonical rational data SHA-256: `{data_sha256}`",
        f"Exact certified upper: `{upper.numerator}/{upper.denominator}`",
        f"Exact target margin: `{margin.numerator}/{margin.denominator}`",
        "-/",
        "",
        "namespace CarmenQExact",
        "",
        "set_option maxHeartbeats 0",
        "set_option maxRecDepth 1000000",
        "",
        "private def q (numerator : Int) (denominator : Nat) : ℚ :=",
        "  (numerator : ℚ) / (denominator : ℚ)",
        "",
    ]
    right_entries = sparse_right_entries(right)
    matrix_rows = {int(row) for row in matrix.indices}
    right_rows = sorted(row for row, _ in right_entries)
    if any(
        row < 0 or row >= matrix.shape[0]
        for row in matrix_rows | set(right_rows)
    ):
        raise ValueError("canonical sparse matrix contains an out-of-range row")
    lines.extend(
        f"private def row{row} : Fin {matrix.shape[0]} := ⟨{row}, by decide⟩"
        for row in right_rows
    )
    lines.append("")
    if len(dual) != matrix.shape[0]:
        raise ValueError("dual length does not match canonical row count")
    if len(objective) != matrix.shape[1]:
        raise ValueError("objective length does not match canonical column count")

    soc_values = ", ".join(str(int(size)) for size in dimensions.soc)
    lines.extend(
        [
            f"def source15818SocSizes : Array Nat := #[{soc_values}]",
            "",
        ]
    )

    lines.append(
        "def source15818Columns : "
        f"Array (Array (RayEntry {matrix.shape[0]})) := #["
    )
    lines.extend(
        f"  source15818Column{column:03d},"
        for column in range(matrix.shape[1])
    )
    lines.extend(["]", ""])

    lines.append(
        f"def source15818RightEntries : Array (RayEntry {matrix.shape[0]}) := #["
    )
    for row, value in right_entries:
        lines.append(
            "  { row := row"
            + str(row)
            + ", value := "
            + lean_q(value)
            + " },"
        )
    lines.extend(["]", ""])

    lines.append("def source15818Objective : Array ℚ := #[")
    lines.extend(
        f"  source15818ObjectiveCoord{column:03d},"
        for column in range(matrix.shape[1])
    )
    lines.extend(["]", ""])


    lines.extend(
        [
            "def source15818Data : CertificateData :=",
            "  { rows := " + str(matrix.shape[0]),
            "    variableCount := " + str(matrix.shape[1]),
            "    zeroDim := " + str(int(dimensions.zero)),
            "    nonnegativeDim := " + str(int(dimensions.nonneg)),
            "    socSizes := source15818SocSizes",
            "    columns := source15818Columns",
            "    rightEntries := source15818RightEntries",
            "    objective := source15818Objective",
            "    dual := source15818Dual",
            "    target := " + lean_q(target) + " }",
        ]
    )
    lines.extend(
        [
            "",
            "end CarmenQExact",
            "",
        ]
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    rendered = "\n".join(lines)
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(rendered, encoding="utf-8", newline="\n")
    temporary.replace(output)


def write_lean_module(
    path: Path,
    imports: list[str],
    title: str,
    body: list[str],
    *,
    large_reduction: bool = False,
) -> Path:
    """Atomically write one deterministic generated Lean module."""

    lines = [*(f"import {module}" for module in imports), "", "/-!", f"# {title}", "", "Generated by `export_exact_socp_certificate_lean.py`.", "-/", "", "namespace CarmenQExact", ""]
    if large_reduction:
        lines.extend(
            [
                "set_option maxHeartbeats 0",
                "set_option maxRecDepth 1000000",
                "",
            ]
        )
    lines.extend(body)
    if lines[-1] != "":
        lines.append("")
    lines.extend(["end CarmenQExact", ""])
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text("\n".join(lines), encoding="utf-8", newline="\n")
    temporary.replace(path)
    return path


def lean_nat_list(values: Iterable[int]) -> str:
    return "[" + ", ".join(str(int(value)) for value in values) + "]"


def integer_chunks(total: int, size: int) -> list[tuple[int, int]]:
    if total < 0 or size <= 0:
        raise ValueError("invalid shard dimensions")
    return [
        (start, min(size, total - start))
        for start in range(0, total, size)
    ]


def stationarity_chunks(total: int, size: int) -> list[tuple[int, int]]:
    """Partition columns, isolating empirically expensive reductions."""

    if total < 0 or size <= 0:
        raise ValueError("invalid stationarity shard dimensions")
    singleton_columns = {
        column
        for column in STATIONARITY_SINGLETON_COLUMNS
        if 0 <= column < total
    }
    chunks: list[tuple[int, int]] = []
    start = 0
    while start < total:
        if start in singleton_columns:
            count = 1
        else:
            next_singleton = min(
                (column for column in singleton_columns if column > start),
                default=total,
            )
            count = min(size, total - start, next_singleton - start)
        if count <= 0:
            raise ValueError("stationarity shard partition did not advance")
        chunks.append((start, count))
        start += count
    flattened = [
        column
        for first, count in chunks
        for column in range(first, first + count)
    ]
    if flattened != list(range(total)):
        raise ValueError("stationarity shard partition is not exact")
    return chunks


def positive_int(raw: str) -> int:
    """Parse a strictly positive integer for an argparse option."""

    value = int(raw)
    if value <= 0:
        raise argparse.ArgumentTypeError("value must be a positive integer")
    return value


def remove_stale_numbered_modules(
    directory: Path, prefix: str, expected_paths: Iterable[Path]
) -> list[Path]:
    """Delete obsolete numbered generated modules for one exact prefix."""
    extension = ".lean"
    expected_names = {path.name for path in expected_paths}
    removed: list[Path] = []
    for candidate in sorted(directory.glob(f"{prefix}*{extension}")):
        name = candidate.name
        digits = name[len(prefix) : -len(extension)]
        is_numbered = bool(digits) and all(
            "0" <= char <= "9" for char in digits
        )
        if candidate.is_file() and is_numbered and name not in expected_names:
            candidate.unlink()
            removed.append(candidate)
    return removed


def proof_dimensions_from_certificate(
    certificate: dict[str, Any],
) -> tuple[int, int, int, int, list[int]]:
    program = require_mapping(certificate.get("canonical_program"), "canonical program")
    rows = require_int(program.get("rows"), "canonical rows")
    variables = require_int(program.get("variables"), "canonical variables")
    zero_dim = require_int(program.get("zero"), "canonical zero dimension")
    nonnegative_dim = require_int(
        program.get("nonnegative"), "canonical nonnegative dimension"
    )
    soc_sizes = [
        require_int(value, f"canonical SOC size {index}")
        for index, value in enumerate(require_list(program.get("soc"), "canonical SOC sizes"))
    ]
    if min(rows, variables, zero_dim, nonnegative_dim) < 0:
        raise ValueError("canonical proof dimensions must be nonnegative")
    if any(size <= 0 for size in soc_sizes):
        raise ValueError("canonical SOC blocks must be nonempty")
    if zero_dim + nonnegative_dim + sum(soc_sizes) != rows:
        raise ValueError("canonical proof cone dimensions do not sum to rows")
    return rows, variables, zero_dim, nonnegative_dim, soc_sizes


def emit_stationarity_modules(
    directory: Path,
    matrix: sp.csc_matrix,
    objective: np.ndarray,
    dual: list[Fraction],
    variable_count: int,
    arithmetic_shard_size: int,
    bridge_shard_size: int,
) -> list[Path]:
    """Emit Data-free arithmetic shards, Data bridges, and the aggregate."""

    if matrix.shape[1] != variable_count or len(objective) != variable_count:
        raise ValueError("stationarity inputs disagree on the variable count")
    if matrix.shape[0] != len(dual):
        raise ValueError("stationarity matrix and dual lengths disagree")

    arithmetic_outputs: list[Path] = []
    arithmetic_chunks = stationarity_chunks(
        variable_count, arithmetic_shard_size
    )
    column_modules: dict[int, str] = {}
    for shard, (first, count) in enumerate(arithmetic_chunks):
        suffix = f"Source15818Stationarity{shard:02d}"
        module_name = f"CarmenQExact.{suffix}"
        shard_rows = sorted(
            {
                int(matrix.indices[pointer])
                for column in range(first, first + count)
                for pointer in range(
                    int(matrix.indptr[column]), int(matrix.indptr[column + 1])
                )
            }
        )
        body: list[str] = [
            "private def q (numerator : Int) (denominator : Nat) : ℚ :=",
            "  (numerator : ℚ) / (denominator : ℚ)",
            "",
        ]
        body.extend(
            f"private def row{row} : Fin {matrix.shape[0]} := ⟨{row}, by decide⟩"
            for row in shard_rows
        )
        body.append("")
        for column in range(first, first + count):
            column_modules[column] = module_name
            body.extend(
                [
                    f"def source15818ObjectiveCoord{column:03d} : ℚ := "
                    f"{lean_q(exact_float(objective[column]))}",
                    "",
                    f"def source15818Column{column:03d} : "
                    f"Array (RayEntry {matrix.shape[0]}) := #[",
                ]
            )
            start = int(matrix.indptr[column])
            stop = int(matrix.indptr[column + 1])
            for pointer in range(start, stop):
                row = int(matrix.indices[pointer])
                coefficient = exact_float(matrix.data[pointer])
                body.append(
                    "  { row := row"
                    + str(row)
                    + ", value := "
                    + lean_q(coefficient)
                    + " },"
                )
            body.extend(
                [
                    "]",
                    "",
                    f"def source15818StationarityValue{column:03d} : ℚ :=",
                    f"  source15818Column{column:03d}.foldl",
                    "    (fun total entry =>",
                    "      total + entry.value * source15818Dual entry.row)",
                    f"    source15818ObjectiveCoord{column:03d}",
                    "",
                    f"theorem source15818StationarityValue{column:03d}Zero :",
                    f"    source15818StationarityValue{column:03d} = 0 := by",
                    "  decide +kernel",
                    "",
                ]
            )
        arithmetic_outputs.append(
            write_lean_module(
                directory / f"{suffix}.lean",
                ["CarmenQExact.Source15818DualData"],
                f"Stationarity arithmetic shard {shard:02d} for source cell 15818",
                body,
                large_reduction=True,
            )
        )

    bridge_outputs: list[Path] = []
    bridge_module_names: list[str] = []
    for shard, (first, count) in enumerate(
        stationarity_chunks(variable_count, bridge_shard_size)
    ):
        suffix = f"Source15818StationarityBridge{shard:02d}"
        module_name = f"CarmenQExact.{suffix}"
        bridge_module_names.append(module_name)
        arithmetic_imports = list(
            dict.fromkeys(
                column_modules[column]
                for column in range(first, first + count)
            )
        )
        body = []
        for column in range(first, first + count):
            body.extend(
                [
                    f"theorem source15818StationarityColumn{column:03d} :",
                    f"    stationarityAt source15818Data {column} = 0 := by",
                    f"  change source15818StationarityValue{column:03d} = 0",
                    f"  exact source15818StationarityValue{column:03d}Zero",
                    "",
                ]
            )
        bridge_outputs.append(
            write_lean_module(
                directory / f"{suffix}.lean",
                ["CarmenQExact.Source15818Data", *arithmetic_imports],
                f"Stationarity bridge shard {shard:02d} for source cell 15818",
                body,
                large_reduction=True,
            )
        )

    aggregate = [
        "theorem source15818Stationarity :",
        "    ∀ column, column < source15818Data.variableCount →",
        "      stationarityAt source15818Data column = 0 := by",
        "  intro column hcolumn",
        f"  change column < {variable_count} at hcolumn",
        "  interval_cases column",
    ]
    aggregate.extend(
        f"  · exact source15818StationarityColumn{column:03d}"
        for column in range(variable_count)
    )
    aggregate.append("")
    aggregate_path = write_lean_module(
        directory / "Source15818Stationarity.lean",
        ["Mathlib.Tactic.IntervalCases", *bridge_module_names],
        "Aggregate stationarity for source cell 15818",
        aggregate,
        large_reduction=True,
    )

    remove_stale_numbered_modules(
        directory, "Source15818Stationarity", arithmetic_outputs
    )
    remove_stale_numbered_modules(
        directory, "Source15818StationarityBridge", bridge_outputs
    )
    return [*arithmetic_outputs, *bridge_outputs, aggregate_path]

def emit_nonnegative_modules(
    directory: Path, zero_dim: int, nonnegative_dim: int
) -> list[Path]:
    outputs: list[Path] = []
    chunks = integer_chunks(nonnegative_dim, NONNEGATIVE_SHARD_SIZE)
    module_names: list[str] = []
    for shard, (offset, count) in enumerate(chunks):
        start = zero_dim + offset
        suffix = f"Source15818Nonnegative{shard:02d}"
        module_names.append(f"CarmenQExact.{suffix}")
        body = [
            f"theorem source15818NonnegativeShard{shard:02d} :",
            "    rationalNonnegativeSegment (rationalDualAt source15818Data.dual)",
            f"      {start} {count} := by",
            "  apply nonnegativeRangeOK_sound",
            "  decide +kernel",
            "",
        ]
        outputs.append(
            write_lean_module(
                directory / f"{suffix}.lean",
                ["CarmenQExact.Source15818Data"],
                f"Nonnegative-cone shard {shard:02d} for source cell 15818",
                body,
                large_reduction=True,
            )
        )

    aggregate = [
        "theorem source15818Nonnegative :",
        "    rationalNonnegativeSegment (rationalDualAt source15818Data.dual)",
        "      source15818Data.zeroDim source15818Data.nonnegativeDim := by",
        "  change rationalNonnegativeSegment (rationalDualAt source15818Data.dual)",
        f"    {zero_dim} {nonnegative_dim}",
    ]
    if not chunks:
        aggregate.append("  intro offset hoffset")
        aggregate.append("  omega")
    else:
        first_offset, first_count = chunks[0]
        aggregate.extend(
            [
                "  have h00 : rationalNonnegativeSegment (rationalDualAt source15818Data.dual)",
                f"      {zero_dim + first_offset} {first_count} :=",
                "    source15818NonnegativeShard00",
            ]
        )
        accumulated = first_count
        previous_name = "h00"
        for shard, (offset, count) in enumerate(chunks[1:], start=1):
            current_name = f"h{shard:02d}"
            aggregate.extend(
                [
                    f"  have {current_name} : rationalNonnegativeSegment (rationalDualAt source15818Data.dual)",
                    f"      {zero_dim} {accumulated + count} := by",
                    "    exact rationalNonnegativeSegment_append",
                    f"      (rationalDualAt source15818Data.dual) {zero_dim} {accumulated} {count}",
                    f"      {previous_name} source15818NonnegativeShard{shard:02d}",
                ]
            )
            accumulated += count
            previous_name = current_name
        aggregate.append(f"  exact {previous_name}")
    aggregate.append("")
    outputs.append(
        write_lean_module(
            directory / "Source15818Nonnegative.lean",
            module_names,
            "Aggregate nonnegative cone for source cell 15818",
            aggregate,
            large_reduction=True,
        )
    )
    return outputs


def soc_chunks_with_starts(
    zero_dim: int, nonnegative_dim: int, soc_sizes: list[int]
) -> list[tuple[int, list[int]]]:
    result: list[tuple[int, list[int]]] = []
    cursor = zero_dim + nonnegative_dim
    for first, count in integer_chunks(len(soc_sizes), SOC_SHARD_SIZE):
        chunk = soc_sizes[first : first + count]
        result.append((cursor, chunk))
        cursor += sum(chunk)
    if cursor != zero_dim + nonnegative_dim + sum(soc_sizes):
        raise ValueError("SOC shard traversal did not consume every block")
    return result


def emit_soc_modules(
    directory: Path,
    rows: int,
    zero_dim: int,
    nonnegative_dim: int,
    soc_sizes: list[int],
) -> list[Path]:
    outputs: list[Path] = []
    chunks = soc_chunks_with_starts(zero_dim, nonnegative_dim, soc_sizes)
    module_names: list[str] = []
    for shard, (start, chunk) in enumerate(chunks):
        suffix = f"Source15818Soc{shard:02d}"
        module_names.append(f"CarmenQExact.{suffix}")
        body = [
            f"theorem source15818SocShard{shard:02d} :",
            "    rationalLorentzBlocks (rationalDualAt source15818Data.dual)",
            f"      {start} {lean_nat_list(chunk)} := by",
            "  apply socBlocksOK_sound",
            "  decide +kernel",
            "",
        ]
        outputs.append(
            write_lean_module(
                directory / f"{suffix}.lean",
                ["CarmenQExact.Source15818Data"],
                f"Lorentz-cone shard {shard:02d} for source cell 15818",
                body,
                large_reduction=True,
            )
        )

    base = zero_dim + nonnegative_dim
    aggregate = [
        "theorem source15818Soc :",
        "    rationalLorentzBlocks (rationalDualAt source15818Data.dual)",
        "      (source15818Data.zeroDim + source15818Data.nonnegativeDim)",
        "      source15818Data.socSizes.toList := by",
        "  change rationalLorentzBlocks (rationalDualAt source15818Data.dual)",
        f"    {base} {lean_nat_list(soc_sizes)}",
    ]
    if not chunks:
        aggregate.append("  trivial")
    else:
        first_start, first_chunk = chunks[0]
        aggregate.extend(
            [
                "  have h00 : rationalLorentzBlocks (rationalDualAt source15818Data.dual)",
                f"      {first_start} {lean_nat_list(first_chunk)} :=",
                "    source15818SocShard00",
            ]
        )
        prefix = list(first_chunk)
        previous_name = "h00"
        for shard, (_, chunk) in enumerate(chunks[1:], start=1):
            current_name = f"h{shard:02d}"
            extended = [*prefix, *chunk]
            aggregate.extend(
                [
                    f"  have {current_name} : rationalLorentzBlocks (rationalDualAt source15818Data.dual)",
                    f"      {base} {lean_nat_list(extended)} := by",
                    "    simpa using",
                    "      (rationalLorentzBlocks_append",
                    f"        (rationalDualAt source15818Data.dual) {base}",
                    f"        {lean_nat_list(prefix)} {lean_nat_list(chunk)}",
                    f"        {previous_name} source15818SocShard{shard:02d})",
                ]
            )
            prefix = extended
            previous_name = current_name
        aggregate.append(f"  exact {previous_name}")
    aggregate.append("")
    outputs.append(
        write_lean_module(
            directory / "Source15818Soc.lean",
            module_names,
            "Aggregate Lorentz cone for source cell 15818",
            aggregate,
            large_reduction=True,
        )
    )
    if base + sum(soc_sizes) != rows:
        raise ValueError("aggregate SOC blocks do not end at the row count")
    return outputs


def emit_proof_modules(
    root: Path,
    matrix: sp.csc_matrix,
    objective: np.ndarray,
    dual: list[Fraction],
    rows: int,
    variable_count: int,
    zero_dim: int,
    nonnegative_dim: int,
    soc_sizes: list[int],
    stationarity_shard_size: int,
    stationarity_bridge_shard_size: int,
) -> list[Path]:
    """Generate every small kernel shard and proposition-level aggregate."""

    directory = root / "formal/CarmenQExact"
    outputs: list[Path] = []
    outputs.append(
        write_lean_module(
            directory / "Source15818Dimensions.lean",
            ["CarmenQExact.Source15818Data"],
            "Dimension check for source cell 15818",
            [
                "theorem source15818DimensionsOK :",
                "    dimensionsOK source15818Data = true := by",
                "  decide +kernel",
                "",
            ],
            large_reduction=True,
        )
    )
    outputs.append(
        write_lean_module(
            directory / "Source15818Indices.lean",
            ["CarmenQExact.Source15818Data"],
            "Top-level index check for source cell 15818",
            [
                "theorem source15818IndicesOK :",
                "    indicesOK source15818Data = true := by",
                "  decide +kernel",
                "",
            ],
            large_reduction=True,
        )
    )
    outputs.extend(
        emit_stationarity_modules(
            directory,
            matrix,
            objective,
            dual,
            variable_count,
            stationarity_shard_size,
            stationarity_bridge_shard_size,
        )
    )
    outputs.extend(emit_nonnegative_modules(directory, zero_dim, nonnegative_dim))
    outputs.extend(
        emit_soc_modules(
            directory, rows, zero_dim, nonnegative_dim, soc_sizes
        )
    )
    outputs.append(
        write_lean_module(
            directory / "Source15818DualCone.lean",
            [
                "CarmenQExact.Source15818Nonnegative",
                "CarmenQExact.Source15818Soc",
            ],
            "Product-cone certificate for source cell 15818",
            [
                "theorem source15818DualCone :",
                "    rationalProductConeDual (rationalDualAt source15818Data.dual)",
                "      source15818Data.zeroDim source15818Data.nonnegativeDim",
                "      source15818Data.socSizes.toList := by",
                "  exact rationalProductConeDual_of_parts",
                "    (rationalDualAt source15818Data.dual) source15818Data.zeroDim",
                "    source15818Data.nonnegativeDim source15818Data.socSizes.toList",
                "    source15818Nonnegative source15818Soc",
                "",
            ],
        )
    )
    outputs.append(
        write_lean_module(
            directory / "Source15818Upper.lean",
            ["CarmenQExact.Source15818Data"],
            "Sparse dual upper bound for source cell 15818",
            [
                "theorem source15818CertifiedUpperLtTarget :",
                "    certifiedUpper source15818Data < source15818Data.target := by",
                "  decide +kernel",
                "",
            ],
            large_reduction=True,
        )
    )
    outputs.append(
        write_lean_module(
            directory / "Source15818Exact.lean",
            [
                "CarmenQExact.Source15818Dimensions",
                "CarmenQExact.Source15818Indices",
                "CarmenQExact.Source15818DualCone",
                "CarmenQExact.Source15818Stationarity",
                "CarmenQExact.Source15818Upper",
            ],
            "Exact certificate for source cell 15818",
            [
                "theorem source15818Exact : CertificateProof source15818Data :=",
                "  { dimensions := source15818DimensionsOK",
                "    indices := source15818IndicesOK",
                "    dualCone := source15818DualCone",
                "    stationarity := source15818Stationarity",
                "    upper := source15818CertifiedUpperLtTarget }",
                "",
            ],
        )
    )
    outputs.append(
        write_lean_module(
            directory / "Source15818.lean",
            [
                "CarmenQExact.EncodedWeakDuality",
                "CarmenQExact.Source15818Exact",
            ],
            "Decoded consequences for source cell 15818",
            [
                "open scoped BigOperators",
                "",
                "theorem source15818DecodedWeakDuality",
                "    (point : Fin source15818Data.variableCount → ℝ)",
                "    (slack : Fin source15818Data.rows → ℝ)",
                "    (feasible : CanonicalFeasible",
                "      (decodedMatrix source15818Data) (decodedRight source15818Data)",
                "      point slack)",
                "    (hslack : ProductConePoint (extendFin slack)",
                "      source15818Data.zeroDim source15818Data.nonnegativeDim",
                "      source15818Data.socSizes.toList) :",
                "    -(∑ column, decodedObjective source15818Data column * point column) ≤",
                "      (certifiedUpper source15818Data : ℝ) :=",
                "  exactCertificate_decoded_weak_duality",
                "    source15818Data source15818Exact point slack feasible hslack",
                "",
                "theorem source15818DecodedStrictTarget",
                "    (point : Fin source15818Data.variableCount → ℝ)",
                "    (slack : Fin source15818Data.rows → ℝ)",
                "    (feasible : CanonicalFeasible",
                "      (decodedMatrix source15818Data) (decodedRight source15818Data)",
                "      point slack)",
                "    (hslack : ProductConePoint (extendFin slack)",
                "      source15818Data.zeroDim source15818Data.nonnegativeDim",
                "      source15818Data.socSizes.toList) :",
                "    -(∑ column, decodedObjective source15818Data column * point column) <",
                "      (source15818Data.target : ℝ) :=",
                "  exactCertificate_decoded_strict_target",
                "    source15818Data source15818Exact point slack feasible hslack",
                "",
            ],
        )
    )
    return outputs


def file_record(path: Path, root: Path) -> dict[str, Any]:
    raw = path.read_bytes()
    return {
        "path": str(path.relative_to(root)),
        "bytes": len(raw),
        "sha256": sha256_bytes(raw),
    }

def main() -> None:
    require_deterministic_hash_seed()
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--proof-only",
        action="store_true",
        help=(
            "reconstruct and exactly audit the certificate, then regenerate "
            "only the sharded proof modules"
        ),
    )
    parser.add_argument(
        "--stationarity-shard-size",
        type=positive_int,
        default=DEFAULT_STATIONARITY_SHARD_SIZE,
        metavar="COLUMNS",
        help=(
            "number of arithmetic stationarity values per Data-free module "
            f"(default: {DEFAULT_STATIONARITY_SHARD_SIZE})"
        ),
    )
    parser.add_argument(
        "--stationarity-bridge-shard-size",
        type=positive_int,
        default=DEFAULT_STATIONARITY_BRIDGE_SHARD_SIZE,
        metavar="COLUMNS",
        help=(
            "number of definitional stationarity bridges per Data-importing "
            f"module (default: {DEFAULT_STATIONARITY_BRIDGE_SHARD_SIZE})"
        ),
    )
    args = parser.parse_args()
    root = project_root()
    certificate_path = root / DEFAULT_CERTIFICATE
    certificate = require_mapping(
        json.loads(certificate_path.read_bytes()), "certificate root"
    )

    output = root / DATA_OUTPUT
    source_meta = require_mapping(certificate.get("source"), "certificate source")
    source_path = resolve_inside_root(root, source_meta.get("path"))
    source_raw = source_path.read_bytes()
    source = require_mapping(json.loads(source_raw), "frontier source")
    cells = require_list(source.get("cells"), "frontier cells")
    if SOURCE_INDEX >= len(cells):
        raise ValueError("source index is outside the frontier cell array")
    stored = require_mapping(cells[SOURCE_INDEX], "stored source cell")

    matrix, right, objective, dimensions = reconstruct_canonical_socp(source, stored)
    _, verified_source_raw, _ = validate_source_and_metadata(
        root, certificate, matrix, right, objective, dimensions
    )
    target, witness, exact = parse_witness(certificate, matrix.shape[0])
    dual = assemble_dual(matrix.shape[0], witness)
    cone_summary = audit_dual_cone(dual, dimensions)
    audit_stationarity(matrix, objective, dual)
    right_entries = sparse_right_entries(right)
    upper, margin = audit_bound(right, right_entries, dual, target, exact)
    data_hash = canonical_data_sha256(matrix, right, objective, dimensions)
    source_hash = canonical_text_sha256(verified_source_raw)

    dual_output = root / DUAL_DATA_OUTPUT
    if not args.proof_only:
        emit_dual_data_module(dual_output, dual, int(matrix.shape[0]))
        emit_data_module(
            output,
            source_hash,
            data_hash,
            matrix,
            right,
            objective,
            dimensions,
            dual,
            target,
            upper,
            margin,
            args.stationarity_shard_size,
        )

    proof_outputs = emit_proof_modules(
        root,
        matrix,
        objective,
        dual,
        int(matrix.shape[0]),
        int(matrix.shape[1]),
        int(dimensions.zero),
        int(dimensions.nonneg),
        [int(size) for size in dimensions.soc],
        args.stationarity_shard_size,
        args.stationarity_bridge_shard_size,
    )
    generated = [*proof_outputs]
    if not args.proof_only:
        generated = [dual_output, output, *generated]
    records = [
        file_record(path, root)
        for path in sorted(generated, key=lambda item: str(item))
    ]
    checked_rows = sorted(
        {int(row) for row in matrix.indices} | {row for row, _ in right_entries}
    )
    arithmetic_chunks = stationarity_chunks(
        int(matrix.shape[1]), args.stationarity_shard_size
    )
    bridge_chunks = stationarity_chunks(
        int(matrix.shape[1]), args.stationarity_bridge_shard_size
    )
    print(
        json.dumps(
            {
                "proof_only": args.proof_only,
                "canonical_data_sha256": data_hash,
                "rows": int(matrix.shape[0]),
                "variables": int(matrix.shape[1]),
                "matrix_nonzeros": int(matrix.nnz),
                "right_nonzeros": len(right_entries),
                "distinct_checked_rows": len(checked_rows),
                "witness_rays": len(witness),
                "dual_domain": f"Fin {matrix.shape[0]} -> Q",
                "dual_chunk_size": DUAL_CHUNK_SIZE,
                "dual_chunks": (
                    len(dual) + DUAL_CHUNK_SIZE - 1
                ) // DUAL_CHUNK_SIZE,
                "upper": [upper.numerator, upper.denominator],
                "margin": [margin.numerator, margin.denominator],
                **cone_summary,
                "stationarity_exact": True,
                "sparse_right_exact": True,
                "strict_target": True,
                "stationarity_shards": len(arithmetic_chunks),
                "stationarity_arithmetic_shards": len(arithmetic_chunks),
                "stationarity_shard_size": args.stationarity_shard_size,
                "stationarity_bridge_shards": len(bridge_chunks),
                "stationarity_bridge_shard_size": (
                    args.stationarity_bridge_shard_size
                ),
                "stationarity_kernel_theorems": int(matrix.shape[1]),
                "nonnegative_shards": len(
                    integer_chunks(
                        int(dimensions.nonneg), NONNEGATIVE_SHARD_SIZE
                    )
                ),
                "soc_shards": len(
                    integer_chunks(len(dimensions.soc), SOC_SHARD_SIZE)
                ),
                "generated_lean_files": len(records),
                "generated_files": records,
            },
            indent=2,
        )
    )

if __name__ == "__main__":
    main()

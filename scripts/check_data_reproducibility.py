"""Compare generated research data across platforms with explicit tolerances."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np


REL_TOL = 1e-12
ABS_TOL = 1e-12


def _number(value: str) -> float | None:
    try:
        return float(value)
    except ValueError:
        return None


def _compare_csv(reference: Path, candidate: Path) -> None:
    with reference.open(encoding="utf-8", newline="") as handle:
        expected = list(csv.DictReader(handle))
    with candidate.open(encoding="utf-8", newline="") as handle:
        observed = list(csv.DictReader(handle))
    if len(expected) != len(observed):
        raise AssertionError(f"{candidate.name}: row count changed")
    if expected and observed and expected[0].keys() != observed[0].keys():
        raise AssertionError(f"{candidate.name}: columns changed")
    for row_index, (left, right) in enumerate(zip(expected, observed, strict=True), 2):
        for column in left:
            left_number = _number(left[column])
            right_number = _number(right[column])
            if left_number is not None and right_number is not None:
                if not math.isclose(
                    left_number, right_number, rel_tol=REL_TOL, abs_tol=ABS_TOL
                ):
                    raise AssertionError(
                        f"{candidate.name}:{row_index}:{column} changed "
                        f"from {left[column]} to {right[column]}"
                    )
            elif left[column] != right[column]:
                raise AssertionError(
                    f"{candidate.name}:{row_index}:{column} changed"
                )


def _compare_json_value(left: Any, right: Any, location: str) -> None:
    if isinstance(left, bool) or isinstance(right, bool):
        if left != right:
            raise AssertionError(f"{location} changed")
    elif isinstance(left, int) or isinstance(right, int):
        # Integer counts and exact rational numerators/denominators are
        # discrete data.  Comparing them exactly also avoids overflowing
        # ``math.isclose`` when a certificate contains very large integers.
        if left != right:
            raise AssertionError(f"{location} changed from {left} to {right}")
    elif isinstance(left, float) and isinstance(right, float):
        if not math.isclose(left, right, rel_tol=REL_TOL, abs_tol=ABS_TOL):
            raise AssertionError(f"{location} changed from {left} to {right}")
    elif isinstance(left, dict) and isinstance(right, dict):
        if left.keys() != right.keys():
            raise AssertionError(f"{location} keys changed")
        for key in left:
            _compare_json_value(left[key], right[key], f"{location}.{key}")
    elif isinstance(left, list) and isinstance(right, list):
        if len(left) != len(right):
            raise AssertionError(f"{location} length changed")
        for index, (left_item, right_item) in enumerate(zip(left, right, strict=True)):
            _compare_json_value(left_item, right_item, f"{location}[{index}]")
    elif left != right:
        raise AssertionError(f"{location} changed")


def _compare_json(reference: Path, candidate: Path) -> None:
    left = json.loads(reference.read_text(encoding="utf-8"))
    right = json.loads(candidate.read_text(encoding="utf-8"))
    _compare_json_value(left, right, candidate.name)


def _compare_npz(reference: Path, candidate: Path) -> None:
    with np.load(reference, allow_pickle=False) as left, np.load(
        candidate, allow_pickle=False
    ) as right:
        if left.files != right.files:
            raise AssertionError(f"{candidate.name}: array keys changed")
        for key in left.files:
            expected = left[key]
            observed = right[key]
            if expected.shape != observed.shape or expected.dtype != observed.dtype:
                raise AssertionError(f"{candidate.name}:{key} structure changed")
            if np.issubdtype(expected.dtype, np.number):
                np.testing.assert_allclose(
                    observed,
                    expected,
                    rtol=REL_TOL,
                    atol=ABS_TOL,
                    err_msg=f"{candidate.name}:{key} changed",
                )
            else:
                np.testing.assert_array_equal(observed, expected)


def compare_directories(reference: Path, candidate: Path) -> None:
    expected_names = {path.name for path in reference.iterdir() if path.is_file()}
    observed_names = {path.name for path in candidate.iterdir() if path.is_file()}
    if expected_names != observed_names:
        raise AssertionError("generated data file set changed")
    for name in sorted(expected_names):
        left = reference / name
        right = candidate / name
        if left.suffix == ".csv":
            _compare_csv(left, right)
        elif left.suffix == ".json":
            _compare_json(left, right)
        elif left.suffix == ".npz":
            _compare_npz(left, right)
        elif left.read_bytes() != right.read_bytes():
            raise AssertionError(f"{name}: bytes changed")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    args = parser.parse_args()
    compare_directories(args.reference, args.candidate)
    print(
        "Data reproducibility check passed "
        f"(rtol={REL_TOL:g}, atol={ABS_TOL:g})."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

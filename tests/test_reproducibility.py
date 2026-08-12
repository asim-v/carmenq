"""Smoke test for deterministic artifact generation."""

from __future__ import annotations

import csv
from pathlib import Path

from reversible_histories.experiments import benchmark_suite


def test_finite_shot_benchmark_is_seed_reproducible(tmp_path: Path) -> None:
    first_data = tmp_path / "first_data"
    first_figures = tmp_path / "first_figures"
    second_data = tmp_path / "second_data"
    second_figures = tmp_path / "second_figures"
    benchmark_suite(first_data, first_figures, seed=17, shots=256)
    benchmark_suite(second_data, second_figures, seed=17, shots=256)

    with (first_data / "benchmark_summary.csv").open(encoding="utf-8") as handle:
        first = list(csv.DictReader(handle))
    with (second_data / "benchmark_summary.csv").open(encoding="utf-8") as handle:
        second = list(csv.DictReader(handle))
    assert first == second
    assert (first_figures / "protocol_controls.pdf").read_bytes() == (
        second_figures / "protocol_controls.pdf"
    ).read_bytes()
    assert (first_figures / "protocol_controls.png").read_bytes() == (
        second_figures / "protocol_controls.png"
    ).read_bytes()

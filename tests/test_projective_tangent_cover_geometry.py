"""Structural checks for the finite projective tangent covers."""

from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
FRONTIER = ROOT / "scratch" / "d2_frontier"
if str(FRONTIER) not in sys.path:
    sys.path.insert(0, str(FRONTIER))

import certify_rank_rank_tangent_cover as rank_rank  # noqa: E402
import certify_remaining_projective_tangent_cover as remaining  # noqa: E402
from merge_rank_rank_tangent_covers import merge_payloads  # noqa: E402
from projective_tangent_interval_certificate import (  # noqa: E402
    Interval,
)
from carmenq.projective_secular import (  # noqa: E402
    endpoint_split_tangent_upper,
    rank_split_tangent_upper,
)


def test_geometric_cell_counts_and_names_are_exact() -> None:
    rank_cells = [
        child
        for leaf in rank_rank.geometric_leaves()
        for child in rank_rank.expand_full_angles(leaf)
    ]
    assert len(rank_rank.geometric_leaves()) == 254
    assert len(rank_cells) == 944
    assert len({cell["certificate_name"] for cell in rank_cells}) == len(rank_cells)

    remaining_cells = remaining.topology_cells()
    counts = Counter(
        (cell["first_kind"], cell["second_kind"])
        for cell in remaining_cells
    )
    assert len(remaining.geometric_leaves()) == 56
    assert counts == {
        ("endpoint", "endpoint"): 56,
        ("endpoint", "rank"): 224,
        ("rank", "endpoint"): 224,
    }
    assert len({cell["certificate_name"] for cell in remaining_cells}) == 504


def test_exact_cell_selection_preserves_catalogue_order() -> None:
    cells = [
        child
        for leaf in rank_rank.geometric_leaves()
        for child in rank_rank.expand_full_angles(leaf)
    ]
    names = [
        "coarse-trace__x89_y01__tangent_a00",
        "angle-exact__x9294_y02_a00",
    ]
    selected = rank_rank.select_named_cells(cells, names)
    assert [cell["certificate_name"] for cell in selected] == names


def test_exact_cell_selection_rejects_unknown_identifiers() -> None:
    with pytest.raises(ValueError, match="unknown certificate names"):
        rank_rank.select_named_cells([], ["does-not-exist"])


def test_merge_prefers_a_complete_cell_over_a_larger_open_checkpoint() -> None:
    cell = next(
        child
        for leaf in rank_rank.geometric_leaves()
        for child in rank_rank.expand_full_angles(leaf)
    )
    base = {
        **{key: value for key, value in cell.items() if key != "bounds"},
        "bounds": {key: list(value) for key, value in cell["bounds"].items()},
    }
    open_record = {
        **base,
        "certificate": {
            "complete": False,
            "boxes_split": 20,
            "open_frontier": [{"upper": 1.01, "box": {}}],
        },
    }
    complete_record = {
        **base,
        "certificate": {
            "complete": True,
            "boxes_split": 30,
            "open_frontier": [],
        },
    }

    def payload(record):
        return {
            "weight": "11/20",
            "level": "7573/10000",
            "cells": [record],
        }

    merged = merge_payloads([payload(open_record), payload(complete_record)])
    assert merged["processed_cell_count"] == 1
    assert merged["complete_cell_count"] == 1
    assert merged["cells"][0]["certificate"]["complete"] is True


def test_angular_cover_encloses_exact_canonical_endpoint() -> None:
    endpoint = Interval.decimal(str(rank_rank.ANGLE_EDGES[-1]))
    assert endpoint.lo * endpoint.lo <= 0.5
    assert endpoint.hi * endpoint.hi >= 0.5


def test_complement_symmetry_exchanges_the_mixed_topologies() -> None:
    x_value = 0.91
    y_value = 0.04
    sine = 0.17
    weight = 0.6
    level = 0.76662
    tangents = (0.61, 1.23, 0.72, 1.11)
    endpoint_rank = sum(
        (
            endpoint_split_tangent_upper(
                x_value, y_value, label, tangents[label], weight, level
            )
            for label in (0, 1)
        )
    ) + sum(
        (
            rank_split_tangent_upper(
                1.0 - y_value,
                1.0 - x_value,
                sine,
                label,
                tangents[2 + label],
                weight,
                level,
            )
            for label in (0, 1)
        )
    )
    complemented_x = 1.0 - y_value
    complemented_y = 1.0 - x_value
    rank_endpoint = sum(
        (
            rank_split_tangent_upper(
                complemented_x,
                complemented_y,
                sine,
                label,
                tangents[2 + label],
                weight,
                level,
            )
            for label in (0, 1)
        )
    ) + sum(
        (
            endpoint_split_tangent_upper(
                1.0 - complemented_y,
                1.0 - complemented_x,
                label,
                tangents[label],
                weight,
                level,
            )
            for label in (0, 1)
        )
    )
    assert endpoint_rank == rank_endpoint

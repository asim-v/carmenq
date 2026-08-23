from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pytest


RESEARCH = Path(__file__).resolve().parents[1] / "scratch" / "d2_frontier"
sys.path.insert(0, str(RESEARCH))

from common_instrument_cells import (  # noqa: E402
    StateCell,
    cover_volume_residual,
    most_excessive_trace_ball_coordinate,
    partition_around_trace_balls,
    partition_one_trace_ball_coordinate,
)


def test_root_cell_covers_ordered_positive_state_coordinates() -> None:
    cell = StateCell.root((0, 1, 2, 3))
    coefficients = np.asarray(
        [
            [0.4, 0.1, -0.2, 0.0],
            [0.3, -0.1, 0.0, 0.2],
            [0.2, 0.0, 0.1, -0.1],
            [0.1, 0.03, 0.02, 0.01],
        ]
    )
    assert cell.contains(coefficients)
    assert cell.bounds[0, 0] == pytest.approx([0.25, 1.0])
    assert cell.bounds[3, 0] == pytest.approx([0.0, 0.25])


def test_trace_ball_partition_preserves_volume_and_contains_center() -> None:
    parent = StateCell.root((0, 1, 2, 3))
    reference = np.asarray(
        [
            [0.32, 0.0, 0.0, 0.0],
            [0.31, 0.0, 0.0, 0.0],
            [0.20, 0.0, 0.0, 0.0],
            [0.17, 0.0, 0.0, 0.0],
        ]
    )
    radii = np.full(4, 0.16)
    sides, center = partition_around_trace_balls(parent, reference, radii)
    assert center is not None
    assert center.contains(reference)
    assert np.all(center.maximum_trace_radii(reference) <= radii + 1e-12)
    assert abs(cover_volume_residual(parent, sides, center)) < 1e-12
    assert 1 <= len(sides) <= 32


def test_partition_covers_random_parent_points() -> None:
    generator = np.random.default_rng(20260823)
    parent = StateCell.root((0, 1, 2, 3))
    reference = parent.center
    radii = np.asarray([0.2, 0.15, 0.12, 0.1])
    sides, center = partition_around_trace_balls(parent, reference, radii)
    children = sides + (() if center is None else (center,))
    for _ in range(500):
        point = generator.uniform(parent.bounds[..., 0], parent.bounds[..., 1])
        assert any(child.contains(point) for child in children)


def test_widest_split_preserves_volume() -> None:
    parent = StateCell.root((0, 1, 2, 3))
    left, right = parent.split_widest()
    assert left.volume + right.volume == pytest.approx(parent.volume)


def test_one_coordinate_partition_has_at_most_three_children() -> None:
    parent = StateCell.root((0, 1, 2, 3))
    reference = np.zeros((4, 4))
    reference[:, 0] = [0.32, 0.31, 0.2, 0.17]
    radii = np.full(4, 0.16)
    selected = most_excessive_trace_ball_coordinate(
        parent, reference, radii
    )
    assert selected is not None
    sides, center, returned = partition_one_trace_ball_coordinate(
        parent, reference, radii, selected
    )
    assert returned == selected
    assert center is not None
    assert len(sides) <= 2
    assert abs(cover_volume_residual(parent, sides, center)) < 1e-12


def test_one_coordinate_partition_stops_inside_inner_box() -> None:
    bounds = np.zeros((4, 4, 2))
    bounds[..., 0] = -0.01
    bounds[..., 1] = 0.01
    cell = StateCell(bounds)
    sides, center, selected = partition_one_trace_ball_coordinate(
        cell, np.zeros((4, 4)), np.full(4, 0.1)
    )
    assert not sides
    assert center is cell
    assert selected is None


def test_flagged_trace_distance_upper_covers_every_sampled_pair() -> None:
    generator = np.random.default_rng(41)
    cell = StateCell.root((0, 1, 2, 3))
    scale = 0.7
    upper = cell.flagged_trace_distance_upper(0, 2, scale)
    for _ in range(500):
        first = generator.uniform(cell.bounds[0, :, 0], cell.bounds[0, :, 1])
        second = generator.uniform(cell.bounds[2, :, 0], cell.bounds[2, :, 1])
        difference = first - scale * second
        exact = max(abs(float(difference[0])), float(np.linalg.norm(difference[1:])))
        assert exact <= upper + 1e-12

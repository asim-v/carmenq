"""Trace-coordinate cells for common-instrument spatial branch-and-cut.

A cell is an axis-aligned box in the sixteen real Pauli coefficients
``(a_z, r_z)`` of four subnormalised qubit states.  Positivity, normalization,
and prefix ordering are imposed by the conic node model; the box itself only
supplies valid coordinate bounds.

``partition_around_trace_balls`` creates a disjoint-up-to-boundaries cover:
one central box lies inside every requested trace-distance ball, while an
ordered family of side boxes covers its complement exactly.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any

import numpy as np


@dataclass(frozen=True)
class StateCell:
    """Axis-aligned bounds with shape ``(4, 4, 2)``."""

    bounds: np.ndarray

    def __post_init__(self) -> None:
        value = np.asarray(self.bounds, dtype=float)
        if value.shape != (4, 4, 2):
            raise ValueError("state-cell bounds must have shape (4, 4, 2)")
        if not np.all(np.isfinite(value)):
            raise ValueError("state-cell bounds must be finite")
        if np.any(value[..., 0] > value[..., 1]):
            raise ValueError("every state-cell lower bound must not exceed its upper bound")
        object.__setattr__(self, "bounds", value.copy())

    @staticmethod
    def root(prefix_order: tuple[int, int, int, int]) -> "StateCell":
        """Return a box covering every ordered normalized state family."""
        if sorted(prefix_order) != [0, 1, 2, 3]:
            raise ValueError("prefix_order must be a permutation of four labels")
        value = np.zeros((4, 4, 2), dtype=float)
        trace_upper_by_rank = (1.0, 0.5, 1.0 / 3.0, 0.25)
        for rank, label in enumerate(prefix_order):
            lower = 0.25 if rank == 0 else 0.0
            upper = trace_upper_by_rank[rank]
            value[label, 0] = (lower, upper)
            value[label, 1:, 0] = -upper
            value[label, 1:, 1] = upper
        return StateCell(value)

    @property
    def widths(self) -> np.ndarray:
        return self.bounds[..., 1] - self.bounds[..., 0]

    @property
    def volume(self) -> float:
        return float(np.prod(self.widths))

    @property
    def center(self) -> np.ndarray:
        return self.bounds.mean(axis=-1)

    def contains(self, coefficients: np.ndarray, tolerance: float = 1e-12) -> bool:
        value = np.asarray(coefficients, dtype=float)
        if value.shape != (4, 4):
            raise ValueError("state coefficients must have shape (4, 4)")
        return bool(
            np.all(value >= self.bounds[..., 0] - tolerance)
            and np.all(value <= self.bounds[..., 1] + tolerance)
        )

    def intersect_coordinate(
        self, label: int, component: int, lower: float, upper: float
    ) -> "StateCell | None":
        value = self.bounds.copy()
        value[label, component, 0] = max(value[label, component, 0], lower)
        value[label, component, 1] = min(value[label, component, 1], upper)
        if value[label, component, 0] > value[label, component, 1] + 1e-15:
            return None
        return StateCell(value)

    def maximum_trace_radii(self, reference: np.ndarray) -> np.ndarray:
        """Maximum qubit trace norm from ``reference`` within the cell."""
        point = np.asarray(reference, dtype=float)
        if point.shape != (4, 4):
            raise ValueError("reference coefficients must have shape (4, 4)")
        radii = []
        for label in range(4):
            deviations = np.maximum(
                np.abs(self.bounds[label, :, 0] - point[label]),
                np.abs(self.bounds[label, :, 1] - point[label]),
            )
            radii.append(max(float(deviations[0]), float(np.linalg.norm(deviations[1:]))))
        return np.asarray(radii)

    def flagged_trace_distance_upper(
        self, first_label: int, second_label: int, scale: float = 1.0
    ) -> float:
        """Upper-bound ``||rho_i - scale*rho_j||_1`` throughout the cell.

        For a Hermitian qubit operator with Pauli coefficients ``(d, v)``,
        the trace norm is ``max(abs(d), ||v||_2)``.  Interval subtraction and
        coordinatewise absolute maxima therefore give a rigorous cellwise
        upper bound without a moment or Jensen relaxation.
        """

        if first_label == second_label:
            raise ValueError("the two input labels must differ")
        if not 0 <= first_label < 4 or not 0 <= second_label < 4:
            raise IndexError("input label outside 0,1,2,3")
        if not math.isfinite(scale) or scale < 0.0:
            raise ValueError("scale must be finite and nonnegative")
        lower = (
            self.bounds[first_label, :, 0]
            - scale * self.bounds[second_label, :, 1]
        )
        upper = (
            self.bounds[first_label, :, 1]
            - scale * self.bounds[second_label, :, 0]
        )
        absolute = np.maximum(np.abs(lower), np.abs(upper))
        return max(float(absolute[0]), float(np.linalg.norm(absolute[1:])))

    def split_widest(self) -> tuple["StateCell", "StateCell"]:
        """Bisect the widest absolute coefficient interval."""
        label, component = np.unravel_index(
            int(np.argmax(self.widths)), self.widths.shape
        )
        lower, upper = self.bounds[label, component]
        midpoint = 0.5 * (lower + upper)
        left = self.intersect_coordinate(label, component, lower, midpoint)
        right = self.intersect_coordinate(label, component, midpoint, upper)
        if left is None or right is None:  # pragma: no cover - finite nonempty box
            raise RuntimeError("failed to bisect a nonempty state cell")
        return left, right

    def to_json(self) -> list[list[list[float]]]:
        return self.bounds.tolist()

    @staticmethod
    def from_json(payload: Any) -> "StateCell":
        return StateCell(np.asarray(payload, dtype=float))


def partition_around_trace_balls(
    cell: StateCell,
    reference: np.ndarray,
    radii: np.ndarray,
) -> tuple[tuple[StateCell, ...], StateCell | None]:
    """Partition ``cell`` into side cells and one witness-valid center.

    For each input label, the scalar coefficient is restricted to distance at
    most ``radius`` and each of the three vector coefficients to distance at
    most ``radius/sqrt(3)``.  Therefore the central box satisfies

    ``max(abs(delta_a), ||delta_r||_2) <= radius``.

    Side cells are emitted sequentially.  Their interiors are disjoint and,
    together with the center, have the same volume as the parent cell.
    """
    point = np.asarray(reference, dtype=float)
    radius_value = np.asarray(radii, dtype=float)
    if point.shape != (4, 4) or radius_value.shape != (4,):
        raise ValueError("expected reference shape (4,4) and four radii")
    if np.any(~np.isfinite(radius_value)) or np.any(radius_value < 0.0):
        raise ValueError("trace radii must be finite and nonnegative")

    remainder: StateCell | None = cell
    sides: list[StateCell] = []
    for label in range(4):
        for component in range(4):
            if remainder is None:
                break
            threshold = (
                radius_value[label]
                if component == 0
                else radius_value[label] / math.sqrt(3.0)
            )
            central_lower = point[label, component] - threshold
            central_upper = point[label, component] + threshold
            lower, upper = remainder.bounds[label, component]
            if lower < central_lower:
                low = remainder.intersect_coordinate(
                    label, component, lower, min(upper, central_lower)
                )
                if low is not None and low.volume > 0.0:
                    sides.append(low)
            if upper > central_upper:
                high = remainder.intersect_coordinate(
                    label, component, max(lower, central_upper), upper
                )
                if high is not None and high.volume > 0.0:
                    sides.append(high)
            remainder = remainder.intersect_coordinate(
                label, component, central_lower, central_upper
            )
    return tuple(sides), remainder


def cover_volume_residual(
    parent: StateCell, sides: tuple[StateCell, ...], center: StateCell | None
) -> float:
    """Return child-volume sum minus parent volume."""
    child_volume = sum(item.volume for item in sides)
    if center is not None:
        child_volume += center.volume
    return float(child_volume - parent.volume)


def most_excessive_trace_ball_coordinate(
    cell: StateCell,
    reference: np.ndarray,
    radii: np.ndarray,
) -> tuple[int, int] | None:
    """Return the coordinate that most exceeds a trace-ball inner box."""
    point = np.asarray(reference, dtype=float)
    radius_value = np.asarray(radii, dtype=float)
    if point.shape != (4, 4) or radius_value.shape != (4,):
        raise ValueError("expected reference shape (4,4) and four radii")
    best: tuple[float, float, int, int] | None = None
    for label in range(4):
        for component in range(4):
            threshold = (
                radius_value[label]
                if component == 0
                else radius_value[label] / math.sqrt(3.0)
            )
            deviation = max(
                abs(float(cell.bounds[label, component, 0] - point[label, component])),
                abs(float(cell.bounds[label, component, 1] - point[label, component])),
            )
            excess = deviation - threshold
            if excess <= 1e-15:
                continue
            ratio = math.inf if threshold <= 1e-15 else deviation / threshold
            candidate = (ratio, excess, label, component)
            if best is None or candidate > best:
                best = candidate
    return None if best is None else (best[2], best[3])


def partition_one_trace_ball_coordinate(
    cell: StateCell,
    reference: np.ndarray,
    radii: np.ndarray,
    coordinate: tuple[int, int] | None = None,
) -> tuple[tuple[StateCell, ...], StateCell | None, tuple[int, int] | None]:
    """Split one excessive coordinate into low, central, and high cells."""
    point = np.asarray(reference, dtype=float)
    radius_value = np.asarray(radii, dtype=float)
    selected = (
        most_excessive_trace_ball_coordinate(cell, point, radius_value)
        if coordinate is None
        else coordinate
    )
    if selected is None:
        return (), cell, None
    label, component = selected
    threshold = (
        radius_value[label]
        if component == 0
        else radius_value[label] / math.sqrt(3.0)
    )
    central_lower = point[label, component] - threshold
    central_upper = point[label, component] + threshold
    lower, upper = cell.bounds[label, component]
    sides = []
    if lower < central_lower:
        low = cell.intersect_coordinate(
            label, component, lower, min(upper, central_lower)
        )
        if low is not None and low.volume > 0.0:
            sides.append(low)
    if upper > central_upper:
        high = cell.intersect_coordinate(
            label, component, max(lower, central_upper), upper
        )
        if high is not None and high.volume > 0.0:
            sides.append(high)
    center = cell.intersect_coordinate(
        label, component, central_lower, central_upper
    )
    return tuple(sides), center, selected

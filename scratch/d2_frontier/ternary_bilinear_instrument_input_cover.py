"""Spatial common-instrument cover of the worst ternary frontier cell.

The norm-tube cover in :mod:`ternary_common_instrument_input_cover` is a
rigorous outer approximation, but it discards the signs and correlations of
all products between an input Pauli coordinate and a shared Choi matrix.
This driver retains every such product with its four McCormick inequalities.

For a fixed input box the model contains four positive Choi matrices whose
sum is trace preserving.  The same matrices generate every path statistic.
The same input box is also coupled to one common effective POVM through a
second family of McCormick products.  Pure-prefix spherical caps and robust
determinant-scaled POVM witnesses further tighten this joint relaxation.
The canonical terminal POVM is enclosed effect by effect, so the remaining
terminal error is

    |q[z,y,t] - q0[z,y,t]| <= e[t] q[z,y],

where ``q0`` is the anchor-POVM prediction and ``q[z,y]`` is the trace of the
conditioned output.  When the input box has zero width, the McCormick
products become exact and the model reduces to one literal shared quantum
instrument (up to the separately enclosed terminal cell).

For a sign-definite input box, Cramer's rule turns positivity of a recovered
POVM effect into a linear inequality whose coefficients are row-replacement
determinants.  Those determinants are multi-affine, so exhaustive vertex
evaluation gives an exact box envelope.  The resulting inequality is valid
for every common POVM in the complete input cell, not merely at the current
optimizer.  Branching targets the bisection that most reduces the closest
uncertified determinant margin.  A hybrid rule limits these targeted splits
to near-robust margins and alternates them with global product-residual
splits, preventing one determinant witness from starving the rest of the
McCormick relaxation.

The JSON checkpoint is resumable.  A node is marked closed only after its
solver-conditional upper bound falls below the requested target.
"""

from __future__ import annotations

import argparse
import heapq
import itertools
import json
import math
from pathlib import Path
from typing import Any

import cvxpy as cp
import numpy as np
from scipy.optimize import minimize_scalar

from ternary_common_povm_input_cover import (
    _compact_result,
    _configuration,
    _node_payload,
    _oracle_keywords,
    _split_coordinate,
    _write_checkpoint,
    localise_candidate_region,
)
from ternary_probability_cone_cover import TernaryConeOracle
from terminal_reconstruction_enclosure import (
    Interval,
    planar_reconstruction,
    reconstruction_intervals,
    terminal_effect_anchor_and_errors,
)


DETERMINANT_NEAR_RELATIVE_GAP = 1.02
ANDO_NEAR_RELATIVE_GAP = 1.4
ANDO_SPLIT_SHORTLIST = 4
MAX_DETERMINANT_BRANCH_STREAK = 1


def box_purity_caps(lower: np.ndarray, upper: np.ndarray) -> np.ndarray:
    """Return spherical caps containing every nonzero vector in each box.

    A cap is stored as ``(normal_x, normal_y, normal_z, cosine)``.  The cone
    ``normal @ r >= cosine * ||r||`` is convex, so checking all eight box
    corners proves that the complete vector box lies in the cap.  For a pure
    subnormalised qubit state ``||r|| = p``; hence
    ``normal @ r >= cosine * p`` is a valid linear outer constraint.  A box
    not contained in one open hemisphere receives the vacuous zero cap.
    """

    lower = np.asarray(lower, dtype=float)
    upper = np.asarray(upper, dtype=float)
    if lower.shape != (4, 4) or upper.shape != (4, 4):
        raise ValueError("input bounds must have shape (4,4)")
    caps = np.zeros((4, 4), dtype=float)
    for z in range(4):
        center = 0.5 * (lower[z, 1:] + upper[z, 1:])
        norm = float(np.linalg.norm(center))
        if norm <= 1e-15:
            continue
        normal = center / norm
        cosine = 1.0
        for bits in range(8):
            corner = np.asarray(
                [
                    upper[z, axis + 1] if bits & (1 << axis) else lower[z, axis + 1]
                    for axis in range(3)
                ]
            )
            corner_norm = float(np.linalg.norm(corner))
            if corner_norm <= 1e-15:
                cosine = 0.0
                break
            cosine = min(cosine, float(normal @ corner / corner_norm))
        if cosine <= 0.0:
            continue
        caps[z, :3] = normal
        caps[z, 3] = max(0.0, float(np.nextafter(cosine, -math.inf)))
    return caps


def determinant_interval(lower: np.ndarray, upper: np.ndarray) -> Interval:
    """Outward-enclose the determinant of a four-by-four interval matrix."""

    lower = np.asarray(lower, dtype=float)
    upper = np.asarray(upper, dtype=float)
    matrix = [
        [
            Interval(float(lower[row, column]), float(upper[row, column]))
            for column in range(4)
        ]
        for row in range(4)
    ]
    result = Interval.point(0.0)
    for permutation in itertools.permutations(range(4)):
        inversions = sum(
            permutation[first] > permutation[second]
            for first in range(4)
            for second in range(first + 1, 4)
        )
        term = Interval.point(-1.0 if inversions % 2 else 1.0)
        for row in range(4):
            term = term * matrix[row][permutation[row]]
        result = result + term
    return result


def determinant_vertex_bounds(lower: np.ndarray, upper: np.ndarray) -> Interval:
    """Return guarded exact extrema of a determinant on a matrix box.

    A determinant is multi-affine in the sixteen matrix entries.  Its minimum
    and maximum on a rectangular box are therefore attained among the
    ``2**16`` vertices.  The exhaustive vectorized evaluation removes the
    dependency inflation of ordinary interval arithmetic, which is decisive
    for detecting sign-definite input-state bases near the frontier.
    """

    lower = np.asarray(lower, dtype=float)
    upper = np.asarray(upper, dtype=float)
    if lower.shape != (4, 4) or upper.shape != (4, 4):
        raise ValueError("input bounds must have shape (4,4)")
    corner_rows = [
        np.asarray(
            [
                [
                    upper[row, coordinate]
                    if bits & (1 << coordinate)
                    else lower[row, coordinate]
                    for coordinate in range(4)
                ]
                for bits in range(16)
            ]
        )
        for row in range(4)
    ]
    selections = np.indices((16, 16, 16, 16)).reshape(4, -1).T
    matrices = np.empty((len(selections), 4, 4), dtype=float)
    for row in range(4):
        matrices[:, row, :] = corner_rows[row][selections[:, row]]

    values = np.zeros(len(selections), dtype=float)
    absolute_sum = np.zeros(len(selections), dtype=float)
    for permutation in itertools.permutations(range(4)):
        inversions = sum(
            permutation[first] > permutation[second]
            for first in range(4)
            for second in range(first + 1, 4)
        )
        term = np.prod(matrices[:, np.arange(4), np.asarray(permutation)], axis=1)
        absolute_sum += np.abs(term)
        values += (-1.0 if inversions % 2 else 1.0) * term
    epsilon = np.finfo(float).eps
    gamma = 128.0 * epsilon / (1.0 - 128.0 * epsilon)
    guard = gamma * float(np.max(absolute_sum)) + 8.0 * np.finfo(float).tiny
    minimum = float(np.nextafter(float(np.min(values)) - guard, -math.inf))
    maximum = float(np.nextafter(float(np.max(values)) + guard, math.inf))
    return Interval(minimum, maximum)


def replacement_determinant_bounds(
    lower: np.ndarray,
    upper: np.ndarray,
    row: int,
    replacement: np.ndarray,
    sign: int,
) -> tuple[float, float]:
    """Return guarded vertex extrema of a row-replacement determinant.

    Replacing row ``z`` by a fixed vector makes the determinant multi-affine
    in the twelve coordinates of the other three rows.  A multi-affine
    polynomial reaches both extrema of a box at a vertex, so the 4096 values
    below are exhaustive rather than sampled.
    """

    if sign not in {-1, 1}:
        raise ValueError("determinant sign must be +1 or -1")
    lower = np.asarray(lower, dtype=float)
    upper = np.asarray(upper, dtype=float)
    replacement = np.asarray(replacement, dtype=float)
    if lower.shape != (4, 4) or upper.shape != (4, 4):
        raise ValueError("input bounds must have shape (4,4)")
    if replacement.shape != (4,):
        raise ValueError("replacement row must have shape (4,)")
    free_rows = [index for index in range(4) if index != row]
    corner_rows: dict[int, np.ndarray] = {}
    for free_row in free_rows:
        corner_rows[free_row] = np.asarray(
            [
                [
                    upper[free_row, coordinate]
                    if bits & (1 << coordinate)
                    else lower[free_row, coordinate]
                    for coordinate in range(4)
                ]
                for bits in range(16)
            ]
        )
    selections = np.indices((16, 16, 16)).reshape(3, -1).T
    matrices = np.empty((len(selections), 4, 4), dtype=float)
    matrices[:, row, :] = replacement
    for position, free_row in enumerate(free_rows):
        matrices[:, free_row, :] = corner_rows[free_row][selections[:, position]]

    values = np.zeros(len(selections), dtype=float)
    absolute_sum = np.zeros(len(selections), dtype=float)
    for permutation in itertools.permutations(range(4)):
        inversions = sum(
            permutation[first] > permutation[second]
            for first in range(4)
            for second in range(first + 1, 4)
        )
        term = np.prod(matrices[:, np.arange(4), np.asarray(permutation)], axis=1)
        absolute_sum += np.abs(term)
        values += (-1.0 if inversions % 2 else 1.0) * term
    values *= float(sign)
    epsilon = np.finfo(float).eps
    gamma = 128.0 * epsilon / (1.0 - 128.0 * epsilon)
    guard = gamma * float(np.max(absolute_sum)) + 8.0 * np.finfo(float).tiny
    minimum = float(np.nextafter(float(np.min(values)) - guard, -math.inf))
    maximum = float(np.nextafter(float(np.max(values)) + guard, math.inf))
    return minimum, maximum


def determinant_povm_witnesses(
    lower: np.ndarray,
    upper: np.ndarray,
    input_pauli: np.ndarray,
    statistics: np.ndarray,
    tolerance: float,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Build robust determinant-scaled POVM witnesses for one input box."""

    determinant_box = determinant_vertex_bounds(lower, upper)
    if determinant_box.lower > 0.0:
        sign = 1
    elif determinant_box.upper < 0.0:
        sign = -1
    else:
        return [], {
            "determinant_interval": [
                float(determinant_box.lower),
                float(determinant_box.upper),
            ],
            "determinant_bounds_method": "exhaustive-multiaffine-vertices",
            "sign_definite": False,
        }
    matrix = np.asarray(input_pauli, dtype=float)
    table = np.asarray(statistics, dtype=float)
    if matrix.shape != (4, 4) or table.shape != (4, 4, 3):
        raise ValueError("invalid determinant-witness point")
    determinant = float(np.linalg.det(matrix))
    if sign * determinant <= 0.0:
        raise RuntimeError("candidate determinant disagrees with interval sign")
    columns = table.reshape(4, 12)
    effects = np.linalg.solve(matrix, columns)
    candidates: list[dict[str, Any]] = []
    negative_effects: list[dict[str, Any]] = []
    for index in range(12):
        vector = effects[1:, index]
        norm = float(np.linalg.norm(vector))
        if norm <= 1e-15:
            if effects[0, index] >= 0.0:
                continue
            direction = np.asarray([1.0, 1.0, 0.0, 0.0])
        else:
            direction = np.concatenate(([1.0], -vector / norm))
        effect_margin = float(direction @ effects[:, index])
        if effect_margin >= -tolerance:
            continue
        upper_coefficients = np.asarray(
            [
                replacement_determinant_bounds(lower, upper, row, direction, sign)[1]
                for row in range(4)
            ]
        )
        robust_lhs = float(upper_coefficients @ columns[:, index])
        exact_lhs = float(abs(determinant) * effect_margin)
        diagnostic = {
            "effect_index": index,
            "path_outcome": list(divmod(index, 3)),
            "direction": direction.tolist(),
            "effect_margin": effect_margin,
            "determinant_scaled_margin": exact_lhs,
            "upper_coefficients": upper_coefficients.tolist(),
            "robust_lhs": robust_lhs,
            "enclosure_gap": robust_lhs - exact_lhs,
            "relative_enclosure_gap": (
                (robust_lhs - exact_lhs) / max(-exact_lhs, tolerance)
            ),
            "robust": bool(robust_lhs < -tolerance),
        }
        negative_effects.append(diagnostic)
        if robust_lhs >= -tolerance:
            continue
        y, t = divmod(index, 3)
        coefficients = np.zeros((4, 4, 3), dtype=float)
        coefficients[:, y, t] = -upper_coefficients
        candidates.append(
            {
                "kind": "positive-effect",
                "coefficients": coefficients.tolist(),
                "bound": 0.0,
                "effect_index": index,
                "path_outcome": [y, t],
                "direction": direction.tolist(),
                "effect_margin": effect_margin,
                "determinant_scaled_margin": exact_lhs,
                "upper_coefficients": upper_coefficients.tolist(),
                "robust_lhs": robust_lhs,
                "violation": -robust_lhs,
                "determinant_sign": sign,
            }
        )
    candidates.sort(key=lambda item: float(item["violation"]), reverse=True)
    negative_effects.sort(
        key=lambda item: (
            float(item["robust_lhs"]),
            float(item["relative_enclosure_gap"]),
        )
    )
    return candidates, {
        "determinant_interval": [
            float(determinant_box.lower),
            float(determinant_box.upper),
        ],
        "determinant_bounds_method": "exhaustive-multiaffine-vertices",
        "candidate_determinant": determinant,
        "sign_definite": True,
        "negative_effect_count": int(
            np.sum(effects[0] - np.linalg.norm(effects[1:], axis=0) < 0.0)
        ),
        "robust_witness_count": len(candidates),
        "minimum_robust_lhs": (
            float(negative_effects[0]["robust_lhs"]) if negative_effects else None
        ),
        "negative_effects": negative_effects,
    }


def planar_ando_direction(
    pulled_effects: np.ndarray,
    reconstruction: np.ndarray,
    grid_size: int = 512,
) -> dict[str, Any]:
    """Return the strongest rank-one planar Ando violation found.

    The columns of pulled_effects are Pauli coordinates of the three pulled
    terminal effects for one instrument outcome. Complete positivity of a
    common qubit instrument requires every positive test functional and
    every planar phase to satisfy the corresponding Ando inequality. For a
    fixed phase, maximisation over the Bloch vector is analytic.

    This numerical search only chooses a separating direction. Rigorous
    validity comes from the interval-enclosed cut built by
    ando_witness_coefficient_bounds.
    """

    effects = np.asarray(pulled_effects, dtype=float)
    reconstruction = np.asarray(reconstruction, dtype=float)
    if effects.shape != (4, 3):
        raise ValueError("pulled effects must have shape (4,3)")
    if reconstruction.shape != (2, 3):
        raise ValueError("reconstruction must have shape (2,3)")
    if grid_size < 16:
        raise ValueError("Ando phase grid must contain at least 16 points")

    total = np.sum(effects, axis=1)
    visible_x = effects @ reconstruction[0]
    visible_y = effects @ reconstruction[1]

    def exposed_violation(phase: float) -> float:
        difference = (
            math.cos(phase) * visible_x
            + math.sin(phase) * visible_y
            - total
        )
        return float(difference[0] + np.linalg.norm(difference[1:]))

    phases = np.linspace(0.0, 2.0 * math.pi, grid_size, endpoint=False)
    values = np.asarray([exposed_violation(phase) for phase in phases])
    index = int(np.argmax(values))
    step = 2.0 * math.pi / grid_size
    centre = float(phases[index])
    refinement = minimize_scalar(
        lambda phase: -exposed_violation(float(phase)),
        bounds=(centre - step, centre + step),
        method="bounded",
        options={"xatol": 1e-14},
    )
    phase = float(refinement.x) % (2.0 * math.pi)
    difference = math.cos(phase) * visible_x + math.sin(phase) * visible_y - total
    vector_norm = float(np.linalg.norm(difference[1:]))
    direction = (
        difference[1:] / vector_norm
        if vector_norm > 1e-15
        else np.zeros(3, dtype=float)
    )
    test = np.concatenate(([1.0], direction))
    total_expectation = float(test @ total)
    visible_expectations = np.asarray(
        [float(test @ visible_x), float(test @ visible_y)]
    )
    exact_margin = float(
        total_expectation
        - math.cos(phase) * visible_expectations[0]
        - math.sin(phase) * visible_expectations[1]
    )
    return {
        "violation": -exact_margin,
        "phase": phase,
        "test_direction": test.tolist(),
        "total_expectation": total_expectation,
        "visible_expectations": visible_expectations.tolist(),
        "exact_margin": exact_margin,
    }


def ando_witness_coefficient_bounds(
    lower: np.ndarray,
    upper: np.ndarray,
    alpha_bounds: tuple[float, float],
    beta_bounds: tuple[float, float],
    test_direction: np.ndarray,
    phase: float,
    sign: int,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    """Enclose every coefficient of a determinant-scaled Ando witness.

    Cramer's rule expresses each tested pulled effect as signed
    row-replacement determinants times observed probabilities. The planar
    complete-positivity inequality couples all three terminal outcomes.
    Interval enclosing both factors over the input and terminal boxes gives
    a robust linear cut because every probability is nonnegative.
    """

    test = np.asarray(test_direction, dtype=float)
    if test.shape != (4,):
        raise ValueError("Ando test direction must have shape (4,)")
    if sign not in {-1, 1}:
        raise ValueError("determinant sign must be +1 or -1")
    replacement = [
        Interval(
            *replacement_determinant_bounds(
                lower,
                upper,
                row,
                test,
                sign,
            )
        )
        for row in range(4)
    ]
    x_intervals, y_intervals = reconstruction_intervals(
        alpha_bounds,
        beta_bounds,
    )
    cosine = Interval.point(math.cos(float(phase)))
    sine = Interval.point(math.sin(float(phase)))
    terminal = [
        Interval.point(1.0) - cosine * x_intervals[t] - sine * y_intervals[t]
        for t in range(3)
    ]
    product = [
        [replacement[z] * terminal[t] for t in range(3)]
        for z in range(4)
    ]
    coefficient_lower = np.asarray(
        [[product[z][t].lower for t in range(3)] for z in range(4)]
    )
    coefficient_upper = np.asarray(
        [[product[z][t].upper for t in range(3)] for z in range(4)]
    )
    return coefficient_lower, coefficient_upper, {
        "replacement_determinant_intervals": [
            [item.lower, item.upper] for item in replacement
        ],
        "terminal_coefficient_intervals": [
            [item.lower, item.upper] for item in terminal
        ],
        "coefficient_intervals": [
            [
                [product[z][t].lower, product[z][t].upper]
                for t in range(3)
            ]
            for z in range(4)
        ],
    }


def determinant_ando_witnesses(
    lower: np.ndarray,
    upper: np.ndarray,
    alpha_bounds: tuple[float, float],
    beta_bounds: tuple[float, float],
    input_pauli: np.ndarray,
    statistics: np.ndarray,
    tolerance: float,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Build robust witnesses for one common planar CP instrument.

    Unlike separate-effect POVM witnesses, each Ando cut couples all three
    terminal outcomes for a fixed path outcome. A negative robust upper
    envelope excludes every common Choi matrix throughout the complete input
    and terminal parameter box.
    """

    determinant_box = determinant_vertex_bounds(lower, upper)
    if determinant_box.lower > 0.0:
        sign = 1
    elif determinant_box.upper < 0.0:
        sign = -1
    else:
        return [], {
            "determinant_interval": [
                float(determinant_box.lower),
                float(determinant_box.upper),
            ],
            "determinant_bounds_method": "exhaustive-multiaffine-vertices",
            "sign_definite": False,
        }
    matrix = np.asarray(input_pauli, dtype=float)
    table = np.asarray(statistics, dtype=float)
    if matrix.shape != (4, 4) or table.shape != (4, 4, 3):
        raise ValueError("invalid Ando-witness point")
    determinant = float(np.linalg.det(matrix))
    if sign * determinant <= 0.0:
        raise RuntimeError("candidate determinant disagrees with interval sign")
    effects = np.linalg.solve(matrix, table.reshape(4, 12)).reshape(4, 4, 3)
    midpoint_reconstruction = planar_reconstruction(
        0.5 * (alpha_bounds[0] + alpha_bounds[1]),
        0.5 * (beta_bounds[0] + beta_bounds[1]),
    )
    candidates: list[dict[str, Any]] = []
    diagnostics: list[dict[str, Any]] = []
    for outcome in range(4):
        direction = planar_ando_direction(
            effects[:, outcome, :],
            midpoint_reconstruction,
        )
        if float(direction["violation"]) <= tolerance:
            continue
        test = np.asarray(direction["test_direction"], dtype=float)
        phase = float(direction["phase"])
        coefficient_lower, coefficient_upper, interval_audit = (
            ando_witness_coefficient_bounds(
                lower,
                upper,
                alpha_bounds,
                beta_bounds,
                test,
                phase,
                sign,
            )
        )
        robust_lhs = float(np.sum(coefficient_upper * table[:, outcome, :]))
        exact_lhs = float(abs(determinant) * float(direction["exact_margin"]))
        diagnostic = {
            "path_outcome": outcome,
            **direction,
            "determinant_scaled_margin": exact_lhs,
            "robust_lhs": robust_lhs,
            "enclosure_gap": robust_lhs - exact_lhs,
            "robust": bool(robust_lhs < -tolerance),
            **interval_audit,
        }
        diagnostics.append(diagnostic)
        if robust_lhs >= -tolerance:
            continue
        coefficients = np.zeros((4, 4, 3), dtype=float)
        coefficients[:, outcome, :] = -coefficient_upper
        candidates.append(
            {
                "kind": "planar-ando",
                "coefficients": coefficients.tolist(),
                "bound": 0.0,
                "path_outcome": outcome,
                "phase": phase,
                "test_direction": test.tolist(),
                "exact_margin": float(direction["exact_margin"]),
                "determinant_scaled_margin": exact_lhs,
                "upper_coefficients": coefficient_upper.tolist(),
                "robust_lhs": robust_lhs,
                "violation": -robust_lhs,
                "determinant_sign": sign,
            }
        )
    candidates.sort(key=lambda item: float(item["violation"]), reverse=True)
    diagnostics.sort(key=lambda item: float(item["robust_lhs"]))
    return candidates, {
        "determinant_interval": [
            float(determinant_box.lower),
            float(determinant_box.upper),
        ],
        "determinant_bounds_method": "exhaustive-multiaffine-vertices",
        "candidate_determinant": determinant,
        "sign_definite": True,
        "violated_direction_count": len(diagnostics),
        "robust_witness_count": len(candidates),
        "minimum_robust_lhs": (
            float(diagnostics[0]["robust_lhs"]) if diagnostics else None
        ),
        "violated_directions": diagnostics,
    }


def determinant_split_scores(
    lower: np.ndarray,
    upper: np.ndarray,
    statistics: np.ndarray,
    determinant_audit: dict[str, Any] | None,
    maximum_effects: int = 2,
) -> np.ndarray:
    """Score bisections by their exact improvement of robust POVM margins.

    For each of the closest nonpositive recovered effects, every candidate
    coordinate is bisected virtually.  On each child, the row-replacement
    determinant coefficients are re-enclosed by exhaustive vertex extrema.
    The score is the relative reduction in the worse child margin.  This is
    only a branching heuristic; all reported upper bounds continue to come
    from the convex oracle and all added cuts remain independently valid.
    """

    scores = np.zeros((4, 4), dtype=float)
    if not determinant_audit or not determinant_audit.get("sign_definite"):
        return scores
    sign = 1 if float(determinant_audit["candidate_determinant"]) > 0.0 else -1
    lower = np.asarray(lower, dtype=float)
    upper = np.asarray(upper, dtype=float)
    columns = np.asarray(statistics, dtype=float).reshape(4, 12)
    effects = [
        item
        for item in determinant_audit.get("negative_effects", [])
        if not bool(item.get("robust", False))
    ]
    if not effects:
        effects = list(determinant_audit.get("negative_effects", []))
    effects = effects[: max(0, int(maximum_effects))]
    for item in effects:
        effect_index = int(item["effect_index"])
        direction = np.asarray(item["direction"], dtype=float)
        probabilities = np.maximum(columns[:, effect_index], 0.0)
        parent_lhs = float(item["robust_lhs"])
        scale = max(
            -float(item["determinant_scaled_margin"]),
            abs(parent_lhs),
            1e-12,
        )
        for row in range(4):
            for coordinate in range(4):
                midpoint = 0.5 * (lower[row, coordinate] + upper[row, coordinate])
                if not lower[row, coordinate] < midpoint < upper[row, coordinate]:
                    continue
                child_lhs: list[float] = []
                for side in range(2):
                    child_lower = lower.copy()
                    child_upper = upper.copy()
                    if side == 0:
                        child_upper[row, coordinate] = midpoint
                    else:
                        child_lower[row, coordinate] = midpoint
                    upper_coefficients = np.asarray(
                        [
                            replacement_determinant_bounds(
                                child_lower,
                                child_upper,
                                replaced_row,
                                direction,
                                sign,
                            )[1]
                            for replaced_row in range(4)
                        ]
                    )
                    child_lhs.append(float(upper_coefficients @ probabilities))
                improvement = parent_lhs - max(child_lhs)
                scores[row, coordinate] += max(0.0, improvement) / scale
    return scores


def ando_input_split_scores(
    lower: np.ndarray,
    upper: np.ndarray,
    alpha_bounds: tuple[float, float],
    beta_bounds: tuple[float, float],
    statistics: np.ndarray,
    ando_audit: dict[str, Any] | None,
    maximum_directions: int = 1,
    maximum_coordinates: int = ANDO_SPLIT_SHORTLIST,
) -> np.ndarray:
    """Score shortlisted input bisections by robust Ando improvement.

    Cofactor sensitivity at the box midpoint cheaply ranks all sixteen input
    coordinates. The strongest candidates are then bisected virtually and
    their rigorous coefficient enclosures recomputed on both children. This
    is only a branching heuristic: every later cut remains independently
    certified on the complete child box where it is applied.
    """

    scores = np.zeros((4, 4), dtype=float)
    if not ando_audit or not ando_audit.get("sign_definite"):
        return scores
    sign = 1 if float(ando_audit["candidate_determinant"]) > 0.0 else -1
    table = np.asarray(statistics, dtype=float)
    midpoint_matrix = 0.5 * (np.asarray(lower) + np.asarray(upper))
    reconstruction = planar_reconstruction(
        0.5 * (alpha_bounds[0] + alpha_bounds[1]),
        0.5 * (beta_bounds[0] + beta_bounds[1]),
    )
    directions = [
        item
        for item in ando_audit.get("violated_directions", [])
        if not bool(item.get("robust", False))
    ][: max(0, int(maximum_directions))]
    for item in directions:
        outcome = int(item["path_outcome"])
        test = np.asarray(item["test_direction"], dtype=float)
        phase = float(item["phase"])
        probabilities = np.maximum(table[:, outcome, :], 0.0)
        terminal_coefficients = (
            1.0
            - math.cos(phase) * reconstruction[0]
            - math.sin(phase) * reconstruction[1]
        )
        sensitivity = np.zeros((4, 4), dtype=float)
        for row in range(4):
            for coordinate in range(4):
                half_width = 0.5 * float(upper[row, coordinate] - lower[row, coordinate])
                if half_width <= 0.0:
                    continue
                for replaced_row in range(4):
                    if replaced_row == row:
                        continue
                    matrix = midpoint_matrix.copy()
                    matrix[replaced_row] = test
                    minor = np.delete(
                        np.delete(matrix, row, axis=0),
                        coordinate,
                        axis=1,
                    )
                    derivative = (
                        sign
                        * (-1.0 if (row + coordinate) % 2 else 1.0)
                        * float(np.linalg.det(minor))
                    )
                    sensitivity[row, coordinate] += half_width * np.sum(
                        np.abs(derivative * terminal_coefficients)
                        * probabilities[replaced_row]
                    )
        ranked = np.argsort(sensitivity.ravel())[::-1]
        selected = [
            int(index)
            for index in ranked[: max(0, int(maximum_coordinates))]
            if sensitivity.ravel()[index] > 0.0
        ]
        parent_lhs = float(item["robust_lhs"])
        scale = max(
            -float(item["determinant_scaled_margin"]),
            abs(parent_lhs),
            1e-12,
        )
        for flat_index in selected:
            row, coordinate = np.unravel_index(flat_index, (4, 4))
            row = int(row)
            coordinate = int(coordinate)
            midpoint = 0.5 * (lower[row, coordinate] + upper[row, coordinate])
            child_lhs: list[float] = []
            for side in range(2):
                child_lower = lower.copy()
                child_upper = upper.copy()
                if side == 0:
                    child_upper[row, coordinate] = midpoint
                else:
                    child_lower[row, coordinate] = midpoint
                _, coefficient_upper, _ = ando_witness_coefficient_bounds(
                    child_lower,
                    child_upper,
                    alpha_bounds,
                    beta_bounds,
                    test,
                    phase,
                    sign,
                )
                child_lhs.append(float(np.sum(coefficient_upper * probabilities)))
            improvement = parent_lhs - max(child_lhs)
            scores[row, coordinate] += max(0.0, improvement) / scale
    return scores


def product_residual_scores(oracle: TernaryConeOracle) -> np.ndarray:
    """Measure where the current McCormick solution departs from a product.

    The score is only a branching heuristic.  It does not enter the bound or
    the certificate: a large value identifies an input coordinate whose
    lifted products are being used most nonphysically by the current optimum.
    """

    if oracle.input_pauli is None or not oracle.common_instrument_products:
        return np.zeros((4, 4), dtype=float)
    inputs = np.asarray(oracle.input_pauli.value, dtype=float)
    scores = np.zeros((4, 4), dtype=float)
    for z in range(4):
        for mu in range(4):
            coordinate = float(inputs[z, mu])
            for y in range(4):
                choi_real, choi_imaginary = oracle.common_instrument_choi[y]
                product_real, product_imaginary = oracle.common_instrument_products[z][
                    mu
                ][y]
                scores[z, mu] += float(
                    np.sum(
                        np.abs(
                            np.asarray(product_real.value)
                            - coordinate * np.asarray(choi_real.value)
                        )
                    )
                    + np.sum(
                        np.abs(
                            np.asarray(product_imaginary.value)
                            - coordinate * np.asarray(choi_imaginary.value)
                        )
                    )
                )
            if oracle.effective_povm is not None and oracle.common_povm_products:
                for index in range(12):
                    effect = np.asarray(oracle.effective_povm.value[index])
                    product = np.asarray(oracle.common_povm_products[z][index].value)
                    scores[z, mu] += abs(float(product[mu] - coordinate * effect[mu]))
    return scores


def _assemble(
    source: dict[str, Any],
    localisation: dict[str, Any],
    terminal_audit: dict[str, Any],
    records: list[dict[str, Any]],
    pending_heap: list[tuple[float, int, dict[str, Any]]],
    unresolved: list[dict[str, Any]],
    next_identifier: int,
    top_solution: dict[str, Any] | None,
    max_nodes: int,
    bound_safety: float,
    minimum_width: float,
    use_top_spectral_cell: bool,
    max_witnesses: int,
    max_new_witnesses_per_node: int,
    witness_tolerance: float,
    use_ando_witnesses: bool,
) -> dict[str, Any]:
    pending = [item[2] for item in sorted(pending_heap)]
    closed = sum(record.get("disposition") == "closed" for record in records)
    split = sum(record.get("disposition") == "split" for record in records)
    frontier_bound = max(
        (float(node["parent_bound"]) for node in pending),
        default=-math.inf,
    )
    statuses_complete = all(
        record.get("status")
        in {
            cp.OPTIMAL,
            cp.OPTIMAL_INACCURATE,
            cp.INFEASIBLE,
            cp.INFEASIBLE_INACCURATE,
        }
        for record in records
    )
    complete = not pending and not unresolved and statuses_complete
    return {
        "support_weight": source["support_weight"],
        "target": source["target"],
        "source": (
            "ternary spatial common-instrument and common-POVM McCormick "
            "cover with determinant-scaled POVM and Ando witnesses"
        ),
        "source_box": source["box"],
        "base_code": source["base_code"],
        "base_plane": source.get("base_plane"),
        "base_sphere": source.get("base_sphere"),
        "top_spectral_cell": bool(use_top_spectral_cell),
        "localisation": localisation,
        "terminal_effect_enclosure": terminal_audit,
        "bound_safety": float(bound_safety),
        "minimum_width": float(minimum_width),
        "max_witnesses": int(max_witnesses),
        "max_new_witnesses_per_node": int(max_new_witnesses_per_node),
        "witness_tolerance": float(witness_tolerance),
        "planar_ando_witnesses": bool(use_ando_witnesses),
        "determinant_bounds_method": "exhaustive-multiaffine-vertices",
        "determinant_near_relative_gap": DETERMINANT_NEAR_RELATIVE_GAP,
        "ando_near_relative_gap": ANDO_NEAR_RELATIVE_GAP,
        "ando_split_shortlist": ANDO_SPLIT_SHORTLIST,
        "maximum_determinant_branch_streak": MAX_DETERMINANT_BRANCH_STREAK,
        "max_nodes": int(max_nodes),
        "solved_nodes": len(records),
        "closed_nodes": int(closed),
        "split_nodes": int(split),
        "pending_nodes": len(pending),
        "unresolved_nodes": len(unresolved),
        "maximum_depth": max((int(record["depth"]) for record in records), default=0),
        "determinant_witness_count": int(
            sum(int(record.get("new_witnesses", 0)) for record in records)
        ),
        "planar_ando_witness_count": int(
            sum(int(record.get("new_ando_witnesses", 0)) for record in records)
        ),
        "maximum_pending_bound": float(frontier_bound),
        "statuses_complete": bool(statuses_complete),
        "complete": bool(complete),
        "next_identifier": int(next_identifier),
        "pending": pending,
        "unresolved": unresolved,
        "records": records,
        "top_solution": top_solution,
        "scope": (
            "one continuous terminal cell and the selected Fourier spectral "
            "cell; all input--Choi and input--POVM products use convergent "
            "McCormick envelopes; pure-prefix caps, determinant-scaled "
            "common-POVM cuts, and optional planar Ando cuts are valid on "
            "complete input and terminal boxes; numerical SDP bounds remain "
            "solver-conditional"
        ),
    }


def cover_candidate_region(
    source: dict[str, Any],
    localisation: dict[str, Any],
    output: Path,
    max_nodes: int,
    bound_safety: float,
    minimum_width: float,
    checkpoint_every: int,
    use_top_spectral_cell: bool,
    max_witnesses: int,
    max_new_witnesses_per_node: int,
    witness_tolerance: float,
    use_ando_witnesses: bool,
    resume: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Cover one localised input region with spatial instrument cells."""

    box, contractions, reconstruction = _configuration(source, use_top_spectral_cell)
    terminal_effects, terminal_errors, _, terminal_audit = (
        terminal_effect_anchor_and_errors(box["terminal_alpha"], box["terminal_beta"])
    )

    lower_parameter = cp.Parameter((4, 4))
    upper_parameter = cp.Parameter((4, 4))
    purity_caps = cp.Parameter((4, 4))
    oracle = TernaryConeOracle(
        **_oracle_keywords(source),
        common_contractions=contractions,
        terminal_reconstruction=reconstruction,
        input_pauli_lower=lower_parameter,
        input_pauli_upper=upper_parameter,
        input_purity_caps=purity_caps,
        common_povm_bilinear=True,
        common_instrument_terminal_effect_anchor=terminal_effects,
        common_instrument_terminal_effect_errors=terminal_errors,
        max_common_instrument_witnesses=max_witnesses,
    )

    target = float(source["target"])
    if resume is None:
        lower = np.asarray(localisation["lower"], dtype=float)
        upper = np.asarray(localisation["upper"], dtype=float)
        root = _node_payload(0, None, 0, lower, upper)
        root["parent_bound"] = math.inf
        root["determinant_witnesses"] = []
        root["determinant_branch_streak"] = 0
        pending: list[tuple[float, int, dict[str, Any]]] = [(-math.inf, 0, root)]
        records: list[dict[str, Any]] = []
        unresolved: list[dict[str, Any]] = []
        next_identifier = 1
        top_solution: dict[str, Any] | None = None
    else:
        pending = []
        for node in resume["pending"]:
            heapq.heappush(
                pending,
                (-float(node["parent_bound"]), int(node["identifier"]), node),
            )
        records = list(resume["records"])
        unresolved = list(resume["unresolved"])
        next_identifier = int(resume["next_identifier"])
        top_solution = resume.get("top_solution")

    while pending and len(records) < max_nodes:
        _, _, node = heapq.heappop(pending)
        lower = np.asarray(node["lower"], dtype=float)
        upper = np.asarray(node["upper"], dtype=float)
        lower_parameter.value = lower
        upper_parameter.value = upper
        purity_caps.value = box_purity_caps(lower, upper)
        witness_records = list(node.get("determinant_witnesses", []))
        new_witnesses = 0
        new_ando_witnesses = 0
        oracle_solves = 0
        determinant_audit: dict[str, Any] | None = None
        ando_audit: dict[str, Any] | None = None
        witness_audit_solve = 0
        while True:
            active_witnesses = tuple(
                (
                    np.asarray(item["coefficients"], dtype=float),
                    float(item["bound"]),
                )
                for item in witness_records
            )
            result = oracle.solve(
                box,
                bound_safety,
                capture=True,
                common_instrument_witnesses=active_witnesses,
            )
            oracle_solves += 1
            if result["status"] not in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}:
                break
            if float(result["bound"]) < target:
                break
            if (
                new_witnesses >= max_new_witnesses_per_node
                or len(witness_records) >= max_witnesses
            ):
                break
            input_pauli = np.column_stack(
                [result["prefix"], result["input_bloch_vectors"]]
            )
            statistics = np.asarray(result["statistics"], dtype=float)
            candidate_witnesses, determinant_audit = determinant_povm_witnesses(
                lower,
                upper,
                input_pauli,
                statistics,
                witness_tolerance,
            )
            if use_ando_witnesses:
                ando_candidates, ando_audit = determinant_ando_witnesses(
                    lower,
                    upper,
                    tuple(map(float, box["terminal_alpha"])),
                    tuple(map(float, box["terminal_beta"])),
                    input_pauli,
                    statistics,
                    witness_tolerance,
                )
                candidate_witnesses.extend(ando_candidates)
            candidate_witnesses.sort(
                key=lambda item: float(item["violation"]),
                reverse=True,
            )
            witness_audit_solve = oracle_solves
            existing = {
                tuple(
                    int(round(value * 1e12))
                    for value in np.asarray(item["coefficients"]).ravel()
                )
                for item in witness_records
            }
            novel = [
                item
                for item in candidate_witnesses
                if tuple(
                    int(round(value * 1e12))
                    for value in np.asarray(item["coefficients"]).ravel()
                )
                not in existing
            ]
            if not novel:
                break
            capacity = min(
                max_new_witnesses_per_node - new_witnesses,
                max_witnesses - len(witness_records),
            )
            additions = novel[:capacity]
            witness_records.extend(additions)
            new_witnesses += len(additions)
            new_ando_witnesses += sum(
                item.get("kind") == "planar-ando" for item in additions
            )
        bound = float(result["bound"])
        widths = upper - lower
        if (
            result["status"] in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}
            and bound >= target
            and witness_audit_solve != oracle_solves
        ):
            input_pauli = np.column_stack(
                [result["prefix"], result["input_bloch_vectors"]]
            )
            statistics = np.asarray(result["statistics"], dtype=float)
            _, determinant_audit = determinant_povm_witnesses(
                lower,
                upper,
                input_pauli,
                statistics,
                witness_tolerance,
            )
            if use_ando_witnesses:
                _, ando_audit = determinant_ando_witnesses(
                    lower,
                    upper,
                    tuple(map(float, box["terminal_alpha"])),
                    tuple(map(float, box["terminal_beta"])),
                    input_pauli,
                    statistics,
                    witness_tolerance,
                )
        residual_scores = (
            product_residual_scores(oracle)
            if result["status"] in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}
            else np.zeros((4, 4), dtype=float)
        )
        negative_effects = (
            determinant_audit.get("negative_effects", []) if determinant_audit else []
        )
        near_determinant_margin = any(
            float(item["relative_enclosure_gap"]) <= DETERMINANT_NEAR_RELATIVE_GAP
            for item in negative_effects
        )
        ando_directions = (
            ando_audit.get("violated_directions", []) if ando_audit else []
        )
        near_ando_margin = any(
            float(item["enclosure_gap"])
            / max(-float(item["determinant_scaled_margin"]), 1e-15)
            <= ANDO_NEAR_RELATIVE_GAP
            for item in ando_directions
        )
        determinant_branch_streak = int(node.get("determinant_branch_streak", 0))
        determinant_scores = (
            determinant_split_scores(
                lower,
                upper,
                np.asarray(result["statistics"], dtype=float),
                determinant_audit,
            )
            if result["status"] in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}
            and bound >= target
            and determinant_branch_streak < MAX_DETERMINANT_BRANCH_STREAK
            and (new_witnesses > 0 or near_determinant_margin)
            else np.zeros((4, 4), dtype=float)
        )
        ando_scores = (
            ando_input_split_scores(
                lower,
                upper,
                tuple(map(float, box["terminal_alpha"])),
                tuple(map(float, box["terminal_beta"])),
                np.asarray(result["statistics"], dtype=float),
                ando_audit,
            )
            if use_ando_witnesses
            and result["status"] in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}
            and bound >= target
            and determinant_branch_streak < MAX_DETERMINANT_BRANCH_STREAK
            and (new_ando_witnesses > 0 or near_ando_margin)
            else np.zeros((4, 4), dtype=float)
        )
        record = {
            **node,
            "determinant_witnesses": witness_records,
            **_compact_result(result),
            "maximum_coordinate_width": float(np.max(widths)),
            "maximum_row_l1_width": float(np.max(np.sum(widths, axis=1))),
            "maximum_product_residual": float(np.max(residual_scores)),
            "maximum_determinant_split_score": float(np.max(determinant_scores)),
            "maximum_ando_split_score": float(np.max(ando_scores)),
            "near_determinant_margin": bool(near_determinant_margin),
            "near_ando_margin": bool(near_ando_margin),
            "determinant_branch_streak": determinant_branch_streak,
            "oracle_solves": oracle_solves,
            "new_witnesses": new_witnesses,
            "new_ando_witnesses": new_ando_witnesses,
            "determinant_audit": determinant_audit,
            "ando_audit": ando_audit,
        }
        records.append(record)
        if top_solution is None or bound > float(top_solution["bound"]):
            top_solution = {
                **record,
                "prefix": result.get("prefix"),
                "input_bloch_vectors": result.get("input_bloch_vectors"),
                "statistics": result.get("statistics"),
                "common_instrument_choi": result.get("common_instrument_choi"),
                "common_instrument_anchor_statistics": result.get(
                    "common_instrument_anchor_statistics"
                ),
            }
        if (
            result["status"] in {cp.INFEASIBLE, cp.INFEASIBLE_INACCURATE}
            or bound < target
        ):
            record["disposition"] = "closed"
        elif result["status"] not in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}:
            record["disposition"] = "solver-unresolved"
            unresolved.append(record)
        else:
            margin_scores = np.maximum(determinant_scores, ando_scores)
            margin_eligible = margin_scores.copy()
            margin_eligible[widths <= minimum_width] = -math.inf
            residual_eligible = residual_scores.copy()
            residual_eligible[widths <= minimum_width] = -math.inf
            use_margin_branch = (
                np.any(np.isfinite(margin_eligible))
                and float(np.max(margin_eligible)) > 1e-12
                and determinant_branch_streak < MAX_DETERMINANT_BRANCH_STREAK
                and (
                    new_witnesses > 0
                    or near_determinant_margin
                    or near_ando_margin
                )
            )
            if use_margin_branch:
                row, coordinate = np.unravel_index(
                    int(np.argmax(margin_eligible)),
                    margin_eligible.shape,
                )
                row = int(row)
                coordinate = int(coordinate)
                record["branching_rule"] = (
                    "ando-margin"
                    if ando_scores[row, coordinate]
                    > determinant_scores[row, coordinate]
                    else "determinant-margin"
                )
            elif (
                np.any(np.isfinite(residual_eligible))
                and float(np.max(residual_eligible)) > 1e-10
            ):
                row, coordinate = np.unravel_index(
                    int(np.argmax(residual_eligible)), residual_eligible.shape
                )
                row = int(row)
                coordinate = int(coordinate)
                record["branching_rule"] = "product-residual"
            else:
                row, coordinate = _split_coordinate(lower, upper)
                record["branching_rule"] = "widest-coordinate"
            width = float(widths[row, coordinate])
            if width <= minimum_width:
                record["disposition"] = "resolution-limit"
                unresolved.append(record)
            else:
                midpoint = 0.5 * (lower[row, coordinate] + upper[row, coordinate])
                record["disposition"] = "split"
                record["split_coordinate"] = [row, coordinate]
                record["split_value"] = float(midpoint)
                for side in range(2):
                    child_lower = lower.copy()
                    child_upper = upper.copy()
                    if side == 0:
                        child_upper[row, coordinate] = midpoint
                    else:
                        child_lower[row, coordinate] = midpoint
                    child = _node_payload(
                        next_identifier,
                        int(node["identifier"]),
                        int(node["depth"]) + 1,
                        child_lower,
                        child_upper,
                    )
                    child["parent_bound"] = bound
                    child["determinant_witnesses"] = witness_records
                    child["determinant_branch_streak"] = (
                        determinant_branch_streak + 1
                        if record["branching_rule"]
                        in {"determinant-margin", "ando-margin"}
                        else 0
                    )
                    heapq.heappush(pending, (-bound, next_identifier, child))
                    next_identifier += 1
        print(
            json.dumps(
                {
                    "solved": len(records),
                    "pending": len(pending),
                    "closed": sum(
                        item.get("disposition") == "closed" for item in records
                    ),
                    "identifier": node["identifier"],
                    "depth": node["depth"],
                    "bound": bound,
                    "maximum_coordinate_width": record["maximum_coordinate_width"],
                    "new_witnesses": new_witnesses,
                    "new_ando_witnesses": new_ando_witnesses,
                    "disposition": record["disposition"],
                }
            ),
            flush=True,
        )
        if checkpoint_every > 0 and len(records) % checkpoint_every == 0:
            _write_checkpoint(
                output,
                _assemble(
                    source,
                    localisation,
                    terminal_audit,
                    records,
                    pending,
                    unresolved,
                    next_identifier,
                    top_solution,
                    max_nodes,
                    bound_safety,
                    minimum_width,
                    use_top_spectral_cell,
                    max_witnesses,
                    max_new_witnesses_per_node,
                    witness_tolerance,
                    use_ando_witnesses,
                ),
            )

    payload = _assemble(
        source,
        localisation,
        terminal_audit,
        records,
        pending,
        unresolved,
        next_identifier,
        top_solution,
        max_nodes,
        bound_safety,
        minimum_width,
        use_top_spectral_cell,
        max_witnesses,
        max_new_witnesses_per_node,
        witness_tolerance,
        use_ando_witnesses,
    )
    _write_checkpoint(output, payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--frontier-json", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--localisation-json", type=Path)
    parser.add_argument("--max-nodes", type=int, default=1000)
    parser.add_argument("--coordinate-safety", type=float, default=2e-6)
    parser.add_argument("--bound-safety", type=float, default=2e-6)
    parser.add_argument("--minimum-width", type=float, default=1e-6)
    parser.add_argument("--checkpoint-every", type=int, default=10)
    parser.add_argument("--max-witnesses", type=int, default=24)
    parser.add_argument("--max-new-witnesses-per-node", type=int, default=4)
    parser.add_argument("--witness-tolerance", type=float, default=2e-9)
    parser.add_argument("--planar-ando-witnesses", action="store_true")
    parser.add_argument("--top-spectral-cell", action="store_true")
    args = parser.parse_args()

    source = json.loads(args.frontier_json.read_text(encoding="utf-8"))
    resume: dict[str, Any] | None = None
    if args.resume:
        resume = json.loads(args.output.read_text(encoding="utf-8"))
        localisation = resume["localisation"]
    elif args.localisation_json is not None:
        localised = json.loads(args.localisation_json.read_text(encoding="utf-8"))
        localisation = localised.get("localisation", localised)
    else:
        localisation = localise_candidate_region(
            source, args.coordinate_safety, args.top_spectral_cell
        )
    payload = cover_candidate_region(
        source,
        localisation,
        args.output,
        args.max_nodes,
        args.bound_safety,
        args.minimum_width,
        args.checkpoint_every,
        args.top_spectral_cell,
        args.max_witnesses,
        args.max_new_witnesses_per_node,
        args.witness_tolerance,
        args.planar_ando_witnesses,
        resume,
    )
    print(
        json.dumps(
            {
                key: payload[key]
                for key in (
                    "solved_nodes",
                    "closed_nodes",
                    "split_nodes",
                    "pending_nodes",
                    "unresolved_nodes",
                    "maximum_pending_bound",
                    "statuses_complete",
                    "complete",
                )
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

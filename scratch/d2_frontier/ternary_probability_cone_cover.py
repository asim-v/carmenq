"""SOCP box cover for a continuum of ternary terminal qubit POVMs.

A rank-one ternary qubit POVM has an inellipse as its normalized probability
range.  Choose the smallest effect weight as the residual coordinate and use
Horwitz reciprocal parameters ``A,B``.  The effect weights are

    w0 = A/(A+B-1),  w1 = B/(A+B-1),
    w2 = (A+B-2)/(A+B-1).

For a terminal syndrome state, let ``q_s`` be its three measurement
probabilities.  Helstrom optimality of this same POVM can be written without
reconstructing Bloch vectors: there is a probability vector ``h`` for the
dual operator such that ``q_s``, ``h``, and ``h-q_s`` all lie in the
homogenized inellipse cone and ``sum_s q_s[s] = sum_t h[t]``.  The normal
Bloch coordinate may be chosen zero, so these cone conditions are equivalent
to positivity of the states, dual operator, and dual slacks.

On a terminal-parameter box, intervals for the three weights give a convergent
outer relaxation of ``q=w*u``.  Selected normalized pulled-effect pairs obey
the independent clean-POVM inellipse constraints.  Every node is therefore
an SOCP outer bound.  Numerical conic bounds still require outward-rounded
dual validation before they are proof-grade.
"""

from __future__ import annotations

import argparse
import heapq
import json
import math
from pathlib import Path
from typing import Any

import cvxpy as cp
import numpy as np

from pairwise_inellipse_box_cover import (
    ANCHOR_LOCATIONS,
    Box,
    box_anchor_relaxations,
    coefficientwise_box_soc_data,
    center_box_anchor_relaxations,
    center_coefficientwise_box_soc_data,
    center_inellipse_soc_data,
    deserialise_box,
    inellipse_soc_data,
    serialise_box,
    split_box,
    write_payload,
)
from pairwise_qubit_helstrom_scip import Column, parse_pair, render_column


OUTCOMES = range(4)
PATHS = tuple((z, y) for z in OUTCOMES for y in OUTCOMES)


def hellinger_hypograph(
    probabilities: list[cp.Expression], constraints: list[cp.Constraint]
) -> cp.Expression:
    """Conic hypograph of ``(sum_i sqrt(p_i))**2 / 16``."""

    cross_terms: list[cp.Expression] = []
    for first in range(len(probabilities)):
        constraints.append(probabilities[first] >= 0.0)
        for second in range(first + 1, len(probabilities)):
            geometric = cp.Variable(nonneg=True)
            constraints.append(
                cp.SOC(
                    probabilities[first] + probabilities[second],
                    cp.hstack(
                        (
                            2.0 * geometric,
                            probabilities[first] - probabilities[second],
                        )
                    ),
                )
            )
            cross_terms.append(geometric)
    return (
        cp.sum(cp.hstack(probabilities))
        + 2.0 * cp.sum(cp.hstack(cross_terms))
    ) / 16.0


TERMINAL_ALPHA = "terminal_alpha"
TERMINAL_BETA = "terminal_beta"
_ORACLE_CACHE: dict[tuple[Any, ...], "TernaryConeOracle"] = {}


def terminal_beta_cap(maximum_weight_floor: float) -> float:
    ratio = (1.0 - maximum_weight_floor) / maximum_weight_floor
    return 1.0 + 2.0 * ratio


def terminal_weights(alpha: float, beta: float) -> tuple[float, float, float]:
    """Return the rank-one effect traces encoded by Horwitz parameters."""

    denominator = alpha + beta - 1.0
    if denominator <= 0.0:
        raise ValueError("Horwitz parameters require alpha+beta>1")
    return (
        alpha / denominator,
        beta / denominator,
        (alpha + beta - 2.0) / denominator,
    )


def projective_comparison_bonus(
    weights: tuple[float, float, float] | np.ndarray,
    priors: tuple[float, float, float] | np.ndarray,
    retained: int,
    complement: int,
) -> float:
    """Upper bound on the ternary audit gain over an aligned projective readout.

    ``retained`` is answered by its rank-one projector, ``complement`` by the
    complementary projector, and the third active ternary label is deleted.
    """

    deleted = 3 - retained - complement
    if retained == complement or min(retained, complement, deleted) < 0:
        raise ValueError("retained and complement must be distinct labels in {0,1,2}")
    value = np.asarray(weights, dtype=float)
    probability = np.asarray(priors, dtype=float)
    return float(
        (1.0 - value[retained]) * probability[complement]
        + value[deleted] * probability[deleted]
    )


def initial_box(
    pair_count: int,
    maximum_weight_floor: float,
    include_priors: bool = False,
) -> Box:
    result: Box = {
        TERMINAL_ALPHA: (1.0, 2.0),
        TERMINAL_BETA: (1.0, terminal_beta_cap(maximum_weight_floor)),
    }
    result.update(
        {
            f"k{pair_index}_{coordinate}": (1.0, 2.0)
            for pair_index in range(pair_count)
            for coordinate in ("alpha", "beta")
        }
    )
    if include_priors:
        result.update(
            {
                "prior_0": (0.25, 1.0),
                "prior_1": (0.0, 0.5),
                "prior_2": (0.0, 1.0 / 3.0),
                "prior_3": (0.0, 0.25),
            }
        )
    return result


def terminal_domain_intersects(
    box: Box, maximum_weight_floor: float
) -> bool:
    al, au = box[TERMINAL_ALPHA]
    bl, bu = box[TERMINAL_BETA]
    if au < bl - 1e-15:
        return False
    # w0=A/(A+B-1) increases with A and decreases with B.
    maximum_w0 = au / (au + bl - 1.0)
    return maximum_w0 >= maximum_weight_floor - 1e-15


def terminal_weight_intervals(box: Box) -> tuple[tuple[float, float], ...]:
    al, au = box[TERMINAL_ALPHA]
    bl, bu = box[TERMINAL_BETA]
    # The monotonicities are exact on A,B >= 1.
    w0 = (al / (al + bu - 1.0), au / (au + bl - 1.0))
    w1 = (bl / (au + bl - 1.0), bu / (al + bu - 1.0))
    w2 = (
        1.0 - 1.0 / (al + bl - 1.0),
        1.0 - 1.0 / (au + bu - 1.0),
    )
    return tuple(
        (max(0.0, lower), min(1.0, upper)) for lower, upper in (w0, w1, w2)
    )


def selected_normalised(
    projective_statistics: cp.Expression,
    probability: cp.Expression,
    column: Column,
    z: int,
) -> cp.Expression:
    kind, y, t = column
    if kind == "b":
        return projective_statistics[z, y, t]
    return probability[z, y] - projective_statistics[z, y, t]


class TernaryConeOracle:
    """DPP SOCP for one chart cell and variable terminal-weight box."""

    def __init__(
        self,
        support_weight: float,
        prefix_order: tuple[int, int, int, int],
        pairs: tuple[tuple[Column, Column], ...],
        coordinate_cases: tuple[str, ...],
        maximum_weight_floor: float = 0.0,
        projective_support_upper: float = 0.76591,
        max_behavior_conditions: int = 0,
        behavior_disjunctions: tuple[
            tuple[tuple[int, tuple[float, float, float, float]], ...], ...
        ] = (),
        disjunction_prior_bounds: tuple[tuple[float, float], ...] | None = None,
        disjunction_source: str = "statistics",
        mip_time_limit: float = 300.0,
        projective_support_lines: tuple[tuple[float, float], ...] = (),
        common_contractions: tuple[dict[str, object], ...] = (),
        terminal_reconstruction: tuple[object, object] | None = None,
    ) -> None:
        self.pairs = pairs
        self.coordinate_cases = coordinate_cases
        self.prefix_order = prefix_order
        constraints: list[cp.Constraint] = []
        self.statistics = cp.Variable((4, 4, 3), nonneg=True)
        self.projective_statistics = cp.Variable((4, 4, 3), nonneg=True)
        self.weights = cp.Variable(3, nonneg=True)
        self.probability = cp.sum(self.statistics, axis=2)
        constraints.append(cp.sum(self.statistics) == 1.0)
        constraints.append(self.projective_statistics <= self.probability[:, :, None])

        # Optional prefix boxes are also used as the pathwise upper endpoints
        # in the q=w*u McCormick envelopes below.  Separate product parameters
        # preserve DPP when both the weight and prefix endpoints vary.
        self.prior_lower = cp.Parameter(4, nonneg=True)
        self.prior_upper = cp.Parameter(4, nonneg=True)
        self.mccormick_lower_cross = cp.Parameter((3, 4), nonneg=True)
        self.mccormick_upper_cross = cp.Parameter((3, 4), nonneg=True)

        # Interval relaxation of q=w*u.  At zero terminal-box width these
        # pairs of inequalities become equality for every path and outcome.
        self.weight_lower = cp.Parameter(3, nonneg=True)
        self.weight_upper = cp.Parameter(3, nonneg=True)
        constraints.extend(
            (
                self.weights >= self.weight_lower,
                self.weights <= self.weight_upper,
                cp.sum(self.weights) == 2.0,
                self.weights[0] >= self.weights[1],
                self.weights[1] >= self.weights[2],
            )
        )
        self.maximum_weight_floor = cp.Parameter(nonneg=True)
        self._maximum_weight_floor = maximum_weight_floor
        constraints.append(self.weights[0] >= self.maximum_weight_floor)

        # The Horwitz parameters obey A=w0/(1-w2) and B=w1/(1-w2).
        # Parameter-box endpoints therefore give exact linear necessary
        # inequalities, avoiding an otherwise very loose independent weight
        # interval relaxation.
        self.terminal_alpha_lower = cp.Parameter(nonneg=True)
        self.terminal_alpha_upper = cp.Parameter(nonneg=True)
        self.terminal_beta_lower = cp.Parameter(nonneg=True)
        self.terminal_beta_upper = cp.Parameter(nonneg=True)
        residual_denominator = 1.0 - self.weights[2]
        constraints.extend(
            (
                self.weights[0]
                >= self.terminal_alpha_lower * residual_denominator,
                self.weights[0]
                <= self.terminal_alpha_upper * residual_denominator,
                self.weights[1]
                >= self.terminal_beta_lower * residual_denominator,
                self.weights[1]
                <= self.terminal_beta_upper * residual_denominator,
            )
        )

        # One global weight vector must serve every path.  The four standard
        # McCormick inequalities outer-convexify q=w*u on each weight box.
        # Prefix ordering supplies constant pathwise upper bounds for u.
        for t in range(3):
            constraints.extend(
                (
                    self.statistics[:, :, t]
                    >= self.weight_lower[t] * self.projective_statistics[:, :, t],
                    self.statistics[:, :, t]
                    <= self.weight_upper[t] * self.projective_statistics[:, :, t],
                )
            )
            for rank, z in enumerate(prefix_order):
                upper = self.prior_upper[z]
                for y in OUTCOMES:
                    normalised = self.projective_statistics[z, y, t]
                    measured = self.statistics[z, y, t]
                    constraints.extend(
                        (
                            measured
                            >= self.weight_upper[t] * normalised
                            + upper * self.weights[t]
                            - self.mccormick_upper_cross[t, z],
                            measured
                            <= self.weight_lower[t] * normalised
                            + upper * self.weights[t]
                            - self.mccormick_lower_cross[t, z],
                        )
                    )

        self.prefix = [cp.sum(self.probability[z, :]) for z in OUTCOMES]
        constraints.extend(
            self.prefix[prefix_order[index]]
            >= self.prefix[prefix_order[index + 1]]
            for index in range(3)
        )
        constraints.append(self.prefix[prefix_order[0]] >= 0.25)
        for rank in range(1, 4):
            constraints.append(
                self.prefix[prefix_order[rank]] <= 1.0 / (rank + 1.0)
            )

        # Optional prefix boxes and full-behaviour Farkas branches.  A branch
        # condition c.B[:,j] <= 0 is relaxed linearly in raw statistics using
        # endpoint bounds on the four positive prefix priors.
        for z in OUTCOMES:
            constraints.extend(
                (
                    self.prefix[z] >= self.prior_lower[z],
                    self.prefix[z] <= self.prior_upper[z],
                )
            )

        # A common quantum instrument is trace-norm contractive even after
        # its classical output flag and the terminal POVM are retained.  For
        # every real coefficient vector c this gives the terminal-geometry-
        # independent necessary condition
        #
        #   sum_{y,t} |sum_z c_z q[z,y,t]|
        #       <= ||sum_z c_z rho[z]||_1.
        #
        # The qubit norm on the right is max(|scalar|, ||Bloch||_2).  A
        # caller supplies one exhaustive spectral branch for each selected
        # coefficient vector.  Vector-active branches use either an exact
        # rotational gauge (``gauge_rank == 0``) or an angular-cap upper
        # envelope.  The measured L1 form is terminal-geometry independent;
        # the optional reconstructed form below retains two visible Bloch
        # coordinates with a box-safe terminal-geometry error budget.
        self.input_vectors = (
            [cp.Variable(3) for _ in OUTCOMES] if common_contractions else []
        )
        for z in OUTCOMES:
            if self.input_vectors:
                constraints.append(cp.SOC(self.prefix[z], self.input_vectors[z]))
        self.common_contraction_values: list[cp.Expression] = []
        self.common_contraction_selectors: list[cp.Variable] = []
        if terminal_reconstruction is not None:
            raw_anchor, raw_errors = terminal_reconstruction
            reconstruction_anchor = (
                raw_anchor
                if isinstance(raw_anchor, cp.Parameter)
                else np.asarray(raw_anchor, dtype=float)
            )
            reconstruction_errors = (
                raw_errors
                if isinstance(raw_errors, cp.Parameter)
                else np.asarray(raw_errors, dtype=float)
            )
            if reconstruction_anchor.shape != (2, 3):
                raise ValueError("terminal reconstruction anchor must have shape (2,3)")
            if reconstruction_errors.shape != (3,):
                raise ValueError("terminal reconstruction errors must be nonnegative")
            if not isinstance(reconstruction_errors, cp.Parameter) and np.any(
                reconstruction_errors < 0.0
            ):
                raise ValueError("terminal reconstruction errors must be nonnegative")
        else:
            reconstruction_anchor = None
            reconstruction_errors = None
        for contraction in common_contractions:
            coefficients = np.asarray(contraction["coefficients"], dtype=float)
            if coefficients.shape != (4,) or np.linalg.norm(coefficients) <= 1e-14:
                raise ValueError("invalid common-instrument contraction coefficients")
            branch = str(contraction["branch"])
            scalar = sum(
                float(coefficients[z]) * self.prefix[z] for z in OUTCOMES
            )
            input_bloch = sum(
                (
                    float(coefficients[z]) * self.input_vectors[z]
                    for z in OUTCOMES
                ),
                cp.Constant(np.zeros(3)),
            )
            signed_statistics = [
                cp.hstack(
                    [
                        sum(
                            float(coefficients[z]) * self.statistics[z, y, t]
                            for z in OUTCOMES
                        )
                        for t in range(3)
                    ]
                )
                for y in OUTCOMES
            ]
            if reconstruction_anchor is None:
                flagged = sum(cp.norm1(vector) for vector in signed_statistics)
            else:
                block_norms: list[cp.Variable] = []
                for y, signed_vector in enumerate(signed_statistics):
                    block_norm = cp.Variable(nonneg=True)
                    trace = cp.sum(signed_vector)
                    # For the true reconstruction R(theta),
                    # ||R(theta)q|| >= ||R0 q|| - sum_t e_t |q_t|.
                    # Positivity of the raw path statistics gives the affine
                    # upper bound |sum_z c_z q_z,t| <= sum_z |c_z| q_z,t.
                    error_budget = sum(
                        abs(float(coefficients[z]))
                        * reconstruction_errors[t]
                        * self.statistics[z, y, t]
                        for z in OUTCOMES
                        for t in range(3)
                    )
                    constraints.extend(
                        (
                            block_norm >= trace,
                            block_norm >= -trace,
                            cp.SOC(
                                block_norm + error_budget,
                                reconstruction_anchor @ signed_vector,
                            ),
                        )
                    )
                    block_norms.append(block_norm)
                flagged = sum(block_norms)
            self.common_contraction_values.append(flagged)
            if branch == "scalar-positive":
                constraints.extend((cp.SOC(scalar, input_bloch), flagged <= scalar))
            elif branch == "scalar-negative":
                constraints.extend((cp.SOC(-scalar, input_bloch), flagged <= -scalar))
            elif branch == "bloch":
                gauge_rank = contraction.get("gauge_rank")
                if gauge_rank == 0:
                    constraints.extend(
                        (
                            input_bloch[0] == 0.0,
                            input_bloch[1] == 0.0,
                            input_bloch[2] >= 0.0,
                            flagged <= input_bloch[2],
                        )
                    )
                else:
                    if gauge_rank == 1:
                        constraints.extend(
                            (input_bloch[1] == 0.0, input_bloch[0] >= 0.0)
                        )
                    cap_data = contraction.get("cap")
                    if isinstance(cap_data, cp.Parameter):
                        if cap_data.shape != (3,):
                            raise ValueError("a scaled contraction cap must have shape (3,)")
                        projection_bound = cap_data @ input_bloch
                    elif cap_data is not None:
                        cap_array = np.asarray(cap_data, dtype=float)
                        if cap_array.shape == (3,):
                            projection_bound = cap_array @ input_bloch
                        elif cap_array.shape == (4,):
                            projection_bound = (
                                cap_array[:3] @ input_bloch / float(cap_array[3])
                            )
                        else:
                            raise ValueError("invalid common-instrument contraction cap")
                    else:
                        raise ValueError("a vector-active contraction needs a cap")
                    constraints.extend(
                        (
                            cp.SOC(projection_bound, input_bloch),
                            flagged <= projection_bound,
                        )
                    )
            elif branch == "spectral-cover":
                raw_caps = tuple(contraction.get("caps", ()))
                if not raw_caps:
                    raise ValueError("a spectral cover requires Bloch caps")
                signs = tuple(contraction.get("scalar_signs", (1, -1)))
                if any(sign not in {1, -1} for sign in signs) or not signs:
                    raise ValueError("spectral-cover scalar signs must be +1 or -1")
                scaled_caps: list[np.ndarray] = []
                for raw_cap in raw_caps:
                    cap_array = np.asarray(raw_cap, dtype=float)
                    if cap_array.shape == (3,):
                        scaled = cap_array
                    elif cap_array.shape == (4,):
                        if cap_array[3] <= 0.0:
                            raise ValueError("spectral-cover cap cosine must be positive")
                        scaled = cap_array[:3] / cap_array[3]
                    else:
                        raise ValueError("invalid spectral-cover cap")
                    scaled_caps.append(scaled)
                coefficient_bound = float(np.max(np.abs(coefficients)))
                cap_scale = max(float(np.linalg.norm(cap)) for cap in scaled_caps)
                big_m = coefficient_bound * (1.0 + max(1.0, cap_scale))
                big_m = float(np.nextafter(big_m, math.inf))
                selectors = cp.Variable(len(signs) + len(scaled_caps), boolean=True)
                self.common_contraction_selectors.append(selectors)
                constraints.append(cp.sum(selectors) == 1.0)
                for index, sign in enumerate(signs):
                    inactive = big_m * (1.0 - selectors[index])
                    signed_scalar = float(sign) * scalar
                    constraints.extend(
                        (
                            cp.SOC(signed_scalar + inactive, input_bloch),
                            flagged <= signed_scalar + inactive,
                        )
                    )
                offset = len(signs)
                for index, cap in enumerate(scaled_caps):
                    inactive = big_m * (1.0 - selectors[offset + index])
                    projection = cap @ input_bloch
                    constraints.extend(
                        (
                            cp.SOC(projection + inactive, input_bloch),
                            flagged <= projection + inactive,
                        )
                    )
            else:
                raise ValueError(f"unknown contraction branch {branch!r}")
        self.max_behavior_conditions = int(max_behavior_conditions)
        self.behavior_disjunctions = behavior_disjunctions
        self.disjunction_prior_bounds = disjunction_prior_bounds
        if disjunction_source not in {"statistics", "projective_statistics"}:
            raise ValueError("unknown behaviour-disjunction source")
        self.disjunction_source = disjunction_source
        self.mip_time_limit = float(mip_time_limit)
        self.branch_coefficients: cp.Parameter | None = None
        self.conditional_behavior: cp.Variable | None = None
        if self.max_behavior_conditions > 0 or (
            self.behavior_disjunctions and self.disjunction_prior_bounds is None
        ):
            self.conditional_behavior = cp.Variable((4, 4, 3), nonneg=True)
            constraints.append(self.conditional_behavior <= 1.0)
            for z in OUTCOMES:
                constraints.append(cp.sum(self.conditional_behavior[z]) == 1.0)
                for y in OUTCOMES:
                    for t in range(3):
                        conditional = self.conditional_behavior[z, y, t]
                        raw = self.statistics[z, y, t]
                        prior = self.prefix[z]
                        # McCormick envelope of raw=prior*conditional for
                        # prior in [l,u] and conditional in [0,1].
                        constraints.extend(
                            (
                                raw >= self.prior_lower[z] * conditional,
                                raw
                                >= self.prior_upper[z] * conditional
                                + prior
                                - self.prior_upper[z],
                                raw <= self.prior_upper[z] * conditional,
                                raw
                                <= self.prior_lower[z] * conditional
                                + prior
                                - self.prior_lower[z],
                            )
                        )
            if self.max_behavior_conditions > 0:
                self.branch_coefficients = cp.Parameter(
                    (self.max_behavior_conditions, 4, 4, 3)
                )
                constraints.extend(
                    cp.sum(
                        cp.multiply(
                            self.branch_coefficients[index], self.conditional_behavior
                        )
                    )
                    <= 0.0
                    for index in range(self.max_behavior_conditions)
                )

        # A nested-ellipsoid Farkas witness supported on two columns says
        # that every physical behaviour obeys at least one of two linear
        # inequalities.  Encode each exact disjunction with one binary
        # selector.  Accumulating clauses lets SCIP share the conic
        # relaxation instead of materialising an exponentially branching
        # witness tree.
        self.disjunction_selectors: list[cp.Variable] = []
        for disjunction in self.behavior_disjunctions:
            if len(disjunction) != 2:
                raise ValueError("behaviour disjunctions must have two branches")
            selector = cp.Variable(boolean=True)
            self.disjunction_selectors.append(selector)
            expressions: list[cp.Expression] = []
            upper_bounds: list[float] = []
            for column, coefficients in disjunction:
                y, t = divmod(column, 3)
                vector = np.asarray(coefficients, dtype=float)
                if self.disjunction_prior_bounds is None:
                    if self.conditional_behavior is None:
                        raise RuntimeError("conditional behaviour was not constructed")
                    expressions.append(vector @ self.conditional_behavior[:, y, t])
                    upper_bounds.append(
                        max(0.0, float(np.maximum(vector, 0.0).sum()))
                    )
                else:
                    if len(self.disjunction_prior_bounds) != 4:
                        raise ValueError("expected four prior intervals")
                    relaxed = np.empty(4, dtype=float)
                    prior_upper = np.empty(4, dtype=float)
                    for z, coefficient in enumerate(vector):
                        lower, upper = self.disjunction_prior_bounds[z]
                        if lower <= 0.0:
                            raise ValueError(
                                "positive prior lower bounds are required for raw disjunctions"
                            )
                        relaxed[z] = coefficient / (
                            upper if coefficient >= 0.0 else lower
                        )
                        prior_upper[z] = upper
                    raw_source = (
                        self.statistics
                        if self.disjunction_source == "statistics"
                        else self.projective_statistics
                    )
                    expressions.append(relaxed @ raw_source[:, y, t])
                    upper_bounds.append(
                        max(
                            0.0,
                            float(np.sum(np.maximum(relaxed, 0.0) * prior_upper)),
                        )
                    )
            constraints.extend(
                (
                    expressions[0] <= upper_bounds[0] * selector,
                    expressions[1] <= upper_bounds[1] * (1.0 - selector),
                )
            )

        # Each entry contains anchor SOC parameters and one coefficientwise
        # lower-envelope SOC.  Index zero is the terminal POVM itself.
        self.soc_parameters: list[
            tuple[
                list[tuple[cp.Parameter, cp.Parameter, cp.Parameter]],
                tuple[cp.Parameter, cp.Parameter, cp.Parameter],
            ]
        ] = []

        def parameter_family() -> tuple[
            list[tuple[cp.Parameter, cp.Parameter, cp.Parameter]],
            tuple[cp.Parameter, cp.Parameter, cp.Parameter],
        ]:
            anchors = [
                (
                    cp.Parameter((2, 2)),
                    cp.Parameter(2),
                    cp.Parameter(nonneg=True),
                )
                for _ in ANCHOR_LOCATIONS
            ]
            lower = (
                cp.Parameter((2, 2)),
                cp.Parameter(2),
                cp.Parameter(nonneg=True),
            )
            return anchors, lower

        terminal_family = parameter_family()
        self.soc_parameters.append(terminal_family)

        def add_cone(
            vector: cp.Expression,
            family: tuple[
                list[tuple[cp.Parameter, cp.Parameter, cp.Parameter]],
                tuple[cp.Parameter, cp.Parameter, cp.Parameter],
            ],
        ) -> None:
            constraints.append(vector >= 0.0)
            scale = cp.sum(vector)
            constraints.extend(
                vector[t] <= self.weight_upper[t] * scale for t in range(3)
            )
            point = cp.hstack([vector[0], vector[1]])
            anchors, lower = family
            for root, offset, radius in (*anchors, lower):
                constraints.append(
                    cp.SOC(radius * scale, root @ point + offset * scale)
                )

        # Every conditioned output is a qubit state measured by the same
        # terminal POVM.  This is stronger than the old componentwise cap and
        # remains a necessary condition for the physical common instrument.
        for z, y in PATHS:
            add_cone(self.statistics[z, y, :], terminal_family)

        terminal_statistics = [
            cp.hstack(
                [
                    sum(
                        self.statistics[z, y, t]
                        for z, y in PATHS
                        if (z ^ y) == syndrome
                    )
                    for t in range(3)
                ]
            )
            for syndrome in OUTCOMES
        ]
        terminal_projective_statistics = [
            cp.hstack(
                [
                    sum(
                        self.projective_statistics[z, y, t]
                        for z, y in PATHS
                        if (z ^ y) == syndrome
                    )
                    for t in range(3)
                ]
            )
            for syndrome in OUTCOMES
        ]
        self.helstrom_probabilities = cp.Variable(3, nonneg=True)
        add_cone(self.helstrom_probabilities, terminal_family)
        for syndrome in OUTCOMES:
            add_cone(
                self.helstrom_probabilities - terminal_statistics[syndrome],
                terminal_family,
            )

        self.audit = sum(terminal_statistics[s][s] for s in range(3))
        constraints.append(self.audit == cp.sum(self.helstrom_probabilities))
        self.cap_weights = cp.Parameter(4, nonneg=True)
        constraints.append(
            self.audit
            <= sum(
                self.cap_weights[index] * self.prefix[prefix_order[index]]
                for index in OUTCOMES
            )
        )
        terminal_prior = [cp.sum(vector) for vector in terminal_statistics]
        constraints.append(
            self.audit
            <= sum(self.weight_upper[t] * terminal_prior[t] for t in range(3))
        )

        for (first, second), coordinate_case in zip(
            pairs, coordinate_cases, strict=True
        ):
            family = parameter_family()
            self.soc_parameters.append(family)
            anchors, lower = family
            for z in OUTCOMES:
                x_value = selected_normalised(
                    self.projective_statistics, self.probability, first, z
                )
                y_value = selected_normalised(
                    self.projective_statistics, self.probability, second, z
                )
                residual = self.prefix[z] - x_value - y_value
                constraints.append(residual >= 0.0)
                if coordinate_case in {"xy", "center"}:
                    point = cp.hstack([x_value, y_value])
                elif coordinate_case == "xr":
                    point = cp.hstack([x_value, residual])
                elif coordinate_case == "yr":
                    point = cp.hstack([y_value, residual])
                else:
                    raise ValueError(f"unknown coordinate case {coordinate_case!r}")
                for root, offset, radius in (*anchors, lower):
                    constraints.append(
                        cp.SOC(
                            radius * self.prefix[z],
                            root @ point + offset * self.prefix[z],
                        )
                    )

        self.returned = hellinger_hypograph(
            [self.probability[z, y] for z, y in PATHS], constraints
        )
        self.score = support_weight * self.audit + (1.0 - support_weight) * self.returned

        # Compare the ternary readout with every aligned binary-projective
        # readout.  If i is retained as Pi_i, j is the complementary answer,
        # and k is deleted, then
        #
        #   A_3-A_2 <= (1-w_i) p_j + w_k p_k.
        #
        # The same leaf and RETURN term form a legal projective protocol, so
        # its support is bounded by the independently covered projective
        # frontier.  Concave McCormick envelopes keep this comparison valid
        # throughout a terminal-weight box.
        comparison_lines = (
            *projective_support_lines,
            (support_weight, projective_support_upper),
        )
        self.exact_projective_constraints: list[
            tuple[float, float, int, int, cp.Constraint]
        ] = []

        # The normalized statistics make every aligned binary-projective
        # replacement available exactly, not merely through a loss bound.
        # Retain projector Pi_i as answer i and use I-Pi_i as answer j:
        #
        #   A_2 = u_{i,i} + p_j - u_{j,i}.
        #
        # This is a legal projective protocol with the same coarse outcomes
        # and hence the same RETURN score.  At nonzero terminal-box width the
        # q=w*u link is outer-convexified, but every physical point still has
        # its exact normalized u and therefore satisfies all these cuts.
        for i in range(3):
            for j in range(4):
                if j == i:
                    continue
                projective_audit = (
                    terminal_projective_statistics[i][i]
                    + terminal_prior[j]
                    - terminal_projective_statistics[j][i]
                )
                for line_weight, line_upper in comparison_lines:
                    constraint = (
                        line_weight * projective_audit
                        + (1.0 - line_weight) * self.returned
                        <= line_upper
                    )
                    constraints.append(constraint)
                    self.exact_projective_constraints.append(
                        (line_weight, line_upper, i, j, constraint)
                    )

        self.comparison_products: list[tuple[cp.Variable, cp.Variable]] = []
        for i in range(3):
            for j in range(3):
                if j == i:
                    continue
                k = 3 - i - j
                first = cp.Variable(nonneg=True)
                second = cp.Variable(nonneg=True)
                self.comparison_products.append((first, second))
                missing_i = 1.0 - self.weights[i]
                missing_lower = 1.0 - self.weight_upper[i]
                missing_upper = 1.0 - self.weight_lower[i]
                constraints.extend(
                    (
                        first <= missing_upper * terminal_prior[j],
                        first
                        <= missing_lower * terminal_prior[j]
                        + missing_i
                        - missing_lower,
                        second <= self.weight_upper[k] * terminal_prior[k],
                        second
                        <= self.weight_lower[k] * terminal_prior[k]
                        + self.weights[k]
                        - self.weight_lower[k],
                    )
                )
                for line_weight, line_upper in comparison_lines:
                    if not 0.0 < line_weight <= 1.0:
                        raise ValueError(
                            "projective support-line weights must lie in (0,1]"
                        )
                    constraints.append(
                        line_weight * self.audit
                        + (1.0 - line_weight) * self.returned
                        <= line_upper + line_weight * (first + second)
                    )
        self.problem = cp.Problem(cp.Maximize(self.score), constraints)
        if not self.problem.is_dpp():
            raise RuntimeError("ternary probability-cone oracle is not DPP")

    @staticmethod
    def assign_soc(
        targets: tuple[cp.Parameter, cp.Parameter, cp.Parameter],
        data: tuple[np.ndarray, np.ndarray, float] | None,
    ) -> None:
        root, offset, radius = targets
        if data is None:
            root.value = np.zeros((2, 2))
            offset.value = np.zeros(2)
            radius.value = 1.0
        else:
            root_value, shift, radius_value = data
            root.value = root_value
            offset.value = root_value @ shift
            radius.value = radius_value

    def assign_family(
        self,
        family_index: int,
        alpha_bounds: tuple[float, float],
        beta_bounds: tuple[float, float],
        center_chart: bool = False,
    ) -> None:
        anchors, lower = self.soc_parameters[family_index]
        values = (
            center_box_anchor_relaxations(alpha_bounds, beta_bounds)
            if center_chart
            else box_anchor_relaxations(alpha_bounds, beta_bounds)
        )
        for (alpha, beta, error), targets in zip(values, anchors, strict=True):
            try:
                data = (
                    center_inellipse_soc_data(alpha, beta, error)
                    if center_chart
                    else inellipse_soc_data(alpha, beta, error)
                )
            except ValueError:
                data = None
            self.assign_soc(targets, data)
        self.assign_soc(
            lower,
            (
                center_coefficientwise_box_soc_data(alpha_bounds, beta_bounds)
                if center_chart
                else coefficientwise_box_soc_data(alpha_bounds, beta_bounds)
            ),
        )

    def solve(
        self,
        box: Box,
        safety: float,
        capture: bool = False,
        behavior_conditions: tuple[
            tuple[int, tuple[float, float, float, float]], ...
        ] = (),
    ) -> dict[str, Any]:
        intervals = terminal_weight_intervals(box)
        self.weight_lower.value = np.asarray([item[0] for item in intervals])
        self.weight_upper.value = np.asarray([item[1] for item in intervals])
        self.maximum_weight_floor.value = self._maximum_weight_floor
        al, au = box[TERMINAL_ALPHA]
        bl, bu = box[TERMINAL_BETA]
        self.terminal_alpha_lower.value = al
        self.terminal_alpha_upper.value = au
        self.terminal_beta_lower.value = bl
        self.terminal_beta_upper.value = bu
        rank_upper = (1.0, 0.5, 1.0 / 3.0, 0.25)
        prior_lower = np.zeros(4, dtype=float)
        prior_upper = np.ones(4, dtype=float)
        for rank, z in enumerate(self.prefix_order):
            prior_upper[z] = rank_upper[rank]
        prior_lower[self.prefix_order[0]] = 0.25
        for z in OUTCOMES:
            if f"prior_{z}" in box:
                prior_lower[z], prior_upper[z] = box[f"prior_{z}"]
        self.prior_lower.value = prior_lower
        self.prior_upper.value = prior_upper
        self.mccormick_lower_cross.value = np.outer(
            np.asarray([item[0] for item in intervals]), prior_upper
        )
        self.mccormick_upper_cross.value = np.outer(
            np.asarray([item[1] for item in intervals]), prior_upper
        )
        if len(behavior_conditions) > self.max_behavior_conditions:
            return {"status": "condition_overflow", "bound": math.inf}
        if self.branch_coefficients is not None:
            coefficients = np.zeros(
                (self.max_behavior_conditions, 4, 4, 3), dtype=float
            )
            for index, (column, vector) in enumerate(behavior_conditions):
                y, t = divmod(column, 3)
                for z, coefficient in enumerate(vector):
                    coefficients[index, z, y, t] = coefficient
            self.branch_coefficients.value = coefficients
        maximum = intervals[0][1]
        self.cap_weights.value = np.asarray(
            [maximum, maximum, max(0.0, 2.0 - 2.0 * maximum), 0.0]
        )
        self.assign_family(
            0, box[TERMINAL_ALPHA], box[TERMINAL_BETA]
        )
        for pair_index in range(len(self.pairs)):
            center_chart = self.coordinate_cases[pair_index] == "center"
            self.assign_family(
                pair_index + 1,
                box[
                    f"k{pair_index}_{'cx' if center_chart else 'alpha'}"
                ],
                box[
                    f"k{pair_index}_{'cy' if center_chart else 'beta'}"
                ],
                center_chart,
            )
        try:
            if self.behavior_disjunctions or self.common_contraction_selectors:
                self.problem.solve(
                    solver="SCIP",
                    verbose=False,
                    scip_params={
                        "limits/time": self.mip_time_limit,
                        "numerics/feastol": 1e-8,
                        "numerics/epsilon": 1e-9,
                        "display/verblevel": 0,
                    },
                    ignore_dpp=False,
                )
            else:
                self.problem.solve(
                    solver="CLARABEL",
                    tol_gap_abs=2e-8,
                    tol_gap_rel=2e-8,
                    tol_feas=2e-8,
                    max_iter=1000,
                    warm_start=True,
                    ignore_dpp=False,
                )
        except cp.SolverError as error:
            return {"status": "solver_error", "error": str(error), "bound": math.inf}
        if self.problem.status in {cp.INFEASIBLE, cp.INFEASIBLE_INACCURATE}:
            return {"status": self.problem.status, "bound": -math.inf}
        if self.problem.status not in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}:
            return {"status": self.problem.status, "bound": math.inf}
        result: dict[str, Any] = {
            "status": self.problem.status,
            "raw_value": float(self.problem.value),
            "bound": float(self.problem.value) + safety,
            "audit": float(self.audit.value),
            "return": float(self.returned.value),
            "terminal_weight_intervals": [list(item) for item in intervals],
            "prior_intervals": np.column_stack(
                [prior_lower, prior_upper]
            ).tolist(),
            "iterations": self.problem.solver_stats.num_iters,
        }
        if capture:
            result.update(
                {
                    "prefix": [float(item.value) for item in self.prefix],
                    "weights": np.asarray(self.weights.value).tolist(),
                    "probability": np.asarray(self.probability.value).tolist(),
                    "statistics": np.asarray(self.statistics.value).tolist(),
                    "projective_statistics": np.asarray(
                        self.projective_statistics.value
                    ).tolist(),
                    "input_bloch_vectors": [
                        np.asarray(vector.value).tolist()
                        for vector in self.input_vectors
                    ],
                    "common_contraction_values": [
                        float(value.value) for value in self.common_contraction_values
                    ],
                    "common_contraction_selectors": [
                        np.asarray(selector.value).tolist()
                        for selector in self.common_contraction_selectors
                    ],
                }
            )
            if self.conditional_behavior is not None:
                result["conditional_behavior"] = np.asarray(
                    self.conditional_behavior.value
                ).tolist()
            result["exact_projective_line_duals"] = [
                {
                    "weight": line_weight,
                    "upper": line_upper,
                    "retained": retained,
                    "complement": complement,
                    "dual": (
                        None
                        if constraint.dual_value is None
                        else float(constraint.dual_value)
                    ),
                }
                for (
                    line_weight,
                    line_upper,
                    retained,
                    complement,
                    constraint,
                ) in self.exact_projective_constraints
            ]
        return result


def solve_node_reusable(
    box: Box,
    support_weight: float,
    prefix_order: tuple[int, int, int, int],
    pairs: tuple[tuple[Column, Column], ...],
    coordinate_cases: tuple[str, ...],
    maximum_weight_floor: float,
    projective_support_upper: float,
    safety: float,
    capture: bool = False,
    max_behavior_conditions: int = 0,
    behavior_conditions: tuple[
        tuple[int, tuple[float, float, float, float]], ...
    ] = (),
) -> dict[str, Any]:
    key = (
        float(support_weight),
        prefix_order,
        pairs,
        coordinate_cases,
        float(maximum_weight_floor),
        float(projective_support_upper),
        int(max_behavior_conditions),
    )
    oracle = _ORACLE_CACHE.get(key)
    if oracle is None:
        oracle = TernaryConeOracle(
            support_weight,
            prefix_order,
            pairs,
            coordinate_cases,
            maximum_weight_floor,
            projective_support_upper,
            max_behavior_conditions,
        )
        _ORACLE_CACHE[key] = oracle
    return oracle.solve(box, safety, capture, behavior_conditions)


def branch_variable(box: Box, root: Box) -> str:
    return max(
        box,
        key=lambda name: (box[name][1] - box[name][0])
        / max(root[name][1] - root[name][0], 1e-15),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lambda", dest="support_weight", type=float, default=0.6)
    parser.add_argument("--maximum-weight-floor", type=float, default=0.88325)
    parser.add_argument("--projective-support-upper", type=float, default=0.76591)
    parser.add_argument("--prefix-order", type=int, nargs=4, default=(0, 1, 2, 3))
    parser.add_argument("--pair", action="append", default=[])
    parser.add_argument(
        "--coordinate-case",
        action="append",
        choices=("xy", "xr", "yr", "center"),
        default=[],
    )
    parser.add_argument("--target", type=float, default=0.76591)
    parser.add_argument("--safety", type=float, default=1e-5)
    parser.add_argument("--max-nodes", type=int, default=101)
    parser.add_argument("--root-box", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    pairs = tuple(parse_pair(text) for text in args.pair)
    cases = tuple(args.coordinate_case)
    if len(pairs) != len(cases):
        raise ValueError("supply one coordinate case for every selected pair")
    root = initial_box(len(pairs), args.maximum_weight_floor)
    if args.root_box is not None:
        payload = json.loads(args.root_box.read_text(encoding="utf-8"))
        if "box" in payload:
            payload = payload["box"]
        root = deserialise_box(payload)
    oracle = TernaryConeOracle(
        args.support_weight,
        tuple(args.prefix_order),
        pairs,
        cases,
        args.maximum_weight_floor,
        args.projective_support_upper,
    )
    result = oracle.solve(root, args.safety)
    solved = 1
    counter = 0
    queue: list[tuple[float, int, Box, dict[str, Any]]] = [
        (-float(result["bound"]), counter, root, result)
    ]
    leaves: list[dict[str, Any]] = []
    while queue:
        negative_bound, node_id, box, result = heapq.heappop(queue)
        if -negative_bound <= args.target:
            leaves.append({"box": serialise_box(box), **result})
            continue
        if solved >= args.max_nodes:
            heapq.heappush(queue, (negative_bound, node_id, box, result))
            break
        name = branch_variable(box, root)
        for child in split_box(box, name):
            if not terminal_domain_intersects(child, args.maximum_weight_floor):
                leaves.append(
                    {"box": serialise_box(child), "bound": -math.inf, "status": "domain_empty"}
                )
                continue
            child_result = oracle.solve(child, args.safety)
            solved += 1
            counter += 1
            bound = float(child_result["bound"])
            if bound <= args.target:
                leaves.append({"box": serialise_box(child), **child_result})
            else:
                heapq.heappush(queue, (-bound, counter, child, child_result))
        current = -queue[0][0] if queue else max(
            (float(item["bound"]) for item in leaves), default=-math.inf
        )
        print(solved, name, current, len(queue), len(leaves), flush=True)

    open_nodes = [
        {"box": serialise_box(box), **result}
        for _, _, box, result in sorted(queue)
    ]
    top = None
    if queue:
        _, _, top_box, _ = queue[0]
        top = {"box": serialise_box(top_box), **oracle.solve(top_box, args.safety, True)}
    payload = {
        "support_weight": args.support_weight,
        "maximum_weight_floor": args.maximum_weight_floor,
        "projective_support_upper": args.projective_support_upper,
        "prefix_order": list(args.prefix_order),
        "pairs": [[render_column(item) for item in pair] for pair in pairs],
        "coordinate_cases": list(cases),
        "target": args.target,
        "safety": args.safety,
        "root_box": serialise_box(root),
        "complete": not open_nodes,
        "solved_nodes": solved,
        "maximum_open_bound": max(
            (float(item["bound"]) for item in open_nodes), default=-math.inf
        ),
        "maximum_leaf_bound": max(
            (float(item["bound"]) for item in leaves), default=-math.inf
        ),
        "open_nodes": open_nodes,
        "leaves": leaves,
        "top_open_solution": top,
        "scope": (
            "one Cartesian product of selected-pair charts over the continuous "
            "sorted ternary terminal-POVM strip"
        ),
        "numerical_status": (
            "finite SOCP outer cover at solver tolerances; dual outward validation pending"
        ),
    }
    write_payload(args.output, payload)
    print(
        json.dumps(
            {key: value for key, value in payload.items() if key not in {"open_nodes", "leaves"}},
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

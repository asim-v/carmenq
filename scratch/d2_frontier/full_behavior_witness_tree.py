"""Disjunctive full-behaviour cover of a fixed ternary frontier cell.

The convex probability/Helstrom outer model can attain a high score with a
``4 x 12`` conditional behaviour that has no qubit factorisation.  A dual
nested-ellipsoid obstruction supported on outcome columns ``J`` says that a
physical behaviour must satisfy at least one disjunct

    c_j . B[:, j] <= 0,  j in J.

For fixed prefix priors this is linear in the joint statistics.  On a prior
box, replacing each reciprocal prior by the endpoint that gives a lower
bound produces a valid linear outer relaxation.  Branching on the disjuncts
and, when necessary, bisecting the four prior intervals gives a convergent
cover without boxing all forty-eight conditional probabilities.

The archived tree is a numerical certificate at the reported conic solver
tolerances.  A publication-grade formal certificate still requires outward
rounding of the SOCP duals and exact PSD validation of every ellipsoid dual.
"""

from __future__ import annotations

import argparse
import heapq
import itertools
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import cvxpy as cp
import numpy as np

from full_behavior_psd_rank_certificate import solve_dual, solve_primal
from pairwise_inellipse_box_cover import (
    OUTCOMES,
    PATHS,
    canonical_three_effect_povm,
    hellinger_hypograph,
    reconstruction_matrix,
)
from terminal_weight_upper import filled_effect_weights


PriorBox = np.ndarray  # shape (4, 2), lower/upper


@dataclass(frozen=True)
class BranchCondition:
    """One exact disjunct c . B[:, column] <= 0."""

    column: int
    coefficients: tuple[float, float, float, float]

    def to_json(self) -> dict[str, Any]:
        return {"column": self.column, "coefficients": list(self.coefficients)}

    @staticmethod
    def from_json(payload: dict[str, Any]) -> "BranchCondition":
        coefficients = tuple(float(x) for x in payload["coefficients"])
        if len(coefficients) != 4:
            raise ValueError("branch condition must have four coefficients")
        return BranchCondition(int(payload["column"]), coefficients)  # type: ignore[arg-type]


class FullBehaviorOuter:
    """Reusable fixed-terminal SOCP with parameterized prior-box cuts."""

    def __init__(
        self,
        terminal: np.ndarray,
        support_weight: float,
        prefix_order: tuple[int, int, int, int],
        max_conditions: int,
    ) -> None:
        self.weights = np.trace(terminal[:3], axis1=1, axis2=2).real
        self.support_weight = float(support_weight)
        self.max_conditions = int(max_conditions)
        constraints: list[cp.Constraint] = []
        self.statistics = cp.Variable((4, 4, 3), nonneg=True)
        self.probability = cp.sum(self.statistics, axis=2)
        constraints.append(cp.sum(self.statistics) == 1.0)
        constraints.extend(
            self.statistics[:, :, t] <= self.weights[t] * self.probability
            for t in range(3)
        )
        self.prefix = cp.hstack(
            [cp.sum(self.probability[z, :]) for z in OUTCOMES]
        )
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

        self.prior_lower = cp.Parameter(4, nonneg=True)
        self.prior_upper = cp.Parameter(4, nonneg=True)
        for z in OUTCOMES:
            constraints.extend(
                (
                    self.prefix[z] >= self.prior_lower[z],
                    self.prefix[z] <= self.prior_upper[z],
                )
            )

        # Each row contains the coefficient of one raw statistic g[z,y,t].
        # Unused rows are zero and impose the harmless inequality 0 <= 0.
        self.branch_coefficients = cp.Parameter(
            (self.max_conditions, 4, 4, 3)
        )
        constraints.extend(
            cp.sum(cp.multiply(self.branch_coefficients[k], self.statistics))
            <= 0.0
            for k in range(self.max_conditions)
        )

        terminal_statistics = [
            [
                sum(
                    self.statistics[z, y, t]
                    for z, y in PATHS
                    if (z ^ y) == syndrome
                )
                for t in range(3)
            ]
            for syndrome in OUTCOMES
        ]
        inverse = reconstruction_matrix(terminal)
        terminal_prior: list[cp.Expression] = []
        terminal_vector: list[cp.Expression] = []
        normal = cp.Variable(4)
        for syndrome in OUTCOMES:
            reconstructed = inverse @ cp.hstack(terminal_statistics[syndrome])
            vector = cp.hstack(
                [reconstructed[1], reconstructed[2], normal[syndrome]]
            )
            constraints.extend(
                (
                    reconstructed[0]
                    == sum(
                        self.probability[z, z ^ syndrome]
                        for z in OUTCOMES
                    ),
                    cp.SOC(reconstructed[0], vector),
                )
            )
            terminal_prior.append(reconstructed[0])
            terminal_vector.append(vector)

        self.audit = sum(terminal_statistics[s][s] for s in range(3))
        cap = filled_effect_weights(float(self.weights.max()))
        constraints.append(
            self.audit
            <= sum(
                cap[index] * self.prefix[prefix_order[index]]
                for index in OUTCOMES
            )
        )
        dual_scalar = cp.Variable(nonneg=True)
        dual_vector = cp.Variable(3)
        constraints.append(cp.SOC(dual_scalar, dual_vector))
        constraints.extend(
            cp.SOC(
                dual_scalar - terminal_prior[s],
                dual_vector - terminal_vector[s],
            )
            for s in OUTCOMES
        )
        constraints.append(self.audit == dual_scalar)
        constraints.append(
            self.audit
            <= sum(self.weights[t] * terminal_prior[t] for t in range(3))
        )
        self.returned = hellinger_hypograph(
            [self.probability[z, y] for z, y in PATHS], constraints
        )
        self.score = (
            support_weight * self.audit
            + (1.0 - support_weight) * self.returned
        )
        self.problem = cp.Problem(cp.Maximize(self.score), constraints)

    @staticmethod
    def relaxed_coefficients(
        condition: BranchCondition, box: PriorBox
    ) -> np.ndarray:
        """Coefficients of a lower bound on c_z g_z / a_z."""

        result = np.zeros((4, 4, 3), dtype=float)
        y, t = divmod(condition.column, 3)
        for z, coefficient in enumerate(condition.coefficients):
            lower, upper = box[z]
            if lower <= 0.0:
                raise ValueError("positive prior lower bounds are required")
            denominator = upper if coefficient >= 0.0 else lower
            result[z, y, t] = coefficient / denominator
        return result

    def solve(
        self,
        box: PriorBox,
        conditions: tuple[BranchCondition, ...],
        safety: float,
    ) -> dict[str, Any]:
        if len(conditions) > self.max_conditions:
            return {"status": "condition_overflow", "bound": math.inf}
        self.prior_lower.value = box[:, 0]
        self.prior_upper.value = box[:, 1]
        coefficients = np.zeros((self.max_conditions, 4, 4, 3), dtype=float)
        for index, condition in enumerate(conditions):
            coefficients[index] = self.relaxed_coefficients(condition, box)
        self.branch_coefficients.value = coefficients
        try:
            self.problem.solve(
                solver="CLARABEL",
                tol_gap_abs=2e-9,
                tol_gap_rel=2e-9,
                tol_feas=2e-9,
                max_iter=1000,
                warm_start=True,
            )
        except cp.SolverError as error:
            return {"status": "solver_error", "error": str(error), "bound": math.inf}
        if self.problem.status in {cp.INFEASIBLE, cp.INFEASIBLE_INACCURATE}:
            return {"status": self.problem.status, "bound": -math.inf}
        if self.problem.status not in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}:
            return {"status": self.problem.status, "bound": math.inf}
        prefix = np.asarray(self.prefix.value, dtype=float)
        statistics = np.asarray(self.statistics.value, dtype=float)
        behavior = statistics.reshape(4, -1) / prefix[:, None]
        exact_conditions = [
            float(
                np.dot(
                    np.asarray(condition.coefficients),
                    behavior[:, condition.column],
                )
            )
            for condition in conditions
        ]
        return {
            "status": self.problem.status,
            "raw_value": float(self.problem.value),
            "bound": float(self.problem.value) + safety,
            "audit": float(self.audit.value),
            "return": float(self.returned.value),
            "prefix": prefix.tolist(),
            "statistics": statistics.tolist(),
            "behavior": behavior.tolist(),
            "exact_condition_values": exact_conditions,
            "iterations": self.problem.solver_stats.num_iters,
            "solve_time": self.problem.solver_stats.solve_time,
        }


def initial_prior_box(order: tuple[int, int, int, int]) -> PriorBox:
    box = np.zeros((4, 2), dtype=float)
    box[:, 1] = 1.0
    box[order[0], 0] = 0.25
    for rank in range(1, 4):
        box[order[rank], 1] = 1.0 / (rank + 1.0)
    return box


def tighten_prior_box(
    oracle: FullBehaviorOuter,
    box: PriorBox,
    conditions: tuple[BranchCondition, ...],
    target: float,
    safety: float,
    steps: int,
) -> tuple[PriorBox, int]:
    """Coordinate hull of the target superlevel set."""

    result = box.copy()
    solves = 0
    for z in OUTCOMES:
        # Lower endpoint: retain only prefix[z] <= midpoint.
        left, right = result[z]
        for _ in range(steps):
            middle = 0.5 * (left + right)
            trial = result.copy()
            trial[z, 1] = middle
            node = oracle.solve(trial, conditions, safety)
            solves += 1
            if float(node["bound"]) >= target:
                right = middle
            else:
                left = middle
        result[z, 0] = left

        # Upper endpoint: retain only prefix[z] >= midpoint.
        left, right = result[z]
        for _ in range(steps):
            middle = 0.5 * (left + right)
            trial = result.copy()
            trial[z, 0] = middle
            node = oracle.solve(trial, conditions, safety)
            solves += 1
            if float(node["bound"]) >= target:
                left = middle
            else:
                right = middle
        result[z, 1] = right
    return result, solves


def dual_is_numerically_valid(dual: dict[str, Any], tolerance: float) -> bool:
    return (
        dual.get("status") in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}
        and float(dual.get("certified_common_margin", -math.inf)) > tolerance
        and float(dual.get("state_dual_min_eigenvalue", -math.inf)) > -tolerance
        and min(dual.get("containment_dual_min_eigenvalues", [-math.inf]))
        > -tolerance
        and float(dual.get("stationarity_frobenius_residual", math.inf))
        < 10.0 * tolerance
    )


PREFERRED_PAIRS = (
    (0, 6),
    (0, 9),
    (0, 10),
    (3, 7),
    (1, 9),
    (0, 11),
    (3, 8),
    (2, 3),
    (2, 9),
    (0, 4),
    (1, 6),
    (9, 11),
)

# Minimal three-column cores repeatedly exposed by the symmetric hard face.
# Singletons and all pairs have been exhaustively ruled out for the reference
# maximizer, so triples are the first genuinely stronger disjunctions.
PREFERRED_TRIPLES = (
    (0, 3, 6),
    (0, 9, 10),
    (3, 6, 8),
    (0, 9, 11),
    (3, 9, 10),
    (1, 6, 8),
    (1, 8, 9),
    (4, 6, 11),
    (0, 1, 10),
    (1, 3, 9),
    (0, 1, 2),
)


def find_pair_witnesses(
    behavior: np.ndarray,
    robust_budget: float,
    tolerance: float,
    exhaustive: bool,
    candidate_columns: tuple[int, ...] | None = None,
) -> list[dict[str, Any]]:
    allowed = (
        set(range(behavior.shape[1]))
        if candidate_columns is None
        else set(candidate_columns)
    )
    pairs = [pair for pair in PREFERRED_PAIRS if set(pair) <= allowed]
    if exhaustive:
        pairs.extend(
            pair
            for pair in itertools.combinations(sorted(allowed), 2)
            if pair not in pairs
        )
    found: list[tuple[float, dict[str, Any]]] = []
    for pair in pairs:
        try:
            dual = solve_dual(
                behavior, list(pair), robust_budget=robust_budget
            )
        except cp.SolverError:
            continue
        if not dual_is_numerically_valid(dual, tolerance):
            continue
        coefficients = np.asarray(dual["halfspace_linear_coefficients"], dtype=float)
        margin = float(dual["certified_common_margin"])
        scale = max(float(np.max(np.sum(np.abs(coefficients), axis=1))), 1e-15)
        quality = margin / scale
        dual["robust_quality"] = quality
        found.append((quality, dual))
    found.sort(key=lambda item: item[0], reverse=True)
    return [item[1] for item in found]


def find_pair_witness(
    behavior: np.ndarray,
    robust_budget: float,
    tolerance: float,
    exhaustive: bool,
    candidate_columns: tuple[int, ...] | None = None,
) -> dict[str, Any] | None:
    """Return the numerically strongest available two-column obstruction."""

    witnesses = find_pair_witnesses(
        behavior,
        robust_budget,
        tolerance,
        exhaustive,
        candidate_columns,
    )
    return witnesses[0] if witnesses else None


def find_small_witness(
    behavior: np.ndarray,
    robust_budget: float,
    tolerance: float,
    exhaustive_pairs: bool,
    exhaustive_triples: bool,
) -> dict[str, Any] | None:
    """Return the strongest validated two- or three-column obstruction."""

    pair = find_pair_witness(
        behavior,
        robust_budget,
        tolerance,
        exhaustive_pairs,
    )
    if pair is not None:
        return pair

    triples = list(PREFERRED_TRIPLES)
    if exhaustive_triples:
        triples.extend(
            triple
            for triple in itertools.combinations(range(behavior.shape[1]), 3)
            if triple not in triples
        )
    found: list[tuple[float, dict[str, Any]]] = []
    for triple in triples:
        try:
            dual = solve_dual(
                behavior,
                list(triple),
                robust_budget=robust_budget,
            )
        except cp.SolverError:
            continue
        if not dual_is_numerically_valid(dual, tolerance):
            continue
        coefficients = np.asarray(
            dual["halfspace_linear_coefficients"], dtype=float
        )
        margin = float(dual["certified_common_margin"])
        scale = max(
            float(np.max(np.sum(np.abs(coefficients), axis=1))), 1e-15
        )
        dual["robust_quality"] = margin / scale
        found.append((margin / scale, dual))
    found.sort(key=lambda item: item[0], reverse=True)
    return found[0][1] if found else None


def condition_key(condition: BranchCondition) -> tuple[int, tuple[int, ...]]:
    vector = np.asarray(condition.coefficients, dtype=float)
    scale = max(float(np.max(np.abs(vector))), 1e-300)
    normalized = vector / scale
    # Positive rescaling leaves the disjunct unchanged.
    return condition.column, tuple(int(round(x * 1e8)) for x in normalized)


def split_prior_box(
    box: PriorBox, coordinate: int, preferred: float | None = None
) -> tuple[PriorBox, PriorBox]:
    lower, upper = box[coordinate]
    midpoint = 0.5 * (lower + upper)
    if preferred is not None:
        # Splitting through the current relaxation maximizer tightens the
        # active reciprocal immediately while retaining two well-sized boxes.
        midpoint = min(
            max(float(preferred), lower + 0.2 * (upper - lower)),
            upper - 0.2 * (upper - lower),
        )
    first, second = box.copy(), box.copy()
    first[coordinate, 1] = midpoint
    second[coordinate, 0] = midpoint
    return first, second


def choose_prior_split(
    box: PriorBox,
    conditions: tuple[BranchCondition, ...],
    node: dict[str, Any],
) -> int:
    widths = box[:, 1] - box[:, 0]
    relative = widths / np.maximum(box[:, 0], 1e-15)
    exact = np.asarray(node.get("exact_condition_values", []), dtype=float)
    if exact.size and float(exact.max()) > 1e-8:
        condition = conditions[int(np.argmax(exact))]
        behavior = np.asarray(node["behavior"], dtype=float)
        contributions = np.zeros(4, dtype=float)
        for z, coefficient in enumerate(condition.coefficients):
            g = behavior[z, condition.column]
            lower, upper = box[z]
            if coefficient >= 0.0:
                contributions[z] = abs(coefficient * g * (1.0 / lower - 1.0 / upper))
            else:
                contributions[z] = abs(coefficient * g * (1.0 / lower - 1.0 / upper))
        if contributions.max() > 0.0:
            return int(np.argmax(contributions))
    return int(np.argmax(relative))


def write_payload(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lambda", dest="weight", type=float, default=0.6)
    parser.add_argument(
        "--fixed-three-povm-weights", type=float, nargs=3, required=True
    )
    parser.add_argument("--prefix-order", type=int, nargs=4, required=True)
    parser.add_argument("--target", type=float, default=0.76591)
    parser.add_argument("--max-nodes", type=int, default=501)
    parser.add_argument("--max-conditions", type=int, default=32)
    parser.add_argument("--safety", type=float, default=2e-7)
    parser.add_argument("--tighten-steps", type=int, default=10)
    parser.add_argument("--branch-tighten-steps", type=int, default=6)
    parser.add_argument("--robust-budget", type=float, default=1000.0)
    parser.add_argument("--witness-tolerance", type=float, default=2e-8)
    parser.add_argument("--exhaustive-pairs", action="store_true")
    parser.add_argument("--exhaustive-triples", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    order = tuple(int(x) for x in args.prefix_order)
    if sorted(order) != list(OUTCOMES):
        raise ValueError("prefix order must be a permutation of 0,1,2,3")
    terminal = canonical_three_effect_povm(
        np.asarray(args.fixed_three_povm_weights, dtype=float)
    )
    oracle = FullBehaviorOuter(
        terminal, args.weight, order, args.max_conditions
    )
    initial = initial_prior_box(order)
    initial, tightening_solves = tighten_prior_box(
        oracle, initial, (), args.target, args.safety, args.tighten_steps
    )
    # The bisection returns excluded-side endpoints.  Keeping them makes the
    # box outward/conservative.  A target superlevel point cannot have a zero
    # prior once every recorded lower endpoint is strictly positive.
    if float(initial[:, 0].min()) <= 0.0:
        raise RuntimeError("superlevel tightening did not prove positive priors")

    root = oracle.solve(initial, (), args.safety)
    queue: list[
        tuple[float, int, PriorBox, tuple[BranchCondition, ...], dict[str, Any]]
    ] = []
    counter = 0
    heapq.heappush(queue, (-float(root["bound"]), counter, initial, (), root))
    leaves: list[dict[str, Any]] = []
    processed = 0
    witnesses = 0
    prior_splits = 0
    no_witness_nodes = 0

    while queue and processed < args.max_nodes:
        _, _, box, conditions, node = heapq.heappop(queue)
        processed += 1
        bound = float(node["bound"])
        if bound < args.target:
            leaves.append(
                {
                    "reason": "score",
                    "bound": bound,
                    "box": box.tolist(),
                    "condition_count": len(conditions),
                }
            )
            continue

        exact_values = np.asarray(node.get("exact_condition_values", []), dtype=float)
        if exact_values.size and float(exact_values.max()) > args.witness_tolerance:
            tightened, extra_solves = tighten_prior_box(
                oracle,
                box,
                conditions,
                args.target,
                args.safety,
                args.branch_tighten_steps,
            )
            tightening_solves += extra_solves
            old_width = np.maximum(box[:, 1] - box[:, 0], 1e-300)
            new_width = np.maximum(tightened[:, 1] - tightened[:, 0], 0.0)
            if float(np.prod(new_width / old_width)) < 0.90:
                tightened_node = oracle.solve(tightened, conditions, args.safety)
                counter += 1
                heapq.heappush(
                    queue,
                    (
                        -float(tightened_node["bound"]),
                        counter,
                        tightened,
                        conditions,
                        tightened_node,
                    ),
                )
                continue
            coordinate = choose_prior_split(box, conditions, node)
            preferred = float(np.asarray(node["prefix"])[coordinate])
            children = split_prior_box(box, coordinate, preferred)
            prior_splits += 1
            for child in children:
                child_node = oracle.solve(child, conditions, args.safety)
                counter += 1
                heapq.heappush(
                    queue,
                    (-float(child_node["bound"]), counter, child, conditions, child_node),
                )
            continue

        behavior = np.asarray(node["behavior"], dtype=float)
        dual = find_small_witness(
            behavior,
            args.robust_budget,
            args.witness_tolerance,
            args.exhaustive_pairs,
            args.exhaustive_triples,
        )
        if dual is None:
            primal = solve_primal(behavior)
            no_witness_nodes += 1
            leaves.append(
                {
                    "reason": "no_small_witness",
                    "bound": bound,
                    "box": box.tolist(),
                    "conditions": [condition.to_json() for condition in conditions],
                    "primal_status": primal["status"],
                    "behavior": behavior.tolist(),
                }
            )
            continue

        witnesses += 1
        child_conditions: list[BranchCondition] = []
        for column, coefficient in zip(
            dual["active_columns"],
            dual["halfspace_linear_coefficients"],
            strict=True,
        ):
            condition = BranchCondition(
                int(column), tuple(float(x) for x in coefficient)
            )
            keys = {condition_key(existing) for existing in conditions}
            if condition_key(condition) in keys:
                # Repeating an already exact-satisfied disjunct cannot refine
                # this branch.  The fallback prior split preserves coverage.
                continue
            child_conditions.append(condition)

        if not child_conditions:
            coordinate = choose_prior_split(box, conditions, node)
            prior_splits += 1
            for child in split_prior_box(box, coordinate):
                child_node = oracle.solve(child, conditions, args.safety)
                counter += 1
                heapq.heappush(
                    queue,
                    (-float(child_node["bound"]), counter, child, conditions, child_node),
                )
            continue

        for condition in child_conditions:
            next_conditions = conditions + (condition,)
            child_node = oracle.solve(box, next_conditions, args.safety)
            counter += 1
            heapq.heappush(
                queue,
                (-float(child_node["bound"]), counter, box.copy(), next_conditions, child_node),
            )

        if processed % 25 == 0:
            print(
                {
                    "processed": processed,
                    "open": len(queue),
                    "leaves": len(leaves),
                    "max_open_bound": -queue[0][0] if queue else -math.inf,
                    "witnesses": witnesses,
                    "prior_splits": prior_splits,
                },
                flush=True,
            )

    open_nodes = [
        {
            "bound": -priority,
            "box": box.tolist(),
            "conditions": [condition.to_json() for condition in conditions],
            "node": {
                key: value
                for key, value in node.items()
                if key not in {"statistics", "behavior"}
            },
        }
        for priority, _, box, conditions, node in sorted(queue)
    ]
    payload = {
        "scope": "fixed ternary terminal POVM and one prefix-prior order",
        "weight": args.weight,
        "terminal_effect_weights": list(args.fixed_three_povm_weights),
        "prefix_order": list(order),
        "target": args.target,
        "safety": args.safety,
        "solver": "CLARABEL",
        "initial_prior_box": initial.tolist(),
        "tightening_solves": tightening_solves,
        "processed_nodes": processed,
        "witness_nodes": witnesses,
        "prior_splits": prior_splits,
        "no_witness_nodes": no_witness_nodes,
        "closed": len(open_nodes) == 0 and no_witness_nodes == 0,
        "max_open_bound": max(
            [float(item["bound"]) for item in open_nodes], default=-math.inf
        ),
        "leaf_count": len(leaves),
        "leaves": leaves,
        "open_nodes": open_nodes,
        "certificate_scope_note": (
            "Numerical conic cover only. Formal use requires outward-rounded "
            "SOCP dual bounds and exact/interval validation of PSD witnesses."
        ),
    }
    write_payload(args.output, payload)
    print(json.dumps({key: payload[key] for key in payload if key not in {"leaves", "open_nodes"}}, indent=2))


if __name__ == "__main__":
    main()

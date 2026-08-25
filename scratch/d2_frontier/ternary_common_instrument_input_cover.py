"""Robust shared-Choi cover of a high-score ternary input region.

For one input box choose anchor states ``rho_z^0`` and an anchor terminal POVM
``E_t^0``.  Four positive Choi matrices ``J_y`` describe one and the same
flagged quantum instrument.  Its anchor probabilities are

    q0[z,y,t] = Tr[J_y (rho_z^0^T tensor E_t^0)].

If ``d_inf[z]`` bounds ``||rho_z-rho_z^0||_inf``, ``e_inf[t]`` bounds
``||E_t-E_t^0||_inf``, and ``u[t]`` bounds ``||E_t||_inf``, then every
physical point in the box obeys

    |q[z,y,t]-q0[z,y,t]|
      <= (d_inf[z] u[t] + ||rho_z^0||_inf e_inf[t]) Tr(J_y).

The complete measured row also obeys the stronger channel-contraction tube

    ||q[z]-q0[z]||_1
      <= d_1[z] + ||rho_z^0||_1 sum_t e_inf[t].

Both inequalities are rigorous outer envelopes.  At zero box width they
reduce to a literal common-instrument SDP, not a collection of independent
effective effects.  This driver localises all relaxed inputs capable of
reaching the target and then branch-and-bounds that region.
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

from ternary_common_povm_input_cover import (
    _compact_result,
    _configuration,
    _node_payload,
    _oracle_keywords,
    _set_common_povm_box,
    _split_coordinate,
    _write_checkpoint,
    localise_candidate_region,
)
from ternary_probability_cone_cover import (
    TernaryConeOracle,
    choi_probability_coefficients,
    terminal_weight_intervals,
)
from terminal_reconstruction_enclosure import (
    terminal_effect_anchor_and_errors,
)


CoefficientParameters = tuple[
    tuple[tuple[cp.Parameter, cp.Parameter], ...], ...
]


def instrument_tube_data(
    lower: np.ndarray,
    upper: np.ndarray,
    terminal_effects: np.ndarray,
    terminal_errors: np.ndarray,
    terminal_norm_upper: np.ndarray,
) -> dict[str, np.ndarray]:
    """Return safe anchor coefficients and robust tube radii for one box."""

    lower = np.asarray(lower, dtype=float)
    upper = np.asarray(upper, dtype=float)
    if lower.shape != (4, 4) or upper.shape != (4, 4):
        raise ValueError("input bounds must have shape (4,4)")
    if np.any(lower > upper):
        raise ValueError("input bounds are reversed")
    anchor = 0.5 * (lower + upper)
    radii = 0.5 * (upper - lower)
    vector_radii = np.linalg.norm(radii[:, 1:], axis=1)
    input_operator_radii = 0.5 * (radii[:, 0] + vector_radii)
    input_trace_radii = np.maximum(radii[:, 0], vector_radii)
    anchor_vector_norms = np.linalg.norm(anchor[:, 1:], axis=1)
    anchor_operator_norms = 0.5 * (
        np.abs(anchor[:, 0]) + anchor_vector_norms
    )
    anchor_trace_norms = np.maximum(
        np.abs(anchor[:, 0]), anchor_vector_norms
    )
    probability_radii = (
        input_operator_radii[:, None] * terminal_norm_upper[None, :]
        + anchor_operator_norms[:, None] * terminal_errors[None, :]
    )
    row_radii = (
        input_trace_radii
        + anchor_trace_norms * float(np.sum(terminal_errors))
    )
    return {
        "anchor": anchor,
        "coordinate_radii": radii,
        "coefficients": choi_probability_coefficients(
            anchor, terminal_effects
        ),
        "probability_radii": probability_radii,
        "row_radii": row_radii,
        "input_operator_radii": input_operator_radii,
        "input_trace_radii": input_trace_radii,
        "anchor_operator_norms": anchor_operator_norms,
        "anchor_trace_norms": anchor_trace_norms,
    }


def robust_witness_error(
    witness: np.ndarray,
    tube: dict[str, np.ndarray],
    terminal_effects: np.ndarray,
    terminal_errors: np.ndarray,
) -> tuple[float, list[dict[str, float]]]:
    """Bound one measured-instrument witness over the complete input box."""

    value = np.asarray(witness, dtype=float)
    if value.shape != (4, 4, 3):
        raise ValueError("a measured-instrument witness must have shape (4,4,3)")
    effect_matrices = [
        np.asarray(
            [
                [row[0] + row[3], row[1] - 1j * row[2]],
                [row[1] + 1j * row[2], row[0] - row[3]],
            ]
        )
        for row in terminal_effects
    ]
    rows: list[dict[str, float]] = []
    total = 0.0
    for z in range(4):
        anchor_observable_norm = max(
            float(
                np.linalg.norm(
                    sum(
                        value[z, y, t] * effect_matrices[t]
                        for t in range(3)
                    ),
                    ord=2,
                )
            )
            for y in range(4)
        )
        terminal_coefficient = max(
            float(
                sum(
                    abs(value[z, y, t]) * terminal_errors[t]
                    for t in range(3)
                )
            )
            for y in range(4)
        )
        input_part = float(tube["input_trace_radii"][z]) * (
            anchor_observable_norm + terminal_coefficient
        )
        terminal_part = float(tube["anchor_trace_norms"][z]) * (
            terminal_coefficient
        )
        row_error = input_part + terminal_part
        total += row_error
        rows.append(
            {
                "input_trace_radius": float(tube["input_trace_radii"][z]),
                "anchor_observable_norm": anchor_observable_norm,
                "terminal_coefficient": terminal_coefficient,
                "anchor_trace_norm": float(tube["anchor_trace_norms"][z]),
                "error": row_error,
            }
        )
    return total, rows


class MeasuredInstrumentProjectionOracle:
    """Project one measured table onto a literal shared-instrument image."""

    def __init__(self, coefficients: CoefficientParameters) -> None:
        self.candidate = cp.Parameter(48)
        self.coefficients = coefficients
        self.choi = [
            (
                cp.Variable((4, 4), symmetric=True),
                cp.Variable((4, 4)),
            )
            for _ in range(4)
        ]
        constraints: list[cp.Constraint] = []
        for real, imaginary in self.choi:
            constraints.extend(
                (
                    imaginary + imaginary.T == 0.0,
                    cp.bmat([[real, -imaginary], [imaginary, real]]) >> 0,
                )
            )
        for i in range(2):
            for j in range(2):
                constraints.extend(
                    (
                        sum(
                            real[2 * i, 2 * j]
                            + real[2 * i + 1, 2 * j + 1]
                            for real, _ in self.choi
                        )
                        == (1.0 if i == j else 0.0),
                        sum(
                            imaginary[2 * i, 2 * j]
                            + imaginary[2 * i + 1, 2 * j + 1]
                            for _, imaginary in self.choi
                        )
                        == 0.0,
                    )
                )
        predictions: list[cp.Expression] = []
        for z in range(4):
            for y in range(4):
                choi_real, choi_imaginary = self.choi[y]
                for t in range(3):
                    coefficient_real, coefficient_imaginary = coefficients[z][t]
                    predictions.append(
                        cp.sum(
                            cp.multiply(coefficient_real, choi_real)
                            - cp.multiply(coefficient_imaginary, choi_imaginary)
                        )
                    )
        self.predictions = cp.hstack(predictions)
        self.projection_problem = cp.Problem(
            cp.Minimize(cp.norm(self.predictions - self.candidate, 2)),
            constraints,
        )
        self.support_coefficients = [
            (cp.Parameter((4, 4)), cp.Parameter((4, 4)))
            for _ in range(4)
        ]
        support_expression = sum(
            cp.sum(
                cp.multiply(self.support_coefficients[y][0], self.choi[y][0])
                - cp.multiply(
                    self.support_coefficients[y][1], self.choi[y][1]
                )
            )
            for y in range(4)
        )
        self.support_problem = cp.Problem(
            cp.Maximize(support_expression), constraints
        )
        if not self.projection_problem.is_dpp() or not self.support_problem.is_dpp():
            raise RuntimeError("measured-instrument projection oracle is not DPP")

    @staticmethod
    def _solve(problem: cp.Problem) -> None:
        problem.solve(
            solver="CLARABEL",
            tol_gap_abs=2e-9,
            tol_gap_rel=2e-9,
            tol_feas=2e-9,
            max_iter=1000,
            warm_start=True,
            ignore_dpp=False,
        )

    def project(
        self,
        statistics: np.ndarray,
        support_safety: float,
        distance_tolerance: float,
    ) -> dict[str, Any]:
        candidate = np.asarray(statistics, dtype=float)
        if candidate.shape != (4, 4, 3):
            raise ValueError("candidate statistics must have shape (4,4,3)")
        self.candidate.value = candidate.reshape(48)
        self._solve(self.projection_problem)
        if self.projection_problem.status not in {
            cp.OPTIMAL,
            cp.OPTIMAL_INACCURATE,
        }:
            return {"status": self.projection_problem.status}
        projected = np.asarray(self.predictions.value, dtype=float)
        residual = candidate.reshape(48) - projected
        distance = float(np.linalg.norm(residual))
        result: dict[str, Any] = {
            "status": self.projection_problem.status,
            "distance": distance,
        }
        if distance <= distance_tolerance:
            result["compatible"] = True
            return result
        witness = residual / distance
        witness_table = witness.reshape(4, 4, 3)
        for y in range(4):
            combined = sum(
                witness_table[z, y, t]
                * (
                    np.asarray(self.coefficients[z][t][0].value)
                    + 1j * np.asarray(self.coefficients[z][t][1].value)
                )
                for z in range(4)
                for t in range(3)
            )
            self.support_coefficients[y][0].value = combined.real
            self.support_coefficients[y][1].value = combined.imag
        self._solve(self.support_problem)
        if self.support_problem.status not in {
            cp.OPTIMAL,
            cp.OPTIMAL_INACCURATE,
        }:
            return {
                "status": self.support_problem.status,
                "distance": distance,
            }
        support = float(self.support_problem.value) + support_safety
        candidate_value = float(witness @ candidate.reshape(48))
        result.update(
            {
                "compatible": False,
                "witness": witness_table,
                "anchor_support": support,
                "candidate_value": candidate_value,
                "anchor_gap": candidate_value - support,
            }
        )
        return result


def _coefficient_parameters() -> CoefficientParameters:
    return tuple(
        tuple(
            (cp.Parameter((4, 4)), cp.Parameter((4, 4)))
            for _ in range(3)
        )
        for _ in range(4)
    )


def _assign_tube(
    coefficients: CoefficientParameters,
    probability_radii: cp.Parameter,
    row_radii: cp.Parameter,
    lower_parameter: cp.Parameter,
    upper_parameter: cp.Parameter,
    lower: np.ndarray,
    upper: np.ndarray,
    terminal_effects: np.ndarray,
    terminal_errors: np.ndarray,
    terminal_norm_upper: np.ndarray,
) -> dict[str, np.ndarray]:
    data = instrument_tube_data(
        lower,
        upper,
        terminal_effects,
        terminal_errors,
        terminal_norm_upper,
    )
    raw = data["coefficients"]
    for z in range(4):
        for t in range(3):
            coefficients[z][t][0].value = raw[z, t].real
            coefficients[z][t][1].value = raw[z, t].imag
    probability_radii.value = data["probability_radii"]
    row_radii.value = data["row_radii"]
    lower_parameter.value = lower
    upper_parameter.value = upper
    return data


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
        "source": "ternary robust shared-Choi input-box cover",
        "source_box": source["box"],
        "base_code": source["base_code"],
        "base_plane": source.get("base_plane"),
        "base_sphere": source.get("base_sphere"),
        "top_spectral_cell": bool(use_top_spectral_cell),
        "localisation": localisation,
        "terminal_effect_enclosure": terminal_audit,
        "bound_safety": float(bound_safety),
        "minimum_width": float(minimum_width),
        "max_nodes": int(max_nodes),
        "solved_nodes": len(records),
        "closed_nodes": int(closed),
        "split_nodes": int(split),
        "pending_nodes": len(pending),
        "unresolved_nodes": len(unresolved),
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
            "cell; one literal shared qubit instrument per robust input box; "
            "numerical SDP bounds remain solver-conditional"
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
    witness_support_safety: float,
    witness_tolerance: float,
    resume: dict[str, Any] | None = None,
) -> dict[str, Any]:
    box, contractions, reconstruction = _configuration(
        source, use_top_spectral_cell
    )
    terminal_effects, terminal_errors, interval_norm_upper, terminal_audit = (
        terminal_effect_anchor_and_errors(
            box["terminal_alpha"], box["terminal_beta"]
        )
    )
    weight_norm_upper = np.asarray(
        [interval[1] for interval in terminal_weight_intervals(box)]
    )
    terminal_norm_upper = np.minimum(
        interval_norm_upper, weight_norm_upper
    )
    terminal_audit["retained_operator_norm_upper"] = terminal_norm_upper.tolist()

    coefficients = _coefficient_parameters()
    probability_radii = cp.Parameter((4, 3), nonneg=True)
    row_radii = cp.Parameter(4, nonneg=True)
    povm_anchor = cp.Parameter((4, 4))
    povm_coordinate_radii = cp.Parameter((4, 4), nonneg=True)
    povm_trace_radii = cp.Parameter(4, nonneg=True)
    lower_parameter = cp.Parameter((4, 4))
    upper_parameter = cp.Parameter((4, 4))
    oracle = TernaryConeOracle(
        **_oracle_keywords(source),
        common_contractions=contractions,
        terminal_reconstruction=reconstruction,
        input_pauli_lower=lower_parameter,
        input_pauli_upper=upper_parameter,
        common_povm_input_anchor=povm_anchor,
        common_povm_input_radii=povm_coordinate_radii,
        common_povm_trace_radii=povm_trace_radii,
        common_instrument_probability_coefficients=coefficients,
        common_instrument_probability_radii=probability_radii,
        common_instrument_row_radii=row_radii,
        max_common_instrument_witnesses=max_witnesses,
    )
    projection_oracle = MeasuredInstrumentProjectionOracle(coefficients)

    target = float(source["target"])
    if resume is None:
        lower = np.asarray(localisation["lower"], dtype=float)
        upper = np.asarray(localisation["upper"], dtype=float)
        root = _node_payload(0, None, 0, lower, upper)
        root["parent_bound"] = math.inf
        root["instrument_witnesses"] = []
        pending: list[tuple[float, int, dict[str, Any]]] = [
            (-math.inf, 0, root)
        ]
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
        tube = _assign_tube(
            coefficients,
            probability_radii,
            row_radii,
            lower_parameter,
            upper_parameter,
            lower,
            upper,
            terminal_effects,
            terminal_errors,
            terminal_norm_upper,
        )
        _set_common_povm_box(
            lower_parameter,
            upper_parameter,
            povm_anchor,
            povm_coordinate_radii,
            povm_trace_radii,
            lower,
            upper,
        )
        witness_records = list(node.get("instrument_witnesses", []))
        oracle_solves = 0
        projection_solves = 0
        new_witnesses = 0
        last_projection: dict[str, Any] | None = None
        witness_error_rows: list[dict[str, float]] | None = None
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
            projection = projection_oracle.project(
                np.asarray(result["statistics"], dtype=float),
                witness_support_safety,
                witness_tolerance,
            )
            projection_solves += 1
            last_projection = {
                key: value
                for key, value in projection.items()
                if key != "witness"
            }
            if projection.get("compatible", True):
                break
            witness = np.asarray(projection["witness"], dtype=float)
            robust_error, witness_error_rows = robust_witness_error(
                witness, tube, terminal_effects, terminal_errors
            )
            robust_bound = float(projection["anchor_support"]) + robust_error
            robust_violation = float(projection["candidate_value"]) - robust_bound
            last_projection.update(
                {
                    "robust_error": robust_error,
                    "robust_bound": robust_bound,
                    "robust_violation": robust_violation,
                    "error_rows": witness_error_rows,
                }
            )
            if robust_violation <= witness_tolerance:
                break
            witness_records.append(
                {
                    "coefficients": witness.tolist(),
                    "bound": robust_bound,
                    "anchor_support": float(projection["anchor_support"]),
                    "robust_error": robust_error,
                    "source_violation": robust_violation,
                }
            )
            new_witnesses += 1
        bound = float(result["bound"])
        record = {
            **node,
            "instrument_witnesses": witness_records,
            **_compact_result(result),
            "maximum_probability_radius": float(
                np.max(tube["probability_radii"])
            ),
            "maximum_row_radius": float(np.max(tube["row_radii"])),
            "oracle_solves": oracle_solves,
            "projection_solves": projection_solves,
            "new_witnesses": new_witnesses,
            "last_projection": last_projection,
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
        if result["status"] in {cp.INFEASIBLE, cp.INFEASIBLE_INACCURATE} or bound < target:
            record["disposition"] = "closed"
        elif result["status"] not in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}:
            record["disposition"] = "solver-unresolved"
            unresolved.append(record)
        else:
            if witness_error_rows is None:
                row, coordinate = _split_coordinate(lower, upper)
            else:
                row = int(
                    np.argmax(
                        [item["error"] for item in witness_error_rows]
                    )
                )
                coordinate = int(np.argmax(upper[row] - lower[row]))
            width = float(upper[row, coordinate] - lower[row, coordinate])
            if width <= minimum_width:
                row, coordinate = _split_coordinate(lower, upper)
                width = float(upper[row, coordinate] - lower[row, coordinate])
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
                    child["instrument_witnesses"] = witness_records
                    heapq.heappush(
                        pending, (-bound, next_identifier, child)
                    )
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
                    "row_radius": record["maximum_row_radius"],
                    "new_witnesses": new_witnesses,
                    "robust_violation": (
                        None
                        if last_projection is None
                        else last_projection.get("robust_violation")
                    ),
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
    parser.add_argument("--max-witnesses", type=int, default=12)
    parser.add_argument("--max-new-witnesses-per-node", type=int, default=2)
    parser.add_argument("--witness-support-safety", type=float, default=2e-7)
    parser.add_argument("--witness-tolerance", type=float, default=2e-7)
    parser.add_argument("--top-spectral-cell", action="store_true")
    args = parser.parse_args()

    source = json.loads(args.frontier_json.read_text(encoding="utf-8"))
    resume: dict[str, Any] | None = None
    if args.resume:
        resume = json.loads(args.output.read_text(encoding="utf-8"))
        localisation = resume["localisation"]
    elif args.localisation_json is not None:
        localised = json.loads(
            args.localisation_json.read_text(encoding="utf-8")
        )
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
        args.witness_support_safety,
        args.witness_tolerance,
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

"""Independent matrix audit of the exact shared-instrument SCIP candidate."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

from common_instrument_exact_scip import _apply_choi
from qubit_discrimination_geometry import discrimination_geometry


OUTCOMES = range(4)
ACTIVE = range(3)
PATHS = tuple((z, y) for z in OUTCOMES for y in OUTCOMES)


def matrix_metrics(
    states: np.ndarray,
    choi: np.ndarray,
    effects: np.ndarray,
    weight: float,
) -> dict[str, Any]:
    outputs = np.asarray(
        [[_apply_choi(choi[y], states[z]) for y in OUTCOMES] for z in OUTCOMES]
    )
    probability = np.trace(outputs, axis1=2, axis2=3).real
    statistics = np.einsum("tij,zyji->zyt", effects, outputs).real
    audit = float(
        sum(
            statistics[z, y, z ^ y]
            for z, y in PATHS
            if (z ^ y) in ACTIVE
        )
    )
    returned = float(np.sqrt(np.maximum(probability, 0.0)).sum() ** 2 / 16.0)
    terminal_states = np.asarray(
        [
            sum(outputs[z, z ^ syndrome] for z in OUTCOMES)
            for syndrome in OUTCOMES
        ]
    )
    geometry = discrimination_geometry(terminal_states)
    return {
        "score": weight * audit + (1.0 - weight) * returned,
        "audit": audit,
        "return": returned,
        "normalisation": float(probability.sum()),
        "statistics_normalisation": float(statistics.sum()),
        "path_probability_consistency": float(
            np.max(np.abs(statistics.sum(axis=2) - probability))
        ),
        "minimum_state_eigenvalue": float(
            min(np.linalg.eigvalsh(item).min() for item in states)
        ),
        "minimum_choi_eigenvalue": float(
            min(np.linalg.eigvalsh(item).min() for item in choi)
        ),
        "minimum_output_eigenvalue": float(
            min(
                np.linalg.eigvalsh(item).min()
                for row in outputs
                for item in row
            )
        ),
        "minimum_effect_eigenvalue": float(
            min(np.linalg.eigvalsh(item).min() for item in effects)
        ),
        "trace_preservation_residual": float(
            np.linalg.norm(
                choi.sum(axis=0)
                .reshape(2, 2, 2, 2)
                .trace(axis1=1, axis2=3)
                - np.eye(2)
            )
        ),
        "povm_completeness_residual": float(
            np.linalg.norm(effects.sum(axis=0) - np.eye(2))
        ),
        "optimal_terminal_audit": float(geometry["optimal_guess_probability"]),
        "fixed_povm_optimality_gap": float(
            geometry["optimal_guess_probability"] - audit
        ),
        "optimal_terminal_effect_traces": geometry["optimal_effect_traces"],
    }


def depolarising_repair(
    states: np.ndarray, safety: float
) -> tuple[np.ndarray, float]:
    """Move slightly into the state cone without changing prefix priors."""

    traces = np.trace(states, axis1=1, axis2=2).real
    minima = np.asarray([np.linalg.eigvalsh(item).min() for item in states])
    epsilon = max(
        [
            0.0,
            *(
                (safety - minimum) / (trace / 2.0 - minimum)
                for minimum, trace in zip(minima, traces, strict=True)
                if minimum < safety
            ),
        ]
    )
    epsilon = min(1.0, epsilon * (1.0 + 1e-6) + 1e-15)
    repaired = (1.0 - epsilon) * states + epsilon * traces[
        :, None, None
    ] * np.eye(2) / 2.0
    return repaired, float(epsilon)


def validate(
    result_path: Path,
    checkpoint_path: Path,
    repaired_output: Path | None,
) -> dict[str, Any]:
    result = json.loads(result_path.read_text(encoding="utf-8"))
    if float(result["weight"]) != 0.55:
        raise RuntimeError("unexpected support weight")
    if result["terminal_effect_weights"] != [0.92, 0.64, 0.44, 0.0]:
        raise RuntimeError("unexpected terminal POVM")
    if result["prefix_order"] != [0, 1, 2, 3]:
        raise RuntimeError("unexpected prefix order")
    if not bool(result["gauge_fix"]):
        raise RuntimeError("canonical exact run must fix the input gauge")

    arrays = np.load(checkpoint_path)
    states = np.asarray(arrays["states"], dtype=complex)
    choi = np.asarray(arrays["choi"], dtype=complex)
    effects = np.asarray(arrays["effects"], dtype=complex)
    if states.shape != (4, 2, 2) or choi.shape != (4, 4, 4):
        raise RuntimeError("invalid exact checkpoint dimensions")
    raw = matrix_metrics(states, choi, effects, 0.55)
    reported = result.get("solution")
    if reported is None:
        raise RuntimeError("canonical exact run contains no feasible point")
    if abs(float(reported["score"]) - float(raw["score"])) > 2e-7:
        raise RuntimeError("independent score disagrees with SCIP report")
    if float(raw["trace_preservation_residual"]) > 1e-10:
        raise RuntimeError("instrument is not trace preserving")
    if float(raw["minimum_choi_eigenvalue"]) < -1e-10:
        raise RuntimeError("a Choi matrix is not positive")
    if float(raw["povm_completeness_residual"]) > 1e-10:
        raise RuntimeError("terminal effects do not form a POVM")

    repaired_states, epsilon = depolarising_repair(states, 1e-12)
    repaired = matrix_metrics(repaired_states, choi, effects, 0.55)
    if float(repaired["minimum_state_eigenvalue"]) < 0.0:
        raise RuntimeError("state-cone repair failed")
    if float(repaired["minimum_output_eigenvalue"]) < -1e-10:
        raise RuntimeError("repaired common instrument has a negative output")
    if repaired_output is not None:
        repaired_output.parent.mkdir(parents=True, exist_ok=True)
        np.savez(
            repaired_output,
            states=repaired_states,
            choi=choi,
            effects=effects,
        )

    return {
        "support_weight": 0.55,
        "terminal_effect_weights": [0.92, 0.64, 0.44, 0.0],
        "prefix_order": [0, 1, 2, 3],
        "exact_formulation": result["formulation"],
        "scip_status": result["status"],
        "scip_primal_bound": float(result["primal_bound"]),
        "scip_dual_bound": float(result["dual_bound"]),
        "scip_gap": float(result["gap"]),
        "scip_nodes": int(result["nodes"]),
        "raw_matrix_audit": raw,
        "state_depolarising_repair": epsilon,
        "repaired_physical_strategy": repaired,
        "previous_checkpoint_score": 0.7225888452252243,
        "improvement_over_previous_checkpoint": float(repaired["score"])
        - 0.7225888452252243,
        "status": (
            "one literal shared instrument verified after an inward state-cone "
            "repair; the global SCIP upper bound remains solver-conditional "
            "and does not close the 0.758 target"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--result",
        type=Path,
        default=Path(__file__).with_name(
            "common_instrument_exact_scip_0123_gauge_seed120s.json"
        ),
    )
    parser.add_argument(
        "--checkpoint",
        type=Path,
        default=Path(__file__).with_name(
            "common_instrument_exact_scip_0123_gauge_seed120s.npz"
        ),
    )
    parser.add_argument("--repaired-output", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = validate(args.result, args.checkpoint, args.repaired_output)
    rendered = json.dumps(payload, indent=2) + "\n"
    print(rendered, end="")
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

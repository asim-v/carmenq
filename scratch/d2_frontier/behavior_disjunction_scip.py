"""Global outer optimisation with qubit-behaviour obstruction cuts.

For a fixed ternary terminal POVM, let ``a[z]`` be the four prefix priors and
let ``B[z,j]`` be the conditional probability of the twelve pulled terminal
outcomes ``j=(y,t)``.  The raw statistic is linked exactly by

    g[z,j] = a[z] B[z,j].

Every common qubit state/effect realisation of ``B`` admits a Bloch ellipsoid
nested between its row tetrahedron and the probability simplex.  A numerical
dual obstruction supported on columns ``j_1,...,j_k`` therefore supplies the
valid disjunction

    c_j . B[:,j] <= 0  for at least one supported column j.

This script inserts small-support obstructions as exact mixed-integer
disjunctions, avoiding reciprocal-prior boxes inside the cut itself.  SCIP
then globally optimises the remaining bilinear model.  The loop is a research
certificate generator: conic witnesses and SCIP bounds remain numerical
until rounded and independently validated.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import subprocess
import sys
from typing import Any

import numpy as np
from pyscipopt import Model, quicksum

OUTCOMES = range(4)
ACTIVE = range(3)
PATHS = tuple((z, y) for z in OUTCOMES for y in OUTCOMES)
PAULIS = (
    np.eye(2, dtype=complex),
    np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex),
    np.array([[0.0, -1j], [1j, 0.0]], dtype=complex),
    np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex),
)


def filled_effect_weights(maximum: float) -> np.ndarray:
    remaining = 2.0
    values = []
    for _ in OUTCOMES:
        value = min(maximum, remaining)
        values.append(value)
        remaining -= value
    if abs(remaining) > 1e-10:
        raise ValueError("effect-norm cap is too small to complete a qubit POVM")
    return np.asarray(values)


def canonical_three_effect_povm(weights: np.ndarray) -> np.ndarray:
    w0, w1, w2 = (float(value) for value in weights)
    cosine = (w2 * w2 - w0 * w0 - w1 * w1) / (2.0 * w0 * w1)
    cosine = float(np.clip(cosine, -1.0, 1.0))
    sine = math.sqrt(max(0.0, 1.0 - cosine * cosine))
    directions = np.asarray(
        [
            [1.0, 0.0, 0.0],
            [cosine, sine, 0.0],
            [-(w0 + w1 * cosine) / w2, -(w1 * sine) / w2, 0.0],
        ]
    )
    effects = np.zeros((4, 2, 2), dtype=complex)
    identity = np.eye(2, dtype=complex)
    for s in ACTIVE:
        effects[s] = 0.5 * weights[s] * (
            identity
            + sum(directions[s, axis] * PAULIS[axis + 1] for axis in range(3))
        )
    return effects


def reconstruction_inverse(effects: np.ndarray) -> np.ndarray:
    traces = np.trace(effects[:3], axis1=1, axis2=2).real
    directions = np.asarray(
        [
            [
                float(np.trace((effects[t] / traces[t]) @ pauli).real)
                for pauli in PAULIS[1:3]
            ]
            for t in ACTIVE
        ]
    )
    matrix = np.asarray(
        [
            [
                0.5 * traces[t],
                0.5 * traces[t] * directions[t, 0],
                0.5 * traces[t] * directions[t, 1],
            ]
            for t in ACTIVE
        ]
    )
    return np.linalg.inv(matrix)


def add_lorentz(model: Model, scalar: object, vector: list[object]) -> None:
    model.addCons(scalar >= 0.0)
    model.addCons(scalar * scalar >= quicksum(item * item for item in vector))


def normalise_cut(cut: dict[str, Any]) -> dict[str, Any]:
    coefficients = np.asarray(cut["coefficients"], dtype=float)
    scale = float(np.max(np.abs(coefficients)))
    if not np.isfinite(scale) or scale <= 0.0:
        raise ValueError("a cut coefficient vector is zero or non-finite")
    return {
        "column": int(cut["column"]),
        "coefficients": (coefficients / scale).tolist(),
    }


def cut_key(cut: dict[str, Any]) -> tuple[int, tuple[int, ...]]:
    item = normalise_cut(cut)
    return int(item["column"]), tuple(
        int(round(float(value) * 1e9)) for value in item["coefficients"]
    )


def find_witness(
    behavior: np.ndarray,
    robust_budget: float,
    tolerance: float,
) -> tuple[dict[str, Any] | None, dict[str, Any]]:
    oracle = Path(__file__).with_name("behavior_witness_oracle.py")
    request = json.dumps(
        {
            "behavior": behavior.tolist(),
            "robust_budget": robust_budget,
            "tolerance": tolerance,
        }
    )
    completed = subprocess.run(
        [sys.executable, str(oracle)],
        input=request,
        text=True,
        capture_output=True,
        check=True,
    )
    payload = json.loads(completed.stdout)
    witness = payload["witness"]
    if witness is not None and len(witness["cuts"]) < 2:
        raise ValueError("a behaviour obstruction needs at least two columns")
    return witness, payload["report"]


def prior_box_from_payload(
    path: Path | None, explicit: list[float] | None = None
) -> np.ndarray:
    if path is not None and explicit is not None:
        raise ValueError("use either --prior-box-source or --prior-box")
    if explicit is not None:
        box = np.asarray(explicit, dtype=float).reshape(4, 2)
        if np.any(~np.isfinite(box)) or np.any(box[:, 0] > box[:, 1]):
            raise ValueError("invalid explicit prior box")
        return box
    if path is None:
        return np.asarray([[0.0, 1.0]] * 4, dtype=float)
    payload = json.loads(path.read_text(encoding="utf-8"))
    box = np.asarray(payload["initial_prior_box"], dtype=float)
    if box.shape != (4, 2):
        raise ValueError("initial_prior_box must have shape (4,2)")
    return box


def build_model(
    weights: np.ndarray,
    support_weight: float,
    prefix_order: tuple[int, int, int, int],
    prior_box: np.ndarray,
    witnesses: list[dict[str, Any]],
    target: float | None,
) -> tuple[Model, dict[str, Any]]:
    effects = canonical_three_effect_povm(weights)
    inverse = reconstruction_inverse(effects)
    traces = np.trace(effects[:3], axis1=1, axis2=2).real
    model = Model("behaviour-disjunction-frontier")

    prior = [
        model.addVar(
            lb=float(prior_box[z, 0]),
            ub=float(prior_box[z, 1]),
            name=f"a_{z}",
        )
        for z in OUTCOMES
    ]
    model.addCons(quicksum(prior) == 1.0)
    for rank in range(3):
        model.addCons(
            prior[prefix_order[rank]] >= prior[prefix_order[rank + 1]]
        )

    behavior: dict[tuple[int, int, int], object] = {}
    raw: dict[tuple[int, int, int], object] = {}
    for z in OUTCOMES:
        for y in OUTCOMES:
            for t in ACTIVE:
                behavior[z, y, t] = model.addVar(
                    lb=0.0, ub=float(traces[t]), name=f"b_{z}_{y}_{t}"
                )
                raw[z, y, t] = model.addVar(
                    lb=0.0, ub=float(prior_box[z, 1] * traces[t]),
                    name=f"g_{z}_{y}_{t}",
                )
                model.addCons(
                    raw[z, y, t] == prior[z] * behavior[z, y, t]
                )
        model.addCons(
            quicksum(behavior[z, y, t] for y in OUTCOMES for t in ACTIVE)
            == 1.0
        )

    conditional_path: dict[tuple[int, int], object] = {}
    probability: dict[tuple[int, int], object] = {}
    for z, y in PATHS:
        conditional_path[z, y] = quicksum(
            behavior[z, y, t] for t in ACTIVE
        )
        probability[z, y] = quicksum(raw[z, y, t] for t in ACTIVE)
        for t in ACTIVE:
            model.addCons(
                behavior[z, y, t]
                <= float(traces[t]) * conditional_path[z, y]
            )
            model.addCons(
                raw[z, y, t] <= float(traces[t]) * probability[z, y]
            )

    terminal_statistics: dict[tuple[int, int], object] = {}
    for syndrome in OUTCOMES:
        for t in ACTIVE:
            terminal_statistics[syndrome, t] = quicksum(
                raw[z, y, t] for z, y in PATHS if (z ^ y) == syndrome
            )

    terminal_prior: list[object] = []
    terminal_vector: list[list[object]] = []
    for syndrome in OUTCOMES:
        reconstructed = [
            quicksum(
                float(inverse[row, t]) * terminal_statistics[syndrome, t]
                for t in ACTIVE
            )
            for row in range(3)
        ]
        normal = model.addVar(lb=-1.0, ub=1.0, name=f"tau_{syndrome}_z")
        terminal_prior.append(reconstructed[0])
        terminal_vector.append([reconstructed[1], reconstructed[2], normal])
        model.addCons(
            reconstructed[0]
            == quicksum(
                probability[z, z ^ syndrome] for z in OUTCOMES
            )
        )
        add_lorentz(model, reconstructed[0], terminal_vector[-1])

    audit = model.addVar(lb=0.0, ub=1.0, name="audit")
    model.addCons(
        audit
        == quicksum(terminal_statistics[t, t] for t in ACTIVE)
    )
    cap = filled_effect_weights(float(traces.max()))
    model.addCons(
        audit
        <= quicksum(
            float(cap[rank]) * prior[prefix_order[rank]]
            for rank in OUTCOMES
        )
    )
    model.addCons(
        audit
        <= quicksum(float(traces[t]) * terminal_prior[t] for t in ACTIVE)
    )

    dual_vector = [
        model.addVar(lb=-1.0, ub=1.0, name=f"dual_{axis}")
        for axis in range(3)
    ]
    add_lorentz(model, audit, dual_vector)
    for syndrome in OUTCOMES:
        add_lorentz(
            model,
            audit - terminal_prior[syndrome],
            [
                dual_vector[axis] - terminal_vector[syndrome][axis]
                for axis in range(3)
            ],
        )

    flat = [probability[z, y] for z, y in PATHS]
    hellinger = []
    for first in range(16):
        for second in range(first + 1, 16):
            item = model.addVar(lb=0.0, ub=0.5, name=f"h_{first}_{second}")
            model.addCons(item * item <= flat[first] * flat[second])
            hellinger.append(item)
    returned = (1.0 + 2.0 * quicksum(hellinger)) / 16.0

    binaries = []
    for index, witness in enumerate(witnesses):
        cuts = [normalise_cut(cut) for cut in witness["cuts"]]
        if len(cuts) < 2:
            raise ValueError("every witness must contain at least two cuts")
        selectors = [
            model.addVar(vtype="B", name=f"witness_{index}_side_{side}")
            for side in range(len(cuts))
        ]
        binaries.extend(selectors)
        model.addCons(quicksum(selectors) == 1.0)
        expressions = []
        upper_bounds = []
        for cut in cuts:
            column = int(cut["column"])
            y, t = divmod(column, 3)
            coefficients = [float(value) for value in cut["coefficients"]]
            expressions.append(
                quicksum(
                    coefficients[z] * behavior[z, y, t] for z in OUTCOMES
                )
            )
            upper_bounds.append(
                float(traces[t])
                * sum(max(coefficient, 0.0) for coefficient in coefficients)
            )
        for expression, upper_bound, selector in zip(
            expressions, upper_bounds, selectors, strict=True
        ):
            # selector=1 activates this member of the valid union; selector=0
            # leaves only its independently computed box upper bound.
            model.addCons(expression <= upper_bound * (1.0 - selector))

    score = model.addVar(lb=0.0, ub=1.0, name="score")
    model.addCons(
        score <= support_weight * audit + (1.0 - support_weight) * returned
    )
    if target is not None:
        model.addCons(score >= target)
    model.setObjective(score, "maximize")
    return model, {
        "prior": prior,
        "behavior": behavior,
        "raw": raw,
        "probability": probability,
        "audit": audit,
        "return": returned,
        "score": score,
        "binaries": binaries,
    }


def extract_solution(
    model: Model, variables: dict[str, Any]
) -> dict[str, Any] | None:
    solution = model.getBestSol()
    if solution is None:
        return None
    behavior = np.asarray(
        [
            [
                model.getSolVal(solution, variables["behavior"][z, y, t])
                for y in OUTCOMES
                for t in ACTIVE
            ]
            for z in OUTCOMES
        ],
        dtype=float,
    )
    probability = np.asarray(
        [
            [
                model.getSolVal(solution, variables["probability"][z, y])
                for y in OUTCOMES
            ]
            for z in OUTCOMES
        ],
        dtype=float,
    )
    return {
        "score": float(model.getSolVal(solution, variables["score"])),
        "audit": float(model.getSolVal(solution, variables["audit"])),
        "return": float(model.getSolVal(solution, variables["return"])),
        "prior": [
            float(model.getSolVal(solution, variable))
            for variable in variables["prior"]
        ],
        "behavior": behavior.tolist(),
        "path_probabilities": probability.tolist(),
    }


def solve_round(
    weights: np.ndarray,
    support_weight: float,
    order: tuple[int, int, int, int],
    prior_box: np.ndarray,
    witnesses: list[dict[str, Any]],
    target: float | None,
    seconds: float,
    gap: float,
) -> tuple[dict[str, Any], dict[str, Any] | None]:
    model, variables = build_model(
        weights, support_weight, order, prior_box, witnesses, target
    )
    model.setRealParam("limits/time", seconds)
    model.setRealParam("limits/gap", gap)
    model.setRealParam("numerics/feastol", 1e-9)
    model.setRealParam("numerics/dualfeastol", 1e-9)
    model.setIntParam("display/verblevel", 2)
    model.optimize()
    row = {
        "status": str(model.getStatus()),
        "primal_bound": float(model.getPrimalbound()),
        "dual_bound": float(model.getDualbound()),
        "gap": float(model.getGap()),
        "nodes": int(model.getNNodes()),
        "solving_time": float(model.getSolvingTime()),
        "witness_count": len(witnesses),
    }
    return row, extract_solution(model, variables)


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
    parser.add_argument("--prior-box-source", type=Path)
    parser.add_argument(
        "--prior-box",
        type=float,
        nargs=8,
        metavar=("L0", "U0", "L1", "U1", "L2", "U2", "L3", "U3"),
        help="four explicit lower/upper prior intervals",
    )
    parser.add_argument("--target", type=float)
    parser.add_argument("--rounds", type=int, default=4)
    parser.add_argument("--seconds-per-round", type=float, default=60.0)
    parser.add_argument("--gap", type=float, default=1e-5)
    parser.add_argument("--robust-budget", type=float, default=1000.0)
    parser.add_argument("--witness-tolerance", type=float, default=2e-8)
    parser.add_argument(
        "--initial-witness-source",
        type=Path,
        help="reuse the witness list from an earlier JSON run",
    )
    parser.add_argument(
        "--initial-witness-limit",
        type=int,
        help="use only the first N witnesses from the source",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    weights = np.asarray(args.fixed_three_povm_weights, dtype=float)
    order = tuple(int(value) for value in args.prefix_order)
    if sorted(order) != list(OUTCOMES):
        raise ValueError("prefix order must be a permutation of 0,1,2,3")
    box = prior_box_from_payload(args.prior_box_source, args.prior_box)
    if args.initial_witness_source is None:
        witnesses: list[dict[str, Any]] = []
    else:
        source = json.loads(
            args.initial_witness_source.read_text(encoding="utf-8")
        )
        if float(source["weight"]) != args.weight:
            raise ValueError("initial witnesses use a different support weight")
        if source["terminal_effect_weights"] != weights.tolist():
            raise ValueError("initial witnesses use a different terminal POVM")
        if source["prefix_order"] != list(order):
            raise ValueError("initial witnesses use a different prefix order")
        witnesses = list(source["witnesses"])
        if args.initial_witness_limit is not None:
            if args.initial_witness_limit < 0:
                raise ValueError("initial witness limit must be nonnegative")
            witnesses = witnesses[: args.initial_witness_limit]
    rounds = []
    seen: set[tuple[tuple[int, tuple[int, ...]], ...]] = {
        tuple(sorted(cut_key(cut) for cut in witness["cuts"]))
        for witness in witnesses
    }
    for round_index in range(args.rounds):
        row, solution = solve_round(
            weights,
            args.weight,
            order,
            box,
            witnesses,
            args.target,
            args.seconds_per_round,
            args.gap,
        )
        row["round"] = round_index + 1
        row["incumbent"] = solution
        rounds.append(row)
        print(json.dumps({k: v for k, v in row.items() if k != "incumbent"}), flush=True)
        payload = {
            "scope": "fixed ternary terminal POVM, fixed prefix order, and recorded prior box",
            "weight": args.weight,
            "terminal_effect_weights": weights.tolist(),
            "prefix_order": list(order),
            "prior_box": box.tolist(),
            "target": args.target,
            "witnesses": witnesses,
            "rounds": rounds,
            "certificate_scope_note": (
                "Numerical SCIP and conic-witness experiment. Publication use "
                "requires outward bounds and exact PSD validation."
            ),
        }
        write_payload(args.output, payload)
        if solution is None:
            break
        witness, witness_report = find_witness(
            np.asarray(solution["behavior"], dtype=float),
            args.robust_budget,
            args.witness_tolerance,
        )
        row["witness_search"] = witness_report
        if witness is None:
            break
        key = tuple(sorted(cut_key(cut) for cut in witness["cuts"]))
        if key in seen:
            row["witness_search"]["duplicate"] = True
            break
        seen.add(key)
        witnesses.append(witness)

    payload = {
        "scope": "fixed ternary terminal POVM, fixed prefix order, and recorded prior box",
        "weight": args.weight,
        "terminal_effect_weights": weights.tolist(),
        "prefix_order": list(order),
        "prior_box": box.tolist(),
        "target": args.target,
        "witnesses": witnesses,
        "rounds": rounds,
        "certificate_scope_note": (
            "Numerical SCIP and conic-witness experiment. Publication use "
            "requires outward bounds and exact PSD validation."
        ),
    }
    write_payload(args.output, payload)
    print(json.dumps({k: v for k, v in payload.items() if k != "rounds"}, indent=2))


if __name__ == "__main__":
    main()

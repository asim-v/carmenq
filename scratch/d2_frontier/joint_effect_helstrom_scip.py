"""Spatial branch-and-bound for the joint-effect/Helstrom outer model.

All qubit states and effects are represented by Bloch four-vectors.  Their
positive-semidefinite constraints are Lorentz-cone quadratic inequalities;
the only nonconvex equalities are the physical state/effect products that
produce path and AUDIT probabilities.  SCIP therefore supplies a global
dual bound for a *fixed* terminal rank-one POVM, subject to its numerical
tolerances.

This is still an outer model: ``G[y,s] <= w_s Q[y]`` is necessary but not
sufficient for all ``G[y,s]`` to be pullbacks of the same terminal POVM by a
common output instrument.  A successful upper certificate is consequently
valid for the Choi programme; a high value need not be physically attainable.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from pyscipopt import Model, quicksum

OUTCOMES = range(4)
PATHS = tuple((z, y) for z in OUTCOMES for y in OUTCOMES)
IDENTITY = np.eye(2, dtype=complex)
PAULIS = (
    np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex),
    np.array([[0.0, -1j], [1j, 0.0]], dtype=complex),
    np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex),
)


def canonical_three_effect_povm(weights: np.ndarray) -> np.ndarray:
    """Return the planar rank-one qubit POVM fixed by three traces."""
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
    for s in range(3):
        effects[s] = 0.5 * weights[s] * (
            IDENTITY
            + sum(directions[s, axis] * PAULIS[axis] for axis in range(3))
        )
    return effects


def bloch(matrix: np.ndarray) -> np.ndarray:
    return np.asarray(
        [
            float(np.trace(matrix).real),
            *(float(np.trace(matrix @ pauli).real) for pauli in PAULIS),
        ]
    )


def add_lorentz(
    model: Model,
    scalar: object,
    vector: tuple[object, object, object] | list[object],
) -> None:
    model.addCons(scalar >= 0.0)
    model.addCons(scalar * scalar >= quicksum(item * item for item in vector))


def build(
    effects: np.ndarray,
    weight: float,
    prefix_order: tuple[int, int, int, int] | None,
    target: float | None,
    fix_rotation_gauge: bool,
    linked_columns: tuple[str, ...] | None,
) -> tuple[Model, dict[str, object]]:
    traces = np.trace(effects, axis1=1, axis2=2).real
    directions = np.zeros((4, 3), dtype=float)
    for s in OUTCOMES:
        if traces[s] > 1e-12:
            directions[s] = bloch(effects[s] / traces[s])[1:]

    model = Model("joint-effect-helstrom")
    variables: dict[str, object] = {}

    state_scalar = []
    state_vector = []
    for z in OUTCOMES:
        scalar = model.addVar(lb=0.0, ub=1.0, name=f"a_{z}")
        vector = tuple(
            model.addVar(lb=-1.0, ub=1.0, name=f"r_{z}_{axis}")
            for axis in range(3)
        )
        add_lorentz(model, scalar, vector)
        state_scalar.append(scalar)
        state_vector.append(vector)
        variables[f"a_{z}"] = scalar
        for axis, item in enumerate(vector):
            variables[f"r_{z}_{axis}"] = item
    model.addCons(quicksum(state_scalar) == 1.0)
    if fix_rotation_gauge:
        # A common SO(3) rotation of all input Bloch vectors and pulled
        # effects is physically irrelevant.  Rotate r_0 to the positive
        # x-axis and r_1 into the upper xy half-plane.  The convention also
        # covers zero/collinear vectors and removes three continuous gauge
        # directions without excluding a physical orbit.
        model.addCons(state_vector[0][1] == 0.0)
        model.addCons(state_vector[0][2] == 0.0)
        model.addCons(state_vector[0][0] >= 0.0)
        model.addCons(state_vector[1][2] == 0.0)
        model.addCons(state_vector[1][1] >= 0.0)

    effect_scalar: dict[tuple[int, int], object] = {}
    effect_vector: dict[tuple[int, int], tuple[object, object, object]] = {}
    for y, s in PATHS:
        upper = 2.0 * float(traces[s])
        scalar = model.addVar(lb=0.0, ub=upper, name=f"g_{y}_{s}_0")
        vector = tuple(
            model.addVar(lb=-upper, ub=upper, name=f"g_{y}_{s}_{axis + 1}")
            for axis in range(3)
        )
        add_lorentz(model, scalar, vector)
        effect_scalar[y, s] = scalar
        effect_vector[y, s] = vector
        variables[f"g_{y}_{s}_0"] = scalar
        for axis, item in enumerate(vector, start=1):
            variables[f"g_{y}_{s}_{axis}"] = item

    model.addCons(quicksum(effect_scalar.values()) == 2.0)
    for axis in range(3):
        model.addCons(
            quicksum(effect_vector[y, s][axis] for y, s in PATHS) == 0.0
        )

    coarse_scalar = {
        y: quicksum(effect_scalar[y, s] for s in OUTCOMES) for y in OUTCOMES
    }
    coarse_vector = {
        y: tuple(
            quicksum(effect_vector[y, s][axis] for s in OUTCOMES)
            for axis in range(3)
        )
        for y in OUTCOMES
    }
    for y, s in PATHS:
        residual_scalar = traces[s] * coarse_scalar[y] - effect_scalar[y, s]
        residual_vector = tuple(
            traces[s] * coarse_vector[y][axis] - effect_vector[y, s][axis]
            for axis in range(3)
        )
        add_lorentz(model, residual_scalar, residual_vector)

    active = tuple(int(s) for s in OUTCOMES if traces[s] > 1e-9)
    linked = (
        {f"b_{3 * y + t}" for y in OUTCOMES for t in active}
        if linked_columns is None
        else set(linked_columns)
    )
    statistics: dict[tuple[int, int, int], object] = {}
    probability: dict[tuple[int, int], object] = {}
    correct: dict[tuple[int, int], object] = {}
    for z, y in PATHS:
        p = model.addVar(lb=0.0, ub=1.0, name=f"p_{z}_{y}")
        model.addCons(
            2.0 * p
            == state_scalar[z] * coarse_scalar[y]
            + quicksum(
                state_vector[z][axis] * coarse_vector[y][axis]
                for axis in range(3)
            )
        )
        for terminal_label in active:
            item = model.addVar(
                lb=0.0,
                ub=float(traces[terminal_label]),
                name=f"q_{z}_{y}_{terminal_label}",
            )
            if f"b_{3 * y + terminal_label}" in linked:
                model.addCons(
                    2.0 * item
                    == state_scalar[z] * effect_scalar[y, terminal_label]
                    + quicksum(
                        state_vector[z][axis]
                        * effect_vector[y, terminal_label][axis]
                        for axis in range(3)
                    )
                )
            model.addCons(item <= traces[terminal_label] * p)
            statistics[z, y, terminal_label] = item
            variables[f"q_{z}_{y}_{terminal_label}"] = item
        model.addCons(p == quicksum(statistics[z, y, t] for t in active))
        s = z ^ y
        d = statistics[z, y, s] if s in active else 0.0
        probability[z, y] = p
        correct[z, y] = d
        variables[f"p_{z}_{y}"] = p

    for name in sorted(linked):
        parts = name.split("_")
        if parts[0] != "d":
            if parts[0] != "b" or len(parts) != 2:
                raise ValueError(f"invalid linked column {name!r}")
            continue
        if len(parts) != 3:
            raise ValueError(f"invalid linked residual {name!r}")
        y, terminal_label = int(parts[1]), int(parts[2])
        if y not in OUTCOMES or terminal_label not in active:
            raise ValueError(f"invalid linked residual {name!r}")
        residual_scalar = (
            traces[terminal_label] * coarse_scalar[y]
            - effect_scalar[y, terminal_label]
        )
        residual_vector = tuple(
            traces[terminal_label] * coarse_vector[y][axis]
            - effect_vector[y, terminal_label][axis]
            for axis in range(3)
        )
        for z in OUTCOMES:
            model.addCons(
                2.0
                * (
                    traces[terminal_label] * probability[z, y]
                    - statistics[z, y, terminal_label]
                )
                == state_scalar[z] * residual_scalar
                + quicksum(
                    state_vector[z][axis] * residual_vector[axis]
                    for axis in range(3)
                )
            )

    model.addCons(quicksum(probability.values()) == 1.0)
    for z in OUTCOMES:
        # This follows algebraically from POVM completeness.  Stating it
        # explicitly tightens node relaxations before all effect products
        # have been resolved.
        model.addCons(
            quicksum(probability[z, y] for y in OUTCOMES) == state_scalar[z]
        )
    audit = model.addVar(lb=0.0, ub=1.0, name="audit")
    model.addCons(audit == quicksum(correct.values()))
    variables["audit"] = audit

    prefix = {
        z: quicksum(probability[z, y] for y in OUTCOMES) for z in OUTCOMES
    }
    if prefix_order is not None:
        for index in range(3):
            model.addCons(prefix[prefix_order[index]] >= prefix[prefix_order[index + 1]])
        maximum = float(traces.max())
        cap = [maximum, maximum, 2.0 - 2.0 * maximum, 0.0]
        model.addCons(
            audit
            <= quicksum(cap[index] * prefix[prefix_order[index]] for index in OUTCOMES)
        )

    syndrome_from_paths = {
        s: quicksum(probability[z, z ^ s] for z in OUTCOMES) for s in OUTCOMES
    }
    terminal_statistics = {
        (s, t): quicksum(
            statistics[z, y, t] for z, y in PATHS if (z ^ y) == s
        )
        for s in OUTCOMES
        for t in active
    }
    reconstruction = np.asarray(
        [
            [
                0.5 * traces[s],
                0.5 * traces[s] * directions[s, 0],
                0.5 * traces[s] * directions[s, 1],
            ]
            for s in active
        ]
    )
    reconstruction_inverse = np.linalg.inv(reconstruction)
    syndrome = {}
    terminal_vector = []
    for s in OUTCOMES:
        reconstructed = tuple(
            quicksum(
                float(reconstruction_inverse[row, column])
                * terminal_statistics[s, active[column]]
                for column in range(3)
            )
            for row in range(3)
        )
        syndrome[s] = reconstructed[0]
        model.addCons(syndrome[s] == syndrome_from_paths[s])
        normal = model.addVar(lb=-1.0, ub=1.0, name=f"tau_{s}_2")
        vector = (reconstructed[1], reconstructed[2], normal)
        add_lorentz(model, syndrome[s], vector)
        terminal_vector.append(vector)
        variables[f"tau_{s}_2"] = normal

    dual_vector = tuple(
        model.addVar(lb=-1.0, ub=1.0, name=f"dual_{axis}")
        for axis in range(3)
    )
    add_lorentz(model, audit, dual_vector)
    for axis, item in enumerate(dual_vector):
        variables[f"dual_{axis}"] = item
    for s in OUTCOMES:
        add_lorentz(
            model,
            audit - syndrome[s],
            tuple(dual_vector[axis] - terminal_vector[s][axis] for axis in range(3)),
        )
    model.addCons(
        audit <= quicksum(traces[s] * syndrome[s] for s in OUTCOMES)
    )

    flat = [probability[z, y] for z, y in PATHS]
    hellinger = []
    for first in range(16):
        for second in range(first + 1, 16):
            item = model.addVar(lb=0.0, ub=0.5, name=f"h_{first}_{second}")
            model.addCons(item * item <= flat[first] * flat[second])
            hellinger.append(item)
            variables[f"h_{first}_{second}"] = item
    returned = (1.0 + 2.0 * quicksum(hellinger)) / 16.0
    score = model.addVar(lb=0.0, ub=1.0, name="score")
    model.addCons(score <= weight * audit + (1.0 - weight) * returned)
    if target is not None:
        # A proof that this constrained model is infeasible certifies the
        # desired strict upper bound and is usually easier than closing the
        # full optimisation gap.
        model.addCons(score >= target)
    model.setObjective(score, "maximize")
    variables.update(
        {
            "score": score,
            "return": returned,
            "terminal_effect_weights": traces,
            "terminal_directions": directions,
        }
    )
    return model, variables


def seed_solution(
    model: Model,
    variables: dict[str, object],
    checkpoint: Path,
    effects: np.ndarray,
    weight: float,
) -> None:
    arrays = np.load(checkpoint)
    states = np.asarray(arrays["states"])
    joint = np.asarray(arrays["joint"])
    coarse = joint.sum(axis=1)
    probabilities = np.einsum("yij,zji->zy", coarse, states).real
    correct = np.asarray(
        [
            [np.trace(joint[y, z ^ y] @ states[z]).real for y in OUTCOMES]
            for z in OUTCOMES
        ]
    )
    audit = float(correct.sum())
    priors = np.asarray(
        [sum(probabilities[z, z ^ s] for z in OUTCOMES) for s in OUTCOMES]
    )
    if "terminal_vectors" not in arrays or "dual_vector" not in arrays:
        return
    terminal_vectors = np.asarray(arrays["terminal_vectors"])
    dual_vector = np.asarray(arrays["dual_vector"])

    values: dict[str, float] = {}
    for z in OUTCOMES:
        vector = bloch(states[z])
        values[f"a_{z}"] = vector[0]
        for axis in range(3):
            values[f"r_{z}_{axis}"] = vector[axis + 1]
    for y, s in PATHS:
        vector = bloch(joint[y, s])
        for axis in range(4):
            values[f"g_{y}_{s}_{axis}"] = vector[axis]
    for z, y in PATHS:
        values[f"p_{z}_{y}"] = probabilities[z, y]
        values[f"d_{z}_{y}"] = correct[z, y]
    values["audit"] = audit
    for s in OUTCOMES:
        for axis in range(3):
            values[f"tau_{s}_{axis}"] = terminal_vectors[s, axis]
    for axis in range(3):
        values[f"dual_{axis}"] = dual_vector[axis]
    flat = probabilities.reshape(16)
    for first in range(16):
        for second in range(first + 1, 16):
            values[f"h_{first}_{second}"] = math.sqrt(
                max(0.0, flat[first] * flat[second])
            )
    returned = float(np.sqrt(np.maximum(probabilities, 0.0)).sum() ** 2 / 16.0)
    values["score"] = weight * audit + (1.0 - weight) * returned

    solution = model.createSol()
    for name, value in values.items():
        if name in variables:
            model.setSolVal(solution, variables[name], float(value))
    model.addSol(solution)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lambda", dest="weight", type=float, default=0.6)
    parser.add_argument("--fixed-three-povm-weights", type=float, nargs=3, required=True)
    parser.add_argument("--prefix-order", type=int, nargs=4)
    parser.add_argument("--seed-npz", type=Path)
    parser.add_argument("--seconds", type=float, default=300.0)
    parser.add_argument("--gap", type=float, default=1e-5)
    parser.add_argument("--target", type=float)
    parser.add_argument("--no-rotation-gauge", action="store_true")
    parser.add_argument(
        "--linked-column",
        action="append",
        default=[],
        help=(
            "link one b_J or d_y_t statistic to its Born product; "
            "when omitted, all twelve b_J columns are linked"
        ),
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    weights = np.asarray(args.fixed_three_povm_weights, dtype=float)
    effects = canonical_three_effect_povm(weights)
    prefix_order = None if args.prefix_order is None else tuple(args.prefix_order)
    if prefix_order is not None and sorted(prefix_order) != list(OUTCOMES):
        raise ValueError("prefix order must be a permutation of 0,1,2,3")
    model, variables = build(
        effects,
        args.weight,
        prefix_order,
        args.target,
        not args.no_rotation_gauge,
        None if not args.linked_column else tuple(dict.fromkeys(args.linked_column)),
    )
    if args.seed_npz is not None:
        seed_solution(model, variables, args.seed_npz, effects, args.weight)
    model.setRealParam("limits/time", args.seconds)
    model.setRealParam("limits/gap", args.gap)
    model.setRealParam("numerics/feastol", 1e-9)
    model.setRealParam("numerics/dualfeastol", 1e-9)
    model.setIntParam("display/verblevel", 2)
    model.optimize()

    payload = {
        "weight": args.weight,
        "terminal_effect_weights": weights.tolist(),
        "prefix_order": None if prefix_order is None else list(prefix_order),
        "target": args.target,
        "rotation_gauge_fixed": not args.no_rotation_gauge,
        "linked_columns": (
            "all_b" if not args.linked_column else list(dict.fromkeys(args.linked_column))
        ),
        "status": str(model.getStatus()),
        "primal_bound": float(model.getPrimalbound()),
        "dual_bound": float(model.getDualbound()),
        "absolute_gap": float(model.getDualbound() - model.getPrimalbound()),
        "relative_gap": float(model.getGap()),
        "nodes": int(model.getNNodes()),
        "solving_time": float(model.getSolvingTime()),
    }
    solution = model.getBestSol()
    if solution is not None:
        payload["incumbent_score"] = float(
            model.getSolVal(solution, variables["score"])
        )
        payload["incumbent_audit"] = float(
            model.getSolVal(solution, variables["audit"])
        )
        payload["incumbent_return"] = float(
            model.getSolVal(solution, variables["return"])
        )
    rendered = json.dumps(payload, indent=2) + "\n"
    print(rendered)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

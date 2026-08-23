"""Global probability bound with exact pairwise-qubit compatibility.

The pulled-statistics probability relaxation can fabricate a different qubit
geometry for every effect column.  A small obstruction is already visible in
pairs of columns.  For two compatible effects ``E,F`` and four prefix labels,
the conditional probability pairs must be the image of four Bloch-ball points
under the same two effects.  Conversely that condition is exact for the pair.

For each requested pair this script introduces an independent planar gauge

    E = c I + u sigma_x,
    F = d I + v_x sigma_x + v_y sigma_y,

with ``E >= 0``, ``F >= 0``, and ``I-E-F >= 0``.  Four subnormalised Bloch
vectors share the already-present prefix traces.  The selected joint
statistics are linked by bilinear Born equalities.  Different pairs need not
share a realization, so the resulting model is still an outer relaxation;
every physical common-instrument strategy is nevertheless feasible.

The alternative inellipse representation is also necessary: preprocessing
cleanness completes every ternary qubit probability ellipse to a rank-one
inellipse.  Three choices of residual coordinate, each with reciprocal
Horwitz parameters in ``[1,2]``, cover that model.  Different selected pairs
still receive independent realizations, so either representation is an outer
relaxation of the common-instrument problem.

The rest of the model keeps the exact fixed-terminal Helstrom KKT system and
the exact Hellinger hypograph.  SCIP's spatial branch-and-bound therefore
gives a numerical global upper bound subject to its recorded tolerances.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from pyscipopt import Model, quicksum

from joint_effect_helstrom_scip import (
    OUTCOMES,
    PATHS,
    add_lorentz,
    bloch,
    canonical_three_effect_povm,
)


Column = tuple[str, int, int]


def parse_column(name: str) -> Column:
    parts = name.split("_")
    if len(parts) == 2 and parts[0] == "b":
        index = int(parts[1])
        if not 0 <= index < 12:
            raise ValueError(f"invalid effect column {name!r}")
        return "b", index // 3, index % 3
    if len(parts) == 3 and parts[0] == "d":
        y, t = int(parts[1]), int(parts[2])
        if not (0 <= y < 4 and 0 <= t < 3):
            raise ValueError(f"invalid residual column {name!r}")
        return "d", y, t
    raise ValueError(f"invalid column {name!r}")


def render_column(column: Column) -> str:
    kind, y, t = column
    return f"b_{3 * y + t}" if kind == "b" else f"d_{y}_{t}"


def parse_pair(text: str) -> tuple[Column, Column]:
    parts = text.split(",")
    if len(parts) != 2:
        raise ValueError("a pair must have the form b_J,b_K or b_J,d_y_t")
    first, second = (parse_column(part.strip()) for part in parts)
    if first[1] == second[1]:
        raise ValueError(
            "the current sufficient compatibility check expects columns "
            "from distinct coarse outcomes"
        )
    return first, second


def selected_statistic(
    statistics: dict[tuple[int, int, int], object],
    probability: dict[tuple[int, int], object],
    weights: np.ndarray,
    column: Column,
    z: int,
) -> object:
    kind, y, t = column
    if kind == "b":
        return statistics[z, y, t]
    return float(weights[t]) * probability[z, y] - statistics[z, y, t]


def add_pair_factorisation(
    model: Model,
    variables: dict[str, object],
    pair_index: int,
    pair: tuple[Column, Column],
    prefix: dict[int, object],
    statistics: dict[tuple[int, int, int], object],
    probability: dict[tuple[int, int], object],
    weights: np.ndarray,
) -> None:
    prefix_name = f"pair_{pair_index}"
    c = model.addVar(lb=0.0, ub=1.0, name=f"{prefix_name}_c")
    d = model.addVar(lb=0.0, ub=1.0, name=f"{prefix_name}_d")
    ux = model.addVar(lb=0.0, ub=1.0, name=f"{prefix_name}_ux")
    vx = model.addVar(lb=-1.0, ub=1.0, name=f"{prefix_name}_vx")
    vy = model.addVar(lb=0.0, ub=1.0, name=f"{prefix_name}_vy")
    # Once these five effect coordinates are fixed, all Born equalities in
    # this pair become linear and the remaining state constraints are conic.
    # Prioritising them makes spatial branching follow that block structure.
    for variable in (c, d, ux, vx, vy):
        model.chgVarBranchPriority(variable, 10_000)

    # Positivity of E, F, and I-E-F in the planar rotational gauge.
    model.addCons(c * c >= ux * ux)
    model.addCons(d * d >= vx * vx + vy * vy)
    model.addCons(ux <= c)
    model.addCons(vx <= d)
    model.addCons(vx >= -d)
    model.addCons(vy <= d)
    residual = 1.0 - c - d
    model.addCons(residual >= 0.0)
    model.addCons(
        residual * residual >= (ux + vx) * (ux + vx) + vy * vy
    )
    model.addCons(ux + vx <= residual)
    model.addCons(ux + vx >= -residual)
    model.addCons(vy <= residual)
    for name, variable in (
        ("c", c),
        ("d", d),
        ("ux", ux),
        ("vx", vx),
        ("vy", vy),
    ):
        variables[f"{prefix_name}_{name}"] = variable

    first, second = pair
    first_cap = float(weights[first[2]])
    second_cap = float(weights[second[2]])
    for z in OUTCOMES:
        rx = model.addVar(lb=-1.0, ub=1.0, name=f"{prefix_name}_r_{z}_x")
        ry = model.addVar(lb=-1.0, ub=1.0, name=f"{prefix_name}_r_{z}_y")
        model.addCons(prefix[z] * prefix[z] >= rx * rx + ry * ry)
        model.addCons(rx <= prefix[z])
        model.addCons(rx >= -prefix[z])
        model.addCons(ry <= prefix[z])
        model.addCons(ry >= -prefix[z])
        first_value = selected_statistic(
            statistics, probability, weights, first, z
        )
        second_value = selected_statistic(
            statistics, probability, weights, second, z
        )
        # Each selected effect is bounded by w_t Q_y.  Dividing it by
        # w_t gives a subeffect of Q_y; columns from distinct y therefore
        # form two outcomes of one three-effect POVM.  The variables above
        # parameterise these *normalised* effects.
        model.addCons(
            first_value == first_cap * (c * prefix[z] + ux * rx)
        )
        model.addCons(
            second_value
            == second_cap * (d * prefix[z] + vx * rx + vy * ry)
        )
        model.addCons(
            first_value / first_cap + second_value / second_cap <= prefix[z]
        )
        variables[f"{prefix_name}_r_{z}_x"] = rx
        variables[f"{prefix_name}_r_{z}_y"] = ry


def add_pair_inellipse(
    model: Model,
    variables: dict[str, object],
    pair_index: int,
    pair: tuple[Column, Column],
    prefix: dict[int, object],
    statistics: dict[tuple[int, int, int], object],
    probability: dict[tuple[int, int], object],
    weights: np.ndarray,
    reciprocal_cap: float,
    coordinate_case: str,
) -> None:
    """Impose Horwitz's homogenised unit-triangle inellipse condition.

    The clean-POVM completion theorem in
    ``notes/clean_povm_inellipse_completion.md`` proves necessity.  One of
    ``xy``, ``xr``, and ``yr`` must hold after relabelling the ternary
    outcomes; an invocation fixes one member of that finite union.
    """

    name = f"inellipse_{pair_index}"
    # Divide Horwitz's equation by w^2 t^2 and use alpha=1/w,
    # beta=1/t.  This prevents w=t=0 from annihilating the whole conic.
    alpha = model.addVar(lb=1.0, ub=reciprocal_cap, name=f"{name}_alpha")
    beta = model.addVar(lb=1.0, ub=reciprocal_cap, name=f"{name}_beta")
    alpha2 = model.addVar(
        lb=1.0, ub=reciprocal_cap**2, name=f"{name}_alpha2"
    )
    beta2 = model.addVar(
        lb=1.0, ub=reciprocal_cap**2, name=f"{name}_beta2"
    )
    alphabeta = model.addVar(
        lb=1.0, ub=reciprocal_cap**2, name=f"{name}_alphabeta"
    )
    cross = model.addVar(
        lb=-2.0 * reciprocal_cap**2,
        ub=8.0 * reciprocal_cap,
        name=f"{name}_cross",
    )
    model.addCons(alpha2 == alpha * alpha)
    model.addCons(beta2 == beta * beta)
    model.addCons(alphabeta == alpha * beta)
    model.addCons(cross == -2.0 * (alphabeta - 2.0 * alpha - 2.0 * beta + 2.0))
    for variable in (alpha, beta, alpha2, beta2, alphabeta, cross):
        model.chgVarBranchPriority(variable, 10_000)
    variables.update(
        {
            f"{name}_alpha": alpha,
            f"{name}_beta": beta,
            f"{name}_alpha2": alpha2,
            f"{name}_beta2": beta2,
            f"{name}_alphabeta": alphabeta,
            f"{name}_cross": cross,
        }
    )

    first, second = pair
    first_cap = float(weights[first[2]])
    second_cap = float(weights[second[2]])
    for z in OUTCOMES:
        x = model.addVar(lb=0.0, ub=1.0, name=f"{name}_x_{z}")
        y = model.addVar(lb=0.0, ub=1.0, name=f"{name}_y_{z}")
        residual = model.addVar(lb=0.0, ub=1.0, name=f"{name}_r_{z}")
        u = model.addVar(lb=0.0, ub=1.0, name=f"{name}_u_{z}")
        v = model.addVar(lb=0.0, ub=1.0, name=f"{name}_v_{z}")
        u2 = model.addVar(lb=0.0, ub=1.0, name=f"{name}_u2_{z}")
        v2 = model.addVar(lb=0.0, ub=1.0, name=f"{name}_v2_{z}")
        uv = model.addVar(lb=0.0, ub=1.0, name=f"{name}_uv_{z}")
        ua = model.addVar(lb=0.0, ub=1.0, name=f"{name}_ua_{z}")
        va = model.addVar(lb=0.0, ub=1.0, name=f"{name}_va_{z}")
        a2 = model.addVar(lb=0.0, ub=1.0, name=f"{name}_a2_{z}")
        first_value = selected_statistic(
            statistics, probability, weights, first, z
        )
        second_value = selected_statistic(
            statistics, probability, weights, second, z
        )
        model.addCons(first_value == first_cap * x)
        model.addCons(second_value == second_cap * y)
        model.addCons(residual == prefix[z] - x - y)
        if coordinate_case == "xy":
            model.addCons(u == x)
            model.addCons(v == y)
        elif coordinate_case == "xr":
            model.addCons(u == x)
            model.addCons(v == residual)
        elif coordinate_case == "yr":
            model.addCons(u == y)
            model.addCons(v == residual)
        else:
            raise ValueError(f"unknown coordinate case {coordinate_case!r}")
        model.addCons(u2 == u * u)
        model.addCons(v2 == v * v)
        model.addCons(uv == u * v)
        model.addCons(ua == u * prefix[z])
        model.addCons(va == v * prefix[z])
        model.addCons(a2 == prefix[z] * prefix[z])
        model.addCons(
            beta2 * u2
            + alpha2 * v2
            + cross * uv
            - 2.0 * beta * ua
            - 2.0 * alpha * va
            + a2
            <= 0.0
        )
        for variable_name, variable in (
            ("x", x),
            ("y", y),
            ("r", residual),
            ("u", u),
            ("v", v),
            ("u2", u2),
            ("v2", v2),
            ("uv", uv),
            ("ua", ua),
            ("va", va),
            ("a2", a2),
        ):
            variables[f"{name}_{variable_name}_{z}"] = variable


def build(
    terminal: np.ndarray,
    support_weight: float,
    prefix_order: tuple[int, int, int, int],
    pairs: tuple[tuple[Column, Column], ...],
    target: float | None,
    known_upper: float | None,
    hellinger_tangents: int,
    pair_geometry: str,
    inellipse_reciprocal_cap: float,
    coordinate_cases: tuple[str, ...],
) -> tuple[Model, dict[str, object]]:
    weights = np.trace(terminal[:3], axis1=1, axis2=2).real
    directions = np.asarray(
        [bloch(terminal[t] / weights[t])[1:] for t in range(3)]
    )
    model = Model("pairwise-qubit-helstrom")
    variables: dict[str, object] = {}

    statistics: dict[tuple[int, int, int], object] = {}
    probability: dict[tuple[int, int], object] = {}
    for z, y in PATHS:
        p = model.addVar(lb=0.0, ub=1.0, name=f"p_{z}_{y}")
        items = []
        for t in range(3):
            item = model.addVar(lb=0.0, ub=float(weights[t]), name=f"q_{z}_{y}_{t}")
            model.addCons(item <= float(weights[t]) * p)
            statistics[z, y, t] = item
            variables[f"q_{z}_{y}_{t}"] = item
            items.append(item)
        model.addCons(p == quicksum(items))
        probability[z, y] = p
        variables[f"p_{z}_{y}"] = p
    model.addCons(quicksum(probability.values()) == 1.0)

    prefix = {
        z: quicksum(probability[z, y] for y in OUTCOMES) for z in OUTCOMES
    }
    for index in range(3):
        model.addCons(
            prefix[prefix_order[index]] >= prefix[prefix_order[index + 1]]
        )
    model.addCons(prefix[prefix_order[0]] >= 0.25)
    for rank in range(1, 4):
        model.addCons(prefix[prefix_order[rank]] <= 1.0 / (rank + 1.0))
    for pair_index, (pair, coordinate_case) in enumerate(
        zip(pairs, coordinate_cases, strict=True)
    ):
        add_pair = (
            add_pair_factorisation
            if pair_geometry == "effects"
            else add_pair_inellipse
        )
        arguments = (
            model,
            variables,
            pair_index,
            pair,
            prefix,
            statistics,
            probability,
            weights,
        )
        if pair_geometry == "effects":
            add_pair(*arguments)
        else:
            add_pair(
                *arguments,
                inellipse_reciprocal_cap,
                coordinate_case,
            )

    terminal_statistics = {
        (syndrome, t): quicksum(
            statistics[z, y, t]
            for z, y in PATHS
            if (z ^ y) == syndrome
        )
        for syndrome in OUTCOMES
        for t in range(3)
    }
    correct = [
        terminal_statistics[syndrome, syndrome]
        for syndrome in range(3)
    ]
    audit = model.addVar(lb=0.0, ub=1.0, name="audit")
    model.addCons(audit == quicksum(correct))
    variables["audit"] = audit

    maximum = float(weights.max())
    cap = [maximum, maximum, 2.0 - 2.0 * maximum, 0.0]
    model.addCons(
        audit
        <= quicksum(
            cap[index] * prefix[prefix_order[index]] for index in OUTCOMES
        )
    )

    reconstruction = np.asarray(
        [
            [
                0.5 * weights[t],
                0.5 * weights[t] * directions[t, 0],
                0.5 * weights[t] * directions[t, 1],
            ]
            for t in range(3)
        ]
    )
    inverse = np.linalg.inv(reconstruction)
    terminal_prior: dict[int, object] = {}
    terminal_vector: dict[int, tuple[object, object, object]] = {}
    for syndrome in OUTCOMES:
        reconstructed = tuple(
            quicksum(
                float(inverse[row, t]) * terminal_statistics[syndrome, t]
                for t in range(3)
            )
            for row in range(3)
        )
        prior_from_paths = quicksum(
            probability[z, z ^ syndrome] for z in OUTCOMES
        )
        model.addCons(reconstructed[0] == prior_from_paths)
        normal = model.addVar(lb=-1.0, ub=1.0, name=f"tau_{syndrome}_z")
        vector = (reconstructed[1], reconstructed[2], normal)
        add_lorentz(model, reconstructed[0], vector)
        terminal_prior[syndrome] = reconstructed[0]
        terminal_vector[syndrome] = vector
        variables[f"tau_{syndrome}_z"] = normal

    dual_vector = tuple(
        model.addVar(lb=-1.0, ub=1.0, name=f"dual_{axis}")
        for axis in range(3)
    )
    add_lorentz(model, audit, dual_vector)
    for syndrome in OUTCOMES:
        add_lorentz(
            model,
            audit - terminal_prior[syndrome],
            tuple(
                dual_vector[axis] - terminal_vector[syndrome][axis]
                for axis in range(3)
            ),
        )
    # Direct terminal-weight cap; redundant physically, useful at the root.
    model.addCons(
        audit
        <= quicksum(
            float(weights[t]) * terminal_prior[t] for t in OUTCOMES if t < 3
        )
    )

    flat = [probability[z, y] for z, y in PATHS]
    hellinger = []
    for first in range(16):
        for second in range(first + 1, 16):
            item = model.addVar(lb=0.0, ub=0.5, name=f"h_{first}_{second}")
            model.addCons(item * item <= flat[first] * flat[second])
            # Tangent majorants of sqrt(p*q).  They are redundant for the
            # exact rotated-cone constraint but make SCIP's root LP much
            # closer to the conic Hellinger hypograph.
            if hellinger_tangents > 0:
                for exponent in np.linspace(
                    -8.0, 8.0, hellinger_tangents
                ):
                    alpha = float(2.0**exponent)
                    model.addCons(
                        2.0 * item
                        <= alpha * flat[first] + flat[second] / alpha
                    )
            hellinger.append(item)
            variables[f"h_{first}_{second}"] = item
    returned = (1.0 + 2.0 * quicksum(hellinger)) / 16.0
    score = model.addVar(lb=0.0, ub=1.0, name="score")
    model.addCons(
        score <= support_weight * audit + (1.0 - support_weight) * returned
    )
    if known_upper is not None:
        model.addCons(score <= known_upper)
    if target is not None:
        model.addCons(score >= target)
    model.setObjective(score, "maximize")
    variables.update(
        {
            "score": score,
            "return": returned,
            "terminal_effect_weights": weights,
        }
    )
    return model, variables


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lambda", dest="weight", type=float, default=0.6)
    parser.add_argument("--fixed-three-povm-weights", type=float, nargs=3, required=True)
    parser.add_argument("--prefix-order", type=int, nargs=4, required=True)
    parser.add_argument(
        "--pair",
        action="append",
        default=[],
        help="compatible effect pair, for example b_0,d_1_0",
    )
    parser.add_argument("--hellinger-tangents", type=int, default=0)
    parser.add_argument(
        "--pair-geometry",
        choices=("effects", "inellipse"),
        default="effects",
    )
    parser.add_argument("--inellipse-reciprocal-cap", type=float, default=2.0)
    parser.add_argument(
        "--coordinate-case",
        action="append",
        choices=("xy", "xr", "yr"),
        default=[],
        help="one inellipse coordinate chart per selected pair",
    )
    parser.add_argument("--target", type=float)
    parser.add_argument(
        "--known-upper",
        type=float,
        help="independently certified redundant upper bound for this weight cell",
    )
    parser.add_argument("--seconds", type=float, default=300.0)
    parser.add_argument("--gap", type=float, default=1e-5)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    order = tuple(args.prefix_order)
    if sorted(order) != list(OUTCOMES):
        raise ValueError("prefix order must be a permutation of 0,1,2,3")
    pair_text = args.pair or [
        "b_0,d_1_0",
        "b_1,b_9",
        "b_3,b_6",
        "b_4,d_0_2",
    ]
    pairs = tuple(parse_pair(text) for text in pair_text)
    coordinate_cases = tuple(args.coordinate_case or ["xy"] * len(pairs))
    if len(coordinate_cases) != len(pairs):
        raise ValueError("supply exactly one --coordinate-case per pair")
    if args.pair_geometry == "inellipse" and args.inellipse_reciprocal_cap < 2.0:
        raise ValueError("the clean-POVM cover requires reciprocal cap at least 2")
    terminal = canonical_three_effect_povm(
        np.asarray(args.fixed_three_povm_weights, dtype=float)
    )
    model, variables = build(
        terminal,
        args.weight,
        order,
        pairs,
        args.target,
        args.known_upper,
        args.hellinger_tangents,
        args.pair_geometry,
        args.inellipse_reciprocal_cap,
        coordinate_cases,
    )
    model.setRealParam("limits/time", args.seconds)
    model.setRealParam("limits/gap", args.gap)
    model.setRealParam("numerics/feastol", 1e-9)
    model.setRealParam("numerics/dualfeastol", 1e-9)
    model.setIntParam("display/verblevel", 2)
    model.optimize()

    payload: dict[str, object] = {
        "weight": args.weight,
        "terminal_effect_weights": args.fixed_three_povm_weights,
        "prefix_order": list(order),
        "pairs": [
            [render_column(first), render_column(second)]
            for first, second in pairs
        ],
        "target": args.target,
        "known_upper": args.known_upper,
        "hellinger_tangents": args.hellinger_tangents,
        "pair_geometry": args.pair_geometry,
        "inellipse_reciprocal_cap": args.inellipse_reciprocal_cap,
        "coordinate_cases": list(coordinate_cases),
        "status": str(model.getStatus()),
        "primal_bound": float(model.getPrimalbound()),
        "dual_bound": float(model.getDualbound()),
        "absolute_gap": float(model.getDualbound() - model.getPrimalbound()),
        "relative_gap": float(model.getGap()),
        "nodes": int(model.getNNodes()),
        "solving_time": float(model.getSolvingTime()),
        "scope": (
            "global pairwise-qubit probability relaxation; distinct pairs "
            "are not required to share one qubit realization"
            if args.pair_geometry == "effects"
            else "one chart cell of the clean-POVM inellipse outer model; "
            "all chart assignments are required for the finite union"
        ),
    }
    solution = model.getBestSol()
    if solution is not None:
        payload.update(
            {
                "incumbent_score": float(
                    model.getSolVal(solution, variables["score"])
                ),
                "incumbent_audit": float(
                    model.getSolVal(solution, variables["audit"])
                ),
                "incumbent_return": float(
                    model.getSolVal(solution, variables["return"])
                ),
            }
        )
    rendered = json.dumps(payload, indent=2) + "\n"
    print(rendered)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

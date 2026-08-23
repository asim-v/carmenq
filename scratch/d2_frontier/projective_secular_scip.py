"""Global secular certificates for the projective two-block topologies.

For a proposed support value ``level``, the Perron rank-one update formula
reduces the support upper bound to

    (1-lambda)/8 * sum_i c_i**2 / (level-lambda*d_i) <= 1.

This model globally maximises the left side with SCIP.  It removes the four
prior amplitudes and the high-degree squared Hellinger objective used by
``projective_topology_scip.py``.  A dual bound at most one certifies the
selected projective topology, subject to SCIP's stated numerical tolerances.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from pyscipopt import Model, quicksum

from projective_topology_scip import KINDS, topology_seed


def build(
    weight: float,
    level: float,
    first_kind: str,
    second_kind: str,
) -> tuple[Model, dict[str, object]]:
    if first_kind not in KINDS or second_kind not in KINDS:
        raise ValueError("unknown topology")
    if not 0.0 < weight < 1.0:
        raise ValueError("weight must lie strictly between zero and one")
    if level <= weight:
        raise ValueError("level must exceed weight so every denominator is positive")

    model = Model(f"secular-{first_kind}-{second_kind}")
    variables: dict[str, object] = {}

    def unit_pair(name: str) -> tuple[object, object]:
        first = model.addVar(lb=0.0, ub=1.0, name=f"{name}_0")
        second = model.addVar(lb=0.0, ub=1.0, name=f"{name}_1")
        model.addCons(first * first + second * second == 1.0)
        variables[f"{name}_0"] = first
        variables[f"{name}_1"] = second
        return first, second

    root_x, root_one_x = unit_pair("x")
    root_y, root_one_y = unit_pair("y")
    model.addCons(root_x >= root_y)
    if first_kind == second_kind:
        # Exchange of the two coarse cosets: (x,y) -> (1-y,1-x).
        model.addCons(root_x * root_x + root_y * root_y <= 1.0)

    ratios = []
    ratio_upper = 2.0 / (level - weight)

    def add_group(
        name: str,
        kind: str,
        high: object,
        low: object,
        high_complement: object,
        low_complement: object,
    ) -> None:
        angle = None
        if kind == "rank":
            sine, cosine = unit_pair(f"{name}_angle")
            model.addCons(sine <= cosine)
            angle = (sine, cosine)

        for label in range(2):
            root_state, root_state_complement = unit_pair(
                f"{name}_state_{label}"
            )
            root_q, root_q_complement = unit_pair(f"{name}_q_{label}")
            model.addCons(
                root_q * root_q
                == high * high * root_state * root_state
                + low * low * root_state_complement * root_state_complement
            )
            model.addCons(
                root_q_complement * root_q_complement
                == high_complement
                * high_complement
                * root_state
                * root_state
                + low_complement
                * low_complement
                * root_state_complement
                * root_state_complement
            )
            if kind == "endpoint":
                decision = root_q * root_q if label == 0 else 0.0
            else:
                assert angle is not None
                sine, cosine = angle
                amplitude = model.addVar(
                    lb=0.0, ub=1.0, name=f"{name}_amplitude_{label}"
                )
                variables[f"{name}_amplitude_{label}"] = amplitude
                if label == 0:
                    expression = (
                        high * cosine * root_state
                        + low * sine * root_state_complement
                    )
                else:
                    expression = (
                        high * sine * root_state
                        + low * cosine * root_state_complement
                    )
                model.addCons(amplitude == expression)
                decision = amplitude * amplitude

            ratio = model.addVar(
                lb=0.0, ub=ratio_upper, name=f"{name}_ratio_{label}"
            )
            variables[f"{name}_ratio_{label}"] = ratio
            model.addCons(
                ratio * (level - weight * decision)
                <= (root_q + root_q_complement)
                * (root_q + root_q_complement)
            )
            ratios.append(ratio)

    add_group(
        "first", first_kind, root_x, root_y, root_one_x, root_one_y
    )
    add_group(
        "second",
        second_kind,
        root_one_y,
        root_one_x,
        root_y,
        root_x,
    )

    secular_sum = quicksum(ratios)
    scaled = (1.0 - weight) * secular_sum / 8.0
    model.setObjective(secular_sum, "maximize")
    variables.update({"secular_sum": secular_sum, "scaled": scaled})
    return model, variables


def seed_solution(
    model: Model,
    variables: dict[str, object],
    row: dict[str, object],
    weight: float,
    level: float,
    first_kind: str,
    second_kind: str,
) -> None:
    values: dict[str, float] = {}
    x_value = float(row["high_eigenvalue"])
    y_value = float(row["low_eigenvalue"])
    if first_kind == second_kind and x_value + y_value > 1.0:
        x_value, y_value = 1.0 - y_value, 1.0 - x_value
        first_coordinates = list(row["second_coordinates"])
        second_coordinates = list(row["first_coordinates"])
    else:
        first_coordinates = list(row["first_coordinates"])
        second_coordinates = list(row["second_coordinates"])
    values.update(
        {
            "x_0": math.sqrt(x_value),
            "x_1": math.sqrt(1.0 - x_value),
            "y_0": math.sqrt(y_value),
            "y_1": math.sqrt(1.0 - y_value),
        }
    )

    def group_values(
        name: str,
        kind: str,
        coordinates: list[float],
        high: float,
        low: float,
    ) -> None:
        cursor = 0
        if kind == "rank":
            angle = (math.pi / 2.0) * float(coordinates[0])
            if angle > math.pi / 4.0:
                angle = math.pi / 2.0 - angle
                coordinates = [coordinates[0], coordinates[2], coordinates[1]]
            values[f"{name}_angle_0"] = math.sin(angle)
            values[f"{name}_angle_1"] = math.cos(angle)
            cursor = 1
        for label in range(2):
            coordinate = float(coordinates[cursor + label])
            root_state = math.sqrt(coordinate)
            root_state_complement = math.sqrt(1.0 - coordinate)
            q_value = low + (high - low) * coordinate
            values[f"{name}_state_{label}_0"] = root_state
            values[f"{name}_state_{label}_1"] = root_state_complement
            values[f"{name}_q_{label}_0"] = math.sqrt(q_value)
            values[f"{name}_q_{label}_1"] = math.sqrt(1.0 - q_value)
            if kind == "endpoint":
                decision = q_value if label == 0 else 0.0
            else:
                sine = values[f"{name}_angle_0"]
                cosine = values[f"{name}_angle_1"]
                if label == 0:
                    amplitude = (
                        math.sqrt(high) * cosine * root_state
                        + math.sqrt(low) * sine * root_state_complement
                    )
                else:
                    amplitude = (
                        math.sqrt(high) * sine * root_state
                        + math.sqrt(low) * cosine * root_state_complement
                    )
                values[f"{name}_amplitude_{label}"] = amplitude
                decision = amplitude * amplitude
            c_squared = (
                math.sqrt(q_value) + math.sqrt(1.0 - q_value)
            ) ** 2
            values[f"{name}_ratio_{label}"] = c_squared / (
                level - weight * decision
            )

    group_values(
        "first", first_kind, first_coordinates, x_value, y_value
    )
    group_values(
        "second",
        second_kind,
        second_coordinates,
        1.0 - y_value,
        1.0 - x_value,
    )

    solution = model.createSol()
    for name, value in values.items():
        model.setSolVal(solution, variables[name], value)
    model.addSol(solution)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lambda", dest="weight", type=float, default=0.6)
    parser.add_argument("--level", type=float, required=True)
    parser.add_argument("--first", choices=KINDS, required=True)
    parser.add_argument("--second", choices=KINDS, required=True)
    parser.add_argument("--seed-json", type=Path)
    parser.add_argument("--seconds", type=float, default=300.0)
    parser.add_argument("--gap", type=float, default=1e-6)
    parser.add_argument(
        "--square-bound",
        nargs=3,
        action="append",
        default=[],
        metavar=("PAIR", "LOW", "HIGH"),
        help=(
            "bound the squared first component of a named unit pair; "
            "for example: --square-bound x 0.94 0.97"
        ),
    )
    parser.add_argument(
        "--linear-bound",
        nargs=3,
        action="append",
        default=[],
        metavar=("VARIABLE", "LOW", "HIGH"),
        help="bound a model variable directly, for example first_angle_0",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    model, variables = build(
        args.weight, args.level, args.first, args.second
    )
    rendered_bounds: list[dict[str, object]] = []
    for pair, lower_text, upper_text in args.square_bound:
        lower = float(lower_text)
        upper = float(upper_text)
        if not 0.0 <= lower <= upper <= 1.0:
            raise ValueError(f"invalid squared bound for {pair!r}")
        variable_name = f"{pair}_0"
        if variable_name not in variables:
            raise ValueError(f"unknown unit pair {pair!r}")
        model.chgVarLb(variables[variable_name], math.sqrt(lower))
        model.chgVarUb(variables[variable_name], math.sqrt(upper))
        rendered_bounds.append(
            {"kind": "squared", "name": pair, "lower": lower, "upper": upper}
        )
    for name, lower_text, upper_text in args.linear_bound:
        lower = float(lower_text)
        upper = float(upper_text)
        if lower > upper or name not in variables:
            raise ValueError(f"invalid direct bound for {name!r}")
        model.chgVarLb(variables[name], lower)
        model.chgVarUb(variables[name], upper)
        rendered_bounds.append(
            {"kind": "linear", "name": name, "lower": lower, "upper": upper}
        )
    seed_row = None
    if args.seed_json is not None:
        seed_row = topology_seed(args.seed_json, args.first, args.second)
        seed_solution(
            model,
            variables,
            seed_row,
            args.weight,
            args.level,
            args.first,
            args.second,
        )
    model.setRealParam("limits/time", args.seconds)
    model.setRealParam("limits/gap", args.gap)
    model.setRealParam("numerics/feastol", 1e-9)
    model.setRealParam("numerics/dualfeastol", 1e-9)
    model.setIntParam("display/verblevel", 2)
    model.optimize()

    solution = model.getBestSol()
    scale = (1.0 - args.weight) / 8.0
    payload: dict[str, object] = {
        "weight": args.weight,
        "level": args.level,
        "first_topology": args.first,
        "second_topology": args.second,
        "bounds": rendered_bounds,
        "status": str(model.getStatus()),
        "seed_support": float(seed_row["score"]) if seed_row else None,
        "primal_secular_sum": float(model.getPrimalbound()),
        "dual_secular_sum": float(model.getDualbound()),
        "primal_scaled": scale * float(model.getPrimalbound()),
        "dual_scaled": scale * float(model.getDualbound()),
        "absolute_scaled_gap": scale
        * float(model.getDualbound() - model.getPrimalbound()),
        "relative_gap": float(model.getGap()),
        "nodes": int(model.getNNodes()),
        "solving_time": float(model.getSolvingTime()),
    }
    if solution is not None:
        payload["incumbent_scaled"] = float(
            model.getSolVal(solution, variables["scaled"])
        )
        unit_pairs = ["x", "y"]
        for group, kind in (
            ("first", args.first), ("second", args.second)
        ):
            if kind == "rank":
                unit_pairs.append(f"{group}_angle")
            unit_pairs.extend(
                [
                    f"{group}_state_0",
                    f"{group}_state_1",
                    f"{group}_q_0",
                    f"{group}_q_1",
                ]
            )
        payload["incumbent_unit_pair_first_squares"] = {
            name: float(model.getSolVal(solution, variables[f"{name}_0"]))
            ** 2
            for name in unit_pairs
        }
    rendered = json.dumps(payload, indent=2) + "\n"
    print(rendered)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

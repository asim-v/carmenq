"""Global polynomial bounds for the four projective binary-line topologies.

The model retains independent eigenvalues of the coarse qubit effect,
independent rank-split angles in its two cosets, four pure prefix states, and
all four Perron prior amplitudes.  SCIP's reported dual bound certifies only
the selected projective topology, up to its numerical tolerances.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from pyscipopt import Model, quicksum


KINDS = ("endpoint", "rank")


def build(
    weight: float, first_kind: str, second_kind: str
) -> tuple[Model, dict[str, object]]:
    if first_kind not in KINDS or second_kind not in KINDS:
        raise ValueError("unknown topology")
    model = Model(f"projective-{first_kind}-{second_kind}")
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
        # Exchanging the two cosets sends (x,y) to (1-y,1-x).
        model.addCons(root_x * root_x + root_y * root_y <= 1.0)

    all_d = []
    all_c = []

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
            all_c.append(root_q + root_q_complement)
            if kind == "endpoint":
                all_d.append(root_q * root_q if label == 0 else 0.0)
            else:
                assert angle is not None
                sine, cosine = angle
                if label == 0:
                    amplitude = (
                        high * cosine * root_state
                        + low * sine * root_state_complement
                    )
                else:
                    amplitude = (
                        high * sine * root_state
                        + low * cosine * root_state_complement
                    )
                all_d.append(amplitude * amplitude)

    add_group(
        "first",
        first_kind,
        root_x,
        root_y,
        root_one_x,
        root_one_y,
    )
    # The second coarse effect is I-E.  Reordering its eigenbasis puts
    # 1-y first and 1-x second; its complementary roots are y and x.
    add_group(
        "second",
        second_kind,
        root_one_y,
        root_one_x,
        root_y,
        root_x,
    )

    priors = [model.addVar(lb=0.0, ub=1.0, name=f"prior_{i}") for i in range(4)]
    variables["priors"] = priors
    model.addCons(quicksum(value * value for value in priors) == 1.0)
    audit = quicksum(priors[i] * priors[i] * all_d[i] for i in range(4))
    hellinger = quicksum(priors[i] * all_c[i] for i in range(4))
    returned = hellinger * hellinger / 8.0
    expression = weight * audit + (1.0 - weight) * returned
    score = model.addVar(lb=0.0, ub=1.0, name="score")
    model.addCons(score <= expression)
    model.setObjective(score, "maximize")
    variables.update(
        {"audit": audit, "return": returned, "score": score}
    )
    return model, variables


def topology_seed(
    path: Path, first_kind: str, second_kind: str
) -> dict[str, object]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    for row in payload["topologies"]:
        if (
            row["first_topology"] == first_kind
            and row["second_topology"] == second_kind
        ):
            return row
    raise ValueError("requested topology is absent from seed JSON")


def seed_solution(
    model: Model,
    variables: dict[str, object],
    row: dict[str, object],
    first_kind: str,
    second_kind: str,
) -> None:
    values: dict[str, float] = {}
    x_value = float(row["high_eigenvalue"])
    y_value = float(row["low_eigenvalue"])
    if first_kind == second_kind and x_value + y_value > 1.0:
        # Use the coset-exchanged representative required by the model.
        x_value, y_value = 1.0 - y_value, 1.0 - x_value
        first_coordinates = list(row["second_coordinates"])
        second_coordinates = list(row["first_coordinates"])
        prior_values = list(row["priors"])[2:] + list(row["priors"])[:2]
    else:
        first_coordinates = list(row["first_coordinates"])
        second_coordinates = list(row["second_coordinates"])
        prior_values = list(row["priors"])
    values.update(
        {
            "x_0": math.sqrt(x_value),
            "x_1": math.sqrt(1.0 - x_value),
            "y_0": math.sqrt(y_value),
            "y_1": math.sqrt(1.0 - y_value),
            "score": float(row["score"]),
        }
    )

    def group_values(
        name: str, kind: str, coordinates: list[float]
    ) -> bool:
        """Populate one group and report whether its two fine labels swapped."""
        cursor = 0
        swapped = False
        if kind == "rank":
            angle = (math.pi / 2.0) * float(coordinates[0])
            if angle > math.pi / 4.0:
                angle = math.pi / 2.0 - angle
                coordinates = [coordinates[0], coordinates[2], coordinates[1]]
                swapped = True
            values[f"{name}_angle_0"] = math.sin(angle)
            values[f"{name}_angle_1"] = math.cos(angle)
            cursor = 1
        for label in range(2):
            coordinate = float(coordinates[cursor + label])
            values[f"{name}_state_{label}_0"] = math.sqrt(coordinate)
            values[f"{name}_state_{label}_1"] = math.sqrt(1.0 - coordinate)
        return swapped

    first_swapped = group_values("first", first_kind, first_coordinates)
    second_swapped = group_values("second", second_kind, second_coordinates)
    if first_swapped:
        prior_values[0], prior_values[1] = prior_values[1], prior_values[0]
    if second_swapped:
        prior_values[2], prior_values[3] = prior_values[3], prior_values[2]

    # Recompute the coarse probabilities from the possibly exchanged seed.
    for name, high, low in (
        ("first", x_value, y_value),
        ("second", 1.0 - y_value, 1.0 - x_value),
    ):
        for label in range(2):
            state_root = values[f"{name}_state_{label}_0"]
            q_value = low + (high - low) * state_root * state_root
            values[f"{name}_q_{label}_0"] = math.sqrt(q_value)
            values[f"{name}_q_{label}_1"] = math.sqrt(1.0 - q_value)

    solution = model.createSol()
    for name, value in values.items():
        model.setSolVal(solution, variables[name], value)
    for index, prior in enumerate(prior_values):
        model.setSolVal(solution, variables["priors"][index], math.sqrt(float(prior)))
    model.addSol(solution)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lambda", dest="weight", type=float, default=0.6)
    parser.add_argument("--first", choices=KINDS, required=True)
    parser.add_argument("--second", choices=KINDS, required=True)
    parser.add_argument("--seed-json", type=Path)
    parser.add_argument("--seconds", type=float, default=300.0)
    parser.add_argument("--gap", type=float, default=1e-6)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    model, variables = build(args.weight, args.first, args.second)
    seed_row = None
    if args.seed_json is not None:
        seed_row = topology_seed(args.seed_json, args.first, args.second)
        seed_solution(
            model, variables, seed_row, args.first, args.second
        )
    model.setRealParam("limits/time", args.seconds)
    model.setRealParam("limits/gap", args.gap)
    model.setRealParam("numerics/feastol", 1e-9)
    model.setRealParam("numerics/dualfeastol", 1e-9)
    model.setIntParam("display/verblevel", 2)
    model.optimize()

    solution = model.getBestSol()
    payload: dict[str, object] = {
        "weight": args.weight,
        "first_topology": args.first,
        "second_topology": args.second,
        "status": str(model.getStatus()),
        "seed_lower_bound": float(seed_row["score"]) if seed_row else None,
        "primal_bound": float(model.getPrimalbound()),
        "dual_bound": float(model.getDualbound()),
        "absolute_gap": float(model.getDualbound() - model.getPrimalbound()),
        "relative_gap": float(model.getGap()),
        "nodes": int(model.getNNodes()),
        "solving_time": float(model.getSolvingTime()),
    }
    if solution is not None:
        payload["audit"] = float(model.getSolVal(solution, variables["audit"]))
        payload["return"] = float(model.getSolVal(solution, variables["return"]))
    rendered = json.dumps(payload, indent=2) + "\n"
    print(rendered)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

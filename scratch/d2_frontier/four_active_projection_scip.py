"""Spatial global certificate for the fully active four-readout relaxation.

This is the polynomial projection relaxation implied by Helstrom
complementarity.  The sorted terminal effects have weights ``w_i``, unit
Bloch supports with projections ``x_i`` on the dual axis, and dual spectral
bias ``t``.  Only the necessary projected closure

    sum_i w_i x_i = 0

is retained; omitting transverse closure enlarges the physical set.

For a fully active readout, positivity of the four syndrome states is exactly

    2 p_i (1+t x_i) >= A (1+2t x_i+t^2),

and dual dominance also gives ``p_i <= A``.  Averaging the exact aligned-
projective comparisons over every retained label except ``j`` yields

    A2_j = A - k_j(A-p_j),   k_j=(1-w_j)/(2-w_j).

The independently certified projective support lines constrain every
``A2_j``.  Four auxiliary square-root variables encode the Hellinger RETURN.
Every physical fully active four-outcome leaf is feasible in this model, so
a spatial dual bound is an upper certificate for that sector.

The certificate remains numerical and conditional on SCIP's recorded
tolerances; it is not an interval-arithmetic proof independent of SCIP.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

from pyscipopt import Model, quicksum


def build(
    support_weight: float,
    maximum_weight_floor: float,
    minimum_active_weight: float,
    projective_lines: tuple[tuple[float, float], ...],
) -> tuple[Model, dict[str, Any]]:
    model = Model("four-active-projection-relaxation")
    weights = [
        model.addVar(
            lb=minimum_active_weight,
            ub=1.0,
            name=f"weight_{index}",
        )
        for index in range(4)
    ]
    projections = [
        model.addVar(lb=-1.0, ub=1.0, name=f"projection_{index}")
        for index in range(4)
    ]
    bias = model.addVar(lb=0.0, ub=1.0, name="dual_bias")
    audit = model.addVar(lb=0.0, ub=1.0, name="audit")
    returned = model.addVar(lb=0.0, ub=1.0, name="return")
    priors = [
        model.addVar(lb=0.0, ub=1.0, name=f"prior_{index}")
        for index in range(4)
    ]
    roots = [
        model.addVar(lb=0.0, ub=1.0, name=f"root_prior_{index}")
        for index in range(4)
    ]
    loss_factors = [
        model.addVar(lb=0.0, ub=0.5, name=f"loss_factor_{index}")
        for index in range(4)
    ]
    transverse_roots = [
        model.addVar(lb=0.0, ub=1.0, name=f"transverse_root_{index}")
        for index in range(4)
    ]
    transverse_weights = [
        model.addVar(lb=0.0, ub=1.0, name=f"transverse_weight_{index}")
        for index in range(4)
    ]

    model.addCons(quicksum(weights) == 2.0)
    model.addCons(weights[0] >= maximum_weight_floor)
    for index in range(3):
        model.addCons(weights[index] >= weights[index + 1])
    model.addCons(
        quicksum(weights[index] * projections[index] for index in range(4))
        == 0.0
    )
    model.addCons(quicksum(priors) == 1.0)
    for index in range(4):
        model.addCons(priors[index] <= audit)
        model.addCons(roots[index] * roots[index] <= priors[index])
        model.addCons(
            2.0
            * priors[index]
            * (1.0 + bias * projections[index])
            >= audit
            * (
                1.0
                + 2.0 * bias * projections[index]
                + bias * bias
            )
        )
        model.addCons(
            loss_factors[index] * (2.0 - weights[index])
            == 1.0 - weights[index]
        )
        model.addCons(
            projections[index] * projections[index]
            + transverse_roots[index] * transverse_roots[index]
            == 1.0
        )
        model.addCons(
            transverse_weights[index]
            == weights[index] * transverse_roots[index]
        )
    transverse_sum = quicksum(transverse_weights)
    for index in range(4):
        # Four planar vectors of these lengths can close iff no length
        # exceeds the sum of the other three.
        model.addCons(2.0 * transverse_weights[index] <= transverse_sum)
    model.addCons(4.0 * returned <= quicksum(roots) * quicksum(roots))

    comparisons = []
    for index in range(4):
        projective_audit = audit - loss_factors[index] * (
            audit - priors[index]
        )
        for line_weight, line_upper in projective_lines:
            constraint = model.addCons(
                line_weight * projective_audit
                + (1.0 - line_weight) * returned
                <= line_upper
            )
            comparisons.append((line_weight, line_upper, index, constraint))

    score = support_weight * audit + (1.0 - support_weight) * returned
    model.setObjective(score, "maximize")
    return model, {
        "weights": weights,
        "projections": projections,
        "bias": bias,
        "audit": audit,
        "return": returned,
        "priors": priors,
        "roots": roots,
        "loss_factors": loss_factors,
        "transverse_roots": transverse_roots,
        "transverse_weights": transverse_weights,
        "score": score,
        "comparisons": comparisons,
    }


def seed_solution(
    model: Model,
    variables: dict[str, Any],
    row: dict[str, Any],
) -> None:
    solution = model.createSol()
    weights = list(map(float, row["weights"]))
    projections = list(map(float, row["bloch_projections"]))
    priors = list(map(float, row["syndrome_priors"]))
    values = []
    for index in range(4):
        values.extend(
            (
                (variables["weights"][index], weights[index]),
                (variables["projections"][index], projections[index]),
                (variables["priors"][index], priors[index]),
                (
                    variables["roots"][index],
                    math.sqrt(max(0.0, priors[index])),
                ),
                (
                    variables["loss_factors"][index],
                    (1.0 - weights[index]) / (2.0 - weights[index]),
                ),
                (
                    variables["transverse_roots"][index],
                    math.sqrt(max(0.0, 1.0 - projections[index] ** 2)),
                ),
                (
                    variables["transverse_weights"][index],
                    weights[index]
                    * math.sqrt(max(0.0, 1.0 - projections[index] ** 2)),
                ),
            )
        )
    values.extend(
        (
            (
                variables["bias"],
                float(row["geometry"]["dual_spectral_bias"]),
            ),
            (variables["audit"], float(row["audit"])),
            (variables["return"], float(row["return"])),
        )
    )
    for variable, value in values:
        model.setSolVal(solution, variable, value)
    model.addSol(solution)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lambda", dest="support_weight", type=float, default=0.6)
    parser.add_argument("--maximum-weight-floor", type=float, default=0.88325)
    parser.add_argument("--minimum-active-weight", type=float, default=0.0003)
    parser.add_argument(
        "--projective-line",
        type=float,
        nargs=2,
        action="append",
        default=None,
        metavar=("WEIGHT", "UPPER"),
    )
    parser.add_argument("--seed-json", type=Path)
    parser.add_argument("--seconds", type=float, default=600.0)
    parser.add_argument("--gap", type=float, default=1e-5)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    raw_lines = (
        ((0.55, 0.7573), (0.6, 0.76591))
        if args.projective_line is None
        else args.projective_line
    )
    lines = tuple(tuple(map(float, line)) for line in raw_lines)
    model, variables = build(
        args.support_weight,
        args.maximum_weight_floor,
        args.minimum_active_weight,
        lines,
    )
    if args.seed_json is not None:
        seed_solution(
            model,
            variables,
            json.loads(args.seed_json.read_text(encoding="utf-8")),
        )
    model.setRealParam("limits/time", args.seconds)
    model.setRealParam("limits/gap", args.gap)
    model.setRealParam("numerics/feastol", 1e-9)
    model.setRealParam("numerics/dualfeastol", 1e-9)
    model.setIntParam("display/verblevel", 3)
    model.optimize()

    solution = model.getBestSol()
    payload: dict[str, Any] = {
        "support_weight": args.support_weight,
        "maximum_weight_floor": args.maximum_weight_floor,
        "minimum_active_weight": args.minimum_active_weight,
        "projective_lines": [list(line) for line in lines],
        "status": str(model.getStatus()),
        "primal_bound": float(model.getPrimalbound()),
        "dual_bound": float(model.getDualbound()),
        "absolute_gap": float(model.getDualbound() - model.getPrimalbound()),
        "relative_gap": float(model.getGap()),
        "nodes": int(model.getNNodes()),
        "solving_time": float(model.getSolvingTime()),
        "certificate_scope": (
            "fully active four-outcome terminal readout with every effect "
            "trace at least minimum_active_weight; projected Helstrom outer "
            "relaxation and averaged exact projective comparisons"
        ),
        "numerical_status": (
            "spatial SCIP dual at recorded tolerances; outward interval "
            "validation pending"
        ),
    }
    if solution is not None:
        payload["incumbent"] = {
            "score": float(model.getSolVal(solution, variables["score"])),
            "weights": [
                float(model.getSolVal(solution, item))
                for item in variables["weights"]
            ],
            "projections": [
                float(model.getSolVal(solution, item))
                for item in variables["projections"]
            ],
            "dual_bias": float(model.getSolVal(solution, variables["bias"])),
            "audit": float(model.getSolVal(solution, variables["audit"])),
            "return": float(model.getSolVal(solution, variables["return"])),
            "priors": [
                float(model.getSolVal(solution, item))
                for item in variables["priors"]
            ],
            "transverse_weights": [
                float(model.getSolVal(solution, item))
                for item in variables["transverse_weights"]
            ],
        }
    rendered = json.dumps(payload, indent=2) + "\n"
    print(rendered, end="")
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

"""Trace-only global relaxation for projective rank/rank certificates.

For a congruence-rank-one split of a coarse qubit effect with ordered
eigenvalues ``a >= b``, the two fine effects have traces ``t`` and
``a+b-t`` with ``(a+b)/2 <= t <= a``.  For any state, its coarse
probability ``q`` and correct probability ``d`` satisfy

    b <= q <= a,   0 <= d <= min(q, t).

Dropping the angle-dependent relation between ``q`` and ``d`` gives a small
outer relaxation of the exact secular program.  It is strong enough to prune
most eigenvalue boxes before running the full spatial branch-and-bound model.
Every reported dual bound is an upper bound on the exact rank/rank topology;
the converse is not true.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from pyscipopt import Model, quicksum

from projective_topology_scip import KINDS


def build(
    weight: float,
    level: float,
    first_kind: str = "rank",
    second_kind: str = "rank",
) -> tuple[Model, dict[str, object]]:
    if not 0.0 < weight < 1.0 or level <= weight:
        raise ValueError("require 0 < lambda < 1 and level > lambda")
    if first_kind not in KINDS or second_kind not in KINDS:
        raise ValueError("unknown topology")
    model = Model(f"projective-trace-{first_kind}-{second_kind}")
    variables: dict[str, object] = {}

    x = model.addVar(lb=0.0, ub=1.0, name="x")
    y = model.addVar(lb=0.0, ub=0.5, name="y")
    variables.update({"x": x, "y": y})
    model.addCons(x >= y)
    model.addCons(x + y <= 1.0)

    ratios = []
    ratio_upper = 2.0 / (level - weight)

    def add_effect(
        name: str,
        low: object,
        high: object,
        trace: object | None,
        endpoint: str | None = None,
    ) -> None:
        q = model.addVar(lb=0.0, ub=1.0, name=f"{name}_q")
        d = model.addVar(lb=0.0, ub=1.0, name=f"{name}_d")
        root_q = model.addVar(lb=0.0, ub=1.0, name=f"{name}_root_q")
        root_one_q = model.addVar(
            lb=0.0, ub=1.0, name=f"{name}_root_one_q"
        )
        ratio = model.addVar(
            lb=0.0, ub=ratio_upper, name=f"{name}_ratio"
        )
        variables.update(
            {
                f"{name}_q": q,
                f"{name}_d": d,
                f"{name}_root_q": root_q,
                f"{name}_root_one_q": root_one_q,
                f"{name}_ratio": ratio,
            }
        )
        model.addCons(q >= low)
        model.addCons(q <= high)
        if endpoint == "active":
            model.addCons(d == q)
        elif endpoint == "null":
            model.addCons(d == 0.0)
        else:
            if trace is None:
                raise ValueError("rank effect requires a trace upper bound")
            model.addCons(d <= q)
            model.addCons(d <= trace)
        model.addCons(root_q * root_q == q)
        model.addCons(root_one_q * root_one_q == 1.0 - q)
        model.addCons(
            ratio * (level - weight * d)
            <= (root_q + root_one_q) * (root_q + root_one_q)
        )
        ratios.append(ratio)

    def add_group(
        name: str, kind: str, low: object, high: object, total: object
    ) -> None:
        if kind == "endpoint":
            add_effect(f"{name}_0", low, high, None, endpoint="active")
            add_effect(f"{name}_1", low, high, None, endpoint="null")
            return
        trace = model.addVar(lb=0.0, ub=1.0, name=f"{name}_trace")
        variables[f"{name}_trace"] = trace
        model.addCons(2.0 * trace >= total)
        model.addCons(trace <= high)
        add_effect(f"{name}_0", low, high, trace)
        add_effect(f"{name}_1", low, high, total - trace)

    add_group("first", first_kind, y, x, x + y)
    add_group("second", second_kind, 1.0 - x, 1.0 - y, 2.0 - x - y)

    secular_sum = quicksum(ratios)
    scaled = (1.0 - weight) * secular_sum / 8.0
    model.setObjective(secular_sum, "maximize")
    variables.update({"secular_sum": secular_sum, "scaled": scaled})
    return model, variables


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lambda", dest="weight", type=float, default=0.6)
    parser.add_argument("--level", type=float, required=True)
    parser.add_argument("--first", choices=KINDS, default="rank")
    parser.add_argument("--second", choices=KINDS, default="rank")
    parser.add_argument("--x", type=float, nargs=2, metavar=("LOW", "HIGH"))
    parser.add_argument("--y", type=float, nargs=2, metavar=("LOW", "HIGH"))
    parser.add_argument(
        "--first-trace", type=float, nargs=2, metavar=("LOW", "HIGH")
    )
    parser.add_argument(
        "--second-trace", type=float, nargs=2, metavar=("LOW", "HIGH")
    )
    parser.add_argument("--seconds", type=float, default=60.0)
    parser.add_argument("--gap", type=float, default=1e-5)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    model, variables = build(
        args.weight, args.level, args.first, args.second
    )
    bounds: dict[str, list[float]] = {}
    for name, supplied in (
        ("x", args.x),
        ("y", args.y),
        ("first_trace", args.first_trace),
        ("second_trace", args.second_trace),
    ):
        if supplied is None:
            continue
        lower, upper = map(float, supplied)
        if not 0.0 <= lower <= upper <= 1.0:
            raise ValueError(f"invalid {name} bounds")
        model.chgVarLb(variables[name], lower)
        model.chgVarUb(variables[name], upper)
        bounds[name] = [lower, upper]

    model.setRealParam("limits/time", args.seconds)
    model.setRealParam("limits/gap", args.gap)
    model.setRealParam("numerics/feastol", 1e-9)
    model.setRealParam("numerics/dualfeastol", 1e-9)
    model.setIntParam("display/verblevel", 2)
    model.optimize()

    scale = (1.0 - args.weight) / 8.0
    payload = {
        "weight": args.weight,
        "level": args.level,
        "first_topology": args.first,
        "second_topology": args.second,
        "bounds": bounds,
        "status": str(model.getStatus()),
        "primal_scaled": scale * float(model.getPrimalbound()),
        "dual_scaled": scale * float(model.getDualbound()),
        "absolute_scaled_gap": scale
        * float(model.getDualbound() - model.getPrimalbound()),
        "relative_gap": float(model.getGap()),
        "nodes": int(model.getNNodes()),
        "solving_time": float(model.getSolvingTime()),
    }
    rendered = json.dumps(payload, indent=2) + "\n"
    print(rendered)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

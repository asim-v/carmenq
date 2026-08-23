"""Validate the finite projective support cover at lambda=0.55.

The archived SCIP jobs were solved at ``OLD_LEVEL``.  Their hierarchy covers
the four projective topologies and, for the rank/rank topology, the complete
reduced domain

    0 <= y <= x <= 1,    x + y <= 1.

The few difficult trace boxes are replaced by exact secular boxes, followed
by a finite partition in ``x``, ``y``, and the two rank-split angles.  This
script verifies that every replacement is a complete partition and retains
only terminal leaves.

The final support value is slightly larger than the level used by the SCIP
runs.  No new optimisation is needed.  For every feasible point and every
decision statistic ``d_i >= 0``, monotonicity gives

    1 / (L1 - lambda*d_i)
        <= (L0/L1) / (L0 - lambda*d_i),       L1 >= L0.

Consequently an archived dual bound ``D0`` at ``L0`` implies the rigorous
post-processing bound ``D1 <= (L0/L1) D0``.  The cover is a *numerical SCIP
certificate*: it is conditional on SCIP's spatial branch-and-bound and the
tolerances recorded in the JSON files, not an interval-arithmetic proof.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parent
OLD_LEVEL = 0.7568534
CERTIFIED_LEVEL = 0.75730
WEIGHT = 0.55
PHYSICAL_LOWER = 0.7568432934790399
TOL = 1e-12

COARSE_EXACT = {
    "x78_y01",
    "x89_y01",
    "x89_y12",
    "x910_y01",
    "x910_y12",
}
COARSE_REFINED = {"x89_y01", "x910_y01"}
FINE_ANGLE_PARENTS = {
    "x8890_y0205",
    "x8890_y0508",
    "x8890_y0810",
    "x9092_y0205",
    "x9092_y0508",
    "x9092_y0810",
    "x9295_y0205",
    "x9295_y0508",
}
PROBE_PARENT = "x9092_y0810"
LOW_ANGLE_INDICES = {(i, j) for i in range(2) for j in range(2)}
ANGLE_EDGES = (0.0, 0.1, 0.3, 0.5, math.sqrt(0.5))


def load_rows(directory: str) -> dict[str, dict[str, object]]:
    path = ROOT / directory
    rows = {
        item.stem: json.loads(item.read_text(encoding="utf-8"))
        for item in path.glob("*.json")
    }
    if not rows:
        raise RuntimeError(f"no certificate files found in {path}")
    return rows


def load_row(filename: str) -> dict[str, object]:
    return json.loads((ROOT / filename).read_text(encoding="utf-8"))


def close(first: float, second: float, *, tol: float = TOL) -> bool:
    return abs(first - second) <= tol


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def bounds_dict(row: dict[str, object]) -> dict[str, tuple[float, float]]:
    raw = row.get("bounds", [])
    if isinstance(raw, dict):
        return {
            str(name): (float(values[0]), float(values[1]))
            for name, values in raw.items()
        }
    result: dict[str, tuple[float, float]] = {}
    for item in raw:
        assert isinstance(item, dict)
        result[str(item["name"])] = (
            float(item["lower"]),
            float(item["upper"]),
        )
    return result


def check_metadata(
    row: dict[str, object],
    *,
    first: str = "rank",
    second: str = "rank",
) -> None:
    require(close(float(row["weight"]), WEIGHT), "wrong lambda in certificate")
    require(close(float(row["level"]), OLD_LEVEL), "wrong archived level")
    require(row["first_topology"] == first, "wrong first topology")
    require(row["second_topology"] == second, "wrong second topology")
    require(math.isfinite(float(row["dual_scaled"])), "non-finite dual bound")


def secular_rescaling_factor(old_level: float, new_level: float) -> float:
    """Return the uniform factor implied by nonnegative decision statistics."""

    if old_level <= WEIGHT or new_level < old_level:
        raise ValueError("require new_level >= old_level > lambda")
    return old_level / new_level


def rescaled_dual(row: dict[str, object]) -> float:
    return secular_rescaling_factor(OLD_LEVEL, CERTIFIED_LEVEL) * float(
        row["dual_scaled"]
    )


def add_leaf(
    leaves: list[dict[str, object]],
    *,
    topology: str,
    stage: str,
    name: str,
    row: dict[str, object],
) -> None:
    first, second = topology.split("/")
    check_metadata(row, first=first, second=second)
    dual = rescaled_dual(row)
    require(dual <= 1.0 + TOL, f"uncertified leaf {stage}/{name}: {dual}")
    leaves.append(
        {
            "topology": topology,
            "stage": stage,
            "name": name,
            "archived_dual_scaled": float(row["dual_scaled"]),
            "rescaled_dual_upper": dual,
            "status": str(row["status"]),
            "nodes": int(row["nodes"]),
            "solving_time": float(row["solving_time"]),
        }
    )


def expected_coarse_names() -> set[str]:
    names: set[str] = set()
    for i in range(10):
        x_lower, x_upper = i / 10.0, (i + 1) / 10.0
        for j in range(5):
            y_lower = j / 10.0
            # Retain every closed box intersecting y <= x and x + y <= 1.
            if x_upper + TOL >= y_lower and x_lower + y_lower <= 1.0 + TOL:
                names.add(f"x{i}{i + 1}_y{j}{j + 1}")
    return names


def expected_fine_bounds() -> set[tuple[float, float, float, float]]:
    boxes: set[tuple[float, float, float, float]] = set()
    for i in range(8):
        x_lower = 0.8 + 0.025 * i
        x_upper = x_lower + 0.025
        for j in range(4):
            y_lower = 0.025 * j
            y_upper = y_lower + 0.025
            if x_lower + y_lower <= 1.0 + TOL:
                boxes.add(tuple(round(v, 12) for v in (
                    x_lower, x_upper, y_lower, y_upper
                )))
    return boxes


def actual_xy_bounds(rows: Iterable[dict[str, object]]) -> set[tuple[float, ...]]:
    result = set()
    for row in rows:
        bounds = bounds_dict(row)
        result.add(tuple(round(v, 12) for v in (*bounds["x"], *bounds["y"])))
    return result


def expected_angle_bounds(i: int, j: int) -> dict[str, tuple[float, float]]:
    return {
        "first_angle_0": (ANGLE_EDGES[i], ANGLE_EDGES[i + 1]),
        "second_angle_0": (ANGLE_EDGES[j], ANGLE_EDGES[j + 1]),
    }


def check_selected_bounds(
    row: dict[str, object],
    expected: dict[str, tuple[float, float]],
    context: str,
) -> None:
    actual = bounds_dict(row)
    for name, pair in expected.items():
        require(name in actual, f"missing {name} bound in {context}")
        require(
            close(actual[name][0], pair[0]) and close(actual[name][1], pair[1]),
            f"wrong {name} interval in {context}: {actual[name]} != {pair}",
        )


def singleton_leaves() -> list[dict[str, object]]:
    leaves: list[dict[str, object]] = []
    for topology, filename in (
        ("endpoint/endpoint", "projective_l055_global_endpoint_endpoint.json"),
        ("endpoint/rank", "projective_l055_global_endpoint_rank.json"),
        ("rank/endpoint", "projective_l055_global_rank_endpoint.json"),
    ):
        add_leaf(
            leaves,
            topology=topology,
            stage="global-exact",
            name=Path(filename).stem,
            row=load_row(filename),
        )
    return leaves


def rank_rank_leaves() -> list[dict[str, object]]:
    leaves: list[dict[str, object]] = []
    coarse_trace = load_rows("projective_l055_trace_coarse")
    coarse_exact = load_rows("projective_l055_exact_coarse")
    require(
        set(coarse_trace) == expected_coarse_names(),
        "coarse trace boxes do not exactly cover the reduced (x,y) domain",
    )
    require(set(coarse_exact) == COARSE_EXACT, "unexpected coarse exact boxes")
    for name, row in coarse_trace.items():
        check_metadata(row)
        if name not in COARSE_EXACT:
            add_leaf(
                leaves, topology="rank/rank", stage="coarse-trace",
                name=name, row=row,
            )
    for name, row in coarse_exact.items():
        check_metadata(row)
        if name not in COARSE_REFINED:
            add_leaf(
                leaves, topology="rank/rank", stage="coarse-exact",
                name=name, row=row,
            )

    fine_trace = load_rows("projective_l055_trace_fine")
    fine_exact = load_rows("projective_l055_exact_fine")
    require(set(fine_trace) == set(fine_exact), "fine trace/exact names differ")
    require(
        actual_xy_bounds(fine_trace.values()) == expected_fine_bounds(),
        "fine boxes do not exactly partition the two unresolved coarse boxes",
    )
    for name, row in fine_exact.items():
        check_metadata(row)
        if name not in FINE_ANGLE_PARENTS:
            add_leaf(
                leaves, topology="rank/rank", stage="fine-exact",
                name=name, row=row,
            )
    require(
        {name for name, row in fine_exact.items() if float(row["dual_scaled"]) > 1.0}
        == FINE_ANGLE_PARENTS,
        "the recorded fine angle parents changed",
    )

    angle_union = load_rows("projective_l055_angle_union")
    require(
        set(angle_union) == {f"a{i}{j}" for i in range(4) for j in range(4)},
        "global angle union is incomplete",
    )
    for i in range(4):
        for j in range(4):
            name = f"a{i}{j}"
            row = angle_union[name]
            check_metadata(row)
            check_selected_bounds(row, expected_angle_bounds(i, j), f"union/{name}")
            check_selected_bounds(
                row,
                {"x": (0.875, 0.95), "y": (0.025, 0.1)},
                f"union/{name}",
            )
            if (i, j) not in LOW_ANGLE_INDICES:
                add_leaf(
                    leaves, topology="rank/rank", stage="angle-union-exact",
                    name=name, row=row,
                )

    local = load_rows("projective_l055_angle_local")
    expected_local = {
        f"{parent}_a{i}{j}"
        for parent in FINE_ANGLE_PARENTS - {PROBE_PARENT}
        for i, j in LOW_ANGLE_INDICES
    }
    require(set(local) == expected_local, "local low-angle boxes are incomplete")

    probe = load_rows("projective_l055_angle_probe")
    require(
        set(probe) == {f"a{i}{j}" for i in range(4) for j in range(4)},
        "physical-box angle probe is incomplete",
    )

    refined_parents: dict[str, dict[str, object]] = {}
    for parent in sorted(FINE_ANGLE_PARENTS):
        source = probe if parent == PROBE_PARENT else local
        parent_xy = bounds_dict(fine_exact[parent])
        for i, j in LOW_ANGLE_INDICES:
            name = f"a{i}{j}" if parent == PROBE_PARENT else f"{parent}_a{i}{j}"
            row = source[name]
            check_metadata(row)
            check_selected_bounds(row, parent_xy, f"low-angle/{name}")
            check_selected_bounds(row, expected_angle_bounds(i, j), f"low-angle/{name}")
            canonical = f"{parent}_a{i}{j}"
            if float(row["dual_scaled"]) <= 1.0:
                add_leaf(
                    leaves, topology="rank/rank", stage="low-angle-exact",
                    name=canonical, row=row,
                )
            else:
                refined_parents[canonical] = row

    refined = load_rows("projective_l055_angle_refined")
    expected_refined = {
        f"{parent}_r{i}{j}"
        for parent in refined_parents
        for i in range(2)
        for j in range(2)
    }
    require(set(refined) == expected_refined, "refined angle partition is incomplete")
    for parent, parent_row in refined_parents.items():
        parent_bounds = bounds_dict(parent_row)
        for i in range(2):
            for j in range(2):
                name = f"{parent}_r{i}{j}"
                row = refined[name]
                check_metadata(row)
                expected = dict(parent_bounds)
                for coordinate, index in (
                    ("first_angle_0", i), ("second_angle_0", j)
                ):
                    lower, upper = parent_bounds[coordinate]
                    middle = (lower + upper) / 2.0
                    expected[coordinate] = (
                        (lower, middle) if index == 0 else (middle, upper)
                    )
                check_selected_bounds(row, expected, f"refined/{name}")
                add_leaf(
                    leaves, topology="rank/rank", stage="refined-angle-exact",
                    name=name, row=row,
                )
    return leaves


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    require(CERTIFIED_LEVEL >= OLD_LEVEL, "rescaling requires L1 >= L0")
    require(OLD_LEVEL > WEIGHT, "all secular denominators must be positive")
    leaves = singleton_leaves() + rank_rank_leaves()
    maximum = max(float(row["rescaled_dual_upper"]) for row in leaves)
    worst = max(leaves, key=lambda row: float(row["rescaled_dual_upper"]))
    require(maximum <= 1.0 + TOL, "projective cover failed after rescaling")

    topologies: dict[str, dict[str, object]] = {}
    for topology in sorted({str(row["topology"]) for row in leaves}):
        selected = [row for row in leaves if row["topology"] == topology]
        topologies[topology] = {
            "leaf_count": len(selected),
            "maximum_rescaled_dual_upper": max(
                float(row["rescaled_dual_upper"]) for row in selected
            ),
            "total_nodes": sum(int(row["nodes"]) for row in selected),
            "total_solving_time": sum(
                float(row["solving_time"]) for row in selected
            ),
        }

    manifest = {
        "weight": WEIGHT,
        "archived_solver_level": OLD_LEVEL,
        "certified_projective_upper": CERTIFIED_LEVEL,
        "explicit_projective_lower": PHYSICAL_LOWER,
        "projective_interval_width": CERTIFIED_LEVEL - PHYSICAL_LOWER,
        "rescaling_factor": secular_rescaling_factor(OLD_LEVEL, CERTIFIED_LEVEL),
        "all_leaf_duals_at_most_one_after_rescaling": True,
        "leaf_count": len(leaves),
        "maximum_rescaled_dual_upper": maximum,
        "worst_leaf": {
            key: worst[key]
            for key in (
                "topology", "stage", "name", "archived_dual_scaled",
                "rescaled_dual_upper",
            )
        },
        "topologies": topologies,
        "certificate_class": (
            "finite numerical SCIP cover with exact analytic monotonic "
            "post-processing; conditional on the recorded solver dual bounds"
        ),
        "scope": (
            "binary-projective terminal sector of the reduced one-way qubit "
            "problem at lambda=0.55"
        ),
    }
    rendered = json.dumps(manifest, indent=2) + "\n"
    print(rendered)
    if args.output is not None:
        output = args.output
        if not output.is_absolute():
            output = ROOT / output
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

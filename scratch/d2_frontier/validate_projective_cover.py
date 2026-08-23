"""Validate the finite SCIP cover of the projective frontier at lambda=0.6.

The cover is hierarchical.  Cheap trace-only relaxations certify most boxes;
an exact secular model replaces every unresolved box, and difficult exact
boxes are partitioned in the two rank-split angles.  This validator checks
that every parent either has dual bound at most one or is completely replaced
by the expected child grid.  It then emits a compact manifest of the leaf
certificates.

The certificate is numerical: its meaning is conditional on SCIP's spatial
branch-and-bound implementation and the tolerances recorded in the individual
JSON files.  It is not an interval-arithmetic proof independent of SCIP.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
ANGLE_ONE = range(4)
ANGLE_TWO = range(3)


def payloads(directory: str) -> dict[str, dict[str, object]]:
    path = ROOT / directory
    rows = {
        item.stem: json.loads(item.read_text(encoding="utf-8"))
        for item in path.glob("*.json")
    }
    if not rows:
        raise RuntimeError(f"no certificate files found in {path}")
    return rows


def certified(row: dict[str, object]) -> bool:
    return float(row["dual_scaled"]) <= 1.0


def add_leaf(
    leaves: list[dict[str, object]],
    topology: str,
    stage: str,
    name: str,
    row: dict[str, object],
) -> None:
    if not certified(row):
        raise RuntimeError(f"uncertified leaf {stage}/{name}")
    leaves.append(
        {
            "topology": topology,
            "stage": stage,
            "name": name,
            "dual_scaled": float(row["dual_scaled"]),
            "status": str(row["status"]),
            "nodes": int(row["nodes"]),
            "solving_time": float(row["solving_time"]),
        }
    )


def rank_rank_leaves() -> list[dict[str, object]]:
    leaves: list[dict[str, object]] = []
    coarse = payloads("projective_trace_boxes")
    coarse_bad = []
    for name, row in coarse.items():
        if certified(row):
            add_leaf(leaves, "rank/rank", "coarse-trace", name, row)
        else:
            coarse_bad.append(name)
    if coarse_bad != ["x910_y01"]:
        raise RuntimeError(f"unexpected unresolved coarse boxes: {coarse_bad}")

    fine_trace = payloads("projective_trace_fine")
    full_fine = payloads("projective_full_fine")
    angle_parents = []
    for name, row in fine_trace.items():
        if certified(row):
            add_leaf(leaves, "rank/rank", "fine-trace", name, row)
            continue
        exact = full_fine.get(name)
        if exact is None:
            raise RuntimeError(f"missing exact fine box {name}")
        if certified(exact):
            add_leaf(leaves, "rank/rank", "fine-exact", name, exact)
        else:
            angle_parents.append(name)

    angle_trace = payloads("projective_trace_angle")
    angle_exact = payloads("projective_full_angle")
    refined_parents = []
    for parent in angle_parents:
        expected = [f"{parent}_a{i}{j}" for i in ANGLE_ONE for j in ANGLE_ONE]
        for name in expected:
            row = angle_trace.get(name)
            if row is None:
                raise RuntimeError(f"missing first-angle trace box {name}")
            if certified(row):
                add_leaf(leaves, "rank/rank", "angle-trace", name, row)
                continue
            exact = angle_exact.get(name)
            if exact is None:
                raise RuntimeError(f"missing first-angle exact box {name}")
            if certified(exact):
                add_leaf(leaves, "rank/rank", "angle-exact", name, exact)
            else:
                if not name.endswith("_a00"):
                    raise RuntimeError(f"unexpected refined angle parent {name}")
                refined_parents.append(parent)

    angle_two_trace = payloads("projective_trace_angle2")
    angle_two_exact = payloads("projective_full_angle2")
    angle_two_long = payloads("projective_full_angle2_long")
    angle_two_extended = payloads("projective_full_angle2_extended")
    for parent in refined_parents:
        expected = [f"{parent}_a{i}{j}" for i in ANGLE_TWO for j in ANGLE_TWO]
        for name in expected:
            row = angle_two_trace.get(name)
            if row is None:
                raise RuntimeError(f"missing refined trace box {name}")
            if certified(row):
                add_leaf(leaves, "rank/rank", "refined-trace", name, row)
                continue
            exact = angle_two_exact.get(name)
            if exact is None:
                raise RuntimeError(f"missing refined exact box {name}")
            if certified(exact):
                add_leaf(leaves, "rank/rank", "refined-exact", name, exact)
                continue
            long_row = angle_two_long.get(name)
            if long_row is None:
                raise RuntimeError(f"missing long exact box {name}")
            if certified(long_row):
                add_leaf(leaves, "rank/rank", "refined-long", name, long_row)
                continue
            extended = angle_two_extended.get(name)
            if extended is None:
                raise RuntimeError(f"missing extended exact box {name}")
            add_leaf(leaves, "rank/rank", "refined-extended", name, extended)
    return leaves


def endpoint_rank_leaves() -> list[dict[str, object]]:
    leaves: list[dict[str, object]] = []
    coarse = payloads("projective_er_trace_boxes")
    bad = []
    for name, row in coarse.items():
        if certified(row):
            add_leaf(leaves, "endpoint/rank", "coarse-trace", name, row)
        else:
            bad.append(name)
    if bad != ["x910_y01"]:
        raise RuntimeError(f"unexpected endpoint/rank coarse boxes: {bad}")
    fine = payloads("projective_er_trace_fine")
    exact = payloads("projective_er_full_fine")
    for name, row in fine.items():
        if certified(row):
            add_leaf(leaves, "endpoint/rank", "fine-trace", name, row)
        else:
            replacement = exact.get(name)
            if replacement is None:
                raise RuntimeError(f"missing endpoint/rank exact box {name}")
            add_leaf(leaves, "endpoint/rank", "fine-exact", name, replacement)
    return leaves


def singleton_leaves() -> list[dict[str, object]]:
    leaves: list[dict[str, object]] = []
    for topology, filename in (
        ("endpoint/endpoint", "projective_trace_endpoint_endpoint_l060.json"),
        ("rank/endpoint", "projective_trace_rank_endpoint_l060.json"),
    ):
        row = json.loads((ROOT / filename).read_text(encoding="utf-8"))
        add_leaf(leaves, topology, "global-trace", Path(filename).stem, row)
    return leaves


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    candidate = json.loads(
        (ROOT / "reduced_four_effect_l060.json").read_text(encoding="utf-8")
    )
    leaves = singleton_leaves() + endpoint_rank_leaves() + rank_rank_leaves()
    levels = {float(row["level"]) for directory in (
        "projective_trace_boxes",
        "projective_er_trace_boxes",
    ) for row in payloads(directory).values()}
    if levels != {0.76591}:
        raise RuntimeError(f"inconsistent support levels: {levels}")
    if any(float(row["dual_scaled"]) > 1.0 for row in leaves):
        raise RuntimeError("cover contains a leaf above the secular threshold")

    by_topology: dict[str, dict[str, object]] = {}
    for topology in sorted({str(row["topology"]) for row in leaves}):
        selected = [row for row in leaves if row["topology"] == topology]
        by_topology[topology] = {
            "leaf_count": len(selected),
            "maximum_dual_scaled": max(
                float(row["dual_scaled"]) for row in selected
            ),
            "total_nodes": sum(int(row["nodes"]) for row in selected),
            "total_solving_time": sum(
                float(row["solving_time"]) for row in selected
            ),
        }
    level = next(iter(levels))
    lower = float(candidate["score"])
    manifest = {
        "weight": 0.6,
        "certified_projective_upper": level,
        "explicit_four_effect_lower": lower,
        "projective_interval_width": level - lower,
        "all_leaf_duals_at_most_one": True,
        "leaf_count": len(leaves),
        "maximum_leaf_dual_scaled": max(
            float(row["dual_scaled"]) for row in leaves
        ),
        "topologies": by_topology,
        "scope": (
            "binary-projective terminal sector of the two-block rank-two "
            "relaxation; numerical SCIP certificate, not the unrestricted "
            "three/four-active terminal POVM sector"
        ),
    }
    rendered = json.dumps(manifest, indent=2) + "\n"
    print(rendered)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

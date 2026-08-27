"""Validate a solver-independent projective tangent cover.

The validator joins the rank/rank cover with the three complementary
topologies, reconstructs every expected geometric cell, and rejects missing,
duplicate, extra, incomplete, or parameter-inconsistent records.  It checks
the finite proof summaries; it does not replay every interval split because
the current artifacts intentionally store compact deterministic summaries.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any

import certify_rank_rank_tangent_cover as rank_rank
import certify_remaining_projective_tangent_cover as remaining
from projective_tangent_interval_certificate import Interval


ROOT = Path(__file__).resolve().parent
REPOSITORY_ROOT = ROOT.parents[1]
KERNEL = "outward-expanded IEEE-754 binary64 intervals"
SUMMARY_SCHEMA = "carmenq.projective-tangent-global-summary.v1"


def expected_rank_rank() -> dict[str, dict[str, Any]]:
    return {
        str(cell["certificate_name"]): cell
        for leaf in rank_rank.geometric_leaves()
        for cell in rank_rank.expand_full_angles(leaf)
    }


def expected_remaining() -> dict[str, dict[str, Any]]:
    return {
        str(cell["certificate_name"]): cell
        for cell in remaining.topology_cells()
    }


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def portable_path(path: Path) -> str:
    try:
        return path.resolve().relative_to(REPOSITORY_ROOT).as_posix()
    except ValueError as error:
        raise RuntimeError(f"artifact lies outside the repository: {path}") from error


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_payload(
    payload: dict[str, Any],
    expected: dict[str, dict[str, Any]],
    name_field: str,
    weight: str,
    requested_level: Fraction,
) -> dict[str, Any]:
    require(payload.get("weight") == weight, "wrong projective weight")
    try:
        certified_level = Fraction(str(payload.get("level")))
    except (ValueError, ZeroDivisionError) as exc:
        raise RuntimeError("malformed projective level") from exc
    require(
        certified_level <= requested_level,
        "projective component is certified only at a weaker level",
    )
    level = str(certified_level)
    require(bool(payload.get("run_complete", True)), "cover run is incomplete")
    require(bool(payload.get("all_cells_complete")), "cover has open cells")
    records = payload.get("cells", [])
    by_name: dict[str, dict[str, Any]] = {}
    for record in records:
        name = str(record[name_field])
        require(name not in by_name, f"duplicate cell {name}")
        by_name[name] = record
    require(set(by_name) == set(expected), "geometric cell set is not exact")
    methods: dict[str, int] = {}
    splits = 0
    for name, record in by_name.items():
        expected_cell = expected[name]
        require(record["bounds"] == {
            key: list(value) for key, value in expected_cell["bounds"].items()
        }, f"bounds changed for {name}")
        certificate = record["certificate"]
        require(bool(certificate["complete"]), f"open certificate {name}")
        require(int(certificate["boxes_remaining"]) == 0, f"open boxes in {name}")
        require(certificate["weight"] == weight, f"wrong weight in {name}")
        require(certificate["level"] == level, f"wrong level in {name}")
        require(certificate["proof_kernel"] == KERNEL, f"wrong kernel in {name}")
        require(not certificate["trusted_optimizers"], f"trusted optimizer in {name}")
        splits += int(certificate["boxes_split"])
        for method, count in certificate["closed_methods"].items():
            methods[method] = methods.get(method, 0) + int(count)
    require(splits == int(payload["total_boxes_split"]), "split total mismatch")
    return {
        "cell_count": len(records),
        "boxes_split": splits,
        "closed_methods": methods,
        "certified_level": level,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--weight", default="0.6")
    parser.add_argument("--level", default="0.766")
    parser.add_argument(
        "--rank-rank",
        type=Path,
        default=ROOT / "rank_rank_tangent_full_l060_L0766.json",
    )
    parser.add_argument(
        "--remaining",
        type=Path,
        default=ROOT / "remaining_projective_tangent_full_l060_L076591.json",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    weight = str(Fraction(args.weight))
    requested_level = Fraction(args.level)
    level = str(requested_level)
    if not Fraction(0) < Fraction(weight) < requested_level < Fraction(1):
        parser.error("require 0 < weight < level < 1")

    rank_payload = json.loads(args.rank_rank.read_text(encoding="utf-8"))
    remaining_payload = json.loads(args.remaining.read_text(encoding="utf-8"))
    require(
        rank_payload.get("schema") == rank_rank.SCHEMA,
        "wrong rank/rank projective schema",
    )
    require(
        remaining_payload.get("schema") == remaining.SCHEMA,
        "wrong remaining-topology projective schema",
    )
    rank_summary = validate_payload(
        rank_payload,
        expected_rank_rank(),
        "certificate_name",
        weight,
        requested_level,
    )
    remaining_summary = validate_payload(
        remaining_payload,
        expected_remaining(),
        "certificate_name",
        weight,
        requested_level,
    )
    angular_endpoint = Interval.decimal(str(rank_rank.ANGLE_EDGES[-1]))
    require(
        angular_endpoint.hi * angular_endpoint.hi >= 0.5,
        "angular partition misses 1/sqrt(2)",
    )
    summary = {
        "schema": SUMMARY_SCHEMA,
        "weight": weight,
        "level": level,
        "certified_projective_upper": float(Fraction(level)),
        "all_four_topologies_complete": True,
        "topology_cell_count": (
            rank_summary["cell_count"] + remaining_summary["cell_count"]
        ),
        "source_artifacts": {
            "rank_rank": {
                "path": portable_path(args.rank_rank),
                "sha256": sha256(args.rank_rank),
            },
            "remaining_topologies": {
                "path": portable_path(args.remaining),
                "sha256": sha256(args.remaining),
            },
        },
        "rank_rank": rank_summary,
        "remaining_topologies": remaining_summary,
        "proof_kernel": KERNEL,
        "solver_independence": (
            "SCIP and archived dual values are not used in any cell inequality"
        ),
        "trust_boundary": (
            "Python IEEE-754 arithmetic and math.nextafter; compact summaries "
            "are deterministically regenerated rather than storing split trees"
        ),
        "scope": (
            "binary-projective terminal sector of the two-block rank-two "
            f"relaxation at lambda={weight}"
        ),
    }
    rendered = json.dumps(summary, indent=2) + "\n"
    print(rendered, end="")
    if args.output is not None:
        output = args.output if args.output.is_absolute() else ROOT / args.output
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

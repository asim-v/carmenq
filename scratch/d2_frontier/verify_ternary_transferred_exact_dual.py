"""Solver-free ternary replay under changed projective support premises.

The source artifacts store cone-feasible dual vectors found for one pair of
named projective lines.  A line change modifies the conic data, so stored
objectives and canonical hashes cannot simply be reused.  This verifier
rebuilds every selected cell with the requested new lines, repairs the stored
dual against the unchanged cones, and recomputes the objective and
stationarity correction in exact rational arithmetic.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any

import ternary_socp_exact_dual_cover as cover


ROOT = Path(__file__).resolve().parent
REPOSITORY_ROOT = ROOT.parents[1]
SCHEMA = "carmenq.ternary-socp-transferred-exact-dual-verification.v1"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def portable_path(path: Path) -> str:
    """Return a repository-relative POSIX path or reject external input."""
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


def verify_artifact(
    artifact_path: Path,
    baseline_target: Fraction,
    baseline_premises: dict[str, str],
    line_055: Fraction,
    line_060: Fraction,
    target: Fraction,
) -> tuple[dict[str, Any], list[str]]:
    artifact = json.loads(artifact_path.read_text(encoding="utf-8"))
    require(
        artifact.get("schema")
        == "carmenq.ternary-socp-clustered-exact-dual-cover.v1",
        f"wrong source schema in {artifact_path}",
    )
    require(artifact.get("support_weight") == "3/5", "wrong support weight")
    require(
        Fraction(artifact["target"]) == baseline_target,
        "wrong baseline ternary target",
    )
    require(
        artifact.get("named_projective_premises") == baseline_premises,
        "wrong baseline projective premises",
    )

    cover.LINE_055_UPPER = line_055
    cover.LINE_060_UPPER = line_060
    certifier = cover.ReusableExactCertifier(target)
    vectors = [
        cover.decode_candidate(candidate) for candidate in artifact["candidates"]
    ]
    keys: list[str] = []
    maximum: Fraction | None = None
    finite_count = 0
    empty_count = 0
    for index, cell in enumerate(artifact["cells"], start=1):
        key = cover.box_key(cell["box"])
        require(key not in keys, f"duplicate cell inside {artifact_path}")
        keys.append(key)
        if cell.get("source_status") == "domain_empty":
            require(cover.exact_domain_empty(cell["box"]), "invalid empty cell")
            require(cell.get("closed") is True, "empty cell is not closed")
            empty_count += 1
            continue
        candidate_id = int(cell["candidate"])
        require(0 <= candidate_id < len(vectors), "candidate index out of range")
        data, enclosure = certifier.data(cell["box"])
        require(
            cover.compact_enclosure(enclosure) == cell["inellipse"],
            "inellipse audit summary mismatch",
        )
        dual, _ = cover.repair_dual_cones(vectors[candidate_id], data["dims"])
        upper, _, _ = certifier.exact_upper(data, dual)
        require(upper <= target, "transferred cell exceeds target")
        maximum = upper if maximum is None else max(maximum, upper)
        finite_count += 1
        if index % 250 == 0:
            print(
                json.dumps(
                    {
                        "artifact": artifact_path.name,
                        "verified": index,
                        "count": len(artifact["cells"]),
                    }
                ),
                flush=True,
            )
    require(
        int(artifact["processed_cell_count"]) == len(keys),
        "source processed-cell count mismatch",
    )
    require(
        int(artifact["closed_cell_count"]) == len(keys),
        "source closed-cell count mismatch",
    )
    require(maximum is not None, "artifact contains no finite cell")
    return (
        {
            "path": portable_path(artifact_path),
            "sha256": sha256(artifact_path),
            "cell_count": len(keys),
            "finite_cell_count": finite_count,
            "domain_empty_cell_count": empty_count,
            "candidate_count": len(vectors),
            "maximum_certified_upper_fraction": [
                maximum.numerator,
                maximum.denominator,
            ],
            "maximum_certified_upper": cover.fraction_decimal(maximum),
            "verified_without_solver": True,
        },
        keys,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("artifacts", nargs="+", type=Path)
    parser.add_argument("--cover", type=Path, default=cover.DEFAULT_COVER)
    parser.add_argument("--baseline-target", default="0.76643")
    parser.add_argument("--baseline-line-055", default="0.7573")
    parser.add_argument("--baseline-line-060", default="0.76591")
    parser.add_argument("--line-055", default="0.7573")
    parser.add_argument("--line-060", default="0.766")
    parser.add_argument("--target", default="0.76652")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--allow-partial", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be positive")

    baseline_target = Fraction(args.baseline_target)
    line_055 = Fraction(args.line_055)
    line_060 = Fraction(args.line_060)
    target = Fraction(args.target)
    baseline_premises = {
        "11/20": str(Fraction(args.baseline_line_055)),
        "3/5": str(Fraction(args.baseline_line_060)),
    }
    source = json.loads(args.cover.read_text(encoding="utf-8"))
    expected = {
        cover.box_key(item["box"]) for item in source["leaves"]
    }
    require(len(expected) == len(source["leaves"]), "source cover has duplicates")

    summaries: list[dict[str, Any]] = []
    union: set[str] = set()
    original_055 = cover.LINE_055_UPPER
    original_060 = cover.LINE_060_UPPER
    try:
        if args.workers == 1:
            completed = [
                verify_artifact(
                    path,
                    baseline_target,
                    baseline_premises,
                    line_055,
                    line_060,
                    target,
                )
                for path in args.artifacts
            ]
        else:
            completed = []
            with ProcessPoolExecutor(max_workers=args.workers) as executor:
                futures = {
                    executor.submit(
                        verify_artifact,
                        path,
                        baseline_target,
                        baseline_premises,
                        line_055,
                        line_060,
                        target,
                    ): path
                    for path in args.artifacts
                }
                for future in as_completed(futures):
                    completed.append(future.result())
        for summary, keys in completed:
            overlap = union.intersection(keys)
            require(not overlap, f"artifacts overlap on {len(overlap)} cells")
            union.update(keys)
            summaries.append(summary)
    finally:
        cover.LINE_055_UPPER = original_055
        cover.LINE_060_UPPER = original_060

    require(union.issubset(expected), "artifact contains an unknown source cell")
    if not args.allow_partial:
        require(union == expected, "artifacts do not cover every source leaf")
    maxima = [
        Fraction(*summary["maximum_certified_upper_fraction"])
        for summary in summaries
    ]
    maximum = max(maxima)
    result = {
        "schema": SCHEMA,
        "baseline_target": str(baseline_target),
        "baseline_projective_premises": baseline_premises,
        "target": str(target),
        "named_projective_premises": {
            "11/20": str(line_055),
            "3/5": str(line_060),
        },
        "source_leaf_count": len(expected),
        "verified_cell_count": len(union),
        "full_source_cover_verified": union == expected,
        "all_cells_at_most_target": maximum <= target,
        "maximum_certified_upper_fraction": [
            maximum.numerator,
            maximum.denominator,
        ],
        "maximum_certified_upper": cover.fraction_decimal(maximum),
        "artifacts": sorted(summaries, key=lambda row: row["path"]),
        "trusted_optimizers": [],
        "verified_without_solver": True,
    }
    rendered = json.dumps(result, indent=2) + "\n"
    print(rendered, end="")
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

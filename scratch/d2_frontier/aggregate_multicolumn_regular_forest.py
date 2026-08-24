"""Aggregate and audit disjoint regular-forest manifests."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from adaptive_multicolumn_regular_forest import open_branch_orbits
from audit_adaptive_multicolumn_certificate import audit_certificate


def aggregate(
    base_cover_path: Path,
    manifest_paths: list[Path],
    target: float,
    require_complete: bool,
) -> dict[str, object]:
    base_cover = json.loads(base_cover_path.read_text(encoding="utf-8"))
    expected = open_branch_orbits(base_cover, target)
    expected_by_key = {
        tuple(item["representative"]): (index, item)
        for index, item in enumerate(expected)
    }

    rows_by_key: dict[tuple[object, ...], dict[str, object]] = {}
    for manifest_path in manifest_paths:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if float(manifest["target"]) != target:
            raise ValueError(f"target mismatch in {manifest_path}")
        for row in manifest["orbits"]:
            key = tuple(row["representative"])
            if key not in expected_by_key:
                raise ValueError(f"unexpected orbit {key} in {manifest_path}")
            if key in rows_by_key:
                raise ValueError(f"duplicate orbit {key}")
            rows_by_key[key] = row

    audited_rows = []
    totals = {
        "expanded_nodes": 0,
        "closed_leaves": 0,
        "source_closed_leaves": 0,
        "infeasible_leaves": 0,
    }
    maximum_depth = 0
    maximum_terminal_bound = float("-inf")
    branching_factors: set[int] = set()
    for key, row in rows_by_key.items():
        if not bool(row["complete"]):
            continue
        certificate_path = Path(str(row["certificate"]))
        certificate = json.loads(certificate_path.read_text(encoding="utf-8"))
        certificate_key = (
            int(certificate["plane_index"]),
            int(certificate["sphere_index"]),
            str(certificate.get("pair_branch", "bloch")),
            int(certificate["pair_cap_index"]),
        )
        if certificate_key != key:
            raise ValueError(f"certificate base key mismatch for {key}")
        audit = audit_certificate(certificate)
        if not bool(audit["certificate_complete"]):
            raise ValueError(f"certificate for {key} is not complete")
        totals["expanded_nodes"] += int(audit["expanded_nodes"])
        totals["closed_leaves"] += int(audit["closed_leaves"])
        totals["source_closed_leaves"] += int(audit["source_closed_nodes"])
        totals["infeasible_leaves"] += int(audit["infeasible_source_nodes"])
        branching_factors.add(int(audit["branches_per_expansion"]))
        maximum_depth = max(maximum_depth, int(audit["maximum_depth"]))
        finite_terminal = []
        for node in certificate["expanded_nodes"]:
            for branch in node["branches"]:
                if branch["child_id"] is not None:
                    continue
                if branch["bound"] is None:
                    totals["infeasible_leaves"] += 1
                else:
                    finite_terminal.append(float(branch["bound"]))
        source_closed = certificate.get(
            "source_closed_nodes", certificate.get("infeasible_source_nodes", [])
        )
        finite_terminal.extend(
            float(node["source_bound"])
            for node in source_closed
            if node.get("source_bound") is not None
        )
        if finite_terminal:
            maximum_terminal_bound = max(maximum_terminal_bound, *finite_terminal)
        expected_index, expected_orbit = expected_by_key[key]
        audited_rows.append(
            {
                "global_index": expected_index,
                **expected_orbit,
                "certificate": certificate_path.as_posix(),
                "expansions": int(audit["expanded_nodes"]),
                "closed_leaves": int(audit["closed_leaves"]),
                "maximum_depth": int(audit["maximum_depth"]),
            }
        )

    missing = [
        {"global_index": index, **item}
        for index, item in enumerate(expected)
        if tuple(item["representative"]) not in rows_by_key
        or not bool(rows_by_key[tuple(item["representative"])]["complete"])
    ]
    complete = not missing
    if require_complete and not complete:
        raise ValueError(f"the regular forest still has {len(missing)} open or missing orbits")
    audited_rows.sort(key=lambda item: int(item["global_index"]))
    if len(branching_factors) > 1:
        raise ValueError("the forest certificates use different branching factors")
    branching_factor = next(iter(branching_factors), None)
    expected_closed_leaves = len(audited_rows)
    if branching_factor is not None:
        expected_closed_leaves += (branching_factor - 1) * totals["expanded_nodes"]
    if totals["closed_leaves"] != expected_closed_leaves:
        raise ValueError("the aggregated forest violates the full-tree leaf identity")
    return {
        "scope": "audited global bbb Fourier/pair multicolumn forest",
        "target": target,
        "solver_conditional": True,
        "complex_conjugation_quotient": True,
        "base_open_branches": sum(len(item["orbit"]) for item in expected),
        "expected_symmetry_orbits": len(expected),
        "audited_complete_orbits": len(audited_rows),
        "missing_or_open_orbits": len(missing),
        "certificate_complete": complete,
        "exhaustive_branches_per_expansion": branching_factor,
        "leaf_identity_verified": True,
        **totals,
        "maximum_depth": maximum_depth,
        "maximum_finite_terminal_bound": (
            maximum_terminal_bound if maximum_terminal_bound != float("-inf") else None
        ),
        "orbits": audited_rows,
        "missing": missing,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("base_cover", type=Path)
    parser.add_argument("manifests", type=Path, nargs="+")
    parser.add_argument("--target", type=float, default=0.758)
    parser.add_argument("--allow-incomplete", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = aggregate(
        args.base_cover,
        args.manifests,
        args.target,
        require_complete=not args.allow_incomplete,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {key: result[key] for key in (
                "certificate_complete",
                "audited_complete_orbits",
                "missing_or_open_orbits",
                "expanded_nodes",
                "closed_leaves",
                "maximum_finite_terminal_bound",
            )},
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

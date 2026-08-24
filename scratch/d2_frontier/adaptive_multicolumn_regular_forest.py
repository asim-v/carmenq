"""Run adaptive multicolumn trees over every open regular base branch."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from adaptive_multicolumn_branch_tree import run_tree
from fourier_behavior_cap_cover import cube_face_caps


BaseKey = tuple[int, int, str, int]


def _reflection_map(normals: list[np.ndarray]) -> dict[int, int]:
    result = {}
    for index, normal in enumerate(normals):
        reflected = normal * np.asarray([1.0, -1.0, 1.0])
        distances = [np.linalg.norm(reflected - candidate) for candidate in normals]
        partner = int(np.argmin(distances))
        if distances[partner] > 1e-12:
            raise ValueError("the cap family is not closed under complex conjugation")
        result[index] = partner
    return result


def open_branch_orbits(
    base_cover: dict[str, object], target: float
) -> list[dict[str, object]]:
    """Quotient open base branches by exact Bloch-y reflection symmetry."""
    grid = int(base_cover["cube_face_grid"])
    sphere_map = _reflection_map(
        [np.asarray(normal, dtype=float) for normal, _ in cube_face_caps(grid)]
    )
    pair_normals = [
        np.eye(3)[axis] * sign
        for axis in range(3)
        for sign in (-1.0, 1.0)
    ]
    pair_map = _reflection_map(pair_normals)

    open_rows: dict[BaseKey, float] = {}
    for cell in base_cover["cells"]:
        for branch in cell["branches"]:
            bound = float(branch["bound"])
            if bound < target:
                continue
            pair_branch = str(branch["pair_branch"])
            pair_cap = -1 if branch["pair_cap"] is None else int(branch["pair_cap"])
            key = (int(cell["plane"]), int(cell["sphere"]), pair_branch, pair_cap)
            open_rows[key] = bound

    representatives: list[dict[str, object]] = []
    visited: set[BaseKey] = set()
    for key, bound in sorted(open_rows.items(), key=lambda item: item[1], reverse=True):
        if key in visited:
            continue
        plane, sphere, pair_branch, pair_cap = key
        partner = (
            plane,
            sphere_map[sphere],
            pair_branch,
            pair_map[pair_cap] if pair_branch == "bloch" else pair_cap,
        )
        if partner not in open_rows:
            raise ValueError(f"open branch {key} has no conjugate partner")
        if abs(open_rows[partner] - bound) > 2e-7:
            raise ValueError(f"conjugate bounds disagree for {key} and {partner}")
        orbit = sorted({key, partner})
        representative = orbit[0]
        visited.update(orbit)
        representatives.append(
            {
                "representative": list(representative),
                "orbit": [list(item) for item in orbit],
                "base_bound": max(open_rows[item] for item in orbit),
            }
        )
    representatives.sort(key=lambda item: float(item["base_bound"]), reverse=True)
    return representatives


def _certificate_name(key: BaseKey) -> str:
    plane, sphere, branch, cap = key
    return f"p{plane:03d}_s{sphere:03d}_{branch}_c{cap}.json"


def _base_key(payload: dict[str, object]) -> BaseKey:
    return (
        int(payload["plane_index"]),
        int(payload["sphere_index"]),
        str(payload.get("pair_branch", "bloch")),
        int(payload["pair_cap_index"]),
    )


def run_forest(args: argparse.Namespace) -> dict[str, object]:
    base_cover = json.loads(args.base_cover.read_text(encoding="utf-8"))
    orbits = open_branch_orbits(base_cover, args.target)
    reused: dict[BaseKey, tuple[Path, dict[str, object]]] = {}
    for path in args.reuse_certificate:
        payload = json.loads(path.read_text(encoding="utf-8"))
        reused[_base_key(payload)] = (path, payload)

    args.certificate_directory.mkdir(parents=True, exist_ok=True)
    rows = []
    selected_pairs = list(enumerate(orbits))[args.start_root :: args.root_stride]
    if args.max_roots is not None:
        selected_pairs = selected_pairs[: args.max_roots]
    selected = [item[1] for item in selected_pairs]
    manifest = {
        "scope": "regular multicolumn forest modulo exact complex conjugation",
        "target": args.target,
        "base_open_branches": sum(len(item["orbit"]) for item in orbits),
        "symmetry_orbits": len(orbits),
        "start_root": args.start_root,
        "root_stride": args.root_stride,
        "selected_orbits": len(selected),
        "processed_orbits": 0,
        "complete_orbits": 0,
        "all_selected_complete": not selected,
        "orbits": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    for local_index, orbit in enumerate(selected):
        global_index = selected_pairs[local_index][0]
        key = tuple(orbit["representative"])
        if key in reused:
            certificate_path, result = reused[key]
        else:
            certificate_path = args.certificate_directory / _certificate_name(key)
            resume = certificate_path.exists()
            result = run_tree(
                None,
                args.target,
                args.max_expansions,
                args.separator_samples,
                args.separator_starts,
                args.seed + global_index * 1009,
                int(base_cover["plane_cells"]),
                int(key[0]),
                int(base_cover["cube_face_grid"]),
                int(key[1]),
                int(key[3]),
                args.contraction_grid,
                certificate_path,
                resume,
                str(key[2]),
            )
        row = {
            **orbit,
            "global_index": global_index,
            "certificate": certificate_path.as_posix(),
            "complete": bool(result["certificate_complete"]),
            "expansions": int(result["expansion_count"]),
            "closed_leaves": int(result["closed_leaf_count"]),
            "open_leaves": int(result["open_leaf_count"]),
            "maximum_open_bound": result["maximum_open_bound"],
        }
        rows.append(row)
        manifest = {
            "scope": "regular multicolumn forest modulo exact complex conjugation",
            "target": args.target,
            "base_open_branches": sum(len(item["orbit"]) for item in orbits),
            "symmetry_orbits": len(orbits),
            "start_root": args.start_root,
            "root_stride": args.root_stride,
            "selected_orbits": len(selected),
            "processed_orbits": len(rows),
            "complete_orbits": sum(bool(item["complete"]) for item in rows),
            "all_selected_complete": all(bool(item["complete"]) for item in rows),
            "orbits": rows,
        }
        args.output.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        print(
            json.dumps(
                {
                    "processed": len(rows),
                    "selected": len(selected),
                    "key": key,
                    "complete": row["complete"],
                    "expansions": row["expansions"],
                }
            ),
            flush=True,
        )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("base_cover", type=Path)
    parser.add_argument("--target", type=float, default=0.758)
    parser.add_argument("--start-root", type=int, default=0)
    parser.add_argument("--root-stride", type=int, default=1)
    parser.add_argument("--max-roots", type=int)
    parser.add_argument("--max-expansions", type=int, default=100)
    parser.add_argument("--separator-samples", type=int, default=8_000)
    parser.add_argument("--separator-starts", type=int, default=8)
    parser.add_argument("--seed", type=int, default=260831)
    parser.add_argument("--contraction-grid", type=int, default=4)
    parser.add_argument("--reuse-certificate", type=Path, action="append", default=[])
    parser.add_argument("--certificate-directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run_forest(args)
    print(
        json.dumps(
            {
                "base_open_branches": result["base_open_branches"],
                "symmetry_orbits": result["symmetry_orbits"],
                "processed_orbits": result["processed_orbits"],
                "complete_orbits": result["complete_orbits"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

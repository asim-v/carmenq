"""Parallel finite-cap cover for the vector-active Fourier branches."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
from pathlib import Path

from fourier_branch_upper import solve_branch


def _solve_task(task: tuple[str, int | None, int | None, str]) -> dict[str, object]:
    code, plane, sphere, output_name = task
    output = Path(output_name)
    if output.exists():
        return json.loads(output.read_text(encoding="utf-8"))
    payload = solve_branch(code, "clarabel", plane, sphere)
    output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not 1 <= args.workers <= 3:
        raise ValueError("workers must lie between one and three")
    directory = args.output.parent
    directory.mkdir(parents=True, exist_ok=True)
    tasks = [
        (
            "bbb",
            plane,
            sphere,
            str(directory / f"fourier_branch_bbb_c{plane}_s{sphere}_l055.json"),
        )
        for plane in range(4)
        for sphere in range(6)
    ]
    rows: list[dict[str, object]] = []
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(_solve_task, task): task for task in tasks}
        for future in as_completed(futures):
            payload = future.result()
            row = {
                "plane_cap": payload["plane_cap"],
                "sphere_cap": payload["sphere_cap"],
                "bound": payload["bound"],
                "status": payload["status"],
            }
            rows.append(row)
            print(json.dumps(row), flush=True)
    rows.sort(key=lambda row: (int(row["plane_cap"]), int(row["sphere_cap"])))
    manifest = {
        "scope": "coarse 4-by-6 cap cover of the bbb spectral branch",
        "target": 0.758,
        "closed": max(float(row["bound"]) for row in rows) < 0.758,
        "maximum_bound": max(float(row["bound"]) for row in rows),
        "cells": rows,
    }
    args.output.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

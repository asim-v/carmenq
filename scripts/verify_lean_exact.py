"""Build the exact Lean certificate in memory-bounded batches.

The concrete certificate is deliberately split into many small modules.  A
plain library build may ask Lake to elaborate too many of those modules at
once, which is unnecessarily memory hungry.  This driver builds the shared
dual witness and independent column atoms first, then assembles the shared
data and checks the remaining proof shards in bounded parallel batches,
and finally builds the aggregate theorem.
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FORMAL = ROOT / "formal"
MODULE_DIR = FORMAL / "CarmenQExact"


def find_lake(explicit: str | None) -> str:
    if explicit:
        return explicit
    discovered = shutil.which("lake")
    if discovered:
        return discovered
    fallback = Path.home() / ".elan" / "bin" / (
        "lake.exe" if sys.platform == "win32" else "lake"
    )
    if fallback.is_file():
        return str(fallback)
    raise SystemExit("lake was not found; install elan or pass --lake PATH")


def generated_modules(prefix: str) -> list[str]:
    pattern = re.compile(rf"{re.escape(prefix)}\d+\.lean$")
    paths = sorted(path for path in MODULE_DIR.glob(f"{prefix}*.lean") if pattern.fullmatch(path.name))
    return [f"CarmenQExact.{path.stem}" for path in paths]


def batches(values: list[str], size: int):
    for start in range(0, len(values), size):
        yield values[start : start + size]


def build(lake: str, modules: list[str]) -> None:
    command = [lake, "build", *modules]
    print("+", " ".join(command), flush=True)
    subprocess.run(command, cwd=FORMAL, check=True)


def trust_scan() -> None:
    forbidden = re.compile(
        r"\b(sorry|sorryAx|admit|axioms?|native_decide|unsafe|opaque|extern|implemented_by)\b"
    )
    offenders: list[str] = []
    # Only first-level non-probe modules enter the production import graph.
    # Experimental probes are intentionally neither scanned nor published.
    production_modules = sorted(
        path
        for path in MODULE_DIR.glob("*.lean")
        if not path.name.startswith("Probe")
    )
    paths = [FORMAL / "CarmenQExact.lean", *production_modules]
    for path in paths:
        for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            # Scan comments too. This is intentionally conservative and prevents
            # comment layout from becoming a lexical bypass for the trust gate.
            if forbidden.search(line):
                offenders.append(f"{path.relative_to(ROOT)}:{line_number}: {line.strip()}")
    if offenders:
        raise SystemExit("forbidden Lean trust escape found:\n" + "\n".join(offenders))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lake", help="path to the Lake executable")
    parser.add_argument(
        "--workers",
        type=int,
        default=3,
        help="maximum lightweight arithmetic modules passed to Lake at once (default: 3)",
    )
    parser.add_argument(
        "--heavy-workers",
        type=int,
        default=1,
        help="maximum Data-importing proof modules passed to Lake at once (default: 1)",
    )
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be positive")
    if args.heavy_workers < 1:
        parser.error("--heavy-workers must be positive")

    lake = find_lake(args.lake)
    trust_scan()

    build(
        lake,
        [
            "CarmenQExact.Checker",
            "CarmenQExact.WeakDuality",
            "CarmenQExact.EncodedWeakDuality",
        ],
    )
    build(lake, ["CarmenQExact.Source15818DualData"])

    for prefix in ("Source15818Stationarity",):
        modules = generated_modules(prefix)
        if not modules:
            raise SystemExit(f"no generated shards found for {prefix}")
        print(f"checking {len(modules)} {prefix} shards", flush=True)
        for group in batches(modules, args.workers):
            build(lake, group)

    build(lake, ["CarmenQExact.Source15818Data"])
    build(
        lake,
        [
            "CarmenQExact.Source15818Dimensions",
            "CarmenQExact.Source15818Indices",
        ],
    )

    heavy_prefixes = (
        "Source15818StationarityBridge",
        "Source15818Nonnegative",
        "Source15818Soc",
    )
    for prefix in heavy_prefixes:
        modules = generated_modules(prefix)
        if not modules:
            raise SystemExit(f"no generated shards found for {prefix}")
        print(f"checking {len(modules)} {prefix} shards", flush=True)
        for group in batches(modules, args.heavy_workers):
            build(lake, group)

    build(
        lake,
        [
            "CarmenQExact.Source15818Stationarity",
            "CarmenQExact.Source15818Nonnegative",
            "CarmenQExact.Source15818Soc",
        ],
    )
    build(
        lake,
        [
            "CarmenQExact.Source15818DualCone",
            "CarmenQExact.Source15818Upper",
        ],
    )
    build(lake, ["CarmenQExact.Source15818Exact"])
    build(lake, ["CarmenQExact"])
    print("exact Lean certificate verified", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Command-line interface for rebuilding the complete simulation artifact."""

from __future__ import annotations

import argparse
from pathlib import Path

from .experiments import generate_all


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Regenerate exact density-matrix tables, deterministic finite-shot "
            "samples, reference states, and publication figures."
        )
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path.cwd(),
        help="Repository root receiving data/ and figures/ (default: current directory).",
    )
    parser.add_argument(
        "--seed", type=int, default=20260812, help="Seed for finite-shot and challenge samples."
    )
    parser.add_argument(
        "--shots", type=int, default=8192, help="Finite-shot samples per benchmark control."
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    metadata = generate_all(args.root, seed=args.seed, shots=args.shots)
    print(
        "Generated all simulation artifacts in "
        f"{Path(args.root).resolve()} (seed={metadata['seed']}, shots={metadata['shots']})."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

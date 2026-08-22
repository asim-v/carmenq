"""Run the complete local reproducibility pipeline.

This command executes the automated tests first and regenerates every numerical
table, reduced-state archive, and figure only after the tests pass.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def run(command: list[str]) -> None:
    print("+", " ".join(command), flush=True)
    subprocess.run(command, cwd=ROOT, check=True)


def main() -> int:
    run([sys.executable, "-m", "pytest", "-q"])
    run(
        [
            sys.executable,
            "scripts/regenerate.py",
            "--seed",
            "20260812",
            "--shots",
            "8192",
        ]
    )
    run(
        [
            sys.executable,
            "scripts/classify_order_sensitive_checks.py",
            "--output",
            "data/order_sensitive_classes.json",
        ]
    )
    run([sys.executable, "scripts/verify_interleaved_counterexample.py"])
    run([sys.executable, "scripts/generate_order_gap_figure.py"])
    print("Reproducibility pipeline completed successfully.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

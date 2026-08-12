"""Build the manuscript and copy it to the stable release filename."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "manuscript" / "main.tex"
OUTPUT_DIR = ROOT / "output" / "pdf"
TECTONIC_OUTPUT = OUTPUT_DIR / "main.pdf"
STABLE_OUTPUT = OUTPUT_DIR / "causal_audit_return_memory.pdf"


def find_tectonic(explicit: str | None) -> str:
    candidates = [explicit, os.environ.get("TECTONIC"), shutil.which("tectonic")]
    for candidate in candidates:
        if candidate and Path(candidate).exists():
            return str(Path(candidate).resolve())
    raise SystemExit(
        "Tectonic was not found. Install it, add it to PATH, set TECTONIC, "
        "or pass --tectonic PATH."
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tectonic", help="Path to the Tectonic executable")
    parser.add_argument(
        "--keep-intermediates",
        action="store_true",
        help="Retain Tectonic auxiliary files for diagnostics",
    )
    args = parser.parse_args()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    command = [
        find_tectonic(args.tectonic),
        str(SOURCE),
        "--outdir",
        str(OUTPUT_DIR),
        "--keep-logs",
    ]
    if args.keep_intermediates:
        command.append("--keep-intermediates")

    print("+", " ".join(command), flush=True)
    subprocess.run(command, cwd=ROOT, check=True)
    shutil.copy2(TECTONIC_OUTPUT, STABLE_OUTPUT)
    print(f"Built {STABLE_OUTPUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

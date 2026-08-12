"""Source-tree wrapper for the installed ``carmenq`` command."""

from __future__ import annotations

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "src"
if str(SOURCE) not in sys.path:
    sys.path.insert(0, str(SOURCE))

from carmenq.audit_return_cli import main  # noqa: E402


if __name__ == "__main__":
    raise SystemExit(main())

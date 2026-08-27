"""Unit tests for the deterministic global-frontier release bundle."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
from zipfile import ZipFile

import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "build_global_frontier_bundle.py"
SPEC = importlib.util.spec_from_file_location("global_frontier_bundle", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
bundle = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(bundle)


def test_archive_metadata_is_reproducible() -> None:
    info = bundle.archive_info("proof/certificate.json")
    assert info.date_time == bundle.FIXED_ZIP_TIME
    assert info.filename == "proof/certificate.json"
    assert info.external_attr == 0o100644 << 16


def test_require_global_manifest_rejects_an_incomplete_file(tmp_path: Path) -> None:
    path = tmp_path / "global.json"
    path.write_text(json.dumps({"schema": bundle.GLOBAL_SCHEMA, "complete": False}))
    with pytest.raises(RuntimeError, match="absent or incomplete"):
        bundle.require_global_manifest(path)


def test_small_bundle_round_trip(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    root = tmp_path / "repo"
    (root / "data").mkdir(parents=True)
    manifest = {
        "schema": bundle.GLOBAL_SCHEMA,
        "complete": True,
        "reported_outward_upper_fraction": [7667, 10000],
        "explicit_physical_lower_fraction": [957373519, 1250000000],
    }
    (root / "data/global_frontier_l060_exact_assembly.json").write_text(
        json.dumps(manifest), encoding="utf-8"
    )
    (root / "payload.txt").write_text("exact certificate\n", encoding="utf-8")
    (root / "pyproject.toml").write_text(
        '[project]\nversion = "2.3.0"\n', encoding="utf-8"
    )
    monkeypatch.setattr(bundle, "ROOT", root)
    output = tmp_path / "proof.zip"
    result = bundle.build_bundle(output, ("payload.txt",))
    assert result["payload_file_count"] == 1
    assert bundle.verify_bundle(output)["valid"] is True
    with ZipFile(output) as archive:
        assert archive.read("payload.txt") == (root / "payload.txt").read_bytes()
    with pytest.raises(RuntimeError, match="not unique"):
        bundle.build_bundle(
            tmp_path / "duplicate.zip", ("payload.txt", "payload.txt")
        )

"""Build and verify the deterministic global-frontier certificate bundle."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
from typing import Any, Iterable
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo


ROOT = Path(__file__).resolve().parents[1]
FIXED_ZIP_TIME = (2026, 8, 27, 0, 0, 0)
BUNDLE_SCHEMA = "carmenq.global-frontier-release-bundle.v1"
GLOBAL_SCHEMA = "carmenq.global-frontier-l060-exact-assembly.v1"
PAYLOAD_PATHS = (
    "data/four_effect_rational_lower_l060.json",
    "data/global_frontier_l060_exact_assembly.json",
    "scratch/d2_frontier/projective_tangent_global_l055_L07573_summary.json",
    "scratch/d2_frontier/rank_rank_tangent_full_l055_L07573.json",
    "scratch/d2_frontier/remaining_projective_tangent_full_l055_L07573.json",
    "scratch/d2_frontier/projective_tangent_global_l060_L0766_summary.json",
    "scratch/d2_frontier/rank_rank_tangent_full_l060_L0766.json",
    "scratch/d2_frontier/remaining_projective_tangent_full_l060_L076591.json",
    "scratch/d2_frontier/low_weight_socp_exact_dual_l060.json",
    "scratch/d2_frontier/low_weight_socp_exact_dual_l060_verified.json",
    "scratch/d2_frontier/ternary_probability_cone_global_0777.json",
    *(
        f"scratch/d2_frontier/ternary_socp_exact_dual_full_l060_shard{i:02d}of08.json"
        for i in range(8)
    ),
    "scratch/d2_frontier/ternary_transferred_exact_dual_full_l060_verified.json",
    "scratch/d2_frontier/four_active_common_bias_fallback_exact_cover_l060.compact.json.gz",
    "scratch/d2_frontier/four_active_common_bias_fallback_exact_cover_l060_verified.json",
)


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def project_version() -> str:
    text = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    match = re.search(r'^version\s*=\s*"([^"]+)"', text, re.MULTILINE)
    if match is None:
        raise RuntimeError("project version is absent from pyproject.toml")
    return match.group(1)


def archive_info(name: str) -> ZipInfo:
    info = ZipInfo(name, date_time=FIXED_ZIP_TIME)
    info.compress_type = ZIP_DEFLATED
    info.create_system = 3
    info.external_attr = 0o100644 << 16
    return info


def bundle_readme(version: str) -> bytes:
    commands = """\
Replay after extracting this archive into a checkout of tag v{version}:

  python scratch/d2_frontier/validate_projective_tangent_global.py --weight .55 --level .7573 --rank-rank scratch/d2_frontier/rank_rank_tangent_full_l055_L07573.json --remaining scratch/d2_frontier/remaining_projective_tangent_full_l055_L07573.json --output projective_tangent_global_l055_L07573_summary.json
  python scratch/d2_frontier/validate_projective_tangent_global.py --weight .6 --level .766 --rank-rank scratch/d2_frontier/rank_rank_tangent_full_l060_L0766.json --remaining scratch/d2_frontier/remaining_projective_tangent_full_l060_L076591.json --output projective_tangent_global_l060_L0766_summary.json
  python scratch/d2_frontier/verify_low_weight_socp_exact_dual.py --certificate scratch/d2_frontier/low_weight_socp_exact_dual_l060.json --output low_weight_socp_exact_dual_l060_verified.json
  python scratch/d2_frontier/verify_ternary_transferred_exact_dual.py scratch/d2_frontier/ternary_socp_exact_dual_full_l060_shard00of08.json scratch/d2_frontier/ternary_socp_exact_dual_full_l060_shard01of08.json scratch/d2_frontier/ternary_socp_exact_dual_full_l060_shard02of08.json scratch/d2_frontier/ternary_socp_exact_dual_full_l060_shard03of08.json scratch/d2_frontier/ternary_socp_exact_dual_full_l060_shard04of08.json scratch/d2_frontier/ternary_socp_exact_dual_full_l060_shard05of08.json scratch/d2_frontier/ternary_socp_exact_dual_full_l060_shard06of08.json scratch/d2_frontier/ternary_socp_exact_dual_full_l060_shard07of08.json --workers 4 --output scratch/d2_frontier/ternary_transferred_exact_dual_full_l060_verified.json
  python scratch/d2_frontier/verify_four_active_mccormick_exact_cover.py scratch/d2_frontier/four_active_common_bias_fallback_exact_cover_l060.compact.json.gz --target .76670 --output scratch/d2_frontier/four_active_common_bias_fallback_exact_cover_l060_verified.json
  python scratch/d2_frontier/verify_global_frontier_l060.py --output ../../data/global_frontier_l060_exact_assembly.json

Every replay above is solver-free. Optimizers were used only to discover candidate
dual vectors. The projective kernel uses outward-expanded binary64 intervals;
the conic replays repair cone membership and evaluate residual corrections in
exact rational arithmetic. See the tagged manuscript for the analytic reductions
and the precise trust boundary.
""".format(version=version)
    return commands.encode("utf-8")


def require_global_manifest(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("schema") != GLOBAL_SCHEMA or payload.get("complete") is not True:
        raise RuntimeError("the global assembly manifest is absent or incomplete")
    if payload.get("reported_outward_upper_fraction") != [7667, 10000]:
        raise RuntimeError("the global upper endpoint changed")
    if payload.get("explicit_physical_lower_fraction") != [957373519, 1250000000]:
        raise RuntimeError("the global lower endpoint changed")
    return payload


def build_bundle(
    output: Path,
    payload_paths: Iterable[str] = PAYLOAD_PATHS,
) -> dict[str, Any]:
    paths = tuple(payload_paths)
    if len(set(paths)) != len(paths):
        raise RuntimeError("bundle payload paths are not unique")
    missing = [name for name in paths if not (ROOT / name).is_file()]
    if missing:
        raise RuntimeError("missing bundle payloads: " + ", ".join(missing))
    require_global_manifest(ROOT / "data/global_frontier_l060_exact_assembly.json")

    records = [
        {
            "path": name,
            "size": (ROOT / name).stat().st_size,
            "sha256": sha256_file(ROOT / name),
        }
        for name in paths
    ]
    version = project_version()
    manifest = {
        "schema": BUNDLE_SCHEMA,
        "version": version,
        "global_interval": {
            "lower": "957373519/1250000000",
            "upper": "7667/10000",
        },
        "payload_file_count": len(records),
        "files": records,
    }
    manifest_bytes = (json.dumps(manifest, indent=2) + "\n").encode("utf-8")
    readme_bytes = bundle_readme(version)
    checksummed = [
        ("BUNDLE_MANIFEST.json", sha256_bytes(manifest_bytes)),
        ("README.txt", sha256_bytes(readme_bytes)),
        *((record["path"], record["sha256"]) for record in records),
    ]
    sums_bytes = "".join(
        f"{digest}  {name}\n" for name, digest in checksummed
    ).encode("ascii")

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    with ZipFile(temporary, "w", compression=ZIP_DEFLATED, compresslevel=9) as archive:
        archive.writestr(archive_info("BUNDLE_MANIFEST.json"), manifest_bytes)
        archive.writestr(archive_info("README.txt"), readme_bytes)
        archive.writestr(archive_info("SHA256SUMS"), sums_bytes)
        for record in records:
            archive.writestr(
                archive_info(record["path"]),
                (ROOT / record["path"]).read_bytes(),
            )
    temporary.replace(output)
    verify_bundle(output)
    return {
        "path": str(output),
        "size": output.stat().st_size,
        "sha256": sha256_file(output),
        "payload_file_count": len(records),
    }


def verify_bundle(path: Path) -> dict[str, Any]:
    with ZipFile(path) as archive:
        manifest = json.loads(archive.read("BUNDLE_MANIFEST.json"))
        if manifest.get("schema") != BUNDLE_SCHEMA:
            raise RuntimeError("wrong bundle schema")
        records = manifest.get("files", [])
        record_paths = [record["path"] for record in records]
        if len(set(record_paths)) != len(record_paths):
            raise RuntimeError("bundle manifest contains duplicate paths")
        members = archive.namelist()
        if manifest.get("payload_file_count") != len(records):
            raise RuntimeError("bundle file count mismatch")
        if len(set(members)) != len(members):
            raise RuntimeError("ZIP contains duplicate members")
        names = set(members)
        expected = {
            "BUNDLE_MANIFEST.json",
            "README.txt",
            "SHA256SUMS",
            *record_paths,
        }
        if names != expected:
            raise RuntimeError("bundle contains missing or unexpected members")
        for record in records:
            data = archive.read(record["path"])
            if len(data) != record["size"]:
                raise RuntimeError(f"size mismatch for {record['path']}")
            if sha256_bytes(data) != record["sha256"]:
                raise RuntimeError(f"hash mismatch for {record['path']}")
        sums = archive.read("SHA256SUMS").decode("ascii").splitlines()
        parsed_sums: dict[str, str] = {}
        for line in sums:
            digest, name = line.split("  ", 1)
            if name in parsed_sums:
                raise RuntimeError(f"duplicate SHA256SUMS entry for {name}")
            parsed_sums[name] = digest
        expected_sums = {
            "BUNDLE_MANIFEST.json": sha256_bytes(
                archive.read("BUNDLE_MANIFEST.json")
            ),
            "README.txt": sha256_bytes(archive.read("README.txt")),
            **{record["path"]: record["sha256"] for record in records},
        }
        if parsed_sums != expected_sums:
            raise RuntimeError("SHA256SUMS entries do not match the manifest")
        for name, digest in parsed_sums.items():
            if sha256_bytes(archive.read(name)) != digest:
                raise RuntimeError(f"SHA256SUMS mismatch for {name}")
    return {
        "path": str(path),
        "version": manifest["version"],
        "payload_file_count": len(records),
        "valid": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT
        / f"output/releases/carmenq-global-frontier-v{project_version()}.zip",
    )
    parser.add_argument("--verify", type=Path)
    args = parser.parse_args()
    result = verify_bundle(args.verify) if args.verify else build_bundle(args.output)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

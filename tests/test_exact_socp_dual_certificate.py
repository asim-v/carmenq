from __future__ import annotations

import hashlib
import importlib
import json
import math
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest
import scipy.sparse as sp


ROOT = Path(__file__).resolve().parents[1]
FRONTIER = ROOT / "scratch" / "d2_frontier"
if str(FRONTIER) not in sys.path:
    sys.path.insert(0, str(FRONTIER))

CERTIFICATE = FRONTIER / "source_15818_exact_socp_certificate.json"


def _module():
    # The helper imports the frontier oracle, whose optional dependency is CVXPY.
    pytest.importorskip("cvxpy")
    return importlib.import_module("exact_socp_dual_certificate")


def _fraction(encoded: list[int]) -> Fraction:
    assert len(encoded) == 2
    return Fraction(int(encoded[0]), int(encoded[1]))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def test_binary64_conversion_and_fraction_encoding_are_exact() -> None:
    certificate = _module()

    assert certificate.exact_float(0.5) == Fraction(1, 2)
    assert certificate.exact_float(-0.0) == 0
    assert certificate.exact_float(0.1) == Fraction(*0.1.as_integer_ratio())
    assert certificate.encode_fraction(Fraction(-3, 7)) == [-3, 7]
    for nonfinite in (math.inf, -math.inf, math.nan):
        with pytest.raises(ValueError, match="finite"):
            certificate.exact_float(nonfinite)


def test_dyadic_soc_ray_is_strictly_inside_the_lorentz_cone() -> None:
    certificate = _module()
    ray = certificate.dyadic_soc_ray(
        np.asarray([0.3, 0.4, -0.2], dtype=float),
        cursor=7,
        label="toy-soc",
        denominator_bits=12,
    )

    assert ray is not None
    assert ray.cone == "soc"
    entries = dict(ray.entries)
    time = entries[7]
    spatial_square = sum(
        (value * value for row, value in entries.items() if row != 7),
        Fraction(0),
    )
    assert time > 0
    assert time * time > spatial_square
    assert all(
        value.denominator & (value.denominator - 1) == 0
        for value in entries.values()
    )
    assert certificate.dyadic_soc_ray(
        np.zeros(3), 0, "zero", denominator_bits=12
    ) is None


def test_exact_audit_accepts_a_toy_dual_certificate() -> None:
    certificate = _module()
    # min -x subject to x + s = 1, s >= 0, hence max x <= 1.
    matrix = sp.csc_matrix(np.asarray([[1.0]]))
    right = np.asarray([1.0])
    objective = np.asarray([-1.0])
    rays = [
        certificate.Ray(
            "nonnegative:0",
            ((0, Fraction(1)),),
            cone="nonnegative",
        )
    ]

    upper, audit = certificate.exact_certificate_audit(
        matrix,
        right,
        objective,
        rays,
        active=[0],
        coefficients=[Fraction(1)],
        target=Fraction(3, 2),
    )

    assert upper == 1
    assert audit == {
        "stationarity_exact": True,
        "conic_coefficients_nonnegative": True,
        "soc_rays_exact": True,
        "strict_target": True,
    }


def test_exact_audit_rejects_negative_conic_coefficients_and_bad_soc_rays() -> None:
    certificate = _module()
    nonnegative_ray = certificate.Ray(
        "nonnegative:0",
        ((0, Fraction(1)),),
        cone="nonnegative",
    )
    with pytest.raises(ValueError, match="negative conic coefficient"):
        certificate.exact_certificate_audit(
            sp.csc_matrix(np.asarray([[1.0]])),
            np.asarray([1.0]),
            np.asarray([1.0]),
            [nonnegative_ray],
            active=[0],
            coefficients=[Fraction(-1)],
            target=Fraction(2),
        )

    bad_soc_ray = certificate.Ray(
        "bad-soc",
        ((0, Fraction(1)), (1, Fraction(2))),
        cone="soc",
    )
    with pytest.raises(ValueError, match="invalid exact SOC ray"):
        certificate.exact_certificate_audit(
            sp.csc_matrix((2, 1)),
            np.zeros(2),
            np.zeros(1),
            [bad_soc_ray],
            active=[0],
            coefficients=[Fraction(0)],
            target=Fraction(1),
        )


def test_exact_audit_rejects_stationarity_and_nonstrict_target() -> None:
    certificate = _module()
    ray = certificate.Ray(
        "nonnegative:0",
        ((0, Fraction(1)),),
        cone="nonnegative",
    )
    matrix = sp.csc_matrix(np.asarray([[1.0]]))

    with pytest.raises(ValueError, match="exact stationarity failed"):
        certificate.exact_certificate_audit(
            matrix,
            np.asarray([1.0]),
            np.asarray([0.0]),
            [ray],
            active=[0],
            coefficients=[Fraction(1)],
            target=Fraction(2),
        )

    with pytest.raises(ValueError, match="does not beat target"):
        certificate.exact_certificate_audit(
            matrix,
            np.asarray([1.0]),
            np.asarray([-1.0]),
            [ray],
            active=[0],
            coefficients=[Fraction(1)],
            target=Fraction(1),
        )


def test_archived_source_15818_certificate_has_pinned_scope_and_integrity() -> None:
    payload = json.loads(CERTIFICATE.read_text(encoding="utf-8"))

    assert payload["schema"] == "carmenq.exact-socp-dual-certificate.v1"
    assert payload["scope"] == "canonical SOCP for stored source spectral cell"
    assert payload["source"] == {
        "path": "scratch\\d2_frontier\\ternary_reconstructed_depth4_g2_top_leaf_bbb_p1_s92_l055.json",
        "sha256": "8d314683c074d0aa59f5cac2677941f908d4d22350ed8187165f1f643a005884",
        "source_index": 15818,
        "source_cell": 608,
        "branches": ["bloch", "bloch", "scalar-negative", "bloch"],
        "caps": [8, 53, None, 8],
    }

    program = payload["canonical_program"]
    assert program["rows"] == 3173
    assert program["variables"] == 274
    assert program["nonzeros"] == 20962
    assert program["zero"] == 6
    assert program["nonnegative"] == 793
    assert len(program["soc"]) == 788
    assert program["zero"] + program["nonnegative"] + sum(program["soc"]) == program["rows"]
    assert program["coefficient_semantics"] == (
        "exact rational values of IEEE-754 binary64 canonical data"
    )

    exact = payload["exact_certificate"]
    target = _fraction(exact["target"])
    upper = _fraction(exact["upper"])
    margin = _fraction(exact["margin"])
    assert target == Fraction(379, 500)
    assert upper == Fraction(
        580812808889032592765723436703891690776891320554313167,
        766247770432944429179173513575154591809369561091801088,
    )
    assert margin == target - upper
    assert upper < target
    assert exact["upper_decimal"] == float(upper)
    assert exact["margin_decimal"] == float(margin)
    assert exact["audit"] == {
        "stationarity_exact": True,
        "conic_coefficients_nonnegative": True,
        "soc_rays_exact": True,
        "strict_target": True,
    }

    recovery = exact["recovery"]
    assert recovery["candidate_columns"] == 956
    assert recovery["active_columns"] == 956
    assert recovery["rank"] == 274
    assert recovery["free_parameters"] == 682
    assert recovery["basis_recovery"] == "FLINT RREF plus positive nullspace LP"

    selected = exact["selected_rays"]
    assert len(selected) == 523
    assert len({item["ray_index"] for item in selected}) == len(selected)
    assert len({item["label"] for item in selected}) == len(selected)
    for item in selected:
        coefficient = _fraction(item["coefficient"])
        if not item["free"]:
            assert coefficient >= 0
        entries = {
            int(entry["row"]): _fraction(entry["value"])
            for entry in item["entries"]
        }
        assert len(entries) == len(item["entries"])
        assert entries
        if item["cone"] == "nonnegative":
            assert all(value >= 0 for value in entries.values())
        elif item["cone"] == "soc":
            first_row = min(entries)
            time = entries[first_row]
            spatial_square = sum(
                (
                    value * value
                    for row, value in entries.items()
                    if row != first_row
                ),
                Fraction(0),
            )
            assert time >= 0
            assert time * time >= spatial_square
        else:
            assert item["cone"] == "free"
            assert item["free"]

    source_path = ROOT / Path(payload["source"]["path"].replace("\\", "/"))
    assert source_path.is_file()
    assert _sha256(source_path) == payload["source"]["sha256"]
    assert payload["epistemic_status"] == (
        "solver-independent exact dual certificate for the serialized canonical "
        "SOCP; upstream physical enclosure semantics not yet formalized"
    )

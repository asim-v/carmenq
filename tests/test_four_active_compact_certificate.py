from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys

import numpy as np
import pytest


RESEARCH = Path(__file__).resolve().parents[1] / "scratch" / "d2_frontier"
sys.path.insert(0, str(RESEARCH))

from compact_four_active_exact_cover import (  # noqa: E402
    COMPACT_SCHEMA,
    compact_artifact,
    read_json,
    write_json,
)
from four_active_mccormick_socp_exact_cover import (  # noqa: E402
    McCormickOracle,
    encode_dual,
    initial_box,
    physical_weight_vertices,
    weight_hull,
)
from verify_four_active_mccormick_exact_cover import verify_order  # noqa: E402


RETAINED = {
    "support_weight": [3, 5],
    "target": [38331, 50000],
    "maximum_weight_floor": [3533, 4000],
    "minimum_active_weight": [3, 10000],
    "projective_lines": [],
    "reserve_perturbations": [],
    "common_bias_coefficient_representatives": [],
    "common_bias_coefficient_orbit_size": 0,
    "weighted_reserve_geometry": "geometry",
    "projective_pair_geometry": "geometry",
    "prefix_order_reduction": "reduction",
    "relaxation": "relaxation",
    "initial_box": {"path": "", "coordinates": []},
    "boxes_split": 1,
    "leaf_count": 2,
    "closed_leaf_count": 1,
    "domain_empty_leaf_count": 1,
}


def source_artifact() -> dict:
    report = {
        "syndrome_permutation": [1, 2, 3],
        "dual_storage_dtype": "f32",
        "dual_zlib_base64": "dual",
        "derived_field": "discard me",
    }
    return {
        **RETAINED,
        "schema": "source.v1",
        "complete": True,
        "all_cells_closed": True,
        "open_boxes": [],
        "leaves": [
            {
                "kind": "closed",
                "box": {"path": "0", "coordinates": []},
                "order_certificates": [report],
                "derived_field": "discard me",
            },
            {
                "kind": "domain-empty",
                "box": {"path": "1", "coordinates": []},
                "derived_field": "discard me",
            },
        ],
    }


def test_compactor_keeps_only_irredundant_leaf_data(tmp_path: Path) -> None:
    compact = compact_artifact(source_artifact())
    assert compact["schema"] == COMPACT_SCHEMA
    assert compact["source_schema"] == "source.v1"
    assert set(compact["leaves"][0]) == {"kind", "box", "order_duals"}
    assert compact["leaves"][0]["order_duals"] == [
        {
            "syndrome_permutation": [1, 2, 3],
            "dual_storage_dtype": "f32",
            "dual_zlib_base64": "dual",
        }
    ]
    assert set(compact["leaves"][1]) == {"kind", "box"}

    destination = tmp_path / "certificate.json.gz"
    write_json(destination, compact)
    assert read_json(destination) == compact


def test_compactor_refuses_an_open_cover() -> None:
    source = source_artifact()
    source["complete"] = False
    with pytest.raises(RuntimeError, match="incomplete"):
        compact_artifact(source)


def test_minimal_compact_dual_is_reconstructed_solver_free() -> None:
    box = initial_box()
    hull = weight_hull(box)
    assert hull is not None
    oracle = McCormickOracle()
    oracle.assign(hull, (1, 2, 3), physical_weight_vertices(box))
    data = oracle.canonical_data()
    report = {
        "syndrome_permutation": [1, 2, 3],
        "dual_storage_dtype": "f32",
        "dual_zlib_base64": encode_dual(
            np.zeros(data["A"].shape[0], dtype=float), "f32"
        ),
    }
    upper = verify_order(
        oracle,
        box,
        hull,
        report,
        Fraction(10000),
        compact=True,
    )
    assert upper <= 10000
    with pytest.raises(RuntimeError, match="missing required fields"):
        verify_order(
            oracle,
            box,
            hull,
            report,
            Fraction(10000),
            compact=False,
        )

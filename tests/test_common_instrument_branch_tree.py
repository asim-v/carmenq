from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pytest


RESEARCH = Path(__file__).resolve().parents[1] / "scratch" / "d2_frontier"
sys.path.insert(0, str(RESEARCH))

from common_instrument_branch_tree import WitnessTemplate  # noqa: E402
from common_instrument_cells import StateCell  # noqa: E402


def test_witness_template_uses_cellwise_trace_radii() -> None:
    reference = np.zeros((4, 4))
    witness = np.zeros((4, 4, 4))
    lipschitz = np.asarray([0.2, 0.3, 0.4, 0.5])
    template = WitnessTemplate(reference, witness, lipschitz, 0.1, 0.28)
    cell = StateCell.root((0, 1, 2, 3))
    cut = template.cut_for_cell(cell)
    assert cut["restrict_to_balls"] is False
    assert np.asarray(cut["input_trace_radii"]).shape == (4,)
    assert template.uniform_radius_budget == pytest.approx(0.2)


def test_witness_template_json_is_complete() -> None:
    template = WitnessTemplate(
        np.zeros((4, 4)), np.zeros((4, 4, 4)), np.ones(4), 0.02, 0.2
    )
    payload = template.to_json()
    assert payload["reference_support"] == pytest.approx(0.02)
    assert payload["uniform_radius_budget"] == pytest.approx(0.05)
    assert np.asarray(payload["witness_bloch"]).shape == (4, 4, 4)
    restored = WitnessTemplate.from_json(payload)
    assert restored.reference_support == pytest.approx(template.reference_support)
    assert restored.source_gap == pytest.approx(template.source_gap)
    assert restored.witness_bloch == pytest.approx(template.witness_bloch)

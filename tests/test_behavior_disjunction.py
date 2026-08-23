from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pytest


RESEARCH = Path(__file__).resolve().parents[1] / "scratch" / "d2_frontier"
sys.path.insert(0, str(RESEARCH))

from behavior_disjunction_scip import (  # noqa: E402
    cut_key,
    normalise_cut,
    prior_box_from_payload,
)


def test_explicit_prior_box_is_self_contained() -> None:
    values = [0.2, 0.4, 0.15, 0.35, 0.1, 0.3, 0.05, 0.2]
    box = prior_box_from_payload(None, values)
    assert box.shape == (4, 2)
    assert box.ravel() == pytest.approx(values)


def test_prior_box_rejects_conflicting_or_inverted_inputs(tmp_path: Path) -> None:
    source = tmp_path / "box.json"
    source.write_text('{"initial_prior_box": [[0, 1], [0, 1], [0, 1], [0, 1]]}')
    with pytest.raises(ValueError, match="either"):
        prior_box_from_payload(source, [0.0, 1.0] * 4)
    with pytest.raises(ValueError, match="invalid"):
        prior_box_from_payload(None, [0.4, 0.2, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0])


def test_cut_normalisation_preserves_halfspace_and_key() -> None:
    cut = {"column": 6, "coefficients": [0.5, -2.0, 1.0, 0.25]}
    normalized = normalise_cut(cut)
    assert np.max(np.abs(normalized["coefficients"])) == pytest.approx(1.0)
    assert normalized["column"] == 6
    assert cut_key(cut) == cut_key(normalized)

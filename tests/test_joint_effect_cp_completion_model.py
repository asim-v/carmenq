from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pytest


RESEARCH = Path(__file__).resolve().parents[1] / "scratch" / "d2_frontier"
sys.path.insert(0, str(RESEARCH))

from joint_effect_helstrom_scip import (  # noqa: E402
    IDENTITY,
    build,
    canonical_three_effect_povm,
    extract_common_instrument,
    seed_from_common_instrument,
)


def test_exact_completion_accepts_and_reconstructs_literal_instrument(
    tmp_path: Path,
) -> None:
    effects = canonical_three_effect_povm(np.asarray([0.92, 0.64, 0.44]))
    states = np.repeat((IDENTITY / 8.0)[None, :, :], 4, axis=0)
    # Four equiprobable completely depolarizing subchannels. Their Choi
    # matrices sum to the Choi matrix of a trace-preserving channel.
    choi = np.repeat((np.eye(4) / 8.0)[None, :, :], 4, axis=0)
    checkpoint = tmp_path / "literal_instrument.npz"
    np.savez(checkpoint, states=states, choi=choi, effects=effects)

    model, variables = build(
        effects,
        weight=0.55,
        prefix_order=(0, 1, 2, 3),
        target=None,
        fix_rotation_gauge=True,
        linked_columns=None,
        require_cp_completion=True,
        ando_direction_count=16,
    )
    model.hideOutput()
    assert seed_from_common_instrument(
        model, variables, checkpoint, effects, weight=0.55
    )
    extracted = extract_common_instrument(model, variables, effects, weight=0.55)
    assert extracted is not None
    report, arrays = extracted
    assert report["score_from_reconstructed_instrument"] == pytest.approx(0.5875)
    assert report["audit_from_reconstructed_instrument"] == pytest.approx(0.25)
    assert report["return_from_reconstructed_instrument"] == pytest.approx(1.0)
    assert report["minimum_choi_eigenvalue"] >= -1e-12
    assert report["trace_preservation_residual"] <= 1e-12
    assert np.linalg.norm(arrays["choi"] - choi) <= 2e-12

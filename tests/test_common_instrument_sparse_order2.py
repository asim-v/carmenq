from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pytest

from carmenq.common_instrument import apply_choi, choi_from_kraus


RESEARCH = Path(__file__).resolve().parents[1] / "scratch" / "d2_frontier"
sys.path.insert(0, str(RESEARCH))

from common_instrument_sparse_order2 import (  # noqa: E402
    PAULIS,
    choi_output_coefficients,
    monomials,
    unique_monomials,
)


def test_sparse_basis_sizes_are_stable() -> None:
    assert len(monomials(range(16), 2)) == 153
    assert len(monomials(range(32), 2)) == 561
    assert unique_monomials([(), (2, 1), (1, 2)]) == [(), (1, 2)]


def test_pauli_choi_convention_matches_direct_application() -> None:
    state = np.asarray([[0.3, 0.07 + 0.02j], [0.07 - 0.02j, 0.2]])
    operator = np.asarray([[0.8, 0.1j], [0.2, 0.5]], dtype=complex)
    choi = choi_from_kraus((operator,))
    state_coefficients = np.asarray(
        [float(np.trace(state @ pauli).real) for pauli in PAULIS]
    )
    choi_coefficients = np.asarray(
        [
            [
                float(np.trace(choi @ np.kron(left, right)).real)
                for right in PAULIS
            ]
            for left in PAULIS
        ]
    )
    expected = apply_choi(choi, state)
    expected_coefficients = np.asarray(
        [float(np.trace(expected @ pauli).real) for pauli in PAULIS]
    )
    assert choi_output_coefficients(
        state_coefficients, choi_coefficients
    ) == pytest.approx(expected_coefficients)

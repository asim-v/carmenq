"""Reproduce the operator-basis obstruction at the Fourier outer optimum."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from audit_common_instrument_candidate import load_reported_family
from carmenq import reconstruct_common_instrument_from_basis
from fourier_behavior_upper import solve_behavior_outer
from fourier_branch_upper import PRIOR_BOX


def reproduce(
    summary_path: Path,
) -> tuple[dict[str, object], np.ndarray]:
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    cell = summary["fully_vectorial_behavior_cover"]["maximum_cell"]
    caps = (
        None,
        tuple(cell["plane_normal"] + [cell["plane_cosine"]]),
        tuple(cell["sphere_normal"] + [cell["sphere_cosine"]]),
    )
    outer = solve_behavior_outer(
        ("bloch", "bloch", "bloch"), PRIOR_BOX, caps
    )
    states, outputs = load_reported_family(outer)
    reconstruction = reconstruct_common_instrument_from_basis(states, outputs)
    eigensystems = [
        np.linalg.eigh(item) for item in reconstruction.signed_choi_numerators
    ]
    numerator_minima = [float(values[0]) for values, _ in eigensystems]
    witnesses = np.asarray([[vectors[:, 0]] for _, vectors in eigensystems])
    return {
        "source": summary_path.name,
        "outer_bound": outer["bound"],
        "outer_status": outer["status"],
        "input_determinant": reconstruction.input_determinant,
        "input_condition_number": reconstruction.input_condition_number,
        "minimum_choi_eigenvalues": reconstruction.minimum_choi_eigenvalues.tolist(),
        "minimum_signed_choi_numerator_eigenvalues": numerator_minima,
        "output_reconstruction_residual": reconstruction.output_residual,
        "trace_preservation_residual": reconstruction.trace_preservation_residual,
        "common_instrument_compatible": reconstruction.is_compatible(),
        "interpretation": (
            "the unique interpolating maps reproduce every conditioned output "
            "and are trace preserving in sum, but all four Choi matrices are non-PSD"
        ),
    }, witnesses


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--summary",
        type=Path,
        default=Path(__file__).with_name("fourier_interior_summary_l055.json"),
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument("--witness-npz", type=Path)
    args = parser.parse_args()
    payload, witnesses = reproduce(args.summary)
    rendered = json.dumps(payload, indent=2) + "\n"
    print(rendered, end="")
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    if args.witness_npz is not None:
        args.witness_npz.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(args.witness_npz, witnesses=witnesses)


if __name__ == "__main__":
    main()

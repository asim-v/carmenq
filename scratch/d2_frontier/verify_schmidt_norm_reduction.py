"""Verify the Schmidt-rank/S(2)-norm representation on a stored leaf."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import torch

from general_two_block_leaf import GeneralTwoBlockLeaf


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("--lambda", dest="weight", type=float, default=0.6)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    model = GeneralTwoBlockLeaf(0)
    model.load_state_dict(
        torch.load(args.checkpoint, map_location="cpu", weights_only=True)
    )
    with torch.no_grad():
        columns = model.columns().numpy()
        effects = model.povm().numpy()

    # Coefficients across (B_L,z) | (B_R,m,y).
    tensor = columns.reshape(4, 4, 2, 4, 4)
    paired = tensor.transpose(0, 3, 1, 2, 4).reshape(16, 32)
    singular = np.linalg.svd(paired, compute_uv=False)

    aligned = np.zeros_like(tensor)
    path_norms = np.zeros((4, 4), dtype=float)
    audit = 0.0
    for z in range(4):
        for y in range(4):
            # paired[(bL,z),(bR,m,y)] with the flattened index conventions
            # inherited from the tensor above.
            block = tensor[:, :, :, z, y]
            norm = float(np.linalg.norm(block))
            path_norms[z, y] = norm
            if norm > 0.0:
                aligned[:, :, :, z, y] = block / (4.0 * norm)
            # The differentiable leaf contracts a transposed physical effect.
            audit += float(
                np.einsum(
                    "lrm,mn,lrn->",
                    block.conj(),
                    effects[z ^ y].T,
                    block,
                ).real
            )

    aligned_paired = aligned.transpose(0, 3, 1, 2, 4).reshape(16, 32)
    overlap = np.vdot(aligned_paired, paired)
    return_from_overlap = float(abs(overlap) ** 2)
    return_from_blocks = float(path_norms.sum() ** 2 / 16.0)
    norm = float(np.linalg.norm(paired) ** 2)
    score = args.weight * audit + (1.0 - args.weight) * return_from_overlap

    payload = {
        "weight": args.weight,
        "normalisation": norm,
        "paired_singular_values": singular[:8].tolist(),
        "paired_schmidt_rank_at_1e-10": int(np.count_nonzero(singular > 1e-10)),
        "aligned_path_block_norms": [
            float(np.linalg.norm(aligned[:, :, :, z, y]))
            for z in range(4)
            for y in range(4)
        ],
        "audit_expectation": audit,
        "return_from_rank_one_overlap": return_from_overlap,
        "return_from_block_norms": return_from_blocks,
        "return_identity_residual": abs(return_from_overlap - return_from_blocks),
        "support_expectation": score,
    }
    rendered = json.dumps(payload, indent=2) + "\n"
    print(rendered)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

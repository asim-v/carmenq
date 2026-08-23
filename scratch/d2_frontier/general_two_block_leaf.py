"""General homogeneous bond-two leaf for the two full-rank block relaxation.

The ordered inputs ``(x1,x2)`` and ``(x3,x4)`` are each granted as one
four-dimensional slot.  Only a qubit crosses the block boundary.  Removing
the two internal one-bit cuts enlarges the physical four-slot class, so its
homogeneous optimum is an upper bound on the target frontier.  The script is
variational and does not certify the global maximum.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import torch


DTYPE = torch.complex128
RDTYPE = torch.float64
BLOCK = 4
TOTAL = 16


def random_complex(generator: torch.Generator, *shape: int) -> torch.Tensor:
    return torch.complex(
        torch.randn(*shape, generator=generator, dtype=RDTYPE),
        torch.randn(*shape, generator=generator, dtype=RDTYPE),
    ) / math.sqrt(2.0)


def inverse_square_root(matrix: torch.Tensor) -> torch.Tensor:
    values, vectors = torch.linalg.eigh(matrix)
    return (
        vectors
        * torch.clamp(values.real, min=1e-12).rsqrt().unsqueeze(-2)
    ) @ vectors.mH


class GeneralTwoBlockLeaf(torch.nn.Module):
    def __init__(self, seed: int) -> None:
        super().__init__()
        generator = torch.Generator().manual_seed(seed)
        self.left = torch.nn.Parameter(random_complex(generator, BLOCK, BLOCK, 2))
        self.right = torch.nn.Parameter(
            random_complex(generator, BLOCK, 2, BLOCK, 2)
        )
        self.raw_povm = torch.nn.Parameter(random_complex(generator, 4, 2, 2))

    @torch.no_grad()
    def load_four_slot_leaf(self, checkpoint: Path) -> None:
        """Embed a four-site bond-two leaf exactly into the block relaxation."""
        from general_single_leaf_bound import GeneralLeaf

        leaf = GeneralLeaf(0)
        leaf.load_state_dict(
            torch.load(checkpoint, map_location="cpu", weights_only=True)
        )
        columns = leaf.columns().detach().reshape((2,) * 9)
        paired = columns.permute(0, 1, 5, 6, 2, 3, 4, 7, 8).reshape(16, 32)
        left_vectors, singular, right_vectors = torch.linalg.svd(
            paired, full_matrices=False
        )
        if singular[2] > 1e-9:
            raise RuntimeError("four-slot checkpoint has middle rank above two")
        left = (left_vectors[:, :2] * singular[:2]).reshape(4, 4, 2)
        right = right_vectors[:2].reshape(2, 4, 2, 4).permute(1, 2, 3, 0)
        self.left.copy_(left)
        self.right.copy_(right)

        effects = leaf.povm().detach()
        values, vectors = torch.linalg.eigh(effects)
        roots = (
            vectors
            * torch.clamp(values.real, min=0.0).sqrt().unsqueeze(-2)
        ) @ vectors.mH
        self.raw_povm.copy_(roots)

    def columns(self) -> torch.Tensor:
        # K[bL,bR,m,z,y] = sum_k R[bR,m,y,k] L[bL,z,k].
        tensor = torch.einsum("rmyk,lzk->lrmzy", self.right, self.left)
        return tensor.reshape(32, TOTAL)

    def povm(self) -> torch.Tensor:
        gram = self.raw_povm.mH @ self.raw_povm
        normalizer = inverse_square_root(gram.sum(dim=0))
        return normalizer[None] @ gram @ normalizer[None]

    def quotient(self, weight: float) -> tuple[torch.Tensor, ...]:
        columns = self.columns()
        total = columns.abs().square().sum()
        terminal = columns.T.reshape(TOTAL, 16, 2)
        states = torch.einsum("xbi,xbj->xij", terminal.conj(), terminal)
        effects = self.povm()
        audit = torch.stack(
            [
                torch.trace(effects[z ^ y] @ states[BLOCK * z + y])
                for z in range(BLOCK)
                for y in range(BLOCK)
            ]
        ).real.sum()
        returned = torch.linalg.svdvals(columns).sum().square() / TOTAL
        score = (weight * audit + (1.0 - weight) * returned) / total
        return score, audit / total, returned / total


def optimise(
    seed: int,
    weight: float,
    steps: int,
    learning_rate: float,
    checkpoint: Path | None = None,
) -> tuple[dict[str, float | int], dict[str, torch.Tensor]]:
    model = GeneralTwoBlockLeaf(seed)
    if checkpoint is not None and seed == 0:
        model.load_four_slot_leaf(checkpoint)
    optimiser = torch.optim.Adam(model.parameters(), lr=learning_rate)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
        optimiser, steps, eta_min=learning_rate * 0.001
    )
    best = (-math.inf, 0.0, 0.0)
    best_state: dict[str, torch.Tensor] = {}
    for iteration in range(steps):
        optimiser.zero_grad(set_to_none=True)
        values = model.quotient(weight)
        (-values[0]).backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 20.0)
        optimiser.step()
        scheduler.step()
        point = tuple(float(value.detach()) for value in values)
        if point[0] > best[0]:
            best = point
            best_state = {
                name: value.detach().cpu().clone()
                for name, value in model.state_dict().items()
            }
        with torch.no_grad():
            scale = model.columns().abs().square().sum().pow(0.25)
            model.left.div_(torch.clamp(scale, min=1e-12))
            model.right.div_(torch.clamp(scale, min=1e-12))
        if iteration % 500 == 0 or iteration + 1 == steps:
            print(seed, iteration, *point, flush=True)
    return (
        {
            "seed": seed,
            "score": best[0],
            "audit": best[1],
            "return": best[2],
        },
        best_state,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lambda", dest="weight", type=float, default=0.6)
    parser.add_argument("--steps", type=int, default=5000)
    parser.add_argument("--restarts", type=int, default=12)
    parser.add_argument("--lr", type=float, default=0.012)
    parser.add_argument("--checkpoint", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rows = [
        optimise(seed, args.weight, args.steps, args.lr, args.checkpoint)
        for seed in range(args.restarts)
    ]
    rows.sort(key=lambda item: item[0]["score"], reverse=True)
    payload = [row[0] for row in rows]
    rendered = json.dumps(payload, indent=2) + "\n"
    print(rendered)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
        torch.save(rows[0][1], args.output.with_suffix(".pt"))


if __name__ == "__main__":
    main()

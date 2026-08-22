"""Relax the first two slots to an arbitrary four-state qubit ensemble.

The last two binary instruments and the common terminal qubit POVM are kept
exact.  This locates whether the middle-to-terminal causal geometry already
certifies the observed frontier or whether constraints from the prefix are
essential.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import torch

from cq_instrument_relaxation import CQInstrument
from general_single_leaf_bound import GeneralLeaf, random_complex
from pauli_complete_general_leaf import row_isometric_cores


def syndrome(word: int) -> int:
    bits = tuple((word >> (3 - index)) & 1 for index in range(4))
    return 2 * (bits[0] ^ bits[2]) + (bits[1] ^ bits[3])


class SuffixRelaxation(torch.nn.Module):
    def __init__(self, seed: int) -> None:
        super().__init__()
        generator = torch.Generator().manual_seed(seed)
        self.state_roots = torch.nn.Parameter(random_complex(generator, 4, 2, 2))
        self.raw3 = torch.nn.Parameter(random_complex(generator, 8, 2))
        self.raw4 = torch.nn.Parameter(random_complex(generator, 8, 2))
        self.raw_povm = torch.nn.Parameter(random_complex(generator, 4, 2, 2))

    @torch.no_grad()
    def load_leaf(self, checkpoint: Path) -> None:
        model = GeneralLeaf(0)
        model.load_state_dict(
            torch.load(checkpoint, map_location="cpu", weights_only=True)
        )
        tensors, _ = row_isometric_cores(model.columns().detach().numpy())

        states = [np.ones((1, 1), dtype=complex)]
        for tensor in tensors[:2]:
            next_states = []
            for state in states:
                for symbol in (0, 1):
                    block = tensor[:, :, symbol, :]
                    output = np.einsum(
                        "bni,ij,bmj->nm", block.conj(), state, block
                    )
                    next_states.append(output)
            states = next_states
        for index, state in enumerate(states):
            values, vectors = np.linalg.eigh(state)
            root = (
                vectors
                * np.sqrt(np.maximum(values.real, 0.0))[None, :]
            ) @ vectors.conj().T
            self.state_roots[index].copy_(torch.from_numpy(root))

        for raw, tensor in zip((self.raw3, self.raw4), tensors[2:], strict=True):
            matrix = tensor.transpose(2, 0, 1, 3).reshape(8, 2)
            raw.copy_(torch.from_numpy(matrix))

        effects = model.povm().detach()
        values, vectors = torch.linalg.eigh(effects)
        roots = (
            vectors * torch.clamp(values.real, min=0.0).sqrt().unsqueeze(-2)
        ) @ vectors.mH
        self.raw_povm.copy_(roots)

    def prefix_states(self) -> torch.Tensor:
        states = self.state_roots.mH @ self.state_roots
        return states / torch.diagonal(states, dim1=-2, dim2=-1).real.sum()

    def instrument(self, raw: torch.Tensor) -> torch.Tensor:
        matrix = raw @ CQInstrument.right_whitener(raw.mH @ raw)
        return matrix.reshape(2, 2, 2, 2)

    def povm(self) -> torch.Tensor:
        gram = self.raw_povm.mH @ self.raw_povm
        normalizer = CQInstrument.right_whitener(gram.sum(dim=0))
        return normalizer.mH[None] @ gram @ normalizer[None]

    def terminal_states(self) -> torch.Tensor:
        states = self.prefix_states()
        for raw in (self.raw3, self.raw4):
            tensor = self.instrument(raw)
            next_states = []
            for state in states:
                branches = torch.einsum(
                    "xbni,ij,xbmj->xnm", tensor.conj(), state, tensor
                )
                next_states.extend(branches.unbind(0))
            states = torch.stack(next_states)
        return states

    def scores(self, weight: float) -> tuple[torch.Tensor, ...]:
        states = self.terminal_states()
        probabilities = torch.diagonal(states, dim1=-2, dim2=-1).real.sum(dim=-1)
        effects = self.povm()
        audit = torch.stack(
            [torch.trace(effects[syndrome(word)] @ states[word]) for word in range(16)]
        ).real.sum()
        returned = torch.sqrt(torch.clamp(probabilities, min=1e-18)).sum().square() / 16
        return weight * audit + (1 - weight) * returned, audit, returned


def optimise(args: argparse.Namespace, seed: int) -> dict[str, float | int]:
    model = SuffixRelaxation(seed)
    if args.checkpoint is not None and seed == 0:
        model.load_leaf(args.checkpoint)
    optimiser = torch.optim.Adam(model.parameters(), lr=args.lr)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
        optimiser, args.steps, eta_min=args.lr * 0.001
    )
    best = (-math.inf, 0.0, 0.0)
    for iteration in range(args.steps):
        optimiser.zero_grad(set_to_none=True)
        values = model.scores(args.weight)
        (-values[0]).backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 20.0)
        optimiser.step()
        scheduler.step()
        point = tuple(float(value.detach()) for value in values)
        if point[0] > best[0]:
            best = point
        if iteration % 400 == 0 or iteration + 1 == args.steps:
            print(seed, iteration, *point, flush=True)
    return {"seed": seed, "score": best[0], "audit": best[1], "return_upper": best[2]}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lambda", dest="weight", type=float, default=0.5)
    parser.add_argument("--steps", type=int, default=3000)
    parser.add_argument("--restarts", type=int, default=6)
    parser.add_argument("--lr", type=float, default=0.008)
    parser.add_argument("--checkpoint", type=Path)
    args = parser.parse_args()
    results = [optimise(args, seed) for seed in range(args.restarts)]
    results.sort(key=lambda item: item["score"], reverse=True)
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()

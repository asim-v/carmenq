"""Homogeneous TT-rank-two upper-relaxation falsifier for one general leaf.

Unlike ``single_leaf_bound.py``, emitted carriers are not assumed QND.  Local
instrument completeness and sibling leaves are deliberately omitted.  A
certified maximum of this quotient would upper-bound every streamed strategy;
the variational maximum is only a diagnostic ceiling.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import torch


DTYPE = torch.complex128
RDTYPE = torch.float64
D = 16
BITS = tuple(tuple((word >> (3 - i)) & 1 for i in range(4)) for word in range(D))
SYNDROME = tuple(2 * (bits[0] ^ bits[2]) + (bits[1] ^ bits[3]) for bits in BITS)


def random_complex(generator: torch.Generator, *shape: int) -> torch.Tensor:
    return torch.complex(
        torch.randn(*shape, generator=generator, dtype=RDTYPE),
        torch.randn(*shape, generator=generator, dtype=RDTYPE),
    ) / math.sqrt(2)


def invsqrt(matrix: torch.Tensor) -> torch.Tensor:
    values, vectors = torch.linalg.eigh(matrix)
    return (vectors * torch.clamp(values.real, min=1e-12).rsqrt().unsqueeze(-2)) @ vectors.mH


def transcript_digits(index: int, arity: int) -> tuple[int, ...]:
    digits = []
    for depth in range(4):
        power = arity ** (3 - depth)
        digits.append(index // power)
        index %= power
    return tuple(digits)


def node_offset(arity: int, depth: int) -> int:
    return (arity**depth - 1) // (arity - 1)


class GeneralLeaf(torch.nn.Module):
    def __init__(self, seed: int) -> None:
        super().__init__()
        generator = torch.Generator().manual_seed(seed)
        self.core1 = torch.nn.Parameter(random_complex(generator, 2, 2, 2))
        self.core2 = torch.nn.Parameter(random_complex(generator, 2, 2, 2, 2))
        self.core3 = torch.nn.Parameter(random_complex(generator, 2, 2, 2, 2))
        self.core4 = torch.nn.Parameter(random_complex(generator, 2, 2, 2, 2))
        self.raw_povm = torch.nn.Parameter(random_complex(generator, 4, 2, 2))

    @torch.no_grad()
    def load_full_arity_leaf(
        self, checkpoint: Path, arity: int, transcript: int
    ) -> None:
        from n4_full_arity_search import FullArityStrategy

        model = FullArityStrategy(arity, seed=0)
        model.load_state_dict(
            torch.load(checkpoint, map_location="cpu", weights_only=True)
        )
        tree = model.tree()
        outcomes = transcript_digits(transcript, arity)
        prefix = 0
        selected = []
        for depth, outcome in enumerate(outcomes):
            node = node_offset(arity, depth) + prefix
            selected.append(tree[node, outcome].detach())
            prefix = arity * prefix + outcome
        # Local tensors are indexed by (B, M_new, A, M_old).
        self.core1.copy_(selected[0][:, :, :, 0].permute(0, 2, 1))
        self.core2.copy_(selected[1].permute(0, 2, 3, 1))
        self.core3.copy_(selected[2].permute(0, 2, 3, 1))
        self.core4.copy_(selected[3].permute(0, 2, 3, 1))
        effects = model.povms()[transcript].detach()
        values, vectors = torch.linalg.eigh(effects)
        roots = (vectors * torch.clamp(values.real, min=0.0).sqrt().unsqueeze(-2)) @ vectors.mH
        self.raw_povm.copy_(roots)

    def columns(self) -> torch.Tensor:
        output = torch.empty((32, D), dtype=DTYPE)
        for word, bits in enumerate(BITS):
            state = self.core1[:, bits[0], :]
            state = torch.einsum("al,blr->abr", state, self.core2[:, bits[1]])
            # The explicit contractions keep carrier indices separate.
            state = torch.einsum("abl,clr->abcr", state, self.core3[:, bits[2]])
            state = torch.einsum("abcl,dlm->abcdm", state, self.core4[:, bits[3]])
            output[:, word] = state.reshape(32)
        return output

    def povm(self) -> torch.Tensor:
        gram = self.raw_povm.mH @ self.raw_povm
        normalizer = invsqrt(gram.sum(dim=0))
        return normalizer[None] @ gram @ normalizer[None]

    def quotient(self, audit_weight: float) -> tuple[torch.Tensor, ...]:
        columns = self.columns()
        total = columns.square().abs().sum() / D
        terminal = columns.T.reshape(D, 16, 2)
        rho_m = torch.einsum("xbi,xbj->xij", terminal.conj(), terminal) / D
        effects = self.povm()
        audit = torch.stack(
            [torch.trace(effects[SYNDROME[word]] @ rho_m[word]) for word in range(D)]
        ).real.sum()
        singular = torch.linalg.svdvals(columns)
        returned = singular.sum().square() / D**2
        score = (audit_weight * audit + (1.0 - audit_weight) * returned) / total
        return score, audit / total, returned / total


def optimize(
    args: argparse.Namespace, seed: int
) -> tuple[dict[str, float | int], dict[str, torch.Tensor]]:
    model = GeneralLeaf(seed)
    if args.checkpoint and seed == 0:
        model.load_full_arity_leaf(args.checkpoint, args.arity, args.transcript)
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
        optimizer, args.steps, eta_min=args.lr * 0.001
    )
    best = (-math.inf, 0.0, 0.0)
    best_state: dict[str, torch.Tensor] = {}
    for iteration in range(args.steps):
        optimizer.zero_grad(set_to_none=True)
        score, audit, returned = model.quotient(args.audit_weight)
        (-score).backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 20.0)
        optimizer.step()
        scheduler.step()
        point = tuple(float(value.detach()) for value in (score, audit, returned))
        if point[0] > best[0]:
            best = point
            best_state = {
                name: value.detach().cpu().clone()
                for name, value in model.state_dict().items()
            }
        with torch.no_grad():
            norm = model.columns().square().abs().sum().pow(1 / 8)
            for core in (model.core1, model.core2, model.core3, model.core4):
                core.div_(norm)
        if iteration % 300 == 0 or iteration + 1 == args.steps:
            print(seed, iteration, *point, flush=True)
    return {
        "seed": seed,
        "score": best[0],
        "audit": best[1],
        "return": best[2],
    }, best_state


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lambda", dest="audit_weight", type=float, default=0.5)
    parser.add_argument("--steps", type=int, default=2500)
    parser.add_argument("--restarts", type=int, default=4)
    parser.add_argument("--lr", type=float, default=0.012)
    parser.add_argument("--checkpoint", type=Path)
    parser.add_argument("--arity", type=int, default=3)
    parser.add_argument("--transcript", type=int, default=46)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    optimized = [optimize(args, seed) for seed in range(args.restarts)]
    optimized.sort(key=lambda item: item[0]["score"], reverse=True)
    results = [item[0] for item in optimized]
    print(json.dumps(results, indent=2))
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
        torch.save(optimized[0][1], args.output.with_suffix(".pt"))


if __name__ == "__main__":
    main()

"""Classical-quantum instrument relaxation of the four-slot leaf frontier.

Pinching ``K^*K`` in the computational basis upper-bounds RETURN by the
classical fidelity of the induced word distribution.  Tracing every emitted
carrier turns a row-canonical Choi MPS into a sequence of binary qubit
instruments.  Optimising that enlarged model gives a diagnostic upper
relaxation; equality at a candidate would isolate the remaining converse to
a sequential cq statement.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import torch

from general_single_leaf_bound import GeneralLeaf, random_complex
from pauli_complete_general_leaf import row_isometric_cores


DTYPE = torch.complex128
RDTYPE = torch.float64


def syndrome(word: int) -> int:
    bits = tuple((word >> (3 - index)) & 1 for index in range(4))
    return 2 * (bits[0] ^ bits[2]) + (bits[1] ^ bits[3])


class CQInstrument(torch.nn.Module):
    def __init__(self, seed: int) -> None:
        super().__init__()
        generator = torch.Generator().manual_seed(seed)
        self.raw1 = torch.nn.Parameter(random_complex(generator, 8, 1))
        self.raw2 = torch.nn.Parameter(random_complex(generator, 8, 2))
        self.raw3 = torch.nn.Parameter(random_complex(generator, 8, 2))
        self.raw4 = torch.nn.Parameter(random_complex(generator, 8, 2))
        self.raw_povm = torch.nn.Parameter(random_complex(generator, 4, 2, 2))

    @torch.no_grad()
    def load_leaf(self, checkpoint: Path) -> None:
        model = GeneralLeaf(0)
        model.load_state_dict(
            torch.load(checkpoint, map_location="cpu", weights_only=True)
        )
        columns = model.columns().numpy()
        tensors, _ = row_isometric_cores(columns)
        for raw, tensor in zip(
            (self.raw1, self.raw2, self.raw3, self.raw4), tensors, strict=True
        ):
            # canonical tensor[b,new,x,old] -> matrix[(x,b,new),old]
            matrix = tensor.transpose(2, 0, 1, 3).reshape(8, tensor.shape[-1])
            raw.copy_(torch.from_numpy(matrix))
        effects = model.povm().detach()
        values, vectors = torch.linalg.eigh(effects)
        roots = (
            vectors
            * torch.clamp(values.real, min=0.0).sqrt().unsqueeze(-2)
        ) @ vectors.mH
        self.raw_povm.copy_(roots)

    @staticmethod
    def right_whitener(gram: torch.Tensor) -> torch.Tensor:
        """Return N with N^* gram N = I using a differentiable Cholesky gauge."""

        lower = torch.linalg.cholesky(gram)
        return torch.linalg.inv(lower.mH)

    @staticmethod
    def canonical(raw: torch.Tensor) -> torch.Tensor:
        matrix = raw @ CQInstrument.right_whitener(raw.mH @ raw)
        return matrix.reshape(2, 2, 2, raw.shape[1])  # x,b,new,old

    def instruments(self) -> tuple[torch.Tensor, ...]:
        return tuple(
            self.canonical(raw) for raw in (self.raw1, self.raw2, self.raw3, self.raw4)
        )

    def povm(self) -> torch.Tensor:
        gram = self.raw_povm.mH @ self.raw_povm
        normalizer = self.right_whitener(gram.sum(dim=0))
        return normalizer.mH[None] @ gram @ normalizer[None]

    def terminal_states(self) -> torch.Tensor:
        states = torch.ones((1, 1, 1), dtype=DTYPE)
        for tensor in self.instruments():
            next_states = []
            for state in states:
                # The surrounding Choi code stores terminal density matrices
                # as amplitude-conjugate times amplitude.  Keep that index
                # convention so a POVM copied from GeneralLeaf is unchanged.
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
        score = weight * audit + (1 - weight) * returned
        return score, audit, returned, probabilities.sum()


def optimise(args: argparse.Namespace, seed: int) -> tuple[dict[str, float | int], dict[str, torch.Tensor]]:
    model = CQInstrument(seed)
    if args.checkpoint is not None and seed == 0:
        model.load_leaf(args.checkpoint)
    optimiser = torch.optim.Adam(model.parameters(), lr=args.lr)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
        optimiser, args.steps, eta_min=args.lr * 0.001
    )
    best = (-math.inf, 0.0, 0.0, 0.0)
    best_state: dict[str, torch.Tensor] = {}
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
            best_state = {
                name: value.detach().cpu().clone()
                for name, value in model.state_dict().items()
            }
        if iteration % 500 == 0 or iteration + 1 == args.steps:
            print(seed, iteration, *point, flush=True)
    return (
        {
            "seed": seed,
            "score": best[0],
            "audit": best[1],
            "return_upper": best[2],
            "normalisation": best[3],
        },
        best_state,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lambda", dest="weight", type=float, default=0.5)
    parser.add_argument("--steps", type=int, default=4000)
    parser.add_argument("--restarts", type=int, default=8)
    parser.add_argument("--lr", type=float, default=0.008)
    parser.add_argument("--checkpoint", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    results = [optimise(args, seed) for seed in range(args.restarts)]
    results.sort(key=lambda item: item[0]["score"], reverse=True)
    payload = [item[0] for item in results]
    print(json.dumps(payload, indent=2))
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        torch.save(results[0][1], args.output.with_suffix(".pt"))


if __name__ == "__main__":
    main()

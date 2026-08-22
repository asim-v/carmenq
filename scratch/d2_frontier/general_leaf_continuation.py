"""Continuation sweep for the exact homogeneous Choi-MPS frontier."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import torch

from general_single_leaf_bound import GeneralLeaf
from n4_full_arity_independent import qubit_dual, syndrome


def normalize(model: GeneralLeaf) -> None:
    with torch.no_grad():
        scale = model.columns().abs().square().sum().pow(1 / 8)
        for core in (model.core1, model.core2, model.core3, model.core4):
            core.div_(torch.clamp(scale, min=1e-12))


def exact_audit(model: GeneralLeaf) -> float:
    with torch.no_grad():
        columns = model.columns().numpy()
    terminal = columns.T.reshape(16, 16, 2)
    states = np.zeros((4, 2, 2), dtype=complex)
    for word in range(16):
        states[syndrome(word)] += terminal[word].conj().T @ terminal[word]
    return qubit_dual(states)


def optimize(
    model: GeneralLeaf,
    weight: float,
    steps: int,
    learning_rate: float,
    lbfgs_steps: int,
) -> tuple[float, float, float]:
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
        optimizer, max(steps, 1), eta_min=learning_rate * 0.001
    )
    best = (-math.inf, 0.0, 0.0)
    best_state: dict[str, torch.Tensor] = {}
    for _ in range(steps):
        optimizer.zero_grad(set_to_none=True)
        point = model.quotient(weight)
        (-point[0]).backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 20.0)
        optimizer.step()
        scheduler.step()
        normalize(model)
        values = tuple(float(value.detach()) for value in point)
        if values[0] > best[0]:
            best = values
            best_state = {
                name: value.detach().clone() for name, value in model.state_dict().items()
            }
    if best_state:
        model.load_state_dict(best_state)

    if lbfgs_steps:
        optimizer_lbfgs = torch.optim.LBFGS(
            model.parameters(),
            lr=0.35,
            max_iter=lbfgs_steps,
            tolerance_grad=1e-12,
            tolerance_change=1e-15,
            history_size=80,
            line_search_fn="strong_wolfe",
        )

        def closure() -> torch.Tensor:
            optimizer_lbfgs.zero_grad(set_to_none=True)
            loss = -model.quotient(weight)[0]
            loss.backward()
            return loss

        optimizer_lbfgs.step(closure)
        normalize(model)

    with torch.no_grad():
        score, audit, returned = model.quotient(weight)
    audit_dual = exact_audit(model)
    returned_value = float(returned)
    exact_score = weight * audit_dual + (1 - weight) * returned_value
    return exact_score, audit_dual, returned_value


def parse_weights(text: str) -> list[float]:
    output = [float(value) for value in text.split(",") if value.strip()]
    if not output or any(not 0 <= value <= 1 for value in output):
        raise ValueError("weights must be a comma-separated list in [0,1]")
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("--weights", required=True)
    parser.add_argument("--steps", type=int, default=1600)
    parser.add_argument("--lr", type=float, default=0.008)
    parser.add_argument("--lbfgs-steps", type=int, default=160)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    model = GeneralLeaf(0)
    model.load_state_dict(torch.load(args.checkpoint, map_location="cpu", weights_only=True))
    args.output.mkdir(parents=True, exist_ok=True)
    rows = []
    for weight in parse_weights(args.weights):
        score, audit, returned = optimize(
            model, weight, args.steps, args.lr, args.lbfgs_steps
        )
        row = {
            "weight": weight,
            "score": score,
            "audit": audit,
            "return": returned,
            "no_record_score": 1 - weight / 2,
            "advantage": score - (1 - weight / 2),
        }
        rows.append(row)
        print(json.dumps(row), flush=True)
        torch.save(model.state_dict(), args.output / f"leaf_l{weight:.6f}.pt")
    (args.output / "frontier.json").write_text(
        json.dumps(rows, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()

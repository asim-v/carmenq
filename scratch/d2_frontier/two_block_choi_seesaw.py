"""Block-convex Choi seesaw for the complete two-block pinching relaxation.

The first block prepares four subnormalised qubit states ``rho[z]``.  The
second block is one four-outcome qubit instrument ``Phi[y]`` represented by
positive 4-by-4 Choi matrices whose sum is trace preserving.  Its sixteen
conditioned outputs are

    sigma[z,y] = Phi[y](rho[z]).

The terminal ensemble for syndrome ``s`` is the sum of ``sigma[z,y]`` over
``z xor y == s``.  RETURN is bounded by the computational-pinching/Hellinger
term.  Consequently this finite model retains the common Stinespring
compatibility missing from the earlier effect-only upper relaxations.

Each of the state, instrument, and terminal-POVM blocks is solved globally by
a conic programme while the other two blocks are held fixed.  The resulting
seesaw is a falsification and stationarity diagnostic, not a global
certificate for the jointly nonconvex model.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import cvxpy as cp
import numpy as np

from qubit_discrimination_geometry import discrimination_geometry


IDENTITY = np.eye(2, dtype=complex)
OUTCOMES = range(4)
PATHS = tuple((z, y) for z in OUTCOMES for y in OUTCOMES)


def inverse_square_root(matrix: np.ndarray) -> np.ndarray:
    values, vectors = np.linalg.eigh(matrix)
    return (vectors * np.maximum(values.real, 1e-14) ** -0.5) @ vectors.conj().T


def random_complex(generator: np.random.Generator, *shape: int) -> np.ndarray:
    return (
        generator.standard_normal(shape) + 1j * generator.standard_normal(shape)
    ) / math.sqrt(2.0)


def canonical_three_effect_povm(weights: np.ndarray) -> np.ndarray:
    """Return the planar rank-one qubit POVM with the requested traces."""
    w0, w1, w2 = (float(value) for value in weights)
    cosine = (w2 * w2 - w0 * w0 - w1 * w1) / (2.0 * w0 * w1)
    cosine = float(np.clip(cosine, -1.0, 1.0))
    sine = math.sqrt(max(0.0, 1.0 - cosine * cosine))
    vectors = np.asarray(
        [
            [1.0, 0.0, 0.0],
            [cosine, sine, 0.0],
            [-(w0 + w1 * cosine) / w2, -(w1 * sine) / w2, 0.0],
        ]
    )
    paulis = (
        np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex),
        np.array([[0.0, -1j], [1j, 0.0]], dtype=complex),
        np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex),
    )
    effects = np.zeros((4, 2, 2), dtype=complex)
    for index in range(3):
        effects[index] = 0.5 * weights[index] * (
            IDENTITY
            + sum(vectors[index, axis] * paulis[axis] for axis in range(3))
        )
    if np.linalg.norm(effects.sum(axis=0) - IDENTITY) > 2e-10:
        raise RuntimeError("canonical three-effect POVM failed completeness")
    return effects


def random_point(
    seed: int,
    povm_arity: int = 4,
    fixed_three_weights: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    generator = np.random.default_rng(seed)

    roots = random_complex(generator, 4, 2, 2)
    states = np.einsum("zai,zaj->zij", roots.conj(), roots)
    states /= np.trace(states, axis1=1, axis2=2).real.sum()

    kraus = random_complex(generator, 4, 2, 2)
    normalizer = inverse_square_root(
        np.einsum("yai,yaj->ij", kraus.conj(), kraus)
    )
    kraus = np.einsum("yab,bi->yai", kraus, normalizer)
    choi = np.stack([choi_from_kraus((operator,)) for operator in kraus])

    if fixed_three_weights is None:
        vectors = random_complex(generator, 4, 2)
        if povm_arity == 3:
            vectors[-1] = 0.0
        effects = np.einsum("si,sj->sij", vectors, vectors.conj())
        normalizer = inverse_square_root(effects.sum(axis=0))
        effects = np.einsum(
            "ai,sij,jb->sab", normalizer, effects, normalizer
        )
    else:
        effects = canonical_three_effect_povm(fixed_three_weights)
    return states, choi, effects


def choi_from_kraus(kraus: tuple[np.ndarray, ...]) -> np.ndarray:
    """Return the input-major Choi matrix of a qubit CP map."""
    result = np.zeros((4, 4), dtype=complex)
    for operator in kraus:
        vector = operator.T.reshape(4)
        result += np.outer(vector, vector.conj())
    return result


def apply_choi(choi: np.ndarray, state: np.ndarray) -> np.ndarray:
    blocks = choi.reshape(2, 2, 2, 2)
    return np.einsum("ij,iajb->ab", state, blocks)


def choi_effect(choi: np.ndarray) -> np.ndarray:
    blocks = choi.reshape(2, 2, 2, 2)
    # H is defined by Tr[H rho] = Tr[Phi(rho)].
    return np.einsum("iaja->ji", blocks)


def pulled_effect(choi: np.ndarray, effect: np.ndarray) -> np.ndarray:
    blocks = choi.reshape(2, 2, 2, 2)
    # H[j,i] = Tr(effect * Phi(|i><j|)).
    return np.einsum("ba,iajb->ji", effect, blocks)


def extract_checkpoint(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Convert a PyTorch two-block checkpoint to physical Choi variables."""
    import torch

    from general_two_block_leaf import GeneralTwoBlockLeaf

    model = GeneralTwoBlockLeaf(0)
    model.load_state_dict(torch.load(path, map_location="cpu", weights_only=True))
    with torch.no_grad():
        columns = model.columns().numpy()
        # The differentiable contraction stores transposed terminal effects.
        effects = model.povm().numpy().transpose(0, 2, 1)

    tensor = columns.reshape(4, 4, 2, 4, 4)
    paired = tensor.transpose(0, 3, 1, 2, 4).reshape(16, 32)
    left, singular, right_h = np.linalg.svd(paired, full_matrices=False)
    prefix = (left[:, :2] * singular[:2]).reshape(4, 4, 2)
    states = np.einsum("bzi,bzj->zij", prefix, prefix.conj())
    right = right_h[:2].reshape(2, 4, 2, 4).transpose(1, 2, 3, 0)

    choi = []
    for y in OUTCOMES:
        operators = tuple(right[r, :, y, :] for r in range(4))
        choi.append(choi_from_kraus(operators))
    return states, np.stack(choi), effects


def terminal_states(states: np.ndarray, choi: np.ndarray) -> np.ndarray:
    terminal = np.zeros((4, 2, 2), dtype=complex)
    for z, y in PATHS:
        terminal[z ^ y] += apply_choi(choi[y], states[z])
    return terminal


def path_probabilities(states: np.ndarray, choi: np.ndarray) -> np.ndarray:
    effects = np.stack([choi_effect(item) for item in choi])
    return np.einsum("yij,zji->zy", effects, states).real


def evaluate(
    states: np.ndarray,
    choi: np.ndarray,
    effects: np.ndarray,
    weight: float,
) -> dict[str, float]:
    terminal = terminal_states(states, choi)
    probabilities = path_probabilities(states, choi)
    audit = float(
        sum(
            np.trace(effects[z ^ y] @ apply_choi(choi[y], states[z])).real
            for z, y in PATHS
        )
    )
    returned = float(np.sqrt(np.maximum(probabilities, 0.0)).sum() ** 2 / 16.0)
    return {
        "score": weight * audit + (1.0 - weight) * returned,
        "audit": audit,
        "return_upper": returned,
        "normalisation": float(probabilities.sum()),
        "terminal_trace": float(np.trace(terminal, axis1=1, axis2=2).real.sum()),
    }


def hellinger_hypograph(
    probabilities: list[cp.Expression], constraints: list[cp.Constraint]
) -> cp.Expression:
    """Return a linear conic representation of ``(sum sqrt(p))**2 / 16``."""
    cross_terms: list[cp.Expression] = []
    for first in range(len(probabilities)):
        constraints.append(probabilities[first] >= 0)
        for second in range(first + 1, len(probabilities)):
            geometric = cp.Variable(nonneg=True)
            constraints.append(
                cp.SOC(
                    probabilities[first] + probabilities[second],
                    cp.hstack(
                        (
                            2.0 * geometric,
                            probabilities[first] - probabilities[second],
                        )
                    ),
                )
            )
            cross_terms.append(geometric)
    return (cp.sum(cp.hstack(probabilities)) + 2.0 * cp.sum(cp.hstack(cross_terms))) / 16.0


def instrument_output_expression(
    choi: cp.Expression, state: np.ndarray
) -> cp.Expression:
    entries = [
        [
            sum(
                state[i, j] * choi[2 * i + a, 2 * j + b]
                for i in range(2)
                for j in range(2)
            )
            for b in range(2)
        ]
        for a in range(2)
    ]
    output = cp.bmat(entries)
    return 0.5 * (output + output.H)


def state_output_expression(
    choi: np.ndarray, state: cp.Expression
) -> cp.Expression:
    blocks = choi.reshape(2, 2, 2, 2)
    entries = [
        [
            sum(
                state[i, j] * blocks[i, a, j, b]
                for i in range(2)
                for j in range(2)
            )
            for b in range(2)
        ]
        for a in range(2)
    ]
    output = cp.bmat(entries)
    return 0.5 * (output + output.H)


def impose_optimal_fixed_povm(
    terminal: list[cp.Expression],
    audit: cp.Expression,
    constraints: list[cp.Constraint],
) -> None:
    """Impose primal-dual equality for the frozen terminal POVM."""
    dual = cp.Variable((2, 2), hermitian=True)
    constraints.extend(dual - terminal[s] >> 0 for s in OUTCOMES)
    constraints.append(cp.real(cp.trace(dual)) == audit)


def solve_problem(problem: cp.Problem) -> None:
    problem.solve(
        solver="CLARABEL",
        tol_gap_abs=2e-9,
        tol_gap_rel=2e-9,
        tol_feas=2e-9,
        max_iter=500,
    )
    if problem.status not in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}:
        raise RuntimeError(f"conic block failed: {problem.status}")


def optimise_instrument(
    states: np.ndarray,
    effects: np.ndarray,
    weight: float,
    require_optimal_povm: bool = False,
) -> np.ndarray:
    variables = [cp.Variable((4, 4), hermitian=True) for _ in OUTCOMES]
    constraints: list[cp.Constraint] = [item >> 0 for item in variables]
    for i in range(2):
        for j in range(2):
            partial_trace = sum(
                variables[y][2 * i, 2 * j]
                + variables[y][2 * i + 1, 2 * j + 1]
                for y in OUTCOMES
            )
            constraints.append(partial_trace == (1.0 if i == j else 0.0))

    probabilities: list[cp.Expression] = []
    audit_terms: list[cp.Expression] = []
    for z, y in PATHS:
        probability_operator = np.kron(states[z].T, IDENTITY)
        audit_operator = np.kron(states[z].T, effects[z ^ y])
        probabilities.append(cp.real(cp.trace(probability_operator @ variables[y])))
        audit_terms.append(cp.real(cp.trace(audit_operator @ variables[y])))
    returned = hellinger_hypograph(probabilities, constraints)
    audit = cp.sum(cp.hstack(audit_terms))
    if require_optimal_povm:
        terminal = [
            sum(
                instrument_output_expression(variables[y], states[z])
                for z, y in PATHS
                if (z ^ y) == s
            )
            for s in OUTCOMES
        ]
        impose_optimal_fixed_povm(terminal, audit, constraints)
    objective = weight * audit + (1.0 - weight) * returned
    problem = cp.Problem(cp.Maximize(objective), constraints)
    solve_problem(problem)
    return np.stack([np.asarray(item.value) for item in variables])


def optimise_states(
    choi: np.ndarray,
    effects: np.ndarray,
    weight: float,
    require_optimal_povm: bool = False,
) -> np.ndarray:
    variables = [cp.Variable((2, 2), hermitian=True) for _ in OUTCOMES]
    constraints: list[cp.Constraint] = [item >> 0 for item in variables]
    constraints.append(
        cp.sum(cp.hstack([cp.real(cp.trace(item)) for item in variables])) == 1.0
    )
    coarse = [choi_effect(item) for item in choi]
    pulled = {
        (y, s): pulled_effect(choi[y], effects[s])
        for y in OUTCOMES
        for s in OUTCOMES
    }
    probabilities: list[cp.Expression] = []
    audit_terms: list[cp.Expression] = []
    for z, y in PATHS:
        probabilities.append(cp.real(cp.trace(coarse[y] @ variables[z])))
        audit_terms.append(
            cp.real(cp.trace(pulled[y, z ^ y] @ variables[z]))
        )
    returned = hellinger_hypograph(probabilities, constraints)
    audit = cp.sum(cp.hstack(audit_terms))
    if require_optimal_povm:
        terminal = [
            sum(
                state_output_expression(choi[y], variables[z])
                for z, y in PATHS
                if (z ^ y) == s
            )
            for s in OUTCOMES
        ]
        impose_optimal_fixed_povm(terminal, audit, constraints)
    objective = weight * audit + (1.0 - weight) * returned
    problem = cp.Problem(cp.Maximize(objective), constraints)
    solve_problem(problem)
    return np.stack([np.asarray(item.value) for item in variables])


def optimise_povm(
    states: np.ndarray,
    choi: np.ndarray,
    trace_floor: float,
) -> np.ndarray:
    terminal = terminal_states(states, choi)
    variables = [cp.Variable((2, 2), hermitian=True) for _ in OUTCOMES]
    constraints: list[cp.Constraint] = [item >> 0 for item in variables]
    constraints.append(sum(variables) == IDENTITY)
    if trace_floor > 0.0:
        constraints.extend(cp.real(cp.trace(item)) >= trace_floor for item in variables)
    objective = cp.sum(
        cp.hstack(
            [cp.real(cp.trace(variables[s] @ terminal[s])) for s in OUTCOMES]
        )
    )
    problem = cp.Problem(cp.Maximize(objective), constraints)
    solve_problem(problem)
    return np.stack([np.asarray(item.value) for item in variables])


def matrix_rank(matrix: np.ndarray, tolerance: float = 2e-7) -> int:
    return int(np.count_nonzero(np.linalg.eigvalsh(matrix).real > tolerance))


def optimise_seed(
    seed: int,
    weight: float,
    rounds: int,
    trace_floor: float,
    checkpoint: Path | None,
    freeze_povm: bool,
    random_povm_arity: int,
    fixed_three_weights: np.ndarray | None,
    require_optimal_fixed_povm: bool,
) -> tuple[dict[str, object], tuple[np.ndarray, np.ndarray, np.ndarray]]:
    if checkpoint is not None and seed == 0:
        states, choi, effects = extract_checkpoint(checkpoint)
    else:
        states, choi, effects = random_point(
            seed + 900_001,
            povm_arity=random_povm_arity,
            fixed_three_weights=fixed_three_weights,
        )

    # A seeded projective checkpoint is infeasible when a positive terminal
    # trace floor is requested.  Project it into the constrained POVM block
    # before recording an incumbent, otherwise the diagnostic would silently
    # report the unconstrained seed as its best point.
    if trace_floor > 0.0:
        if not freeze_povm:
            effects = optimise_povm(states, choi, trace_floor)

    history = []
    best = evaluate(states, choi, effects, weight)
    if require_optimal_fixed_povm:
        best["score"] = -math.inf
    best_arrays = (states.copy(), choi.copy(), effects.copy())
    for round_index in range(rounds):
        choi = optimise_instrument(
            states, effects, weight, require_optimal_fixed_povm
        )
        states = optimise_states(
            choi, effects, weight, require_optimal_fixed_povm
        )
        if not freeze_povm:
            effects = optimise_povm(states, choi, trace_floor)
        point = evaluate(states, choi, effects, weight)
        point["round"] = round_index + 1
        history.append(point)
        if point["score"] > best["score"]:
            best = point.copy()
            best_arrays = (states.copy(), choi.copy(), effects.copy())
        print(seed, point, flush=True)

    states, choi, effects = best_arrays
    geometry = discrimination_geometry(terminal_states(states, choi))
    row: dict[str, object] = {
        "seed": seed,
        "weight": weight,
        "povm_frozen": freeze_povm,
        "random_povm_arity": random_povm_arity,
        "fixed_three_povm_weights": (
            None if fixed_three_weights is None else fixed_three_weights.tolist()
        ),
        "optimal_fixed_povm_required": require_optimal_fixed_povm,
        **best,
        "state_ranks": [matrix_rank(item) for item in states],
        "choi_ranks": [matrix_rank(item) for item in choi],
        "effect_ranks": [matrix_rank(item) for item in effects],
        "effect_traces": [float(np.trace(item).real) for item in effects],
        "exact_audit_socp": geometry["optimal_guess_probability"],
        "active_indices": geometry["active_indices"],
        "constraint_slacks": geometry["constraint_slacks"],
        "kkt_weights": geometry["kkt_weights"],
        "history": history,
    }
    return row, best_arrays


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lambda", dest="weight", type=float, default=0.6)
    parser.add_argument("--rounds", type=int, default=8)
    parser.add_argument("--restarts", type=int, default=8)
    parser.add_argument("--seed-offset", type=int, default=0)
    parser.add_argument("--povm-trace-floor", type=float, default=0.0)
    parser.add_argument("--checkpoint", type=Path)
    parser.add_argument("--freeze-povm", action="store_true")
    parser.add_argument(
        "--random-povm-arity", type=int, choices=(3, 4), default=4
    )
    parser.add_argument("--fixed-three-povm-weights", type=float, nargs=3)
    parser.add_argument("--require-optimal-fixed-povm", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if not 0.0 <= args.povm_trace_floor <= 0.5:
        raise ValueError("povm trace floor must lie in [0, 0.5]")
    fixed_three_weights = (
        None
        if args.fixed_three_povm_weights is None
        else np.asarray(args.fixed_three_povm_weights, dtype=float)
    )
    if fixed_three_weights is not None:
        if (
            np.any(fixed_three_weights <= 0.0)
            or np.any(fixed_three_weights > 1.0)
            or abs(float(fixed_three_weights.sum()) - 2.0) > 1e-9
            or 2.0 * float(fixed_three_weights.max()) > 2.0 + 1e-9
        ):
            raise ValueError(
                "three-effect weights must be positive, at most one, and sum to two"
            )
        if args.checkpoint is not None:
            raise ValueError("fixed three-effect weights cannot seed a checkpoint")
    if args.require_optimal_fixed_povm and not args.freeze_povm:
        raise ValueError("optimal-fixed-POVM constraints require --freeze-povm")

    results = [
        optimise_seed(
            seed,
            args.weight,
            args.rounds,
            args.povm_trace_floor,
            args.checkpoint,
            args.freeze_povm,
            args.random_povm_arity,
            fixed_three_weights,
            args.require_optimal_fixed_povm,
        )
        for seed in range(args.seed_offset, args.seed_offset + args.restarts)
    ]
    results.sort(key=lambda item: float(item[0]["score"]), reverse=True)
    rows = [item[0] for item in results]
    rendered = json.dumps(rows, indent=2) + "\n"
    print(rendered)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
        states, choi, effects = results[0][1]
        np.savez_compressed(
            args.output.with_suffix(".npz"),
            states=states,
            choi=choi,
            effects=effects,
        )


if __name__ == "__main__":
    main()

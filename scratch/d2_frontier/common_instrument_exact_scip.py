"""Exact spatial model for one shared qubit instrument.

This is the nonconvex counterpart of ``common_instrument_sparse_order2.py``.
There are four subnormalised input states ``rho[z]`` and exactly four Choi
matrices ``J[y]``.  Every conditioned output is tied to those same variables,

    sigma[z,y] = Phi[J[y]](rho[z]),

and ``sum_y Phi[J[y]]`` is trace preserving.  Choi positivity is imposed by
an explicit complex Cholesky factorisation.  Consequently, every feasible
point is a literal common quantum instrument rather than a pseudo-moment
mixture of input-dependent instruments.

SCIP spatially relaxes the quadratic equalities and reports global primal and
dual bounds.  Those bounds remain solver-conditional until independently
outward-rounded; the formulation itself contains no moment relaxation.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from pyscipopt import Model, quicksum

from joint_effect_helstrom_scip import canonical_three_effect_povm


OUTCOMES = range(4)
ACTIVE = range(3)
PATHS = tuple((z, y) for z in OUTCOMES for y in OUTCOMES)
PAULIS = (
    np.eye(2, dtype=complex),
    np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex),
    np.array([[0.0, -1j], [1j, 0.0]], dtype=complex),
    np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex),
)
TRANSPOSE_SIGN = np.asarray([1.0, 1.0, -1.0, 1.0])
CHOI_BASIS = np.asarray(
    [[np.kron(left, right) for right in PAULIS] for left in PAULIS]
)


def add_lorentz(model: Model, scalar: object, vector: list[object]) -> None:
    """Impose ``scalar >= ||vector||_2`` as quadratic constraints."""

    model.addCons(scalar >= 0.0)
    model.addCons(scalar * scalar >= quicksum(item * item for item in vector))


def matrix_entry_from_pauli(
    coefficients: dict[tuple[int, int], object],
    row: int,
    column: int,
) -> tuple[object, object]:
    """Return real and imaginary parts of one Choi matrix entry."""

    real = quicksum(
        float(CHOI_BASIS[mu, nu, row, column].real)
        * coefficients[mu, nu]
        / 4.0
        for mu in OUTCOMES
        for nu in OUTCOMES
    )
    imaginary = quicksum(
        float(CHOI_BASIS[mu, nu, row, column].imag)
        * coefficients[mu, nu]
        / 4.0
        for mu in OUTCOMES
        for nu in OUTCOMES
    )
    return real, imaginary


def add_complex_cholesky(
    model: Model,
    coefficients: dict[tuple[int, int], object],
    label: str,
) -> dict[tuple[int, int, str], object]:
    """Impose positivity of a 4x4 Hermitian matrix by ``J = L L*``."""

    root_two = math.sqrt(2.0)
    factor: dict[tuple[int, int, str], object] = {}
    for row in OUTCOMES:
        for column in range(row + 1):
            factor[row, column, "real"] = model.addVar(
                lb=0.0 if row == column else -root_two,
                ub=root_two,
                name=f"{label}_L_{row}_{column}_re",
            )
            if row != column:
                factor[row, column, "imag"] = model.addVar(
                    lb=-root_two,
                    ub=root_two,
                    name=f"{label}_L_{row}_{column}_im",
                )

    def component(row: int, column: int, part: str) -> object:
        if part == "imag" and row == column:
            return 0.0
        return factor[row, column, part]

    for row in OUTCOMES:
        for column in range(row + 1):
            target_real, target_imaginary = matrix_entry_from_pauli(
                coefficients, row, column
            )
            product_real = quicksum(
                component(row, inner, "real")
                * component(column, inner, "real")
                + component(row, inner, "imag")
                * component(column, inner, "imag")
                for inner in range(column + 1)
            )
            model.addCons(target_real == product_real)
            if row != column:
                product_imaginary = quicksum(
                    component(row, inner, "imag")
                    * component(column, inner, "real")
                    - component(row, inner, "real")
                    * component(column, inner, "imag")
                    for inner in range(column + 1)
                )
                model.addCons(target_imaginary == product_imaginary)
    return factor


def pauli_coefficients(matrix: np.ndarray) -> np.ndarray:
    """Return ``Tr(matrix sigma_mu)`` for a 2x2 Hermitian matrix."""

    return np.asarray(
        [float(np.trace(matrix @ pauli).real) for pauli in PAULIS]
    )


def choi_pauli_coefficients(matrix: np.ndarray) -> np.ndarray:
    """Return ``Tr(J (sigma_mu tensor sigma_nu))``."""

    return np.asarray(
        [
            [
                float(np.trace(matrix @ CHOI_BASIS[mu, nu]).real)
                for nu in OUTCOMES
            ]
            for mu in OUTCOMES
        ]
    )


def output_coefficients(state: np.ndarray, choi: np.ndarray) -> np.ndarray:
    """Apply an input-major Choi matrix in Pauli coordinates."""

    return 0.5 * np.einsum(
        "m,m,mn->n", TRANSPOSE_SIGN, state, choi
    )


def gauge_rotation(state_values: np.ndarray) -> np.ndarray:
    """Return an SO(3) rotation putting states 0 and 1 in canonical gauge."""

    first = np.asarray(state_values[0, 1:], dtype=float)
    second = np.asarray(state_values[1, 1:], dtype=float)
    first_norm = float(np.linalg.norm(first))
    if first_norm <= 1e-12:
        return np.eye(3)
    x_axis = first / first_norm
    transverse = second - float(np.dot(second, x_axis)) * x_axis
    transverse_norm = float(np.linalg.norm(transverse))
    if transverse_norm <= 1e-12:
        trial = np.asarray([0.0, 0.0, 1.0])
        if abs(float(np.dot(trial, x_axis))) > 0.9:
            trial = np.asarray([0.0, 1.0, 0.0])
        transverse = trial - float(np.dot(trial, x_axis)) * x_axis
        transverse_norm = float(np.linalg.norm(transverse))
    y_axis = transverse / transverse_norm
    z_axis = np.cross(x_axis, y_axis)
    return np.vstack([x_axis, y_axis, z_axis])


def rotate_input_gauge(
    state_values: np.ndarray, choi_values: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Rotate input states and the Choi input leg without changing outputs."""

    rotation = gauge_rotation(state_values)
    rotated_states = np.array(state_values, copy=True)
    rotated_states[:, 1:] = state_values[:, 1:] @ rotation.T
    sign = np.diag(TRANSPOSE_SIGN[1:])
    choi_rotation = sign @ rotation @ sign
    rotated_choi = np.array(choi_values, copy=True)
    rotated_choi[:, 1:, :] = np.einsum(
        "ab,ybn->yan", choi_rotation, choi_values[:, 1:, :]
    )
    return rotated_states, rotated_choi


def build_model(
    weights: np.ndarray,
    support_weight: float,
    prefix_order: tuple[int, int, int, int],
    target: float | None,
    gauge_fix: bool,
) -> tuple[Model, dict[str, Any]]:
    """Build the exact common-instrument QCQP."""

    effects = canonical_three_effect_povm(weights)
    traces = np.trace(effects, axis1=1, axis2=2).real
    directions = np.zeros((4, 3), dtype=float)
    for terminal in ACTIVE:
        directions[terminal] = pauli_coefficients(
            effects[terminal] / traces[terminal]
        )[1:]

    model = Model("exact-common-instrument")
    variables: dict[str, Any] = {}

    state: dict[tuple[int, int], object] = {}
    for rank, z in enumerate(prefix_order):
        trace_upper = 1.0 / (rank + 1.0) if rank else 1.0
        state[z, 0] = model.addVar(
            lb=0.25 if rank == 0 else 0.0,
            ub=trace_upper,
            name=f"rho_{z}_0",
        )
        for mu in range(1, 4):
            state[z, mu] = model.addVar(
                lb=-trace_upper, ub=trace_upper, name=f"rho_{z}_{mu}"
            )
        add_lorentz(
            model, state[z, 0], [state[z, mu] for mu in range(1, 4)]
        )
    model.addCons(quicksum(state[z, 0] for z in OUTCOMES) == 1.0)
    for rank in range(3):
        model.addCons(
            state[prefix_order[rank], 0]
            >= state[prefix_order[rank + 1], 0]
        )
    if gauge_fix:
        # Quotient the simultaneous input-unitary orbit.  This convention is
        # also valid on the zero and collinear strata.
        model.addCons(state[0, 2] == 0.0)
        model.addCons(state[0, 3] == 0.0)
        model.addCons(state[0, 1] >= 0.0)
        model.addCons(state[1, 3] == 0.0)
        model.addCons(state[1, 2] >= 0.0)

    choi: dict[tuple[int, int, int], object] = {}
    factors: dict[int, dict[tuple[int, int, str], object]] = {}
    for y in OUTCOMES:
        local: dict[tuple[int, int], object] = {}
        for mu in OUTCOMES:
            for nu in OUTCOMES:
                lower, upper = ((0.0, 2.0) if (mu, nu) == (0, 0) else (-2.0, 2.0))
                item = model.addVar(
                    lb=lower, ub=upper, name=f"J_{y}_{mu}_{nu}"
                )
                choi[y, mu, nu] = item
                local[mu, nu] = item
                if (mu, nu) != (0, 0):
                    model.addCons(item <= choi[y, 0, 0])
                    model.addCons(item >= -choi[y, 0, 0])
        factors[y] = add_complex_cholesky(model, local, f"J_{y}")
        # The partial input trace of each subchannel is a positive qubit
        # effect.  It follows from Choi positivity but materially tightens
        # spatial node relaxations.
        add_lorentz(
            model,
            choi[y, 0, 0],
            [choi[y, mu, 0] for mu in range(1, 4)],
        )

    # One and the same four-outcome instrument is trace preserving.  These
    # equations are shared globally, rather than imposed independently for
    # each input or terminal statistic.
    for mu in OUTCOMES:
        model.addCons(
            quicksum(choi[y, mu, 0] for y in OUTCOMES)
            == (2.0 if mu == 0 else 0.0)
        )

    output: dict[tuple[int, int, int], object] = {}
    probability: dict[tuple[int, int], object] = {}
    statistics: dict[tuple[int, int, int], object] = {}
    for z, y in PATHS:
        for nu in OUTCOMES:
            item = model.addVar(
                lb=0.0 if nu == 0 else -1.0,
                ub=1.0,
                name=f"sigma_{z}_{y}_{nu}",
            )
            output[z, y, nu] = item
            model.addCons(
                2.0 * item
                == quicksum(
                    float(TRANSPOSE_SIGN[mu])
                    * state[z, mu]
                    * choi[y, mu, nu]
                    for mu in OUTCOMES
                )
            )
        probability[z, y] = output[z, y, 0]
        add_lorentz(
            model,
            probability[z, y],
            [output[z, y, nu] for nu in range(1, 4)],
        )
        model.addCons(probability[z, y] <= state[z, 0])
        for terminal in ACTIVE:
            item = model.addVar(
                lb=0.0,
                ub=float(traces[terminal]),
                name=f"q_{z}_{y}_{terminal}",
            )
            statistics[z, y, terminal] = item
            model.addCons(
                2.0 * item
                == float(traces[terminal])
                * (
                    probability[z, y]
                    + quicksum(
                        float(directions[terminal, axis])
                        * output[z, y, axis + 1]
                        for axis in range(3)
                    )
                )
            )
            model.addCons(
                item <= float(traces[terminal]) * probability[z, y]
            )
        model.addCons(
            probability[z, y]
            == quicksum(statistics[z, y, terminal] for terminal in ACTIVE)
        )

    for z in OUTCOMES:
        model.addCons(
            quicksum(probability[z, y] for y in OUTCOMES) == state[z, 0]
        )
    model.addCons(quicksum(probability.values()) == 1.0)

    syndrome_prior = {
        syndrome: quicksum(
            probability[z, z ^ syndrome] for z in OUTCOMES
        )
        for syndrome in OUTCOMES
    }
    syndrome_vector = {
        syndrome: [
            quicksum(
                output[z, z ^ syndrome, axis + 1] for z in OUTCOMES
            )
            for axis in range(3)
        ]
        for syndrome in OUTCOMES
    }
    for syndrome in OUTCOMES:
        add_lorentz(model, syndrome_prior[syndrome], syndrome_vector[syndrome])

    audit = model.addVar(lb=0.0, ub=1.0, name="audit")
    model.addCons(
        audit
        == quicksum(
            statistics[z, y, z ^ y]
            for z, y in PATHS
            if (z ^ y) in ACTIVE
        )
    )
    cap = [float(traces.max()), float(traces.max()), 2.0 - 2.0 * float(traces.max()), 0.0]
    model.addCons(
        audit
        <= quicksum(
            cap[rank] * state[prefix_order[rank], 0]
            for rank in OUTCOMES
        )
    )
    model.addCons(
        audit
        <= quicksum(
            float(traces[terminal]) * syndrome_prior[terminal]
            for terminal in ACTIVE
        )
    )

    # The fixed terminal POVM is required to be Helstrom optimal in this
    # sector.  A dual qubit operator with trace equal to the achieved success
    # enforces the matching upper certificate.
    dual_vector = [
        model.addVar(lb=-1.0, ub=1.0, name=f"dual_{axis}")
        for axis in range(3)
    ]
    add_lorentz(model, audit, dual_vector)
    for syndrome in OUTCOMES:
        add_lorentz(
            model,
            audit - syndrome_prior[syndrome],
            [
                dual_vector[axis] - syndrome_vector[syndrome][axis]
                for axis in range(3)
            ],
        )

    flat = [probability[z, y] for z, y in PATHS]
    hellinger: dict[tuple[int, int], object] = {}
    for first in range(16):
        for second in range(first + 1, 16):
            item = model.addVar(
                lb=0.0, ub=0.5, name=f"h_{first}_{second}"
            )
            model.addCons(item * item <= flat[first] * flat[second])
            hellinger[first, second] = item
    returned = (1.0 + 2.0 * quicksum(hellinger.values())) / 16.0
    score = model.addVar(lb=0.0, ub=1.0, name="score")
    model.addCons(
        score <= support_weight * audit + (1.0 - support_weight) * returned
    )
    if target is not None:
        model.addCons(score >= target)
    model.setObjective(score, "maximize")

    variables.update(
        {
            "state": state,
            "choi": choi,
            "factors": factors,
            "output": output,
            "probability": probability,
            "statistics": statistics,
            "syndrome_prior": syndrome_prior,
            "syndrome_vector": syndrome_vector,
            "dual_vector": dual_vector,
            "audit": audit,
            "hellinger": hellinger,
            "return": returned,
            "score": score,
            "effects": effects,
        }
    )
    return model, variables


def seed_from_checkpoint(
    model: Model,
    variables: dict[str, Any],
    checkpoint: Path,
    support_weight: float,
    gauge_fix: bool = False,
) -> bool:
    """Submit a physical state--Choi checkpoint as a complete SCIP start."""

    arrays = np.load(checkpoint)
    states = np.asarray(arrays["states"], dtype=complex)
    choi_matrices = np.asarray(arrays["choi"], dtype=complex)
    effects = np.asarray(variables["effects"], dtype=complex)
    if states.shape != (4, 2, 2) or choi_matrices.shape != (4, 4, 4):
        raise ValueError("checkpoint must contain states (4,2,2) and choi (4,4,4)")

    state_values = np.asarray([pauli_coefficients(item) for item in states])
    choi_values = np.asarray(
        [choi_pauli_coefficients(item) for item in choi_matrices]
    )
    if gauge_fix:
        state_values, choi_values = rotate_input_gauge(
            state_values, choi_values
        )
        states = np.asarray(
            [
                sum(state_values[z, mu] * PAULIS[mu] for mu in OUTCOMES)
                / 2.0
                for z in OUTCOMES
            ]
        )
        choi_matrices = np.asarray(
            [
                sum(
                    choi_values[y, mu, nu] * CHOI_BASIS[mu, nu]
                    for mu in OUTCOMES
                    for nu in OUTCOMES
                )
                / 4.0
                for y in OUTCOMES
            ]
        )
    output_values = np.asarray(
        [
            [output_coefficients(state_values[z], choi_values[y]) for y in OUTCOMES]
            for z in OUTCOMES
        ]
    )
    statistics_values = np.asarray(
        [
            [
                [
                    float(
                        0.5
                        * np.dot(
                            pauli_coefficients(effects[terminal]),
                            output_values[z, y],
                        )
                    )
                    for terminal in ACTIVE
                ]
                for y in OUTCOMES
            ]
            for z in OUTCOMES
        ]
    )
    probabilities = output_values[:, :, 0]
    audit_value = float(
        sum(
            statistics_values[z, y, z ^ y]
            for z, y in PATHS
            if (z ^ y) in ACTIVE
        )
    )
    root_sum = float(np.sqrt(np.maximum(probabilities, 0.0)).sum())
    return_value = root_sum * root_sum / 16.0

    # At an optimum using the fixed POVM, sum_s E_s tau_s is the Helstrom
    # dual.  Symmetrisation only removes checkpoint roundoff.
    terminal_states = np.asarray(
        [
            sum(
                _apply_choi(choi_matrices[z ^ syndrome], states[z])
                for z in OUTCOMES
            )
            for syndrome in OUTCOMES
        ]
    )
    dual_matrix = sum(
        effects[terminal] @ terminal_states[terminal]
        for terminal in ACTIVE
    )
    dual_matrix = 0.5 * (dual_matrix + dual_matrix.conj().T)
    dual_coefficients = pauli_coefficients(dual_matrix)

    solution = model.createSol()
    for z in OUTCOMES:
        for mu in OUTCOMES:
            model.setSolVal(
                solution, variables["state"][z, mu], float(state_values[z, mu])
            )
    for y in OUTCOMES:
        matrix = 0.5 * (choi_matrices[y] + choi_matrices[y].conj().T)
        factor = np.linalg.cholesky(matrix)
        for mu in OUTCOMES:
            for nu in OUTCOMES:
                model.setSolVal(
                    solution,
                    variables["choi"][y, mu, nu],
                    float(choi_values[y, mu, nu]),
                )
        for row in OUTCOMES:
            for column in range(row + 1):
                model.setSolVal(
                    solution,
                    variables["factors"][y][row, column, "real"],
                    float(factor[row, column].real),
                )
                if row != column:
                    model.setSolVal(
                        solution,
                        variables["factors"][y][row, column, "imag"],
                        float(factor[row, column].imag),
                    )
    for z, y in PATHS:
        for nu in OUTCOMES:
            model.setSolVal(
                solution,
                variables["output"][z, y, nu],
                float(output_values[z, y, nu]),
            )
        for terminal in ACTIVE:
            model.setSolVal(
                solution,
                variables["statistics"][z, y, terminal],
                float(statistics_values[z, y, terminal]),
            )
    for axis in range(3):
        model.setSolVal(
            solution,
            variables["dual_vector"][axis],
            float(dual_coefficients[axis + 1]),
        )
    model.setSolVal(solution, variables["audit"], audit_value)
    for first in range(16):
        z1, y1 = divmod(first, 4)
        for second in range(first + 1, 16):
            z2, y2 = divmod(second, 4)
            model.setSolVal(
                solution,
                variables["hellinger"][first, second],
                float(math.sqrt(probabilities[z1, y1] * probabilities[z2, y2])),
            )
    model.setSolVal(
        solution,
        variables["score"],
        support_weight * audit_value + (1.0 - support_weight) * return_value,
    )
    return bool(model.addSol(solution, free=True))


def _apply_choi(choi: np.ndarray, state: np.ndarray) -> np.ndarray:
    blocks = choi.reshape(2, 2, 2, 2)
    output = np.einsum("ij,iajb->ab", state, blocks)
    return 0.5 * (output + output.conj().T)


def extract_solution(model: Model, variables: dict[str, Any]) -> dict[str, Any] | None:
    solution = model.getBestSol()
    if solution is None:
        return None

    def get(variable: object) -> float:
        return float(model.getSolVal(solution, variable))

    state = np.asarray(
        [[get(variables["state"][z, mu]) for mu in OUTCOMES] for z in OUTCOMES]
    )
    choi = np.asarray(
        [
            [
                [get(variables["choi"][y, mu, nu]) for nu in OUTCOMES]
                for mu in OUTCOMES
            ]
            for y in OUTCOMES
        ]
    )
    output = np.asarray(
        [
            [
                [get(variables["output"][z, y, nu]) for nu in OUTCOMES]
                for y in OUTCOMES
            ]
            for z in OUTCOMES
        ]
    )
    probabilities = output[:, :, 0]
    audit = get(variables["audit"])
    returned = float(np.sqrt(np.maximum(probabilities, 0.0)).sum() ** 2 / 16.0)
    dual = np.asarray(
        [audit, *(get(item) for item in variables["dual_vector"])]
    )
    syndrome = np.asarray(
        [
            [
                sum(output[z, z ^ s, nu] for z in OUTCOMES)
                for nu in OUTCOMES
            ]
            for s in OUTCOMES
        ]
    )
    helstrom_slacks = [
        dual[0]
        - syndrome[s, 0]
        - float(np.linalg.norm(dual[1:] - syndrome[s, 1:]))
        for s in OUTCOMES
    ]
    return {
        "score": get(variables["score"]),
        "objective_from_reported": float(
            model.getObjVal()
        ),
        "audit": audit,
        "return": returned,
        "prefix_priors": state[:, 0].tolist(),
        "state_pauli_coefficients": state.tolist(),
        "choi_pauli_coefficients": choi.tolist(),
        "path_probabilities": probabilities.tolist(),
        "helstrom_dual_pauli_coefficients": dual.tolist(),
        "minimum_helstrom_lorentz_slack": float(min(helstrom_slacks)),
        "minimum_state_lorentz_slack": float(
            min(state[z, 0] - np.linalg.norm(state[z, 1:]) for z in OUTCOMES)
        ),
        "trace_preservation_residual": float(
            max(
                abs(choi[:, mu, 0].sum() - (2.0 if mu == 0 else 0.0))
                for mu in OUTCOMES
            )
        ),
    }


def save_solution_arrays(path: Path, payload: dict[str, Any], effects: np.ndarray) -> None:
    """Archive a feasible point in the package's state/Choi checkpoint format."""

    state = np.asarray(payload["state_pauli_coefficients"], dtype=float)
    choi = np.asarray(payload["choi_pauli_coefficients"], dtype=float)
    states = np.asarray(
        [
            sum(state[z, mu] * PAULIS[mu] for mu in OUTCOMES) / 2.0
            for z in OUTCOMES
        ]
    )
    choi_matrices = np.asarray(
        [
            sum(
                choi[y, mu, nu] * CHOI_BASIS[mu, nu]
                for mu in OUTCOMES
                for nu in OUTCOMES
            )
            / 4.0
            for y in OUTCOMES
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(path, states=states, choi=choi_matrices, effects=effects)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lambda", dest="weight", type=float, default=0.55)
    parser.add_argument(
        "--fixed-three-povm-weights", type=float, nargs=3, required=True
    )
    parser.add_argument("--prefix-order", type=int, nargs=4, required=True)
    parser.add_argument("--target", type=float)
    parser.add_argument("--seconds", type=float, default=60.0)
    parser.add_argument("--gap", type=float, default=1e-4)
    parser.add_argument("--gauge-fix", action="store_true")
    parser.add_argument("--checkpoint", type=Path)
    parser.add_argument("--solution-npz", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    order = tuple(int(item) for item in args.prefix_order)
    if sorted(order) != list(OUTCOMES):
        raise ValueError("prefix order must be a permutation of 0,1,2,3")
    weights = np.asarray(args.fixed_three_povm_weights, dtype=float)
    model, variables = build_model(
        weights, args.weight, order, args.target, args.gauge_fix
    )
    seed_accepted = None
    if args.checkpoint is not None:
        seed_accepted = seed_from_checkpoint(
            model, variables, args.checkpoint, args.weight, args.gauge_fix
        )
    model.setRealParam("limits/time", args.seconds)
    model.setRealParam("limits/gap", args.gap)
    model.setRealParam("numerics/feastol", 1e-9)
    model.setRealParam("numerics/dualfeastol", 1e-9)
    model.setIntParam("display/verblevel", 2)
    model.optimize()
    extracted = extract_solution(model, variables)
    payload = {
        "scope": "fixed terminal POVM and fixed prefix order",
        "formulation": "exact shared Choi variables with complex Cholesky positivity",
        "weight": args.weight,
        "terminal_effect_weights": [*weights.tolist(), 0.0],
        "prefix_order": list(order),
        "target": args.target,
        "gauge_fix": args.gauge_fix,
        "checkpoint": None if args.checkpoint is None else str(args.checkpoint),
        "seed_accepted": seed_accepted,
        "status": str(model.getStatus()),
        "primal_bound": float(model.getPrimalbound()),
        "dual_bound": float(model.getDualbound()),
        "gap": float(model.getGap()),
        "nodes": int(model.getNNodes()),
        "solving_time": float(model.getSolvingTime()),
        "solution": extracted,
        "certificate_scope_note": (
            "Exact nonconvex formulation; SCIP bounds remain numerical and "
            "solver-conditional until outward validation."
        ),
    }
    rendered = json.dumps(payload, indent=2) + "\n"
    print(rendered, end="")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered, encoding="utf-8")
    if args.solution_npz is not None:
        if extracted is None:
            raise RuntimeError("cannot write a checkpoint without a feasible solution")
        save_solution_arrays(
            args.solution_npz, extracted, np.asarray(variables["effects"])
        )


if __name__ == "__main__":
    main()

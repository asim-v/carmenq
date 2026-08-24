"""Small SOCP outer bound based on common-instrument Fourier contraction.

This relaxation discards the Choi matrix completely.  It keeps only input
Bloch vectors, flagged conditioned outputs, the fixed Helstrom face, and the
trace-norm contractions forced by one common instrument.  Consequently every
physical point is feasible, while angular cap covers make the reverse-convex
Bloch-dominated trace norm finite and certifiable.
"""

from __future__ import annotations

import numpy as np
import cvxpy as cp

from terminal_weight_upper import filled_effect_weights
from two_block_choi_seesaw import (
    OUTCOMES,
    PATHS,
    canonical_three_effect_povm,
    hellinger_hypograph,
)


class InfeasibleBehaviorOuter(RuntimeError):
    """The selected conic behavior cell is solver-conditionally infeasible."""


CHARACTERS = np.asarray(
    [
        [1.0, 1.0, -1.0, -1.0],
        [1.0, -1.0, 1.0, -1.0],
        [1.0, -1.0, -1.0, 1.0],
    ]
)


def solve_behavior_outer(
    branches: tuple[str, str, str],
    prior_box: np.ndarray,
    caps: tuple[object | None, ...] = (
        None,
        None,
        None,
    ),
    solver: str = "clarabel",
    build_only: bool = False,
    pairwise_contractions: tuple[dict[str, object], ...] = (),
) -> dict[str, object]:
    if len(branches) != 3 or len(caps) != 3:
        raise ValueError("three Fourier branches and three caps are required")
    effects = canonical_three_effect_povm(np.asarray([0.92, 0.64, 0.44]))
    traces = np.trace(effects, axis1=1, axis2=2).real
    directions = np.zeros((4, 3), dtype=float)
    active = tuple(int(s) for s in OUTCOMES if traces[s] > 1e-9)
    paulis = (
        np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex),
        np.array([[0.0, -1j], [1j, 0.0]], dtype=complex),
        np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex),
    )
    for s in active:
        projector = effects[s] / traces[s]
        directions[s] = [
            float(np.trace(projector @ pauli).real) for pauli in paulis
        ]

    probability = cp.Variable((4, 4), nonneg=True)
    output = [[cp.Variable(3) for _ in OUTCOMES] for _ in OUTCOMES]
    input_vector = [cp.Variable(3) for _ in OUTCOMES]
    constraints: list[cp.Constraint] = [cp.sum(probability) == 1.0]
    prior = [cp.sum(probability[z, :]) for z in OUTCOMES]
    for z in OUTCOMES:
        constraints.extend(
            (
                prior[z] >= float(prior_box[z, 0]),
                prior[z] <= float(prior_box[z, 1]),
                cp.SOC(prior[z], input_vector[z]),
            )
        )
        for y in OUTCOMES:
            constraints.append(cp.SOC(probability[z, y], output[z][y]))
    constraints.extend(prior[index] >= prior[index + 1] for index in range(3))

    syndrome_prior = [
        sum(probability[z, z ^ s] for z in OUTCOMES) for s in OUTCOMES
    ]
    terminal_vector = [
        sum((output[z][z ^ s] for z in OUTCOMES), cp.Constant(np.zeros(3)))
        for s in OUTCOMES
    ]
    for s in OUTCOMES:
        constraints.append(cp.SOC(syndrome_prior[s], terminal_vector[s]))
    audit = sum(
        0.5
        * traces[s]
        * (syndrome_prior[s] + directions[s] @ terminal_vector[s])
        for s in active
    )
    cap = filled_effect_weights(float(traces.max()))
    constraints.append(audit <= sum(cap[z] * prior[z] for z in OUTCOMES))
    dual_trace = cp.Variable(nonneg=True)
    dual_vector = cp.Variable(3)
    constraints.append(cp.SOC(dual_trace, dual_vector))
    for s in OUTCOMES:
        constraints.append(
            cp.SOC(
                dual_trace - syndrome_prior[s],
                dual_vector - terminal_vector[s],
            )
        )
    constraints.append(audit == dual_trace)
    for s in active:
        constraints.append(
            terminal_vector[s]
            == dual_vector + (dual_trace - syndrome_prior[s]) * directions[s]
        )

    flagged_values: list[cp.Expression] = []
    input_fourier_vectors: list[cp.Expression] = []
    bloch_seen = 0
    for character_index, (character, branch) in enumerate(
        zip(CHARACTERS, branches)
    ):
        scalar = sum(character[z] * prior[z] for z in OUTCOMES)
        vector = sum(
            (character[z] * input_vector[z] for z in OUTCOMES),
            cp.Constant(np.zeros(3)),
        )
        input_fourier_vectors.append(vector)
        block_norms = []
        for y in OUTCOMES:
            block_scalar = sum(
                character[z] * probability[z, y] for z in OUTCOMES
            )
            block_vector = sum(
                (character[z] * output[z][y] for z in OUTCOMES),
                cp.Constant(np.zeros(3)),
            )
            block_norm = cp.Variable(nonneg=True)
            constraints.extend(
                (
                    block_norm >= block_scalar,
                    block_norm >= -block_scalar,
                    cp.SOC(block_norm, block_vector),
                )
            )
            block_norms.append(block_norm)
        flagged = sum(block_norms)
        flagged_values.append(flagged)
        if branch == "scalar-positive":
            constraints.extend((cp.SOC(scalar, vector), flagged <= scalar))
        elif branch == "scalar-negative":
            constraints.extend((cp.SOC(-scalar, vector), flagged <= -scalar))
        elif branch == "bloch":
            if bloch_seen == 0:
                constraints.extend(
                    (
                        vector[0] == 0.0,
                        vector[1] == 0.0,
                        vector[2] >= 0.0,
                        flagged <= vector[2],
                    )
                )
            elif bloch_seen == 1:
                constraints.extend((vector[1] == 0.0, vector[0] >= 0.0))
            cap_data = caps[character_index]
            if cap_data is not None:
                if isinstance(cap_data, cp.Parameter):
                    if cap_data.shape != (3,):
                        raise ValueError("a scaled cap parameter must have shape (3,)")
                    projection_bound = cap_data @ vector
                else:
                    cap_array = np.asarray(cap_data, dtype=float)
                    normal = cap_array[:3]
                    cosine = float(cap_array[3])
                    projection_bound = normal @ vector / cosine
                constraints.extend(
                    (
                        cp.SOC(projection_bound, vector),
                        flagged <= projection_bound,
                    )
                )
            bloch_seen += 1
        else:
            raise ValueError(f"unknown Fourier branch {branch}")

    pairwise_flagged: list[cp.Expression] = []
    for contraction in pairwise_contractions:
        if "coefficients" in contraction:
            coefficients = np.asarray(contraction["coefficients"], dtype=float)
            if coefficients.shape != (4,) or np.linalg.norm(coefficients) <= 1e-14:
                raise ValueError("invalid general contraction coefficients")
        else:
            first_z, second_z = tuple(contraction["pair"])
            scale = float(contraction["scale"])
            if (
                first_z not in OUTCOMES
                or second_z not in OUTCOMES
                or first_z == second_z
                or not np.isfinite(scale)
                or scale < 0.0
            ):
                raise ValueError("invalid pairwise contraction")
            coefficients = np.zeros(4)
            coefficients[first_z] = 1.0
            coefficients[second_z] = -scale
        branch = str(contraction["branch"])
        scalar = sum(float(coefficients[z]) * prior[z] for z in OUTCOMES)
        vector = sum(
            (float(coefficients[z]) * input_vector[z] for z in OUTCOMES),
            cp.Constant(np.zeros(3)),
        )
        block_norms = []
        for y in OUTCOMES:
            block_scalar = sum(
                float(coefficients[z]) * probability[z, y] for z in OUTCOMES
            )
            block_vector = sum(
                (float(coefficients[z]) * output[z][y] for z in OUTCOMES),
                cp.Constant(np.zeros(3)),
            )
            block_norm = cp.Variable(nonneg=True)
            constraints.extend(
                (
                    block_norm >= block_scalar,
                    block_norm >= -block_scalar,
                    cp.SOC(block_norm, block_vector),
                )
            )
            block_norms.append(block_norm)
        flagged = sum(block_norms)
        pairwise_flagged.append(flagged)
        if branch == "scalar-positive":
            constraints.extend((cp.SOC(scalar, vector), flagged <= scalar))
        elif branch == "scalar-negative":
            constraints.extend((cp.SOC(-scalar, vector), flagged <= -scalar))
        elif branch == "bloch":
            cap_data = contraction.get("cap")
            if isinstance(cap_data, cp.Parameter):
                if cap_data.shape != (3,):
                    raise ValueError("a scaled pair cap parameter has shape (3,)")
                projection_bound = cap_data @ vector
            elif cap_data is not None:
                cap_array = np.asarray(cap_data, dtype=float)
                projection_bound = cap_array[:3] @ vector / float(cap_array[3])
            else:
                raise ValueError("a vector-active pairwise contraction needs a cap")
            constraints.extend(
                (
                    cp.SOC(projection_bound, vector),
                    flagged <= projection_bound,
                )
            )
        else:
            raise ValueError(f"unknown pairwise branch {branch}")

    returned = hellinger_hypograph(
        [probability[z, y] for z, y in PATHS], constraints
    )
    problem = cp.Problem(cp.Maximize(0.55 * audit + 0.45 * returned), constraints)
    if build_only:
        return {
            "problem": problem,
            "prior_expressions": prior,
            "audit_expression": audit,
            "return_expression": returned,
            "flagged_expressions": flagged_values,
            "input_fourier_expressions": input_fourier_vectors,
            "pairwise_flagged_expressions": pairwise_flagged,
        }
    solver_used = solver
    if solver == "clarabel":
        try:
            problem.solve(
                solver="CLARABEL",
                tol_gap_abs=1e-9,
                tol_gap_rel=1e-9,
                tol_feas=1e-9,
                max_iter=1000,
            )
        except cp.SolverError:
            problem.solve(solver="SCS", eps=1e-7, max_iters=200_000)
            solver_used = "scs_after_clarabel_error"
    elif solver == "scs":
        problem.solve(solver="SCS", eps=1e-7, max_iters=200_000)
    else:
        raise ValueError("solver must be clarabel or scs")
    if problem.status in {cp.INFEASIBLE, cp.INFEASIBLE_INACCURATE}:
        raise InfeasibleBehaviorOuter(problem.status)
    if problem.status not in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}:
        raise RuntimeError(f"Fourier behavior outer failed: {problem.status}")
    return {
        "bound": float(problem.value),
        "status": problem.status,
        "branches": list(branches),
        "caps": [
            None
            if cap is None
            else "parameter"
            if isinstance(cap, cp.Parameter)
            else list(cap)
            for cap in caps
        ],
        "prior": [float(item.value) for item in prior],
        "audit": float(audit.value),
        "return": float(returned.value),
        "flagged_norms": [float(item.value) for item in flagged_values],
        "input_fourier_vectors": [
            np.asarray(item.value, dtype=float).tolist()
            for item in input_fourier_vectors
        ],
        "pairwise_flagged_norms": [
            float(item.value) for item in pairwise_flagged
        ],
        "path_probabilities": np.asarray(probability.value, dtype=float).tolist(),
        "conditioned_output_bloch_vectors": [
            [np.asarray(output[z][y].value, dtype=float).tolist() for y in OUTCOMES]
            for z in OUTCOMES
        ],
        "prefix_bloch_coefficients": [
            [float(prior[z].value), *np.asarray(input_vector[z].value).tolist()]
            for z in OUTCOMES
        ],
        "solver": solver_used,
        "iterations": problem.solver_stats.num_iters,
        "solve_time": problem.solver_stats.solve_time,
    }

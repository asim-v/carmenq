"""Exact and finite-statistics tests for the audit--return benchmark."""

from __future__ import annotations

import json

import numpy as np
import pytest

from reversible_histories.audit_return import (
    BenchmarkCounts,
    PhenomenologicalNoise,
    certify_classical_memory,
    classical_memory_bound,
    classical_memory_frontier,
    collective_classical_record_bound,
    plan_experiment,
    return_curve,
    simulate_counts,
)
from reversible_histories.audit_return_cli import main as cli_main


def test_balanced_classical_bound_and_collective_comparator() -> None:
    for n_steps in range(2, 9):
        point = classical_memory_frontier(n_steps, 0.5)
        assert np.isclose(point.score, 0.75, atol=1e-12)
        assert np.isclose(point.local_strength, 0.0, atol=1e-12)
    assert np.isclose(
        collective_classical_record_bound(0.5),
        (1.0 + 1.0 / np.sqrt(2.0)) / 2.0,
    )


def test_one_and_two_slot_closed_forms() -> None:
    for weight in np.linspace(0.0, 1.0, 101):
        expected_one = 0.5 + 0.5 * np.sqrt(weight**2 + (1.0 - weight) ** 2)
        assert np.isclose(classical_memory_bound(1, float(weight)), expected_one)
        expected_two = (
            1.0 - weight / 2.0
            if weight <= 0.5
            else weight / 2.0 + weight**2 / (3.0 * weight - 1.0)
        )
        assert np.isclose(classical_memory_bound(2, float(weight)), expected_two)


def test_first_order_transition_for_three_steps() -> None:
    below = classical_memory_frontier(3, 0.624)
    above = classical_memory_frontier(3, 0.626)
    assert below.strategy == "no_record"
    assert np.isclose(below.local_strength, 0.0)
    assert above.strategy == "equal_weak_measurement"
    assert above.local_strength > 0.97


def test_return_curve_endpoints() -> None:
    assert np.isclose(return_curve(0.0), 1.0)
    assert np.isclose(return_curve(1.0), 0.5)
    with pytest.raises(ValueError):
        return_curve(1.1)


def test_perfect_coherent_data_certifies_but_classical_boundary_does_not() -> None:
    perfect = certify_classical_memory(
        BenchmarkCounts(5_000, 5_000, 5_000, 5_000),
        n_steps=4,
        alpha=0.01,
    )
    assert perfect.certified
    assert perfect.margin > 0.2

    boundary = certify_classical_memory(
        BenchmarkCounts(2_500, 5_000, 5_000, 5_000),
        n_steps=4,
        alpha=0.01,
    )
    assert np.isclose(boundary.observed_score, 0.75)
    assert not boundary.certified


def test_systematic_penalty_and_null_slack_are_conservative() -> None:
    counts = BenchmarkCounts(900, 1_000, 900, 1_000)
    nominal = certify_classical_memory(counts, 3, alpha=0.01)
    robust = certify_classical_memory(
        counts,
        3,
        alpha=0.01,
        audit_systematic=0.02,
        return_systematic=0.03,
        null_slack=0.01,
    )
    assert robust.margin < nominal.margin
    assert np.isclose(robust.systematic_penalty, 0.025)


def test_power_plan_guarantees_declared_gap() -> None:
    plan = plan_experiment(
        n_steps=8,
        audit_probability=0.97,
        return_fidelity=0.95,
        alpha=0.01,
        beta=0.1,
        audit_systematic=0.005,
        return_systematic=0.005,
        null_slack=0.005,
    )
    assert plan.feasible
    assert plan.total_trials > 0
    assert plan.alpha_radius + plan.beta_radius < plan.adjusted_gap

    impossible = plan_experiment(3, 0.5, 1.0, alpha=0.01, beta=0.1)
    assert not impossible.feasible


def test_seeded_simulation_is_reproducible() -> None:
    first = simulate_counts(0.91, 0.87, 1_000, 1_200, seed=17)
    second = simulate_counts(0.91, 0.87, 1_000, 1_200, seed=17)
    assert first == second


def test_phenomenological_noise_is_monotone() -> None:
    model = PhenomenologicalNoise()
    points = [model.point(n_steps) for n_steps in range(1, 20)]
    assert np.all(np.diff([point[0] for point in points]) < 0.0)
    assert np.all(np.diff([point[1] for point in points]) < 0.0)


def test_cli_emits_strict_json_for_bounds_and_infeasible_plans(capsys) -> None:
    assert cli_main(["bound", "--steps", "8"]) == 0
    bound = json.loads(capsys.readouterr().out)
    assert np.isclose(bound["streaming_classical_memory"]["score"], 0.75)

    assert (
        cli_main(
            [
                "plan",
                "--steps",
                "8",
                "--audit-probability",
                "0.5",
                "--return-fidelity",
                "1.0",
            ]
        )
        == 0
    )
    plan = json.loads(capsys.readouterr().out)
    assert plan["feasible"] is False
    assert plan["alpha_radius"] is None

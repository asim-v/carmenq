"""Statistical benchmark for causal audit--return quantum memory.

The benchmark is a resource-restricted test.  Fresh halves of EPR pairs are
sent through a sequential device and every returned carrier is sequestered.
After the complete prefix and terminal-memory commitment, the verifier chooses
between predicting the parity of later Z measurements (AUDIT) and returning
every EPR pair while visibly resetting that memory port (RETURN).

This module implements the exact support function for the adaptive
classical-memory null, a fixed-sample finite-statistics certificate, power
planning, and a transparent phenomenological noise forecast.  It does not
claim device-independent certification.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from math import ceil, exp, log, sqrt

import numpy as np
from scipy.optimize import brentq


def _probability(value: float, name: str) -> float:
    value = float(value)
    if not 0.0 <= value <= 1.0:
        raise ValueError(f"{name} must lie in [0, 1]")
    return value


def return_curve(strength: float | np.ndarray) -> float | np.ndarray:
    """One-slot return fidelity at classical audit strength ``strength``."""
    values = np.asarray(strength, dtype=float)
    if np.any((values < 0.0) | (values > 1.0)):
        raise ValueError("strength must lie in [0, 1]")
    result = (1.0 + np.sqrt(np.maximum(0.0, 1.0 - values**2))) / 2.0
    if result.ndim == 0:
        return float(result)
    return result


@dataclass(frozen=True)
class FrontierPoint:
    """One exposed point on an audit--return support curve."""

    audit_probability: float
    return_fidelity: float
    score: float
    local_strength: float
    strategy: str


def score(
    audit_probability: float, return_fidelity: float, audit_weight: float = 0.5
) -> float:
    """Weighted audit--return score."""
    audit_probability = _probability(audit_probability, "audit_probability")
    return_fidelity = _probability(return_fidelity, "return_fidelity")
    audit_weight = _probability(audit_weight, "audit_weight")
    return audit_weight * audit_probability + (1.0 - audit_weight) * return_fidelity


def _streaming_point(n_steps: int, audit_weight: float, y: float) -> FrontierPoint:
    strength = 2.0 * y / (1.0 + y * y)
    audit = (1.0 + strength**n_steps) / 2.0
    returned = (1.0 / (1.0 + y * y)) ** n_steps
    value = score(audit, returned, audit_weight)
    if y == 0.0:
        strategy = "no_record"
    elif y == 1.0:
        strategy = "projective_each_slot"
    else:
        strategy = "equal_weak_measurement"
    return FrontierPoint(audit, returned, value, strength, strategy)


def classical_memory_frontier(
    n_steps: int, audit_weight: float = 0.5
) -> FrontierPoint:
    """Exact support point for an adaptive classical-memory streaming comb.

    The allowed null has arbitrary local instruments, disposable within-slot
    ancillas, unlimited adaptive classical memory, and transcript-conditioned
    joint recovery, but no coherent quantum state or shared entanglement
    persisting between slots.  Previous returned carriers are inaccessible.
    """
    if int(n_steps) != n_steps or n_steps < 1:
        raise ValueError("n_steps must be a positive integer")
    n_steps = int(n_steps)
    audit_weight = _probability(audit_weight, "audit_weight")

    if n_steps == 1:
        norm = sqrt(audit_weight**2 + (1.0 - audit_weight) ** 2)
        if norm == 0.0:  # Unreachable, but keeps the endpoint explicit.
            strength = 0.0
        else:
            strength = audit_weight / norm
        audit = (1.0 + strength) / 2.0
        returned = float(return_curve(strength))
        return FrontierPoint(
            audit, returned, score(audit, returned, audit_weight), strength, "one_slot"
        )

    candidates = [
        _streaming_point(n_steps, audit_weight, 0.0),
        _streaming_point(n_steps, audit_weight, 1.0),
    ]

    if n_steps == 2:
        if audit_weight > 0.5:
            y = sqrt(max(0.0, 2.0 - 1.0 / audit_weight))
            candidates.append(_streaming_point(n_steps, audit_weight, y))
    elif 0.0 < audit_weight < 1.0:
        turning_point = sqrt((n_steps - 2.0) / n_steps)

        def stationary(y: float) -> float:
            return (
                audit_weight
                * 2.0 ** (n_steps - 2)
                * y ** (n_steps - 2)
                * (1.0 - y * y)
                - (1.0 - audit_weight)
            )

        at_turn = stationary(turning_point)
        if at_turn >= 0.0:
            # The larger root is the only nonzero local maximum.
            high_root = brentq(stationary, turning_point, 1.0, xtol=1e-14)
            candidates.append(_streaming_point(n_steps, audit_weight, high_root))

    # Stable tie breaking prefers the least disturbing strategy.
    return max(candidates, key=lambda item: (round(item.score, 14), -item.local_strength))


def classical_memory_bound(n_steps: int, audit_weight: float = 0.5) -> float:
    """Exact classical-memory support value."""
    return classical_memory_frontier(n_steps, audit_weight).score


def collective_classical_record_bound(audit_weight: float = 0.5) -> float:
    """Support for a joint instrument retaining only a classical outcome."""
    audit_weight = _probability(audit_weight, "audit_weight")
    return 0.5 + 0.5 * sqrt(audit_weight**2 + (1.0 - audit_weight) ** 2)


@dataclass(frozen=True)
class BenchmarkCounts:
    """Fixed-sample binary outcomes from the two late challenge branches."""

    audit_successes: int
    audit_trials: int
    return_successes: int
    return_trials: int

    def __post_init__(self) -> None:
        for successes, trials, label in (
            (self.audit_successes, self.audit_trials, "audit"),
            (self.return_successes, self.return_trials, "return"),
        ):
            if int(successes) != successes or int(trials) != trials:
                raise ValueError(f"{label} counts must be integers")
            if trials <= 0:
                raise ValueError(f"{label}_trials must be positive")
            if not 0 <= successes <= trials:
                raise ValueError(f"{label}_successes must lie in [0, trials]")

    @property
    def audit_probability(self) -> float:
        return self.audit_successes / self.audit_trials

    @property
    def return_fidelity(self) -> float:
        return self.return_successes / self.return_trials


@dataclass(frozen=True)
class CertificationResult:
    """A conservative one-sided fixed-sample certification result."""

    n_steps: int
    audit_weight: float
    alpha: float
    observed_audit_probability: float
    observed_return_fidelity: float
    observed_score: float
    systematic_penalty: float
    adjusted_score: float
    confidence_radius: float
    lower_confidence_score: float
    classical_bound: float
    null_slack: float
    certification_threshold: float
    certified: bool
    margin: float
    p_value_upper_bound: float
    effective_shots: float

    def to_dict(self) -> dict[str, float | int | bool]:
        return asdict(self)


def weighted_hoeffding_radius(
    alpha: float,
    audit_trials: int,
    return_trials: int,
    audit_weight: float = 0.5,
) -> float:
    """One-sided Hoeffding radius for a weighted sum of two sample means."""
    alpha = float(alpha)
    if not 0.0 < alpha < 1.0:
        raise ValueError("alpha must lie in (0, 1)")
    if audit_trials <= 0 or return_trials <= 0:
        raise ValueError("Both trial counts must be positive")
    audit_weight = _probability(audit_weight, "audit_weight")
    denominator = (
        audit_weight**2 / audit_trials
        + (1.0 - audit_weight) ** 2 / return_trials
    )
    return sqrt(0.5 * denominator * log(1.0 / alpha))


def certify_classical_memory(
    counts: BenchmarkCounts,
    n_steps: int,
    audit_weight: float = 0.5,
    alpha: float = 0.01,
    *,
    audit_systematic: float = 0.0,
    return_systematic: float = 0.0,
    null_slack: float = 0.0,
) -> CertificationResult:
    """Test the adaptive classical-memory null with fixed sample sizes.

    ``audit_systematic`` and ``return_systematic`` are declared upper bounds
    on upward biases of the corresponding observed probabilities.  They are
    subtracted before inference.  ``null_slack`` enlarges the null for source,
    isolation, or model uncertainty.  These quantities must be established by
    calibration or a separate proof; fitting them to the observed data would
    invalidate the stated confidence level.
    """
    audit_weight = _probability(audit_weight, "audit_weight")
    audit_systematic = _probability(audit_systematic, "audit_systematic")
    return_systematic = _probability(return_systematic, "return_systematic")
    null_slack = _probability(null_slack, "null_slack")
    radius = weighted_hoeffding_radius(
        alpha, counts.audit_trials, counts.return_trials, audit_weight
    )
    observed = score(
        counts.audit_probability, counts.return_fidelity, audit_weight
    )
    systematic_penalty = (
        audit_weight * audit_systematic
        + (1.0 - audit_weight) * return_systematic
    )
    adjusted = observed - systematic_penalty
    lower = adjusted - radius
    bound = classical_memory_bound(n_steps, audit_weight)
    threshold = min(1.0, bound + null_slack)
    margin = lower - threshold
    denominator = (
        audit_weight**2 / counts.audit_trials
        + (1.0 - audit_weight) ** 2 / counts.return_trials
    )
    raw_margin = max(0.0, adjusted - threshold)
    p_upper = min(1.0, exp(-2.0 * raw_margin**2 / denominator))
    return CertificationResult(
        n_steps=int(n_steps),
        audit_weight=audit_weight,
        alpha=float(alpha),
        observed_audit_probability=counts.audit_probability,
        observed_return_fidelity=counts.return_fidelity,
        observed_score=observed,
        systematic_penalty=systematic_penalty,
        adjusted_score=adjusted,
        confidence_radius=radius,
        lower_confidence_score=lower,
        classical_bound=bound,
        null_slack=null_slack,
        certification_threshold=threshold,
        certified=bool(margin > 0.0),
        margin=margin,
        p_value_upper_bound=p_upper,
        effective_shots=1.0 / denominator,
    )


@dataclass(frozen=True)
class PowerPlan:
    """A Hoeffding-guaranteed fixed-sample plan for a declared alternative."""

    n_steps: int
    audit_weight: float
    audit_probability: float
    return_fidelity: float
    alpha: float
    beta: float
    alternative_score: float
    adjusted_gap: float
    audit_trials: int
    return_trials: int
    total_trials: int
    alpha_radius: float | None
    beta_radius: float | None
    feasible: bool

    def to_dict(self) -> dict[str, float | int | bool | None]:
        return asdict(self)


def plan_experiment(
    n_steps: int,
    audit_probability: float,
    return_fidelity: float,
    audit_weight: float = 0.5,
    alpha: float = 0.01,
    beta: float = 0.1,
    *,
    audit_systematic: float = 0.0,
    return_systematic: float = 0.0,
    null_slack: float = 0.0,
    max_total_trials: int = 100_000_000,
) -> PowerPlan:
    """Plan shots with false-positive ``alpha`` and power at least ``1-beta``.

    The guarantee follows from two one-sided Hoeffding bounds and assumes
    independent Bernoulli trials, a fixed sample size, and preregistered
    ``n_steps`` and ``audit_weight``.  Trial allocation is asymptotically
    optimal: audit and return fractions approach ``audit_weight`` and
    ``1-audit_weight``.
    """
    audit_probability = _probability(audit_probability, "audit_probability")
    return_fidelity = _probability(return_fidelity, "return_fidelity")
    audit_weight = _probability(audit_weight, "audit_weight")
    if not 0.0 < audit_weight < 1.0:
        raise ValueError("Power planning requires 0 < audit_weight < 1")
    if not 0.0 < alpha < 1.0:
        raise ValueError("alpha must lie in (0, 1)")
    if not 0.0 < beta < 1.0:
        raise ValueError("beta must lie in (0, 1)")
    if int(max_total_trials) != max_total_trials or max_total_trials < 2:
        raise ValueError("max_total_trials must be an integer of at least two")
    audit_systematic = _probability(audit_systematic, "audit_systematic")
    return_systematic = _probability(return_systematic, "return_systematic")
    null_slack = _probability(null_slack, "null_slack")
    alternative = score(audit_probability, return_fidelity, audit_weight)
    penalty = (
        audit_weight * audit_systematic
        + (1.0 - audit_weight) * return_systematic
    )
    threshold = classical_memory_bound(n_steps, audit_weight) + null_slack
    gap = alternative - penalty - threshold
    if gap <= 0.0:
        return PowerPlan(
            int(n_steps),
            audit_weight,
            audit_probability,
            return_fidelity,
            float(alpha),
            float(beta),
            alternative,
            gap,
            0,
            0,
            0,
            None,
            None,
            False,
        )

    leading = (
        sqrt(log(1.0 / alpha)) + sqrt(log(1.0 / beta))
    ) ** 2 / (2.0 * gap**2)
    total = max(2, int(ceil(leading)))
    while total <= max_total_trials:
        audit_trials = max(1, int(round(audit_weight * total)))
        return_trials = max(1, total - audit_trials)
        alpha_radius = weighted_hoeffding_radius(
            alpha, audit_trials, return_trials, audit_weight
        )
        beta_radius = weighted_hoeffding_radius(
            beta, audit_trials, return_trials, audit_weight
        )
        if alpha_radius + beta_radius < gap:
            return PowerPlan(
                int(n_steps),
                audit_weight,
                audit_probability,
                return_fidelity,
                float(alpha),
                float(beta),
                alternative,
                gap,
                audit_trials,
                return_trials,
                audit_trials + return_trials,
                alpha_radius,
                beta_radius,
                True,
            )
        total += 1
    return PowerPlan(
        int(n_steps),
        audit_weight,
        audit_probability,
        return_fidelity,
        float(alpha),
        float(beta),
        alternative,
        gap,
        0,
        0,
        0,
        None,
        None,
        False,
    )


@dataclass(frozen=True)
class PhenomenologicalNoise:
    """Explicit, replaceable forecast model rather than hardware calibration.

    Audit contrast is ``audit_base_contrast * audit_step_contrast**n`` and
    RETURN fidelity is ``return_base_fidelity * return_step_fidelity**n``.
    This simple multiplicative law is useful for planning and sensitivity
    analysis; experimental claims must use measured branch probabilities.
    """

    audit_base_contrast: float = 0.995
    audit_step_contrast: float = 0.998
    return_base_fidelity: float = 0.995
    return_step_fidelity: float = 0.997

    def __post_init__(self) -> None:
        for name in (
            "audit_base_contrast",
            "audit_step_contrast",
            "return_base_fidelity",
            "return_step_fidelity",
        ):
            _probability(getattr(self, name), name)

    def point(self, n_steps: int) -> tuple[float, float]:
        if int(n_steps) != n_steps or n_steps < 1:
            raise ValueError("n_steps must be a positive integer")
        audit = 0.5 * (
            1.0
            + self.audit_base_contrast * self.audit_step_contrast ** int(n_steps)
        )
        returned = self.return_base_fidelity * self.return_step_fidelity ** int(
            n_steps
        )
        return float(audit), float(returned)


def simulate_counts(
    audit_probability: float,
    return_fidelity: float,
    audit_trials: int,
    return_trials: int,
    seed: int = 20260812,
) -> BenchmarkCounts:
    """Draw a deterministic-seed synthetic fixed-sample data set."""
    audit_probability = _probability(audit_probability, "audit_probability")
    return_fidelity = _probability(return_fidelity, "return_fidelity")
    if audit_trials <= 0 or return_trials <= 0:
        raise ValueError("Both trial counts must be positive")
    rng = np.random.default_rng(seed)
    return BenchmarkCounts(
        int(rng.binomial(audit_trials, audit_probability)),
        int(audit_trials),
        int(rng.binomial(return_trials, return_fidelity)),
        int(return_trials),
    )

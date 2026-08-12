"""Reference simulations for reversible quantum-history interferometry."""

from .audit_return import (
    BenchmarkCounts,
    CertificationResult,
    FrontierPoint,
    PhenomenologicalNoise,
    PowerPlan,
    certify_classical_memory,
    classical_memory_bound,
    classical_memory_frontier,
    collective_classical_record_bound,
    plan_experiment,
    return_curve,
    score,
    simulate_counts,
    weighted_hoeffding_radius,
)
from .protocol import (
    NoiseModel,
    ProtocolConfig,
    ProtocolResult,
    conditional_record_information,
    environment_conditional_information,
    run_protocol,
)

__all__ = [
    "BenchmarkCounts",
    "CertificationResult",
    "FrontierPoint",
    "NoiseModel",
    "PhenomenologicalNoise",
    "PowerPlan",
    "ProtocolConfig",
    "ProtocolResult",
    "certify_classical_memory",
    "classical_memory_bound",
    "classical_memory_frontier",
    "collective_classical_record_bound",
    "conditional_record_information",
    "environment_conditional_information",
    "plan_experiment",
    "return_curve",
    "run_protocol",
    "score",
    "simulate_counts",
    "weighted_hoeffding_radius",
]

__version__ = "1.0.0"

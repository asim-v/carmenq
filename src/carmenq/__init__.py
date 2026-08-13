"""CARMEN-Q: causal audit-return quantum-memory benchmarks.

The top-level API exposes concise names for routine use while retaining the
long-form scientific function names used in the accompanying manuscript.
"""

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

# Concise public aliases. The descriptive originals remain part of the API.
streaming_bound = classical_memory_bound
collective_bound = collective_classical_record_bound
certify = certify_classical_memory
plan = plan_experiment

__all__ = [
    "BenchmarkCounts",
    "CertificationResult",
    "FrontierPoint",
    "NoiseModel",
    "PhenomenologicalNoise",
    "PowerPlan",
    "ProtocolConfig",
    "ProtocolResult",
    "certify",
    "certify_classical_memory",
    "classical_memory_bound",
    "classical_memory_frontier",
    "collective_bound",
    "collective_classical_record_bound",
    "conditional_record_information",
    "environment_conditional_information",
    "plan",
    "plan_experiment",
    "return_curve",
    "run_protocol",
    "score",
    "simulate_counts",
    "streaming_bound",
    "weighted_hoeffding_radius",
]

__version__ = "2.0.2"

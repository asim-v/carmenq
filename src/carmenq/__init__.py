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
from .order_sensitive import (
    GROUPED_CHECK_MATRIX,
    INTERLEAVED_CHECK_MATRIX,
    INTERLEAVED_PERFECT_AUDIT_ENDPOINT,
    GroupedFrontierPoint,
    PerfectAuditEndpoint,
    full_crossing_cuts,
    gf2_rank,
    grouped_frontier,
    rank_two_static_qubit_support,
    trellis_connectivity_profile,
    trellis_connectivity_tau,
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
    "GROUPED_CHECK_MATRIX",
    "GroupedFrontierPoint",
    "INTERLEAVED_CHECK_MATRIX",
    "INTERLEAVED_PERFECT_AUDIT_ENDPOINT",
    "NoiseModel",
    "PhenomenologicalNoise",
    "PowerPlan",
    "PerfectAuditEndpoint",
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
    "full_crossing_cuts",
    "gf2_rank",
    "grouped_frontier",
    "plan",
    "plan_experiment",
    "return_curve",
    "rank_two_static_qubit_support",
    "run_protocol",
    "score",
    "simulate_counts",
    "streaming_bound",
    "trellis_connectivity_profile",
    "trellis_connectivity_tau",
    "weighted_hoeffding_radius",
]

__version__ = "2.1.0"

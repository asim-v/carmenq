"""Minimal CARMEN-Q example."""

from carmenq import BenchmarkCounts, certify, collective_bound, streaming_bound


n_steps = 8
weight = 0.5

print(f"Streaming classical ceiling: {streaming_bound(n_steps, weight):.6f}")
print(f"Collective classical ceiling: {collective_bound(weight):.6f}")

result = certify(
    BenchmarkCounts(
        audit_successes=9700,
        audit_trials=10000,
        return_successes=9500,
        return_trials=10000,
    ),
    n_steps=n_steps,
    audit_weight=weight,
    alpha=0.01,
    audit_systematic=0.005,
    return_systematic=0.005,
    null_slack=0.005,
)

print(f"Certified: {result.certified}")
print(f"Conservative margin: {result.margin:.6f}")

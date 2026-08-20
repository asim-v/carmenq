# Python API

CARMEN-Q exposes a compact top-level API for the common audit-return workflow. The same objects remain available under their long-form scientific names so published analyses can state exactly which bound or test they use.

| Task | Concise API | Scientific name |
|---|---|---|
| Streaming classical ceiling | `streaming_bound(n_steps, audit_weight)` | `classical_memory_bound(...)` |
| Collective classical ceiling | `collective_bound(audit_weight)` | `collective_classical_record_bound(...)` |
| Fixed-sample certificate | `certify(counts, n_steps, ...)` | `certify_classical_memory(...)` |
| Power planning | `plan(n_steps, audit_probability, return_fidelity, ...)` | `plan_experiment(...)` |
| Binary cut rank | `gf2_rank(matrix)` | `carmenq.order_sensitive.gf2_rank(...)` |
| Ordered cut profile | `trellis_connectivity_profile(matrix)` | same |
| Exact grouped point | `grouped_frontier(audit_weight)` | same |
| Interleaved analytic lower bound | `interleaved_candidate_lower_bound(audit_weight)` | same |
| Verified balanced counterexample | `INTERLEAVED_BALANCED_COUNTEREXAMPLE` | same |

```python
from carmenq import BenchmarkCounts, certify, streaming_bound

ceiling = streaming_bound(n_steps=8, audit_weight=0.5)
result = certify(
    BenchmarkCounts(9700, 10000, 9500, 10000),
    n_steps=8,
    audit_weight=0.5,
    alpha=0.01,
)
```

`certify` returns an immutable `CertificationResult` dataclass containing the observed score, systematic penalty, confidence radius, enlarged null threshold, decision, and conservative margin. It can be serialized with `result.to_dict()`.

The lower-level modules `carmenq.protocol`, `carmenq.linalg`, and `carmenq.experiments` implement the exact density-matrix reference protocol and deterministic reproduction pipeline. They are public for research use, but the top-level API is the stability boundary for ordinary users.

The `carmenq.order_sensitive` module contains the separate rank-two syndrome result. `GROUPED_CHECK_MATRIX` and `INTERLEAVED_CHECK_MATRIX` define the canonical coordinate orders. `grouped_frontier` evaluates the exact attainable grouped boundary, while `INTERLEAVED_PERFECT_AUDIT_ENDPOINT` records the proved and attained interleaved endpoint. `interleaved_candidate_scores(q, v)` evaluates an exact two-parameter streamed construction, and `interleaved_candidate_lower_bound` optimizes that construction deterministically. Its result keeps `support_is_globally_optimal=False`: a stored finite-outcome non-QND instrument now strictly exceeds this family, and the arbitrary-instrument interior frontier remains open. Run `python scripts/verify_interleaved_counterexample.py` to reproduce that correction independently.

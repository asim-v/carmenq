# Python API

CARMEN-Q exposes a compact top-level API for the common audit-return workflow. The same objects remain available under their long-form scientific names so published analyses can state exactly which bound or test they use.

| Task | Concise API | Scientific name |
|---|---|---|
| Streaming classical ceiling | `streaming_bound(n_steps, audit_weight)` | `classical_memory_bound(...)` |
| Collective classical ceiling | `collective_bound(audit_weight)` | `collective_classical_record_bound(...)` |
| Fixed-sample certificate | `certify(counts, n_steps, ...)` | `certify_classical_memory(...)` |
| Power planning | `plan(n_steps, audit_probability, return_fidelity, ...)` | `plan_experiment(...)` |

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

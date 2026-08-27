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
| Full-crossing square-law bound | `full_crossing_perfect_audit_return_bound(rank, dimension, alphabet_size=2)` | same |
| Full-rank-block temporal bound | `full_rank_block_perfect_audit_return_bound(rank, dimension, block_count, alphabet_size=2)` | same |
| Approximate-AUDIT temporal bound | `full_rank_block_approximate_audit_return_bound(P_A, rank, dimension, block_count, alphabet_size=2)` | same |
| Ordered binary block packing | `full_rank_block_packing_number(matrix)` | same |
| Ordered binary endpoint bound | `ordered_check_perfect_audit_return_bound(matrix, dimension)` | same |
| Exact grouped point | `grouped_frontier(audit_weight)` | same |
| Interleaved analytic lower bound | `interleaved_candidate_lower_bound(audit_weight)` | same |
| Compact MPS lower bound | `interleaved_compact_lower_bound(audit_weight)` | same |
| Interleaved RETURN certificate | `interleaved_return_upper_bound(P_A)` | same |
| Interleaved support certificate | `interleaved_support_upper_bound(audit_weight)` | same |
| Certified interior interval at weight 0.6 | `INTERLEAVED_L060_CERTIFIED_INTERVAL` | same |
| Verified balanced counterexample | `INTERLEAVED_BALANCED_COUNTEREXAMPLE` | same |
| Flagged common-instrument cuts | `scan_flagged_trace_norm_cuts(states, outputs, scales=None)` | same |
| Exact fixed-input instrument projection | `project_to_common_instrument(states, outputs)` | same |
| Exact four-input basis reconstruction | `reconstruct_common_instrument_from_basis(states, outputs)` | same |
| Exact four-input effective POVM reconstruction | `reconstruct_effective_povm_from_basis(states, probabilities)` | same |

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

`scan_flagged_trace_norm_cuts` evaluates the necessary inequalities
`sum_y ||sigma[z,y] - t sigma[z',y]||_1 <= ||rho[z] - t rho[z']||_1` for one
flagged channel shared by every input. `project_to_common_instrument` goes
further: for fixed prefix states it projects the complete output family onto
the exact Choi-compatible set and, when the distance is nonzero, verifies a
linear separating witness with a second support SDP. The projection requires
the optional `frontier` dependencies. Its numerical witness is
solver-conditional rather than an interval-arithmetic certificate.

When exactly four prefix states span the Hermitian qubit operators,
`reconstruct_common_instrument_from_basis` needs no optimizer. It reconstructs
the unique Pauli-transfer and Choi matrices compatible with all sixteen
conditioned outputs. The family comes from one instrument precisely when the
four Choi matrices are positive semidefinite and their sum is trace preserving.
The result also contains determinant-scaled Choi numerators, which expose the
criterion as polynomial matrix inequalities for global-optimization work.
Singular input families should continue to use `project_to_common_instrument`.

If only the measured probability table is needed,
`reconstruct_effective_povm_from_basis` reconstructs the unique common POVM
seen by the four input states. The returned `BasisPovmReconstruction` contains
all effect matrices, their minimum eigenvalues, completeness and interpolation
residuals, and determinant-scaled positive-effect numerators. Every common
instrument followed by a terminal measurement must pass this test. Passing it
does not by itself prove that the POVM factors through the specified terminal
measurement and instrument.

The `carmenq.order_sensitive` module contains the separate syndrome-order result. `GROUPED_CHECK_MATRIX` and `INTERLEAVED_CHECK_MATRIX` define the canonical coordinate orders. `grouped_frontier` evaluates the exact attainable grouped boundary, while `INTERLEAVED_PERFECT_AUDIT_ENDPOINT` records the proved and attained interleaved endpoint. `full_rank_block_perfect_audit_return_bound` evaluates the finite-field theorem `F_R <= min(1,(d/q**r)**m)` for `m` consecutive full-rank blocks; `full_crossing_perfect_audit_return_bound` is its two-block specialization, and `q` must be a prime power. The robust function `full_rank_block_approximate_audit_return_bound` applies the causal-list/Ky-Fan theorem when `P_A < 1`. For the canonical interleaved matrix, `interleaved_return_upper_bound` and `interleaved_support_upper_bound` give its RETURN and support certificates; the latter is strictly below the grouped/static support for `3/7 < audit_weight < 1`. For a supplied binary matrix, `full_rank_block_packing_number` computes the maximum valid exponent `mu(H)` and `ordered_check_perfect_audit_return_bound` evaluates the resulting matrix-level endpoint bound. These are upper bounds and do not assert attainability for every matrix.

`interleaved_candidate_scores(q, v)` evaluates the original two-parameter streamed construction, and `interleaved_candidate_lower_bound` optimizes it deterministically. `interleaved_compact_lower_bound` evaluates the stronger three-effect bond-two Choi-MPS construction; at balanced weight it reaches `0.759802783851444`. Both return objects keep `support_is_globally_optimal=False`. The compact construction is physically complete after local Pauli completion and matches unrestricted variational searches, but the exact global MPS maximum remains open. Run `python scripts/verify_interleaved_counterexample.py` for the historical ternary counterexample and `python scripts/verify_compact_interleaved_candidate.py` for the current MPS construction.

`INTERLEAVED_L060_CERTIFIED_INTERVAL` exposes the published, outward-rounded
support enclosure at audit weight `0.6`. Its lower endpoint is backed by a
rational physical four-effect witness, and its upper endpoint is the maximum
of an exhaustive five-sector certificate. `lower_fraction`,
`upper_fraction`, and `width_fraction` preserve the exact declared
rationals. The object deliberately records `exact_optimum_known=False`: the
interval is rigorous within the documented proof trust boundary, but neither
endpoint is claimed to equal the unknown optimum.

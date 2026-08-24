<p align="center">
  <img src="https://raw.githubusercontent.com/asim-v/carmenq/main/assets/carmenq-logo.svg" width="560" alt="CARMEN-Q">
</p>

<p align="center">
  <strong>Exact causal audit-return bounds, statistical certificates, and reproducible quantum-memory benchmarks.</strong>
</p>

<p align="center">
  <a href="https://github.com/asim-v/carmenq/releases"><img alt="GitHub release" src="https://img.shields.io/github/v/release/asim-v/carmenq?color=102A43"></a>
  <img alt="Python 3.10+" src="https://img.shields.io/badge/python-3.10%2B-18A999">
  <a href="LICENSE"><img alt="MIT license" src="https://img.shields.io/badge/code-MIT-F2A93B"></a>
  <a href="LICENSE-CONTENT.md"><img alt="CC BY 4.0 content" src="https://img.shields.io/badge/content-CC%20BY%204.0-526777"></a>
</p>

CARMEN-Q stands for **Causal Audit-Return Memory Evaluation and Numerics for Quantum Processes**. It asks a precise operational question: can a streamed device retain a temporal predicate while preserving the ability to return every input coherently, beyond what unlimited adaptive classical memory can achieve?

The library implements the exact classical-memory frontier derived in the accompanying manuscript, a collective-access comparator, fixed-sample certification with systematic allowances, power planning, and an explicit density-matrix reference protocol. At balanced branch weight, the streaming classical ceiling is `0.75` for every stream length `n >= 2`; one persistent coherent qubit attains `1.0` in the ideal circuit.

## Install

The package is installable directly from the public release:

```bash
python -m pip install "carmenq @ git+https://github.com/asim-v/carmenq.git@lit/bounded-coherent-memory-review"
```

PyPI packaging is ready, but this repository does not claim that the `carmenq` name has already been uploaded to PyPI.
The recorded interior-frontier certificates additionally require the research
solver stack: `python -m pip install -e ".[frontier]"`.

## Use it in Python

```python
from carmenq import BenchmarkCounts, certify, collective_bound, streaming_bound

n = 8
print(streaming_bound(n, 0.5))  # 0.75
print(collective_bound(0.5))    # 0.853553...

result = certify(
    BenchmarkCounts(9700, 10000, 9500, 10000),
    n_steps=n,
    audit_weight=0.5,
    alpha=0.01,
    audit_systematic=0.005,
    return_systematic=0.005,
    null_slack=0.005,
)
print(result.certified, result.margin)
```

The same workflow is available from the terminal:

```bash
carmenq bound --steps 8
carmenq plan --steps 8 --forecast-model
python examples/quickstart.py
python examples/order_gap.py
```

The [Python API guide](docs/python-api.md) explains the concise and long-form scientific interfaces. The [benchmark specification](docs/audit-return-benchmark-v0.1.md) states the trusted access model, while the [preregistration example](docs/audit-return-preregistration.example.json) records the experimental assumptions that must be fixed before data are observed.

For frontier research, the package can also test whether every conditioned
qubit output comes from one shared quantum instrument. The inexpensive scan
checks flagged trace-norm data processing; the exact fixed-input projection
uses the optional solver stack and returns a separating witness when the
family is incompatible.

```python
from carmenq import project_to_common_instrument, scan_flagged_trace_norm_cuts

cuts = scan_flagged_trace_norm_cuts(prefix_states, conditioned_outputs)
projection = project_to_common_instrument(prefix_states, conditioned_outputs)
print(min(cut.slack for cut in cuts), projection.separation_gap)
```

The validated `lambda=0.55` strengthening experiment and its explicit local
versus global boundary are documented in
[`notes/common_instrument_strengthening_l055.md`](notes/common_instrument_strengthening_l055.md).
The latest exact theory isolates why a naive mixed-state hierarchy stalls and
replaces its concave RETURN extension by a block-coherence functional with the
same optimum as the deterministic pure-leaf problem; see
[`notes/coherence_preserving_convexification.md`](notes/coherence_preserving_convexification.md).

## Explore the order-sensitive theorem

Version 2.2 adds a robust bounded-memory result that depends on temporal order rather than only on terminal memory dimension. Two rank-two check matrices represent the same code up to a coordinate permutation. The grouped order reaches return fidelity `0.5` at perfect syndrome audit; the interleaved order has the exact and attained ceiling `0.25` under the declared one-qubit streaming interface. A causal list-decoding theorem now certifies that their support functions remain strictly separated throughout the interior interval `3/7 < audit_weight < 1`.

```python
from carmenq import (
    GROUPED_CHECK_MATRIX,
    INTERLEAVED_BALANCED_COUNTEREXAMPLE,
    INTERLEAVED_ORDER_GAP_WEIGHT_THRESHOLD,
    INTERLEAVED_PERFECT_AUDIT_ENDPOINT,
    full_crossing_perfect_audit_return_bound,
    full_rank_block_approximate_audit_return_bound,
    full_rank_block_packing_number,
    full_rank_block_perfect_audit_return_bound,
    grouped_frontier,
    interleaved_best_known_lower_bound,
    interleaved_candidate_lower_bound,
    interleaved_compact_lower_bound,
    interleaved_four_effect_lower_bound,
    interleaved_return_upper_bound,
    interleaved_support_upper_bound,
    ordered_check_perfect_audit_return_bound,
    rank_two_static_qubit_support,
    trellis_connectivity_tau,
)

print(trellis_connectivity_tau(GROUPED_CHECK_MATRIX))  # 1
print(grouped_frontier(1.0).return_fidelity)            # 0.5
print(INTERLEAVED_PERFECT_AUDIT_ENDPOINT.maximum_return_fidelity)  # 0.25
print(full_crossing_perfect_audit_return_bound(2, 2))  # 0.25
print(full_rank_block_perfect_audit_return_bound(2, 2, 3))  # 0.125
print(full_rank_block_packing_number(GROUPED_CHECK_MATRIX))  # 1
print(ordered_check_perfect_audit_return_bound(GROUPED_CHECK_MATRIX, 2))  # 0.5
print(interleaved_return_upper_bound(0.99))  # rigorous approximate-AUDIT bound
print(interleaved_support_upper_bound(0.5))  # 0.8415063509... upper certificate
print(rank_two_static_qubit_support(0.5))     # 0.8535533905... grouped value
print(INTERLEAVED_ORDER_GAP_WEIGHT_THRESHOLD)  # 3/7
print(interleaved_candidate_lower_bound(0.5).support_value)  # restricted family
print(INTERLEAVED_BALANCED_COUNTEREXAMPLE.support_value)  # 0.7594489703...
print(interleaved_compact_lower_bound(0.5).support_value)  # 0.7598027839...
print(interleaved_four_effect_lower_bound(0.6).support_value)  # 0.7658988153...
print(interleaved_best_known_lower_bound(0.6).strategy)  # four_effect_mps
```

The original two-parameter result is an exact physical construction, but it is
not the strongest known lower bound. The stored complete finite-outcome
non-QND instrument reaches `0.7594489703` at balanced weight and falsifies the
earlier frontier conjecture. A later three-effect bond-two Choi-MPS
construction reaches `0.7598027839`. Its local Pauli completion is a legal
four-outcome instrument at every slot, and the independent NumPy verifier
checks the complete tensor construction. A distinct symmetric four-effect
phase becomes stronger at larger AUDIT weights; at weight `0.6` it reaches
`0.7658988153`, exceeding the three-effect value by `0.0101928807`. This is
an exact physical lower bound. A complete exhaustion of projective, ternary,
and four-active terminal readouts now places the two-block rank-two relaxation
at that same support direction in
`[0.7658988152646944, 0.76662]`. The upper endpoint is a finite
solver-conditional numerical enclosure, not yet a solver-independent interval
proof or a claim that the four-effect construction is exactly optimal.

The arbitrary adaptive support problem now has an exact compact variational
form: maximise one normalised bond-two Choi-MPS leaf. A local Weyl/Pauli
completion proves that every feasible leaf is attainable by a complete
streamed instrument, so unbounded outcome trees are unnecessary. This closes
the outcome/completeness part of the problem. At `lambda=0.6`, the nonconvex
two-block maximum has additionally been enclosed to relative width `0.0942%`
by independent projective covers, a complete 12,008-leaf ternary SOCP cover,
and a spatial four-active bound. The entire support curve and exact equality
at that point remain open. See
[`notes/mps_leaf_completion.md`](notes/mps_leaf_completion.md) and
[`notes/interleaved_interior_frontier_l060.md`](notes/interleaved_interior_frontier_l060.md).

The endpoint theorem now extends to every finite field. If the ordered check
matrix splits into `m` consecutive blocks that each retain syndrome rank `r`,
perfect AUDIT with coherent dimension `d` across the block boundaries implies
`F_R <= min(1, (d / q**r)**m)`. Every additional full-rank temporal block adds
one factor `d/q**r`. The law is tight for repeated identity blocks when
`d=q**k`; the earlier interleaved square law is its `m=2` case. See the
[dimension-bound proof](notes/full_crossing_dimension_bound.md).
For a concrete ordered binary matrix, `full_rank_block_packing_number(H)`
computes the maximum exponent `mu(H)` greedily, and
`ordered_check_perfect_audit_return_bound(H, d)` evaluates
`min(1, (d/2**rank(H))**mu(H))`.

This exponent gives an exact asymptotic order effect. With `m` copies of each
basis column and `d=q**k`, batching equal columns has perfect-AUDIT optimum
`d/q**r`, while cycling through all basis columns `m` times has optimum
`(d/q**r)**m`. The matrices have the same column multiset; their exact return
fidelities differ by `(q**r/d)**(m-1)`.

The endpoint law also has a general approximate-AUDIT extension. Causal
dimension-to-list conversion bounds the summed Ky Fan spectral tail by
`m * D * (1 - P_A)` and turns it into an explicit RETURN certificate through
`full_rank_block_approximate_audit_return_bound`. In the canonical interleaved
instance this yields
`interleaved_support_upper_bound(lambda) < rank_two_static_qubit_support(lambda)`
for every `3/7 < lambda < 1`. This resolves the existence of a genuine
interior order effect, including balanced weight. It does not close the
remaining gap between the rigorous upper certificate and the best complete
finite-outcome lower strategy.

The [canonical result note](notes/order_sensitive_memory_result.md) states the
model, proof map, corrected numerical status, robust interval, novelty
boundary, and open problems. Run
`python scripts/classify_order_sensitive_checks.py` to reproduce the finite
classification of four-slot orders.
Run `python scripts/verify_interleaved_candidate.py` to contract the complete
two-parameter instrument independently and compare all 64 terminal vectors
with the public closed formulas.
Run `python scripts/verify_interleaved_counterexample.py` to contract the
stored 81-branch instrument, audit every local completeness relation, and
verify its strict excess over that construction without PyTorch.
Run `python scripts/verify_compact_interleaved_candidate.py` to construct the
stronger three-effect MPS, verify all temporal ranks and row-isometry
conditions, and check its local Pauli completion independently.
Run `python scripts/verify_four_effect_interleaved_candidate.py --lambda 0.6`
to perform the corresponding independent construction and checks for the new
four-effect phase.

## Reproduce the paper

Clone the repository and run:

Building the manuscript requires [Tectonic](https://tectonic-typesetting.github.io/), tested here with version 0.17.0.

```bash
python -m pip install -e ".[dev,reproducibility]"
python -m pytest -q
python scripts/run_all.py
python scripts/build_pdf.py
python scripts/build_order_pdf.py
```

The pipeline regenerates the numerical tables and publication figures under `data/` and `figures/`. Continuous integration compares cross-platform numerical outputs at declared relative and absolute tolerances of `1e-12`, while same-platform tests retain byte-level determinism. The visually verified manuscript is available as [PDF](output/pdf/CARMEN-Q-paper.pdf), with LaTeX sources in `manuscript/` and audited references in `references/`. The package source lives in `src/carmenq/`; `tests/` contains the analytic, statistical, protocol, and reproducibility checks.

The focused temporal-order result has its own visually verified
[PDF](output/pdf/CARMEN-Q-order-paper.pdf) and two-column LaTeX source in
`manuscript-order/`. It proves the linear rank-tail theorem, the certified
interior order gap, and the exact homogeneous Choi-MPS reduction. The compact
three-effect curve in that manuscript is a verified physical lower bound. The
accompanying research notes and machine-readable certificates additionally
verify the stronger four-effect lower bound and enclose, at `lambda=0.6`, the
streamed support and its two-block relaxation between
`0.7658988152646944` and `0.76662`. This newer enclosure is not yet part of
the PDF; exact equality and the complete support curve remain explicitly open.

## Scientific boundary

CARMEN-Q is a trusted-interface resource witness. A positive score rejects the declared adaptive classical-memory model under the stated source, timing, sequestration, and measurement assumptions. It does not identify a unique microscopic implementation, prove deletion of inaccessible environmental records, or select an interpretation of quantum mechanics.

The project began as an investigation of reversible quantum histories. Its originality audit showed that compute-phase-uncompute circuits alone are established quantum computing. CARMEN-Q retains the part that survived that audit: an exact same-task separation between streamed classical memory, collective classical recording, and coherent temporal memory.

The [focused literature review](notes/literature_review_bounded_coherent_memory.md) records the originality gate. Its conclusion is deliberately narrow: bounded coherent memory, late choice, syndrome accumulation, entanglement recovery, and trellis order dependence are prior art. What survives is the exact operational consequence of permuting one code under a fixed coherent-memory constraint, the perfect-AUDIT temporal product law for consecutive full-rank blocks, and its linear-tail approximate-AUDIT extension. A broad interval of interior order dependence is proved, and one interior support direction is now completely enclosed numerically; the exact support curve and the optimal robustness coefficient remain open.

## Citation and authorship

The manuscript and software are by **Javier Emilio Bazán Sánchez**, Facultad de Ciencias, Universidad Nacional Autónoma de México, `bazan@ciencias.unam.mx`. Machine-readable citation metadata are in [`CITATION.cff`](CITATION.cff).

Code is distributed under the MIT License. Unless noted otherwise, the manuscript, figures, and data are distributed under CC BY 4.0.

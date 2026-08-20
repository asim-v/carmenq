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
```

The [Python API guide](docs/python-api.md) explains the concise and long-form scientific interfaces. The [benchmark specification](docs/audit-return-benchmark-v0.1.md) states the trusted access model, while the [preregistration example](docs/audit-return-preregistration.example.json) records the experimental assumptions that must be fixed before data are observed.

## Explore the order-sensitive theorem

Version 2.1 adds a bounded-memory result that depends on temporal order rather than only on terminal memory dimension. Two rank-two check matrices represent the same code up to a coordinate permutation. The grouped order reaches return fidelity `0.5` at perfect syndrome audit; the interleaved order has the exact and attained ceiling `0.25` under the declared one-qubit streaming interface.

```python
from carmenq import (
    GROUPED_CHECK_MATRIX,
    INTERLEAVED_PERFECT_AUDIT_ENDPOINT,
    grouped_frontier,
    trellis_connectivity_tau,
)

print(trellis_connectivity_tau(GROUPED_CHECK_MATRIX))  # 1
print(grouped_frontier(1.0).return_fidelity)            # 0.5
print(INTERLEAVED_PERFECT_AUDIT_ENDPOINT.maximum_return_fidelity)  # 0.25
```

The [canonical result note](notes/order_sensitive_memory_result.md) states the model, proof map, robust interval, novelty boundary, and open problems. Run `python scripts/classify_order_sensitive_checks.py` to reproduce the finite classification of four-slot orders.

## Reproduce the paper

Clone the repository and run:

Building the manuscript requires [Tectonic](https://tectonic-typesetting.github.io/), tested here with version 0.17.0.

```bash
python -m pip install -e ".[dev,reproducibility]"
python -m pytest -q
python scripts/run_all.py
python scripts/build_pdf.py
```

The pipeline regenerates the numerical tables and publication figures under `data/` and `figures/`. Continuous integration compares cross-platform numerical outputs at declared relative and absolute tolerances of `1e-12`, while same-platform tests retain byte-level determinism. The visually verified manuscript is available as [PDF](output/pdf/CARMEN-Q-paper.pdf), with LaTeX sources in `manuscript/` and audited references in `references/`. The package source lives in `src/carmenq/`; `tests/` contains the analytic, statistical, protocol, and reproducibility checks.

## Scientific boundary

CARMEN-Q is a trusted-interface resource witness. A positive score rejects the declared adaptive classical-memory model under the stated source, timing, sequestration, and measurement assumptions. It does not identify a unique microscopic implementation, prove deletion of inaccessible environmental records, or select an interpretation of quantum mechanics.

The project began as an investigation of reversible quantum histories. Its originality audit showed that compute-phase-uncompute circuits alone are established quantum computing. CARMEN-Q retains the part that survived that audit: an exact same-task separation between streamed classical memory, collective classical recording, and coherent temporal memory.

The [focused literature review](notes/literature_review_bounded_coherent_memory.md) records the originality gate. Its conclusion is deliberately narrow: bounded coherent memory, late choice, syndrome accumulation, entanglement recovery, and trellis order dependence are prior art. What survives is the exact operational consequence of permuting one code under a fixed coherent-memory constraint. The full interleaved interior frontier remains open.

## Citation and authorship

The manuscript and software are by **Javier Emilio Bazán Sánchez**, Facultad de Ciencias, Universidad Nacional Autónoma de México, `bazan@ciencias.unam.mx`. Machine-readable citation metadata are in [`CITATION.cff`](CITATION.cff).

Code is distributed under the MIT License. Unless noted otherwise, the manuscript, figures, and data are distributed under CC BY 4.0.

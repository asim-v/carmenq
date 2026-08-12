# Causal Audit--Return Quantum Memory

**An exact classical-memory frontier and a coherent one-qubit separation**

This repository studies one operational question:

> Can a streamed quantum device retain a temporal predicate while preserving
> the ability to return every input coherently, beyond what unlimited adaptive
> classical memory can achieve?

The project is formulated within standard quantum mechanics. It does not attempt to prove a preferred interpretation, communicate with macroscopic branches, or attribute consciousness to small quantum registers.

[Read the current manuscript (PDF)](output/pdf/causal_audit_return_memory.pdf)

## Version 1.1 deliverables

Version `v1.1` is a reproducible theory-and-benchmark research artifact containing:

- the exact adaptive classical-memory streaming support;
- a strict collective classical-record comparator;
- a one-qubit coherent strategy attaining perfect AUDIT and RETURN;
- a proof of the first-order optimal-strategy transition for `n >= 3`;
- conservative finite-sample certification and power planning;
- a hardware-facing protocol and preregistration template;
- adversarial numerical validation, deterministic data, figures, and tests;
- a 19-page English LaTeX manuscript and visually verified final PDF.

## The central limitation

The closed circuit

```text
prepare -> U_history -> predicate phase -> U_history^dagger -> interfere
```

is an established instance of reversible computation and phase kickback. Viewed only at its endpoints, it is equivalent to a directly compiled effective unitary. A candidate contribution begins only when the protocol declares testable access, locality, time, or resource restrictions and requires multitime evidence of causal memory. This is a kill criterion, not a footnote.

## Current development: causal audit--return benchmark

The originality audit led to a narrower and experimentally useful task. Fresh
entangled carriers pass through a sequential device and are immediately
sequestered. Only after the stream is committed, the verifier asks the device
either to predict a temporal parity from its committed terminal memory
(**AUDIT**) or to restore every entangled pair and visibly reset that memory
(**RETURN**).

For an adaptive device with unlimited classical memory but no quantum state
persisting between slots, the exact weighted null is implemented in the
package. At the balanced operating point it is `0.75` for every stream length
`n >= 2`; a coherent parity-memory qubit attains `1`. The repository now
includes fixed-sample certification, systematic-error allowances, conservative
power planning, deterministic frontier/forecast figures, and a preregistration
template. This is a trusted-interface quantum-memory witness--not a claim about
Everett or consciousness.

```bash
audit-return-benchmark bound --steps 8
audit-return-benchmark plan --steps 8 --forecast-model
```

See the [benchmark specification](docs/audit-return-benchmark-v0.1.md) and the
[example preregistration](docs/audit-return-preregistration.example.json).

## Reproduce the results

Python 3.10 or later is required.

```bash
python -m pip install -e ".[dev,reproducibility]"
python -m pytest
python scripts/run_all.py
```

The `reproducibility` extra pins the numerical and plotting versions used for
the committed release artifacts. For normal library use, `python -m pip
install -e .` keeps the broader supported dependency ranges.

Generated data are written to `data/` and figures to `figures/`. Every output records its parameters and deterministic seed.

Compile the paper with Tectonic (available on `PATH`, through the `TECTONIC`
environment variable, or with `--tectonic PATH`):

```bash
python scripts/build_pdf.py
```

## Repository map

- `manuscript/`: self-contained LaTeX preprint.
- `output/pdf/`: visually verified final PDF.
- `src/`: exact frontier, statistical planner, and reference simulators.
- `scripts/`: data, figure, and manuscript reproduction commands.
- `tests/`: algebraic and numerical checks.
- `figures/`: reproducible paper figures in vector and raster formats.
- `data/`: numerical results, reduced states, and metadata.
- `docs/`: versioned protocol specification.
- `notes/`: formal theory, simulation design, and the dated novelty audit.
- `references/`: audited BibTeX library.
- `CLAIMS.md`: claim ledger with evidence status.
- `CHANGELOG.md`: versioned release history.

## Epistemic firewall

A positive result rejects the declared adaptive classical-memory null under the
trusted source, timing, sequestration, and measurement assumptions. It does
not identify a unique implementation or establish an interpretation of quantum
mechanics. A negative result may reflect ordinary noise, leakage, drift, or
insufficient power.

## Licence and citation

Code is released under the MIT License. Unless noted otherwise, the manuscript,
data, and figures are released under CC BY 4.0; see `LICENSE-CONTENT.md`.
Citation metadata are in `CITATION.cff`.

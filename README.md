# Zero-Record Quantum Histories

**Causal and interferometric certification of reversible quantum histories**

This repository studies one operational question:

> Can a quantum process be certified to have executed a nontrivial causal history, retain only a selected predicate about it, and finish without a physically distinguishable transcript of the complete history?

The project is formulated within standard quantum mechanics. It does not attempt to prove a preferred interpretation, communicate with macroscopic branches, or attribute consciousness to small quantum registers.

[Read the final manuscript (PDF)](output/pdf/zero_record_quantum_histories.pdf)

## Version 1 deliverables

Version `v1.0` is a reproducible theory-and-simulation research artifact containing:

- an operational definition of **zero-record certification**;
- formal results about phase kickback, no-transcription, and visibility;
- an explicit separation between a trivially compilable closed circuit and a causal task with restricted temporal access;
- ideal and noisy simulations of reversible histories;
- incoherent, retained-memory, absent-predicate, and direct-phase controls;
- a multistep reversible transducer with post-memory challenge parameters,
  plus a separate analytic late-challenge task;
- automated tests, generated data, and publication figures;
- a LaTeX manuscript and visually verified final PDF.

## The central limitation

The closed circuit

```text
prepare -> U_history -> predicate phase -> U_history^dagger -> interfere
```

is an established instance of reversible computation and phase kickback. Viewed only at its endpoints, it is equivalent to a directly compiled effective unitary. A candidate contribution begins only when the protocol declares testable access, locality, time, or resource restrictions and requires multitime evidence of causal memory. This is a kill criterion, not a footnote.

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
- `src/`: density-matrix reference simulator.
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

A positive result establishes interference and coherent control under the declared model and restrictions. It does not establish a many-worlds ontology. An anomalously negative result does not establish objective collapse before noise, drift, leakage, preparation, measurement, and imperfect reversal have been excluded.

## Licence and citation

Code is released under the MIT License. Unless noted otherwise, the manuscript,
data, and figures are released under CC BY 4.0; see `LICENSE-CONTENT.md`.
Citation metadata are in `CITATION.cff`.

# Version 1 specification: zero-record certification

## 1. Scientific object

The basic object is not an ontological "branch" but a multitime quantum process with explicitly declared physical records. The joint objective is to certify:

1. that an internal memory causally mediates between early interventions and late responses;
2. that an allowed property of the process changes a final phase or statistic;
3. that distinguishable history records are reset or decoupled;
4. that interferometric coherence is recovered; and
5. that the data exclude incoherent controls within the experimental model.

## 2. Systems and access

This specification uses

\[
\mathcal H=\mathcal H_B\otimes\mathcal H_W\otimes\mathcal H_M
\otimes\mathcal H_G\otimes\mathcal H_A\otimes\mathcal H_E,
\]

where `B` labels alternatives, `W` is the controlled world, `M` is memory, `G` is reversible garbage, `A` is the readout register, and `E` contains uncontrolled degrees of freedom.

A history labelled by \(b\) is a sequence of channels or isometries

\[
\mathbf U_b=(U_{b,T},\ldots,U_{b,1})
\]

interleaved with intervention ports. The full process can be represented as a comb or process tensor. Without ports or access restrictions, its endpoint action cannot certify its internal structure.

## 3. Models that must not be conflated

### 3.1 Known closed-circuit model

The experimenter knows and controls all of \(U_H\). For a clean predicate \(p\),

\[
U_H^\dagger O_p U_H |b\rangle|0\rangle
=(-1)^{p_b}|b\rangle|0\rangle.
\]

This model validates physical control, noise sensitivity, and scaling. It does not show that an internal trajectory was necessary: the effective unitary can be compiled directly.

### 3.2 Causal-certification model

The verifier has only declared ports, times, and admissible operations. At least one challenge \(x_t\) is chosen after its proposed causal parent exists. Statistics are compared against alternative classes: no memory, bounded memory dimension, nonadaptive policy, or direct phase insertion without authorised access.

Every conclusion is relative to that comparison class and those access assumptions.

## 4. Operational definition

A protocol family is a **zero-record certification** with parameters

\[
(c,s,\varepsilon,\eta,v;\mathfrak A,\mathfrak C)
\]

for valid processes \(\mathfrak A\) and shortcuts \(\mathfrak C\) if:

- **completeness:** a valid process accepts with probability at least \(c\);
- **soundness:** every process in \(\mathfrak C\) accepts with probability at most \(s<c\);
- **zero record:** for a declared non-injective predicate and history distribution, all accessible residual transcript systems \(T_f\) obey

  \[
  I(H:T_f\mid p(H))\leq\varepsilon;
  \]

- **reset:** \(F(\rho_{WMG}^{f},|w_0,0,0\rangle)\geq1-\eta\); and
- **recoherence:** normalised coherence or visibility obeys \(V\geq v\).

This is not automatically a cryptographic zero-knowledge proof. That term requires a simulator and a complete adversarial model. Here, zero record is a physical final-decoupling property.

## 5. Reference instances

### ZR-0: phase survival

Two alternatives, one memory bit, one binary predicate, and exact uncomputation. This is an algebraic test and positive control, not a novelty claim.

### ZR-1: partial record

The environment retains states \(|e_0\rangle,|e_1\rangle\) with overlap \(\gamma\). For balanced amplitudes and pure records,

\[
V=|\gamma|,\qquad D=\sqrt{1-|\gamma|^2}.
\]

This instance validates leakage, visibility, and transcript-residue calculations.

### ZR-2: multistep reversible agent

The implemented reversible transducer stores a two-bit world symbol, computes challenge-parameter-dependent decisions into a work bit, conditionally changes the world, and uncomputes the work bit. The final internal predicate depends on those actions. Every transition is a reversible permutation and is inverted after phase kickback.

### ZR-3: late challenge

In the analytic task, an independent challenge \(x\) is selected after an intermediate memory is created. The valid response depends on both \(x\) and that memory. The challenge may remain at the end if it is independent of the history label; history-dependent responses and memories must be uncomputed.

The present simulator has configured challenge parameters in the post-memory gate sequence, but no independent \(X\) register and no ZR-3 game score. The random-access bound applies to a classical readout game, not directly to the coherent phase endpoint. Excluding a precompiled phase requires an additional inequality or task and is not claimed in version 1.

### ZR-4: four-history coarse-graining

Four fine-grained labels share a non-injective predicate such as parity. This makes the conditional leakage \(I(H:R_f\mid p(H))\) nonvacuous. Exact reset gives zero leakage, whereas retained or partial within-class records give positive leakage.

## 6. Mandatory metrics

- control visibility \(V\);
- \(l_1\) coherence or selected off-diagonal elements;
- reset fidelity \(F_{\mathrm{reset}}\);
- predicate fidelity \(F_p\);
- residual mutual information \(I(H:T_f\mid p(H))\), with the tested transcript system and prior stated explicitly;
- depth, gates, ancillas, and history steps;
- parameters of every noise channel; and
- numerical error or confidence intervals when applicable.

## 7. Mandatory controls

Every empirical claim must be compared with:

1. ideal coherent evolution;
2. a classical mixture of labels;
3. correlated memory that is not uncomputed;
4. an unrecovered history-dependent environment;
5. perturbed reversal;
6. an absent predicate;
7. a directly compiled phase; and
8. a memoryless or nonadaptive policy.

The direct-phase control is expected to match the known closed model. Equality is the endpoint non-certifiability diagnosis.

## 8. Version 1 success criteria

Version 1 succeeds if it:

- reproduces ZR-0 exactly;
- recovers \(V=|\gamma|\) for ZR-1 within numerical tolerance;
- demonstrates memory participation and reset in ZR-2, and states the ZR-3 classical bound analytically;
- shows monotonic visibility loss under the declared dephasing and leakage models;
- makes direct-phase endpoint equivalence explicit; and
- does not present an unproved causal or cryptographic advantage as a theorem.

## 9. Kill criteria

The programme does not justify a new contribution if, after an access model is imposed, all proposed tasks:

- reduce to a known Hadamard test or phase oracle;
- fail to exclude memoryless models;
- yield no new bound, resource separation, or operational primitive;
- depend only on renaming ancillas as observers; or
- require interpretive language to produce distinct predictions.

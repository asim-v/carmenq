# Simulation design and interpretation

## Scope

This repository contains a small, explicit density-matrix reference model for
reversible quantum-history interferometry.  It demonstrates phase kickback,
history-dependent reversible computation, uncomputation, controlled loss of
which-history information, and final predicate readout.  It is not evidence
for a particular interpretation of quantum mechanics and it does not certify
that a physical device executed an irreducible internal history.  In
particular, the included direct-phase control is endpoint-equivalent to the
ideal trusted circuit.

The implementation uses only NumPy, SciPy, and Matplotlib.  No quantum SDK or
opaque circuit optimizer is required.

## Registers

The complete Hilbert space is

\[
\mathcal H = \mathcal H_B\otimes\mathcal H_W\otimes\mathcal H_M
\otimes\mathcal H_G\otimes\mathcal H_A,
\]

with dimensions

\[
4\times4\times4\times2\times2=256.
\]

- `B`: four history labels, encoded as two logical bits;
- `W`: a two-bit internal world state;
- `M`: a two-bit reversible memory;
- `G`: one decision/work bit;
- `A`: a phase-kickback ancilla initialized in \(|-\rangle\).

Four histories are used rather than two.  The allowed predicate is two-bit
parity,

\[
p(b)=\operatorname{popcount}(b)\bmod 2,
\qquad p=(0,1,1,0).
\]

Each predicate class therefore contains two labels.  A residual record can
reveal one full bit about `H` even after parity is given, so
\(\chi(H:R\mid P)\) is a nontrivial diagnostic.  With only two labels and a
predicate that identifies the label, the same conditional quantity would be
zero by construction.

## Reversible history and late challenge

The initial coherent state is

\[
|\Psi_0\rangle=
\frac{1}{2}\sum_{b=0}^{3}|b\rangle_B|0\rangle_W|0\rangle_M
|0\rangle_G|-\rangle_A.
\]

The forward history first prepares `W += B (mod 4)` and copies the world into
memory by `M += W (mod 4)`.  After this memory exists in the gate sequence,
the simulator applies a configured challenge parameter \(c_r\) in round \(r\).
These parameters are known when the circuit is constructed: the simulator does
not implement an independent late source, hide earlier access to \(c_r\), or
score the analytic ZR-3 random-access game.  The reversible transducer computes

\[
d_r=m_{r\bmod2}\oplus c_r
\]

into `G`, conditionally flips one world bit, and uncomputes `G`.  Repetition
produces a multistage, challenge-responsive internal trajectory without
leaving work garbage.  The predicate oracle reads final internal state and
uses

\[
p_{\mathrm{internal}}
=\operatorname{parity}(W)\oplus\bigoplus_r d_r
=\operatorname{parity}(B)
\]

on the intended causal path.  Controlled `X` on the ancilla implements

\[
|x\rangle|-\rangle\longmapsto(-1)^{p_{\mathrm{internal}}(x)}
|x\rangle|-\rangle.
\]

The agent, memory, and world gates are then inverted.  Applying
\(H\otimes H\) to `B` maps the parity phase pattern
\((1,-1,-1,1)/2\) to \(|11\rangle\), so the ideal predicate readout is
deterministic.

The `action_bypassed` control keeps the oracle definition fixed while removing
the internal actions.  Its failure shows that the simulated oracle is wired to
the declared internal trajectory.  This is a white-box causal check, not a
device-independent certificate.  The `direct_phase` control applies parity
directly to `B` and reaches the same endpoint; excluding that shortcut in a
future experiment requires access constraints, randomized interventions, an
unknown process, or a cryptographic soundness model.

## Reported metrics

The reduced state of the history label immediately before readout is
\(\rho_B\).  The normalized \(l_1\) coherence is

\[
V_{l_1}=\frac{1}{K-1}\sum_{b\ne b'}|\rho_{bb'}|,
\qquad K=4.
\]

It equals one for a maximally coherent branch state and zero for a branch
mixture.  It is a multipath coherence measure; it should not be confused with
every possible operational definition of two-path fringe visibility.

Reset fidelity is the final population of the target subspace

\[
F_{\mathrm{reset}}=
\Pr[W=0,M=0,G=0],
\]

with `B` and `A` unrestricted.  Predicate fidelity is
\(F_p=\Pr[B=3]\) after the two-qubit Hadamard readout.  The normalized target
contrast is

\[
C_p=\frac{F_p-1/4}{1-1/4}.
\]

For branch-conditioned final record states
\(\rho_{WMG}^{(b)}\), residual transcript information is reported as

\[
R_{\mathrm{transcript}}
=\chi(H:WMG\mid P).
\]

This conditional Holevo quantity is an upper bound on accessible classical
information in the modeled record ensemble.  It is zero after ideal
uncomputation and one bit when `M` retains the complete four-valued label.

The benchmark table also reports the sum of the separately modeled internal-
record and environment-record conditional Holevo quantities.  That column is
only a bookkeeping sum of marginals.  It is not asserted to bound the Holevo
information of an arbitrary correlated joint record, which must be evaluated
from its joint states.

## Environment leakage

Uncontrolled environment records are modeled by pure states satisfying

\[
\langle e_{b'}|e_b\rangle=\eta\quad(b\ne b'),
\qquad 0\le\eta\le1.
\]

After tracing out the environment, each off-diagonal branch block is multiplied
by \(\eta\).  In the otherwise ideal protocol this gives exactly

\[
V_{l_1}=\eta,
\qquad F_p=\frac14+\frac34\eta.
\]

The complementary pure-record ensemble has
\(\chi(H:E\mid P)=h_2((1+\eta)/2)\): zero at \(\eta=1\) and one bit at
\(\eta=0\).  This constant-overlap model is analytically transparent, but it
does not represent every structured environment.

## Noise and inversion models

After each addressed logical operation, the exact density matrix can undergo:

- qudit dephasing, multiplying local off-diagonal terms by \(1-p\);
- qudit depolarization,
  \((1-p)\rho+p(I/d)\otimes\operatorname{Tr}_{\mathrm{target}}\rho\);
- qudit relaxation, taking every excited level toward level zero with
  probability \(p\).

These are comparative Markovian models, not platform calibrations.  In
particular, the four-level relaxation channel is not identical to independent
amplitude damping on two hardware qubits.

Imperfect inversion is modeled as a stochastic skipped inverse operation,

\[
\mathcal E_{\mathrm{inv}}(\rho)
=(1-\epsilon)U^\dagger\rho U+\epsilon\rho.
\]

This deliberately simple model separates inverse reliability from the three
local channels.  Coherent over-rotation, leakage outside the logical Hilbert
space, correlated noise, non-Markovian environments, measurement error, and
error mitigation are not modeled.

## Reproducibility

Run from the repository root:

```powershell
python scripts/regenerate.py --seed 20260812 --shots 8192
python -m pytest
```

The exact density-matrix results are deterministic.  The seed controls only
finite-shot multinomial samples and the challenge strings in the depth sweep.
Generated CSV tables, reduced reference states, version metadata, and both PDF
and PNG figures are written to `data/` and `figures/`.

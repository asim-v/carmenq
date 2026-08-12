# Originality gate: frozen scientific core

Cutoff: 2026-08-12

Status: **GO for a cautious theorem-level research claim; NO-GO for the original
reversible-history framing as the novelty claim.** This is a research-priority
assessment, not proof that no equivalent unpublished or differently named
result exists.

## 1. Frozen contribution

Working technical description:

> **Causal audit--return duality for a reversible streaming predicate.**

The verifier sends `n` fresh halves of EPR pairs through a sequential device,
one at a time. Every returned carrier is immediately sequestered. Only after
the complete prefix has been committed does the verifier choose between:

1. **AUDIT:** measure the reference halves in `Z` and ask the committed
   classical transcript for the parity of the resulting temporal bit string;
2. **RETURN:** use the transcript to control an arbitrary joint recovery and
   test whether every carrier is returned maximally entangled with its
   reference.

The null class has arbitrary within-slot instruments, finite outcomes,
disposable local quantum ancillas, classical randomness, unlimited adaptive
classical memory, and transcript-conditioned joint recovery. It has no quantum
state or pre-shared entanglement persisting between slots and no later access
to sequestered outputs.

For

\[
 f(t)=\frac{1+\sqrt{1-t^2}}{2},
\]

the exact support function of that full classical-memory class is

\[
 \boxed{
 \sup\{\lambda P_A+(1-\lambda)F_R\}
 =\max_{0\le t\le1}
 \left[
 \frac{\lambda}{2}(1+t^n)+(1-\lambda)f(t)^n
 \right].
 }
\]

The bound covers adaptive and non-QND local instruments. It is attained by
independent binary symmetric weak Lueders-type instruments. The proof has four
nontrivial parts:

1. an arbitrary-instrument recovery lemma reducing each transcript leaf to
   its basis likelihoods;
2. causal factorization of those likelihoods in the absence of coherent
   inter-slot memory;
3. multiplication of conditional parity bias and return fidelity along each
   classical path, making adaptation unable to improve a linear support;
4. a convexity argument proving that equal local strengths are globally
   optimal, rather than merely a symmetric ansatz.

## 2. New-looking consequences

For `n=2`, the optimum changes continuously at `lambda=1/2`. For every
`n>=3`, it undergoes a first-order transition. Let `z_n` be the unique root

\[
 (1-z_n)(1+z_n)^{n-1}=1,
 \qquad z_n\in((n-2)/n,1).
\]

Then

\[
 \lambda_c(n)=
 \frac{1}{1+2^{n-2}z_n^{(n-2)/2}(1-z_n)}.
\]

At the transition the optimal classical strategy jumps from recording
nothing to an almost projective measurement. Numerically,

| `n` | `lambda_c` | measurement-strength jump |
|---:|---:|---:|
| 3 | 0.6247789017 | 0.9717365435 |
| 4 | 0.6495455168 | 0.9961752259 |
| 5 | 0.6588954571 | 0.9992936175 |

Moreover `lambda_c(n)` approaches `2/3` from below.

A collective instrument that processes all carriers at once but retains only
a classical outcome has the exact one-effective-bit support

\[
 \sup S_{\rm coll}^{\rm classical}
 =\frac12+\frac12\sqrt{\lambda^2+(1-\lambda)^2},
\]

strictly above the streaming classical-memory value for `n>=2` and interior
`lambda`. This comparator is **not** an unrestricted coherent device.

One persistent coherent qubit closes the gap completely. Sequential CNOTs
accumulate parity coherently; AUDIT measures the accumulator, while RETURN
uses the sequestered carriers to uncompute it. This strategy attains

\[
 (P_A,F_R)=(1,1).
\]

Thus the result is an exact same-task separation among causal classical
memory, collective classical recording, and one qubit of coherent temporal
memory.

## 3. What is not claimed as new

The following ingredients are occupied and must be presented as antecedents:

- one-shot information--disturbance and weak-measurement curves;
- parity accumulation with an ancilla;
- multiplicative bias for estimating parity from independent records;
- collective advantages for parity discrimination;
- read-versus-return quantum seals;
- late Ways-versus-Phases choices;
- process tensors, quantum-memory witnesses, classical-memory hierarchies,
  sequential QRACs, and bounded-memory testers;
- EPR/process-fidelity tests, uncomputation, erasure, and reversible computing.

The old headline "histories can be queried and recombined" is therefore not a
defensible originality claim. The candidate contribution is only the precise
streaming parity audit--global-return frontier, its causal tensorization,
first-order phase diagram, and strict coherent-memory separation.

## 4. Priority verdict

An adversarial search across quantum information--disturbance, quantum seals,
parity cryptography, process combs, temporal-memory witnesses, sequential
QRACs, bounded-memory testers, and weak parity measurement found no primary
source stating the full theorem above under the same access model.

The most dangerous parity collision is Bennett--Mor--Smolin (1996): it proves
that collective measurements can extract more information about the parity of
product nonorthogonal states than separate measurements and contains the same
multiplicative-bias motif. It does not pose the late AUDIT/RETURN alternative,
score joint EPR return, characterize the causally factorized
zero-coherent-memory instrument class, or derive this support and transition.

The most dangerous architectural collisions are quantum seals (read or return
for tamper detection), multipath Ways/Phases games, the formal hierarchy of
multi-time classical memory, quantum-memory witnesses, sequential QRAC
witnesses, and constrained-separability descriptions of classically adaptive
testers. They occupy the ingredients and vocabulary, but the exact conjunction
and formula were not located by the cutoff.

Closest primary sources checked include:

- Bennett, Mor, and Smolin, *Physical Review A* 54, 2675 (1996),
  <https://doi.org/10.1103/PhysRevA.54.2675>;
- Banaszek, *Physical Review Letters* 86, 1366 (2001),
  <https://doi.org/10.1103/PhysRevLett.86.1366>;
- Bisio *et al.*, *Physical Review A* 85, 032333 (2012),
  <https://doi.org/10.1103/PhysRevA.85.032333>;
- Bagan *et al.*, *Physical Review Letters* 120, 050402 (2018),
  <https://doi.org/10.1103/PhysRevLett.120.050402>;
- Kimmel and Kolkowitz, *Physical Review A* 100, 052326 (2019),
  <https://doi.org/10.1103/PhysRevA.100.052326>;
- Giarmatzi and Costa, *Quantum* 5, 440 (2021),
  <https://doi.org/10.22331/q-2021-04-26-440>;
- Taranto *et al.*, *Quantum* 8, 1328 (2024),
  <https://doi.org/10.22331/q-2024-05-02-1328>;
- Roy *et al.*, *Physical Review A* 110, 012608 (2024),
  <https://doi.org/10.1103/PhysRevA.110.012608>;
- Vieira, Ku, and Budroni, *Physical Review Research* 7, 043281
  (2025), <https://doi.org/10.1103/r8lf-bb4p>;
- Ohst *et al.*, *Quantum* 10, 1988 (2026),
  <https://doi.org/10.22331/q-2026-01-28-1988>.

Safe priority language:

> To our knowledge, for the precisely specified streaming class, this work
> gives the exact parity-audit versus global-entanglement-return support,
> including arbitrary adaptive non-QND instruments, and proves a strict
> same-task separation from collective classical-record strategies and from a
> one-qubit coherent accumulator. The claimed contribution is the causal
> tensorization and phase diagram, not weak measurement, parity processing,
> delayed choice, or quantum-memory witnessing in general.

## 5. Validation record

The proof and scope audit are in
`notes/originality_gate_parity_theory.md`. Expensive adversarial searches are
implemented in `scratch/originality_gate/streaming_parity.py` with stored output
in `streaming_parity.json`. They found:

- maximum adaptive-tree excess over the theorem: `6.52e-14`;
- maximum asymmetric-product excess: `3.53e-14`;
- positive collective-classical gap in every tested interior direction;
- no numerical counterexample for `n=2,...,5`.

The independent fast validator
`scratch/originality_gate/validate_streaming_parity.py` additionally checks:

- one- and two-slot closed forms;
- the transition equations for `n=3,...,12`;
- random product-recovery inequalities;
- random collective likelihood tables against the classical-record bound;
- the exact coherent-accumulator circuit;
- the stored asymmetric and adaptive searches.

Its frozen result is `streaming_parity_validation.json`: `PASS`.

## 6. Hard scope conditions and kill gates

The theorem ceases to apply if a device can retain an uncharged coherent system
between slots, share entanglement among slot ancillas, revisit earlier outputs,
use the returned carriers in AUDIT, exploit correlated verifier inputs, or
postselect recovery. "Classical memory" must mean the explicit classical-path
factorization, not merely that an intermediate channel is
entanglement-breaking.

Visible transcript reset is not global erasure. If RETURN also requires reset
of every purification of a classical record, the null class changes. This must
be specified rather than hidden in the word "erase."

Priority is killed if expert review produces an earlier theorem with all of:

1. fresh sequential entangled inputs and immediate output sequestration;
2. arbitrary adaptive instruments with no coherent inter-slot state;
3. late transcript-only global-parity prediction versus joint EPR return;
4. the same exact support over arbitrary non-QND instruments;
5. the `n>=3` first-order transition or an immediately specializing general
   theorem.

## 7. Decision

The originality gate is passed at a cautious theorem level. The scientific
project should now be rebuilt around **causal audit--return duality for
reversible streaming predicates**, with the parity theorem as its first exact
result. Everettian language, observers, alternate worlds, and generic
uncomputation should remain motivation or interpretation at most, never the
novelty claim.

# Originality gate: numerical audit of delayed-choice reveal-or-erase games

Date: 2026-08-12
Status: adversarial scratch audit; not a release claim
Code: `scratch/originality_gate/audit_delayed_choice.py`
Deterministic output: `scratch/originality_gate/results.json`
Seed: `20260812`

## Executive verdict

The most natural delayed-choice game is **unsound**. If its erase branch only
checks that a known coherent state, such as \(|+_d\rangle\), reappears, an
entanglement-breaking device can score perfectly: measure the history label,
store it classically, answer every late reveal query, and in the erase branch
discard everything and prepare a fresh \(|+_d\rangle\). Thus

\[
P_{\rm reveal}=1,\qquad V_{\rm erase}=1
\]

without a coherent reversible prefix. Auditing a named memory reset does not
repair this attack if an unobserved environment may receive the discarded
record.

A verifier-held entangled reference repairs the fresh-state attack. If erase is
scored by entanglement fidelity, every entanglement-breaking prefix and all of
its subsequent processing obey

\[
F_{\rm erase}\leq \frac1d.
\]

However, that repaired test certifies only that quantum information survived
somewhere. A coherent SWAP of the input into a private memory scores perfectly
on both reveal and erase without making a second record or implementing the
advertised world-to-memory history. If the reveal responder may access the
emitted world register, even the identity channel is a perfect memoryless
strategy. Therefore a delayed choice is not a causal restriction by itself;
the locations and post-prefix access rights are indispensable parts of the
game.

The most defensible next candidate is a **two-location coherent
broadcast-and-erase game**: after the common prefix, the world and memory ports
are separated and queried simultaneously in reveal trials, while erase trials
reunite them and test entanglement recovery. This blocks both the identity and
SWAP shortcuts under enforced noncommunication. It operationally tests a
reversible redundant record of one classical observable. Its ingredients are,
however, very close to coherent fanout, quantum nondemolition measurement,
which-path complementarity, and quantum erasure. It is a promising benchmark,
not yet a demonstrated new primitive.

## 1. Access model under audit

The common-prefix idea has four times:

1. The verifier prepares a system \(B\), possibly entangled with an inaccessible
   reference \(R\).
2. A device applies one fixed prefix before learning the branch choice. It emits
   a world port \(Q\) and may retain a memory port \(M\).
3. The verifier selects `reveal` or `erase` after the prefix.
4. In a reveal trial it asks for a past fact. In an erase trial it requests
   reversal and scores recovered coherence.

This description is incomplete until it says:

- whether the reveal responder can access \(Q\);
- whether \(Q\) is returned in the erase branch;
- whether the verifier audits a reference system or only a known state;
- whether private environments and fresh ancillas are allowed;
- whether aborts and losses count as failures;
- which ports may communicate after the choice;
- whether a reset is checked globally or only on named registers; and
- whether the recovery is verifier-controlled, fixed in advance, or supplied by
  the untrusted device.

Every positive statement below is relative to its stated answers.

## 2. Falsifying strategies

### 2.1 Fresh-state attack against a known-state erase test

Let the common input be \(|+_d\rangle=d^{-1/2}\sum_h|h\rangle\). An
entanglement-breaking prefix measures in the \(h\) basis and stores the outcome.
On reveal it returns any requested function of \(h\). On erase it traces out
the record and prepares a new \(|+_d\rangle\). Both tests accept with
probability one.

This remains possible if the verifier checks that a named memory is \(|0\rangle\)
at the end. A CPTP map can reset that memory while exporting its entropy to an
unobserved environment. Requiring the device to call its operation an inverse
does not impose an operational constraint.

**Consequence:** ordinary visibility of one known input cannot distinguish
recoherence from replacement. Random phases help only to the extent that they
are secret from the replacement device. A process-level or reference-assisted
test is cleaner.

### 2.2 Identity attack when the world remains accessible

Suppose the device emits \(Q=B\) unchanged and keeps no memory. If, after the
late choice, its reveal response can still measure \(Q\), it reads the requested
computational-basis fact perfectly. In erase trials it does nothing and passes
an identity-channel or Bell-recovery test perfectly.

**Consequence:** the analytic memoryless value \(P_{\rm reveal}=1/2\) is valid
only if every system correlated with the fact, especially \(Q\), is withheld
from the response module. “The challenge arrived late” is insufficient.

### 2.3 SWAP attack against a one-location memory witness

The isometry

\[
|h\rangle_B|0\rangle_M\longmapsto
|0\rangle_Q|h\rangle_M
\]

keeps the complete quantum input in \(M\). The device measures \(M\) in reveal
trials and swaps it back in erase trials. With a verifier-held reference it has

\[
P_{\rm reveal}=F_{\rm erase}=1.
\]

Yet no second record was created, and no persistent world port carried the past
fact. This is a valid quantum memory, not evidence of an observation-like
world-to-memory correlation.

**Consequence:** the repaired one-location game may certify quantum storage or
non-entanglement-breaking behavior, but it cannot identify the intended
history implementation.

### 2.4 Direct-phase and branch-dependent replacement shortcuts

If erase asks only for a particular final phase pattern, a direct unitary can
apply that phase without executing the prefix. If the untrusted device controls
post-choice processing, it may implement unrelated reveal and erase channels
after a minimal common prefix. A meaningful game must make its score depend on
information that was unavailable after the prefix except through a constrained
persistent system.

Entanglement fidelity removes a pure-state replacement and a fixed direct-phase
shortcut, but it does not identify where the quantum state was stored. Hidden
coherent environments count as quantum memory. Branch-dependent loss is another
attack: conditional fidelity must not be reported without charging every abort
as failure.

## 3. Reference-assisted repair

The verifier prepares

\[
|\Phi_d\rangle_{RB}=\frac1{\sqrt d}\sum_{h=0}^{d-1}|h\rangle_R|h\rangle_B
\]

and keeps \(R\). In reveal trials it measures \(R\) in the computational basis,
thereby defining a uniformly random fact \(H\), and sends a late query to the
allowed response port. In erase trials it asks for a recovered output \(B'\)
and performs the Bell projection

\[
F_{\rm erase}=\langle\Phi_d|\rho_{RB'}|\Phi_d\rangle.
\]

No channel that first becomes entanglement-breaking can exceed \(1/d\), because
later channels preserve separability and every separable state has overlap at
most \(1/d\) with \(|\Phi_d\rangle\). The exact matrix audit recovered

\[
F_{\rm EB}^{\max}=1/2\quad(d=2),\qquad
F_{\rm EB}^{\max}=1/4\quad(d=4).
\]

Fifty thousand random product states in each dimension respected these bounds;
the script also inserts an attaining product state, giving the ceiling to
machine precision.

The values of explicit four-history strategies are:

| Strategy | reveal | known-state erase | Bell erase | world reveal | memory reveal |
|---|---:|---:|---:|---:|---:|
| EB measure, copy, reprepare | 1 | 1 | 0.25 | 1 | 1 |
| Identity, world accessible | 1 | 1 | 1 | 1 | 0.5 |
| Coherent SWAP | 1 | 1 | 1 | 0.5 | 1 |
| Coherent fanout and inverse | 1 | 1 | 1 | 1 | 1 |

The last three columns show why a single scalar score is not sufficient to
identify a record-creating implementation.

## 4. Candidate bounds for restricted comparison classes

### 4.1 Strictly memoryless response

Assume the prefix emits \(Q\), retains no system correlated with \(H\), and the
reveal responder is denied \(Q\). Then its answer is independent of a uniform
queried bit:

\[
P_{\rm reveal}^{\rm memless}=1/2.
\]

It may preserve \(Q\) perfectly, so \(F_{\rm erase}=1\). The equal-weight score
\((P_{\rm reveal}+F_{\rm erase})/2\) is at most \(3/4\). This bound fails
immediately if \(Q\), a correlated environment, or pre-choice challenge
information reaches the responder.

### 4.2 Entanglement-breaking prefix

For a complete entanglement-breaking prefix, including all outputs that may be
used in recovery,

\[
F_{\rm erase}\le 1/d.
\]

It can reveal a computational-basis fact perfectly. Therefore the equal-weight
ceiling is \(3/4\) for \(d=2\) and \(5/8\) for \(d=4\). This is a clean class
separation, but it is fundamentally an entanglement-preservation test.

### 4.3 Coherent reversible prefix

The coherent fanout

\[
|h\rangle_Q|0\rangle_M\mapsto|h\rangle_Q|h\rangle_M
\]

allows perfect reveal by measuring \(M\) and perfect erase by applying the
inverse before the Bell test. Hence \((P_{\rm reveal},F_{\rm erase})=(1,1)\).
The SWAP strategy reaches the same pair, which is why this pair does not by
itself certify fanout or an internal history.

### 4.4 Classical retained flag with a coherent public system: qubit

A more informative intermediate class permits a coherent world output \(Q\)
but restricts the private memory to a classical instrument flag \(K\). For a
two-outcome Lüders instrument with effect \(E\), reveal guesses the past
computational-basis bit from \(K\), while erase uses the optimal flag-dependent
polar correction on \(Q\).

For the symmetric unsharp family with strength \(t\),

\[
P_{\rm reveal}=\frac{1+t}{2},\qquad
F_{\rm erase}=\frac{1+\sqrt{1-t^2}}{2}.
\]

It saturates

\[
(2P_{\rm reveal}-1)^2+(2F_{\rm erase}-1)^2=1
\]

and has equal-weight optimum

\[
\frac{P_{\rm reveal}+F_{\rm erase}}2
\le \frac12+\frac{\sqrt2}{4}
=0.8535533906.
\]

A deterministic search over 100,000 random two-outcome qubit effects found
maximum circle left-hand side `0.999999999995` and best equal-weight score
`0.8530696704`, within `4.84e-4` of the analytic symmetric optimum. No
counterexample was found in this restricted family. The inequality is a form
of standard information--disturbance/complementarity, not a new historical
bound.

### 4.5 Four histories and a late random-access query

For \(H=(H_0,H_1)\), the reveal branch requests one uniformly selected bit. A
diagonal classical-flag instrument is described by \(q(k|h)\). Its optimal
reveal probability and flag-corrected entanglement fidelity are

\[
P_{\rm reveal}
=\frac18\sum_{k,x}\max_{y\in\{0,1\}}
\sum_{h:h_x=y}q(k|h),
\]

\[
F_{\rm erase}
=\frac1{16}\sum_k\left(\sum_h\sqrt{q(k|h)}\right)^2.
\]

A 48,398-evaluation differential-evolution search over four outcomes found

\[
P_{\rm reveal}=0.7417088,qquad
F_{\rm erase}=0.8817810,qquad
\frac{P+F}{2}=0.8117448933.
\]

The optimizer reached its iteration cap, so this is a lower bound on the
comparison-class optimum, not a certified global bound. It independently
approached an exactly optimized XOR-covariant family
\(q(k|h)=r(k\oplus h)\), whose error probabilities are

\[
(a,b,b,c)=
(0.5670715054,0.1746458477,0.1746458477,0.0836367992).
\]

Writing \(v=(\sqrt a,\sqrt{2b},\sqrt c)\), the covariant equal-weight objective
is \(v^TAv\) on the unit sphere, with

\[
A=\frac12\operatorname{diag}(1,1/2,0)
+\frac18
\begin{pmatrix}1\\ \sqrt2\\ 1\end{pmatrix}
\begin{pmatrix}1&\sqrt2&1\end{pmatrix}.
\]

Its exact numerical eigenvalue is `0.8117449009`, at
\(P=0.7417173531\), \(F=0.8817724488\). A product of independent equal-strength
unsharp measurements is slightly worse (`0.8105544165`).

This is the most concrete possible mathematical subproblem found in the audit:
prove or refute that the XOR-covariant value is optimal over a carefully chosen
class of classical-memory combs. Even if proved, it would still be primarily a
new random-access information--disturbance inequality. Historical language
would not be part of the theorem.

## 5. Stronger two-location candidate

To reject both identity and SWAP shortcuts, split \(Q\) and \(M\) after the
prefix. In a reveal trial, noncommunicating responders at both locations must
answer the same late computational-basis query. In an erase trial, the systems
are reunited and a global decoder is tested against the verifier-held
reference.

Report a triple rather than one average:

\[
(P_Q,P_M,F_{\rm erase}).
\]

The relevant exact points for a qubit are

\[
\begin{array}{c|ccc}
\text{strategy}&P_Q&P_M&F_{\rm erase}\\\hline
\text{identity}&1&1/2&1\\
\text{SWAP}&1/2&1&1\\
\text{EB classical copy}&1&1&\le1/2\\
\text{coherent fanout}&1&1&1.
\end{array}
\]

Thus only coherent fanout occupies the ideal corner under the declared
dimensions and noncommunication. Twenty thousand Haar-random
qubit-to-two-qubit isometries were globally reversible by construction. The
largest sampled minimum local Helstrom success was `0.9917720`; exact coherent
fanout attains one. This numerical search is a sanity check, not a soundness
proof.

Perfect local distinguishability at both ports plus global reversibility
forces an encoding equivalent, on the tested basis, to coherently broadcasting
a classical observable. A robust version could ask how proximity to the ideal
triple constrains distance to such an isometry. That robust self-testing
question is more specific and potentially more valuable than the original
one-location game.

Risks remain:

- it resembles standard coherent copying, nondemolition measurement, GHZ
  encoding, redundant records, and quantum erasure;
- dimensions and port isolation must be enforced or bounded;
- shared hidden systems can relocate the memory;
- simultaneous local tests and global erase occur on different trials and need
  an i.i.d. or sequential-adversary model;
- a full device-independent version may reduce to existing nonlocal games or
  quantum-memory certification.

## 6. Repaired nondemolition AUDIT/ECHO game

The strongest one-step repair uses one common channel

\[
\mathcal N:A\longrightarrow BM
\]

and a verifier-held reference. The verifier always prepares
\(|\Phi_d\rangle_{RA}\), so the device sees the same maximally mixed input
before the choice.

- In `AUDIT`, the verifier measures \(R\) in a fixed computational basis and
  obtains the hidden label \(j\). It jointly checks \(B=j\) and \(M=j\).
- In `ECHO`, a decoder on \(BM\) must recover
  \(|\Phi_d\rangle_{RA'}|0\rangle_M\).

The coherent premeasurement isometry

\[
V|j\rangle_A=|j\rangle_B|j\rangle_M
\]

and its inverse attain `(AUDIT,ECHO)=(1,1)`. Unlike the one-location game,
identity and SWAP each fail AUDIT on one output. For a qubit they score
`AUDIT=1/2`, `ECHO=1`.

### 6.1 Exact ideal-point characterization

If AUDIT is exactly one, every Kraus operator of \(\mathcal N\) must satisfy

\[
K_\alpha|j\rangle=c_{\alpha j}|jj\rangle.
\]

Consequently the channel is a Schur-dephased coherent copy,

\[
\mathcal N(|j\rangle\langle k|)
=g_{jk}|jj\rangle\langle kk|,
\qquad
g_{jk}=\sum_\alpha c_{\alpha j}c_{\alpha k}^*.
\]

The Gram matrix \(G=(g_{jk})\) represents information exported to an
environment. Perfect ECHO requires that environment to be branch independent,
so \(G\) has rank one with \(|g_{jk}|=1\). The remaining freedom is a diagonal
phase on the copied basis, a branch-independent ancillary state, and arbitrary
action outside the tested input/code subspace. A matching decoder removes the
phase.

Thus, under fixed dimensions and full unconditional scoring, the exact ideal
corner does certify a coherent nondemolition premeasurement of the tested
classical observable, not merely storage. It still does not certify a unique
Hamiltonian, gate sequence, physical carrier, or interpretation.

### 6.2 Exact and random channel audit

The independent script `scratch/originality_gate/nondemolition_game.py` used
normalized Choi matrices and unitary decoders. Qubit results were:

| Channel | AUDIT | ECHO |
|---|---:|---:|
| coherent copy | 1 | 1 |
| coherent copy with arbitrary diagonal phases | 1 | 1 with matched decoder |
| identity embedding | 0.5 | 1 |
| SWAP embedding | 0.5 | 1 |
| measure-and-prepare copy | 1 | 0.5 |
| copy with environment overlap 0.5 | 1 | 0.75 |
| economic phase-covariant cloner isometry | 0.5 | 1 |
| universal symmetric cloner | 0.6667 | between 0.4523 and 0.5 for reset-unitary decoders |

The universal-cloner lower value comes from 100,000 random decoders. The upper
value is the largest normalized-Choi eigenvalue and applies only because a
unitary decoder merely rotates the pure Bell-reset target. Neither number is an
optimum for arbitrary irreversible CPTP recovery. The example is included to
seek an accidental perfect strategy; none was found, and it must not be used as
a general cloner bound.

For uniform real environment overlap \(\gamma\), the exact family obeyed

\[
P_{\rm AUDIT}=1,
\qquad
P_{\rm ECHO}=\frac{1+\gamma}{2}
\]

at all eleven sampled values from zero to one. Measure-and-prepare is the
\(\gamma=0\) endpoint. Fifty thousand Haar-random qubit-to-two-qubit
isometries had mean AUDIT `0.25003` and best sampled AUDIT `0.91751`; every
isometry has perfect ECHO with its own inverse. No non-copy isometry reached
the ideal pair.

The entanglement-breaking and identity/SWAP values are not a bound for every
notion of “classical memory.” A weak nondemolition instrument may leave a
coherent public \(B\) while retaining only a classical flag in \(M\). In the
symmetric qubit family its joint AUDIT probability and recoverable ECHO
fidelity are exactly the pair from Section 4.4,

\[
P_{\rm AUDIT}=\frac{1+t}{2},\qquad
P_{\rm ECHO}=\frac{1+\sqrt{1-t^2}}2.
\]

It reaches equal scores `0.8535533906` at \(t=1/\sqrt2\). Thus a proposed
classical null of `1/2` would be false; the comparison class must distinguish a
fully entanglement-breaking prefix from a coherent public system with only a
classical private record.

There is no constant gap against the set of all “not exactly honest” channels:
an arbitrarily weak mixture of a coherent copy with noise approaches `(1,1)`.
A useful theorem must therefore be robust, for example lower-bounding diamond
distance to the coherent-copy equivalence class from observed score deficits.

### 6.3 Remaining loopholes and novelty status

- If the device chooses its own decoder, diagonal phases and more general
  consistent code-space gauges are unavoidable. A verifier-fixed decoder tests
  a stricter implementation but becomes calibration-dependent.
- Extra hidden dimensions permit branch-independent gauge systems. A
  branch-dependent hidden environment is penalized by ECHO only if it is not
  returned to the decoder.
- Postselection can fake a high conditional ECHO score; aborts and losses must
  be failures.
- A dimension-unbounded or communicating implementation can relocate ports;
  physical separation or dimension witnesses remain necessary.
- The game certifies a coherent premeasurement channel at its declared ports.
  It cannot certify phenomenological observation or a unique internal story.

This repaired game is substantially better than the original proposal. The
specific robust joint certification problem may be worth pursuing, but the
ideal construction itself is a coherent copy followed by reversal and is close
to established nondemolition measurement and quantum-instrument theory.

## 7. Distributed records and access structures

For pure branch-conditioned environment states
\(|e_0\rangle_{S\bar S}\) and
\(|e_1\rangle_{S\bar S}\), suppose only fragments in \(S\) are accessible to
a branch-controlled inverse. Uhlmann's theorem gives the maximum recovered
control visibility

\[
V_{\max}(S)
=\max_{U_S}\left|
\langle e_0|(U_S\otimes I_{\bar S})|e_1\rangle
\right|
=F_{\rm root}\!\left(
\rho^{0}_{\bar S},\rho^{1}_{\bar S}
\right).
\]

Perfect recovery is therefore controlled by **decoupling of the inaccessible
complement**, not by how much information each accessible or inaccessible
fragment carries separately. The deterministic script
`scratch/originality_gate/distributed_records.py` enumerated every subset for
four exact encodings.

| Record family | fragments | minimum controlled fragments for perfect recovery | access behavior |
|---|---:|---:|---|
| product repetition \(|0^n\rangle,|1^n\rangle\) | 4 | 4 | every uncontrolled fragment retains the full bit |
| phase-secret GHZ \((|0^n\rangle\pm|1^n\rangle)/\sqrt2\) | 4 | 1 | every proper complement is decoupled |
| four-qubit line graph, \(|G\rangle,Z_1|G\rangle\) | 4 | 1, but subset-dependent | one-fragment sets `2` and `3` recover; `0` and `1` do not |
| five-qubit perfect-code logical states | 5 | 3 | exact `(3,5)` threshold |

For the five-qubit code, every controlled subset of size at most two had
\(V_{\max}=0\) and inaccessible Holevo information one bit; every subset of
size at least three had \(V_{\max}=1\) and inaccessible Holevo information
zero. The graph-state example shows that cardinality alone need not determine
recovery: geometry and stabilizer support matter.

### 7.1 Synergistic leakage counterexample

Let

\[
|e_0\rangle=|\Phi^+\rangle,
\qquad
|e_1\rangle=|\Phi^-\rangle.
\]

Each fragment alone is maximally mixed in both branches:

\[
\chi(H:E_1)=\chi(H:E_2)=0.
\]

Jointly the states are orthogonal:

\[
\chi(H:E_1E_2)=1,
\qquad V(S=\varnothing)=0.
\]

Access to either fragment permits a controlled local \(Z\) that maps one Bell
record to the other, giving \(V_{\max}=1\). The four-fragment GHZ-phase instance
is stronger: **every proper marginal** has zero leakage, while the complete
uncontrolled coalition has one bit and destroys all visibility.

This falsifies any recoherence criterion based on summing fragmentwise marginal
leakages. The correct object is the joint complementary state of the entire
uncontrolled coalition. It also supplies a substantive benchmark family:
recoverability access structures can range from Darwinist repetition, through
secret-shared thresholds, to geometry-dependent stabilizer patterns.

### 7.2 Novelty assessment

The equality between optimal local overlap and complement fidelity is an
application of Uhlmann's theorem, while threshold behavior is closely related
to quantum secret sharing, erasure correction, and the cleaning lemma. Those
facts are not new. A potentially new contribution would have to be narrower:

- an experimentally usable **recoherence access-structure tomography**;
- robust bounds from incomplete fragment coalitions;
- a comparison between redundant classical records and synergistic quantum
  records under the same inverse-control budget; or
- a causal-agent benchmark whose record geometry changes the minimum inverse
  coalition in a measurable way.

The synergistic example is nevertheless important even if not novel: it is a
hard correctness requirement for every future leakage ledger.

## 8. Random-cut sequential extension

Consider a reversible automaton with prefix isometries \(V_t\). After each
prefix the verifier secretly chooses either:

- a cutwise AUDIT of the declared state, transcript, and action; or
- ECHO, applying an inverse and Bell-testing every streamed input while checking
  resets.

Choosing the cut after the prefix prevents a device from learning which test it
faces before producing that boundary state. Shielding future streamed inputs
also prevents literal precomputation of functions that depend on them.

The scratch calculation `sequential_gauges.py` applied independent diagonal
phase gauges to five four-dimensional cuts. Every cut retained
`AUDIT=1` and `ECHO=1` with a matched inverse. An ungauged fixed decoder obtained
ECHO values between `0.0070` and `0.4315`. A compiled transition and a decomposed
transition with the same product had exactly zero Frobenius distance at the
tested endpoint.

### What cutwise tests can certify

Under dimensions, isolated ports, unpredictable streamed challenges, and
unconditional scoring, exact tests can certify:

- the classical transition table at every audited boundary;
- coherent recoverability of each tested prefix;
- causal dependence on inputs that were physically unavailable earlier;
- redundant records when separated ports both pass an audit; and
- absence of branch information from the declared uncontrolled complement, to
  the accuracy of the ECHO test.

### What remains gauge or implementation-dependent

They cannot distinguish:

1. cut-dependent diagonal phases removed by the matching decoder;
2. consistent hidden basis changes
   \(U'_t=G_tU_tG_{t-1}^\dagger\);
3. arbitrary action on states and subspaces never challenged;
4. a lookup-table transition from the advertised internal algorithm;
5. a compiled unitary from its chosen gate decomposition;
6. teleportation or measurement-based realizations from direct unitary
   realizations of the same comb; or
7. Stinespring dilations related by environment isometries.

If the device supplies both prefix and decoder, the test characterizes their
composition and boundary comb only up to these gauges. A random cut creates
more operational boundaries and can certify a longer causal input/output
structure, but no set of boundary tests reveals uninstrumented microhistory.

The sequential direction becomes scientifically distinctive only if it proves
a robust comb-level rigidity or memory-dimension result beyond independently
applying the one-step test at each cut.

## 9. Spectral classical-QND frontier

The most successful mathematical extension found in this audit is an exact
support function for an \(n\)-bit late random-access AUDIT versus ECHO tradeoff,
within the declared class of quantum-nondemolition instruments whose retained
outcome is classical.

Let \(h\in\{0,1\}^n\) be uniform and let \(q(k|h)\) be the classical outcome
distribution. After a uniformly random coordinate \(X\) is disclosed, the
optimal answer associated with \(k\) is a guess word \(g(k)\). Then

\[
P=\frac1{n2^n}\sum_{k,h}q(k|h)
\bigl[n-w(h,g(k))\bigr],
\]

while optimal flag-conditioned phase correction has the efficient-instrument
upper bound

\[
F\le\frac1{4^n}\sum_k
\left(\sum_h\sqrt{q(k|h)}\right)^2.
\]

For support weight \(0\le\lambda\le1\), the tested formula is

\[
\sup_q\{\lambda P+(1-\lambda)F\}
=\lambda_{\max}(A_{n,\lambda}),
\]

where

\[
(A_{n,\lambda})_{wz}
=\lambda\left(1-\frac wn\right)\delta_{wz}
+\frac{1-\lambda}{2^n}
\sqrt{\binom nw\binom nz},
\qquad 0\le w,z\le n.
\]

The Perron eigenvector gives an explicit translated instrument. If \(v_w\) is
the normalized shell vector, lift it as

\[
a(h)=v_{w(h)}/\sqrt{\binom n{w(h)}}
\]

and set

\[
q(g|h)=a(g\oplus h)^2.
\]

No additional \(1/2^n\) factor appears: translation already makes every input
column sum to one.

### 9.1 Stress-test results

`spectral_qnd_audit.py` produced these checks:

- Every explicit attainer on a seven-point \(\lambda\) grid for
  \(n=1,2,3,4\) normalized within `6.7e-16` and attained the eigenvalue within
  floating-point error.
- Unsymmetrized softmax optimizations over all \(q(k|h)\), at
  \(\lambda=0.25,0.5,0.75\), never exceeded the spectral value. Some runs,
  especially audit-heavy \(n=3,4\), remained below it by up to `0.0141`, which
  is consistent with local optimization failure and is not evidence of a
  stricter bound.
- Independent differential evolution for \(n=2,\lambda=1/2\) found
  `0.8117448881`, versus spectral `0.8117449009`.
- The largest root of
  \(\mu^3-10\mu^2+24\mu-8=0\), divided by eight, agreed with the spectral
  value to `1.1e-15`.
- Two thousand random hidden-multi-Kraus diagonal instruments at each of
  (n=1,2), including arbitrary phases and a common optimized correction,
  produced no violation of the efficient Lüders upper bound. The maximum ratio
  was `1+1e-15`; destructive hidden-environment interference reduced it as far
  as `0.342` for \(n=2\).

These calculations did not find a counterexample. They support the stated
spectral theorem but do not establish literature priority.

### 9.2 Orbital extensions

The same radial construction was tested on small transitive spaces:

- q-ary Hamming spaces \(H(n,q)\) for `(q,n)=(3,1),(3,2),(3,3),(4,2)`;
- Johnson spaces \(J(v,k)\) for `(v,k)=(4,2),(5,2),(6,3)`.

For Hamming reward \(1-w/n\), shell valencies are
\(\binom nw(q-1)^w\). For Johnson normalized-intersection reward, valencies are
\(\binom kw\binom{v-k}w\). The explicit radial instruments agreed with their
shell spectral values within `8.9e-16` over three support directions.
Unsymmetrized optimizations never exceeded them. This supports a broader
association-scheme or transitive-history-space formulation; that formulation
should be checked against known symmetrization and spherical-function results
before being called new.

### 9.3 q-ary Hamming asymptotics

A stable secular solver using binomial log-probabilities tested \(q=2,3,4\) up
to \(n=4096\). It confirms

\[
s_n(\lambda)\longrightarrow
\max\left\{\lambda,
1-\lambda\left(1-\frac1q\right)\right\},
\]

with crossing

\[
\lambda_c=\frac{q}{2q-1}.
\]

At the crossing, the finite-size excess has the observed limits

\[
n\,[s_n(\lambda_c)-\lambda_c]
\longrightarrow
\begin{cases}
1/3,&q=2,\\
1/5,&q=3,\\
1/7,&q=4,
\end{cases}
\]

strongly suggesting \(1/(2q-1)\). The eigenvector AUDIT value at
\(\lambda_c\) tends \(1/q\). The two limiting support lines have slopes
\(1\) and \(-(1-1/q)\), corresponding to limiting AUDIT values \(1\) and
\(1/q\). The derivative therefore develops a jump: a first-order
localization transition between an Echo-dominated binomial shell distribution
and a perfect-guess shell localized at \(w=0\). Away from the crossing, grid
errors decrease monotonically with \(n\). Double precision remained stable for
the support through \(n=4096\).

## 10. Ordered histories on rooted trees

To make temporal order explicit, let histories be depth-\(n\) leaves of a
rooted \(q\)-ary tree. For a guessed leaf \(g\), define reward

\[
r(g,h)=\frac{\operatorname{LCP}(g,h)}n.
\]

The orbit labelled by longest-common-prefix length \(\ell\) has valency

\[
v_\ell=(q-1)q^{n-\ell-1}\quad(\ell<n),
\qquad v_n=1.
\]

`tree_lcp_audit.py` diagonalized the corresponding shell matrices for
\(q=2,3,4,8\) through depth 256. At \(\lambda=1/2\), the support converges to
\(1/2\) from above; AUDIT tends zero like \(1/[(q-1)n]\), ECHO tends one, and
the perfect-leaf mass vanishes. Numerically located midpoint transitions tend
to \(1/2\) from above. For example, at depth 256:

| \(q\) | support at \(1/2\) | AUDIT | ECHO | midpoint transition |
|---:|---:|---:|---:|---:|
| 2 | 0.501969 | 0.0040 | 1.0000 | 0.500986 |
| 3 | 0.500982 | 0.0020 | 1.0000 | 0.500492 |
| 4 | 0.500654 | 0.0013 | 1.0000 | 0.500327 |
| 8 | 0.500280 | 0.0006 | 1.0000 | 0.500140 |

The fixed-\(\lambda\) expansions supplied by the shell moments also survived.
For \(\lambda<1/2\), errors after retaining the stated \(1/n^2\) terms decay
at the expected next order. At \(q=2,\lambda=0.25\), for example,
`n^3` times the support error tends about `0.167`; at
\(\lambda=0.4\) it tends about `1.07`. For \(\lambda>1/2\), ratios of the
exponentially small corrections approach

\[
\frac{s-\lambda}{q^{-n}}
\to\frac{\lambda(1-\lambda)}{2\lambda-1},
\]

\[
\frac{1-P}{q^{-n}}
\to\frac{(1-\lambda)^2}{(2\lambda-1)^2},
\qquad
\frac{F}{q^{-n}}
\to\frac{\lambda^2}{(2\lambda-1)^2}.
\]

Float64 loses the first two ratios once \(q^{-n}\) falls below resolution, but
the moderate-depth approach is clear and \(F/q^{-n}\) remains stable longer.

### 10.1 Inconsistent late-cut guesses

A late cut \(T\) could allow a separate guessed prefix for every \(T\), not one
consistent leaf. `tree_inconsistent_profiles.py` enumerated all binary decision
profiles at depths two, three, and four: respectively 8, 64, and 1024 profiles,
versus 4, 8, and 16 consistent leaves. It tested uniform, front-loaded, and
back-loaded cut distributions at eleven support weights.

No inconsistent profile improved the consistent-leaf support; the maximum and
minimum numerical improvement were both exactly zero. Many inconsistent
profiles are suboptimal or tie a consistent optimum.

There is also a general proof of equality of the optimum values. Under a
uniform leaf, the correctness indicator \(I_t\) for a declared depth-\(t\)
prefix always has marginal probability \(q^{-t}\), whether or not the guesses
at different depths are compatible. A consistent leaf makes these indicators
nested and therefore comonotone. Consequently
\(R_{\rm leaf}=\sum_t\pi_t I_t\) dominates the reward of every inconsistent
profile in convex order for any nonnegative cut weights \(\pi_t\). The largest
eigenvalue is the unique root above the diagonal edge of

\[
1=(1-\lambda)\,\mathbb E[(s-\lambda R)^{-1}].
\]

The integrand is increasing and convex in \(R\), so the comonotone reward has
the largest root. Thus late, independently decoded cuts do not raise the QND
support for any depth or alphabet. They can still produce additional tied or
suboptimal strategies, so this equality of optimum values is not a rigidity
statement.

The tree theorem is the candidate most naturally about ordered histories, but
its mathematical core is still a group/orbit-symmetrized
information--disturbance problem. Temporal interpretation alone does not
establish novelty.

## 11. Other computationally testable originality candidates

### Candidate A: complementary-basis reveal plus Bell erase

Let the late reveal challenge select an \(X\) or \(Z\) fact. A coherent qubit
memory answers either basis perfectly and can be returned for perfect Bell
recovery. A measure-first classical strategy has optimal BB84/Breidbart success

\[
P_{\rm reveal}^{\rm classical}=cos^2(\pi/8)=0.8535533906,
\]

and an entanglement-breaking prefix has Bell fidelity at most \(1/2\).

**Assessment:** computationally clean but not original enough. It is the
standard distinction between classical and quantum memory under delayed basis
information. Zero-record terminology adds no operational prediction.

### Candidate B: robust coherent broadcast-and-erase self-test

Use the triple \((P_Q,P_M,F_{\rm erase})\) and derive a robust lower bound on
the distance from every swap, identity, EB, and bounded-classical-memory
strategy. Numerically optimize low-dimensional combs, then seek an SDP/NPA or
dimension-bounded relaxation.

**Assessment:** the strongest candidate from this audit. It has explicit
anti-shortcut ports and an experimentally interpretable target. Its novelty
depends on whether an equivalent robust instrument/self-testing theorem already
exists. Literature clearance is mandatory before investment.

### Candidate C: zero-record oracle-query separation

The script simulated Bernstein--Vazirani phase queries for every secret through
eight bits. One query returned the exact secret with unit fidelity, while the
phase ancilla remained unchanged. This already realizes a clean predicate
readout without a work transcript.

**Assessment:** reject as a novelty direction. Standard phase-oracle algorithms
already have this structure. Relabeling their query superposition as histories
does not create a new separation. A viable query-complexity result must add a
constraint not satisfied by ordinary clean quantum algorithms, such as a
physically persistent unknown process, causal ports, or composable transcript
privacy against an external verifier.

### Candidate D: coarse-graining versus recovery

Numerical optimization could compare recoverability for nested predicates, but
the conjecture is undefined until the tasks share the same input ensemble,
allowed controls, output encoding, and recovery class. With task-dependent
oracles or gate costs, a finer predicate can be native while a coarse function
requires extra computation; with unrestricted recovery, forgetting more cannot
hurt. Either outcome can be manufactured by changing the task.

**Assessment:** do not simulate or claim monotonicity until a common resource
theory is fixed. The current formulation is too plastic to survive an
originality gate.

## 12. Recommendation and kill conditions

Do not promote the naive common-prefix reveal-or-erase game. It is falsified by
a one-line entanglement-breaking strategy.

Do not describe the Bell-reference one-location repair as certification of a
history. At most it certifies preservation of quantum information under strict
access assumptions; SWAP is a perfect counterexample to the stronger reading.

Continue only with one of these sharply delimited goals:

1. Prove a robust AUDIT/ECHO rigidity bound for the nondemolition premeasurement
   equivalence class, or determine that it follows immediately from existing
   instrument-reversal results.
2. Prove a nontrivial four-history information--disturbance bound for a declared
   classical-memory comb class, ideally resolving whether `0.8117449009` is the
   global equal-weight value.
3. Develop access-structure recoherence benchmarks only if they add a robust or
   experimentally economical result beyond Uhlmann recovery and quantum secret
   sharing.

Abandon the originality claim if either reduces directly to an existing
which-path duality, quantum nondemolition measurement, quantum-memory witness,
or known self-testing result without a new bound or access model.

## 13. Reproduction

The environment did not contain CVXPY or another SDP package, so no SDP claim is
made. Exact matrix calculations, analytic eigenvalue optimization, SciPy global
optimization, and deterministic random searches were used.

Run:

```powershell
python scratch/originality_gate/audit_delayed_choice.py
python scratch/originality_gate/nondemolition_game.py
python scratch/originality_gate/distributed_records.py
python scratch/originality_gate/sequential_gauges.py
python scratch/originality_gate/spectral_qnd_audit.py
python scratch/originality_gate/orbital_qnd_audit.py
python scratch/originality_gate/hamming_asymptotics.py
python scratch/originality_gate/tree_lcp_audit.py
python scratch/originality_gate/tree_inconsistent_profiles.py
```

The programs rewrite only files below `scratch/originality_gate/`. They do not
modify the release simulator, manuscript, figures, data, or PDF.

## 14. Online causal classical-memory audit

The collective spectral theorem permits one joint diagonal instrument on an
entire history. A stricter temporal null was therefore tested: event qudits
arrive sequentially, each returned carrier is immediately sequestered, and the
device may retain unlimited classical state but no quantum system between
slots. The late branch either asks for a uniformly random past symbol or
returns all carriers for flag-conditioned joint Bell recovery.

For local dimension \(Q\), define

\[
g_Q(p)=\frac1Q
\left(\sqrt p+\sqrt{(Q-1)(1-p)}\right)^2.
\]

The exact candidate, subsequently proved in the theory audit, is

\[
C_{n,Q}(\lambda)=
\max_{1/Q\le p\le1}
\left[\lambda p+(1-\lambda)g_Q(p)^n\right].
\]

It is attained by the same symmetric weak Lüders instrument independently at
every slot. The audit actively searched for improvements from nonidentical
strengths, measuring only a subset of slots, time sharing, adaptive classical
policies, larger outcome alphabets, shared classical randomness, and non-QND
local effects.

### 14.1 Binary exact comparisons

For \(Q=2\), adaptive binary and ternary trees were optimized through five
slots. An independent differential-evolution search over the complete
three-node binary tree at \(n=2\) also reached the formula. No adaptive search
exceeded it. The strict gaps to the unrestricted collective Hamming support
\(G_{n,2}\) were:

| \(n\) | \(G-C\), \(\lambda=0.25\) | \(G-C\), \(\lambda=0.50\) | \(G-C\), \(\lambda=0.75\) |
|---:|---:|---:|---:|
| 2 | 0.0000342132 | 0.0011904845 | 0.0033456107 |
| 3 | 0.0000213139 | 0.0010673790 | 0.0077738007 |
| 4 | 0.0000136723 | 0.0007793253 | 0.0123299960 |
| 5 | 0.0000093736 | 0.0005598070 | 0.0139616410 |

For example, at \(n=2,\lambda=1/2\),

\[
C_{2,2}=0.8105544164707321,
\qquad
G_{2,2}=0.8117449009293669.
\]

The gap is small but reproducibly nonzero. Optimizing a fixed measured subset
selected all \(n\) slots, and unrestricted nonidentical local strengths
collapsed to equal strengths.

### 14.2 q-ary and non-QND checks

Adaptive q-ary trees were also tested for \(Q=3\), \(n=2,3,4\), and \(Q=4\),
\(n=2,3\), at the same three support weights. No run exceeded the i.i.d.
formula, and every formula value was strictly below its collective spectral
counterpart. Differential evolution over nonidentical strengths equalized all
slots to within `1.83e-6` in the worst run.

To probe the step that extends the theorem beyond diagonal local channels,
20,000 random non-diagonal positive effects were sampled for each
\(Q=2,3,4,8\). None violated the local effect-fidelity envelope. The largest
ratios of actual value to envelope were respectively `0.999973`, `0.990550`,
`0.952760`, and `0.849358`.

Classical randomization cannot improve the support because it only forms a
convex mixture of refined deterministic paths. The numerical results therefore
found no cheating strategy based on adaptation, enlarged classical state,
time sharing, correlated classical coins, non-QND effects, or joint final
recovery. They support the exact product-versus-collective theorem; they do not
establish literature priority or guarantee experimental enforcement of the
sequestering and no-side-channel assumptions.

The full theorem and proof were independently preserved in
`notes/originality_gate_online_backup.md` because concurrent theory-note edits
created a file-version race.

Additional reproduction commands:

```powershell
python scratch/originality_gate/online_causal_qnd.py
python scratch/originality_gate/online_causal_qary.py
python scratch/originality_gate/asymptotic_formula_checks.py
```

## 15. Streaming parity: a finite-automaton-specific variant

A more history-specific late query asks for the parity of \(n\) sequential
binary events. With no cross-slot quantum memory, a product weak instrument
with local biases \(t_i\in[0,1]\) has

\[
P_{\rm parity}=\frac{1+\prod_i t_i}{2},
\qquad
F_{\rm echo}=\prod_i
\frac{1+\sqrt{1-t_i^2}}2.
\]

For a support direction \(\lambda\), the online product value is therefore

\[
C_n^{\rm parity}(\lambda)
=\max_{0\le t_i\le1}
\left[
\frac\lambda2\left(1+\prod_i t_i\right)
+(1-\lambda)\prod_i\frac{1+\sqrt{1-t_i^2}}2
\right].
\]

Every fully interior stationary point has equal \(t_i\): differentiating with
respect to each bias makes \((1-\sqrt{1-t_i^2})/\sqrt{1-t_i^2}\) independent
of \(i\). If any \(t_i=0\), parity information vanishes and setting all biases
to zero maximizes ECHO. A point with \(t_i=1\) cannot be optimal for interior
\(\lambda\), because decreasing it produces a square-root ECHO gain against a
linear AUDIT loss. Thus the product optimization reduces to one common bias
\(t\), including the no-record endpoint.

The collective classical-QND parity frontier instead depends only on the two
parity sectors and is independent of history length:

\[
G^{\rm parity}(\lambda)
=\frac{1+\sqrt{\lambda^2+(1-\lambda)^2}}2.
\]

### 15.1 Adversarial results

`streaming_parity.py` optimized unrestricted nonidentical biases and complete
adaptive binary outcome trees for \(n=2,3,4,5\). Across six support weights,
the largest asymmetric-product excess was `3.6e-14`, and the largest adaptive
excess was `6.5e-14`. No counterexample to the symmetric product value was
found.

Representative results are:

| \(n\) | \(\lambda\) | online score | optimal \(t\) | collective minus online |
|---:|---:|---:|---:|---:|
| 2 | 0.50 | 0.7500000000 | 0 | 0.1035533906 |
| 2 | 0.75 | 0.8250000000 | 0.979796 | 0.0702847075 |
| 3 | 0.50 | 0.7500000000 | 0 | 0.1035533906 |
| 3 | 0.75 | 0.7859840863 | 0.994819 | 0.1093006212 |
| 4 | 0.75 | 0.7671219899 | 0.998843 | 0.1281627176 |
| 5 | 0.75 | 0.7582575144 | 0.999740 | 0.1370271932 |

The product null has a notable phase structure. At \(n=2\), the nonzero-bias
solution activates continuously at \(\lambda=1/2\). For \(n>2\), weak local
records initially lose ECHO at order \(t^2\) but acquire parity information
only at order \(t^n\); the no-record endpoint therefore remains locally stable
until a strong-record solution overtakes it discontinuously. Numerically, the
switch weights and biases were:

| \(n\) | switch weight | bias at switch |
|---:|---:|---:|
| 2 | 0.5000000 | 0 |
| 3 | 0.6247789 | 0.971737 |
| 4 | 0.6495455 | 0.996176 |
| 5 | 0.6588955 | 0.999294 |

The sequence suggests a limiting switch at \(2/3\). This discontinuity is a
property of the optimized support problem, not a thermodynamic phase
transition.

An honest device with one coherent parity accumulator can update
\(M\leftarrow M\oplus h_t\) at every slot, answer the parity in AUDIT, or use
the returned carriers to invert all updates in ECHO. It reaches
\((P_{\rm parity},F_{\rm echo})=(1,1)\). This makes the parity variant a sharper
finite-automaton benchmark than independent random-coordinate retrieval,
although its relation to existing quantum-memory and automata witnesses still
requires literature clearance.

Reproduce with:

```powershell
python scratch/originality_gate/streaming_parity.py
```

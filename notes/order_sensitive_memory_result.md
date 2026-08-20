# Temporal Order Changes a Bounded-Memory Information--Recovery Frontier

**Canonical result note**<br>
**Status:** exact endpoint theorem, exact rank-two order dichotomy, quantitative
near-endpoint separation, and numerical interior lower bounds<br>
**Literature cutoff:** 19 August 2026<br>
**Scope:** binary rank-two checks, one persistent coherent qubit, and the
trusted streaming interface below

## 1. Result in one paragraph

Temporal order can change how much unknown quantum information remains
recoverable after a streamed device has learned a classical syndrome, even
when the code, terminal coherent dimension, and static information--disturbance
bound are unchanged. Consider the same rank-two binary check under the grouped
order \(A,A,B,B\) and the interleaved order \(A,B,A,B\). With one persistent
coherent qubit, unlimited genuinely classical memory, and a choice made only
after the stream between complete-syndrome **AUDIT** and all-carrier EPR
**RETURN**, the grouped order attains the full static qubit frontier. In
particular, perfect AUDIT is compatible with \(F_{\rm R}=1/2\). For the
interleaved order, every adaptive non-QND strategy instead satisfies the exact
and attainable endpoint

\[
P_{\rm A}=1\quad\Longrightarrow\quad F_{\rm R}\leq\frac14.
\]

The separation is caused solely by a coordinate permutation. A robust theorem
also certifies a strict support-function gap for
\(0.997339868377\ldots\leq\lambda<1\). The complete interleaved interior
frontier remains open.

## 2. Frozen operational model

At slot \(i\in\{1,\ldots,n\}\), a verifier prepares

\[
|\Phi^+\rangle_{R_iA_i}
=\frac{|00\rangle+|11\rangle}{\sqrt2}
\]

and sends \(A_i\) to a streaming device. The device may apply an arbitrary
adaptive instrument to \(A_i\), fresh local ancillas, and a persistent memory
\(M_{i-1}\). It emits \(B_i\), which is immediately sequestered, and carries
\(M_i\) to the next slot with

\[
\dim M_i\leq2.
\]

The device may retain an unrestricted finite classical transcript. That
transcript must be genuinely dephased: its coherent purification is not a free
side channel. No other coherent degree of freedom, discarded environment, or
previous carrier may cross a temporal cut uncharged.

Only after the last slot does the verifier choose the task.

In **AUDIT**, the references are measured in the computational basis, giving
a uniform word \(X\in\mathbb F_2^n\). Using only the classical transcript and
terminal qubit, not the sequestered carriers, the device must output the full
syndrome \(HX\in\mathbb F_2^2\). Its unconditional success probability is
\(P_{\rm A}\).

In **RETURN**, all carriers are supplied to a transcript-conditioned decoder.
The decoder may act jointly on \(B_1\cdots B_nM_n\), but has no reference
access. The verifier tests restoration of all \(n\) EPR pairs and reset of the
charged memory. The unconditional entanglement fidelity, with failure and no
postselection hidden, is \(F_{\rm R}\).

For audit weight \(0\leq\lambda\leq1\), define

\[
\beta^{\rm stream}_{H,2}(\lambda)
=\sup_{\mathcal S}
\left[\lambda P_{\rm A}(\mathcal S)
+(1-\lambda)F_{\rm R}(\mathcal S)\right],
\]

where the supremum includes all finite-outcome adaptive non-QND streamed
strategies satisfying this interface.

## 3. The canonical order-sensitive pair

The two four-slot matrices are

\[
H_{\rm G}=
\begin{pmatrix}
1&1&0&0\\
0&0&1&1
\end{pmatrix},
\qquad
H_{\rm I}=
\begin{pmatrix}
1&0&1&0\\
0&1&0&1
\end{pmatrix}.
\]

They are the same rank-two code up to a coordinate permutation. Both obey the
same universal terminal-qubit relaxation

\[
B_{4,2}(\lambda)
=\frac{1+\sqrt{\lambda^2+(1-\lambda)^2}}2.
\tag{1}
\]

This is a static ceiling. It does not assert that every temporal ordering can
realize the corresponding terminal instrument through a one-qubit bond.

### Theorem 1: grouped order attains the entire static boundary

For \(H_{\rm G}\), every point

\[
P_{\rm A}=\frac{1+t}{2},
\qquad
F_{\rm R}=\frac{1+\sqrt{1-t^2}}2,
\qquad 0\leq t\leq1,
\tag{2}
\]

is attainable with one persistent qubit. Consequently,

\[
\beta^{\rm stream}_{H_{\rm G},2}(\lambda)=B_{4,2}(\lambda)
\quad\text{for every }\lambda.
\tag{3}
\]

The construction accumulates the first parity, weakly records it with
strength \(t\), and then hands the same qubit over to the second parity. The
classical weak-measurement flag and terminal qubit identify the two-bit
syndrome in AUDIT; reversing the four controlled updates resets the qubit in
RETURN. At perfect AUDIT, equation (2) gives \(F_{\rm R}=1/2\).

### Theorem 2: exact interleaved perfect-AUDIT endpoint

For \(H_{\rm I}\), every admissible streamed strategy satisfies

\[
P_{\rm A}=1\quad\Longrightarrow\quad F_{\rm R}\leq\frac14.
\tag{4}
\]

The bound is attained. Hence

\[
\max\{F_{\rm R}:P_{\rm A}=1\}=\frac14.
\tag{5}
\]

An attainer coherently stores
\(S_1=X_1\oplus X_3\), projectively records \(X_2\) and \(X_4\), and obtains
\(S_2=X_2\oplus X_4\) from their classical parity. On RETURN it undoes the two
coherent controlled operations. The two projective records each contribute
optimal EPR recovery fidelity \(1/2\), giving \(F_{\rm R}=1/4\).

Theorem 2 is not restricted to this construction class. It quantifies over
arbitrary adaptive non-QND instruments and arbitrary transcript-conditioned
joint RETURN decoders allowed by the frozen model.

## 4. Why the converse works

The exact proof has five steps.

First, reveal every Kraus refinement as additional classical transcript. This
is a relaxation: both decoders may ignore the extra flag, while no coherent
cut dimension is added. Every complete refined transcript \(c\) is therefore
a single sequential Kraus leaf \(K_c\), with basis likelihood

\[
p_c(x)=\|K_c|x\rangle\|^2.
\]

Second, perfect AUDIT forces the terminal-qubit states belonging to different
supported syndromes in a leaf to have orthogonal supports. A qubit can support
at most two of the four syndromes.

Third, cut the interleaved stream after slot two. For prefix \(z\) and suffix
\(y\), causality and sequestration give

\[
K_c|z,y\rangle=(I_{B_1B_2}\otimes L_y)|\psi_z\rangle,
\qquad
L_y:M_2\longrightarrow B_3B_4M_4.
\tag{6}
\]

If a fixed suffix supports both syndrome labels available to the leaf, their
terminal supports are orthogonal. Thus \(L_y\) has rank two and, because its
domain is a qubit, is injective. Prefix vectors rejected after the cut must
therefore already have been zero. If no suffix supports both labels, at most
one word survives for each suffix. Either way,

\[
|\operatorname{supp}p_c|\leq4.
\tag{7}
\]

Fourth, optimal flagged polar recovery followed by computational-basis
pinching gives, with \(D=16\),

\[
F_{{\rm R},c}
\leq\frac1{D^2}
\left(\sum_x\sqrt{p_c(x)}\right)^2.
\tag{8}
\]

Cauchy--Schwarz, equation (7), and instrument completeness yield

\[
F_{\rm R}
\leq\frac4{256}\sum_{c,x}p_c(x)
=\frac14.
\]

Finally, the explicit construction above proves equality. The proof neither
diagonalizes the instrument nor assumes QND, Clifford, stabilizer, Lüders, or
covariant dynamics. Independent red-team analysis found no escape through
prefix--carrier entanglement, adaptive flags, Kraus refinement, or a general
joint recovery.

## 5. General rank-two order theorem

For an ordered rank-two binary matrix, define the standard trellis
connectivity

\[
\tau(H)=\max_{1\leq i<n}
\left[
\operatorname{rank}H_{\leq i}
+\operatorname{rank}H_{>i}
-\operatorname{rank}H
\right].
\tag{9}
\]

The invariant and its coding-theoretic meaning are prior art; they are not a
new definition. For rank two, \(\tau(H)=2\) exactly when some cut is **full
crossing**:

\[
\operatorname{rank}H_{\leq i}
=\operatorname{rank}H_{>i}=2.
\]

### Theorem 3: full-crossing endpoint and width-one converse

Let one qubit be the only coherent system crossing every cut.

If \(H\) has a full-crossing cut, then every perfect-AUDIT strategy obeys

\[
F_{\rm R}\leq\frac14.
\tag{10}
\]

For a length-\(n\) leaf, the proof generalizes equation (7) to support at
most \(2^{n-2}\) of the \(2^n\) computational words; equation (8) then gives
the same ratio \(1/4\). This is an upper theorem. Attainability at \(1/4\)
may still depend on the ordered columns.

Conversely, if \(\tau(H)\leq1\), a one-qubit accumulate--weakly-record--handoff
strategy attains the complete static boundary (2). Thus, for every binary
rank-two check order,

\[
\begin{array}{ll}
\tau(H)\leq1 &: \text{the complete static qubit frontier is attainable},\\
\tau(H)=2 &: \text{perfect AUDIT forces }F_{\rm R}\leq1/4.
\end{array}
\tag{11}
\]

This dichotomy holds for arbitrary length, including repeated and zero
columns. It is not a claimed formula \(F_{\rm R}=2^{-\tau(H)}\) beyond this
rank-two theorem.

## 6. Quantitative robustness

Let a full-crossing cut occur after slot \(i\), put

\[
J=2^{n-i},
\qquad
C_i=4\sqrt J,
\qquad
\epsilon=1-P_{\rm A}.
\]

### Theorem 4: robust full-crossing bound

Every admissible strategy satisfies

\[
\boxed{
F_{\rm R}
\leq\frac14
+\frac{C_i}{2}\sqrt\epsilon
+\frac{\sqrt{3C_i}}2\epsilon^{1/4}.}
\tag{12}
\]

For the four-slot interleaved instance, a sharpened leaf estimate gives
\(C_i=7\), hence

\[
F_{\rm R}
\leq\frac14
+\frac72\sqrt{1-P_{\rm A}}
+\frac{\sqrt{21}}2(1-P_{\rm A})^{1/4}.
\tag{13}
\]

The proof puts each leaf in a right-canonical cut gauge, shows that the two
most likely syndromes contain all but the AUDIT-error mass, replaces exact
injectivity by a smallest-singular-value bound, and obtains a set of exactly
\(D/4\) basis words outside which the leaf likelihood is at most
\(C_i\sqrt{m_ce_c}\). Applying equation (8), summing completeness, and using
Cauchy--Schwarz and Hölder gives equation (12). A separate audit checked the
rank-deficient gauge, sign case in the singular-value lemma, exact cardinality
needed by the pinching step, constants, and endpoint domain.

Combining equation (13) with the static ceiling proves the explicit strict
order gap

\[
\beta^{\rm stream}_{H_{\rm I},2}(\lambda)
<B_{4,2}(\lambda)
\quad\text{for}\quad
0.997339868377\ldots\leq\lambda<1.
\tag{14}
\]

The endpoint \(\lambda=1\) is deliberately excluded: both support functions
then equal one because only \(P_{\rm A}\) is weighted. The exact constrained
endpoint (5) supplies the secondary separation at \(P_{\rm A}=1\).

Equation (12) is non-sharp. The legal weak-record family

\[
P_{\rm A}=\frac{1+t^2}{2},
\qquad
F_{\rm R}=\left(\frac{1+\sqrt{1-t^2}}2\right)^2
\]

has, for \(\epsilon=1-P_{\rm A}\),

\[
F_{\rm R}=\frac14+\sqrt{\frac\epsilon2}+\frac\epsilon2.
\tag{15}
\]

Therefore no valid endpoint correction can be \(o(\sqrt\epsilon)\). The
\(\epsilon^{1/4}\) term in equation (12) is a proof loss, not evidence for the
true scaling.

## 7. Numerical status of the unknown interior

Numerics are lower-bound searches only. An unrestricted complex binary-tree
optimization over adaptive non-QND node instruments and
transcript-conditioned qubit POVMs found, at \(\lambda=1/2\),

\[
(P_{\rm A},F_{\rm R})
\approx(0.6446434022,0.8662314902),
\qquad
\lambda P_{\rm A}+(1-\lambda)F_{\rm R}
\approx0.7554374462.
\tag{16}
\]

This exceeds the best balanced score \(0.75\) in the simple weak-record
family, so that family is not the full frontier. It remains well below the
static value \(B_{4,2}(1/2)=0.8535533906\), but the missing interval is not a
certified gap at balanced weight.

High-AUDIT runs approach the exact endpoint: approximately
\((0.9982256,0.2808604)\) at \(\lambda=0.9\) and
\((0.9999871,0.2525485)\) at \(\lambda=0.99\). These values are feasible lower
bounds, not optimality certificates. No closed form, exact support function,
or full-interior dual certificate is currently known for \(H_{\rm I}\).

## 8. What is occupied and what may be new

The broad ingredients are occupied. Quantum comb memory cost with free
classical assistance is established [@bisio2012memory]. Sequential quantum
generation already relates persistent ancilla dimension and temporal order to
Schmidt/MPS width [@schoen2005sequential; @li2022emitters]. Trellis
connectivity, coordinate-order dependence, and matroid pathwidth are classical
coding prior art [@forney1994trellis; @mceliece1996bcjr;
@kashyap2008pathwidth]. Quantum trellises and coherent code decoding are also
established [@ollivier2006trellises; @piveteau2022message;
@piveteau2025belief]. Bounded-coherent-memory testers and discrimination
hierarchies are active frameworks [@ohst2026memory; @zonnios2026bounded].
Finally, the numerical value \(1/4\) occurs in spatial nondestructive
discrimination [@bilash2024nondestructive; @lim2025local]. None of those
headlines, invariants, or values should be claimed as new.

The defensible candidate contribution is the conjunction:

> Within the primary literature located through 19 August 2026, this appears
> to be the first exact separation in a late-choice complete-syndrome
> AUDIT/all-carrier EPR-RETURN game produced solely by permuting the temporal
> coordinates of one rank-two linear code under a one-qubit coherent-memory
> constraint. A full-crossing order forces \(F_{\rm R}\leq1/4\) at perfect
> audit, while a noncrossing order attains \(F_{\rm R}=1/2\).

The qualifiers “appears,” the literature date, exact interface, rank, memory
dimension, and perfect-AUDIT condition are essential. The priority search is a
focused collision audit, not proof of priority. Its residual collision risk is
low to moderate, especially through an unrecognized corollary of a quantum-comb
rank theorem, tensor-network cut lemma, or recoverability inequality.

## 9. Limitations and falsifiers

The result does not solve the complete interleaved tradeoff curve. It proves
one exact endpoint and a conservative strict interval near it. The robust
exponent is non-sharp, and the general full-crossing \(1/4\) upper bound need
not always be attainable.

The theorem is limited to two binary checks and one persistent qubit. It does
not establish a general \(2^{-\tau}\) law, a higher-rank characterization, a
device-independent memory witness, or the physical dimension of an
uncharacterized laboratory system.

The access assumptions are structural, not cosmetic. The proof no longer
applies if AUDIT receives emitted carriers, if a carrier can re-enter before
commitment, if the classical transcript keeps an uncharged coherent
purification, if the process or decoder accesses the references, or if RETURN
is postselected. An experimental violation would exclude only the frozen
trusted-interface model after leakage and statistical allowances; it would not
identify a unique hidden mechanism.

Nothing here tests a quantum interpretation, macroscopic branching,
consciousness, objective collapse, or non-unitarity. The result is an
operational theorem about temporal coherent-memory constraints.

## 10. Reproducibility and source map

The detailed exact proof is in **notes/order_gap_analytic.md**; the quantitative
proof is in **notes/order_gap_robust.md**. Both proofs were independently
audited before consolidation. The literature and collision audit is recorded
in **notes/order_gap_priority.md** and
**notes/literature_review_bounded_coherent_memory.md**, with primary-source
metadata in **references/library.bib**.

The public Python surface is **carmenq.order_sensitive** in
**src/carmenq/order_sensitive.py**. It provides the two canonical matrices,
binary-rank and trellis-connectivity utilities, full-crossing-cut detection,
the static support ceiling, the exact grouped frontier, and explicit metadata
for the interleaved perfect-AUDIT endpoint. It intentionally does not label a
numerical interior curve as exact. Regression tests are in
**tests/test_order_sensitive.py** and **tests/test_public_api.py**.

For independent computation,
**scripts/classify_order_sensitive_checks.py** enumerates the four-column
order classes. It finds 78 full-rank nonzero-column sequences and nine classes
after quotienting by row-basis changes and time reversal: four have
\(\tau=1\) and five have \(\tau=2\). This is a structural check, not a
replacement for the analytic theorem. Exploratory variational searches used
to falsify simple interior ansatzes are not part of the stable package and do
not enter any proved claim.

## 11. Next theorem, not a claim already made

The highest-value unresolved target is the exact function

\[
\beta^{\rm stream}_{H_{\rm I},2}(\lambda),
\qquad 0<\lambda<1,
\]

including attaining strategies and analytic or computer-assisted upper
certificates. A sharp \(O(\sqrt{1-P_{\rm A}})\) robust converse would be a
second meaningful advance. Higher-rank work should begin only after checking
whether a proposed generalization contains information beyond standard
trellis width and comb memory cost.

# Priority and Collision Audit: Temporal-Order AUDIT/RETURN Gap

**Search date and cutoff:** 20 August 2026

**Scope:** focused primary-literature search, not a proof of priority

**Object audited:** `notes/order_gap_analytic.md`,
`notes/order_gap_robust.md`, `notes/order_sensitive_memory_result.md`, and
`notes/full_crossing_dimension_bound.md`

## Bottom line

No exact collision was located for the following conjunction: a fixed rank-two
binary syndrome code; two streams differing only by a column permutation; one
persistent coherent qubit with an unrestricted classical transcript; sequestered
carrier outputs; a choice made only after the stream between perfect syndrome
AUDIT and all-carrier EPR RETURN; and the exact separation

\[
F_{\rm R}^{\rm grouped}=\tfrac12,
\qquad
F_{\rm R}^{\rm interleaved}=\tfrac14
\quad\text{at }P_{\rm A}=1.
\]

This narrow result therefore survives the present priority attack. Its broad
ingredients do not. Temporal order changing quantum memory cost, trellis
connectivity, quantum processing of linear codes, bounded-coherent-memory
strategies, and a nondestructive-discrimination value of \(1/4\) are all occupied.
The publishable object, if the proof survives independent review, is the *exact
late-choice information--recoverability separation under a column permutation*,
not any one of those ingredients.

The endpoint proof now has a finite-field temporal generalization. If the
ordered check matrix splits into \(m\) consecutive nonempty blocks that each
retain syndrome rank \(r\), \(N=q^r\), and the charged coherent dimension at
boundary \(j\) is \(d_j\), perfect AUDIT implies

\[
F_{\rm R}\leq\prod_{j=1}^m\min\{1,d_j/N\}.
\]

For uniform dimension \(d\), this becomes \(F_{\rm R}\leq(d/N)^m\), capped
at one. The uniform law is tight on the repeated-identity family for
\(d=q^k\). No exact
primary-source collision was located for this temporal power law in the
declared late AUDIT/all-carrier RETURN interface. Its proof is elementary
enough that the risk of it being an unstated corollary of broader comb or
tensor-network machinery remains material.

If \(\mu(H)\) denotes the maximum number of consecutive full-rank blocks in
an ordered partition, the uniform statement is
\(F_{\rm R}\leq(d/q^r)^{\mu(H)}\), capped at one. The maximum is found by
closing each block at its earliest full-rank endpoint. Trellis sectionalization
itself is established [@lafourcade1996sectionalization]; neither the partition
nor its greedy computation is assigned novelty here.

The matrix-level form supplies an exact exponential order separation. With
the same multiset containing \(m\) copies of each basis column, a batched order
has \(\mu=1\) and a cycled order has \(\mu=m\). For \(d=q^k\),
\(1\leq k<r\), the exact perfect-AUDIT return optima differ by
\((q^r/d)^{m-1}\). No searched source combined trellis sectionalization with
this late-choice recovery exponent; the claim remains scoped to the frozen
interface.

## What was tested

The two matrices

\[
H_{\rm G}=\begin{pmatrix}1&1&0&0\\0&0&1&1\end{pmatrix},
\qquad
H_{\rm I}=\begin{pmatrix}1&0&1&0\\0&1&0&1\end{pmatrix}
\]

represent the same code up to coordinate order. The claimed theorem quantifies
over arbitrary adaptive instruments, not only QND syndrome accumulators. The
order dichotomy in the notes is specifically the two-check statement

\[
\tau(H)=2,\ d=2,\ P_{\rm A}=1
\quad\Longrightarrow\quad F_{\rm R}\leq\tfrac14,
\qquad
\tau(H)=\max_i\bigl(
\operatorname{rank}H_{\leq i}+\operatorname{rank}H_{>i}
-\operatorname{rank}H\bigr).
\]

The search did **not** treat a general formula \(F_{\rm R}=2^{-\tau(H)}\) as
proved. The separate higher-rank theorem uses the stronger hypothesis of a
partition into full-rank blocks and gives a product of boundary-dimension
factors; it is not a formula in \(\tau(H)\) alone.

## Collision map

| Primary line of work | What it already occupies | Why it does not imply the audited theorem |
|---|---|---|
| Classical minimal trellises, sectionalization, and matroid pathwidth [@forney1994trellis; @mceliece1996bcjr; @lafourcade1996sectionalization; @kashyap2008pathwidth; @sheshadri2026trellis] | The cut expression \(\operatorname{rank}H_L+\operatorname{rank}H_R-\operatorname{rank}H\), coordinate-order dependence, and consecutive trellis partitioning are established. | These works contain neither quantum disturbance nor a late audit-versus-EPR-return objective. **The invariants and partitions themselves are not new.** |
| Sequential quantum-state generation [@schoen2005sequential; @li2022emitters] | Persistent ancilla dimension is governed by MPS/Schmidt width. Li--Economou--Barnes explicitly show that emission order changes the minimum number of coherent emitters, even from linear scaling to two emitters for one graph family. | This is fixed target-state generation. There is no classical syndrome transcript, late task choice, or recovery of unknown input entanglement. |
| Graph-state scheduling [@elman2025scheduling] | Quantum scheduling cost is exactly tied to path decompositions; order--width connections in a quantum protocol are therefore occupied. | The resource is the number of active graph-state qubits, not disturbance of streamed EPR halves after extracting a syndrome. |
| Quantum trellises and quantum code decoding [@ollivier2006trellises; @piveteau2022message; @piveteau2025belief] | Trellises for stabilizer codes are old. BPQM gives coherent message passing for classical linear codes, and the 2025 preprint gives optimal quantum decoders for every code with an efficient trellis, including deferred measurement and uncomputation. | Those decoders discriminate classical codewords encoded into quantum channel outputs. They may use the received quantum block/message registers and are not constrained to return every input carrier entangled with a reference. There is no late RETURN branch. |
| Quantum factor-graph processing [@mandal2026abelian] | Check, equality, homomorphism, marginalization, and automorphism factors over finite abelian groups already have quantum message-update rules. | It does not impose a one-qubit temporal bond or optimize an information--entanglement tradeoff. A claim of the “first quantum treatment of linear checks” would be false. |
| Memory cost of combs and bounded-memory testers [@bisio2012memory; @ohst2026memory; @zonnios2026bounded] | Global coherent-memory cost with free classical memory, and testers retaining a full classical record plus bounded coherent memory, are established frameworks. | Their tasks are protocol realization or process discrimination. None of the inspected sources uses a delayed syndrome audit versus all-output entanglement recovery, or proves the column-order endpoint. |
| Sequential and nondestructive discrimination [@bergou2013sequential; @bilash2024nondestructive; @lim2025local] | Sequential observers extracting information while leaving a state for a later observer are established. Bilash *et al.* prove a no-resource \(1/4\) bound for nondestructive Bell discrimination; for \(K\) equiprobable maximally entangled states, Lim--Hhan--Kwon derive a tight local information--disturbance relation whose perfect-guessing fidelity is \(1/K\), hence \(1/4\) for \(K=4\). | These restrictions are sequential-observer or spatial separable/LOCC processing, and the instrument must discriminate and preserve in the same task. The present restriction is a temporal one-qubit bond followed by an otherwise joint RETURN decoder. The grouped protocol's achievable \((P_{\rm A},F_{\rm R})=(1,1/2)\) rules out a direct identification with the four-state local ceiling. **The numerical value \(1/4\) is not itself novel.** |
| Quantum branching programs and ordered quantum decision diagrams [@ablayev2005branching; @khadiev2022reordering] | Width is quantum memory and variable order can change width dramatically. | Their input variables are classical controls; there are no incoming systems entangled with references and consequently no RETURN fidelity to protect. |
| Post-measurement information [@ballester2008postmeasurement] | Measurements selected before later classical task information are a mature discrimination model. | The later information changes how a classical guess is decoded; it is not a choice between reading a destructive transcript and physically reversing a streamed interaction. |

## The closest numerical collision does not furnish a reduction

The Bilash *et al.* and Lim--Hhan--Kwon results are the most dangerous
superficial collisions because both produce a no-resource value \(1/4\) in
nondestructive discrimination of four maximally entangled alternatives
[@bilash2024nondestructive; @lim2025local]. A direct reduction would need
to map every legal temporal-bond strategy, including its joint terminal RETURN
recovery, into their separable/LOCC nondestructive discriminator while preserving
both scores. That cannot hold as stated: the legal grouped strategy has perfect
audit and return fidelity \(1/2\), whereas their four-state, no-resource bound
would cap the corresponding nondestructive fidelity at \(1/4\). The coincidence
should be cited and explained, but it does not kill the theorem.

The closest structural collision is instead Li--Economou--Barnes. They prove
that the minimum persistent coherent resource in sequential graph-state
generation is the maximum cut entanglement and that changing the output order can
change this resource drastically [@li2022emitters]. This kills any broad framing
such as “quantum memory has never before been linked to temporal order or a cut
width.” What remains different here is that a cut obstruction becomes an exact
penalty in *recoverable input entanglement conditioned on a later operational
choice*.

Piveteau--Renes is the closest code-theoretic collision. Their 2025 preprint
turns efficient classical trellises into efficient optimal quantum decoders and
explicitly uses coherent decoding/uncomputation [@piveteau2025belief]. It kills
claims of a first quantum trellis decoder, first coherent linear-code processor,
or first trellis/interference connection. It does not return unknown channel
inputs or compare coordinate orders under a fixed coherent bond.

## Priority-safe and unsafe wording

A defensible provisional claim is:

> Within the primary literature located through 20 August 2026, this appears to
> be the first exact separation in a late-choice syndrome-AUDIT/all-carrier-RETURN
> game produced solely by permuting the temporal coordinates of one rank-two
> linear code under a one-qubit coherent-memory constraint. A full-crossing cut
> forces \(F_{\rm R}\leq1/4\) at perfect audit, while a noncrossing order attains
> \(F_{\rm R}=1/2\).

For the higher-rank extension, the safe sentence is:

> In the same declared interface, \(m\) consecutive full-rank syndrome blocks
> give the perfect-AUDIT temporal law
> \(F_{\rm R}\leq\prod_j\min\{1,d_j/q^r\}\); its uniform-dimension form is
> attained by the repeated-identity family when \(d=q^k\).

The corresponding consequence can be stated without priority language: two
column permutations with the same basis-column multiset have exact
perfect-AUDIT optima \(d/q^r\) and \((d/q^r)^m\), respectively, for
\(d=q^k\) and \(1\leq k<r\).

The qualifiers “appears,” the date, the exact model, and “at perfect audit” are
essential. The following claims are not supportable:

- first demonstration that temporal order changes quantum memory cost;
- first quantum connection to trellises, pathwidth, matroids, or linear codes;
- a new definition of \(\tau(H)\) or a new classical trellis-width formula;
- first coherent or reversible syndrome computation;
- first \(1/4\) information--disturbance or nondestructive-discrimination bound;
- a solved interleaved tradeoff curve away from \(P_{\rm A}=1\);
- a general \(2^{-\tau}\) law beyond what is actually proved.

## Residual novelty risk

The exact-collision risk after this search is **low to moderate**, not zero. The
remaining danger is not an obvious same-statement paper; it is that the endpoint
may be recognized as a short corollary of an existing quantum-comb rank theorem,
a tensor-network bond lemma, or a recoverability inequality once the game is
translated into that language. Bisio *et al.* already warn that memory costs must
be optimized globally, and the sequential-generation literature already turns
cut Schmidt rank into persistent memory. A paper must therefore foreground the
arbitrary-instrument support lemma and the late-choice operational score, and must
explicitly show why an MPS argument restricted to canonical coherent accumulators
would not suffice.

Literature absence does not validate the proof. The arbitrary
non-QND/adaptive converse and the robust extension have undergone independent
adversarial checks, with no unresolved theorem-level defect. Before a formal
submission, the model should also be sent to experts in quantum combs,
BPQM/code trellises, and nondestructive discrimination with a request for
counterexamples or reductions.

## Search record and limitations

The search covered publisher pages and arXiv primary records through the cutoff,
using combinations of: bounded coherent memory, quantum comb memory cost,
sequential nondemolition discrimination, EPR/entanglement recovery, syndrome and
parity measurement, trellis/pathwidth and coordinate order, quantum belief
propagation, quantum branching programs/OBDD reordering, hybrid quantum-classical
automata, sequential graph-state generation, and late/post-measurement choice.
Exact-phrase and conjunction searches for a syndrome audit versus EPR return, a
same-code column permutation, and a one-qubit streamed recovery game produced no
matching primary source. Search-engine indexing is incomplete, terminology may
differ, and unpublished work or results embedded in appendices can still collide.

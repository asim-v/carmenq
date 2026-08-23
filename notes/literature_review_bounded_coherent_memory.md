# Focused Literature Review: Causal Audit--Return Frontiers with Bounded Coherent Memory

**Cutoff:** 21 August 2026<br>
**Document type:** focused scoping review and adversarial novelty audit; not a PRISMA systematic review<br>
**Language:** English<br>
**Evidence base:** primary papers, publisher records, and first-party preprint archives<br>
**Bibliography:** `references/library.bib`<br>
**Project author:** Javier Emilio Bazán Sánchez, Facultad de Ciencias, Universidad Nacional Autónoma de México, `bazan@ciencias.unam.mx`

---

## 1. Executive verdict

The proposed pivot is scientifically coherent, but most of its vocabulary and most of its ingredients are already occupied. Quantum combs already assign a memory cost to a multi-step protocol. Process tensors and matrix-product representations already identify a temporal bond dimension with an effective quantum memory. Bounded-memory channel discrimination already has dimension-indexed optimization hierarchies. Random-access codes and post-measurement-information games already ask what can be learned from a bounded quantum register after a late query. Information--disturbance theory already relates classical information gained by an instrument to entanglement that can subsequently be recovered. Nondestructive discrimination already connects stabilizer-like labels, recovery fidelity, and an entanglement resource. Classical coding theory already connects a parity-check matrix to weight enumerators and minimal trellis state spaces.

Consequently, none of the following is a defensible novelty claim: introducing a coherent-memory dimension $d$; accumulating a linear syndrome in $k$ qubits; obtaining the perfect-return endpoint $P_{\mathrm A}\leq d/2^k$; expressing an independent-noise calculation through the row-code weight enumerator of $H$; constructing a numerical semidefinite hierarchy; or choosing late between reading information and testing recovery. Each is known, standard, or a direct specialization of known machinery.

One narrower problem survived the search:

> For one precisely specified streaming device class, determine the exact support function of a post-commitment game that asks either for a complete linear syndrome or for transcript-conditioned restoration of every streamed EPR pair, when the device may retain unlimited classical state but only a $d$-dimensional coherent system across temporal cuts.

Writing $H\in\mathbb F_2^{k\times n}$ for a full-rank check matrix, the target quantity is

\[
\beta_{H,d}(\lambda)
=
\sup_{\mathcal S_d}
\left\{
\lambda P_{\mathrm A}(\mathcal S_d)
+(1-\lambda)F_{\mathrm R}(\mathcal S_d)
\right\},
\qquad 0\leq\lambda\leq1,
\]

where the supremum covers arbitrary adaptive, non-QND slot instruments satisfying the declared coherent-bond constraint. The search found no primary source stating this full optimization or an exact intermediate-$d$ solution for it. That absence is not evidence of priority. The problem lies close enough to recent work by Ohst *et al.*, Hsieh *et al.*, Lim *et al.*, and Zonnios--Binder that a result is publishable only if it is more than an application of their frameworks.

The derivation completed after this search changes that conditional verdict. A
full interior frontier was not obtained, but an exact order-sensitive endpoint
and a quantitative robust neighbourhood were proved over the full declared
causal class: two temporal permutations of the same rank-two code have
perfect-AUDIT return optima (1/2) and (1/4), respectively. The endpoint proof
then generalized to a finite-field temporal product law. If \(m\) consecutive
blocks each have full syndrome rank \(r\), perfect AUDIT forces
\(F_{\mathrm R}\leq\prod_j\min\{1,d_j/q^r\}\), where \(d_j\) is the coherent
dimension crossing boundary \(j\); its uniform-dimension form is tight on
repeated identity blocks. These results clear the minimum bar for a focused
theorem paper, while the surrounding architecture remains prior art. A generic
continuity bound, monotonicity statement, numerical seesaw, or SDP formulation
alone would still not have done so.

## 2. The candidate task, frozen before the search

Novelty cannot be assessed against a moving protocol. The following interface is the object searched in this review.

At slot $i\in\{1,\ldots,n\}$, the verifier prepares

\[
|\Phi^+\rangle_{R_iA_i}
=\frac{|00\rangle+|11\rangle}{\sqrt2}
\]

and supplies only $A_i$ to the device. The device applies an arbitrary instrument to $A_i$, its persistent coherent memory $M_{i-1}$, and fresh local ancillas. It emits a carrier $B_i$, which the verifier immediately sequesters, and updates a persistent memory $M_i$ satisfying

\[
\dim M_i\leq d
\]

at every temporal cut. The device may retain an arbitrarily large finite classical transcript $C_i$, including adaptive outcomes and private randomness. No other coherent system, purification, or entangled common cause may cross a cut. Every fresh ancilla must be returned, discarded into a declared inaccessible sink that is charged as coherent memory, or irreversibly dephased into the classical transcript. This accounting rule is essential: without it, $d$ is not a physical constraint.

Only after the final carrier has been sequestered does the verifier sample the branch.

In **AUDIT**, the verifier measures $R_1\cdots R_n$ in the computational basis, obtaining a uniform string $X\in\mathbb F_2^n$. Without access to the sequestered carriers, the device must output the complete rank-$k$ syndrome

\[
S=HX\in\mathbb F_2^k.
\]

The score $P_{\mathrm A}$ is the unconditional probability of returning all $k$ bits correctly.

In **RETURN**, the verifier makes the carriers $B_1\cdots B_n$ available to a declared transcript-conditioned decoder. The decoder may act jointly on those carriers and $M_n$, but it receives no reference system $R_i$. The verifier tests restoration of

\[
|\Phi^+\rangle_{R_1\widehat B_1}\otimes\cdots\otimes
|\Phi^+\rangle_{R_n\widehat B_n}
\]

together with reset of the charged persistent memory. The return score $F_{\mathrm R}$ is unconditional entanglement fidelity; heralding and postselection are included in the score rather than hidden in a conditional fidelity.

This is a trusted-interface resource test. It does not identify a unique internal circuit. A violation of a $d$-dimensional bound excludes the declared model, but an unaccounted coherent side channel, prior entanglement across slots, leakage from the sequestration bank, or decoder access to the references could produce the same violation.

## 3. Search method and evidentiary standard

The search covered combinations of *bounded coherent memory*, *quantum comb memory cost*, *process-tensor bond dimension*, *classically adaptive tester*, *post-measurement information*, *random access code*, *syndrome*, *linear sketch*, *nondestructive discrimination*, *information gain and recovery*, *entanglement fidelity*, *stabilizer*, *trellis complexity*, and *quantum branching program*. Searches were run against arXiv, APS journals, Quantum, IEEE publication records, IOP, Nature-family journals, and the reference lists of the closest papers. The cutoff is 21 August 2026.

Peer-reviewed results and preprints are distinguished below. The search is a focused novelty audit, not a database-complete systematic review. It did not include a formal Scopus or Web of Science export, citation-count screening, dual independent reviewers, or author contact. “No exact collision located” therefore means only that no collision appeared in the searched primary corpus. It must never be rewritten as “the first” without a further professional priority search.

## 4. What the nearest literatures already establish

### 4.1 Multi-time protocols already have memory costs and dimension hierarchies

Bisio *et al.* define the memory cost of a quantum strategy globally, including implementations assisted by free classical memory. Their central warning is directly relevant: memory cannot in general be minimized independently at each step, because compatibility across cuts matters [@bisio2012memory]. Thus the very act of optimizing a comb under a bound $\dim M_i\leq d$ is established formal territory.

Taranto *et al.* characterize inequivalent forms of classical memory in multi-time processes [@taranto2024hierarchy]. Vieira *et al.* construct semidefinite hierarchies that lower-bound an effective temporal environment dimension from observed correlations [@vieira2024dimension]. Roy *et al.* use sequential random-access codes to witness quantum memory in a non-Markovian process [@roy2024semidevice]. These papers prevent any broad claim that a temporal task, a dimension witness, or a separation between classical and quantum memory is new.

Ohst *et al.* are the closest peer-reviewed methodological collision. They formulate limited-memory channel discrimination as constrained separability, include adaptive multi-use protocols with classical memory, and obtain both hierarchy relations and explicit dimension-dependent examples [@ohst2026memory]. Their clock--shift ensemble gives an exact success law of the form $\min\{1,d_E/d\}$, demonstrating that linear dimension ratios in discrimination are already standard. A proposed audit--return paper must therefore do more than encode its strategy set as a constrained-separability problem or report numerically separated values.

Zonnios and Binder, in a June 2026 preprint, define recurrent process testers that retain the full classical outcome record and a coherent memory of prescribed dimension. Their hierarchy is monotone in that dimension and complete at finite time once the memory is sufficiently large [@zonnios2026bounded]. It does not contain the present syndrome-versus-return game, but it occupies the headline “distinguishing temporal quantum resources with bounded coherent memory.”

The tensor-network language is equally established. Sequential interaction with a $d$-dimensional ancilla generates matrix-product states of bond dimension at most $d$ [@schoen2005sequential]. Process tensors admit matrix-product representations whose bond dimension quantifies effective memory [@guomodi2020tensor; @guo2022memorycomplexity]. Bond-dimension witnesses are known independently of process theory [@navascues2018bond]. Calling $d$ a temporal bond dimension is useful notation, not a contribution. The proposed resource must nevertheless be defined as the Hilbert-space dimension of the coherent system crossing a cut: a purified MPS bond $d$ and a density-operator MPO bond that can scale as $d^2$ are not interchangeable.

### 4.2 Guessing versus recovery is an established information--disturbance problem

Banaszek's fidelity balance and Barnum's analysis of least-disturbing square-root dynamics already cover the basic optimization pattern “gain classical information while preserving or recovering quantum information” [@banaszek2001fidelity; @barnum2002information]. Berta, Coles, and Wehner give an exact operational equality between a guessing probability and recoverable entanglement fidelity in a complementary-measurement setting [@berta2014guessing]. Puzzuoli and Watrous show that ancillary dimension can be essential for optimal channel discrimination and connect entanglement preservation to reversibility [@puzzuoli2017ancilla]. These results make a generic $P_{\mathrm A}$--$F_{\mathrm R}$ inequality unsurprising.

Khandelwal and Tavakoli close another tempting novelty route.  They give a
complete characterization of projective-instrument simulability for qubit
inputs and exhibit a strict nonprojective advantage in a one-shot
information--disturbance trade-off [@khandelwal2025instruments].  Therefore,
neither “nonprojective instruments improve information--disturbance” nor “an
SDP distinguishes projective from nonprojective qubit instruments” is new.
Their task averages hemisphere discrimination and output-state fidelity for
one travelling qubit.  It has no sequestered-carrier EPR return, ordered
parity check, or coherent-dimension constraint across several temporal cuts.
The present residual question is consequently about the causal composition
of instruments, not the existence of a nonprojective advantage by itself.

Hsieh *et al.* sharpen the collision. Their 2026 resource theory of interactive quantum instruments assigns a robustness to the ability of an instrument to produce a classical outcome while retaining a nontrivial quantum output. The same robustness has exact operational interpretations in maximally-entangled-state preservation, average state preservation, and recovery of classical information generated by measuring half of a maximally entangled state [@hsieh2026interactive]. This paper makes “classical flag plus entanglement return” a standard instrument resource.

The overlap is exact on the RETURN coordinate. Collapse the committed prefix into a flagged instrument $\mathcal E=\{\mathcal E_c:A^n\to B^nM\}$, set $D=2^n$, and absorb any declared memory reset into the decoder. The optimized RETURN score is then Hsieh *et al.*'s maximally-entangled recovery fidelity, for which their Result 2 gives

\[
1+R(\mathcal E)=D^2F_{\mathrm R}(\mathcal E).
\]

This identity rewrites the candidate objective as a syndrome-guessing term plus a robustness term, but it does not solve the optimization: Hsieh *et al.* impose neither the streamed $d$-dimensional realization nor AUDIT's restriction to the terminal memory and transcript while the carriers remain sequestered. The possible contribution is therefore a constrained joint frontier, not a new recovery measure.

Lim, Hhan, and Kwon study nondestructive local discrimination of entangled states. For $K$ equiprobable maximally entangled states, their information--disturbance relation yields a guessing-plus-fidelity ceiling saturated by random guessing without an entanglement resource; preshared entanglement enables perfect nondestructive discrimination. They also give adaptive stabilizer-based protocols and entanglement-cost bounds [@lim2025local]. The spatial LOCC restriction is not the temporal bond restriction used here, but the resemblance is structural rather than cosmetic. A reduction between the two tasks could eliminate much of the claimed novelty, and it must be tested explicitly.

### 4.3 Late information and bounded storage already support functional queries

König, Maurer, and Renner compare classical and quantum storage of a random string when a predicate is selected later [@konig2005power]. Ballester, Wehner, and Winter allow unlimited classical information together with a bounded quantum register, reveal side information after the memory bound applies, and optimize the later computation of a function $f(X)$ [@ballester2008postmeasurement]. In their two-basis construction, one qubit can suffice to compute any Boolean function, a warning against informal claims that the number of possible queries alone lower-bounds memory.

Shah proves the exact dimension-only message-discrimination lemma needed by the AUDIT coordinate: for arbitrary priors, a $d$-dimensional quantum encoding cannot achieve a guessing probability larger than the sum of the $d$ largest prior masses, and a classical $d$-level message attains the bound [@shah2025qudits]. Thus that step is explicit recent prior art rather than an original memory theorem.

Doriguello and Montanaro define random-access codes for Boolean functions and connect achievable performance to Fourier analysis and noise stability [@doriguello2021boolean]. Roy *et al.* then use a sequential RAC operationally as a temporal quantum-memory witness [@roy2024semidevice]. If the proposed AUDIT branch were changed from returning the complete syndrome to answering one syndrome bit selected after commitment, it would move even closer to this occupied RAC literature. The complete-syndrome version is therefore the cleaner originality target.

Mohan, Tavakoli, and Brunner derive a tight two-score frontier for sequential
\(2\to1\) QRACs and use equality to self-test a qubit instrument
[@mohan2019sequential].  Their polar-decomposition, extremal-instrument, and
Bloch-sphere proof technology is the closest clear template for the new
two-parameter interleaved construction below.  It is not an exact reduction:
their downstream score is a second random-access decoding probability on one
travelling qubit, not transcript-conditioned recovery of four sequestered EPR
carriers after a common late-choice prefix.

Late choice itself adds no novelty. Delayed-choice wave--particle games, post-measurement information, quantum seals, and audit-versus-recovery formulations already make the decision after a common prefix [@bagan2018duality; @kimmel2019seals]. The contribution, if any, must lie in the exact bounded-memory support for the same causal device.

### 4.4 Linear syndromes bring established coding and streaming structure

At the classical level, $Hx$ is an $\mathbb F_2$-linear sketch. Linear sketches are a standard small-space streaming and communication primitive, with a mature theory of when they optimally represent Boolean functions [@kannan2018linear]. The exact syndrome accumulator therefore cannot serve as the interdisciplinary novelty claim.

For an independent binary error $Z$ with correlation parameter $t\in[0,1]$, a symmetric bitwise record gives

\[
\Pr[HZ=0]
=2^{-k}\sum_{u\in\mathbb F_2^k}
t^{\operatorname{wt}(H^{\mathsf T}u)}.
\]

The sum is the weight enumerator of the row code of $H$, evaluated at a particular point. This is a direct character/Fourier calculation and belongs to the MacWilliams weight-enumerator tradition [@macwilliams1963weights]. It is a useful reduction for constructing examples, but it cannot be advertised as a new connection between quantum memory and coding theory.

The order dependence of sequential code processing also has a classical antecedent. Minimal conventional trellises realize a linear block code with a time-dependent state space, their state-complexity profile depends on the coordinate ordering, and optimizing consecutive trellis sectionalizations is established [@forney1994trellis; @mceliece1996bcjr; @kashyap2008pathwidth; @lafourcade1996sectionalization]. Quantum and classical branching programs likewise use width as a memory measure [@ablayev2005branching]. A theorem that merely identifies a syndrome accumulator, a consecutive column partition, or its greedy computation with a code trellis or branching program would therefore be expository.

For the present terminal task there is an additional negative result: if the device must output the entire syndrome, two classical prefixes are equivalent for every common suffix precisely when their partial syndromes agree. The deterministic width after cut $i$ is therefore $2^{\operatorname{rank}(H_{\leq i})}$, reaching $2^k$ at the end for full-rank $H$. Ordinary trellis merging does not compress exact terminal syndrome output. Trellis or pathwidth could become central only after changing the task to intermediate/random-cut audits, a decision problem such as $Hx=0$, or partial row queries.

The open opportunity is subtler: the audit score is a decoding objective, the return score charges irreversible information extraction, and a coherent bond of dimension $d$ can carry a superposition of partial-syndrome states. The post-review derivation recorded in Section 8 establishes that the attainable region can depend on the ordered columns of $H$, even when rank and abstract code are unchanged. The result is an endpoint and robust-neighbourhood theorem rather than a solution of the complete interior frontier.

### 4.5 The 2026 frontier is moving quickly

Three July 2026 preprints substantially raise the priority risk. Arunachalam and Schatzki prove tight sample-complexity tradeoffs for stabilizer testing and learning when only $q$ coherent qubits may persist between sequential measurements [@arunachalam2026stabilizer]. Bravo-Prieto, Gong, and Mele prove a query-complexity advantage for coherent-memory process tomography even against adaptive protocols with unlimited classical computation and fresh ancillas [@bravoprieto2026tomography]. Sundar and Elliott identify an adaptive agent's memory with a temporal matrix-product bond and derive a fidelity certificate for truncating that dimension [@sundar2026reduction]. None states the present audit--return frontier, but together they occupy the general message that bounded coherent memory is an operational resource with certifiable performance consequences.

These works are preprints at the cutoff and must be labelled as such. They are nevertheless close enough that any submission should be re-searched immediately before posting.

### 4.6 The terminal-readout geometry is prior art; its temporal use is the candidate contribution

The final interior enclosure uses standard minimum-error discrimination more
deeply than the earlier draft acknowledged. Bae and Hwang formulate qubit
discrimination through the Helstrom/KKT dual and complementary states
[@bae2013qubitdiscrimination]. Weir, Barnett, and Croke make explicit that an
optimal nontrivial qubit readout consists of weighted rank-one projectors and
that an optimal readout with at most four outcomes always exists
[@weir2017qubitdiscrimination]. Rouhbakhsh and Ghoreishi subsequently give a
constructive Bloch-vector treatment that includes arbitrary priors and genuine
three- and four-outcome solutions [@rouhbakhsh2023qubitdiscrimination]. The
underlying extremal-POVM decomposition and rank restrictions are older still
[@darianolo2005extremal].

Consequently, neither the active-set classification nor reconstruction of a
syndrome state from a Helstrom operator should be advertised as new quantum
state-discrimination theory. The literature search through 23 August 2026 did
not locate the geometry-free weighted projective comparison used here, but it
is an elementary consequence of those standard KKT identities and is not a
safe standalone priority claim. The defensible contribution is narrower: the
comparison becomes a certified cut linking every genuinely multi-outcome
terminal readout to independent support lines while preserving the temporal
RETURN score. Together with the common-instrument constraints, probability-
cone cover, and small-effect deletion argument, it exhausts all terminal
arities for one bounded-memory support direction. The originality resides in
that constrained temporal synthesis and enclosure, if it survives review, not
in Helstrom geometry itself.

## 5. Claim--collision matrix

| Primary source | What is already established | Clause absent from that source | Consequence for the proposed claim |
|---|---|---|---|
| Bisio *et al.* 2012 [@bisio2012memory] | Global quantum memory cost of combs with classical assistance | No syndrome audit or recovery score | Memory cost itself is occupied |
| Taranto *et al.* 2024 [@taranto2024hierarchy] | Structural hierarchy of classical multi-time memory | No bounded coherent dimension frontier | “Classical versus quantum temporal memory” is too broad |
| Vieira *et al.* 2024 [@vieira2024dimension] | SDP dimension witnesses for temporal environments | No same-prefix entanglement return | A numerical dimension witness is insufficient |
| Ohst *et al.* 2026 [@ohst2026memory] | Limited-memory discrimination, classical adaptation, constrained-separability hierarchy | Different task and no audit--return support | Generic SDP/hierarchy contribution is occupied |
| Zonnios--Binder 2026, preprint [@zonnios2026bounded] | Complete finite-time hierarchy indexed by coherent memory dimension | Autonomous discrimination rather than syndrome/recovery | The dimension-indexed tester architecture is occupied |
| Banaszek 2001; Barnum 2002 [@banaszek2001fidelity; @barnum2002information] | Optimal information--disturbance and square-root dynamics | No causal multi-slot dimension constraint | One-slot tradeoff machinery is prior art |
| Khandelwal--Tavakoli 2025 [@khandelwal2025instruments] | Complete qubit projective-instrument simulability test and a nonprojective information--disturbance advantage | One-shot hemisphere task; no ordered code, bounded temporal bond, or all-carrier EPR return | Nonprojectivity and its SDP characterization cannot be the headline |
| Berta--Coles--Wehner 2014 [@berta2014guessing] | Exact guessing/recoverable-entanglement relation | Complementary measurement task, no streaming bond | Guessing and return are not a new pairing |
| Hsieh *et al.* 2026 [@hsieh2026interactive] | Instrument resource with exact classical-recovery and entanglement-preservation meanings | No bounded temporal bond or syndrome family | “Interactive instrument” framing is occupied |
| Lim--Hhan--Kwon 2025 [@lim2025local] | Nondestructive MES discrimination, tight tradeoff, stabilizer protocol, entanglement cost | Spatial LOCC resource instead of temporal memory | A possible reduction is a hard novelty risk |
| König *et al.* 2005; Ballester *et al.* 2008 [@konig2005power; @ballester2008postmeasurement] | Late functional query with bounded quantum memory and free classical data | No obligation to restore EPR carriers | Late-query bounded storage is occupied |
| Doriguello--Montanaro 2021 [@doriguello2021boolean] | Boolean-function RACs and noise-stability analysis | No return branch or causal comb | Boolean/noise-stability generalization is occupied |
| Mohan--Tavakoli--Brunner 2019 [@mohan2019sequential] | Tight sequential-QRAC score frontier and qubit-instrument self-test | Second QRAC score rather than all-carrier EPR return; no four-slot local-completion constraint | Sequential weak-instrument/Bloch machinery is occupied, but does not solve the present support |
| MacWilliams 1963; McEliece 1996 [@macwilliams1963weights; @mceliece1996bcjr] | Weight enumerators and minimal trellis state complexity | No quantum recovery objective | Code enumerators and trellis memory are ingredients, not novelty |
| Schön *et al.* 2005; Guo *et al.* 2020 [@schoen2005sequential; @guomodi2020tensor] | Sequential ancilla dimension equals an MPS/process bond resource | No audit score | The tensor-network connection is occupied |
| Arunachalam--Schatzki 2026, preprint [@arunachalam2026stabilizer] | Tight limited-coherent-memory separations for sequential stabilizer tasks | Multiple-copy testing/learning, no return | Stabilizer plus bounded memory is not a headline |

No row in the table contains all of the candidate's clauses. Several rows, however, provide machinery that might imply its easy bounds. That is why the residual claim is narrow.

## 6. Subclaims that should be treated as known or routine

### 6.1 Perfect-return endpoint

Let $K=2^k$. At $F_{\mathrm R}=1$, all classical side information that survives outside the charged recovery system must be independent of the reference label. The AUDIT information can therefore reside only in a $d$-dimensional quantum register. For $K$ equiprobable labels, the elementary dimension bound gives

\[
P_{\mathrm A}\leq \min\left\{1,\frac dK\right\}.
\]

When $d=2^q\leq K$, coherently storing $q$ independent syndrome bits and guessing the other $k-q$ bits attains $d/K$; a $K$-dimensional coherent syndrome accumulator attains $(P_{\mathrm A},F_{\mathrm R})=(1,1)$. For arbitrary non-power-of-two $d$, the dimension upper bound need not be causally tight. In every case this endpoint follows from standard state-discrimination, conditional recovery, and min-entropy arguments [@gregoratti2003lostfound; @konig2009operational]. It is a required lemma, not the paper's main result.

### 6.2 Full-memory construction

For columns $h_i\in\mathbb F_2^k$ of $H$, a $k$-qubit register can update

\[
|m\rangle_M|x_i\rangle_{A_i}
\longmapsto
|m\oplus h_i x_i\rangle_M|x_i\rangle_{B_i}.
\]

After all slots, measuring $M$ reveals $HX$. If RETURN is selected, applying the updates in reverse with the sequestered carriers resets $M$ and restores the EPR pairs. This is reversible linear computation and stabilizer-syndrome extraction [@bennett1973logical]. It proves completeness but is not original.

Three baseline points should be used to reject incorrect candidate frontiers. Doing nothing gives $(P_{\mathrm A},F_{\mathrm R})=(1/K,1)$. Coherently retaining $q$ independent syndrome bits in $d=2^q$ dimensions gives $(2^q/K,1)$. Projectively measuring and classically recording the full syndrome gives $(1,1/K)$, because the measurement destroys coherence between the $K$ syndrome sectors. Any proposed upper bound below one of these points is false.

### 6.3 Independent weak-record curve

If every slot is subjected to the same binary weak record with classical correlation $t$, the audit success is the row-code weight-enumerator expression above. If the least-disturbing local instrument has return factor $f(t)$, the product strategy has

\[
F_{\mathrm R}=f(t)^n.
\]

This gives a useful lower bound for $d=1$, and in special symmetric cases it may be optimal. The weight enumerator, square-root/Lüders update, and product law are known pieces. Global optimality over adaptive non-QND instruments is the only theorem-level content such a curve could carry.

### 6.4 Generic robust dimension bound

A statement of the schematic form

\[
P_{\mathrm A}
\leq
\frac d{2^k}+c\sqrt{1-F_{\mathrm R}}
\]

is likely obtainable by combining approximate complementary-channel decoupling, a fidelity--trace-distance inequality, and the $d/K$ discrimination bound [@kretschmann2008information]. Unless its constant is optimal and it is tight for a meaningful family, this is not enough for a paper. The square-root scaling is already the natural continuity scale of information--disturbance theory.

### 6.5 SDP formulation

Finite $n$ and $d$ permit a rank- or Schmidt-number-constrained comb formulation, outer relaxations, and seesaw lower bounds. Ohst *et al.* already provide the relevant constrained-separability methodology, while bounded-bond tensor-network optimization is widespread. Fixed-dimension temporal strategy sets can also be nonconvex [@mao2022dimension], so classical flags and shared randomization must be represented explicitly rather than silently convexifying the model. A solver is valuable infrastructure, but “the problem can be written as an SDP hierarchy” is not the intended contribution.

## 7. The residual research gap

The strongest defensible target is now:

> **Exact causal syndrome audit--return theorem.** For a declared ordered full-rank matrix $H$, coherent bond dimension $d$, and branch weight $\lambda$, derive $\beta_{H,d}(\lambda)$ over every adaptive finite-outcome streaming instrument allowed by the interface, including non-QND operations and transcript-conditioned joint recovery. Characterize optimal strategies and equality cases, and prove at least one strict intermediate-dimension separation that is robust to experimental error.

Five features must remain joined. The same common-prefix comb must face both branches. The output carriers must be inaccessible during AUDIT. Every coherent system crossing a temporal cut must count toward $d$. The return test must include failures and residual memory rather than postselecting. Finally, the optimization must cover the full declared strategy class rather than an iid, diagonal, Clifford, covariant, or stabilizer ansatz.

The most valuable result would show that matrices with the same rank $k$ but different ordered code structure have different intermediate-$d$ frontiers. That would establish that the theorem is not merely the endpoint $d/2^k$ in disguise. It would also make the coding connection substantive: the row-code weight distribution would control one classical slice, while a coherent temporal bond would introduce a genuinely different optimization related to, but not identical with, trellis state complexity.

The operational payoff would be a calibrated dimension witness. If an experiment attains

\[
\lambda\widehat P_{\mathrm A}+(1-\lambda)\widehat F_{\mathrm R}
>
\beta_{H,d}(\lambda)+\Delta_{\rm stat}+\Delta_{\rm sys},
\]

then, under the audited interface, no realization with coherent bond dimension at most $d$ is compatible with the data. This is more informative than the current binary statement “some coherent cross-slot resource was present.”

## 8. Post-review theorem resolution

The proposed triangle gate was solved exactly, but it failed the stronger originality test. For

\[
H_\triangle=
\begin{pmatrix}
1&1&0\\
0&1&1
\end{pmatrix},
\qquad d=2,
\]

the support function reduces to the static four-label/qubit information--recovery curve. The proof is valid for arbitrary adaptive non-QND instruments, but it loses the ordered structure of $H$. It is therefore a useful baseline rather than the flagship contribution.

The first genuinely order-sensitive instance compares

\[
H_{\rm G}=\begin{pmatrix}1&1&0&0\\0&0&1&1\end{pmatrix},
\qquad
H_{\rm I}=\begin{pmatrix}1&0&1&0\\0&1&0&1\end{pmatrix}.
\]

These matrices differ only by a coordinate permutation. With one persistent coherent qubit, an unrestricted genuinely classical transcript, immediate sequestration of emitted carriers, and the same late AUDIT/RETURN choice, the grouped order attains the complete static boundary. In particular, at perfect AUDIT it attains $F_{\mathrm R}=1/2$. For the interleaved order, an arbitrary-instrument converse proves

\[
P_{\mathrm A}=1
\quad\Longrightarrow\quad
F_{\mathrm R}\leq\frac14,
\]

and the bound is attained. The proof refines each transcript to a single Kraus leaf, uses the full-crossing cut to limit every perfect-audit leaf to one quarter of the computational words, and then applies flagged polar recovery and pinching. It makes no QND, Clifford, covariance, or finite-transcript-alphabet assumption.

The endpoint now has a linear-tail robust version. For \(m\) consecutive
full-rank blocks, if \(t_c\) is the spectral mass of a refined leaf below the
perfect-AUDIT rank and \(D\) is the input dimension, then

\[
\sum_ct_c\leq mD(1-P_{\mathrm A}).
\]

For the four-slot interleaved stream this yields

\[
F_{\mathrm R}\leq
\frac14+\frac12\theta
+\frac{\sqrt3}{2}\sqrt{\theta(1-\theta)},
\qquad
\theta=\min\{2(1-P_{\mathrm A}),3/4\}.
\]

The square-root order is sharp. Maximising the bound at fixed audit weight
gives

\[
\beta^{\rm stream}_{H_{\rm I},2}(\lambda)
\leq\frac12+\frac\lambda4
+\frac14\sqrt{7\lambda^2-10\lambda+4},
\]

which certifies a strict gap below the grouped/static support throughout
\(3/7<\lambda<1\), including balanced weight.

The exact interior frontier remains open.  For \(q,v\in[0,1]\), a complete
four-slot QND instrument achieves

\[
P_{\rm can}=\frac12+qv\sqrt{1-v^2}-q(1-q)v^2,
\]

\[
F_{\rm can}=\frac14\left[
\sqrt{1-(1-q^2)v^2}
+v\bigl(1-q+2\sqrt{q(1-q)}\bigr)
\right]^2.
\]

Optimizing this explicit family gives a first-order onset within that family at
\(\lambda=0.477812793357157\ldots\) and the balanced lower bound
\(0.755437446228747\ldots\).  Unrestricted binary-tree searches and
three-/four-outcome QND controls found no improvement, but those ansatzes did
not cover arbitrary finite-outcome non-QND instruments.  A subsequent
complete ternary-outcome search produced

\[
(P_{\rm A},F_{\rm R})
=(0.625754561820\ldots,0.893143378814\ldots)
\]

at balanced weight, with support
\(0.759448970317\ldots\).  The stored instrument independently verifies a
strict excess of \(0.004011524089\ldots\), so the two-parameter frontier
conjecture is false.  This correction does not alter either proved endpoint
theorem.  It instead sharpens the remaining problem: local finite-outcome
completion and non-QND polar geometry cannot be discarded in an interior
converse.

A subsequent local Weyl-completion theorem removes the unbounded adaptive
tree from that remaining problem exactly: the full support equals the maximum
of one Hilbert--Schmidt-normalised bond-two Choi MPS.  At balanced weight, a
three-effect physical leaf reaches
\(0.759802783851444\ldots\).  That curve is itself not the complete attainable
interior.  A symmetric four-effect physical leaf reaches
\(0.765898815264694\ldots\) at \(\lambda=0.6\), exceeding the three-effect
value by \(0.010192880679\ldots\).  Its explicit tensor has diagonal Choi Gram
matrix and exact local Pauli completion.  The revised three-/four-effect
envelope is the strongest known lower bound. A later exhaustive terminal-
readout calculation at \(\lambda=0.6\) gives the complete solver-conditional
enclosure

\[
0.7658988152646944
\leq\beta_{4\mathrm s}(0.6)
\leq\beta_{2\mathrm b}(0.6)
\leq0.76662.
\]

The reduction covers projective, ternary, and four-active readouts rather than
assuming terminal projectivity. Its logical sector exhaustion is exact, but
the upper endpoint still depends on recorded CLARABEL and SCIP tolerances; it
is not an outward interval proof, a closed form, or an equality result for the
explicit four-effect leaf. This distinction is important for the originality
claim: compact variational reduction is proved and one interior direction is
tightly enclosed; the complete support curve is not.

The structural condition is the standard trellis-connectivity quantity

\[
\tau(H)=\max_i\bigl(
\operatorname{rank}H_{\leq i}
+\operatorname{rank}H_{>i}
-\operatorname{rank}H\bigr).
\]

For rank-two checks, $\tau\leq1$ admits a one-qubit handoff construction that reaches the static boundary, whereas $\tau=2$ supplies a full-crossing cut and the perfect-audit $1/4$ upper bound. The invariant and its order dependence are established trellis theory [@forney1994trellis; @mceliece1996bcjr; @kashyap2008pathwidth]. The new object is the late-choice information--recoverability consequence under the declared interface, not $\tau$ itself.

The endpoint mechanism extends beyond binary rank two. Let
\(H\in\mathbb F_q^{r\times n}\) split into \(m\) consecutive nonempty blocks,
each of rank \(r\), and let \(d_j\) be the coherent dimension crossing the
boundary after block \(j\). For arbitrary adaptive non-QND instruments with a
genuinely classical transcript,

\[
P_{\mathrm A}=1
\quad\Longrightarrow\quad
F_{\mathrm R}\leq
\prod_{j=1}^m\min\left\{1,\frac{d_j}{q^r}\right\}.
\]

The proof bounds every refined leaf to at most \(d_j\) cumulative syndrome
labels at boundary \(j\), counts the resulting block-syndrome tuples, and
converts computational support into flagged recovery fidelity. For uniform
\(d=q^k\), repeated identity blocks attain
\(F_{\mathrm R}=(d/q^r)^m\). Thus the result is a temporal product law, not a
new definition of memory or trellis width.

Let \(\mu(H)\) be the maximum number of consecutive full-rank blocks in a
partition of the ordered columns. The earliest-full-rank greedy partition is
optimal, and the uniform-dimension theorem becomes

\[
F_{\mathrm R}\leq
\min\left\{1,\left(\frac{d}{q^r}\right)^{\mu(H)}\right\}.
\]

This descriptor gives \(\mu(H_{\rm G})=1\) and
\(\mu(H_{\rm I})=2\), so the original order separation is the first two values
of the same law. Trellis sectionalization is prior art
[@lafourcade1996sectionalization]; no novelty is assigned to partitioning an
ordered code. The candidate contribution is the exact EPR-return exponent in
the frozen late-choice interface.

This yields an exact asymptotic order separation using the same column
multiset. Put \(m\) copies of each standard column into either a batched order
or the cycled order \([I_r\mid\cdots\mid I_r]\). For \(d=q^k\),
\(1\leq k<r\), both perfect-AUDIT bounds are attainable and

\[
F_{\mathrm R}^{\star}(H_{\rm batched})=\frac d{q^r},
\qquad
F_{\mathrm R}^{\star}(H_{\rm cycled})
=\left(\frac d{q^r}\right)^m.
\]

The ratio \((q^r/d)^{m-1}\) grows exponentially. This is the strongest
resolved contribution: it promotes the four-slot example from an isolated
constant gap to an exact family while keeping the code columns, memory bound,
and late tasks fixed.

## 9. Updated collision assessment

The exact result survives the focused priority search, but only as a narrow conjunction. Sequential quantum-state generation already links emission order to persistent coherent resources [@li2022emitters]. Quantum trellises, coherent message passing, and uncomputation for linear codes are also occupied [@ollivier2006trellises; @piveteau2022message; @piveteau2025belief]. Lim--Hhan--Kwon obtain a numerical value $1/4$ in a different spatial nondestructive-discrimination model [@lim2025local]. None of these sources located through the cutoff contains the same-code column permutation, one-qubit temporal bond, late complete-syndrome AUDIT, all-carrier EPR RETURN, and exact grouped $1/2$ versus interleaved $1/4$ endpoint. Nor did the search locate the finite-field product law for repeated full-rank temporal blocks in this interface.

The same focused update located no primary source combining causal
dimension-to-list conversion, full-rank temporal block intersections, a Ky
Fan tail bound, and all-carrier entanglement recovery. Bounded-quantum-storage
and postmeasurement-information bounds are nearby prior art
[@konig2005power; @ballester2008postmeasurement; @ohst2026memory], but they do not state the
rank-tail or support certificate in this interface. This is not a claim of
absolute priority. An equivalent statement may exist under quantum-comb rank,
recoverability, list-decoding, or tensor-network language. The exact-collision
risk is assessed as low to moderate, and the search should be repeated before
submission.

## 10. Remaining kill criteria and research sequence

The result should not be promoted as a full frontier: the grouped curve, the
interleaved perfect-AUDIT endpoint, and the interval \(3/7<\lambda<1\) are
proved as upper statements over the full strategy class, but the latter upper
bound does not match the best complete lower strategy. A general
\(2^{-\tau}\) law and a theorem for arbitrary partial-rank block sequences
have not been established.

The next mathematical target is the optimal coefficient in the sharp
\(\sqrt{1-P_{\mathrm A}}\) law, a solver-independent closure of the
\(\lambda=0.6\) enclosure, or the exact interleaved support curve. A
legal weak-measurement family proves that a smaller asymptotic order is
impossible. A useful experimental result must also charge every coherent side
channel, verify carrier sequestration, include all RETURN failures, and
conclude only incompatibility with the declared one-qubit streaming model.

## 11. Safe novelty language

The following sentence is supportable at this cutoff:

> Within the primary literature located through 21 August 2026, this appears to be the first exact late-choice syndrome-AUDIT/all-carrier-RETURN separation produced solely by permuting the temporal coordinates of one rank-two linear code under a one-qubit coherent-memory constraint: the noncrossing order attains $F_{\mathrm R}=1/2$ at perfect audit, whereas a full-crossing order is bounded by, and attains, $F_{\mathrm R}=1/4$.

For the higher-rank endpoint, the safe statement is: in the same declared
interface, \(m\) consecutive full-rank syndrome blocks impose
\(F_{\mathrm R}\leq\prod_j\min\{1,d_j/q^r\}\) at perfect AUDIT; the
uniform-dimension law is tight on repeated identity blocks for \(d=q^k\).
For approximate AUDIT, the safe statement is that the summed spectral mass
below the perfect-AUDIT rank is at most \(mD(1-P_{\mathrm A})\), which implies
the explicit recovery bound stated above. Neither coefficient optimality nor
the exact interior support curve is claimed. Separately, at
\(\lambda=0.6\) the two-block rank-two support is enclosed in
\([0.7658988152646944,0.76662]\) by a complete finite sector exhaustion
conditional on the recorded conic and spatial-solver tolerances.
Equivalently, a uniform bound \(d\) gives
\(F_{\mathrm R}\leq(d/q^r)^{\mu(H)}\), capped at one.

The qualifiers “appears,” the date, the exact interface, and “at perfect audit” are essential. Unsafe formulations include “the first order dependence of quantum memory,” “a new trellis invariant,” “the first quantum processing of linear checks,” “the first $1/4$ information--disturbance relation,” and “the complete frontier is solved.”

## 12. Final go/no-go decision

**GO for a focused theorem paper; NO for the original broad framing.** The order-sensitive endpoint, its robust extension, the finite-field temporal product law, and the exact exponential same-columns order separation are mathematically nontrivial and were not reduced by the closest located literature. They establish a real operational consequence of temporal code order and repeated full-rank temporal fragmentation under bounded coherent memory.

The paper should be built around this result, with the triangle solution presented as a negative originality gate, the \(\lambda=0.6\) enclosure presented as a tightly scoped computational theorem, and the full interleaved support curve listed as open. The broader reversible-histories narrative can remain motivation, but it is not the contribution.

---

### Search update policy

Because four of the closest items appeared in 2026, repeat the title/abstract and citation-neighbourhood search immediately before any preprint submission. In particular monitor bounded coherent-memory process discrimination, interactive-instrument resource theories, limited-memory stabilizer learning, nondestructive discrimination with restricted entanglement, and dimension-reduced adaptive agents.

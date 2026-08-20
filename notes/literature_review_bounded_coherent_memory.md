# Focused Literature Review: Causal Audit--Return Frontiers with Bounded Coherent Memory

**Cutoff:** 19 August 2026<br>
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

The present verdict is therefore **amber-green for the exact theorem and red for the surrounding architecture**. A closed analytic frontier for at least one nontrivial intermediate dimension, with a strict dimension separation and equality cases over the full causal strategy class, would be a real contribution. A generic continuity bound, a monotonicity statement, a numerical seesaw, or an SDP formulation alone would not clear the novelty bar.

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

The search covered combinations of *bounded coherent memory*, *quantum comb memory cost*, *process-tensor bond dimension*, *classically adaptive tester*, *post-measurement information*, *random access code*, *syndrome*, *linear sketch*, *nondestructive discrimination*, *information gain and recovery*, *entanglement fidelity*, *stabilizer*, *trellis complexity*, and *quantum branching program*. Searches were run against arXiv, APS journals, Quantum, IEEE publication records, IOP, Nature-family journals, and the reference lists of the closest papers. The cutoff is 19 August 2026.

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

Hsieh *et al.* sharpen the collision. Their 2026 resource theory of interactive quantum instruments assigns a robustness to the ability of an instrument to produce a classical outcome while retaining a nontrivial quantum output. The same robustness has exact operational interpretations in maximally-entangled-state preservation, average state preservation, and recovery of classical information generated by measuring half of a maximally entangled state [@hsieh2026interactive]. This paper makes “classical flag plus entanglement return” a standard instrument resource.

The overlap is exact on the RETURN coordinate. Collapse the committed prefix into a flagged instrument $\mathcal E=\{\mathcal E_c:A^n\to B^nM\}$, set $D=2^n$, and absorb any declared memory reset into the decoder. The optimized RETURN score is then Hsieh *et al.*'s maximally-entangled recovery fidelity, for which their Result 2 gives

\[
1+R(\mathcal E)=D^2F_{\mathrm R}(\mathcal E).
\]

This identity rewrites the candidate objective as a syndrome-guessing term plus a robustness term, but it does not solve the optimization: Hsieh *et al.* impose neither the streamed $d$-dimensional realization nor AUDIT's restriction to the terminal memory and transcript while the carriers remain sequestered. The possible contribution is therefore a constrained joint frontier, not a new recovery measure.

Lim, Hhan, and Kwon study nondestructive local discrimination of entangled states. For $K$ equiprobable maximally entangled states, their information--disturbance relation yields a guessing-plus-fidelity ceiling saturated by random guessing without an entanglement resource; preshared entanglement enables perfect nondestructive discrimination. They also give adaptive stabilizer-based protocols and entanglement-cost bounds [@lim2025local]. The spatial LOCC restriction is not the temporal bond restriction used here, but the resemblance is structural rather than cosmetic. A reduction between the two tasks could eliminate much of the claimed novelty, and it must be tested explicitly.

### 4.3 Late information and bounded storage already support functional queries

König, Maurer, and Renner compare classical and quantum storage of a random string when a predicate is selected later [@konig2005power]. Ballester, Wehner, and Winter allow unlimited classical information together with a bounded quantum register, reveal side information after the memory bound applies, and optimize the later computation of a function $f(X)$ [@ballester2008postmeasurement]. In their two-basis construction, one qubit can suffice to compute any Boolean function, a warning against informal claims that the number of possible queries alone lower-bounds memory.

Doriguello and Montanaro define random-access codes for Boolean functions and connect achievable performance to Fourier analysis and noise stability [@doriguello2021boolean]. Roy *et al.* then use a sequential RAC operationally as a temporal quantum-memory witness [@roy2024semidevice]. If the proposed AUDIT branch were changed from returning the complete syndrome to answering one syndrome bit selected after commitment, it would move even closer to this occupied RAC literature. The complete-syndrome version is therefore the cleaner originality target.

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

The order dependence of sequential code processing also has a classical antecedent. Minimal conventional trellises realize a linear block code with a time-dependent state space, and their state-complexity profile depends on the coordinate ordering [@forney1994trellis; @mceliece1996bcjr; @kashyap2008pathwidth]. Quantum and classical branching programs likewise use width as a memory measure [@ablayev2005branching]. A theorem that merely identifies a syndrome accumulator with a code trellis or branching program would therefore be expository.

For the present terminal task there is an additional negative result: if the device must output the entire syndrome, two classical prefixes are equivalent for every common suffix precisely when their partial syndromes agree. The deterministic width after cut $i$ is therefore $2^{\operatorname{rank}(H_{\leq i})}$, reaching $2^k$ at the end for full-rank $H$. Ordinary trellis merging does not compress exact terminal syndrome output. Trellis or pathwidth could become central only after changing the task to intermediate/random-cut audits, a decision problem such as $Hx=0$, or partial row queries.

The open opportunity is subtler: the audit score is a decoding objective, the return score charges irreversible information extraction, and a coherent bond of dimension $d$ can carry a superposition of partial-syndrome states. The exact joint optimum could depend on the ordered columns of $H$ in a way not captured by the elementary terminal-width argument. That dependence has not been established by the sources located here and must be demonstrated rather than assumed.

### 4.5 The 2026 frontier is moving quickly

Three July 2026 preprints substantially raise the priority risk. Arunachalam and Schatzki prove tight sample-complexity tradeoffs for stabilizer testing and learning when only $q$ coherent qubits may persist between sequential measurements [@arunachalam2026stabilizer]. Bravo-Prieto, Gong, and Mele prove a query-complexity advantage for coherent-memory process tomography even against adaptive protocols with unlimited classical computation and fresh ancillas [@bravoprieto2026tomography]. Sundar and Elliott identify an adaptive agent's memory with a temporal matrix-product bond and derive a fidelity certificate for truncating that dimension [@sundar2026reduction]. None states the present audit--return frontier, but together they occupy the general message that bounded coherent memory is an operational resource with certifiable performance consequences.

These works are preprints at the cutoff and must be labelled as such. They are nevertheless close enough that any submission should be re-searched immediately before posting.

## 5. Claim--collision matrix

| Primary source | What is already established | Clause absent from that source | Consequence for the proposed claim |
|---|---|---|---|
| Bisio *et al.* 2012 [@bisio2012memory] | Global quantum memory cost of combs with classical assistance | No syndrome audit or recovery score | Memory cost itself is occupied |
| Taranto *et al.* 2024 [@taranto2024hierarchy] | Structural hierarchy of classical multi-time memory | No bounded coherent dimension frontier | “Classical versus quantum temporal memory” is too broad |
| Vieira *et al.* 2024 [@vieira2024dimension] | SDP dimension witnesses for temporal environments | No same-prefix entanglement return | A numerical dimension witness is insufficient |
| Ohst *et al.* 2026 [@ohst2026memory] | Limited-memory discrimination, classical adaptation, constrained-separability hierarchy | Different task and no audit--return support | Generic SDP/hierarchy contribution is occupied |
| Zonnios--Binder 2026, preprint [@zonnios2026bounded] | Complete finite-time hierarchy indexed by coherent memory dimension | Autonomous discrimination rather than syndrome/recovery | The dimension-indexed tester architecture is occupied |
| Banaszek 2001; Barnum 2002 [@banaszek2001fidelity; @barnum2002information] | Optimal information--disturbance and square-root dynamics | No causal multi-slot dimension constraint | One-slot tradeoff machinery is prior art |
| Berta--Coles--Wehner 2014 [@berta2014guessing] | Exact guessing/recoverable-entanglement relation | Complementary measurement task, no streaming bond | Guessing and return are not a new pairing |
| Hsieh *et al.* 2026 [@hsieh2026interactive] | Instrument resource with exact classical-recovery and entanglement-preservation meanings | No bounded temporal bond or syndrome family | “Interactive instrument” framing is occupied |
| Lim--Hhan--Kwon 2025 [@lim2025local] | Nondestructive MES discrimination, tight tradeoff, stabilizer protocol, entanglement cost | Spatial LOCC resource instead of temporal memory | A possible reduction is a hard novelty risk |
| König *et al.* 2005; Ballester *et al.* 2008 [@konig2005power; @ballester2008postmeasurement] | Late functional query with bounded quantum memory and free classical data | No obligation to restore EPR carriers | Late-query bounded storage is occupied |
| Doriguello--Montanaro 2021 [@doriguello2021boolean] | Boolean-function RACs and noise-stability analysis | No return branch or causal comb | Boolean/noise-stability generalization is occupied |
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

## 8. Minimum theorem that would change the publication verdict

The smallest useful instance is

\[
n=3,
\qquad
k=2,
\qquad
H_\triangle=
\begin{pmatrix}
1&1&0\\
0&1&1
\end{pmatrix},
\qquad
d\in\{1,2,4\}.
\]

Here $d=4$ gives the full coherent construction, while $d=2$ is a genuinely intermediate resource. A publishable core would provide an exact analytic expression or a finite algebraic characterization of $\beta_{H_\triangle,2}(\lambda)$, prove it against all adaptive non-QND strategies, and exhibit an interval of $\lambda$ on which it is strictly separated from both the $d=1$ region and trivial classical time-sharing with the $d=4$ strategy. A dual certificate or a symmetry reduction with fully justified extremality would be substantially stronger than floating-point optimization.

The next step would compare two rank-two matrices with different ordered column structure. If the frontiers coincide for every such example and depend only on $k$ and $d$, that is also informative: it suggests a representation-theoretic theorem rather than a code-specific one. Either outcome is scientifically cleaner than assuming code dependence from the outset.

A robust corollary should translate a finite score gap into a lower bound on coherent dimension after explicit allowances for source infidelity, sequestration leakage, decoder error, and classical-transcript side channels. It should state only incompatibility with the declared $d$-bounded model, not a unique microscopic memory size.

## 9. Hard kill criteria

The pivot should be abandoned as the flagship contribution if the exact optimization reduces directly to Hsieh *et al.*'s instrument robustness, Lim *et al.*'s nondestructive discrimination tradeoff, or Ohst *et al.*'s published constrained-separability example after relabelling. It should also be killed if the only analytic statement is the perfect-return ratio $d/2^k$, if intermediate dimensions give only the convex hull of known endpoints, or if the proof covers only diagonal/iid instruments while the claim says arbitrary adaptive instruments.

A purely numerical hierarchy is appropriate as a tool, not as the central theorem. Likewise, the project should not claim a coding-theory contribution if $H$ appears only through a standard weight enumerator, or a tensor-network contribution if $d$ appears only as an ordinary MPS bond dimension. Finally, any experimental claim fails if an uncharged coherent subsystem can cross a cut, if the device can access $B_i$ during AUDIT, if separate devices implement the two branches, or if high conditional return fidelity is bought with unreported failure probability.

## 10. Recommended research sequence

The current parity theorem should remain frozen as the $k=1,d=1$ baseline. Work on the manuscript should pause while the new theorem is tested. First, formulate the $n=3,k=2,d=2$ comb exactly and derive primal and dual finite-dimensional programs. Second, solve it numerically at high precision for a grid of $\lambda$, including unrestricted transcript-conditioned recovery, and use the optimizer to conjecture an analytic strategy. Third, search for a reduction to the four closest frameworks: constrained-separability memory bounds, interactive-instrument robustness, nondestructive MES discrimination, and bounded-storage functional decoding. Only if no reduction closes the problem should effort move to an analytic proof and a robust statistical witness.

If that instance has a strict frontier, the natural software contribution is not another generic simulator. It is a reproducible bounded-bond optimizer that exports primal strategies, certified upper bounds, symmetry reductions, and dimension-witness thresholds. A later circuit layer can compile the honest accumulator and the verifier's decoder. Hardware forecasts should wait until the trusted coherent storage required for the sequestered $B_i$ carriers is counted explicitly.

## 11. Safe novelty language

The following sentence is supportable at this cutoff:

> Within the primary-source corpus searched through 19 August 2026, no work was located that derives the exact complete-syndrome audit versus all-carrier entanglement-return support function for a single adaptive streaming process with unlimited classical transcript and a bounded coherent temporal bond. The candidate contribution is this exact intermediate-dimension frontier and its robust dimension-witness consequence, not bounded memory, late choice, syndrome accumulation, weight enumerators, entanglement recovery, or semidefinite optimization separately.

Unsafe formulations include “the first coherent-memory hierarchy,” “the first reversible syndrome observer,” “a new relation between coding and quantum memory,” and “certification of $d$ quantum memory levels” without the access assumptions. The search result should be reported as a bounded absence statement, never as proof that no equivalent theorem exists.

## 12. Final go/no-go decision

**GO, but only to the $n=3,k=2,d=2$ theorem gate.** The problem is sufficiently differentiated to justify a focused derivation and adversarial numerical search. It is not yet sufficiently differentiated to justify rewriting the public paper around it.

The decision after that gate is binary. An exact strict intermediate-dimension frontier, with a proof over the full causal class, is a plausible solid publication and a useful functional quantum-memory benchmark. Failure to beat endpoint interpolation, or discovery that the frontier is a direct corollary of the closest papers, means the pivot should be documented as a negative result and not marketed as a new research program.

---

### Search update policy

Because four of the closest items appeared in 2026, repeat the title/abstract and citation-neighbourhood search immediately before any preprint submission. In particular monitor bounded coherent-memory process discrimination, interactive-instrument resource theories, limited-memory stabilizer learning, nondestructive discrimination with restricted entanglement, and dimension-reduced adaptive agents.

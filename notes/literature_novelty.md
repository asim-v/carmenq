# Literature and Novelty Audit

## Causally certified, reversibly erased quantum histories

**Cutoff:** 12 August 2026

**Document type:** adversarial scoping review and novelty audit; not a PRISMA systematic review

**Language:** English

**Evidence base:** primary papers, publisher records, DOI metadata, and first-party preprint archives

**Bibliography:** `references/library.bib`

---

## 1. Executive finding

The elementary protocol

\[
U_H \longrightarrow O_p \longrightarrow U_H^\dagger
\]

is not new. Its ingredients are already covered by reversible computation, phase kickback, quantum erasure, measurement reversal, and ordinary coherent circuit synthesis [@bennett1973logical; @bennett1989timespace; @cleve1998quantum; @herzog1995complementarity; @katz2008reversal]. Calling the intermediate components “histories,” “agents,” or “branches” does not change the implemented channel.

Several closer territories are also occupied:

- consistent/decoherent histories and redundant records already connect temporal alternatives to physical records, while a March 2026 peer-reviewed model separates approximate histories from objectivity [@griffiths1984consistent; @paz1993environment; @riedel2016objective; @ferte2026histories];
- process tensors and quantum combs already represent operational multi-time processes with interventions and memory [@chiribella2009networks; @pollock2018processes; @white2022tomography];
- temporal correlations already witness a minimum effective memory dimension, and device-independent quantum-memory certification is an active line [@vieira2024dimension; @sekatski2023memory; @santos2026deviceindependent];
- Wigner-friend analogues already include explicit quantum erasure of the friend’s memory [@elouard2021erasing];
- 2026 preprints already use “inter-branch communication” language and run related uncomputation circuits on superconducting hardware [@violaris2026branches; @altman2026wigner; @cogburn2026interbranch];
- quantum zero knowledge already covers QMA, classical-verifier quantum-prover protocols, non-destructive proofs on quantum states, and security against a verifier receiving a superposition of transcripts [@broadbent2016qma; @vidick2020classical; @colisson2025nondestructive; @coladangelo2026mpc];
- computationally constrained erasure and thermodynamic models of measurement already give nontrivial complexity/work and record-formation tradeoffs [@munson2025complexity; @latune2025measurement]; and
- a preprint posted only eight days before this cutoff claims an exact QEC–Quantum Darwinism tradeoff and a model-independent no-go theorem [@maity2026tradeoff].

The strongest remaining research candidate is therefore **not** “interference between reversible branches.” It is a joint operational task:

> Certify, under an explicit intervention and access model, that a multi-time quantum process required nontrivial causal memory; reveal only an allowed predicate of that process; decouple every other accessible final record from the hidden history; and recover interference in the same logical degree of freedom.

This audit found adjacent results for every clause separately, but did **not** locate a primary source that formulates and proves all clauses as one cryptographic/causal/interferometric primitive. That is a bounded search result, not a proof of absence or priority. A publishable manuscript must use wording such as “we did not identify” or “to our knowledge under the search described here,” never “this has never been considered.”

### Novelty verdict

| Proposed contribution | Verdict at cutoff |
|---|---|
| Compute a predicate, phase-kick it, uncompute, and interfere | Occupied; textbook-level quantum computation |
| Erase a path/friend memory and restore interference | Occupied experimentally and theoretically |
| Model a multi-time process with memory | Occupied by quantum combs/process tensors |
| Witness or certify quantum memory from temporal statistics | Occupied and rapidly advancing |
| Build a reversible finite-state “agent” | Occupied in adjacent quantum-agent/computational-mechanics work unless a new task or bound is supplied |
| Call coherent components “branches” and benchmark hardware | Directly occupied by 2026 preprints; interpretive language is not novelty |
| Relate record redundancy to recoverability/QEC | Occupied in models; a very recent 2026 preprint claims a general tradeoff |
| Join causal-memory certification, predicate-only disclosure, physical final-record decoupling, and recovered interference under anti-shortcut constraints | **Apparently underexplored; candidate only, requiring a formal priority search and a nontrivial theorem** |

---

## 2. Epistemic labels and search method

### 2.1 Evidence labels

The audit uses four labels:

- **Peer reviewed (PR):** version of record verified on a publisher or DOI page.
- **Canonical authority (CA):** a peer-reviewed foundational paper or field-defining review used to establish terminology and known results.
- **Preprint (PP):** first-party arXiv record; not treated as settled or as equivalent to peer review.
- **Audit inference (AI):** a conclusion drawn by comparing sources. It must not be cited as if it were a theorem in any one source.

### 2.2 Sources searched

Searches were run through current web indexes, then resolved to primary or authoritative destinations:

- APS journals and DOI pages;
- Nature, Nature Communications, Nature Physics, and npj Quantum Information;
- Quantum journal article pages;
- SIAM, ACM, IEEE, IBM Research, Springer, and institutional publication records;
- arXiv abstracts for 2026 preprints;
- PhilArchive only for a documented terminology collision, not as substantive physics evidence.

Representative query families included combinations of:

- `reversible computation`, `uncomputation`, `phase kickback`, `pebble game`;
- `decoherent histories`, `records`, `recoherence`, `quantum eraser`, `measurement reversal`;
- `process tensor`, `multi-round process matrix`, `causal discovery`, `memory dimension witness`;
- `device-independent quantum memory`, `temporal correlations`, `interventions`;
- `Wigner friend memory erasure`, `inter-branch communication circuit`;
- `zero knowledge QMA`, `non-destructive quantum zero knowledge`, `superposition transcript`;
- `QUALM`, `quantum algorithmic measurement`, `quantum adaptive agent`, `causal states`;
- `complexity-constrained erasure`, `thermodynamic cost measurement`;
- `Quantum Darwinism error correction tradeoff`, `redundant records`;
- `collapse model test`, `CSL`, `nanoparticle matter-wave interferometry`; and
- exact-phrase probes for `zero-record`, `erased history`, and combinations of causal certification, recoherence, and zero knowledge.

### 2.3 Limitations

This is a high-resolution scoping audit, not a complete priority search. It does not include exhaustive citation-network traversal, subscription database exports, conference proceedings in every subfield, non-English literature, patents, theses, or all versions posted after the cutoff. Negative search results are therefore weak evidence.

Before a first-submission novelty claim, perform all of the following:

1. reproducible searches in Web of Science, Scopus, INSPIRE, dblp, ACM DL, IEEE Xplore, and Google Scholar;
2. backward and forward citation chasing from the closest papers in Sections 4–10;
3. an author-name/topic search around process tensors, quantum verification, and coherent cryptography;
4. explicit searches for unpublished workshop talks and theses;
5. a dated search log and exported result set; and
6. external review by at least one expert in process tensors and one in quantum cryptography.

---

## 3. Adjacency map

The map below separates what each field already supplies from what would have to be added.

| Field | Established capability | Territory already occupied | Missing element for the joint task |
|---|---|---|---|
| Reversible computation | Simulate irreversible computations reversibly; trade time for workspace | Compute–copy output–uncompute; reversible pebbling; measurement-assisted qubit recycling | No certification that a particular internal causal history occurred |
| Phase algorithms | Convert a Boolean predicate into a relative phase | The core (U_H O_p U_H^\dagger) pattern | Anti-shortcut condition excluding a directly compiled phase oracle |
| Decoherence and records | Quantify loss of reduced-state coherence through environmental correlations | Record formation, einselection, information–disturbance, approximate recovery | Joint certification of record formation and later final-record decoupling |
| Consistent histories | Assign probabilities to compatible sets of temporal alternatives | Decoherence functional, class operators, redundant records of histories | Interactive operational proof about an erased history, rather than a descriptive formalism |
| Quantum eraser / reversal | Restore conditional interference when path information is made unavailable; reverse some partial measurements | Photonic erasers, superconducting weak-measurement reversal, trapped-ion measurement undoing | Deterministic unitary reset plus causal-memory certification and predicate-only disclosure |
| Process tensors / combs | Represent all multi-time statistics under interventions | Memory, non-Markovianity, process tomography, multi-round local memory | A final interference/privacy requirement and adversarial shortcut model |
| Causal discovery / memory witnesses | Infer causal order or lower-bound effective memory from temporal data | Quantum causal discovery, dimension witnesses, DI memory proposals and experiments | Couple the witness to subsequent coherent erasure in a single defined task |
| Wigner-friend analogues | Treat a small “friend” or lab as a controlled quantum system | Local-friendliness tests; explicit memory erasure; branch-message circuits | Functional-agent complexity certified independently of interpretive labels |
| Quantum zero knowledge | Prove statements while limiting verifier knowledge | QMA ZK, classical-verifier ZK, non-destructive state proofs, superposition-secure transcripts | Physical decoupling of final history records plus interferometric readout; these are not implied by simulator-based ZK |
| QUALMs / complexity of experiments | Compare coherent experimental access models with classical or incoherent ones | Provable resource separations for experimental tasks | A separation tailored to history predicates and causal-memory access |
| Computational mechanics / agents | Quantify minimal predictive memory; construct quantum memory-efficient agents | Quantum causal states and adaptive agents | Recoherence cost as a function of irreducible causal memory and policy depth |
| Quantum thermodynamics | Bound work and entropy costs of measurement/reset under resource constraints | Complexity entropy; finite-time measurement costs; reversible limits | A new work–visibility–causal-complexity theorem, not a naive Landauer restatement |
| Quantum Darwinism / QEC | Quantify redundant classical records and encoded quantum information | Darwinism–encoding transitions; decoherent-histories/objectivity separation; recent claimed QEC–Darwinism no-go tradeoff | Only a more restricted operational theorem may remain; general rhetoric is occupied |
| Collapse-model tests | Predict model-specific loss/heating/radiation and constrain parameters | CSL/DP bounds, optomechanical proposals, x-ray searches, matter-wave interferometry | A pre-registered parameter-estimation experiment at appropriate mass/geometry; small logical circuits are not competitive tests |

The key structure is not a simple bridge between two fields. It is a five-way junction:

```text
multi-time causal certification
              |
              v
reversible agent/history --> predicate-only disclosure --> physical record decoupling
              |                                      |
              +--------------> recovered interference+
                                     |
                         complexity/access separation
```

Each arrow has prior art. The candidate contribution is a rigorous definition and theorem for the whole junction.

---

## 4. Reversible computation, phase kickback, and the endpoint-equivalence problem

Landauer tied logically irreversible operations to heat generation, while Bennett established that general computation can be embedded in logically reversible dynamics and later quantified time/space tradeoffs [@landauer1961irreversibility; @bennett1973logical; @bennett1989timespace]. Quantum algorithms routinely use reversible evaluation followed by a phase kickback and uncomputation [@cleve1998quantum]. Modern pebbling work further shows that the workspace cost of coherent classical subroutines is itself a developed resource question, including measurement-assisted “spooky” pebbling [@kornerup2025spooky].

Consequently, the following statement is established and not novel:

\[
U_f\,|x\rangle|0\rangle|{-}\rangle
\mapsto
|x\rangle|f(x)\rangle(-1)^{f(x)}|{-}\rangle
\mapsto
(-1)^{f(x)}|x\rangle|0\rangle|{-}\rangle.
\]

Calling (f(x)) a property of a “history,” and calling the work register an “observer,” does not alter this equivalence.

### 4.1 Endpoint-equivalence obstruction

If an experimenter has unrestricted knowledge and control of a circuit (U), then any endpoint transformation can in principle be represented by a directly synthesized circuit implementing the same channel. Output statistics alone cannot identify the implementation history.

This is an identifiability statement, not a claim that circuit synthesis is always efficient. A scientifically meaningful “genuine history” criterion therefore requires a declared access model, such as:

- unknown but repeatable dynamics used as an oracle;
- locality restrictions on allowed gates;
- challenges selected only after an intermediate memory has formed;
- a black-box prover not trusted to implement a claimed internal circuit;
- a lower bound on memory dimension from interventional data; or
- a QUALM-style separation between coherent access to one evolving process and independent samples.

Without at least one such restriction, a “direct phase” control is not merely an experimental baseline: it proves that the proposed history is not operationally identifiable.

### 4.2 What could still be new

A lower bound of the following form could be new if it is not reducible to existing pebbling, query-complexity, or process-discrimination results:

\[
\text{cost}_{\rm shortcut}(P,\mathcal A)
\;\geq\;
g(d_{\rm mem},T,\ell),
\]

where (mathcal A) is the allowed access model, (d_{\rm mem}) is a causally certified memory dimension, (T) is temporal depth, and (ell) is a locality parameter. The theorem must compare against an optimized direct implementation, not an intentionally inefficient straw man.

---

## 5. Decoherence, records, histories, erasure, and recovery

Decoherence theory already explains how correlations with inaccessible degrees of freedom suppress reduced-state interference [@zurek2003decoherence]. For two controlled alternatives correlated with environmental states,

\[
|\Psi\rangle
=\alpha|0\rangle|\Phi_0\rangle|e_0\rangle
+\beta|1\rangle|\Phi_1\rangle|e_1\rangle,
\]

the off-diagonal term is proportional to (\langle e_1|e_0\rangle). The visibility–which-way tradeoff is not new [@englert1996fringe], nor are the standard fidelity/trace-distance and information–disturbance relations [@uhlmann1976transition; @fuchs1999cryptographic; @kretschmann2008information].

Quantum eraser experiments show that appropriate conditioning or basis choice can recover interference when distinguishing information is made unavailable [@herzog1995complementarity; @kim2000delayed]. These experiments must not be described as retrocausal rewriting. Moreover, postselected quantum erasure and deterministic unitary uncomputation are operationally different and must be reported separately.

Measurement reversibility has its own literature. Logical and probabilistic reversibility conditions were formalized before modern quantum processors [@ueda1996logical], partial-measurement reversal was demonstrated in a superconducting qubit [@katz2008reversal], and measurement undoing was demonstrated with trapped ions [@schindler2013undoing]. The interpretation of reversibility relative to accessible records has also been analysed explicitly [@zurek2018reversibility]. A 2026 peer-reviewed paper studies foundational implications of information erasure on measurement statistics [@montina2026erasure].

Consistent histories supply class operators, decoherence conditions, and probability rules for compatible temporal alternatives [@griffiths1984consistent]. Decoherence and histories were explicitly joined by Paz and Zurek [@paz1993environment]. Riedel, Zurek, and Zwolak then connected objective pasts to redundant records of consistent histories [@riedel2016objective]. A March 2026 peer-reviewed model finds approximate decoherent histories in both an apparatus phase and a scrambler phase, while redundant pointer-state records distinguish the apparatus/objectivity regime [@ferte2026histories]. Thus “histories plus records,” and even the distinction between histories and objectivity, is occupied territory.

The January 2026 preprint by Strasberg and collaborators is especially close. It studies approximate decoherence, records, and recoherence in isolated systems, including Petz-recovery structure and long histories [@strasberg2026records]. It is provisional, but any general record–recoherence theorem must be compared against it explicitly.

### 5.1 Correct interpretation of “zero record”

Unitary evolution does not literally destroy arbitrary quantum information. A defensible requirement is **decoupling from the forbidden history variable**, conditional on the allowed predicate. Introduce a reference (H) that labels or purifies the history alternatives and let (R_f) contain every final accessible register other than the authorized readout (P). A candidate privacy condition is

\[
I(H:R_f\mid P)\leq \varepsilon,
\]

or a composable trace-distance/diamond-norm analogue. This says that the remaining accessible records reveal no appreciable information beyond (P). It does not say that information vanished from the global wavefunction, and it does not permit ignoring uncontrolled radiation, controllers, calibration logs, measurement records, or error-correction syndromes correlated with (H).

### 5.2 Occupied versus open

- **Occupied:** path information suppresses interference; erasure can conditionally recover it; partial measurements can be reversed; redundant records define an objective past; approximate recoherence can be studied with recovery maps.
- **Candidate:** prove a composable bound that simultaneously includes causally certified internal memory, authorized predicate leakage, final forbidden-record decoupling, and recovered visibility under a declared access model.

---

## 6. Multi-time processes, causal discovery, and memory witnesses

Quantum combs and process tensors are the natural operational language for the project. Quantum networks represent higher-order transformations and memory channels [@chiribella2009networks]. Process tensors capture multi-time statistics under arbitrary allowed interventions and provide an operational account of non-Markovian quantum processes [@pollock2018processes]. Process-tensor tomography has been developed and demonstrated on quantum processors [@white2022tomography].

The process framework prevents a major conceptual mistake: a final state or channel does not determine the sequence of interventions and memory mechanisms that produced it. A history claim must be relative to an intervention set (\mathfrak I).

Giarmatzi and Costa give a quantum causal-discovery algorithm that distinguishes causal order, Markovian structure, and latent memory at the level of a process matrix [@giarmatzi2018causal]. Multi-round process matrices explicitly allow local parties to retain memory across several rounds [@hoffreumon2021multiround]. Instrument-specific memory effects and process recovery have already been experimentally demonstrated [@guo2021memory].

Memory certification is now particularly close to the proposed core:

- temporal correlations can lower-bound the minimum dimension of an effective environment acting as memory [@vieira2024dimension];
- a device-independent certification framework for quantum memories was published in 2023 [@sekatski2023memory];
- randomized measurements can efficiently certify temporal quantum correlations [@liu2025temporal]; and
- a January 2026 preprint reports a two-time, causal-inequality-based, device-independent approach and a trapped-ion proof of principle certifying 35 ms of qubit memory [@santos2026deviceindependent].

The comparison class must also state whether its memory is classical or quantum. Quantum random-access codes already show that one qubit can outperform one classical bit on selected retrieval tasks [@ambainis2002dense].

### 6.1 Consequence for novelty

“Certify that memory existed” is no longer an open headline. The possible opening is the conjunction:

1. lower-bound the causal memory or exclude a memoryless process class;
2. use the same hidden multi-time process to compute a predicate;
3. decouple the final apparatus from all other history information; and
4. recover an interferometric signal whose visibility is benchmarked against a matched identity echo.

Whether all four can be certified in one shot is itself nontrivial. Ordinary intermediate readout changes the process and leaves a transcript. A practical protocol may require coherent challenges, randomized experimental contexts across repeated trials, or a cryptographic verifier. The manuscript must state which certification is single-shot and which is inferred from an ensemble of matched runs.

---

## 7. Wigner-friend analogues: useful model, weak novelty anchor

Extended Wigner-friend work formalizes nested agents and observer-dependent records [@frauchiger2018quantum]. Proof-of-principle experiments use small photonic degrees of freedom as “observers,” not conscious humans [@proietti2019observer; @bong2020strong]. The word “observer” in these experiments denotes a controlled information-processing subsystem.

Elouard and collaborators are a direct antecedent: their paper explicitly analyzes “quantum erasing the memory of Wigner’s friend” [@elouard2021erasing]. This blocks any novelty claim based merely on a friend acquiring a record and an external controller erasing it.

The 2026 literature makes an interpretation-led framing even riskier:

- Jones and Mueller argue, in a peer-reviewed June 2026 paper, that key Wigner-friend structures can have classical analogues involving duplicated agents and should be understood beyond quantum foundations [@jones2026significance].
- Violaris proposes “communication across multiverse branches” but requires the sender not to retain memory of the message [@violaris2026branches].
- Altman compiles a related five-qubit Wigner-friend-style circuit on IBM hardware and explicitly frames it as an operational benchmark, not an interpretation test [@altman2026wigner].
- Cogburn benchmarks the related message-transfer primitive across superconducting architectures and scales message registers up to 32 qubits [@cogburn2026interbranch].

The last three are preprints, not settled results. They nevertheless establish a serious priority risk for any claim framed as reversible inter-branch communication or a Wigner-friend circuit benchmark.

### Recommendation

Use Wigner-friend systems only as one implementation family. Define “agent” functionally and quantify it through causal-state or process-memory complexity. Do not use an Everettian interpretation in the title, abstract, theorem names, or primary claims.

---

## 8. Quantum zero knowledge and the danger of a false analogy

Zero knowledge is an attractive connection because the desired output is “learn a proposition without learning the witness.” The literature is mature:

- Watrous established zero knowledge against quantum attacks [@watrous2009zeroknowledge].
- Broadbent, Ji, Song, and Watrous showed that every QMA problem has a quantum zero-knowledge proof system under stated cryptographic assumptions [@broadbent2016qma].
- Vidick and Zhang constructed classical-verifier, quantum-prover zero-knowledge arguments for quantum computations [@vidick2020classical].
- Colisson, Grosshans, and Kashefi gave non-interactive, non-destructive zero-knowledge proofs about quantum states [@colisson2025nondestructive].
- Coladangelo and collaborators generalized MPC-in-the-head to quantum computation and, in a peer-reviewed July 2026 paper, obtained protocols secure even when a verifier can obtain a superposition of transcripts [@coladangelo2026mpc].

Therefore neither “zero knowledge about a quantum computation,” “preserve the quantum state,” nor “security against coherent transcript access” can be claimed as new.

### 8.1 Zero knowledge is not zero physical record

Cryptographic zero knowledge is normally simulator-based: the verifier’s view can be simulated without the witness, subject to computational or statistical definitions. A transcript may physically exist and still be zero knowledge. Conversely, a final register may be decoupled from a history but the protocol may fail cryptographic soundness or simulation-based privacy.

The proposed physical condition

\[
I(H:R_f\mid P)\le \varepsilon
\]

is therefore not automatically a zero-knowledge property. A rigorous paper must either:

1. avoid the term “zero knowledge” and call the condition final-record decoupling; or
2. define an interactive proof with completeness, soundness, and a simulator, then separately prove physical decoupling and recoherence.

### 8.2 Candidate cryptographic bridge

A nontrivial bridge could use challenges chosen after relevant memory is formed:

1. a prover/process evolves coherently and commits to an intermediate causal state;
2. a verifier supplies a later challenge (x_t), possibly coherently;
3. a response depends on the earlier memory and later challenge;
4. an accept bit is phase-kicked into an authorized control;
5. response, proof, memory, and work registers are uncomputed; and
6. the final interference estimates acceptance while forbidden transcript registers satisfy a decoupling bound.

This resembles an interactive proof but adds a physical recoherence condition. It is promising only if a cheating direct-phase strategy is excluded within a formal access/locality model.

---

## 9. QUALMs, computational mechanics, and operational observer complexity

The QUALM framework studies physical experiments from a computational-complexity perspective and proves that coherent access can offer large resource advantages for some experimental tasks [@aharonov2022qualm]. This is a more defensible comparison class than “quantum superposition versus classical random sampling,” which is usually too weak.

Computational mechanics provides a separate way to define memory complexity. Quantum models can reproduce stochastic processes with less stored information than optimal classical models [@gu2012complexity], and quantum adaptive agents with efficient long-term memories have been characterized [@elliott2022agents]. Thus “quantum agent with memory” is already an established connection.

### 9.1 Recommended observer metrics

Do not define observer complexity only by raw qubit count. At minimum, report:

- minimum compatible memory dimension under the declared instruments;
- classical and quantum statistical memory costs;
- number of causal states or effective process-tensor bond dimension;
- policy depth and number of adaptive decisions;
- temporal horizon over which prior inputs affect later actions;
- redundancy of internally stored records;
- reversible workspace and non-Clifford cost; and
- recovery-channel complexity at a target visibility.

A candidate project-specific quantity is

\[
\mathcal C_{\rm CR}(d,v,\varepsilon;\mathfrak I)
=
\min_{\Pi}
\left\{
\operatorname{cost}(\Pi):
d_{\rm causal}(\Pi)\ge d,
V(\Pi)\ge v,
I(H:R_f\mid P)\le\varepsilon
\right\},
\]

where (\mathfrak I) fixes the allowed interventions. This expression is only a research definition. It becomes a contribution only if the authors prove operational properties, monotonicity or bounds, and separation from process-tensor rank, ordinary circuit depth, and existing statistical complexity.

### 9.2 Strong target theorem

The most valuable result would be a resource separation:

> Under a specified oracle or multi-time access model, coherent interaction with one persistent process estimates a history predicate using polynomial resources, whereas every classical-sampling, memoryless, or endpoint-only strategy needs asymptotically more resources.

The comparison must be against the best allowed classical/incoherent strategy, not merely against dephasing the control qubit.

---

## 10. Thermodynamics and Quantum Darwinism: valuable secondary axes, poor novelty headlines

### 10.1 Thermodynamic cost

Landauer's original result concerns logically irreversible operations [@landauer1961irreversibility]. Complexity-constrained quantum thermodynamics relates the minimum work of erasure to computational limitations through complexity entropy [@munson2025complexity]. A microscopic thermodynamic model of quantum measurement already treats correlation formation, redundant environmental records, apparatus reset, reversible limits, and finite-time work costs [@latune2025measurement].

It would therefore be incorrect to claim that uncomputing (m) coherent bits necessarily costs (m k_B T\ln 2). Ideal reversible uncomputation is not blind many-to-one erasure. A thermodynamic work package must distinguish:

- logically reversible uncomputation;
- irreversible reset of unknown or uncontrolled registers;
- finite-time dissipation;
- control-field and calibration cost;
- refrigeration overhead;
- error correction and syndrome handling; and
- entropy exported into inaccessible degrees of freedom.

A possible new result is a constrained Pareto surface among causal memory, reset error, visibility, work, time, and control complexity. It is not a primary novelty claim until shown to go beyond [@munson2025complexity] and [@latune2025measurement].

### 10.2 Quantum Darwinism and QEC

Quantum Darwinism formalizes the redundant proliferation of classically accessible information [@blumekohout2006darwinism]. Redundant records have already been linked to consistent histories [@riedel2016objective]. Solvable random-circuit models already exhibit transitions between a Darwinistic phase and a quantum-encoding phase [@ferte2024transition]. A March 2026 peer-reviewed model distinguishes the emergence of approximate decoherent histories from the emergence of redundant objective records [@ferte2026histories], while a May 2026 peer-reviewed paper extends pointer-state and objectivity analysis to noncommuting evolutions [@chisholm2026noncommuting].

Most importantly, Maity, Onggadinata, and Koh posted arXiv:2608.03944 on 4 August 2026, claiming an exact QEC–Quantum Darwinism tradeoff in a solvable model and a model-independent no-go theorem [@maity2026tradeoff]. At this cutoff it is an extremely recent, unreviewed preprint. Its claims must be checked line by line, but its existence makes “objectivity versus recoverability” an unsafe headline for novelty.

Recommended use:

- treat record redundancy as a scaling parameter in simulations;
- reproduce known Darwinism–encoding behavior as validation;
- compare any theorem explicitly to the August 2026 preprint; and
- keep this as a secondary work package unless a sharper causal/interactive restriction yields a distinct bound.

---

## 11. Collapse models and tests of non-unitarity

Objective-collapse models make quantitative, model-dependent predictions. The field has mature reviews and experimental constraints [@bassi2013collapse]. Optomechanical and near-field interferometric proposals quantify sensitivity to CSL-like modifications under realistic noise [@nimmrichter2014optomechanical; @gasbarri2021nearfield]. The Majorana Demonstrator found no spontaneous-radiation signal and set strong CSL constraints in a specified parameter region [@arnquist2022majorana]. In March 2026, XENONnT reported stronger x-ray constraints on Markovian CSL and Diósi–Penrose parameters [@aprile2026collapse]. In January 2026, nanoparticle matter-wave interferometry extended coherent matter-wave experiments to massive sodium nanoparticles [@pedalino2026nanoparticle].

These results do not imply that a small logical reversible-history circuit is a competitive collapse-model test. Any non-unitarity claim must:

1. choose a named model such as mass-proportional CSL;
2. predefine the parameter region and predicted scaling with mass, geometry, and superposition separation;
3. calibrate ordinary dephasing, relaxation, leakage, control error, and non-Markovian noise;
4. normalize against a matched identity echo of the same depth;
5. avoid postselection or report success probability separately;
6. fit competing environmental and collapse models; and
7. report exclusion regions, not “unexplained decoherence.”

Collapse tests are therefore a long-range extension, not evidence for the core protocol and not an interpretation test.

---

## 12. Candidate joint primitive

This section states the narrow research object that survived the audit. It is deliberately presented as a **candidate definition**, not an established novelty claim.

### 12.1 Task instance

An instance contains:

- a persistent multi-time quantum process (\Upsilon_{T:0});
- an allowed instrument family (\mathfrak I=\{\mathcal I^{(x_t)}_t\}_{t=0}^{T-1});
- a hidden history reference (H);
- an internal memory (M), world/workspace (W,G), history control (B), and authorized output (P);
- a challenge schedule (X=(x_1,\ldots,x_T)), at least part of which is unavailable before relevant internal states form;
- a predicate (p(H,X)); and
- a model class (\mathcal M_{<d}) representing shortcut processes with memory dimension below (d), missing causal links, or forbidden direct access.

### 12.2 Desired properties

1. **Completeness.** A valid process and honest reversible strategy cause acceptance/interference with probability at least (c).
2. **Soundness.** Every strategy in the excluded shortcut class succeeds with probability at most (s<c).
3. **Causal authenticity.** Interventional statistics reject (\mathcal M_{<d}) at a declared confidence level or violate a valid memory/causal witness.
4. **Predicate correctness.** The authorized readout estimates (p(H,X)) with error at most (\delta_p).
5. **Final-record decoupling.** Every accessible forbidden record (R_f) obeys a composable bound such as (I(H:R_f\mid P)\le\varepsilon).
6. **Reset fidelity.** The channel on (W,M,G) is close to a branch-independent reset channel, preferably in entanglement fidelity or diamond distance, not only for one basis state.
7. **Recoherence.** The history control retains visibility (V\ge v), measured in at least two complementary equatorial bases and normalized to a matched echo.
8. **No postselection ambiguity.** If success is heralded, report both conditional performance and the unconditional success probability.
9. **Interpretive neutrality.** All formal claims are stated in channel/process language.

### 12.3 Why the conjunction is nontrivial

- A direct phase oracle can imitate endpoint statistics unless the access model blocks it.
- A classical interventional measurement can certify memory but leave a transcript and destroy the coherence needed later.
- A coherent challenge can preserve reversibility but makes the verifier part of the quantum process, complicating soundness and zero-knowledge definitions.
- Simulator-based zero knowledge does not imply physical record decoupling.
- Final-record decoupling does not imply a cryptographic proof of knowledge.
- Reset on one initial state does not certify channel-level reversal.
- Perfect visibility can coexist with a trivial internal implementation.

These tensions are the project’s real theoretical content.

### 12.4 Conservative novelty language

Permitted draft wording:

> We define a task that jointly imposes (i) an interventional lower bound on causal memory, (ii) predicate-only authorized disclosure, (iii) final decoupling of all other accessible history records, and (iv) recovered interference. In the primary literature located by our search through 12 August 2026, these conditions appear separately but we did not identify a protocol analyzing their conjunction under an explicit anti-shortcut access model.

Forbidden wording:

- “No one has ever connected these fields.”
- “The first machine to communicate between branches.”
- “A proof of many worlds.”
- “A zero-knowledge proof” before a simulator-based definition and theorem exist.
- “Information is destroyed” when only selected registers are reset or decoupled.

---

## 13. Simulation and experimental implications of the audit

A statevector demonstration of (U_H O_p U_H^\dagger) is necessary for validation but insufficient for novelty. The reference simulation suite should include the following layers.

### 13.1 Required ideal baselines

1. **Direct compiled phase:** implements the same endpoint predicate without the claimed history.
2. **Coherent multi-time history:** full world, memory, policy, predicate, and inverse.
3. **Memoryless process:** optimized process in (\mathcal M_{<d}).
4. **Classical mixture:** same populations, duration, and marginal noise as the coherent run.
5. **Retained memory:** skip memory reset.
6. **Retained garbage:** reset visible registers but leave one correlated work ancilla.
7. **Matched identity echo:** same gate depth/connectivity with no predicate.
8. **Challenge-ablation model:** reveal all challenges in advance to quantify the shortcut.

### 13.2 Required noise models

- local dephasing and amplitude damping;
- coherent over-rotation and inverse-gate correlation;
- branch-correlated leakage into an explicit environment;
- non-Markovian noise represented as a small process tensor;
- crosstalk and spectator memory;
- readout error and state-preparation error;
- randomized compiler or routing variation; and
- optional model-specific CSL channel only at scales where a meaningful bound can be stated.

### 13.3 Required reported quantities

- raw and echo-normalized (X/Y) visibility;
- predicate accuracy;
- state reset fidelity and channel/entanglement fidelity;
- residual mutual information or trace-distance leakage;
- minimum compatible memory dimension from a temporal witness;
- process-tensor distance from a memoryless model;
- logical qubits, physical qubits, depth, two-qubit gates, non-Clifford gates, and ancillas;
- postselection/heralding probability;
- confidence intervals and shot complexity; and
- optimized shortcut cost under the declared access model.

### 13.4 The first meaningful numerical result

A credible first result would not be “the predicate is recovered perfectly in an ideal simulator.” It would be one of:

- a separation between the valid process and every bounded-memory shortcut in a specified model;
- a tight relation among certified memory dimension, residual leakage, and echo-normalized visibility;
- an impossibility theorem showing that a chosen form of causal certification cannot coexist with same-shot transcript erasure;
- a sample-complexity advantage from coherent multi-time access; or
- a verified protocol where the challenges prevent precompilation and all correlated work registers are audited.

---

## 14. Claim-to-source ledger

The ledger records what may be safely asserted and its evidential basis.

| Claim | Status | Primary/authority source | Manuscript use |
|---|---|---|---|
| Logical irreversibility is associated with thermodynamic cost | Established | [@landauer1961irreversibility] | Distinguish blind erasure from unitary uncomputation |
| General computation can be embedded reversibly | Established | [@bennett1973logical] | Background, not novelty |
| Reversible simulation has time/workspace tradeoffs | Established | [@bennett1989timespace; @kornerup2025spooky] | Resource model |
| A predicate can be moved into phase and work registers uncomputed | Standard quantum algorithmic technique | [@cleve1998quantum] | Lemma only; explicitly non-novel |
| Environmental records suppress reduced-state coherence | Established | [@zurek2003decoherence] | Physical model |
| Which-way distinguishability constrains fringe visibility | Established | [@englert1996fringe] | Baseline bound |
| Root fidelity bounds purification overlap and is related to trace distance | Established | [@uhlmann1976transition; @fuchs1999cryptographic] | Mixed-record recovery bounds |
| Erasing path distinguishability can recover conditional interference | Experimentally established | [@herzog1995complementarity; @kim2000delayed] | Control, not novelty |
| Some partial quantum measurements can be reversed conditionally | Theoretical/experimental | [@ueda1996logical; @katz2008reversal; @schindler2013undoing] | Distinguish from deterministic uncomputation |
| Consistent histories supply probabilities for compatible temporal sequences | Established formalism | [@griffiths1984consistent] | Definitions |
| Decoherence and consistency of histories are linked | Established | [@paz1993environment] | Background |
| Redundant records can select an objective past among histories | Established model/result | [@riedel2016objective] | Darwinism adjacency |
| Approximate records and recoherence in long isolated histories are under active study | 2026 preprint, provisional | [@strasberg2026records] | Closest-prior-work discussion |
| Quantum combs/process tensors represent multi-time processes with memory | Established | [@chiribella2009networks; @pollock2018processes] | Core formalism |
| Non-Markovian process tensors can be reconstructed experimentally | Established | [@white2022tomography] | Simulation/experimental method |
| Interventions allow quantum causal discovery and latent-memory detection | Established | [@giarmatzi2018causal] | Causal-authenticity definition |
| Multi-round process matrices include local memory across rounds | Established | [@hoffreumon2021multiround] | Multi-agent protocol formalism |
| Temporal correlations can lower-bound effective memory dimension | Established | [@vieira2024dimension] | Observer-complexity metric |
| Device-independent quantum-memory certification is an existing program | Established theory; 2026 experiment still preprint | [@sekatski2023memory; @santos2026deviceindependent] | Direct adjacency; no novelty claim |
| A qubit can outperform one classical bit in a random-access retrieval task | Established | [@ambainis2002dense] | Null classes must distinguish classical and quantum memory |
| Wigner-friend analogues use small physical systems, not conscious observers | Established proof-of-principle practice | [@proietti2019observer; @bong2020strong] | Epistemic firewall |
| Quantum erasure of a Wigner-friend memory has been explicitly proposed | Established | [@elouard2021erasing] | Direct antecedent |
| Inter-branch message circuits have 2026 superconducting-hardware preprints | Provisional | [@altman2026wigner; @cogburn2026interbranch] | Priority-risk disclosure |
| These circuits do not discriminate interpretations | Explicitly acknowledged in direct preprint | [@altman2026wigner] | Interpretive neutrality |
| Quantum zero knowledge against quantum verifiers exists | Established | [@watrous2009zeroknowledge] | Cryptographic background |
| QMA has quantum zero-knowledge proof systems under stated assumptions | Established | [@broadbent2016qma] | Direct adjacency |
| Classical verifiers can have ZK arguments for quantum computations | Established | [@vidick2020classical] | Direct adjacency |
| Non-destructive ZK proofs about quantum states exist | Established 2025 | [@colisson2025nondestructive] | Prevent overclaim |
| ZK security against superpositions of transcripts exists under assumptions | Peer-reviewed 2026 | [@coladangelo2026mpc] | Prevent overclaim; coherent-verifier model |
| Coherent experimental access can yield complexity advantages | Established in QUALM tasks | [@aharonov2022qualm] | Model for a strong comparison |
| Quantum models can reduce predictive memory cost | Established | [@gu2012complexity] | Observer metric |
| Quantum adaptive agents with efficient memory already exist as a formal line | Established | [@elliott2022agents] | Prevent “first quantum agent” claim |
| Complexity constraints change minimum erasure work | Established 2025 | [@munson2025complexity] | Thermodynamic adjacency |
| Measurement, redundant records, reset, and finite-time cost have a microscopic thermodynamic model | Established 2025 | [@latune2025measurement] | Prevent naive thermodynamic novelty claim |
| Quantum Darwinism quantifies redundant classical records | Established | [@blumekohout2006darwinism] | Redundancy metric |
| Models exhibit Darwinism-to-encoding transitions | Established 2024 | [@ferte2024transition] | Scaling baseline |
| Approximate decoherent histories need not coincide with redundant objective records | Peer-reviewed 2026 model | [@ferte2026histories] | Prevent conflating history consistency with objectivity |
| A general QEC–Darwinism tradeoff is claimed | Very recent preprint, provisional | [@maity2026tradeoff] | High-priority comparison; do not treat as settled |
| Collapse models are quantitatively testable and constrained | Established field | [@bassi2013collapse; @arnquist2022majorana; @aprile2026collapse] | Long-range extension only |
| Nanoparticle matter-wave interference advanced in 2026 | Peer-reviewed experiment | [@pedalino2026nanoparticle] | Technology calibration; not an observer experiment |
| No endpoint statistic alone proves a specific internal implementation | Audit inference from operational equivalence and process formalism | [@cleve1998quantum; @pollock2018processes] | State as inference and prove within chosen model |
| No located paper combines all four joint requirements | Search result, not a literature fact | This dated audit | Use only with caveat and updated search |

---

## 15. What is occupied, what is plausibly open, and what would falsify the framing

### 15.1 Occupied

The following cannot serve as the principal contribution:

1. the algebraic survival of a predicate phase after uncomputation;
2. path-memory erasure and interference recovery;
3. a one-bit friend/observer implemented as an ancilla;
4. a generic visibility–distinguishability inequality;
5. process-tensor modeling of a multi-step memory;
6. witnessing that an environment has minimum dimension (d);
7. quantum zero knowledge for a quantum computation;
8. non-destructive proof about a quantum state;
9. “superposition of transcripts” as a new adversarial idea;
10. coherent experimental access as a general source of advantage;
11. quantum agents with compressed memories;
12. work cost of erasure under complexity constraints;
13. objectivity-versus-encoding rhetoric; or
14. inter-branch message transfer as a circuit/hardware benchmark.

### 15.2 Plausibly open joint task

The following conjunction remains a defensible target **only pending a complete priority search**:

> A protocol with interventional causal-memory soundness, an authorized predicate phase, composable final-record decoupling, and recovered interference, together with an explicit access model that makes a direct compiled phase an invalid or resource-expensive shortcut.

The strongest versions would add either:

- a provable coherent-versus-incoherent resource separation;
- a lower bound connecting irreducible causal memory to recovery complexity;
- a no-go theorem for simultaneous same-shot certification and record erasure; or
- a device-independent or semi-device-independent statement with realistic assumptions.

### 15.3 Kill criteria

Abandon the new-field framing if any of the following occurs:

1. a prior paper is found with the same joint task and security/causal definitions;
2. the “agent” can be removed with no change in allowed resources or certified statistics;
3. all challenges are known early enough that the result can be compiled into a direct phase at equal cost;
4. the memory witness certifies only a generic environment, not the claimed agent/policy structure;
5. final-record privacy is merely ordinary cryptographic ZK with no independent physical decoupling statement;
6. physical decoupling is only reset fidelity on one input state;
7. the observed visibility is fully explained by an ordinary matched echo;
8. postselection is essential but its success probability erases the claimed advantage;
9. the resource measure reduces to process-tensor bond dimension, circuit depth, or a known statistical-complexity quantity without a new theorem; or
10. the result’s only distinctive feature is interpretation language.

---

## 16. 2026 watchlist and provisionality notice

The cutoff falls in an unusually active year. These sources require special handling.

| Date | Source | Status at cutoff | Relevance |
|---|---|---|---|
| 13 Jan 2026 | Violaris, arXiv:2601.08102 | Preprint | Direct “inter-branch communication” and memory-erasure framing |
| 20 Jan 2026 | Santos et al., arXiv:2601.14191 | Preprint | Device-independent temporal memory certification; trapped-ion proof of principle |
| 21 Jan 2026 | Pedalino et al., Nature 649 | Peer reviewed | Nanoparticle matter-wave interferometry; scale calibration |
| 22 Jan 2026 | Altman, arXiv:2601.16004 | Preprint | Wigner-friend circuit benchmark on superconducting hardware |
| 27 Jan 2026 | Strasberg et al., arXiv:2601.19703 | Preprint | Approximate records, decoherence, recoherence, recovery |
| 27 Jan 2026 | Cogburn, arXiv:2601.19762 | Preprint | Multi-architecture inter-branch message benchmark |
| 5 Mar 2026 | Ferté, Farci, and Cao, Phys. Rev. Lett. 136 | Peer reviewed | Decoherent histories versus redundant objective records |
| 23 Mar 2026 | XENON Collaboration, Phys. Rev. Lett. 136 | Peer reviewed | Current spontaneous-radiation constraints on collapse models |
| 30 Apr 2026 | Montina and Wolf, Phys. Rev. A 113 | Peer reviewed | Information erasure upon measurement |
| 22 May 2026 | Chisholm et al., APS Open Science 1 | Peer reviewed | Quantum Darwinism for noncommuting evolutions |
| 4 Jun 2026 | Wu, PhilArchive | Independent manuscript; not peer reviewed | Exact `zero-record` terminology collision only |
| 30 Jun 2026 | Jones and Mueller, Quantum 10 | Peer reviewed | Wigner-friend structure beyond specifically quantum narratives |
| 15 Jul 2026 | Coladangelo et al., Quantum 10 | Peer reviewed | Superposition-secure quantum zero knowledge |
| 4 Aug 2026 | Maity et al., arXiv:2608.03944 | Very recent preprint | Claimed exact QEC–Quantum Darwinism tradeoff/no-go theorem |

The preprints above may be revised, rejected, split, or published after this cutoff. Conversely, papers published after 12 August 2026 may alter the novelty assessment. Every manuscript release should regenerate this table.

The exact expression “zero-record” is also not safely claimable as a coinage: an independent June 2026 PhilArchive manuscript uses “zero-record-holonomy” [@wu2026zerorecord]. That manuscript is not used here as scientific support, only as evidence of a terminology collision.

---

## 17. Recommended manuscript positioning

### 17.1 Primary question

> Under explicit temporal-intervention, access, and resource constraints, what is the minimum cost of certifying a quantum process with nontrivial causal memory, extracting only an authorized predicate, decoupling all other accessible history records, and recovering logical interference?

### 17.2 First defensible paper

A coherent first paper should contain:

1. an endpoint-equivalence/no-certification proposition;
2. a formal process-tensor and adversarial access model;
3. the joint task definition from Section 12;
4. either a soundness theorem, resource separation, or impossibility result;
5. a minimal multi-round reversible agent with late challenges;
6. ideal and noisy simulations including all Section 13 baselines;
7. an explicit comparison to quantum erasure, Wigner-friend memory erasure, memory witnesses, quantum ZK, QUALMs, and the 2026 direct preprints;
8. a claim ledger separating operational results from interpretations; and
9. a dated limitations and priority-search statement.

### 17.3 Safe one-sentence description

> We study whether a multi-time quantum process can be causally certified as using internal memory, queried for a restricted predicate, and then coherently reset so that no other accessible record of the hidden history remains while interference is recovered.

This description is technically ambitious without depending on a public-facing name, a many-worlds ontology, or a claim to observe an alternative universe.

---

## 18. Audit conclusion

The project is physically coherent at small circuit scale, but almost every intuitive ingredient is prior art. The novelty opportunity lies in a stringent conjunction, not in a new metaphor:

\[
\boxed{
\text{causal-memory certification}
+\text{predicate-only disclosure}
+\text{final-record decoupling}
+\text{recovered interference}
+\text{anti-shortcut access model}
}
\]

Even this conjunction must currently be described as **apparently underexplored**, not unprecedented. The 2026 literature sharply increases priority risk: direct branch-message circuits, record/recoherence theory, device-independent memory certification, superposition-secure zero knowledge, and a claimed QEC–Darwinism tradeoff are all active within months or days of the cutoff.

The project earns a distinct research identity only if it produces a theorem or protocol whose causal soundness and privacy/recoherence conditions cannot be compiled away into an ordinary phase oracle. If that target fails, the honest fallback is still useful: a rigorous benchmark suite for reversible multi-time memories and their loss of recoverability under controlled leakage.

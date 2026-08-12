# Originality gate: post-commitment audit-or-recovery protocols

**Search cutoff:** 12 August 2026
**Purpose:** adversarial novelty assessment, not a publication claim
**Status:** working research note; primary sources were used for technical conclusions
**Scope:** protocols in which a common multi-time quantum prefix creates a physical memory and a challenge sampled only after that prefix asks the device either to reveal a randomly selected past fact or to undo the prefix and certify recovered coherence plus removal of the distinguishing record.

## 1. Bottom-line verdict

The broad idea is **not original**. At least six mature literatures occupy most of its verbal description:

1. random delayed-choice quantum erasers already make a late choice between path information and conditional interference;
2. certified deletion, secure software leasing, and blind delegation with certified deletion already formalize a late use-or-delete choice, including computation on encoded data before certified deletion;
3. weak-measurement protocols already certify entanglement and then probabilistically recover it by reversal measurements;
4. mirror randomized benchmarking, Loschmidt echoes, and process-tensor witnesses already certify reversibility, temporal correlations, and quantum memory;
5. quantum secret sharing, private/correctable-channel duality, and operator-algebra quantum error correction already connect fragment access, absent information in a complement, and recovery; 2026 preprints now make the quantum-Darwinism/QEC connection explicit;
6. coherent classical communication, stabilizer/process verification, and quantum-instrument self-testing already identify and certify the coherent-copy isometry underlying the repaired two-port proposal;
7. operational wave--particle duality games already use a random late choice between extracting which-path information and demonstrating phase coherence, with tight multipath trade-off regions.

The closest cryptographic collision is particularly important. Bartusek *et al.* give a maliciously secure blind-delegation protocol in which an evaluator computes on protected data and later proves information-theoretic deletion of the protected software/data [Software with Certified Deletion, EUROCRYPT 2024](https://doi.org/10.1007/978-3-031-58737-5_4). Earlier work explicitly discusses coherent evaluation followed by reversal and deletion [Cryptography with Certified Deletion, CRYPTO 2023](https://doi.org/10.1007/978-3-031-38554-4_7). Therefore, **compute, uncompute, and certify deletion** cannot be claimed as a new combination.

The closest physical collision is Kim *et al.*, who use weak measurements to certify entanglement and reversal measurements to probabilistically recover the original entanglement [Science Advances 9, eadi5261 (2023)](https://doi.org/10.1126/sciadv.adi5261). Therefore, **certify a quantum property and subsequently recover entanglement** is also occupied.

The repaired one-step candidate is more strongly occupied than the multi-time flagship. Its ideal map \(V_{\mathrm{copy}}|j\rangle=|j\rangle_B|j\rangle_M\) is Harrow's coherent-bit (cobit) isometry [PRL 92, 097902 (2004)](https://doi.org/10.1103/PhysRevLett.92.097902). On half of a maximally entangled input, its Choi state is a GHZ state. A late computational-basis copy audit and a reverse-and-Bell echo are therefore complementary GHZ/stabilizer or channel-fidelity tests. With a trusted target decoder and a joint Bell-plus-reset check, the echo alone already tests the target Choi state, so the audit adds no ideal rigidity. **That repaired candidate is red unless it is reformulated in an untrusted-decoder model and yields a genuinely new robust theorem that is not a restatement of channel, stabilizer, gate, or instrument verification.**

The static diagonal/QND variant is also occupied at the architectural level. Bagan *et al.* define an \(N\)-path game in which a house flips a coin after a common path--detector interaction and requests either path discrimination (Ways) or phase discrimination (Phases), and derive a tight arbitrary-\(N\) region [PRL 120, 050402 (2018)](https://doi.org/10.1103/PhysRevLett.120.050402). Hillery later studies weaker games that request only partial path or phase information [J. Phys. A 54, 495301 (2021)](https://doi.org/10.1088/1751-8121/ac367d). Therefore, a late **audit-or-coherence** game is not itself new. Narrower mathematical gaps may remain for the exact frontiers between either a randomly requested Boolean coordinate or a committed random-prefix score and flag-conditioned entanglement return; Section 4.10 records why these are amber rather than green.

Within the primary-source corpus queried for this audit, I did **not** locate a protocol satisfying all of the following clauses together:

- one committed, multi-time prefix rather than two independently prepared experiments;
- verifier-chosen inputs or events injected throughout the prefix, so the history is not a precompiled label;
- a challenge generated only after the prefix is complete;
- an audit branch that queries a randomly selected earlier event and yields a quantitative causal-memory or memory-dimension certificate;
- an alternative recovery branch that implements a coherent inverse and certifies entanglement/interference recovery;
- a conditional decoupling statement for every residual record accessible to the device, allowing only a declared predicate or public transcript;
- one soundness theorem against direct-label, precomputation, measure-and-reprepare, postselection, and separate-device strategies.

That conjunction is a **plausible residual research task, not an established absence result**. Database search cannot prove that no equivalent construction exists, and several adjacent fields use incompatible terminology. A priority claim would still require expert review by researchers in quantum cryptography, process tensors, and quantum verification.

The project should therefore pass the originality gate only under a narrower description such as:

> A post-commitment, two-challenge game for a single multi-time quantum process, jointly certifying queryable causal memory in one branch and coherent recoverability with conditional record decoupling in the other.

It should not be described as the first reveal-or-erase experiment, the first reversible observer, the first erasure of a quantum record, the first computation followed by certified deletion, or the first recovery of entanglement after certification.

## 2. The exact object that would need to be new

Let a verifier interact sequentially with a device through a common prefix process \(\mathcal P\). The verifier samples fresh inputs \(X_1,\ldots,X_T\), and the device evolves an internal system and memory. Only after the last prefix interaction, the verifier samples a challenge \(C\).

### Audit challenge

For \(C=(A,J)\), with \(J\) sampled after the prefix, the device must answer a fact \(f_J(X_{1:T},Y_{1:T})\) about an earlier event. Acceptance must imply more than correlation with a branch-control bit. Under an explicit interface and dimension model, the score should lower-bound a memory resource: for example, effective environment dimension, process-tensor memory, or the success probability of a temporal random-access code.

### Recovery challenge

For \(C=E\), the device must implement or cooperate in implementing a coherent inverse. The verifier should test at least two logically distinct properties:

1. **recovery:** entanglement fidelity, an interference visibility, or a self-tested identity-channel fidelity relative to a verifier-held reference;
2. **record removal:** the residual device register \(Z\) is close to a state independent of the detailed history \(H\), conditional on the explicitly allowed output \(P\) and public challenge transcript \(T_{\mathrm{pub}}\):

\[
  \rho_{HZ\mid P,T_{\mathrm{pub}}}
  \approx
  \rho_{H\mid P,T_{\mathrm{pub}}}\otimes
  \rho_{Z\mid P,T_{\mathrm{pub}}}.
\]

Unconditional decoupling is usually the wrong target because the public predicate, challenge, acceptance bit, and verifier transcript can legitimately remain correlated with the run.

### Cross-branch soundness

The novelty-bearing theorem cannot be two independent completeness statements. It must show that the *same committed process* could not know in advance which branch will be requested, and that high average score excludes a specified comparison class. Candidate comparison classes include:

- a direct phase gate or direct access to the history label;
- a device storing only a precomputed predicate;
- a classical finite-state memory of bounded dimension;
- a measure-and-reprepare channel;
- a strategy that leaves a hidden transcript outside the declared inverse;
- a postselected weak-measurement strategy whose success probability is not charged;
- two separately calibrated devices, one optimized for audit and one for echo;
- a simulator given the entire future challenge in advance.

The challenge bit is not cosmetic. Without post-prefix unpredictability and a commitment to one physical process, the experiment reduces to running one memory test and one echo benchmark side by side.

## 3. Collision matrix

Codes: **P** = common prefix before choice; **A** = late fact/path audit; **E** = deletion or inverse; **R** = physical coherence/entanglement recovery; **D** = residual-state deletion/decoupling guarantee; **M** = causal or dimension-bounded memory soundness. “Partial” means that the source realizes a weaker or differently modelled version.

| Primary source | P | A | E | R | D | M | Originality consequence |
|---|---:|---:|---:|---:|---:|---:|---|
| [Kim *et al.*, delayed-choice quantum eraser (2000)](https://doi.org/10.1103/PhysRevLett.84.1) | yes | path | basis erasure | conditional interference | no | no | Occupies the delayed reveal/erase headline. |
| [Scarcelli, Zhou, and Shih (2007)](https://doi.org/10.1140/epjd/e2007-00164-y) | yes | random late path readout | basis erasure | coincidence interference | no | no | An explicit random late choice between reading path information and erasing it. |
| [Ma *et al.* (2013)](https://doi.org/10.1073/pnas.1213201110) | yes | path information | eraser basis | conditional interference | no | no | The choice is causally disconnected; lateness alone is not novel. |
| [Elouard *et al.* (2021)](https://doi.org/10.22331/q-2021-07-08-498) | yes | friend/lab context | friend-memory erasure | interferometric analogue | no | no | Explicitly occupies “erasing the memory of Wigner’s friend.” |
| [Katz *et al.* (2008)](https://doi.org/10.1103/PhysRevLett.101.200401) | yes | partial measurement record | conditional uncollapse | state recovery | no | no | Measurement reversal is established and probabilistic success must be charged. |
| [Kim *et al.* (2023)](https://doi.org/10.1126/sciadv.adi5261) | sequential | entanglement certification | reversal measurement | probabilistic entanglement recovery | no | no | Occupies certification followed by recovery and its trade-off curves. |
| [Broadbent and Islam (2020)](https://doi.org/10.1007/978-3-030-64381-2_4) | ciphertext | decrypt/use | certified deletion | no | plaintext hiding | no | Exact logical use-or-delete primitive; deletion is cryptographic, not recoherence. |
| [Poremba (2023)](https://doi.org/10.4230/LIPIcs.ITCS.2023.90) | encrypted input | computed output | FHE proof of deletion | no | plaintext hiding | no | Computation plus certified deletion is occupied. |
| [Bartusek *et al.* (2024)](https://doi.org/10.1007/978-3-031-58737-5_4) | delegated computation | repeated functionality | maliciously secure deletion | no | information-theoretic deletion | no | The strongest collision with compute-then-delete. |
| [Çakan, Goyal, and Raizes (2024 preprint)](https://arxiv.org/abs/2411.05176) | cryptographic object | verification/use | certified deniability | no | simulation-based “no trace” | no | Blocks novelty claims based only on leaving no usable residual evidence. |
| [Mahadev (2018)](https://doi.org/10.1109/FOCS.2018.00033) | cryptographic commitment | test/standard-basis branch | no | Hadamard-basis computational branch, not recovery | computational soundness | no | Occupies post-commitment random test-versus-function challenge architecture. |
| [Gunn *et al.* (2025)](https://doi.org/10.1145/3717823.3718264) | basis-independent classical commitment to a quantum state | late standard-basis opening | destructive opening | late Hadamard-basis opening, not recovery | computational binding | no | Occupies a single commitment opened later in either complementary basis. |
| [Harrow (2004)](https://doi.org/10.1103/PhysRevLett.92.097902) | coherent-copy isometry | either output carries the basis label | inverse isometry | preserves superpositions/creates entanglement | no | no | The target \(|j\rangle\mapsto|j\rangle|j\rangle\) is the established cobit primitive. |
| [Hofmann (2005)](https://doi.org/10.1103/PhysRevLett.94.160504) | one target process | classical truth table in one basis | no | complementary-basis truth table | no | process-fidelity bound | Two complementary classical fidelities bound the quantum process fidelity; this is the closest score-level collision. |
| [Okamoto *et al.* (2005)](https://doi.org/10.1103/PhysRevLett.95.210506) | optical CNOT | computational truth table | no | complementary truth table | no | process/entangling fidelity | Experimentally estimates CNOT process fidelity and entangling capability from two classical truth tables. |
| [Mohan, Tavakoli, and Brunner (2019)](https://doi.org/10.1088/1367-2630/ab3773) | one prepare-transform-measure device | first sequential QRAC score | no | second preserved-state QRAC score | no | instrument rigidity | Gives an optimal joint frontier and self-tests an information-gain/disturbance instrument. |
| [Bagan *et al.* (2018)](https://doi.org/10.1103/PhysRevLett.120.050402) | common path--detector state | late full-path discrimination | no | late phase discrimination / coherence score | no | no | Exact late Ways-or-Phases architecture and tight arbitrary-\(N\) region; kills novelty of the static game itself. |
| [Hillery (2021)](https://doi.org/10.1088/1751-8121/ac367d) | common path--detector state | partial/coarse path answer | no | partial/coarse phase answer | no | no | Partial-information duality games are occupied, although not the exact random-coordinate Bayes score. |
| [Wagner *et al.* (2020)](https://doi.org/10.22331/q-2020-03-19-243) | unknown measurement instrument | classical outcome | no | post-measurement entanglement tests | no | instrument rigidity | Robustly self-tests quantum instruments from pre/post-measurement correlations. |
| [Dangniam, Han, and Zhu (2020)](https://doi.org/10.1103/PhysRevResearch.2.043323) | target stabilizer state | \(Z\)-type stabilizers | no | \(X\)-type coherence stabilizer | no | state fidelity | The qubit copy-audit/echo observables are a GHZ stabilizer verification test. |
| [Hsieh *et al.* (2026)](https://doi.org/10.1103/dx6d-4kjy) | interactive instrument | recovery of generated classical information | outcome-conditioned recovery | maximally-entangled-state restoration | no | resource characterization | Peer-reviewed 2026 result makes classical recovery and entanglement restoration operational meanings of one instrument resource. |
| [Sarkar (2026)](https://doi.org/10.1103/m1tx-9mx1) | network-tested operation | no | no | arbitrary unitary self-test | no | operation rigidity | Peer-reviewed 2026 theorem self-tests any unitary in a quantum network. |
| [Kimmel and Kolkowitz (2019)](https://doi.org/10.1103/PhysRevA.100.052326) | sealed state | read message | return for tamper test | state-return test, not inverse | disturbance bound | no | Read-versus-return verification and its no-go bounds are occupied. |
| [Giarmatzi and Costa (2021)](https://doi.org/10.22331/q-2021-04-26-440) | multi-time process | temporal witness | no | no | no | quantum-vs-classical | Generic quantum-memory witnessing is occupied. |
| [Vieira *et al.* (2024)](https://doi.org/10.22331/q-2024-01-10-1224) | multi-time process | random temporal statistics | no | no | no | environment dimension | Random past/future correlations can already lower-bound effective memory dimension. |
| [White *et al.* (2025)](https://doi.org/10.22331/q-2025-04-08-1695) | multi-time process | unitary-only witness | no | terminal measurement | no | temporal entanglement | Even unitary control plus one final measurement can witness multi-time properties. |
| [Proctor *et al.* (2022)](https://doi.org/10.1103/PhysRevLett.129.150502) | random circuit | no | inverse mirror | circuit polarization/fidelity | no | no | Run-and-invert as a scalable benchmark is occupied. |
| [Laeuchli and Trujillo-Rasua (2024)](https://doi.org/10.1007/s11128-024-04421-x) | memory state | random-address checksum challenge | no | entanglement-assisted attestation | no | classical memory | Random challenge attestation is adjacent, with strong access-model assumptions. |

No row in this table, by itself, instantiates all six columns. Several rows occupy five or six words of the proposed title, however, so novelty can reside only in a formal joint task and theorem.

## 4. Findings by literature sector

### 4.1 Delayed-choice erasers: the public story is occupied

Kim *et al.* implemented a delayed-choice eraser in which the which-path or both-path measurement of an entangled partner is selected after registration of the signal photon [PRL 84, 1 (2000)](https://doi.org/10.1103/PhysRevLett.84.1). Scarcelli, Zhou, and Shih made a random delayed choice after signal detection between reading complete path information and erasing it [EPJ D 44, 167–173 (2007)](https://doi.org/10.1140/epjd/e2007-00164-y), [arXiv:quant-ph/0512207](https://arxiv.org/abs/quant-ph/0512207). Ma *et al.* placed the eraser choice and the interferometer events in causally disconnected regions [PNAS 110, 1221–1226 (2013)](https://doi.org/10.1073/pnas.1213201110).

These experiments sort joint outcomes according to complementary measurements. They do not implement a deterministic inverse of a multi-step memory-bearing computation, certify that all work registers returned to a reference state, or prove a causal-memory lower bound. That distinction is technically meaningful, but it must be stated. “The later choice reveals a past fact or restores interference” is already a description of prior quantum erasers.

A November 2025 preprint proposes a memory-delayed eraser discriminator without postselection [Ohwada, arXiv:2511.22827](https://arxiv.org/abs/2511.22827). It was not located as a peer-reviewed result in this audit and its claims should be treated as provisional, but it increases priority risk for any novelty claim based on a quantum memory making the choice genuinely late.

### 4.2 Friend-memory erasure and measurement reversal: record creation and undo are occupied

Elouard *et al.* explicitly model the quantum erasure of Wigner’s friend’s memory and identify incompatible experimental contexts [Quantum 5, 498 (2021)](https://doi.org/10.22331/q-2021-07-08-498). This blocks novelty based merely on treating a memory-bearing observer as a controlled subsystem and undoing its record.

Partial quantum measurements can be conditionally reversed. Katz *et al.* experimentally erased the information extracted by a weak partial-collapse measurement and recovered a superconducting-qubit state [PRL 101, 200401 (2008)](https://doi.org/10.1103/PhysRevLett.101.200401). Kim *et al.* went closer to the proposed certification language: weak measurements certify entanglement, after which reversal measurements probabilistically recover the initial entanglement [Science Advances 9, eadi5261 (2023)](https://doi.org/10.1126/sciadv.adi5261), [arXiv:2305.06852](https://arxiv.org/abs/2305.06852).

The residual distinction is not “measurement plus reversal” or “certification plus recovery.” It is a post-prefix random choice between a destructive causal-history audit and an alternative coherent inverse, with a single adversarial soundness model and no uncharged postselection.

### 4.3 Certified deletion and leasing: the closest logical collision

Broadbent and Islam formalized encryption with certified deletion: a receiver may retain a quantum ciphertext for later decryption or produce a classical certificate after which the plaintext remains hidden even if the key is later released [TCC 2020](https://doi.org/10.1007/978-3-030-64381-2_4), [arXiv:1910.03551](https://arxiv.org/abs/1910.03551). Kundu and Tan gave a composable device-independent construction based on nonlocal games [Quantum 7, 1047 (2023)](https://doi.org/10.22331/q-2023-07-06-1047).

Computation before deletion is no longer a gap. Poremba combined fully homomorphic encryption with a proof of deletion [ITCS 2023](https://doi.org/10.4230/LIPIcs.ITCS.2023.90). Bartusek and Khurana developed a general certified-deletion framework [CRYPTO 2023](https://doi.org/10.1007/978-3-031-38554-4_7), and Bartusek *et al.* subsequently supplied the first fully maliciously secure blind delegation with certified deletion, correcting weaker prior security models [EUROCRYPT 2024](https://doi.org/10.1007/978-3-031-58737-5_4), [IACR ePrint 2023/265](https://eprint.iacr.org/2023/265).

Secure software leasing likewise lets a user evaluate a quantum-encoded program and later return it in a way that prevents future use [Broadbent *et al.*, arXiv:2101.12739](https://arxiv.org/abs/2101.12739); [Kitagawa, Nishimaki, and Yamakawa, arXiv:2010.11186](https://arxiv.org/abs/2010.11186). Certified deniability goes beyond inability to recover a secret: the adversary’s residual state after accepted deletion must be simulable without having received the deletable object [Çakan, Goyal, and Raizes, arXiv:2411.05176](https://arxiv.org/abs/2411.05176). As of the cutoff it was a preprint accepted for a QCrypt 2026 talk, not a located journal or archival proceedings article.

Cryptographic deletion and physical recoherence are different resources:

- a deletion certificate normally proves an inability to recover a protected message within a security game;
- an echo or Bell-pair test proves that a specified channel approximately preserved or restored quantum correlations;
- a physical decoupling statement concerns the full residual state under a declared system boundary;
- none of these, alone, proves that a queried fact was generated by a nontrivial causal history rather than encoded directly.

The new task would need to relate these notions quantitatively. Merely relabeling certified deletion as record erasure would not be new.

### 4.4 Test-versus-compute and random-access challenges: late challenge architecture is occupied

Mahadev’s measurement protocol makes a quantum prover produce a cryptographic commitment before a classical verifier selects a test or Hadamard round; it enforces standard- or Hadamard-basis measurement behavior under the Learning with Errors assumption [FOCS 2018](https://doi.org/10.1109/FOCS.2018.00033), [arXiv:1804.01082](https://arxiv.org/abs/1804.01082). Brakerski *et al.* use the same trapdoor-claw-free paradigm for a preimage-versus-equation challenge that certifies quantumness and randomness of one untrusted device [FOCS 2018](https://doi.org/10.1109/FOCS.2018.00038), [JACM 68, 31 (2021)](https://doi.org/10.1145/3441309). Thus, a commitment followed by an unpredictable classical-versus-phase-sensitive challenge is standard architecture in quantum verification.

Gunn *et al.* sharpen the collision: their basis-independent classical commitment fixes a quantum state before the receiver chooses whether to open any qubit in the standard or Hadamard basis [STOC 2025](https://doi.org/10.1145/3717823.3718264), [arXiv:2404.14438v2](https://arxiv.org/abs/2404.14438). This matters because Mahadev's original key generation can depend on the requested measurement basis, whereas the later commitment construction explicitly removes that syntactic escape.

The repaired physical game is **not literally an unkeyed instance** of those cryptographic protocols. It gives the verifier a quantum reference and trusted quantum tests, seeks information-theoretic channel fidelity, and has no classical commitment string, trapdoor, or computational assumption. Conversely, the cryptographic protocols give a classical verifier computational binding/extraction statements, not a physical inverse, Bell-pair return, or reset of the prover's memory. Nevertheless, the commit-then-complementary-challenge skeleton and the claim that both branches bind one underlying quantum object are occupied. The physical protocol must earn novelty from a new channel/instrument rigidity statement, not from the late challenge itself.

Random-access codes already ask for a bit selected only after an encoding has been prepared. Ambainis *et al.* give quantum encodings from which any requested source bit can be recovered with bounded success and prove size limitations [JACM 49, 496–511 (2002)](https://doi.org/10.1145/581771.581773), [arXiv:quant-ph/9804043](https://arxiv.org/abs/quant-ph/9804043). A historical-fact audit is therefore a temporal random-access task unless it exploits and certifies additional causal structure.

A genuine contribution would have to combine post-commitment challenge soundness with a physical recovery/decoupling certificate, not merely borrow the narrative of cut-and-choose.

### 4.5 Process tensors and temporal memory witnesses: causal memory is measurable already

Giarmatzi and Costa map quantum memory in a non-Markovian process to an entanglement-witness problem [Quantum 5, 440 (2021)](https://doi.org/10.22331/q-2021-04-26-440). Taranto, Pollock, and Modi define operational memory strength for quantum stochastic processes and bound process recoverability from it [npj Quantum Information 7, 149 (2021)](https://doi.org/10.1038/s41534-021-00481-4). Any claim of a generic new “memory versus recoverability” inequality must be compared directly with this result.

Vieira *et al.* use temporal correlations and a convergent SDP hierarchy to lower-bound the effective environment dimension [Quantum 8, 1224 (2024)](https://doi.org/10.22331/q-2024-01-10-1224). White *et al.* show that unitary sequences followed only by a terminal measurement can witness genuine multi-time entanglement and bound non-Markovianity and related quantities [Quantum 9, 1695 (2025)](https://doi.org/10.22331/q-2025-04-08-1695). This is a close collision with any claim that reversible controls and one final measurement are a new way to interrogate a multi-time process.

Rosset, Buscemi, and Liang provide faithful verification games for non-entanglement-breaking quantum memories with minimal trust assumptions [PRX 8, 021033 (2018)](https://doi.org/10.1103/PhysRevX.8.021033). Sekatski *et al.* develop robust device-independent self-testing of a qubit identity channel [PRL 131, 170802 (2023)](https://doi.org/10.1103/PhysRevLett.131.170802). These supply mature choices for the recovery branch.

There is a subtle modelling hazard. Vieira, Ku, and Budroni show that entanglement-breaking channels can still generate genuinely nonclassical multi-time correlations and therefore are not, without further assumptions, synonymous with classical temporal memory [Physical Review Research 7, 043281 (2025)](https://doi.org/10.1103/r8lf-bb4p). A soundness theorem must define the classical comparison class at the process level, not identify “classical memory” with an entanglement-breaking channel at one time step.

Ohst *et al.* characterize the value of bounded classical and quantum memory in adaptive channel-discrimination strategies [Quantum 10, 1988 (2026)](https://doi.org/10.22331/q-2026-01-28-1988). This is a peer-reviewed 2026 result and directly threatens novelty based on memory-dimension hierarchies for adaptive tasks.

### 4.6 Echo and randomized benchmarking: inversion as a certificate is occupied

Randomized mirror circuits execute a circuit followed by a structured inverse and infer scalable error metrics. Proctor *et al.* introduced scalable mirror-circuit randomized benchmarking [PRL 129, 150502 (2022)](https://doi.org/10.1103/PhysRevLett.129.150502); Hines *et al.* demonstrated scalable randomized benchmarking of universal gate sets [PRX 13, 041030 (2023)](https://doi.org/10.1103/PhysRevX.13.041030). Peters *et al.* used Loschmidt echoes as a program-specific reversibility diagnostic on superconducting processors [PRX Quantum 3, 040333 (2022)](https://doi.org/10.1103/PRXQuantum.3.040333). Shaffer *et al.* developed time-reversal, multi-basis, and randomized inversion protocols for analog simulators [npj Quantum Information 7, 46 (2021)](https://doi.org/10.1038/s41534-021-00380-8).

Consequently, an inverse circuit followed by return probability is not an originality anchor. The recovery branch should instead certify a resource not implied by a basis-state echo alone, ideally entanglement fidelity relative to a reference plus an explicit conditional record-decoupling bound.

A June 2026 preprint augments mirror benchmarking with injected Bell pairs and tracks per-edge mutual information while retaining a global infidelity estimate [Pettugani *et al.*, arXiv:2606.20123](https://arxiv.org/abs/2606.20123). It is provisional, but it narrows the gap around combining mirror circuits with entanglement/correlation diagnostics.

### 4.7 Quantum seals and remote memory attestation: read-versus-return has bounds

Quantum sealing asks that a message be readable while later testing whether it was read. Kimmel and Kolkowitz derive no-go bounds between readability and the probability of detecting that the seal was broken [PRA 100, 052326 (2019)](https://doi.org/10.1103/PhysRevA.100.052326). This is structurally close to “answer a query or return an undisturbed quantum object” and should be treated as a baseline, not merely an analogy.

Quantum-assisted remote memory attestation uses unpredictable challenges and entanglement to authenticate a prover and checksum classical memory [Laeuchli and Trujillo-Rasua, Quantum Information Processing 23, 208 (2024)](https://doi.org/10.1007/s11128-024-04421-x). Their follow-up identifies correctness and security flaws in another proposed quantum-attestation scheme and proves limitations for attesting quantum memory, including the problem that an unknown memory cannot simply be copied for verification [IET Quantum Communication (2025)](https://doi.org/10.1049/qtc2.70019), [arXiv:2503.04311](https://arxiv.org/abs/2503.04311).

This literature offers a useful adversarial lesson: a memory challenge is meaningful only after the trust boundary, communication abilities, timing, device isolation, and ability to send quantum states to collaborators are formalized.

### 4.8 Coalition control and record-access hypergraphs: the exact core is occupied

The proposed coalition formulation is useful notation, but its exact mathematical core is not a defensible novelty claim. Let the environment fragments be indexed by \(V=[n]\). If a family \(\mathcal I\subseteq 2^V\) contains every coalition whose reduced state has any dependence on the history label, then \(\mathcal I\) is upward closed by data processing. Let \(\mathcal H_{\min}\) be its inclusion-minimal members. An uncontrolled set \(U=V\setminus C\) is label-decoupled exactly when it contains no member of \(\mathcal H_{\min}\). Set-theoretically,

\[
  U\text{ is label-decoupled}
  \quad\Longleftrightarrow\quad
  C\cap A\neq\varnothing\ \text{for every }A\in\mathcal H_{\min}.
\]

Thus the minimum cardinality of a controlling coalition is the transversal number \(\tau(\mathcal H_{\min})\). This is the ordinary blocker/hitting-set dual of an upward-closed access structure. Calling it a “recoherence number” may be convenient, but the equality itself is a combinatorial restatement once the informative family is defined.

The quantum-information content behind the restatement is also established:

- Cleve, Gottesman, and Lo prove the complement duality of authorized and unauthorized sets for pure-state quantum secret sharing and connect threshold sharing to erasure correction [PRL 83, 648–651 (1999)](https://doi.org/10.1103/PhysRevLett.83.648).
- Gottesman develops general monotone quantum access structures [PRA 61, 042311 (2000)](https://doi.org/10.1103/PhysRevA.61.042311).
- Kretschmann, Kribs, and Spekkens prove that a subsystem is private for a channel precisely when it is correctable for a complementary channel, including approximate diamond-norm versions [PRA 78, 032330 (2008)](https://doi.org/10.1103/PhysRevA.78.032330).
- Ouyang *et al.* prove an approximate quantum-secret-sharing equivalence between reconstructability by a coalition and small leakage through the complementary channel [PRA 108, 012425 (2023)](https://doi.org/10.1103/PhysRevA.108.012425).
- Girard, Cheng, and Cao explicitly identify Darwinist records with locally recoverable commuting logical subalgebras of quantum codes and classify which logical operators are supported on which fragments [arXiv:2606.06588v2](https://arxiv.org/abs/2606.06588). This June 2026 preprint is provisional, but it is the closest collision with the proposed Darwinism–QEC access-structure bridge.

Two older physical literatures further occupy “how much of the environment must be controlled.” Miatto *et al.* optimize the coherence recoverable on a qubit when only a fraction of its environment can be measured, finding a half-environment transition for Haar-random pure states [PRA 92, 062331 (2015)](https://doi.org/10.1103/PhysRevA.92.062331). Environment-assisted restoration with measurement and weak-measurement reversal is also explicit [Wang, Zhao, and Yu, PRA 89, 042320 (2014)](https://doi.org/10.1103/PhysRevA.89.042320). Localizable coherence formalizes the coherence obtainable in one subsystem by measuring or discarding the rest [Hamma, Styliaris, and Zanardi, Physics Letters A 397, 127264 (2021)](https://doi.org/10.1016/j.physleta.2021.127264).

The threshold examples follow, but only under carefully stated models:

- For a **perfect classical** \((k,n)\) record-sharing structure, coalitions of size below \(k\) have zero information. Requiring \(|V\setminus C|<k\) gives \(|C|\ge n-k+1\). This is an immediate access-structure corollary.
- For an ideal spectrum-broadcast/repetition-code record, every singleton is informative, so the minimal hyperedges are \(\{1\},\ldots,\{n\}\) and \(\tau=n\). Exact unconditional interference requires controlling every record-bearing fragment. Spectrum-broadcast states and redundant environmental records are established constructions [Korbicz, Horodecki, and Horodecki, PRL 112, 120402 (2014)](https://doi.org/10.1103/PhysRevLett.112.120402); [Blume-Kohout and Zurek, PRA 73, 062310 (2006)](https://doi.org/10.1103/PhysRevA.73.062310).
- For a **pure quantum** threshold scheme, complement duality forces \(n=2k-1\), so the same expression becomes \(n-k+1=k\). Mixed-state threshold schemes hide a purifying share outside the listed \(n\) systems; omitting that share from the control boundary invalidates an exact recoherence conclusion.

There are three important failure modes in the proposed iff statement.

1. “Unauthorized” must mean **forbidden/decoupled**, not merely unable to reconstruct the record perfectly. Ramp or noisy coalitions can retain partial which-history information and still suppress exact interference.
2. Equality of the complement’s reduced states for the computational history labels is only classical privacy. Exact recovery of an *arbitrary coherent logical input* requires the complementary channel to be constant on the full relevant logical algebra, including off-diagonal operators. This is the private/correctable-channel condition, not merely a classical secret-sharing condition.
3. The allowed control class matters. Unitary inversion, measurement plus feed-forward, heralded erasure, and postselected localization have different control costs. The half-environment result for random states already shows that measurement-assisted coherence recovery need not obey a naïve “control every microscopic degree of freedom” count.

Recent work makes this territory especially high-risk. Girard *et al.* connect Quantum Darwinism directly to operator-algebra QEC [arXiv:2606.06588](https://arxiv.org/abs/2606.06588). Maity, Onggadinata, and Koh give a model-specific exact Darwinism–logical-fidelity tradeoff plus an information-theoretic no-go bound [arXiv:2608.03944](https://arxiv.org/abs/2608.03944). Torvinen, Keski-Vakkuri, and Pranzini study Petz recovery of einselected information from environment fragments [arXiv:2605.06848](https://arxiv.org/abs/2605.06848). Strasberg *et al.* analyze approximate records and recoherence in decoherent histories [arXiv:2601.19703](https://arxiv.org/abs/2601.19703). All four are 2026 preprints at the cutoff and must be treated as provisional; together they eliminate any safe claim that a Quantum-Darwinism/QEC/recovery connection is itself new.

**Residual opportunity.** A potentially original result would be an *approximate, weighted coalition-control frontier* for a declared operation class: given noisy or ramp conditional fragment states, heterogeneous control costs, and correlated fragments, optimize the achievable entanglement fidelity and residual history leakage. It would need a theorem or algorithm stronger than applying approximate secret-sharing/complementary-channel recovery coalition by coalition, and it should reproduce the perfect-threshold, spectrum-broadcast, and random-state half-environment cases as baselines.

**Hard kill criterion.** Kill this candidate as a foundational theorem if its proof is only (i) private/correctable complementarity or Uhlmann/Stinespring recovery, followed by (ii) the elementary hitting-set identity. Also kill it if it treats “not authorized” as “zero information,” certifies only reset of named fragments while ignoring a purification, or guarantees recovery only for one known initial superposition rather than a reference-entangled logical input.

### 4.9 Two-port copy audit versus reverse-and-Bell echo: the one-step candidate is occupied

#### Exact repaired candidate

The repaired proposal gives a device an orthogonal label \(j\) and asks it to realize a nondemolition coherent-copy map

\[
  V_{\mathrm{copy}}:\ |j\rangle_Q\longmapsto |j\rangle_B|j\rangle_M .
\]

Only after this prefix, a random challenge requests either:

1. a **copy audit**, in which both output ports \(B\) and \(M\) must reproduce \(j\); or
2. an **echo**, in which a decoder reverses the prefix, the input is returned entangled with a verifier-held reference, and \(M\) is reset or at least decoupled.

Identity and SWAP strategies send the label to only one audited port and therefore fail the two-port audit. A computational-basis measure-and-prepare channel can pass that audit, but it is entanglement breaking; no decoder can make its Choi state maximally entangled. For a \(d\)-dimensional input, a separable normalized Choi state has overlap at most \(1/d\) with a maximally entangled target, so the echo rejects the ideal entanglement-breaking cheat.

The narrative sounds distinctive. The mathematical object is not. Harrow defined exactly this basis-copying isometry as a coherent classical bit, or **cobit**, two decades ago [PRL 92, 097902 (2004)](https://doi.org/10.1103/PhysRevLett.92.097902). The task is therefore certification of a cobit/premeasurement isometry, not a new kind of quantum observer.

#### Source replacement exposes a GHZ/Choi test

Prepare

\[
  |\Phi_d\rangle_{RQ}
  =\frac{1}{\sqrt d}\sum_{j=0}^{d-1}|j\rangle_R|j\rangle_Q
\]

before the challenge and apply the common prefix only to \(Q\). Measuring \(R\) in the computational basis *after* the prefix generates the hidden audit label \(j\); leaving \(R\) coherent supports the echo. The post-prefix state of the honest device is

\[
  |\Gamma_V\rangle_{RBM}
  =(\mathbb I_R\otimes V_{\mathrm{copy}})|\Phi_d\rangle
  =\frac{1}{\sqrt d}\sum_j |j,j,j\rangle .
\]

For \(d=2\), this is the three-qubit GHZ state. The audit accepts the \(Z\)-type correlation subspace

\[
  \Pi_{\mathrm{audit}}
  =|000\rangle\!\langle000|+|111\rangle\!\langle111|,
\]

equivalently checking two independent \(Z\)-parity stabilizers. The coherent part of the echo checks the \(X_R X_B X_M\) stabilizer. With the ideal inverse CNOT from \(B\) to \(M\), a Bell test on \(RB\) plus a reset test on \(M\) pulls back to the full GHZ target. Random \(Z\)-correlation and \(X\)-coherence settings are standard stabilizer verification; optimal Pauli verification of GHZ and other stabilizer states is already developed [Dangniam, Han, and Zhu, Physical Review Research 2, 043323 (2020)](https://doi.org/10.1103/PhysRevResearch.2.043323).

This equivalence has a sharp consequence. If the source, target decoder, and Bell-plus-reset measurement are trusted, the echo probability is already the Choi-state fidelity with \(V_{\mathrm{copy}}\). At unit score it identifies the target channel on the tested input space; the audit is redundant for rigidity. If the echo omits the reset test, the audit can add a missing \(Z\)-parity constraint, but the combined test is still an ordinary GHZ stabilizer witness.

#### It is also a restricted CNOT process-fidelity test

With a blank target \(M=|0\rangle\), \(V_{\mathrm{copy}}\) is the restriction of CNOT to that input subspace. Hofmann proved that two complementary classical fidelities \(F_1,F_2\) bound target process fidelity as

\[
  F_1+F_2-1\le F_{\mathrm{proc}}\le\min(F_1,F_2)
\]

[PRL 94, 160504 (2005)](https://doi.org/10.1103/PhysRevLett.94.160504). Okamoto *et al.* then estimated an optical CNOT's process fidelity and entangling capability from two complementary truth tables [PRL 95, 210506 (2005)](https://doi.org/10.1103/PhysRevLett.95.210506). Modern quantum-process verification explicitly reduces process tests to verification of Choi states or prepare-and-measure output tests [Liu *et al.*, PRA 101, 042315 (2020)](https://doi.org/10.1103/PhysRevA.101.042315); efficient local verification protocols cover generalized CNOT and other Clifford gates [Zhu and Zhang, PRA 101, 042316 (2020)](https://doi.org/10.1103/PhysRevA.101.042316).

The proposed audit and echo are not numerically identical to Hofmann's full two-basis CNOT truth tables: they test the blank-target isometric subspace, and one branch uses a reference-entangled input. But this is a restriction and repackaging of process/Choi verification, not a new certification principle.

#### Instrument and self-testing collisions

The nearest results under weaker trust assumptions are also mature:

| Primary source | What it certifies | Relation to the repaired game |
|---|---|---|
| [Mohan, Tavakoli, and Brunner, NJP 21, 083034 (2019)](https://doi.org/10.1088/1367-2630/ab3773) | An optimal pair of sequential random-access-code scores, including a self-test of a qubit measurement instrument and a tight information-gain/disturbance frontier | Already supplies a two-score instrument-rigidity theorem in a prepare-transform-measure model. |
| [Wagner *et al.*, Quantum 4, 243 (2020)](https://doi.org/10.22331/q-2020-03-19-243) | Robust device-independent characterization of an instrument from the input state, outcome probabilities, and conditional post-measurement states | Occupies certification of pre/post-measurement behavior, although its explicit output register is classical rather than a coherent second port. |
| [Sekatski *et al.*, PRL 121, 180505 (2018)](https://doi.org/10.1103/PhysRevLett.121.180505) | Device-independent certification of coherent storage, processing, and transfer operations | General operation-certification baseline. |
| [Hsieh *et al.*, PRA 113, 062445 (2026)](https://doi.org/10.1103/dx6d-4kjy) | A resource theory of interactive instruments whose robustness equals, among other tasks, optimal maximally-entangled-state restoration and classical-information recovery | The closest resource-level collision between the audit and echo capabilities. This is peer reviewed, published 17 June 2026. |
| [Sarkar, PRL 137, 030802 (2026)](https://doi.org/10.1103/m1tx-9mx1) | Device-independent self-testing of any unitary gate in a network with independent sources | Removes any broad 2026 claim that coherent inverse operations or CNOT-class gates lack device-independent certification. Published 15 July 2026. |
| [Paul, Roy, and Pan, arXiv:2604.19911](https://arxiv.org/abs/2604.19911) | Semi-device-independent self-testing of unitary operations via a prepare-measure random-access game | Provisional 21 April 2026 preprint; directly adjacent to an untrusted-operation, random-access formulation. |

Hsieh *et al.* do not present the same two quantum output ports, late challenge, or cobit rigidity score. Sarkar tests unitary gates, whereas \(V_{\mathrm{copy}}\) is an input-output isometry unless a blank ancilla and its full unitary extension are included. These are real modelling differences, but they leave very little room for a novelty claim based only on the one-step functionality.

#### Comparison with Mahadev-style test/Hadamard rounds

The cryptographic family has the same high-level commitment pattern:

- Mahadev commits a quantum prover before a random test-versus-Hadamard challenge and proves computational soundness from LWE [FOCS 2018](https://doi.org/10.1109/FOCS.2018.00033).
- Brakerski, Christiano, Mahadev, Vazirani, and Vidick use a preimage-versus-equation challenge to certify quantumness of one untrusted device [JACM 68, 31 (2021)](https://doi.org/10.1145/3441309).
- Gunn, Kalai, Natarajan, and Villányi give a basis-independent classical commitment that can later be opened in either the standard or Hadamard basis [STOC 2025](https://doi.org/10.1145/3717823.3718264).

The physical cobit test is not merely an unkeyed special case in a theorem-preserving sense. It has a quantum verifier/reference and information-theoretic channel metrics; the cryptographic tests have a classical verifier, classical commitment, computational binding, and no physical Bell-return or memory-reset guarantee. What is already occupied is the **architectural claim** that a late complementary challenge binds one earlier quantum commitment. What remains to be proved, if anything, is a distinct physical rigidity theorem under a clearly weaker trust model.

#### Trust-model verdict

| Model | What the pair establishes | Originality status |
|---|---|---|
| Trusted source, fixed target decoder, trusted Bell-plus-reset readout | Choi/GHZ fidelity of a known cobit target; the full echo already contains the audit constraints | **Killed:** standard state/process verification. |
| Trusted source and finite ports, unknown prefix, decoder optimized by the verifier | High echo means the prefix is correctable; perfect two-port audit fixes its basis action | **Essentially occupied:** the exact statement follows from elementary Stinespring/Choi arguments; only a sharp robust bound might be publishable. |
| Prefix and decoder supplied by one untrusted device | Echo alone permits identity/SWAP or another reversible encoding; audit selects a basis-copy encoding | **Possible narrow gap**, but only with fixed dimensions, one committed device, explicit no-leakage boundary, and a new robust self-test beyond existing instrument/operation certification. |
| Unbounded hidden environment or undeclared side channel | Named ports may pass while a transcript leaks elsewhere | **No global erasure claim is sound** without cryptographic or physical assumptions. |

At perfect score in the finite-dimensional trusted-port model, the rigidity proof is short. A Stinespring dilation maps each basis input to \(|j,j\rangle|e_j\rangle\) because the two-port audit is perfect. Perfect entanglement recovery forces the complementary channel to be input independent, hence all \(|e_j\rangle\) coincide up to an irrelevant phase and the off-diagonal coherences are preserved. This is the cobit isometry. A robust version obtained only by continuity of Stinespring dilation, Fuchs--van de Graaf inequalities, or existing process-fidelity bounds is not a new foundational theorem.

**Exact-collision finding.** I did not locate a primary paper using the exact prose “sample a late bit; either measure two coherent-copy output ports in \(Z\), or apply a decoder and test Bell return plus memory reset.” That wording-level absence does not rescue originality because the trusted score is unitarily equivalent to GHZ/Choi verification and the target is the standard cobit/CNOT isometry.

**Hard kill criterion.** Kill the one-step claim if (i) the decoder and target Bell/reset measurement are trusted; (ii) source replacement turns the two scores into ordinary GHZ stabilizers; (iii) the theorem is a direct Choi-fidelity, Hofmann-fidelity, or Stinespring-continuity bound; (iv) “late choice” means only choosing a complementary measurement setting on one static Choi state; or (v) hidden junk is excluded by assumption but the conclusion says that no record exists globally.

**Only defensible residual direction.** Move the random cut inside a genuinely multi-time comb. Fresh verifier inputs must enter at multiple times, the cut \(J\) must be sampled after a common sequential prefix, the audit must query a causal fact not reducible to the initial label, and the alternative branch must reverse that same comb. A publishable result would require a nonfactorizing robust bound against direct-label, precompiled-predicate, separate-device, and hidden-transcript strategies. If that sequential game reduces under source replacement to verification of one static graph/Choi state or to independent per-step Hofmann tests, it is killed as well.

### 4.10 Coarse-grained QND audit versus entanglement return: only exact payoff-specific frontiers remain open

#### Precise restricted model

The strongest static candidate found during this audit is no longer a claim about a history-bearing observer. It is a finite-dimensional information--coherence optimization. Let \(H\) be uniform on \(\{0,1\}^{n}\), let \(d=2^n\), and restrict the device to a canonical diagonal/QND instrument

\[
  M_k=\sum_{h\in\{0,1\}^n}\sqrt{q(k\mid h)}\,|h\rangle\!\langle h|,
  \qquad \sum_k q(k\mid h)=1.
\]

The only accessible record is the classical flag \(K=k\); the public \(d\)-level system remains coherent. A late uniformly random coordinate \(X\) defines the audit score

\[
 P_{\mathrm{coord}}
 =\frac{1}{nd}\sum_{k,x}
   \max_{b\in\{0,1\}}
   \sum_{h:h_x=b}q(k\mid h).
\]

For the alternative challenge, a flag-conditioned decoder tries to return half of a maximally entangled state. In the real, nonnegative, no-hidden-junk model above, its entanglement fidelity is

\[
 F_{\mathrm{return}}
 =\frac{1}{d^2}\sum_k
   \left(\sum_h\sqrt{q(k\mid h)}\right)^2.
\]

This formula must not be silently extended to arbitrary instruments. Extra outcome-dependent phases can require a correcting diagonal unitary; extra Kraus multiplicity or inaccessible junk can lower the attainable return; and allowing a coherent rather than classical \(K\) changes the access model.

#### The exact architecture collides with multipath duality games

Bagan, Calsamiglia, Bergou, and Hillery already formulate a common \(N\)-path state followed by a random late **Ways** or **Phases** challenge and derive a tight arbitrary-\(N\) discrimination region [PRL 120, 050402 (2018)](https://doi.org/10.1103/PhysRevLett.120.050402), [arXiv:1708.03968](https://arxiv.org/abs/1708.03968). Their phase success is an operational form of the path state's \(\ell_1\) coherence. For the canonical diagonal instrument above, let

\[
 \rho_{hh'}=\frac{1}{d}\sum_k
 \sqrt{q(k\mid h)q(k\mid h')}.
\]

All entries are nonnegative, and therefore

\[
 F_{\mathrm{return}}
 =\frac1d+\frac{d-1}{d}\,C_{\ell_1}^{\mathrm{norm}}(\rho),
 \qquad
 C_{\ell_1}^{\mathrm{norm}}(\rho)
 =\frac{1}{d-1}\sum_{h\ne h'}|\rho_{hh'}|.
\]

Thus the proposed echo score is not a new wave quantity: in this restricted model it is affine in the same normalized coherence used in multipath duality. The difference is solely on the audit axis. Bagan *et al.* ask for the entire path label; the proposed score asks for one randomly selected Boolean coordinate. Hillery subsequently studies partial Ways/Phases games, including set-valued answers and fixed partitions, and proves a mutual-information duality relation [J. Phys. A 54, 495301 (2021)](https://doi.org/10.1088/1751-8121/ac367d), [arXiv:2106.01514](https://arxiv.org/abs/2106.01514). That is a direct collision with the phrase "partial path information," although it is not the same Bayes random-coordinate functional.

The general finite-group extension of Bagan *et al.* trades group asymmetry against discrimination of representation subspaces [J. Phys. A 51, 414015 (2018)](https://doi.org/10.1088/1751-8121/aabb21), [arXiv:1803.04079](https://arxiv.org/abs/1803.04079). It reinforces that symmetry and group-labelled alternatives are occupied ingredients. It does not state the Boolean-coordinate frontier below.

#### Nearest exact trade-off and recovery results

| Primary source | Exact result already occupied | Why it is not yet the same theorem |
|---|---|---|
| [Banaszek, PRL 86, 1366 (2001)](https://doi.org/10.1103/PhysRevLett.86.1366) | Tight arbitrary-dimensional estimation-fidelity versus operation-fidelity frontier for a Haar-uniform unknown pure state | The information score is universal state estimation, not a random commuting coordinate of a known orthogonal label. |
| [Barnum, arXiv:quant-ph/0205155 (2002)](https://arxiv.org/abs/quant-ph/0205155) | General measurement model, convex information--disturbance frontier, square-root least-disturbing dynamics, and covariance reduction for the uniform ensemble | It supplies much of the optimization methodology, but not this Boolean payoff or Hamming-shell spectrum. This item remained a preprint in the corpus located here. |
| [Mišta and Filip, PRA 72, 034307 (2005)](https://doi.org/10.1103/PhysRevA.72.034307) | A QND measurement realizes the optimal Banaszek fidelity trade-off | QND implementation is occupied; the utility and ensemble differ. |
| [Cheong and Lee, PRL 109, 150402 (2012)](https://doi.org/10.1103/PhysRevLett.109.150402) | Tight arbitrary-dimensional estimation-gain versus probabilistic exact-reversal bound | Its reversal score is success probability governed by the smallest singular values, not deterministic flag-assisted entanglement fidelity. |
| [Lee, Kim, and Nha, Quantum 5, 414 (2021)](https://doi.org/10.22331/q-2021-03-17-414) | Complete trade-offs among information gain, operation disturbance, and probabilistic reversibility | Same general sector, different information ensemble and recovery functional. |
| [Berta, Coles, and Wehner, PRA 90, 062127 (2014)](https://doi.org/10.1103/PhysRevA.90.062127) | An exact equality between guessing complementary measurement outcomes and recoverable entanglement fidelity | This is the most important score-level warning. It concerns noncommuting MUB choices and quantum side information, not a classical flag guessing random commuting Boolean coordinates while a separate coherent output is recovered. No exact variable substitution was located. |
| [Gregoratti and Werner, J. Mod. Opt. 50, 915 (2003)](https://doi.org/10.1080/0950034021000058021) and [Memarzadeh, Macchiavello, and Mancini, NJP 13, 103031 (2011)](https://doi.org/10.1088/1367-2630/13/10/103031) | Environment measurement, flag-conditioned correction, and entanglement-fidelity optimization | They occupy environment-assisted recovery, but do not jointly optimize a random-coordinate audit payoff. |
| [Liu *et al.*, PRA 111, 062215 (2025)](https://doi.org/10.1103/PhysRevA.111.062215) | A resource-theoretic bound trading classical distinguishability extracted from orthogonal pure states against coherence retained after discrimination | It is a close general information/coherence collision, but its resources and perfect-discrimination setting are not an arbitrary orbital Bayes payoff versus flag-conditioned entanglement return. |
| [Mohan, Tavakoli, and Brunner, NJP 21, 083034 (2019)](https://doi.org/10.1088/1367-2630/ab3773) and [Wei *et al.*, arXiv:2103.03075 (2021)](https://arxiv.org/abs/2103.03075) | Tight sequential \(2\to1\) and \(3\to1\) QRAC score pairs with instrument characterization/self-testing | These optimize two sequential decoding scores for qubit communication, not the diagonal classical-record score versus reference-entangled channel return. The \(3\to1\) item was located as a preprint, not an archival journal article. |
| [Asadian, Gams, and Sponar, PR Research 8, L012011 (2026)](https://doi.org/10.1103/llgb-gql9) | A covariant correlation--disturbance relation for sequential \(n\)-outcome measurements in arbitrary dimension | Peer reviewed and highly adjacent, but its correlation and disturbance measures are not \(P_{\mathrm{coord}}\) and \(F_{\mathrm{return}}\). |
| [Zhang *et al.*, accepted PRA, 14 July 2026](https://doi.org/10.1103/yxq6-n2nr), [arXiv:2608.06726](https://arxiv.org/abs/2608.06726) | Fully device-independent unsharp-instrument characterization through an entanglement-assisted sequential QRAC | Current and directly relevant, but still a sequential unsharp-measurement score rather than the arbitrary-\(n\) diagonal Boolean frontier. It was an accepted manuscript, not yet a version-of-record article at the cutoff, so details remain provisional. |
| [Krajenbrink *et al.*, arXiv:2606.04843 (2026)](https://arxiv.org/abs/2606.04843) | Extends decoded quantum interferometry from Hamming space to translation association schemes and analyzes shell amplitudes through a finite spectral problem | This very recent preprint occupies the association-scheme/shell-eigenproblem methodology, but not a measurement audit, environment flag, or entanglement-recovery frontier. |

#### Hamming reduction: promising calculation, not yet a global theorem

For the XOR-covariant subclass \(q(k\mid h)=r(k\oplus h)\), relabel the outcome as the coordinatewise guess and let

\[
 p_m=\sum_{e:|e|=m}r(e),\qquad
 z_m=\sqrt{p_m},\qquad
 u_m=\sqrt{\binom nm}.
\]

After averaging within each Hamming shell,

\[
 P_{\mathrm{coord}}=\sum_{m=0}^{n}\left(1-\frac mn\right)z_m^2,
 \qquad
 F_{\mathrm{return}}=\frac1{2^n}
 \left(\sum_{m=0}^{n}u_m z_m\right)^2,
 \qquad \sum_m z_m^2=1.
\]

Consequently the covariant-subclass support function is

\[
 \max_{\|z\|=1}
 \left[\lambda P_{\mathrm{coord}}+(1-\lambda)F_{\mathrm{return}}\right]
 =\lambda_{\max}\!\left[
 \lambda\,\operatorname{diag}\!\left(1,1-\frac1n,\ldots,0\right)
 +\frac{1-\lambda}{2^n}uu^{\mathsf T}
 \right].
\]

For \(n=2\) and \(\lambda=1/2\), this is the \(3\times3\) diagonal-plus-rank-one matrix

\[
 \frac12\operatorname{diag}(1,1/2,0)
 +\frac18(1,\sqrt2,1)^{\mathsf T}(1,\sqrt2,1),
\]

whose largest eigenvalue agrees with the independently obtained numerical value \(0.8117449009\). This agreement validates the covariant calculation, not the missing global reduction.

A publishable theorem would have to optimize over **all** finite-output channels \(q(k\mid h)\), not assume the answer. In particular it must prove that outcome refinement and group averaging can preserve the Bayes coordinate score while not reducing the recovery score; that every flag can be labelled by an optimal Hamming centre; that stabilizer/permutation averaging legitimately reduces it to Hamming shells; and that no noncovariant instrument, hidden Kraus refinement, or more general diagonal quantum instrument beats the matrix value. The max in \(P_{\mathrm{coord}}\) and the square roots in \(F_{\mathrm{return}}\) make an unqualified "by symmetry" step inadequate.

#### Association-scheme strengthening

There is a plausible broader mathematical statement. Let a finite group act transitively on a history space \(\Omega\), let an audit guess \(g\in\Omega\), and suppose the payoff \(w(g,h)=w_r\) depends only on the orbital/relation \(r\) containing \((g,h)\). If an optimal covariant reduction is proved for the same canonical diagonal instrument class, and if \(v_r\) is the valency of relation \(r\), the orbit-amplitude vector would give the candidate support matrix

\[
 \lambda\,\operatorname{diag}(w_r)
 +\frac{1-\lambda}{|\Omega|}uu^{\mathsf T},
 \qquad u_r=\sqrt{v_r}.
\]

The Boolean Hamming scheme, \(q\)-ary Hamming schemes, and Johnson schemes would then be corollaries. The orbit/valency and Bose--Mesner machinery itself is classical [Bose and Mesner, Ann. Math. Statist. 30, 21--38 (1959)](https://doi.org/10.1214/aoms/1177706356); using it is not a novelty claim. Moreover, a general coherent configuration need not be commutative or multiplicity free, so a universal scalar diagonal-plus-rank-one reduction may fail and require matrix blocks.

The searches through 12 August 2026 did not locate a primary source proving this exact association-scheme information-versus-entanglement-return support theorem. Searches combined "association scheme," "coherent configuration," "Hamming scheme," "group-covariant instrument," "information disturbance," "environment-assisted recovery," "entanglement fidelity," "random access code," and "wave-particle duality." The closest group and partial-information papers are those above. This is evidence of an unlocated conjunction, not proof of priority.

#### \(q\)-ary large-\(n\) spectral transition

The \(q\)-ary Hamming family supplies a sharper theorem candidate. Its shell valencies and normalized weights are

\[
 v_w=\binom nw(q-1)^w,\qquad
 \mu_n(w)=\frac{v_w}{q^n}
 =\operatorname{Bin}\!\left(n,\frac{q-1}{q}\right)(w).
\]

For the covariant support matrix, any top eigenvalue \(s\) above the diagonal edge \(\lambda\) obeys the rank-one secular equation

\[
 1=(1-\lambda)\sum_{w=0}^{n}
 \frac{\mu_n(w)}
 {s-\lambda(1-w/n)}.
\]

Concentration of \(\mu_n\) at \(w/n=(q-1)/q\) suggests, and a proof should make uniform away from the critical point, the limiting support

\[
 \lim_{n\to\infty}s_n(\lambda)
 =\max\!\left\{
 \lambda,\,
 1-\lambda\!\left(1-\frac1q\right)
 \right\},
 \qquad
 \lambda_c=\frac{q}{2q-1}.
\]

Below \(\lambda_c\), the candidate optimizer remains distributed over the typical Hamming shell, with \(F_{\mathrm{return}}\to1\) and \(P_{\mathrm{coord}}\to1/q\). Above \(\lambda_c\), it localizes at the zero-error shell, \(P_{\mathrm{coord}}\to1\), and the return fidelity vanishes asymptotically. The behavior exactly at \(\lambda_c\), finite-size scaling, convergence rates, and optimizer localization still require proof; they should not be inferred from the pointwise secular equation alone.

No primary source located in the targeted searches states this Hamming/QND/RAC spectral transition or the threshold \(q/(2q-1)\). Searches included "phase transition wave-particle duality," "Hamming random access disturbance," "weak measurement many bits," "multipath large-\(N\) critical," and rank-one perturbations. However, the qualitative asymptotic geometry is strongly Bagan-adjacent: Bagan *et al.* already show that their tight normalized full-path/coherence region approaches the triangle \(x+y\leq1\) as the number of paths tends to infinity [PRL 120, 050402 (2018)](https://doi.org/10.1103/PhysRevLett.120.050402). In raw coordinate-success units, the two endpoints \((P,F)=(1/q,1)\) and \((1,0)\) alone produce the same maximum of two affine functions and the same \(q\)-dependent crossing after rescaling.

Accordingly, the limiting kink by itself is **amber-red**, not a strong standalone novelty claim. It becomes potentially substantive only if the work proves the full finite-\(n\) frontier first and then derives nontrivial finite-size scaling, an eigenvector-localization law, critical exponents/window, or a universality theorem across association schemes that is not merely the support function of a limiting line segment.

#### Rooted-tree random-prefix variant: new payoff candidate, old geometry

A temporal-looking alternative labels the \(q^n\) leaves by strings \(h\in[q]^n\). The safest static game requires the apparatus to commit a complete candidate transcript \(g\in[q]^n\) **before** a uniformly random cut \(T\in\{1,\ldots,n\}\) is revealed. It wins exactly when \(g_{1:T}=h_{1:T}\). If

\[
 \ell(g,h)=\max\{t:g_{1:t}=h_{1:t}\}
\]

is the longest-common-prefix length, averaging over the late cut gives the orbital payoff

\[
 w(g,h)=\frac{\ell(g,h)}{n}.
\]

Around any fixed leaf, the shell valencies are

\[
 v_\ell=(q-1)q^{\,n-\ell-1}\quad(0\leq\ell<n),
 \qquad v_n=1.
\]

For the tree-automorphism-covariant subclass of the canonical QND model, write \(p_\ell\) for the total conditional mass in shell \(\ell\), \(z_\ell=\sqrt{p_\ell}\), and \(u_\ell=\sqrt{v_\ell}\). The candidate scores then reduce to

\[
 P_{\mathrm{prefix}}=\sum_{\ell=0}^{n}\frac{\ell}{n}z_\ell^2,
 \qquad
 F_{\mathrm{return}}=\frac1{q^n}
 \left(\sum_{\ell=0}^{n}u_\ell z_\ell\right)^2,
\]

and the covariant support function is

\[
 \lambda_{\max}\!\left[
 \lambda\,\operatorname{diag}\!\left(0,\frac1n,\ldots,\frac{n-1}{n},1\right)
 \;+\;\frac{1-\lambda}{q^n}uu^{\mathsf T}
 \right].
\]

This calculation is another instance of the association-scheme ansatz above, not yet a proof of global optimality over all declared instruments.

The normalized shell law is geometric rather than binomial:

\[
 \mu_n(\ell)=\frac{v_\ell}{q^n}
 =
 \begin{cases}
  (q-1)q^{-\ell-1},&0\leq\ell<n,\\
  q^{-n},&\ell=n.
 \end{cases}
\]

Hence a covariant top eigenvalue \(s>\lambda\) obeys

\[
 1=(1-\lambda)\sum_{\ell=0}^{n}
 \frac{\mu_n(\ell)}{s-\lambda\ell/n}.
\]

For fixed \(\lambda\) away from \(1/2\), this suggests the candidate asymptotics

\[
 \lim_{n\to\infty}s_n(\lambda)=\max\{1-\lambda,\lambda\},
 \qquad \lambda_c=\frac12.
\]

On the distributed side \(\lambda<1/2\), a secular-equation expansion gives

\[
 s_n(\lambda)
 =1-\lambda+\frac{\lambda}{(q-1)n}+O(n^{-2}),
\]

while for fixed \(\lambda>1/2\), localization at the exact-transcript shell suggests

\[
 s_n(\lambda)
 =\lambda+
 \frac{\lambda(1-\lambda)}{2\lambda-1}\,q^{-n}
 +o(q^{-n}).
\]

These expansions are candidate consequences of the covariant matrix, not established global instrument bounds; uniform control in the \(O(1/n)\) critical window remains to be proved. The limiting kink alone is again just the support function of the two endpoints \((P,F)\to(0,1)\) and \((1,0)\), so it is weaker originality evidence than the exact finite-\(n\) curve.

The commitment clause is essential. If the device sees \(T\) and may then choose a separate prefix \(a_T\), the natural score is instead

\[
 P_{\mathrm{free\ cut}}
 =\frac1{nq^n}\sum_{k,T}
 \max_{a\in[q]^T}
 \sum_{h:h_{1:T}=a}q(k\mid h).
\]

The maximizing prefixes for different \(T\) need not be mutually consistent, so they need not be prefixes of any single \(g\). Therefore \(P_{\mathrm{free\ cut}}\geq P_{\mathrm{prefix}}\), and the LCP spectral formula does **not** establish the frontier for the freer decoder. Any theorem or experiment must say which game is meant.

The tree geometry itself is occupied mathematics. Up to reversing coordinate order, the LCP relation is the one-chain kernel/Niederreiter--Rosenbloom--Tsfasman relation used to build ordered Hamming schemes [Martin and Stinson, Canadian J. Math. 51, 326--346 (1999)](https://doi.org/10.4153/CJM-1999-017-5). Equivalently, its \(n+1\) relations arise from an iterated wreath product of one-class schemes; the relevant Bose--Mesner and Terwilliger structures are established [Bhattacharyya, Song, and Tanaka, J. Algebraic Combin. 31, 455--466 (2010)](https://doi.org/10.1007/s10801-009-0196-x), [Song and Xu, arXiv:1008.2228 (2010)](https://arxiv.org/abs/1008.2228). Ordered-Hamming schemes have also appeared in quantum-walk and spin-network models [Miki, Tsujimoto, and Vinet, arXiv:1712.09200](https://arxiv.org/abs/1712.09200). None of those mathematical or quantum-walk sources supplies the proposed information--return frontier.

The closest physical papers remain partial-path complementarity rather than this exact nested-prefix score. Hillery permits set-valued partial answers in late Ways/Phases games [J. Phys. A 54, 495301 (2021)](https://doi.org/10.1088/1751-8121/ac367d). Banerjee *et al.* derive a wave--particle relation when only an incomplete set of path pointers is accessible [arXiv:2108.05849 (2021)](https://arxiv.org/abs/2108.05849); this item remained a preprint in the corpus located here. Wu and Wang discuss complementarity for hierarchically arranged bipartite systems, including a binary-tree representation [Entropy 22, 813 (2020)](https://doi.org/10.3390/e22080813), but not LCP discrimination, a late random prefix, or entanglement-return optimization.

Targeted searches through 12 August 2026 did not locate a primary source giving the same committed-transcript random-prefix payoff together with optimal flag-conditioned entanglement return. This is a report about the searched corpus, not an absence or priority theorem. The nearest sources already occupy partial path information, hierarchical complementarity, the underlying tree association scheme, and the generic information--disturbance architecture.

Most importantly, a static orthogonal input \(|h\rangle\) does not certify that \(h_1,\ldots,h_n\) occurred as a causal history. The rooted-tree score becomes temporally meaningful only if fresh verifier-controlled symbols enter a genuine sequential comb, the device commits before the random cut, and soundness excludes a direct transcript port or precompiled leaf label. Even then, simply storing the entire transcript may be an honest optimal strategy unless a memory, communication, or control constraint makes temporal structure operationally relevant.

#### Originality verdict and kill conditions

- **Red:** the late Ways/Phases or audit/coherence architecture; the \(n=1\) two-path case; full-label path discrimination; QND realization; entanglement fidelity as the recovery score; and symmetry reduction as a generic method.
- **Amber:** the exact all-\(n\) Boolean random-coordinate frontier over the declared canonical diagonal instrument class. No exact primary-source collision was located, but it is one coarse-grained-payoff specialization away from multipath duality and may have a short proof from known ingredients.
- **Amber, potentially stronger:** a rigorously stated transitive association-scheme theorem with several nontrivial families, equality cases, asymptotics, and a proof that covers all allowed instruments. Its physics interpretation should be "coarse-grained multipath information versus coherent return," not a new reversible-history principle.
- **Amber-red in isolation:** the \(q\)-ary large-\(n\) kink at \(\lambda_c=q/(2q-1)\). Its exact Hamming spectral realization was not located, but its limiting two-endpoint geometry is an affine reparametrization of the large-\(N\) triangular duality region already identified by Bagan *et al.*
- **Amber:** the exact finite-\(n\) committed-prefix frontier. The LCP payoff and rank-one matrix were not located in the searched quantum literature, but the tree scheme is standard and the calculation is a direct specialization of the proposed association-scheme reduction. The unrestricted cut-dependent decoder remains a separate unsolved optimization.

Kill the candidate if any of the following occurs:

1. the claimed frontier is an exact substitution into Bagan *et al.*, Hillery, Berta--Coles--Wehner, or a general cost-sensitive discrimination theorem;
2. a noncovariant \(q(k\mid h)\), an outcome refinement, or a permitted diagonal instrument exceeds the Hamming eigenvalue;
3. the proof establishes only the XOR-covariant subclass or only \(n=2\) numerics;
4. the "optimal recovery" formula ignores inaccessible Kraus labels, phases, residual environments, or the declared decoder access;
5. the association-scheme statement merely restates orbit averaging and Rayleigh--Ritz without a new global optimization or a useful tight inequality;
6. the result is presented as certification of a temporal history even though the entire experiment is one static QND interaction with a pre-existing orthogonal label;
7. the large-\(n\) result consists only of the maximum of the two endpoint scores, without a new finite-size or localization theorem;
8. the random-prefix theorem permits cut-dependent, mutually inconsistent answers but proves only the committed-full-transcript LCP formula; or
9. the rooted-tree result is only the association-scheme corollary, without a global instrument proof, a distinct inequality, or genuine sequential soundness.

**Gate recommendation.** This narrow theorem is worth proving before investing in a large experiment. Treat it as a standalone combinatorial quantum-information result. If global covariance fails or the result collapses to an existing multipath duality theorem, kill the static candidate and retain only the genuinely sequential random-cut comb direction.

#### Formal priority verdict for the structured-payoff candidates

No source located in the searches through 12 August 2026 states the following theorem at the advertised level of generality: an arbitrary finite transitive history space \(\Omega\), an arbitrary invariant/orbital audit payoff, a classical flag produced by a diagonal QND instrument, optimal flag-conditioned entanglement recovery, and a global reduction of the joint support function to the valency matrix

\[
 \lambda\operatorname{diag}(w_r)
 +(1-\lambda)|\Omega|^{-1}uu^{\mathsf T}.
\]

This must be reported as **not located**, not as nonexistent. The constituent ingredients are close: Barnum supplies covariant information--disturbance optimization for a different ensemble and scores; Bagan *et al.* supply a tight finite-group Ways/Phases architecture; Berta--Coles--Wehner identify a guessing probability with recoverable entanglement fidelity in a different complementary-observable setting; Gregoratti--Werner and Memarzadeh *et al.* supply environment-flag-conditioned recovery; Liu *et al.* trade extracted distinguishability against retained coherence; and Krajenbrink *et al.* use association-scheme shell spectral reductions in a different interferometric algorithm. A new theorem must do more than place those tools side by side.

| Candidate claim | Formal verdict at the cutoff | Subjective chance that expert review finds an exact or essentially equivalent prior theorem |
|---|---|---:|
| Arbitrary transitive/orbital payoff versus flag-conditioned entanglement return | **Amber.** No exact source located; components and likely proof methods are occupied. | 25--40% |
| Exact finite-\(n\) Boolean/q-ary Hamming random-coordinate frontier over all declared instruments | **Amber.** Covariant diagonal-plus-rank-one curve is clear; the global reduction is missing. | 35--55% |
| Hamming large-\(n\) support \(\max\{\lambda,1-\lambda(1-1/q)\}\) and crossing \(q/(2q-1)\), claimed alone | **Amber-red.** The precise Hamming derivation was not located, but its limiting geometry is Bagan-adjacent and may be an implicit rescaling. | 55--75% |
| Exact finite-\(n\) rooted-tree committed-prefix frontier over all declared instruments | **Amber.** The payoff-specific quantum curve was not located; its LCP/wreath-product geometry is standard. | 20--35% |
| Rooted-tree large-\(n\) support \(\max\{\lambda,1-\lambda\}\) and crossing \(1/2\), claimed alone | **Amber-red.** It is the generic two-endpoint limiting triangle; only finite-size or critical-window structure could carry novelty. | 65--85% |
| Genuine sequential-comb version with fresh verifier inputs, post-prefix transcript commitment, and one nonfactorizing audit/return soundness theorem | **Amber-green as a research question, not as a priority claim.** No exact conjunction was located, but cryptographic test/Hadamard games and temporal-memory witnesses are close. | 20--35% |

These percentages are planning priors, not bibliometric measurements. “Essentially equivalent” includes a general theorem whose specialization yields the proposed formula even if it uses different terminology.

**Safe novelty sentence now:**

> Within a declared diagonal-QND instrument model, we formulate payoff-specific audit--return frontiers for Hamming-coordinate and committed-prefix rooted-tree queries; the spectral reductions are established only for covariant strategies, while global optimality and genuine temporal soundness remain open.

If the global optimization is later proved, “formulate” may become “prove an exact finite-\(n\) trade-off,” but “first audit-or-coherence game,” “first hierarchical complementarity relation,” and “first association-scheme quantum reduction” would still be unsafe.

### 4.11 Online, classical-memory-only version: the resource class is occupied; the joint frontier was not located

The online restriction is materially different from the static diagonal-QND model in Section 4.10, but it does **not** create a new resource class. Quantum combs with bounded coherent memory and free or bounded classical memory, multi-time processes whose memory carrier transmits only classical information, classically adaptive testers, and sequential RAC witnesses are all established subjects.

#### Access model that must be fixed before any theorem

A clean zero-coherent-memory interface is the following. At round \(t\), the verifier supplies a fresh public system \(S_t\), possibly entangled with a retained reference \(R_t\). Conditional on a classical state \(c_{t-1}\), the device applies an instrument

\[
 \left\{\mathcal J^{(t)}_{z_t\mid c_{t-1}}:
 \mathsf L(S_t)\rightarrow\mathsf L(S'_t)\right\}_{z_t},
\]

returns \(S'_t\) to the verifier, and updates \(c_t\) from \((c_{t-1},z_t)\). No quantum system passes inside the device from round \(t\) to \(t+1\). The verifier sequesters every returned \(S'_t\) until a challenge is sampled only after round \(n\). The two challenges are then:

- **audit:** reveal a random coordinate, or a prefix fact, using only the committed classical state/transcript;
- **return:** apply a precisely charged decoder to all sequestered public outputs, optionally conditioned on the classical transcript, and score global entanglement fidelity together with decoupling/reset of the transcript.

Three distinctions are essential.

1. **No quantum memory inside the device is not no quantum memory in the experiment.** The sequestered public outputs can preserve coherence. The claim must identify the storage boundary rather than call the whole experiment classical-memory-only.
2. **The decoder is a resource.** If the trusted verifier performs an arbitrary joint decoder, that operation is outside the device's prefix memory budget. If the device performs it, the protocol must say when the device receives the outputs back and what post-challenge quantum workspace is allowed.
3. **Entanglement breaking is not by itself a safe synonym for classical multi-time memory.** Taranto *et al.* distinguish several structurally different classical-memory notions, and Vieira, Ku, and Budroni show that an entanglement-breaking step can still coexist with nonclassical multi-time correlations. The admissible conditional-instrument factorization, hidden environments, and side channels must therefore be stated directly.

For independent QND-labelled inputs and one fine-grained classical transcript \(\mathbf z\), the zero-coherent-memory restriction gives a causal product constraint. In the simplest single-Kraus diagonal realization,

\[
 m_{\mathbf z}(\mathbf h)
 =\prod_{t=1}^{n}
 m^{(t)}_{z_t\mid \mathbf z_{<t}}(h_t),
 \qquad
 q(\mathbf z\mid\mathbf h)=|m_{\mathbf z}(\mathbf h)|^2.
\]

Thus the online feasible set is not the arbitrary global kernel \(q(k\mid h)\) of the static calculation. It is a causally factored, adaptively chosen local instrument, with a classical hidden-state restriction if \(|c_t|\) is bounded. Multiple Kraus operators and unobserved outcomes enlarge the notation but do not remove the need for an explicit causal factorization. This is the most plausible source of a genuinely new theorem: a tight audit--return frontier for that causal subset, or a strict separation from collective static instruments and from combs with one or more coherent memory qubits.

#### Exact nearest primary sources

| Primary source | What is already occupied | What it does **not** establish for this candidate |
|---|---|---|
| [Kretschmann and Werner, *PRA* 72, 062323 (2005)](https://doi.org/10.1103/PhysRevA.72.062323) | General causal channels with memory, represented by concatenated memory channels. | No late audit/return score and no zero-coherent-memory history-query frontier. |
| [Bisio *et al.*, *PRA* 85, 032333 (2012)](https://doi.org/10.1103/PhysRevA.85.032333) | Quantum-strategy memory cost with assistance from classical memory; importantly, memory cost is a global property and cannot generally be optimized one time step at a time. | No random past-fact audit or entanglement-return alternative. This source blocks an unsupported claim that a per-round optimization automatically proves the global bound. |
| [Taranto *et al.*, *Quantum* 8, 1328 (2024)](https://doi.org/10.22331/q-2024-05-02-1328) | A hierarchy of multi-time processes with classical memory. Their classical-memory process has a sequence of conditional instruments, while a memoryless process has independent channels. | No joint information-readout versus reversal score for a verifier-sequestered stream. |
| [Roy *et al.*, *PRA* 110, 012608 (2024)](https://doi.org/10.1103/PhysRevA.110.012608) | Sequential QRAC scores as a semi-device-independent witness separating classical-memory/Markovian processes from quantum-memory processes; includes joint intermediate/later success regions and unsharp instruments. | The later score is another RAC decoder's success, not global entanglement recovery of all earlier public systems. |
| [Nakahira and Kato, *PRA* 103, 062606 (2021), with erratum](https://doi.org/10.1103/PhysRevA.103.062606) | Convex optimization, duality, and symmetry for general adaptive multi-time process-discrimination objectives and constraints. | Entanglement-return fidelity of the same committed prefix is not supplied as the candidate's second counterfactual branch. A claimed theorem must show why it is not a direct instance of this general optimization framework. |
| [Ohst *et al.*, *Quantum* 10, 1988 (28 January 2026)](https://doi.org/10.22331/q-2026-01-28-1988) | The sharpest formal collision. Adaptive discrimination with bounded quantum memory and classical memory is expressed through constrained separability; the intermediate instrument outcomes encode classical memory. The paper explicitly studies classically adaptive protocols without coherent transfer between channel uses and finds non-hierarchical relations with parallel strategies. | Its objective is channel discrimination, not a late alternative between a temporal audit and returned entanglement/reset. The paper's constrained-separability machinery is nevertheless an immediate baseline for any proposed SDP. |
| [Zonnios and Binder, arXiv:2606.19511 (17 June 2026)](https://arxiv.org/abs/2606.19511) | **Provisional 2026 preprint.** Machines for autonomous distinction repeatedly apply an instrument, retain the full classical outcome record, and carry coherent memory of dimension \(d_A\); setting \(d_A=1\) gives the zero-coherent-memory endpoint by inference. | It studies recurrent process discrimination and a coherent-memory hierarchy, not audit versus recovery. Because it is a recent preprint, use it as a priority warning rather than as definitive authority for an unreviewed theorem. |
| [Mohan, Tavakoli, and Brunner, *NJP* 21, 083034 (2019)](https://doi.org/10.1088/1367-2630/ab3773) | Tight sequential-QRAC information--disturbance curves and semi-device-independent self-testing of instruments with classical and quantum outputs. | One system is passed through sequential decoders; there is no multi-round sequestered-output return test. |
| [Chiribella and Goswami, *PRX Quantum* 6, 020335 (2025)](https://doi.org/10.1103/PRXQuantum.6.020335) | Maximum and minimum quantum causal effects connect classical signalling capacity and approximate quantum invertibility, with a duality/monogamy relation between recoverability toward one output and causal influence toward another. | The information score is a channel-level causal-effect measure, not a random stored-coordinate/prefix audit under a zero-coherent-memory controller. Any new bound must nevertheless be compared to this causal information--recovery duality. |
| [Lim, Hhan, and Kwon, *QST* 10, 025048 (2025)](https://doi.org/10.1088/2058-9565/adc034) | Information-gain/non-disturbance bounds for local discrimination of entangled states, plus adaptive nondestructive strategies assisted by preshared entanglement. | The ensemble, locality constraint, and nondestructive score differ from a random past-coordinate audit followed counterfactually by global entanglement return. |
| [Wiesner and Crutchfield, *Physica D* 237, 1173--1195 (2008)](https://doi.org/10.1016/j.physd.2008.01.021) | Stochastic and quantum finite-state transducers for repeatedly measured input--output processes and their process-language hierarchy. | No late counterfactual audit/return test or entanglement-fidelity frontier. The finite-state vocabulary itself is occupied. |
| [Gregoratti and Werner, *J. Mod. Opt.* 50, 915--933 (2003)](https://doi.org/10.1080/0950034021000058021) | Classical environment outcomes can condition quantum correction. | No causal online audit score; it supplies the correct flag-conditioned recovery baseline. |

These sources occupy almost every ingredient. The literature search through 12 August 2026 did **not** locate one primary paper that combines all of the following clauses: fresh sequential inputs; public outputs sequestered after each round; no coherent memory retained by the device; one committed classical state/transcript; a post-prefix random coordinate or prefix audit; an alternative joint entanglement-return/reset test; and a tight adversarial frontier for the same prefix instrument. That statement is a search result, not proof of absence.

#### Collision and residual-novelty verdict

| Claim | Verdict |
|---|---|
| "We introduce online quantum instruments with classical memory but no quantum memory." | **Red.** Bisio, Taranto *et al.*, Ohst *et al.*, and the 2026 MAD preprint already provide overlapping formalisms. |
| "Sequential random-access success certifies quantum rather than classical temporal memory." | **Red.** Roy *et al.* directly do this, building on sequential QRAC work. |
| "Local/adaptive information extraction trades against nondisturbance." | **Red.** This is a mature information--disturbance and nondestructive-discrimination theme. |
| "A late coin commits one prefix before either score is selected." | **Red as novelty by itself.** The coin is operationally useful against branch-specific devices, but it does not create a new feasible region if both scores are merely two functionals of an already specified instrument. |
| Exact optimum for a random past-coordinate audit versus transcript-conditioned global entanglement return over the causally factored zero-coherent-memory class | **Amber.** No exact collision located; close general optimization machinery exists. |
| A strict, quantitative separation among (i) collective static diagonal instruments, (ii) online classical-memory instruments, and (iii) online instruments with \(d_A>1\) coherent memory, under the same audit--return game | **Amber-green as a research target.** It would connect comb memory cost to information--recovery frontiers in a way not located here. |
| A finite-state law relating classical-state size, history length, audit score, and return fidelity | **Amber.** Defensible only if it is not an ordinary classical finite-state/RAC bound pasted onto an independent channel-recovery inequality. |

**Safe residual novelty sentence:**

> We study a late two-score game over classically adaptive, zero-coherent-memory quantum instruments, with public outputs sequestered between rounds; the proposed contribution is a tight joint audit--return frontier, or a separation from collective and coherent-memory strategies, rather than the instrument class or delayed choice itself.

**Hard kill criteria for the online version.** Kill or sharply reframe it if:

1. the verifier supplies a directly readable classical history label, making the audit a trivial classical transcript task;
2. the public systems are not sequestered, so the device can retain or revisit an undeclared quantum side channel;
3. the trusted joint decoder's quantum memory and computation are attributed to the prefix device;
4. the return score omits the persistent classical transcript or other side information from its decoupling condition;
5. "classical memory" is imposed only by an entanglement-breaking link, without the stronger declared conditional-instrument/hidden-state model;
6. the optimum is obtained by multiplying or intersecting known one-round information--disturbance bounds, contrary to the global-memory warning of Bisio *et al.*;
7. the causal factorization gives the same optimum as the unrestricted static kernel for every stated payoff, with no separation, finite-state law, or new equality condition;
8. the proposed SDP is simply Ohst *et al.*'s constrained-separability program with a relabelled linear discrimination payoff and an independent fidelity calculation; or
9. the audit and return quantities specialize directly to Chiribella--Goswami maximum/minimum causal effects and their existing duality; or
10. the freer random-cut decoder is allowed to give mutually inconsistent prefixes but the proof covers only a committed full transcript.

**Priority estimate.** Conditional on this precise access model, the subjective chance that expert review finds an already published *exact* audit--return frontier is roughly 20--40%. The chance that the first proposed proof reduces to known comb optimization, sequential QRAC, or one-round information--disturbance machinery is higher, roughly 50--70%. These are planning priors, not bibliometric measurements. The online version should therefore begin with a two- or three-round separation theorem, not with a broad foundational claim.

## 5. Fundamental constraints that any new theorem must respect

### 5.1 Information versus recoverability is already controlled by channel theory

Kretschmann, Schlingemann, and Werner relate information leaked to an environment to correctability of the main channel through continuity of Stinespring dilations [IEEE Transactions on Information Theory 54, 1708–1717 (2008)](https://doi.org/10.1109/TIT.2008.917696). Bény and Oreshkov give necessary and sufficient conditions for approximate quantum error correction and express optimal entanglement recovery through a dual optimization on the environment [PRL 104, 120501 (2010)](https://doi.org/10.1103/PhysRevLett.104.120501).

Therefore, a proposed record–recovery bound is not new merely because its variables are named “history,” “memory,” and “recoherence.” It must exploit the multi-time challenge structure, restricted instruments, causal-memory dimension, computational restrictions, or a new operational quantity that is not an immediate corollary of complementary-channel recovery.

### 5.2 One run cannot deliver both destructive evidence and perfect recovery

The protocol can randomly test two counterfactual capabilities across repeated trials. It cannot both read an orthogonal, externally retained history transcript and obtain perfect interference between the corresponding alternatives in the same trial. The claim must be “the committed device succeeds under either challenge with the stated average score,” not “one execution simultaneously reveals the history and erases every record.”

### 5.3 Global absence of records is not device-independently observable without assumptions

An experiment observes selected ports. A dishonest or uncontrolled device can leak a transcript into an unobserved degree of freedom. A decoupling conclusion therefore requires one of:

- a trusted finite-dimensional dilation and a complete system boundary;
- a cryptographic residual-state security definition;
- a self-testing theorem with explicit side-information assumptions;
- an energy, dimension, communication, or isolation bound that rules out hidden storage.

Resetting the named work qubits to zero is not a certificate that no external record remains.

### 5.4 A fact query does not by itself certify a causal history

If the circuit designer knows the full history string, a device can directly prepare its answer or copy a branch label. To certify temporal memory, the verifier must introduce unpredictable interventions during the prefix and later ask a randomly chosen function of them. Even then, the theorem certifies a resource relative to an access model; it does not establish a unique internal narrative.

### 5.5 Postselection must be priced

Quantum erasers and measurement-reversal experiments often recover interference or entanglement conditionally. Any score must include the success probability. A high conditional fidelity with vanishing acceptance probability is not a useful recovery certificate and can make a classical or measure-and-prepare cheat look impressive.

## 6. Residual flagship candidate

### Post-commitment temporal random-access-or-entanglement-return game

This is the strongest surviving formulation.

1. The verifier holds a reference \(R\) entangled with a token \(Q\) supplied to the device.
2. Across \(T\) ordered interactions, the verifier supplies fresh random inputs \(X_t\). The device must update an internal memory and use it causally in later operations.
3. After the final interaction, a fresh challenge is sampled.
4. In an **audit** round, a fresh index or predicate \(J\) is requested. The score is a temporal random-access or process-witness score that lower-bounds memory dimension or excludes a declared classical process class.
5. In a **return** round, the device must undo its history-bearing dynamics and return \(Q\). The verifier tests entanglement fidelity of \(RQ\); a security definition bounds history information in the device’s residual state conditioned on allowed outputs.
6. The security theorem bounds the weighted success of any device outside the target resource class, including devices with challenge foreknowledge, direct history-label ports, or hidden external storage.

What could be genuinely new is not either branch. It is a **single-comb completeness and soundness theorem** connecting a temporal random-access lower bound to an entanglement-return/decoupling guarantee after an unpredictable challenge.

A useful target theorem would have the following form, with all norms and device classes made explicit:

> If one committed \(T\)-step comb attains audit score \(s_A\) and return entanglement fidelity \(F_E\), then (i) every compatible realization requires at least \(d(s_A)\) effective causal-memory dimension, while (ii) accepted return rounds leave at most \(\delta(F_E)\) information about the detailed history outside the allowed predicate. Conversely, an honest reversible realization attains both bounds up to noise \(\varepsilon\).

This theorem is only interesting if \(d\) and \(\delta\) are jointly tighter than applying an existing temporal dimension witness and an existing information–disturbance theorem independently.

## 7. Five alternative cross-disciplinary novelty candidates

### Candidate A — Cryptographic audit-or-coherence receipt

**Nearest prior art.** Mahadev’s post-commitment test architecture; Broadbent–Islam certified deletion; maliciously secure blind delegation with certified deletion; certified deniability; device-independent memory verification.

**Occupied part.** Test-versus-function challenges, compute-then-delete, public deletion certificates, and residual-state simulation are all occupied.

**Potentially new result.** A malicious-server protocol in which the late branch is either a random query about a verifier-driven multi-time computation or a *coherence receipt*: a reference-entanglement test proving approximate inverse-channel behavior, accompanied by composable conditional deletion of every nonallowed historical predicate. The theorem would relate cryptographic soundness to entanglement fidelity rather than merely to complementary-basis deletion.

**Hard kill criterion.** Kill if the coherence receipt is reducible to an existing certified-deletion verification measurement; if it requires a fully trusted server dilation; if the computation output can be prepared classically without executing the committed history; or if the claimed residual security is already exactly certified deniability under different notation.

### Candidate B — Instrument-restricted process-tensor audit/echo frontier

**Nearest prior art.** White *et al.* on unitary-only multi-time witnesses; Taranto–Pollock–Modi on memory strength and recovery; Vieira *et al.* on environment dimension; mirror randomized benchmarking.

**Occupied part.** Unitary-only terminal witnesses, multi-time memory bounds, process recoverability, and mirror inversion are occupied individually.

**Potentially new result.** An SDP hierarchy for the achievable pair \((s_A,F_E)\) when the *same* finite-memory process tensor is subjected after commitment to either a random temporal witness or an inverse/entanglement-fidelity test. A nontrivial separation between bounded classical-memory combs, bounded quantum-memory combs, and unrestricted combs would be the result.

**Hard kill criterion.** Kill if the feasible region factorizes into the Cartesian product of an existing memory-witness bound and an existing recovery bound; if different process ensembles are used in the two branches; or if the separation disappears once entanglement-breaking channels with nonclassical multi-time correlations are admitted.

### Candidate C — Reversible temporal random-access code

**Nearest prior art.** Quantum random-access codes, temporal dimension witnesses, quantum seals, secure software/key leasing, and memory self-testing.

**Occupied part.** Randomly querying one of many encoded facts and returning a protected quantum object are occupied.

**Potentially new result.** A code generated sequentially from verifier-supplied events that permits either one late random fact query or high-fidelity entanglement return after a coherent decoder, with a tight rate region among history length, memory dimension, audit probability, return fidelity, and residual leakage.

**Hard kill criterion.** Kill if the sequential encoder can be compressed to an ordinary static QRAC with no operational loss; if the return branch is simply the identity channel; if the “history” is available at a direct input port; or if quantum-seal bounds already give the claimed rate region after variable substitution.

### Candidate D — Reversible remote history attestation

**Nearest prior art.** Quantum-assisted remote memory attestation, its 2025 security critique, quantum sealing, quantum authentication, and summoning/return tasks.

**Occupied part.** Random memory-address challenges, entanglement-assisted authentication, and read-versus-return verification are occupied.

**Potentially new result.** A verifier injects unpredictable events into a remote reversible computation, then asks either for a checksum/predicate of a random temporal slice or for coherent return of an entangled memory token. A proof would explicitly account for communication, timing, collusion, and hidden storage.

**Hard kill criterion.** Kill if security requires the trusted hardware the protocol claims to avoid; if quantum memory must be copied or hashed in violation of the 2025 attestation limitations; if communication grows linearly with the entire memory and removes any advantage; or if a proxy can answer one branch while a second device answers the other.

### Candidate E — No-postselection causal certification-and-recovery benchmark

**Nearest prior art.** Kim *et al.* on entanglement certification followed by probabilistic recovery; weak-measurement uncollapse; process-tensor witnesses; mirror benchmarking.

**Occupied part.** Sequential certification and probabilistic entanglement recovery are directly occupied.

**Potentially new result.** A deterministic or success-probability-accounted benchmark in which a random late challenge selects either a causal-history witness or a recovery test, and a theorem gives the optimal audit–recovery frontier for a multi-step adaptive memory under realistic noise.

**Hard kill criterion.** Kill if the optimal frontier follows directly from known weak-measurement information–disturbance curves; if postselection is hidden in heralding or data filtering; if adaptivity does not change the observable region; or if the same score is achievable by a single partial measurement with no multi-time memory.

## 8. Minimum viable originality experiment

Before expanding the manuscript, the project should implement a small adversarial benchmark, not a larger narrative simulation.

1. Specify a two- or three-step verifier-driven prefix with fresh random events.
2. Define one audit score and one entanglement-return score, including failure events.
3. State the trusted ports, maximum dimensions, and side-information model.
4. Optimize the joint score over at least these cheats: direct-label, precomputed predicate, bounded classical automaton, entanglement-breaking stepwise channel, hidden one-qubit transcript, postselected weak measurement, and separate audit/echo devices.
5. Compare the resulting bound with a true reversible-memory strategy.
6. Demonstrate a strict gap that disappears if the post-prefix challenge is revealed early.
7. Check whether the gap is merely an instance of a known RAC, seal, certified-deletion, complementary-channel recovery, or process-witness inequality.

If no strict gap survives these baselines, the framing should be killed or reduced to a pedagogical synthesis. If a gap survives but is an immediate numerical intersection of known witnesses, it may be a useful benchmark but not a new foundational theorem.

## 9. Safe and unsafe novelty language

### Safe now

- “We formulate a joint audit-or-recovery task at the intersection of temporal memory certification and coherent recovery.”
- “In the primary-source corpus searched through 12 August 2026, we did not locate a single protocol containing all of our formal clauses.”
- “The proposed contribution, if established, is a cross-branch soundness theorem for one committed multi-time process.”
- “The work combines existing primitives; originality depends on a nonreducible joint bound.”

### Unsafe now

- “the first delayed-choice reveal-or-erase protocol”;
- “the first experiment to erase an observer’s memory”;
- “the first protocol to compute and then certify deletion”;
- “the first recovery of entanglement after certification”;
- “the first witness of quantum causal memory”;
- “the first mirror/echo benchmark with entanglement information”;
- “proof that a real history occurred internally” without a comparison class and access model;
- “proof that no record exists anywhere” from reset of named registers.

## 10. 2025–2026 provisionality and priority watchlist

| Date/status at cutoff | Primary item | Relevance | Treatment |
|---|---|---|---|
| 28 Nov 2025, arXiv preprint | [Ohwada, arXiv:2511.22827](https://arxiv.org/abs/2511.22827) | Memory-delayed eraser proposal claiming a marginal-statistics discriminator without postselection | Provisional; do not rely on its claims, but monitor for priority. |
| 20 Jan 2026, arXiv preprint | [Santos *et al.*, arXiv:2601.14191](https://arxiv.org/abs/2601.14191) | Device-independent temporal-correlation certification of quantum memory; trapped-ion proof of principle | Provisional until peer review; directly relevant to audit soundness. |
| 27 Jan 2026, arXiv preprint | [Strasberg *et al.*, arXiv:2601.19703](https://arxiv.org/abs/2601.19703) | Approximate records and recoherence in isolated-system decoherent histories | Provisional; directly occupies records–recoherence language. |
| 28 Jan 2026, peer reviewed | [Ohst *et al.*, Quantum 10, 1988](https://doi.org/10.22331/q-2026-01-28-1988) | Classical/quantum memory dimension in adaptive discrimination | Treat as established primary literature. |
| 13 Jan 2026, peer reviewed | [Asadian, Gams, and Sponar, PR Research 8, L012011](https://doi.org/10.1103/llgb-gql9) | Covariant correlation--disturbance trade-off for sequential \(n\)-outcome measurements in arbitrary dimension | Treat as established, directly adjacent measurement trade-off literature; its scores differ from the candidate's. |
| 15 Jan 2026, arXiv preprint | [Dey and Safavi-Naini, arXiv:2601.10542](https://arxiv.org/abs/2601.10542) | Hybrid encryption with certified deletion | Provisional; reinforces that deletion primitives are still advancing. |
| 8 Mar 2026, arXiv preprint | [Murshid, Sarkar, and Mandal, arXiv:2603.07646](https://arxiv.org/abs/2603.07646) | Publicly verifiable certified deletion in registered ABE | Provisional; not a physical recoherence certificate. |
| 7 May 2026, arXiv preprint; v2 3 Jul | [Torvinen, Keski-Vakkuri, and Pranzini, arXiv:2605.06848](https://arxiv.org/abs/2605.06848) | Petz recovery of einselected information from Darwinist environment fragments | Provisional; it concerns state reconstruction, not the proposed coalition inverse, but narrows the framing. |
| 4 Jun 2026, arXiv preprint; v2 24 Jun | [Girard, Cheng, and Cao, arXiv:2606.06588](https://arxiv.org/abs/2606.06588) | Identifies Darwinist objectivity with local recoverability of operator-algebra QEC codes and tracks fragment-supported logical algebras | Provisional and the sharpest collision with the access-structure/hypergraph proposal. |
| 3 Jun 2026, arXiv preprint | [Krajenbrink *et al.*, arXiv:2606.04843](https://arxiv.org/abs/2606.04843) | Extends decoded quantum interferometry beyond Hamming space using translation association schemes and shell-level spectral reductions | Provisional; adjacent mathematical methodology, not an information--recovery measurement theorem. |
| 17 Jun 2026, arXiv preprint | [Zonnios and Binder, arXiv:2606.19511](https://arxiv.org/abs/2606.19511) | Recurrent process testers that retain a full classical outcome record and a bounded coherent memory; the \(d_A=1\) endpoint is directly adjacent to the online no-coherent-memory model | Provisional; process discrimination only, but a strong priority warning against claiming the memory architecture. |
| QCrypt 2026 accepted talk; underlying 7 Nov 2024 preprint | [Çakan, Goyal, and Raizes, arXiv:2411.05176](https://arxiv.org/abs/2411.05176) | Simulation-based deletion leaving no usable evidence | Do not call it peer-reviewed archival publication on the basis of conference acceptance alone. |
| 18 Jun 2026, arXiv preprint | [Pettugani *et al.*, arXiv:2606.20123](https://arxiv.org/abs/2606.20123) | Bell-pair injection plus mirror randomized benchmarking and mutual-information tracking | Provisional; directly narrows the recovery-benchmark gap. |
| 4 Aug 2026, arXiv preprint | [Maity, Onggadinata, and Koh, arXiv:2608.03944](https://arxiv.org/abs/2608.03944) | Exact model tradeoff between Darwinistic redundancy and post-recovery logical fidelity, plus a general information bound | Provisional and only eight days old at the cutoff; independently verify before relying on its theorem. |
| Accepted 14 Jul 2026; arXiv posted 7 Aug | [Zhang *et al.*, accepted PRA](https://doi.org/10.1103/yxq6-n2nr), [arXiv:2608.06726](https://arxiv.org/abs/2608.06726) | Fully device-independent unsharp-instrument characterization using an entanglement-assisted sequential QRAC | Accepted primary work, but no version-of-record article was available at the cutoff; theorem details and pagination remain provisional. |

No 2026 preprint in this table should support a definitive theorem claim. It may establish priority risk or motivate a comparison, with its status and version date stated explicitly.

## 11. Claim-to-source ledger

| Technical claim used in this audit | Primary source |
|---|---|
| A random late choice can select complete path readout or erasure with conditional interference | [Scarcelli, Zhou, and Shih (2007)](https://doi.org/10.1140/epjd/e2007-00164-y) |
| The eraser choice can be causally disconnected from interferometer events | [Ma *et al.* (2013)](https://doi.org/10.1073/pnas.1213201110) |
| A Wigner-friend memory-erasure context has been explicitly analyzed | [Elouard *et al.* (2021)](https://doi.org/10.22331/q-2021-07-08-498) |
| Partial measurement can be conditionally reversed in experiment | [Katz *et al.* (2008)](https://doi.org/10.1103/PhysRevLett.101.200401) |
| Entanglement can be weakly certified and probabilistically recovered by reversal | [Kim *et al.* (2023)](https://doi.org/10.1126/sciadv.adi5261) |
| Quantum ciphertext deletion can be classically certified | [Broadbent and Islam (2020)](https://doi.org/10.1007/978-3-030-64381-2_4) |
| Device-independent composable certified deletion exists | [Kundu and Tan (2023)](https://doi.org/10.22331/q-2023-07-06-1047) |
| FHE computation can be combined with proof of deletion | [Poremba (2023)](https://doi.org/10.4230/LIPIcs.ITCS.2023.90) |
| Fully malicious blind delegation with certified deletion exists | [Bartusek *et al.* (2024)](https://doi.org/10.1007/978-3-031-58737-5_4) |
| Post-deletion residual state can be required to have a simulation-based no-trace property | [Çakan, Goyal, and Raizes, preprint](https://arxiv.org/abs/2411.05176) |
| A quantum prover can be cryptographically committed before a random basis/test challenge | [Mahadev (2018)](https://doi.org/10.1109/FOCS.2018.00033) |
| Any requested bit can be recovered from a quantum random-access encoding, subject to lower bounds | [Ambainis *et al.* (2002)](https://doi.org/10.1145/581771.581773) |
| Quantum memory in a multi-time process can be witnessed by process entanglement | [Giarmatzi and Costa (2021)](https://doi.org/10.22331/q-2021-04-26-440) |
| Multi-time processes that transmit only classical information forward form a conditional-instrument hierarchy distinct from both quantum-memory and memoryless processes | [Taranto *et al.* (2024)](https://doi.org/10.22331/q-2024-05-02-1328) |
| The quantum-memory cost of a multi-step strategy can be defined with classical-memory assistance and generally cannot be optimized independently at each step | [Bisio *et al.* (2012)](https://doi.org/10.1103/PhysRevA.85.032333) |
| Sequential QRAC statistics can separate classical-memory processes from quantum-memory processes under stated semi-device-independent assumptions | [Roy *et al.* (2024)](https://doi.org/10.1103/PhysRevA.110.012608) |
| Adaptive channel discrimination with bounded coherent and classical memory can be treated through constrained-separability optimization | [Ohst *et al.* (2026)](https://doi.org/10.22331/q-2026-01-28-1988) |
| Recurrent process testers with a full classical outcome record and bounded coherent memory were proposed in a June 2026 preprint | [Zonnios and Binder (2026 preprint)](https://arxiv.org/abs/2606.19511) |
| Temporal correlations can lower-bound effective environment dimension via an SDP hierarchy | [Vieira *et al.* (2024)](https://doi.org/10.22331/q-2024-01-10-1224) |
| Multi-time entanglement and non-Markovianity can be bounded with unitary controls plus terminal measurement | [White *et al.* (2025)](https://doi.org/10.22331/q-2025-04-08-1695) |
| Non-Markovian memory strength bounds process recoverability | [Taranto, Pollock, and Modi (2021)](https://doi.org/10.1038/s41534-021-00481-4) |
| An entanglement-breaking step can still support nonclassical multi-time correlations | [Vieira, Ku, and Budroni (2025)](https://doi.org/10.1103/r8lf-bb4p) |
| Identity-channel/quantum-memory behavior can be robustly self-tested | [Sekatski *et al.* (2023)](https://doi.org/10.1103/PhysRevLett.131.170802) |
| Mirror circuits give a scalable run-and-invert benchmark | [Proctor *et al.* (2022)](https://doi.org/10.1103/PhysRevLett.129.150502) |
| Readability and later detection of reading obey no-go bounds for quantum seals | [Kimmel and Kolkowitz (2019)](https://doi.org/10.1103/PhysRevA.100.052326) |
| Information leakage to the complementary channel controls recoverability | [Kretschmann, Schlingemann, and Werner (2008)](https://doi.org/10.1109/TIT.2008.917696) |
| Approximate error correction admits an environment-dual characterization of optimal recovery | [Bény and Oreshkov (2010)](https://doi.org/10.1103/PhysRevLett.104.120501) |
| Pure-state quantum secret sharing has authorized/unauthorized complement duality and is equivalent to erasure correction | [Cleve, Gottesman, and Lo (1999)](https://doi.org/10.1103/PhysRevLett.83.648) |
| Privacy for a channel is correctability for its complementary channel, also approximately | [Kretschmann, Kribs, and Spekkens (2008)](https://doi.org/10.1103/PhysRevA.78.032330) |
| Approximate QSS reconstructability is equivalent to small complementary-channel leakage under the paper's capacity metric | [Ouyang *et al.* (2023)](https://doi.org/10.1103/PhysRevA.108.012425) |
| Coherence recovery with access to only part of an environment has already been optimized for random pure states | [Miatto *et al.* (2015)](https://doi.org/10.1103/PhysRevA.92.062331) |
| Spectrum-broadcast structures encode redundant distinguishable pointer records in environment fragments | [Korbicz, Horodecki, and Horodecki (2014)](https://doi.org/10.1103/PhysRevLett.112.120402) |
| Darwinist objectivity has been mapped to locally recoverable logical subalgebras of QEC codes | [Girard, Cheng, and Cao, 2026 preprint](https://arxiv.org/abs/2606.06588) |
| A late coin flip between full path discrimination and phase discrimination has a tight arbitrary-\(N\) operational duality region | [Bagan *et al.* (2018)](https://doi.org/10.1103/PhysRevLett.120.050402) |
| Partial/set-valued Ways and Phases games and a mutual-information duality bound have been studied | [Hillery (2021)](https://doi.org/10.1088/1751-8121/ac367d) |
| Universal pure-state estimation and operation fidelity obey a tight arbitrary-dimensional trade-off | [Banaszek (2001)](https://doi.org/10.1103/PhysRevLett.86.1366) |
| Guessing complementary measurement outcomes can equal recoverable entanglement fidelity | [Berta, Coles, and Wehner (2014)](https://doi.org/10.1103/PhysRevA.90.062127) |
| Classical information obtained from an environment can condition quantum correction | [Gregoratti and Werner (2003)](https://doi.org/10.1080/0950034021000058021) |
| Partial environment access and conditional feedback have been optimized using entanglement fidelity | [Memarzadeh, Macchiavello, and Mancini (2011)](https://doi.org/10.1088/1367-2630/13/10/103031) |
| Sequential QRAC score pairs can yield tight information--disturbance frontiers and instrument self-tests | [Mohan, Tavakoli, and Brunner (2019)](https://doi.org/10.1088/1367-2630/ab3773) |
| Maximum/minimum quantum causal effects connect classical causal influence, approximate invertibility, and a causal monogamy relation | [Chiribella and Goswami (2025)](https://doi.org/10.1103/PRXQuantum.6.020335) |
| Local information gain and nondestructive preservation of entangled states obey explicit trade-offs, including adaptive entanglement-assisted strategies | [Lim, Hhan, and Kwon (2025)](https://doi.org/10.1088/2058-9565/adc034) |
| Finite-state stochastic and quantum transducers for repeatedly measured input--output processes predate the candidate | [Wiesner and Crutchfield (2008)](https://doi.org/10.1016/j.physd.2008.01.021) |
| A 2026 peer-reviewed covariant correlation--disturbance relation covers sequential \(n\)-outcome measurements in arbitrary dimension | [Asadian, Gams, and Sponar (2026)](https://doi.org/10.1103/llgb-gql9) |
| Extracted classical distinguishability and retained coherence obey a resource-theoretic trade-off for orthogonal pure-state ensembles | [Liu *et al.* (2025)](https://doi.org/10.1103/PhysRevA.111.062215) |
| The single-chain ordered/NRT relation underlying an LCP shell decomposition is established association-scheme mathematics | [Martin and Stinson (1999)](https://doi.org/10.4153/CJM-1999-017-5) |
| Iterated wreath products of one-class schemes and their Terwilliger algebras are established | [Bhattacharyya, Song, and Tanaka (2010)](https://doi.org/10.1007/s10801-009-0196-x) |
| Partial Ways/Phases games allow set-valued path or phase answers but do not use the committed random-prefix/return score | [Hillery (2021)](https://doi.org/10.1088/1751-8121/ac367d) |
| Wave--particle complementarity has been extended to hierarchical bipartite systems, including a binary-tree representation | [Wu and Wang (2020)](https://doi.org/10.3390/e22080813) |
| Translation association schemes and shell-level eigenproblems entered decoded quantum interferometry in a June 2026 preprint | [Krajenbrink *et al.* (2026 preprint)](https://arxiv.org/abs/2606.04843) |

## 12. Search method and limits

The audit searched publisher pages, DOI records, arXiv, IACR ePrint, Quantum, APS journals, Nature-family journals, PNAS/PMC, ACM/IEEE records, and citation trails. Query families included combinations of:

- delayed-choice quantum eraser, random reveal/erase, quantum memory eraser;
- weak-measurement reversal, entanglement certification recovery, uncollapse;
- certified deletion, blind delegation, secure software leasing, certified deniability;
- test-or-compute, test/Hadamard round, cut-and-choose, quantum commitment challenge;
- random-access code, dimension witness, temporal witness, causal memory;
- process tensor, unitary-only multi-time witness, memory strength, recoverability;
- mirror randomized benchmarking, Loschmidt echo, Bell-pair injection;
- quantum seal, read disturbance, state return, remote memory attestation;
- information–disturbance, complementary channel, entanglement-breaking recovery;
- multipath Ways/Phases games, partial path information, \(\ell_1\) coherence, and QND fidelity trade-offs;
- coherence/distinguishability resource theories, co-bits, retained coherence after discrimination, and decoded quantum interferometry beyond Hamming space;
- random-coordinate and sequential QRAC disturbance, recoverable entanglement fidelity, and environment-assisted correction;
- Boolean and \(q\)-ary Hamming schemes, association schemes, coherent configurations, group-covariant instruments, and rank-one spectral transitions.
- rooted-tree automorphisms, longest-common-prefix and ultrametric orbitals, ordered/NRT Hamming schemes, iterated wreath products, random-prefix discrimination, nested partitions, and partial/hierarchical path information.
- online and classically adaptive quantum instruments, zero coherent memory, comb/tester memory cost, causal hidden-state factorizations, recurrent testers, sequential QRAC memory witnesses, local nondestructive discrimination, and transcript-conditioned entanglement recovery.

Searches were run through 12 August 2026. Technical conclusions above point to primary papers. Secondary indexes were used only to locate metadata when a publisher page was difficult to query. The search is adversarial but not a registered systematic review: terminology drift, unpublished manuscripts, and nonindexed conference material remain possible. The correct epistemic statement is “no exact conjunction was located in this search,” not “the conjunction has never been proposed.”

## 13. Gate decision

**Current decision: amber, with a narrow go condition.**

Proceed only if the next theoretical unit is the formal two-challenge comb game and it produces a nonfactorizing soundness or rate-region result. Stop or reframe if any of the following occurs:

1. the audit branch is just a static QRAC or a known temporal dimension witness;
2. the recovery branch is just mirror benchmarking, a quantum eraser, or weak-measurement uncollapse;
3. deletion means only certified deletion or certified deniability under renamed registers;
4. the joint score is the product or intersection of two independent known tests;
5. a direct-label or precompiled-predicate strategy passes;
6. hidden transcripts are excluded only by assumption but the claim is phrased as physical proof of global erasure;
7. the history inputs are known before commitment;
8. the only novelty is interpretive language about branches or observers.

**Separate static-mathematics gate:** red for the late audit/coherence game, amber for the exact arbitrary-\(n\) random-coordinate frontier, amber for the finite-\(n\) committed-prefix frontier, and amber-red for the large-\(n\) \(q\)-ary kink alone. The tree/ultrametric geometry itself is red as a novelty claim. The static line should proceed only if it proves global optimality beyond the covariant ansatz and delivers a distinct tight inequality, finite-size/localization result, or genuine sequential soundness not implied by known multipath duality and association-scheme machinery.

**Separate online-memory gate:** red for introducing classical-memory-only combs, classically adaptive testers, or sequential QRAC memory witnesses; amber for the exact causally factored audit--return frontier; and amber-green only for a demonstrated strict separation from both unrestricted collective instruments and bounded-coherent-memory combs. The minimal publishable unit is a two- or three-round theorem with all storage and decoder ports charged explicitly.

The strongest defensible research question is no longer whether a reversible history can be queried and erased. It is whether an unpredictable post-prefix challenge can turn *the same verifier-driven multi-time process* into either a causal-memory certificate or a coherence-return certificate, with a joint adversarial theorem that is not reducible to certified deletion, quantum sealing, temporal dimension witnessing, or complementary-channel recovery.

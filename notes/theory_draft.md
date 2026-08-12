# Theoretical Core: Causal Certification of Reversible Quantum Histories

**Working document:** v0.2
**Cutoff date:** 12 August 2026
**Provisional program name:** **Causal Certification of Reversible Quantum Histories (CCRQH)**
**Proposed primitive:** **Zero-Record Causal Proof (ZRCP)**

This document intentionally uses no public-facing name and assumes no interpretation of quantum mechanics. Its subject is an operational task: certify nontrivial multi-time causal structure, extract a limited function of a history, and recover coherence after decoupling every record that could identify the fine-grained history.

---

## 0. Epistemic labels

Every formal claim is marked with one of the following labels.

- **[K — known]:** a standard result or a direct consequence of published work. No novelty is claimed.
- **[D — derived here]:** a statement proved in this document from the declared model. It may organize the project usefully, but priority is not claimed until a systematic prior-art search is complete.
- **[P — proposal]:** a working definition, quantity, or protocol. Its value and novelty remain research hypotheses.
- **[C — conjecture]:** an unproved statement that requires a proof, counterexample, or reduction to prior literature.
- **[A — assumption]:** an operational or trust restriction without which the indicated conclusion is not identifiable.

The words *history*, *memory*, and *agent* are functional. They imply neither consciousness nor collapse, and they carry no branch ontology. All predictions below belong to standard quantum mechanics.

---

## 1. Technical identity and falsifiable thesis

### 1.1 Scope

CCRQH studies families of multi-time quantum processes with four approximate properties:

1. a coherent logical variable controls alternatives with distinguishable internal dynamics;
2. the alternatives contain causal memory that can be certified under a declared intervention model;
3. a limited function of the process is transferred to a phase or acceptance bit;
4. degrees of freedom containing any additional history information are reset or decoupled before final interference.

The central question is:

> Under which assumptions about access, locality, and resources can an experiment certify that a process had nontrivial internal causal memory while its final accessible state retains no transcript of the fine-grained history beyond an authorized output?

### 1.2 Claims outside scope

The protocol does not select a quantum interpretation, recover a naturally decohered macroscopic alternative, or copy arbitrary outcomes from incompatible alternatives. An ideal implementation is ordinary unitary evolution. Any original content must reside in the **certification task**, **access model**, **resource bound**, or **complexity separation**, not in renaming phase kickback.

### 1.3 Minimum publishable unit

A defensible contribution must provide at least one of the following:

- a ZRCP definition not operationally reducible to ordinary circuit fidelity;
- a multi-time witness excluding an explicit class of processes without the claimed causal memory;
- a resource bound involving memory, depth, locality, or access;
- a separation between coherent process access and independent classical samples;
- a no-go theorem delimiting which properties can survive without a transcript.

---

## 2. Kinematic model

All Hilbert spaces are finite-dimensional. Continuous-variable extensions are outside the first model.

### 2.1 Registers

Let

\[
\mathcal H_{\mathrm{all}}
=\mathcal H_B\otimes\mathcal H_S\otimes\mathcal H_M
\otimes\mathcal H_G\otimes\mathcal H_X\otimes\mathcal H_Q
\otimes\mathcal H_E .
\]

- \(B\): coherent control or logical alternative label;
- \(S\): accessible system, world, or intervention interface;
- \(M\): internal memory persistent across time steps;
- \(G\): reversible work registers and computational garbage;
- \(X\): verifier challenge, classical or coherent;
- \(Q\): response, acceptance, or phase-readout qubit;
- \(E\): uncontrolled or unrecovered degrees of freedom.

Write \(R=SMG\) for internal registers intended to be reset, and \(L=RE\) for all systems that can retain path information outside \(B\). Which physical degrees belong to \(E\) is part of the device model and must be audited experimentally.

### 2.2 Channels, instruments, and coherent records

A channel is completely positive and trace preserving (CPTP). An instrument
\(\{\mathcal A^{a|x}\}_a\) is a family of completely positive, trace-nonincreasing maps whose sum is CPTP. Here \(x\) is an intervention and \(a\) its outcome.

Every CPTP map has a Stinespring dilation. A measurement or decision can therefore be represented coherently by an isometry,

\[
V_x|\psi\rangle_S|0\rangle_M|0\rangle_G
=\sum_a K_{a|x}|\psi\rangle_S|a\rangle_M|g_{a,x}(\psi)\rangle_G ,
\]

provided its outcome record and dilation environment are retained. **[K]** A measurement whose outcome is sent to an external classical controller is not reversible from \(SMG\) unless that external record is independent of the alternative or is itself included in the inverse.

### 2.3 Controlled histories

Relative to a declared basis \(\{|b\rangle_B\}_{b=0}^{K-1}\), an ideal controlled history is an isometry

\[
V_H=\sum_{b=0}^{K-1}|b\rangle\!\langle b|_B\otimes V_b,
\qquad
V_b:\mathcal H_{R,\mathrm{in}}\rightarrow
\mathcal H_{R,\mathrm{out}}\otimes\mathcal H_E .
\]

For a pure reference input,

\[
V_H\left(\sum_b\alpha_b|b\rangle_B\right)|0,e_0\rangle_{RE}
=\sum_b\alpha_b|b\rangle_B|r_b\rangle_R|e_b\rangle_E .
\]

The word *alternative* is always relative to this basis, the chosen interface, and the declared coarse-graining.

### 2.4 Class operators and histories

For a property sequence \(h=(a_1,\ldots,a_T)\), one may define

\[
C_h=P^{(T)}_{a_T}U_{T:T-1}\cdots
P^{(1)}_{a_1}U_{1:0}
\]

and the decoherence functional

\[
D(h,h')=\operatorname{Tr}\!\left(C_h\rho_0 C_{h'}^\dagger\right).
\]

**[K]** Only when the relevant off-diagonal terms \(D(h,h')\) vanish approximately can the family be treated as classical mutually exclusive histories with ordinary probabilities. CCRQH deliberately recovers some interference terms. It must therefore not claim simultaneously that all paths are a decoherent consistent family. Class operators are optional notation; the primary object below is an operational multi-time process.

---

## 3. Multi-time processes and causal memory

### 3.1 Process tensor

Let \(t_0<t_1<\cdots<t_T\). At slot \(j\), an experimenter inserts
\(\mathcal A_j^{a_j|x_j}:\mathsf L(S_j^{\mathrm i})\to
\mathsf L(S_j^{\mathrm o})\). A multi-time quantum process is a multilinear map

\[
\mathbf T_{T:0}[\mathcal A_{T-1},\ldots,\mathcal A_0]=\rho_T
\]

that returns a valid state for all completely positive interventions and obeys temporal causality. Its Choi representation is a positive operator
\(\Upsilon_{T:0}\ge0\). With a fixed Choi convention,

\[
p(\mathbf a|\mathbf x)
=\operatorname{Tr}\!\left[
\Upsilon_{T:0}
\left(
M_{T-1}^{a_{T-1}|x_{T-1}}\otimes\cdots\otimes
M_0^{a_0|x_0}
\right)^{\mathsf T}
\right].
\]

The transpose is convention-dependent. Causality imposes recursive trace constraints of the schematic form

\[
\operatorname{Tr}_{S_j^{\mathrm o}}\Upsilon_{j:0}
=I_{S_j^{\mathrm i}}\otimes\Upsilon_{j-1:0},
\]

with the appropriate initial normalization. **[K]** This is the established quantum-comb/process-tensor framework.

### 3.2 Realization memory

A process is Markovian relative to the declared temporal partition when its Choi operator factorizes into independent step channels with the appropriate link convention,

\[
\Upsilon^{\mathrm{Markov}}_{T:0}
=\rho_0\otimes J(\Lambda_{1:0})\otimes\cdots
\otimes J(\Lambda_{T:T-1}).
\]

A realization with memory dimension \(d_M\) consists of sequential channels

\[
\Lambda_j:S_j^{\mathrm o}\otimes M_j
\longrightarrow S_{j+1}^{\mathrm i}\otimes M_{j+1},
\qquad \dim M_j\le d_M .
\]

For an intervention set \(\mathfrak I\), define the realization-memory dimension

\[
d_{\min}(\Upsilon;\mathfrak I)
=\min\{d_M:\text{a realization with memory }d_M
\text{ reproduces all statistics in }\mathfrak I\}.
\]

**[P]** This quantity is relative to the interface, intervention set, and trust model. It does not identify a unique microscopic architecture. Existing temporal dimension witnesses already bound closely related quantities; a novelty claim would require proving exactly what differs here.

### 3.3 Reversible agent

A length-\(T\) reversible agent is an isometric realization with explicit persistent memory,

\[
V_j^{(x_j)}:
S_j^{\mathrm i}\otimes M_j\otimes G_j
\longrightarrow
S_j^{\mathrm o}\otimes M_{j+1}\otimes G_{j+1},
\]

admitting a controlled inverse after retaining every challenge \(x_j\) and every datum needed for reversibility. At least one action must change a later input to the same agent; otherwise *agent* may be only a new name for a detector or ancilla.

A coherent adaptive update can be written

\[
|s_j\rangle_S|m_j\rangle_M|0\rangle_A
\mapsto
|s_j\rangle_S|m_{j+1}(s_j,m_j)\rangle_M
|\pi(s_j,m_j)\rangle_A ,
\]

followed by world dynamics controlled by \(A\). Sufficient information must remain to invert any logically noninjective update.

### 3.4 Relative causal certification

Let \(\mathfrak C_0\) be a declared null class: for example Markovian processes, realizations with \(d_M<d\), circuits lacking the causal edge \(M_j\to S_{j+1}\), or circuits forbidden from receiving a late input through a given port. A behavior \(p(\mathbf a|\mathbf x)\) certifies structure outside \(\mathfrak C_0\) when

\[
p\notin\mathcal P(\mathfrak C_0,\mathfrak I),
\]

where \(\mathcal P(\mathfrak C_0,\mathfrak I)\) is the set of behaviors produced by the null class under the allowed interventions.

If this set is convex and closed and \(p\) lies outside it, finite-dimensional separation gives a linear witness

\[
W[p]=\sum_{\mathbf a,\mathbf x}c_{\mathbf a,\mathbf x}
p(\mathbf a|\mathbf x)
\]

such that

\[
W[p]>\beta_0:=
\sup_{q\in\mathcal P(\mathfrak C_0,\mathfrak I)}W[q].
\]

**[K]** This is convex separation, not a new physical theorem. The scientific work is to choose a meaningful null class, compute a robust bound, and avoid circular assumptions.

---

## 4. Fine histories, predicates, and physical transcripts

### 4.1 Fine-grained history variable

Let \(H\) be a classical label for a fine-grained history in a finite set
\(\mathscr H\). Operationally, \(H\) may encode a sequence of world states, memory updates, challenges, actions, and work-register values at the declared resolution. It is a bookkeeping variable for an ensemble of interventions; it is not assumed to be a simultaneously measurable hidden trajectory.

For two alternatives with residual states \(\rho_T^0,\rho_T^1\), define

\[
D_T=\frac12\|\rho_T^0-\rho_T^1\|_1.
\]

For equal priors the optimal guessing probability is
\(P_{\mathrm{guess}}=(1+D_T)/2\). Orthogonal supports give a perfect record; identical states give no accessible record.

### 4.2 Noninjective authorized coarse-graining

Let

\[
p:\mathscr H\to\mathscr Z
\]

be the authorized predicate or coarse-graining. Its fibers are

\[
\mathscr H_z=\{h\in\mathscr H:p(h)=z\}.
\]

**Nonvacuity requirement [P].** A privacy claim about information beyond \(p(H)\) is meaningful only if the tested distribution has at least two distinguishable fine histories with nonzero probability in at least one fiber:

\[
\exists z,\quad
|\operatorname{supp}(H|p(H)=z)|\ge2.
\]

Preferably every tested fiber should satisfy this condition. If \(p\) is injective on the support—for example two histories with \(p(h_b)=b\)—then

\[
I(H:T\mid p(H))=0
\]

for every transcript \(T\), even one containing a complete copy of \(H\). The conditional mutual information is then vacuous. A credible benchmark therefore needs \(K\ge4\), a separate microhistory variable \(\mu\), or another noninjective coarse-graining such as parity:

\[
H=(b,\mu),\qquad p(H)=b\bmod2
\quad\text{or}\quad p(H)=f(\mu)
\]

with several values of \((b,\mu)\) per predicate value.

### 4.3 Residual transcript

Let \(T\) contain every final accessible register except the intended interferometric control and the authorized output \(Z\). For a fine history \(h\), let its residual state be \(\tau_T^h\). Define within-fiber leakage

\[
\epsilon_{\mathrm{fib}}
=\max_z\max_{h,h'\in\mathscr H_z}
\frac12\|\tau_T^h-\tau_T^{h'}\|_1 .
\]

This tests whether histories sharing the same authorized output remain distinguishable. It is meaningful only under the nonvacuity requirement.

An absolute zero-record condition is stronger: there must exist a history-independent \(\sigma_T\) such that

\[
\epsilon_{\mathrm{ZR}}
=\inf_{\sigma_T}\max_{h\in\mathscr H}
\frac12\|\tau_T^h-\sigma_T\|_1
\le\epsilon .
\]

Absolute zero record permits the authorized output only in the explicitly excluded register \(Z\). A scheme that meets only the first condition must be called **predicate-conditioned zero record**.

### 4.4 Information-theoretic privacy and its traps

For a noninjective predicate and a declared prior, one may report

\[
I(H:T\mid Z),\qquad Z=p(H).
\]

This is a useful diagnostic, but it is prior-dependent, noncomposable by itself, and vacuous for injective \(p\). It must therefore be accompanied by fiber cardinalities and a worst-case trace-distance or channel criterion.

Let \(\mathcal T:H\to TZ\) be the induced transcript channel on an orthogonal encoding of histories and let
\(\mathcal P:H\to Z\) return only \(p(h)\). Physical transcript privacy holds with error \(\epsilon_{\mathrm{priv}}\) if a simulator channel
\(\mathcal S:Z\to TZ\) exists such that

\[
\frac12\|\mathcal T-\mathcal S\circ\mathcal P\|_\diamond
\le\epsilon_{\mathrm{priv}}.
\]

**[P]** This is inspired by simulation-based cryptography but is not automatically cryptographic zero knowledge. It omits computational efficiency, malicious auxiliary inputs, concurrency, and full composability. Unless those are proved, the correct phrase is *physical transcript privacy*.

### 4.5 Coherent history privacy

Classical labels do not test preservation of superpositions. Let
\(\mathcal N:B\to BT\) be the actual pre-readout channel and let the ideal be

\[
\mathcal N_{\mathrm{id}}(\rho_B)
=U_p\rho_BU_p^\dagger\otimes\sigma_T,
\qquad
U_p=\sum_b e^{i\phi p_b}|b\rangle\!\langle b|.
\]

Define

\[
\epsilon_{\diamond}
=\frac12\|\mathcal N-\mathcal N_{\mathrm{id}}\|_\diamond .
\]

This tests inputs entangled with a reference and simultaneously captures phase action, transcript decoupling, and logical recovery. It is the preferred ideal metric. Large devices may require lower bounds from entanglement fidelity, complementary-basis tests, or restricted witnesses, with the unverified remainder stated explicitly.

---

## 5. Visibility, leakage, and reset

### 5.1 Visibility

For two balanced alternatives, inserting an analysis phase \(\theta\) and recombining gives

\[
P_0(\theta)=\frac12[1+V\cos(\theta+\varphi)].
\]

If \(\rho_B\) has populations \(1/2\),

\[
V=2|\langle0|\rho_B|1\rangle|.
\]

For unequal populations, use normalized visibility

\[
V_N=
\frac{|(\rho_B)_{01}|}
{\sqrt{(\rho_B)_{00}(\rho_B)_{11}}}
\]

when the denominator is nonzero.

### 5.2 Pure leakage lemma **[K]**

**Lemma 1.** Let

\[
|\Psi\rangle_{BL}
=\frac{|0\rangle|\ell_0\rangle+
e^{i\phi}|1\rangle|\ell_1\rangle}{\sqrt2},
\]

with normalized leakage states. Then

\[
V=|\langle\ell_1|\ell_0\rangle|.
\]

If \(L\) perfectly distinguishes the alternatives, \(V=0\). If \(V=1\), the two leakage states differ only by a phase and no measurement on \(L\) reveals the alternative.

**Proof.** Tracing out \(L\) makes the off-diagonal entry of
\(\rho_B\) equal to
\(e^{-i\phi}\langle\ell_1|\ell_0\rangle/2\). The visibility formula follows. Pure states are perfectly distinguishable exactly when orthogonal, and unit overlap means the same ray. \(\square\)

### 5.3 Distinguishability–visibility bound **[K]**

For mixed residual states define root fidelity

\[
f(\rho,\sigma)=\|\sqrt\rho\sqrt\sigma\|_1.
\]

The recovered visibility cannot exceed the appropriate fidelity of the unresolved leakage states. The Fuchs–van de Graaf inequalities give

\[
D_L\le\sqrt{1-f(\rho_L^0,\rho_L^1)^2}.
\]

Hence, whenever \(V\le f\),

\[
D_L^2+V^2\le1.
\]

This is known complementarity/information disturbance. A project-specific bound must add nontrivial constraints involving memory, locality, depth, or allowed recovery maps.

### 5.4 State reset fidelity

For a fixed input and intended-reset registers \(R=SMG\), define squared fidelity

\[
F_{\mathrm{reset}}^{\mathrm{state}}
=F(\rho_R^{\mathrm{fin}},|0\rangle\!\langle0|_R),
\qquad
F(\rho,\sigma)=\|\sqrt\rho\sqrt\sigma\|_1^2.
\]

It must be reported per challenge and logical input, not only as an average. It is necessary but not sufficient for channel recovery or coherence.

### 5.5 Perfect visible reset need not imply coherence **[D]**

**Proposition 1.** \(F_{\mathrm{reset}}^{\mathrm{state}}=1\) on \(R\) does not imply nonzero visibility.

**Proof by counterexample.** Let

\[
|\Psi\rangle_{BRE}
=\frac{|0\rangle_B|0\rangle_R|e_0\rangle_E+
|1\rangle_B|0\rangle_R|e_1\rangle_E}{\sqrt2},
\qquad
\langle e_0|e_1\rangle=0.
\]

\(R\) is exactly reset, but \(E\) contains a perfect transcript. Lemma 1 gives \(V=0\). \(\square\)

Thus “all measured ancillas returned zero” is not a recoherence certificate. One must test complementary observables of \(B\), use a matched identity echo, and bound complementary-channel leakage.

### 5.6 Joint fidelity implies a contrast bound **[D]**

**Proposition 2.** Suppose

\[
F\!\left(
\rho_{BL},
|\psi_\phi\rangle\!\langle\psi_\phi|_B\otimes\sigma_L
\right)\ge1-\varepsilon,
\qquad
|\psi_\phi\rangle=(|0\rangle+e^{i\phi}|1\rangle)/\sqrt2 .
\]

Then, for

\[
X_\phi=e^{-i\phi}|0\rangle\!\langle1|
+e^{i\phi}|1\rangle\!\langle0|,
\]

\[
\operatorname{Tr}[(X_\phi\otimes I)\rho_{BL}]
\ge1-2\sqrt\varepsilon .
\]

**Proof.** Fuchs–van de Graaf bounds the trace distance by
\(\sqrt\varepsilon\). Since
\(\|X_\phi\otimes I\|_\infty=1\), the expectation differs from its ideal value one by at most twice the trace distance. \(\square\)

This proposition assumes a joint bound including every relevant leakage degree; a marginal reset fidelity does not satisfy the premise.

---

## 6. Minimal phase protocol

### 6.1 Circuit

Prepare

\[
|\Psi_0\rangle
=\left(\sum_b\alpha_b|b\rangle_B\right)
|0\rangle_R|-\rangle_Q,
\qquad
|-\rangle=(|0\rangle-|1\rangle)/\sqrt2.
\]

Let

\[
U_H|b\rangle_B|0\rangle_R
=|b\rangle_B|h_b\rangle_R
\]

and let a clean reversible predicate circuit satisfy

\[
C_p|h_b\rangle_R|q\rangle_Q
=|h_b\rangle_R|q\oplus p(h_b)\rangle_Q ,
\]

including all predicate garbage within \(R\).

### 6.2 Phase-survival lemma **[K]**

**Lemma 2.**

\[
U_H^\dagger C_pU_H
\left(\sum_b\alpha_b|b\rangle_B\right)|0\rangle_R|-\rangle_Q
=
\left(\sum_b\alpha_b(-1)^{p(h_b)}|b\rangle_B\right)
|0\rangle_R|-\rangle_Q .
\]

**Proof.** On each \(b\) term, \(C_p\) applies
\(X^{p(h_b)}\) to \(|-\rangle\), and
\(X|-\rangle=-|-\rangle\). Then
\(U_H^\dagger|b,h_b\rangle=|b,0\rangle\). Linearity completes the proof. \(\square\)

This is standard phase kickback plus uncomputation. It is not a novel theorem.

### 6.3 Readout

For \(B=(|0\rangle+|1\rangle)/\sqrt2\), the final state is, up to global phase,

\[
\frac{|0\rangle+
(-1)^{p(h_0)\oplus p(h_1)}|1\rangle}{\sqrt2}.
\]

A Hadamard reads the parity in the ideal case; it does not return both predicate values separately. With visibility \(V\),

\[
P(B_{\mathrm{out}}=0)
=\frac12\left[1+
V(-1)^{p(h_0)\oplus p(h_1)}\right]
\]

after phase calibration.

### 6.4 Endpoint-equivalence theorem **[D; standard channel consequence]**

**Theorem 1.** Experiments that only prepare \(B\), apply the complete device, and measure \(B\) at the endpoint cannot certify a particular internal history, memory, or agent.

**Proof.** The full device induces a CPTP channel \(\mathcal E_B\) from the input of \(B\) to its output. Every accessible statistic is

\[
p(y|\rho,x)=\operatorname{Tr}[M_y^x\mathcal E_B(\rho)].
\]

An alternative box implementing \(\mathcal E_B\) directly, without the proposed internal realization, reproduces all such statistics. In the ideal unitary case of Lemma 2 it suffices to apply

\[
U_{\mathrm{eff}}
=\sum_b(-1)^{p(h_b)}|b\rangle\!\langle b|.
\]

Therefore the internal realization is not identifiable from endpoint data. \(\square\)

**Corollary.** More internal steps do not create a new observable unless access changes. Excluding the direct phase shortcut requires intermediate interventions, a black-box promise, locality or latency restrictions, a commitment mechanism, or a proved resource separation.

---

## 7. No-transcription results

### 7.1 Exact no-transcription theorem **[K/D]**

**Theorem 2.** In a purified two-alternative interferometer, a persistent external transcript that perfectly identifies the alternative is incompatible with nonzero visibility. Perfect visibility is incompatible with any distinguishable external alternative record.

**Proof.** A perfect transcript means orthogonal record supports. Including all transcript degrees in a purification gives
\(\langle\ell_1|\ell_0\rangle=0\), so Lemma 1 yields \(V=0\). Conversely, \(V=1\) implies the same leakage ray and therefore the same reduced state for every external subsystem; every external measurement has identical statistics. \(\square\)

The qualitative statement is established in complementarity and channel theory. Its role here is to fix the meaning of zero record and prevent the contradictory claim that a complete memory can persist while the same alternatives interfere perfectly.

### 7.2 Approximate form

If the residual system guesses the alternative with
\(P_{\mathrm{guess}}=\tfrac12(1+D_L)\), then

\[
V\le\sqrt{1-D_L^2}
=2\sqrt{P_{\mathrm{guess}}(1-P_{\mathrm{guess}})}.
\]

A guessing advantage
\(P_{\mathrm{guess}}=\tfrac12+\delta\) therefore implies

\[
V\le\sqrt{1-4\delta^2}.
\]

The bound is indifferent to whether the record is called memory, observation, work, or environment.

### 7.3 Authorized predicate versus a copied transcript

The authorized output \(Z=p(h)\) may be recorded classically after the relevant interference, or in a register that does not distinguish the alternatives being recombined. If \(p(h_0)\ne p(h_1)\) is copied into orthogonal states before recombination, the copy is path information and destroys their interference. A relative phase avoids a readable pre-recombination transcript.

For privacy, the benchmark must additionally contain fine histories sharing the same predicate. A two-path demonstration with \(p(h_b)=b\) can test phase survival and no-transcription between paths, but it cannot test privacy *beyond the predicate*, because that predicate already reveals the full path label.

---

## 8. Late-challenge protocol

The purpose of a late challenge is not to make an implementation magically identifiable. It defines a setting in which answers cannot be precomputed under explicit restrictions.

### 8.1 Timeline

Let \(t_c\) be a commitment point. For a classical challenge \(x\sim q(x)\):

1. **Preparation:** initialize \(\sum_b\alpha_b|b\rangle_B\), \(R\), and \(Q\).
2. **Prehistory:** before \(x\) is available,

   \[
   U_<|b\rangle|0\rangle_R
   =|b\rangle|r_b^<\rangle_R .
   \]

3. **Late challenge:** after \(t_c\), sample \(x\) independently and deliver it only through the declared interface.
4. **Continuation and response:** use pre-existing memory,

   \[
   U_{>,x}|r_b^<\rangle_R|0\rangle_Q
   =|r_{b,x}^>\rangle_R|a(b,x)\rangle_Q .
   \]

5. **Phase transfer:** apply \(Z_Q\), or equivalent phase kickback, to obtain \((-1)^{a(b,x)}\).
6. **Recovery:** apply \(U_{>,x}^\dagger\) and \(U_<^\dagger\). Retain \(x\) until this finishes; \(x\) is branch-independent.
7. **Recombination:** measure complementary bases of \(B\) and estimate relative phase for each context \(x\).

Ideally,

\[
|\Psi_f(x)\rangle
=\sum_b\alpha_b(-1)^{a(b,x)}
|b\rangle_B|0\rangle_R|x\rangle_X .
\]

For two alternatives the endpoint reveals
\(a(0,x)\oplus a(1,x)\).

### 8.2 Coherent challenge

Prepare

\[
|\chi\rangle_X=\sum_x\sqrt{q_x}|x\rangle_X
\]

and control \(U_{>,x}\) by \(X\). The ideal final state is

\[
\sum_{b,x}\alpha_b\sqrt{q_x}(-1)^{a(b,x)}
|b\rangle_B|x\rangle_X|0\rangle_R .
\]

Tracing \(X\) gives a \(B\)-coherence factor

\[
\Gamma=\sum_xq_x(-1)^{a(0,x)\oplus a(1,x)}.
\]

Measuring \(X\) computationally estimates conditioned answers; complementary measurements of \(X\) test challenge coherence. This variant increases control demands and is not cryptographic zero knowledge without a complete adversarial analysis.

### 8.3 A late challenge alone does not exclude direct phase

After receiving \(x\), a box may still apply

\[
U_{\mathrm{eff}}(x)
=\sum_b(-1)^{a(b,x)}|b\rangle\!\langle b|.
\]

The challenge excludes shortcuts only together with at least one explicit condition:

- **[A1: causal interface]** after \(t_c\), \(x\) reaches only a local port lacking all information needed to evaluate \(a\);
- **[A2: commitment]** a binding state or commitment created before \(x\) fixes part of the response;
- **[A3: black-box instance]** the prehistory queried an unknown physical instance not reproducible by the verifier;
- **[A4: resource limit]** direct computation of every valid response exceeds a proved depth, query, communication, or memory bound;
- **[A5: multi-time interventions]** an informationally sufficient intervention family separates the history class from shortcuts;
- **[A6: independence]** the challenge source is independent of the device before \(t_c\), or residual dependence is bounded.

Without one of these, “the challenge proves the history occurred” is false.

### 8.4 Two-layer certification

Intermediate measurements used for process-tensor reconstruction generally disturb the run whose coherence one wishes to recover. The first benchmark should therefore separate:

1. a **causal layer**, using many instrumented runs to estimate
   \(p(\mathbf a|\mathbf x)\), a witness \(W\), or a bound on \(d_{\min}\);
2. an **echo layer**, using matched runs without invasive readout—or with fully coherent dilations—to estimate phase, visibility, reset, and leakage.

Inferring that both layers characterize the same implementation requires randomized contexts, drift bounds, and stability tests. A single-run version requires coherent interventions and control of every dilation.

---

## 9. Definition of a Zero-Record Causal Proof

### 9.1 Protocol object

**Definition 1 [P].** Relative to a security parameter \(\lambda\), valid causal class \(\mathfrak C_1\), null class
\(\mathfrak C_0\), intervention family \(\mathfrak I\), noninjective predicate \(p\), and access model \(\mathfrak A\), a ZRCP is a family

\[
\Pi_\lambda=(\mathrm{Prep},\mathrm{Commit},\mathrm{Challenge},
\mathrm{Respond},\mathrm{Phase},\mathrm{Recover},\mathrm{Read})
\]

that produces a decision \(Z\), a final verifier state \(T\), and recovery metrics. The specification must state:

- trusted devices and calibrations;
- auxiliary inputs and possible reference entanglement;
- messages and accessible ports at every time;
- the authorized output and the noninjective fine-history coarse-graining;
- the shortcut class included in \(\mathfrak C_0\);
- whether guarantees are information-theoretic or computational.

### 9.2 Completeness

For every valid \(\Theta\in\mathfrak C_1\) satisfying the statement,

\[
\Pr_\Pi[Z=1|\Theta]\ge c,
\qquad
\epsilon_\diamond(\Theta)\le\epsilon_{\mathrm{rec}},
\qquad
\epsilon_{\mathrm{priv}}(\Theta)
\le\epsilon_{\mathrm{priv}}^{\max}.
\]

\(c\) is completeness. A noisy experiment must report confidence intervals rather than only an ideal probability.

### 9.3 Causal soundness

For every \(\Theta_0\in\mathfrak C_0\) compatible with the access model,

\[
\Pr_\Pi[Z=1|\Theta_0]\le s.
\]

\(s\) is soundness. The gap \(c-s>0\) has meaning only relative to
\(\mathfrak C_0\). There is no absolute soundness against an unrestricted box: Theorem 1 allows it to implement the effective channel.

### 9.4 Physical privacy

A simulator given only the public challenge and authorized coarse-grained output must reproduce the residual verifier channel:

\[
\frac12\|\mathcal T_\Theta-
\mathcal S\circ\mathcal P_{p,\Theta}\|_\diamond
\le\epsilon_{\mathrm{priv}}.
\]

The tested support must satisfy the nonvacuity requirement. The report must provide fiber sizes and at least one worst-case within-fiber distinguishability. If the guarantee covers only an honest verifier, call it *honest physical transcript privacy*, not zero knowledge.

### 9.5 Recoverability

At minimum,

\[
F_{\mathrm{reset}}^{\mathrm{state}}\ge1-\epsilon_R,
\qquad
V_{\mathrm{excess}}\ge v,
\qquad
\epsilon_{\mathrm{ZR}}\le\epsilon_Z,
\]

where

\[
V_{\mathrm{excess}}=\frac{V_\Pi}{V_{\mathrm{identity\ echo}}}
\]

uses the same schedule, compilation, and two-qubit-gate pattern. Whenever feasible, require the channel criterion
\(\epsilon_\diamond\le\epsilon_{\mathrm{rec}}\).

### 9.6 Causal authenticity

Acceptance must depend on a multi-time behavior violating a null-class bound,

\[
W[p_{\mathrm{obs}}]>
\beta_0+\Delta_{\mathrm{stat}}+\Delta_{\mathrm{sys}} .
\]

\(\Delta_{\mathrm{stat}}\) is statistical uncertainty and
\(\Delta_{\mathrm{sys}}\) a declared systematic-error budget. An endpoint phase test alone does not meet this clause.

### 9.7 Relation to cryptographic proofs

A ZRCP borrows the vocabulary of completeness, soundness, and simulation, but it does not inherit quantum zero-knowledge theorems. A cryptographic version must additionally define languages, witnesses, malicious polynomial-time strategies, auxiliary inputs, negligible error, composition, simulator efficiency, and any computational assumptions.

---

## 10. Conditional certification results

### 10.1 Witness separation **[D/K]**

**Proposition 3.** Let \(\mathcal P_0\) be the convex closed set of null-class behaviors under \(\mathfrak I\). If
\(p_*\notin\mathcal P_0\), then a linear witness \(W\) and
\(\delta>0\) exist such that

\[
W[p_*]\ge\sup_{q\in\mathcal P_0}W[q]+\delta .
\]

**Proof.** Apply strict separation of a point from a closed convex set in finite dimension. \(\square\)

This does not guarantee an efficient or device-independent witness, nor does it show that \(p_*\) has a zero-record realization.

### 10.2 Combining causal and echo layers **[D]**

**Proposition 4.** Suppose:

1. causal runs estimate \(W\) with total error below \(\delta/2\);
2. the null class obeys \(W\le\beta_0\), while the implementation obeys \(W\ge\beta_0+\delta\);
3. echo runs satisfy \(\epsilon_\diamond\le\varepsilon\);
4. a stability test bounds operational distance between the two run families by \(\eta\).

Then the experiment certifies, at the declared confidence level, (i) causal structure outside \(\mathfrak C_0\) in the causal layer and (ii) recovery within \(\varepsilon+\eta\) for an implementation compatible with the echo layer.

**Proof.** The witness gap and error prove (i). The triangle inequality between causal-layer, echo-layer, and ideal process gives (ii). \(\square\)

The stability premise is substantive. The proposition does not imply that one individual run was both invasively measured and recohered.

### 10.3 No unrestricted soundness **[D]**

**Theorem 3.** No ZRCP has a positive soundness gap against the class of all physically allowed combs with unrestricted access to every input, challenge, and output, if the class is not required to realize the specified internal history.

**Proof.** Compose the honest operations, conditional on every verifier message, into an effective quantum comb. An unrestricted dishonest realization can implement the same comb through a different Stinespring realization without the designated internal semantics. Every verifier-accessible behavior is identical, so the dishonest acceptance probability equals the honest one and \(s\ge c\). \(\square\)

Nontrivial soundness therefore requires physical interfaces, locality, resource limits, commitment, or computational assumptions.

---

## 11. Excluding the direct-phase shortcut

The direct-phase model is the central null hypothesis, not an optional control.

### 11.1 Minimal shortcut class

For a challenge \(x\), let \(\mathfrak C_{\mathrm{phase}}\) contain circuits that receive \((b,x)\), or coherent equivalents, and apply
\(e^{i\phi f(b,x)}\) without using the internal causal edge being certified. With unrestricted resources this class contains the honest effective channel and cannot be excluded.

### 11.2 Legitimate routes to a separation

1. **Oracle query:** the history accesses an unknown physical black box. Compare coherent and independent-sample query complexities.
2. **Space-time locality:** required information is distributed, and a response deadline prevents a local shortcut from collecting it under a precise architecture or relativistic model.
3. **Memory dimension:** temporal correlations violate a bound for realizations with \(d_M<d\).
4. **Cryptographic commitment:** a device is bound before receiving \(x\); soundness then rests on explicit cryptographic assumptions.
5. **Circuit complexity:** prove that direct synthesis of every
   \(U_{\mathrm{eff}}(x)\) in the allowed gate set costs more than history execution plus recovery.
6. **Multi-time tomography or witnessing:** characterize enough slots of a comb to rule out the null causal factorization.

### 11.3 Matched controls

Every experimental family should include:

- a directly compiled phase with matched duration and entangling-gate pattern;
- the best known optimized direct phase, even if shallower;
- history with memory replaced by a fresh state at each step;
- the action-to-future-world edge cut;
- challenge delivered before versus after commitment;
- deliberately challenge-correlated runs to calibrate independence sensitivity;
- retained memory, retained garbage, and unreversed environment;
- an identity echo with the same native gates;
- unrecorded random \(Z\) dephasing of \(B\) to form an incoherent mixture while preserving populations;
- both \(X_B\) and \(Y_B\) readout, not only a computational probability.

### 11.4 A possible null-class witness

For a concrete benchmark with inputs \(\mathbf x\) and outputs
\(\mathbf a\), define a score

\[
W=\sum_{\mathbf a,\mathbf x}
c_{\mathbf a,\mathbf x}p(\mathbf a|\mathbf x).
\]

The project must calculate

\[
\beta(d,\mathcal G,\mathcal L)
=\sup_{\Theta\in\mathfrak C_0(d,\mathcal G,\mathcal L)}W[\Theta],
\]

where \(d\) bounds memory, \(\mathcal G\) restricts gates or queries, and \(\mathcal L\) specifies causal communication. Omitting these arguments would make the bound appear more device-independent than it is.

An SDP hierarchy may bound \(\beta\) when the process dimension and intervention operators are controlled. A fully device-independent claim needs a causal inequality or self-testing statement and substantially stronger analysis.

---

## 12. Noise and elementary scaling laws

### 12.1 Independent leakage per step **[D under stated model]**

Assume step \(j\) creates branch-conditioned pure leakage with overlap

\[
\gamma_j=\langle e_{1,j}|e_{0,j}\rangle
\]

and different steps leak into independent factors. Then

\[
V_T=V_0\prod_{j=1}^{T}|\gamma_j|.
\]

If \(|\gamma_j|=1-\ell_j\), with \(\ell_j\ll1\),

\[
-\log(V_T/V_0)
=-\sum_j\log(1-\ell_j)
\simeq\sum_j\ell_j.
\]

Thus log-visibility loss is approximately additive under independent weak leakage. Correlated environments invalidate the product and should be represented by a process tensor.

### 12.2 Effective dephasing **[K]**

The channel

\[
\mathcal D_\lambda(\rho)
=(1-\lambda)\rho+\lambda Z\rho Z
\]

multiplies coherence by \(1-2\lambda\). After \(n\) independent applications,

\[
V_n=V_0|1-2\lambda|^n.
\]

This is a baseline control, not an observer theory.

### 12.3 Coherent inversion error **[D]**

Let the implemented inverse be

\[
\widetilde U_H^\dagger=U_H^\dagger e^{-i\epsilon K}.
\]

For internal state \(|h_b\rangle\), the echo amplitude is

\[
A_b=\langle h_b|e^{-i\epsilon K}|h_b\rangle
=1-i\epsilon\langle K\rangle_b
-\frac{\epsilon^2}{2}\langle K^2\rangle_b+O(\epsilon^3),
\]

and

\[
|A_b|^2=
1-\epsilon^2\operatorname{Var}_b(K)+O(\epsilon^3).
\]

Coherent errors can shift phase at first order while reducing return probability only at second order. Both interferometric quadratures are therefore required.

### 12.4 Partial record leakage model

Introduce an environment qubit with

\[
|e_0\rangle=|0\rangle,\qquad
|e_1\rangle=\cos\vartheta|0\rangle+\sin\vartheta|1\rangle .
\]

Then

\[
V=|\cos\vartheta|,
\qquad
D=|\sin\vartheta|,
\qquad
D^2+V^2=1.
\]

This exact one-parameter model is the canonical simulation for memory retained continuously between no record and a perfect transcript.

### 12.5 Complexity is not a fundamental penalty without restrictions

For every exactly known finite unitary, \(U^\dagger U=I\), independently of how many steps one calls observations. Standard ideal quantum theory has no special loss term for “agent complexity.” Complexity matters through:

- the number of error locations;
- circuit depth relative to coherence time;
- leakage into uncontrolled degrees;
- inverse-synthesis difficulty;
- thermodynamic or control restrictions;
- redundancy of records;
- characterization and certification complexity.

A universal law \(V=V(\text{cognitive complexity})\) has no basis without additional physics.

---

## 13. Candidate resource quantities

These definitions are proposals and must be checked against prior literature before being named as new measures.

### 13.1 Recoherable causal complexity **[P]**

For target visibility \(v\), privacy error \(\epsilon\), and allowed controls
\(\mathfrak U\), define

\[
\mathsf{RCC}_{\mathfrak U}(v,\epsilon)
=\sup_{\Theta}\left\{
\log d_{\min}(\Theta;\mathfrak I):
V_{\mathrm{excess}}(\Theta)\ge v,
\epsilon_{\mathrm{priv}}(\Theta)\le\epsilon
\right\}.
\]

It asks for the largest certifiable causal memory that remains recoverable in a resource class. Without dimension, depth, or control restrictions the supremum can be trivial or unbounded.

### 13.2 Zero-record proof cost **[P]**

\[
\mathsf C_{\mathrm{ZR}}(p,\varepsilon)
=\min_{\Pi}
\left\{
(N_q,D,N_a,N_{\mathrm{shots}}):
\Pi\text{ proves }p,
\epsilon_\diamond,\epsilon_{\mathrm{priv}}\le\varepsilon
\right\}.
\]

The tuple contains logical qubits, depth, ancillas, and sample count. A Pareto frontier is preferable to arbitrary scalar weights.

### 13.3 Coherent-access advantage **[P]**

For a task \(\mathcal T\), compare:

- \(Q_{\mathrm{coh}}(\mathcal T)\): calls to the same process instance with coherent control across time and an inverse interface;
- \(Q_{\mathrm{samp}}(\mathcal T)\): independent preparations, classical intermediate outcomes, and restarts.

A proved separation in \(Q_{\mathrm{samp}}/Q_{\mathrm{coh}}\) is a stronger criterion than simply comparing a superposition with a dephased mixture.

### 13.4 Nonvacuous privacy capacity **[P]**

For each output \(z\), define the effective tested fiber size

\[
n_z^{\mathrm{eff}}=
\exp H_{\min}(H|Z=z),
\]

or, in a purely combinatorial benchmark, use
\(|\operatorname{supp}(H|Z=z)|\). Report

\[
\mathsf P_{\mathrm{NV}}
=\min_{z:q_z>0}\log n_z^{\mathrm{eff}}
\]

alongside leakage. A privacy error is not scientifically impressive if
\(\mathsf P_{\mathrm{NV}}=0\), because no hidden fine-history choice existed after conditioning on the authorized predicate.

---

## 14. Research conjectures

### Conjecture A: zero-record causal separation **[C]**

There exists an explicit family of oracle processes
\(\{\Theta_n\}\) and noninjective predicates \(p_n\) for which a coherent ZRCP uses \(\operatorname{poly}(n)\) calls and leaves
\(\epsilon_{\mathrm{priv}}=\operatorname{negl}(n)\), whereas every independent-sample protocol with the same completeness and soundness needs \(2^{\Omega(n)}\) calls.

**Status.** Plausible by analogy with quantum algorithmic measurement, but unproved under the zero-record condition. The first task is to attempt a reduction to existing QUALM separations; a direct reduction would be useful but not a new separation.

### Conjecture B: local recovery cost of causal memory **[C]**

Under fixed local geometry, bounded inverse depth, and nonzero local noise, jointly certifying memory dimension \(d_M\) and recovering visibility \(v\) requires cost growing with a nonconstant function of
\(\log d_M\), causal depth, and record redundancy.

**Caveat.** Without locality, noise, or a depth bound, the fundamental version is false: apply the exact global inverse.

### Conjecture C: late commitment versus precompiled phase **[C]**

There exists an interactive game with a post-commitment challenge and explicit communication/locality restrictions in which memory-bearing processes achieve value \(c\), null processes without the designated memory edge satisfy \(s<c\), and the honest game admits a coherent recoverable implementation.

**Status.** Unproved. The game, both values, and signaling loopholes must be established before using the phrase “certified genuine history.”

### Conjecture D: composable physical transcript privacy **[C]**

Diamond-norm closeness of the complementary transcript channel to a channel depending only on a noninjective \(p(h)\), together with approximate logical recovery, yields a useful composable zero-record notion under sequential composition.

**Status.** This likely overlaps continuity of Stinespring dilations, approximate quantum error correction, and simulation-based security. A literature reduction must precede any novelty claim.

### Conjecture E: coarse-graining versus recoverability **[C]**

For a fixed physical process and allowed recovery class, coarser authorized predicates admit weakly better optimal zero-record recovery than finer predicates:

\[
p_1=g\circ p_2
\quad\Longrightarrow\quad
\epsilon_{\mathrm{rec}}^*(p_1)
\le \epsilon_{\mathrm{rec}}^*(p_2)
\]

under a carefully specified common task. The intuitive reason is that a coarser output permits more histories to be erased together. The claim is not yet rigorous because task-dependent phase encodings and verifier power can reverse naive orderings.

---

## 15. Fundamental and practical limits

### 15.1 Established boundaries

1. **No realization identification from endpoints:** Theorem 1.
2. **No perfect transcript with perfect interference:** Theorem 2.
3. **No unlimited readout:** finite quantum systems and finite queries yield bounded accessible classical information; extracting many predicates requires new preparations or calls.
4. **Interpretive neutrality:** interpretations sharing unitary dynamics and the Born rule predict the same channel statistics here.
5. **Inaccessible environment:** coherence recovery requires controlling degrees that retain distinguishable information, or implementing an equivalent correction channel. There is no generic operational access to arbitrarily dispersed natural records.

### 15.2 Limits of the proposed formalism

- \(d_{\min}\) is interface-relative; memory can sometimes be reassigned between accessible system and environment.
- Full process-tensor tomography scales poorly. Witnesses, low-rank assumptions, or tensor-network structure will be needed.
- An imperfectly independent challenge source weakens soundness.
- Explicit inversion approximately doubles logical depth and can magnify coherent calibration errors.
- Absolute zero record cannot audit every degree of freedom in the universe; experiments bound it within a physical environment model.
- A causal witness in one batch and recoherence in another need stability assumptions.
- This draft provides neither a computational advantage theorem, a complete cryptographic protocol, nor fully device-independent certification.

### 15.3 What could count as evidence of nonunitarity

Unexplained visibility loss is insufficient. A claim would require:

1. a calibrated environmental model with uncertainties;
2. a preregistered prediction from a parameterized nonunitary model;
3. comparison against alternative noise hypotheses;
4. multiplatform replication;
5. blind analysis and selection controls;
6. predeclared significance or Bayesian model-comparison criteria.

The first diagnosis of an anomaly must be an incomplete device model.

---

## 16. Kill criteria

The claim of a new research primitive must be abandoned or reduced to benchmarking/pedagogy if any of the following holds:

1. under its own access model, the ZRCP reduces to an effective diagonal gate with no resource gap;
2. no noncircular, physically motivated null class can be stated;
3. the causal witness is only a renamed standard non-Markovianity, fidelity, or OTOC measure with no new task;
4. privacy merely means the experimenter chose not to read an ancilla that still stores the history;
5. the authorized output is copied before interference and destroys the coherence claimed to be recovered;
6. reset fidelity is reported only on visible ancillas without a complementary-basis or leakage test;
7. causal structure and recoherence are measured in different circuits without a stability argument;
8. the late challenge is available before commitment, correlated with the device, or exposed through a side channel;
9. an optimized direct phase matches every claimed resource and no black-box task remains;
10. the advantage disappears against standard coherent oracle access;
11. the tested predicate is injective on history support, making privacy beyond the predicate vacuous;
12. prior literature already contains the same definition, game, and bound;
13. interpretive language is the only remaining distinction.

---

## 17. Required theorem and simulation program

### 17.1 Symbolic tests

- verify Lemma 2 for \(K=2\) and \(K=4\);
- display the effective-phase equivalence explicitly;
- construct Proposition 1's counterexample;
- sample residual states and verify the appropriate
  \(D^2+V^2\le1\) relation;
- distinguish visible reset from complementary-channel decoupling;
- compare absolute zero record with within-fiber privacy;
- include at least four fine histories with a two-valued noninjective predicate.

### 17.2 Circuit ladder

1. **P0 — phase baseline:** two alternatives, Boolean predicate, phase kickback, and direct-phase control.
2. **P1 — partial record:** one memory qubit retained, rotated, or erased; map the \(D\)-\(V\) curve.
3. **P2 — reversible automaton:** at least three steps, persistent memory, and an action changing a later world state.
4. **P3 — noninjective privacy:** at least four fine histories, parity or threshold output, and within-fiber transcript tests.
5. **P4 — classical late challenge:** randomized challenge contexts and causal witness against a declared memoryless model.
6. **P5 — coherent challenge:** joint \(BX\) phase readout.
7. **P6 — dimension restriction:** compare behaviors against processes with bounded \(d_M\).

### 17.3 Noise suite

Model separately:

- dephasing of \(B\);
- one- and two-qubit depolarization;
- amplitude damping;
- coherent overrotation in the inverse;
- branch-dependent crosstalk;
- adjustable leakage to an unreversed ancilla;
- drift between causal and echo batches;
- readout error with documented mitigation.

Report:

\[
V,\ \varphi,\ F_{\mathrm{reset}},\
\epsilon_{\mathrm{fib}},\ \epsilon_{\mathrm{ZR}},\
F_p,\ W,\ N_{\mathrm{shots}},
\]

with confidence intervals and a matched identity-echo normalization.

### 17.4 First-prototype success condition

The first prototype need not claim quantum advantage. It must show jointly:

1. a multi-time statistic incompatible with an explicit null class;
2. causal dependence of an internal action on a later observation;
3. transfer of a noninjective predicate to phase;
4. recovery with visibility significantly above retained-memory and retained-garbage controls;
5. residual within-fiber indistinguishability within an experimental budget;
6. endpoint equivalence with a direct phase, reported openly as a control.

---

## 18. Claim ledger

| Claim | Status | Authorized conclusion |
|---|---|---|
| Phase kickback survives uncomputation | [K] | A reversible predicate can remain as phase. |
| Visibility equals pure leakage overlap | [K] | Distinguishable traces suppress coherence. |
| Perfect transcript conflicts with perfect interference | [K/D] | Exact persistent path information forbids recoherence. |
| Endpoints do not certify internal history | [D; standard channel theory] | Additional access or assumptions are necessary. |
| Visible ancilla reset does not imply recoherence | [D] | Test the control coherence and complementary leakage. |
| Conditional privacy is vacuous for injective predicates | [D; elementary information theory] | Use noninjective coarse-graining and multiple fine histories per fiber. |
| ZRCP definition | [P] | Research framework; novelty not established. |
| Late challenge plus causal witness | [P] | Candidate architecture; assumptions remain essential. |
| Exponential coherent-versus-sample separation | [C] | Research target, not a result. |
| Universal complexity–recoherence law | Not established | False without physical resource restrictions. |
| Test of a unique interpretation | Unauthorized | Unitary data remain interpretively underdetermined. |

---

## 19. Minimal reference backbone

These primary sources support components of the framework; none by itself establishes the proposed combination as novel.

1. C. H. Bennett, “Logical Reversibility of Computation,” *IBM Journal of Research and Development* 17, 525–532 (1973). <https://doi.org/10.1147/rd.176.0525>
2. B.-G. Englert, “Fringe Visibility and Which-Way Information: An Inequality,” *Physical Review Letters* 77, 2154 (1996). <https://doi.org/10.1103/PhysRevLett.77.2154>
3. D. Kretschmann, D. Schlingemann, and R. F. Werner, “The Information-Disturbance Tradeoff and the Continuity of Stinespring's Representation,” *IEEE Transactions on Information Theory* 54, 1708–1717 (2008). <https://doi.org/10.1109/TIT.2008.917696>
4. C. Bény and O. Oreshkov, “Approximate simulation of quantum channels,” *Physical Review A* 84, 022333 (2011). <https://doi.org/10.1103/PhysRevA.84.022333>
5. G. Chiribella, G. M. D'Ariano, and P. Perinotti, “Theoretical framework for quantum networks,” *Physical Review A* 80, 022339 (2009). <https://doi.org/10.1103/PhysRevA.80.022339>
6. F. A. Pollock et al., “Non-Markovian quantum processes: Complete framework and efficient characterization,” *Physical Review A* 97, 012127 (2018). <https://doi.org/10.1103/PhysRevA.97.012127>
7. C. Giarmatzi and F. Costa, “A quantum causal discovery algorithm,” *npj Quantum Information* 4, 17 (2018). <https://doi.org/10.1038/s41534-018-0062-6>
8. C. Giarmatzi and F. Costa, “Witnessing quantum memory in non-Markovian processes,” *Quantum* 5, 440 (2021). <https://doi.org/10.22331/q-2021-04-26-440>
9. L. B. Vieira, S. Milz, G. Vitagliano, and C. Budroni, “Witnessing environment dimension through temporal correlations,” *Quantum* 8, 1224 (2024). <https://doi.org/10.22331/q-2024-01-10-1224>
10. A. Aharonov, J. Cotler, and X.-L. Qi, “Quantum algorithmic measurement,” *Nature Communications* 13, 887 (2022). <https://doi.org/10.1038/s41467-022-28401-2>
11. J. Watrous, “Zero-knowledge against quantum attacks,” *SIAM Journal on Computing* 39, 25–58 (2009). <https://doi.org/10.1137/060670997>
12. A. Broadbent, Z. Ji, F. Song, and J. Watrous, “Zero-knowledge proof systems for QMA,” *SIAM Journal on Computing* 49, 245–283 (2020). <https://doi.org/10.1137/16M1101220>
13. T. Vidick and T. Zhang, “Classical zero-knowledge arguments for quantum computations,” *Quantum* 4, 266 (2020). <https://doi.org/10.22331/q-2020-05-14-266>
14. A. Coladangelo et al., “MPC in the Quantum Head (or: Superposition-Secure (Quantum) Zero-Knowledge),” *Quantum* 10, 2161 (2026). <https://doi.org/10.22331/q-2026-07-15-2161>
15. L. S. V. Santos et al., “Device-independent quantum memory certification in two-point measurement experiments,” arXiv:2601.14191 (2026), preprint at the cutoff date. <https://arxiv.org/abs/2601.14191>
16. C. Elouard et al., “Quantum erasing the memory of Wigner's friend,” *Quantum* 5, 498 (2021). <https://doi.org/10.22331/q-2021-07-08-498>
17. W. H. Zurek, “Quantum reversibility is relative, or does a quantum measurement reset initial conditions?” *Philosophical Transactions of the Royal Society A* 376, 20170315 (2018). <https://doi.org/10.1098/rsta.2017.0315>

---

## 20. Canonical statement

> Causal Certification of Reversible Quantum Histories studies multi-time processes with internal memory whose causal structure can be distinguished, under an explicit intervention and access model, from processes lacking that memory. A Zero-Record Causal Proof transfers only an authorized noninjective coarse-graining to an interferometric output, recovers the logical channel, and leaves all remaining registers approximately simulable from that coarse-grained output alone. An endpoint measurement never certifies the internal realization; authenticity must come from interventions, commitments, physical interfaces, or resource separations, while absence of transcript is quantified by decoupling and complementary-channel information.

This is a research agenda, not yet an advantage theorem or a priority claim. Its first mandatory result is to locate precisely where standard phase kickback ends and a genuinely distinct causal-certification task begins.

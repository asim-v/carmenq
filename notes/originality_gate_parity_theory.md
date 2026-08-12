# Streaming parity audit--return theorem (scratch originality gate)

Status: exact mathematical result under an explicit **classical-memory comb** null model. This note is for adversarial novelty review only. It is not yet a manuscript claim. An adversarial search through 2026-08-12 found no exact primary-source collision, but this establishes a cautious "to our knowledge" claim rather than proof of priority.

## 1. Operational model

Fix `n >= 1`. For every slot `i`, the verifier prepares

\[
|\Phi\rangle_{R_iA_i}=\frac{|00\rangle+|11\rangle}{\sqrt 2}
\]

and sends `A_i` to the tested device. The device must return an output system `B_i` before `A_{i+1}` is supplied; returned systems are sequestered. After the last return, a hidden coin selects one of two tests.

* **AUDIT.** The verifier measures every `R_i` in the computational basis, obtaining independent uniform bits `X_i`, and asks the device's committed classical transcript to predict
  \(Z=X_1\oplus\cdots\oplus X_n\). The returned `B_i` are unavailable in this branch.
* **RETURN.** The verifier uses the committed classical transcript as a control for an otherwise arbitrary joint recovery channel on `B_1\cdots B_n`, and accepts with the projector onto \(|\Phi\rangle^{\otimes n}\) (and, if desired, checks that the *visible* controller is placed in a standard reset state).

Write the two acceptance probabilities as `P_A` and `F_R`. For a weight `0 <= lambda <= 1`, the support score is

\[
S_{n,\lambda}=\lambda P_A+(1-\lambda)F_R.
\]

### Classical-memory comb null

The null may have unlimited classical memory, unlimited shared classical randomness, adaptive control, arbitrary finite outcome alphabets, and arbitrary disposable quantum ancillas within a slot. It may implement a completely general (not necessarily QND) instrument `A_i -> B_i` in each slot. It may not carry a coherent quantum system from one slot to the next, and it may not begin with entanglement shared among different slots. Conditional on a fully refined classical history, its leaf channel therefore factorizes over slots. Previous outputs cannot be touched after commitment.

Shared randomness can be included in the refined transcript. Revealing that refinement can only help both the optimal parity decoder and the transcript-conditioned recovery, so this gives a valid upper bound even when some randomness was initially hidden.

This is a resource-restricted theorem, not a device-independent theorem. Operationally certifying that an implementation has no hidden coherent memory is a separate problem.

## 2. The exact theorem

Define

\[
f(t)=\frac{1+\sqrt{1-t^2}}2,\qquad 0\leq t\leq1.
\]

**Theorem (exact streaming parity frontier).** Over all classical-memory combs above, including arbitrary non-QND local instruments and arbitrary transcript-conditioned joint recovery,

\[
\boxed{
\sup S_{n,\lambda}
=
\max_{0\leq t_1,\ldots,t_n\leq1}
\left[
\frac{\lambda}{2}\left(1+\prod_{i=1}^n t_i\right)
+(1-\lambda)\prod_{i=1}^n f(t_i)
\right].
}
\tag{1}
\]

The maximum is attained by independent binary-symmetric diagonal instruments. Moreover, an optimizer of (1) can always be chosen with equal strengths,

\[
t_1=\cdots=t_n=t.
\tag{2}
\]

Thus an equivalent one-variable expression is

\[
\boxed{
\sup S_{n,\lambda}
=\max_{0\leq t\leq1}
\left[
\frac{\lambda}{2}(1+t^n)+(1-\lambda)f(t)^n
\right].
}
\tag{3}
\]

The reduction is exact: adaptivity, outcome-dependent measurement strengths, non-QND dynamics, extra within-slot ancillas, and a joint rather than product recovery do not improve the score.

## 3. Proof against arbitrary local instruments

### 3.1 A recovery lemma that removes the non-QND freedom

Let `N_k : A -> B` be any trace-nonincreasing CP map associated with a complete classical transcript `k`, where `dim(A)=D`. Define its basis likelihoods

\[
p(k|x)=\operatorname{Tr}N_k(|x\rangle\!\langle x|),
\qquad x\in\{0,1\}^n.
\]

For any recovery channel `R_k : B -> A`, its contribution to maximally-entangled-state fidelity obeys

\[
F_k
:=
\langle\Phi_D|({\rm id}\otimes R_k\circ N_k)
(\Phi_D)|\Phi_D\rangle
\leq
\frac1{D^2}\left(\sum_x\sqrt{p(k|x)}\right)^2.
\tag{4}
\]

To prove (4), put `Lambda_k=R_k o N_k` and expand

\[
D^2F_k
=\sum_{x,y}\langle x|\Lambda_k(|x\rangle\!\langle y|)|y\rangle.
\]

Positivity of the Choi matrix gives, for every `x,y`,

\[
\left|\langle x|\Lambda_k(|x\rangle\!\langle y|)|y\rangle\right|
\leq
\sqrt{
\langle x|\Lambda_k(|x\rangle\!\langle x|)|x\rangle
\langle y|\Lambda_k(|y\rangle\!\langle y|)|y\rangle}
\leq\sqrt{p(k|x)p(k|y)}.
\]

Summing proves the claim. It already allows `R_k` to be an arbitrary joint channel on all returned systems.

The bound is tight for every likelihood table: take the positive diagonal QND Kraus operator

\[
K_k=\sum_x\sqrt{p(k|x)}|x\rangle\!\langle x|
\]

and use the identity recovery. Consequently, once the classical likelihoods are fixed, arbitrary non-QND output dynamics can only decrease RETURN performance.

### 3.2 Causality factorizes each transcript likelihood

For a refined leaf transcript `k=(k_1,...,k_n)`, absence of inter-slot quantum memory implies

\[
p(k|x)=\prod_{i=1}^n p_i(k_i|x_i,k_{<i}).
\tag{5}
\]

At a prefix `h=k_{<i}`, set

\[
r_i(k_i|h)=\frac{p_i(k_i|0,h)+p_i(k_i|1,h)}2
\]

and, whenever `r_i>0`, define a signed posterior bias by

\[
s_i\delta_i
=\frac{p_i(k_i|0,h)-p_i(k_i|1,h)}
{p_i(k_i|0,h)+p_i(k_i|1,h)},
\qquad s_i\in\{-1,+1\},\quad0\leq\delta_i\leq1.
\tag{6}
\]

Because every fresh bit is uniform and independent of the earlier transcript, `r_i` is exactly the conditional probability of observing `k_i` given `h`. Along a realized transcript, Bayes' rule and (5) make the posterior over `x_1,...,x_n` a product distribution. Its parity bias is therefore

\[
\left|\mathbb E[(-1)^{X_1+\cdots+X_n}|k]\right|
=\prod_i\delta_i.
\]

Hence the optimal transcript-only parity decoder has

\[
P_A=\frac12\left(1+\mathbb E_k\prod_i\delta_i\right).
\tag{7}
\]

Applying (4), using (5), and summing the input strings separately at every slot gives

\[
\begin{aligned}
F_R
&\leq4^{-n}\sum_k
\left(\sum_x\sqrt{p(k|x)}\right)^2\\
&=\sum_k\prod_i
\frac{(\sqrt{p_i(k_i|0,h)}+\sqrt{p_i(k_i|1,h)})^2}{4}\\
&=\mathbb E_k\prod_i
\frac{1+\sqrt{1-\delta_i^2}}2
=\mathbb E_k\prod_i f(\delta_i).
\end{aligned}
\tag{8}
\]

Equations (7)--(8) are the key causal tensorization. They are valid even when the strengths `delta_i` are random and chosen adaptively from the entire earlier transcript.

### 3.3 Adaptivity cannot improve a linear support score

Combining (7)--(8),

\[
S_{n,\lambda}
\leq\frac\lambda2+
\mathbb E_k\left[
\frac\lambda2\prod_i\delta_i
+(1-\lambda)\prod_i f(\delta_i)
\right].
\]

Every realized path supplies one point `(delta_1,...,delta_n)` in `[0,1]^n`. An expectation of the bracket cannot exceed its largest pointwise value. This proves the upper bound in (1), without assuming symmetric outcomes or nonadaptive control.

For attainability, at slot `i` use outcomes `k_i in {0,1}` and Kraus operators

\[
K^{(i)}_0=
\begin{pmatrix}
\sqrt{(1+t_i)/2}&0\\0&\sqrt{(1-t_i)/2}
\end{pmatrix},
\qquad
K^{(i)}_1=
\begin{pmatrix}
\sqrt{(1-t_i)/2}&0\\0&\sqrt{(1+t_i)/2}
\end{pmatrix}.
\tag{9}
\]

The AUDIT guess is `k_1 xor ... xor k_n`. It has bias `product_i t_i`. Identity recovery attains `product_i f(t_i)`. Thus every value on the right side of (1), and in particular its maximum, is physically attained by an allowed null strategy.

## 4. Why equal strengths are globally optimal

Set

\[
t_i=\frac{2y_i}{1+y_i^2},\qquad
f(t_i)=\frac1{1+y_i^2},\qquad0\leq y_i\leq1.
\tag{10}
\]

Apart from the additive `lambda/2`, the objective becomes

\[
G(y_1,...,y_n)=
\frac{(1-\lambda)+\lambda2^{n-1}\prod_i y_i}
{\prod_i(1+y_i^2)}.
\tag{11}
\]

For fixed nonzero `Y=product_i y_i`, write `z_i=log y_i`. The function

\[
\phi(z)=\log(1+e^{2z})
\]

is strictly convex. Jensen's inequality says that, at fixed `sum_i z_i=log Y`, the denominator in (11) is minimized exactly when all `y_i` are equal. If `Y=0`, the numerator is `1-lambda` and the denominator is minimized by setting every `y_i=0`. This proves (2)--(3), including all boundary cases.

This is stronger than merely proving that a symmetric stationary point exists: no asymmetric, sparse, time-shared, or adaptive strength profile can win a support optimization.

## 5. Closed-form optimization and phase transition

Put `b=1-lambda`, `c=lambda 2^{n-1}` and let the common parameter in (10) be `y`. The nonconstant part of the support is

\[
g_n(y)=\frac{b+c y^n}{(1+y^2)^n}.
\tag{12}
\]

### One slot

For `n=1`, this is the usual binary information--recovery curve:

\[
\sup S_{1,\lambda}
=\frac12+\frac12\sqrt{\lambda^2+(1-\lambda)^2},
\qquad
t_*=\frac{\lambda}{\sqrt{\lambda^2+(1-\lambda)^2}}.
\tag{13}
\]

### Two slots

For `n=2`, the transition is continuous at `lambda_c=1/2`:

\[
\sup S_{2,\lambda}=
\begin{cases}
1-\lambda/2,&0\leq\lambda\leq1/2,\\[2mm]
\lambda/2+\dfrac{\lambda^2}{3\lambda-1},&1/2<\lambda\leq1.
\end{cases}
\tag{14}
\]

Above threshold,

\[
y_*^2=2-\frac1\lambda,
\qquad
t_*=\frac{2\sqrt{\lambda(2\lambda-1)}}{3\lambda-1}.
\]

### Three or more slots: a first-order all-or-nothing transition

For `n >= 3`, differentiation gives

\[
g_n'(y)=
\frac{ny}{(1+y^2)^{n+1}}
\left[c y^{n-2}(1-y^2)-2b\right].
\tag{15}
\]

The no-measurement endpoint `y=0` is always a local maximum. When a second local maximum exists, it is the larger solution of

\[
\lambda2^{n-2}y^{n-2}(1-y^2)=1-\lambda.
\tag{16}
\]

The global switch occurs at a nonzero `z_n=y_c^2` uniquely determined by

\[
\boxed{(1-z_n)(1+z_n)^{n-1}=1,\qquad
z_n\in\left(\frac{n-2}{n},1\right).}
\tag{17}
\]

The critical support weight is

\[
\boxed{
\lambda_c(n)=
\frac1{1+2^{n-2}z_n^{(n-2)/2}(1-z_n)}.
}
\tag{18}
\]

For `lambda < lambda_c`, the unique support optimum is no measurement (`t=0`) and

\[
\sup S_{n,\lambda}=1-\lambda/2.
\tag{19}
\]

At `lambda=lambda_c`, this endpoint ties a finite-strength strategy. For `lambda>lambda_c`, the unique nontrivial optimum is the high-`y` solution of (16), and

\[
\sup S_{n,\lambda}
=\frac\lambda2+
\frac{\lambda2^{n-2}y^{n-2}}{(1+y^2)^{n-1}}.
\tag{20}
\]

Thus the optimal local measurement strength jumps discontinuously at the transition for every `n >= 3`. For example, at `n=3`,

\[
z_3=\frac{\sqrt5-1}{2},
\qquad
\lambda_c(3)\approx0.6248,
\qquad
t_c\approx0.9717.
\]

As `n -> infinity`, put `a=2^{1-n}`. Equation (17) gives

\[
1-z_n=a+\frac{n-1}{2}a^2+O(n^2a^3),
\qquad
\lambda_c(n)=\frac23-\frac{a}{9}+O(na^2).
\tag{21}
\]

The optimal null therefore changes from learning nothing to measuring almost every bit projectively; the intermediate weak-measurement branch is not on the supported upper hull. Shared-random time sharing supplies the straight segment at the coexistence tangent but cannot raise a linear support value.

## 6. Collective classical-record comparator

Now allow all `n` carriers to participate in one collective instrument, while still requiring that only a classical outcome survives the commit point. The device can weakly measure the single global parity observable. With parity projectors `Pi_0,Pi_1`, use

\[
L_y=
\sqrt{\frac{1+t}{2}}\,\Pi_y
+\sqrt{\frac{1-t}{2}}\,\Pi_{1-y}.
\tag{22}
\]

This gives

\[
P_A^{\rm coll}=\frac{1+t}{2},
\qquad
F_R^{\rm coll}=f(t),
\]

and hence

\[
\boxed{
\sup S^{\rm coll}_{n,\lambda}
=\frac12+\frac12\sqrt{\lambda^2+(1-\lambda)^2},
}
\tag{23}
\]

independently of `n`. This is exactly the one-slot curve because collective access turns parity into one effective binary observable.

This value is optimal over arbitrary collective instruments with a classical surviving record, not merely over the displayed weak parity family. To see this, let `D=2^n` and `d=D/2`. For outcome `k`, write

\[
a_k=\sum_{x:\,\operatorname{par}(x)=0}p(k|x),
\qquad
b_k=\sum_{x:\,\operatorname{par}(x)=1}p(k|x).
\]

The optimal parity decoder gives

\[
P_A=\frac12+\frac1{2D}\sum_k|a_k-b_k|.
\]

The recovery lemma (4), followed by Cauchy--Schwarz within each parity sector, gives

\[
F_R\leq \frac1{D^2}\sum_k
\left(\sqrt{d a_k}+\sqrt{d b_k}\right)^2.
\]

Put `A_k=a_k/d`, `B_k=b_k/d`,

\[
T=\frac12\lVert A-B\rVert_1,
\qquad
C=\sum_k\sqrt{A_kB_k}.
\]

Then `P_A=(1+T)/2`, while `F_R <= (1+C)/2`. The classical fidelity--total-variation relation `T^2+C^2 <= 1` yields `F_R <= f(T)`. Optimizing the weighted score gives (23), and the weak parity instrument (22) attains equality.

For `n >= 2` and `0 < lambda < 1`, (23) is strictly larger than the streaming classical-memory value. One quick proof is

\[
\prod_i f(t_i)\leq f\!\left(\prod_i t_i\right),
\]

followed by concavity of `f` for randomized/adaptive profiles. Equality relevant to a support optimum occurs only at the trivial endpoints. In particular, a strong collective parity measurement retains `F_R=1/2`, whereas measuring all slots strongly and classically in the stream leaves only `F_R=2^{-n}`.

This strict gap between a streaming classical-memory comb and a collective classical-record instrument is the causal content of the result. A global `2^n`-level QND information--disturbance theorem does not capture it.

The qualifier **classical-record** is essential. A genuinely unrestricted collective device may retain coherent quantum memory and, like the coherent accumulator below, can attain `(P_A,F_R)=(1,1)`.

## 7. Honest coherent one-qubit accumulator

A device with one coherent qubit of inter-slot memory attains the algebraic point

\[
(P_A,F_R)=(1,1).
\]

Initialize `M=|0>`. At slot `i`, apply `CNOT_{A_i -> M}` and immediately return `A_i`. After the stream, `M` coherently stores the parity.

* In AUDIT, measure `M` in `Z`; the result is exactly the parity of the later `Z` outcomes on the references.
* In RETURN, hand `M` to the verifier/recovery station and apply the same CNOTs from the sequestered carriers to `M`. This uncomputes `M` to `|0>` and restores all `n` EPR pairs exactly.

The coherent accumulator, not the word "observer" or an Everettian interpretation, is the operational resource. The theorem therefore supplies an exact witness separating one persistent coherent qubit from unlimited adaptive classical memory in a genuinely streamed history.

## 8. Adversarial audit and scope limitations

1. **Reset ambiguity is material.** Equations (1)--(3) allow a classical transcript to control recovery and then be irreversibly reset, with the erasure environment untested. Checking only the visible controller does not certify global erasure. If RETURN instead demands a pure reset of the transcript *and every purifying environment*, any genuinely classical measurement record is disallowed and the null collapses to a different, much stronger model. The protocol must state which reset notion is tested.

2. **Immediate sequestration is essential.** If the device can delay a carrier, touch an earlier carrier again, or receive all carriers at once, it can implement (22); the streaming theorem no longer applies.

3. **No coherent inter-slot side channel is essential.** A hidden qubit, pre-shared entanglement across slot ancillas, a quantum clock, or a coherent control register invalidates factorization (5). That is precisely the resource the test is intended to witness, but enforcing the laboratory boundary is nontrivial.

4. **AUDIT must use the transcript only.** Giving the device returned carriers in AUDIT reduces the problem to a different discrimination task.

5. **The fresh inputs must be independent.** Correlated source states generally destroy the posterior-product identity behind (7).

6. **The RETURN score here is global entanglement fidelity.** Average single-pair fidelity, postselected recovery, or accepting only a random pair produces different tensorization laws and different bounds.

7. **No interpretational conclusion follows.** The separation is ordinary unitary quantum information: persistent coherent memory beats a classical-memory comb. It neither proves branching nor favors a quantum interpretation.

8. **Priority remains qualified.** The proof combines a flagged entanglement-fidelity bound, causal likelihood factorization, parity-bias multiplication, and a new-looking support optimization with a first-order transition. Each ingredient is close to established information--disturbance, parity-discrimination, sequential random-access-code, quantum-comb, quantum-seal, and memory-witness techniques. A search through 2026-08-12 found no source stating the full theorem under this access model. The defensible candidate claim is therefore the *exact streamed parity audit--global-return frontier and its phase diagram, to our knowledge*, not the general idea of witnessing quantum memory.

## 9. Cautious theorem-level novelty sentence after literature clearance

> For a sequential instrument with no coherent inter-round memory, we derive the exact tradeoff between late parity prediction from its classical transcript and joint entanglement recovery of all immediately returned probes. Arbitrary adaptive non-QND instruments reduce to a product-strength frontier, whose supported optimum undergoes a first-order transition for three or more time steps, while one coherent accumulator qubit attains perfect late audit and perfect reversible return.

This sentence must remain qualified as “to our knowledge.” Literature search can lower collision risk but cannot prove the absence of an equivalent result under different terminology.

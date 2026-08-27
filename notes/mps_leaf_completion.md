# Exact homogeneous-leaf reduction by local Weyl completion

**Date:** 22 August 2026

**Status:** theorem proved and construction verified numerically

> **Numerical update (22 August 2026).** The homogeneous-leaf reduction below
> is unchanged. The balanced three-effect benchmark remains valid, but a new
> four-effect leaf is stronger at larger AUDIT weights. At \(\lambda=0.6\) it
> reaches \(0.765898815264694\ldots\), disproving the former three-effect
> full-frontier conjecture. See `notes/interleaved_four_effect_frontier.md`.

## 1. Result

Consider the late-choice AUDIT--RETURN game with $n$ streamed
$q$-dimensional carriers, genuinely classical finite transcript, immediate
output sequestration, and coherent-memory profile
$\mathbf d=(d_1,\ldots,d_n)$. Let $D=q^n$. After refining every transcript to one
Kraus operator, write

\[
K_c:A_1\cdots A_n\longrightarrow B_1\cdots B_nM_n,
\qquad G_c=K_c^\dagger K_c,
\qquad t_c=\operatorname{Tr}G_c.
\tag{1}
\]

The vectorisation of $K_c$ is an open-boundary matrix-product state whose
bond at cut $i$ is at most $d_i$. For a Hilbert--Schmidt-normalised leaf
$K$, define

\[
\begin{aligned}
A_H(K)&=
\max_{\{E_s\}}
\sum_{x\in\mathbb F_q^n}
\operatorname{Tr}\!\left[
E_{Hx}\operatorname{Tr}_{B^n}
K|x\rangle\!\langle x|K^\dagger
\right],\\
R(K)&=\frac1D
\left(\operatorname{Tr}\sqrt{K^\dagger K}\right)^2,
\qquad \operatorname{Tr}K^\dagger K=1.
\end{aligned}
\tag{2}
\]

Here the POVM acts only on the terminal memory. Then the exact streamed
support function is

\[
\boxed{
\beta^{\rm stream}_{H,\mathbf d}(\lambda)
=\max_{\substack{\operatorname{Tr}K^\dagger K=1\\
\operatorname{TT-rank}_i(|K\rangle\!\rangle)\le d_i}}
\left[\lambda A_H(K)+(1-\lambda)R(K)\right].}
\tag{3}
\]

Equation (3) is an exact compact variational characterisation. It contains no
unbounded outcome alphabet, no adaptive tree, and no local-completeness
variables. It does not imply that the nonconvex maximum is available in
closed form.

## 2. Upper direction: every strategy is an average of leaves

For a refined physical strategy, instrument completeness gives

\[
\sum_cG_c=I_D,
\qquad
\sum_ct_c=D.
\tag{4}
\]

For $t_c>0$, put $\widehat K_c=K_c/\sqrt{t_c}$. The conditional optimal
AUDIT and RETURN scores of the leaf are respectively
$A_H(\widehat K_c)$ and $R(\widehat K_c)$. The exact flagged polar
recovery formula is used for RETURN. Consequently,

\[
P_{\rm A}=\sum_c\frac{t_c}{D}A_H(\widehat K_c),
\qquad
F_{\rm R}=\sum_c\frac{t_c}{D}R(\widehat K_c).
\tag{5}
\]

The coefficients in (5) are a probability distribution. Every refined leaf
inherits the physical coherent-memory profile, so its vectorisation has the
declared TT ranks. Taking the largest homogeneous support quotient proves the
upper inequality in (3). This argument is independent of the number of
leaves.

## 3. Lower direction: completing an arbitrary Choi MPS

Let $K$ be any normalised operator on the right-hand side of (3). Put its
vectorisation in sequential row-isometric gauge. Thus it has local tensors

\[
T_i:A_i\otimes M_{i-1}\longrightarrow B_i\otimes M_i
\tag{6}
\]

whose matrix elements reproduce $K$, and whose canonical condition is

\[
\operatorname{Tr}_{A_i}(T_i^\dagger T_i)=I_{M_{i-1}}.
\tag{7}
\]

Equation (7) is the row-isometry condition of a finite open-boundary MPS; it
follows constructively from successive Schmidt decompositions. Zero Schmidt
subspaces may be removed and smaller bonds may be padded to $d_i$.

Let

\[
W_{ab}=X^aZ^b,
\qquad a,b\in\mathbb F_q,
\tag{8}
\]

be the Heisenberg--Weyl unitary error basis. Its one-design identity is

\[
\sum_{a,b}
(W_{ab}^\dagger\otimes I)Y(W_{ab}\otimes I)
=qI_A\otimes\operatorname{Tr}_A Y.
\tag{9}
\]

At slot $i$, define $q^2$ local Kraus operators

\[
L^{(i)}_{ab}
=\frac1{\sqrt q}
T_i(W_{ab}\otimes I_{M_{i-1}}).
\tag{10}
\]

Equations (7)--(9) give local instrument completeness exactly:

\[
\sum_{a,b}(L^{(i)}_{ab})^\dagger L^{(i)}_{ab}
=I_{A_iM_{i-1}}.
\tag{11}
\]

The outcome $(a,b)$ is dephased and stored classically. No outcome-dependent
future tensor is required. For a complete Weyl transcript
$\boldsymbol\omega$, the global Kraus operator is

\[
K_{\boldsymbol\omega}
=q^{-n/2}K W_{\boldsymbol\omega},
\qquad
W_{\boldsymbol\omega}=\bigotimes_iW_{a_i b_i}.
\tag{12}
\]

There are (q^{2n}=D^2) leaves. Each has

\[
\operatorname{Tr}K_{\boldsymbol\omega}^\dagger
K_{\boldsymbol\omega}=D^{-1},
\tag{13}
\]

so every transcript has probability (D^{-2}) on the maximally mixed input.
Right-unitary invariance makes the conditional RETURN score equal to $R(K)$.
Moreover, $W_{ab}|x\rangle$ is a phase times $|x+a\rangle$. Relabelling the
terminal AUDIT answer by the known syndrome shift $H(a_1,\ldots,a_n)$ makes
every conditional AUDIT score equal to $A_H(K)$. Summing all transcripts
proves the lower inequality in (3).

## 4. Consequences

The theorem supplies an exact outcome-cardinality bound: $q^2$ outcomes per
slot suffice for every support direction. Arbitrary transcript adaptivity may
be useful for finding a compact implementation, but is unnecessary for
attaining the optimum. The complete strategy produced by (10) respects every
memory cut and immediate sequestration.

The reduction is specific to games whose leaf score is invariant under the
chosen monomial unitary error basis up to a known classical relabelling. The
RETURN score has this invariance, and a linear syndrome AUDIT has it because
Weyl shifts translate the syndrome. It should not be asserted for an
arbitrary terminal task without checking that covariance.

The construction also shows that a leafwise counterexample to a proposed
homogeneous inequality is automatically a complete-strategy counterexample.
In particular, the previously tested coefficient (3/2) in the interleaved
linear-tail conjecture is false globally, not merely after postselection.

## 5. Four-slot interleaved benchmark

For $q=2$, $n=4$, and $d_i=2$, a general Choi-MPS leaf optimised at balanced weight
has

\[
A_{H_{\rm I}}(K)=0.620085075585902\ldots,
\qquad
R(K)=0.899520492116986\ldots,
\tag{14}
\]

and hence

\[
\boxed{S_{1/2}(K)=0.759802783851444\ldots.}
\tag{15}
\]

The local Pauli completion has four outcomes at each of four slots and 256
complete transcripts. Direct contraction gives total probability
$0.999999999999992$, maximum local-completeness residual
$1.05\times10^{-15}$, and reproduces (14)--(15) within
$7\times10^{-15}$. Thus (15) is a physical lower bound, improving the
earlier complete ternary strategy $0.759448970317260\ldots$.

Continuation of the same complete three-effect MPS branch gives a coexistence crossing
with the no-record point near

\[
\lambda_c\simeq0.441437845,
\tag{16}
\]

and approaches the exact perfect-AUDIT endpoint
$(P_{\rm A},F_{\rm R})=(1,1/4)$. The digits in (15)--(16) are variational
results, not globally certified upper bounds. What is now exact is the
reduction (3): closing the frontier means certifying the compact MPS maximum,
not searching over unbounded adaptive instruments. The later four-effect
construction changes the candidate maximiser, not this reduction theorem.

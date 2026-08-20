# A quantitative robust bound at a full-crossing cut

**Status (19 August 2026).** This note proves a non-sharp but fully
quantitative version of the crossed-pair endpoint theorem. It covers
arbitrary adaptive non-QND instruments, an unrestricted finite classical
transcript, and transcript-conditioned AUDIT and RETURN decoders. For the
four-slot interleaved instance,

\[
F_{\rm R}\leq \frac14+\frac72\sqrt{1-P_{\rm A}}
 +\frac{\sqrt{21}}2(1-P_{\rm A})^{1/4}.                       \tag{1}
\]

It follows that the streamed support function is strictly below the common
static qubit ceiling whenever

\[
0.997339869\leq\lambda<1.                                    \tag{2}
\]

The threshold is conservative. Its purpose is to replace the previous
compactness-only neighbourhood by an explicit certified interval. The
exponent \(1/4\) in (1) is a proof loss; a legal weak-measurement family
shows that the optimal exponent can be no larger than \(1/2\).

## 1. Interface and leaf scores

Let \(H\in\mathbb F_2^{2\times n}\) have rank two. Suppose that a cut after
slot \(i\) is full crossing:

\[
H=(H_L\mid H_R),\qquad
\operatorname{rank}H_L=\operatorname{rank}H_R=2.              \tag{3}
\]

Only a qubit \(M\) crosses the cut coherently. Earlier output carriers are
sequestered, and the classical transcript is genuinely dephased. AUDIT sees
the transcript and terminal qubit, but not the emitted carriers. RETURN may
act jointly on every returned carrier and \(M\), conditioned on the
transcript. Put \(D=2^n\).

Refine every outcome to one Kraus label and reveal that refinement to both
decoders. This is a relaxation: it can only increase both scores and does
not enlarge the coherent memory. For a complete refined leaf \(c\), write

\[
p_c(x)=\lVert K_c|x\rangle\rVert^2,\qquad
m_c=\sum_xp_c(x).                                             \tag{4}
\]

Let \(\rho_{c,x}\) be the corresponding unnormalised terminal-qubit state
after all output carriers are traced out. If
\(\{E_{c,s}:s\in\mathbb F_2^2\}\) is the leaf's AUDIT POVM, define

\[
a_c=\sum_x\operatorname{Tr}(E_{c,Hx}\rho_{c,x}),\qquad
e_c=m_c-a_c.                                                  \tag{5}
\]

Trace preservation and the uniform computational-basis ensemble give

\[
\sum_cm_c=D,\qquad
\sum_ce_c=D(1-P_{\rm A}).                                     \tag{6}
\]

Every estimate below is homogeneous in a leaf, so zero-mass leaves may be
discarded.

## 2. Right-canonical form at the cut

Fix one leaf and suppress \(c\). Write \(x=(z,y)\), where
\(z\in\mathbb F_2^i\) and \(y\in\mathbb F_2^{n-i}\). Causality and immediate
sequestration imply

\[
K|z,y\rangle=(I_{B_L}\otimes L_y)|\psi_z\rangle,              \tag{7}
\]

with \(|\psi_z\rangle\in B_L\otimes M_i\) and
\(L_y:M_i\to B_R\otimes M_n\). No QND assumption is present.

One may choose the gauge

\[
\sum_yL_y^\dagger L_y=I,\qquad
\sum_z\lVert\psi_z\rVert^2=m.                                 \tag{8}
\]

Here \(I\) means the identity on the effective support, whose dimension is
at most two. To obtain the gauge, let \(Q=\sum_yL_y^\dagger L_y\), absorb
\(Q^{1/2}\) into every prefix vector, and absorb its inverse on
\(\operatorname{supp}Q\) into every continuation. Components in \(\ker Q\)
never contribute to the complete leaf and can be removed. In this gauge each
\(L_y\) is a contraction and

\[
\sum_yp(z,y)=\lVert\psi_z\rVert^2.                             \tag{9}
\]

Dividing all prefix vectors by \(\sqrt m\) reduces the leafwise proof to
\(m=1\).

## 3. Two elementary lemmas

### 3.1 Two syndromes carry almost all likelihood

For \(s\in\mathbb F_2^2\), put

\[
q_s=\sum_{x:Hx=s}p(x).
\]

Because \(E_s\geq0\),

\[
a\leq\sum_sq_s\lVert E_s\rVert_\infty.
\]

Moreover,

\[
0\leq\lVert E_s\rVert_\infty\leq1,\qquad
\sum_s\lVert E_s\rVert_\infty
\leq\sum_s\operatorname{Tr}E_s=2.
\]

The resulting two-unit linear programme is maximised by the two largest
\(q_s\). If \(T\) is that pair, then

\[
\sum_{s\notin T}q_s\leq m-a=e.                                \tag{10}
\]

Every two-point subset of \(\mathbb F_2^2\) is an affine line, so
\(T=\{s:u\cdot s=b\}\) for some nonzero \(u\).

### 3.2 Approximate discrimination forces approximate injectivity

The next lemma is the quantitative replacement for the rank-one
contradiction in the exact endpoint proof.

> **Lemma 1.** Let \(L:\mathbb C^2\to\mathcal Y\otimes\mathbb C^2\) be a
> contraction with smaller singular value \(\sigma\). Let two possibly
> spectator-entangled inputs have output probabilities \(p_0,p_1\).
> Suppose two distinct terminal-qubit POVM effects identify them with error
> probabilities \(\epsilon_0,\epsilon_1\). If
> \(r=\min(p_0,p_1)\), then
> \[
> \sigma\geq
> \frac{(\sqrt r-\sqrt{\epsilon_0+\epsilon_1})_+}{1+\sqrt2}.
>                                                               \tag{11}
> \]

Let \(L_0\) be the best rank-one approximation to \(L\), so
\(\lVert L-L_0\rVert=\sigma\). Denote the actual and rank-one outputs by
\(v_j\) and \(v_j^{(0)}\). The input norms are at most one, hence

\[
\lVert v_j-v_j^{(0)}\rVert\leq\sigma,\qquad
\lVert v_j^{(0)}\rVert\geq\sqrt r-\sigma.                     \tag{12}
\]

If \(\sigma\geq\sqrt r\), (11) is immediate because its right-hand side is
at most \(\sqrt r/(1+\sqrt2)\). Hence assume
\(\sigma<\sqrt r\) in what follows, so the lower bound in (12) is positive.

Since \(L_0\) has rank one, even a spectator-entangled input gives
\(v_j^{(0)}=|\chi_j\rangle\otimes|w\rangle\), with the same normalised
future vector \(|w\rangle\). Let \(\alpha_j\) be the probability of the
declared label on \(|w\rangle\). The two distinct POVM effects obey
\(\alpha_0+\alpha_1\leq1\). Applying the reverse triangle inequality after
the square root of the incorrect-label effect gives

\[
\sqrt{1-\alpha_j}\,(\sqrt r-\sigma)
\leq\sqrt{\epsilon_j}+\sigma.                                 \tag{13}
\]

Since \((1-\alpha_0)+(1-\alpha_1)\geq1\), summing the squares of (13) yields

\[
(\sqrt r-\sigma)^2
\leq(\sqrt{\epsilon_0+\epsilon_1}+\sqrt2\,\sigma)^2,
\]

which proves (11). The proof is unchanged if an artificial orthogonal label
combines several prefix vectors into either input ensemble.

## 4. A robust \(D/4\)-word lemma

Let

\[
J=2^{n-i}
\]

be the number of suffix words at the full-crossing cut. For a fixed suffix
\(y\), the prefixes whose syndrome lies in \(T\) split into two equal
classes, one for each label in \(T\). Let their total probabilities be
\(A_y,B_y\), and define

\[
r_y=\min(A_y,B_y),\qquad R=\sum_yr_y.                          \tag{14}
\]

There are two useful candidate subsets of input words.

For each suffix, first retain the larger of its two syndrome classes. This
keeps at most \(2^{i-2}J=D/4\) words and discards mass at most

\[
\delta\leq e+R.                                               \tag{15}
\]

For the second candidate, choose \(y_*\) with
\(r_{y_*}\geq R/J\). Combine all compatible prefix vectors for either
syndrome label into two orthogonal direct sums and apply Lemma 1 to
\(L_{y_*}\). Their combined AUDIT error is at most \(e\), so

\[
\sigma_{\min}(L_{y_*})
\geq\frac{(\sqrt{R/J}-\sqrt e)_+}{1+\sqrt2}.                  \tag{16}
\]

Retain the affine half of the prefixes compatible with \(T\) at \(y_*\), and
the affine half of the suffixes with matching \(uH_R\)-parity. This again
contains \(D/4\) words. Every rejected prefix has syndrome outside \(T\) at
\(y_*\), so its total probability at that suffix is at most \(e\).
Equations (9) and (16) imply

\[
\delta\leq e+\frac{e}{\sigma_{\min}(L_{y_*})^2}.               \tag{17}
\]

Define

\[
C_i=4\sqrt J=4\,2^{(n-i)/2}.                                  \tag{18}
\]

If \(e\geq C_i^{-2}\), the trivial estimate
\(\delta\leq1\leq C_i\sqrt e\) suffices. Suppose
\(e<C_i^{-2}\). If \(R<12Je\), (15) gives

\[
\delta<(1+12J)e\leq C_i\sqrt e.
\]

If \(R\geq12Je\), then

\[
\sqrt{R/J}-\sqrt e
\geq\frac{1-1/\sqrt{12}}{\sqrt J}\sqrt R.
\]

Using the better of (15) and (17),

\[
\begin{aligned}
\delta
&\leq e+\min\left\{R,
\frac{(1+\sqrt2)^2J}{(1-1/\sqrt{12})^2}\frac eR\right\}\\
&\leq e+
\frac{(1+\sqrt2)\sqrt J}{1-1/\sqrt{12}}\sqrt e
\leq C_i\sqrt e.                                             \tag{19}
\end{aligned}
\]

The final inequality follows from \(e<C_i^{-2}\) and \(J\geq1\).
Restoring the leaf normalisation proves the following statement.

> **Lemma 2 (robust full-crossing support).** Every refined leaf admits
> \(S_c\subset\mathbb F_2^n\), with \(|S_c|\leq D/4\), such that
> \[
> \sum_{x\notin S_c}p_c(x)\leq C_i\sqrt{m_ce_c}.               \tag{20}
> \]

For the interleaved four-slot instance, \(J=4\). Repeating the same two-case
calculation with thresholds \(e=1/49\) and \(R=48e\) improves \(C_i=8\) to

\[
C_i=7.                                                        \tag{21}
\]

In the small-\(R\) case,
\(e+R<49e\leq7\sqrt e\). In the large-\(R\) case, the coefficient of
\(\sqrt e\) is

\[
\frac17+\frac{4\sqrt3(1+\sqrt2)}{2\sqrt3-1}
=6.930792\ldots<7.                                            \tag{22}
\]

No diagonality, covariance, Clifford structure, or finite a priori bound on
the number of leaves has been used.

## 5. Robust RETURN theorem

For a refined leaf, let \(G_c=K_c^\dagger K_c\). Optimal flagged polar
recovery and computational-basis pinching give

\[
F_{{\rm R},c}\leq\frac1{D^2}
\left(\sum_x\sqrt{p_c(x)}\right)^2.                            \tag{23}
\]

If \(S_c\) is supplied by Lemma 2, pad it arbitrarily to exactly \(D/4\)
words; this cannot increase its discarded mass \(\delta_c\).
Cauchy--Schwarz on the padded set and its \(3D/4\)-word complement gives

\[
\begin{aligned}
\left(\sum_x\sqrt{p_c(x)}\right)^2
&\leq
\left(\sqrt{\frac D4(m_c-\delta_c)}
      +\sqrt{\frac{3D}4\delta_c}\right)^2\\
&\leq\frac D4m_c+\frac D2\delta_c
      +\frac D2\sqrt{3m_c\delta_c}.                           \tag{24}
\end{aligned}
\]

Substitute (20), sum over all leaves, and use (6) and Hölder's inequality.
With \(\epsilon=1-P_{\rm A}\),

\[
\sum_c\sqrt{m_ce_c}\leq D\sqrt\epsilon,\qquad
\sum_cm_c^{3/4}e_c^{1/4}\leq D\epsilon^{1/4}.                 \tag{25}
\]

Therefore

\[
\boxed{
F_{\rm R}\leq\frac14+\frac{C_i}{2}\sqrt\epsilon
+\frac{\sqrt{3C_i}}2\epsilon^{1/4}.}                          \tag{26}
\]

The proof permits arbitrarily many refined leaves because (25) has no
alphabet-size factor. If the physical strategy was initially unrefined,
revealing its Kraus labels increases \(P_{\rm A}\) and \(F_{\rm R}\). The
right-hand side of (26) decreases with \(P_{\rm A}\), so (26) also holds for
the original unrefined scores.

For \(H_{\rm I}\), equations (21) and (26) give (1).

## 6. Explicit interval of strict order dependence

For the four-slot problem, set

\[
b=\frac72,\qquad c=\frac{\sqrt{21}}2,\qquad
q=1-\lambda,\qquad R_\lambda=\frac{\lambda}{q}.
\]

Writing \(\epsilon=u^4\), equation (1) gives

\[
\lambda P_{\rm A}+qF_{\rm R}
\leq\lambda+\frac q4
+q\bigl(cu+bu^2-R_\lambda u^4\bigr).                          \tag{27}
\]

The static qubit ceiling is

\[
B_{4,2}(\lambda)=\frac{1+\sqrt{\lambda^2+q^2}}2,
\]

and

\[
B_{4,2}(\lambda)-\left(\lambda+\frac q4\right)
>\frac q4\qquad(0<q<1).                                      \tag{28}
\]

Define

\[
u_*=\frac{\sqrt{301}-3\sqrt{21}}{28},\qquad
R_*=\frac{7u_*+\sqrt{21}/2}{4u_*^3}
=374.921248196\ldots .                                       \tag{29}
\]

The one-variable identity

\[
R_*=\sup_{u>0}\frac{bu^2+cu-1/4}{u^4}                         \tag{30}
\]

follows by differentiation. The unique positive stationary point solves
\(2bu^2+3cu-1=0\), which gives \(u_*\). Consequently, if
\(\lambda/(1-\lambda)\geq R_*\), the last parenthesis in (27) is at most
\(1/4\). Equation (28) is strict, so

\[
\boxed{
\beta^{\rm stream}_{H_{\rm I},2}(\lambda)
<B_{4,2}(\lambda)
\quad\text{for}\quad
0.997339868377\ldots
=\frac{R_*}{1+R_*}\leq\lambda<1.}                            \tag{31}
\]

The grouped ordering attains the right-hand side for every \(\lambda\), so
(31) is an explicit temporal-order gap.

## 7. The complementary width-one construction

The full-crossing condition has a clean converse at the level of static
attainability for every rank-two binary check matrix. Define

\[
\tau(H)=\max_{1\leq i<n}
\left(\operatorname{rank}H_{\leq i}
+\operatorname{rank}H_{>i}-2\right).                          \tag{32}
\]

> **Proposition 3 (rank-two order dichotomy).** If \(\tau(H)\leq1\), a
> streamed qubit strategy attains the complete static boundary
> \[
> P_{\rm A}=\frac{1+t}{2},\qquad
> F_{\rm R}=\frac{1+\sqrt{1-t^2}}2
> \quad(0\leq t\leq1).                                        \tag{33}
> \]
> If \(\tau(H)=2\), a full-crossing cut exists. Perfect AUDIT then obeys
> \(F_{\rm R}\leq1/4\), and the robust estimate (26) holds with the
> constant of any full-crossing cut.

Let \(j\) be the first slot with
\(\operatorname{rank}H_{\leq j}=2\). A single column has rank at most one,
so \(j\geq2\). Since \(\tau(H)\leq1\),

\[
\operatorname{rank}H_{>j}\leq1.                               \tag{34}
\]

Choose a nonzero row coefficient whose character \(a\) annihilates
\(H_{>j}\). It is nonzero on the complete word because
\(H_{\leq j}\) has rank two. Let \(\ell\leq j\) be the last nonzero
coordinate of \(a\), and choose an independent character \(b\). If
\(\ell<j\), the restricted row code through \(\ell\) has rank one and
\(a_{\leq\ell}\neq0\). Of the two coefficient vectors independent of \(a\),
exactly one has the same restriction as \(a\); choose it for \(b\). Thus
\(b_{\leq\ell}=a_{\leq\ell}\) whenever the handoff precedes the rank-two
slot.

For every cut at which \(a\) remains live, \(i<j\), and minimality of \(j\)
gives

\[
\operatorname{rank}_{\mathbb F_2}
\{a_{\leq i},b_{\leq i}\}\leq1.                               \tag{35}
\]

A qubit can therefore accumulate the unique live prefix form. If \(a\)
completes before \(j\), the chosen \(b\) already has the same prefix value.
If \(a\) completes at \(j\), the old live form and current carrier span the
new rank-two prefix. A reversible two-bit linear transformation can expose
\(a\), weakly read it, and then place \(b\) in memory while the complementary
form leaves in the carrier. Afterwards only \(b\) is accumulated. Zero
columns do nothing, and the argument also covers \(j=n\).

At the terminal cut, the flag weakly identifies \(a\) and the qubit stores
\(b\). The pair identifies the complete syndrome. Reversing the linear
stream on RETURN gives (33). Hence, for rank-two checks, \(\tau\leq1\)
means that every static support direction is attainable, whereas
\(\tau=2\) forces a strict perfect-AUDIT coherence loss. This all-length
dichotomy is unaffected by repeated or zero columns.

## 8. Sharpness and the remaining gap

The legal strategy that stores one syndrome row coherently and makes equal
local weak measurements of the other row has

\[
P_{\rm A}=\frac{1+t^2}{2},\qquad
F_{\rm R}=\left(\frac{1+\sqrt{1-t^2}}2\right)^2.
\]

If \(\epsilon=1-P_{\rm A}\), then

\[
F_{\rm R}=\frac14+\sqrt{\frac\epsilon2}+\frac\epsilon2.        \tag{36}
\]

Thus no endpoint theorem can replace the correction in (1) by
\(o(\sqrt\epsilon)\). The expected sharp robust order is
\(\sqrt\epsilon\), with any universal leading constant at least
\(1/\sqrt2\). Equation (1) is weaker because Lemma 2 first bounds the
discarded likelihood by \(O(\sqrt e)\), and the RETURN pinching bound takes
another square root of that tail. A sharp square-root law will require
either a linear robust support statement after summing completeness across
leaves, or a direct recovery argument that avoids the second square-root
loss. The exact interior frontier remains open.

## 9. Self-audit

The proof uses the frozen interface in four essential places. Immediate
sequestration gives the spectator identity in (7). A genuinely classical
transcript permits Kraus refinement without an uncharged coherent
purification. AUDIT's lack of carrier access reduces its decoder to a qubit
POVM, which is needed in (10) and Lemma 1. Unconditional RETURN, with no
reference access or postselection, is needed for the flagged polar bound
(23).

The right-canonical gauge is a factorisation identity, not a physical
post-processing step, so it does not assume that a normalized leaf is itself
a trace-nonincreasing operation. If the effective support in (8) is
one-dimensional, Lemma 1 has \(\sigma=0\) and excludes the large-\(R\)
case directly; division by \(\sigma_{\min}\) is used only when (16) makes it
strictly positive. The artificial direct sum used before (16) only labels
different prefix basis words in the proof; it is not a side channel given to
the device. Prefix entanglement with \(B_L\) is harmless because both \(L_y\)
and its rank-one approximation act only on \(M_i\).

The constants can be checked independently from the two scalar inequalities

\[
\frac17+\frac{4\sqrt3(1+\sqrt2)}{2\sqrt3-1}<7,\qquad
2\left(\frac72\right)u_*^2
+3\left(\frac{\sqrt{21}}2\right)u_*=1.
\]

The first certifies \(C_i=7\); the second certifies the stationary point in
(30). Substitution gives \(R_*=374.921248196\ldots\) and
\(R_*/(1+R_*)=0.997339868377\ldots\). No numerical optimization enters the
theorem.

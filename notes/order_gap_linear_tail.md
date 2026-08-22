# Linear rank-tail robustness from causal list decoding

**Date:** 21 August 2026
**Status:** proved; algebraically and numerically audited
**Scope:** arbitrary adaptive non-QND instruments, genuinely classical
transcript, immediate output sequestration, terminal syndrome AUDIT, and
unconditional all-carrier RETURN

## 1. Main result

Let \(H\in\mathbb F_q^{r\times n}\) be partitioned into \(m\) consecutive
full-rank blocks,

\[
H=[H^{(1)}|\cdots|H^{(m)}],\qquad \operatorname{rank}H^{(j)}=r.
\tag{1}
\]

Put \(N=q^r\), \(D=q^n\), and let the coherent memory at the \(j\)-th
block boundary have dimension at most \(d_j\leq N\). Define

\[
k=D\prod_{j=1}^m\frac{d_j}{N}.
\tag{2}
\]

The integer \(k=q^{n-mr}\prod_jd_j\) is the exact support bound at perfect
AUDIT. For a complete refined transcript \(c\), let \(K_c\) be its
sequential Kraus operator, \(G_c=K_c^\dagger K_c\), and

\[
t_c=\sum_{\ell>k}\lambda_\ell(G_c),
\tag{3}
\]

where the eigenvalues are in nonincreasing order.  If
\(\epsilon=1-P_{\rm A}\), then

\[
\boxed{\sum_c t_c\leq mD\epsilon.}
\tag{4}
\]

Thus the approximate theorem has a *linear* spectral tail.  It does not
pass through a leafwise approximate-injectivity estimate and therefore does
not incur the square-root loss in the earlier robust support lemma.

Let

\[
\alpha=\frac{k}{D}=\prod_{j=1}^m\frac{d_j}{N},\qquad
\theta=\min\{m\epsilon,1-\alpha\}.
\tag{5}
\]

Optimal flagged polar recovery then obeys

\[
\boxed{
F_{\rm R}\leq
\alpha+(1-2\alpha)\theta
+2\sqrt{\alpha(1-\alpha)\theta(1-\theta)}.}
\tag{6}
\]

The expression is valid for every \(0<\alpha<1\). As a function of the
actual tail fraction it increases up to \(1-\alpha\), precisely the universal
tail ceiling.

For the four-slot interleaved binary benchmark,

\[
m=2,\qquad D=16,\qquad N=4,\qquad d_1=d_2=2,
\qquad \alpha=\frac14.
\tag{7}
\]

Consequently, with

\[
\theta=\min\{2(1-P_{\rm A}),3/4\},
\tag{8}
\]

every legal strategy satisfies

\[
\boxed{
F_{\rm R}\leq
\frac14+\frac12\theta
+\frac{\sqrt3}{2}\sqrt{\theta(1-\theta)}.}
\tag{9}
\]

Near perfect AUDIT this is

\[
F_{\rm R}
\leq \frac14+(1-P_{\rm A})
+\sqrt{\frac32(1-P_{\rm A})[1-2(1-P_{\rm A})]}.
\tag{10}
\]

It has the optimal square-root exponent.  The leading coefficient
\(\sqrt{3/2}\) is not claimed sharp; the legal weak-record family requires
at least \(1/\sqrt2\).

## 2. Dimension-to-list lemma

The proof uses the following elementary observation.

> **Lemma 1 (quantum dimension gives a classical list).**  Let
> \(\{\rho_x\}_{x\in\mathcal X}\) be subnormalised states on a
> \(d\)-dimensional system and let \(\{E_x\}_{x\in\mathcal X}\) be a POVM.
> If \(p_x=\operatorname{Tr}\rho_x\), then
> \[
> \sum_x\operatorname{Tr}(E_x\rho_x)
> \leq \max_{L\subseteq\mathcal X,\ |L|\leq d}\sum_{x\in L}p_x.
> \tag{11}
> \]

Indeed,

\[
\operatorname{Tr}(E_x\rho_x)
\leq p_x\lVert E_x\rVert_\infty,
\tag{12}
\]

while \(0\leq\lVert E_x\rVert_\infty\leq1\) and

\[
\sum_x\lVert E_x\rVert_\infty
\leq\sum_x\operatorname{Tr}E_x=d.
\tag{13}
\]

Maximising the resulting fractional \(d\)-unit linear programme selects the
\(d\) largest \(p_x\). The lemma remains valid when some states have zero
mass and when fewer than \(d\) labels are present.

## 3. A high-probability list at every full-rank boundary

Refine every physical instrument outcome to one Kraus label and reveal the
refinement to both late decoders. This is a relaxation and preserves every
coherent-memory bound. Write a word as a sequence of block inputs
\(x=(x_1,\ldots,x_m)\), and define the block and cumulative syndromes

\[
v_j=H^{(j)}x_j,\qquad
u_j=v_1+\cdots+v_j\in\mathbb F_q^r.
\tag{14}
\]

Let \(a_j\) be the refined classical transcript through boundary \(j\).
Conditioned on \(a_j\) and \(u_j=u\), trace the sequestered prefix carriers
and all within-fibre variables. This gives a subnormalised state
\(\rho_{a_j,u}\) on the \(d_j\)-dimensional cut memory.

Fix the future block inputs and compose the suffix instrument with the
terminal AUDIT POVM. Relabel an AUDIT answer \(u_m\) as the prefix answer

\[
u_j=u_m-(v_{j+1}+\cdots+v_m).
\tag{15}
\]

Because future inputs are independent and uniform, averaging these pulled-
back effects over all future inputs gives a POVM on the cut memory.  It guesses
\(u_j\) with exactly the contribution of the original AUDIT decoder. This
remains true with prefix--carrier entanglement: the future comb acts only on
the cut memory, and the earlier carriers are spectators that are traced in
AUDIT.

Apply Lemma 1 for every prefix transcript \(a_j\). There is a list

\[
L_j(a_j)\subseteq\mathbb F_q^r,\qquad |L_j(a_j)|\leq d_j,
\tag{16}
\]

such that, in the actual experiment,

\[
\Pr[U_j\in L_j(A_j)]\geq P_{\rm A}.
\tag{17}
\]

No joint measurability of the boundary POVMs for different \(j\) is needed:
the lists are classical objects used only in the proof.

## 4. Intersecting the causal lists

For a complete transcript \(c\), let \(a_j(c)\) denote its prefix through
boundary \(j\), and define

\[
S_c=\{x:u_j(x)\in L_j(a_j(c))\text{ for all }j\}.
\tag{18}
\]

The map between cumulative labels and block labels,

\[
(u_1,\ldots,u_m)\longleftrightarrow(v_1,\ldots,v_m),
\quad v_1=u_1,\quad v_j=u_j-u_{j-1},
\tag{19}
\]

is bijective. Full rank of every block gives a fibre of size
\(q^{n_j-r}\). Therefore

\[
|S_c|\leq q^{n-mr}\prod_jd_j=k.
\tag{20}
\]

Let \(E_j\) be the event in equation (17). The union bound gives

\[
\Pr\!\left[\bigcap_{j=1}^m E_j\right]
\geq1-\sum_j\Pr[E_j^c]
\geq1-m\epsilon.
\tag{21}
\]

In likelihood notation, \(p_c(x)=\lVert K_c|x\rangle\rVert^2\), equations
(18)--(21) are exactly

\[
\sum_c\sum_{x\notin S_c}p_c(x)\leq mD\epsilon.
\tag{22}
\]

For \(G_c=K_c^\dagger K_c\), the diagonal entries in the computational
basis are \(p_c(x)\). Schur--Horn majorisation implies

\[
\sum_{\ell=1}^{k}\lambda_\ell(G_c)
\geq\sum_{x\in S_c}p_c(x).
\tag{23}
\]

Subtracting from \(m_c=\operatorname{Tr}G_c\) and summing equation (22)
proves the linear tail bound (4).

The argument also covers \(m\epsilon\geq1\), where equation (4) follows
trivially from \(\sum_ct_c\leq\sum_c\operatorname{Tr}G_c=D\).

## 5. From spectral tail to RETURN

Put \(m_c=\operatorname{Tr}G_c\). Cauchy--Schwarz on the first \(k\) and
last \(D-k\) singular values gives

\[
\operatorname{Tr}\sqrt{G_c}
\leq\sqrt{k(m_c-t_c)}+\sqrt{(D-k)t_c}.
\tag{24}
\]

Let \(T=\sum_ct_c=\theta'D\). Instrument completeness gives
\(\sum_cm_c=D\). Summing the square of equation (24) and applying
Cauchy--Schwarz once more yields

\[
F_{\rm R}
\leq
\alpha+(1-2\alpha)\theta'
+2\sqrt{\alpha(1-\alpha)\theta'(1-\theta')}.
\tag{25}
\]

Ky Fan anti-norm superadditivity and \(\sum_cG_c=I_D\) also give

\[
T=\sum_c\sum_{\ell>k}\lambda_\ell(G_c)
\leq\sum_{\ell>k}\lambda_\ell(I_D)=D-k,
\tag{26}
\]

so \(0\leq\theta'\leq\min\{m\epsilon,1-\alpha\}\). The right-hand side of
equation (25) is the squared overlap of the unit vectors

\[
(\sqrt{\alpha},\sqrt{1-\alpha})
\quad\text{and}\quad
(\sqrt{1-\theta'},\sqrt{\theta'}).
\tag{26a}
\]

It increases on \(0\leq\theta'\leq1-\alpha\). Substitution of equation (5)
therefore proves equation (6).

The flagged polar-recovery formula already optimises an arbitrary
transcript-conditioned joint RETURN decoder on every emitted carrier and the
terminal memory.  Dropping a demanded memory reset only relaxes the task.

## 6. Explicit interleaved support certificate

For equation (9), write

\[
\theta=\sin^2\phi,qquad 0\leq\phi\leq\frac\pi3.
\tag{27}
\]

Then

\[
P_{\rm A}\leq1-\frac{\theta}{2},\qquad
F_{\rm R}\leq
\left(\frac12\cos\phi+\frac{\sqrt3}{2}\sin\phi\right)^2.
\tag{28}
\]

For audit weight \(\lambda\) and \(q=1-\lambda\), maximising the resulting
quadratic form gives

\[
\boxed{
\beta^{\rm stream}_{H_{\rm I},2}(\lambda)
\leq U_{\rm I}(\lambda)
=\frac12+\frac\lambda4
+\frac14\sqrt{7\lambda^2-10\lambda+4}.}
\tag{29}
\]

The Perron eigenvector of the positive \(2\times2\) quadratic form lies in
the interval in equation (27); at \(\lambda=0\) it reaches the endpoint
\(\theta=3/4\), and thereafter \(\theta\) decreases.

The grouped ordering attains the static qubit value

\[
B_{4,2}(\lambda)
=\frac{1+\sqrt{\lambda^2+(1-\lambda)^2}}2.
\tag{30}
\]

Direct algebra gives

\[
U_{\rm I}(\lambda)<B_{4,2}(\lambda)
\quad\text{for}\quad \frac37<\lambda<1,
\tag{31}
\]

with equality at \(\lambda=0,3/7,1\) where relevant.  One way to check all
possible crossings is to eliminate the two square roots; the resulting
polynomial factors as

\[
-4\lambda^2(\lambda-1)(7\lambda-3).
\tag{31a}
\]

Checking the sign before squaring then gives equation (31).  Therefore temporal
order is now certified to change the support function throughout the broad
interior interval

\[
\boxed{\frac37<\lambda<1,}
\tag{32}
\]

including balanced weight.  Equation (32) replaces the previous conservative
near-endpoint interval \(0.997339868\ldots\leq\lambda<1\).

This is an upper certificate, not the exact interleaved frontier.  At
\(\lambda=1/2\), it gives

\[
\beta^{\rm stream}_{H_{\rm I},2}(1/2)
\leq\frac58+\frac{\sqrt3}{8}
=0.841506350946\ldots,
\tag{33}
\]

whereas the best verified complete strategy currently gives
\(0.759448970317\ldots\). Closing that residual interval remains a distinct
optimisation problem.

## 7. Assumption audit

Immediate sequestration is used when a future AUDIT strategy is pulled back
to a POVM on the cut memory: past carriers cannot participate.  A genuinely
classical transcript is used to condition separate lists without smuggling a
coherent purification across a boundary.  Full rank of each block is used
only in the uniform fibre count (20).  AUDIT must lack carrier and reference
access.  RETURN must be unconditional and have no reference access for the
flagged polar bound.

Kraus refinement can only help both tasks, and all inequalities are
independent of the number of refined leaves.  The proof never assumes QND,
Clifford, diagonal, Lüders, covariant, projective, or binary-outcome local
instruments.

## 8. Independent checks and rejected strengthening

The public verifier evaluates the causal lists and spectral tail of the stored
complete ternary-outcome strategy without using the proof formulas. It finds

\[
P_{\rm A}=0.625754561820\ldots,qquad
\frac1D\sum_ct_c=0.452260604\ldots,
\]

while equation (4) permits
\(2(1-P_{\rm A})=0.748490876\ldots\). The best four-word coordinate mass and
the spectral top-four mass coincide for this artifact to numerical precision,
so the Schur--Horn step is also directly visible. Symbolic elimination and a
dense floating-point grid independently reproduce equation (29), its balanced
value, and the only interior crossing at \(\lambda=3/7\).

A tempting stronger conjecture replaces the coefficient \(m=2\) in equation
(4) by \(3/2\). It holds on all complete binary, ternary, and quaternary trees
tested so far, but it is not proved. More importantly, it fails for an
admissible postselected tensor-train leaf once sibling completeness is
removed. Any proof of the \(3/2\) coefficient would therefore have to use the
shared local completeness relations between leaves. Neither that conjecture
nor its implied all-weight order gap is claimed here.

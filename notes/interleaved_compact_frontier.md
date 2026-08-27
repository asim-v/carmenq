# Compact variational frontier for the four-slot interleaved game

**Date:** 22 August 2026
**Status:** exact physical construction and exact homogeneous-leaf reduction;
superseded as a full-frontier candidate by a stronger four-effect branch

> **Correction (22 August 2026).** The three-effect family remains a valid
> physical lower bound and is exposed near balanced weight, but it is not the
> full attainable curve. A symmetric four-effect bond-two MPS reaches
> \(0.765898815264694\ldots\) at \(\lambda=0.6\), compared with
> \(0.755705934586018\ldots\) here. See
> `notes/interleaved_four_effect_frontier.md`. Statements below describing 3E
> as the frontier candidate should be read as historical unless explicitly
> restricted to the three-effect family.

## 1. What changed

The local Weyl-completion theorem reduces the full adaptive streamed support
function exactly to one Hilbert--Schmidt-normalised Choi MPS,

\[
\beta^{\rm stream}_{H_{\rm I},(2,2,2,2)}(\lambda)
=\max_{\operatorname{TT-rank}\leq2}
\bigl[\lambda A_{H_{\rm I}}(K)+(1-\lambda)R(K)\bigr].
\tag{1}
\]

The proof and the precise covariance assumptions are in
`notes/mps_leaf_completion.md`.  Equation (1), not the numerical maximiser
reported below, is an exact theorem.

At balanced weight an unrestricted complex Choi-MPS optimisation gives

\[
(P_{\rm A},F_{\rm R})=
(0.620085075585902\ldots,0.899520492116986\ldots),
\tag{2}
\]

and hence

\[
S_{1/2}=0.759802783851444\ldots .
\tag{3}
\]

Local Pauli completion produces a legal four-outcome instrument at every
slot, with 256 equiprobable covariance-related transcripts.  Direct
contraction reproduces (2)--(3) to $7\times10^{-15}$.  Thus (3) is a
rigorous physical lower bound, not a postselected leaf value.

## 2. Pinched cq relaxation

Put a normalised MPS in row-isometric gauge.  For every computational input
symbol $x_i$, trace the emitted carrier from its local tensor.  This gives a
binary quantum instrument on the charged qubit.  Let
$p_x=\lVert K\lvert x\rangle\rVert^2$.  Pinching $K^\dagger K$ in the computational
basis gives

\[
R(K)\leq \frac1{16}\left(\sum_x\sqrt{p_x}\right)^2.
\tag{4}
\]

The candidate saturates (4): its Gram matrix is diagonal to numerical
precision below $10^{-15}$.  Eight multistarts of the enlarged cq-instrument
problem at $\lambda=1/2$ found (3) three times, a secondary stationary point
at $0.757115893155365\ldots$ three times, and the no-record point twice.  No
larger cq value was found at balanced weight.  This is local diagnostic
evidence, not a global certificate; at larger support weight the same
relaxation finds the stronger four-effect phase.

Relaxing the first two slots still leaves (3) unchanged.  Consequently the
active geometry is already visible at the middle cut: an arbitrary
four-state qubit ensemble is followed by the two-slot suffix.

## 3. Three-effect construction

The active solution has a particularly small exact realisation.  At the
middle cut use four subnormalised pure qubit states with priors

\[
(a_0,a_\varnothing,a_+,a_-),\qquad
a_0+a_\varnothing+a_++a_-=1,
\tag{5}
\]

and a rank-one qubit POVM with three nonzero effects,

\[
G_0=tP_{n_0},\qquad
G_\varnothing=0,\qquad
G_+=gP_{n_+},\qquad
G_-=gP_{n_-},
\tag{6}
\]

where

\[
g=1-\frac t2,\qquad
n_0\!\cdot n_\pm=-\frac{t}{2-t},\qquad
t\in[0,1].
\tag{7}
\]

Equations (7) are exactly the completeness conditions for (6).  The first
signal state points along $n_0$.  The two symmetric signal states have Bloch
coordinate $r$ along $n_0$ and transverse components aligned with $n_+$ and
$n_-$.  Their correct-effect probability is

\[
h(t,r)=\frac g2\left[
1-\frac{tr}{2-t}
+\sqrt{1-\frac{t^2}{(2-t)^2}}\sqrt{1-r^2}
\right].
\tag{8}
\]

The coarse third-slot probability for these states is

\[
e(t,r)=\frac{t(1+r)}2.
\tag{9}
\]

Define the Hellinger factors

\[
c_0(t)=\sqrt t+\sqrt{1-t},\qquad
c_1(t,r)=\sqrt{e(t,r)}+\sqrt{1-e(t,r)}.
\tag{10}
\]

The null AUDIT state is chosen only to maximise recoverability.  Its factor is
therefore

\[
c_\varnothing(t)=
\begin{cases}
\sqrt2,&t\geq1/2,\\
c_0(t),&t<1/2.
\end{cases}
\tag{11}

The kink at $t=1/2$ is real.  Above it the null state can make the two coarse
outcomes exactly equiprobable; below it the largest achievable first-outcome
probability is $t$, so the state aligns with $n_0$.

For fixed $t,r,\lambda$, optimising all priors in (5) is a three-dimensional
eigenvalue problem.  Put

\[
d=(t,0,h(t,r)),\qquad
c=(c_0(t),c_\varnothing(t),\sqrt2c_1(t,r))
\tag{12}
\]

and

\[
M_\lambda(t,r)=
\lambda\operatorname{diag}d+
\frac{1-\lambda}{8}cc^{\mathsf T}.
\tag{13}
\]

If $u\geq0$ is the normalised Perron eigenvector of (13), then

\[
(a_0,a_\varnothing,a_+,a_-)
=(u_0^2,u_1^2,u_2^2/2,u_2^2/2).
\tag{14}

The exact score of this physical family is

\[
\boxed{
L_{\rm 3E}(\lambda)=
\max_{0\leq t\leq1,\,-1\leq r\leq1}
\lambda_{\max}\!\left(M_\lambda(t,r)\right).}
\tag{15}
\]

The actual lower bound also includes the no-record strategy:

\[
L_{\rm I}(\lambda)=
\max\left\{1-\frac\lambda2,L_{\rm 3E}(\lambda)\right\}.
\tag{16}

## 4. Why (15) is physically attainable

Any four pure middle-cut states form a two-site MPS with bond at most two.
The first two emitted qubits may record the corresponding input bits, so the
four columns are orthogonal.

At slot three, realise each rank-one effect $G_{ij}$ as a Kraus map whose
output memory is $\lvert j\rangle$.  Choose mutually orthogonal states of the
emitted qubit and memory for the three nonzero effects.  Because the effects
sum to the identity, these maps form a complete binary-input row isometry.
At slot four use

\[
T_{x_4}=2^{-1/2}\lvert x_4\rangle_{B_4}\otimes X^{x_4}_{M}.
\tag{17}

Measure $M$ as the second syndrome bit and always report zero for the first
syndrome bit.  The slot-four word bit is uniform, and all input columns remain
orthogonal.  Therefore (4) is saturated and (13)--(16) are exact physical
scores.  Applying the local Pauli completion then gives a complete streamed
instrument with four outcomes per slot.

## 5. Numerical frontier and phases

The coexistence equation between (15) and the no-record line gives

\[
\lambda_c=0.441437845098\ldots .
\tag{18}

At coexistence,

\[
t=0.556247377052\ldots,\qquad
r=0.028845281145\ldots .
\tag{19}

The active branch has two smooth sectors separated when $t=1/2$, at
$\lambda\simeq0.476920825$.  Below that point $c_\varnothing=\sqrt2$; above it
$c_\varnothing=c_0$.  At balanced weight,

\[
\begin{aligned}
t&=0.45807398\ldots,&r&=0.01637352\ldots,\\
(a_0,a_\varnothing,a_+,a_-)
&=(0.19921398,0.09721290,0.35178656,0.35178656),
\end{aligned}
\tag{20}

and (15) gives (2)--(3).

For $\varepsilon=1-\lambda\downarrow0$, the construction has

\[
t=\frac{\varepsilon^2}{8}+O(\varepsilon^3),\qquad
r=\frac{\varepsilon^2}{16}+O(\varepsilon^3),
\tag{21}

and

\[
L_{\rm 3E}(1-\varepsilon)
=1-\frac34\varepsilon+\frac18\varepsilon^2+O(\varepsilon^3).
\tag{22}

This approaches the exact endpoint $(P_{\rm A},F_{\rm R})=(1,1/4)$.

## 6. Exact versus conjectural statements

The following statements are proved:

1. the arbitrary adaptive problem equals the compact Choi-MPS quotient (1);
2. every normalised feasible leaf can be completed locally by Pauli/Weyl
   covariance;
3. equations (15)--(17) define a legal streamed strategy; and
4. the numerical values in (2)--(3) are rigorous lower bounds once the stored
   parameters are rounded with an explicit residual allowance.

The following historical conjecture is false:

\[
\beta^{\rm stream}_{H_{\rm I},2}(\lambda)=L_{\rm I}(\lambda)
\quad\text{for every }\lambda.
\tag{23}

At balanced weight, unrestricted Choi-MPS and pinched cq-instrument searches
agree with (15) to double precision.  They do not agree with it throughout
the interior: the physical four-effect phase is stronger beyond the observed
first-order crossing near $\lambda=0.515250589$.  Equation (1) remains an
exact variational characterisation and (15) remains an exact attainable
three-effect curve, but equation (23) is disproved.  The curated exploratory
optimisers, example commands, and representative outputs are documented in
`scratch/d2_frontier/README.md`.

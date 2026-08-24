# Schmidt-number form of the complete two-block relaxation

**Date:** 22 August 2026
**Status:** exact reformulation of the pinched two-block outer programme;
global optimisation over its terminal variables remains open

## 1. One rank-two vector replaces the state--instrument product

Vectorise a normalised two-block leaf across its middle cut as

\[
|K\rangle\!\rangle\in
\mathcal H_L\otimes\mathcal H_R,
\qquad
\dim\mathcal H_L=16,\quad
\dim\mathcal H_R=32.
\tag{1}
\]

Here \(L=(B_L,z)\) and \(R=(B_R,m,y)\).  The one-qubit bond condition is
exactly

\[
\operatorname{SR}(|K\rangle\!\rangle)\leq2,
\qquad
\langle\!\langle K|K\rangle\!\rangle=1.
\tag{2}
\]

For a terminal POVM \(P=\{P_s\}_{s\in\mathbb F_2^2}\), let
\(D(P)\) be the positive block-diagonal operator that applies
\(P_{z\oplus y}\) to the terminal-memory factor in path \((z,y)\) and acts
as the identity on both emitted registers.  Up to the harmless transpose
fixed by vectorisation convention,

\[
A_P(K)=\langle\!\langle K|D(P)|K\rangle\!\rangle.
\tag{3}
\]

## 2. Pinched RETURN is a rank-one perturbation

Let \(\Pi_i\), \(i=1,\ldots,16\), project onto one input-path block
\((z,y)\), and define

\[
\mathcal V=
\left\{|v\rangle:
\Pi_i|v\rangle=|v_i\rangle,
\ \lVert v_i\rVert_2=\frac14\text{ for every }i
\right\}.
\tag{4}
\]

Blockwise Cauchy--Schwarz gives the exact identity

\[
\frac1{16}
\left(\sum_i\lVert\Pi_iK\rVert_2\right)^2
=\max_{v\in\mathcal V}
|\langle v|K\rangle\!\rangle|^2.
\tag{5}
\]

The maximiser simply aligns each \(v_i\) with \(\Pi_iK\).  Equation (5) is
the computational-pinching RETURN term, not an additional relaxation.

## 3. Exact \(S(2)\)-norm programme

For a positive bipartite operator \(X\), define

\[
\lVert X\rVert_{S(2)}
:=\max_{\substack{\lVert\psi\rVert=1\\
\operatorname{SR}(\psi)\leq2}}
\langle\psi|X|\psi\rangle.
\tag{6}
\]

Maximising equations (3) and (5) over the same leaf and commuting the compact
maxima yields

\[
\boxed{
\overline\beta_{\rm Choi}(\lambda)
=\max_{\substack{P_s\succeq0,\ \sum_sP_s=I_2\\v\in\mathcal V}}
\left\lVert
\lambda D(P)+(1-\lambda)|v\rangle\!\langle v|
\right\rVert_{S(2)}.}
\tag{7}
\]

Thus the common state--instrument compatibility in the Choi programme is
equivalently a Schmidt-number-two constraint on one positive operator norm.
This does not by itself make the problem easy: computing an \(S(k)\) norm is
equivalent to testing \(k\)-block positivity in general.  The relevant norm
and SDP upper bounds were developed by Johnston and Kribs,
[arXiv:0909.3907](https://arxiv.org/abs/0909.3907) and
[arXiv:1006.0898](https://arxiv.org/abs/1006.0898).  Complete hierarchies for
rank-constrained semidefinite optimisation give a second route; see Yu *et
al.*, [arXiv:2012.00554](https://arxiv.org/abs/2012.00554).

Equation (7) is nevertheless useful.  It replaces the bilinear product of
four prefix states and four Choi matrices by a standard, named convex cone:
the cone of states of Schmidt number at most two.  A dual certificate at
level \(L\) is a proof that

\[
L I-\lambda D(P)-(1-\lambda)|v\rangle\!\langle v|
\tag{8}
\]

is two-block-positive for every terminal \(P\) and every equal-block-norm
\(v\).  The projective certificate already establishes this after restricting
\(P\) to two complementary rank-one outcomes.  The remaining interior gate
is equation (8) for extremal three- and four-effect qubit POVMs.

There is also an exact convexification that does not reward classical path
mixing. `notes/coherence_preserving_convexification.md` replaces the
pure-state Hellinger expression by the sum of trace norms of the off-diagonal
path blocks. Maximising that convex block-coherence extension over the mixed
Schmidt-number-two cone is exactly equal to the present pure-state programme.

## 4. Scope and numerical check

`scratch/d2_frontier/verify_schmidt_norm_reduction.py` reconstructs the
four-effect checkpoint, builds its aligned \(v\), and verifies independently
that its Schmidt rank is two, all sixteen blocks of \(v\) have norm \(1/4\),
and the two expectation values in equations (3) and (5) reproduce AUDIT and
pinched RETURN.

This equivalence identifies the correct certification technology; it is not
itself a certificate of equation (8).  In particular, a low-level SDP value
without an archived dual witness must not be reported as a completed
unrestricted frontier.

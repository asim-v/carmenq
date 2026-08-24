# Coherence-preserving convexification of pinched RETURN

**Date:** 23 August 2026
**Status:** exact reformulation theorem; certifying its global maximum remains
the open interior problem

## 1. Why the obvious convexification fails

Let \(\{\Pi_i\}_{i=1}^N\) be the orthogonal projectors onto the path blocks
of a normalized pure leaf \(|K\rangle\!\rangle\). With

\[
p_i=\lVert\Pi_iK\rVert_2^2,
\]

pinched RETURN is

\[
R(K)=\frac1N\left(\sum_i\sqrt{p_i}\right)^2.
\tag{1}
\]

Applying equation (1) to the diagonal of a mixed Schmidt-number relaxation
is not a harmless extension: its concavity rewards classical path mixing.
The exact \(0.7902656097\ldots\) barrier derived in
[`convexification_barrier.md`](convexification_barrier.md) is attained by a
separable path-diagonal state, so no tighter description of the mixed
Schmidt-number cone can remove it.

## 2. The block-coherence extension

For an arbitrary positive operator \(\rho\), define

\[
\mathcal R_{\rm coh}(\rho)
=\frac1N\left[
\operatorname{Tr}\rho
+\sum_{i\ne j}\lVert\Pi_i\rho\Pi_j\rVert_1
\right].
\tag{2}
\]

This is normalized block \(l_1\)-coherence. Ordinary \(l_1\)-coherence was
introduced as a resource measure by Baumgratz, Cramer, and Plenio
([arXiv:1311.0275](https://arxiv.org/abs/1311.0275)); block coherence and its
free operations are developed by Dey *et al.*
([arXiv:1908.01882](https://arxiv.org/abs/1908.01882)). The established
resource measure is not claimed as new. The problem-specific contribution is
to use it as the exact mixed-state extension of this temporal RETURN score.

## 3. Exact convexification theorem

> **Theorem 1 (coherence-preserving convexification).** Let
> \(\mathcal P_k\) be the normalized pure states of Schmidt rank at most
> \(k\) across the middle cut, and let
> \[
> \mathcal S_k=\operatorname{conv}
> \{|\psi\rangle\!\langle\psi|:|\psi\rangle\in\mathcal P_k\}
> \]
> be the states of Schmidt number at most \(k\). For every Hermitian audit
> operator \(D\) and \(0\leq\lambda\leq1\),
> \[
> \max_{|\psi\rangle\in\mathcal P_k}
> \left[\lambda\langle\psi|D|\psi\rangle
> +(1-\lambda)R(\psi)\right]
> =
> \max_{\rho\in\mathcal S_k}
> \left[\lambda\operatorname{Tr}(D\rho)
> +(1-\lambda)\mathcal R_{\rm coh}(\rho)\right].
> \tag{3}
> \]

For a pure state, the \((i,j)\) block is

\[
\Pi_i|\psi\rangle\!\langle\psi|\Pi_j
=|\psi_i\rangle\!\langle\psi_j|,
\]

whose trace norm is
\(\lVert\psi_i\rVert_2\lVert\psi_j\rVert_2\). Consequently,

\[
\mathcal R_{\rm coh}(|\psi\rangle\!\langle\psi|)
=\frac1N\left(\sum_i\lVert\psi_i\rVert_2\right)^2
=R(\psi).
\tag{4}
\]

Equation (2) is convex because every trace norm is convex. If
\(\rho=\sum_aq_a|\psi_a\rangle\!\langle\psi_a|\in\mathcal S_k\), then

\[
\begin{split}
&\lambda\operatorname{Tr}(D\rho)
+(1-\lambda)\mathcal R_{\rm coh}(\rho)\\
&\quad\leq\sum_aq_a\left[
\lambda\langle\psi_a|D|\psi_a\rangle
+(1-\lambda)R(\psi_a)\right],
\end{split}
\]

which is no greater than the pure-state maximum. The reverse inequality
follows because every member of \(\mathcal P_k\) is feasible in
\(\mathcal S_k\). This proves equation (3).

The theorem is elementary but decisive: purity can be convexified without
changing the optimum, provided the off-diagonal path coherences are retained.
A path-diagonal classical mixture now receives only
\(\mathcal R_{\rm coh}=1/N\), not the spurious Hellinger bonus.

## 4. Exact polar-witness program

For each unordered pair \(i<j\), let \(U_{ij}:\mathcal H_j\to\mathcal H_i\)
be a contraction and form the Hermitian block operator

\[
W(U)_{ij}=U_{ij},\qquad
W(U)_{ji}=U_{ij}^{\dagger},\qquad
W(U)_{ii}=0.
\tag{5}
\]

The trace-norm variational identity gives

\[
\mathcal R_{\rm coh}(\rho)
=\max_{\lVert U_{ij}\rVert_\infty\leq1}
\frac1N\operatorname{Tr}[(I+W(U))\rho].
\tag{6}
\]

For a Hermitian operator \(X\), write

\[
h_k(X)=\max_{\substack{\lVert\psi\rVert=1\\
\operatorname{SR}(\psi)\leq k}}
\langle\psi|X|\psi\rangle.
\tag{7}
\]

Compactness lets the two maxima commute. Combining equations (3) and (6)
therefore gives a second exact form of the frontier:

\[
\boxed{
\beta_k(\lambda,D)
=\max_{\lVert U_{ij}\rVert_\infty\leq1}
h_k\!\left(
\lambda D+\frac{1-\lambda}{N}[I+W(U)]
\right).}
\tag{8}
\]

For positive arguments, \(h_k\) is the usual \(S(k)\)-norm; in general it
is its one-sided Hermitian support version. Equation (8) is not a numerical
upper bound by itself. It is an exact certification target: a spatial cover
now ranges over bounded contraction matrices, while every fixed cell reduces
to a linear Schmidt-number support problem. It is equivalent to the aligned
rank-one program in [`two_block_schmidt_norm.md`](two_block_schmidt_norm.md),
but exposes pairwise path coherence rather than square-root probabilities.

## 5. Reproducible checks and remaining gate

[`coherence_polar_program.py`](../scratch/d2_frontier/coherence_polar_program.py)
constructs the optimal polar contractions for arbitrary unequal path-block
sizes. The tests verify that all blocks have operator norm at most one, the
witness is Hermitian, and equation (6) agrees with equation (2) on random
mixed states. A separate test confirms that a path-diagonal mixture receives
exactly \(1/N\).

This closes a formulation error, not the target inequality. The unresolved
task is to certify equation (8) below the desired target for every admissible
terminal POVM and every contraction family. A convergent route must either
cover those contractions with rigorous remainder bounds or derive a
problem-specific two-block-positive majorant of \(I+W(U)\).

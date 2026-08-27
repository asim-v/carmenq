# Complete finite Choi programme for the interior two-block frontier

**Date:** 27 August 2026
**Status:** exact finite outer formulation plus a solver-independent finite
replay enclosure at \(\lambda=0.6\); exact equality, a fully kernel-formalised
physical reduction, and the complete support curve remain open

## 1. The common-instrument constraint

After the middle cut, every two-block rank-two leaf can be written with four
subnormalised prefix states

\[
\rho_z\succeq0,\qquad
\sum_{z\in\mathbb F_2^2}\operatorname{Tr}\rho_z=1,
\tag{1}
\]

and one four-outcome qubit instrument
\(\{\Phi_y:y\in\mathbb F_2^2\}\).  In input-major Choi convention let

\[
J_y\succeq0,\qquad
\sum_y\operatorname{Tr}_{\rm out}J_y=I_2.
\tag{2}
\]

The sixteen conditioned terminal states and their path probabilities are

\[
\sigma_{zy}
=\operatorname{Tr}_{\rm in}
 \left[(\rho_z^{\mathsf T}\otimes I_2)J_y\right],
\qquad
p_{zy}=\operatorname{Tr}\sigma_{zy}.
\tag{3}
\]

The syndrome ensemble seen by the terminal qubit is

\[
\tau_s=\sum_{z\oplus y=s}\sigma_{zy}.
\tag{4}
\]

Equations (1)--(4) preserve the constraint that all sixteen paths arise from
the same Stinespring-compatible continuation.  Replacing the four maps by
unrelated pulled-back effects loses this constraint and gives the much looser
static value \(0.8605551275\ldots\) at \(\lambda=0.6\).

Conversely, every feasible collection in equations (1)--(2) has a two-block
realisation.  Purify each \(\rho_z\) into the left emitted register and use a
Kraus representation of each \(J_y\) in the four-dimensional right emitted
register.  Thus no additional adaptive-tree variables are missing from this
description.

## 2. Exact AUDIT and the RETURN outer term

The terminal POVM can be eliminated exactly:

\[
A(J,\rho)
=\max_{\{P_s\}}\sum_s\operatorname{Tr}(P_s\tau_s)
=\min_{Y\succeq\tau_s}\operatorname{Tr}Y.
\tag{5}
\]

For qubits, writing
\(\tau_s=(\pi_sI+r_s\cdot\sigma)/2\) turns equation (5) into

\[
A(J,\rho)
=\min_{c\in\mathbb R^3}\max_s
 \left(\pi_s+\lVert c-r_s\rVert_2\right).
\tag{6}
\]

Computational pinching of the input Gram matrix gives

\[
R(K)\leq R_{\rm pin}(J,\rho)
:=\frac1{16}
\left(\sum_{z,y}\sqrt{p_{zy}}\right)^2.
\tag{7}
\]

Therefore the exact two-block support is bounded by the finite programme

\[
\boxed{
\beta_{2\rm b}(\lambda)
\leq
\overline\beta_{\rm Choi}(\lambda)
:=\max_{\rho,J}
\left[\lambda A(J,\rho)
+(1-\lambda)R_{\rm pin}(J,\rho)\right]
}
\tag{8}
\]

subject only to equations (1)--(3).  The inequality in (8) is solely the
pinching step.  It is an equality for the explicit 3E and 4E leaves because
their sixteen input columns are orthogonal.

The four-effect checkpoint at \(\lambda=0.6\) evaluates in this formulation
to

\[
(A,R_{\rm pin},L)
=(0.869930022902\ldots,
  0.609852003743\ldots,
  0.765898815238\ldots),
\tag{9}
\]

within \(3\times10^{-11}\) of the independent physical contraction.

## 3. Why this is the right finite certification target

For fixed prefix states and terminal POVM, equation (8) is a conic
optimisation over the four Choi matrices.  For fixed instrument and POVM it
is a conic optimisation over the four prefix states.  For fixed states and
instrument, equation (5) is an SDP over the terminal POVM.  The Hellinger
term has an exact second-order-cone hypograph: for every two paths introduce
\(h_{ij}\geq0\) and impose

\[
\left\lVert(2h_{ij},p_i-p_j)\right\rVert_2
\leq p_i+p_j.
\tag{10}
\]

Since \(\sum_i p_i=1\), maximising

\[
\frac{1+2\sum_{i<j}h_{ij}}{16}
\tag{11}
\]

recovers equation (7) exactly.  Hence each block update is globally solved;
only the coupling among the three blocks is nonconvex.

Two standard extremality reductions keep the remaining hierarchy finite.
For a fixed leaf an optimal terminal POVM can be chosen extremal, so on a
qubit it has at most four rank-one effects.  For fixed coarse effects
\(Q_y=\Phi_y^*(I)\), prefix states, and terminal POVM, RETURN is independent
of the post-filter channels and AUDIT is linear in them.  An extreme qubit
channel suffices for each outcome; the usual linear-independence criterion
then bounds its Kraus rank by two.  These observations do not prove terminal
projectivity, but they remove unbounded POVM and Kraus alphabets from the
interior problem.  General extremal-instrument criteria are given by
D'Ariano, Perinotti, and Sedlak,
[arXiv:1101.4889](https://arxiv.org/abs/1101.4889).

## 4. Exhaustive arity and weight split

The weighted-ball KKT conditions split every extremal terminal readout into
two-, three-, or four-active sectors. The projective sector is covered in all
four rank/topology classes by outward-expanded interval arithmetic:

\[
\beta_{\rm projective}(0.6)\leq0.76600.
\tag{12}
\]

For the remaining sectors, the weight of the largest effect gives the first
exact partition. If \(w_{\max}\leq0.88325\), an exact-residual replay of 576
capped-weight SOCP cells yields \(0.765893818\). Above that threshold, the
three-active sector is covered by 12,008 transferred ternary cells at
\(0.76652\). The four-active sector is split once more at
\(w_{\min}=0.0003\): 90 spatial leaves, replayed in all six affine orders,
give \(0.76670\) when the smallest effect is retained. When it is smaller,
deleting it reduces AUDIT by at most \(0.0003\), hence

\[
0.76652+\frac35(0.0003)=0.76670.
\tag{13}
\]

A rational four-effect physical strategy supplies the other side. Five
rational half-angle coordinates make every state, effect, and prior exact;
192-bit dyadic square-root floors bound the RETURN radicals from below.
Consequently,

\[
\boxed{
\frac{957373519}{1250000000}
\leq\beta^{\rm stream}_{H_{\rm I},2}(0.6)
\leq\overline\beta_{2\rm b}(0.6)
\leq\frac{7667}{10000}}
\tag{14}
\]

or, in outward decimals,
\([0.7658988152,0.76670]\). The width is \(0.0008011848\).

The upper replay calls no optimiser and does not trust the optimisers used to
propose dual vectors. Its finite proof boundary consists of directed
binary64 interval kernels for the projective cells, exact rational residual
repair for the conic cells, and exact final assembly. Python matrix
canonicalisation and the analytic reductions from the physical problem to
these finite cells remain trusted; the result is therefore solver independent
in verification, but not an end-to-end kernel-formalised theorem. The full
sector ledger is in `notes/interleaved_interior_frontier_l060.md`.

## 5. Reproducible block-convex evidence

`scratch/d2_frontier/two_block_choi_seesaw.py` implements equations
(1)--(11) directly in complex CVXPY.  It independently reconstructs a
two-block tensor checkpoint, solves every block with CLARABEL, and classifies
the final ensemble with the weighted-ball SOCP.

The corrected random rank-one POVM probe gives the same qualitative split.
At \(\lambda=0.6\), nine of twelve deterministic seeds remained four-active
and converged to the reversible point
\((A,R_{\rm pin})=(1/2,1)\), of support \(0.7\).  The other three became
two-active and converged to 3E or 4E; the largest value was
\(0.765898815264695\ldots\).  No three-active basin appeared.  This is useful
falsification evidence, but it is not equation (13).

For other support directions, a natural rigorous continuation is a
dimension-constrained moment/SOS relaxation of equations (1)--(11), split by
the three- and four-active KKT systems. Finite-dimensional moment hierarchies
are established tools rather than a new claim of this project; see Navascues
and Vertesi, [arXiv:1412.0924](https://arxiv.org/abs/1412.0924), and Navascues
*et al.*, [arXiv:1507.07521](https://arxiv.org/abs/1507.07521). A solver value
is not a certificate until its dual witness, residuals, numerical precision,
and a solver-independent replay are archived.

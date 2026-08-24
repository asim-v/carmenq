# Planar instrument completion through numerical radius

**Date:** 23 August 2026

**Status:** exact specialization of established operator-system and numerical-
radius results; new benchmark-specific reduction and incompatibility witness,
not yet a closed global frontier bound

## 1. The missing compatibility condition

Fix a nondegenerate rank-one qubit POVM with three active effects in one
Bloch plane,

\[
P_s=\frac{w_s}{2}(I+n_{sx}X+n_{sy}Y),\qquad s=0,1,2.
\tag{1}
\]

For one outcome \(y\) of a putative instrument, let

\[
H_{ys}=\Phi_y^*(P_s),\qquad F_y=\sum_sH_{ys}=\Phi_y^*(I).
\tag{2}
\]

The three effects in equation (1) form a basis of the planar operator system
\(\mathcal S=\operatorname{span}\{I,X,Y\}\). Thus the observed pullbacks
uniquely determine

\[
F_y=\Phi_y^*(I),\qquad
B_{yx}=\Phi_y^*(X),\qquad
B_{yy}=\Phi_y^*(Y).
\tag{3}
\]

Positivity of every \(H_{ys}\), the bounds \(H_{ys}\preceq w_sF_y\), and
normalization \(\sum_yF_y=I\) are necessary, but they do not imply that the
three columns arise from one completely positive map. The missing question is
whether a Hermitian fourth pullback \(B_{yz}=\Phi_y^*(Z)\) exists.

## 2. Exact Choi completion

The Choi matrix of the adjoint map has the block form

\[
C_y(B_{yz})=\\frac12
\begin{bmatrix}
F_y+B_{yz} & B_{yx}+iB_{yy}\\
B_{yx}-iB_{yy} & F_y-B_{yz}
\end{bmatrix}.
\tag{4}
\]

Consequently, the planar data have a completely positive extension if and
only if some Hermitian \(B_{yz}\) makes equation (4) positive semidefinite.
This is the finite-dimensional operator-system extension problem. General
semidefinite formulations for partially specified quantum operations are
given by Heinosaari, Jivulescu, Reeb, and Wolf
([arXiv:1205.0641](https://arxiv.org/abs/1205.0641)).

When \(F_y\) is positive definite, congruence by
\(\operatorname{diag}(F_y^{-1/2},F_y^{-1/2})\) turns equation (4) into an
Ando block completion. Define

\[
T_y=F_y^{-1/2}(B_{yx}+iB_{yy})F_y^{-1/2}.
\tag{5}
\]

Ando's theorem gives the exact scalar criterion

\[
\boxed{\quad C_y(B_{yz})\succeq0\text{ for some }B_{yz}
\quad\Longleftrightarrow\quad w(T_y)\leq1,\quad}
\tag{6}
\]

where \(w(T)=\max_{\lVert v\rVert=1}|v^\dagger Tv|\) is the numerical
radius. For singular \(F_y\), the same statement holds on its support,
together with the necessary support inclusion for \(B_{yx}\) and
\(B_{yy}\).

To see equation (6), set

\[
X=\frac12F_y^{-1/2}(F_y+B_{yz})F_y^{-1/2}.
\]

The normalized Choi matrix is positive precisely when

\[
\begin{bmatrix}
X&T_y/2\\T_y^\dagger/2&I-X
\end{bmatrix}\succeq0
\]

for some \(0\preceq X\preceq I\). This is Ando's positive-block
characterization of the numerical-radius unit ball. The theorem originates
in Ando's 1973 paper
[*Structure of operators with numerical radius one*](https://acta.bibl.u-szeged.hu/14381/1/math_034_011-015.pdf);
Bhatia and Jain give a modern block-matrix treatment and explicitly discuss
its completely positive map interpretation
([arXiv:2307.16014](https://arxiv.org/abs/2307.16014)).

## 3. Exact instrument statement

> **Corollary.** A family \(\{H_{ys}\}_{y,s}\) for the fixed planar POVM is
> produced by one qubit instrument \(\{\Phi_y\}_y\) if and only if:
>
> 1. every outcome family has the support property and \(w(T_y)\leq1\); and
> 2. \(\sum_yF_y=I\).

Necessity follows from the positive Choi matrices of the physical
subchannels. For sufficiency, equation (6) supplies one positive Choi
completion per \(y\); the second condition makes their sum trace preserving.
The resulting completely positive maps therefore form a quantum instrument
and reproduce all observed pullbacks in equation (2).

This corollary is the precise way to require that all conditioned columns
come from one common instrument without retaining unobserved output-channel
coordinates. It also identifies the obstruction with one numerical-radius
number per instrument outcome.

## 4. Explicit failure of columnwise positivity

The archived outer strategy satisfies, to numerical precision,

\[
H_{ys}\succeq0,\qquad H_{ys}\preceq w_sF_y,qquad \sum_yF_y=I.
\]

Nevertheless, its four Ando radii are

\[
(0.9999999994,\ 1.7244406626,\ 1.5335799867,\ 0.9999999996).
\tag{7}
\]

The second and third outcomes therefore have no common CP completion. This is
a direct multicolumn incompatibility witness, not a failure of an optimizer
to find hidden Choi variables. More explicitly, the radius maximizer supplies
a pure input vector \(v_y\) for which

\[
\sqrt{
  \langle B_{yx}\rangle_{v_y}^{\,2}
  +\langle B_{yy}\rangle_{v_y}^{\,2}}
>
\langle F_y\rangle_{v_y}.
\tag{8}
\]

Equation (8) says that the conditioned planar output would lie outside its
own positive Bloch disk. It is a scalar separating certificate that can be
checked directly from the archived matrices, without constructing a missing
Choi coordinate or trusting a nonlinear solve.

[`planar_cp_completion.py`](../scratch/d2_frontier/planar_cp_completion.py)
reconstructs equations (3) and (5), evaluates numerical radius through its
Hermitian support formula, and emits the four-outcome report. Unit tests check
identity and random Kraus maps, known numerical-radius examples, and an
explicit family that passes every individual effect inequality while
violating equation (6).

The reproducibility bundle contains the
[`outer checkpoint`](../scratch/d2_frontier/joint_effect_outer_092_064_044_l055_4r.npz),
its [`radius report`](../scratch/d2_frontier/planar_cp_completion_outer_l055.json),
the [`sampled outer run`](../scratch/d2_frontier/joint_effect_pathpos_ando64_092_064_044_l055_30s.json),
and the [`exact-completion run`](../scratch/d2_frontier/joint_effect_cpcomplete_linked_092_064_044_l055_60s.json).
The scalar witness report is regenerated with

```text
python scratch/d2_frontier/planar_cp_completion.py scratch/d2_frontier/joint_effect_outer_092_064_044_l055_4r.npz
```

## 5. Novelty boundary and next use

The equivalence between numerical-radius contraction and positive block
completion is established mathematics, as is completely positive extension
from operator systems. The specialization to a planar POVM is therefore not
claimed as a new standalone theorem. The contribution candidate is narrower:
using this equivalence as the exact common-instrument closure for the temporal
AUDIT--RETURN frontier, together with a concrete incompatibility witness and
a reduced global optimization model. A targeted literature search found the
abstract CP-extension and numerical-radius machinery, but not this
planar-pullback use in a temporal AUDIT--RETURN frontier; that negative search
is evidence for positioning, not a proof of priority.

The reduced SCIP formulation reconstructs only the three observed pullbacks,
adds one missing \(Z\) operator per outcome, and links that same operator to
the terminal Helstrom certificate. A 30-second outer pilot using 64 sampled
forms of equation (8), together with positivity of every reconstructed path
output, returned the solver-conditional dual \(0.7699011264\). This remains
above \(0.758\), and finite sampling is only a necessary outer closure. The
exact Choi-complete model is therefore retained as the proof target. An exact
60-second linked run returned the solver-conditional dual $0.7699427861$
without a feasible target incumbent; it is a progress record,
not an infeasibility certificate. The next analytic opportunity is to combine
the radius constraint with the Hellinger RETURN term before spatial branching,
rather than treating the Choi completion and state--effect products as
unrelated quadratics.

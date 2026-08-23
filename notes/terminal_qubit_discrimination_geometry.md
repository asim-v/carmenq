# Terminal-qubit discrimination as a weighted smallest-ball problem

**Date:** 22 August 2026
**Status:** exact reduction; the three- and four-active sectors are now
numerically covered at \(\lambda=0.6\), while the complete support curve and
solver-independent validation remain open

## 1. Exact elimination of the terminal POVM

Let the four subnormalised terminal-memory states associated with the four
syndrome labels be

\[
\tau_s=\frac12\bigl(p_s I+r_s\mathbin{\cdot}\sigma\bigr),
\qquad s\in\mathbb F_2^2,
\]

where \(p_s=\operatorname{Tr}\tau_s\), \(r_s\in\mathbb R^3\), and
\(\lVert r_s\rVert_2\le p_s\).  Minimum-error AUDIT is the semidefinite
program

\[
P_{\rm guess}=\max_{\{E_s\}}\sum_s\operatorname{Tr}(E_s\tau_s)
=\min_{Y\succeq\tau_s}\operatorname{Tr}Y.
\tag{1}
\]

Write the dual variable as

\[
Y=\frac12\bigl(tI+y\mathbin{\cdot}\sigma\bigr).
\]

For a two-by-two Hermitian matrix, positivity is equivalent to its scalar
part dominating the Euclidean norm of its Bloch part.  Therefore

\[
Y\succeq\tau_s
\quad\Longleftrightarrow\quad
t-p_s\geq\lVert y-r_s\rVert_2,
\]

and equation (1) becomes the additively weighted Euclidean one-centre problem

\[
\boxed{
P_{\rm guess}
=\min_{y\in\mathbb R^3}\max_s
\left(p_s+\lVert y-r_s\rVert_2\right).}
\tag{2}
\]

Thus the terminal POVM can be removed exactly and replaced by a
second-order-cone problem.  This is the qubit minimum-error discrimination
geometry in a form adapted to the temporal frontier.  It is consistent with
the general geometric solution of qubit discrimination; see Bae and Hwang,
[arXiv:1204.2313](https://arxiv.org/abs/1204.2313).

## 2. Active constraints classify the missing sectors

Let \((t_\star,y_\star)\) solve equation (2).  A label is active when

\[
t_\star-p_s=\lVert y_\star-r_s\rVert_2.
\tag{3}
\]

The Karush--Kuhn--Tucker conditions give nonnegative weights \(\mu_s\),
supported on active labels, such that

\[
\sum_s\mu_s=1,
\qquad
\sum_s\mu_s
\frac{r_s-y_\star}{\lVert r_s-y_\star\rVert_2}=0.
\tag{4}
\]

Complementary slackness reconstructs an optimal POVM effect as

\[
E_s=\mu_s\left(I+n_s\mathbin{\cdot}\sigma\right),
\qquad
n_s=\frac{r_s-y_\star}{\lVert r_s-y_\star\rVert_2}.
\tag{5}
\]

Caratheodory's theorem limits a nondegenerate active set in three-dimensional
Bloch space to four labels.  The cases are now geometrically distinct:

* two active labels force antipodal \(n_s\) and hence a binary projective
  POVM;
* three active labels give a genuine three-effect qubit POVM; and
* four active labels give a genuine four-effect qubit POVM.

Consequently, the unproved terminal-projectivity claim is equivalent to a
precise statement about the weighted smallest ball:

> Among support-maximising two-block rank-two leaves, one can choose a
> terminal syndrome ensemble whose weighted one-centre has at most two active
> constraints.

This formulation removes ambiguity about what still has to be shown.  It also
prevents ordinary POVM extremality from being mistaken for a proof: extremal
qubit POVMs with three or four rank-one effects correspond exactly to the
three- and four-active cases.

## 3. The exposed four-effect leaf is strictly two-active

`scratch/d2_frontier/analyze_two_block_leaf.py` now solves equation (2) with
an independent SOCP after reconstructing the terminal ensemble from the
unrestricted complex checkpoint.  At \(\lambda=0.6\) it returns

\[
P_{\rm guess}=0.869930022901752\ldots,
\]

in agreement with the directly contracted AUDIT score.  Only labels zero and
three are active.  The two inactive weighted-ball constraints have slacks

\[
0.7821363837\ldots,
\qquad
0.7813522117\ldots .
\]

The KKT weights are \(1/2,0,0,1/2\) to numerical precision.  This is stronger
than observing two tiny parametrised POVM effects: it proves that the exact
optimal discrimination problem for the reconstructed ensemble lies well
inside the two-active region.  Small perturbations of this leaf therefore do
not open a genuine nonprojective readout sector.

## 4. A tempting dimension bound is valid but too loose

For arbitrary syndrome priors \(p_s\), discrimination of freely chosen states
in dimension two obeys

\[
P_{\rm guess}\leq p_{(1)}+p_{(2)},
\tag{6}
\]

where the right side is the sum of the two largest priors.  Indeed,
\(\tau_s\preceq p_sI\), so

\[
\sum_s\operatorname{Tr}(E_s\tau_s)
\leq\sum_sp_s\operatorname{Tr}E_s,
\qquad
\sum_s\operatorname{Tr}E_s=2,
\]

and the resulting linear programme assigns unit trace to the two largest
weights.  A binary projective encoding attains the bound when terminal states
may be chosen freely.

This does **not** close the temporal problem.  Optimising equation (6) together
with the exact central-qubit path probabilities but dropping the common
Stinespring channel gives

\[
\lambda P_{\rm guess}+(1-\lambda)R
\leq 0.860555127546400\ldots
\quad(\lambda=0.6),
\]

well above the physical value \(0.7658988153\ldots\).  The same number is
obtained by the independent compatible-effect relaxation.  The missing
physics is therefore not terminal dimension alone.  It is the requirement
that all sixteen conditioned outputs arise from the same four continuation
maps \(L_y\), or equivalently from one common Stinespring-compatible channel
geometry.

## 5. Exact boundary after the projective certificate

The binary-projective sector has a complete numerical spatial
branch-and-bound cover at \(\lambda=0.6\):

\[
0.765898815264694\ldots
\leq\beta_{\rm projective}(0.6)
\leq0.76591.
\tag{7}
\]

The lower value is a verified four-slot physical leaf.  The upper value is a
finite cover of all four extreme split topologies, validated by
`scratch/d2_frontier/validate_projective_cover.py`; its width is
\(1.1185\times10^{-5}\).  This is a SCIP numerical certificate subject to the
recorded solver tolerances, not an interval proof independent of the solver.

After equation (7), the unrestricted converse consisted only of weighted-ball
solutions with three or four active constraints.  The earlier trace-barrier
interpretation was not evidence for such a solution: once the terminal POVM
is eliminated exactly, its intermediate ensembles can already be two-active
even when all four parametrised effect traces are positive.

That logical gap is now covered directly at \(\lambda=0.6\), without proving
or assuming terminal projectivity.  Fully active complementary slackness
reconstructs

\[
\tau_i=Y-(A-p_i)\Pi_i^\perp,
\qquad
p_i\geq A\frac{1+2tx_i+t^2}{2(1+tx_i)},
\]

while POVM closure gives both a scalar longitudinal constraint and an exact
transverse polygon criterion.  Replacing the readout by aligned binary
projectors also yields a geometry-free averaged comparison with the certified
projective support lines.  A complete ternary SOCP cover, a four-active
spatial SCIP dual, and deletion of an effect below trace \(0.0003\) combine to
give

\[
0.7658988152646944
\leq\beta_{2\mathrm b}(0.6)
\leq0.76662.
\]

The result is a complete finite solver-conditional numerical enclosure, not a
formal interval proof or an equality theorem for the exposed 4E leaf.  The
full derivation and certificate audit are in
`notes/interleaved_interior_frontier_l060.md`.  Extending the enclosure across
all support directions, or replacing its numerical duals by analytic or
outward-rounded certificates, is the remaining theorem gate.

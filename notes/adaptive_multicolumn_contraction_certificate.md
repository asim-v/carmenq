# Adaptive multicolumn contraction certificate

## Result and scope

The most permissive cell in the existing Fourier/pair outer cover of the
fixed \(\lambda=0.55\) interior benchmark is now excluded at score \(0.758\).
The calculation uses 48 adaptive contractions, 98 exhaustive spectral
branches at every expansion, and no unclosed leaf.  The resulting tree has
4,657 terminal leaves: 2,640 are conically infeasible and the largest finite
terminal upper bound is \(0.7579901084\).  This is a
solver-conditional certificate for one cell, not yet a closure of the full
regular interior.  In particular, the coefficient directions selected in
this tree cannot silently be imposed on the other open cells without giving
those cells their own exhaustive trace-norm branch covers.

The base cell is planar Fourier sector 4 of 8, spherical cube-face cell 18 at
grid 4, and Bloch cap 3 for the \((2,3)\) pair contraction.  Before the new
cuts its conic upper bound is \(0.7633263647\), the largest value among the
768 cells in the earlier Fourier/pair calculation.  The adaptive tree has
maximum depth 3.  Its complete branch table, parent links, solver statuses,
and numerical upper bounds are stored in
`scratch/d2_frontier/adaptive_multicolumn_worstcell_l055_auditable.json`.

## The common-instrument inequality

Let \(A_z\) be any Hermitian input family and let
\(\{\Phi_y\}_y\) be one quantum instrument.  For every real coefficient
vector \(c\), define

\[
 X_c=\sum_z c_z A_z,
 \qquad
 X_{c,y}=\sum_z c_z\Phi_y(A_z).
\]

The flagged map \(\Gamma(X)=\bigoplus_y\Phi_y(X)\) is a channel.  Trace-norm
data processing therefore gives the necessary condition

\[
  \sum_y\lVert X_{c,y}\rVert_1
  =\lVert\Gamma(X_c)\rVert_1
  \leq\lVert X_c\rVert_1. \tag{1}
\]

Unlike the character-only Fourier cuts, (1) holds for every real linear
combination of the four columns.  It is precisely where the requirement that
all conditioned states arise from a *single* instrument enters.  The
inequality itself is standard trace-distance data processing; the substantive
work here is its adaptive use inside a finite, reproducible outer
certificate, not a claim that data processing is new.

For a Hermitian qubit operator
\(X=(dI+\mathbf r\cdot\boldsymbol\sigma)/2\),

\[
 \lVert X\rVert_1=\max\{|d|,\lVert\mathbf r\rVert_2\}.
\]

Consequently the right side of (1) has two scalar-dominated branches and one
vector-dominated branch.  The latter is covered by 96 normalized cube-face
caps.  If a cap with centre \(\mathbf n\) and covering cosine \(\gamma\)
contains the direction of \(\mathbf r\), then

\[
 \lVert\mathbf r\rVert_2
 \leq \frac{\mathbf n\cdot\mathbf r}{\gamma}.
\]

Thus the scalar-positive, scalar-negative, and 96 Bloch-cap problems form an
exhaustive 98-way disjunction.  Each child is a convex conic outer problem.
The cap inequality is deliberately one-sided, so every child remains an
upper relaxation: closing all children is valid, while leaving one open says
nothing about physical attainability.

## Adaptive separation

At an open node, the conic maximizer supplies four prefix states and sixteen
conditioned output blocks.  A deterministic sphere search followed by local
refinement selects a coefficient vector \(c\) with positive violation

\[
 \sum_y\left\lVert\sum_z c_z\sigma_{zy}\right\rVert_1
 -\left\lVert\sum_z c_z\rho_z\right\rVert_1>0.
\]

The numerical search chooses a useful direction but is not trusted as a
proof.  Once chosen, the stored vector defines an ordinary checkable instance
of (1), and all 98 spectral branches are solved.  Any child still at or above
\(0.758\) is expanded again.  The final tree contains every branch bound, not
only a count of the closed leaves.

The certificate is conditional on CLARABEL's floating-point optimality and
infeasibility decisions at the recorded tolerances.  It is not an interval
proof.  The geometric exhaustiveness of each 98-way split and the tree
topology are independent of the separator optimizer and can be checked
without solving an optimization problem.

## Reproduction and audit

Regenerate the tree from the compact seed with:

```bash
python scratch/d2_frontier/adaptive_multicolumn_branch_tree.py \
  scratch/d2_frontier/adaptive_multicolumn_seed_l055.json \
  --max-expansions 50 --separator-samples 20000 --separator-starts 12 \
  --output scratch/d2_frontier/adaptive_multicolumn_worstcell_l055_auditable.json
```

The fast audit checks the exhaustive branch labels, closure decisions,
parent-child links, leaf counts, and stored maxima:

```bash
python scratch/d2_frontier/audit_adaptive_multicolumn_certificate.py \
  scratch/d2_frontier/adaptive_multicolumn_worstcell_l055_auditable.json
```

Adding `--recompute` rebuilds and re-solves every one of the 4,704 conic
subproblems before comparing their bounds with the stored certificate.

## What remains

The earlier Fourier/pair cover has 306 Fourier cells whose maximum pair
branch is at or above \(0.758\).  The present tree closes the single largest
base branch within that collection.  A full regular-interior certificate must
enumerate every still-open pair branch in those cells and run the same
adaptive process until the global queue is empty.  If arbitrary
multicolumn trace-norm contraction stalls at a positive-map relaxation, the
remaining regular case will require the determinant-scaled Choi positivity
condition rather than more scalar coefficient cuts.

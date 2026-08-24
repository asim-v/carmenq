# Adaptive multicolumn contraction certificate

## Result and scope

The complete Fourier/pair outer cover of the fixed \(\lambda=0.55\) benchmark
is now excluded at score \(0.758\).  Of the 5,376 base branches, 4,670 were
already below target.  The remaining 706 Bloch branches form 353 exact
complex-conjugation orbits, and every representative now has its own adaptive
tree with 98 exhaustive spectral branches at each expansion.

The audited global forest contains 2,698 expansions and 262,059 closed leaves,
including 151,733 conically infeasible leaves and 236 nodes closed by a fresh
source solve.  Its maximum depth is 7 and its largest finite terminal upper
bound is \(0.7579983961\).  There are no missing or open orbits.  The independent
full-tree identity

\[
 262059=353+(98-1)2698
\]

also holds, so the aggregate cannot have silently dropped a branch.

The original worst cell remains a useful compact example.  Its calculation
uses 48 adaptive contractions and has 4,657 terminal leaves: 2,640 are
conically infeasible and the largest finite terminal upper bound is
\(0.7579901084\).

That base cell is planar Fourier sector 4 of 8, spherical cube-face cell 18 at
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

### Positive-map completeness on an operator basis

Suppose the four inputs \(\{A_z\}_{z=0}^3\) form a real basis of the Hermitian
qubit operators.  The conditioned outputs then define a unique
Hermiticity-preserving linear interpolation

\[
 \Gamma(A_z)=\bigoplus_y\sigma_{zy}.
\]

If the behaviour constraints preserve trace on the four basis elements, the
extension is trace preserving everywhere.  In this regular case, the family
of inequalities (1) for every real \(c\) is equivalent to positivity of
\(\Gamma\).  This is the finite-basis specialization of the established
characterization that a trace-preserving linear map is positive exactly when
its induced Schatten-1 norm is one; see Pérez-García, Wolf, Petz, and Ruskai,
[*J. Math. Phys.* **47**, 083506
(2006)](https://doi.org/10.1063/1.2218675), especially their discussion after
Theorem II.1.  The equivalence itself is therefore prior art.

One direction is the usual trace-norm contraction of a positive
trace-preserving map on Hermitian inputs.  Conversely, take any \(X\succeq0\)
and expand it in the input basis.  Contractivity and trace preservation imply

\[
 \lVert\Gamma(X)\rVert_1
 \leq \lVert X\rVert_1
 =\operatorname{Tr}X
 =\operatorname{Tr}\Gamma(X).
\]

For every Hermitian operator \(Y\),
\(\lVert Y\rVert_1\geq|\operatorname{Tr}Y|\), with equality at positive
trace only when \(Y\succeq0\).  Hence \(\Gamma(X)\succeq0\), proving the
converse.

Accordingly, an exhaustive all-coefficient tree imposes one common *positive*
flagged instrument exactly at the level of the interpolating input-output
data.  It does not impose complete positivity.  The latter remains the
degree-four Choi matrix condition described in
`operator_basis_instrument_criterion.md`.  This theorem both explains why the
multicolumn hierarchy is substantially stronger than a few pairwise cuts and
gives a sharp kill criterion: if the tree converges to a positive but non-CP
point above target, scalar trace-norm directions cannot remove it.

The gap is concrete.  Put the four input states at tetrahedral Bloch vectors
and send them through matrix transposition into one flag.  Transposition is
positive and trace preserving, so every real-coefficient contraction is
saturated.  Its Choi matrix is the swap operator and has a negative
eigenvalue, so the map is not completely positive.  The regression test
`test_all_contractions_deliberately_retain_the_positive_non_cp_transpose`
keeps this limitation visible in the implementation.

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
infeasibility decisions at the recorded tolerances, with SCS fallback after a
CLARABEL solver exception.  A solver exception by itself is stored as an open
infinite bound and never closes a leaf.  This is not an interval proof.  The
geometric exhaustiveness of each 98-way split and the tree topology are
independent of the separator optimizer and can be checked without solving an
optimization problem.

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

The global aggregate is reproduced and structurally audited with:

```bash
python scratch/d2_frontier/aggregate_multicolumn_regular_forest.py \
  scratch/d2_frontier/fourier_behavior_pair_cover_p8_g4_l055_auditable.json \
  scratch/d2_frontier/regular_multicolumn_forest_top3_l055.json \
  scratch/d2_frontier/regular_multicolumn_forest_batch0_l055.json \
  scratch/d2_frontier/regular_multicolumn_forest_batch1_l055.json \
  scratch/d2_frontier/regular_multicolumn_forest_batch2_l055.json \
  --output scratch/d2_frontier/regular_multicolumn_forest_complete_l055.json
```

This command audits all 353 component trees, rejects duplicate or missing
orbits, checks the 98 branch labels at every expansion, verifies all
parent--child links and closure decisions, and enforces the global leaf
identity.  Replaying every conic solve is deliberately a separate, much more
expensive operation.

## What remains

Together with the eleven non-fully-vectorial spectral regimes, whose largest
bound is \(0.75115255\), this forest closes the fixed terminal-POVM/prior-box
benchmark at target \(0.758\), conditional on the recorded conic solves.  It
does **not** close the unrestricted interior frontier: the terminal POVM has
been fixed to weights \((0.92,0.64,0.44,0)\), only \(\lambda=0.55\) is covered,
and no interval-verified numerical enclosure is claimed.  The next global
task is a certified outer cover of terminal-measurement geometry (and then
parameter continuation), not additional splitting of these 353 cells.

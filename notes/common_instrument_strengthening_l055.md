# Common-instrument strengthening at the interior direction lambda = 0.55

**Date:** 23 August 2026

**Status:** exact fixed-input oracle, sparse order-two hierarchy, and a literal
shared-instrument spatial model implemented; the fixed interior upper bound
remains open

## 1. What is being imposed

Let four subnormalised qubit states \(\rho_z\) enter a four-outcome quantum
instrument.  The sixteen conditioned outputs are not independent states.  One
collection of completely positive maps must generate all of them:

\[
 \sigma_{zy}=\Phi_y(\rho_z),\qquad
 \sum_y\Phi_y\text{ trace preserving}.
\]

Equivalently, the flagged channel

\[
 \Lambda(\rho)=\bigoplus_y\Phi_y(\rho)
\]

is one CPTP map shared by every input label.  In input-major Choi convention,
the exact fixed-input feasibility system is

\[
 J_y\succeq0,\qquad
 \sum_y\operatorname{Tr}_{\rm out}J_y=I,
\]

\[
 \sigma_{zy}=\operatorname{Tr}_{\rm in}
 [ (\rho_z^{\mathsf T}\otimes I)J_y].
\]

The new public module `carmenq.common_instrument` implements this system as an
SDP projection.  A nonzero residual produces a normalized separating witness
\(W_{zy}\), and a second support SDP checks the separation independently.

## 2. Inexpensive necessary cuts

Trace-norm contractivity of the flagged channel gives, for every input pair
and every \(t\geq0\),

\[
 \sum_y\lVert\sigma_{zy}-t\sigma_{z'y}\rVert_1
 \leq
 \lVert\rho_z-t\rho_{z'}\rVert_1.
\tag{1}
\]

The facially reduced Choi moment programme now accepts an arbitrary finite
scale grid.  For a physical rank-one moment point, its quadratic constraints
are consequences of equation (1).  They therefore retain every physical
common-instrument realization.

The benchmark fixes the ternary terminal POVM with effect traces

\[
 (0.92,0.64,0.44,0),
\]

representative of the difficult terminal region in the completed
\(\lambda=0.55\) coarse cover, and fixes prefix-prior order \((0,1,2,3)\).
The single-scale and five-scale first-level bounds are

\[
 U_{t=1}=0.7633285377447268,
\]

\[
 U_{\{1/4,1/2,1,2,4\}}=0.7633274596929631.
\]

The improvement, about \(1.08\times10^{-6}\), is negligible.  This negative
result is informative: merely sampling more Alberti--Uhlmann-type comparisons
does not repair a high-rank moment relaxation.

## 3. Exact audit of the relaxed maximizer

The reported first-moment prefix states and conditioned outputs were extracted
and audited independently.  Every one of the 30 pair/scale comparisons in the
five-scale grid is violated.  The worst one is \((z,z',t)=(0,1,1)\):

\[
 \lVert\rho_0-\rho_1\rVert_1
 =1.9747\times10^{-5},
\]

while

\[
 \sum_y\lVert\sigma_{0y}-\sigma_{1y}\rVert_1
 =0.5729649212.
\]

The exact Choi projection gives

\[
 d_{\rm Choi}=0.2009299661,
 \qquad
 g_W=0.2009297479.
\]

Thus the first moments of the relaxed maximizer are not close to any common
instrument.  The moment constraints survive only because the 76-dimensional
moment matrix has full numerical rank and stores artificial variance in its
second moments.

## 4. A witness that remains valid when inputs vary

For the separating witness define

\[
 M_z=\max_y\lVert W_{zy}\rVert_\infty.
\]

If \(h_W(\bar\rho)\) is the exact support value at reference inputs
\(\bar\rho_z\), flagged-channel contractivity gives the robust bound

\[
 \sum_{z,y}\operatorname{Tr}[W_{zy}\Phi_y(\rho_z)]
 \leq
 h_W(\bar\rho)
 +\sum_z M_z\lVert\rho_z-\bar\rho_z\rVert_1.
\tag{2}
\]

For the audited maximizer,

\[
 M=(0.34955,0.35035,0.25335,0.23197),
\]

and equal trace-radius balls preserve positive separation up to radius
\(0.16953\).  Imposing equation (2) on the conservative radius-\(0.1\) branch
reduces the local moment bound to

\[
 U_{\rm local}=0.7342319934643531.
\]

The drop of \(0.02909547\) shows that the compatibility oracle removes the
spurious basin rather than merely perturbing it.  This is a local upper bound:
it does not constrain prefix states outside the four trace-distance balls and
must not be compared as a global converse with the known physical lower bound.

## 5. Consequence for the complete frontier

The next certification target is no longer an undirected higher moment level.
It is a spatial common-instrument branch-and-cut algorithm:

1. cover the prefix-state domain by trace-distance cells;
2. solve the current moment upper model on each cell;
3. audit the cell maximizer by exact Choi projection;
4. add equation (2) whenever its robust radius covers the cell;
5. split only cells whose strengthened upper bound remains above the target.

Every terminal POVM geometry and prefix-prior order still require an outer
cover.  Until that tree is complete, the result is a strong obstruction and a
working separation mechanism, not the full \(\lambda=0.55\) theorem.

## 6. Reproduction

The principal artifacts are:

- `src/carmenq/common_instrument.py`: public trace-norm and exact Choi tools;
- `scratch/d2_frontier/audit_common_instrument_candidate.py`: independent
  candidate audit and witness archive;
- `scratch/d2_frontier/choi_moment_reduced_upper.py`: multi-scale and robust
  branch-local witness constraints;
- `scratch/d2_frontier/common_instrument_exact_scip.py`: exact nonconvex model
  with one shared collection of Choi matrices;
- `scratch/d2_frontier/validate_exact_common_instrument.py`: independent
  matrix audit and inward repair of the physical candidate;
- `scratch/d2_frontier/validate_common_instrument_strengthening.py`: numerical
  chain validator.

Run the validator from the repository root:

```text
python scratch/d2_frontier/validate_common_instrument_strengthening.py
```

## 7. The obstruction is genuinely multicolumn

The first-moment maximizer was also tested without referring to Choi
coordinates.  Its complete \(4\times12\) terminal behaviour cannot be produced
by a common four-preparation qubit model.  Exhaustive checks of all twelve
single columns and all sixty-six pairs found no obstruction, whereas a
three-column dual does.  Thus the minimum validated support has cardinality
three.  This rules out the tempting explanation that one anomalous probability
or one incompatible pair is responsible for the gap.

The selected-column clique hierarchy consequently couples three different
continuation branches.  At SCS tolerance \(10^{-4}\), a \(273\times273\)
state--effect cross bridge gave

\[
 U_{\rm selected}=0.7634298037.
\]

The number is above the first-level Choi bound because the selected behaviour
model leaves the other nine columns free; it is not a nested comparison with
the full common-instrument model.  Its role is diagnostic: the obstruction
requires joint dimension geometry across at least three outcomes.

## 8. Sparse order-two common-instrument hierarchy

To impose one instrument more directly, the new programme
`common_instrument_sparse_order2.py` uses eighty physical coordinates: sixteen
for the four input states and sixty-four for four Choi matrices. Every pair
\((z,y)\) has the cross basis

\[
 \{1,\rho_z,J_y,\rho_zJ_y\}.
\]

For every outcome \(y\), a bridge contains all four state blocks and the same
\(J_y\). Choi positivity is localised against all state coordinates, while
state positivity is localised against the relevant Choi coordinates.  The
trace-preservation equations are multiplied by every state monomial through
degree two.

A further \(169\times169\) instrument bridge contains all degree-two state
monomials and all partial-trace Choi coordinates across the four outcomes.
Inside this positive matrix,

\[
 L_\alpha=\sum_y j_{y,\alpha0}-2\delta_{\alpha0}
\]

is forced into the moment kernel.  This strengthens trace preservation from a
first-moment equality to a shared multibranch identity.  The fixed rank-one
terminal POVM is imposed using the exact exposed Helstrom face

\[
 Y-\tau_s=\kappa_s(I-\Pi_s),\qquad \kappa_s\geq0,
\]

for the three active outcomes.  This removes the singular complementary-
slackness cone that contaminated earlier numerical solves.

The resulting model contains 64,509 scalar moments. At SCS tolerance
\(10^{-5}\), it reports

\[
 U_{\rm sparse\ O2}=0.7633280830,
\]

with minimum moment eigenvalue \(-1.8\times10^{-10}\), minimum localiser
eigenvalue \(-6.4\times10^{-8}\), and independent Helstrom residual
\(7.0\times10^{-15}\). The value differs from the first-level PPT bound by
only \(7.0\times10^{-7}\), within the numerical scale of the two solves.

This is a negative but decisive result: merely adding the natural sparse
order-two state--instrument moments does not close the interior point.

## 9. Why the hierarchy stalls

At the order-two optimum the first moments of the four input Bloch vectors are
nearly zero, and \(\rho_0\) and \(\rho_1\) differ in trace norm by only
\(5.3\times10^{-6}\). Nevertheless, their conditioned first-moment outputs
differ by \(0.54434\) in flagged trace norm. The exact fixed-input Choi
projection remains far away:

\[
 d_{\rm Choi}=0.19979246,\qquad g_W=0.19979224.
\]

The contradiction is carried by covariance moments such as
\(\mathbb E[\rho J]-\mathbb E[\rho]\mathbb E[J]\). Positive moment matrices
describe pseudo-distributions over state--instrument parameters; sharing the
same symbolic \(J_y\) does not force that pseudo-distribution to be a single
rank-one evaluation.  No convex constraint of this finite order can simply
declare all coordinate variances zero without losing the upper-bound
property.

Consequently, the remaining gap is now localized.  It is not missing Choi
positivity, trace preservation, Helstrom optimality, pairwise data processing,
or the obvious order-two cross moments.  It is the flatness/determinism gap of
the polynomial relaxation.

## 10. Deterministic spatial remedy

For a state-coordinate cell \(C\), interval arithmetic gives the valid bound

\[
 R_{zz'}(t;C)=
 \max\!\left\{
   \max_C|a_z-ta_{z'}|,
   \sqrt{\sum_{k=1}^3\max_C|r_{zk}-tr_{z'k}|^2}
 \right\}.
\]

The branch model now imposes

\[
 \sum_y\|\sigma_{zy}-t\sigma_{z'y}\|_1\le R_{zz'}(t;C)
\]

directly.  This upper bound cannot be inflated by second-moment variance and
converges to the exact input trace distance as the cell diameter tends to
zero.  It complements the robust Choi separating cuts: the witness controls a
neighbourhood of a full four-state tuple, whereas the contraction cut controls
large regions described only by pairwise differences.

The current theorem-level statement remains deliberately limited.  The fixed
interior benchmark has a working convergent branch mechanism and a diagnosed
relaxation gap, but its complete spatial cover has not yet been closed.  No
claim about the entire \(\lambda=0.55\) frontier follows until that fixed-POVM
tree and the outer terminal-geometry cover are both complete.

A twelve-node pilot quantifies the practical limitation of the present
axis-aligned cover.  It generated twelve independent Choi witnesses, closed no
cell, and left twenty-five cells open.  The maximum inherited upper bound moved
only from \(0.7633296353\) to \(0.7633295671\). The method is convergent in the
zero-diameter limit, but this parametrisation is not competitive enough to be
the primary completion route without a symmetry quotient, diagonal
difference cells, or a mixed-integer representation of the behaviour
disjunctions.

The expanded numerical chain is checked by

```text
python scratch/d2_frontier/validate_common_instrument_hierarchy.py
```

## 11. Literal shared-instrument formulation

The hierarchy diagnosis motivates a nonconvex model with no pseudo-moment
interpretation. The programme `common_instrument_exact_scip.py` contains four
physical state variables and exactly four shared Choi variables. Positivity is
represented by complex Cholesky factors,

\[
 J_y=L_yL_y^\dagger,
\]

and every output coordinate obeys the same bilinear evaluation rule,

\[
 s_{zy,\nu}=\frac12\sum_{\mu=0}^3
 \eta_\mu r_{z,\mu}j_{y,\mu\nu},
 \qquad \eta=(1,1,-1,1).
\]

The four global equations

\[
 \sum_y j_{y,\mu0}=2\delta_{\mu0}
\]

make the collection trace preserving. There is no outcome-, column-, or
input-dependent copy of \(J_y\). Thus every exactly feasible point is a
literal instrument satisfying
\(\sigma_{zy}=\Phi_y(\rho_z)\) simultaneously for all sixteen pairs.

A simultaneous input-unitary gauge fixes three continuous directions while
leaving every conditioned output invariant. The invariance was checked to
\(4.2\times10^{-17}\) on the archived seed. In 120 seconds, SCIP explored 56
nodes and reported

\[
 L_{\rm SCIP\ incumbent}=0.7239479693,
 \qquad
 U_{\rm exact\ QCQP}=0.7699031152.
\]

The raw state matrices had a minimum eigenvalue of
\(-9.1\times10^{-10}\), consistent with the solver feasibility tolerance. An
independent inward repair mixes only \(7.7\times10^{-9}\) of the maximally
mixed state into each prefix state, preserves all priors, and gives the
explicit physical value

\[
 L_{\rm repaired}=0.7239479195.
\]

For that repaired strategy the minimum state, Choi, and output eigenvalues are
respectively \(1.0\times10^{-12}\), \(3.5\times10^{-11}\), and
\(1.9\times10^{-11}\); the trace-preservation residual is
\(3.3\times10^{-16}\). Its fixed terminal POVM is only
\(1.93\times10^{-9}\) below the independently recomputed Helstrom optimum.
The physical value improves the preceding archived checkpoint by
\(0.00135907\).

This cleanly separates two statements. The repaired lower bound is an
explicit shared-instrument construction. The SCIP dual is only a
solver-conditional global upper bound and remains above the target
\(0.758\); it is not a proof that the target is attainable.

## 12. Mixed-integer behaviour disjunctions

The minimum three-column ellipsoid obstruction supplies a valid union: for a
physical qubit behaviour, at least one supported half-space must be satisfied.
Each archived witness therefore becomes one small mixed-integer disjunction,
rather than an invalid simultaneous collection of all its half-spaces. With
four active witnesses, a 300-second spatial solve gave

\[
 U_{\rm disj}=0.7699022234.
\]

The relaxation incumbent remained at \(0.7633299692\) and generated a fifth
incompatibility witness, so it is not a physical lower bound. The independent
exact-Choi and behaviour-disjunction upper bounds differ by only
\(8.9\times10^{-7}\). This agreement is useful validation, but both are still
time-limited numerical bounds and neither closes \(0.758\).

## 13. Literature positioning and novelty gate

The fixed-input feasibility SDP is an instance of the established completely
positive interpolation problem.  Ambrozie and Gheondea give necessary and
sufficient affine-positive criteria for mapping a finite matrix family to
prescribed outputs and explicitly allow trace-preserving restrictions
([arXiv:1308.0667](https://arxiv.org/abs/1308.0667)).  The projection oracle and
its witness are useful implementations for this benchmark, but their existence
is not a new theorem.

Likewise, convergent SDP hierarchies for constrained bilinear optimization and
finite de Finetti approximations to convex hulls of product quantum channels
already exist; see Berta, Borderi, Fawzi, and Scholz
([arXiv:1810.12197](https://arxiv.org/abs/1810.12197)).  Dimension-constrained
moment hierarchies and sequential prepare-and-measure hierarchies are also
established nearby frameworks
([arXiv:1308.3410](https://arxiv.org/abs/1308.3410),
[arXiv:2409.17185](https://arxiv.org/abs/2409.17185)).  Channel-compatibility
witnesses have a separate mature formulation in the quantum channel marginal
problem ([arXiv:2102.10926](https://arxiv.org/abs/2102.10926)).

Therefore neither “use a Choi matrix,” “add a moment hierarchy,” “factor a PSD
matrix,” nor “extract a separating witness” passes the novelty gate. A
defensible contribution would have to be the problem-specific result that
those tools enable: for example, a
closed interior bound, a provably convergent and quantitatively superior
state-cell certificate, a new finite-level flatness criterion tailored to a
shared instrument, or a nontrivial theorem connecting minimum multicolumn
behaviour obstructions to branch complexity.  The current work supplies the
diagnosis and machinery for such a result, but does not yet claim one.

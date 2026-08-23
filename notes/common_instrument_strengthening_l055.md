# Common-instrument strengthening at the interior direction lambda = 0.55

**Date:** 23 August 2026

**Status:** exact fixed-input compatibility oracle and validated branch-local
strengthening; a global covering tree remains open

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
- `scratch/d2_frontier/validate_common_instrument_strengthening.py`: numerical
  chain validator.

Run the validator from the repository root:

```text
python scratch/d2_frontier/validate_common_instrument_strengthening.py
```

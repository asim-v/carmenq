# Terminal-reconstructed common-instrument frontier at lambda = 0.55

## Result in one sentence

The common-instrument contractions can now be imposed throughout a variable
ternary terminal cell while retaining the two Bloch coordinates reconstructed
by that POVM.  On the previously worst cell this improves the open `pbb`
Fourier family from `0.7631780712` to `0.7613156805`; the `bbb` family remains
at `0.7666498770`, and one adaptive coefficient lowers its selected worst
angular cell to `0.7665057550`, still above the target `0.758`.

This is a solver-conditional strengthening and a concrete bridge to the global
common-instrument problem.  It is not a local or global closure.

## Why the earlier measured contraction was insufficient

Let `rho_z` be the four subnormalised qubit inputs and let a common instrument
produce subnormalised outputs

\[
  \sigma_{zy}=\Phi_y(\rho_z).
\]

For any real coefficient vector `c`, complete positivity and trace
preservation of the flagged map imply

\[
  \sum_y\left\|\sum_z c_z\sigma_{zy}\right\|_1
  \leq
  \left\|\sum_z c_z\rho_z\right\|_1.
\]

Measuring every output with the terminal POVM gives signed statistics

\[
  q_{yt}(c)=\sum_z c_z\,\mathrm{tr}(P_t\sigma_{zy})
\]

and hence the weaker classical contraction

\[
  \sum_{y,t}|q_{yt}(c)|
  \leq
  \left\|\sum_z c_z\rho_z\right\|_1.
\]

That inequality is valid throughout the terminal geometry, but the local
pilots show that it throws away too much operator information.  The worst
`pbb` cell remains at `0.7631780712`; adaptive coefficients expose large local
violations, but the relaxation can rotate the unknown input Bloch vectors and
move the violation to another coefficient.

## The visible-Bloch reconstruction

For a nondegenerate planar rank-one ternary POVM, its three probabilities
determine the trace and two Bloch coordinates of any Hermitian qubit operator.
In the canonical chart, with reciprocal Horwitz parameters `A,B`, write

\[
  D=A+B-1,\qquad
  w_0=\frac{A}{D},\quad
  w_1=\frac{B}{D},\quad
  w_2=\frac{A+B-2}{D},
\]

and

\[
  \gamma
  =1-\frac{2}{A}-\frac{2}{B}+\frac{2}{AB},
  \qquad s=\sqrt{1-\gamma^2}.
\]

If `q=(q_0,q_1,q_2)` is the terminal probability vector of an operator with
trace `d`, then `d=q_0+q_1+q_2` and

\[
  \begin{pmatrix}x\\y\end{pmatrix}=R(A,B)q,
\]

where

\[
  R_x=\left(\frac{2}{w_0}-1,-1,-1\right)
\]

and

\[
  R_y=\frac{1}{s}\left(
    -1-\gamma R_{x0},
    \frac{2}{w_1}-1+\gamma,
    -1+\gamma
  \right).
\]

For a Hermitian qubit operator with Bloch coefficients `(x,y,z)`,

\[
  \|X\|_1=\max\{|d|,\sqrt{x^2+y^2+z^2}\}
  \geq \max\{|d|,\sqrt{x^2+y^2}\}.
\]

Therefore the reconstructed planar norm is a lower bound on every output-block
trace norm.  It can replace the classical `L1` term in the left side of the
common-instrument contraction without assuming that the physical output is
itself planar.

## A box-safe version

For a terminal cell `Theta`, choose a midpoint reconstruction `R_0` and
outward-rounded column errors

\[
  e_t\geq\sup_{\theta\in\Theta}
  \|(R(\theta)-R_0)_{:,t}\|_2.
\]

Then

\[
  \|R(\theta)q_y(c)\|_2
  \geq
  \|R_0q_y(c)\|_2
  -\sum_t e_t|q_{yt}(c)|.
\]

Because every raw path statistic is nonnegative,

\[
  |q_{yt}(c)|
  \leq \sum_z |c_z|q_{zyt},
\]

so the error term has an affine upper bound.  One second-order cone per output
flag imposes

\[
  \ell_y+\sum_{z,t}|c_z|e_tq_{zyt}
  \geq \|R_0q_y(c)\|_2,
  \qquad
  \ell_y\geq|d_y(c)|.
\]

The common-instrument condition is then

\[
  \sum_y\ell_y
  \leq
  \left\|\sum_z c_z\rho_z\right\|_1.
\]

The implementation encloses every scalar coefficient of `R(A,B)` with
elementary interval arithmetic and applies `nextafter` outward after every
operation.  This is stronger than dense sampling and remains auditable without
an interval package.  Cells whose sine interval reaches zero are explicitly
rejected; the projective boundary must be covered separately.

## Spectral cover on the input side

For a qubit Hermitian combination with scalar coefficient `a(c)` and Bloch
vector `v(c)`,

\[
  \left\|\sum_zc_z\rho_z\right\|_1
  =\max\{|a(c)|,\|v(c)\|_2\}.
\]

The reverse-convex maximum is covered by three branches: scalar positive,
scalar negative, and Bloch active.  Global rotational symmetry fixes the first
active Bloch vector to `+z` and the second to the `xz` half-plane; proved plane
and cube-face caps cover the remaining directions.  The present pilot uses the
three nontrivial `Z_2^2` Fourier coefficient vectors.

## Numerical evidence archived in the repository

The parent terminal cell is

\[
  A\in[1.923828125,1.92578125],\qquad
  B\in[1.1453718354430378,1.149525316455696].
\]

Its earlier probability-cone bound was `0.7666498763533095`.  The following
artifacts isolate the strengthening:

- `ternary_common_instrument_top_leaf_pbb_p16_l055.json`: measured-only `pbb`,
  maximum `0.7631780711822648`;
- `ternary_common_instrument_top_leaf_coarse_l055.json`: measured-only full
  coarse spectral cover, whose `bbb` maximum is `0.7666498773672448`;
- `ternary_reconstructed_fourier_top_leaf_pbb_p16_l055.json`: reconstructed
  `pbb`, maximum `0.7613156804571753`;
- `ternary_reconstructed_fourier_top_leaf_bbb_p4_g4_l055.json`: reconstructed
  `bbb` with a 384-cell angular cover, maximum `0.766649877013737`;
- `ternary_reconstructed_multicolumn_top_leaf_bbb_p1_s92_l055.json`: one
  exhaustive adaptive contraction on the selected worst `bbb` cell, maximum
  `0.7665057550449533`.

Measured-only and reconstructed adaptive pilots repeatedly moved the violation
to a new real coefficient after the previous one was imposed.  In the
reconstructed tree, the second expansion closes all 98 children of the chosen
orientation, but 95 sibling orientations inherited from the first expansion
remain open.  This negative result isolates the remaining obstacle as branch
compression rather than absence of a separating inequality.

## Relation to prior work and novelty boundary

Trace-norm or base-norm contractivity under positive trace-preserving maps is
standard.  Testing-region criteria for state convertibility include the
Alberti--Uhlmann line and its multi-state extensions.  Positive-semidefinite
rank and nested-spectrahedron geometry likewise provide the established
language for common qubit behaviours.  Relevant starting points include:

- Reeb, Kastoryano, and Wolf, *Hilbert's projective metric in quantum
  information theory*, arXiv:1102.5170;
- Dall'Arno, Buscemi, and Scarani, *Extension of the Alberti-Uhlmann criterion
  beyond qubit dichotomies*, arXiv:1910.04294;
- Fawzi, Gouveia, Parrilo, Robinson, and Thomas, *Positive semidefinite rank*,
  arXiv:1407.4095;
- Kubjas, Robeva, and Robinson, *Positive semidefinite rank and nested
  spectrahedra*, arXiv:1512.08766.

Accordingly, neither trace-norm contraction nor planar tomography is claimed
as new.  The candidate contribution is the combination of:

1. an adaptive all-coefficient common-instrument certificate;
2. terminal-geometry interval reconstruction inside that certificate;
3. exact spectral/angular covers of the reverse-convex input norm; and
4. integration with the audit/return upper-bound problem.

Publication-level novelty still requires a completed global certificate and a
literature comparison against quantum statistical comparison, qubit behaviour
factorisation, and robust tomography bounds.

## Remaining step

Even the worst local terminal cell remains open.  The immediate task is to
compress the product of adaptive spectral branches.  A one-hot MISOCP encoding
has been implemented, but the local SCIP binding currently collides with a
second OpenMP runtime; the unsafe duplicate-runtime override is deliberately
not used.  The safe fallback is the explicit SOCP tree or a separate clean SCIP
environment.

After local closure, the previous terminal cover still has 5,550 finite leaves
above `0.758`.  Of those, 5,363 admit a nondegenerate reconstruction enclosure
as currently partitioned; 187 touch a projective boundary and must be closed by
further terminal subdivision plus the independent projective envelope.  A
global certificate must reuse parameterised reconstruction SOCPs across those
leaves, combine all spectral/angular branches, and count every inherited and
newly split leaf.

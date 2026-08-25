# Common-effective-POVM frontier at lambda = 0.55

## Result

The adaptive trace-norm frontier has reached its kill criterion.  On the
selected worst continuous terminal/Fourier cell, adding a fourth globally
valid common-instrument contraction required 21,190 SOCP solves, increased the
number of open spectral cells from 815 to 2,216, and did not improve the
maximum bound beyond numerical noise:

\[
  0.7635145903\longrightarrow 0.7635145905.
\]

The obstruction is not close to a physical common instrument.  Reconstructing
the unique effective twelve-outcome POVM from the depth-four maximiser gives
ten nonpositive effects.  The worst positivity margin is

\[
  -0.1532642720,
\]

while the input operator basis has determinant `0.0089984985` and condition
number `3.2406541`.  The violation is therefore large and is not caused by an
ill-conditioned basis.

This suggests a sharper route: impose the finite common-effective-POVM
criterion jointly, rather than discover its separating hyperplanes one at a
time.

## Exact common-effective-POVM criterion

Let four subnormalised qubit states be written in Pauli coordinates as

\[
  \rho_z=\frac12\left(p_z I+\mathbf r_z\cdot\boldsymbol\sigma\right),
  \qquad z=0,1,2,3,
\]

and define the operator-basis matrix

\[
  R=
  \begin{pmatrix}
    p_0 & \mathbf r_0^{\mathsf T}\\
    p_1 & \mathbf r_1^{\mathsf T}\\
    p_2 & \mathbf r_2^{\mathsf T}\\
    p_3 & \mathbf r_3^{\mathsf T}
  \end{pmatrix}.
\]

Flatten the measured path statistics into a `4 x 12` matrix `Q`, with

\[
  Q_{z,(y,t)}=q_{zyt}.
\]

Suppose first that `det R` is nonzero.  There is then exactly one Hermitian
family of input effects `F_yt` reproducing these probabilities.  If

\[
  F_k=a_{0k}I+\mathbf a_k\cdot\boldsymbol\sigma,
\]

its coordinates are the columns of

\[
  A=R^{-1}Q.
\]

The statistics admit a common qubit POVM if and only if

\[
  a_{0k}\geq\|\mathbf a_k\|_2
  \qquad\text{for every }k=(y,t),
\]

and the effects sum to the identity.  In the present model the second
condition is automatic.  Row normalisation gives

\[
  Q\mathbf 1=(p_0,p_1,p_2,p_3)^{\mathsf T}=R e_0,
\]

so invertibility implies

\[
  A\mathbf 1=e_0.
\]

Thus twelve Lorentz-cone tests are necessary and sufficient for the common
effective POVM.

## Why every common instrument must pass this test

If `{Phi_y}` is a quantum instrument and `{E_t}` is the terminal POVM, define

\[
  F_{yt}=\Phi_y^*(E_t).
\]

Complete positivity makes every `F_yt` positive.  Trace preservation of the
flagged instrument gives

\[
  \sum_{y,t}F_{yt}
  =\sum_y\Phi_y^*\!\left(\sum_tE_t\right)
  =\sum_y\Phi_y^*(I)=I,
\]

and

\[
  q_{zyt}=\operatorname{Tr}(\rho_zF_{yt}).
\]

Consequently the common-effective-POVM criterion is an exact necessary
condition for the common instrument.  It is not sufficient for the more
specific sequential factorisation through the prescribed terminal POVM: a
POVM can be common without admitting that particular instrument/terminal
decomposition.  This asymmetry is useful for an upper certificate.  If the
weaker common-POVM set is already below the target, every physical common
instrument is excluded as well.

## Determinant-scaled polynomial SOC

The inverse can be eliminated.  Put

\[
  \delta=\det R,
  \qquad
  N=\operatorname{adj}(R)Q=\delta A.
\]

On a fixed determinant-sign branch `s = sign(delta)`, positivity is exactly

\[
  sN_{0k}\geq\|N_{1:3,k}\|_2,
  \qquad k=1,\ldots,12.
\]

Every entry of `N` has degree four in the joint input/statistics variables:
degree three from the adjugate and degree one from `Q`.  The condition is a
finite polynomial second-order-cone system.  It replaces the unbounded list
of real-coefficient contractions on the nonsingular stratum.  The singular
stratum `delta = 0` must remain a separate operator-system extension problem.

At fixed input states the criterion is an ordinary SOCP.  With both states and
statistics variable it is nonconvex and should be imposed by a determinant-
sign moment localiser, interval branch-and-bound, or spatial product
relaxation.  The important change is that positivity is now joint and finite;
angular spectral products are no longer required.

## Numerical audit of the depth-four maximiser

The archived maximiser has score bound `0.7635145905202629`.  The exact basis
audit reports:

| quantity | value |
|---|---:|
| `det R` | `0.00899849846918432` |
| `cond R` | `3.2406541491246794` |
| interpolation residual | `3.47e-17` |
| POVM completeness residual | `1.56e-16` |
| nonpositive effects | `10 / 12` |
| worst positivity margin | `-0.15326427204020443` |

The completeness and interpolation residuals show that reconstruction itself
is numerically exact at double precision.  The large negative margins show
that the high point violates common-measurement structure, before the stronger
question of complete positivity of each instrument subchannel is even asked.

## What the adaptive frontier established

The contraction tree was not wasted.  It established that the earlier
terminal reconstruction inequalities detect genuine common-map violations,
and it supplied a controlled scalability experiment:

| stage | open cells | maximum bound |
|---|---:|---:|
| first adaptive cover | `22` | `0.7665057558` |
| shared second separator | `246` | `0.7650422590` |
| third separator, coarse | `826` | `0.7647160479` |
| third separator, nested grid-4 refinement | `815` | `0.7635145903` |
| fourth separator, coarse | `2,216` | `0.7635145905` |

The last row is the decisive negative result: local separator generation has
entered a branch-proliferation regime.  A fifth blind separator is not the
next unit of progress.

## Novelty boundary

Operator-basis process and detector tomography are standard, as is quantum
statistical comparison.  The matrix inversion and Lorentz-cone positivity
test are therefore not claimed as new.  Relevant foundations include Ziman,
Plesch, and Buzek on superoperator reconstruction
([arXiv:quant-ph/0406088](https://arxiv.org/abs/quant-ph/0406088)), Buscemi on
comparison of quantum statistical models
([arXiv:1004.3794](https://arxiv.org/abs/1004.3794)), and Dall'Arno, Buscemi,
and Scarani on extensions of Alberti--Uhlmann comparison
([arXiv:1910.04294](https://arxiv.org/abs/1910.04294)).

The candidate contribution is the integration of this exact common-POVM
localiser with a continuous terminal-geometry upper certificate and, if
needed, its subsequent strengthening to the exact flagged-instrument Choi
criterion.  Publication value depends on closing a nontrivial region or
proving a useful impossibility/scaling theorem with that integrated method.

## Reproducibility

The arbitrary-depth angular engine is
`scratch/d2_frontier/ternary_extend_separator_frontier.py`.  The exact basis
audit is `scratch/d2_frontier/common_effective_povm_audit.py`.  Its numerical
output is `scratch/d2_frontier/common_effective_povm_audit_depth4_top_l055.json`.
The solver-independent reconstruction is also exposed as the stable library
function `carmenq.reconstruct_effective_povm_from_basis`, returning a
`BasisPovmReconstruction` with the unique effects, spectral margins,
determinant-scaled numerators, and residuals.
Unit tests include a physical POVM, a nonpositive unique reconstruction with
nonnegative sampled statistics, nested cube-face refinement, and legacy-to-
arbitrary-depth frontier conversion.

The next implementation target is the determinant-scaled twelve-SOC
localiser on one fixed-sign nonsingular terminal cell, followed by a separate
singular-stratum cover.  Only if that remains above `0.758` should the model
pay the higher cost of the full four Choi-matrix inequalities for the flagged
instrument.

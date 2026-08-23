# Complete sector exhaustion at the interior direction \(\lambda=0.6\)

**Date:** 23 August 2026
**Status:** complete solver-conditional numerical enclosure for the two-block
rank-two relaxation at one support direction; an exact closed form and
solver-independent interval certification remain open.

## 1. Result and precise scope

Let \(\beta_{4\mathrm s}(\lambda)\) be the support of the original four-slot
streamed problem and let \(\beta_{2\mathrm b}(\lambda)\) be the support of the
finite two-block rank-two relaxation defined in
`notes/two_block_choi_program.md`. At \(\lambda=0.6\), the explicit
four-effect streamed construction and a complete exhaustion of the relaxed
terminal-readout sectors give the stronger chain

\[
\boxed{
0.7658988152646944
\;\leq\;
\beta_{4\mathrm s}(0.6)
\;\leq\;
\beta_{2\mathrm b}(0.6)
\;\leq\;
0.76662 .}
\tag{1}
\]

The reported width is \(7.2118473530558\times10^{-4}\), or approximately
\(0.0942\%\) of the physical lower bound.  The unrounded assembly of the
sector bounds is \(0.76661\); equation (1) reports one additional decimal
safety step.

The lower endpoint is an exact physical construction evaluated in double
precision and independently checked as a complete four-slot instrument.  The
upper endpoint is a finite *solver-conditional numerical enclosure*: its
logical reductions are analytic, while its terminal boxes use CLARABEL values
with an explicit safety margin and its nonconvex spatial bounds use SCIP dual
bounds at recorded tolerances.  It is not an outward-rounded interval proof.
Nor does (1) prove that the explicit four-effect construction is the exact
optimum.

## 2. Fully active Helstrom reconstruction

Consider an optimal rank-one terminal POVM

\[
E_i=w_i\Pi_i,
\qquad
\sum_i E_i=I_2,
\qquad
\Pi_i=\frac{I+n_i\cdot\sigma}{2},
\tag{2}
\]

for subnormalised syndrome states \(\tau_i\), with
\(p_i=\operatorname{Tr}\tau_i\).  Write the Helstrom dual as

\[
Y=\frac A2\left(I+t\,\hat z\cdot\sigma\right),
\qquad A=\operatorname{Tr}Y,
\qquad 0\leq t\leq1,
\tag{3}
\]

and set \(x_i=\hat z\cdot n_i\).  If every displayed effect is active,
complementary slackness and two-dimensional positivity imply the exact
identity

\[
Y-\tau_i=(A-p_i)\Pi_i^\perp,
\qquad
\tau_i=Y-(A-p_i)\Pi_i^\perp.
\tag{4}
\]

Thus \(p_i\leq A\), and positivity of \(\tau_i\) is equivalent to

\[
p_i\geq A f_t(x_i),
\qquad
f_t(x)=\frac{1+2tx+t^2}{2(1+tx)}.
\tag{5}
\]

At the only singular corner, \((t,x)=(1,-1)\), equation (5) is understood
with the exact boundary value \(f_1(-1)=0\): the support expectation of
\(Y\) then vanishes and positivity requires no prior reserve.

POVM completeness supplies

\[
\sum_iw_i=2,
\qquad
\sum_iw_i n_i=0,
\qquad
\sum_iw_i x_i=0.
\tag{6}
\]

After choosing the longitudinal coordinates, transverse closure is no longer
an opaque vector constraint.  The lengths

\[
a_i=w_i\sqrt{1-x_i^2}
\tag{7}
\]

form a closed planar polygon exactly when

\[
2a_i\leq\sum_j a_j
\quad\text{for every }i.
\tag{8}
\]

For four effects this is the closed-quadrilateral criterion.  Equations
(5)--(8) give a low-dimensional necessary representation of the entire
four-active Helstrom geometry without reconstructing four Bloch vectors.

## 3. Exact projective comparisons

Fix two distinct active labels \(i\) and \(j\).  Retain \(\Pi_i\) as the
answer \(i\) and use \(I-\Pi_i\) as the answer \(j\).  Equation (4) gives the
binary-projective AUDIT score exactly:

\[
A_2^{i\leftarrow j}
=\operatorname{Tr}(\Pi_i\tau_i)
+p_j-\operatorname{Tr}(\Pi_i\tau_j)
=\frac12\left[(1-n_i\cdot n_j)A
+(1+n_i\cdot n_j)p_j\right].
\tag{9}
\]

This replacement does not alter RETURN.  Averaging equation (9) over
\(i\neq j\) with weights \(w_i/(2-w_j)\), and using vector closure in
equation (6), removes the remaining geometry:

\[
\boxed{
\overline A_{2,j}
=A-\frac{1-w_j}{2-w_j}(A-p_j).}
\tag{10}
\]

At least one of the projective comparisons is no smaller than its weighted
average.  Consequently, every independently certified projective support
line yields a simultaneous linear restriction on a fully active readout.
The enclosure uses the certified lines

\[
\beta_{\rm proj}(0.55)\leq0.75730,
\qquad
\beta_{\rm proj}(0.60)\leq0.76591.
\tag{11}
\]

The first is a 128-leaf rescaled secular cover.  Its worst rescaled dual is
\(0.9999892110586747\).  The second is the independent 312-leaf projective
cover already used for the physical support direction.

## 4. RETURN cap and finite readout arity

If \(r_{zy}\) are the refined classical path probabilities and \(p_s\) the
four syndrome priors, the pinching/Hellinger step gives

\[
R\leq\frac1{16}\left(\sum_{z,y}\sqrt{r_{zy}}\right)^2
\leq\frac14\left(\sum_s\sqrt{p_s}\right)^2.
\tag{12}
\]

The terminal discrimination objective is linear in the POVM.  An optimal
extreme qubit POVM may therefore be chosen.  The extremality rank bound
\(\sum_i\operatorname{rank}(E_i)^2\leq4\) leaves only the following
nontrivial possibilities: a binary projective readout, a ternary rank-one
readout, or a four-outcome rank-one readout.  Hence the sector split below is
exhaustive; it is not an ansatz over a selected family of measurements.
The one-outcome case is included as a degenerate binary-projective readout by
assigning both projector outcomes to the same answer.

There is also a useful deletion inequality.  If one rank-one effect has
\(w_k=\operatorname{Tr}E_k<\delta\), merge its outcome into any retained
answer.  RETURN is unchanged and the AUDIT loss is at most

\[
\operatorname{Tr}(E_k\tau_k)
\leq\lVert E_k\rVert_\infty p_k
\leq w_k<\delta.
\tag{13}
\]

Thus a four-active support obeys \(S_4\leq S_{\leq3}+\lambda\delta\).

## 5. Exhaustion of all terminal sectors

At \(\lambda=0.6\), use the thresholds

\[
w_{\max}=0.88325,
\qquad
\delta=0.0003.
\tag{14}
\]

The complete partition and recorded upper bounds are

| terminal sector | method | support upper bound |
|---|---|---:|
| at most two active effects | four-topology spatial projective cover | \(0.76591\) |
| three/four active, \(w_{\max}\leq0.88325\) | analytic prior-reserve cap | \(0.7658931806287275\) |
| three active, \(w_{\max}>0.88325\) | complete terminal-weight SOCP cover | \(0.76643\) |
| four active, \(w_{\max}>0.88325\), \(w_{\min}\geq0.0003\) | projected Helstrom spatial SCIP relaxation | \(0.7663946336432972\) |
| four active, \(w_{\min}<0.0003\) | deletion into the ternary sector | \(0.76643+0.6(0.0003)=0.76661\) |

The ternary cover contains 12,008 terminal boxes and 24,002 solved nodes.  Its
largest stored leaf bound is \(0.7664281458427126\), below its declared target
\(0.76643\), with no open boxes.  The fully active spatial run explored 15,361
nodes in 60 seconds and returned primal/dual values
\(0.7660017525973855/0.7663946336432972\).  The deletion row is therefore the
largest assembled upper bound and produces equation (1).

## 6. Reproduction and audit chain

The combined validator reconstructs the deterministic dyadic split tree of
the ternary cover, rejects missing, overlapping, or extra leaves, checks both
projective manifests, checks the active-sector artifacts, and recomputes the
final arithmetic:

The recorded environment used NumPy 2.2.6, SciPy 1.15.3, CVXPY 1.7.5,
CLARABEL 0.11.1, PySCIPOpt 6.2.1, and SCIP 10.0.2. It can be recreated with
`python -m pip install -e ".[dev,reproducibility]"`.

```bash
python scratch/d2_frontier/validate_projective_line_l055.py
python scratch/d2_frontier/validate_projective_cover.py
python scratch/d2_frontier/validate_interior_frontier.py
python -m pytest \
  tests/test_active_readout_geometry.py \
  tests/test_projective_secular_rescaling.py \
  tests/test_ternary_probability_cone.py \
  tests/test_inellipse_geometry.py -q
```

The decisive machine-readable manifest is
`scratch/d2_frontier/interior_frontier_l060_certificate.json`.  It references
the physical lower-bound artifact, the two projective covers, the complete
ternary cover, the prior-reserve cap, and the fully active SCIP run.

## 7. What is closed, and what is not

The logical gap left by the earlier “terminal projectivity” conjecture is
closed numerically at \(\lambda=0.6\): genuine ternary and four-active
readouts have now been bounded directly, so no projectivity assumption enters
equation (1).  This is the complete *interior direction* requested here.

Three stronger claims remain unwarranted.  The calculation does not yet give
the complete support curve for every \(\lambda\); it does not identify the
exact maximiser at \(\lambda=0.6\); and it is not a solver-independent formal
proof.  The next mathematically meaningful frontier is therefore either an
outward interval reconstruction of the recorded conic and spatial duals, or
an analytic inequality that collapses the upper endpoint of equation (1) to
the explicit four-effect value.

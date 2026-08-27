# Complete sector exhaustion at the interior direction \(\lambda=0.6\)

**Date:** 27 August 2026
**Status:** complete finite computer-assisted enclosure for the two-block
rank-two relaxation at one support direction. Verification of the conic
sectors calls no optimizer; an exact maximizer and the full support curve
remain open.

## 1. Result and precise scope

Let \(\beta_{4\mathrm s}(\lambda)\) be the support of the original four-slot
streamed problem and let \(\beta_{2\mathrm b}(\lambda)\) be the support of the
finite two-block rank-two relaxation defined in
`notes/two_block_choi_program.md`. At \(\lambda=0.6\), the explicit
four-effect streamed construction and a complete exhaustion of the relaxed
terminal-readout sectors give the stronger chain

\[
\boxed{
0.7658988152
\;\leq\;
\beta_{4\mathrm s}(0.6)
\;\leq\;
\beta_{2\mathrm b}(0.6)
\;\leq\;
0.76670 .}
\tag{1}
\]

The exact rational width is \(0.0008011848\), approximately \(0.10461\%\) of
the lower endpoint. The lower bound uses a fixed physical four-effect Choi-MPS
whose amplitudes are rational unit-circle parametrizations; four radicals in
its RETURN score are rounded downward at 192 dyadic bits. The computed
rational score is \(0.7658988152646940319\ldots\), so the displayed lower
endpoint is a strict outward truncation.

The upper endpoint combines directed interval arithmetic for the projective
sector with exact-rational residual dual replay for the conic sectors. Solvers
are used only to propose dual vectors during certificate generation. The
result is not a closed form or an end-to-end kernel-formalized proof: Python
matrix canonicalization, the interval implementation, and the documented
analytic reductions remain in the trust boundary.

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
\beta_{\rm proj}(0.60)\leq0.76600.
\tag{11}
\]

For each weight, the outward interval replay covers all 1,448 cells obtained
from the rank/rank, rank/endpoint, endpoint/rank, and endpoint/endpoint
topologies. Archived SCIP trees supply only the geometric partition; no
archived primal or dual value enters a cell inequality.

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
| at most two active effects | outward interval replay of four projective topologies | \(0.76600\) |
| three/four active, \(w_{\max}\leq0.88325\) | 576-cell exact dual replay | \(0.765893817258146\) |
| three active, \(w_{\max}>0.88325\) | 12,008-cell transferred exact dual replay | \(0.76652\) |
| four active, \(w_{\max}>0.88325\), \(w_{\min}\geq0.0003\) | 90-leaf McCormick--SOCP exact dual replay, six orders per leaf | \(0.76670\) |
| four active, \(w_{\min}<0.0003\) | deletion into the ternary sector | \(0.76652+0.6(0.0003)=0.76670\) |

The exact maximum reconstructed over all ternary cells is
\(0.7665158446239910257\ldots\). The compact four-active tree has 88 closed
and two domain-empty leaves; its independently replayed maximum is
\(0.7666751671519558614\ldots\). Both values are below their displayed
outward targets. The deletion row equals the largest reported sector bound
and therefore fixes equation (1).

## 6. Reproduction and audit chain

The recorded environment used NumPy 2.2.6, SciPy 1.15.3, CVXPY 1.7.5,
CLARABEL 0.11.1, PySCIPOpt 6.2.1, and SCIP 10.0.2. It can be recreated with
`python -m pip install -e ".[dev,reproducibility]"`.

```bash
python scripts/verify_four_effect_rational_lower.py
python scratch/d2_frontier/verify_low_weight_socp_exact_dual.py
python scratch/d2_frontier/verify_four_active_mccormick_exact_cover.py \
  scratch/d2_frontier/four_active_common_bias_fallback_exact_cover_l060.compact.json.gz \
  --target 0.76670
python scratch/d2_frontier/verify_global_frontier_l060.py
python -m pytest -q \
  tests/test_four_effect_rational_lower.py \
  tests/test_projective_tangent_interval_certificate.py \
  tests/test_ternary_projective_line_sensitivity.py \
  tests/test_four_active_compact_certificate.py \
  tests/test_global_frontier_assembly.py
```

The full ternary replay accepts the eight source shards as positional
arguments; it reconstructs all 12,008 selected cones under the new projective
premises. The public release bundle contains those shards, the compact
four-active tree, both projective covers, and SHA-256 checksums. The decisive
small manifest is `data/global_frontier_l060_exact_assembly.json`.

## 7. What is closed, and what is not

The earlier terminal-projectivity gap is closed at \(\lambda=0.6\): genuine
ternary and four-active readouts are bounded directly, and no projectivity
assumption enters equation (1). The result is solver independent in the
specific sense that replaying every stored dual and assembling the endpoint
requires no optimizer.

Three stronger claims remain unwarranted. The certificate does not give the
complete support curve, identify the exact maximizer at \(\lambda=0.6\), or
formalize every numerical kernel in Lean. The next mathematical frontier is
an analytic inequality or a sharper certified cover that collapses the
remaining width to the rational four-effect witness. A proof of equality
would be stronger than merely adding more digits to the present enclosure.

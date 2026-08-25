# State--Choi PSD/PPT cluster completion

Status: complete solver-conditional cover of one fixed base angular cell,
25 August 2026. This is not a solver-independent theorem and does not close
the full terminal leaf or terminal strip.

## What changed

The previous common-instrument relaxation lifted products
$P_{z\mu y}=x_{z\mu}J_y$ between an input Pauli coordinate and each outcome
Choi matrix. Entrywise McCormick envelopes, product trace preservation, and
PSD interval sandwiches made that relaxation convergent, but a difficult
scalar/Bloch branch still required 797 spatial nodes and ended only
$3.44\times10^{-6}$ below the target.

The stronger model retains the joint cone structure of all four Pauli
coordinates. For the subnormalised input state

\[
 \rho(x)=\frac12
 \begin{pmatrix}
 x_0+x_3 & x_1-i x_2\\
 x_1+i x_2 & x_0-x_3
 \end{pmatrix}
\]

and $J_y\succeq0$, an exact lifted point satisfies both

\[
 K_{zy}=
 \frac12
 \begin{pmatrix}
 P_{z0y}+P_{z3y} & P_{z1y}-iP_{z2y}\\
 P_{z1y}+iP_{z2y} & P_{z0y}-P_{z3y}
 \end{pmatrix}
 =\rho(x_z)\otimes J_y\succeq0
\]

and its input partial transpose

\[
 K_{zy}^{\Gamma}=
 \frac12
 \begin{pmatrix}
 P_{z0y}+P_{z3y} & P_{z1y}+iP_{z2y}\\
 P_{z1y}-iP_{z2y} & P_{z0y}-P_{z3y}
 \end{pmatrix}
 =\rho(x_z)^T\otimes J_y\succeq0.
\]

Both constraints are affine LMIs in the lifted products. Their validity is
immediate because the exact matrix is a product, hence separable across the
input-state/Choi split and positive under partial transpose. The implementation
uses realified $16\times16$ LMIs and adds both constraints for every input and
instrument outcome.

These tensor localizers are not a new generic optimization construction. The
project had already tested the state--Choi and PPT ideas in a global moment
relaxation, where the state--Choi version still gave `0.7633276193`. The new
technical step is to place the same exact cone information directly in the
localized input--Choi product model, together with shared trace rules, PSD
sandwiches, and an adaptive angular cover. That combination changes the
computational frontier.

## Quantitative ablation

On the former maximum `bbbb` spectral cell, the earlier product localizers
needed 49 spatial nodes and gave `0.7579083037`. Adding the state--Choi PSD
localizer closed the localized box at its root with `0.7563939115`. Adding the
PPT companion lowered the same root to `0.7497073209`.

The harder `bb-b` cell is more diagnostic:

| Relaxation | Solved spatial nodes | Upper bound |
|---|---:|---:|
| PSD sandwiches + product trace | 797 | `0.7579965603` |
| plus state--Choi PSD | 35 | `0.7579103294` |
| plus state--Choi PPT | 1 | `0.7529043854` |

Thus PSD alone reduced the tree by 95.6%, and the PPT companion removed the
spatial tree entirely. The exact-physical-product and DPP construction tests
cover both tensor signs.

## Angular cluster lemma

Solving all 2,216 open spectral cells independently would still waste most of
the gain. Suppose child cap $i$ has unit centre $n_i$ and angular radius
$\alpha_i$. Choose a unit parent centre $n$ and set

\[
 \beta=\max_i\{\arccos(n\mathbin{\cdot}n_i)+\alpha_i\}.
\]

The spherical triangle inequality shows that every direction in every child
cap is contained in the parent cap $(n,\cos\beta)$. Whenever
$\cos\beta>0$, this gives the valid convex constraint

\[
 \lVert r\rVert_2\leq \frac{n\mathbin{\cdot}r}{\cos\beta}.
\]

The implementation rounds the parent cosine downward. It begins with clusters
indexed by branch pattern and cube-face tuple. If a cluster root is not below
target, it splits the cap coordinate with the largest parent angular radius,
using the two farthest child centres as seeds. Only a discrete singleton would
be allowed to open a spatial tree; none was needed in the completed run.

Two probes illustrate the compression. One parent cap containing eight former
`bbbb` cells closed in one node at `0.7571985468`. A broad `bb-b` parent
containing caps 20 and 22 closed in one node at `0.7529210003`. A complete cube
face was too coarse (`0.7641378701`), and the adaptive tree split it rather than
spending a large spatial-node budget.

## Completed selected-cell result

The source spectral cover contains 21,190 cells. Of those, 2,216 had finite
bounds at or above `0.758`. The adaptive cover produced:

| Quantity | Value |
|---|---:|
| Initial face clusters | 66 |
| Total adaptive cluster nodes | 382 |
| Angular split nodes | 158 |
| Closed clusters | 224 |
| Source-open cells covered | 2,216 / 2,216 |
| Pending or unresolved clusters | 0 |
| Cap-containment audits | 1,282 |
| Largest closed cluster | 56 source-open cells |
| Maximum new cluster bound | `0.7579919715540884` |

All 224 closed clusters terminated at their product root; no spatial leaf
cover was required. The 382 cluster nodes incurred 12,988 solver calls: one
base solve, 32 coordinate-support solves, and one tensor-product solve per
node. This is about 5.8 times fewer cluster evaluations than the 2,216-cell
singleton route, and the number of terminal closures fell by a factor of 9.9.

Combining the new cover with the cells already below target in the source gives

\[
 U_{\mathrm{base\ angular\ cell}}
 =0.7579979090029844 < 0.758.
\]

The margin is `2.0909970e-6`. The limiting component is source cell 15,818,
which was already closed by the earlier spectral relaxation; the worst newly
closed cluster has margin `8.0284459e-6`.

## Numerical boundary

This is a complete combinatorial cover but remains solver-conditional. The
12,988 recorded solves comprise 12,599 `optimal`, 388 `optimal_inaccurate`,
and one `infeasible` status. Two child bounds exceed their parent bounds by
more than `1e-8`; the largest inversion is `9.62e-5`. Such inversions cannot
occur for exact optima and expose ordinary floating-point solver error. The
stored bounds include a `2e-6` additive safety allowance, but that does not
replace a checked dual or interval certificate.

Accordingly, the defensible result is:

> Conditional on the recorded CLARABEL solutions and safety convention, one
> fixed terminal box and one fixed base Fourier angular cell are completely
> below `0.758` after imposing common-instrument state--Choi PSD/PPT product
> localizers and an exhaustive adaptive cap-cluster cover.

It does not close the full terminal leaf, the 384-cell base Fourier cover, the
5,581-leaf terminal partition, or the global physical optimization problem.

## Reproduction

From `scratch/d2_frontier`, run the resumable cover:

```bash
python spectral_cap_cluster_cover.py \
  --frontier-json ternary_reconstructed_depth4_g2_top_leaf_bbb_p1_s92_l055.json \
  --output spectral_cap_cluster_cover_l055.json \
  --max-cluster-nodes 500 --leaf-max-nodes 100
```

Then generate and validate the compact audit:

```bash
python summarize_spectral_cap_cluster_cover.py \
  --source ternary_reconstructed_depth4_g2_top_leaf_bbb_p1_s92_l055.json \
  --checkpoint spectral_cap_cluster_cover_l055.json \
  --output spectral_cap_cluster_cover_l055_summary.json
```

The raw 7.8 MB checkpoint is regenerable and ignored by Git. The committed
summary records its SHA-256 digest, the complete pattern counts, limiting
cells, worst clusters, and numerical monotonicity audit.

The immediate next step is not another conceptual localizer. It is numerical
certification: re-solve source cell 15,818 and the ten worst cluster roots with
an independent solver/tighter tolerances, extract dual feasible points, and
replace the narrow floating-point margins by checked upper bounds. Only then
should the same cluster engine be expanded across the remaining base angular
cells of the terminal leaf.

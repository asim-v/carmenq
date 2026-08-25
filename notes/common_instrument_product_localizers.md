# Common-instrument conic-RLT localizers

Status: validated localized result, 25 August 2026. This note records a
mathematically sound strengthening and a solver-conditional numerical cover.
It does not claim a global frontier theorem or priority for the underlying
optimization construction.

## The obstruction

The localized ternary frontier model contains products between one input-state
Pauli coordinate and the Choi matrices of a single outcome-resolved quantum
instrument. Entrywise McCormick envelopes make the spatial relaxation
convergent under subdivision, but they forget two pieces of operator structure:
positive-semidefinite order and the fact that every outcome belongs to the same
trace-preserving instrument. The determinant and planar-Ando cuts developed
earlier constrain induced effects, yet a 1,000-node Ando run still had maximum
pending upper bound `0.7633726473409502`, above target `0.758`.

The missing connection is between global nonconvex optimization and quantum
process geometry. Reformulation--linearization is most often applied to scalar
products. Here the second factor is not an arbitrary vector of entries: it lies
in the positive-semidefinite Choi cone and its outcome sum has a fixed partial
trace. Preserving those facts after lifting produces a much stronger relaxation.

## Validity lemma

Let $x\in[\ell,u]$, let $J_y\succeq0$ be the Choi matrices of the outcomes
of one instrument, and suppose

\[
  \operatorname{Tr}_{\mathrm{out}}\sum_y J_y=I.
\]

Introduce exact product matrices $P_y=xJ_y$. Every physical point then
satisfies

\[
  P_y-\ell J_y=(x-\ell)J_y\succeq0,
  \qquad
  uJ_y-P_y=(u-x)J_y\succeq0,
\]

and

\[
  \operatorname{Tr}_{\mathrm{out}}\sum_y P_y=xI.
\]

The proof is immediate: the two scalar coefficients in the first line are
nonnegative, and the second identity follows by multiplying the common
trace-preservation equation by $x$. Consequently these constraints cannot
remove any exact physical strategy. In the complex implementation, each
Hermitian matrix $A$ is constrained through its realification

\[
  \mathcal R(A)=
  \begin{pmatrix}
  \operatorname{Re}A&-\operatorname{Im}A\\
  \operatorname{Im}A& \operatorname{Re}A
  \end{pmatrix},
\]

for which $A\succeq0$ if and only if $\mathcal R(A)\succeq0$.

These are cone-valued RLT constraints. The general optimization idea is known:
Anstreicher derives stronger lifted constraints from Kronecker products of
positive-semidefinite affine matrix pencils, explicitly as a generalization of
RLT and SOC-RLT
([SIAM Journal on Optimization, 2017](https://doi.org/10.1137/16M1078859)).
Conic strengthenings of RLT have also been studied computationally across
polynomial optimization
([González-Rodríguez et al., 2022](https://arxiv.org/abs/2208.05608)). The
present lemma is the particularly cheap scalar-cone specialization obtained by
multiplying $x-\ell\ge0$ and $u-x\ge0$ by the Choi PSD constraint.

## Why the two constraints reinforce each other

The trace identity alone is almost redundant: independent entrywise envelopes
can redistribute lifted mass among outcomes while preserving their sum. The
PSD sandwiches alone constrain each outcome in Löwner order but still allow
such a redistribution. Together they require every lifted outcome to lie in
the same matrix-order interval while their common partial trace is exactly
normalized. The observed gain is therefore a coupling effect, not the sum of
two independently strong cuts.

This also explains why the construction is specific enough to be useful.
Quantum bilinear optimization and convergent SDP hierarchies are established
subjects
([Berta, Fawzi, and Scholz, 2016](https://authors.library.caltech.edu/records/xe6ym-mqt54);
[constrained bilinear SDP hierarchies, 2021](https://doi.org/10.1007/s10107-021-01650-1)).
Quantum-instrument compatibility likewise has an independent formal literature
([Mitra and Farkas, 2022](https://doi.org/10.1103/PhysRevA.105.052202)). What the
present computation contributes is a small spatial SDP relaxation that keeps
the common-instrument normalization directly in the product lift, rather than
replacing the whole problem by a much larger generic moment hierarchy.

## Controlled ablation

The 20 leading cells of the prior Ando checkpoint were solved with identical
boxes, objectives, inherited witnesses, safety margin, and solver settings.
Two baseline cells were already infeasible, leaving 18 finite comparisons.

| Added constraints | Mean finite bound improvement | Strictly improved finite cells | Finite bounds below `0.758` | Other outcomes |
|---|---:|---:|---:|---|
| POVM product sums | `8.57e-8` | 1/18 | 0 | 2 infeasible |
| Choi product trace | `1.32e-10` | 0/18 | 0 | 2 infeasible |
| Choi PSD sandwiches | `1.3279e-2` | 18/18 | 14 | 2 infeasible |
| PSD sandwiches + product trace | `2.6589e-2` | 18/18 | 18 | 1 infeasible, 1 solver error |

The one solver error is not counted as a closure. The strict JSON artifact
`product_localizer_ablation_top20_l055.json` records every status and uses
`null` plus an explicit extended-real class instead of nonstandard `Infinity`
or `NaN` values.

The stronger evidence comes from complete branch-and-bound runs. With PSD
sandwiches but without the product-trace identities, 100 solved nodes left 95
pending cells and maximum pending upper bound `0.7631986663448878`. Enabling
both constraints closed the same localized region in 49 solved nodes: 24 split
nodes, 25 closed leaves, no pending or unresolved nodes, and maximum leaf upper
bound

\[
  0.7579083037237451 < 0.758.
\]

The numerical margin is `9.16962762549e-5`. The maximizing leaf had solver
status `optimal`; the full tree contains 34 `optimal` and 15
`optimal_inaccurate` solves. Every stored bound includes an additive `2e-6`
safety margin, but this remains a CLARABEL double-precision result rather than
an interval-arithmetic certificate. The validated compact artifact
`product_localizer_cover_l055_summary.json` contains the raw checkpoint hashes,
the complete leaf list, an independent parent/child audit, and the exact
comparison against the sandwich-only control.

## Epistemic status and originality boundary

The result is real but localized. It proves, conditional on the recorded
solver bounds, that one continuous terminal cell paired with the selected
Fourier spectral cell lies below the target. It does not yet cover all terminal
and spectral cells needed for a global statement.

The generic cone-RLT inequality is not new and must not be advertised as such.
A defensible contribution claim would instead be: a Choi-specialized conic-RLT
relaxation with shared product trace normalization, integrated into a spatial
cover for common quantum instruments, together with a quantitative example in
which the two individually weak/insufficient ingredients become target-closing
when combined. Priority for that exact combination still requires a broader
systematic search before publication.

## Reproduction

Run the component ablation from `scratch/d2_frontier`:

```bash
python product_localizer_ablation.py \
  --frontier-json ternary_reconstructed_depth4_g2_top_leaf_bbb_p1_s92_l055.json \
  --checkpoint ternary_exactdet_ando_guided_instrument_topcell_pilot100_l055.json \
  --limit 20 --output product_localizer_ablation_top20_l055.json
```

Run the target-closing cover:

```bash
python ternary_bilinear_instrument_input_cover.py \
  --frontier-json ternary_reconstructed_depth4_g2_top_leaf_bbb_p1_s92_l055.json \
  --localisation-json ternary_exactdet_hybrid_instrument_topcell_pilot_l055.json \
  --output ternary_exactdet_ando_matrixlocalizer_instrument_topcell_pilot100_l055.json \
  --top-spectral-cell --planar-ando-witnesses \
  --common-instrument-product-trace-rules \
  --common-instrument-product-psd-sandwiches \
  --max-nodes 100 --checkpoint-every 10
```

The raw production and sandwich-only checkpoints are regenerable and ignored by
Git. Their SHA-256 digests are committed in the compact summary. The next
scientific step is to apply the same localizers to the full terminal/spectral
partition, then replace the numerical SDP leaf bounds by independently checked
dual or interval certificates where the final global margin is small.

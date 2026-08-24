# Singular common-instrument stratum

## Result in one sentence

For the fixed interior benchmark at \(\lambda=0.55\), the complete singular
operator-basis stratum \(\det R=0\) has been excluded at score \(0.758\) by an
exact four-chart cover and 149 solver-infeasible terminal cells. This is a
solver-conditional numerical certificate, not an interval or rational proof.

## Why a separate treatment is necessary

For four linearly independent qubit inputs, the conditioned maps are uniquely
reconstructed by

\[
L_y=Q_y^{\mathsf T}(R^{\mathsf T})^{-1}.
\]

This gives the determinant-scaled Choi criterion described in
`operator_basis_instrument_criterion.md`. It deliberately says nothing at
\(\det R=0\), where the inverse does not exist and the output family is an
operator-system extension problem. Bounding \(|\det R|\) by a small number did
not repair the relaxation: the direct \(\det R=0\) solve retained a dual bound
near \(0.76337\). The singular geometry itself had to be exposed.

## Exact four-chart cover

Let the rows of \(R\) be the Pauli coordinates of the four subnormalised input
states. If \(R\) is singular, there is a nonzero real vector \(c\) satisfying

\[
c^{\mathsf T}R=0.
\]

Choose an index \(k\) for which \(|c_k|\) is maximal and rescale so that
\(c_k=1\). Every other coefficient then lies in \([-1,1]\). Consequently,

\[
R_k+\sum_{z\ne k}c_zR_z=0,
\qquad (c_z)_{z\ne k}\in[-1,1]^3,
\]

and the four values \(k=0,1,2,3\) cover every singular \(R\). Conversely, the
displayed relation makes \(R\) singular, so no nonsingular point is inserted.

If all conditioned outputs come from a common linear instrument, the same
relation must hold after every outcome map:

\[
Q^{(y)}_k+\sum_{z\ne k}c_zQ^{(y)}_z=0
\quad\text{for every }y.
\]

The model imposes these input and output identities together with the literal
shared-instrument Choi completion. An axis-aligned bisection of \([-1,1]^3\)
is an exact partition: it introduces no discretisation error and leaves every
unresolved cell explicitly present in the checkpoint.

## Fixed benchmark certificate

The benchmark uses terminal traces \((0.92,0.64,0.44,0)\), prefix order
\(0\ge1\ge2\ge3\), and the prior box

\[
\begin{split}
a_0&\in[0.296875,0.42596435546875],\\
a_1&\in[0.224609375,0.34832000732421875],\\
a_2&\in[0.15234375,0.258392333984375],\\
a_3&\in[0.1083984375,0.201324462890625].
\end{split}
\]

Each cell was solved as a target-feasibility problem at score \(0.758\), with
SCIP primal and dual feasibility tolerances \(10^{-9}\) and ten seconds per
node. A terminal cell was counted as closed only when SCIP returned
`infeasible`; time-limited nodes were bisected.

| pivot | solved nodes | infeasible terminal cells | maximum depth | solver time (s) | open volume |
|---:|---:|---:|---:|---:|---:|
| 0 | 9 | 5 | 4 | 61 | 0 |
| 1 | 15 | 8 | 6 | 98 | 0 |
| 2 | 101 | 51 | 10 | 717 | 0 |
| 3 | 169 | 85 | 11 | 1245 | 0 |
| **total** | **294** | **149** | **11** | **2121** | **0** |

No cell produced a feasible target incumbent. All four checkpoints report a
complete cover with zero open volume. The compact machine-readable summary is
`scratch/d2_frontier/basis_null_chart_summary_l055.json`; the full cell trees
are the four `basis_null_chart_cover_p*_l055.json` files.

## What this does and does not establish

The result closes the lower-dimensional singular stratum inside one fixed
terminal-POVM cell. Together with the determinant-sign symmetry, the only
remaining part of this interior benchmark is one regular orientation,
\(\det R>0\). Its exact determinant-scaled Choi condition is known, but the
current degree-four spatial relaxation still reports a root dual near
\(0.7634\), above \(0.758\).

Several valid regular-stratum strengthenings were tested. A determinant shell,
explicit bounded inverse, lifted adjugate, finite Choi-witness net, and one
adapted flagged trace-norm contraction did not close the gap. The scalar
branch of the adapted contraction did lower the dual to \(0.761919\), while
the vector branch remained near \(0.7634\). These are negative numerical
results and are not presented as theorems.

The present contribution is therefore precise: it supplies an exact compact
charting of the singular compatibility problem and a complete
solver-conditional exclusion of that charted domain. It does not yet close
the full regular interior, the terminal-POVM cover, or the global frontier.

## Reproduction

For one chart, run for example

```text
python scratch/d2_frontier/basis_null_chart_cover.py \
  --pivot 2 --seconds-per-node 10 --max-nodes 128 --max-depth 12 \
  --output scratch/d2_frontier/basis_null_chart_cover_p2_l055.json
```

The checkpoint is resumable with `--resume` and a larger `--max-nodes`. Unit
tests verify max-entry normalisation, rejection of nonsingular matrices, and
the exact volume identity for every binary cell partition.

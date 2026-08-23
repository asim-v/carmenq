# A prior-reserve cap for fully active terminal qubit readouts

**Date:** 23 August 2026
**Status:** analytic necessary bound incorporated into a complete
solver-conditional four-active sector enclosure at \(\lambda=0.6\);
solver-independent interval validation remains open.

Let an optimal minimum-error terminal measurement have `k=3` or `k=4`
nonzero rank-one effects

\[
E_s=w_s\Pi_s,\qquad \sum_sE_s=I_2.
\]

Write the Helstrom dual as `Y`, put `A=Tr Y`, and let its eigenvalue bias be
`t in [0,1]`.  If `x_s` is the projection of the Bloch vector of `Pi_s` on
the Bloch axis of `Y`, complementary slackness fixes the row and column of
the syndrome state `tau_s` supported on `Pi_s`.  Positivity of the remaining
two-by-two matrix then gives

\[
p_s=\operatorname{Tr}\tau_s
\ge A f_t(x_s),\qquad
f_t(x)=\frac{1+2tx+t^2}{2(1+tx)}.
\tag{1}
\]

The singular corner is assigned its exact boundary value
\(f_1(-1)=0\), because the corresponding support expectation of the
rank-one Helstrom operator vanishes.

POVM completeness implies

\[
\sum_sw_s=2,\qquad \sum_sw_sx_s=0.
\tag{2}
\]

Since the syndrome priors sum to one, equations (1)--(2) prove

\[
\boxed{
A\le
\left[
\min_{0\le t\le1}
\min_{-1\le x_s\le1,\ \sum_sw_sx_s=0}
\sum_sf_t(x_s)
\right]^{-1}.}
\tag{3}
\]

For fixed `t`, the function `f_t` is increasing and concave.  A minimum over
the box section in equation (3) is therefore attained at a vertex: all but
one of the `x_s` equal `+1` or `-1`.  The inner optimisation is consequently
a finite enumeration, followed by a one-dimensional minimisation in `t`.
`scratch/d2_frontier/active_readout_audit_cap.py` implements exactly this
reduction.

For the balanced four-effect weights `(1/2,1/2,1/2,1/2)`, equation (3) gives
`A <= 1/2`.  Hence at support weight `lambda=0.6` this whole terminal geometry
has support at most `0.6(1/2)+0.4=0.7`, even if RETURN is granted its algebraic
maximum.  At weights `(0.9,0.5,0.4,0.2)` the cap is approximately `0.588304`,
again enough to lie below `0.76591` with perfect RETURN.  Highly skewed
near-projective weights remain the difficult region; the cap is a pruning
theorem, not by itself the complete four-active converse.

The minimisation in equation (3) relaxes the transverse components of the
Bloch-vector closure and is therefore safely one-sided.  Any later use in a
global certificate must round the scalar minima outward or cover the
one-dimensional branches with interval arithmetic.

## 2. Role in the completed sector split

For \(\lambda=0.6\), the recorded scalar calculation bounds every fully
active readout with maximum effect trace at most \(0.88325\) by

\[
S\leq0.7658931806287275.
\]

Above that threshold, exact Helstrom reconstruction, transverse polygon
closure, and averaged projective comparisons give a compact spatial model.
SCIP bounds the region in which all four effect traces are at least
\(0.0003\) by \(0.7663946336432972\).  If the smallest trace lies below that
threshold, deleting that effect costs at most
\(0.6\times0.0003=0.00018\), so the complete ternary bound \(0.76643\) yields
\(0.76661\).  Together with the projective and ternary covers, this exhausts
all terminal arities and gives the reported outward decimal bound
\(\beta_{2\mathrm b}(0.6)\leq0.76662\).

The derivation and all artifact links are collected in
`notes/interleaved_interior_frontier_l060.md`.

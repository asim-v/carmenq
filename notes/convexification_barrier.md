# A convexification barrier for the two-block leaf program

**Date:** 23 August 2026
**Status:** exact elementary theorem and reproducible numerical
specialization; it diagnoses, but does not close, the nonconvex frontier

## 1. Statement

Fix a terminal POVM \(P=\{P_s\}_{s=1}^q\), let
\(w_s=\lVert P_s\rVert_\infty\), and suppose every syndrome occurs on \(m\)
path blocks. Put \(N=qm\). Consider the relaxation obtained after tracing
the emitted registers, replacing one pure Schmidt-rank-two leaf by an
arbitrary mixed state of Schmidt number at most two, and evaluating RETURN
from the mixed state's path probabilities:

\[
R(p)=\frac1N\left(\sum_{i=1}^N\sqrt{p_i}\right)^2.
\]

Then the relaxed support is exactly

\[
\boxed{
\beta_{\rm mix}(\lambda)
=\lambda_{\max}\!\left[
\lambda\,\operatorname{diag}(\widetilde w)
+\frac{1-\lambda}{N}\mathbf1\mathbf1^{\mathsf T}
\right],}
\tag{1}
\]

where every \(w_s\) is repeated \(m\) times in \(\widetilde w\). The
optimum is attained by a separable state diagonal in the path labels.
Therefore the value in equation (1) is unavoidable for every convex outer
approximation that contains the Schmidt-number-two cone and uses this same
concave path-probability RETURN formula.

## 2. Proof

For an arbitrary relaxed state, let \(p_i\) be its path probabilities. On
path \(i\), whose syndrome is \(s(i)\),

\[
\operatorname{Tr}(P_{s(i)}\rho_i)\leq w_{s(i)}p_i.
\]

Hence its score is no greater than

\[
\lambda\sum_i\widetilde w_i p_i
+\frac{1-\lambda}{N}\left(\sum_i\sqrt{p_i}\right)^2.
\tag{2}
\]

Set \(x_i=\sqrt{p_i}\). The constraints become \(x_i\geq0\) and
\(\lVert x\rVert_2=1\), while equation (2) is the Rayleigh quotient of the
matrix in equation (1). That matrix is entrywise nonnegative, so a leading
eigenvector can be chosen nonnegative by Perron--Frobenius. Its largest
eigenvalue is therefore attainable under the constraint \(x_i\geq0\).

For the reverse inequality, choose the path-diagonal state

\[
\rho=\sum_i p_i\,|z_i\rangle\!\langle z_i|
\otimes|y_i\rangle\!\langle y_i|
\otimes|\phi_{s(i)}\rangle\!\langle\phi_{s(i)}|,
\]

where \(P_s|\phi_s\rangle=w_s|\phi_s\rangle\). This state is separable
across the middle cut, has Schmidt number one, attains the AUDIT term in
equation (2), and has the prescribed path probabilities. Taking
\(p_i=x_i^2\) from the Perron eigenvector proves equality.

## 3. The fixed interior benchmark

For

\[
\lambda=0.55,\qquad
w=(0.92,0.64,0.44,0),\qquad m=4,
\]

the exact value is

\[
\boxed{\beta_{\rm mix}=0.790265609741724\ldots .}
\tag{3}
\]

This lies well above the target \(0.758\). Partial-transpose trace norm,
realignment, symmetric extensions, or even exact membership in the convex
Schmidt-number-two cone cannot repair the gap: the optimizer in the proof is
already separable. The indispensable missing condition is that all terms
come from one deterministic pure leaf before RETURN is evaluated.

This theorem explains the flatness failure seen in the state--instrument
moment hierarchy. It does not prove the desired upper bound for the pure leaf
problem. The exact replacement is the block-coherence functional in
[`coherence_preserving_convexification.md`](coherence_preserving_convexification.md).

## 4. Reproduction

Run

```text
python scratch/d2_frontier/convexification_barrier.py \
  --lambda 0.55 --effect-norms 0.92 0.64 0.44 0 \
  --output scratch/d2_frontier/convexification_barrier_l055.json
```

The script archives the Perron distribution, its AUDIT and RETURN values,
and the eigen-equation residual. The tests compare one hundred deterministic
random path distributions against equation (1).

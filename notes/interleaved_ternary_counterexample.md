# A Finite-Outcome Counterexample to the Two-Parameter Interleaved Candidate

**Date:** 20 August 2026<br>
**Status:** explicit reproducible lower bound; complete interior frontier still
open

## 1. What was falsified

The two-parameter QND construction in
`notes/interleaved_interior_candidate.md` gives, at balanced support weight,

\[
\beta_{\rm can}(1/2)=0.755437446228747\ldots .
\]

That construction remains a valid physical lower bound, but it is not the
unrestricted streamed frontier.  A complete adaptive instrument with three
fine-grained outcomes per temporal node gives

\[
P_{\rm A}=0.625754561820388\ldots,
\qquad
F_{\rm R}=0.893143378814133\ldots,
\]

and therefore

\[
\boxed{
\frac{P_{\rm A}+F_{\rm R}}2
=0.759448970317260\ldots
>\beta_{\rm can}(1/2).
}
\tag{1}
\]

The strict excess is (0.004011524088513\ldots), far larger than the
reported numerical residuals.  Consequently, the conjectured equality

\[
\beta^{\rm stream}_{H_{\rm I},2}(\lambda)=\beta_{\rm can}(\lambda)
\]

is false.  None of the exact endpoint or robust full-crossing theorems depend
on that conjecture.

## 2. Frozen interface and construction class

The counterexample uses the same four-slot interleaved check matrix

\[
H_{\rm I}=\begin{pmatrix}1&0&1&0\\0&1&0&1\end{pmatrix}.
\]

At every slot, a complete local instrument maps the fresh carrier and one
persistent qubit to one emitted carrier, the updated qubit, and a classical
outcome.  Emitted carriers are sequestered immediately.  After the fourth
slot, the classical transcript and terminal qubit are available to AUDIT;
RETURN may additionally act jointly on all sequestered carriers.  The stored
strategy has three local outcomes, 81 possible four-symbol transcripts, and
16 transcripts of probability greater than (10^{-8}).

The phrase “three outcomes” describes a compact representation, not a new
physical resource assumption.  A finite-outcome instrument can be decomposed
into binary sub-instruments inside a slot.  The scientific point is that the
earlier one-bit-per-slot variational ansatz did not cover the declared
finite-outcome model.

The optimizer was allowed arbitrary complex non-QND local maps.  The resulting
leaf effects (K_c^\dagger K_c) are diagonal to numerical precision, while
their polar isometries are genuinely non-QND.  Thus the gain is not explained
by a more complicated classical likelihood table alone: coherent output
geometry matters.

## 3. Framework-neutral verification

The file `data/interleaved_ternary_counterexample.npz` stores the local Kraus
operators and AUDIT effects after normalization.  It contains no optimizer
state and requires no PyTorch deserialization.  The independent command

```text
python scripts/verify_interleaved_counterexample.py
```

contracts all 81 transcript branches with NumPy, checks every local
completeness relation, checks positivity and completeness of every terminal
POVM, evaluates AUDIT directly, and evaluates exact flagged RETURN from the
singular values of each branch.  The frozen artifact gives

\[
\max_v\|V_v^\dagger V_v-I_4\|_F
=1.58\times10^{-15},
\]

\[
\max_c\left\|\sum_sQ_{s|c}-I_2\right\|_F
=1.68\times10^{-15},
\]

and total transcript probability (1-1.1\times10^{-16}).  The smallest
computed POVM eigenvalue is (-1.7\times10^{-16}), consistent with floating
point roundoff.

## 4. Adversarial checks

Five independent ternary searches were run for 1,800 optimization steps at
(lambda=1/2).  Two random seeds converged to the same nontrivial score within
(2.5\times10^{-11}); the others converged to the exact no-record point.  An
independent NumPy contraction reproduced the stored POVM score. A later
multi-start qubit-dual calculation gives the same AUDIT value at the displayed
precision; the earlier report of a slightly larger dual value was a local
optimizer failure and is withdrawn.

The ternary strategy was embedded into complete four- and five-outcome trees
and reoptimized.  Neither embedding improved equation (1) at the displayed
precision.  This is a falsification check, not an outcome-cardinality theorem.

Every active refined leaf has the same normalized singular spectrum and the
same homogeneous score to numerical precision.  The 16 leaves collapse into
four likelihood orbits of total probability (1/4), related by flips of the
first and third input bits.  This covariance is evidence that the improvement
is structural rather than a single postselected anomaly.

## 5. What is known after the counterexample

The exact arbitrary-instrument statements remain:

* the grouped order attains the static rank-two frontier;
* for the interleaved order, (P_{\rm A}=1) implies
  (F_{\rm R}\leq1/4), and the bound is attained;
* the linear-tail theorem certifies a strict grouped-versus-interleaved support
  gap for \(3/7<\lambda<1\), including balanced weight.

For the balanced interior, the explicit lower bound is now equation (1).  A
general single-leaf tensor-train relaxation reaches

\[
0.759802783851444\ldots,
\tag{2}
\]

with ((P_{\rm A},F_{\rm R})=(0.620085075586\ldots,
0.899520492117\ldots)).  Equation (2) is not a certified upper bound because
its nonconvex variational maximum has not been globally certified.  It does,
however, leave only (3.54\times10^{-4}) between the best complete strategy
and the best located postselected-leaf relaxation.

Continuation searches also reveal a weak-information branch absent from the
two-parameter family. A homogeneous active leaf located at
\((P_{\rm A},F_{\rm R})\approx(0.57217925,0.94295569)\) crosses the no-record
support at \(\lambda\approx0.44143895\). Complete-instrument continuation is
consistent with that onset, but neither the decimal nor a closed form is
proved. The earlier guess \(4/9\) is withdrawn.

## 6. Revised theorem gate

The complete interior problem is not solved.  A correct solution now requires
one of the following:

1. a causal upper certificate matching the finite-outcome construction;
2. a stronger complete strategy, possibly obtained by completing the optimal
   single-leaf tensor; or
3. a proof that a declared separability/comb hierarchy converges tightly
   enough to separate the two values in equations (1) and (2).

Until then, equation (1) is the public achievable benchmark, not a claimed
frontier.  This negative result is scientifically material: it identifies the
precise hidden restriction in the previous ansatz and prevents a false
optimality claim from entering the manuscript.

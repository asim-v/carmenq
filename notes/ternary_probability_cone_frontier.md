# The ternary probability-cone reduction

**Date:** 23 August 2026
**Status:** exact reduction plus a completed finite numerical cover of the
full sorted ternary terminal-weight strip at \(\lambda=0.6\); outward-rounded
dual validation remains open.

## 1. Why this reduction matters

The earlier fixed-readout models reconstructed Bloch coordinates by inverting
a terminal measurement matrix.  That representation becomes ill-conditioned
at the binary-projective boundary and obscures which constraints are genuinely
necessary.  A ternary rank-one qubit POVM admits a simpler description entirely
in its own probability coordinates.  In those coordinates, positivity of the
four syndrome operators and optimality of the terminal discrimination
measurement are Lorentz-cone conditions.  No inverse matrix and no hidden
normal Bloch coordinate are required.

This reduction is useful independently of the present audit--return problem:
it converts the KKT system for minimum-error discrimination by a fixed planar
ternary qubit POVM into a collection of copies of one homogenized ellipse cone.

## 2. Rank-one ternary POVMs and Horwitz parameters

Let

\[
E_t=w_t\Pi_t,\qquad t=0,1,2,
\]

be a rank-one qubit POVM, with the effects relabelled so that
\(w_0\geq w_1\geq w_2\).  Write its reciprocal inellipse parameters as
\(A,B\in[1,2]\).  The effect traces are

\[
w_0=\frac{A}{A+B-1},\qquad
w_1=\frac{B}{A+B-1},\qquad
w_2=\frac{A+B-2}{A+B-1}.
\tag{1}
\]

For a positive, possibly subnormalised operator \(X\), define its terminal
probability vector

\[
q_t=\operatorname{Tr}(E_tX),\qquad p=\sum_tq_t=\operatorname{Tr}X.
\]

The vector \(q\) belongs to the homogenized inellipse cone precisely when it
is nonnegative and

\[
B^2q_0^2+A^2q_1^2
-2(AB-2A-2B+2)q_0q_1
-2Bq_0p-2Aq_1p+p^2\leq0.
\tag{2}
\]

Equation (2) is one second-order-cone constraint whenever the ellipse is
nondegenerate.  Its degenerate limits give the corresponding line-segment
probability cones.

## 3. Probability-cone positivity lemma

Let \(\mathcal C_E\) denote the cone defined above.  Then

\[
\mathcal C_E
=\left\{
(\operatorname{Tr}E_0X,\operatorname{Tr}E_1X,
\operatorname{Tr}E_2X):X\succeq0
\right\}.
\tag{3}
\]

To see this, choose the Bloch plane spanned by the three effect directions.
Removing the component of \(X\) normal to that plane is a positive dephasing
map and leaves all three terminal probabilities unchanged.  Positivity of the
remaining planar operator is exactly the disk whose affine image is the
inellipse.  Thus membership in \(\mathcal C_E\) is neither merely necessary nor
a relaxation: it is equivalent to the existence of a positive qubit operator
with those probabilities.

## 4. Helstrom optimality without Bloch reconstruction

Let \(q_s\in\mathcal C_E\) be the terminal probability vector of syndrome
operator \(\tau_s\), for \(s=0,1,2,3\).  The POVM \(E\), with a zero fourth
effect, is minimum-error optimal for this ensemble if and only if there is a
vector \(h\) such that

\[
h\in\mathcal C_E,\qquad h-q_s\in\mathcal C_E\quad\text{for every }s,
\tag{4}
\]

and

\[
\sum_{s=0}^{2}q_{s,s}=\sum_{t=0}^{2}h_t.
\tag{5}
\]

Indeed, \(h\) is the terminal-probability representation of a Helstrom dual
operator \(Y\succeq0\), while \(h-q_s\in\mathcal C_E\) is equivalent to
\(Y-\tau_s\succeq0\).  The left side of equation (5) is the success
probability of \(E\), and the right side is \(\operatorname{Tr}Y\).  Weak
duality plus equality therefore proves optimality.  Conversely, the usual
Helstrom dual supplies
such an \(h\); planar dephasing preserves every inequality and terminal
probability.

This is the central simplification used by
`scratch/d2_frontier/ternary_probability_cone_cover.py`.

## 5. A projective comparison lemma

Fix distinct labels \(i,j,k\in\{0,1,2\}\).  Compare the ternary POVM with the
binary-projective readout

\[
Q_i=\Pi_i,\qquad Q_j=I-\Pi_i,\qquad Q_k=0.
\]

If \(p_s=\operatorname{Tr}\tau_s\), direct use of
\(E_j+E_k=I-w_i\Pi_i\) gives

\[
A_3-A_2
\leq (1-w_i)p_j+w_kp_k.
\tag{6}
\]

The RETURN score is unchanged by this terminal replacement.  Hence, if
\(\beta_{\rm proj}(\lambda)\) is any valid upper bound for the binary-projective
sector,

\[
S_3\leq\beta_{\rm proj}(\lambda)
+\lambda\big[(1-w_i)p_j+w_kp_k\big].
\tag{7}
\]

All six ordered choices of \((i,j)\) are valid simultaneously.  The continuous
SOCP cover imposes their concave McCormick envelopes on every terminal-weight
box.  At \(w_i=w_j=1,w_k=0\), equation (7) reduces exactly to the independently
covered projective frontier and removes the spurious degeneracy that appears
in a probability-only relaxation.

## 6. Continuous box relaxation

For every path, write \(q_t=w_tu_t\), where \(u_t\) is the probability of the
rank-one projector \(\Pi_t\).  A terminal \((A,B)\) box supplies exact intervals
for \(w_t\).  One global weight vector is shared by all paths, obeys
\(\sum_tw_t=2\), and is coupled to the Horwitz box through the linear relations

\[
A_{\min}(1-w_2)\leq w_0\leq A_{\max}(1-w_2),
\]

and the analogous inequalities for \((B,w_1)\).  Four McCormick inequalities
outer-convexify every product \(q_t=w_tu_t\).  These constraints converge to
equality as the terminal box shrinks.

The remaining common-qubit requirement is imposed through clean-POVM
inellipse constraints on selected compatible effect pairs.  Each parameter
box is an SOCP outer relaxation, so a finite cover below a target is a global
certificate for the stated outer model rather than a collection of local
optimisations.

## 7. Completed computational cover

The fixed-weight adaptive cover at
\((w_0,w_1,w_2)=(0.99,0.80,0.21)\) was the first difficult local test.  The
final continuous calculation supersedes it by covering the entire sorted
ternary weight strip.  The single-member archive
`scratch/d2_frontier/continuous_terminal_projective_l055cert_complete.json.tar.gz`
records 12,008 terminal boxes, 24,002 solved SOCP nodes, no open boxes, and
maximum leaf bound

\[
0.7664281458427126<0.76643.
\]

The cover incorporates the exact projective comparison in equation (7) and
independent projective support lines at \(\lambda=0.55\) and \(0.6\).  It
therefore remains stable at the degenerate projective edge instead of relying
on an ill-conditioned inverse-POVM chart.

The word *numerical* is essential.  CLARABEL feasibility and dual values were
computed at explicit tolerances and every node receives a conservative score
margin, but the conic duals have not yet been outward-rounded with interval
arithmetic.  The result is strong computational evidence and a replayable
solver certificate, not yet a solver-independent formal proof.

## 8. What remains

The ternary sector itself is now closed at this support direction, subject to
the numerical qualifier above.  A separate projected Helstrom spatial bound
and a small-effect deletion argument also close the four-active sector, giving
the combined enclosure

\[
0.7658988152646944
\leq\beta_{2\mathrm b}(0.6)
\leq0.76662.
\]

The complete exhaustion is documented in
`notes/interleaved_interior_frontier_l060.md`.  What remains here is
solver-independent outward validation of the conic duals, not an uncovered
terminal-weight region.

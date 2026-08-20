# Analytic Candidate for the Interleaved Interior Frontier

**Date:** 20 August 2026<br>
**Status:** exact achievable family and exact reduced optimization; global
arbitrary-instrument converse still open

## 1. Problem

For

\[
H_{\rm I}=
\begin{pmatrix}
1&0&1&0\\
0&1&0&1
\end{pmatrix},
\]

the exact streamed support function with one persistent coherent qubit is

\[
\beta^{\rm stream}_{H_{\rm I},2}(\lambda)
=\sup_{\mathcal S}
\bigl[\lambda P_{\rm A}(\mathcal S)
+(1-\lambda)F_{\rm R}(\mathcal S)\bigr].
\]

The perfect-AUDIT endpoint and a robust interval next to it are already
proved for arbitrary adaptive non-QND instruments.  This note isolates an
analytic two-parameter family that exactly reconstructs every nontrivial
variational optimum found so far.  It does not silently promote the numerical
agreement to a global theorem.

## 2. Canonical causal instrument

Only slots 1 and 3 produce classical flags, denoted by \(c_1,c_3\).  Slots 2
and 4 perform controlled qubit updates.  Fix \(q,v\in[0,1]\), and define

\[
p=1-(1-q)v^2,
\qquad
u=\frac{qv^2}{p}.
\tag{1}
\]

The slot-1 outcome guesses \(x_1\).  Conditional on a correct outcome, which
has probability \(p\), the persistent qubit is prepared in

\[
|\psi\rangle=\sqrt u\,|0\rangle+\sqrt{1-u}\,|1\rangle.
\]

Conditional on a wrong outcome, it is prepared in \(|0\rangle\).  Slot 2
applies \(Z^{x_2}\).  At slot 3, after the known correction determined by
\(x_3\) and the classical flags, the correct-effect operator is

\[
E_{\rm c}=q|0\rangle\!\langle0|+|1\rangle\!\langle1|,
\]

and the wrong-effect operator is

\[
E_{\rm w}=(1-q)|0\rangle\!\langle0|.
\]

Their positive square roots are used as Kraus operators.  Slot 4 applies the
final \(Z^{x_4}\) update.  Bit-flip covariance supplies the other outcome
branches, so the construction is a complete QND instrument rather than a
postselected filter.

## 3. Exact score formulas

For record errors \(e_1=x_1\oplus c_1\) and
\(e_3=x_3\oplus c_3\), every leaf has the likelihood table

\[
\begin{array}{c|cc}
 &e_3=0&e_3=1\\\hline
e_1=0&A&B\\
e_1=1&B&C
\end{array},
\tag{2}
\]

where

\[
A=1-(1-q^2)v^2,
\qquad
B=q(1-q)v^2,
\qquad
C=(1-q)^2v^2.
\tag{3}
\]

Thus \(A+2B+C=1\).  The terminal qubit carries useful information about the
second syndrome bit only in the favoured \(A\) sector.  Its unnormalised
Helstrom bias is

\[
D=\frac{2B}{C}\sqrt{AC-B^2}.
\tag{4}
\]

Equation (4) follows directly by evaluating the two pure states obtained from
\(E_{\rm c}^{1/2}Z^{x_2}|\psi\rangle\).  The remaining sectors either give an
incorrect first syndrome bit or no second-bit bias.  Consequently,

\[
P(q,v)=\frac12+qv\sqrt{1-v^2}-q(1-q)v^2,
\tag{5}
\]

while optimal flagged polar recovery gives

\[
F(q,v)=\frac14\left[
\sqrt{1-(1-q^2)v^2}
+v\bigl(1-q+2\sqrt{q(1-q)}\bigr)
\right]^2.
\tag{6}
\]

The same family can be written using \(q\leq p\leq1\):

\[
A=p(1+q)-q,
\quad B=q(1-p),
\quad C=(1-q)(1-p),
\]

\[
P=\frac12-q(1-p)
+\frac{q}{1-q}\sqrt{(1-p)(p-q)},
\tag{7}
\]

\[
F=\frac14\left[
\sqrt{p(1+q)-q}
+2\sqrt{q(1-p)}
+\sqrt{(1-q)(1-p)}
\right]^2.
\tag{8}
\]

## 4. Candidate support function and phase transition

The exact support within this family is the two-variable maximisation

\[
\beta_{\rm can}(\lambda)
=\max\left\{
1-\frac\lambda2,
\max_{0\leq q,v\leq1}
\bigl[\lambda P(q,v)+(1-\lambda)F(q,v)\bigr]
\right\}.
\tag{9}
\]

The first term is the distinct no-record strategy
\((P_{\rm A},F_{\rm R})=(1/2,1)\); it is not a singular point of equations
(5)--(6).  The two branches coexist at

\[
\lambda_{\rm c}=0.477812793357157\ldots,
\tag{10}
\]

where the nontrivial stationary point is

\[
q_{\rm c}=0.576397215294735\ldots,
\qquad
v_{\rm c}=0.812012949137950\ldots,
\]

\[
(P_{\rm c},F_{\rm c})
=(0.612174911465904\ldots,
0.897357485763005\ldots).
\tag{11}
\]

The jump from \((1/2,1)\) to equation (11) is first order.  On the smooth
branch, \((q,v)\) obeys

\[
\det\frac{\partial(P,F)}{\partial(q,v)}=0,
\tag{12}
\]

and the support weight is recovered from

\[
\frac{\lambda}{1-\lambda}
=-\frac{\partial_qF}{\partial_qP}
=-\frac{\partial_vF}{\partial_vP}.
\tag{13}
\]

At balanced weight,

\[
q=0.616895603071868\ldots,
\qquad
v=0.800317703643181\ldots,
\]

and

\[
\beta_{\rm can}(1/2)
=0.755437446228747\ldots.
\tag{14}
\]

As \(h=1-\lambda\downarrow0\), the optimizer satisfies

\[
1-q=\frac{h^2}{8}+O(h^3),
\qquad
v=\frac1{\sqrt2}+O(h^2),
\]

and

\[
\beta_{\rm can}(1-h)
=1-\frac{3h}{4}+\frac{h^2}{8}+O(h^3).
\tag{15}
\]

This approaches the exact interleaved endpoint
\((P_{\rm A},F_{\rm R})=(1,1/4)\) with the expected quadratic support gain
over the endpoint chord.

## 5. Verification and failed converse relaxation

For every stored nontrivial QND checkpoint from
\(\lambda=.48\) through \(.80\), equations (3), (5), and (6) reproduce the
direct tensor contraction to floating-point precision.  The Jacobian in
equation (12) is between \(10^{-16}\) and \(10^{-13}\) on the well-converged
checkpoints.  The unrestricted adaptive binary-tree searches at
\(\lambda=.5,.9,.99\) collapse to the same QND transducer.

A tempting converse is false.  If local instrument completeness is dropped
and a single arbitrary TT-rank-two leaf is optimised homogeneously, its score
at balanced weight exceeds equation (14).  Such a leaf need not admit a
causal completion with the same one-qubit bond.  Therefore TT rank alone does
not prove the frontier: the missing theorem must use local completeness at
the intermediate cuts.  This failed relaxation is useful evidence that the
order effect is genuinely about causal instrument completion, not merely a
static tensor-rank constraint.

## 6. Exact status and theorem gate

Equations (1)--(15) are exact for an explicit physical family.  They provide
the strongest current lower bound and replace the previous 320-real-parameter
description by a transparent two-parameter instrument.  They do **not** yet
establish

\[
\beta^{\rm stream}_{H_{\rm I},2}(\lambda)=\beta_{\rm can}(\lambda)
\]

for arbitrary finite-outcome adaptive non-QND instruments.  A complete
solution requires either:

1. a dynamic converse using local completeness and the qubit bond at every
   temporal cut; or
2. a counterexample that exceeds equation (9).

The current extensive searches support the equality, but numerical agreement
is not a substitute for that theorem.

## 7. Closest proof technology

Mohan, Tavakoli, and Brunner, *New Journal of Physics* **21**, 083034
(2019), DOI
[10.1088/1367-2630/ab3773](https://doi.org/10.1088/1367-2630/ab3773),
prove a tight two-score frontier for a sequential \(2\to1\) QRAC.  Their
argument uses polar decomposition, extremality of qubit instruments, an
eigenvalue optimisation over the downstream measurement, and a Bloch-sphere
reduction to antipodal pure preparations.  Those ingredients are directly
relevant to the missing converse here.

The result is not an exact reduction of the present problem.  Their second
score is another random-access decoding probability on the same travelling
qubit.  Here the second score is transcript-conditioned recovery of four EPR
pairs; previous outputs are sequestered; the same prefix must support a late
AUDIT/RETURN choice; and local completeness must hold across four temporal
slots.  In particular, substituting their sequential-QRAC curve for
equation (9) would be unjustified.  The appropriate next proof attempt is to
adapt their extremal-instrument/Bloch reduction while retaining the full
flagged entanglement-recovery functional.

# Exact projective sector of the two-block qubit frontier

**Date:** 22 August 2026
**Status:** theorem and finite numerical certificate for the
binary-projective terminal sector; the unrestricted readout is separately
enclosed, solver-conditionally, at \(\lambda=0.6\)

## 1. Why the two-block model is an upper relaxation

Group the four interleaved inputs into the two full-rank blocks

\[
z=(x_1,x_2),\qquad y=(x_3,x_4),\qquad z,y\in\mathbb F_2^2.
\]

Every normalised four-site Choi MPS with bond at most two has Schmidt rank at
most two across the paired cut

\[
(A_1B_1A_2B_2)\mid(A_3B_3A_4B_4M).
\]

Dropping the two internal one-bit cuts therefore enlarges the feasible set to
all two-block leaves of paired Schmidt rank at most two.  If
\(\beta_{4\rm s}(\lambda)\) denotes the target four-slot support and
\(\beta_{2\rm b}(\lambda)\) the relaxed two-block support, then

\[
\beta_{4\rm s}(\lambda)\leq\beta_{2\rm b}(\lambda).
\tag{1}
\]

The relaxation is useful because the four-effect leaf remains feasible and
has exactly the same value.  Thus a sharp converse for the two-block model
would also close the original four-slot frontier.

Put a two-block leaf in right-canonical form.  There are prefix vectors
\(|\psi_z\rangle\in B_L\otimes\mathbb C^2\) and continuation maps

\[
L_y:\mathbb C^2\longrightarrow B_R\otimes\mathbb C^2
\]

such that

\[
K|z,y\rangle=(I_{B_L}\otimes L_y)|\psi_z\rangle,
\qquad
\sum_yL_y^\dagger L_y=I_2.
\tag{2}
\]

Let

\[
\rho_z=\operatorname{Tr}_{B_L}|\psi_z\rangle\!\langle\psi_z|,
\quad
Q_y=L_y^\dagger L_y,
\quad
p_{zy}=\operatorname{Tr}(Q_y\rho_z).
\tag{3}
\]

The normalisation is \(\sum_z\operatorname{Tr}\rho_z=1\).  Computational
pinching of \(K^\dagger K\) gives the universal RETURN upper bound

\[
R(K)\leq\frac1{16}
\left(\sum_{z,y}\sqrt{p_{zy}}\right)^2.
\tag{4}
\]

No QND, real-amplitude, orthogonal-column, or measure-and-prepare assumption
has entered equations (1)--(4).

## 2. Binary terminal readout produces one affine line

Assume now that the terminal AUDIT POVM has exactly two nonzero effects and
that they are complementary rank-one projectors.  Its supported syndrome
labels form a two-point affine line \(T=\{t_0,t_1\}\).  Translate the answer
labels so that \(t_0=0\) and write \(t=t_1-t_0\neq0\).  If the terminal
projectors are \(P_0,P_1\), define

\[
G_{y,b}=L_y^\dagger(I_{B_R}\otimes P_b)L_y,
\qquad b\in\{0,1\}.
\tag{5}
\]

They obey

\[
G_{y,0}+G_{y,1}=Q_y,
\qquad
\sum_yQ_y=I_2.
\tag{6}
\]

For prefix label \(z\), the two successful suffixes are \(y=z\) and
\(y=z+t\).  Hence its correct-decision effect is

\[
H_z=G_{z,0}+G_{z+t,1}.
\tag{7}
\]

The four effects \(H_z\) form a POVM.  More is true.  Let \(C\) be either
coset of the subgroup \(\{0,t\}\).  Since addition by \(t\) permutes the two
elements of \(C\),

\[
\sum_{z\in C}H_z
=\sum_{y\in C}(G_{y,0}+G_{y,1})
=\sum_{y\in C}Q_y.
\tag{8}
\]

Calling the first coarse effect \(E\), equation (8) proves that every legal
binary-projective leaf induces exactly the binary-line geometry

\[
H_0+H_1=E,\qquad H_2+H_3=I-E
\tag{9}
\]

after relabelling within the two cosets.

## 3. The RETURN reduction is sharp

Apply Cauchy--Schwarz separately to the two suffixes in each coset.  For
every prefix state,

\[
\sum_{y\in C}\sqrt{\operatorname{Tr}(Q_y\rho_z)}
\leq
\sqrt{2\operatorname{Tr}(E\rho_z)},
\tag{10}
\]

and the complementary coset gives the same expression with \(I-E\).  Write
\(a_z=\operatorname{Tr}\rho_z\), \(\widehat\rho_z=\rho_z/a_z\), and use

\[
q_z=
\begin{cases}
\operatorname{Tr}(E\widehat\rho_z),&z\in C,\\
\operatorname{Tr}((I-E)\widehat\rho_z),&z\notin C,
\end{cases}
\qquad
c_z=\sqrt{q_z}+\sqrt{1-q_z}.
\tag{11}
\]

Equations (4) and (10) imply

\[
R(K)\leq\frac18\left(\sum_z\sqrt{a_z}\,c_z\right)^2.
\tag{12}
\]

The AUDIT score is

\[
A(K)=\sum_za_zd_z,
\qquad
d_z=\operatorname{Tr}(H_z\widehat\rho_z).
\tag{13}
\]

For fixed effects and normalised states, put \(v_z=\sqrt{a_z}\).  The exact
prior optimisation is the Perron eigenvalue

\[
\max_{\|v\|_2=1}
v^{\mathsf T}\left[
\lambda\operatorname{diag}(d_z)
+\frac{1-\lambda}{8}cc^{\mathsf T}
\right]v.
\tag{14}
\]

The two inequalities used above are sharp for the explicit 3E and 4E leaves.
Split each relevant coarse outcome into equiprobable suffixes, use square-root
Kraus operators for the fine rank-one effects, and use the emitted carriers to
orthogonalise the input columns.  The endpoint cases follow by deleting or
coalescing a fine effect.  For an arbitrary scalar topology, equations
(12)--(18) are used conservatively as an upper reduction; equality is needed
only at the exposed lower-bound construction.  This distinction avoids
claiming that every abstract scalar point has already been lifted to a
four-slot leaf.

## 4. Finite scalar topology theorem

Diagonalise

\[
E=\operatorname{diag}(x,y),\qquad1\geq x\geq y\geq0.
\tag{15}
\]

For a trial support value \(L>\lambda\max_zd_z\), the rank-one determinant
lemma turns equation (14) into the secular condition

\[
\frac{1-\lambda}{8}
\sum_z\frac{c_z^2}{L-\lambda d_z}\leq1.
\tag{16}
\]

Fix \(E\) and one binary split \((G,E-G)\).  After maximising each summand
over its state, the left side of equation (16) is a convex function of
\(G\in[0,E]\): the reciprocal of the positive affine denominator is convex,
and a pointwise supremum preserves convexity.  A maximum is consequently
attained at an extreme point of the operator interval.  Those points are

\[
G=E^{1/2}PE^{1/2},
\tag{17}
\]

where \(P\) is a projection on the support of \(E\).  In dimension two there
are only two inequivalent split types: an endpoint \((E,0)\), or a pair of
congruence-rank-one effects obtained from a rank-one \(P\).

At fixed coarse expectation \(q\), every summand increases with the correct
expectation \(d\).  A linear functional on a qubit Bloch slice is maximised
at a pure state.  Thus mixed prefix states are unnecessary.  After equation
(15), one projector angle and one state coordinate per active effect are
enough.  The complete sector is the maximum over four finite scalar
topologies:

\[
(\mathrm{endpoint},\mathrm{endpoint}),\quad
(\mathrm{endpoint},\mathrm{rank}),\quad
(\mathrm{rank},\mathrm{endpoint}),\quad
(\mathrm{rank},\mathrm{rank}).
\tag{18}
\]

The implementation in `scratch/d2_frontier/oneway_exact_topologies.py`
evaluates equation (18).  The exposed numerical winners found so far are
the endpoint/rank 3E phase and the symmetric rank/rank 4E phase, in addition
to the no-record point.

At \(\lambda=0.6\), the remaining scalar boxes have now also been covered by
SCIP spatial branch and bound.  The secular formulation removes the four
prior amplitudes and asks directly whether

\[
\frac{1-\lambda}{8}
\sum_z\frac{c_z^2}{L-\lambda d_z}\leq1.
\tag{19}
\]

At the fixed test level \(L=0.76591\), a trace-only outer relaxation certifies
most boxes and the full polynomial model certifies every unresolved child.
The machine-checked manifest reports 312 leaf boxes over the four topologies,
with maximum scaled dual bound \(0.9999905164\).  Therefore

\[
0.765898815264694\ldots
\leq\beta_{\rm projective}(0.6)
\leq0.76591.
\tag{20}
\]

The interval width is \(1.1185\times10^{-5}\).  Run
`scratch/d2_frontier/validate_projective_cover.py` to validate the hierarchy
and regenerate its compact summary.  This is a global numerical certificate
conditional on SCIP and its recorded tolerances; it is not a directed-rounding
interval proof independent of the solver.

## 5. What the new reduction does and does not close

At \(\lambda=0.6\), unrestricted complex optimisation of the two-block leaf
returns

\[
(A,R,L)=(0.8699300229\ldots,
0.6098520038\ldots,
0.7658988153\ldots),
\]

the same support as the explicit four-site 4E construction.  Its terminal
POVM has two unit-trace projectors and two effects below \(10^{-102}\).
At \(\lambda=0.5\), seeding the two-block relaxation with the four-site 3E
leaf reproduces \(0.759802783851444\ldots\); at \(\lambda=0.9\), random
two-block restarts reproduce \(0.9272008537\ldots\), within optimisation
tolerance of 4E.

Fixed symmetric nonprojective tests do not reveal a competing phase at
\(\lambda=0.6\): a trine readout reaches \(0.6972087779\ldots\), and a
tetrahedral readout reaches the no-record value \(0.7\).  A corrected trace
barrier diagnostic shows why nonzero parametrised effect traces are not, by
themselves, evidence for a genuine nonprojective optimum: after eliminating
the readout by the exact SOCP, the barrier stages can already be two-active.
The earlier trace-only interpretation is therefore withdrawn.

A more targeted random-POVM probe first optimises a leaf against a frozen
extremal rank-one POVM, replaces that readout by the exact SOCP optimum, and
then releases all parameters.  After correcting the physical/transposed POVM
convention at reinjection, twelve deterministic seeds leave nine trajectories
four-active, all converging to the no-record score \(0.7\).  The three
trajectories that exceeded that basin became two-active
and converged to 3E or 4E.  The largest score was
\(0.765898815264695\ldots\), and no three-active basin appeared.  This is a
targeted falsification result, not a global upper bound.

The terminal POVM can be eliminated exactly for any reconstructed leaf.  If

\[
\tau_s=\tfrac12(p_sI+r_s\cdot\sigma),
\]

then qubit minimum-error discrimination is the additively weighted
smallest-ball problem

\[
P_{\rm guess}
=\min_{y\in\mathbb R^3}\max_s
\left(p_s+\lVert y-r_s\rVert_2\right).
\tag{21}
\]

Two active ball constraints are equivalent to a binary projective optimum;
three or four active constraints are exactly the genuine nonprojective
sectors.  For the unrestricted complex 4E checkpoint, the independent SOCP
has only labels zero and three active.  The other two constraint slacks are
\(0.7821\ldots\) and \(0.7814\ldots\), so the candidate lies strictly inside
the projective region rather than on an arity-degenerate boundary.  The proof
and diagnostics are in `notes/terminal_qubit_discrimination_geometry.md`.

These results originally narrowed the unrestricted converse to one tempting
claim:

> A support-maximising two-block rank-two leaf admits an optimal terminal
> AUDIT POVM with two projective outcomes.

That claim is not proved, and it is not assumed in the final enclosure.
General extremal qubit POVMs can have three or four rank-one effects, so
ordinary POVM extremality does not establish it.  Instead, the three- and
four-active weighted-ball sectors have now been bounded directly.  Combining
them with equation (20) gives

\[
0.7658988152646944
\leq\beta_{2\mathrm b}(0.6)
\leq0.76662.
\]

The stronger upper endpoint is solver-conditional and does not prove exact
terminal projectivity or equality with the physical four-effect leaf.  See
`notes/interleaved_interior_frontier_l060.md` for the exhaustive arity split.

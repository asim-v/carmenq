# A Causal Order Gap Beyond the Triangle Instance

**Status:** rigorous grouped attainer, rigorous canonical-bond obstruction,
and an exact arbitrary-instrument theorem at the perfect-AUDIT endpoint for
the interleaved order. A causal-list extension now proves a linear rank-tail
bound and a strict support gap for \(3/7<\lambda<1\). The complete interleaved
support function remains open.

Consider the two rank-two, four-slot check matrices

\[
H_{\rm G}=
\begin{pmatrix}
1&1&0&0\\
0&0&1&1
\end{pmatrix},
\qquad
H_{\rm I}=
\begin{pmatrix}
1&0&1&0\\
0&1&0&1
\end{pmatrix}.
\tag{1}
\]

They have the same rank and the same abstract row code up to a coordinate
permutation. They differ only in temporal order. The persistent coherent
dimension is \(d=2\).

The universal terminal-dimension argument gives both problems the same static
ceiling

\[
B_{4,2}(\lambda)
=\frac{1+\sqrt{\lambda^2+(1-\lambda)^2}}2.
\tag{2}
\]

The grouped order attains (2). At perfect AUDIT its RETURN value is therefore
\(1/2\). In contrast, the exact perfect-AUDIT RETURN value of the interleaved
order is only \(1/4\). This is the first instance in the present search at
which the ordered columns, rather than only \(K=4\) and \(d=2\), provably
enter.

## 1. Grouped order saturates the static ceiling

Write

\[
S_1=X_1\oplus X_2,\qquad S_2=X_3\oplus X_4.
\]

One qubit \(M\) initially in \(|0\rangle\) implements the complete static
boundary. Accumulate \(S_1\) with CNOTs in slots 1 and 2, apply the binary weak
\(Z_M\) instrument of strength \(t\) after the second CNOT, and keep its
outcome \(c\) classically. Continue applying CNOTs from slots 3 and 4. At the
terminal cut,

\[
M=S_1\oplus S_2.
\]

AUDIT uses \(c\) to estimate \(S_1\) and \(M\) to recover \(S_2\). RETURN
applies the four CNOTs in reverse order and resets \(M\) exactly. Therefore

\[
P_{\rm A}=\frac{1+t}{2},
\qquad
F_{\rm R}=\frac{1+\sqrt{1-t^2}}2,
\]

which attains (2) at
\(t=\lambda/\sqrt{\lambda^2+(1-\lambda)^2}\). Immediate sequestration is
respected; no coherent system other than \(M\) crosses a cut.

## 2. Canonical equality tensors and their live-form width

Every two-versus-two partition of the four syndromes is an affine level set
of a nonzero row-code character. A canonical equality branch weakly records
one character \(a\in\operatorname{row}(H)\setminus\{0\}\) and retains an
independent character \(b\) in the terminal qubit:

\[
K_{a,b,c}|x\rangle
=\ell_c(a\cdot x)\,|x\rangle_B\,|b\cdot x\rangle_M,
\tag{3}
\]

where both values of \(\ell_c\) are nonzero and unequal at an interior point.

There is a closed Schmidt-rank calculation for (3). At a temporal cut
\(x=(x_L,x_R)\), write

\[
\ell_c(z)=u_c+v_c(-1)^z,\qquad u_cv_c\neq0,
\]

and expand the terminal-bit Kronecker delta in its two Fourier characters.
The coefficient matrix across the cut becomes a sum of four rank-one terms
indexed by \((p,q)\in\mathbb F_2^2\). Their left characters are

\[
p\,a_L+q\,b_L,
\]

and their right characters, including the terminal memory bit, are

\[
(p\,a_R+q\,b_R,\ q_M).
\]

If \(a_R\neq0\), the four right characters are distinct. Fourier characters
are linearly independent, so

\[
\operatorname{SchmidtRank}_{L|R}(K_{a,b,c})
=2^{\operatorname{rank}_{\mathbb F_2}\{a_L,b_L\}}.
\tag{4}
\]

This suggests the ordered live-form width

\[
w(a,b)
=\max_{i:\,a_{>i}\neq0}
\operatorname{rank}_{\mathbb F_2}
\{a_{\leq i},b_{\leq i}\}.
\tag{5}
\]

The canonical branch has coherent bond at least \(2^{w(a,b)}\). Conversely,
when \(w(a,b)=1\), the usual accumulate--weakly-read--handoff construction
gives a bond-two realization: before \(a\) is complete, the two live prefix
forms are the same or one is zero; after \(a\) is recorded, the qubit may be
continued into \(b\).

For \(H_{\rm G}\), choose

\[
a=1100,\qquad b=1111.
\]

Then \(w(a,b)=1\), reproducing the construction above. For the triangle
instance, \(a=110\) and \(b=101\) likewise give \(w=1\).

For \(H_{\rm I}\), the middle cut has

\[
H_{{\rm I},L}=H_{{\rm I},R}=I_2.
\]

Hence the restrictions of every independent pair of row-code characters
remain independent on the first two slots, while every nonzero \(a\) remains
live on the last two slots. Thus

\[
w(a,b)=2
\quad\text{for all independent }a,b\in\operatorname{row}(H_{\rm I}),
\tag{6}
\]

and every canonical equality branch has Schmidt rank four at the middle cut.
An exhaustive direct calculation over the six ordered independent pairs gives
the same result.

Equation (6) is a rigorous obstruction for the complete canonical
linear/Lüders equality family. By itself it is not yet a converse for every
non-QND instrument: a general polar isometry could encode the two favoured
syndromes differently. The next section gives a converse at one important
endpoint without assuming the canonical form.

## 3. Exact interleaved perfect-AUDIT endpoint

> **Theorem (crossed-pair endpoint).** For every streamed strategy for
> \(H_{\rm I}\) with a two-dimensional coherent system at every cut,
> genuinely classical transcript, sequestered earlier outputs, and an
> arbitrary transcript-conditioned RETURN decoder,
> \[
> P_{\rm A}=1\quad\Longrightarrow\quad F_{\rm R}\leq\frac14.
> \tag{7}
> \]
> The bound is attained. Hence
> \[
> \max\{F_{\rm R}:P_{\rm A}=1\}=\frac14.
> \tag{8}
> \]

The theorem covers adaptive non-QND instruments. It neither assumes a
Clifford realization nor diagonal or Lüders branch operators.

### 3.1 Fine branches and the four-word support lemma

Refine every physical outcome to its Kraus label and reveal the refinement to
both decoders. This can only increase both optimized scores, preserves the
streamed realization, and turns every complete transcript \(c\) into a single
sequential Kraus branch. Define its basis likelihood

\[
p(c\mid x)=\lVert K_c|x\rangle\rVert^2,
\qquad x\in\mathbb F_2^4.
\tag{9}
\]

If the unrefined strategy has \(P_{\rm A}=1\), so does the refinement. Zero
average AUDIT error then implies that, for each \(c\), the terminal-memory
states associated with different supported syndromes are perfectly
distinguishable. A qubit can support at most two such syndromes. A one-syndrome
branch already has at most four supported words, because every syndrome fibre
has size four. It remains to consider a branch whose supported syndrome set is
a pair

\[
T=\{s:u\mathbin{\cdot}s=b\},
\qquad 0\neq u\in\mathbb F_2^2.
\tag{10}
\]

Split the branch after slot 2. Write the basis word as \(x=(z,y)\), where
\(z=(x_1,x_2)\) is the prefix and \(y=(x_3,x_4)\) the suffix. Once the prefix
outcomes in \(c\) are fixed, the prefix produces a possibly zero pure vector
\(|\psi_z\rangle\in B_1B_2\otimes M_2\). For a fixed basis suffix \(y\), the
rest of the same complete branch is a linear map

\[
I_{B_1B_2}\otimes L_y,
\qquad L_y:M_2\longrightarrow B_3B_4\otimes M_4.
\tag{11}
\]

Immediate sequestration is exactly what makes \(B_1B_2\) a spectator in
(11). Since both halves of the interleaved check restrict to the identity,
the words compatible with (10) at fixed \(y\) have prefix in the affine line

\[
L(y)=\{z:u\mathbin{\cdot}z=b\oplus u\mathbin{\cdot}y\},
\tag{12}
\]

which contains two prefixes.

Suppose first that, for some \(y\), both prefixes in \(L(y)\) have positive
branch likelihood. They lead to the two distinct syndromes in \(T\), so their
terminal qubit supports are orthogonal. Therefore the continuation map
\(L_y\) has rank two and is injective. Indeed, a rank-one map has the form
\(|\omega\rangle\langle\phi|\); even for a prefix entangled with
\(B_1B_2\), every nonzero output then has the same terminal-memory support and
cannot encode two orthogonal labels.

The two prefixes outside \(L(y)\) lead to syndromes outside \(T\), hence their
full branch vectors under \(L_y\) vanish. Injectivity of \(L_y\) implies that
their prefix vectors \(|\psi_z\rangle\) were already zero. Thus this complete
branch can occur only for the two prefix values in \(L(y)\). Among suffixes,
only the two values with the same \(u\)-parity as \(y\) keep those prefixes in
the supported syndrome pair. The branch consequently supports at most
\(2\times2=4\) basis words.

If no suffix has both compatible prefixes with positive likelihood, then at
most one word is supported for each of the four suffix values, again giving
at most four. We have proved the branchwise statement

\[
\bigl|\operatorname{supp}p(c\mid\cdot)\bigr|\leq4
\quad\text{for every refined transcript }c.
\tag{13}
\]

The injectivity step remains valid for arbitrary prefix entanglement: writing
a Schmidt decomposition of \(|\psi_z\rangle\) across \(B_1B_2|M_2\), an
injective map on \(M_2\) cannot annihilate a nonzero Schmidt component. The
complete transcript fixes every adaptive path, so no QND assumption entered.

### 3.2 RETURN bound

For \(D=16\), the flagged polar-recovery bound and computational-basis
pinching give, for each refined leaf,

\[
F_{{\rm R},c}
\leq\frac1{16^2}
\left(\sum_x\sqrt{p(c\mid x)}\right)^2.
\tag{14}
\]

Cauchy--Schwarz and (13) imply

\[
F_{{\rm R},c}
\leq\frac4{256}\sum_xp(c\mid x).
\]

Summing over \(c\) and using trace preservation,

\[
F_{\rm R}
\leq\frac4{256}\sum_x\sum_cp(c\mid x)
=\frac4{256}\,16
=\frac14.
\tag{15}
\]

The polar-recovery lemma already optimizes over arbitrary joint decoders on
all returned carriers and the terminal memory; a demanded memory reset can
only lower the value. Thus (15) is an unconditional RETURN upper bound.

### 3.3 Attainer

Initialize \(M=0\). Apply \(\operatorname{CNOT}_{A_1\to M}\), projectively
measure \(A_2\) in the computational basis and record \(c_2\), apply
\(\operatorname{CNOT}_{A_3\to M}\), and projectively measure \(A_4\), recording
\(c_4\). At the terminal cut,

\[
M=X_1\oplus X_3=S_1,
\qquad c_2\oplus c_4=S_2,
\]

so AUDIT is perfect. In RETURN, apply
\(\operatorname{CNOT}_{B_3\to M}\) and then
\(\operatorname{CNOT}_{B_1\to M}\), which exactly resets \(M\) and restores
pairs 1 and 3. Each of the two flagged projective measurements has optimal
EPR recovery fidelity \(1/2\); their product gives \(F_{\rm R}=1/4\).

### 3.4 A strict support gap near the endpoint

The closure issue caused by an unbounded classical alphabet can be handled
directly at the level of refined leaves. Put \(D=16\). For every nonzero leaf
define

\[
w_c=D^{-1}\operatorname{Tr}(K_c^\dagger K_c),
\qquad \widehat K_c=K_c/\sqrt{w_c}.
\]

Then \(\sum_cw_c=1\) and
\(\operatorname{Tr}(\widehat K_c^\dagger\widehat K_c)=D\). Regard the strategy
as the probability measure

\[
\mu=\sum_cw_c\,
\delta_{(\widehat K_c,E_c,\mathcal D_c)},
\]

where \(E_c\) is its terminal qubit POVM and \(\mathcal D_c\) its RETURN
decoder. The normalized streamed leaf tensors form a compact set: their norm
is fixed, and every temporal matricization has rank at most two, a closed
determinantal condition. The sets of qubit POVMs and finite-dimensional CPTP
decoders are compact as well. Hence any sequence of strategies has a
weakly-convergent subsequence of these probability measures. The constraint

\[
\int \widehat K^\dagger\widehat K\,d\mu=I
\]

and both scores are continuous under this convergence.

For a normalized leaf, let \(a\in[0,1]\) and \(f\in[0,1]\) denote its AUDIT
and RETURN contributions per unit weight. If a limiting strategy has
\(P_{\rm A}=\int a\,d\mu=1\), then \(a=1\) almost everywhere. The four-word
support lemma therefore applies almost everywhere, and (14) gives
\(f\leq1/4\) almost everywhere. Thus every operational limit still obeys
\(F_{\rm R}=\int f\,d\mu\leq1/4\).

Let \(\widetilde{\mathcal R}_{\rm I}\subset[0,1]^2\) be the score image of
all probability measures on this compact leaf space that satisfy the displayed
closed completeness constraint. It is a compact relaxation containing the
closure of the genuinely streamed achievable region, even though a generic
measure in the relaxation need not retain all shared-prefix compatibility
conditions of an adaptive instrument. The preceding leafwise argument proves
that the relaxation itself still has
\(F_{\rm R}\leq1/4\) on its perfect-AUDIT slice. This avoids any assumption
that the unbounded classical transcript has a preassigned finite alphabet.

The static dimension relaxation has, for every \(0<\lambda<1\), the unique
exposed boundary point

\[
z_\lambda=
\left(
\frac{1+t_\lambda}{2},
\frac{1+\sqrt{1-t_\lambda^2}}2
\right),
\qquad
t_\lambda=
\frac{\lambda}{\sqrt{\lambda^2+(1-\lambda)^2}},
\tag{16}
\]

and \(z_\lambda\to(1,1/2)\) as \(\lambda\uparrow1\). The endpoint bound says
that the compact relaxation \(\widetilde{\mathcal R}_{\rm I}\) does not
contain \((1,1/2)\); hence a whole neighbourhood of that point is disjoint
from \(\widetilde{\mathcal R}_{\rm I}\). Therefore there exists
\(\lambda_0<1\) such that
\(z_\lambda\notin\widetilde{\mathcal R}_{\rm I}\) whenever
\(\lambda_0<\lambda<1\). Since \(z_\lambda\) is the unique point attaining
the universal static support value, and the homogeneous static proof applies
unchanged to the leaf-measure relaxation, compactness gives the strict
inequality

\[
\beta^{\rm stream}_{H_{\rm I},2}(\lambda)
<B_{4,2}(\lambda),
\qquad \lambda_0<\lambda<1.
\tag{17}
\]

Thus the endpoint theorem already proves a nonempty interval of genuine
causal order dependence. The later proof in `order_gap_linear_tail.md`
certifies the much larger interval \(3/7<\lambda<1\), although the exact full
interior frontier remains unknown.

### 3.5 General full-crossing cut criterion

The support argument is not peculiar to four slots. Let \(H\) have two
independent rows and split an arbitrary \(n\)-slot order at a cut \(i\):

\[
H=(H_L\mid H_R),
\qquad
\operatorname{rank}H_L=\operatorname{rank}H_R=2.
\tag{18}
\]

For a perfect-AUDIT refined leaf supported on a syndrome pair
\(T=\{s:u\cdot s=b\}\), fix a suffix \(y\). If the leaf is nonzero on at least
one prefix for each of the two syndromes in \(T\), the continuation must again
have rank two and hence be injective. All prefix vectors outside the affine
half-space

\[
u\cdot H_Lz=b\oplus u\cdot H_Ry
\]

were therefore already zero. This leaves at most \(2^{i-1}\) prefixes, and
each can be paired with at most \(2^{n-i-1}\) suffixes because
\(u\cdot H_R\neq0\). If no suffix sees both supported syndrome labels, then
each suffix supports prefixes from at most one fibre of the surjective map
\(H_L\), again at most \(2^{i-2}\) prefixes per suffix. In either case, and
also for a one-syndrome leaf,

\[
|\operatorname{supp}p(c\mid\cdot)|
\leq2^{n-2}.
\]

The flagged recovery argument consequently gives the general implication

\[
P_{\rm A}=1
\quad\Longrightarrow\quad
F_{\rm R}\leq\frac{2^{n-2}}{2^n}=\frac14
\tag{19}
\]

for every two-check stream with a full-crossing cut and a persistent qubit.
This is an upper theorem; attainability of \(1/4\) can still depend on the
ordered columns. The interleaved four-slot matrix attains it by Section 3.3.

A convenient order invariant is the maximum syndrome-trellis cut dimension

\[
\tau(H)=\max_i\left(
\operatorname{rank}H_{\leq i}
+\operatorname{rank}H_{>i}
-\operatorname{rank}H
\right).
\tag{20}
\]

For two checks, condition (18) is equivalent to \(\tau(H)=2\). The
interleaved order has \(\tau=2\), whereas the grouped order and the triangle
instance both have \(\tau=1\). Thus \(\tau\) cleanly separates this exact
order gap and supplies the promised connection to ordered linear-code or
trellis width. Whether \(\tau\), together with finer branch data,
characterizes the complete interior frontier remains open.

## 4. What numerics currently say

A simple legal interleaved lower strategy stores one row coherently and weakly
measures the two separated carrier bits belonging to the other row. With
equal local strength,

\[
P_{\rm A}=\frac{1+t^2}{2},\qquad
F_{\rm R}=\left(\frac{1+\sqrt{1-t^2}}2\right)^2.
\tag{21}
\]

This family is not optimal. An exploratory unrestricted binary-tree search,
which optimises arbitrary complex non-QND node instruments and
transcript-conditioned qubit POVMs, found at \(\lambda=1/2\)

\[
P_{\rm A}=0.6446433867,\qquad
F_{\rm R}=0.8662315060,\qquad
\mathcal S=0.7554374463.
\tag{22}
\]

This exceeds the best value \(0.75\) of the simple family at balanced weight,
but remains far below the static ceiling

\[
B_{4,2}(1/2)=0.8535533906.
\]

The computation is a lower-bound search, not a certificate. It is valuable
mainly because it falsifies the tempting conjecture that independent weak
records plus one exact row give the complete interleaved frontier.

Independent high-AUDIT runs are consistent with the exact endpoint. At
\(\lambda=0.9\) they give approximately
\((P_{\rm A},F_{\rm R})=(0.998226,0.280860)\), and at \(\lambda=0.99\) they
give \((0.999987,0.252549)\). These are lower bounds, not ingredients of the
proof.

## 5. Remaining theorem

The remaining exact target is the full support function

\[
\beta^{\rm stream}_{H_{\rm I},2}(\lambda)
\quad(0<\lambda<1),
\tag{23}
\]

together with attaining strategies and analytic or computer-assisted upper
certificates. Equation (17) already certifies strict separation from the
static ceiling on a nonempty interval. The endpoint theorem fixes one boundary
value exactly, while the live-form calculation identifies the obstruction in
the canonical family. The linear-tail companion gives a quantitative gap for
\(3/7<\lambda<1\). None of these results supplies the exact interior support.

A complete proof will probably need one of the following:

1. a rigidity theorem promoting equality in the terminal dimension/recovery
   bounds to a branchwise live-form lower bound;
2. a two-cut constrained-separability dual certificate specialised to the
   middle identity-rank cut; or
3. a dynamic-programming recursion over qubit instruments that preserves the
   full classical transcript without silently convexifying the fixed
   dimension.

This order comparison is a stronger originality target than the triangle
theorem. The grouped and interleaved matrices have identical rank and differ
only by temporal permutation, so a strict certified gap would make causal
order and intermediate memory cuts operationally indispensable rather than
decorative.

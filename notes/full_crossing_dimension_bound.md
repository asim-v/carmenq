# A temporal product law and exponential order separation

**Date:** 20 August 2026<br>
**Status:** proved endpoint theorem with a saturating family; priority language
remains conservative<br>
**Scope:** finite-field linear syndromes, perfect AUDIT, bounded coherent
memory, free genuinely classical transcript, and unconditional EPR RETURN

## 1. Result

Let \(H\in\mathbb F_q^{r\times n}\) have rank \(r\). Partition its ordered
columns into \(m\) consecutive nonempty blocks,

\[
H=[H^{(1)}\mid H^{(2)}\mid\cdots\mid H^{(m)}],
\qquad
\operatorname{rank}H^{(j)}=r
\quad(1\le j\le m).
\tag{1}
\]

Write \(N=q^r\) for the number of syndrome labels. A streamed device may
retain arbitrary finite classical information. Let \(d_j\) bound the coherent
dimension crossing the boundary after block \(j\), including the terminal
AUDIT memory at \(j=m\). Earlier output carriers are sequestered. If the device
identifies the full syndrome with certainty, then its optimized unconditional
EPR-return fidelity satisfies the product law

\[
\boxed{
F_{\rm R}
\le
\prod_{j=1}^m\min\left\{1,\frac{d_j}{N}\right\}.}
\tag{2}
\]

In particular, if every boundary has coherent dimension at most \(d\),

\[
F_{\rm R}\leq\min\left\{1,\left(\frac dN\right)^m\right\}.
\tag{3}
\]

For an arbitrary ordered full-rank check matrix, define its full-rank block
packing number

\[
\mu(H)=\max\left\{m:
H=[H^{(1)}\mid\cdots\mid H^{(m)}],\quad
\operatorname{rank}H^{(j)}=r\ \text{for every }j
\right\}.
\tag{3a}
\]

Choosing a maximizing partition in equation (3) gives the matrix-level
corollary

\[
\boxed{
P_{\rm A}=1
\quad\Longrightarrow\quad
F_{\rm R}\leq
\min\left\{1,\left(\frac d{q^r}\right)^{\mu(H)}\right\}.}
\tag{3b}
\]

The value \(\mu(H)\) is computed greedily: close each block at the earliest
column where it reaches rank \(r\), and append any final rank-deficient tail
to the last block. If \(t_j\) is the \(j\)-th greedy endpoint and \(s_j\) is
the endpoint in any competing partition, induction gives \(t_j\leq s_j\).
Thus no other consecutive partition contains more full-rank blocks.

Equation (3) is tight on the canonical repeated-identity family

\[
H=[I_r\mid I_r\mid\cdots\mid I_r]
\tag{4}
\]

whenever \(d=q^k\). Coherently retaining \(k\) syndrome coordinates and
projectively recording the remaining coordinates in every block gives

\[
F_{\rm R}
=q^{-m(r-k)}
=\left(\frac{d}{N}\right)^m.
\tag{5}
\]

The cases have direct meanings:

\[
\begin{array}{c|c}
m=1 & F_{\rm R}\le d/N\quad\text{(static dimension limit)},\\
m=2 & F_{\rm R}\le(d/N)^2\quad\text{(full-crossing square law)},\\
m>2 & \text{each full-rank temporal block adds one factor }d/N.
\end{array}
\tag{6}
\]

### Exponential separation under a column permutation

The product law yields an exact family in which temporal order alone creates
an exponential gap. Let \(e_1,\ldots,e_r\) be the standard columns and compare

\[
H_{\rm batched}
=
[\underbrace{e_1\mid\cdots\mid e_1}_{m}\mid\cdots\mid
  \underbrace{e_r\mid\cdots\mid e_r}_{m}],
\qquad
H_{\rm cycled}
=
[I_r\mid\cdots\mid I_r].
\tag{6a}
\]

The matrices contain the same multiset of columns and differ only by a
permutation. They satisfy

\[
\mu(H_{\rm batched})=1,
\qquad
\mu(H_{\rm cycled})=m.
\tag{6b}
\]

For the batched order, full rank first appears only when the stream reaches
the \(e_r\) batch, after which the remaining columns have rank one. For the
cycled order, every identity block has full rank, and no rank-\(r\) block can
use fewer than \(r\) columns. These observations prove equation (6b).

Let \(d=q^k\) with \(1\leq k<r\). At perfect AUDIT, their exact optimized
RETURN fidelities are

\[
F_{\rm R}^{\star}(H_{\rm batched})
=q^{-(r-k)}=\frac dN,
\qquad
F_{\rm R}^{\star}(H_{\rm cycled})
=q^{-m(r-k)}=\left(\frac dN\right)^m.
\tag{6c}
\]

For the batched order, one memory qudit sequentially accumulates and is
projectively read for each of the first \(r-k\) coordinate syndromes; it is
reset after each read. The final \(k\) coordinate syndromes are retained
coherently in the \(k\) memory qudits. Every measured syndrome contributes
one factor \(1/q\), so the matrix-level upper bound is attained. The cycled
value is attained by the repeated-identity construction of Section 4.
Consequently,

\[
\frac{F_{\rm R}^{\star}(H_{\rm batched})}
     {F_{\rm R}^{\star}(H_{\rm cycled})}
=
\left(\frac Nd\right)^{m-1}
=q^{(r-k)(m-1)}.
\tag{6d}
\]

Thus a coordinate permutation of one linear check can change recoverable EPR
fidelity by an exponential factor under the same coherent-memory limit.

For the binary rank-two interleaved benchmark, \(\mu(H)=2\), \(N=4\), and
\(d=2\), so equation (3) recovers the exact endpoint

\[
P_{\rm A}=1\quad\Longrightarrow\quad F_{\rm R}\le\frac14.
\tag{7}
\]

Equivalently, a perfect-AUDIT experiment supplies coherent-dimension
witnesses

\[
\left(\prod_{j=1}^m d_j\right)^{1/m}
\ge N F_{\rm R}^{1/m},
\qquad
d\ge N F_{\rm R}^{1/m}\ \text{in the uniform case}.
\tag{8}
\]

The grouped permutation has \(\mu(H)=1\), while the interleaved permutation
has \(\mu(H)=2\). Repeated identity blocks have arbitrary \(\mu(H)=m\). The
contribution is not a new definition of quantum memory, a new recovery
formula, or a claim that the classical sectionalization idea is new. It is the
multiplicative recovery consequence of repeatedly forcing a full linear
syndrome through a bounded coherent bond.

## 2. Frozen interface

Let \(q\) be a prime power. The verifier streams halves of \(n\) maximally
entangled \(q\)-level pairs.
At each slot the device acts on the fresh carrier, a charged coherent memory,
fresh disposable ancillas, and its classical transcript. The emitted carrier
is immediately sequestered and cannot be used at a later slot. The classical
transcript is genuinely dephased; a coherent purification is not free.
Disposable ancillas begin in fixed product states: no uncharged quantum side
channel or entangled common cause crosses slots.

After the stream, a hidden choice selects one of two tasks. In **AUDIT**, the
references are measured in the computational basis, producing uniform
\(X\in\mathbb F_q^n\), and the device must output \(HX\) using only its
transcript and terminal coherent memory. In **RETURN**, a
transcript-conditioned decoder receives all emitted carriers and the terminal
memory, and attempts to restore all EPR pairs. Failure probability and memory
reset are included; there is no postselection.

The bounded-memory comb model itself is established
[@bisio2012memory; @ohst2026memory], as are flagged recovery and general
information--disturbance formulations
[@gregoratti2003lostfound; @hsieh2026interactive]. Equation (2) concerns the
additional restriction imposed by the ordered full-rank block sequence.

## 3. Proof

Put \(D=q^n\). Reveal every local Kraus refinement as extra classical
transcript. This can only help both late decoders and does not enlarge the
coherent memory. A complete refined transcript \(c\) is therefore represented
by one sequential Kraus leaf \(K_c\).

For a basis word \(x\), define

\[
p_c(x)=\lVert K_c|x\rangle\rVert^2.
\tag{9}
\]

### 3.1 At most \(d_j\) cumulative labels cross boundary \(j\)

At the terminal boundary, let \(U_{c,m}\subseteq\mathbb F_q^r\) be the set of
complete syndrome labels represented by nonzero words in the leaf. Perfect
AUDIT makes their terminal-memory supports mutually orthogonal. Therefore

\[
|U_{c,m}|\le \min\{d_m,N\}.
\tag{10}
\]

Now fix an earlier block boundary \(j<m\). Write a word as \(x=(z,y)\), where
\(z\) contains the first \(j\) blocks. Immediate sequestration gives

\[
K_c|z,y\rangle
=
(I_{B_{\le j}}\otimes L_{c,y})|\psi_{c,z}\rangle,
\tag{11}
\]

where \(L_{c,y}\) acts only on the memory crossing that boundary.

Let \(U_{c,j}\) contain the cumulative prefix labels

\[
u=H^{(1)}z_1+\cdots+H^{(j)}z_j
\tag{12}
\]

that occur in some nonzero complete word. For \(u\in U_{c,j}\), let
\(R_u\) be the span of every cut-memory Schmidt vector appearing in prefix
states with cumulative label \(u\).

Choose one supported suffix \(y_u\), and let \(s_u\) be the resulting complete
syndrome. If \(Q_{s_u}\) projects onto its terminal AUDIT support, then

\[
A_u=(I\otimes Q_{s_u})L_{c,y_u}
\tag{13}
\]

is nonzero on \(R_u\). It annihilates every \(R_{u'}\) with \(u'\ne u\):
under the same suffix, that prefix either gives zero or ends in the orthogonal
terminal support belonging to the distinct label \(s_u+(u'-u)\).

Choose \(r_u\in R_u\) with \(A_ur_u\ne0\). Applying \(A_u\) to a linear
relation among the \(r_{u'}\) isolates its \(u\)-th coefficient. Hence these
vectors are linearly independent and

\[
|U_{c,j}|\le\min\{d_j,N\}
\qquad(1\le j\le m).
\tag{14}
\]

This is the causal step: each cumulative prefix-syndrome class needs a private
direction in the coherent memory.

### 3.2 Counting block-syndrome tuples

Let \(v_j=H^{(j)}x_j\) be the syndrome contribution of block \(j\), and define
the cumulative values

\[
u_j=v_1+\cdots+v_j.
\tag{15}
\]

The map

\[
(v_1,\ldots,v_m)
\longleftrightarrow
(u_1,\ldots,u_m)
\tag{16}
\]

is bijective, since \(v_1=u_1\) and \(v_j=u_j-u_{j-1}\). Every supported
tuple has \(u_j\in U_{c,j}\); equation (14) therefore leaves at most
\(\prod_j\min\{d_j,N\}\) supported block-syndrome tuples.

If block \(j\) contains \(n_j\) input symbols, its full-rank map has uniform
fiber size \(q^{n_j-r}\). Thus each supported tuple represents exactly
\(q^{n-mr}\) computational words, and

\[
|\operatorname{supp}p_c|
\le q^{n-mr}\prod_{j=1}^m\min\{d_j,N\}
=D\prod_{j=1}^m\min\left\{1,\frac{d_j}{N}\right\}.
\tag{17}
\]

### 3.3 From support to EPR return

For a refined leaf, optimal flagged recovery is bounded by

\[
F_{{\rm R},c}
\le
\frac{1}{D^2}
\left(\operatorname{Tr}\sqrt{K_c^\dagger K_c}\right)^2.
\tag{18}
\]

Schur--Horn majorization and concavity of the square root imply

\[
\operatorname{Tr}\sqrt{K_c^\dagger K_c}
\le\sum_x\sqrt{p_c(x)}.
\tag{19}
\]

If the support contains at most \(\alpha D\) words, Cauchy--Schwarz gives

\[
F_{{\rm R},c}
\le
\frac{\alpha}{D}\sum_xp_c(x).
\tag{20}
\]

The reset requirement may be dropped for an upper bound; retaining it can
only decrease RETURN. Instrument completeness gives
\(\sum_c p_c(x)=1\) for every word. Summing equation (20) over all transcripts
yields \(F_{\rm R}\le\alpha\), which is equation (2).

## 4. Tight repeated-identity construction

Take equation (4), let \(d=q^k\), and use a memory of \(k\) qudits. In every
block, coherently add the first \(k\) input symbols into memory and
projectively record the remaining \(r-k\) input symbols.

AUDIT reads the \(k\) coherently accumulated syndrome coordinates and computes
the others from the classical records. RETURN reverses every coherent update.
Each of the \(m(r-k)\) projectively recorded carriers contributes flagged EPR
fidelity \(1/q\), giving equation (5).

Thus neither the exponent \(m\) nor its coefficient can be improved under the
declared assumptions. The construction proves family-level tightness, not
that every full-rank block sequence attains equation (3).

## 5. Interpretation and limits

The theorem turns temporal fragmentation into a quantitative resource law.
For fixed syndrome size \(N\) and coherent dimension \(d<N\), perfect classical
readout becomes exponentially more destructive in \(\mu(H)\). Conversely,
observed perfect-AUDIT return fidelity lower-bounds the coherent bond dimension
through equation (8).

No claim is made about the unresolved interior \(P_{\rm A}<1\), approximate
AUDIT, matrices that cannot be partitioned into full-rank blocks, nonclassical
uncharged side channels, or access to sequestered carriers during AUDIT.

## 6. Priority-safe claim

Classical trellis connectivity, coordinate-order dependence, and trellis
sectionalization are standard [@forney1994trellis;
@lafourcade1996sectionalization]. Nondestructive discrimination and bounded
quantum-memory hierarchies are also established
[@lim2025local; @ohst2026memory]. A focused search found no primary source
stating equations (2)--(3) for the late syndrome-AUDIT versus all-carrier
EPR-RETURN game. Absence from that search is not proof of priority. The safe
claim is:

> For the declared streamed interface, \(m\) consecutive full-rank linear
> syndrome blocks impose the perfect-AUDIT product law
> \(F_{\rm R}\le\prod_j\min\{1,d_j/q^r\}\); the uniform-dimension power law is
> tight on the repeated-identity family for \(d=q^k\). Equivalently, a uniform
> bound \(d\) gives \(F_{\rm R}\le(d/q^r)^{\mu(H)}\), capped at one.

The theorem should receive independent external checking before stronger
priority language is used.

# Common-instrument Fourier contraction

## Status

This note records a rigorous necessary condition for the fixed interior
frontier and the numerical attempt to turn it into a complete certificate. It
does **not** close the target value \(0.758\). Its main value is that it
replaces the vague requirement that all conditioned states “come from one
instrument” by explicit trace-norm inequalities and identifies exactly why
the three canonical inequalities are not sufficient.

## Harmonic contraction theorem

Let \(G\) be a finite abelian group and let \(\{\Phi_y\}_{y\in G}\) be a
quantum instrument, so each \(\Phi_y\) is completely positive and the flagged
channel

\[
  \Gamma(X)=\bigoplus_{y\in G}\Phi_y(X)
\]

is trace preserving. For a Hermitian family \(\{A_z\}_{z\in G}\), define the
convolutional readout

\[
  T_s=\sum_{z\in G}\Phi_{z-s}(A_z).
\]

For every character \(\chi\in\widehat G\), write
\(\widehat A_\chi=\sum_z\chi(z)A_z\). Character multiplicativity gives

\[
  \widehat T_\chi
  =\sum_y\chi(-y)\Phi_y(\widehat A_\chi).
\]

The triangle inequality and trace-norm contractivity of the flagged CPTP
channel therefore imply

\[
  \|\widehat T_\chi\|_1
  \leq \sum_y\|\Phi_y(\widehat A_\chi)\|_1
  =\|\Gamma(\widehat A_\chi)\|_1
  \leq\|\widehat A_\chi\|_1.
\]

The middle, flagged inequality is stronger than contraction of the terminal
Fourier component because it forbids cancellations between instrument
outcomes. The proof uses a single common instrument. Applying a different
channel to every input would not justify it.

For a Hermitian qubit operator
\(X=(dI+\mathbf r\cdot\boldsymbol\sigma)/2\),

\[
  \|X\|_1=\max\{|d|,\|\mathbf r\|_2\}.
\]

Thus each character has three regimes: positive scalar, negative scalar, or
Bloch-vector dominated. Scalar regimes are second-order-cone representable.
The vector regime is reverse convex, but it has a finite conic outer cover. If
\(\mathbf r\) lies in the spherical cap
\(\mathbf n\cdot\mathbf r\geq c\|\mathbf r\|_2\), then

\[
  \sum_y\|\Phi_y(X)\|_1
  \leq \|\mathbf r\|_2
  \leq \frac{\mathbf n\cdot\mathbf r}{c}.
\]

A collection of caps covering the sphere therefore converts the exact
disjunction into finitely many convex outer problems. For
\(G=\mathbb Z_2^2\), Parseval supplies the additional joint inequality

\[
 \sum_{\chi\neq 1}\|\mathbf r_\chi\|_2^2
 =4\sum_z\|\mathbf r_z\|_2^2
  -\left\|\sum_z\mathbf r_z\right\|_2^2.
\]

On a prior cell \(a_z\in[\ell_z,u_z]\), positivity and the secant bound on
\(a_z^2\) turn the right-hand side into an affine upper bound. These are
necessary conditions at every physical point, not heuristic penalties.

## Gauge reduction and cap covers

Simultaneously rotating all four input states and precomposing the common
instrument by the inverse unitary leaves every conditioned output unchanged.
The first vector-active Fourier component can therefore be aligned with
\(+z\). The residual rotation about \(z\) puts the second active component in
the \(xz\) plane with nonnegative \(x\). Only the second planar angle and, if
present, the third spherical direction need finite covers.

The implemented spherical cover uses normalized cube-face charts. On a face,
write a direction as
\(f(u,v)=(1,u,v)/\|(1,u,v)\|\). The derivative norm of \(f\) is at most one.
A square of half-width \(1/N\) is therefore contained in a cap with chord
radius at most \(\sqrt2/N\), hence covering cosine at least \(1-1/N^2\). This
supplies a proved cover rather than a sampled set of directions.

## Interior benchmark

The benchmark fixes terminal effect weights \((0.92,0.64,0.44,0)\), prefix
order \(a_0\geq a_1\geq a_2\geq a_3\), and \(\lambda=0.55\). An independently
generated outer superlevel box is

\[
\begin{aligned}
a_0&\in[0.296875,0.42596435546875],\\
a_1&\in[0.224609375,0.34832000732421875],\\
a_2&\in[0.15234375,0.258392333984375],\\
a_3&\in[0.1083984375,0.201324462890625].
\end{aligned}
\]

Order makes the scalar parts of the first two nontrivial characters
nonnegative, so the exact spectral cover has twelve regimes rather than
twenty-seven. In the scalar–vector–vector regime, the facially reduced Choi
moment outer gives \(0.75115255\); an independent SCS run gives
\(0.75116220\). Both are below \(0.758\), but they cover only that regime.

The fully vector-active regime remains the obstruction. A reduced SOCP with
16 planar sectors and 384 proved spherical caps solved 6,144 cells and
returned \(0.76293706\). At its worst point, the three character contractions
are nearly saturated, yet the flagged channel expands
\(\rho_2-\rho_3\): its output norm is \(0.32304\) while its input norm is
\(0.22977\). Adding that pairwise contraction moves the violation first to
\(\rho_1-2\rho_2\), then to \(\rho_0-2\rho_3\). Hence the three group
characters are necessary but not sufficient to enforce a common instrument.

The exact operator-basis reconstruction makes the failure stronger. The four
input states at the worst cap have determinant \(0.0084665\) and condition
number \(3.55\), so this is not a singular-input artifact. The unique four
interpolating subchannels reproduce all outputs to \(5.2\times10^{-17}\) and
are trace preserving in sum to \(1.8\times10^{-16}\), yet their minimum Choi
eigenvalues are \((-0.206,-0.218,-0.271,-0.275)\). The relaxation has retained
positivity and several data-processing shadows while violating complete
positivity in every outcome.

This is a useful negative result about the relaxation, not a physical
counterexample. The exact nondegenerate replacement is the degree-four
polynomial Choi matrix inequality derived in
`operator_basis_instrument_criterion.md`. For a planar terminal POVM, its
positive and negative determinant branches are related by an exact
orientation-reversing symmetry, so only one sign and the singular stratum need
independent treatment. A physical constrained seesaw establishes both sign
branches at score \(0.72288131\).

The remaining mathematical task is to exploit the polynomial matrix
inequality with a strong matrix localizer or verified spatial cover and to
handle the singular input stratum separately. Unassisted scalar coefficient
tests cannot certify complete positivity; their failure is now diagnosed,
not merely observed.

## Reproducibility map

The large moment implementation is in
`scratch/d2_frontier/choi_moment_reduced_upper.py`. The small parameterized
SOCP is in `scratch/d2_frontier/fourier_behavior_upper.py`; its cap-cover
driver is `scratch/d2_frontier/fourier_behavior_cap_cover.py`. The exact basis
audit is reproduced by `reproduce_operator_basis_obstruction.py`, and
`make_planar_conjugate_seed.py` constructs the opposite determinant-sign
strategy. The compact canonical result is
`scratch/d2_frontier/fourier_interior_summary_l055.json`; the full 6,144-cell
JSON is deliberately regenerated rather than committed. These files report
solver-conditional conic bounds. No claim of interval-verified numerics is
made.

# CarmenQ exact certificates

This Lean project is the trusted checking layer for exact rational conic
certificates. Numerical solvers, rational-recovery programs, CVXPY, and the
Python exporter are discovery and serialization tools; none is trusted by the
final Lean arithmetic proof.

## Source cell 15818

`CarmenQExact/Source15818Data.lean` assembles a self-contained generated module
graph into one `CertificateData` value. The stationarity modules own the sparse
matrix columns and objective coordinates, while `Source15818DualData.lean`
owns the exact dual function; the assembler also supplies the right-hand side,
cone dimensions, and target. Each binary64 coefficient is interpreted as the
exact rational returned by `float.as_integer_ratio()`. Thus the Lean theorem is
about the rational program defined by `source15818Data`, not about a decimal
approximation and not directly about the upstream physical model.

The generated theorem

```lean
source15818Exact : CertificateProof source15818Data
```

kernel-checks the cone dimensions, encoded-array indices, product-dual-cone
membership, all 274 stationarity equations, and the strict exact bound

$$
\operatorname{certifiedUpper}=
\frac{1161625617776817804636132569745984328794723998120992987}
     {1532495540865888858358347027150309183618739122183602176}
<\frac{379}{500}.
$$

`EncodedWeakDuality.lean` gives the sparse arrays their mathematical meaning.
Applying `exactCertificate_decoded_weak_duality` or
`exactCertificate_decoded_strict_target` to `source15818Exact` proves the
corresponding upper bound, respectively strict target inequality, for every
point and slack feasible for the decoded `CertificateData` program and its
encoded product cone. The public instance module exports those applications
directly as `source15818DecodedWeakDuality` and
`source15818DecodedStrictTarget`.

## Serialization boundary

`scratch/d2_frontier/source_15818_exact_socp_certificate.json` records source
identity, rational ray coefficients, and audit metadata, but the JSON alone
does not contain the complete $A,b,c$ arrays. The generated Source15818 module
graph contains the complete sparse rational $A,b,c$ data
and dual witness, with `Source15818Data.lean` assembling their named atoms into
the single object checked by Lean. The Python exporter rebuilds that graph from
the pinned source and checks it before emission; Lean then checks the emitted
object independently.

## Local semantic bridge

`CarmenQExact/Horwitz.lean` proves the exact monotonicities of all three
Horwitz weights and the theorem `source15818_horwitz_outer_enclosure`. For
every real `(alpha, beta)` in the stored source-15818 box, the theorem places
all three mathematical weights inside the six exact rational endpoints used
by the corrected canonical program. The Python builder obtains those binary64
endpoints by directed outward rounding; exact tests compare their rational
values without tolerances.

`TerminalReconstruction.lean` proves the stored coefficient intervals and
three-column error budget. `SpectralCaps.lean` proves all five source-15818 cap
containments, while `SpectralBranches.lean` proves the branch split is
exhaustive and packages the source branch into SOC premises.
`SourceConstraints.lean` proves the ordered-simplex, exact and concave
McCormick, Hellinger-SOC, Horwitz, inellipse tangent-remainder, and generic
outer-relaxation identities. `PhysicalBridge.lean` then states
`source15818_conditional_physical_frontier`: an explicit
`Source15818Embedding` into the decoded SOCP inherits the kernel-checked strict
target.

The physical theorem remains deliberately conditional. Its structure carries
the projective support bounds at weights `0.55` and `0.60` as named premises,
because the directed-interval projective covers are replayed in Python rather
than imported as Lean objects. The Python canonicalization audit maps all
1,142 explicit and 164 implicit source constraints to all 1,306 canonical
constraints and recompiles the 274 affine columns, but that executable audit
is not a Lean formalization of CVXPY. No solver claim is imported as an axiom.

## Global-frontier lemmas

`ProjectiveTangent.lean` proves both the affine Hellinger tangent used to
eliminate the remaining state variable and the scalar Rayleigh perturbation
lemma used near a zero coarse eigenvalue. The latter retains the certified
face quotient and therefore avoids an unnecessary second inverse spectral-gap
factor.

`CommonBias.lean` proves the support-cut bridge for a single shared Helstrom
bias. It also checks nonnegativity and instantiates the three rational
directions whose complete label orbits generate the 48 production cuts.
`SupportInterpolation.lean` records the general convex interpolation rule;
it deliberately assumes no numerical endpoint. `GlobalFrontier.lean` checks
the exact five-sector case split at `lambda = 3/5`, including

```text
76652/100000 + (3/5)(3/10000) = 76670/100000.
```

The same module checks the rational physical-witness endpoint and the exact
width of the resulting enclosure:

```text
0.7658988152 <= 0.76670,
0.76670 - 0.7658988152 = 0.0008011848.
```

These modules formalize the analytic bridges and final rational arithmetic.
They do not claim that Lean has replayed the interval trees, reconstructed the
CVXPY cone matrices, or formalized the physical reduction into those cones.

## Verification

Run the memory-bounded verifier from the repository root:

```console
python scripts/verify_lean_exact.py --workers 3 --heavy-workers 1
```

Use `--lake PATH` when `lake` is not on `PATH`. The driver rejects Lean trust
escapes, builds the checker and weak-duality layers, verifies the shared dual
and 274 column atoms before assembling the data, checks the bridge and cone
shards in bounded batches, then builds the aggregate source-15818 theorem and
the top-level `CarmenQExact` module. The conservative scan covers every production
Lean source file, including comments and the root module.

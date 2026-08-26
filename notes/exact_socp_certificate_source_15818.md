# Exact SOCP certificate for source cell 15818

## Status in one sentence

`source15818Exact : CertificateProof source15818Data` is a kernel-checked Lean
theorem for the rational conic program assembled as `source15818Data`. Its
generated module graph interprets every binary64 coefficient of the
reconstructed canonical source-cell 15818 program as its exact rational value
and proves a strict objective bound below (379/500=0.758).
`source15818DecodedWeakDuality` gives the exact upper bound for every feasible
point of the decoded program, and `source15818DecodedStrictTarget` gives the
strict target inequality.
`source15818_conditional_physical_frontier` now assembles the corresponding
physical statement when supplied with an explicit feasible embedding and two
named projective-line premises. Those premises are not hidden: their current
artifacts are numerical SCIP covers. This is therefore not an unconditional
proof of the full physical frontier theorem.

This distinction is part of the result, not a disclaimer to be removed later.
The formal certificate closes the floating-point arithmetic and conic
weak-duality gap for one literal `CertificateData` program. It does not by
itself establish that every upstream reduction, cap, enclosure, source-cell
construction, or canonicalizer transformation faithfully contains every
physical strategy that the proposed theorem is meant to cover.

## Exact statement certified

The stored source record is identified by all of the following data:

| field | value |
|---|---|
| frontier artifact | `scratch/d2_frontier/ternary_reconstructed_depth4_g2_top_leaf_bbb_p1_s92_l055.json` |
| SHA-256 | `8d314683c074d0aa59f5cac2677941f908d4d22350ed8187165f1f643a005884` |
| source index | `15818` |
| source cell | `608` |
| branch pattern | `bloch, bloch, scalar-negative, bloch` |
| cap indices | `8, 53, null, 8` |

With the repository's source-to-oracle code and canonicalization environment,
that record produces a cone program with 274 variables, 3,173 rows, and 20,962
nonzero matrix coefficients. Its cone is a product of a six-dimensional zero
cone, a 793-dimensional nonnegative cone, and 788 Lorentz cones. There are no
positive-semidefinite, exponential, or power-cone blocks in this particular
canonical program.

The exporter converts each binary64 entry of this canonical (A,b,c) to the
exact rational returned by `float.as_integer_ratio()` and writes the complete
sparse program into a generated Source15818 module graph.
`formal/CarmenQExact/Source15818Data.lean` assembles those named atoms into one
`CertificateData` value. The Lean theorem quantifies over the program decoded
from that value; it does not rerun or assume the correctness of the Python
reconstruction or CVXPY canonicalizer.

Write the Clarabel-form primal as

$$
  \min_x c^\mathsf{T}x
  \quad\text{subject to}\quad
  Ax+s=b,\qquad s\in K.
$$

The project score for this program is $-c^\mathsf{T}x$. The archived dual
vector is represented as an exact rational conic combination of sparse rays,

$$
  z=\sum_j \alpha_j r_j.
$$

Coefficients of zero-cone dual rays are unrestricted. Coefficients of
nonnegative and Lorentz-cone rays are nonnegative, and every Lorentz ray is
checked by the exact inequality $t^2\geq\lVert u\rVert_2^2$, with $t\geq0$.
The exact audit verifies

$$
  c+A^\mathsf{T}z=0.
$$

Therefore, for every feasible (x,s), self-duality of the nonnegative and
Lorentz cones gives

$$
  -c^\mathsf{T}x
  =z^\mathsf{T}Ax
  =z^\mathsf{T}(b-s)
  \leq b^\mathsf{T}z.
$$

The exact archived value is

$$
  U=b^\mathsf{T}z
  =\frac{
    1161625617776817804636132569745984328794723998120992987
  }{
    1532495540865888858358347027150309183618739122183602176
  }
  \approx 0.7579960833820616.
$$

Consequently,

$$
  U < \frac{379}{500},
  \qquad
  \frac{379}{500}-U
  =\frac{
    750274940743749936809604243754048535032061772182801
  }{
    191561942608236107294793378393788647952342390272950272000
  }
  \approx 3.916617938450015\times10^{-6}.
$$

All coefficients used in these equalities are exact rationals. A canonical
binary64 coefficient is interpreted by `float.as_integer_ratio()`, so neither
the exporter nor Lean silently treats a decimal printout as exact.

The source builder now rounds every proof-producing terminal interval endpoint
outward. `source15818_horwitz_outer_enclosure` independently proves in Lean
that the six stored rational endpoints enclose the three exact Horwitz weights
for every real parameter pair in the entire source-15818 terminal box.
Completed-square SOC radii are additionally advanced by 32 binary64 ULPs.
The independent `Fraction` audit checks the resulting perturbation rather than
assuming it harmless: all 29 anchor relaxations are outer, with minimum exact
margin

$$
\frac{13479967549230269}
{81129638414606681695789005144064}>0,
$$

and the coefficientwise lower completion has exact excess

$$
-\frac{187039510876414829}
{324518553658426726783156020576256}<0.
$$

The Lean bridge layer proves the Horwitz bounds, terminal reconstruction
interval arithmetic and three-column error budget, all five source-cell
spectral-cap containments, and exhaustiveness of the scalar-positive,
scalar-negative, and Bloch-dominant branches. It also proves the generic
ordered-simplex, McCormick, Hellinger-SOC, Horwitz inellipse, tangent-remainder,
and coefficientwise-relaxation identities used by the source model.

The concrete Lean proof checks the dimension sum; the objective and column
array sizes; membership of the witness in the product dual cone; all 274 exact
stationarity equations; and the displayed strict rational bound. The decoded
weak-duality layer then proves, for every point and slack satisfying the
literal sparse equation (Ax+s=b) and the encoded product-cone condition, that
its score is at most (U) and therefore strictly below (379/500).

## Trust chain

The result has two deliberately separate trust chains. For the theorem about
the literal encoded program, the generated Source15818 module graph contains
complete sparse rational (A,b,c), cone dimensions, target, and dual data;
`Source15818Data.lean` assembles them into the checked object. The generated
`source15818Exact` proof is checked by Lean's kernel and establishes dimensions,
indices, product-cone membership, 274 stationarity identities, and the strict
upper bound. `EncodedWeakDuality.lean` decodes the sparse arrays, proves that
the checked fold expressions equal the corresponding finite sums, and derives
weak duality and the strict target for every feasible point of that decoded
program.

Clarabel, HiGHS, FLINT, CVXPY, and the Python exporter are not premises of this
literal-data theorem. They discover and emit a candidate. If their proposed
rational witness fails any encoded equality or cone inequality, the kernel
proof cannot be completed. The formal proof uses `by decide +kernel`; it does
not use `native_decide`, `sorry`, `admit`, a solver oracle, or a custom axiom.

A different chain supplies provenance for why this particular literal program
is scientifically relevant. The certificate binds the stored source artifact
by SHA-256 and record identifiers. The exporter reruns
`build_localisation_oracle`, `set_cell_caps`, and CVXPY canonicalization,
interprets every resulting binary64 coefficient exactly, verifies the source
metadata and canonical digest, and emits the formal object only after its own
independent rational audit succeeds.

The separate canonicalization audit follows the actual reduction chain
`FlipObjective -> Dcp2Cone -> CvxAttr2Constr -> ConeMatrixStuffing -> CLARABEL`.
It maps every one of the 1,142 explicit source constraints plus 164 implicit
domain constraints to the 1,306 canonical constraints and all 3,173 rows. It
then evaluates the canonical expressions at zero and on every one of the 274
coordinate basis vectors, recovering the archived right-hand side and
objective exactly and the sparse matrix to maximum absolute discrepancy
`1.1102230246251565e-16`. The canonical-data SHA-256 is
`0861e28c987a2fdf03864ec8f753f70698e8cd3e8ba3b241ba715379acf0f1cf`;
the provenance-manifest SHA-256 is
`6f1ec1704c952520b3f677d8ad2388bd5730d3080e9ad98c8f0243e390b14064`.
This is strong executable evidence, not a kernel proof of CVXPY's
canonicalizer.

### Important serialization boundary

The JSON certificate is not self-contained: it stores source identity,
canonical dimensions, selected rays, rational ray coefficients, the objective
upper and margin, and audit flags, but it omits the full canonical (A,b,c).
Consequently, the JSON alone cannot support a fresh 274-coordinate
stationarity check without regenerating the canonical program.

The formal object has a different boundary. The generated Source15818 module
graph transitively encodes all 20,962 sparse matrix coefficients, all nonzero
right-hand-side entries, the complete objective, cone layout, exact functional
dual witness, and target. `Source15818Data.lean` is the small assembler for that
graph. Together with the checker and weak-duality modules, it is a self-contained
formal artifact: a third party can check the theorem without trusting a
numerical solver or regenerating (A,b,c) from the JSON.

Self-contained arithmetic is not the complete semantics. Lean now proves the
local cap, reconstruction, branch, and generic source-constraint lemmas and
the theorem for the data written in `CertificateData`. The executable audit
links every source constraint to canonical rows and recompiles their affine
action. Lean still does not verify CVXPY's reduction implementation, nor does
the repository contain exact witnesses for the two projective support covers
that delimit the source box. `PhysicalBridge.lean` therefore exposes those
support bounds and the physical-to-canonical embedding as fields, so no
unproved numerical statement enters the kernel as an axiom.

## What the certificate does not prove

The exact result should not be cited as any of the following:

- a kernel proof of CVXPY's canonicalization or of the recorded SCIP dual
  bounds for the projective support lines at weights `0.55` and `0.60`;
- an unconditional proof that every physically admissible common instrument
  has the `Source15818Embedding` required by the conditional Lean theorem;
- a proof for all source cells, an entire
  terminal leaf or strip, the whole interior frontier, or the global benchmark;
- an exact symbolic theorem for ideal physical constants—the certified
  coefficients are the exact rational values of binary64 canonical data;
- a proof that (U) is the exact optimum—the certificate supplies a rigorous
  feasible-dual upper bound, not dual attainment or a matching primal value;
- a state--Choi/PPT or semidefinite certificate. Such constraints may motivate
  upstream localizers and caps, but this selected canonical cell contains only
  zero, nonnegative, and second-order cones;
- a proof that a common quantum instrument is the only physical explanation of
  an experimental observation, or any claim about quantum interpretations,
  branches, observers, or consciousness.

The correct citation is: **a kernel-checked exact strict upper bound, with
decoded weak duality for every feasible point, for the rational canonical SOCP
literally encoded in `source15818Data`; kernel-checked local cap,
reconstruction, branch, and source-constraint lemmas; exact-rational enclosure
and provenance audits; and an end-to-end physical theorem whose projective
support and embedding assumptions are explicit. The unconditional global
frontier remains open.**

## Archived evidence and verification

The proof-related files are:

- `scratch/d2_frontier/source_15818_exact_socp_certificate.json`: source
  identity, sparse exact ray decomposition, rational upper, and audit metadata;
- `scratch/d2_frontier/export_exact_socp_certificate_lean.py`: deterministic
  reconstruction, independent exact audit, and Lean-data/proof generator;
- `formal/CarmenQExact/Horwitz.lean`: exact parameter monotonicities and the
  source-15818 theorem enclosing all three terminal weights in their directed
  binary64 intervals;
- `formal/CarmenQExact/TerminalReconstruction.lean`, `SpectralCaps.lean`, and
  `SpectralBranches.lean`: exact reconstruction enclosures, cap containment,
  and exhaustive branch selection;
- `formal/CarmenQExact/SourceConstraints.lean`: exact generic simplex,
  McCormick, Hellinger-SOC, inellipse, and outer-relaxation lemmas;
- `formal/CarmenQExact/PhysicalBridge.lean`: the data-carrying conditional
  physical theorem with both projective-line assumptions exposed;
- `formal/CarmenQExact/Source15818DualData.lean` and the numbered stationarity
  modules: the shared exact dual, sparse column, and objective atoms;
- `formal/CarmenQExact/Source15818Data.lean`: the assembler for the complete
  literal rational canonical program and dual witness;
- `formal/CarmenQExact/Checker.lean`, `WeakDuality.lean`, and
  `EncodedWeakDuality.lean`: exact checking, cone weak duality, and decoding of
  the sparse formal data;
- `formal/CarmenQExact/Source15818Exact.lean`: the assembled concrete
  `CertificateProof`;
- `formal/CarmenQExact/Source15818.lean`: the public decoded weak-duality and
  strict-target corollaries for this concrete program;
- `scripts/verify_lean_exact.py`: memory-bounded kernel verification and trust
  scan;
- `scratch/d2_frontier/exact_socp_dual_certificate.py` and
  `exact_socp_dual_certificate_flint_positive.py`: exact auditing and the
  rational recovery route used to discover the archived witness;
- `scratch/d2_frontier/audit_source15818_enclosures.py`: exact-rational audit
  of all anchor and coefficientwise completed-square SOCs;
- `scratch/d2_frontier/audit_source15818_canonicalization.py`: complete
  source-constraint provenance and independent affine-expression evaluation;
- `tests/test_exact_socp_dual_certificate.py`: toy audit/rejection tests,
  archive invariants, stored-ray cone checks, and source-hash verification.
- `tests/test_source15818_semantic_audits.py`: regression tests for both
  semantic audits and their pinned digests.

To verify the formal result from the repository root, run:

```console
python scripts/verify_lean_exact.py --workers 3 --heavy-workers 1
```

The driver first rejects explicit trust escapes throughout the production Lean
sources. It then builds the checker and decoded weak-duality theory, the shared
dual and 274 stationarity/column atoms, the data assembler, the bounded data
bridges and cone shards, and finally `Source15818Exact`, its public decoded
corollaries, and the top-level package. Pass `--lake PATH` if `lake` is not
available on `PATH`; lower either worker count on a memory-constrained machine.

The focused Python tests remain useful:

```console
python -m pytest tests/test_exact_socp_dual_certificate.py \
  tests/test_source15818_semantic_audits.py -q
```

Those tests deliberately do not reconstruct the 274-dimensional stationarity
identity from the JSON alone, because the JSON omits (A,b,c). That limitation
belongs to the JSON test path, not to the Lean artifact: the generated
Source15818 module graph contains the full arrays, and the kernel checks all
274 identities.

To reconstruct the canonical source program, repeat every exact exporter audit,
and deterministically regenerate the Lean modules, run:

```console
python scratch/d2_frontier/export_exact_socp_certificate_lean.py
```

To rediscover a rational ray decomposition rather than consume the archived
one, install the pinned `reproducibility` and `exact` dependencies and write to
a new JSON file:

```console
python scratch/d2_frontier/exact_socp_dual_certificate_flint_positive.py \
  --frontier-json scratch/d2_frontier/ternary_reconstructed_depth4_g2_top_leaf_bbb_p1_s92_l055.json \
  --source-index 15818 \
  --output scratch/d2_frontier/source_15818_exact_socp_certificate.regenerated.json
```

Discovery uses floating-point solvers, but emission occurs only after the
independent rational audit succeeds. A rediscovered certificate can use a
different redundant ray decomposition and still prove the same strict target;
compare such certificates semantically rather than requiring identical JSON.

## Next theorem gate

The next meaningful proof obligation is no longer another generic local lemma
or a tighter solver tolerance. It is an exact certificate format for the 128
projective leaves at weight `0.55` and the corresponding cover at weight
`0.60`, followed by a kernel-checked union of all required source cells. The
present result is an exact local computational theorem plus a transparent
conditional physical bridge; closing those external support covers is the
remaining gate to an unconditional frontier statement.

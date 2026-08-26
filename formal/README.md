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
\frac{2323251235560640692765186324170448172025567705991486731}
     {3064991081731777716716694054300618367237478244367204352}
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

This does not formalize the semantic bridge from the physical model through
the source construction, enclosure arguments, and CVXPY canonicalizer to the
literal `CertificateData`. Consequently, this certificate is an exact theorem
for one encoded canonical program, not a proof of the full physical frontier.

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

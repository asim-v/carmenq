# Interleaved-frontier research scripts

This directory contains both exploratory optimisers and the curated finite
certificates used to challenge and bound the compact three-effect construction
in the focused temporal-order manuscript. Exploratory outputs remain research
diagnostics; the certificate scripts state their exact scope and numerical
solver assumptions explicitly. The dependency-light construction check is
`../../scripts/verify_compact_interleaved_candidate.py`; the public lower-bound
implementation is in `../../src/carmenq/order_sensitive.py`.

The curated search chain is:

- `general_single_leaf_bound.py`: unrestricted complex bond-two Choi-MPS
  multistart search;
- `general_leaf_continuation.py`: continuation of one Choi-MPS checkpoint over
  support directions;
- `pauli_complete_general_leaf.py`: independent row-gauge and local Pauli
  completion checks for a stored leaf;
- `cq_instrument_relaxation.py`: enlarged pinched classical-quantum instrument
  relaxation;
- `suffix_cq_relaxation.py`: further relaxation of the first two slots to an
  arbitrary four-state qubit ensemble; and
- `reduced_three_effect_frontier.py`: direct two-parameter optimisation of the
  compact physical family;
- `reduced_four_effect_frontier.py`: the stronger symmetric four-effect
  physical phase;
- `verify_four_effect_candidate.py`: independent NumPy contraction of that
  phase;
- `oneway_exact_topologies.py`: direct optimisation of all four extreme
  binary-split topologies; and
- `four_effect_scip_certificate.py`: spatial global certificate for the
  symmetric four-parameter family;
- `projective_secular_scip.py`: global secular model for any one of the four
  projective topologies;
- `projective_trace_relaxation_scip.py`: cheap trace-only outer relaxation
  used to prune the finite projective cover;
- `validate_projective_cover.py`: coverage checker and compact manifest for
  the complete projective-sector certificate at `lambda=0.6`;
- `validate_projective_line_l055.py`: independent 128-leaf secular-cover
  validator for the auxiliary projective line at `lambda=0.55`;
- `continuous_terminal_projective_envelope_cover.py`: resumable dyadic SOCP
  cover of the full sorted ternary terminal-weight strip;
- `four_active_projection_scip.py`: projected Helstrom spatial relaxation for
  the genuinely four-active region;
- `audit_common_instrument_candidate.py`: exact fixed-input Choi projection,
  flagged trace-norm audit, and separating-witness archive for a relaxed
  first-moment family;
- `validate_common_instrument_strengthening.py`: validates the multi-scale
  and branch-local common-instrument experiment at `lambda=0.55`;
- `common_instrument_sparse_order2.py`: sparse second-order state--Choi
  hierarchy with shared outcome and trace-preservation bridges;
- `validate_common_instrument_hierarchy.py`: checks the multicolumn
  obstruction, sparse Choi residuals, exact first-moment incompatibility, and
  the independent deterministic formulations;
- `common_instrument_branch_tree.py`: resumable state-cell cover using robust
  Choi witnesses and deterministic cellwise trace-contraction bounds;
- `behavior_disjunction_scip.py`: exact mixed-integer union of small-support
  qubit-behaviour obstruction clauses;
- `common_instrument_exact_scip.py`: nonconvex spatial model whose sixteen
  outputs all use one literal collection of positive Choi matrices;
- `ternary_common_povm_input_cover.py`: localises the complete input basis of
  a selected ternary frontier cell and spatially refines one shared effective
  POVM;
- `ternary_common_instrument_input_cover.py`: robust shared-Choi tube over an
  input box, with inherited fixed-input separating witnesses;
- `ternary_bilinear_instrument_input_cover.py`: convergent input--Choi and
  input--POVM McCormick model with pure-prefix caps, exact multi-affine
  determinant signs, robust Cramer witnesses, optional planar Ando cuts, and
  optional common-instrument conic-RLT product localizers;
- `product_localizer_ablation.py`: identical-cell ablation of scalar product
  sums, Choi product traces, and matrix-order product sandwiches;
- `summarize_product_localizer_cover.py`: strict-JSON validation of the
  target-closing localizer tree and its sandwich-only control;
- `summarize_determinant_povm_cover.py`: validates a large checkpoint's tree
  accounting and emits the compact committed determinant-witness audit;
- `summarize_ando_instrument_cover.py`: validates the common-instrument
  witnesses and emits the compact equal-budget comparison against the
  determinant-only baseline;
- `validate_exact_common_instrument.py`: direct matrix audit and inward repair
  of the exact model's archived physical strategy;
- `active_readout_geometry_probe.py` and
  `four_active_geometry_support_probe.py`: exact Helstrom reconstruction,
  projective-comparison, polygon-closure, and prior-reserve helpers;
- `validate_interior_frontier.py`: combined sector-exhaustion validator for
  the complete numerical enclosure at `lambda=0.6`;
- `general_two_block_leaf.py` and `analyze_two_block_leaf.py`: unrestricted
  two-block relaxation and gauge-invariant diagnostics, including exact
  terminal qubit discrimination as a weighted smallest-ball SOCP;
- `verify_schmidt_norm_reduction.py`: reconstructs one two-block checkpoint as
  a single Schmidt-rank-two vector and checks the rank-one RETURN identity;
- `convexification_barrier.py`: exact diagonal-plus-rank-one obstruction to
  mixed Schmidt-number relaxations that retain concave Hellinger RETURN;
- `coherence_polar_program.py`: optimal path-pair polar contractions for the
  block-coherence extension that preserves the pure-leaf optimum;
- `planar_cp_completion.py`: exact Ando numerical-radius test for whether
  three planar pulled-effect columns share one completely positive map;
- `joint_effect_helstrom_seesaw.py`: generates the Helstrom-constrained
  columnwise outer checkpoints audited by the planar completion test;
- `qubit_discrimination_geometry.py`: lightweight implementation of that
  weighted-ball Helstrom dual, shared by the analyzer and Choi seesaw;
- `two_block_choi_seesaw.py`: finite common-instrument formulation with four
  qubit Choi matrices and globally solved state/instrument/POVM blocks;
- `two_block_random_povm_probe.py`, `two_block_rank_one_povm.py`, and
  `two_block_fixed_povm.py`: nonprojective-sector falsification searches (the
  random probe handles the physical/transposed POVM convention explicitly);
  and
- `syndrome_prior_upper.py`: a deliberately loose terminal-dimension upper
  bound.

The first five scripts, the general two-block optimiser, and its checkpoint
analyzer require PyTorch in addition to the package dependencies. The Choi
seesaw itself needs PyTorch only when a `.pt` checkpoint is supplied.
Run them from this directory so their local research imports resolve:

```bash
python -m pip install torch
python general_single_leaf_bound.py --lambda 0.5 --restarts 8 \
  --output general_leaf_reproduction.json
python pauli_complete_general_leaf.py general_leaf_reproduction.pt --lambda 0.5
python cq_instrument_relaxation.py --lambda 0.5 --restarts 8 \
  --checkpoint general_leaf_reproduction.pt --output cq_reproduction.json
python suffix_cq_relaxation.py --lambda 0.5 --restarts 6 \
  --checkpoint general_leaf_reproduction.pt
python reduced_three_effect_frontier.py --sweep \
  --output reduced_three_effect_reproduction.json
python reduced_four_effect_frontier.py --lambda 0.6 \
  --output reduced_four_effect_reproduction.json
python verify_four_effect_candidate.py --lambda 0.6 \
  --output four_effect_verification.json
python validate_projective_cover.py \
  --output projective_cover_l060_summary.json
python validate_projective_line_l055.py
python validate_interior_frontier.py
python general_two_block_leaf.py --lambda 0.6 --steps 5000 --restarts 12 \
  --output general_two_block_reproduction_l060.json
python analyze_two_block_leaf.py general_two_block_reproduction_l060.pt \
  --lambda 0.6
python two_block_choi_seesaw.py --lambda 0.6 --rounds 8 --restarts 4 \
  --output two_block_choi_reproduction_l060.json
python ternary_bilinear_instrument_input_cover.py \
  --frontier-json ternary_reconstructed_depth4_g2_top_leaf_bbb_p1_s92_l055.json \
  --output ternary_exactdet_hybrid_instrument_topcell_pilot_l055.json \
  --top-spectral-cell --max-nodes 1000 --checkpoint-every 20
python summarize_determinant_povm_cover.py \
  ternary_exactdet_hybrid_instrument_topcell_pilot_l055.json \
  --output determinant_povm_cover_l055_summary.json

python ternary_bilinear_instrument_input_cover.py \
  --frontier-json ternary_reconstructed_depth4_g2_top_leaf_bbb_p1_s92_l055.json \
  --localisation-json ternary_exactdet_hybrid_instrument_topcell_pilot_l055.json \
  --output ternary_exactdet_ando_guided_instrument_topcell_pilot100_l055.json \
  --top-spectral-cell --planar-ando-witnesses \
  --max-nodes 1000 --checkpoint-every 25

python summarize_ando_instrument_cover.py \
  ternary_exactdet_ando_guided_instrument_topcell_pilot100_l055.json \
  --baseline-summary determinant_povm_cover_l055_summary.json \
  --output ando_instrument_cover_l055_summary.json

python product_localizer_ablation.py \
  --frontier-json ternary_reconstructed_depth4_g2_top_leaf_bbb_p1_s92_l055.json \
  --checkpoint ternary_exactdet_ando_guided_instrument_topcell_pilot100_l055.json \
  --limit 20 --output product_localizer_ablation_top20_l055.json

python ternary_bilinear_instrument_input_cover.py \
  --frontier-json ternary_reconstructed_depth4_g2_top_leaf_bbb_p1_s92_l055.json \
  --localisation-json ternary_exactdet_hybrid_instrument_topcell_pilot_l055.json \
  --output _ternary_exactdet_ando_sandwichonly_instrument_topcell_pilot100_l055.json \
  --top-spectral-cell --planar-ando-witnesses \
  --common-instrument-product-psd-sandwiches \
  --max-nodes 100 --checkpoint-every 10

python ternary_bilinear_instrument_input_cover.py \
  --frontier-json ternary_reconstructed_depth4_g2_top_leaf_bbb_p1_s92_l055.json \
  --localisation-json ternary_exactdet_hybrid_instrument_topcell_pilot_l055.json \
  --output ternary_exactdet_ando_matrixlocalizer_instrument_topcell_pilot100_l055.json \
  --top-spectral-cell --planar-ando-witnesses \
  --common-instrument-product-trace-rules \
  --common-instrument-product-psd-sandwiches \
  --max-nodes 100 --checkpoint-every 10

python summarize_product_localizer_cover.py \
  ternary_exactdet_ando_matrixlocalizer_instrument_topcell_pilot100_l055.json \
  --baseline-summary ando_instrument_cover_l055_summary.json \
  --sandwich-checkpoint _ternary_exactdet_ando_sandwichonly_instrument_topcell_pilot100_l055.json \
  --output product_localizer_cover_l055_summary.json
```

The committed JSON files record representative double-precision searches.
Binary PyTorch checkpoints are omitted because they are
framework-version-specific and can be regenerated. The four-effect branch
strictly falsifies the earlier three-effect full-frontier conjecture at
`lambda=0.6`.

The final certificate no longer trusts the archived optimum values. Directed
projective interval covers establish `beta_projective(0.60) <= 0.76600` and
the auxiliary line `beta_projective(0.55) <= 0.7573`. Exact-residual conic
replay then bounds the capped-weight sector by `0.765893818`, all 12,008
ternary cells by `0.76652`, and the 90-leaf four-active tree by `0.76670` in
all six affine orders. Deleting an effect below `0.0003` gives the remaining
bound `0.76652 + 0.6(0.0003) = 0.76670`. A separate rational physical witness
certifies `0.7658988152` from below, so

```text
0.7658988152 <= beta_stream(0.6) <= beta_2b(0.6) <= 0.76670.
```

The conic replays call no optimiser and evaluate repaired residuals with exact
rational arithmetic. The projective kernels use outward-expanded binary64
intervals. Python matrix canonicalisation and the analytic physical reduction
remain in the documented trust boundary; this is solver independent in
verification, not end-to-end kernel formalisation. Assemble the final manifest
after the component files are present:

```bash
python verify_global_frontier_l060.py \
  --output ../../data/global_frontier_l060_exact_assembly.json
```

See `../../notes/interleaved_interior_frontier_l060.md` and the checksummed
global manifest in `../../data/global_frontier_l060_exact_assembly.json`.

The determinant-witness pilot at `lambda=0.55` is a strengthening experiment,
not another completed sector certificate.  Its 1,000-node run leaves 935
cells pending at maximum upper bound `0.7633741152384063`, above target
`0.758`.  The full checkpoint is regenerated by the commands above; only its
11 KB validated summary is committed.  See
`../../notes/determinant_scaled_common_povm_witnesses.md`.

The Ando-strengthened pilot converts exact planar common-instrument completion
conditions into box-valid determinant-scaled probability cuts. The equal-budget
run generates eight such cuts and lowers the maximum pending bound from
`0.7633741152384063` to `0.7633726473409502`. This is a strict gain, but it
does not close target `0.758` and the strengthened tree closes fewer nodes.
The compact audit is `ando_instrument_cover_l055_summary.json`; see
`../../notes/determinant_scaled_ando_instrument_witnesses.md`.

The common-instrument product localizers resolve that selected localized cell.
The matrix inequalities retain the Choi order of each lifted scalar--matrix
product, while exact partial-trace identities couple all outcome products to
one trace-preserving instrument.  A 20-cell ablation finds no material gain
from the trace rule alone, mean improvement `0.0132790` from the PSD sandwich,
and mean improvement `0.0265894` when both are imposed.  The sandwich-only
100-node control remains open at `0.7631986663448878`; the combined tree closes
in 49 nodes with solver-conditional cover upper bound `0.7579083037237451`,
below target `0.758`.  This is not a global frontier certificate.  See
`../../notes/common_instrument_product_localizers.md` and the strict compact
audit `product_localizer_cover_l055_summary.json`.

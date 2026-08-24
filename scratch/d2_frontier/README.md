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
```

The committed JSON files record representative double-precision runs. Binary
PyTorch checkpoints are deliberately omitted from the curated release because
they are framework-version-specific and can be regenerated. The four-effect
branch strictly falsifies the earlier three-effect full-frontier conjecture at
`lambda=0.6`. The projective terminal sector is globally enclosed in
`[0.765898815264694, 0.76591]` by the recorded SCIP cover. Equal values across
the unrestricted MPS and two-block searches initially left the genuine
three-/four-active terminal-POVM sectors open. Those sectors are now exhausted
directly at `lambda=0.6`, producing the complete solver-conditional enclosure
`[0.7658988152646944, 0.76662]`. The qualifier is essential: the logical
sector reduction is exact, but the upper endpoint uses CLARABEL values with a
safety margin and SCIP spatial duals at recorded tolerances rather than
solver-independent interval arithmetic. See
`../../notes/interleaved_interior_frontier_l060.md` and
`interior_frontier_l060_certificate.json`.

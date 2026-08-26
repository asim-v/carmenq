import CarmenQExact.Source15818
import CarmenQExact.SourceConstraints
import CarmenQExact.SpectralBranches
import CarmenQExact.TerminalReconstruction

namespace CarmenQExact

noncomputable section

open scoped BigOperators

/-!
# Conditional physical bridge for source cell 15818

The kernel-checked dual theorem concerns a literal canonical SOCP.  This file
states the missing semantic interface without adding unchecked declarations: a caller must
provide an explicit embedding of each physical strategy into that program.
The two projective support lines are carried as named premises because their
current repository covers are numerical SCIP certificates, not exact proofs.
-/

def WeightedSupportBound {Strategy : Type}
    (weight upper : ℝ) (audit returned : Strategy → ℝ) : Prop :=
  ∀ strategy,
    weight * audit strategy + (1 - weight) * returned strategy ≤ upper

structure ProjectiveFrontierData where
  Strategy : Type
  audit : Strategy → ℝ
  returned : Strategy → ℝ

structure Source15818ProjectivePremises
    (data : ProjectiveFrontierData) : Prop where
  line055 : WeightedSupportBound (11 / 20) (7573 / 10000)
    data.audit data.returned
  line060 : WeightedSupportBound (3 / 5) (76591 / 100000)
    data.audit data.returned

def source15818CanonicalScore
    (point : Fin source15818Data.variableCount → ℝ) : ℝ :=
  -(∑ column, decodedObjective source15818Data column * point column)

structure Source15818Embedding (score : ℝ) where
  point : Fin source15818Data.variableCount → ℝ
  slack : Fin source15818Data.rows → ℝ
  feasible : CanonicalFeasible
    (decodedMatrix source15818Data) (decodedRight source15818Data)
    point slack
  cone : ProductConePoint (extendFin slack)
    source15818Data.zeroDim source15818Data.nonnegativeDim
    source15818Data.socSizes.toList
  score_eq : score = source15818CanonicalScore point

theorem Source15818Embedding.strict_target {score : ℝ}
    (embedding : Source15818Embedding score) :
    score < (source15818Data.target : ℝ) := by
  rw [embedding.score_eq]
  exact source15818DecodedStrictTarget
    embedding.point embedding.slack embedding.feasible embedding.cone

structure PhysicalFrontierData where
  Strategy : Type
  score : Strategy → ℝ

structure Source15818ConditionalBridge
    (physical : PhysicalFrontierData)
    (projective : ProjectiveFrontierData) where
  projectivePremises : Source15818ProjectivePremises projective
  embed : ∀ strategy, Source15818Embedding (physical.score strategy)

theorem source15818_conditional_physical_frontier
    (physical : PhysicalFrontierData)
    (projective : ProjectiveFrontierData)
    (bridge : Source15818ConditionalBridge physical projective) :
    ∀ strategy, physical.score strategy < (source15818Data.target : ℝ) := by
  intro strategy
  exact (bridge.embed strategy).strict_target

end

end CarmenQExact

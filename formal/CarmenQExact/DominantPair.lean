import Mathlib

namespace CarmenQExact

/-!
# Pairwise closure bound for a four-effect qubit POVM

For any selected pair of rank-one effects, full Bloch closure says that its
weighted Bloch-vector resultant is the negative resultant of the complementary
pair.  The latter has norm at most its total weight.  The theorem below turns
that triangle bound into the exact loss factor used by every ordered-pair
projective comparison.  The historical theorem names are retained.
-/

theorem dominantPairNumeratorIdentity (w₀ w₁ : ℝ) :
    (2 - w₀ - w₁) ^ 2 - (w₀ - w₁) ^ 2 =
      4 * (1 - w₀) * (1 - w₁) := by
  ring

theorem dominantPairClosureBound
    (w₀ w₁ d : ℝ)
    (hw₀ : 0 < w₀) (hw₁ : 0 < w₁)
    (hclosure :
      w₀ ^ 2 + w₁ ^ 2 + 2 * w₀ * w₁ * d ≤
        (2 - w₀ - w₁) ^ 2) :
    (1 + d) / 2 ≤ (1 - w₀) * (1 - w₁) / (w₀ * w₁) := by
  have hproduct : 0 < w₀ * w₁ := mul_pos hw₀ hw₁
  apply (le_div_iff₀ hproduct).2
  have hid := dominantPairNumeratorIdentity w₀ w₁
  nlinarith

theorem dominantPairProjectiveLoss
    (A p actual upper : ℝ)
    (hgap : 0 ≤ A - p) (hfactor : actual ≤ upper) :
    A - upper * (A - p) ≤ A - actual * (A - p) := by
  nlinarith

end CarmenQExact

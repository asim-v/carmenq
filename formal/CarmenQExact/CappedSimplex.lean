import Mathlib

namespace CarmenQExact

/-!
# Four-coordinate capped-simplex support function

The low-terminal-weight certificate replaces an unknown four-effect weight
vector by `(u, u, 2 - 2*u, 0)` after sorting the nonnegative marginal on
which it acts.  The theorem below is the exact majorisation step.  It uses
only the coordinate caps, nonnegativity, total weight two, and the declared
order of the marginal; no optimiser or floating-point computation enters.
-/

theorem cappedSimplexFourSupport
    (q₀ q₁ q₂ q₃ w₀ w₁ w₂ w₃ u : ℝ)
    (hq₀₁ : q₁ ≤ q₀) (hq₁₂ : q₂ ≤ q₁) (hq₂₃ : q₃ ≤ q₂)
    (hw₀ : w₀ ≤ u) (hw₁ : w₁ ≤ u) (hw₃ : 0 ≤ w₃)
    (hsum : w₀ + w₁ + w₂ + w₃ = 2) :
    w₀ * q₀ + w₁ * q₁ + w₂ * q₂ + w₃ * q₃ ≤
      u * q₀ + u * q₁ + (2 - 2 * u) * q₂ := by
  have hfirst : 0 ≤ (u - w₀) * (q₀ - q₁) :=
    mul_nonneg (sub_nonneg.mpr hw₀) (sub_nonneg.mpr hq₀₁)
  have hsecond : 0 ≤ (2 * u - w₀ - w₁) * (q₁ - q₂) := by
    apply mul_nonneg
    · linarith
    · exact sub_nonneg.mpr hq₁₂
  have hthird : 0 ≤ w₃ * (q₂ - q₃) :=
    mul_nonneg hw₃ (sub_nonneg.mpr hq₂₃)
  calc
    w₀ * q₀ + w₁ * q₁ + w₂ * q₂ + w₃ * q₃
        ≤ w₀ * q₀ + w₁ * q₁ + w₂ * q₂ + w₃ * q₃
          + (u - w₀) * (q₀ - q₁)
          + (2 * u - w₀ - w₁) * (q₁ - q₂)
          + w₃ * (q₂ - q₃) := by linarith
    _ = u * q₀ + u * q₁ + (2 - 2 * u) * q₂ := by
      have hw₂ : w₂ = 2 - w₀ - w₁ - w₃ := by
        linarith
      rw [hw₂]
      ring

theorem capTailAt88325 :
    (2 : ℚ) - 2 * (3533 / 4000 : ℚ) = 467 / 2000 := by
  norm_num

end CarmenQExact

import CarmenQExact.SpectralCaps

namespace CarmenQExact

noncomputable section

/-! Spectral branch exhaustivity and the exact branch predicate for source
cell 15818.  For a Hermitian qubit contraction `s I + v·sigma`, its operator
norm is `max |s| ‖v‖`.  The three cases below are therefore exhaustive and
are precisely the scalar-positive, scalar-negative, and Bloch-dominant cases
used by the branch tree.
-/

def ScalarPositiveBranch (scalar : ℝ) (bloch : Bloch3) : Prop :=
  ‖bloch‖ ≤ scalar

def ScalarNegativeBranch (scalar : ℝ) (bloch : Bloch3) : Prop :=
  ‖bloch‖ ≤ -scalar

def BlochDominantBranch (scalar : ℝ) (bloch : Bloch3) : Prop :=
  |scalar| ≤ ‖bloch‖

theorem spectral_branch_exhaustive (scalar : ℝ) (bloch : Bloch3) :
    ScalarPositiveBranch scalar bloch ∨
      ScalarNegativeBranch scalar bloch ∨
      BlochDominantBranch scalar bloch := by
  by_cases hpositive : ‖bloch‖ ≤ scalar
  · exact Or.inl hpositive
  by_cases hnegative : ‖bloch‖ ≤ -scalar
  · exact Or.inr (Or.inl hnegative)
  · right
    right
    unfold BlochDominantBranch
    rw [abs_le]
    constructor <;> linarith

def PhysicalFlagBound (scalar flagged : ℝ) (bloch : Bloch3) : Prop :=
  flagged ≤ max |scalar| ‖bloch‖

theorem physical_flag_scalar_positive {scalar flagged : ℝ} {bloch : Bloch3}
    (hphysical : PhysicalFlagBound scalar flagged bloch)
    (hbranch : ScalarPositiveBranch scalar bloch) :
    flagged ≤ scalar := by
  have hscalar : 0 ≤ scalar := (norm_nonneg bloch).trans hbranch
  simpa [PhysicalFlagBound, abs_of_nonneg hscalar, max_eq_left hbranch] using hphysical

theorem physical_flag_scalar_negative {scalar flagged : ℝ} {bloch : Bloch3}
    (hphysical : PhysicalFlagBound scalar flagged bloch)
    (hbranch : ScalarNegativeBranch scalar bloch) :
    flagged ≤ -scalar := by
  have hscalar : scalar ≤ 0 := by
    unfold ScalarNegativeBranch at hbranch
    linarith [norm_nonneg bloch]
  simpa [PhysicalFlagBound, abs_of_nonpos hscalar, max_eq_left hbranch] using hphysical

theorem physical_flag_bloch {scalar flagged : ℝ} {bloch : Bloch3}
    (hphysical : PhysicalFlagBound scalar flagged bloch)
    (hbranch : BlochDominantBranch scalar bloch) :
    flagged ≤ ‖bloch‖ := by
  simpa [PhysicalFlagBound, max_eq_right hbranch] using hphysical

theorem bloch_cap_flag {scalar flagged : ℝ} {bloch cap : Bloch3}
    (hphysical : PhysicalFlagBound scalar flagged bloch)
    (hbranch : BlochDominantBranch scalar bloch)
    (hcap : InScaledCap cap bloch) :
    flagged ≤ inner ℝ cap bloch :=
  (physical_flag_bloch hphysical hbranch).trans hcap

theorem norm_positive_z_axis {z : ℝ} (hz : 0 ≤ z) :
    ‖(!₂[0, 0, z] : Bloch3)‖ = z := by
  have hv : (!₂[0, 0, z] : Bloch3) = z • (!₂[0, 0, 1] : Bloch3) := by
    ext i
    fin_cases i <;> simp
  have hunit : ‖(!₂[0, 0, 1] : Bloch3)‖ = 1 := by
    rw [EuclideanSpace.norm_eq]
    norm_num [Fin.sum_univ_succ]
  rw [hv, norm_smul, Real.norm_eq_abs, abs_of_nonneg hz, hunit, mul_one]

structure Source15818SpectralPremises
    (scalar₀ scalar₁ scalar₂ scalar₃ scalar₄ scalar₅ scalar₆ : ℝ)
    (bloch₀ bloch₁ bloch₂ bloch₃ bloch₄ bloch₅ bloch₆ : Bloch3)
    (flagged₀ flagged₁ flagged₂ flagged₃ flagged₄ flagged₅ flagged₆ : ℝ) : Prop where
  angular : Source15818AngularCells bloch₁ bloch₂ bloch₃ bloch₄ bloch₆
  base0Axis : ∃ z : ℝ, 0 ≤ z ∧ bloch₀ = !₂[0, 0, z]
  dominant0 : BlochDominantBranch scalar₀ bloch₀
  dominant1 : BlochDominantBranch scalar₁ bloch₁
  dominant2 : BlochDominantBranch scalar₂ bloch₂
  dominant3 : BlochDominantBranch scalar₃ bloch₃
  dominant4 : BlochDominantBranch scalar₄ bloch₄
  negative5 : ScalarNegativeBranch scalar₅ bloch₅
  dominant6 : BlochDominantBranch scalar₆ bloch₆
  physical0 : PhysicalFlagBound scalar₀ flagged₀ bloch₀
  physical1 : PhysicalFlagBound scalar₁ flagged₁ bloch₁
  physical2 : PhysicalFlagBound scalar₂ flagged₂ bloch₂
  physical3 : PhysicalFlagBound scalar₃ flagged₃ bloch₃
  physical4 : PhysicalFlagBound scalar₄ flagged₄ bloch₄
  physical5 : PhysicalFlagBound scalar₅ flagged₅ bloch₅
  physical6 : PhysicalFlagBound scalar₆ flagged₆ bloch₆

structure Source15818SpectralSOCP
    (scalar₅ : ℝ)
    (bloch₀ bloch₁ bloch₂ bloch₃ bloch₄ bloch₅ bloch₆ : Bloch3)
    (flagged₀ flagged₁ flagged₂ flagged₃ flagged₄ flagged₅ flagged₆ : ℝ) : Prop where
  base0 : ∃ z : ℝ, 0 ≤ z ∧ bloch₀ = !₂[0, 0, z] ∧ flagged₀ ≤ z
  planeLorentz : InScaledCap source15818PlaneCap bloch₁
  planeFlag : flagged₁ ≤ inner ℝ source15818PlaneCap bloch₁
  sphereLorentz : InScaledCap source15818BaseSphereCap bloch₂
  sphereFlag : flagged₂ ≤ inner ℝ source15818BaseSphereCap bloch₂
  separator0Lorentz : InScaledCap source15818Separator0Cap bloch₃
  separator0Flag : flagged₃ ≤ inner ℝ source15818Separator0Cap bloch₃
  separator1Lorentz : InScaledCap source15818Separator1Cap bloch₄
  separator1Flag : flagged₄ ≤ inner ℝ source15818Separator1Cap bloch₄
  separator2Lorentz : ‖bloch₅‖ ≤ -scalar₅
  separator2Flag : flagged₅ ≤ -scalar₅
  separator3Lorentz : InScaledCap source15818Separator3Cap bloch₆
  separator3Flag : flagged₆ ≤ inner ℝ source15818Separator3Cap bloch₆

theorem source15818_spectral_cell_to_socp
    {scalar₀ scalar₁ scalar₂ scalar₃ scalar₄ scalar₅ scalar₆ : ℝ}
    {bloch₀ bloch₁ bloch₂ bloch₃ bloch₄ bloch₅ bloch₆ : Bloch3}
    {flagged₀ flagged₁ flagged₂ flagged₃ flagged₄ flagged₅ flagged₆ : ℝ}
    (h : Source15818SpectralPremises
      scalar₀ scalar₁ scalar₂ scalar₃ scalar₄ scalar₅ scalar₆
      bloch₀ bloch₁ bloch₂ bloch₃ bloch₄ bloch₅ bloch₆
      flagged₀ flagged₁ flagged₂ flagged₃ flagged₄ flagged₅ flagged₆) :
    Source15818SpectralSOCP scalar₅
      bloch₀ bloch₁ bloch₂ bloch₃ bloch₄ bloch₅ bloch₆
      flagged₀ flagged₁ flagged₂ flagged₃ flagged₄ flagged₅ flagged₆ := by
  have hcaps := source15818_all_angular_caps_contain h.angular
  rcases h.base0Axis with ⟨z, hz, hbloch₀⟩
  have hflag₀ : flagged₀ ≤ z := by
    calc
      flagged₀ ≤ ‖bloch₀‖ := physical_flag_bloch h.physical0 h.dominant0
      _ = z := by rw [hbloch₀, norm_positive_z_axis hz]
  exact {
    base0 := ⟨z, hz, hbloch₀, hflag₀⟩
    planeLorentz := hcaps.plane
    planeFlag := bloch_cap_flag h.physical1 h.dominant1 hcaps.plane
    sphereLorentz := hcaps.sphere
    sphereFlag := bloch_cap_flag h.physical2 h.dominant2 hcaps.sphere
    separator0Lorentz := hcaps.sep0
    separator0Flag := bloch_cap_flag h.physical3 h.dominant3 hcaps.sep0
    separator1Lorentz := hcaps.sep1
    separator1Flag := bloch_cap_flag h.physical4 h.dominant4 hcaps.sep1
    separator2Lorentz := h.negative5
    separator2Flag := physical_flag_scalar_negative h.physical5 h.negative5
    separator3Lorentz := hcaps.sep3
    separator3Flag := bloch_cap_flag h.physical6 h.dominant6 hcaps.sep3
  }

end

end CarmenQExact

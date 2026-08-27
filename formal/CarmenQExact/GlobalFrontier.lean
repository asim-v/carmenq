import Mathlib

namespace CarmenQExact

/-!
# Exact assembly of the lambda = 3/5 five-sector frontier

The numerical proof kernels certify five exhaustive terminal-readout sectors.
This file checks only their final rational arithmetic and the case split.  It
does not import a numerical artifact or claim to formalise the interval/SOCP
kernels themselves.
-/

noncomputable section

def projective055Level : ℝ := 7573 / 10000
def projective060Level : ℝ := 766 / 1000
def lowWeightLevel : ℝ := 76591 / 100000
def ternaryLevel : ℝ := 76652 / 100000
def fourActiveLevel : ℝ := 76670 / 100000
def deletionLevel : ℝ := 76670 / 100000
def globalFrontierLevel : ℝ := 76670 / 100000
def physicalLowerLevel : ℝ := 7658988152 / 10000000000
def certifiedIntervalWidth : ℝ := 1001481 / 1250000000

theorem deletionLevel_exact :
    (76652 : ℝ) / 100000 + (3 / 5) * (3 / 10000) = deletionLevel := by
  norm_num [deletionLevel]

theorem physicalLower_below_global :
    physicalLowerLevel ≤ globalFrontierLevel := by
  norm_num [physicalLowerLevel, globalFrontierLevel]

theorem certifiedIntervalWidth_exact :
    globalFrontierLevel - physicalLowerLevel = certifiedIntervalWidth := by
  norm_num [globalFrontierLevel, physicalLowerLevel, certifiedIntervalWidth]

theorem projective055_below_global :
    projective055Level ≤ globalFrontierLevel := by
  norm_num [projective055Level, globalFrontierLevel]

theorem projective060_below_global :
    projective060Level ≤ globalFrontierLevel := by
  norm_num [projective060Level, globalFrontierLevel]

theorem fourActive_below_global : fourActiveLevel ≤ globalFrontierLevel := by
  norm_num [fourActiveLevel, globalFrontierLevel]

theorem lowWeight_below_global : lowWeightLevel ≤ globalFrontierLevel := by
  norm_num [lowWeightLevel, globalFrontierLevel]

theorem ternary_below_global : ternaryLevel ≤ globalFrontierLevel := by
  norm_num [ternaryLevel, globalFrontierLevel]

theorem deletion_below_global : deletionLevel ≤ globalFrontierLevel := by
  norm_num [deletionLevel, globalFrontierLevel]

theorem terminalSectorPartition
    (activeCount : ℕ)
    (maximumWeight minimumWeight : ℝ)
    (harity : activeCount ≤ 4) :
    activeCount ≤ 2 ∨
      maximumWeight ≤ 3533 / 4000 ∨
      activeCount = 3 ∨
      (activeCount = 4 ∧ 3 / 10000 ≤ minimumWeight) ∨
      (activeCount = 4 ∧ minimumWeight < 3 / 10000) := by
  by_cases hprojective : activeCount ≤ 2
  · exact Or.inl hprojective
  by_cases hlowWeight : maximumWeight ≤ (3533 : ℝ) / 4000
  · exact Or.inr (Or.inl hlowWeight)
  have harityCases : activeCount = 3 ∨ activeCount = 4 := by omega
  rcases harityCases with hternary | hfour
  · exact Or.inr (Or.inr (Or.inl hternary))
  · by_cases hminimum : (3 : ℝ) / 10000 ≤ minimumWeight
    · exact Or.inr (Or.inr (Or.inr (Or.inl ⟨hfour, hminimum⟩)))
    · exact Or.inr (
        Or.inr (
          Or.inr (
            Or.inr ⟨hfour, lt_of_not_ge hminimum⟩)))

theorem globalFrontierFromFiveSectors
    (score projective lowWeight ternary fourActive smallEffect : ℝ)
    (hcover :
      score ≤ projective ∨
      score ≤ lowWeight ∨
      score ≤ ternary ∨
      score ≤ fourActive ∨
      score ≤ smallEffect)
    (hprojective : projective ≤ projective060Level)
    (hlowWeight : lowWeight ≤ lowWeightLevel)
    (hternary : ternary ≤ ternaryLevel)
    (hfourActive : fourActive ≤ fourActiveLevel)
    (hsmallEffect : smallEffect ≤ deletionLevel) :
    score ≤ globalFrontierLevel := by
  rcases hcover with h | h | h | h | h
  · exact h.trans (hprojective.trans projective060_below_global)
  · exact h.trans (hlowWeight.trans lowWeight_below_global)
  · exact h.trans (hternary.trans ternary_below_global)
  · exact h.trans (hfourActive.trans fourActive_below_global)
  · exact h.trans (hsmallEffect.trans deletion_below_global)

theorem fiveSectorMaximum_exact :
    max projective060Level
      (max lowWeightLevel
        (max ternaryLevel (max fourActiveLevel deletionLevel))) =
      globalFrontierLevel := by
  norm_num [projective060Level, lowWeightLevel, ternaryLevel, fourActiveLevel,
    deletionLevel, globalFrontierLevel, max_def]

end

end CarmenQExact

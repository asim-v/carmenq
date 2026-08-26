import Mathlib

namespace CarmenQExact

noncomputable section

/-! Exact monotonicity and interval enclosure for the Horwitz
parameterisation of a three-outcome rank-one qubit POVM.  This is the first
source-to-canonical bridge used by source cell 15818: it proves that the
stored binary64 endpoints are an outer enclosure, rather than relying on
nearest-rounding arithmetic. -/

def horwitzDenominator (alpha beta : ℝ) : ℝ := alpha + beta - 1

def horwitzW0 (alpha beta : ℝ) : ℝ :=
  alpha / horwitzDenominator alpha beta

def horwitzW1 (alpha beta : ℝ) : ℝ :=
  beta / horwitzDenominator alpha beta

def horwitzW2 (alpha beta : ℝ) : ℝ :=
  1 - 1 / horwitzDenominator alpha beta

theorem horwitzW0_mono_alpha
    (alpha₁ alpha₂ beta : ℝ)
    (halpha : 0 < alpha₁) (hbeta : 1 ≤ beta)
    (horder : alpha₁ ≤ alpha₂) :
    horwitzW0 alpha₁ beta ≤ horwitzW0 alpha₂ beta := by
  have hdenom₁ : 0 < horwitzDenominator alpha₁ beta := by
    simp only [horwitzDenominator]
    linarith
  have hdenom₂ : 0 < horwitzDenominator alpha₂ beta := by
    simp only [horwitzDenominator]
    linarith
  unfold horwitzW0
  rw [div_le_div_iff₀ hdenom₁ hdenom₂]
  simp only [horwitzDenominator]
  nlinarith

theorem horwitzW0_antitone_beta
    (alpha beta₁ beta₂ : ℝ)
    (halpha : 0 < alpha) (hbeta : 1 ≤ beta₁)
    (horder : beta₁ ≤ beta₂) :
    horwitzW0 alpha beta₂ ≤ horwitzW0 alpha beta₁ := by
  have hdenom₁ : 0 < horwitzDenominator alpha beta₁ := by
    simp only [horwitzDenominator]
    linarith
  have hdenom₂ : 0 < horwitzDenominator alpha beta₂ := by
    simp only [horwitzDenominator]
    linarith
  unfold horwitzW0
  rw [div_le_div_iff₀ hdenom₂ hdenom₁]
  simp only [horwitzDenominator]
  nlinarith

theorem horwitzW1_antitone_alpha
    (alpha₁ alpha₂ beta : ℝ)
    (halpha : 1 ≤ alpha₁) (hbeta : 0 < beta)
    (horder : alpha₁ ≤ alpha₂) :
    horwitzW1 alpha₂ beta ≤ horwitzW1 alpha₁ beta := by
  have hdenom₁ : 0 < horwitzDenominator alpha₁ beta := by
    simp only [horwitzDenominator]
    linarith
  have hdenom₂ : 0 < horwitzDenominator alpha₂ beta := by
    simp only [horwitzDenominator]
    linarith
  unfold horwitzW1
  rw [div_le_div_iff₀ hdenom₂ hdenom₁]
  simp only [horwitzDenominator]
  nlinarith

theorem horwitzW1_mono_beta
    (alpha beta₁ beta₂ : ℝ)
    (halpha : 1 ≤ alpha) (hbeta : 0 < beta₁)
    (horder : beta₁ ≤ beta₂) :
    horwitzW1 alpha beta₁ ≤ horwitzW1 alpha beta₂ := by
  have hdenom₁ : 0 < horwitzDenominator alpha beta₁ := by
    simp only [horwitzDenominator]
    linarith
  have hdenom₂ : 0 < horwitzDenominator alpha beta₂ := by
    simp only [horwitzDenominator]
    linarith
  unfold horwitzW1
  rw [div_le_div_iff₀ hdenom₁ hdenom₂]
  simp only [horwitzDenominator]
  nlinarith

theorem horwitzW2_mono_denominator
    (alpha₁ beta₁ alpha₂ beta₂ : ℝ)
    (hpositive : 0 < horwitzDenominator alpha₁ beta₁)
    (horder : horwitzDenominator alpha₁ beta₁ ≤
      horwitzDenominator alpha₂ beta₂) :
    horwitzW2 alpha₁ beta₁ ≤ horwitzW2 alpha₂ beta₂ := by
  have hpositive₂ : 0 < horwitzDenominator alpha₂ beta₂ :=
    hpositive.trans_le horder
  have hreciprocal :
      1 / horwitzDenominator alpha₂ beta₂ ≤
        1 / horwitzDenominator alpha₁ beta₁ := by
    rw [div_le_div_iff₀ hpositive₂ hpositive]
    linarith
  unfold horwitzW2
  linarith

structure HorwitzBoxBounds
    (alphaLower alphaUpper betaLower betaUpper alpha beta : ℝ) : Prop where
  w0Lower : horwitzW0 alphaLower betaUpper ≤ horwitzW0 alpha beta
  w0Upper : horwitzW0 alpha beta ≤ horwitzW0 alphaUpper betaLower
  w1Lower : horwitzW1 alphaUpper betaLower ≤ horwitzW1 alpha beta
  w1Upper : horwitzW1 alpha beta ≤ horwitzW1 alphaLower betaUpper
  w2Lower : horwitzW2 alphaLower betaLower ≤ horwitzW2 alpha beta
  w2Upper : horwitzW2 alpha beta ≤ horwitzW2 alphaUpper betaUpper

theorem horwitz_box_bounds
    (alphaLower alphaUpper betaLower betaUpper alpha beta : ℝ)
    (halphaOne : 1 ≤ alphaLower) (hbetaOne : 1 ≤ betaLower)
    (halphaBox : alphaLower ≤ alpha ∧ alpha ≤ alphaUpper)
    (hbetaBox : betaLower ≤ beta ∧ beta ≤ betaUpper) :
    HorwitzBoxBounds
      alphaLower alphaUpper betaLower betaUpper alpha beta := by
  have halphaPositive : 0 < alphaLower := lt_of_lt_of_le zero_lt_one halphaOne
  have hbetaPositive : 0 < betaLower := lt_of_lt_of_le zero_lt_one hbetaOne
  have halphaCurrent : 1 ≤ alpha := halphaOne.trans halphaBox.1
  have hbetaCurrent : 1 ≤ beta := hbetaOne.trans hbetaBox.1
  have halphaUpper : 1 ≤ alphaUpper := halphaCurrent.trans halphaBox.2
  have hbetaUpper : 1 ≤ betaUpper := hbetaCurrent.trans hbetaBox.2
  have halphaCurrentPositive : 0 < alpha :=
    lt_of_lt_of_le zero_lt_one halphaCurrent
  have hbetaCurrentPositive : 0 < beta :=
    lt_of_lt_of_le zero_lt_one hbetaCurrent
  have halphaUpperPositive : 0 < alphaUpper :=
    lt_of_lt_of_le zero_lt_one halphaUpper
  refine {
    w0Lower := ?_
    w0Upper := ?_
    w1Lower := ?_
    w1Upper := ?_
    w2Lower := ?_
    w2Upper := ?_
  }
  · exact (horwitzW0_mono_alpha alphaLower alpha betaUpper
      halphaPositive hbetaUpper halphaBox.1).trans
        (horwitzW0_antitone_beta alpha beta betaUpper
          halphaCurrentPositive hbetaCurrent hbetaBox.2)
  · exact (horwitzW0_mono_alpha alpha alphaUpper beta
      halphaCurrentPositive hbetaCurrent halphaBox.2).trans
        (horwitzW0_antitone_beta alphaUpper betaLower beta
          halphaUpperPositive hbetaOne hbetaBox.1)
  · exact (horwitzW1_antitone_alpha alpha alphaUpper betaLower
      halphaCurrent hbetaPositive halphaBox.2).trans
        (horwitzW1_mono_beta alpha betaLower beta
          halphaCurrent hbetaPositive hbetaBox.1)
  · exact (horwitzW1_mono_beta alpha beta betaUpper
      halphaCurrent hbetaCurrentPositive hbetaBox.2).trans
        (horwitzW1_antitone_alpha alphaLower alpha betaUpper
          halphaOne (lt_of_lt_of_le zero_lt_one hbetaUpper) halphaBox.1)
  · apply horwitzW2_mono_denominator alphaLower betaLower alpha beta
    · simp only [horwitzDenominator]
      linarith
    · simp only [horwitzDenominator]
      linarith
  · apply horwitzW2_mono_denominator alpha beta alphaUpper betaUpper
    · simp only [horwitzDenominator]
      linarith
    · simp only [horwitzDenominator]
      linarith

structure HorwitzStoredEnclosure
    (alpha beta w0Lower w0Upper w1Lower w1Upper w2Lower w2Upper : ℝ) : Prop where
  w0 : w0Lower ≤ horwitzW0 alpha beta ∧ horwitzW0 alpha beta ≤ w0Upper
  w1 : w1Lower ≤ horwitzW1 alpha beta ∧ horwitzW1 alpha beta ≤ w1Upper
  w2 : w2Lower ≤ horwitzW2 alpha beta ∧ horwitzW2 alpha beta ≤ w2Upper

def source15818AlphaLower : ℚ := 985 / 512
def source15818AlphaUpper : ℚ := 493 / 256
def source15818BetaLower : ℚ := 2579148085650963 / 2251799813685248
def source15818BetaUpper : ℚ := 647125223355353 / 562949953421312

def source15818W0Lower : ℚ := 261175671186813 / 281474976710656
def source15818W0Upper : ℚ := 8374994374731957 / 9007199254740992
def source15818W1Lower : ℚ := 77829474067439 / 140737488355328
def source15818W1Upper : ℚ := 4993843966331245 / 9007199254740992
def source15818W2Lower : ℚ := 4654212870180819 / 9007199254740992
def source15818W2Upper : ℚ := 4667021567102341 / 9007199254740992

theorem source15818_horwitz_outer_enclosure
    (alpha beta : ℝ)
    (halpha : (source15818AlphaLower : ℝ) ≤ alpha ∧
      alpha ≤ (source15818AlphaUpper : ℝ))
    (hbeta : (source15818BetaLower : ℝ) ≤ beta ∧
      beta ≤ (source15818BetaUpper : ℝ)) :
    HorwitzStoredEnclosure alpha beta
      source15818W0Lower source15818W0Upper
      source15818W1Lower source15818W1Upper
      source15818W2Lower source15818W2Upper := by
  have hbox := horwitz_box_bounds
    (source15818AlphaLower : ℝ) (source15818AlphaUpper : ℝ)
    (source15818BetaLower : ℝ) (source15818BetaUpper : ℝ)
    alpha beta
    (by norm_num [source15818AlphaLower])
    (by norm_num [source15818BetaLower]) halpha hbeta
  refine ⟨⟨?_, ?_⟩, ⟨?_, ?_⟩, ⟨?_, ?_⟩⟩
  · exact (by
      norm_num [source15818W0Lower, source15818AlphaLower,
        source15818BetaUpper, horwitzW0, horwitzDenominator] :
        (source15818W0Lower : ℝ) ≤
          horwitzW0 source15818AlphaLower source15818BetaUpper).trans
        hbox.w0Lower
  · exact hbox.w0Upper.trans (by
      norm_num [source15818W0Upper, source15818AlphaUpper,
        source15818BetaLower, horwitzW0, horwitzDenominator] :
        horwitzW0 source15818AlphaUpper source15818BetaLower ≤
          (source15818W0Upper : ℝ))
  · exact (by
      norm_num [source15818W1Lower, source15818AlphaUpper,
        source15818BetaLower, horwitzW1, horwitzDenominator] :
        (source15818W1Lower : ℝ) ≤
          horwitzW1 source15818AlphaUpper source15818BetaLower).trans
        hbox.w1Lower
  · exact hbox.w1Upper.trans (by
      norm_num [source15818W1Upper, source15818AlphaLower,
        source15818BetaUpper, horwitzW1, horwitzDenominator] :
        horwitzW1 source15818AlphaLower source15818BetaUpper ≤
          (source15818W1Upper : ℝ))
  · exact (by
      norm_num [source15818W2Lower, source15818AlphaLower,
        source15818BetaLower, horwitzW2, horwitzDenominator] :
        (source15818W2Lower : ℝ) ≤
          horwitzW2 source15818AlphaLower source15818BetaLower).trans
        hbox.w2Lower
  · exact hbox.w2Upper.trans (by
      norm_num [source15818W2Upper, source15818AlphaUpper,
        source15818BetaUpper, horwitzW2, horwitzDenominator] :
        horwitzW2 source15818AlphaUpper source15818BetaUpper ≤
          (source15818W2Upper : ℝ))

end

end CarmenQExact

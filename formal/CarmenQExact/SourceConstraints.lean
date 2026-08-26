import CarmenQExact.Horwitz

namespace CarmenQExact

noncomputable section

/-!
Exact, solver-independent lemmas used to embed a physical point in the
source SOCP.  These statements deliberately stop before any numerical cover
claim: a covered projective support line must enter the final theorem as an
explicit premise until its own exact certificate has been checked.
-/

structure OrderedSimplex4 (p₀ p₁ p₂ p₃ : ℝ) : Prop where
  nonneg₀ : 0 ≤ p₀
  nonneg₁ : 0 ≤ p₁
  nonneg₂ : 0 ≤ p₂
  nonneg₃ : 0 ≤ p₃
  total : p₀ + p₁ + p₂ + p₃ = 1
  order01 : p₁ ≤ p₀
  order12 : p₂ ≤ p₁
  order23 : p₃ ≤ p₂

theorem OrderedSimplex4.rank_bounds {p₀ p₁ p₂ p₃ : ℝ}
    (h : OrderedSimplex4 p₀ p₁ p₂ p₃) :
    1 / 4 ≤ p₀ ∧ p₁ ≤ 1 / 2 ∧ p₂ ≤ 1 / 3 ∧ p₃ ≤ 1 / 4 := by
  constructor
  · linarith [h.total, h.order01, h.order12, h.order23]
  constructor
  · linarith [h.total, h.nonneg₂, h.nonneg₃, h.order01]
  constructor
  · linarith [h.total, h.nonneg₃, h.order01, h.order12]
  · linarith [h.total, h.order01, h.order12, h.order23]

/-! The four McCormick inequalities used for `q = w u`. -/

theorem exact_product_in_mccormick_box
    (w u wLower wUpper uUpper : ℝ)
    (hwLower : wLower ≤ w) (hwUpper : w ≤ wUpper)
    (huLower : 0 ≤ u) (huUpper : u ≤ uUpper) :
    wLower * u ≤ w * u ∧
      w * u ≤ wUpper * u ∧
      wUpper * u + uUpper * w - wUpper * uUpper ≤ w * u ∧
      w * u ≤ wLower * u + uUpper * w - wLower * uUpper := by
  have hwuLower : 0 ≤ (w - wLower) * u :=
    mul_nonneg (sub_nonneg.mpr hwLower) huLower
  have hwuUpper : 0 ≤ (wUpper - w) * u :=
    mul_nonneg (sub_nonneg.mpr hwUpper) huLower
  have hcornerUpper : 0 ≤ (wUpper - w) * (uUpper - u) :=
    mul_nonneg (sub_nonneg.mpr hwUpper) (sub_nonneg.mpr huUpper)
  have hcornerLower : 0 ≤ (w - wLower) * (uUpper - u) :=
    mul_nonneg (sub_nonneg.mpr hwLower) (sub_nonneg.mpr huUpper)
  refine ⟨?_, ?_, ?_, ?_⟩ <;> nlinarith

/-! Concave McCormick upper envelopes used for the comparison loss. -/

theorem exact_product_below_concave_envelope
    (x p lower upper : ℝ)
    (hxLower : lower ≤ x) (hxUpper : x ≤ upper)
    (hpLower : 0 ≤ p) (hpUpper : p ≤ 1) :
    x * p ≤ upper * p ∧ x * p ≤ lower * p + x - lower := by
  have hfirst : 0 ≤ (upper - x) * p :=
    mul_nonneg (sub_nonneg.mpr hxUpper) hpLower
  have hsecond : 0 ≤ (x - lower) * (1 - p) :=
    mul_nonneg (sub_nonneg.mpr hxLower) (sub_nonneg.mpr hpUpper)
  constructor <;> nlinarith

/-! The rotated two-dimensional SOC used by the Hellinger hypograph. -/

def Lorentz3 (t x y : ℝ) : Prop := 0 ≤ t ∧ x ^ 2 + y ^ 2 ≤ t ^ 2

theorem hellinger_soc_identity (p q : ℝ) (hp : 0 ≤ p) (hq : 0 ≤ q) :
    Lorentz3 (p + q) (2 * √(p * q)) (p - q) := by
  have hpq : 0 ≤ p * q := mul_nonneg hp hq
  have hsqrt : (√(p * q)) ^ 2 = p * q := Real.sq_sqrt hpq
  constructor
  · positivity
  · nlinarith [hsqrt]

theorem hellinger_auxiliary_is_exact
    (p q g : ℝ) (hp : 0 ≤ p) (hq : 0 ≤ q)
    (hg : g = √(p * q)) : Lorentz3 (p + q) (2 * g) (p - q) := by
  subst g
  exact hellinger_soc_identity p q hp hq

/-! Exact identities behind the Horwitz weight constraints. -/

theorem horwitz_weight_sum (alpha beta : ℝ)
    (hdenom : horwitzDenominator alpha beta ≠ 0) :
    horwitzW0 alpha beta + horwitzW1 alpha beta +
      horwitzW2 alpha beta = 2 := by
  unfold horwitzW0 horwitzW1 horwitzW2
  simp only [horwitzDenominator] at hdenom ⊢
  field_simp [hdenom]
  ring

theorem horwitz_w0_residual (alpha beta : ℝ)
    (hdenom : horwitzDenominator alpha beta ≠ 0) :
    horwitzW0 alpha beta = alpha * (1 - horwitzW2 alpha beta) := by
  unfold horwitzW0 horwitzW2
  simp only [horwitzDenominator] at hdenom ⊢
  field_simp [hdenom]
  ring

theorem horwitz_w1_residual (alpha beta : ℝ)
    (hdenom : horwitzDenominator alpha beta ≠ 0) :
    horwitzW1 alpha beta = beta * (1 - horwitzW2 alpha beta) := by
  unfold horwitzW1 horwitzW2
  simp only [horwitzDenominator] at hdenom ⊢
  field_simp [hdenom]
  ring

theorem horwitz_weight_order (alpha beta : ℝ)
    (hdenom : 0 < horwitzDenominator alpha beta)
    (hbetaAlpha : beta ≤ alpha) (halphaTwo : alpha ≤ 2) :
    horwitzW2 alpha beta ≤ horwitzW1 alpha beta ∧
      horwitzW1 alpha beta ≤ horwitzW0 alpha beta := by
  have hdenomNe : horwitzDenominator alpha beta ≠ 0 := ne_of_gt hdenom
  constructor
  · unfold horwitzW1 horwitzW2
    rw [le_div_iff₀ hdenom]
    simp only [horwitzDenominator] at hdenom hdenomNe ⊢
    field_simp [hdenomNe]
    linarith
  · unfold horwitzW0 horwitzW1
    exact (div_le_div_iff_of_pos_right hdenom).mpr hbetaAlpha

theorem source15818_weight_floor
    (alpha beta : ℝ)
    (halpha : (source15818AlphaLower : ℝ) ≤ alpha ∧
      alpha ≤ (source15818AlphaUpper : ℝ))
    (hbeta : (source15818BetaLower : ℝ) ≤ beta ∧
      beta ≤ (source15818BetaUpper : ℝ)) :
    (79 / 100 : ℝ) ≤ horwitzW0 alpha beta := by
  have h := source15818_horwitz_outer_enclosure alpha beta halpha hbeta
  exact (by norm_num [source15818W0Lower] :
    (79 / 100 : ℝ) ≤ source15818W0Lower).trans h.w0.1

/-! Polynomial form of a Horwitz inellipse and two generic outer-relaxation
principles.  These isolate exactly what an interval generator must certify. -/

def horwitzCross (alpha beta : ℝ) : ℝ :=
  -2 * (alpha * beta - 2 * alpha - 2 * beta + 2)

def horwitzInellipsePolynomial
    (alpha beta x y scale : ℝ) : ℝ :=
  beta ^ 2 * x ^ 2 + alpha ^ 2 * y ^ 2 +
    horwitzCross alpha beta * x * y -
    2 * beta * x * scale - 2 * alpha * y * scale + scale ^ 2

def horwitzAlphaDerivative
    (alpha beta x y scale : ℝ) : ℝ :=
  2 * alpha * y ^ 2 + (4 - 2 * beta) * x * y - 2 * y * scale

def horwitzBetaDerivative
    (alpha beta x y scale : ℝ) : ℝ :=
  2 * beta * x ^ 2 + (4 - 2 * alpha) * x * y - 2 * x * scale

theorem horwitz_tangent_remainder
    (alpha beta deltaAlpha deltaBeta x y scale : ℝ) :
    horwitzInellipsePolynomial
        (alpha + deltaAlpha) (beta + deltaBeta) x y scale -
      horwitzInellipsePolynomial alpha beta x y scale -
      horwitzAlphaDerivative alpha beta x y scale * deltaAlpha -
      horwitzBetaDerivative alpha beta x y scale * deltaBeta =
        (deltaBeta * x - deltaAlpha * y) ^ 2 := by
  simp only [horwitzInellipsePolynomial, horwitzCross,
    horwitzAlphaDerivative, horwitzBetaDerivative]
  ring

theorem horwitz_tangent_lower_bound
    (alpha beta deltaAlpha deltaBeta x y scale : ℝ) :
    horwitzInellipsePolynomial alpha beta x y scale +
        horwitzAlphaDerivative alpha beta x y scale * deltaAlpha +
        horwitzBetaDerivative alpha beta x y scale * deltaBeta ≤
      horwitzInellipsePolynomial
        (alpha + deltaAlpha) (beta + deltaBeta) x y scale := by
  have hsquare : 0 ≤ (deltaBeta * x - deltaAlpha * y) ^ 2 := sq_nonneg _
  nlinarith [horwitz_tangent_remainder
    alpha beta deltaAlpha deltaBeta x y scale]

theorem horwitz_tangent_outer_relaxation
    (alpha beta deltaAlpha deltaBeta x y scale error : ℝ)
    (htrue : horwitzInellipsePolynomial
      (alpha + deltaAlpha) (beta + deltaBeta) x y scale ≤ 0)
    (herror :
      -(horwitzAlphaDerivative alpha beta x y scale * deltaAlpha +
        horwitzBetaDerivative alpha beta x y scale * deltaBeta) ≤
          error * scale ^ 2) :
    horwitzInellipsePolynomial alpha beta x y scale -
        error * scale ^ 2 ≤ 0 := by
  have htangent := horwitz_tangent_lower_bound
    alpha beta deltaAlpha deltaBeta x y scale
  linarith

theorem anchor_error_outer_relaxation
    (qTrue qAnchor error scale : ℝ)
    (htrue : qTrue ≤ 0)
    (herror : qAnchor - qTrue ≤ error * scale ^ 2) :
    qAnchor - error * scale ^ 2 ≤ 0 := by
  linarith

theorem coefficientwise_quadratic_outer_relaxation
    (a b c d e f aLower bLower cLower dLower eLower fLower
      x y scale : ℝ)
    (hx : 0 ≤ x) (hy : 0 ≤ y) (hscale : 0 ≤ scale)
    (ha : aLower ≤ a) (hb : bLower ≤ b) (hc : cLower ≤ c)
    (hd : dLower ≤ d) (he : eLower ≤ e) (hf : fLower ≤ f)
    (htrue : a * x ^ 2 + b * y ^ 2 + c * x * y +
      d * x * scale + e * y * scale + f * scale ^ 2 ≤ 0) :
    aLower * x ^ 2 + bLower * y ^ 2 + cLower * x * y +
      dLower * x * scale + eLower * y * scale +
      fLower * scale ^ 2 ≤ 0 := by
  have hxSq : 0 ≤ x ^ 2 := sq_nonneg x
  have hySq : 0 ≤ y ^ 2 := sq_nonneg y
  have hsSq : 0 ≤ scale ^ 2 := sq_nonneg scale
  have hxy : 0 ≤ x * y := mul_nonneg hx hy
  have hxs : 0 ≤ x * scale := mul_nonneg hx hscale
  have hys : 0 ≤ y * scale := mul_nonneg hy hscale
  have hA : 0 ≤ (a - aLower) * x ^ 2 :=
    mul_nonneg (sub_nonneg.mpr ha) hxSq
  have hB : 0 ≤ (b - bLower) * y ^ 2 :=
    mul_nonneg (sub_nonneg.mpr hb) hySq
  have hC : 0 ≤ (c - cLower) * (x * y) :=
    mul_nonneg (sub_nonneg.mpr hc) hxy
  have hD : 0 ≤ (d - dLower) * (x * scale) :=
    mul_nonneg (sub_nonneg.mpr hd) hxs
  have hE : 0 ≤ (e - eLower) * (y * scale) :=
    mul_nonneg (sub_nonneg.mpr he) hys
  have hF : 0 ≤ (f - fLower) * scale ^ 2 :=
    mul_nonneg (sub_nonneg.mpr hf) hsSq
  nlinarith

def AffineSoc2
    (radius scale u₀ u₁ : ℝ) : Prop :=
  0 ≤ radius * scale ∧ u₀ ^ 2 + u₁ ^ 2 ≤ (radius * scale) ^ 2

theorem affine_soc_implies_completed_quadratic
    (radius scale u₀ u₁ : ℝ)
    (h : AffineSoc2 radius scale u₀ u₁) :
    u₀ ^ 2 + u₁ ^ 2 - (radius * scale) ^ 2 ≤ 0 := by
  exact sub_nonpos.mpr h.2

end

end CarmenQExact

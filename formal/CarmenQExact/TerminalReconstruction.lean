import CarmenQExact.Horwitz

namespace CarmenQExact

noncomputable section

/-! Exact enclosure of the terminal reconstruction used by source cell 15818.

Every decimal below denotes the exact rational value of the binary64 number
stored in the source artifact.  The proof reasons over arbitrary real
`alpha,beta` in the source box.  It therefore does not trust the Python
`nextafter` implementation: the endpoint comparisons are discharged again in
Lean over exact rationals.
-/

structure InInterval (x lower upper : ℝ) : Prop where
  lower_le : lower ≤ x
  le_upper : x ≤ upper

theorem InInterval.weaken {x lower upper lower' upper' : ℝ}
    (h : InInterval x lower upper) (hlower : lower' ≤ lower)
    (hupper : upper ≤ upper') : InInterval x lower' upper' :=
  ⟨hlower.trans h.lower_le, h.le_upper.trans hupper⟩

theorem InInterval.point (x : ℝ) : InInterval x x x := ⟨le_rfl, le_rfl⟩

theorem InInterval.add {x y xl xu yl yu : ℝ}
    (hx : InInterval x xl xu) (hy : InInterval y yl yu) :
    InInterval (x + y) (xl + yl) (xu + yu) :=
  ⟨add_le_add hx.lower_le hy.lower_le, add_le_add hx.le_upper hy.le_upper⟩

theorem InInterval.neg {x xl xu : ℝ} (hx : InInterval x xl xu) :
    InInterval (-x) (-xu) (-xl) :=
  ⟨neg_le_neg hx.le_upper, neg_le_neg hx.lower_le⟩

theorem InInterval.sub {x y xl xu yl yu : ℝ}
    (hx : InInterval x xl xu) (hy : InInterval y yl yu) :
    InInterval (x - y) (xl - yu) (xu - yl) := by
  simpa [sub_eq_add_neg] using hx.add hy.neg

theorem InInterval.smul_nonneg {x xl xu r : ℝ}
    (hx : InInterval x xl xu) (hr : 0 ≤ r) :
    InInterval (r * x) (r * xl) (r * xu) :=
  ⟨mul_le_mul_of_nonneg_left hx.lower_le hr,
    mul_le_mul_of_nonneg_left hx.le_upper hr⟩

theorem InInterval.mul_nonneg {x y xl xu yl yu : ℝ}
    (hx : InInterval x xl xu) (hy : InInterval y yl yu)
    (hxl : 0 ≤ xl) (hyl : 0 ≤ yl) :
    InInterval (x * y) (xl * yl) (xu * yu) := by
  have hx0 : 0 ≤ x := hxl.trans hx.lower_le
  have hy0 : 0 ≤ y := hyl.trans hy.lower_le
  exact ⟨mul_le_mul hx.lower_le hy.lower_le hyl hx0,
    mul_le_mul hx.le_upper hy.le_upper hy0 (hx0.trans hx.le_upper)⟩

theorem InInterval.mul_nonpos_nonneg {x y xl xu yl yu : ℝ}
    (hx : InInterval x xl xu) (hy : InInterval y yl yu)
    (hxu : xu ≤ 0) (hyl : 0 ≤ yl) :
    InInterval (x * y) (xl * yu) (xu * yl) := by
  have hy0 : 0 ≤ y := hyl.trans hy.lower_le
  have hxl : xl ≤ 0 := hx.lower_le.trans (hx.le_upper.trans hxu)
  constructor
  · calc
      xl * yu ≤ xl * y := mul_le_mul_of_nonpos_left hy.le_upper hxl
      _ ≤ x * y := mul_le_mul_of_nonneg_right hx.lower_le hy0
  · calc
      x * y ≤ xu * y := mul_le_mul_of_nonneg_right hx.le_upper hy0
      _ ≤ xu * yl := mul_le_mul_of_nonpos_left hy.lower_le hxu

theorem InInterval.reciprocal_pos {x xl xu : ℝ}
    (hx : InInterval x xl xu) (hxl : 0 < xl) :
    InInterval (1 / x) (1 / xu) (1 / xl) := by
  have hx0 : 0 < x := hxl.trans_le hx.lower_le
  exact ⟨one_div_le_one_div_of_le hx0 hx.le_upper,
    one_div_le_one_div_of_le hxl hx.lower_le⟩

theorem InInterval.div_nonneg {x y xl xu yl yu : ℝ}
    (hx : InInterval x xl xu) (hy : InInterval y yl yu)
    (hxl : 0 ≤ xl) (hyl : 0 < yl) :
    InInterval (x / y) (xl / yu) (xu / yl) := by
  have hyu : 0 < yu := hyl.trans_le (hy.lower_le.trans hy.le_upper)
  simpa [div_eq_mul_inv] using
    hx.mul_nonneg (hy.reciprocal_pos hyl) hxl
      (le_of_lt (one_div_pos.mpr hyu))

theorem InInterval.div_nonpos {x y xl xu yl yu : ℝ}
    (hx : InInterval x xl xu) (hy : InInterval y yl yu)
    (hxu : xu ≤ 0) (hyl : 0 < yl) :
    InInterval (x / y) (xl / yl) (xu / yu) := by
  have hyu : 0 < yu := hyl.trans_le (hy.lower_le.trans hy.le_upper)
  simpa [div_eq_mul_inv] using
    hx.mul_nonpos_nonneg (hy.reciprocal_pos hyl) hxu
      (le_of_lt (one_div_pos.mpr hyu))

theorem InInterval.sq_nonpos {x xl xu : ℝ}
    (hx : InInterval x xl xu) (hxu : xu ≤ 0) :
    InInterval (x ^ 2) (xu ^ 2) (xl ^ 2) := by
  constructor
  · have hleft : 0 ≤ -xu := by linarith
    have hright : 0 ≤ -x := by linarith [hx.le_upper]
    have hsquare := (sq_le_sq₀ hleft hright).2 (by linarith [hx.le_upper])
    simpa using hsquare
  · have hleft : 0 ≤ -x := by linarith [hx.le_upper]
    have hright : 0 ≤ -xl := by linarith [hx.lower_le, hx.le_upper]
    have hsquare := (sq_le_sq₀ hleft hright).2 (by linarith [hx.lower_le])
    simpa using hsquare

theorem InInterval.sqrt {x xl xu sl su : ℝ}
    (hx : InInterval x xl xu) (hxl : 0 ≤ xl)
    (hsl : 0 ≤ sl) (hsu : 0 ≤ su)
    (hslSq : sl ^ 2 ≤ xl) (hsuSq : xu ≤ su ^ 2) :
    InInterval (√x) sl su := by
  have hx0 : 0 ≤ x := hxl.trans hx.lower_le
  constructor
  · exact (Real.le_sqrt hsl hx0).2 (hslSq.trans hx.lower_le)
  · exact Real.sqrt_le_iff.mpr ⟨hsu, hx.le_upper.trans hsuSq⟩

def reconstructionCosine (alpha beta : ℝ) : ℝ :=
  1 - 2 / alpha - 2 / beta + 2 / (alpha * beta)

def reconstructionSine (alpha beta : ℝ) : ℝ :=
  √(1 - reconstructionCosine alpha beta ^ 2)

def reconstructionX0 (alpha beta : ℝ) : ℝ :=
  1 + 2 * (beta - 1) / alpha

def reconstructionY0 (alpha beta : ℝ) : ℝ :=
  (-1 - reconstructionCosine alpha beta * reconstructionX0 alpha beta) /
    reconstructionSine alpha beta

def reconstructionY1 (alpha beta : ℝ) : ℝ :=
  (1 + 2 * (alpha - 1) / beta + reconstructionCosine alpha beta) /
    reconstructionSine alpha beta

def reconstructionY2 (alpha beta : ℝ) : ℝ :=
  (-1 + reconstructionCosine alpha beta) / reconstructionSine alpha beta

theorem reconstructionX0_eq_horwitz (alpha beta : ℝ)
    (halpha : alpha ≠ 0) (_hdenom : horwitzDenominator alpha beta ≠ 0) :
    reconstructionX0 alpha beta = 2 / horwitzW0 alpha beta - 1 := by
  rw [reconstructionX0, horwitzW0, horwitzDenominator]
  field_simp
  ring

theorem reconstructionY1Numerator_eq_horwitz (alpha beta : ℝ)
    (hbeta : beta ≠ 0) (_hdenom : horwitzDenominator alpha beta ≠ 0) :
    1 + 2 * (alpha - 1) / beta + reconstructionCosine alpha beta =
      2 / horwitzW1 alpha beta - 1 + reconstructionCosine alpha beta := by
  rw [horwitzW1, horwitzDenominator]
  field_simp
  ring

def source15818CosineLower : ℚ := -1986765117239993 / 2251799813685248
def source15818CosineUpper : ℚ := -1960733567034699 / 2251799813685248
def source15818SineLower : ℚ := 8479119840478813 / 18014398509481984
def source15818SineUpper : ℚ := 8858446924887201 / 18014398509481984

def source15818X0Lower : ℚ := 1293598198356577 / 1125899906842624
def source15818X0Upper : ℚ := 5212808775141477 / 4503599627370496
def source15818X12Lower : ℚ := -4503599627370497 / 4503599627370496
def source15818X12Upper : ℚ := -9007199254740991 / 9007199254740992

def source15818Y0Lower : ℚ := 4069749381756277 / 4611686018427387904
def source15818Y0Upper : ℚ := 6503975756596965 / 144115188075855872
def source15818Y1Lower : ℚ := 3933081649912597 / 1125899906842624
def source15818Y1Upper : ℚ := 4193422069918103 / 1125899906842624
def source15818Y2Lower : ℚ := -9005085340283451 / 2251799813685248
def source15818Y2Upper : ℚ := -2141635427130777 / 562949953421312

structure ReconstructionCoefficientEnclosure (alpha beta : ℝ) : Prop where
  cosine : InInterval (reconstructionCosine alpha beta)
    source15818CosineLower source15818CosineUpper
  sine : InInterval (reconstructionSine alpha beta)
    source15818SineLower source15818SineUpper
  x0 : InInterval (reconstructionX0 alpha beta)
    source15818X0Lower source15818X0Upper
  x1 : InInterval (-1 : ℝ) source15818X12Lower source15818X12Upper
  x2 : InInterval (-1 : ℝ) source15818X12Lower source15818X12Upper
  y0 : InInterval (reconstructionY0 alpha beta)
    source15818Y0Lower source15818Y0Upper
  y1 : InInterval (reconstructionY1 alpha beta)
    source15818Y1Lower source15818Y1Upper
  y2 : InInterval (reconstructionY2 alpha beta)
    source15818Y2Lower source15818Y2Upper

theorem source15818_reconstruction_coefficient_enclosure
    (alpha beta : ℝ)
    (halpha : (source15818AlphaLower : ℝ) ≤ alpha ∧
      alpha ≤ (source15818AlphaUpper : ℝ))
    (hbeta : (source15818BetaLower : ℝ) ≤ beta ∧
      beta ≤ (source15818BetaUpper : ℝ)) :
    ReconstructionCoefficientEnclosure alpha beta := by
  let ha : InInterval alpha source15818AlphaLower source15818AlphaUpper :=
    ⟨halpha.1, halpha.2⟩
  let hb : InInterval beta source15818BetaLower source15818BetaUpper :=
    ⟨hbeta.1, hbeta.2⟩
  have haPos : (0 : ℝ) < source15818AlphaLower := by
    norm_num [source15818AlphaLower]
  have hbPos : (0 : ℝ) < source15818BetaLower := by
    norm_num [source15818BetaLower]
  have haUpperPos : (0 : ℝ) < source15818AlphaUpper := by
    norm_num [source15818AlphaUpper]
  have hbUpperPos : (0 : ℝ) < source15818BetaUpper := by
    norm_num [source15818BetaUpper]
  have haRec := ha.reciprocal_pos haPos
  have hbRec := hb.reciprocal_pos hbPos
  have habRec := haRec.mul_nonneg hbRec
    (le_of_lt (one_div_pos.mpr haUpperPos))
    (le_of_lt (one_div_pos.mpr hbUpperPos))
  have htwoA := haRec.smul_nonneg (r := (2 : ℝ)) (by norm_num)
  have htwoB := hbRec.smul_nonneg (r := (2 : ℝ)) (by norm_num)
  have htwoAB := habRec.smul_nonneg (r := (2 : ℝ)) (by norm_num)
  have hcosRaw := ((InInterval.point (1 : ℝ)).sub htwoA).sub htwoB |>.add htwoAB
  have hcos : InInterval (reconstructionCosine alpha beta)
      source15818CosineLower source15818CosineUpper := by
    apply InInterval.weaken
    · simpa [reconstructionCosine, div_eq_mul_inv, mul_comm, mul_left_comm,
        mul_assoc] using hcosRaw
    · norm_num [source15818AlphaLower, source15818AlphaUpper,
        source15818BetaLower, source15818BetaUpper, source15818CosineLower]
    · norm_num [source15818AlphaLower, source15818AlphaUpper,
        source15818BetaLower, source15818BetaUpper, source15818CosineUpper]
  have hcosSq := hcos.sq_nonpos (by norm_num [source15818CosineUpper])
  have hsineSq := (InInterval.point (1 : ℝ)).sub hcosSq
  have hsine : InInterval (reconstructionSine alpha beta)
      source15818SineLower source15818SineUpper := by
    rw [reconstructionSine]
    apply hsineSq.sqrt
    · norm_num [source15818CosineLower]
    · norm_num [source15818SineLower]
    · norm_num [source15818SineUpper]
    · norm_num [source15818SineLower, source15818CosineLower]
    · norm_num [source15818SineUpper, source15818CosineUpper]
  have hbetaMinus := hb.sub (InInterval.point (1 : ℝ))
  have hbetaFrac := hbetaMinus.div_nonneg ha
    (by norm_num [source15818BetaLower]) haPos
  have hxRaw := (InInterval.point (1 : ℝ)).add
    (hbetaFrac.smul_nonneg (r := (2 : ℝ)) (by norm_num))
  have hx : InInterval (reconstructionX0 alpha beta)
      source15818X0Lower source15818X0Upper := by
    apply InInterval.weaken
    · simpa [reconstructionX0, div_eq_mul_inv, mul_assoc] using hxRaw
    · norm_num [source15818AlphaLower, source15818AlphaUpper,
        source15818BetaLower, source15818BetaUpper, source15818X0Lower]
    · norm_num [source15818AlphaLower, source15818AlphaUpper,
        source15818BetaLower, source15818BetaUpper, source15818X0Upper]
  have hcosX := hcos.mul_nonpos_nonneg hx
    (by norm_num [source15818CosineUpper])
    (by norm_num [source15818X0Lower])
  have hnum0 := (InInterval.point (-1 : ℝ)).sub hcosX
  have hy0Raw := hnum0.div_nonneg hsine
    (by norm_num [source15818CosineUpper, source15818X0Lower])
    (by norm_num [source15818SineLower])
  have hy0 : InInterval (reconstructionY0 alpha beta)
      source15818Y0Lower source15818Y0Upper := by
    rw [reconstructionY0]
    apply hy0Raw.weaken
    · norm_num [source15818CosineLower, source15818CosineUpper,
        source15818SineLower, source15818SineUpper,
        source15818X0Lower, source15818X0Upper, source15818Y0Lower]
    · norm_num [source15818CosineLower, source15818CosineUpper,
        source15818SineLower, source15818SineUpper,
        source15818X0Lower, source15818X0Upper, source15818Y0Upper]
  have halphaMinus := ha.sub (InInterval.point (1 : ℝ))
  have halphaFrac := halphaMinus.div_nonneg hb
    (by norm_num [source15818AlphaLower]) hbPos
  have hnum1 := ((InInterval.point (1 : ℝ)).add
    (halphaFrac.smul_nonneg (r := (2 : ℝ)) (by norm_num))).add hcos
  have hy1Raw := hnum1.div_nonneg hsine
    (by norm_num [source15818AlphaLower, source15818AlphaUpper,
        source15818BetaLower, source15818BetaUpper, source15818CosineLower])
    (by norm_num [source15818SineLower])
  have hy1 : InInterval (reconstructionY1 alpha beta)
      source15818Y1Lower source15818Y1Upper := by
    apply InInterval.weaken
    · simpa [reconstructionY1, div_eq_mul_inv, mul_assoc] using hy1Raw
    · norm_num [source15818AlphaLower, source15818AlphaUpper,
        source15818BetaLower, source15818BetaUpper,
        source15818CosineLower, source15818CosineUpper,
        source15818SineLower, source15818SineUpper, source15818Y1Lower]
    · norm_num [source15818AlphaLower, source15818AlphaUpper,
        source15818BetaLower, source15818BetaUpper,
        source15818CosineLower, source15818CosineUpper,
        source15818SineLower, source15818SineUpper, source15818Y1Upper]
  have hnum2 := (InInterval.point (-1 : ℝ)).add hcos
  have hy2Raw := hnum2.div_nonpos hsine
    (by norm_num [source15818CosineUpper])
    (by norm_num [source15818SineLower])
  have hy2 : InInterval (reconstructionY2 alpha beta)
      source15818Y2Lower source15818Y2Upper := by
    rw [reconstructionY2]
    apply hy2Raw.weaken
    · norm_num [source15818CosineLower, source15818CosineUpper,
        source15818SineLower, source15818SineUpper, source15818Y2Lower]
    · norm_num [source15818CosineLower, source15818CosineUpper,
        source15818SineLower, source15818SineUpper, source15818Y2Upper]
  exact {
    cosine := hcos
    sine := hsine
    x0 := hx
    x1 := by constructor <;> norm_num [source15818X12Lower, source15818X12Upper]
    x2 := by constructor <;> norm_num [source15818X12Lower, source15818X12Upper]
    y0 := hy0
    y1 := hy1
    y2 := hy2
  }

abbrev Reconstruction2 := EuclideanSpace ℝ (Fin 2)

theorem reconstruction2_sub (a b c d : ℝ) :
    (!₂[a, b] : Reconstruction2) - !₂[c, d] = !₂[a - c, b - d] := by
  ext i
  fin_cases i <;> rfl

def reconstructionColumn0 (alpha beta : ℝ) : Reconstruction2 :=
  !₂[reconstructionX0 alpha beta, reconstructionY0 alpha beta]

def reconstructionColumn1 (alpha beta : ℝ) : Reconstruction2 :=
  !₂[-1, reconstructionY1 alpha beta]

def reconstructionColumn2 (alpha beta : ℝ) : Reconstruction2 :=
  !₂[-1, reconstructionY2 alpha beta]

def source15818AnchorColumn0 : Reconstruction2 :=
  !₂[2596795519490791 / 2251799813685248,
      6472328410994315 / 288230376151711744]

def source15818AnchorColumn1 : Reconstruction2 :=
  !₂[-1, 8118128379027461 / 2251799813685248]

def source15818AnchorColumn2 : Reconstruction2 :=
  !₂[-1, -8778208511218709 / 2251799813685248]

def source15818ReconstructionError0 : ℚ := 103911665412485 / 4503599627370496
def source15818ReconstructionError1 : ℚ := 8598904345879841 / 72057594037927936
def source15818ReconstructionError2 : ℚ := 7260058530071745 / 72057594037927936

theorem abs_sub_le_of_interval {x lower upper anchor radius : ℝ}
    (hx : InInterval x lower upper)
    (hlower : anchor - radius ≤ lower) (hupper : upper ≤ anchor + radius) :
    |x - anchor| ≤ radius := by
  rw [abs_le]
  constructor <;> linarith [hx.lower_le, hx.le_upper]

theorem norm_pair_le {x y dx dy error : ℝ}
    (hx : |x| ≤ dx) (hy : |y| ≤ dy)
    (herror : 0 ≤ error) (hsquares : dx ^ 2 + dy ^ 2 ≤ error ^ 2) :
    ‖(!₂[x, y] : Reconstruction2)‖ ≤ error := by
  apply (sq_le_sq₀ (norm_nonneg _) herror).mp
  rw [EuclideanSpace.real_norm_sq_eq]
  have hdx : 0 ≤ dx := (abs_nonneg x).trans hx
  have hdy : 0 ≤ dy := (abs_nonneg y).trans hy
  have hxSq : x ^ 2 ≤ dx ^ 2 := by
    rw [sq_le_sq, abs_of_nonneg hdx]
    exact hx
  have hySq : y ^ 2 ≤ dy ^ 2 := by
    rw [sq_le_sq, abs_of_nonneg hdy]
    exact hy
  simp only [Fin.sum_univ_two, Matrix.cons_val_zero, Matrix.cons_val_one]
  linarith

structure ReconstructionColumnErrors (alpha beta : ℝ) : Prop where
  column0 : ‖reconstructionColumn0 alpha beta - source15818AnchorColumn0‖ ≤
    source15818ReconstructionError0
  column1 : ‖reconstructionColumn1 alpha beta - source15818AnchorColumn1‖ ≤
    source15818ReconstructionError1
  column2 : ‖reconstructionColumn2 alpha beta - source15818AnchorColumn2‖ ≤
    source15818ReconstructionError2

theorem source15818_reconstruction_column_errors
    (alpha beta : ℝ)
    (halpha : (source15818AlphaLower : ℝ) ≤ alpha ∧
      alpha ≤ (source15818AlphaUpper : ℝ))
    (hbeta : (source15818BetaLower : ℝ) ≤ beta ∧
      beta ≤ (source15818BetaUpper : ℝ)) :
    ReconstructionColumnErrors alpha beta := by
  have h := source15818_reconstruction_coefficient_enclosure alpha beta halpha hbeta
  have hx0 : |reconstructionX0 alpha beta -
      2596795519490791 / 2251799813685248| ≤
      19217736159895 / 4503599627370496 := by
    apply abs_sub_le_of_interval h.x0 <;>
      norm_num [source15818X0Lower, source15818X0Upper]
  have hy0 : |reconstructionY0 alpha beta -
      6472328410994315 / 288230376151711744| ≤
      6535623102199615 / 288230376151711744 := by
    apply abs_sub_le_of_interval h.y0 <;>
      norm_num [source15818Y0Lower, source15818Y0Upper]
  have hx1 : |(-1 : ℝ) - (-1)| ≤ 1 / 4503599627370496 := by norm_num
  have hy1 : |reconstructionY1 alpha beta -
      8118128379027461 / 2251799813685248| ≤
      268715760808745 / 2251799813685248 := by
    apply abs_sub_le_of_interval h.y1 <;>
      norm_num [source15818Y1Lower, source15818Y1Upper]
  have hx2 : |(-1 : ℝ) - (-1)| ≤ 1 / 4503599627370496 := by norm_num
  have hy2 : |reconstructionY2 alpha beta -
      (-8778208511218709 / 2251799813685248)| ≤
      113438414532371 / 1125899906842624 := by
    apply abs_sub_le_of_interval h.y2 <;>
      norm_num [source15818Y2Lower, source15818Y2Upper]
  constructor
  · have hpair : ‖(!₂[reconstructionX0 alpha beta -
        2596795519490791 / 2251799813685248,
        reconstructionY0 alpha beta -
          6472328410994315 / 288230376151711744] : Reconstruction2)‖ ≤
          (source15818ReconstructionError0 : ℝ) := by
      apply norm_pair_le hx0 hy0
      · norm_num [source15818ReconstructionError0]
      · norm_num [source15818ReconstructionError0]
    rw [reconstructionColumn0, source15818AnchorColumn0, reconstruction2_sub]
    exact hpair
  · have hpair : ‖(!₂[(-1 : ℝ) - (-1), reconstructionY1 alpha beta -
        8118128379027461 / 2251799813685248] : Reconstruction2)‖ ≤
          (source15818ReconstructionError1 : ℝ) := by
      apply norm_pair_le hx1 hy1
      · norm_num [source15818ReconstructionError1]
      · norm_num [source15818ReconstructionError1]
    rw [reconstructionColumn1, source15818AnchorColumn1, reconstruction2_sub]
    exact hpair
  · have hpair : ‖(!₂[(-1 : ℝ) - (-1), reconstructionY2 alpha beta -
        (-8778208511218709 / 2251799813685248)] : Reconstruction2)‖ ≤
          (source15818ReconstructionError2 : ℝ) := by
      apply norm_pair_le hx2 hy2
      · norm_num [source15818ReconstructionError2]
      · norm_num [source15818ReconstructionError2]
    rw [reconstructionColumn2, source15818AnchorColumn2, reconstruction2_sub]
    exact hpair

theorem three_column_error_budget
    (q₀ q₁ q₂ : ℝ) (d₀ d₁ d₂ : Reconstruction2)
    (e₀ e₁ e₂ : ℝ)
    (hd₀ : ‖d₀‖ ≤ e₀) (hd₁ : ‖d₁‖ ≤ e₁) (hd₂ : ‖d₂‖ ≤ e₂) :
    ‖q₀ • d₀ + q₁ • d₁ + q₂ • d₂‖ ≤
      |q₀| * e₀ + |q₁| * e₁ + |q₂| * e₂ := by
  calc
    ‖q₀ • d₀ + q₁ • d₁ + q₂ • d₂‖ ≤
        ‖q₀ • d₀ + q₁ • d₁‖ + ‖q₂ • d₂‖ := norm_add_le _ _
    _ ≤ (‖q₀ • d₀‖ + ‖q₁ • d₁‖) + ‖q₂ • d₂‖ :=
      add_le_add (norm_add_le _ _) le_rfl
    _ = (|q₀| * ‖d₀‖ + |q₁| * ‖d₁‖) + |q₂| * ‖d₂‖ := by
      simp only [norm_smul, Real.norm_eq_abs]
    _ ≤ (|q₀| * e₀ + |q₁| * e₁) + |q₂| * e₂ := by
      exact add_le_add
        (add_le_add
          (mul_le_mul_of_nonneg_left hd₀ (abs_nonneg q₀))
          (mul_le_mul_of_nonneg_left hd₁ (abs_nonneg q₁)))
        (mul_le_mul_of_nonneg_left hd₂ (abs_nonneg q₂))
    _ = |q₀| * e₀ + |q₁| * e₁ + |q₂| * e₂ := rfl

theorem source15818_reconstruction_error_budget
    (alpha beta q₀ q₁ q₂ : ℝ)
    (halpha : (source15818AlphaLower : ℝ) ≤ alpha ∧
      alpha ≤ (source15818AlphaUpper : ℝ))
    (hbeta : (source15818BetaLower : ℝ) ≤ beta ∧
      beta ≤ (source15818BetaUpper : ℝ)) :
    ‖q₀ • (reconstructionColumn0 alpha beta - source15818AnchorColumn0) +
        q₁ • (reconstructionColumn1 alpha beta - source15818AnchorColumn1) +
        q₂ • (reconstructionColumn2 alpha beta - source15818AnchorColumn2)‖ ≤
      |q₀| * source15818ReconstructionError0 +
        |q₁| * source15818ReconstructionError1 +
        |q₂| * source15818ReconstructionError2 := by
  have h := source15818_reconstruction_column_errors alpha beta halpha hbeta
  exact three_column_error_budget q₀ q₁ q₂ _ _ _ _ _ _
    h.column0 h.column1 h.column2

end

end CarmenQExact

import Mathlib

namespace CarmenQExact

noncomputable section

/-!
# Exact tangent majorant for a binary Hellinger term

The projective secular certificate removes its remaining pure-state variables
with the scalar inequality proved below.  `hellinger_amgm_scaled` is Young's
inequality applied to `t * sqrt q` and `sqrt (1 - q)`.  Dividing by positive
`t` yields the affine tangent majorant used by the interval branch-and-bound
checker.

The statements are entirely solver-independent.  In particular, the tangent
parameter need only be positive; a numerical routine may choose it for
tightness without participating in soundness.
-/

theorem hellinger_amgm_scaled
    (q t : ℝ) (hq0 : 0 ≤ q) (hq1 : q ≤ 1) :
    2 * t * √(q * (1 - q)) ≤ t ^ 2 * q + (1 - q) := by
  have h1q : 0 ≤ 1 - q := sub_nonneg.mpr hq1
  have hsq_q : (√q) ^ 2 = q := Real.sq_sqrt hq0
  have hsq_1q : (√(1 - q)) ^ 2 = 1 - q := Real.sq_sqrt h1q
  have hsqrt_mul : √(q * (1 - q)) = √q * √(1 - q) := by
    rw [Real.sqrt_mul hq0]
  have hsquare : 0 ≤ (t * √q - √(1 - q)) ^ 2 := sq_nonneg _
  rw [hsqrt_mul]
  nlinarith

theorem hellinger_tangent
    (q t : ℝ) (hq0 : 0 ≤ q) (hq1 : q ≤ 1) (ht : 0 < t) :
    1 + 2 * √(q * (1 - q)) ≤
      1 + 1 / t + (t - 1 / t) * q := by
  have hscaled := hellinger_amgm_scaled q t hq0 hq1
  have htne : t ≠ 0 := ne_of_gt ht
  have hdiv : 2 * √(q * (1 - q)) ≤
      (t ^ 2 * q + (1 - q)) / t := by
    apply (le_div_iff₀ ht).2
    nlinarith
  have hrhs : 1 / t + (t - 1 / t) * q =
      (t ^ 2 * q + (1 - q)) / t := by
    field_simp [htne]
    ring
  calc
    1 + 2 * √(q * (1 - q)) ≤
        1 + (t ^ 2 * q + (1 - q)) / t := by linarith
    _ = 1 + 1 / t + (t - 1 / t) * q := by rw [← hrhs]; ring

/-!
The low-eigenvalue interval kernel compares two generalized Rayleigh
quotients.  The following scalar lemma is the algebraic core of the sharper
perturbation estimate: using the already certified face quotient avoids the
extra inverse spectral-gap factor incurred by replacing it with `‖N₀‖ / d`.
-/

theorem rayleighPerturbation
    (n n0 d d0 deltaN deltaD lambda0 dmin : ℝ)
    (hdmin : 0 < dmin)
    (hdenominator : dmin ≤ d)
    (hnumerator : n ≤ n0 + deltaN)
    (hface : n0 ≤ lambda0 * d0)
    (hdenominatorChange : d0 ≤ d + deltaD)
    (hlambda : 0 ≤ lambda0)
    (hdeltaN : 0 ≤ deltaN)
    (hdeltaD : 0 ≤ deltaD) :
    n / d ≤ lambda0 + (deltaN + lambda0 * deltaD) / dmin := by
  have hd : 0 < d := lt_of_lt_of_le hdmin hdenominator
  have hface' : lambda0 * d0 ≤ lambda0 * (d + deltaD) :=
    mul_le_mul_of_nonneg_left hdenominatorChange hlambda
  have hnum : n ≤ lambda0 * d + (deltaN + lambda0 * deltaD) := by
    calc
      n ≤ n0 + deltaN := hnumerator
      _ ≤ lambda0 * d0 + deltaN := by
        simpa [add_comm] using add_le_add_right hface deltaN
      _ ≤ lambda0 * (d + deltaD) + deltaN :=
        by simpa [add_comm] using add_le_add_right hface' deltaN
      _ = lambda0 * d + (deltaN + lambda0 * deltaD) := by ring
  have hcorrection : 0 ≤ deltaN + lambda0 * deltaD := by positivity
  have hquotient :
      n / d ≤ lambda0 + (deltaN + lambda0 * deltaD) / d := by
    apply (div_le_iff₀ hd).2
    calc
      n ≤ lambda0 * d + (deltaN + lambda0 * deltaD) := hnum
      _ = (lambda0 + (deltaN + lambda0 * deltaD) / d) * d := by
        field_simp [ne_of_gt hd]
  have hgap :
      (deltaN + lambda0 * deltaD) / d ≤
        (deltaN + lambda0 * deltaD) / dmin := by
    exact div_le_div_of_nonneg_left hcorrection hdmin hdenominator
  linarith

end

end CarmenQExact

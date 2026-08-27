import Mathlib

namespace CarmenQExact

/-!
# Convex interpolation of support bounds

For a fixed strategy, the weighted AUDIT--RETURN score is affine in the
AUDIT weight.  Multiplying two valid endpoint inequalities by nonnegative
barycentric coefficients and adding them therefore supplies a valid bound at
every intermediate weight.  This file records that elementary bridge without
assuming that any particular endpoint bound is available.
-/

theorem supportInterpolation
    (theta weightZero weightOne upperZero upperOne audit returned : ℝ)
    (hthetaZero : 0 ≤ theta)
    (hthetaOne : theta ≤ 1)
    (hzero :
      weightZero * audit + (1 - weightZero) * returned ≤ upperZero)
    (hone :
      weightOne * audit + (1 - weightOne) * returned ≤ upperOne) :
    ((1 - theta) * weightZero + theta * weightOne) * audit
        + (1 - ((1 - theta) * weightZero + theta * weightOne)) * returned
      ≤ (1 - theta) * upperZero + theta * upperOne := by
  have hzeroScaled :=
    mul_le_mul_of_nonneg_left hzero (sub_nonneg.mpr hthetaOne)
  have honeScaled := mul_le_mul_of_nonneg_left hone hthetaZero
  calc
    ((1 - theta) * weightZero + theta * weightOne) * audit
          + (1 - ((1 - theta) * weightZero + theta * weightOne)) * returned
        = (1 - theta)
              * (weightZero * audit + (1 - weightZero) * returned)
            + theta * (weightOne * audit + (1 - weightOne) * returned) := by
              ring
    _ ≤ (1 - theta) * upperZero + theta * upperOne :=
      add_le_add hzeroScaled honeScaled

end CarmenQExact

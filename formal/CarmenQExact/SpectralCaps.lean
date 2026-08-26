import Mathlib

namespace CarmenQExact

noncomputable section

/-! Exact geometric bridge for every angular cap used by source cell 15818.

The Python cover represents a cube-face cell by the conic hull of its four
unnormalised corner rays.  The SOCP receives a binary64 vector `a` and imposes
the Lorentz inequality `‖v‖ ≤ ⟨a,v⟩`.  This file checks the actual binary64
coefficients exactly, not their intended irrational centres, and proves that
each complete source cell lies inside the corresponding Lorentz cone.
-/

abbrev Bloch3 := EuclideanSpace ℝ (Fin 3)

def InScaledCap (a v : Bloch3) : Prop :=
  ‖v‖ ≤ inner ℝ a v

theorem inScaledCap_of_sq (a v : Bloch3)
    (hnonneg : 0 ≤ inner ℝ a v)
    (hsq : ‖v‖ ^ 2 ≤ inner ℝ a v ^ 2) :
    InScaledCap a v := by
  unfold InScaledCap
  nlinarith [norm_nonneg v]

theorem inScaledCap_add {a v w : Bloch3}
    (hv : InScaledCap a v) (hw : InScaledCap a w) :
    InScaledCap a (v + w) := by
  unfold InScaledCap at *
  calc
    ‖v + w‖ ≤ ‖v‖ + ‖w‖ := norm_add_le v w
    _ ≤ inner ℝ a v + inner ℝ a w := add_le_add hv hw
    _ = inner ℝ a (v + w) := by rw [inner_add_right]

theorem inScaledCap_smul {a v : Bloch3} {r : ℝ}
    (hr : 0 ≤ r) (hv : InScaledCap a v) :
    InScaledCap a (r • v) := by
  unfold InScaledCap at *
  rw [norm_smul, Real.norm_eq_abs, abs_of_nonneg hr, inner_smul_right]
  exact mul_le_mul_of_nonneg_left hv hr

def InFourRayCone (c₀ c₁ c₂ c₃ v : Bloch3) : Prop :=
  ∃ r₀ r₁ r₂ r₃ : ℝ,
    0 ≤ r₀ ∧ 0 ≤ r₁ ∧ 0 ≤ r₂ ∧ 0 ≤ r₃ ∧
      v = r₀ • c₀ + r₁ • c₁ + r₂ • c₂ + r₃ • c₃

theorem fourRayCone_subset_scaledCap {a c₀ c₁ c₂ c₃ v : Bloch3}
    (hc₀ : InScaledCap a c₀) (hc₁ : InScaledCap a c₁)
    (hc₂ : InScaledCap a c₂) (hc₃ : InScaledCap a c₃)
    (hv : InFourRayCone c₀ c₁ c₂ c₃ v) :
    InScaledCap a v := by
  rcases hv with ⟨r₀, r₁, r₂, r₃, hr₀, hr₁, hr₂, hr₃, rfl⟩
  exact inScaledCap_add
    (inScaledCap_add
      (inScaledCap_add (inScaledCap_smul hr₀ hc₀) (inScaledCap_smul hr₁ hc₁))
      (inScaledCap_smul hr₂ hc₂))
    (inScaledCap_smul hr₃ hc₃)

def InTwoRayCone (c₀ c₁ v : Bloch3) : Prop :=
  ∃ r₀ r₁ : ℝ, 0 ≤ r₀ ∧ 0 ≤ r₁ ∧ v = r₀ • c₀ + r₁ • c₁

theorem twoRayCone_subset_scaledCap {a c₀ c₁ v : Bloch3}
    (hc₀ : InScaledCap a c₀) (hc₁ : InScaledCap a c₁)
    (hv : InTwoRayCone c₀ c₁ v) :
    InScaledCap a v := by
  rcases hv with ⟨r₀, r₁, hr₀, hr₁, rfl⟩
  exact inScaledCap_add (inScaledCap_smul hr₀ hc₀) (inScaledCap_smul hr₁ hc₁)

private theorem exact_cap_ray (a v : Bloch3)
    (hnonneg : 0 ≤ inner ℝ a v)
    (hsq : (∑ i, (v i) ^ 2) ≤ inner ℝ a v ^ 2) :
    InScaledCap a v := by
  apply inScaledCap_of_sq a v hnonneg
  rw [EuclideanSpace.real_norm_sq_eq]
  exact hsq

def source15818PlaneCap : Bloch3 :=
  !₂[1, 0, -7461808180621107 / 18014398509481984]

def source15818BaseSphereCap : Bloch3 :=
  !₂[77236096414237 / 140737488355328,
      -77236096414237 / 140737488355328,
      6590813560681557 / 9007199254740992]

def source15818Separator0Cap : Bloch3 :=
  !₂[-7536883653373163 / 9007199254740992,
      -7536883653373163 / 36028797018963968,
      353291421251867 / 562949953421312]

def source15818Separator1Cap : Bloch3 :=
  !₂[-4529103411573883 / 18014398509481984,
      4529103411573883 / 4503599627370496,
      -4529103411573883 / 18014398509481984]

def source15818Separator3Cap : Bloch3 :=
  !₂[4902898263487189 / 9007199254740992,
      -4902898263487189 / 4503599627370496,
      4902898263487189 / 9007199254740992]

def source15818PlaneCell (v : Bloch3) : Prop :=
  InTwoRayCone (!₂[1, 0, 0]) (!₂[1, 0, -1]) v

def source15818BaseSphereCell (v : Bloch3) : Prop :=
  InFourRayCone
    (!₂[1 / 2, -1, 1]) (!₂[1 / 2, -1 / 2, 1])
    (!₂[1, -1, 1]) (!₂[1, -1 / 2, 1]) v

def source15818Separator0Cell (v : Bloch3) : Prop :=
  InFourRayCone
    (!₂[-1, 0, 1]) (!₂[-1, 0, 1 / 2])
    (!₂[-1, -1 / 2, 1]) (!₂[-1, -1 / 2, 1 / 2]) v

def source15818Separator1Cell (v : Bloch3) : Prop :=
  InFourRayCone
    (!₂[-1 / 2, 1, -1 / 2]) (!₂[-1 / 2, 1, 0])
    (!₂[0, 1, -1 / 2]) (!₂[0, 1, 0]) v

def source15818Separator3Cell (v : Bloch3) : Prop :=
  InFourRayCone
    (!₂[1, -1, 1]) (!₂[1, -1, 0])
    (!₂[0, -1, 1]) (!₂[0, -1, 0]) v

private theorem plane_ray0 :
    InScaledCap source15818PlaneCap (!₂[1, 0, 0]) := by
  apply exact_cap_ray
  · norm_num [source15818PlaneCap, PiLp.inner_apply, Fin.sum_univ_succ]
  · norm_num [source15818PlaneCap, PiLp.inner_apply, Fin.sum_univ_succ]

private theorem plane_ray1 :
    InScaledCap source15818PlaneCap (!₂[1, 0, -1]) := by
  apply exact_cap_ray
  · norm_num [source15818PlaneCap, PiLp.inner_apply, Fin.sum_univ_succ]
  · norm_num [source15818PlaneCap, PiLp.inner_apply, Fin.sum_univ_succ]

theorem source15818_plane_cap_contains {v : Bloch3}
    (hv : source15818PlaneCell v) : InScaledCap source15818PlaneCap v :=
  twoRayCone_subset_scaledCap plane_ray0 plane_ray1 hv

private theorem baseSphere_ray0 :
    InScaledCap source15818BaseSphereCap (!₂[1 / 2, -1, 1]) := by
  apply exact_cap_ray <;>
    norm_num [source15818BaseSphereCap, PiLp.inner_apply, Fin.sum_univ_succ]

private theorem baseSphere_ray1 :
    InScaledCap source15818BaseSphereCap (!₂[1 / 2, -1 / 2, 1]) := by
  apply exact_cap_ray <;>
    norm_num [source15818BaseSphereCap, PiLp.inner_apply, Fin.sum_univ_succ]

private theorem baseSphere_ray2 :
    InScaledCap source15818BaseSphereCap (!₂[1, -1, 1]) := by
  apply exact_cap_ray <;>
    norm_num [source15818BaseSphereCap, PiLp.inner_apply, Fin.sum_univ_succ]

private theorem baseSphere_ray3 :
    InScaledCap source15818BaseSphereCap (!₂[1, -1 / 2, 1]) := by
  apply exact_cap_ray <;>
    norm_num [source15818BaseSphereCap, PiLp.inner_apply, Fin.sum_univ_succ]

theorem source15818_base_sphere_cap_contains {v : Bloch3}
    (hv : source15818BaseSphereCell v) :
    InScaledCap source15818BaseSphereCap v :=
  fourRayCone_subset_scaledCap baseSphere_ray0 baseSphere_ray1
    baseSphere_ray2 baseSphere_ray3 hv

private theorem separator0_ray0 :
    InScaledCap source15818Separator0Cap (!₂[-1, 0, 1]) := by
  apply exact_cap_ray <;>
    norm_num [source15818Separator0Cap, PiLp.inner_apply, Fin.sum_univ_succ]

private theorem separator0_ray1 :
    InScaledCap source15818Separator0Cap (!₂[-1, 0, 1 / 2]) := by
  apply exact_cap_ray <;>
    norm_num [source15818Separator0Cap, PiLp.inner_apply, Fin.sum_univ_succ]

private theorem separator0_ray2 :
    InScaledCap source15818Separator0Cap (!₂[-1, -1 / 2, 1]) := by
  apply exact_cap_ray <;>
    norm_num [source15818Separator0Cap, PiLp.inner_apply, Fin.sum_univ_succ]

private theorem separator0_ray3 :
    InScaledCap source15818Separator0Cap (!₂[-1, -1 / 2, 1 / 2]) := by
  apply exact_cap_ray <;>
    norm_num [source15818Separator0Cap, PiLp.inner_apply, Fin.sum_univ_succ]

theorem source15818_separator0_cap_contains {v : Bloch3}
    (hv : source15818Separator0Cell v) :
    InScaledCap source15818Separator0Cap v :=
  fourRayCone_subset_scaledCap separator0_ray0 separator0_ray1
    separator0_ray2 separator0_ray3 hv

private theorem separator1_ray0 :
    InScaledCap source15818Separator1Cap (!₂[-1 / 2, 1, -1 / 2]) := by
  apply exact_cap_ray <;>
    norm_num [source15818Separator1Cap, PiLp.inner_apply, Fin.sum_univ_succ]

private theorem separator1_ray1 :
    InScaledCap source15818Separator1Cap (!₂[-1 / 2, 1, 0]) := by
  apply exact_cap_ray <;>
    norm_num [source15818Separator1Cap, PiLp.inner_apply, Fin.sum_univ_succ]

private theorem separator1_ray2 :
    InScaledCap source15818Separator1Cap (!₂[0, 1, -1 / 2]) := by
  apply exact_cap_ray <;>
    norm_num [source15818Separator1Cap, PiLp.inner_apply, Fin.sum_univ_succ]

private theorem separator1_ray3 :
    InScaledCap source15818Separator1Cap (!₂[0, 1, 0]) := by
  apply exact_cap_ray <;>
    norm_num [source15818Separator1Cap, PiLp.inner_apply, Fin.sum_univ_succ]

theorem source15818_separator1_cap_contains {v : Bloch3}
    (hv : source15818Separator1Cell v) :
    InScaledCap source15818Separator1Cap v :=
  fourRayCone_subset_scaledCap separator1_ray0 separator1_ray1
    separator1_ray2 separator1_ray3 hv

private theorem separator3_ray0 :
    InScaledCap source15818Separator3Cap (!₂[1, -1, 1]) := by
  apply exact_cap_ray <;>
    norm_num [source15818Separator3Cap, PiLp.inner_apply, Fin.sum_univ_succ]

private theorem separator3_ray1 :
    InScaledCap source15818Separator3Cap (!₂[1, -1, 0]) := by
  apply exact_cap_ray <;>
    norm_num [source15818Separator3Cap, PiLp.inner_apply, Fin.sum_univ_succ]

private theorem separator3_ray2 :
    InScaledCap source15818Separator3Cap (!₂[0, -1, 1]) := by
  apply exact_cap_ray <;>
    norm_num [source15818Separator3Cap, PiLp.inner_apply, Fin.sum_univ_succ]

private theorem separator3_ray3 :
    InScaledCap source15818Separator3Cap (!₂[0, -1, 0]) := by
  apply exact_cap_ray <;>
    norm_num [source15818Separator3Cap, PiLp.inner_apply, Fin.sum_univ_succ]

theorem source15818_separator3_cap_contains {v : Bloch3}
    (hv : source15818Separator3Cell v) :
    InScaledCap source15818Separator3Cap v :=
  fourRayCone_subset_scaledCap separator3_ray0 separator3_ray1
    separator3_ray2 separator3_ray3 hv

structure Source15818AngularCells
    (basePlane baseSphere separator0 separator1 separator3 : Bloch3) : Prop where
  plane : source15818PlaneCell basePlane
  sphere : source15818BaseSphereCell baseSphere
  sep0 : source15818Separator0Cell separator0
  sep1 : source15818Separator1Cell separator1
  sep3 : source15818Separator3Cell separator3

structure Source15818LorentzCaps
    (basePlane baseSphere separator0 separator1 separator3 : Bloch3) : Prop where
  plane : InScaledCap source15818PlaneCap basePlane
  sphere : InScaledCap source15818BaseSphereCap baseSphere
  sep0 : InScaledCap source15818Separator0Cap separator0
  sep1 : InScaledCap source15818Separator1Cap separator1
  sep3 : InScaledCap source15818Separator3Cap separator3

theorem source15818_all_angular_caps_contain
    {basePlane baseSphere separator0 separator1 separator3 : Bloch3}
    (h : Source15818AngularCells basePlane baseSphere separator0 separator1 separator3) :
    Source15818LorentzCaps basePlane baseSphere separator0 separator1 separator3 :=
  ⟨source15818_plane_cap_contains h.plane,
    source15818_base_sphere_cap_contains h.sphere,
    source15818_separator0_cap_contains h.sep0,
    source15818_separator1_cap_contains h.sep1,
    source15818_separator3_cap_contains h.sep3⟩

end

end CarmenQExact

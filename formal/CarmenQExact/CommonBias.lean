import Mathlib

namespace CarmenQExact

/-!
# Common-Helstrom-bias support cuts

For a fully active terminal readout, positivity gives componentwise reserve
bounds `A * f i ≤ p i`.  A scalar minimisation over the *single shared*
Helstrom bias and longitudinal POVM closure supplies a support value
`m ≤ ∑ i, c i * f i`.  The theorem below is the exact algebraic bridge used
by every rational all-prior cut in the four-active certificate.  It makes
explicit why four independently chosen reserve vectors would be invalid.
-/

theorem commonBiasSupportCut
    (A m : ℝ)
    (p reserve coefficient : Fin 4 → ℝ)
    (hA : 0 ≤ A)
    (hcoefficient : ∀ i, 0 ≤ coefficient i)
    (hprior : ∀ i, A * reserve i ≤ p i)
    (hsupport : m ≤ ∑ i, coefficient i * reserve i) :
    m * A ≤ ∑ i, coefficient i * p i := by
  calc
    m * A ≤ (∑ i, coefficient i * reserve i) * A :=
      mul_le_mul_of_nonneg_right hsupport hA
    _ = ∑ i, coefficient i * (A * reserve i) := by
      rw [Finset.sum_mul]
      apply Finset.sum_congr rfl
      intro i _
      ring
    _ ≤ ∑ i, coefficient i * p i := by
      apply Finset.sum_le_sum
      intro i _
      exact mul_le_mul_of_nonneg_left (hprior i) (hcoefficient i)

noncomputable def commonBiasDirectionOne : Fin 4 → ℝ :=
  ![1, 1 / 3, 2 / 5, 2 / 5]

noncomputable def commonBiasDirectionTwo : Fin 4 → ℝ :=
  ![1, 1 / 2, 1 / 3, 1 / 3]

noncomputable def commonBiasDirectionThree : Fin 4 → ℝ :=
  ![1, 10 / 29, 1 / 2, 16 / 27]

theorem commonBiasDirectionOne_nonnegative :
    ∀ i, 0 ≤ commonBiasDirectionOne i := by
  intro i
  fin_cases i <;> norm_num [commonBiasDirectionOne]

theorem commonBiasDirectionTwo_nonnegative :
    ∀ i, 0 ≤ commonBiasDirectionTwo i := by
  intro i
  fin_cases i <;> norm_num [commonBiasDirectionTwo]

theorem commonBiasDirectionThree_nonnegative :
    ∀ i, 0 ≤ commonBiasDirectionThree i := by
  intro i
  fin_cases i <;> norm_num [commonBiasDirectionThree]

theorem commonBiasDirectionOneCut
    (A m : ℝ) (p reserve : Fin 4 → ℝ)
    (hA : 0 ≤ A)
    (hprior : ∀ i, A * reserve i ≤ p i)
    (hsupport : m ≤ ∑ i, commonBiasDirectionOne i * reserve i) :
    m * A ≤ ∑ i, commonBiasDirectionOne i * p i :=
  commonBiasSupportCut A m p reserve commonBiasDirectionOne hA
    commonBiasDirectionOne_nonnegative hprior hsupport

theorem commonBiasDirectionTwoCut
    (A m : ℝ) (p reserve : Fin 4 → ℝ)
    (hA : 0 ≤ A)
    (hprior : ∀ i, A * reserve i ≤ p i)
    (hsupport : m ≤ ∑ i, commonBiasDirectionTwo i * reserve i) :
    m * A ≤ ∑ i, commonBiasDirectionTwo i * p i :=
  commonBiasSupportCut A m p reserve commonBiasDirectionTwo hA
    commonBiasDirectionTwo_nonnegative hprior hsupport

theorem commonBiasDirectionThreeCut
    (A m : ℝ) (p reserve : Fin 4 → ℝ)
    (hA : 0 ≤ A)
    (hprior : ∀ i, A * reserve i ≤ p i)
    (hsupport : m ≤ ∑ i, commonBiasDirectionThree i * reserve i) :
    m * A ≤ ∑ i, commonBiasDirectionThree i * p i :=
  commonBiasSupportCut A m p reserve commonBiasDirectionThree hA
    commonBiasDirectionThree_nonnegative hprior hsupport

end CarmenQExact

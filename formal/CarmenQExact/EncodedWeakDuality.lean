import CarmenQExact.WeakDuality

open scoped BigOperators

namespace CarmenQExact

/-!
# Weak duality for the encoded sparse canonical program

`CertificateData` stores the canonical matrix by sparse columns and stores the
right-hand side by sparse row entries.  This module gives those arrays a
mathematical meaning and connects the exact checker predicates to the generic
weak-duality theorem.  Repeated row entries are interpreted additively; no
ordering or uniqueness assumption is required.
-/

private def sparseValueAt {rows : Nat}
    (entries : Array (RayEntry rows)) (row : Fin rows) : ℚ :=
  (entries.toList.map fun entry =>
    if entry.row = row then entry.value else 0).sum

def decodedMatrixQ (data : CertificateData)
    (row : Fin data.rows) (column : Fin data.variableCount) : ℚ :=
  sparseValueAt (data.columns.getD column.val #[]) row

def decodedRightQ (data : CertificateData) (row : Fin data.rows) : ℚ :=
  sparseValueAt data.rightEntries row

def decodedObjectiveQ (data : CertificateData)
    (column : Fin data.variableCount) : ℚ :=
  data.objective.getD column.val 0

def decodedDualQ (data : CertificateData) (row : Fin data.rows) : ℚ :=
  data.dual row

def decodedMatrix (data : CertificateData)
    (row : Fin data.rows) (column : Fin data.variableCount) : ℝ :=
  (decodedMatrixQ data row column : ℝ)

def decodedRight (data : CertificateData) (row : Fin data.rows) : ℝ :=
  (decodedRightQ data row : ℝ)

def decodedObjective (data : CertificateData)
    (column : Fin data.variableCount) : ℝ :=
  (decodedObjectiveQ data column : ℝ)

private theorem foldl_entry_dot_eq_acc_add_sum
    {rows : Nat} (entries : List (RayEntry rows))
    (dual : Fin rows → ℚ) (accumulator : ℚ) :
    entries.foldl
        (fun total entry => total + entry.value * dual entry.row)
        accumulator =
      accumulator +
        (entries.map fun entry => entry.value * dual entry.row).sum := by
  induction entries generalizing accumulator with
  | nil => simp
  | cons entry rest inductionHypothesis =>
      simp only [List.foldl_cons, List.map_cons, List.sum_cons]
      rw [inductionHypothesis]
      ring

private theorem list_entry_dot_eq_sum_sparseValue
    {rows : Nat} (entries : List (RayEntry rows))
    (dual : Fin rows → ℚ) :
    (entries.map fun entry => entry.value * dual entry.row).sum =
      ∑ row : Fin rows,
        (entries.map fun entry =>
          if entry.row = row then entry.value else 0).sum * dual row := by
  induction entries with
  | nil => simp
  | cons entry rest inductionHypothesis =>
      simp only [List.map_cons, List.sum_cons]
      simp_rw [add_mul]
      rw [Finset.sum_add_distrib]
      rw [← inductionHypothesis]
      have hsingle :
          (∑ row : Fin rows,
            (if entry.row = row then entry.value else 0) * dual row) =
            entry.value * dual entry.row := by
        simp
      rw [hsingle]

private theorem array_entry_dot_eq_sum_sparseValue
    {rows : Nat} (entries : Array (RayEntry rows))
    (dual : Fin rows → ℚ) :
    entries.foldl
        (fun total entry => total + entry.value * dual entry.row) 0 =
      ∑ row : Fin rows, sparseValueAt entries row * dual row := by
  rw [← Array.foldl_toList]
  rw [foldl_entry_dot_eq_acc_add_sum entries.toList dual 0]
  simp only [zero_add]
  exact list_entry_dot_eq_sum_sparseValue entries.toList dual

private theorem array_entry_dot_eq_acc_add_sum_sparseValue
    {rows : Nat} (entries : Array (RayEntry rows))
    (dual : Fin rows → ℚ) (accumulator : ℚ) :
    entries.foldl
        (fun total entry => total + entry.value * dual entry.row)
        accumulator =
      accumulator +
        ∑ row : Fin rows, sparseValueAt entries row * dual row := by
  rw [← Array.foldl_toList]
  rw [foldl_entry_dot_eq_acc_add_sum entries.toList dual accumulator]
  rw [list_entry_dot_eq_sum_sparseValue entries.toList dual]
  simp only [sparseValueAt]

theorem stationarityAt_eq_decoded_sum
    (data : CertificateData) (column : Fin data.variableCount) :
    stationarityAt data column.val =
      decodedObjectiveQ data column +
        ∑ row : Fin data.rows,
          decodedMatrixQ data row column * decodedDualQ data row := by
  unfold stationarityAt decodedObjectiveQ decodedMatrixQ decodedDualQ
  exact array_entry_dot_eq_acc_add_sum_sparseValue
    (data.columns.getD column.val #[])
    data.dual
    (data.objective.getD column.val 0)

theorem certifiedUpper_eq_decoded_dotQ (data : CertificateData) :
    certifiedUpper data =
      ∑ row : Fin data.rows,
        decodedRightQ data row * decodedDualQ data row := by
  unfold certifiedUpper decodedRightQ decodedDualQ
  exact array_entry_dot_eq_sum_sparseValue
    data.rightEntries data.dual

theorem stationarityAt_cast_eq_decoded_sum
    (data : CertificateData) (column : Fin data.variableCount) :
    (stationarityAt data column.val : ℝ) =
      (∑ row : Fin data.rows,
        decodedMatrix data row column * certificateDual data row) +
          decodedObjective data column := by
  have rationalIdentity := stationarityAt_eq_decoded_sum data column
  have realIdentity := congrArg (fun value : ℚ => (value : ℝ)) rationalIdentity
  simpa [decodedMatrix, decodedObjective, certificateDual,
    decodedDualQ, add_comm] using realIdentity

theorem certifiedUpper_cast_eq_decoded_dot (data : CertificateData) :
    (certifiedUpper data : ℝ) =
      dot (decodedRight data) (certificateDual data) := by
  have rationalIdentity := certifiedUpper_eq_decoded_dotQ data
  have realIdentity := congrArg (fun value : ℚ => (value : ℝ)) rationalIdentity
  simpa [dot, decodedRight, certificateDual,
    decodedDualQ] using realIdentity

theorem exactCertificate_dimensions_eq
    (data : CertificateData) (proof : ExactCertificate data) :
    data.zeroDim + data.nonnegativeDim + data.socSizes.toList.sum =
      data.rows := by
  simpa [dimensionsOK] using proof.dimensions

theorem exactCertificate_dualCone_to_real
    (data : CertificateData) (proof : ExactCertificate data) :
    ProductConeDual (rationalFunctionAsReal (rationalDualAt data.dual))
      data.zeroDim data.nonnegativeDim data.socSizes.toList := by
  exact rationalProductConeDual_to_real
    (rationalDualAt data.dual) data.zeroDim data.nonnegativeDim data.socSizes.toList
    proof.dualCone

theorem extendFin_certificateDual_eq_rationalFunctionAsReal
    (data : CertificateData) :
    extendFin (certificateDual data) =
      rationalFunctionAsReal (rationalDualAt data.dual) := by
  funext index
  unfold extendFin certificateDual rationalFunctionAsReal rationalDualAt
  by_cases inRows : index < data.rows
  · simp [inRows]

  · simp [inRows]
theorem exactCertificate_decoded_stationary
    (data : CertificateData) (proof : ExactCertificate data) :
    ∀ column : Fin data.variableCount,
      (∑ row : Fin data.rows,
        decodedMatrix data row column * certificateDual data row) +
          decodedObjective data column = 0 := by
  intro column
  rw [← stationarityAt_cast_eq_decoded_sum data column]
  exact_mod_cast proof.stationarity column.val column.isLt

theorem exactCertificate_decoded_weak_duality
    (data : CertificateData)
    (proof : ExactCertificate data)
    (point : Fin data.variableCount → ℝ)
    (slack : Fin data.rows → ℝ)
    (feasible : CanonicalFeasible
      (decodedMatrix data) (decodedRight data) point slack)
    (hslack : ProductConePoint (extendFin slack)
      data.zeroDim data.nonnegativeDim data.socSizes.toList) :
    -(∑ column : Fin data.variableCount,
        decodedObjective data column * point column) ≤
      (certifiedUpper data : ℝ) := by
  have bound := canonical_productCone_weak_duality
    (decodedMatrix data) (decodedRight data) (decodedObjective data)
    point slack (certificateDual data)
    data.zeroDim data.nonnegativeDim data.socSizes.toList
    feasible (exactCertificate_decoded_stationary data proof)
    (exactCertificate_dimensions_eq data proof) hslack
    (by
      rw [extendFin_certificateDual_eq_rationalFunctionAsReal data]
      exact exactCertificate_dualCone_to_real data proof)
  rw [← certifiedUpper_cast_eq_decoded_dot] at bound
  exact bound

theorem exactCertificate_decoded_strict_target
    (data : CertificateData)
    (proof : ExactCertificate data)
    (point : Fin data.variableCount → ℝ)
    (slack : Fin data.rows → ℝ)
    (feasible : CanonicalFeasible
      (decodedMatrix data) (decodedRight data) point slack)
    (hslack : ProductConePoint (extendFin slack)
      data.zeroDim data.nonnegativeDim data.socSizes.toList) :
    -(∑ column : Fin data.variableCount,
        decodedObjective data column * point column) <
      (data.target : ℝ) := by
  apply lt_of_le_of_lt
    (exactCertificate_decoded_weak_duality
      data proof point slack feasible hslack)
  exact_mod_cast proof.upper

end CarmenQExact


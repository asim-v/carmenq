import Mathlib.Data.Rat.Defs
import Mathlib.Algebra.BigOperators.Group.Finset.Basic

open scoped BigOperators

namespace CarmenQExact

structure RayEntry (rowCount : Nat) where
  row : Fin rowCount
  value : ℚ
deriving DecidableEq, Repr


structure CertificateData where
  rows : Nat
  variableCount : Nat
  zeroDim : Nat
  nonnegativeDim : Nat
  socSizes : Array Nat
  columns : Array (Array (RayEntry rows))
  rightEntries : Array (RayEntry rows)
  objective : Array ℚ
  dual : Fin rows → ℚ
  target : ℚ
def rationalDualAt {rows : Nat}
    (dual : Fin rows → ℚ) (index : Nat) : ℚ :=
  if h : index < rows then dual ⟨index, h⟩ else 0


def stationarityAt (data : CertificateData) (column : Nat) : ℚ :=
  (data.columns.getD column #[]).foldl
    (fun total entry =>
      total + entry.value * data.dual entry.row)
    (data.objective.getD column 0)

def certifiedUpper (data : CertificateData) : ℚ :=
  data.rightEntries.foldl
    (fun total entry =>
      total + entry.value * data.dual entry.row)
    0

def rationalSpatialSquare
    (dual : Nat → ℚ) (start size : Nat) : ℚ :=
  ∑ offset ∈ Finset.range (size - 1),
    let value := dual (start + offset + 1)
    value * value

def rationalLorentzBlock
    (dual : Nat → ℚ) (start size : Nat) : Prop :=
  0 < size ∧
    0 ≤ dual start ∧
    rationalSpatialSquare dual start size ≤
      dual start * dual start

def rationalLorentzBlocks
    (dual : Nat → ℚ) : Nat → List Nat → Prop
  | _, [] => True
  | start, size :: rest =>
      rationalLorentzBlock dual start size ∧
        rationalLorentzBlocks dual (start + size) rest

def rationalNonnegativeSegment
    (dual : Nat → ℚ) (start size : Nat) : Prop :=
  ∀ offset, offset < size →
    0 ≤ dual (start + offset)

def rationalProductConeDual
    (dual : Nat → ℚ) (zeroDim nonnegativeDim : Nat)
    (socSizes : List Nat) : Prop :=
  rationalNonnegativeSegment dual zeroDim nonnegativeDim ∧
    rationalLorentzBlocks dual
      (zeroDim + nonnegativeDim) socSizes

def socBlockOK (dual : Nat → ℚ) (start size : Nat) : Bool :=
  if size = 0 then
    false
  else
    let time := dual start
    let spatialSquare := rationalSpatialSquare dual start size
    decide (0 ≤ time ∧ spatialSquare ≤ time * time)

def socBlocksOK (dual : Nat → ℚ) : Nat → List Nat → Bool
  | _, [] => true
  | start, size :: rest =>
      socBlockOK dual start size && socBlocksOK dual (start + size) rest
theorem socBlocksOK_append
    (dual : Nat → ℚ) (start : Nat) (first second : List Nat) :
    socBlocksOK dual start (first ++ second) =
      (socBlocksOK dual start first &&
        socBlocksOK dual (start + first.sum) second) := by
  induction first generalizing start with
  | nil => simp [socBlocksOK]
  | cons size rest ih =>
      simp only [List.cons_append, List.sum_cons, socBlocksOK]
      rw [ih]
      simp only [Bool.and_assoc, Nat.add_assoc]

theorem socBlocksOK_append_of_true
    (dual : Nat → ℚ) (start : Nat) (first second : List Nat)
    (hfirst : socBlocksOK dual start first = true)
    (hsecond :
      socBlocksOK dual (start + first.sum) second = true) :
    socBlocksOK dual start (first ++ second) = true := by
  rw [socBlocksOK_append, Bool.and_eq_true]
  exact ⟨hfirst, hsecond⟩

theorem socBlockOK_sound
    (dual : Nat → ℚ) (start size : Nat)
    (h : socBlockOK dual start size = true) :
    rationalLorentzBlock dual start size := by
  unfold socBlockOK at h
  by_cases hzero : size = 0
  · simp [hzero] at h
  · simp only [hzero, ↓reduceIte] at h
    have hcone := of_decide_eq_true h
    exact ⟨Nat.pos_of_ne_zero hzero, hcone⟩

theorem socBlocksOK_sound
    (dual : Nat → ℚ) (start : Nat) (sizes : List Nat)
    (h : socBlocksOK dual start sizes = true) :
    rationalLorentzBlocks dual start sizes := by
  induction sizes generalizing start with
  | nil => trivial
  | cons size rest ih =>
      simp only [socBlocksOK, Bool.and_eq_true] at h
      exact ⟨socBlockOK_sound dual start size h.1,
        ih (start + size) h.2⟩


def dimensionsOK (data : CertificateData) : Bool :=
  data.zeroDim + data.nonnegativeDim + data.socSizes.toList.sum = data.rows

def indicesOK (data : CertificateData) : Bool :=
  data.objective.size = data.variableCount &&
  data.columns.size = data.variableCount

def nonnegativeRangeOK
    (dual : Nat → ℚ) (start size : Nat) : Bool :=
  (List.range size).all
    (fun offset => decide (0 ≤ dual (start + offset)))

theorem nonnegativeRangeOK_sound
    (dual : Nat → ℚ) (start size : Nat)
    (h : nonnegativeRangeOK dual start size = true) :
    rationalNonnegativeSegment dual start size := by
  intro offset hoffset
  unfold nonnegativeRangeOK at h
  simp only [List.all_eq_true] at h
  exact of_decide_eq_true
    (h offset (List.mem_range.mpr hoffset))

theorem rationalNonnegativeSegment_append
    (dual : Nat → ℚ) (start left right : Nat)
    (hleft : rationalNonnegativeSegment dual start left)
    (hright :
      rationalNonnegativeSegment dual (start + left) right) :
    rationalNonnegativeSegment dual start (left + right) := by
  intro offset hoffset
  by_cases hsplit : offset < left
  · exact hleft offset hsplit
  · have hleft_le : left ≤ offset := Nat.le_of_not_gt hsplit
    let shifted := offset - left
    have hshifted : shifted < right := by omega
    have hindex :
        start + offset = (start + left) + shifted := by omega
    rw [hindex]
    exact hright shifted hshifted

theorem rationalLorentzBlocks_append
    (dual : Nat → ℚ) (start : Nat) (left right : List Nat)
    (hleft : rationalLorentzBlocks dual start left)
    (hright :
      rationalLorentzBlocks dual (start + left.sum) right) :
    rationalLorentzBlocks dual start (left ++ right) := by
  induction left generalizing start with
  | nil => simpa using hright
  | cons size rest ih =>
      rcases hleft with ⟨hhead, htail⟩
      constructor
      · exact hhead
      · apply ih (start := start + size) htail
        simpa [List.sum_cons, Nat.add_assoc, Nat.add_left_comm,
          Nat.add_comm] using hright

theorem rationalProductConeDual_of_parts
    (dual : Nat → ℚ) (zeroDim nonnegativeDim : Nat)
    (socSizes : List Nat)
    (hnonnegative :
      rationalNonnegativeSegment dual zeroDim nonnegativeDim)
    (hsoc :
      rationalLorentzBlocks dual
        (zeroDim + nonnegativeDim) socSizes) :
    rationalProductConeDual dual zeroDim nonnegativeDim socSizes :=
  ⟨hnonnegative, hsoc⟩

def dualConeOK (data : CertificateData) : Bool :=
  nonnegativeRangeOK (rationalDualAt data.dual) data.zeroDim data.nonnegativeDim &&
    socBlocksOK (rationalDualAt data.dual)
      (data.zeroDim + data.nonnegativeDim) data.socSizes.toList

theorem dualConeOK_of_parts
    (data : CertificateData)
    (hnonnegative :
      nonnegativeRangeOK (rationalDualAt data.dual)
        data.zeroDim data.nonnegativeDim = true)
    (hsoc :
      socBlocksOK (rationalDualAt data.dual)
        (data.zeroDim + data.nonnegativeDim)
        data.socSizes.toList = true) :
    dualConeOK data = true := by
  unfold dualConeOK
  rw [Bool.and_eq_true]
  exact ⟨hnonnegative, hsoc⟩

theorem dualConeOK_sound (data : CertificateData)
    (h : dualConeOK data = true) :
    rationalProductConeDual (rationalDualAt data.dual)
      data.zeroDim data.nonnegativeDim data.socSizes.toList := by
  unfold dualConeOK at h
  rw [Bool.and_eq_true] at h
  constructor
  · intro offset hoffset
    have hall := h.1
    unfold nonnegativeRangeOK at hall
    simp only [List.all_eq_true] at hall
    exact of_decide_eq_true
      (hall offset (List.mem_range.mpr hoffset))
  · exact socBlocksOK_sound (rationalDualAt data.dual)
      (data.zeroDim + data.nonnegativeDim)
      data.socSizes.toList h.2

def stationarityOK (data : CertificateData) : Bool :=
  (List.range data.variableCount).all
    (fun column => decide (stationarityAt data column = 0))
theorem stationarityOK_of_columns
    (data : CertificateData)
    (hcolumns : ∀ column, column < data.variableCount →
      stationarityAt data column = 0) :
    stationarityOK data = true := by
  unfold stationarityOK
  simp only [List.all_eq_true]
  intro column hmember
  exact decide_eq_true
    (hcolumns column (List.mem_range.mp hmember))

theorem stationarityOK_sound
    (data : CertificateData) (h : stationarityOK data = true) :
    ∀ column, column < data.variableCount →
      stationarityAt data column = 0 := by
  intro column hcolumn
  unfold stationarityOK at h
  simp only [List.all_eq_true] at h
  exact of_decide_eq_true
    (h column (List.mem_range.mpr hcolumn))

structure CertificateProof (data : CertificateData) : Prop where
  dimensions : dimensionsOK data = true
  indices : indicesOK data = true
  dualCone :
    rationalProductConeDual (rationalDualAt data.dual)
      data.zeroDim data.nonnegativeDim data.socSizes.toList
  stationarity :
    ∀ column, column < data.variableCount →
      stationarityAt data column = 0
  upper : certifiedUpper data < data.target

abbrev ExactCertificate (data : CertificateData) : Prop :=
  CertificateProof data


def checkCertificate (data : CertificateData) : Bool :=
  dimensionsOK data &&
  indicesOK data &&
  dualConeOK data &&
  stationarityOK data &&
  decide (certifiedUpper data < data.target)

theorem checked_dimensionsOK (data : CertificateData)
    (h : checkCertificate data = true) :
    dimensionsOK data = true := by
  unfold checkCertificate at h
  repeat' rw [Bool.and_eq_true] at h
  exact h.1.1.1.1

theorem checked_indicesOK (data : CertificateData)
    (h : checkCertificate data = true) :
    indicesOK data = true := by
  unfold checkCertificate at h
  repeat' rw [Bool.and_eq_true] at h
  exact h.1.1.1.2

theorem checked_dualConeOK (data : CertificateData)
    (h : checkCertificate data = true) :
    dualConeOK data = true := by
  unfold checkCertificate at h
  repeat' rw [Bool.and_eq_true] at h
  exact h.1.1.2

theorem checked_stationarityOK (data : CertificateData)
    (h : checkCertificate data = true) :
    stationarityOK data = true := by
  unfold checkCertificate at h
  repeat' rw [Bool.and_eq_true] at h
  exact h.1.2

theorem checked_upper_lt_target (data : CertificateData)
    (h : checkCertificate data = true) : certifiedUpper data < data.target := by
  unfold checkCertificate at h
  repeat' rw [Bool.and_eq_true] at h
  exact of_decide_eq_true h.2

theorem certificateProof_of_check (data : CertificateData)
    (h : checkCertificate data = true) : CertificateProof data := by
  refine ⟨checked_dimensionsOK data h, checked_indicesOK data h,
    ?_, ?_, checked_upper_lt_target data h⟩
  · exact dualConeOK_sound data (checked_dualConeOK data h)
  · exact stationarityOK_sound data (checked_stationarityOK data h)

end CarmenQExact

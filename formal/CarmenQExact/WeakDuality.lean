import CarmenQExact.Checker
import Mathlib.Analysis.Real.Sqrt
import Mathlib.Data.Rat.BigOperators

open scoped BigOperators

namespace CarmenQExact

def dot {n : Nat} (left right : Fin n → ℝ) : ℝ :=
  ∑ index, left index * right index

structure CanonicalFeasible {rows variableCount : Nat}
    (matrix : Fin rows → Fin variableCount → ℝ)
    (right : Fin rows → ℝ)
    (point : Fin variableCount → ℝ)
    (slack : Fin rows → ℝ) : Prop where
  equation : ∀ row,
    (∑ column, matrix row column * point column) + slack row = right row

theorem canonical_weak_duality
    {rows variableCount : Nat}
    (matrix : Fin rows → Fin variableCount → ℝ)
    (right : Fin rows → ℝ)
    (objective : Fin variableCount → ℝ)
    (point : Fin variableCount → ℝ)
    (slack dual : Fin rows → ℝ)
    (feasible : CanonicalFeasible matrix right point slack)
    (stationary : ∀ column,
      (∑ row, matrix row column * dual row) + objective column = 0)
    (conePair : 0 ≤ dot dual slack) :
    -(∑ column, objective column * point column) ≤ dot right dual := by
  have hstationary : ∀ column,
      -objective column = ∑ row, matrix row column * dual row := by
    intro column
    linarith [stationary column]
  have hcross :
      ∑ column, point column * (∑ row, matrix row column * dual row) =
        ∑ row, (∑ column, matrix row column * point column) * dual row := by
    simp_rw [Finset.mul_sum, Finset.sum_mul]
    rw [Finset.sum_comm]
    apply Finset.sum_congr rfl
    intro row _
    apply Finset.sum_congr rfl
    intro column _
    ring
  calc
    -(∑ column, objective column * point column) =
        ∑ column, point column * (-objective column) := by
          rw [← Finset.sum_neg_distrib]
          apply Finset.sum_congr rfl
          intro column _
          ring
    _ = ∑ column, point column * (∑ row, matrix row column * dual row) := by
          apply Finset.sum_congr rfl
          intro column _
          rw [hstationary column]
    _ = ∑ row, (∑ column, matrix row column * point column) * dual row := hcross
    _ ≤ ∑ row, ((∑ column, matrix row column * point column) + slack row) * dual row := by
          have hpair : 0 ≤ ∑ row, slack row * dual row := by
            simpa [dot, mul_comm] using conePair
          have hdecomp :
              ∑ row,
                  ((∑ column, matrix row column * point column) + slack row) *
                    dual row =
                (∑ row, (∑ column, matrix row column * point column) * dual row) +
                  ∑ row, slack row * dual row := by
            rw [← Finset.sum_add_distrib]
            apply Finset.sum_congr rfl
            intro row _
            ring
          rw [hdecomp]
          exact le_add_of_nonneg_right hpair
    _ = dot right dual := by
          unfold dot
          apply Finset.sum_congr rfl
          intro row _
          rw [feasible.equation row]

/-! ### Exact product-cone pairing

The canonical certificate format orders its rows as a zero cone, a
nonnegative cone, and then a list of Lorentz blocks.  The following
definitions use `Nat`-indexed functions because this is exactly the layout
emitted by the certificate exporter. -/

def segmentDot (left right : Nat → ℝ) (start size : Nat) : ℝ :=
  ∑ offset ∈ Finset.range size, left (start + offset) * right (start + offset)

def lorentzBlock (values : Nat → ℝ) (start size : Nat) : Prop :=
  0 < size ∧
    0 ≤ values start ∧
    (∑ offset ∈ Finset.range (size - 1),
      values (start + offset + 1) ^ 2) ≤ values start ^ 2

theorem lorentzBlock_pair_nonnegative
    (left right : Nat → ℝ) (start size : Nat)
    (hleft : lorentzBlock left start size)
    (hright : lorentzBlock right start size) :
    0 ≤ segmentDot left right start size := by
  rcases hleft with ⟨hsize, hleftTime, hleftSquare⟩
  rcases hright with ⟨_, hrightTime, hrightSquare⟩
  have hsizeEq : size = (size - 1) + 1 := by omega
  have hleftRoot :
      Real.sqrt (∑ offset ∈ Finset.range (size - 1),
        left (start + offset + 1) ^ 2) ≤ left start :=
    Real.sqrt_le_iff.mpr ⟨hleftTime, hleftSquare⟩
  have hrightRoot :
      Real.sqrt (∑ offset ∈ Finset.range (size - 1),
        right (start + offset + 1) ^ 2) ≤ right start :=
    Real.sqrt_le_iff.mpr ⟨hrightTime, hrightSquare⟩
  have hcs := Real.sum_mul_le_sqrt_mul_sqrt
    (Finset.range (size - 1))
    (fun offset => -left (start + offset + 1))
    (fun offset => right (start + offset + 1))
  have hrootProduct :
      Real.sqrt (∑ offset ∈ Finset.range (size - 1),
        left (start + offset + 1) ^ 2) *
          Real.sqrt (∑ offset ∈ Finset.range (size - 1),
            right (start + offset + 1) ^ 2) ≤
        left start * right start := by
    exact mul_le_mul hleftRoot hrightRoot (Real.sqrt_nonneg _) hleftTime
  have hspatial :
      -(∑ offset ∈ Finset.range (size - 1),
        left (start + offset + 1) * right (start + offset + 1)) ≤
        left start * right start := by
    have hcs' :
        -(∑ offset ∈ Finset.range (size - 1),
          left (start + offset + 1) * right (start + offset + 1)) ≤
          Real.sqrt (∑ offset ∈ Finset.range (size - 1),
            left (start + offset + 1) ^ 2) *
            Real.sqrt (∑ offset ∈ Finset.range (size - 1),
              right (start + offset + 1) ^ 2) := by
      simpa only [neg_mul, neg_sq, Finset.sum_neg_distrib] using hcs
    exact hcs'.trans hrootProduct
  have hspatial' :
      -(∑ offset ∈ Finset.range (size - 1),
        left (start + (offset + 1)) * right (start + (offset + 1))) ≤
        left start * right start := by
    simpa only [Nat.add_assoc] using hspatial
  rw [segmentDot, hsizeEq, Finset.sum_range_succ']
  simp only [Nat.add_zero]
  linarith [hspatial']

def nonnegativeSegment (values : Nat → ℝ) (start size : Nat) : Prop :=
  ∀ offset, offset < size → 0 ≤ values (start + offset)

def zeroSegment (values : Nat → ℝ) (size : Nat) : Prop :=
  ∀ offset, offset < size → values offset = 0

def lorentzBlocks (values : Nat → ℝ) : Nat → List Nat → Prop
  | _, [] => True
  | start, size :: rest =>
      lorentzBlock values start size ∧ lorentzBlocks values (start + size) rest

def lorentzBlocksDot (left right : Nat → ℝ) : Nat → List Nat → ℝ
  | _, [] => 0
  | start, size :: rest =>
      segmentDot left right start size +
        lorentzBlocksDot left right (start + size) rest

theorem lorentzBlocks_pair_nonnegative
    (left right : Nat → ℝ) (start : Nat) (sizes : List Nat)
    (hleft : lorentzBlocks left start sizes)
    (hright : lorentzBlocks right start sizes) :
    0 ≤ lorentzBlocksDot left right start sizes := by
  induction sizes generalizing start with
  | nil => simp [lorentzBlocksDot]
  | cons size rest ih =>
      rcases hleft with ⟨hleftHead, hleftTail⟩
      rcases hright with ⟨hrightHead, hrightTail⟩
      simp only [lorentzBlocksDot]
      exact add_nonneg
        (lorentzBlock_pair_nonnegative left right start size hleftHead hrightHead)
        (ih (start + size) hleftTail hrightTail)

theorem nonnegativeSegment_pair_nonnegative
    (left right : Nat → ℝ) (start size : Nat)
    (hleft : nonnegativeSegment left start size)
    (hright : nonnegativeSegment right start size) :
    0 ≤ segmentDot left right start size := by
  apply Finset.sum_nonneg
  intro offset hoffset
  have hoffset' : offset < size := Finset.mem_range.mp hoffset
  exact mul_nonneg (hleft offset hoffset') (hright offset hoffset')

theorem segmentDot_add
    (left right : Nat → ℝ) (start first second : Nat) :
    segmentDot left right start (first + second) =
      segmentDot left right start first +
        segmentDot left right (start + first) second := by
  simp only [segmentDot, Finset.sum_range_add, Nat.add_assoc]

theorem zeroSegment_pair_zero
    (left right : Nat → ℝ) (size : Nat)
    (hleft : zeroSegment left size) :
    segmentDot left right 0 size = 0 := by
  unfold segmentDot
  apply Finset.sum_eq_zero
  intro offset hoffset
  simp only [Nat.zero_add]
  rw [hleft offset (Finset.mem_range.mp hoffset), zero_mul]

theorem lorentzBlocksDot_eq_segmentDot
    (left right : Nat → ℝ) (start : Nat) (sizes : List Nat) :
    lorentzBlocksDot left right start sizes =
      segmentDot left right start sizes.sum := by
  induction sizes generalizing start with
  | nil => simp [lorentzBlocksDot, segmentDot]
  | cons size rest ih =>
      rw [lorentzBlocksDot, List.sum_cons, ih]
      exact (segmentDot_add left right start size rest.sum).symm

structure ProductConePoint
    (values : Nat → ℝ) (zeroDim nonnegativeDim : Nat)
    (socSizes : List Nat) : Prop where
  zero : zeroSegment values zeroDim
  nonnegative : nonnegativeSegment values zeroDim nonnegativeDim
  lorentz :
    lorentzBlocks values (zeroDim + nonnegativeDim) socSizes

structure ProductConeDual
    (values : Nat → ℝ) (zeroDim nonnegativeDim : Nat)
    (socSizes : List Nat) : Prop where
  nonnegative : nonnegativeSegment values zeroDim nonnegativeDim
  lorentz :
    lorentzBlocks values (zeroDim + nonnegativeDim) socSizes

def productConeDot
    (left right : Nat → ℝ) (zeroDim nonnegativeDim : Nat)
    (socSizes : List Nat) : ℝ :=
  segmentDot left right 0 zeroDim +
    segmentDot left right zeroDim nonnegativeDim +
    lorentzBlocksDot left right (zeroDim + nonnegativeDim) socSizes

theorem productCone_pair_nonnegative
    (left right : Nat → ℝ) (zeroDim nonnegativeDim : Nat)
    (socSizes : List Nat)
    (hleft : ProductConePoint left zeroDim nonnegativeDim socSizes)
    (hright : ProductConeDual right zeroDim nonnegativeDim socSizes) :
    0 ≤ productConeDot left right zeroDim nonnegativeDim socSizes := by
  rw [productConeDot, zeroSegment_pair_zero left right zeroDim hleft.zero]
  simp only [zero_add]
  exact add_nonneg
    (nonnegativeSegment_pair_nonnegative
      left right zeroDim nonnegativeDim hleft.nonnegative hright.nonnegative)
    (lorentzBlocks_pair_nonnegative
      left right (zeroDim + nonnegativeDim) socSizes hleft.lorentz hright.lorentz)

theorem productConeDot_eq_segmentDot
    (left right : Nat → ℝ) (zeroDim nonnegativeDim : Nat)
    (socSizes : List Nat) :
    productConeDot left right zeroDim nonnegativeDim socSizes =
      segmentDot left right 0
        (zeroDim + nonnegativeDim + socSizes.sum) := by
  have hprefix :
      segmentDot left right 0 (zeroDim + nonnegativeDim) =
        segmentDot left right 0 zeroDim +
          segmentDot left right zeroDim nonnegativeDim := by
    simpa using segmentDot_add left right 0 zeroDim nonnegativeDim
  have htotal :
      segmentDot left right 0
          (zeroDim + nonnegativeDim + socSizes.sum) =
        segmentDot left right 0 (zeroDim + nonnegativeDim) +
          segmentDot left right (zeroDim + nonnegativeDim) socSizes.sum := by
    simpa using
      segmentDot_add left right 0 (zeroDim + nonnegativeDim) socSizes.sum
  unfold productConeDot
  rw [lorentzBlocksDot_eq_segmentDot]
  calc
    segmentDot left right 0 zeroDim +
          segmentDot left right zeroDim nonnegativeDim +
        segmentDot left right (zeroDim + nonnegativeDim) socSizes.sum =
      segmentDot left right 0 (zeroDim + nonnegativeDim) +
        segmentDot left right (zeroDim + nonnegativeDim) socSizes.sum := by
          rw [hprefix]
    _ = segmentDot left right 0
          (zeroDim + nonnegativeDim + socSizes.sum) := htotal.symm

def extendFin {n : Nat} (values : Fin n → ℝ) (index : Nat) : ℝ :=
  if h : index < n then values ⟨index, h⟩ else 0

theorem segmentDot_extendFin
    {n : Nat} (left right : Fin n → ℝ) :
    segmentDot (extendFin left) (extendFin right) 0 n =
      dot left right := by
  unfold segmentDot dot
  simp only [Nat.zero_add]
  rw [← Fin.sum_univ_eq_sum_range]
  simp [extendFin]

theorem fin_productCone_pair_nonnegative
    {rows : Nat} (slack dual : Fin rows → ℝ)
    (zeroDim nonnegativeDim : Nat) (socSizes : List Nat)
    (hdimensions :
      zeroDim + nonnegativeDim + socSizes.sum = rows)
    (hslack :
      ProductConePoint (extendFin slack)
        zeroDim nonnegativeDim socSizes)
    (hdual :
      ProductConeDual (extendFin dual)
        zeroDim nonnegativeDim socSizes) :
    0 ≤ dot dual slack := by
  have hpair := productCone_pair_nonnegative
    (extendFin slack) (extendFin dual)
    zeroDim nonnegativeDim socSizes hslack hdual
  rw [productConeDot_eq_segmentDot, hdimensions,
    segmentDot_extendFin] at hpair
  simpa [dot, mul_comm] using hpair

theorem canonical_productCone_weak_duality
    {rows variableCount : Nat}
    (matrix : Fin rows → Fin variableCount → ℝ)
    (right : Fin rows → ℝ)
    (objective : Fin variableCount → ℝ)
    (point : Fin variableCount → ℝ)
    (slack dual : Fin rows → ℝ)
    (zeroDim nonnegativeDim : Nat) (socSizes : List Nat)
    (feasible : CanonicalFeasible matrix right point slack)
    (stationary : ∀ column,
      (∑ row, matrix row column * dual row) + objective column = 0)
    (hdimensions :
      zeroDim + nonnegativeDim + socSizes.sum = rows)
    (hslack :
      ProductConePoint (extendFin slack)
        zeroDim nonnegativeDim socSizes)
    (hdual :
      ProductConeDual (extendFin dual)
        zeroDim nonnegativeDim socSizes) :
    -(∑ column, objective column * point column) ≤ dot right dual := by
  apply canonical_weak_duality
    matrix right objective point slack dual feasible stationary
  exact fin_productCone_pair_nonnegative
    slack dual zeroDim nonnegativeDim socSizes hdimensions hslack hdual

def rationalFunctionAsReal (values : Nat → ℚ) (index : Nat) : ℝ :=
  (values index : ℝ)

theorem rationalLorentzBlock_to_real
    (dual : Nat → ℚ) (start size : Nat)
    (h : rationalLorentzBlock dual start size) :
    lorentzBlock (rationalFunctionAsReal dual) start size := by
  rcases h with ⟨hsize, htime, hsquare⟩
  refine ⟨hsize, ?_, ?_⟩
  · change (0 : ℝ) ≤ (dual start : ℝ)
    exact_mod_cast htime
  · unfold rationalFunctionAsReal
    unfold rationalSpatialSquare at hsquare
    have hsquare' :
        (∑ offset ∈ Finset.range (size - 1),
          (dual (start + offset + 1)) ^ 2) ≤
          (dual start) ^ 2 := by
      simpa [pow_two] using hsquare
    have hcast :
        ((∑ offset ∈ Finset.range (size - 1),
          (dual (start + offset + 1)) ^ 2 : ℚ) : ℝ) ≤
          (((dual start) ^ 2 : ℚ) : ℝ) :=
      (Rat.cast_le (K := ℝ)).2 hsquare'
    simpa using hcast

theorem rationalLorentzBlocks_to_real
    (dual : Nat → ℚ) (start : Nat) (sizes : List Nat)
    (h : rationalLorentzBlocks dual start sizes) :
    lorentzBlocks (rationalFunctionAsReal dual) start sizes := by
  induction sizes generalizing start with
  | nil => trivial
  | cons size rest ih =>
      rcases h with ⟨hhead, htail⟩
      exact ⟨rationalLorentzBlock_to_real dual start size hhead,
        ih (start + size) htail⟩

theorem rationalNonnegativeSegment_to_real
    (dual : Nat → ℚ) (start size : Nat)
    (h : rationalNonnegativeSegment dual start size) :
    nonnegativeSegment (rationalFunctionAsReal dual) start size := by
  intro offset hoffset
  change (0 : ℝ) ≤ (dual (start + offset) : ℝ)
  exact_mod_cast h offset hoffset

theorem rationalProductConeDual_to_real
    (dual : Nat → ℚ) (zeroDim nonnegativeDim : Nat)
    (socSizes : List Nat)
    (h : rationalProductConeDual dual
      zeroDim nonnegativeDim socSizes) :
    ProductConeDual (rationalFunctionAsReal dual)
      zeroDim nonnegativeDim socSizes := by
  exact ⟨rationalNonnegativeSegment_to_real
      dual zeroDim nonnegativeDim h.1,
    rationalLorentzBlocks_to_real dual
      (zeroDim + nonnegativeDim) socSizes h.2⟩

theorem checked_dualCone_to_real
    (data : CertificateData) (h : checkCertificate data = true) :
    ProductConeDual (rationalFunctionAsReal (rationalDualAt data.dual))
      data.zeroDim data.nonnegativeDim data.socSizes.toList := by
  exact rationalProductConeDual_to_real
    (rationalDualAt data.dual) data.zeroDim data.nonnegativeDim data.socSizes.toList
    (dualConeOK_sound data (checked_dualConeOK data h))

theorem checked_dimensions_eq
    (data : CertificateData) (h : checkCertificate data = true) :
    data.zeroDim + data.nonnegativeDim + data.socSizes.toList.sum =
      data.rows := by
  have hdimensions := checked_dimensionsOK data h
  simpa [dimensionsOK] using hdimensions

theorem checked_productCone_pair_nonnegative
    (data : CertificateData) (slack : Nat → ℝ)
    (h : checkCertificate data = true)
    (hslack : ProductConePoint slack
      data.zeroDim data.nonnegativeDim data.socSizes.toList) :
    0 ≤ productConeDot slack (rationalFunctionAsReal (rationalDualAt data.dual))
      data.zeroDim data.nonnegativeDim data.socSizes.toList := by
  exact productCone_pair_nonnegative
    slack (rationalFunctionAsReal (rationalDualAt data.dual))
    data.zeroDim data.nonnegativeDim data.socSizes.toList
    hslack (checked_dualCone_to_real data h)

theorem checked_full_pair_nonnegative
    (data : CertificateData) (slack : Nat → ℝ)
    (h : checkCertificate data = true)
    (hslack : ProductConePoint slack
      data.zeroDim data.nonnegativeDim data.socSizes.toList) :
    0 ≤ segmentDot slack (rationalFunctionAsReal (rationalDualAt data.dual)) 0 data.rows := by
  have hpair := checked_productCone_pair_nonnegative data slack h hslack
  rw [productConeDot_eq_segmentDot, checked_dimensions_eq data h] at hpair
  exact hpair

def certificateDual (data : CertificateData) (row : Fin data.rows) : ℝ :=
  (data.dual row : ℝ)

theorem segmentDot_extendFin_certificateDual
    (data : CertificateData) (slack : Fin data.rows → ℝ) :
    segmentDot (extendFin slack) (rationalFunctionAsReal (rationalDualAt data.dual))
        0 data.rows =
      dot slack (certificateDual data) := by
  unfold segmentDot dot certificateDual rationalFunctionAsReal rationalDualAt
  simp only [Nat.zero_add]
  rw [← Fin.sum_univ_eq_sum_range]
  simp [extendFin]

theorem checked_fin_pair_nonnegative
    (data : CertificateData) (slack : Fin data.rows → ℝ)
    (h : checkCertificate data = true)
    (hslack : ProductConePoint (extendFin slack)
      data.zeroDim data.nonnegativeDim data.socSizes.toList) :
    0 ≤ dot (certificateDual data) slack := by
  have hpair := checked_full_pair_nonnegative
    data (extendFin slack) h hslack
  rw [segmentDot_extendFin_certificateDual] at hpair
  simpa [dot, mul_comm] using hpair

theorem checked_dual_weak_duality
    (data : CertificateData) {variableCount : Nat}
    (matrix : Fin data.rows → Fin variableCount → ℝ)
    (right : Fin data.rows → ℝ)
    (objective : Fin variableCount → ℝ)
    (point : Fin variableCount → ℝ)
    (slack : Fin data.rows → ℝ)
    (h : checkCertificate data = true)
    (feasible : CanonicalFeasible matrix right point slack)
    (stationary : ∀ column,
      (∑ row, matrix row column * certificateDual data row) +
        objective column = 0)
    (hslack : ProductConePoint (extendFin slack)
      data.zeroDim data.nonnegativeDim data.socSizes.toList) :
    -(∑ column, objective column * point column) ≤
      dot right (certificateDual data) := by
  apply canonical_weak_duality
    matrix right objective point slack (certificateDual data)
    feasible stationary
  exact checked_fin_pair_nonnegative data slack h hslack

theorem certificateProof_dualCone_to_real
    (data : CertificateData) (proof : CertificateProof data) :
    ProductConeDual (rationalFunctionAsReal (rationalDualAt data.dual))
      data.zeroDim data.nonnegativeDim data.socSizes.toList := by
  exact rationalProductConeDual_to_real
    (rationalDualAt data.dual) data.zeroDim data.nonnegativeDim data.socSizes.toList
    proof.dualCone

theorem certificateProof_dimensions_eq
    (data : CertificateData) (proof : CertificateProof data) :
    data.zeroDim + data.nonnegativeDim + data.socSizes.toList.sum =
      data.rows := by
  simpa [dimensionsOK] using proof.dimensions

theorem certificateProof_productCone_pair_nonnegative
    (data : CertificateData) (proof : CertificateProof data)
    (slack : Nat → ℝ)
    (hslack : ProductConePoint slack
      data.zeroDim data.nonnegativeDim data.socSizes.toList) :
    0 ≤ productConeDot slack (rationalFunctionAsReal (rationalDualAt data.dual))
      data.zeroDim data.nonnegativeDim data.socSizes.toList := by
  exact productCone_pair_nonnegative
    slack (rationalFunctionAsReal (rationalDualAt data.dual))
    data.zeroDim data.nonnegativeDim data.socSizes.toList
    hslack (certificateProof_dualCone_to_real data proof)

theorem certificateProof_full_pair_nonnegative
    (data : CertificateData) (proof : CertificateProof data)
    (slack : Nat → ℝ)
    (hslack : ProductConePoint slack
      data.zeroDim data.nonnegativeDim data.socSizes.toList) :
    0 ≤ segmentDot slack (rationalFunctionAsReal (rationalDualAt data.dual)) 0 data.rows := by
  have hpair :=
    certificateProof_productCone_pair_nonnegative data proof slack hslack
  rw [productConeDot_eq_segmentDot,
    certificateProof_dimensions_eq data proof] at hpair
  exact hpair

theorem certificateProof_fin_pair_nonnegative
    (data : CertificateData) (proof : CertificateProof data)
    (slack : Fin data.rows → ℝ)
    (hslack : ProductConePoint (extendFin slack)
      data.zeroDim data.nonnegativeDim data.socSizes.toList) :
    0 ≤ dot (certificateDual data) slack := by
  have hpair := certificateProof_full_pair_nonnegative
    data proof (extendFin slack) hslack
  rw [segmentDot_extendFin_certificateDual] at hpair
  simpa [dot, mul_comm] using hpair

theorem certificateProof_weak_duality
    (data : CertificateData) (proof : CertificateProof data)
    {variableCount : Nat}
    (matrix : Fin data.rows → Fin variableCount → ℝ)
    (right : Fin data.rows → ℝ)
    (objective : Fin variableCount → ℝ)
    (point : Fin variableCount → ℝ)
    (slack : Fin data.rows → ℝ)
    (feasible : CanonicalFeasible matrix right point slack)
    (stationary : ∀ column,
      (∑ row, matrix row column * certificateDual data row) +
        objective column = 0)
    (hslack : ProductConePoint (extendFin slack)
      data.zeroDim data.nonnegativeDim data.socSizes.toList) :
    -(∑ column, objective column * point column) ≤
      dot right (certificateDual data) := by
  apply canonical_weak_duality
    matrix right objective point slack (certificateDual data)
    feasible stationary
  exact certificateProof_fin_pair_nonnegative data proof slack hslack

theorem certificateProof_strict_upper
    (data : CertificateData) (proof : CertificateProof data)
    {variableCount : Nat}
    (matrix : Fin data.rows → Fin variableCount → ℝ)
    (right : Fin data.rows → ℝ)
    (objective : Fin variableCount → ℝ)
    (point : Fin variableCount → ℝ)
    (slack : Fin data.rows → ℝ)
    (feasible : CanonicalFeasible matrix right point slack)
    (stationary : ∀ column,
      (∑ row, matrix row column * certificateDual data row) +
        objective column = 0)
    (hslack : ProductConePoint (extendFin slack)
      data.zeroDim data.nonnegativeDim data.socSizes.toList)
    (hright :
      dot right (certificateDual data) =
        (certifiedUpper data : ℝ)) :
    -(∑ column, objective column * point column) <
      (data.target : ℝ) := by
  have hweak := certificateProof_weak_duality
    data proof matrix right objective point slack feasible stationary hslack
  rw [hright] at hweak
  have hupper : (certifiedUpper data : ℝ) < (data.target : ℝ) := by
    exact_mod_cast proof.upper
  exact hweak.trans_lt hupper

end CarmenQExact

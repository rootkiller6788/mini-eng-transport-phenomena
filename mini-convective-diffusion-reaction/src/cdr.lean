/-
  Formalization: Convection-Diffusion-Reaction Transport Phenomena
  Reference: Bird, Stewart, Lightfoot "Transport Phenomena" (2007)

  This file formalizes key conservation laws and dimensionless number
  relationships for CDR systems in Lean 4.

  L4: Conservation laws formalized
  L1: Dimensionless number definitions
  L2: Transport regime relationships

  NOTE: Float arithmetic proof is restricted in Lean 4 (Float is not Ring).
  Theorems use structural/inductive reasoning or specific constant evaluation.
-/

-- ===========================================================================
-- L1: Dimensionless Number Definitions
-- ===========================================================================

/- Peclet Number: Pe = u * L / D -/
def PecletNumber : Type := Float

def peclet (u L D : Float) : PecletNumber :=
  if D == 0.0 then 0.0 else u * L / D

/- Damkohler Number (Type I): Da_I = k * tau -/
def DamkohlerI : Type := Float

def damkohlerI (k tau : Float) : DamkohlerI := k * tau

/- Damkohler Number (Type II): Da_II = k * L^2 / D -/
def DamkohlerII : Type := Float

def damkohlerII (k L D : Float) : DamkohlerII :=
  if D == 0.0 then 0.0 else k * L * L / D

/- Thiele Modulus: phi = L_c * sqrt(k / D_eff) -/
def ThieleModulus : Type := Float

def thiele (Lc k Deff : Float) : ThieleModulus :=
  if Deff <= 0.0 || Lc <= 0.0 then 0.0
  else Lc * Float.sqrt (k / Deff)

/- Schmidt Number: Sc = mu / (rho * D) -/
def SchmidtNumber : Type := Float

def schmidt (mu rho D : Float) : SchmidtNumber :=
  if rho == 0.0 || D == 0.0 then 0.0 else mu / (rho * D)

/- Reynolds Number: Re = rho * u * L / mu -/
def ReynoldsNumber : Type := Float

def reynolds (rho u L mu : Float) : ReynoldsNumber :=
  if mu == 0.0 then 0.0 else rho * u * L / mu

/- Sherwood Number: Sh = k_c * L / D -/
def SherwoodNumber : Type := Float

def sherwood (kc L D : Float) : SherwoodNumber :=
  if D == 0.0 then 0.0 else kc * L / D

-- ===========================================================================
-- L2: Transport Regime Classification
-- ===========================================================================

inductive TransportRegime
  | diffusionDominated
  | convectionDominated
  | reactionDominated
  | mixedRegime
  | dispersionRegime
  | massTransferLimited
  deriving BEq, Repr, Inhabited

theorem diffusion_dominated_distinct_from_convection :
  TransportRegime.diffusionDominated != TransportRegime.convectionDominated := by
  intro h
  injection h

theorem regimes_distinct :
  List.Chain' (. != .) [
    TransportRegime.diffusionDominated,
    TransportRegime.convectionDominated,
    TransportRegime.reactionDominated,
    TransportRegime.mixedRegime,
    TransportRegime.dispersionRegime,
    TransportRegime.massTransferLimited
  ] := by
  decide

-- ===========================================================================
-- L4: Conservation Laws
-- ===========================================================================

/- Ficks First Law (1D): J = -D * dC/dx -/
def FickFirstLaw (D dCdx : Float) : Float := -D * dCdx

theorem fick_zero_gradient (D : Float) : FickFirstLaw D 0.0 = 0.0 := by
  unfold FickFirstLaw
  simp

/- First-Order Reaction Rate Law: r = k * C_A -/
def firstOrderRate (k C_A : Float) : Float := k * C_A

theorem first_order_zero_rate (C_A : Float) : firstOrderRate 0.0 C_A = 0.0 := by
  unfold firstOrderRate
  simp

/- Arrhenius Law: k(T) = A * exp(-Ea / (R * T)) where R = 8.314 J/(mol*K) -/
def arrheniusRate (A Ea T : Float) : Float :=
  if T <= 0.0 then 0.0
  else
    let R : Float := 8.314
    A * Float.exp (-Ea / (R * T))

theorem arrhenius_zero_activation (A T : Float) (hT : T > 0) :
    arrheniusRate A 0.0 T = A := by
  unfold arrheniusRate
  simp [hT]

theorem arrhenius_zero_prefactor (Ea T : Float) :
    arrheniusRate 0.0 Ea T = 0.0 := by
  unfold arrheniusRate
  by_cases hT : T <= 0.0
  . simp [hT]
  . simp [hT]

-- ===========================================================================
-- L4: Effectiveness Factor
-- ===========================================================================

/- Effectiveness factor for slab geometry (1st order):
  eta = tanh(phi) / phi for phi > 0, eta = 1 for phi = 0 -/
def effectivenessFactorSlab (phi : Float) : Float :=
  if phi <= 0.0 then 1.0
  else Float.tanh phi / phi

theorem effectiveness_no_limitation : effectivenessFactorSlab 0.0 = 1.0 := by
  unfold effectivenessFactorSlab
  simp

theorem effectiveness_bounded_at_zero : effectivenessFactorSlab 0.0 <= 1.0 := by
  rw [effectiveness_no_limitation]
  rfl

-- ===========================================================================
-- L5: Reactor Design Equations
-- ===========================================================================

/- CSTR Design Equation (1st-order): C_A = C_A0 / (1 + k * tau) -/
def cstrOutletConcentration (C_A0 k tau : Float) : Float :=
  C_A0 / (1.0 + k * tau)

theorem cstr_no_reaction (C_A0 tau : Float) :
    cstrOutletConcentration C_A0 0.0 tau = C_A0 := by
  unfold cstrOutletConcentration
  simp

theorem cstr_zero_tau (C_A0 k : Float) :
    cstrOutletConcentration C_A0 k 0.0 = C_A0 := by
  unfold cstrOutletConcentration
  simp

/- PFR Design Equation (1st-order): C_A(L) = C_A0 * exp(-k * tau) -/
def pfrOutletConcentration (C_A0 k tau : Float) : Float :=
  C_A0 * Float.exp (-k * tau)

theorem pfr_no_reaction (C_A0 tau : Float) :
    pfrOutletConcentration C_A0 0.0 tau = C_A0 := by
  unfold pfrOutletConcentration
  simp

theorem pfr_differs_from_cstr_specific :
    pfrOutletConcentration 1.0 1.0 1.0 !=
    cstrOutletConcentration 1.0 1.0 1.0 := by
  unfold pfrOutletConcentration cstrOutletConcentration
  native_decide

/- Conversion: X = 1 - C_A / C_A0 -/
def conversion (C_A C_A0 : Float) : Float :=
  if C_A0 == 0.0 then 0.0 else 1.0 - C_A / C_A0

theorem conversion_no_reaction (C_A0 : Float) (h : C_A0 != 0.0) :
    conversion C_A0 C_A0 = 0.0 := by
  unfold conversion
  simp [h]

theorem conversion_complete (C_A0 : Float) (h : C_A0 != 0.0) :
    conversion 0.0 C_A0 = 1.0 := by
  unfold conversion
  simp [h]

-- ===========================================================================
-- L6: Reactor Volume Ratio
-- ===========================================================================

/- V_CSTR / V_PFR for first-order reaction
  ratio = X / ((1-X) * ln(1/(1-X))) -/
def cstrPfrVolumeRatio (X : Float) : Float :=
  if X <= 0.0 || X >= 1.0 then 1.0
  else
    let one_minus_X := 1.0 - X
    X / (one_minus_X * Float.log (1.0 / one_minus_X))

theorem volume_ratio_above_one_at_half :
    cstrPfrVolumeRatio 0.5 > 1.0 := by
  unfold cstrPfrVolumeRatio
  native_decide

-- ===========================================================================
-- L5: Two-Film Theory
-- ===========================================================================

/- Overall gas-side mass transfer coefficient: 1/K_OG = 1/k_G + H/(k_L * P) -/
def overallMassTransferCoeff (k_G k_L H P : Float) : Float :=
  if k_G <= 0.0 || k_L <= 0.0 || P <= 0.0 then 0.0
  else 1.0 / (1.0 / k_G + H / (k_L * P))

theorem overall_approx_gas_side (k_G H P : Float) (hG : k_G > 0) (hP : P > 0) :
    overallMassTransferCoeff k_G 1e10 H P = k_G := by
  unfold overallMassTransferCoeff
  have hkL : (1e10 : Float) > 0.0 := by native_decide
  simp [hG, hkL, hP]
  native_decide

theorem overall_zero_resistances (H P : Float) :
    overallMassTransferCoeff 0.0 0.0 H P = 0.0 := by
  unfold overallMassTransferCoeff
  simp

-- ===========================================================================
-- L3: Engineering Quantities -- Numerical Ranges
-- ===========================================================================

structure QuantityRange where
  min : Float
  max : Float
  unit : String
  description : String

def gasDiffusivityRange : QuantityRange :=
  { min := 1.0e-5, max := 1.0e-4, unit := "m2/s",
    description := "Typical binary gas diffusion coefficient at STP" }

def liquidDiffusivityRange : QuantityRange :=
  { min := 1.0e-10, max := 1.0e-9, unit := "m2/s",
    description := "Typical liquid-phase diffusion coefficient" }

def gasMassTransferCoeffRange : QuantityRange :=
  { min := 1.0e-3, max := 1.0e-1, unit := "m/s",
    description := "Typical gas-phase mass transfer coefficient" }

def inRange (v : Float) (r : QuantityRange) : Bool :=
  r.min <= v && v <= r.max

theorem inRange_min_true (r : QuantityRange) : inRange r.min r := by
  unfold inRange
  simp

theorem inRange_max_true (r : QuantityRange) : inRange r.max r := by
  unfold inRange
  simp

-- ===========================================================================
-- L5: Half-Life Analysis
-- ===========================================================================

/- Half-life for first-order reaction: t_1/2 = ln(2) / k -/
def halfLifeFirstOrder (k : Float) : Float :=
  if k <= 0.0 then Float.infinity
  else Float.log 2.0 / k

theorem half_life_zero_rate : halfLifeFirstOrder 0.0 = Float.infinity := by
  unfold halfLifeFirstOrder
  simp

/-
  Theorem: For k = ln(2), half-life is exactly 1.
  Uses Float.div which simplifies since Float.log 2.0 / Float.log 2.0 = 1.0
  when the denominator is nonzero.
-/
theorem half_life_ln2 : halfLifeFirstOrder (Float.log 2.0) = 1.0 := by
  unfold halfLifeFirstOrder
  have hpos : Float.log 2.0 > 0.0 := by
    have h2gt1 : (2.0 : Float) > 1.0 := by native_decide
    exact Float.log_pos h2gt1
  simp [hpos]

-- ===========================================================================
-- L3: Catalyst Geometry Enum
-- ===========================================================================

inductive CatalystGeometry
  | slab
  | cylinder
  | sphere
  | arbitrary3D
  deriving BEq, Repr, Inhabited

theorem catalyst_geometry_cardinality :
  (List.map (fun _ => ()) [
    CatalystGeometry.slab,
    CatalystGeometry.cylinder,
    CatalystGeometry.sphere,
    CatalystGeometry.arbitrary3D
  ]).length = 4 := by
  native_decide

theorem catalyst_geometries_distinct :
  CatalystGeometry.slab != CatalystGeometry.sphere := by
  intro h
  injection h

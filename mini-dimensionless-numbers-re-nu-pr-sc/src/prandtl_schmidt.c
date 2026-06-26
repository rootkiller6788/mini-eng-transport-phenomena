/*
 * prandtl_schmidt.c — Prandtl, Schmidt, and Lewis Number Computations
 *
 * Reference: Bird, Stewart, Lightfoot (2007), Incropera & DeWitt (2007),
 *            Cussler (2009), Wilke-Chang (1955)
 *
 * Knowledge: L1-L5 complete
 */

#include "../include/prandtl_schmidt.h"
#include <math.h>

/* ==========================================================================
 * L1: Prandtl Number
 * ========================================================================== */

double prandtl_from_material(double cp, double mu, double k)
{
    /* Pr = cp.mu/k */
    if (k <= 0.0 || cp < 0.0 || mu < 0.0) return -1.0;
    if (k < 1e-15) return -1.0;
    return cp * mu / k;
}

double prandtl_from_diffusivities(double nu, double alpha)
{
    /* Pr = nu/alpha */
    if (alpha <= 0.0 || nu < 0.0) return -1.0;
    if (alpha < 1e-20) return -1.0;
    return nu / alpha;
}

double prandtl_thermal_diffusivity(double nu, double Pr)
{
    /* alpha = nu/Pr */
    if (Pr <= 0.0 || nu < 0.0) return -1.0;
    return nu / Pr;
}

double prandtl_kinematic_viscosity(double Pr, double alpha)
{
    /* nu = Pr.alpha */
    if (alpha < 0.0 || Pr < 0.0) return -1.0;
    return Pr * alpha;
}

double prandtl_air(double T)
{
    /*
     * Pr(T) for dry air at 1 atm, 200 K < T < 2000 K.
     *
     * Data from Incropera Table A.4:
     *   T [K]  | Pr
     *   250     | 0.724
     *   300     | 0.707
     *   400     | 0.690
     *   500     | 0.684
     *   600     | 0.680
     *   1000    | 0.726
     *   1500    | 0.738
     *   2000    | 0.746
     *
     * For gases, Pr is weakly T-dependent (primarily through cp(T)).
     * The curve is slightly U-shaped with a minimum near 500 K.
     *
     * Piecewise linear fit for engineering use (±3% accuracy):
     */
    if (T < 200.0 || T > 2000.0) return -1.0;

    if (T < 400.0) {
        /* 200-400 K: decreasing */
        return 0.742 - 0.00013 * (T - 200.0);
    } else if (T < 600.0) {
        /* 400-600 K: near minimum */
        return 0.690 - 0.000025 * (T - 400.0);
    } else if (T < 1000.0) {
        /* 600-1000 K: increasing */
        return 0.680 + 0.000115 * (T - 600.0);
    } else {
        /* 1000-2000 K: slowly increasing */
        return 0.726 + 0.00002 * (T - 1000.0);
    }
}

double prandtl_water(double T)
{
    /*
     * Pr(T) for saturated liquid water, 0 degC < T < 200 degC.
     *
     * Data from NIST / Incropera Table A.6:
     *   T [ degC] |  Pr
     *   0       |  13.5
     *   20      |  7.01
     *   40      |  4.34
     *   60      |  3.00
     *   80      |  2.22
     *   100     |  1.76
     *   150     |  1.16
     *   200     |  0.94
     *
     * Water's Pr decreases strongly with T because mu(T) decreases
     * much faster than k(T) and cp(T) change.
     *
     * Fit: ln(Pr) ≈ a + b.ln(T+273) + c.[ln(T+273)]^2
     * For engineering use, piecewise linear in ln(Pr) vs T:
     */
    if (T < 0.0 || T > 200.0) return -1.0;

    /* Simplified fit with exponential decay + constant floor */
    if (T < 75.0) {
        return 13.5 * exp(-0.0223 * T) + 1.8 * (1.0 - exp(-0.0223 * T));
    } else {
        /* High-temperature asymptote approaching Pr ≈ 1 */
        return 2.0 * exp(-0.006 * (T - 75.0)) + 1.0 * (1.0 - exp(-0.006 * (T - 75.0)));
    }
}

double prandtl_engine_oil(double T)
{
    /*
     * Pr(T) for unused engine oil (SAE 10W-30 typical).
     *
     * Characterised by very high Pr at low temperatures due to
     * extremely high viscosity: Pr(0 degC) ~ 10⁴, Pr(100 degC) ~ 200,
     * Pr(150 degC) ~ 80.
     *
     * This strong temperature dependence drives the design of
     * oil coolers in automotive applications. The viscosity
     * decrease with temperature (~ Arrhenius-like) dominates Pr.
     *
     * Simplified fit for 0 < T < 200 degC:
     */
    if (T < 0.0 || T > 200.0) return -1.0;

    /* Exponential decay: Pr ~ exp(-b.T) */
    double Pr_0 = 15000.0;  /* Pr at T ≈ 0 degC */
    double decay_rate = 0.035; /* per  degC */

    return Pr_0 * exp(-decay_rate * T) + 80.0 * (1.0 - exp(-decay_rate * T));
}

/* ==========================================================================
 * L1-L2: Schmidt Number
 * ========================================================================== */

double schmidt_basic(double nu, double D_AB)
{
    /* Sc = nu/D_AB */
    if (D_AB <= 0.0 || nu < 0.0) return -1.0;
    if (D_AB < 1e-20) return -1.0;
    return nu / D_AB;
}

double schmidt_from_material(double mu, double rho, double D_AB)
{
    /* Sc = mu/(rho.D_AB) */
    if (rho <= 0.0 || D_AB <= 0.0 || mu < 0.0) return -1.0;
    if (D_AB < 1e-20 || rho < 1e-20) return -1.0;
    return mu / (rho * D_AB);
}

double schmidt_mass_diffusivity(double nu, double Sc)
{
    /* D_AB = nu/Sc */
    if (Sc <= 0.0 || nu < 0.0) return -1.0;
    return nu / Sc;
}

double schmidt_gas_typical(void)
{
    /*
     * Typical Schmidt number for gases at 1 atm, 300 K.
     *
     * For most gas pairs: Sc ≈ 0.5-2.0.
     *   CO2 in air:  Sc ≈ 1.0
     *   H2O in air:  Sc ≈ 0.6
     *   O2 in N2:    Sc ≈ 0.74
     *   He in air:   Sc ≈ 0.2
     *
     * Since nu ~ D_AB for gases (both ~ T^(3/2)/p via Chapman-Enskog),
     * Sc is typically O(1).
     */
    return 1.0;
}

double schmidt_liquid_typical(int solute_type)
{
    /*
     * Typical Sc for solutes in liquids at 300 K.
     *
     * In liquids, nu is much larger than D_AB (by 10^2-10⁴×):
     *   Small molecules in water (O2, CO2):  Sc ≈ 200-500
     *   Large molecules/proteins in water:   Sc ≈ 10^3-10⁵
     *   Polymers in organic solvents:        Sc ≈ 10⁴-10⁶
     *
     * The large Sc means the concentration boundary layer is
     * much thinner than the velocity boundary layer:
     *   delta_C/delta ≈ Sc^(-1/3)
     *
     * For Sc = 1000: delta_C ≈ delta/10 (thin mass transfer BL).
     */
    switch (solute_type) {
        case 0: return 300.0;    /* small molecule in water */
        case 1: return 5000.0;   /* protein in water */
        case 2: return 50000.0;  /* polymer in organic solvent */
        default: return -1.0;
    }
}

/* ==========================================================================
 * L1-L2: Lewis Number
 * ========================================================================== */

double lewis_basic(double alpha, double D_AB)
{
    /* Le = alpha/D_AB */
    if (D_AB <= 0.0 || alpha < 0.0) return -1.0;
    if (D_AB < 1e-20) return -1.0;
    return alpha / D_AB;
}

double lewis_from_sc_pr(double Sc, double Pr)
{
    /* Le = Sc/Pr */
    if (Pr <= 0.0 || Sc < 0.0) return -1.0;
    return Sc / Pr;
}

double lewis_typical(int category)
{
    /*
     * Typical Lewis numbers:
     *
     * Gases:       Le ≈ Sc/Pr ≈ 1.0/0.7 ≈ 1.4
     *   Heat and mass diffuse at similar rates.
     *   The psychrometric ratio (h/h_m)/cp ≈ Le^(2/3) ≈ 1.
     *
     * Water:       Le ≈ 500/7 ≈ 70
     *   Heat diffuses much faster than mass.
     *   Temperature field develops much faster than concentration.
     *
     * Liquid metals: Le ≈ 1/0.01 ≈ 100
     *   Pr << 1 (thermal diffuses fast), Sc ~ 1 (mass diffuses slow).
     *   Very different from gases!
     */
    switch (category) {
        case 0: return 1.4;    /* gas mixture (air-water vapor) */
        case 1: return 70.0;   /* aqueous solution */
        case 2: return 100.0;  /* liquid metal */
        default: return -1.0;
    }
}

/* ==========================================================================
 * L3: Diffusivity Models
 * ========================================================================== */

double thermal_diffusivity_definition(double k, double rho, double cp)
{
    /* alpha = k/(rho.cp) */
    if (rho <= 0.0 || cp <= 0.0 || k < 0.0) return -1.0;
    return k / (rho * cp);
}

double mass_diffusivity_gas_estimate(double T, double p, double M_A, double M_B)
{
    /*
     * Binary gas diffusivity — simplified Chapman-Enskog form.
     *
     * D_AB = (1.86×10⁻⁷).T^(1.75) / p  .  (1/M_A + 1/M_B)^(1/2)
     *
     * Units: D in [m^2/s], T in [K], p in [Pa], M in [kg/kmol].
     *
     * For air (M ≈ 29 kg/kmol) at 300 K, 1 atm (101325 Pa):
     *   D_AB ≈ 2.0×10⁻⁵ m^2/s (typical gas diffusivity)
     *
     * Temperature dependence: D ~ T^1.75 (faster than nu ~ T^1.5),
     * so Sc increases slightly with T for gases.
     *
     * Pressure dependence: D ~ 1/p. At 10 atm, D_AB ~ 2×10⁻⁶ m^2/s.
     */
    if (T <= 0.0 || p <= 0.0 || M_A <= 0.0 || M_B <= 0.0) return -1.0;

    double reduced_mass_factor = sqrt(1.0 / M_A + 1.0 / M_B);
    return 1.86e-7 * pow(T, 1.75) / p * reduced_mass_factor;
}

double mass_diffusivity_liquid_wilke_chang(double T, double mu_B, double V_A,
                                           double M_B, double phi_B)
{
    /*
     * Wilke-Chang (1955) correlation for liquid diffusivity.
     *
     * D_AB = (7.4×10⁻¹^2).(φ_B.M_B)^(1/2).T / (mu_B.V_A^0.6)
     *
     * Units: T [K], mu_B [cP], V_A [cm^3/mol], M_B [g/mol]
     *        Returns D_AB [m^2/s].
     *
     * Typical values for V_A (molar volume at boiling point):
     *   O2: 25.6, N2: 31.2, CO2: 34.0, H2O: 18.9, benzene: 96.0
     *
     * Association parameter φ_B:
     *   Water: 2.6, Methanol: 1.9, Ethanol: 1.5,
     *   Benzene, heptane, ether: 1.0
     *
     * For a small solute in water at 300 K:
     *   mu_B ≈ 0.86 cP, V_A ≈ 30 cm^3/mol, M_B = 18 g/mol, φ_B = 2.6
     *   -> D_AB ≈ 1.2×10⁻⁹ m^2/s
     *
     * Compare: D for gases ~ 2×10⁻⁵ m^2/s — liquids are ~10⁴× slower!
     */
    if (T <= 0.0 || mu_B <= 0.0 || V_A <= 0.0 || M_B <= 0.0 || phi_B <= 0.0)
        return -1.0;

    double numerator = 7.4e-12 * sqrt(phi_B * M_B) * T;
    double denominator = mu_B * pow(V_A, 0.6);

    if (denominator < 1e-20) return -1.0;
    return numerator / denominator;
}

/* ==========================================================================
 * L2: Boundary Layer Ratios
 * ========================================================================== */

double thermal_to_velocity_bl_ratio(double Pr)
{
    /*
     * delta_T/delta ≈ Pr^(-1/3) for Pr > 0.5 (Pohlhausen similarity)
     *
     * For Pr -> 0 (liquid metals):
     *   delta_T/delta ≈ 1.325.Pr^(-1/2)
     *
     * Physical interpretation:
     *   Pr << 1: Heat diffuses much faster than momentum
     *            -> thermal BL is thicker (delta_T >> delta)
     *            Example: mercury fills the whole pipe
     *
     *   Pr ~ 1: Similar thicknesses (gases)
     *
     *   Pr >> 1: Momentum diffuses much faster than heat
     *            -> velocity BL is thicker (delta >> delta_T)
     *            Example: viscous oil, delta_T confined to thin layer near wall
     */
    if (Pr <= 0.0) return -1.0;

    if (Pr < 0.05) {
        return 1.325 / sqrt(Pr);
    }
    return 1.0 / cbrt(Pr);
}

double concentration_to_velocity_bl_ratio(double Sc)
{
    /*
     * delta_C/delta ≈ Sc^(-1/3) for Sc > 0.5
     *
     * Same scaling as thermal BL (mathematically identical).
     * Mass and heat transfer analogous when replacing Pr -> Sc, Nu -> Sh.
     */
    if (Sc <= 0.0) return -1.0;
    if (Sc < 0.05) {
        return 1.325 / sqrt(Sc);
    }
    return 1.0 / cbrt(Sc);
}

int triple_analogy_bl_ratios(double Pr, double Sc,
                             double *ratio_T, double *ratio_C)
{
    if (!ratio_T || !ratio_C) return -1;
    if (Pr <= 0.0 || Sc <= 0.0) return -1;

    *ratio_T = thermal_to_velocity_bl_ratio(Pr);
    *ratio_C = concentration_to_velocity_bl_ratio(Sc);
    return 0;
}

/* ==========================================================================
 * L5: Combined Transport Analogies
 * ========================================================================== */

int colburn_analogy_verify(double St_H, double Pr, double St_M, double Sc,
                           double f, double *j_H, double *j_D)
{
    /*
     * Chilton-Colburn analogy verification:
     *   j_H = St_H × Pr^(2/3)
     *   j_D = St_M × Sc^(2/3)
     *
     * The analogy predicts: j_H ≈ j_D ≈ f/8.
     *
     * Engineering significance:
     *   If you know the friction factor (easy to measure with DeltaP),
     *   you can predict heat transfer (harder to measure) AND
     *   mass transfer (even harder to measure) for the same geometry.
     */
    if (!j_H || !j_D) return -1;
    if (St_H < 0.0 || Pr <= 0.0 || St_M < 0.0 || Sc <= 0.0 || f <= 0.0)
        return -1;

    *j_H = St_H * pow(Pr, 2.0 / 3.0);
    *j_D = St_M * pow(Sc, 2.0 / 3.0);
    return 0;
}

double reynolds_analogy_stanton(double f)
{
    /*
     * Simplest Reynolds analogy (1874): St = f/8.
     *
     * Valid only for Pr ≈ 1 (gases) because it assumes
     * identical turbulent transport of momentum and heat.
     *
     * The factor 8 comes from the Darcy friction factor convention.
     * In terms of Fanning friction factor Cf = f/4: St = Cf/2.
     */
    if (f <= 0.0) return -1.0;
    return f / 8.0;
}

double reynolds_analogy_nusselt(double f, double Re, double Pr)
{
    /*
     * From St = f/8 and St = Nu/(Re.Pr):
     *   Nu = (f/8).Re.Pr
     *
     * Good approximation for turbulent gas flows.
     */
    if (f <= 0.0 || Re <= 0.0 || Pr <= 0.0) return -1.0;
    return (f / 8.0) * Re * Pr;
}

double chilton_colburn_sherwood(double f, double Re, double Sc)
{
    /*
     * Chilton-Colburn analogy for mass transfer:
     *   j_D = f/8  ->  St_M × Sc^(2/3) = f/8
     *   St_M = Sh/(Re.Sc)
     *   -> Sh = (f/8).Re.Sc^(1/3)
     *
     * Valid for: 0.6 < Sc < 3000, turbulent flow.
     *
     * This is the mass transfer analog of the heat transfer Colburn
     * analogy. Widely used for estimating mass transfer coefficients
     * in packed beds, wetted-wall columns, and heat exchangers.
     */
    if (f <= 0.0 || Re <= 0.0 || Sc <= 0.0) return -1.0;
    return (f / 8.0) * Re * pow(Sc, 1.0 / 3.0);
}

double heat_mass_analogy_nusselt_to_sherwood(double Nu, double Pr, double Sc)
{
    /*
     * Sh = Nu.(Sc/Pr)^(1/3)
     *
     * This is the direct heat-to-mass transfer analogy:
     * same geometry, same flow -> same functional form.
     *
     * If Nu = C.Re^m.Pr^(1/3), then Sh = C.Re^m.Sc^(1/3).
     * Therefore: Sh/Nu = (Sc/Pr)^(1/3) = Le^(1/3).
     *
     * This assumes the same C and m for heat and mass transfer
     * (i.e., same flow configuration with heat and mass analogies).
     */
    if (Nu < 0.0 || Pr <= 0.0 || Sc <= 0.0) return -1.0;
    return Nu * pow(Sc / Pr, 1.0 / 3.0);
}

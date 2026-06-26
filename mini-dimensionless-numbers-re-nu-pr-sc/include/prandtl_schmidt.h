/*
 * prandtl_schmidt.h — Prandtl, Schmidt, and Lewis Numbers
 *
 * Reference: Bird, Stewart, Lightfoot (2007) "Transport Phenomena", Ch. 9-10, 17-24
 *            Incropera & DeWitt (2007) "Fundamentals of Heat and Mass Transfer"
 *            Cussler (2009) "Diffusion: Mass Transfer in Fluid Systems"
 *
 * These numbers relate the three transport diffusivities:
 *   nu  (momentum diffusivity, m^2/s)
 *   alpha  (thermal diffusivity, m^2/s)
 *   D  (mass diffusivity, m^2/s)
 *
 * Pr = nu/alpha  — Prandtl number: momentum vs. thermal diffusion
 * Sc = nu/D  — Schmidt number:  momentum vs. mass diffusion
 * Le = alpha/D  — Lewis number:    thermal vs. mass diffusion
 *
 * Knowledge Coverage: L1 Definitions, L2 Transport Analogy,
 *                     L3 Property Ranges, L4 Species Conservation,
 *                     L5 Correlation Methods
 */

#ifndef PRANDTL_SCHMIDT_H
#define PRANDTL_SCHMIDT_H

#ifdef __cplusplus
extern "C" {
#endif

/* ==========================================================================
 * L1: Prandtl Number — Momentum diffusivity / Thermal diffusivity
 * ========================================================================== */

/**
 * prandtl_from_material — Pr = cp.mu/k
 *
 * @param cp  specific heat at constant pressure [J/(kg.K)]
 * @param mu  dynamic viscosity [Pa.s]
 * @param k   thermal conductivity [W/(m.K)]
 * @return    Prandtl number [-]
 */
double prandtl_from_material(double cp, double mu, double k);

/**
 * prandtl_from_diffusivities — Pr = nu/alpha
 *
 * @param nu    kinematic viscosity [m^2/s]
 * @param alpha thermal diffusivity [m^2/s]
 * @return      Prandtl number [-]
 */
double prandtl_from_diffusivities(double nu, double alpha);

/**
 * prandtl_thermal_diffusivity — alpha = nu/Pr
 */
double prandtl_thermal_diffusivity(double nu, double Pr);

/**
 * prandtl_kinematic_viscosity — nu = Pr.alpha
 */
double prandtl_kinematic_viscosity(double Pr, double alpha);

/**
 * prandtl_temperature_dependence_air — Pr(T) for air.
 *
 * Approximate fit for air at 1 atm, 200 K < T < 2000 K.
 * Pr(T) ≈ 0.69 + 0.0001.(T-300) [weak T-dependence for gases]
 *
 * @param T  temperature [K]
 * @return   Prandtl number [-], -1 if T out of range
 *
 * Reference: Incropera & DeWitt (2007), Table A.4.
 */
double prandtl_air(double T);

/**
 * prandtl_temperature_dependence_water — Pr(T) for saturated water.
 *
 * Strong T-dependence: Pr ~ 13 at 0 degC, ~2 at 100 degC, ~1 at 200 degC.
 * Fit: Pr(T) ≈ a + b.T + c.T^2 (piecewise polynomial)
 *
 * @param T  temperature [ degC], 0 < T < 200
 * @return   Prandtl number [-], -1 if T out of range
 *
 * Reference: NIST / Incropera Appendix A.
 */
double prandtl_water(double T);

/**
 * prandtl_engine_oil — Pr(T) for unused engine oil.
 *
 * Pr varies from ~10⁴ at 0 degC to ~100 at 150 degC.
 *
 * @param T  temperature [ degC], 0 < T < 200
 * @return   Prandtl number [-], -1 if out of range
 */
double prandtl_engine_oil(double T);

/* ==========================================================================
 * L1-L2: Schmidt Number — Momentum diffusivity / Mass diffusivity
 * ========================================================================== */

/**
 * schmidt_basic — Sc = nu/D_AB
 *
 * @param nu   kinematic viscosity [m^2/s]
 * @param D_AB binary mass diffusivity [m^2/s]
 * @return     Schmidt number [-]
 */
double schmidt_basic(double nu, double D_AB);

/**
 * schmidt_from_material — Sc = mu/(rho.D_AB)
 *
 * @param mu   dynamic viscosity [Pa.s]
 * @param rho  density [kg/m^3]
 * @param D_AB binary mass diffusivity [m^2/s]
 * @return     Schmidt number [-]
 */
double schmidt_from_material(double mu, double rho, double D_AB);

/**
 * schmidt_mass_diffusivity — recover D_AB from Sc: D = nu/Sc
 *
 * @param nu  kinematic viscosity [m^2/s]
 * @param Sc  Schmidt number [-]
 * @return    D_AB [m^2/s]
 */
double schmidt_mass_diffusivity(double nu, double Sc);

/**
 * schmidt_gas_estimate — typical Sc values for gases.
 *
 * For most gas pairs at 1 atm, 300 K: Sc ≈ 0.5-2.0.
 * Typical: CO2 in air -> Sc ≈ 1.0; H2O in air -> Sc ≈ 0.6.
 *
 * @return typical Sc = 1.0 for gases
 */
double schmidt_gas_typical(void);

/**
 * schmidt_liquid_estimate — typical Sc for solutes in liquids.
 *
 * For dilute solutions in water: Sc ≈ 10^2-10^3.
 * For polymers in solvents: Sc ≈ 10⁴-10⁶.
 *
 * @param solute_type  0 = small molecule in water, 1 = protein in water,
 *                     2 = polymer in organic solvent
 * @return             typical Sc estimate [-]
 */
double schmidt_liquid_typical(int solute_type);

/* ==========================================================================
 * L1-L2: Lewis Number — Thermal diffusivity / Mass diffusivity
 * ========================================================================== */

/**
 * lewis_basic — Le = alpha/D_AB = Sc/Pr
 *
 * @param alpha thermal diffusivity [m^2/s]
 * @param D_AB  mass diffusivity [m^2/s]
 * @return      Lewis number [-]
 */
double lewis_basic(double alpha, double D_AB);

/**
 * lewis_from_sc_pr — Le = Sc/Pr
 *
 * @param Sc  Schmidt number [-]
 * @param Pr  Prandtl number [-]
 * @return    Lewis number [-]
 */
double lewis_from_sc_pr(double Sc, double Pr);

/**
 * lewis_typical — return typical Le for engineering categories.
 *
 * @param category  0 = gas mixture (Le ~ 1), 1 = aqueous solution (Le ~ 100),
 *                  2 = liquid metal (Le ~ 10⁴)
 * @return          typical Lewis number [-]
 *
 *   Gases:      Le ≈ Sc/Pr = 1.0/0.7 ≈ 1.4 (heat and mass at similar rates)
 *   Water:      Le ≈ 500/7 ≈ 70 (heat diffuses much faster)
 *   Liquid metals: Le ≈ 1/0.01 ≈ 100 (mass diffuses faster, but both fast)
 */
double lewis_typical(int category);

/* ==========================================================================
 * L3: Diffusivity Models for Property Estimation
 * ========================================================================== */

/**
 * thermal_diffusivity_definition — alpha = k/(rho.cp)
 *
 * @param k    thermal conductivity [W/(m.K)]
 * @param rho  density [kg/m^3]
 * @param cp   specific heat [J/(kg.K)]
 * @return     alpha [m^2/s]
 */
double thermal_diffusivity_definition(double k, double rho, double cp);

/**
 * mass_diffusivity_gas_chapman_enskog — binary gas diffusivity.
 *
 * D_AB = (1.86×10⁻⁷.T^(3/2))/(p.sigma^2_AB.Omega_D) . (1/M_A + 1/M_B)^(1/2)
 *
 * Simplified engineering form at 1 atm, 300 K:
 *   D_AB ≈ (1.0×10⁻⁷).T^(1.75)/p
 *
 * @param T    temperature [K]
 * @param p    pressure [Pa]
 * @param M_A  molecular weight of A [kg/kmol]
 * @param M_B  molecular weight of B [kg/kmol]
 * @return     D_AB [m^2/s], -1 on invalid inputs
 *
 * Reference: Bird, Stewart, Lightfoot (2007), Eq. 17.3-12.
 * Accuracy: ±10% for nonpolar gas pairs.
 */
double mass_diffusivity_gas_estimate(double T, double p, double M_A, double M_B);

/**
 * mass_diffusivity_liquid_wilke_chang — liquid diffusivity.
 *
 * D_AB = (7.4×10⁻¹^2).(φ_B.M_B)^(1/2).T / (mu_B.V_A^0.6)
 *
 * @param T    temperature [K]
 * @param mu_B solvent viscosity [cP = mPa.s]
 * @param V_A  solute molar volume at boiling point [cm^3/mol]
 * @param M_B  solvent molecular weight [g/mol]
 * @param phi_B solvent association parameter (1.0 for nonpolar, 1.9 for ethanol,
 *              2.6 for water, 1.5 for methanol)
 * @return     D_AB [m^2/s], -1 on invalid inputs
 *
 * Reference: Wilke-Chang (1955).
 * Accuracy: ±20% for dilute solutions.
 */
double mass_diffusivity_liquid_wilke_chang(double T, double mu_B, double V_A,
                                           double M_B, double phi_B);

/* ==========================================================================
 * L2: Boundary Layer Ratios
 * ========================================================================== */

/**
 * thermal_to_velocity_bl_ratio_laminar — delta_T/delta for laminar BL.
 *
 * delta_T/delta ≈ Pr^(-1/3) for Pr > 0.5 (Pohlhausen solution)
 *
 * For liquid metals (Pr << 1): delta_T/delta ≈ 1.325.Pr^(-1/2)
 *
 * @param Pr  Prandtl number [-]
 * @return    delta_T/delta [-]
 */
double thermal_to_velocity_bl_ratio(double Pr);

/**
 * concentration_to_velocity_bl_ratio_laminar — delta_C/delta for laminar BL.
 *
 * delta_C/delta ≈ Sc^(-1/3) for Sc > 0.5
 *
 * @param Sc  Schmidt number [-]
 * @return    delta_C/delta [-]
 */
double concentration_to_velocity_bl_ratio(double Sc);

/**
 * triple_analogy_bl_ratios — compute all three BL thickness ratios.
 *
 * Simultaneously compute delta_T/delta and delta_C/delta from Pr and Sc.
 *
 * @param Pr       Prandtl number [-]
 * @param Sc       Schmidt number [-]
 * @param[out] ratio_T thermal/velocity ratio
 * @param[out] ratio_C concentration/velocity ratio
 * @return         0 on success, -1 on invalid inputs
 */
int triple_analogy_bl_ratios(double Pr, double Sc,
                             double *ratio_T, double *ratio_C);

/* ==========================================================================
 * L5: Combined Transport Analogies (Colburn, Chilton-Colburn)
 * ========================================================================== */

/**
 * colburn_analogy_verify — verify the Chilton-Colburn analogy.
 *
 * j_H = j_D = f/8
 * where j_H = St_H.Pr^(2/3), j_D = St_M.Sc^(2/3)
 *
 * @param St_H  Stanton number for heat [-]
 * @param Pr    Prandtl number [-]
 * @param St_M  Stanton number for mass [-]
 * @param Sc    Schmidt number [-]
 * @param f     Darcy friction factor [-]
 * @param[out] j_H  computed j-factor for heat
 * @param[out] j_D  computed j-factor for mass
 * @return      0 on success, -1 on invalid inputs
 *
 * Returns the absolute difference |j_H - j_D| via output params.
 * The Chilton-Colburn analogy holds when j_H ≈ j_D ≈ f/8.
 */
int colburn_analogy_verify(double St_H, double Pr, double St_M, double Sc,
                           double f, double *j_H, double *j_D);

/**
 * reynolds_analogy_heat — simplest analogy: St = f/8.
 *
 * @param f  Darcy friction factor [-]
 * @return   St_H = f/8
 *
 * Valid for: Pr ≈ 1 (gases), turbulent flow.
 * This is the original Reynolds analogy (1874).
 */
double reynolds_analogy_stanton(double f);

/**
 * reynolds_analogy_nusselt — Nu from Reynolds analogy.
 *
 * Nu = (f/8).Re.Pr
 *
 * @param f   friction factor [-]
 * @param Re  Reynolds number [-]
 * @param Pr  Prandtl number [-]
 * @return    Nusselt number [-]
 *
 * Good estimate for gases (Pr ~ 0.7) in turbulent flow.
 */
double reynolds_analogy_nusselt(double f, double Re, double Pr);

/**
 * chilton_colburn_mass_transfer — Sh from j_D = f/8.
 *
 * Sh = (f/8).Re.Sc^(1/3)
 *
 * @param f   friction factor [-]
 * @param Re  Reynolds number [-]
 * @param Sc  Schmidt number [-]
 * @return    Sherwood number [-]
 *
 * Valid for: 0.6 < Sc < 3000, turbulent flow.
 */
double chilton_colburn_sherwood(double f, double Re, double Sc);

/**
 * heat_mass_analogy_nusselt — predict Sh from Nu using Pr/Sc analogy.
 *
 * Sh = Nu.(Sc/Pr)^(1/3)
 *
 * @param Nu  Nusselt number [-]
 * @param Pr  Prandtl number [-]
 * @param Sc  Schmidt number [-]
 * @return    predicted Sherwood number [-]
 *
 * This is the heat-mass transfer analogy: same geometry,
 * same flow -> same functional form, different property exponent.
 */
double heat_mass_analogy_nusselt_to_sherwood(double Nu, double Pr, double Sc);

#ifdef __cplusplus
}
#endif

#endif /* PRANDTL_SCHMIDT_H */

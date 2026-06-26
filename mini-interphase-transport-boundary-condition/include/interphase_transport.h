/**
 * interphase_transport.h
 * =======================
 * Transport coefficient calculation and flux models across phase interfaces.
 *
 * Reference: Bird, Stewart, Lightfoot (2007) "Transport Phenomena"
 *   Chapter 14 - Interphase Transport in Nonisothermal Systems
 *   Chapter 21 - Macroscopic Balances for Isothermal Systems
 *   Chapter 22 - Interphase Transport in Isothermal Systems
 *
 * Treybal, R.E. (1980) "Mass-Transfer Operations" (3rd ed.)
 * Cussler, E.L. (2009) "Diffusion: Mass Transfer in Fluid Systems" (3rd ed.)
 *
 * Knowledge coverage: L2 Core Concepts, L3 Engineering Quantities, L5 Engineering Methods
 */

#ifndef INTERPHASE_TRANSPORT_H
#define INTERPHASE_TRANSPORT_H

#include "interphase_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * L5: Mass Transfer Coefficient Models
 * ============================================================================ */

/**
 * Compute mass transfer coefficient from two-film theory (Whitman, 1923).
 * 1/K_overall = 1/k_A + H/(k_B) for gas-side basis.
 * 1/K_overall = 1/(H*k_A) + 1/k_B for liquid-side basis.
 * Assumes linear equilibrium y = m*x (Henry's law form).
 *
 * @param k_a       Individual mass transfer coefficient phase A (m/s), > 0
 * @param k_b       Individual mass transfer coefficient phase B (m/s), > 0
 * @param m_partition Equilibrium slope m = dy/dx (-), > 0
 * @param film      Output two-film model with overall coefficients
 * Complexity: O(1)
 *
 * Theorem: 1/K_og = 1/k_g + m/k_l  (gas-film basis).
 * Theorem: 1/K_ol = 1/(m*k_g) + 1/k_l  (liquid-film basis).
 * Theorem: For m << 1 (highly soluble gas), K_og ~ k_g (gas-film controlled).
 * Theorem: For m >> 1 (low solubility), K_ol ~ k_l (liquid-film controlled).
 */
void two_film_mass_transfer(double k_a, double k_b, double m_partition,
                             TwoFilmModel *film);

/**
 * Compute mass transfer coefficient from penetration theory (Higbie, 1935).
 * k_L(t) = sqrt(D / (pi * t)), instantaneous coefficient.
 * k_L_avg = (1/t_e) * integral_0^{t_e} k_L(t) dt = 2*sqrt(D/(pi*t_e)).
 *
 * @param diffusivity    Molecular diffusivity D (m^2/s), > 0
 * @param exposure_time  Contact time t_e (s), > 0
 * @param model          Output penetration model with computed coefficients
 * Complexity: O(1)
 *
 * Theorem: k_avg = 2 * sqrt(D / (pi * t_e)).
 * Theorem: Compared to film theory, penetration gives k ~ D^{0.5}
 *          while film theory gives k ~ D^{1.0} (different scaling laws).
 */
void penetration_mass_transfer(double diffusivity, double exposure_time,
                                PenetrationModel *model);

/**
 * Compute mass transfer coefficient from surface renewal theory
 * (Danckwerts, 1951). Elements at the surface are randomly replaced
 * with age distribution s * exp(-s * t), giving k = sqrt(D * s).
 *
 * @param diffusivity   Molecular diffusivity D (m^2/s), > 0
 * @param renewal_rate  Surface renewal rate s (1/s), > 0
 * @param model         Output surface renewal model
 * Complexity: O(1)
 *
 * Theorem: k = sqrt(D * s).
 * Theorem: Mean element age = 1/s.
 * Theorem: For mass transfer with 1st order reaction, k = sqrt(D*(s + k1)).
 */
void surface_renewal_mass_transfer(double diffusivity, double renewal_rate,
                                    SurfaceRenewalModel *model);

/**
 * Compare mass transfer theories for a given system.
 * Returns a dimensionless ratio of predictions.
 *
 * @param film_k        Film theory prediction k_f (m/s)
 * @param penetration_k Penetration theory prediction k_p (m/s)
 * @param renewal_k     Surface renewal theory prediction k_s (m/s)
 * @param ratios        Output array: [k_p/k_f, k_s/k_f, k_s/k_p]
 * Complexity: O(1)
 *
 * Theorem: For long contact times, film and penetration converge.
 * Theorem: Surface renewal interpolates between film (s->0, randomness low)
 *          and penetration (s*t_e = 1, single-age ensemble).
 */
void compare_mass_transfer_theories(double film_k, double penetration_k,
                                     double renewal_k, double ratios[3]);

/* ============================================================================
 * L5: Heat Transfer Coefficient Correlations at Interfaces
 * ============================================================================ */

/**
 * Compute Nusselt number for laminar forced convection over a flat plate.
 * Nu_x = 0.332 * Re_x^(1/2) * Pr^(1/3)  for Pr >= 0.6 (Pohlhausen, 1921).
 * Average: Nu_L = 0.664 * Re_L^(1/2) * Pr^(1/3).
 *
 * @param Re_x    Local Reynolds number rho*U*x/mu, > 0
 * @param Pr      Prandtl number, > 0
 * @param Nu_x    Output local Nusselt number
 * @param Nu_avg  Output average Nusselt number over length L
 * Complexity: O(1)
 *
 * Theorem: Nu ~ Re^(1/2) Pr^(1/3) for laminar flow.
 * Theorem: Transition at Re_crit ~ 5e5 yields turbulent Nu ~ Re^(4/5) Pr^(1/3).
 */
void nusselt_flat_plate_laminar(double Re_x, double Pr,
                                 double *Nu_x, double *Nu_avg);

/**
 * Compute Nusselt number for turbulent forced convection over a flat plate.
 * Nu_x = 0.0296 * Re_x^(4/5) * Pr^(1/3) (local, Colburn analogy).
 * Nu_L = 0.037 * Re_L^(4/5) * Pr^(1/3) (average, for Re up to ~10^7).
 *
 * @param Re_x    Local Reynolds number, > 5e5 (turbulent regime)
 * @param Pr      Prandtl number, > 0
 * @param Nu_x    Output local Nusselt number
 * @param Nu_avg  Output average Nusselt number
 * Complexity: O(1)
 */
void nusselt_flat_plate_turbulent(double Re_x, double Pr,
                                   double *Nu_x, double *Nu_avg);

/**
 * Compute Nusselt number for laminar flow inside a circular tube.
 * Constant wall temperature: Nu = 3.66 (fully developed).
 * Constant wall heat flux: Nu = 4.36 (fully developed).
 * Laminar entry region: Nu = 1.86 * (Re*Pr*D/L)^(1/3) * (mu/mu_w)^0.14.
 *
 * @param Re_D    Reynolds number based on diameter, > 0
 * @param Pr      Prandtl number, > 0
 * @param L_over_D Length-to-diameter ratio, > 1
 * @param mu_ratio Bulk-to-wall viscosity ratio, ~1 for moderate heating
 * @param is_const_T 1 for constant Tw BC, 0 for constant qw BC
 * @param Nu      Output Nusselt number
 * Complexity: O(1)
 */
void nusselt_pipe_laminar(double Re_D, double Pr, double L_over_D,
                           double mu_ratio, int is_const_T, double *Nu);

/**
 * Compute Nusselt number for turbulent flow in a smooth circular tube.
 * Dittus-Boelter (1930): Nu = 0.023 * Re^(4/5) * Pr^n
 *   n = 0.4 for heating, n = 0.3 for cooling.
 * Gnielinski (1976): Nu = (f/8)*(Re-1000)*Pr/(1 + 12.7*sqrt(f/8)*(Pr^(2/3)-1)).
 *
 * @param Re_D      Reynolds number, > 2300
 * @param Pr        Prandtl number, > 0.5
 * @param heating   1 for heating (n=0.4), 0 for cooling (n=0.3)
 * @param Nu_dittus Output Dittus-Boelter Nusselt number
 * @param Nu_gnielinski Output Gnielinski Nusselt number
 * Complexity: O(1)
 */
void nusselt_pipe_turbulent(double Re_D, double Pr, int heating,
                             double *Nu_dittus, double *Nu_gnielinski);

/**
 * Compute Nusselt number for natural convection from a vertical plate.
 * Churchill-Chu correlation: Nu = {0.825 + 0.387*Ra^(1/6)/[1+(0.492/Pr)^(9/16)]^(8/27)}^2.
 * Valid for all Ra (laminar, transition, turbulent).
 *
 * @param Ra      Rayleigh number Gr*Pr, > 0
 * @param Pr      Prandtl number, > 0
 * @param Nu      Output Nusselt number
 * Complexity: O(1)
 */
void nusselt_natural_convection_vertical(double Ra, double Pr, double *Nu);

/**
 * Compute Nusselt number for condensation on a vertical surface.
 * Nusselt (1916) film condensation theory:
 * Nu_L = 0.943 * [g*rho_L*(rho_L-rho_v)*h_fg*L^3 / (mu_L*k_L*(T_sat-T_w))]^(1/4).
 *
 * @param g         Gravitational acceleration (m/s^2), typically 9.81
 * @param rho_L     Liquid density (kg/m^3)
 * @param rho_v     Vapor density (kg/m^3)
 * @param h_fg      Latent heat (J/kg)
 * @param L         Plate length (m)
 * @param mu_L      Liquid viscosity (Pa*s)
 * @param k_L       Liquid thermal conductivity (W/(m*K))
 * @param delta_T   T_sat - T_w (K), > 0
 * @param Nu        Output average Nusselt number
 * Complexity: O(1)
 */
void nusselt_film_condensation(double g, double rho_L, double rho_v,
                                double h_fg, double L, double mu_L,
                                double k_L, double delta_T, double *Nu);

/**
 * Compute Nusselt number for nucleate pool boiling.
 * Rohsenow (1952) correlation:
 * q = mu_L * h_fg * sqrt(g*(rho_L-rho_v)/sigma) * (cp_L*delta_T/(C_sf*h_fg*Pr_L^n))^3.
 * where C_sf = surface-fluid combination constant.
 *
 * @param q_flux    Heat flux (W/m^2), > 0
 * @param rho_L     Liquid density (kg/m^3)
 * @param rho_v     Vapor density (kg/m^3)
 * @param h_fg      Latent heat (J/kg)
 * @param cp_L      Liquid specific heat (J/(kg*K))
 * @param mu_L      Liquid viscosity (Pa*s)
 * @param Pr_L      Liquid Prandtl number
 * @param sigma     Surface tension (N/m)
 * @param C_sf      Rohsenow surface-fluid constant (~0.013 for water-Pt)
 * @param n_exp     Pr exponent (1.0 for water, 1.7 for others)
 * @param delta_T   Output wall superheat T_w - T_sat (K)
 * Complexity: O(1)
 */
void rohsenow_pool_boiling(double q_flux, double rho_L, double rho_v,
                            double h_fg, double cp_L, double mu_L,
                            double Pr_L, double sigma, double C_sf,
                            double n_exp, double *delta_T);

/* ============================================================================
 * L5: Sherwood Number Correlations for Mass Transfer at Interfaces
 * ============================================================================ */

/**
 * Compute Sherwood number for gas absorption in a falling liquid film.
 * Laminar film, short contact time: Sh = 2/sqrt(pi) * sqrt(Re * Sc * d/L).
 * Long contact time (penetration): Sh = 3.41 (fully developed concentration).
 *
 * @param Re_film Film Reynolds number 4*Gamma/mu, > 0
 * @param Sc      Schmidt number, > 0
 * @param L_over_d Length-to-film-thickness ratio, > 0
 * @param Sh      Output Sherwood number
 * Complexity: O(1)
 */
void sherwood_falling_film(double Re_film, double Sc, double L_over_d,
                            double *Sh);

/**
 * Compute Sherwood number for mass transfer from a single sphere.
 * Froessling (1938): Sh = 2.0 + 0.552 * Re^(1/2) * Sc^(1/3).
 * The constant 2.0 is the pure diffusion limit (creeping flow).
 * Valid for Re < 1000, 0.6 < Sc < 2.7.
 *
 * @param Re_sphere Particle Reynolds number, > 0
 * @param Sc        Schmidt number, > 0
 * @param Sh        Output Sherwood number
 * Complexity: O(1)
 *
 * Theorem: As Re -> 0, Sh -> 2 (pure molecular diffusion around sphere).
 */
void sherwood_single_sphere(double Re_sphere, double Sc, double *Sh);

/**
 * Compute Sherwood number for mass transfer in packed beds.
 * Wakao & Funazkri (1978): Sh = 2.0 + 1.1 * Re^(0.6) * Sc^(1/3).
 * Valid for 3 < Re < 10000, gas-solid systems.
 *
 * @param Re_p    Particle Reynolds number rho*u*d_p/mu, > 0
 * @param Sc      Schmidt number, > 0
 * @param epsilon Bed void fraction, (0, 1)
 * @param Sh      Output Sherwood number
 * Complexity: O(1)
 */
void sherwood_packed_bed(double Re_p, double Sc, double epsilon, double *Sh);

/**
 * Compute Sherwood number for mass transfer in bubble columns.
 * Akita & Yoshida (1974): Sh = 0.5 * Sc^(1/2) * Bo^(3/8) * Ga^(1/4).
 * where Bo = g*d_b^2*rho_L/sigma, Ga = g*d_b^3/nu_L^2.
 *
 * @param d_b     Bubble diameter (m), > 0
 * @param g       Gravitational acceleration (m/s^2)
 * @param rho_L   Liquid density (kg/m^3)
 * @param sigma   Surface tension (N/m)
 * @param nu_L    Liquid kinematic viscosity (m^2/s)
 * @param D_L     Liquid phase diffusivity (m^2/s)
 * @param Sh      Output Sherwood number
 * Complexity: O(1)
 */
void sherwood_bubble_column(double d_b, double g, double rho_L,
                             double sigma, double nu_L, double D_L,
                             double *Sh);

/* ============================================================================
 * L2: Interfacial Resistance and Enhancement
 * ============================================================================ */

/**
 * Compute the overall mass transfer resistance with chemical reaction
 * enhancement. E = k_with_rxn / k_without_rxn.
 * For instantaneous reaction: E_inf = 1 + D_B*C_B / (D_A*C_Ai) (film theory).
 *
 * @param D_A       Diffusivity of transferring species A (m^2/s)
 * @param D_B       Diffusivity of reactant B (m^2/s)
 * @param C_B_bulk  Bulk concentration of reactant B (mol/m^3)
 * @param C_Ai      Interfacial concentration of A (mol/m^3)
 * @param nu        Stoichiometric coefficient of B per A
 * @param E_inf     Output enhancement factor for instantaneous reaction
 * Complexity: O(1)
 *
 * Theorem: E_inf >= 1 always (reaction enhances mass transfer).
 * Theorem: For slow reaction (Ha < 0.3), E ~ 1.
 * Theorem: For fast pseudo-1st-order (Ha > 3), E = Ha.
 */
void enhancement_factor_instantaneous(double D_A, double D_B,
                                       double C_B_bulk, double C_Ai,
                                       double nu, double *E_inf);

/**
 * Compute Hatta number (Ha) for mass transfer with chemical reaction.
 * Ha = sqrt(k_rxn * D_A) / k_L  for first-order reaction.
 * Ha >> 1: reaction in film (fast); Ha << 1: reaction in bulk (slow).
 *
 * @param k_rxn    First-order reaction rate constant (1/s), >= 0
 * @param D_A      Diffusivity of reactant A (m^2/s), > 0
 * @param k_L      Mass transfer coefficient without reaction (m/s), > 0
 * @param Ha       Output Hatta number, dimensionless
 * Complexity: O(1)
 *
 * Theorem: Ha < 0.3 -> slow reaction, no enhancement.
 * Theorem: 0.3 < Ha < 3 -> intermediate regime.
 * Theorem: Ha > 3 -> fast reaction, E = Ha (significant enhancement).
 */
void hatta_number(double k_rxn, double D_A, double k_L, double *Ha);

/**
 * Compute overall mass transfer coefficient with chemical reaction.
 * K_overall_with_rxn = K_overall * E.
 *
 * @param K_overall  Overall mass transfer coefficient without reaction (m/s)
 * @param E          Enhancement factor (>= 1)
 * @param K_with_rxn Output enhanced overall coefficient (m/s)
 * Complexity: O(1)
 */
void overall_coefficient_with_reaction(double K_overall, double E,
                                        double *K_with_rxn);

/* ============================================================================
 * L3: Engineering Quantity - Interfacial Area Estimation
 * ============================================================================ */

/**
 * Compute specific interfacial area for a packed column.
 * Onda et al. (1968) correlation:
 * a/a_t = 1 - exp[-1.45*(sigma_c/sigma)^0.75 * Re_L^0.1 * Fr_L^{-0.05} * We_L^0.2].
 *
 * @param a_t       Total packing surface area per unit volume (m^2/m^3)
 * @param sigma_c   Critical surface tension of packing material (N/m)
 * @param sigma_L   Liquid surface tension (N/m)
 * @param Re_L      Liquid Reynolds number rho_L*u_L/(a_t*mu_L)
 * @param Fr_L      Liquid Froude number a_t*u_L^2/g
 * @param We_L      Liquid Weber number rho_L*u_L^2/(sigma_L*a_t)
 * @param a_eff     Output effective interfacial area (m^2/m^3)
 * Complexity: O(1)
 */
void interfacial_area_packed_column(double a_t, double sigma_c,
                                     double sigma_L, double Re_L,
                                     double Fr_L, double We_L, double *a_eff);

/* ============================================================================
 * L6/L7: Application Functions (defined in interphase_applications.c)
 * ============================================================================ */

/* Falling film gas absorption */
void falling_film_absorption(double Q_L, double L, double W,
                              double p_CO2, double H_CO2,
                              double C_in, double k_L, double D_CO2,
                              double rho_L, double mu_L,
                              double *N_abs, double *C_out);

/* Evaporative cooling / wet-bulb temperature */
void evaporative_cooling(double T_air, double RH, double p_atm,
                          double h, double k_m, double h_fg,
                          double cp_air, double Le,
                          double *T_wet_bulb, double *evap_flux);

/* Distillation tray efficiency (Murphree) */
void murphree_tray_efficiency(double k_g, double k_l, double m,
                               double a, double h_f,
                               double G, double L,
                               double *E_MV);

/* Membrane contactor sizing */
void membrane_contactor_sizing(double Q_L, double C_in, double C_out_target,
                                double C_star, double K_overall,
                                double d_fiber, double packing_frac,
                                double *L_module, double *A_total);

/* Condensation on horizontal tube (Nusselt) */
void condensation_horizontal_tube(double D_tube, double T_sat, double T_wall,
                                   double rho_L, double rho_v,
                                   double k_L, double h_fg, double mu_L,
                                   int N_tubes,
                                   double *h_avg, double *Q_cond);

/* Wastewater aeration design */
void wastewater_aeration_design(double Q_wastewater, double BOD_in,
                                 double BOD_out, double T,
                                 double alpha_factor, double beta_factor,
                                 double *V_tank, double *Q_air,
                                 double *P_blower);

/* PWR steam generator sizing */
void pwr_steam_generator_sizing(double P_thermal,
                                 double T_primary_in, double T_primary_out,
                                 double T_sat_secondary,
                                 double cp_primary, double U_overall,
                                 double D_tube, double L_tube,
                                 double *N_tubes, double *A_ht,
                                 double *m_dot_primary);

/* PEM fuel cell water flux */
void pem_fuel_cell_water_flux(double i_current, double T_cell,
                               double lambda_anode, double lambda_cathode,
                               double t_membrane, double rho_dry,
                               double EW,
                               double *J_net, double *drag_flux,
                               double *diff_flux);

#ifdef __cplusplus
}
#endif

#endif /* INTERPHASE_TRANSPORT_H */

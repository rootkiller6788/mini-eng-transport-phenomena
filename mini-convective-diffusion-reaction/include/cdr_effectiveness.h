/**
 * @file cdr_effectiveness.h
 * @brief Catalyst Effectiveness Factor and Intraparticle Diffusion-Reaction
 *
 * Effectiveness factor (eta) for various catalyst geometries,
 * Weisz-Prater criterion for diffusion limitations,
 * and Thiele modulus approach for porous catalyst analysis.
 *
 * Reference: Bird, Stewart, Lightfoot "Transport Phenomena" (2007), Ch. 18, 21
 * Reference: Fogler "Elements of Chemical Reaction Engineering" (2016), Ch. 12, 15
 * Reference: Froment, Bischoff, De Wilde "Chemical Reactor Analysis & Design" (2011)
 *
 * L1: Effectiveness factor, Thiele modulus, internal diffusion limitation
 * L2: Diffusion-reaction coupling in porous catalysts
 * L4: Mass balance with reaction inside catalyst pellet
 * L5: Weisz-Prater criterion, Mears criterion for external transport
 * L6: Catalyst pellet design, washcoat optimization
 */

#ifndef CDR_EFFECTIVENESS_H
#define CDR_EFFECTIVENESS_H

#include "cdr_core.h"

/* ---------------------------------------------------------------------------
 * L1: Effectiveness Factor Definition
 * ------------------------------------------------------------------------- */

/**
 * @brief Effectiveness factor = (actual rate) / (rate at surface conditions).
 *
 * eta = r_actual / r(C_s, T_s)
 *
 * For eta ~ 1: No diffusion limitation (kinetic control)
 * For eta << 1: Strong diffusion limitation (diffusion control)
 */
typedef struct {
    double eta;                 /**< Effectiveness factor [-] */
    double eta_internal;        /**< Internal (pore diffusion) factor */
    double eta_external;        /**< External (film) factor */
    double eta_overall;         /**< Overall (internal * external) */
} EffectivenessFactors;

/**
 * @brief Diffusion limitation diagnostics.
 */
typedef struct {
    double Weisz_Prater;   /**< Weisz-Prater parameter: C_WP = eta * phi^2 */
    double Mears_external; /**< Mears criterion for external transport */
    double Carberry;       /**< Carberry number: Ca = r_obs / (k_c * a * C_bulk) */
    int    is_internal_limited;  /**< 1 if internal diffusion is limiting */
    int    is_external_limited;  /**< 1 if external mass transfer is limiting */
} DiffusionLimitationDiagnosis;

/* ---------------------------------------------------------------------------
 * L4: Effectiveness Factor for Different Geometries
 * ------------------------------------------------------------------------- */

/**
 * @brief Effectiveness factor for 1st-order reaction in a slab.
 *
 * eta = tanh(phi) / phi
 *
 * where phi = L * sqrt(k/D_eff)
 *
 * @param phi  Thiele modulus [-]
 * @return Effectiveness factor eta [0-1]
 *
 * Complexity: O(1)
 */
double cdr_effectiveness_slab(double phi);

/**
 * @brief Effectiveness factor for 1st-order reaction in an infinite cylinder.
 *
 * eta = 2 * I_1(phi) / (phi * I_0(phi))
 *
 * where I_0, I_1 are modified Bessel functions of the first kind.
 *
 * Uses a rational approximation for computational speed.
 *
 * @param phi  Thiele modulus (based on R/2) [-]
 * @return Effectiveness factor eta [0-1]
 */
double cdr_effectiveness_cylinder(double phi);

/**
 * @brief Effectiveness factor for 1st-order reaction in a sphere.
 *
 * eta = (3/phi) * (1/tanh(phi) - 1/phi)
 *
 * @param phi  Thiele modulus (based on R/3) [-]
 * @return Effectiveness factor eta [0-1]
 */
double cdr_effectiveness_sphere(double phi);

/**
 * @brief Generalized effectiveness factor (any geometry, 1st-order).
 *
 * Uses the generalized Thiele modulus:
 *   Phi_gen = L_c * sqrt(k/D_eff) where L_c = V_p/S_p
 *
 * For Phi_gen > 5: eta ~ 1/Phi_gen (all geometries)
 * For Phi_gen < 0.4: eta ~ 1 (all geometries)
 *
 * @param phi      Generalized Thiele modulus [-]
 * @param geometry Catalyst geometry
 * @return Effectiveness factor eta [0-1]
 */
double cdr_effectiveness_general(double phi, CatalystGeometry geometry);

/**
 * @brief Effectiveness factor for nth-order irreversible reaction.
 *
 * Uses generalized Thiele modulus approach:
 *   phi_n = L_c * sqrt( (n+1)*k_n*C_s^(n-1) / (2*D_eff) )
 *
 * Then eta ~ tanh(phi_n)/phi_n (approximate, works for n != 1 as well)
 *
 * @param k_n    nth-order rate constant [(m^3/mol)^(n-1)/s]
 * @param n      Reaction order
 * @param C_s    Surface concentration [mol/m^3]
 * @param D_eff  Effective diffusivity [m^2/s]
 * @param L_c    Characteristic length V_p/S_p [m]
 * @return Effectiveness factor eta [0-1]
 */
double cdr_effectiveness_nth_order(double k_n, double n, double C_s,
                                    double D_eff, double L_c);

/**
 * @brief Effectiveness factor for Michaelis-Menten kinetics.
 *
 * Uses the general modulus approach with linearized kinetics.
 *
 * phi_MM = L_c * sqrt( V_max / (D_eff * K_m) ) / sqrt(1 + S_s/K_m)  (approximate)
 *
 * @param V_max  Maximum rate [mol/(m^3.s)]
 * @param K_m    Michaelis constant [mol/m^3]
 * @param S_s    Substrate concentration at surface [mol/m^3]
 * @param D_eff  Effective diffusivity [m^2/s]
 * @param L_c    Characteristic length [m]
 * @return Effectiveness factor eta [0-1]
 */
double cdr_effectiveness_michaelis_menten(double V_max, double K_m,
                                           double S_s, double D_eff,
                                           double L_c);

/* ---------------------------------------------------------------------------
 * L5: Diffusion Limitation Diagnostics
 * ------------------------------------------------------------------------- */

/**
 * @brief Weisz-Prater criterion for internal diffusion limitation.
 *
 * C_WP = r_obs * L_c^2 / (D_eff * C_s)
 *
 * C_WP << 1 : No internal diffusion limitation (eta ~ 1)
 * C_WP >> 1 : Strong internal diffusion limitation
 *
 * Practical threshold: C_WP > 1 ? investigate internal limitation.
 *
 * @param r_obs   Observed (measured) reaction rate [mol/(m^3.s)]
 * @param L_c     Characteristic length V_p/S_p [m]
 * @param D_eff   Effective diffusivity [m^2/s]
 * @param C_s     Surface concentration [mol/m^3]
 * @return Weisz-Prater parameter [-]
 */
double cdr_weisz_prater_criterion(double r_obs, double L_c,
                                   double D_eff, double C_s);

/**
 * @brief Mears criterion for external mass transfer limitation.
 *
 * C_M = r_obs * L_c * n / (k_c * C_bulk)
 *
 * C_M < 0.15 : External mass transfer negligible (< 5% error)
 * C_M > 0.15 : External mass transfer significant
 *
 * @param r_obs    Observed rate [mol/(m^3.s)]
 * @param L_c      Characteristic length [m]
 * @param n        Reaction order
 * @param k_c      Mass transfer coefficient [m/s]
 * @param C_bulk   Bulk concentration [mol/m^3]
 * @return Mears parameter for external transport [-]
 */
double cdr_mears_external_criterion(double r_obs, double L_c, double n,
                                     double k_c, double C_bulk);

/**
 * @brief Full diffusion limitation diagnosis.
 *
 * Combines Weisz-Prater and Mears criteria into a comprehensive assessment.
 *
 * @param r_obs    Observed rate [mol/(m^3.s)]
 * @param L_c      Characteristic length [m]
 * @param D_eff    Effective diffusivity [m^2/s]
 * @param C_s      Surface concentration [mol/m^3]
 * @param n        Reaction order
 * @param k_c      Mass transfer coefficient [m/s]
 * @param C_bulk   Bulk concentration [mol/m^3]
 * @param diag     Output: populated diagnosis structure
 */
void cdr_diagnose_limitations(double r_obs, double L_c, double D_eff,
                               double C_s, double n, double k_c,
                               double C_bulk, DiffusionLimitationDiagnosis *diag);

/* ---------------------------------------------------------------------------
 * L5: Observed Rate and Intrinsic Kinetics Recovery
 * ------------------------------------------------------------------------- */

/**
 * @brief Observed (apparent) rate from intrinsic kinetics and effectiveness.
 *
 * r_obs = eta * r_intrinsic(C_s)
 *
 * @param r_intrinsic  Intrinsic rate at surface conditions [mol/(m^3.s)]
 * @param eta          Effectiveness factor [-]
 * @return Observed (diffusion-disguised) rate [mol/(m^3.s)]
 */
double cdr_observed_rate(double r_intrinsic, double eta);

/**
 * @brief Recover intrinsic activation energy from observed (diffusion-limited).
 *
 * In strong pore diffusion limitation (eta << 1):
 *   E_obs = (E_true + E_diff) / 2  [for 1st order in sphere]
 *
 * where E_diff is activation energy for diffusion (usually ~10-20 kJ/mol).
 *
 * @param E_obs      Observed (apparent) activation energy [J/mol]
 * @param E_diff     Activation energy for diffusion [J/mol]
 * @return Estimated true activation energy [J/mol]
 */
double cdr_intrinsic_activation_energy(double E_obs, double E_diff);

/**
 * @brief Apparent reaction order under strong diffusion limitation.
 *
 * For true nth-order and strong diffusion limitation:
 *   n_app = (n + 1) / 2
 *
 * @param n_true  True reaction order
 * @return Apparent reaction order under diffusion limitation
 */
double cdr_apparent_reaction_order(double n_true);

/**
 * @brief Compute the concentration profile inside a spherical catalyst pellet.
 *
 * d/dr (r^2 * D_eff * dC/dr) = r^2 * k * C   (1st order)
 *
 * Analytical solution: C(r) = C_s * (R/r) * sinh(phi*r/R) / sinh(phi)
 *
 * @param r       Radial position [m] (0 <= r <= R)
 * @param R       Pellet radius [m]
 * @param C_s     Surface concentration [mol/m^3]
 * @param phi     Thiele modulus for sphere [-]
 * @return Concentration at r [mol/m^3]
 */
double cdr_pellet_concentration_profile(double r, double R, double C_s,
                                         double phi);

/**
 * @brief Compute optimal catalyst pellet size for given constraints.
 *
 * Balances pressure drop vs. effectiveness. For 1st-order sphere:
 * Minimizes cost per unit conversion.
 *
 * The optimal phi is approximately:
 *   phi_opt = 2.0 to 3.0 (for first-order, balancing kinetics and diffusion)
 *
 * @param k        Rate constant [1/s]
 * @param D_eff    Effective diffusivity [m^2/s]
 * @param target_eta Minimum acceptable effectiveness factor
 * @return Optimal pellet radius [m]
 */
double cdr_optimal_pellet_radius(double k, double D_eff, double target_eta);

/**
 * @brief Washcoat effectiveness factor for monolith catalysts.
 *
 * For a thin washcoat layer (slab geometry):
 *   eta = tanh(phi_w) / phi_w
 *   where phi_w = L_w * sqrt(k/D_eff,w)
 *
 * @param L_w     Washcoat thickness [m]
 * @param k       First-order rate constant [1/s] (per washcoat volume)
 * @param D_eff_w Effective diffusivity in washcoat [m^2/s]
 * @return Washcoat effectiveness factor [-]
 */
double cdr_washcoat_effectiveness(double L_w, double k, double D_eff_w);

#endif /* CDR_EFFECTIVENESS_H */

/**
 * @file transport_coefficients.c
 * @brief Implementation of transport coefficient computations
 *
 * Implements Sutherland's law, Chapman-Enskog kinetic theory,
 * Eucken formula, Wilke-Chang correlation, and mixture rules.
 *
 * Reference: Bird, Stewart, Lightfoot (2007), Transport Phenomena.
 *            Reid, Prausnitz, Poling (1987), Properties of Gases & Liquids.
 */

#include "transport_coefficients.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

/* ─── Physical Constants ────────────────────────────────────────── */

#define KB 1.380649e-23   /**< Boltzmann constant [J/K] */
#define R_GAS 8.314462618  /**< Universal gas constant [J/(mol·K)] */
#define PI 3.14159265358979323846
#define SQRT_PI 1.77245385090551602729

/* ─── L2: Sutherland's Law ──────────────────────────────────────── */

/**
 * Sutherland's law for gas viscosity.
 *
 * Derived from elementary kinetic theory with a temperature-dependent
 * collision cross-section. The Sutherland constant S accounts for the
 * attractive intermolecular potential.
 *
 * Typical values: Air: μ₀=1.716e-5 Pa·s, T₀=273.15 K, S=110.4 K
 *
 * Complexity: O(1)
 */
double sutherland_viscosity(double T, const SutherlandParams *params)
{
    if (!params || T <= 0.0) return -1.0;
    double ratio = T / params->T0;
    return params->mu0 * ratio * sqrt(ratio)
           * (params->T0 + params->S) / (T + params->S);
}

/* ─── L2: Chapman-Enskog Kinetic Theory ────────────────────────── */

/**
 * Chapman-Enskog viscosity for monatomic gases.
 *
 * μ = (5/16) · sqrt(π·M·kB·T) / (π·σ²·Ω_μ)
 *
 * This is the rigorous first-order solution of the Boltzmann equation
 * for dilute monatomic gases. The collision integral Ω_μ accounts for
 * the intermolecular potential (Lennard-Jones 12-6).
 *
 * For T* = kT/ε:
 *   Ω_μ ≈ 1.16145·(T*)^(-0.14874) + 0.52487·exp(-0.77320·T*)
 *          + 2.16178·exp(-2.43787·T*)  (Neufeld et al., 1972)
 *
 * Complexity: O(1)
 */
static double collision_integral_mu(double T_star)
{
    /* Neufeld, Janzen, Aziz (1972) empirical fit */
    if (T_star <= 0.0) return 2.0;  /* limit for very low T */
    return 1.16145 * pow(T_star, -0.14874)
         + 0.52487 * exp(-0.77320 * T_star)
         + 2.16178 * exp(-2.43787 * T_star);
}

double chapman_enskog_viscosity(double T, const KineticTheoryParams *params)
{
    if (!params || T <= 0.0 || params->sigma <= 0.0) return -1.0;
    double T_star = T / params->epsilon_kB;
    double omega = collision_integral_mu(T_star);
    double numerator = (5.0 / 16.0) * sqrt(PI * params->M * KB * T);
    double denominator = PI * params->sigma * params->sigma * omega;
    return numerator / denominator;
}

/* ─── L2: Power-Law Viscosity ──────────────────────────────────── */

double power_law_viscosity(double T, const PowerLawViscosityParams *params)
{
    if (!params || T <= 0.0 || params->T_ref <= 0.0) return -1.0;
    return params->mu_ref * pow(T / params->T_ref, params->n);
}

/* ─── L2: Eucken Thermal Conductivity ──────────────────────────── */

/**
 * Eucken formula for polyatomic gas thermal conductivity.
 *
 * k = μ · (Cp + 5/4 · R/M)     for monatomic: Cp = 5/2·R/M → k = (15/4)·(R/M)·μ
 * k = μ · (Cp + 9/4 · R/M)     modified Eucken for polyatomic gases
 *
 * The correction accounts for internal degrees of freedom
 * (rotation, vibration) transferring energy.
 *
 * Complexity: O(1)
 */
double eucken_thermal_conductivity(double mu, double cp, double M)
{
    if (M <= 0.0 || mu < 0.0) return -1.0;
    double R_over_M = R_GAS / M;
    /* Modified Eucken for polyatomic gases */
    return mu * (cp + 1.25 * R_over_M);
}

/* ─── L2: Chapman-Enskog Diffusivity ───────────────────────────── */

/**
 * Collision integral for diffusion (differs from viscosity).
 *
 * Ω_D ≈ 1.06036·(T*)^(-0.15610) + 0.19300·exp(-0.47635·T*)
 *        + 1.03587·exp(-1.52996·T*) + 1.76474·exp(-3.89411·T*)
 *
 * Complexity: O(1)
 */
static double collision_integral_D(double T_star)
{
    if (T_star <= 0.0) return 2.0;
    return 1.06036 * pow(T_star, -0.15610)
         + 0.19300 * exp(-0.47635 * T_star)
         + 1.03587 * exp(-1.52996 * T_star)
         + 1.76474 * exp(-3.89411 * T_star);
}

double chapman_enskog_diffusivity(double T, double P,
                                  const KineticTheoryParams *gasA,
                                  const KineticTheoryParams *gasB)
{
    if (!gasA || !gasB || T <= 0.0 || P <= 0.0) return -1.0;

    /* Combining rules for Lennard-Jones parameters */
    double sigma_AB = 0.5 * (gasA->sigma + gasB->sigma);
    double eps_AB = sqrt(gasA->epsilon_kB * gasB->epsilon_kB);
    double M_AB = 2.0 / (1.0 / gasA->M + 1.0 / gasB->M);

    double T_star = T / eps_AB;
    double omega_D = collision_integral_D(T_star);

    /* Full Chapman-Enskog binary diffusion coefficient */
    double num = (3.0 / 16.0) * sqrt(2.0 * PI * KB * KB * KB * T * T * T / M_AB);
    double den = P * PI * sigma_AB * sigma_AB * omega_D;
    return num / den;
}

/* ─── L2: Wilke-Chang Liquid Diffusivity ───────────────────────── */

/**
 * Wilke-Chang correlation for dilute liquid diffusion.
 *
 * D_AB = 7.4×10⁻⁸ · (φ·M_B)^(1/2) · T / (μ_B · V_A^(0.6))
 *
 * φ = association factor: 2.6 (water), 1.9 (methanol), 1.5 (ethanol), 1.0 (unassociated)
 * V_A = molar volume of solute at normal boiling point [cm³/mol]
 * μ_B = solvent viscosity [cP]
 *
 * Result in cm²/s, converted to m²/s by ×10⁻⁴.
 *
 * Complexity: O(1)
 */
double wilke_chang_diffusivity(double T, double mu_B, double phi_MB,
                               double V_A)
{
    if (T <= 0.0 || mu_B <= 0.0 || V_A <= 0.0) return -1.0;
    /* μ_B expected in Pa·s, convert to cP: 1 Pa·s = 1000 cP */
    double mu_cP = mu_B * 1000.0;
    /* D in cm²/s */
    double D_cm2s = 7.4e-8 * sqrt(phi_MB) * T / (mu_cP * pow(V_A, 0.6));
    /* Convert to m²/s */
    return D_cm2s * 1.0e-4;
}

/* ─── L2: Fuller Diffusivity ───────────────────────────────────── */

/**
 * Fuller-Schettler-Giddings correlation for gas-phase diffusion.
 *
 * D_AB [cm²/s] = 1.0×10⁻³ · T^(1.75) · sqrt(1/M_A + 1/M_B)
 *                 / [P_atm · (Σv_A^(1/3) + Σv_B^(1/3))²]
 *
 * Atomic diffusion volumes Σv:
 *   C=15.9, H=2.31, O=6.11, N=4.54, Cl=21.0, S=22.9, Air=19.7
 *
 * Complexity: O(1)
 */
double fuller_diffusivity(double T, double P,
                          double M_A, double M_B,
                          double sum_v_A, double sum_v_B)
{
    if (T <= 0.0 || P <= 0.0 || M_A <= 0.0 || M_B <= 0.0) return -1.0;
    double P_atm = P / 101325.0;  /* Convert Pa to atm */
    double term1 = 1.0e-3 * pow(T, 1.75) * sqrt(1.0 / M_A + 1.0 / M_B);
    double sum_v = pow(sum_v_A, 1.0 / 3.0) + pow(sum_v_B, 1.0 / 3.0);
    double D_cm2s = term1 / (P_atm * sum_v * sum_v);
    return D_cm2s * 1.0e-4;  /* cm²/s → m²/s */
}

/* ─── L1: Diffusivity Definitions ───────────────────────────────── */

double thermal_diffusivity(double k, double rho, double cp)
{
    if (rho <= 0.0 || cp <= 0.0) return -1.0;
    return k / (rho * cp);
}

double momentum_diffusivity(double mu, double rho)
{
    if (rho <= 0.0) return -1.0;
    return mu / rho;
}

/* ─── L3: Engineering Benchmarks ────────────────────────────────── */

/**
 * Air viscosity: Sutherland's law with standard parameters.
 * μ₀ = 1.716×10⁻⁵ Pa·s at T₀ = 273.15 K, S = 110.4 K
 */
double air_viscosity(double T)
{
    SutherlandParams air = {1.716e-5, 273.15, 110.4};
    return sutherland_viscosity(T, &air);
}

/**
 * Water viscosity: empirical fit.
 * μ(T) = A · 10^(B/(T-C))
 * A = 2.414×10⁻⁵ Pa·s, B = 247.8 K, C = 140 K  (Vogel-Fulcher-Tammann)
 *
 * For simplicity and typical accuracy (273-373K), use power-law fit:
 * μ ≈ 1.0×10⁻³ · (T/293.15)^(-1.5)   [Pa·s]
 */
double water_viscosity(double T)
{
    if (T <= 0.0) return -1.0;
    /* VFT parameters optimized for water in range 273-373K */
    double A = 2.414e-5;
    double B = 247.8;
    double C = 140.0;
    return A * pow(10.0, B / (T - C));
}

/**
 * Air thermal conductivity: empirical fit.
 * k ≈ 1.52×10⁻⁴ · T^(0.5) / (1 + 245/T · 10^(-12/T))
 * Approximate: k_air(T) ≈ k_0 · (T/T_0)^(0.82)
 */
double air_thermal_conductivity(double T)
{
    if (T <= 0.0) return -1.0;
    /* Power-law fit: k ≈ 0.026 W/(m·K) at 300K */
    double k0 = 0.02624;
    double T0 = 300.0;
    return k0 * pow(T / T0, 0.82);
}

/**
 * Water thermal conductivity: empirical fit.
 * k ≈ 0.6 at 300K, increases slowly with T (to ~0.68 at 373K).
 * k(T) = 0.6065 + 0.00164·(T-273.15) - 7.3×10⁻⁶·(T-273.15)²  [W/(m·K)]
 */
double water_thermal_conductivity(double T)
{
    if (T <= 0.0) return -1.0;
    double T_C = T - 273.15;  /* Convert to Celsius */
    return 0.6065 + 0.00164 * T_C - 7.3e-6 * T_C * T_C;
}

/**
 * Water vapor in air diffusivity.
 * Uses Fuller correlation.
 * M in g/mol: M_H2O = 18.0, M_air = 29.0,
 * Σv_H2O = 12.7, Σv_air = 20.1
 *
 * Benchmark: D ≈ 2.6×10⁻⁵ m²/s at 298K, 1 atm.
 */
double water_vapor_air_diffusivity(double T, double P)
{
    /* Fuller correlation expects M in g/mol, returns D in cm²/s.
     * We convert to m²/s by multiplying by 1e-4.
     * D = 1.0e-3 * T^1.75 * sqrt(1/M_A + 1/M_B) / [P_atm * (Σv_A^{1/3} + Σv_B^{1/3})^2]
     * where M in g/mol, Σv are atomic diffusion volumes (dimensionless).
     */
    double M_A = 18.0;     /* H2O g/mol */
    double M_B = 29.0;     /* air g/mol */
    double Sv_A = 12.7;    /* water */
    double Sv_B = 20.1;    /* air */
    double P_atm = P / 101325.0;

    double term1 = 1.0e-3 * pow(T, 1.75) * sqrt(1.0/M_A + 1.0/M_B);
    double sum_v_cuberoot = pow(Sv_A, 1.0/3.0) + pow(Sv_B, 1.0/3.0);
    double D_cm2s = term1 / (P_atm * sum_v_cuberoot * sum_v_cuberoot);
    return D_cm2s * 1.0e-4;  /* cm²/s → m²/s */
}

/* ─── L4: Conservation Laws ──────────────────────────────────────── */

double newtons_law_shear_stress(double mu, double dv_dy)
{
    /* τ_yx = -μ · dv_x/dy, return magnitude */
    return mu * fabs(dv_dy);
}

double fouriers_law_heat_flux(double k, double dT_dy)
{
    /* q_y = -k · dT/dy, return magnitude */
    return k * fabs(dT_dy);
}

double ficks_law_mass_flux(double D_AB, double dC_dy)
{
    /* J_Ay = -D_AB · dC_A/dy, return magnitude */
    return D_AB * fabs(dC_dy);
}

double generalized_transport_flux(double transport_coeff,
                                  double driving_force_gradient)
{
    return transport_coeff * fabs(driving_force_gradient);
}

/* ─── Mixture Rules ─────────────────────────────────────────────── */

/**
 * Wilke interaction parameter Φ_ij.
 *
 * Φ_ij = [1 + (μ_i/μ_j)^(1/2)·(M_j/M_i)^(1/4)]²
 *        / [8·(1 + M_i/M_j)]^(1/2)
 *
 * This quantifies the deviation from linear mixing.
 * For identical species, Φ_ii = 1.
 *
 * Complexity: O(1)
 */
double wilke_interaction_parameter(double mu_i, double mu_j,
                                   double M_i, double M_j)
{
    if (M_i <= 0.0 || M_j <= 0.0) return -1.0;
    double term1 = 1.0 + sqrt(mu_i / mu_j) * pow(M_j / M_i, 0.25);
    double term2 = sqrt(8.0 * (1.0 + M_i / M_j));
    return (term1 * term1) / term2;
}

/**
 * Wilke's mixing rule for gas viscosity.
 *
 * μ_mix = Σ [x_i · μ_i / Σ_j (x_j · Φ_ij)]
 *
 * Derived from the Stefan-Maxwell equations for multicomponent
 * gas mixtures in the low-density limit.
 *
 * Complexity: O(n²)
 */
double wilke_mixture_viscosity(size_t n, const double *x, const double *mu,
                               const double *M)
{
    if (!x || !mu || !M || n == 0) return -1.0;
    double mu_mix = 0.0;
    for (size_t i = 0; i < n; i++) {
        double denominator = 0.0;
        for (size_t j = 0; j < n; j++) {
            double phi = wilke_interaction_parameter(mu[i], mu[j], M[i], M[j]);
            denominator += x[j] * phi;
        }
        if (denominator > 0.0) {
            mu_mix += x[i] * mu[i] / denominator;
        }
    }
    return mu_mix;
}

/**
 * Mason-Saxena mixing rule for gas mixture thermal conductivity.
 *
 * k_mix = Σ [x_i · k_i / Σ_j (x_j · Φ_ij)]  (same form as Wilke)
 * but Φ_ij = ε · Φ_ij_Wilke where ε ≈ 1.065 for most gases.
 *
 * Complexity: O(n²)
 */
double mason_saxena_mixture_conductivity(size_t n, const double *x,
                                         const double *k, const double *M)
{
    if (!x || !k || !M || n == 0) return -1.0;
    double k_mix = 0.0;

    /* For thermal conductivity, Φ_ij uses the same form as viscosity
     * but with an empirical factor ε ≈ 1.065 (Mason & Saxena, 1958) */
    for (size_t i = 0; i < n; i++) {
        double denominator = 0.0;
        for (size_t j = 0; j < n; j++) {
            /* Approximate: use viscosity ratio proxy from k-ratio
             * (both scale similarly with sqrt(T) and 1/σ²) */
            double mu_ratio = k[i] / k[j];
            double phi = 1.065 * wilke_interaction_parameter(mu_ratio, 1.0,
                                                              M[i], M[j]);
            denominator += x[j] * phi;
        }
        if (denominator > 0.0) {
            k_mix += x[i] * k[i] / denominator;
        }
    }
    return k_mix;
}

/* ─── Transport State Utilities ─────────────────────────────────── */

void compute_transport_state(const PowerLawViscosityParams *visc_params,
                             double k, double D, double rho,
                             double cp, double T, double P,
                             TransportState *state)
{
    if (!state) return;
    state->mu = power_law_viscosity(T, visc_params);
    state->k = k;
    state->D = D;
    state->rho = rho;
    state->cp = cp;
    state->T = T;
    state->P = P;
}

int validate_transport_state(const TransportState *state)
{
    if (!state) return -1;
    if (state->mu <= 0.0)  return 1;
    if (state->k  <= 0.0)  return 2;
    if (state->D  <= 0.0)  return 3;
    if (state->rho <= 0.0) return 4;
    if (state->cp  <= 0.0) return 5;
    if (state->T   <= 0.0) return 6;
    if (state->P   <= 0.0) return 7;
    return 0;
}

void print_transport_summary(const TransportState *state)
{
    if (!state) {
        printf("TransportState: NULL\n");
        return;
    }
    printf("━ Transport Properties ━━━━━━━━━━━━━━━━━━━━━━━━\n");
    printf("  Temperature:  %.2f K (%.2f °C)\n",
           state->T, state->T - 273.15);
    printf("  Pressure:     %.3f kPa (%.4f atm)\n",
           state->P / 1000.0, state->P / 101325.0);
    printf("  Density:      %.4f kg/m³\n", state->rho);
    printf("  Cp:           %.1f J/(kg·K)\n", state->cp);
    printf("  Viscosity μ:  %.3e Pa·s\n", state->mu);
    printf("  Therm. Cond k: %.4f W/(m·K)\n", state->k);
    printf("  Diffusivity D: %.3e m²/s\n", state->D);
    printf("  Kinematic ν:  %.3e m²/s\n",
           momentum_diffusivity(state->mu, state->rho));
    printf("  Thermal α:    %.3e m²/s\n",
           thermal_diffusivity(state->k, state->rho, state->cp));
    printf("  Prandtl Pr:   %.3f\n",
           state->mu * state->cp / state->k);
    printf("  Schmidt Sc:   %.3f\n",
           state->mu / (state->rho * state->D));
    printf("  Lewis Le:     %.3f\n",
           state->k / (state->rho * state->cp * state->D));
    printf("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");
}

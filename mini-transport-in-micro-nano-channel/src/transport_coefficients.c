#include "micro_nano_transport.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * transport_coefficients.c - Transport Coefficients for Micro/Nano Systems
 *
 * Computes viscosity, thermal conductivity, and mass diffusivity
 * from kinetic theory and empirical correlations. These coefficients
 * are essential inputs to all continuum transport equations.
 *
 * Knowledge: L1 (definitions), L3 (engineering quantities), L5 (computation)
 * Reference: Chapman & Cowling (1970), Hirschfelder et al. (1954)
 */

/*
 * compute_kinetic_fluid_props - L3 Engineering Quantities
 *
 * Computes fluid transport properties from hard-sphere kinetic theory.
 *
 * Dynamic viscosity (Chapman-Enskog for hard spheres):
 *   eta = (5/16) * sqrt(pi*m*kB*T) / (pi*d^2)
 *
 * Thermal conductivity for monatomic gas (modified Eucken):
 *   k = (15/4) * (R/M) * eta
 *
 * Self-diffusion coefficient:
 *   D = (3/8) * sqrt(pi*(kB*T)^3/m) / (pi*d^2*p)
 *
 * These are the simplest models that capture the correct
 * temperature and pressure scaling.
 */
void compute_kinetic_fluid_props(double T, double p, double m_mol,
                                  double d_col, FluidProperties *out)
{
    if (out == NULL || T <= 0.0 || p <= 0.0 || m_mol <= 0.0
        || d_col <= 0.0) return;

    double n = p / (KB * T);
    double d2 = d_col * d_col;
    double m = m_mol;

    /* Viscosity: eta ~ sqrt(T), independent of p (Maxwell 1860) */
    double eta = (5.0 / 16.0) * sqrt(M_PI * m * KB * T) / (M_PI * d2);
    out->dynamic_viscosity = eta;
    out->molecular_mass = m;
    out->collision_diameter = d_col;
    out->temperature = T;

    /* Density from ideal gas law */
    out->density = p * m / (KB * T);
    out->kinematic_viscosity = eta / out->density;

    /* Thermal conductivity for monatomic gas */
    double cv = 1.5 * R_GAS / (m * NA); /* per unit mass */
    double gamma = 5.0 / 3.0;
    out->thermal_conductivity = (15.0 / 4.0) * (R_GAS / (m * NA)) * eta;
    out->specific_heat_cp = gamma * R_GAS / (m * NA);

    /* Thermal diffusivity */
    out->thermal_diffusivity = out->thermal_conductivity
                               / (out->density * out->specific_heat_cp);

    /* Mass diffusivity (self-diffusion) */
    out->mass_diffusivity = (3.0 / 8.0)
        * sqrt(M_PI * pow(KB * T, 3.0) / m) / (M_PI * d2 * p);

    /* Dimensionless numbers */
    out->prandtl_number = out->kinematic_viscosity / out->thermal_diffusivity;
    out->schmidt_number = out->kinematic_viscosity / out->mass_diffusivity;
    out->relative_permittivity = 1.0;
}

/*
 * compute_viscosity_sutherland - L3
 *
 * Sutherland's law for gas viscosity temperature dependence:
 *
 *   eta(T) = eta_0 * (T/T_0)^(3/2) * (T_0 + S) / (T + S)
 *
 * where S is the Sutherland constant (characteristic temperature).
 * This is more accurate than the hard-sphere model for real gases.
 *
 * Typical values for air: eta_0 = 1.716e-5 Pa*s, T_0 = 273.15 K, S = 110.4 K
 */
double compute_viscosity_sutherland(double T, double eta_0, double T_0,
                                     double S)
{
    if (T <= 0.0 || T_0 <= 0.0 || eta_0 <= 0.0) return -1.0;
    return eta_0 * pow(T / T_0, 1.5) * (T_0 + S) / (T + S);
}

/*
 * compute_viscosity_power_law - L3
 *
 * Simple power-law temperature scaling for viscosity:
 *
 *   eta(T) = eta_0 * (T/T_0)^omega
 *
 * where omega is the viscosity temperature exponent.
 * omega = 0.5 for hard spheres, 0.74 for VHS model of air.
 */
double compute_viscosity_power_law(double T, double eta_0, double T_0,
                                    double omega)
{
    if (T <= 0.0 || T_0 <= 0.0 || eta_0 <= 0.0) return -1.0;
    return eta_0 * pow(T / T_0, omega);
}

/*
 * compute_thermal_conductivity_eucken - L3
 *
 * Modified Eucken correlation for polyatomic gas thermal conductivity:
 *
 *   k = eta * (c_p + 1.25 * R/M)
 *
 * where c_p is the specific heat at constant pressure.
 * For monatomic gases: c_p = (5/2)*(R/M), giving k = (15/4)*(R/M)*eta.
 */
double compute_thermal_conductivity_eucken(double eta, double cp,
                                            double molar_mass)
{
    if (eta <= 0.0 || cp <= 0.0 || molar_mass <= 0.0) return -1.0;
    return eta * (cp + 1.25 * R_GAS / molar_mass);
}

/*
 * compute_binary_diffusion_coefficient - L3
 *
 * Binary diffusion coefficient from Chapman-Enskog:
 *
 *   D_AB = (3/16) * sqrt(2*pi*(kB*T)^3*(1/m_A+1/m_B)) / (p*pi*sigma_AB^2*Omega_D)
 *
 * Simplified form with Omega_D ~ 1 (hard sphere approximation):
 *   D_AB = 1.858e-7 * T^(3/2) * sqrt(1/M_A + 1/M_B) / (p * sigma_AB^2)
 *
 * where D_AB is in m^2/s, p in atm, sigma in Angstroms.
 */
double compute_binary_diffusion_coefficient(double T, double p_atm,
    double M_A, double M_B, double sigma_AB_angstrom)
{
    if (T <= 0.0 || p_atm <= 0.0 || M_A <= 0.0 || M_B <= 0.0
        || sigma_AB_angstrom <= 0.0) return -1.0;

    double sigma_m = sigma_AB_angstrom * 1.0e-10;
    double p_Pa = p_atm * 101325.0;
    double reduced_mass = (M_A * M_B) / (M_A + M_B);
    double factor = 1.858e-7;

    return factor * pow(T, 1.5) * sqrt(1.0 / reduced_mass)
           / (p_atm * sigma_AB_angstrom * sigma_AB_angstrom)
           * 1.0e-4; /* convert cm^2/s to m^2/s */
}

/*
 * compute_transport_cross_section - L3
 *
 * Effective collision cross-section for transport properties:
 *
 *   sigma_transport = pi * d^2 * Omega(T*)
 *
 * where Omega(T*) is the collision integral, a function of
 * the reduced temperature T* = kB*T/epsilon (epsilon = well depth).
 *
 * Omega ~ 1 at high T (hard sphere limit)
 * Omega > 1 at low T (long-range attraction increases cross-section)
 */
double compute_transport_cross_section(double d_col, double T,
                                        double epsilon_kB)
{
    if (d_col <= 0.0 || T <= 0.0 || epsilon_kB <= 0.0) return -1.0;

    double T_star = T / epsilon_kB;
    double Omega;
    if (T_star < 1.0) {
        Omega = 1.0 / T_star * 2.0;
    } else if (T_star < 10.0) {
        Omega = 1.0 + 0.2 / sqrt(T_star);
    } else {
        Omega = 1.0;
    }

    return M_PI * d_col * d_col * Omega;
}

/*
 * compute_prandtl_number_estimate - L3
 *
 * Prandtl number estimation from kinetic theory:
 *
 *   Pr = c_p * eta / k ~ (4*gamma) / (9*gamma - 5)
 *
 * For monatomic gas (gamma=5/3): Pr = 2/3 = 0.667
 * For diatomic gas (gamma=7/5): Pr ~ 0.72
 * For polyatomic (gamma~1.3): Pr ~ 0.78
 *
 * The Eucken formula gives the same result.
 */
double compute_prandtl_number_estimate(double gamma)
{
    if (gamma <= 1.0) return -1.0;
    return (4.0 * gamma) / (9.0 * gamma - 5.0);
}

/*
 * compute_schmidt_number_estimate - L3
 *
 * Schmidt number estimation:
 *
 *   Sc = nu / D = eta / (rho * D)
 *
 * For gases at STP, Sc is typically 0.7-2.0.
 * Sc = 1 implies equal momentum and mass diffusion rates
 *     (Reynolds analogy holds).
 */
double compute_schmidt_number_estimate(double eta, double rho, double D)
{
    if (rho <= 0.0 || D <= 0.0) return -1.0;
    return eta / (rho * D);
}

/*
 * compute_lewis_number - L3
 *
 * Lewis number:
 *
 *   Le = alpha / D = Sc / Pr
 *
 * Compares thermal and mass diffusion rates.
 * Le < 1: mass diffuses faster than heat (most gases)
 * Le > 1: heat diffuses faster than mass
 * Le = 1: equal rates (common simplifying assumption)
 */
double compute_lewis_number(double alpha, double D)
{
    if (D <= 0.0) return -1.0;
    return alpha / D;
}

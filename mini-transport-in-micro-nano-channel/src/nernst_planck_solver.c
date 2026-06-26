#include "micro_nano_transport.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * nernst_planck_solver.c - Nernst-Planck Equation Solvers
 *
 * The Nernst-Planck equation describes the flux of ionic species i
 * under combined concentration gradient (diffusion) and electric
 * potential gradient (migration):
 *
 *   J_i = -D_i * [dc_i/dx + (z_i*e/(kB*T)) * c_i * d(phi)/dx]
 *
 * Coupled with the Poisson equation, this forms the Poisson-Nernst-Planck
 * (PNP) system, the fundamental model for ionic transport.
 *
 * Knowledge: L4 (Nernst-Planck equation), L5 (PNP solvers)
 * Reference: Newman & Thomas-Alyea (2004) Electrochemical Systems
 *            Zheng & Wei (2011) Microfluid. Nanofluid. 10, 211
 */

/*
 * compute_nernst_planck_steady_1d - L4-L5
 *
 * Steady-state 1D Nernst-Planck flux for a single ion species:
 *
 *   J = -D * [dc/dx + (z*e/(kB*T)) * c * d(phi)/dx]
 *
 * Using the Goldman-Hodgkin-Katz constant field approximation
 * (d(phi)/dx = -V_applied/L = constant):
 *
 *   J = D * z*e*V / (kB*T*L) * [c_0*exp(z*e*V/(kB*T)) - c_L] / [exp(z*e*V/(kB*T)) - 1]
 *
 * This is the formula used to model ion transport through
 * membrane channels and nanopores.
 *
 * Complexity: O(n) for profile computation
 */
void compute_nernst_planck_steady_1d(const double *c, const double *phi,
    double D, double z, double T, double L, int n, double *flux)
{
    if (c == NULL || phi == NULL || flux == NULL || n < 2
        || D <= 0.0 || T <= 0.0 || L <= 0.0) return;

    double beta = E_CHARGE / (KB * T);
    double V = phi[n-1] - phi[0];
    double c0 = c[0];
    double cL = c[n-1];

    /* Goldman-Hodgkin-Katz flux */
    double arg = fabs(z) * beta * V;
    if (fabs(arg) < 1.0e-10) {
        /* Near-equilibrium: flux = D * (c0 - cL) / L (pure diffusion) */
        *flux = D * (c0 - cL) / L;
    } else {
        double exp_arg = exp(z * beta * V);
        if (fabs(exp_arg - 1.0) < 1.0e-15) {
            *flux = D * (c0 - cL) / L;
        } else {
            *flux = D * z * beta * V / L
                    * (c0 * exp_arg - cL) / (exp_arg - 1.0);
        }
    }
}

/*
 * compute_diffusion_potential - L4
 *
 * Diffusion (liquid junction) potential between two solutions
 * of different concentrations:
 *
 *   V_diff = (kB*T/e) * (D_+ - D_-) / (D_+ + D_-) * ln(c_hi / c_lo)
 *
 * This potential arises when cations and anions have different
 * diffusion coefficients, creating a charge separation.
 *
 * For KCl (D_+/D_- ~ 1): V_diff ~ 0 (used in salt bridges)
 * For HCl (D_+/D_- ~ 5): V_diff ~ 38 mV for 10x concentration ratio
 */
void compute_diffusion_potential(double c_hi, double c_lo, double z,
    double D_plus, double D_minus, double T, double *V_diff)
{
    if (V_diff == NULL || c_hi <= 0.0 || c_lo <= 0.0
        || T <= 0.0 || (D_plus + D_minus) <= 0.0) {
        if (V_diff) *V_diff = 0.0;
        return;
    }

    double mobility_diff = (D_plus - D_minus) / (D_plus + D_minus);
    double conc_ratio = c_hi / c_lo;

    *V_diff = (KB * T / (fabs(z) * E_CHARGE))
              * mobility_diff * log(conc_ratio);
}

/*
 * compute_nernst_equilibrium_potential - L4
 *
 * Nernst equilibrium potential for a single ion species
 * across a selectively permeable membrane:
 *
 *   E_eq = (RT/(zF)) * ln(c_out / c_in)
 *
 * where F = NA * e is Faraday's constant.
 *
 * This is the voltage required to maintain zero net flux
 * when there is a concentration difference across the membrane.
 *
 * For K+ at 37C with [K+]_in/[K+]_out = 140/5 mM:
 *   E_K ~ -89 mV (the resting potential of neurons)
 */
double compute_nernst_equilibrium_potential(double c_in, double c_out,
    double z, double T)
{
    if (c_in <= 0.0 || c_out <= 0.0 || T <= 0.0 || z == 0.0) return 0.0;
    double RT_over_zF = KB * T / (z * E_CHARGE);
    return RT_over_zF * log(c_out / c_in);
}

/*
 * compute_ghk_current_voltage - L5
 *
 * Goldman-Hodgkin-Katz (GHK) current equation for a membrane
 * permeable to multiple ion species:
 *
 *   I = V * F^2/(R*T) * SUM[ P_i * z_i^2 * (c_i_in - c_i_out*exp(-z_i*F*V/(R*T)))
 *                            / (1 - exp(-z_i*F*V/(R*T))) ]
 *
 * where P_i = D_i / L is the permeability of species i.
 *
 * This is the basis of the Hodgkin-Huxley model of neural
 * action potentials, and also applies to synthetic nanopores.
 *
 * Complexity: O(n_species)
 */
double compute_ghk_current_voltage(int n_species, const double *P,
    const int *z, const double *c_in, const double *c_out,
    double V, double T)
{
    if (P == NULL || z == NULL || c_in == NULL || c_out == NULL
        || n_species <= 0 || T <= 0.0) return 0.0;

    double I_total = 0.0;
    double F = NA * E_CHARGE;
    double F_over_RT = F / (R_GAS * T);

    for (int i = 0; i < n_species; i++) {
        double zi = (double)z[i];
        double exp_term = exp(-zi * F_over_RT * V);
        double numerator = P[i] * zi * zi
            * (c_in[i] - c_out[i] * exp_term);
        double denominator = 1.0 - exp_term;

        if (fabs(denominator) < 1.0e-15) {
            /* L'Hopital: I_i = P_i * z_i * F * c_in */
            I_total += P[i] * zi * F * c_in[i];
        } else {
            I_total += numerator / denominator;
        }
    }

    return V * F * F_over_RT * I_total;
}

/*
 * compute_ion_flux_fickian - L4
 *
 * Pure Fickian diffusion flux (no electric field):
 *
 *   J_diff = -D * (c_L - c_0) / L
 *
 * This is the limiting case when the electric field vanishes
 * (e.g., for neutral species or at very high salt concentration
 * where migration effects are screened).
 */
double compute_ion_flux_fickian(double D, double c0, double cL, double L)
{
    if (D <= 0.0 || L <= 0.0) return 0.0;
    return -D * (cL - c0) / L;
}

/*
 * compute_ion_flux_migrative - L4
 *
 * Pure migrative flux (no concentration gradient):
 *
 *   J_mig = -(z*e*D/(kB*T)) * c * d(phi)/dx
 *
 * This is the electrophoretic contribution, dominant when
 * concentration gradients are negligible (e.g., in buffered
 * solutions with strong electric fields).
 */
double compute_ion_flux_migrative(double D, double z, double c,
    double T, double E_field)
{
    if (D <= 0.0 || T <= 0.0) return 0.0;
    return -(z * E_CHARGE * D / (KB * T)) * c * E_field;
}

/*
 * compute_transference_number - L3
 *
 * Transference (transport) number for ion species i:
 *
 *   t_i = |z_i| * u_i * c_i / SUM(|z_j| * u_j * c_j)
 *
 * where u_i = |z_i|*e*D_i/(kB*T) is the ion mobility
 * (from the Einstein relation).
 *
 * t_i represents the fraction of total current carried by
 * species i. In a 1:1 electrolyte:
 *   t_+ ~ 0.4, t_- ~ 0.6 for NaCl (Cl- is more mobile)
 *   t_+ ~ 0.5, t_- ~ 0.5 for KCl (equal mobility)
 */
double compute_transference_number(int n_species, const double *D,
    const int *z, const double *c, int species_idx)
{
    if (D == NULL || z == NULL || c == NULL || n_species <= 0
        || species_idx < 0 || species_idx >= n_species) return -1.0;

    double numerator = fabs((double)z[species_idx]) * D[species_idx]
                       * c[species_idx];
    double denominator = 0.0;
    for (int i = 0; i < n_species; i++) {
        denominator += fabs((double)z[i]) * D[i] * c[i];
    }

    if (denominator <= 0.0) return 0.0;
    return numerator / denominator;
}

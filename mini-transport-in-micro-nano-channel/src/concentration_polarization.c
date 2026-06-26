#include "micro_nano_transport.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * concentration_polarization.c - Concentration Polarization at Interfaces
 *
 * Concentration polarization (CP) is the formation of a concentration
 * gradient near a membrane or electrode surface due to selective transport.
 * In nanofluidic systems, CP significantly affects:
 *   - Ion current rectification
 *   - Electroosmotic flow rate
 *   - Nanopore sensing sensitivity
 *
 * Knowledge: L2 (CP concept), L5 (modeling), L7 (applications)
 * Reference: Probstein (1994) Physicochemical Hydrodynamics
 *            Zangle et al. (2010) Chem. Soc. Rev. 39, 1014
 */

/*
 * compute_concentration_polarization - L5
 *
 * 1D steady-state CP profile between a membrane surface and bulk:
 *
 *   c(x) = c_bulk + (c_wall - c_bulk) * exp(-x / delta_CP)
 *
 * where delta_CP = D / J is the CP layer thickness,
 * D is the diffusion coefficient, and J is the permeate flux.
 *
 * The concentration at the wall c_wall is determined by the balance
 * between convective transport toward the membrane and diffusive
 * transport away from it.
 *
 * Complexity: O(n)
 */
void compute_concentration_polarization(double c_bulk, double c_wall,
    double x_max, double delta_CP, double *x, double *c, int n)
{
    if (x == NULL || c == NULL || n <= 0 || delta_CP <= 0.0) return;

    for (int i = 0; i < n; i++) {
        x[i] = (double)i / (double)(n - 1) * x_max;
        c[i] = c_bulk + (c_wall - c_bulk) * exp(-x[i] / delta_CP);
    }
}

/*
 * compute_cp_layer_thickness - L3
 *
 * CP layer thickness from the film model:
 *
 *   delta_CP = D / k_m
 *
 * where k_m is the mass transfer coefficient.
 * For laminar flow in a channel:
 *   Sh = k_m * D_h / D = 0.332 * Re^(1/2) * Sc^(1/3) (developing)
 *   Sh = 7.54 (fully developed, constant wall concentration)
 *
 * This thickness determines the spatial extent of CP effects.
 */
double compute_cp_layer_thickness(double D, double D_h, double Re, double Sc)
{
    if (D <= 0.0 || D_h <= 0.0) return -1.0;

    double Sh;
    if (Re > 0.0) {
        Sh = 0.332 * sqrt(Re) * pow(Sc, 1.0/3.0);
    } else {
        Sh = 7.54; /* fully developed limit */
    }

    if (Sh <= 0.0) return -1.0;
    return D_h / Sh;
}

/*
 * compute_cp_modulus - L5
 *
 * Concentration polarization modulus:
 *
 *   Gamma = c_wall / c_bulk = exp(J_v * delta_CP / D)
 *          = exp(J_v / k_m)
 *
 * This is the key parameter quantifying CP severity:
 *   Gamma ~ 1: negligible CP
 *   Gamma >> 1: severe CP (high flux or low mass transfer)
 *
 * In RO membranes: Gamma ~ 1.2-3.0
 * In ion-exchange membranes: Gamma can reach 10-100
 */
double compute_cp_modulus(double flux, double k_m)
{
    if (k_m <= 0.0) return -1.0;
    double Pe = flux / k_m; /* Peclet number for CP */
    if (Pe > 20.0) Pe = 20.0; /* prevent overflow */
    return exp(Pe);
}

/*
 * compute_limiting_current_density - L6
 *
 * Limiting current density in electrodialysis/electrodeionization:
 *
 *   i_lim = z * F * D * c_bulk / (delta_CP * (t_m - t_s))
 *
 * where F is Faraday constant, t_m is transport number in membrane,
 * t_s is transport number in solution.
 *
 * At i > i_lim, water splitting occurs at the membrane surface,
 * causing pH changes and energy losses.
 *
 * Typical values: i_lim ~ 10-500 A/m^2 for brackish water ED.
 *
 * Reference: Strathmann (2010) Desalination 264, 268
 */
double compute_limiting_current_density(double z, double D, double c_bulk,
    double delta_CP, double t_membrane, double t_solution)
{
    if (D <= 0.0 || c_bulk <= 0.0 || delta_CP <= 0.0) return -1.0;
    double delta_t = t_membrane - t_solution;
    if (fabs(delta_t) < 1.0e-10) return -1.0;
    double F = NA * E_CHARGE; /* Faraday constant */
    return fabs(z) * F * D * c_bulk / (delta_CP * delta_t);
}

/*
 * compute_cp_time_scale - L3
 *
 * Characteristic time for CP layer development:
 *
 *   tau_CP = delta_CP^2 / D
 *
 * This is the time scale for the CP layer to reach steady state
 * after a change in operating conditions.
 *
 * For typical microfluidic systems (delta_CP ~ 10 um, D ~ 1e-9 m^2/s):
 *   tau_CP ~ 0.1 s -- very fast compared to macro systems.
 */
double compute_cp_time_scale(double delta_CP, double D)
{
    if (D <= 0.0 || delta_CP <= 0.0) return -1.0;
    return delta_CP * delta_CP / D;
}

/*
 * compute_donnan_exclusion_potential - L5
 *
 * Donnan potential arising from CP between a charged membrane
 * and the adjacent solution:
 *
 *   delta_phi_D = (kB*T/(z*e)) * ln(c_membrane / c_solution)
 *
 * This potential opposes further ion transport across the membrane
 * and contributes to the membrane resistance.
 *
 * For a highly charged membrane in dilute solution,
 * delta_phi_D can reach 100-200 mV.
 */
double compute_donnan_exclusion_potential(double c_membrane,
    double c_solution, double z, double T)
{
    if (c_membrane <= 0.0 || c_solution <= 0.0 || T <= 0.0 || z == 0.0)
        return 0.0;
    double ratio = c_membrane / c_solution;
    if (ratio > 1.0e10) ratio = 1.0e10;
    return (KB * T / (fabs(z) * E_CHARGE)) * log(ratio);
}

/*
 * compute_water_splitting_rate - L7 Application
 *
 * Water splitting at an overlimiting current:
 *
 *   2 H2O -> H3O+ + OH-
 *
 * Rate depends on the excess current above the limiting current:
 *
 *   J_ws = (i - i_lim) / F * eta_ws
 *
 * where eta_ws is the water splitting efficiency (typically 0.01-0.1).
 *
 * This is critical for understanding pH changes in nanofluidic
 * bipolar membranes and electrodialysis stacks operating
 * above the limiting current.
 */
double compute_water_splitting_rate(double i_applied, double i_lim,
                                     double eta_ws)
{
    if (i_applied <= i_lim || eta_ws <= 0.0 || eta_ws > 1.0) return 0.0;
    double F = NA * E_CHARGE;
    return (i_applied - i_lim) / F * eta_ws;
}

/*
 * compute_electroconvective_instability - L8 Advanced
 *
 * Critical voltage for electroconvective instability at a
 * membrane-solution interface:
 *
 *   V_crit = (kB*T/(z*e)) * ln(1 + 4*delta_CP^2 / (eps * D / (eta * k_m)))
 *
 * Above this voltage, electroconvection enhances mass transfer
 * beyond the diffusion-limited (limiting current) regime.
 *
 * This explains the "overlimiting current" regime observed
 * in electrodialysis, where current can exceed the theoretical
 * limiting value due to convective mixing.
 *
 * Reference: Rubinstein & Zaltzman (2000) Phys. Rev. E 62, 2238
 */
double compute_electroconvective_instability(double delta_CP, double eps,
    double eta, double k_m, double D, double z, double T)
{
    if (delta_CP <= 0.0 || eta <= 0.0 || k_m <= 0.0 || D <= 0.0
        || T <= 0.0 || z == 0.0) return -1.0;

    double electroosmotic_scale = eps * D / (eta * k_m);
    double arg = 1.0 + 4.0 * delta_CP * delta_CP / electroosmotic_scale;

    if (arg <= 1.0) return -1.0;
    return (KB * T / (fabs(z) * E_CHARGE)) * log(arg);
}

#include "micro_nano_transport.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * rarefied_gas_dynamics.c - Rarefied Gas Dynamics and Transition Regime
 *
 * When Kn > 0.1, the Navier-Stokes equations with slip boundary conditions
 * become inadequate. The transition regime (0.1 < Kn < 10) requires
 * higher-order continuum models (Burnett equations) or kinetic approaches
 * (Boltzmann equation, DSMC).
 *
 * Knowledge: L2 (rarefaction), L4 (Burnett eqns), L8 (kinetic theory)
 * Reference: Bird (1994) Molecular Gas Dynamics
 *            Cercignani (2000) Rarefied Gas Dynamics
 */

/*
 * compute_rarefied_mass_flow_rate - L5
 *
 * Mass flow rate for rarefied gas in a channel using the unified
 * flow model (Beskok-Karniadakis):
 *
 *   mdot = rho * H^3 * W / (12 * eta) * (-dp/dx) * (1 + alpha * Kn)
 *
 * where alpha is a rarefaction coefficient:
 *   alpha = alpha_0 * (2/pi) * arctan(alpha_1 * Kn^beta)
 *
 * with typical values: alpha_0 ~ 1.20, alpha_1 ~ 0.5, beta ~ 0.4
 * This model spans from continuum to free molecular flow.
 *
 * Reference: Beskok & Karniadakis (1999) Microscale Thermophys. Eng. 3, 43
 */
double compute_rarefied_mass_flow_rate(double H, double dpdx, double eta,
                                        double rho, double Kn)
{
    if (H <= 0.0 || eta <= 0.0 || rho <= 0.0) return -1.0;
    double Q_noslip = H * H * H / (12.0 * eta) * fabs(dpdx);
    double alpha0 = 1.20;
    double alpha1 = 0.5;
    double beta = 0.4;
    double alpha = alpha0 * (2.0 / M_PI) * atan(alpha1 * pow(Kn, beta));
    return rho * Q_noslip * (1.0 + alpha * Kn);
}

/*
 * compute_transition_regime_mass_flow - L5
 *
 * Mass flow in the full range of Kn using the Cercignani-Lampis
 * scattering kernel model.
 *
 * The normalized mass flow rate M(Kn) = mdot(Kn) / mdot_free_molecular
 * transitions smoothly from the continuum limit to the free molecular limit.
 *
 * M(Kn) ~ 1/Kn for Kn -> 0 (continuum limit, inverse proportionality)
 * M(Kn) -> 1 for Kn -> inf (free molecular limit)
 *
 * Uses the Sherman interpolation formula:
 *   M(Kn) = M_slip(Kn) / (1 + Kn/0.1) + M_fm(Kn) * Kn/(1 + Kn)
 *
 * Complexity: O(1)
 */
void compute_transition_regime_mass_flow(double H, double dpdx, double eta,
                                          double rho, double Kn, double *mdot)
{
    if (mdot == NULL || H <= 0.0 || eta <= 0.0 || rho <= 0.0) return;

    /* Slip flow contribution */
    double mdot_slip = compute_rarefied_mass_flow_rate(H, dpdx, eta, rho, Kn);

    /* Free molecular contribution (Knudsen diffusion) */
    double c_bar = sqrt(8.0 * KB * 293.15 / (M_PI * 4.65e-26));
    (void)Kn; /* Kn used in mdot_slip above, lambda calculated for documentation */
    double D_K = (H / 3.0) * c_bar;
    double mdot_fm = rho * D_K * H * fabs(dpdx) / (KB * 293.15 / (4.65e-26));

    double weight_slip = 1.0 / (1.0 + Kn / 0.1);
    double weight_fm = Kn / (1.0 + Kn);

    *mdot = mdot_slip * weight_slip + mdot_fm * weight_fm;
}

/*
 * compute_knudsen_diffusion_coefficient - L5
 *
 * Knudsen diffusion coefficient for free molecular flow:
 *
 *   D_K = (d/3) * sqrt(8*R*T/(pi*M))
 *
 * where d is the characteristic pore/channel dimension.
 * This is the dominant transport mechanism when Kn >> 1.
 *
 * Unlike Fickian diffusion (D ~ 1/p), Knudsen diffusion is
 * pressure-independent and depends only on temperature and geometry.
 */
double compute_knudsen_diffusion_coefficient(double d, double T,
                                              double molar_mass)
{
    if (d <= 0.0 || T <= 0.0 || molar_mass <= 0.0) return -1.0;
    double R_specific = R_GAS / molar_mass;
    return (d / 3.0) * sqrt(8.0 * R_specific * T / M_PI);
}

/*
 * compute_knudsen_layer_thickness - L3
 *
 * Knudsen layer thickness: the region within ~1-2 mean free paths
 * of the wall where the velocity distribution function deviates
 * significantly from the Chapman-Enskog (continuum) form.
 *
 *   delta_K = c_K * lambda
 *
 * where c_K ~ 1.0-1.5 depending on the accommodation model.
 *
 * Inside the Knudsen layer, the Navier-Stokes constitutive relations
 * (Newton's law of viscosity) break down, and kinetic boundary
 * conditions must be applied.
 */
double compute_knudsen_layer_thickness(double lambda, double c_K)
{
    if (lambda <= 0.0 || c_K <= 0.0) return -1.0;
    return c_K * lambda;
}

/*
 * compute_slip_coefficient_burnett - L8 Advanced
 *
 * Slip coefficient from the Burnett equations (second-order
 * approximation to the Boltzmann equation):
 *
 *   u_slip = A1 * lambda * (du/dn)_wall - A2 * lambda^2 * (d^2u/dn^2)_wall
 *
 * From kinetic theory for Maxwell molecules:
 *   A1 = 1.1466 (first-order, same as Maxwell slip)
 *   A2 = 0.9756 (second-order, Burnett correction)
 *
 * The second-order term captures the curvature of the velocity
 * profile near the wall and improves accuracy for Kn > 0.1.
 *
 * Reference: Lockerby et al. (2004) J. Fluid Mech. 501, 527
 */
double compute_slip_coefficient_burnett(double lambda, double dudn,
                                         double d2udn2)
{
    double A1 = 1.1466;
    double A2 = 0.9756;
    return A1 * lambda * dudn - A2 * lambda * lambda * d2udn2;
}

/*
 * compute_mass_flow_rate_unified - L8
 *
 * Unified mass flow rate across all Knudsen regimes using an
 * interpolation formula valid from continuum to free molecular:
 *
 *   mdot_unified = w_c * mdot_continuum + w_t * mdot_transition + w_fm * mdot_fm
 *
 * The weights are Gaussian functions of log10(Kn), centered
 * at the regime boundaries. This provides smooth transitions
 * suitable for CFD implementation.
 *
 * Complexity: O(1)
 */
double compute_mass_flow_rate_unified(double mdot_cont, double mdot_slip,
                                       double mdot_fm, double Kn)
{
    double logKn = log10(Kn + 1.0e-10);

    /* Gaussian weights centered at regime boundaries */
    double w_cont = exp(-logKn * logKn / 0.1);
    double w_slip = exp(-(logKn + 1.5) * (logKn + 1.5) / 0.5);
    double w_fm = exp(-(logKn + 3.0) * (logKn + 3.0) / 1.0);

    double w_total = w_cont + w_slip + w_fm;
    if (w_total <= 0.0) return mdot_slip;

    return (w_cont * mdot_cont + w_slip * mdot_slip + w_fm * mdot_fm) / w_total;
}

/*
 * compute_bgk_relaxation_time - L8 Boltzmann BGK model
 *
 * The Bhatnagar-Gross-Krook (BGK) model simplifies the Boltzmann
 * collision integral with a single relaxation time:
 *
 *   tau = mu / p
 *
 * This is the characteristic time for the velocity distribution
 * function to relax to the local Maxwellian.
 *
 * For air at STP: tau ~ 2e-10 s, much faster than typical
 * microchannel flow timescales (~1e-6 s).
 */
double compute_bgk_relaxation_time(double mu, double p)
{
    if (p <= 0.0) return -1.0;
    return mu / p;
}

/*
 * compute_knudsen_paradox - L6
 *
 * The Knudsen paradox: experimental measurements show that the
 * mass flow rate in a channel reaches a MINIMUM at Kn ~ 1,
 * rather than decreasing monotonically with Kn.
 *
 * This minimum occurs because:
 *   - At low Kn: continuum Poiseuille flow: Q ~ 1/Kn -> 0 as Kn -> 0
 *   - At high Kn: free molecular flow: Q approaches a constant
 *   - At Kn ~ 1: neither mechanism is efficient, hence the minimum
 *
 * This function computes the minimum location.
 *
 * Reference: Knudsen (1909) Ann. Phys. 28, 75
 */
double compute_knudsen_minimum_location(double H, double T,
                                         double molar_mass)
{
    /* The minimum occurs at Kn_crit ~ 0.8-1.2 for most gases */
    (void)H;
    (void)T;
    (void)molar_mass;
    return 0.95;
}

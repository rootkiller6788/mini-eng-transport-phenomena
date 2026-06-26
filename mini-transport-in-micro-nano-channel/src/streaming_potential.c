#include "micro_nano_transport.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * streaming_potential.c - Streaming Potential and Current
 *
 * When a pressure gradient drives flow through a charged channel,
 * the net motion of counter-ions in the EDL generates a convective
 * current. This charge separation produces a potential difference
 * (streaming potential) that opposes further charge transport.
 *
 * Knowledge: L4 (streaming current/potential), L5 (computation)
 * Reference: Hunter (1981) Zeta Potential in Colloid Science
 */

void compute_streaming_potential(double eps, double zeta, double eta,
    double sigma_bulk, double dpdx, double *V_str)
{
    if (V_str == NULL || eta <= 0.0 || sigma_bulk <= 0.0) {
        if (V_str) *V_str = 0.0;
        return;
    }
    *V_str = (eps * zeta / (eta * sigma_bulk)) * fabs(dpdx);
}

double compute_streaming_current(double eps, double zeta, double eta,
    double A, double dpdx, double L, double H, double kappa)
{
    if (eta <= 0.0 || L <= 0.0 || H <= 0.0 || kappa <= 0.0) return -1.0;
    double u_HS = -(eps * zeta / eta);
    double flow_rate = u_HS * A * (1.0 - 2.0 / (kappa * H));
    double Q = flow_rate * fabs(dpdx) * H * H / (12.0 * eta);
    return eps * zeta * A * dpdx / (eta * L) * (1.0 - tanh(kappa * H / 2.0) / (kappa * H / 2.0));
}

double compute_helmholtz_smoluchowski_streaming(double eps, double zeta,
    double eta, double dp, double L, double sigma)
{
    if (eta <= 0.0 || L <= 0.0 || sigma <= 0.0) return -1.0;
    double delta_V = (eps * zeta / (eta * sigma)) * dp;
    return delta_V;
}

double compute_streaming_current_full(double eps, double zeta, double eta,
    double H, double W, double dpdx, double kappa)
{
    if (eta <= 0.0 || H <= 0.0 || W <= 0.0 || kappa <= 0.0) return -1.0;
    double I_str = 0.0;
    double u_HS = -(eps * zeta / eta);
    double half_H = 0.5 * H;
    double cosh_kH = cosh(kappa * half_H);
    for (int i = 0; i < 100; i++) {
        double y = -half_H + (double)i / 99.0 * H;
        double rho_e = -eps * kappa * kappa * zeta * cosh(kappa * y) / cosh_kH;
        double u = u_HS * (1.0 - cosh(kappa * y) / cosh_kH);
        I_str += rho_e * u * W * (H / 100.0);
    }
    return I_str * fabs(dpdx) / (12.0 * eta / (H * H));
}

double compute_surface_conductance(double sigma_s, double mu_counterion,
    double H)
{
    if (H <= 0.0) return -1.0;
    return fabs(sigma_s) * mu_counterion / H;
}

void compute_zeta_from_streaming(double V_str, double dp, double eps,
    double eta, double sigma, double *zeta)
{
    if (zeta == NULL) return;
    if (eta <= 0.0 || sigma <= 0.0 || dp == 0.0) {
        *zeta = 0.0;
        return;
    }
    *zeta = V_str * eta * sigma / (eps * dp);
}

double compute_electroosmotic_pressure(double eps, double zeta,
    double V_applied, double H)
{
    if (H <= 0.0) return -1.0;
    return 12.0 * eps * zeta * V_applied / (H * H);
}

double compute_streaming_potential_efficiency(double V_str, double I_str,
    double dp, double Q)
{
    if (dp <= 0.0 || Q <= 0.0) return -1.0;
    double P_mech = dp * Q;
    double P_elec = V_str * I_str;
    if (P_mech <= 0.0) return -1.0;
    return P_elec / P_mech;
}

double compute_zeta_potential_estimation(double pH, double pI,
    double max_zeta, double slope)
{
    return max_zeta * tanh(slope * (pH - pI));
}

double compute_debye_length_simple(double c_bulk, double z_val,
    double eps_r, double T)
{
    if (c_bulk <= 0.0 || T <= 0.0 || eps_r <= 0.0 || z_val == 0.0)
        return -1.0;
    double eps = eps_r * EPSILON0;
    double kappa_sq = 2.0 * z_val * z_val * E_CHARGE * E_CHARGE
                      * c_bulk * NA / (eps * KB * T);
    return 1.0 / sqrt(kappa_sq);
}

/* Streaming potential for parallel-plate channel with full EDL solution */
void compute_streaming_potential_full(double eps, double zeta, double eta,
    double sigma_bulk, double dpdx, double H, double kappa, double L, double *V_str)
{
    if (V_str == NULL || eta <= 0.0 || sigma_bulk <= 0.0 || kappa <= 0.0
        || H <= 0.0 || L <= 0.0) {
        if (V_str) *V_str = 0.0;
        return;
    }
    double kH = kappa * H;
    double half_kH = 0.5 * kH;
    double shape_factor;
    if (half_kH > 20.0) {
        shape_factor = 1.0 - 2.0 / kH;
    } else {
        shape_factor = 1.0 - tanh(half_kH) / half_kH;
    }
    double I_str = (eps * eps * zeta * zeta / (eta * L)) * dpdx * shape_factor;
    double R_bulk = L / (sigma_bulk * H);
    *V_str = I_str * R_bulk;
}

/* Electro-viscous effect: apparent viscosity increase due to EDL */
double compute_electroviscous_effect(double eps, double zeta, double eta,
    double sigma_bulk, double H, double kappa)
{
    if (eta <= 0.0 || sigma_bulk <= 0.0 || kappa <= 0.0 || H <= 0.0)
        return -1.0;
    double kH = kappa * H;
    double half_kH = 0.5 * kH;
    double G = (eps * zeta / eta) * (eps * zeta / (sigma_bulk * H * H));
    double F;
    if (half_kH > 20.0) {
        F = 3.0 * (1.0 - 1.0 / kH) / (kH * kH);
    } else {
        F = 3.0 / (kH * kH) * (1.0 - tanh(half_kH) / half_kH);
    }
    return eta * (1.0 + G * F);
}

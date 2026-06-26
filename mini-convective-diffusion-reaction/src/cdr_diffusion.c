/**
 * @file cdr_diffusion.c
 * @brief Implementation of diffusion transport laws and models.
 *
 * Fick's first and second laws, effective diffusivity models,
 * analytical solutions for transient and steady-state diffusion.
 *
 * Reference: Bird, Stewart, Lightfoot "Transport Phenomena" (2007), Ch. 17-18
 * Reference: Cussler "Diffusion: Mass Transfer in Fluid Systems" (2009)
 */

#include "cdr_diffusion.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef R_GAS
#define R_GAS 8.314462618  /* J/(mol.K) */
#endif

/* ---------------------------------------------------------------------------
 * L2: Effective Diffusivity Models
 * ------------------------------------------------------------------------- */

double cdr_diffusion_effective(double D_bulk, double porosity, double tortuosity)
{
    /* D_eff = (epsilon / tau) * D_bulk */
    if (porosity <= 0.0 || porosity > 1.0 || tortuosity < 1.0 || D_bulk < 0.0) {
        return 0.0;
    }
    return D_bulk * porosity / tortuosity;
}

double cdr_diffusion_knudsen(double pore_diameter, double temperature,
                             double molecular_weight)
{
    /*
     * D_K = (d_pore / 3) * sqrt(8*R*T / (pi * M))
     *
     * molecular_weight in [kg/mol], R in [J/(mol.K)] = [kg.m^2/(s^2.mol.K)]
     */
    if (pore_diameter <= 0.0 || temperature <= 0.0 || molecular_weight <= 0.0) {
        return 0.0;
    }

    double mean_speed_factor = sqrt(8.0 * R_GAS * temperature /
                                    (M_PI * molecular_weight));
    return (pore_diameter / 3.0) * mean_speed_factor;
}

double cdr_diffusion_bosanquet(double D_bulk, double D_K)
{
    /* 1/D_eff = 1/D_bulk + 1/D_K */
    if (D_bulk <= 0.0 || D_K <= 0.0) {
        return 0.0;
    }
    return 1.0 / (1.0 / D_bulk + 1.0 / D_K);
}

double cdr_diffusion_fuller(double T, double P,
                            double M_A, double M_B,
                            double V_A, double V_B)
{
    /*
     * Fuller-Schettler-Giddings correlation for binary gas diffusion.
     *
     * D_AB = 1.0e-7 * T^1.75 * sqrt(1/M_A + 1/M_B) /
     *        (P_atm * (V_A^(1/3) + V_B^(1/3))^2)   [m^2/s]
     *
     * T in K, P in atm, M in g/mol, V in cm^3/mol.
     */
    if (T <= 0.0 || P <= 0.0 || M_A <= 0.0 || M_B <= 0.0 ||
        V_A <= 0.0 || V_B <= 0.0) {
        return 0.0;
    }

    double T_factor = pow(T, 1.75);
    double M_factor = sqrt(1.0 / M_A + 1.0 / M_B);
    double V_factor = pow(pow(V_A, 1.0/3.0) + pow(V_B, 1.0/3.0), 2.0);

    return 1.0e-7 * T_factor * M_factor / (P * V_factor);
}

double cdr_diffusion_wilke_chang(double T, double mu_B, double M_B,
                                  double V_A, double phi_B)
{
    /*
     * Wilke-Chang correlation for liquid diffusion.
     *
     * D_AB = 7.4e-12 * (phi_B * M_B)^0.5 * T / (mu_B * V_A^0.6)
     *
     * T in K, mu_B in cP, M_B in g/mol, V_A in cm^3/mol.
     * Result in m^2/s.
     */
    if (T <= 0.0 || mu_B <= 0.0 || M_B <= 0.0 || V_A <= 0.0 || phi_B <= 0.0) {
        return 0.0;
    }

    double numerator = 7.4e-12 * sqrt(phi_B * M_B) * T;
    double denominator = mu_B * pow(V_A, 0.6);

    return numerator / denominator;
}

/* ---------------------------------------------------------------------------
 * L4: Fick's Law Implementations
 * ------------------------------------------------------------------------- */

double cdr_fick_first_law_1d(double D, double C1, double C2, double dx)
{
    /*
     * J = -D * (C2 - C1) / dx = D * (C1 - C2) / dx
     *
     * Returns positive flux when C1 > C2 (flow from 1 to 2).
     */
    if (dx <= 0.0) {
        return 0.0;
    }
    return D * (C1 - C2) / dx;
}

void cdr_fick_first_law_3d(double D, double dCdx, double dCdy, double dCdz,
                            DiffusionFlux *flux)
{
    /* J = -D * grad(C) = [-D*dCdx, -D*dCdy, -D*dCdz] */
    if (flux == NULL) {
        return;
    }

    flux->J_diff_x = -D * dCdx;
    flux->J_diff_y = -D * dCdy;
    flux->J_diff_z = -D * dCdz;

    flux->magnitude = sqrt(flux->J_diff_x * flux->J_diff_x +
                           flux->J_diff_y * flux->J_diff_y +
                           flux->J_diff_z * flux->J_diff_z);
    flux->J_diff = flux->magnitude;
}

double cdr_fick_second_law_1d_step(ConcentrationProfile1D *profile,
                                    double D, double dt)
{
    /*
     * dC/dt = D * d^2C/dx^2
     *
     * FTCS (Forward-Time Centered-Space) scheme:
     * C_i^{n+1} = C_i^n + D*dt/dx^2 * (C_{i+1}^n - 2*C_i^n + C_{i-1}^n)
     *
     * Stability: D*dt/dx^2 <= 0.5
     */
    if (profile == NULL || profile->C == NULL || profile->n_points < 3) {
        return 0.0;
    }

    double dx = profile->dx;
    if (dx <= 0.0) {
        return 0.0;
    }

    double r = D * dt / (dx * dx);

    /* Check stability */
    if (r > 0.5) {
        /* Return the maximum stable dt */
        return 0.5 * dx * dx / D;
    }

    /* Allocate temporary array for old values */
    double *C_old = (double *)malloc(profile->n_points * sizeof(double));
    if (C_old == NULL) {
        return 0.0;
    }
    memcpy(C_old, profile->C, profile->n_points * sizeof(double));

    /* FTCS update for interior points */
    for (size_t i = 1; i < profile->n_points - 1; i++) {
        profile->C[i] = C_old[i] + r * (C_old[i+1] - 2.0 * C_old[i] + C_old[i-1]);
    }

    /* Boundary conditions: zero-flux at both ends (Neumann) */
    profile->C[0] = profile->C[1];
    profile->C[profile->n_points - 1] = profile->C[profile->n_points - 2];

    profile->t += dt;
    free(C_old);

    return dt;  /* Step was stable */
}

/* ---------------------------------------------------------------------------
 * L5: Analytical Solutions
 * ------------------------------------------------------------------------- */

/* Complementary error function approximation (Abramowitz & Stegun 7.1.26) */
static double erfc_approx(double x)
{
    /*
     * erf(x) = 2/sqrt(pi) * integral_0^x exp(-t^2) dt
     *
     * Using Horner's method for rational approximation.
     * Maximum error: 1.5e-7
     */
    const double p  =  0.3275911;
    const double a1 =  0.254829592;
    const double a2 = -0.284496736;
    const double a3 =  1.421413741;
    const double a4 = -1.453152027;
    const double a5 =  1.061405429;

    int sign = (x >= 0.0) ? 1 : -1;
    x = fabs(x);

    /* Exact at zero: erfc(0) = 1.0 */
    if (x == 0.0) {
        return (sign < 0) ? 2.0 : 1.0;
    }

    double t = 1.0 / (1.0 + p * x);
    double erfc_val = ((((a5 * t + a4) * t + a3) * t + a2) * t + a1) * t *
                       exp(-x * x);

    if (sign < 0) {
        return 2.0 - erfc_val;
    }
    return erfc_val;
}

double cdr_diffusion_semi_infinite(double C0, double Cs, double D,
                                    double x, double t)
{
    /*
     * C(x,t) = C_s + (C_0 - C_s) * erf(x / (2*sqrt(D*t)))
     *
     * Boundary: C(0,t) = C_s, C(inf,t) = C_0
     */
    if (t <= 0.0 || D <= 0.0) {
        return C0;
    }

    double z = x / (2.0 * sqrt(D * t));
    double erf_z = 1.0 - erfc_approx(z);  /* erf(z) = 1 - erfc(z) */

    return Cs + (C0 - Cs) * erf_z;
}

double cdr_diffusion_penetration_depth(double D, double t)
{
    /*
     * delta_99 = 4 * sqrt(D * t)
     *
     * Where C(delta, t) ~ 0.99 * (C0 - Cs) + Cs
     */
    if (D <= 0.0 || t < 0.0) {
        return 0.0;
    }
    return 4.0 * sqrt(D * t);
}

double cdr_diffusion_timescale(double L, double D)
{
    /* tau_diff = L^2 / D */
    if (D <= 0.0 || L <= 0.0) {
        return 0.0;
    }
    return L * L / D;
}

double cdr_diffusion_membrane_flux(double D, double C_high, double C_low,
                                    double L)
{
    /* J_ss = D * (C_high - C_low) / L (steady-state flux through membrane) */
    if (L <= 0.0) {
        return 0.0;
    }
    return D * (C_high - C_low) / L;
}

double cdr_diffusion_cylindrical_ss(double r, double R_inner, double R_outer,
                                     double C_inner, double C_outer)
{
    /*
     * C(r) = C_inner + (C_outer - C_inner) * ln(r/R_inner) / ln(R_outer/R_inner)
     *
     * Valid for R_inner <= r <= R_outer.
     */
    if (R_inner <= 0.0 || R_outer <= R_inner || r < R_inner || r > R_outer) {
        return C_inner;  /* Return inner value if r is invalid */
    }
    if (r == R_inner) {
        return C_inner;
    }
    if (r == R_outer) {
        return C_outer;
    }

    double log_ratio = log(R_outer / R_inner);
    if (log_ratio == 0.0) {
        return C_inner;
    }

    return C_inner + (C_outer - C_inner) * log(r / R_inner) / log_ratio;
}

double cdr_diffusion_spherical_ss(double r, double R_inner, double R_outer,
                                   double C_inner, double C_outer)
{
    /*
     * C(r) = [C_inner*R_inner*(R_outer - r) + C_outer*R_outer*(r - R_inner)] /
     *        [r * (R_outer - R_inner)]
     *
     * More directly: C(r) = A/r + B where A and B come from BCs.
     *
     * C(r) = C_inner + (C_outer - C_inner) *
     *        (1/R_inner - 1/r) / (1/R_inner - 1/R_outer)
     */
    if (R_inner <= 0.0 || R_outer <= R_inner || r < R_inner || r > R_outer) {
        return C_inner;
    }
    if (r == R_inner) {
        return C_inner;
    }
    if (r == R_outer) {
        return C_outer;
    }

    double inv_r_inner = 1.0 / R_inner;
    double inv_r_outer = 1.0 / R_outer;
    double inv_r = 1.0 / r;

    double denom = inv_r_inner - inv_r_outer;
    if (fabs(denom) < 1e-15) {
        return C_inner;
    }

    double fraction = (inv_r_inner - inv_r) / denom;

    return C_inner + (C_outer - C_inner) * fraction;
}

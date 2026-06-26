#include "micro_nano_transport.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * slip_flow_solver.c - Pressure-Driven Slip Flow
 * Knowledge: L4 (NS+slip BC), L5 (analytical solutions)
 * Reference: Maxwell (1879), Karniadakis et al. (2005)
 */

double compute_maxwell_slip_length(double sigma, double mean_free_path)
{
    if (sigma <= 0.0 || sigma > 1.0 || mean_free_path <= 0.0) return -1.0;
    return ((2.0 - sigma) / sigma) * mean_free_path;
}

double compute_temperature_jump_length(double alpha_T, double mfp,
                                        double gamma, double Pr)
{
    if (alpha_T <= 0.0 || alpha_T > 1.0 || mfp <= 0.0
        || gamma <= 1.0 || Pr <= 0.0) return -1.0;
    return ((2.0 - alpha_T) / alpha_T) * (2.0 * gamma / ((gamma + 1.0) * Pr)) * mfp;
}

void solve_pressure_driven_slip_flow(double H, double dpdx, double eta,
                                      double slip_len, int n,
                                      double *y, double *u)
{
    if (y == NULL || u == NULL || n <= 0 || H <= 0.0 || eta <= 0.0) return;
    double half_H = 0.5 * H;
    double coeff = 1.0 / (2.0 * eta);
    double dpdx_abs = fabs(dpdx);
    for (int i = 0; i < n; i++) {
        double t = (double)i / (double)(n - 1);
        y[i] = -half_H + t * H;
        double y_sq = y[i] * y[i];
        u[i] = coeff * dpdx_abs * (half_H * half_H - y_sq + slip_len * H);
    }
}

double compute_slip_enhanced_flow_rate(double H, double W, double dpdx,
                                        double eta, double slip_len)
{
    if (H <= 0.0 || W <= 0.0 || eta <= 0.0) return -1.0;
    double Q_noslip = (H * H * H * W / (12.0 * eta)) * fabs(dpdx);
    double enhancement = 1.0 + 6.0 * slip_len / H;
    return Q_noslip * enhancement;
}

double compute_poiseuille_number_slip(double H, double slip_len)
{
    if (H <= 0.0) return -1.0;
    return 96.0 / (1.0 + 6.0 * slip_len / H);
}

double compute_enhancement_factor_slip(double H, double slip_len)
{
    if (H <= 0.0) return -1.0;
    return 1.0 + 6.0 * slip_len / H;
}

void solve_second_order_slip_flow(double H, double dpdx, double eta,
                                    double lambda, double A1, double A2,
                                    int n, double *y, double *u)
{
    if (y == NULL || u == NULL || n <= 0 || H <= 0.0 || eta <= 0.0
        || lambda < 0.0) return;
    double half_H = 0.5 * H;
    double coeff = 1.0 / (2.0 * eta);
    double dpdx_abs = fabs(dpdx);
    for (int i = 0; i < n; i++) {
        double t = (double)i / (double)(n - 1);
        y[i] = -half_H + t * H;
        double second_order_term = 2.0 * A1 * lambda * half_H
                                   + 4.0 * A2 * lambda * lambda;
        u[i] = coeff * dpdx_abs * (half_H * half_H - y[i]*y[i] + second_order_term);
    }
}

double compute_wall_shear_stress_slip(double H, double dpdx)
{
    if (H <= 0.0) return -1.0;
    return 0.5 * H * fabs(dpdx);
}

double compute_centerline_to_slip_ratio(double H, double slip_len)
{
    if (H <= 0.0) return -1.0;
    if (slip_len <= 0.0) return 1.0e15;
    return 1.0 + (H * H) / (4.0 * slip_len * H);
}

double compute_slip_correction_mass_flow_gas(double H, double W, double L,
                                               double p_in, double p_out,
                                               double eta, double R_specific,
                                               double T, double lambda_ref,
                                               double p_ref)
{
    if (H <= 0.0 || W <= 0.0 || L <= 0.0 || p_in <= 0.0 || p_out <= 0.0
        || eta <= 0.0 || R_specific <= 0.0 || T <= 0.0 || p_ref <= 0.0)
        return -1.0;
    double p_avg = 0.5 * (p_in + p_out);
    double Kn_avg = (lambda_ref * p_ref) / (p_avg * H);
    double slip_factor = 1.0 + 6.0 * Kn_avg;
    double geo_factor = H * H * H * W / (12.0 * eta * R_specific * T * L);
    double pressure_factor = (p_in * p_in - p_out * p_out) / 2.0;
    return geo_factor * slip_factor * pressure_factor;
}

void compute_velocity_profile_stats(double *y, double *u, int n,
                                     double width, double *u_max,
                                     double *u_avg, double *Q)
{
    if (y == NULL || u == NULL || n < 2 || u_max == NULL
        || u_avg == NULL || Q == NULL) return;
    *u_max = u[0];
    for (int i = 1; i < n; i++) {
        if (u[i] > *u_max) *u_max = u[i];
    }
    double sum = 0.0;
    for (int i = 0; i < n - 1; i++) {
        double dy = y[i+1] - y[i];
        sum += 0.5 * (u[i] + u[i+1]) * dy;
    }
    double total_y = y[n-1] - y[0];
    *u_avg = sum / total_y;
    *Q = sum * width;
}

double compute_maxwell_slip_velocity(double H, double dpdx, double eta,
                                      double slip_len)
{
    if (H <= 0.0 || eta <= 0.0) return -1.0;
    return slip_len * H * fabs(dpdx) / (2.0 * eta);
}

#include "micro_nano_transport.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * electroosmotic_flow.c - Electroosmotic Flow in Micro/Nano Channels
 * Knowledge: L2 (EOF concept), L4 (NS + electric body force), L5 (analytical)
 * Reference: Probstein (1994), Li (2004)
 */

double compute_helmholtz_smoluchowski_velocity(double eps, double zeta,
                                                double eta, double E)
{
    if (eta <= 0.0) return 0.0;
    return -(eps * zeta / eta) * E;
}

double compute_eof_mobility(double eps, double zeta, double eta)
{
    if (eta <= 0.0) return 0.0;
    return -(eps * zeta) / eta;
}

void solve_electroosmotic_flow(double H, double eps, double zeta,
                                double eta, double Ex, double kappa,
                                int n, double *y, double *u)
{
    if (y == NULL || u == NULL || n <= 0 || H <= 0.0
        || eta <= 0.0 || kappa <= 0.0) return;

    double u_HS = -(eps * zeta / eta) * Ex;
    double half_H = 0.5 * H;
    double cosh_kH_half = cosh(kappa * half_H);

    for (int i = 0; i < n; i++) {
        double t = (double)i / (double)(n - 1);
        y[i] = -half_H + t * H;
        if (cosh_kH_half > 1.0e10) {
            u[i] = u_HS;
        } else {
            u[i] = u_HS * (1.0 - cosh(kappa * y[i]) / cosh_kH_half);
        }
    }
}

double compute_eof_flow_rate(double H, double W, double eps, double zeta,
                              double eta, double Ex, double kappa)
{
    if (H <= 0.0 || eta <= 0.0 || kappa <= 0.0) return -1.0;

    double u_HS = -(eps * zeta / eta) * Ex;
    double kH = kappa * H;
    double half_kH = 0.5 * kH;
    double shape_factor;

    if (half_kH > 20.0) {
        shape_factor = 1.0 - 2.0 / kH;
    } else {
        shape_factor = 1.0 - tanh(half_kH) / half_kH;
    }

    return u_HS * H * W * shape_factor;
}

double compute_eof_with_backpressure(double H, double W, double dpdx,
                                      double eps, double zeta,
                                      double eta, double Ex, double kappa)
{
    if (H <= 0.0 || W <= 0.0 || eta <= 0.0 || kappa <= 0.0) return -1.0;

    double Q_eof = compute_eof_flow_rate(H, W, eps, zeta, eta, Ex, kappa);
    double Q_pressure = -(H * H * H * W / (12.0 * eta)) * dpdx;

    return Q_eof + Q_pressure;
}

double compute_max_eof_backpressure(double H, double L, double eps,
                                     double zeta, double eta, double Ex)
{
    if (H <= 0.0 || L <= 0.0 || eta <= 0.0) return -1.0;

    double u_HS = -(eps * zeta / eta) * Ex;
    return (12.0 * eta * L / (H * H)) * u_HS;
}

void compute_eof_profile_full(double H, double eps, double zeta,
                               double eta, double Ex, double kappa,
                               double dpdx, int n,
                               double *y, double *u)
{
    if (y == NULL || u == NULL || n <= 0 || H <= 0.0
        || eta <= 0.0 || kappa <= 0.0) return;

    double u_HS = -(eps * zeta / eta) * Ex;
    double half_H = 0.5 * H;
    double cosh_kH_half = cosh(kappa * half_H);

    for (int i = 0; i < n; i++) {
        double t = (double)i / (double)(n - 1);
        y[i] = -half_H + t * H;
        double psi_y = zeta * cosh(kappa * y[i]) / cosh_kH_half;
        double u_eof = u_HS * (1.0 - psi_y / zeta);
        double u_p = (1.0 / (2.0 * eta)) * fabs(dpdx)
                     * (half_H * half_H - y[i] * y[i]);
        u[i] = u_eof + u_p;
    }
}

double compute_edl_thickness_ratio(double kappa, double H)
{
    if (kappa <= 0.0 || H <= 0.0) return -1.0;
    return 2.0 / (kappa * H);
}

double compute_eof_power_consumption(double sigma_bulk, double A,
                                      double L, double V)
{
    if (L <= 0.0 || A <= 0.0) return -1.0;
    double R = L / (sigma_bulk * A);
    if (R <= 0.0) return -1.0;
    return V * V / R;
}

double compute_eof_efficiency(double Q, double delta_p,
                               double current, double voltage)
{
    if (current <= 0.0 || voltage <= 0.0) return -1.0;
    double P_hyd = Q * delta_p;
    double P_elec = current * voltage;
    if (P_elec <= 0.0) return -1.0;
    return P_hyd / P_elec;
}

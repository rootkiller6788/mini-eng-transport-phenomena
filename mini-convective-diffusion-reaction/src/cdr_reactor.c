/**
 * @file cdr_reactor.c
 * @brief Implementation of ideal and non-ideal reactor models.
 *
 * CSTR, PFR, Batch reactor design equations, RTD analysis,
 * dispersion model, and reactor network optimization.
 *
 * Reference: Levenspiel "Chemical Reaction Engineering" (1999), Ch. 5-15
 * Reference: Fogler "Elements of Chemical Reaction Engineering" (2016)
 */

#include "cdr_reactor.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---------------------------------------------------------------------------
 * Helper: factorial for RTD tanks-in-series
 * ------------------------------------------------------------------------- */

static double factorial_double(int n)
{
    /* Stirling-approximation-based factorial for large n */
    if (n < 0) return 1.0;
    if (n <= 1) return 1.0;

    double result = 1.0;
    for (int i = 2; i <= n; i++) {
        result *= (double)i;
    }
    return result;
}

/* ---------------------------------------------------------------------------
 * L4: CSTR Design Equations
 * ------------------------------------------------------------------------- */

double cdr_cstr_volume_first_order(double C_A0, double Q, double k,
                                    double X_target)
{
    /*
     * CSTR design equation (1st-order, constant density):
     * V = Q * X / (k * (1 - X))
     *
     * Derivation:
     *   C_A = C_A0 * (1 - X)
     *   V = Q * (C_A0 - C_A) / (k * C_A) = Q * X / (k*(1-X))
     */
    (void)C_A0;  /* Not needed for 1st-order (volume independent of C_A0) */

    if (k <= 0.0 || Q <= 0.0 || X_target >= 1.0 || X_target < 0.0) {
        return 0.0;
    }
    if (X_target <= 0.0) {
        return 0.0;
    }

    return Q * X_target / (k * (1.0 - X_target));
}

double cdr_cstr_outlet_first_order(double C_A0, double tau, double k)
{
    /*
     * C_A = C_A0 / (1 + k*tau)
     */
    if (C_A0 < 0.0 || tau < 0.0 || k < 0.0) {
        return 0.0;
    }
    return C_A0 / (1.0 + k * tau);
}

double cdr_cstr_volume_second_order(double C_A0, double Q, double k,
                                     double X_target)
{
    /*
     * CSTR for 2nd-order 2A -> products, constant density:
     *
     * C_A = C_A0 * (1 - X)
     * V = Q * (C_A0 - C_A) / (k * C_A^2)
     *   = Q * X / (k * C_A0 * (1 - X)^2)
     */
    if (C_A0 <= 0.0 || k <= 0.0 || Q <= 0.0 ||
        X_target >= 1.0 || X_target < 0.0) {
        return 0.0;
    }
    if (X_target <= 0.0) {
        return 0.0;
    }

    double one_minus_X = 1.0 - X_target;
    return Q * X_target / (k * C_A0 * one_minus_X * one_minus_X);
}

/* ---------------------------------------------------------------------------
 * L4: PFR Design Equations
 * ------------------------------------------------------------------------- */

double cdr_pfr_volume_first_order(double C_A0, double Q, double k,
                                   double X_target)
{
    /*
     * PFR design (1st-order, constant density):
     * V = Q * ln(1/(1-X)) / k
     *
     * Derivation: dX/dV = k*(1-X)/Q → integral dX/(1-X) = k/Q * integral dV
     */
    (void)C_A0;  /* Not needed for 1st-order (volume independent of C_A0) */

    if (k <= 0.0 || Q <= 0.0 || X_target >= 1.0 || X_target < 0.0) {
        return 0.0;
    }
    if (X_target <= 0.0) {
        return 0.0;
    }

    return Q * (-log(1.0 - X_target)) / k;
}

double cdr_pfr_volume_second_order(double C_A0, double Q, double k,
                                    double X_target)
{
    /*
     * PFR for 2nd-order 2A -> products, constant density:
     *
     * dX/dV = k * C_A0*(1-X)^2 / Q
     * V = Q * X / (k * C_A0 * (1-X))
     */
    if (C_A0 <= 0.0 || k <= 0.0 || Q <= 0.0 ||
        X_target >= 1.0 || X_target < 0.0) {
        return 0.0;
    }
    if (X_target <= 0.0) {
        return 0.0;
    }

    return Q * X_target / (k * C_A0 * (1.0 - X_target));
}

double cdr_pfr_profile_first_order(double C_A0, double k, double tau, double z_L)
{
    /*
     * C_A(z/L) = C_A0 * exp(-k * tau * z/L)
     *
     * tau = V/Q is space time, z/L is fractional position.
     */
    if (C_A0 < 0.0 || k < 0.0 || tau < 0.0 || z_L < 0.0 || z_L > 1.0) {
        return 0.0;
    }
    return C_A0 * exp(-k * tau * z_L);
}

/* ---------------------------------------------------------------------------
 * L4: Batch Reactor Design
 * ------------------------------------------------------------------------- */

double cdr_batch_time_first_order(double k, double X)
{
    /*
     * t_batch = -ln(1 - X) / k   (isothermal, constant volume)
     */
    if (k <= 0.0 || X >= 1.0 || X < 0.0) {
        return INFINITY;
    }
    if (X <= 0.0) {
        return 0.0;
    }
    return -log(1.0 - X) / k;
}

double cdr_batch_time_second_order(double k, double C_A0, double X)
{
    /*
     * 2A -> products: t = X / (k * C_A0 * (1 - X))
     */
    if (k <= 0.0 || C_A0 <= 0.0 || X >= 1.0 || X < 0.0) {
        return INFINITY;
    }
    if (X <= 0.0) {
        return 0.0;
    }
    return X / (k * C_A0 * (1.0 - X));
}

/* ---------------------------------------------------------------------------
 * L5: Reactor Comparison
 * ------------------------------------------------------------------------- */

double cdr_cstr_pfr_volume_ratio(double X)
{
    /*
     * V_CSTR / V_PFR = X/(1-X) / (-ln(1-X))   (1st-order)
     *
     * Ratio >= 1 for positive-order kinetics.
     * Ratio -> 1 as X -> 0.
     * Ratio -> infinity as X -> 1.
     */
    if (X <= 0.0 || X >= 1.0) {
        return 1.0;
    }

    double V_CSTR = X / (1.0 - X);
    double V_PFR  = -log(1.0 - X);

    return V_CSTR / V_PFR;
}

double cdr_pbr_catalyst_mass_first_order(double Q, double k_prime,
                                          double X_target)
{
    /*
     * W = Q * ln(1/(1-X)) / k'
     *
     * where k' has units [m^3/(kg_cat.s)].
     *
     * Analogous to PFR volume but with mass-based rate constant.
     */
    if (k_prime <= 0.0 || Q <= 0.0 || X_target >= 1.0 || X_target < 0.0) {
        return 0.0;
    }
    if (X_target <= 0.0) {
        return 0.0;
    }

    return Q * (-log(1.0 - X_target)) / k_prime;
}

/* ---------------------------------------------------------------------------
 * L5: RTD Analysis
 * ------------------------------------------------------------------------- */

double cdr_rtd_cstr_E(double t, double tau)
{
    /*
     * E(t) = (1/tau) * exp(-t/tau)
     */
    if (tau <= 0.0 || t < 0.0) {
        return 0.0;
    }
    return exp(-t / tau) / tau;
}

double cdr_rtd_tanks_in_series_E(double t, double tau, int N)
{
    /*
     * E(t) = (N^N / ((N-1)! * tau^N)) * t^(N-1) * exp(-N*t/tau)
     *
     * N CSTRs in series, total tau.
     */
    if (tau <= 0.0 || t < 0.0 || N < 1) {
        return 0.0;
    }

    double theta = t / tau;
    double N_double = (double)N;

    double prefactor = pow(N_double, N_double) / factorial_double(N - 1);
    double E = prefactor * pow(theta, N - 1) * exp(-N_double * theta) / tau;

    return E;
}

double cdr_rtd_dispersion_E(double t, double tau, double Pe)
{
    /*
     * E_theta(theta) = sqrt(Pe/(4*pi*theta)) * exp( -Pe*(1-theta)^2 / (4*theta) )
     *
     * where theta = t/tau.
     * This is the open-open vessel dispersion model RTD.
     * Valid for Pe > 10 (approx).
     */
    if (tau <= 0.0 || t <= 0.0 || Pe <= 0.0) {
        return 0.0;
    }

    double theta = t / tau;
    double pre = sqrt(Pe / (4.0 * M_PI * theta));
    double exponent = -Pe * (1.0 - theta) * (1.0 - theta) / (4.0 * theta);

    return pre * exp(exponent) / tau;
}

double cdr_rtd_mean_time(const double *time, const double *E, size_t n_points)
{
    /*
     * t_bar = integral t*E(t) dt / integral E(t) dt
     *
     * Trapezoidal integration.
     */
    if (time == NULL || E == NULL || n_points < 2) {
        return 0.0;
    }

    double sum_E  = 0.0;
    double sum_tE = 0.0;

    for (size_t i = 1; i < n_points; i++) {
        double dt = time[i] - time[i-1];
        double avg_E = (E[i] + E[i-1]) / 2.0;
        double avg_tE = (time[i]*E[i] + time[i-1]*E[i-1]) / 2.0;
        sum_E  += avg_E * dt;
        sum_tE += avg_tE * dt;
    }

    if (sum_E <= 0.0) {
        return 0.0;
    }

    return sum_tE / sum_E;
}

double cdr_rtd_variance(const double *time, const double *E,
                         size_t n_points, double t_bar)
{
    /*
     * sigma^2 = integral (t - t_bar)^2 * E(t) dt
     */
    if (time == NULL || E == NULL || n_points < 2) {
        return 0.0;
    }

    double sum_variance = 0.0;
    double sum_E = 0.0;

    for (size_t i = 1; i < n_points; i++) {
        double dt = time[i] - time[i-1];
        double diff1 = time[i-1] - t_bar;
        double diff2 = time[i] - t_bar;
        double f1 = E[i-1] * diff1 * diff1;
        double f2 = E[i] * diff2 * diff2;
        double avg = (f1 + f2) / 2.0;
        double avg_E_segment = (E[i-1] + E[i]) / 2.0;

        sum_variance += avg * dt;
        sum_E += avg_E_segment * dt;
    }

    if (sum_E <= 0.0) {
        return 0.0;
    }

    return sum_variance;  /* Already normalized if E is normalized */
}

double cdr_rtd_peclet_from_variance(double variance, double tau)
{
    /*
     * For open-open dispersion model:
     * sigma_theta^2 = sigma^2/tau^2 = 2/Pe + 8/Pe^2
     *
     * Solve for Pe: Pe = (-2 + sqrt(4 + 32*sigma_theta^2)) / (2*sigma_theta^2)
     *
     * Simplified: Pe ? 2 / sigma_theta^2 (for Pe >> 1)
     */
    if (tau <= 0.0 || variance < 0.0) {
        return 0.0;
    }

    double sigma_theta_sq = variance / (tau * tau);

    if (sigma_theta_sq < 1e-10) {
        return 1e10;  /* Essentially plug flow */
    }

    /* Solve quadratic: sigma^2 * Pe^2 - 2*Pe - 8 = 0 */
    double a = sigma_theta_sq;
    double b = -2.0;
    double c = -8.0;

    double discriminant = b * b - 4.0 * a * c;
    if (discriminant < 0.0) {
        return 0.0;
    }

    double Pe = (-b + sqrt(discriminant)) / (2.0 * a);
    return Pe;
}

double cdr_rtd_N_tanks_from_variance(double tau, double variance)
{
    /*
     * N = tau^2 / sigma^2   (tanks-in-series model)
     */
    if (variance <= 0.0 || tau <= 0.0) {
        return 1e10;  /* Plug flow */
    }

    /* Approximate: N >= 1, integer */
    double N = tau * tau / variance;
    return (N < 1.0) ? 1.0 : N;
}

double cdr_segregation_model_conversion(const double *time, const double *E,
                                         size_t n_points, double k)
{
    /*
     * X_bar = integral_0^inf (1 - exp(-k*t)) * E(t) dt
     *
     * Uses trapezoidal integration for first-order kinetics.
     */
    if (time == NULL || E == NULL || n_points < 2 || k <= 0.0) {
        return 0.0;
    }

    double sum = 0.0;

    for (size_t i = 1; i < n_points; i++) {
        double dt = time[i] - time[i-1];
        double x_batch1 = 1.0 - exp(-k * time[i-1]);
        double x_batch2 = 1.0 - exp(-k * time[i]);
        double f1 = E[i-1] * x_batch1;
        double f2 = E[i] * x_batch2;
        sum += (f1 + f2) / 2.0 * dt;
    }

    return sum;
}

/* ---------------------------------------------------------------------------
 * L5: Dispersion Model
 * ------------------------------------------------------------------------- */

double cdr_dispersion_outlet_concentration(double C_A0, double Pe, double Da)
{
    /*
     * Analytical solution for 1st-order reaction in dispersion model
     * with Danckwerts boundary conditions.
     *
     * C_out/C_A0 = 4*q*exp(Pe/2) / ((1+q)^2*exp(q*Pe/2) - (1-q)^2*exp(-q*Pe/2))
     *
     * where q = sqrt(1 + 4*Da/Pe)
     */
    if (C_A0 < 0.0 || Pe <= 0.0 || Da < 0.0) {
        return 0.0;
    }
    if (Da <= 0.0) {
        return C_A0;  /* No reaction */
    }

    double q = sqrt(1.0 + 4.0 * Da / Pe);
    double exp_Pe_half = exp(Pe / 2.0);
    double exp_qPe_half = exp(q * Pe / 2.0);
    double exp_minus_qPe_half = exp(-q * Pe / 2.0);

    double numerator = 4.0 * q * exp_Pe_half;
    double one_plus_q_sq = (1.0 + q) * (1.0 + q);
    double one_minus_q_sq = (1.0 - q) * (1.0 - q);

    double denominator = one_plus_q_sq * exp_qPe_half -
                         one_minus_q_sq * exp_minus_qPe_half;

    if (fabs(denominator) < 1e-15) {
        return 0.0;
    }

    return C_A0 * numerator / denominator;
}

double cdr_dispersion_conversion(double Pe, double Da)
{
    /*
     * X = 1 - C_out/C_A0
     */
    double C_ratio = cdr_dispersion_outlet_concentration(1.0, Pe, Da);
    return 1.0 - C_ratio;
}

/* ---------------------------------------------------------------------------
 * L5: Reactor Network Optimization
 * ------------------------------------------------------------------------- */

double cdr_cstr_series_conversion(double k, double total_tau, int N)
{
    /*
     * N equal-volume CSTRs in series:
     * X_N = 1 - 1/(1 + k*tau_i)^N
     *
     * where tau_i = total_tau / N
     */
    if (k <= 0.0 || total_tau <= 0.0 || N < 1) {
        return 0.0;
    }

    double tau_i = total_tau / (double)N;
    double X_N = 1.0 - 1.0 / pow(1.0 + k * tau_i, (double)N);

    return X_N;
}

double cdr_pfr_recycle_volume_factor(double X_target, double k,
                                      double tau, double R)
{
    /*
     * PFR with recycle.
     *
     * For first-order reaction:
     * V/Q0 = (R+1)/k * ln( (1 + R*(1-X)) / (R*(1-X)) )  [when R*(1-X) > 0]
     *
     * Without recycle: V_PFR/Q0 = ln(1/(1-X))/k
     *
     * This function returns the volume factor (V_recycle / V_no_recycle).
     * The ratio is independent of tau (both scale with Q), so tau is unused.
     *
     * R = recycle ratio = Q_recycle / Q_fresh
     */
    (void)tau;  /* Volume ratio is independent of space time */

    if (k <= 0.0 || X_target >= 1.0 || X_target <= 0.0 || R < 0.0) {
        return 1.0;
    }
    if (R == 0.0) {
        return 1.0;  /* No recycle, same as PFR */
    }

    /* Denominator for recycle case: R*(1-X) must be > 0 */
    double denom = R * (1.0 - X_target);
    if (denom <= 0.0) {
        return 1e10;  /* Impractical, requires infinite volume */
    }

    double V_recycle_factor = (R + 1.0) * log((1.0 + denom) / denom) / k;
    double V_no_recycle_factor = (-log(1.0 - X_target)) / k;

    if (V_no_recycle_factor <= 0.0) {
        return 1e10;
    }

    return V_recycle_factor / V_no_recycle_factor;
}

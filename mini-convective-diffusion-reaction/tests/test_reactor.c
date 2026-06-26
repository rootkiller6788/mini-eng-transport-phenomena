/**
 * @file test_reactor.c
 * @brief Tests for reactor models, RTD, and dispersion.
 */

#include "../include/cdr_reactor.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TOL 1e-8

static void test_cstr_design(void)
{
    /* CSTR 1st-order */
    double V = cdr_cstr_volume_first_order(1.0, 0.001, 0.1, 0.9);
    assert(V > 0.0);

    /* CSTR outlet concentration */
    double C_out = cdr_cstr_outlet_first_order(1.0, 100.0, 0.01);
    assert(fabs(C_out - 0.5) < TOL);

    /* CSTR 2nd-order */
    double V2 = cdr_cstr_volume_second_order(1.0, 0.001, 0.1, 0.9);
    assert(V2 > 0.0);

    /* Edge cases */
    assert(cdr_cstr_volume_first_order(1.0, 1.0, 0.0, 0.9) == 0.0);
    assert(cdr_cstr_volume_first_order(1.0, 1.0, 1.0, 1.0) == 0.0);

    printf("  [PASS] test_cstr_design\n");
}

static void test_pfr_design(void)
{
    /* PFR 1st-order: V = Q*ln(1/(1-X))/k */
    double V = cdr_pfr_volume_first_order(1.0, 0.001, 0.1, 0.9);
    double expected = 0.001 * log(10.0) / 0.1;
    assert(fabs(V - expected) < TOL);

    /* PFR 2nd-order */
    double V2 = cdr_pfr_volume_second_order(1.0, 0.001, 0.1, 0.9);
    assert(V2 > 0.0);

    /* PFR should be smaller than CSTR for positive order */
    double V_cstr = cdr_cstr_volume_first_order(1.0, 0.001, 0.1, 0.9);
    assert(V < V_cstr);

    /* PFR concentration profile */
    double C_mid = cdr_pfr_profile_first_order(1.0, 0.1, 100.0, 0.5);
    double C_out = cdr_pfr_profile_first_order(1.0, 0.1, 100.0, 1.0);
    assert(C_mid > C_out);  /* Concentration decreases along reactor */

    printf("  [PASS] test_pfr_design\n");
}

static void test_batch_reactor(void)
{
    /* Batch 1st-order */
    double t = cdr_batch_time_first_order(0.1, 0.9);
    double expected = -log(0.1) / 0.1;
    assert(fabs(t - expected) < TOL);

    /* Batch 2nd-order */
    double t2 = cdr_batch_time_second_order(0.1, 1.0, 0.9);
    assert(t2 > 0.0);

    assert(isinf(cdr_batch_time_first_order(0.0, 0.9)));
    printf("  [PASS] test_batch_reactor\n");
}

static void test_cstr_pfr_comparison(void)
{
    /* V_CSTR / V_PFR >= 1 for positive-order reactions */
    double ratio = cdr_cstr_pfr_volume_ratio(0.5);
    assert(ratio > 1.0);

    double ratio2 = cdr_cstr_pfr_volume_ratio(0.9);
    assert(ratio2 > ratio);  /* Ratio increases with conversion */

    /* As X -> 0, ratio -> 1 */
    double ratio3 = cdr_cstr_pfr_volume_ratio(0.001);
    assert(fabs(ratio3 - 1.0) < 0.01);

    printf("  [PASS] test_cstr_pfr_comparison\n");
}

static void test_pbr(void)
{
    double W = cdr_pbr_catalyst_mass_first_order(0.001, 0.001, 0.9);
    assert(W > 0.0);

    assert(cdr_pbr_catalyst_mass_first_order(1.0, 0.0, 0.9) == 0.0);
    printf("  [PASS] test_pbr\n");
}

static void test_rtd_cstr(void)
{
    /* E(t) must integrate to 1 (approximate) */
    double tau = 10.0;
    double dt = 0.1;
    double integral = 0.0;
    for (double t = 0.0; t < 100.0; t += dt) {
        integral += cdr_rtd_cstr_E(t + dt/2.0, tau) * dt;
    }
    assert(fabs(integral - 1.0) < 0.01);

    /* E(0) = 1/tau */
    assert(fabs(cdr_rtd_cstr_E(0.0, tau) - 0.1) < TOL);

    printf("  [PASS] test_rtd_cstr\n");
}

static void test_rtd_tanks_in_series(void)
{
    /* N tanks in series, should approach PFR as N -> large */
    double tau = 10.0;
    double E_peak_N1 = cdr_rtd_tanks_in_series_E(9.0, tau, 1);
    double E_peak_N10 = cdr_rtd_tanks_in_series_E(9.0, tau, 10);

    /* More tanks ? sharper distribution near tau */
    assert(E_peak_N10 > E_peak_N1);

    /* At t=tau (=10), distribution peaks for N>1 */
    /* For N=1 (CSTR), peak is at t=0 */
    assert(cdr_rtd_tanks_in_series_E(0.0, tau, 1) >
           cdr_rtd_tanks_in_series_E(tau, tau, 1));

    printf("  [PASS] test_rtd_tanks_in_series\n");
}

static void test_rtd_dispersion(void)
{
    double tau = 10.0;

    /* Large Pe ? near plug flow */
    double E_hi = cdr_rtd_dispersion_E(10.0, tau, 100.0);
    assert(E_hi > 0.0);

    /* Low Pe ? closer to CSTR behavior */
    double E_lo = cdr_rtd_dispersion_E(10.0, tau, 1.0);
    /* Both should be positive */
    assert(E_lo > 0.0);

    printf("  [PASS] test_rtd_dispersion\n");
}

static void test_rtd_statistics(void)
{
    /* Generate CSTR RTD data and compute statistics */
    double tau = 10.0;
    size_t n = 1000;
    double *time = (double *)malloc(n * sizeof(double));
    double *E = (double *)malloc(n * sizeof(double));

    for (size_t i = 0; i < n; i++) {
        time[i] = i * 0.1;
        E[i] = cdr_rtd_cstr_E(time[i], tau);
    }

    double t_bar = cdr_rtd_mean_time(time, E, n);
    assert(fabs(t_bar - tau) / tau < 0.05);

    double var = cdr_rtd_variance(time, E, n, t_bar);
    double expected_var = tau * tau;
    assert(fabs(var - expected_var) / expected_var < 0.1);

    /* Pe from variance (should be ~0 for pure CSTR ? small Pe) */
    double Pe = cdr_rtd_peclet_from_variance(var, tau);
    assert(Pe > 0.0);

    /* N tanks from variance (pure CSTR ? N=1) */
    double N = cdr_rtd_N_tanks_from_variance(tau, var);
    assert(fabs(N - 1.0) < 0.2);

    free(time);
    free(E);
    printf("  [PASS] test_rtd_statistics\n");
}

static void test_segregation_model(void)
{
    /* For CSTR RTD with first-order reaction */
    double tau = 10.0;
    double k = 0.1;
    size_t n = 1000;
    double *time = (double *)malloc(n * sizeof(double));
    double *E = (double *)malloc(n * sizeof(double));

    for (size_t i = 0; i < n; i++) {
        time[i] = i * 0.1;
        E[i] = cdr_rtd_cstr_E(time[i], tau);
    }

    double X_seg = cdr_segregation_model_conversion(time, E, n, k);

    /* CSTR analytical: X = k*tau/(1+k*tau) = 1.0/2.0 = 0.5 */
    double X_analytical = k * tau / (1.0 + k * tau);
    assert(fabs(X_seg - X_analytical) < 0.05);

    free(time);
    free(E);
    printf("  [PASS] test_segregation_model\n");
}

static void test_dispersion_model(void)
{
    /* Large Pe, moderate Da: should approach PFR */
    double Cout1 = cdr_dispersion_outlet_concentration(1.0, 100.0, 2.302585);
    /* For large Pe: X ~ 1-exp(-Da) = 1-exp(-2.3026) = 0.9 */
    double expected_cout = exp(-2.302585);
    assert(fabs(Cout1 - expected_cout) < 0.05);

    /* Conversion */
    double X = cdr_dispersion_conversion(100.0, 2.302585);
    assert(fabs(X - 0.9) < 0.05);

    printf("  [PASS] test_dispersion_model\n");
}

static void test_cstr_series(void)
{
    /* N CSTRs in series: as N increases, approach PFR performance */
    double X1 = cdr_cstr_series_conversion(0.1, 10.0, 1);
    double X5 = cdr_cstr_series_conversion(0.1, 10.0, 5);
    double X10 = cdr_cstr_series_conversion(0.1, 10.0, 10);

    assert(X5 > X1);
    assert(X10 > X5);

    /* Limit: as N -> infinity, X -> 1-exp(-k*tau) = 1-exp(-1) = 0.6321 */
    double X_PFR_limit = 1.0 - exp(-1.0);
    assert(X10 < X_PFR_limit + 0.01);

    printf("  [PASS] test_cstr_series\n");
}

static void test_pfr_recycle(void)
{
    /* With recycle, the reactor behaves more like CSTR */
    double factor = cdr_pfr_recycle_volume_factor(0.9, 0.1, 23.02585, 1.0);
    assert(factor > 1.0);  /* Recycle requires more volume */

    /* R=0 is same as no recycle */
    double factor_R0 = cdr_pfr_recycle_volume_factor(0.9, 0.1, 23.02585, 0.0);
    assert(fabs(factor_R0 - 1.0) < TOL);

    printf("  [PASS] test_pfr_recycle\n");
}

int main(void)
{
    printf("CDR Reactor Tests:\n");
    test_cstr_design();
    test_pfr_design();
    test_batch_reactor();
    test_cstr_pfr_comparison();
    test_pbr();
    test_rtd_cstr();
    test_rtd_tanks_in_series();
    test_rtd_dispersion();
    test_rtd_statistics();
    test_segregation_model();
    test_dispersion_model();
    test_cstr_series();
    test_pfr_recycle();
    printf("All reactor tests passed.\n");
    return 0;
}

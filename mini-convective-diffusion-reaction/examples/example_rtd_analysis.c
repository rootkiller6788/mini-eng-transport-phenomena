/**
 * @file example_rtd_analysis.c
 * @brief Example: Residence Time Distribution analysis for reactor diagnosis.
 *
 * Problem: A chemical reactor is suspected of having dead zones and bypassing.
 * A tracer pulse injection experiment is simulated, and RTD analysis is used
 * to diagnose the non-ideal behavior.
 *
 * This example computes E(t) for CSTR, N-tanks-in-series, and dispersion
 * models, then shows how to estimate model parameters from RTD data.
 *
 * Application: Industrial reactor troubleshooting, FDA process validation,
 * wastewater treatment (Detroit WWTP uses RTD for chlorine contactor design).
 *
 * Reference: Levenspiel (1999), Ch. 11-13
 * Reference: Fogler (2016), Ch. 16 ? RTD and Diagnosing Reactor Ills
 */

#include "../include/cdr_reactor.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main(void)
{
    printf("??????????????????????????????????????????????????\n");
    printf("  RTD Analysis ? Reactor Diagnosis\n");
    printf("??????????????????????????????????????????????????\n\n");

    double tau = 60.0;  /* Design mean residence time [s] */
    double Q   = 0.01;  /* Volumetric flow rate [m?/s] */
    double V   = tau * Q; /* Reactor volume = 0.6 m? */

    printf("Reactor specifications:\n");
    printf("  Design ? = %.0f s, V = %.3f m?, Q = %.3f m?/s\n\n", tau, V, Q);

    /* Generate RTD data for N=3 tanks-in-series (moderate backmixing) */
    int N_model = 3;
    size_t n_pts = 200;
    double *time = (double *)malloc(n_pts * sizeof(double));
    double *E_data = (double *)malloc(n_pts * sizeof(double));

    printf("Simulated RTD data (N=%d tanks-in-series model):\n", N_model);
    printf("  Time [s]   E(t) [1/s]    F(t)\n");
    printf("  ???????????????????????????????\n");

    double F_cum = 0.0;
    double dt = tau * 5.0 / (double)n_pts;

    for (size_t i = 0; i < n_pts; i++) {
        time[i] = i * dt;
        E_data[i] = cdr_rtd_tanks_in_series_E(time[i], tau, N_model);

        /* Trapezoidal integration for F(t) */
        if (i > 0) {
            F_cum += (E_data[i] + E_data[i-1]) / 2.0 * dt;
        }

        /* Print selected points */
        if (i % 20 == 0 || i == n_pts - 1) {
            printf("  %8.1f   %.6f      %.4f\n", time[i], E_data[i], F_cum);
        }
    }

    /* Compute RTD statistics from the simulated data */
    double t_bar = cdr_rtd_mean_time(time, E_data, n_pts);
    double variance = cdr_rtd_variance(time, E_data, n_pts, t_bar);

    printf("\nRTD Statistics from simulated data:\n");
    printf("  Mean residence time t_bar = %.2f s (design: %.0f s)\n", t_bar, tau);
    printf("  Variance sigma^2 = %.2f s^2\n", variance);
    printf("  sigma^2/tau^2 = %.4f\n\n", variance / (tau * tau));

    /* Diagnostic: Does the reactor have dead volume? */
    double V_active = Q * t_bar;
    double V_dead = V - V_active;
    printf("Reactor Diagnosis:\n");
    printf("  Active volume  = %.4f m? (%.1f%% of design)\n",
           V_active, V_active / V * 100);
    if (V_dead > V * 0.05) {
        printf("  Dead volume    = %.4f m? (%.1f%%)  ?? Possible dead zones!\n",
               V_dead, V_dead / V * 100);
    } else {
        printf("  Dead volume    = %.4f m? (negligible)\n", V_dead);
    }

    /* Estimate model parameters from RTD statistics */
    double N_est = cdr_rtd_N_tanks_from_variance(t_bar, variance);
    double Pe_est = cdr_rtd_peclet_from_variance(variance, t_bar);

    printf("\nModel parameter estimation:\n");
    printf("  Tanks-in-series: N ? %.1f (true: %d)\n", N_est, N_model);
    printf("  Dispersion:      Pe ? %.1f (D/uL = %.4f)\n", Pe_est, 1.0/Pe_est);

    /* Conversion prediction using segregation model */
    double k = 0.02;  /* 1/s */
    double X_cstr_ideal = k * tau / (1.0 + k * tau);
    double X_pfr_ideal = 1.0 - exp(-k * tau);
    double X_seg = cdr_segregation_model_conversion(time, E_data, n_pts, k);
    double X_tanks = 1.0 - 1.0 / pow(1.0 + k * tau / N_model, N_model);

    printf("\nConversion prediction (k=%.3f 1/s):\n", k);
    printf("  Ideal CSTR:         X = %.4f (%.1f%%)\n",
           X_cstr_ideal, X_cstr_ideal * 100);
    printf("  Ideal PFR:          X = %.4f (%.1f%%)\n",
           X_pfr_ideal, X_pfr_ideal * 100);
    printf("  N=%d tanks-in-series: X = %.4f (%.1f%%)\n",
           N_model, X_tanks, X_tanks * 100);
    printf("  Segregation model:  X = %.4f (%.1f%%)\n",
           X_seg, X_seg * 100);

    printf("\n  Note: For 1st-order kinetics, segregation = max. mixedness.\n");
    printf("  The segregation model result matches the tanks-in-series model\n");
    printf("  because both represent the same micromixing state for linear kinetics.\n");

    free(time);
    free(E_data);

    printf("\n??????????????????????????????????????????????????\n");
    printf("  RTD analysis complete.\n");
    printf("??????????????????????????????????????????????????\n");

    return 0;
}

/**
 * @file demo_conversion.c
 * @brief Interactive demo: conversion vs. space time for different reactor types.
 *
 * Demonstrates CSTR, PFR, and dispersion model behavior.
 */

#include "../include/cdr_reactor.h"
#include <stdio.h>
#include <math.h>

int main(void)
{
    printf("╔════════════════════════════════════════════╗\n");
    printf("║  Reactor Conversion Demo                  ║\n");
    printf("║  A → B  (1st-order, k=0.1 s⁻¹)           ║\n");
    printf("╚════════════════════════════════════════════╝\n\n");

    double k = 0.1;  /* 1/s */

    printf("Conversion vs. Space Time:\n");
    printf("  τ [s]    X_CSTR    X_PFR     X_disp(Pe=5)  X_disp(Pe=20)\n");
    printf("  ─────────────────────────────────────────────────────────\n");

    double tau_values[] = {1, 2, 5, 10, 20, 30, 50, 100};
    int n_tau = 8;

    for (int i = 0; i < n_tau; i++) {
        double tau = tau_values[i];
        double Da = k * tau;

        double X_cstr = Da / (1.0 + Da);
        double X_pfr  = 1.0 - exp(-Da);
        double X_d5   = cdr_dispersion_conversion(5.0, Da);
        double X_d20  = cdr_dispersion_conversion(20.0, Da);

        printf("  %6.0f   %7.4f   %7.4f   %8.4f      %8.4f\n",
               tau, X_cstr, X_pfr, X_d5, X_d20);
    }

    printf("\nObservations:\n");
    printf("  • PFR always achieves higher conversion than CSTR.\n");
    printf("  • Dispersion model (Pe=5) lies between PFR and CSTR.\n");
    printf("  • As Pe → ∞ (plug flow), dispersion → PFR.\n");
    printf("  • As Pe → 0 (complete mixing), dispersion → CSTR.\n");

    printf("\nVolume ratio V_CSTR/V_PFR for different conversions:\n");
    printf("  X        V_CSTR/V_PFR\n");
    printf("  ─────────────────────\n");

    double X_values[] = {0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 0.99};
    int n_X = 7;
    for (int i = 0; i < n_X; i++) {
        double X = X_values[i];
        double ratio = cdr_cstr_pfr_volume_ratio(X);
        printf("  %.2f     %.3f\n", X, ratio);
    }

    printf("\nKey insight: The CSTR penalty grows dramatically at high conversion.\n");
    printf("At 99%% conversion, the CSTR requires ~21.5× the PFR volume.\n");

    return 0;
}

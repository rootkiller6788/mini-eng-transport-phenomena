/**
 * @file example_dispersion.c
 * @brief Example: Axial dispersion model for a tubular reactor.
 *
 * Problem: A first-order liquid-phase reaction is carried out in a tubular
 * reactor with significant axial dispersion. Determine the conversion
 * as a function of the dispersion number (D/(uL) = 1/Pe) and compare
 * with ideal PFR and CSTR predictions.
 *
 * Application: This models real tubular reactors where backmixing occurs
 * due to turbulence, packing, or low length-to-diameter ratios.
 * Examples include industrial polymerization reactors and wastewater
 * treatment UV disinfection channels (EPA design guidelines).
 *
 * Reference: Levenspiel (1999), Ch. 13 ? The Dispersion Model
 * Reference: EPA "UV Disinfection Guidance Manual" (2006) ? dispersion in open channels
 */

#include "../include/cdr_reactor.h"
#include <stdio.h>
#include <math.h>

int main(void)
{
    printf("???????????????????????????????????????????????????\n");
    printf("  Axial Dispersion Model ? Tubular Reactor\n");
    printf("  A ? B (1st-order, liquid phase)\n");
    printf("???????????????????????????????????????????????????\n\n");

    double k    = 0.05;    /* 1/s */
    double tau  = 100.0;   /* s (space time) */
    double Da   = k * tau; /* Damkohler-I = 5.0 */

    /* Ideal limits */
    double X_cstr = Da / (1.0 + Da);
    double X_pfr  = 1.0 - exp(-Da);

    printf("Parameters:\n");
    printf("  k = %.3f s??, ? = %.0f s, Da_I = %.2f\n", k, tau, Da);
    printf("  CSTR limit: X = %.4f (%.1f%%)\n", X_cstr, X_cstr * 100);
    printf("  PFR  limit: X = %.4f (%.1f%%)\n\n", X_pfr, X_pfr * 100);

    /* Dispersion model for various Peclet numbers */
    printf("Dispersion model results:\n");
    printf("  Pe       D/(uL)    X         C_out/C_A0\n");
    printf("  ????????????????????????????????????????\n");

    double Pe_values[] = {0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0, 100.0, 1000.0};
    int n_Pe = 9;

    for (int i = 0; i < n_Pe; i++) {
        double Pe = Pe_values[i];
        double X_disp = cdr_dispersion_conversion(Pe, Da);
        double C_out_ratio = 1.0 - X_disp;

        printf("  %7.1f  %.4f    %.4f    %.4f",
               Pe, 1.0/Pe, X_disp, C_out_ratio);

        if (Pe < 1.0) {
            printf("   ? near CSTR");
        } else if (Pe > 50.0) {
            printf("   ? near PFR");
        }
        printf("\n");
    }

    /* RTD analysis for the dispersion model */
    printf("\nRTD E(?) for different dispersion numbers:\n");
    printf("  ?=t/?    Pe=1 (D/uL=1)   Pe=10          Pe=100\n");
    printf("  ????????????????????????????????????????????????\n");

    double theta_values[] = {0.1, 0.3, 0.5, 0.7, 0.9, 1.0, 1.2, 1.5, 2.0, 3.0};
    int n_theta = 10;

    for (int i = 0; i < n_theta; i++) {
        double theta = theta_values[i];
        double t = theta * tau;

        double E1   = cdr_rtd_dispersion_E(t, tau, 1.0);
        double E10  = cdr_rtd_dispersion_E(t, tau, 10.0);
        double E100 = cdr_rtd_dispersion_E(t, tau, 100.0);

        printf("  %.1f     %.4f          %.4f          %.4f\n",
               theta, E1 * tau, E10 * tau, E100 * tau);
    }
    /* Note: E_?(?) = ? * E(t), so we multiply by ? to get dimensionless */

    /* RTD statistics: estimate Pe from variance */
    printf("\nRTD statistics demonstration:\n");
    double tau_test = 10.0;
    printf("  For a CSTR (?=%.0f s):\n", tau_test);
    printf("    Ideal variance ?? = ?? = %.0f s?\n", tau_test * tau_test);

    /* Compare tanks-in-series with dispersion */
    printf("\nTanks-in-series vs. Dispersion equivalence:\n");
    printf("  N tanks    Equivalent Pe ? 2*(N-1)\n");
    printf("  ?????????????????????????????????\n");
    for (int N = 1; N <= 10; N++) {
        double N_d = (double)N;
        double Pe_eq = 2.0 * (N_d - 1.0);  /* Approximate equivalence for N>=2 */
        if (N == 1) Pe_eq = 0.0;
        printf("  %2d         %.0f\n", N, Pe_eq);
    }

    /* Segregation model ? late mixing vs. early mixing */
    printf("\nSegregation Model (late mixing) for a real reactor:\n");
    printf("  (Mixedness does NOT matter for 1st-order kinetics)\n");
    printf("  For linear kinetics, segregation = maximum mixedness.\n");

    printf("\n???????????????????????????????????????????????????\n");
    printf("  Analysis complete.\n");
    printf("???????????????????????????????????????????????????\n");

    return 0;
}

/**
 * @file example_cstr_design.c
 * @brief Example: CSTR and PFR sizing for a first-order liquid-phase reaction.
 *
 * Problem: A first-order irreversible reaction A ? B (k=0.5 min??) is carried
 * out in liquid phase. The feed is C_A0=2.0 mol/L at Q=10 L/min. Compare the
 * reactor volumes needed for 90% conversion in a CSTR vs. a PFR.
 *
 * Reference: Levenspiel (1999), Example 5.1-5.3
 */

#include "../include/cdr_reactor.h"
#include <stdio.h>
#include <math.h>

int main(void)
{
    double C_A0 = 2.0;       /* mol/L */
    double Q   = 10.0;       /* L/min */
    double k   = 0.5;        /* 1/min */
    double X_target = 0.90;  /* 90% conversion */

    printf("??????????????????????????????????????????\n");
    printf("  Reactor Sizing: A ? B (1st order)\n");
    printf("  C_A0 = %.1f mol/L, Q = %.1f L/min\n", C_A0, Q);
    printf("  k = %.2f min??, Target X = %.0f%%\n", k, X_target * 100);
    printf("??????????????????????????????????????????\n\n");

    /* CSTR design */
    double V_cstr = cdr_cstr_volume_first_order(C_A0, Q, k, X_target);
    double tau_cstr = V_cstr / Q;
    double C_A_cstr = cdr_cstr_outlet_first_order(C_A0, tau_cstr, k);
    double X_cstr_check = 1.0 - C_A_cstr / C_A0;

    printf("CSTR Design:\n");
    printf("  Space time  ? = V/Q = %.2f min\n", tau_cstr);
    printf("  Volume      V = %.2f L\n", V_cstr);
    printf("  Outlet C_A  = %.4f mol/L\n", C_A_cstr);
    printf("  Conversion  X = %.4f (check: %.4f)\n\n", X_cstr_check,
           k * tau_cstr / (1.0 + k * tau_cstr));

    /* PFR design */
    double V_pfr = cdr_pfr_volume_first_order(C_A0, Q, k, X_target);
    double tau_pfr = V_pfr / Q;

    printf("PFR Design:\n");
    printf("  Space time  ? = V/Q = %.2f min\n", tau_pfr);
    printf("  Volume      V = %.2f L\n", V_pfr);
    printf("  Outlet C_A  = %.4f mol/L\n",
           C_A0 * exp(-k * tau_pfr));
    printf("  Conversion  X = %.4f\n\n",
           1.0 - exp(-k * tau_pfr));

    /* Comparison */
    double ratio = cdr_cstr_pfr_volume_ratio(X_target);
    printf("??????????????????????????????????????????\n");
    printf("  V_CSTR / V_PFR = %.2f\n", ratio);
    printf("  The CSTR requires %.1f? more volume\n", ratio);
    printf("  than the PFR for 90%% conversion.\n");
    printf("??????????????????????????????????????????\n\n");

    /* Show concentration profile along PFR */
    printf("PFR axial concentration profile:\n");
    printf("  z/L    C_A(z/L) [mol/L]\n");
    printf("  ---------------------------\n");
    for (int i = 0; i <= 10; i++) {
        double z_L = i / 10.0;
        double C = cdr_pfr_profile_first_order(C_A0, k, tau_pfr, z_L);
        printf("  %.1f    %.4f\n", z_L, C);
    }

    /* Demonstrate CSTR-in-series improvement */
    printf("\nCSTRs in series (equal total ?):\n");
    for (int N = 1; N <= 6; N++) {
        double X_N = cdr_cstr_series_conversion(k, tau_cstr, N);
        printf("  N=%d: X = %.4f\n", N, X_N);
    }
    printf("  PFR limit: X = %.4f\n", 1.0 - exp(-k * tau_cstr));

    return 0;
}

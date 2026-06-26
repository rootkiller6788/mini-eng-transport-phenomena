/**
 * example_condenser.c
 * ====================
 * Example: Power plant condenser design using Nusselt film condensation theory.
 *
 * Demonstrates L6 problem solving for heat exchanger design with
 * interphase transport (steam condensation on horizontal tubes).
 *
 * Scenario: Design the condenser for a 500 MWe steam turbine.
 * The condenser must reject ~1000 MW_th of waste heat from the
 * Rankine cycle to cooling water.
 *
 * Reference: Incropera & DeWitt (2007) "Fundamentals of Heat and Mass Transfer"
 *            Nusselt (1916) - Film condensation theory
 */

#include "interphase_types.h"
#include "interphase_transport.h"
#include <stdio.h>
#include <math.h>

static void condenser_example(void) {
    printf("=== Power Plant Condenser Design Example ===\n\n");

    /* Design parameters */
    double P_thermal_reject = 1000.0e6;  /* 1000 MW_th to reject */
    double T_sat = 310.0;   /* K (~37C condensing at ~6 kPa) */
    double T_cw_in = 293.15;  /* Cooling water inlet (20C) */
    double T_cw_out = 303.15; /* Cooling water outlet (30C) */

    /* Tube geometry */
    double D_tube = 0.0254;   /* 1 inch OD */
    double L_tube = 15.0;     /* 15 m tube length */
    int N_tubes_per_column = 20;  /* tubes in vertical column */

    /* Properties (water/steam at condensing conditions) */
    double rho_L = 993.0;     /* kg/m^3 (water at 37C) */
    double rho_v = 0.044;     /* kg/m^3 (steam at 37C, 6 kPa) */
    double k_L   = 0.628;     /* W/(m*K) */
    double h_fg  = 2.41e6;    /* J/kg (latent heat at 37C) */
    double mu_L  = 6.9e-4;    /* Pa*s */
    double cp_cw = 4181.0;    /* J/(kg*K) cooling water */

    /* Condensation HTC (Nusselt for horizontal tube) */
    double h_cond, Q_per_tube;
    condensation_horizontal_tube(D_tube, T_sat,
                                  (T_sat - 5.0),  /* T_wall ~ 5K subcooling */
                                  rho_L, rho_v, k_L, h_fg, mu_L,
                                  N_tubes_per_column,
                                  &h_cond, &Q_per_tube);

    printf("Condensation heat transfer:\n");
    printf("  Nusselt HTC:       %.0f W/(m^2*K)\n", h_cond);
    printf("  Heat per tube:     %.0f W/m\n", Q_per_tube);
    printf("  Heat per tube:     %.1f kW\n", Q_per_tube / 1000.0);
    printf("\n");

    /* Required number of tubes */
    double Q_per_tube_total = Q_per_tube * L_tube;
    int N_tubes = (int) ceil(P_thermal_reject / Q_per_tube_total);
    double A_total = N_tubes * M_PI * D_tube * L_tube;

    printf("Condenser sizing:\n");
    printf("  Heat rejection:    %.0f MW_th\n", P_thermal_reject / 1.0e6);
    printf("  Tubes required:    %d\n", N_tubes);
    printf("  Total area:        %.0f m^2\n", A_total);
    printf("  Area/MW_th:        %.1f m^2/MW\n", A_total / (P_thermal_reject / 1.0e6));
    printf("\n");

    /* Cooling water flow rate */
    double m_dot_cw = P_thermal_reject / (cp_cw * (T_cw_out - T_cw_in));
    double Q_cw = m_dot_cw / rho_L;

    printf("Cooling water requirements:\n");
    printf("  Mass flow rate:    %.0f kg/s\n", m_dot_cw);
    printf("  Volumetric flow:   %.2f m^3/s\n", Q_cw);
    printf("                      (%.0f m^3/h)\n", Q_cw * 3600.0);
    printf("\n");

    /* Overall HTC estimate */
    double U_estimate = 2500.0;  /* typical condenser U ~ 2000-3000 W/(m^2*K) */
    double LMTD = (T_sat - T_cw_in);  /* simplified: condensing side constant T */

    /* Verify: Q = U * A * LMTD */
    double Q_check = U_estimate * A_total * LMTD;
    printf("Design verification (U = %.0f W/(m^2*K)):\n", U_estimate);
    printf("  Q = U*A*LMTD = %.0f MW_th\n", Q_check / 1.0e6);
    printf("  Target:          %.0f MW_th\n", P_thermal_reject / 1.0e6);
    printf("  Margin:          %.1f %%\n",
           100.0 * (Q_check - P_thermal_reject) / P_thermal_reject);

    printf("\nExample completed successfully.\n");
}

int main(void) {
    condenser_example();
    return 0;
}

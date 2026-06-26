/*
 * example_flat_plate_heat_transfer.c — Flat Plate Convective Heat Transfer
 *
 * Engineering Problem (L6): Compute heat transfer from a heated flat plate
 * in cross-flow — the canonical external convection problem.
 *
 * This models applications ranging from solar panel cooling to
 * aircraft wing de-icing to semiconductor wafer cooling.
 */

#include "../include/nusselt_number.h"
#include "../include/dimensionless_numbers.h"
#include "../include/transport_correlations.h"
#include <stdio.h>
#include <math.h>

int main(void)
{
    printf("=== Flat Plate Forced Convection Analysis ===\n");
    printf("Reference: Incropera & DeWitt (2007), Ch. 7\n\n");

    /* Problem: Solar panel in wind — 1 m long, 30 degC above ambient */
    printf("--- Problem 1: Solar Panel (Photovoltaic) Cooling ---\n");

    double T_plate = 333.0;   /* 60 degC — hot solar panel */
    double T_air   = 303.0;   /* 30 degC — summer day */
    double U_wind  = 5.0;     /* 5 m/s — moderate breeze */
    double L_panel = 1.0;     /* 1 m panel length in wind direction */
    double W_panel = 1.6;     /* 1.6 m wide */

    /* Air properties at film temperature T_f = (333+303)/2 = 318 K */
    double rho_air = 1.10;    /* kg/m^3 */
    double mu_air  = 1.93e-5; /* Pa.s */
    double k_air   = 0.0276;  /* W/(m.K) */
    double cp_air  = 1007.0;  /* J/(kg.K) */
    double Pr_air  = 0.704;
    double nu_air  = mu_air / rho_air;

    /* Reynolds number */
    double Re_L = rho_air * U_wind * L_panel / mu_air;
    printf("  Re_L = %.0f (%.0f%% of critical)\n",
           Re_L, Re_L / 500000.0 * 100.0);

    /* Laminar heat transfer */
    double Nu_lam = nu_flat_plate_laminar_avg(Re_L, Pr_air);
    double h_lam  = h_from_nusselt(Nu_lam, L_panel, k_air);
    double q_lam  = h_lam * (T_plate - T_air);
    double Q_lam  = q_lam * L_panel * W_panel;

    printf("  Laminar: Nu_L = %.1f, h = %.2f W/m^2K, Q = %.1f W\n",
           Nu_lam, h_lam, Q_lam);

    /* Check transition location */
    double x_crit = 500000.0 * nu_air / U_wind;
    printf("  Transition at x_crit = %.2f m (%.0f%% of plate)\n",
           x_crit, x_crit / L_panel * 100.0);

    /* Mixed laminar-turbulent */
    if (Re_L > 500000.0) {
        double Nu_mixed = nu_flat_plate_mixed(Re_L, Pr_air, 500000.0);
        double h_mixed  = h_from_nusselt(Nu_mixed, L_panel, k_air);
        double Q_mixed  = h_mixed * (T_plate - T_air) * L_panel * W_panel;
        printf("  Mixed:   Nu_L = %.1f, h = %.2f W/m^2K, Q = %.1f W\n",
               Nu_mixed, h_mixed, Q_mixed);
    }

    /* Air mass flow over panel */
    double m_dot_air = rho_air * U_wind * L_panel * W_panel;
    double delta_T_air = Q_lam / (m_dot_air * cp_air);
    printf("  Air temperature rise: DeltaT_air ≈ %.2f K\n", delta_T_air);
    printf("  Solar panel spacing > 0.3 m recommended for natural ventilation.\n\n");

    /* Problem 2: Aircraft wing anti-icing */
    printf("--- Problem 2: Boeing 737 Wing Anti-Icing ---\n");

    double T_wing   = 280.0;   /* 7 degC — heated wing surface */
    double T_amb    = 253.0;   /* -20 degC — cruise altitude */
    double U_flight = 230.0;   /* 230 m/s ~ Mach 0.78 at altitude */
    double L_wing   = 3.5;     /* chord length */
    double W_wing   = 15.0;    /* span of heated section */

    /* Air at 10 km altitude (~250 K) */
    double rho_alt = 0.41;    /* kg/m^3 */
    double mu_alt  = 1.42e-5; /* Pa.s */
    double k_alt   = 0.022;   /* W/(m.K) */
    double Pr_alt  = 0.72;

    double Re_wing = rho_alt * U_flight * L_wing / mu_alt;
    printf("  Re_chord = %.2e (fully turbulent)\n", Re_wing);

    double Nu_wing = nu_flat_plate_turbulent_avg(Re_wing, Pr_alt);
    double h_wing  = h_from_nusselt(Nu_wing, L_wing, k_alt);
    double q_wing  = h_wing * (T_wing - T_amb);
    double Q_wing  = q_wing * L_wing * W_wing / 1000.0; /* kW */

    printf("  Nu_L = %.0f, h = %.1f W/m^2K\n", Nu_wing, h_wing);
    printf("  Heat required per wing: Q = %.1f kW\n", Q_wing);
    printf("  Engine bleed air: ~200 degC, sufficient for anti-icing.\n");
    printf("  Note: Actual system uses piccolo tube + hot air injection.\n\n");

    /* Problem 3: iPhone heat dissipation */
    printf("--- Problem 3: iPhone 15 Thermal Management ---\n");

    double T_phone  = 313.0;   /* 40 degC — surface temperature */
    double T_room   = 298.0;   /* 25 degC — room temperature */
    double U_conv   = 0.3;     /* 0.3 m/s — natural convection + slight air motion */
    double L_phone  = 0.15;    /* 150 mm phone length */

    /* Air properties at 30 degC */
    double rho_r = 1.16;
    double mu_r  = 1.86e-5;
    double k_r   = 0.026;
    double Pr_r  = 0.71;

    double Re_phone = rho_r * U_conv * L_phone / mu_r;
    printf("  Re = %.0f (laminar)\n", Re_phone);

    /* But natural convection dominates: Ra-based */
    double g = 9.81;
    double beta = 1.0 / 303.0;
    double nu_r = mu_r / rho_r;
    double alpha_r = k_r / (rho_r * 1007.0);
    double Gr = g * beta * (T_phone - T_room) * L_phone * L_phone * L_phone
                / (nu_r * nu_r);
    double Ra = Gr * Pr_r;
    printf("  Gr = %.2e, Ra = %.2e\n", Gr, Ra);

    double Nu_nat = nu_natural_vertical_plate_churchill_chu(Ra, Pr_r);
    double h_nat  = h_from_nusselt(Nu_nat, L_phone, k_r);
    double A_phone = L_phone * 0.075; /* 150×75 mm */
    double Q_phone = h_nat * A_phone * (T_phone - T_room);

    printf("  Natural convection: Nu = %.1f, h = %.2f W/m^2K\n", Nu_nat, h_nat);
    printf("  Heat dissipated: Q = %.2f W\n", Q_phone);
    printf("  iPhone 15 TDP ~ 6-8 W; requires additional conduction to frame.\n");
    printf("  Apple uses graphite sheets + aluminum frame as heat spreader.\n");

    printf("\n=== Analysis Complete ===\n");
    printf("Key insight: Nu = f(Re, Pr, geometry) — the same dimensionless\n");
    printf("correlation predicts heat transfer from mm-scale chips to km-scale aircraft.\n");

    return 0;
}

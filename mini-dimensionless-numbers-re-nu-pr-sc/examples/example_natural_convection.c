/*
 * example_natural_convection.c — Natural Convection Problems
 *
 * Engineering Problems (L6): Natural convection governs passive cooling,
 * building ventilation, electronics thermal management, and geophysical flows.
 *
 * Demonstrates the use of Gr, Ra, and Nu for buoyancy-driven problems.
 */

#include "../include/nusselt_number.h"
#include "../include/transport_correlations.h"
#include "../include/dimensionless_numbers.h"
#include <stdio.h>
#include <math.h>

int main(void)
{
    printf("=== Natural Convection Analysis ===\n");
    printf("Reference: Churchill & Chu (1975), Bejan (2013)\n\n");

    /* Problem 1: Room radiator (vertical panel heater) */
    printf("--- Problem 1: Household Radiator ---\n");

    double T_rad  = 333.0;    /* 60 degC — radiator surface */
    double T_room = 293.0;    /* 20 degC — room temperature */
    double H_rad  = 0.6;      /* 0.6 m height */
    double W_rad  = 1.2;      /* 1.2 m width */

    /* Air properties at film temperature */
    double g = 9.81;
    double beta = 1.0 / 313.0;  /* ideal gas */
    double nu_a = 1.7e-5;       /* m^2/s at ~40 degC */
    double alpha_a = 2.4e-5;    /* m^2/s */
    double k_a = 0.027;         /* W/(m.K) */
    double Pr_a = nu_a / alpha_a;

    double delta_T = T_rad - T_room;
    double Gr_H = g * beta * delta_T * H_rad * H_rad * H_rad / (nu_a * nu_a);
    double Ra_H = Gr_H * Pr_a;

    printf("  DeltaT = %.0f K\n", delta_T);
    printf("  Gr_H = %.2e, Ra_H = %.2e\n", Gr_H, Ra_H);

    double Nu_CC = nu_natural_vertical_plate_churchill_chu(Ra_H, Pr_a);
    double h_rad = h_from_nusselt(Nu_CC, H_rad, k_a);
    double Q_rad = h_rad * H_rad * W_rad * delta_T;

    printf("  Nu_H = %.1f (Churchill-Chu)\n", Nu_CC);
    printf("  h = %.2f W/m^2K\n", h_rad);
    printf("  Q = %.1f W\n", Q_rad);
    printf("  (Typical radiator: 500-2000 W — this is panel-only, no fins)\n\n");

    /* Problem 2: DC motor natural cooling (Toyota Prius MG2) */
    printf("--- Problem 2: Toyota Prius Motor-Generator Cooling ---\n");

    double T_amb  = 313.0;    /* 40 degC — under-hood ambient */
    double P_loss = 500.0;    /* 500 W — electrical losses at ~90% efficiency */
    double D_motor = 0.25;    /* 250 mm motor diameter */
    double L_motor = 0.35;    /* 350 mm motor length */

    double T_surf, Ra_D, Nu_D;
    dc_motor_natural_cooling(T_amb, P_loss, D_motor, L_motor,
                             k_a, nu_a, alpha_a, &T_surf, &Ra_D, &Nu_D);

    printf("  Motor housing: D=%.0f mm, L=%.0f mm\n", D_motor*1000, L_motor*1000);
    printf("  Power loss: %.0f W\n", P_loss);
    printf("  Ra_D = %.2e, Nu_D = %.1f\n", Ra_D, Nu_D);
    printf("  Surface temperature: %.1f  degC\n", T_surf - 273.15);
    printf("  Toyota uses oil spray + water jacket for > 50 kW peak.\n\n");

    /* Problem 3: Enclosure — double-pane window */
    printf("--- Problem 3: Double-Pane Window Heat Loss ---\n");

    double T_indoor  = 295.0;   /* 22 degC indoor */
    double T_outdoor = 273.0;   /* 0 degC outdoor */
    double gap       = 0.012;   /* 12 mm air gap */
    double H_window  = 1.5;     /* 1.5 m window height */
    double k_argon   = 0.016;   /* argon fill (better than air 0.026) */

    double delta_T_window = T_indoor - T_outdoor;
    double beta_w = 1.0 / 288.0;
    double Ra_gap = g * beta_w * delta_T_window * gap * gap * gap
                    / (nu_a * alpha_a);
    double aspect_w = H_window / gap;

    double Nu_gap, h_gap, q_gap;
    natural_convection_enclosure(k_argon, gap, H_window,
                                 T_indoor, T_outdoor,
                                 Ra_gap, Pr_a, 90.0,
                                 &Nu_gap, &h_gap, &q_gap);

    printf("  Air gap: %.0f mm, DeltaT = %.0f K\n", gap*1000, delta_T_window);
    printf("  Ra_gap = %.2e -> %s\n", Ra_gap,
           (Ra_gap > 1708.0) ? "convection active" : "conduction only");
    printf("  Nu = %.2f (Nu=1 = pure conduction)\n", Nu_gap);
    printf("  Heat loss: q\" = %.1f W/m^2\n", q_gap);
    printf("  For 1.5×1.0 m window: Q = %.1f W\n", q_gap * H_window * 1.0);
    printf("  Triple-pane with low-E coating reduces to < 0.5 W/m^2K.\n");
    printf("  NASA uses multi-layer insulation (MLI) for spacecraft — same physics.\n\n");

    /* Problem 4: Mars rover thermal — natural convection in thin CO2 */
    printf("--- Problem 4: Mars Rover Heat Rejection ---\n");

    /* Mars atmosphere: 95% CO2, 600 Pa, ~210 K */
    double T_mars_amb = 210.0;    /* -63 degC, equatorial daytime */
    double T_rover    = 273.0;    /* 0 degC — minimum electronics temp */
    double P_rhg      = 2.0;      /* 2 W — single Radioisotope Heater Unit */

    double T_final, Q_lost;
    mars_rover_night_survival(P_rhg, 1.0, T_rover + 20.0, 50400.0,
                               20000.0, &T_final, &Q_lost);

    printf("  Mars night: 14 hours (~50400 s)\n");
    printf("  Heater power: %.0f W (RHU — Pu-238)\n", P_rhg);
    printf("  Thin atmosphere (600 Pa): h ≈ 1-2 W/m^2K\n");
    printf("  Final temperature: %.1f K (%.0f degC)\n", T_final, T_final - 273.15);
    printf("  Total heat lost: %.0f kJ\n", Q_lost / 1000.0);
    printf("  Perseverance uses 8 RHUs + MMRTG waste heat.\n\n");

    /* Summary table */
    printf("=== Natural Convection Regime Map ===\n");
    printf("  Ra          | Convection Type\n");
    printf("  ------------|----------------\n");
    printf("  < 1708      | Conduction only (enclosure)\n");
    printf("  10^3 - 10⁸   | Laminar natural convection\n");
    printf("  10⁸ - 10⁹   | Transitional\n");
    printf("  10⁹ - 10¹^2  | Turbulent natural convection\n");
    printf("  > 10¹^2      | Fully turbulent (geophysical)\n\n");

    printf("Key insight: Ra = Gr.Pr = (gbetaDeltaT L^3/nu^2).(nu/alpha) = gbetaDeltaT L^3/(nualpha).\n");
    printf("  Ra > Ra_crit -> natural convection begins spontaneously.\n");
    printf("  Nu ~ Ra^(1/4) (laminar), Nu ~ Ra^(1/3) (turbulent).\n");

    return 0;
}

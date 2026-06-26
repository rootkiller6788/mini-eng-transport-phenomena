/**
 * example_droplet_evaporation.c
 * ==============================
 * Example: Droplet evaporation with interphase transport.
 *
 * Demonstrates L6 problem solving for spray drying / fuel spray
 * evaporation. A spherical droplet evaporates in hot gas,
 * governed by simultaneous heat and mass transfer at the
 * liquid-vapor interface.
 *
 * This is the d^2-law (Godsave, 1953; Spalding, 1953):
 *   d^2(t) = d_0^2 - K * t
 * where K is the evaporation constant.
 *
 * Reference: Sirignano (2010) "Fluid Dynamics and Transport of Droplets and Sprays"
 *            (2nd ed., Cambridge University Press)
 */

#include "interphase_types.h"
#include "interphase_transport.h"
#include "interphase_boundary_conditions.h"
#include "interphase_jump_conditions.h"
#include <stdio.h>
#include <math.h>

static void droplet_evaporation_example(void) {
    printf("=== Droplet Evaporation (d^2-Law) Example ===\n\n");

    /* Initial conditions */
    double d0 = 100.0e-6;   /* initial droplet diameter (m), 100 um */
    double T_drop = 300.0;  /* initial droplet temperature (K) */
    double T_gas  = 600.0;  /* hot gas temperature (K) */

    /* Gas properties (air at 600K) */
    double rho_g = 0.6;     /* kg/m^3 */
    double mu_g  = 3.0e-5;  /* Pa*s */
    double k_g   = 0.045;   /* W/(m*K) */
    double cp_g  = 1050.0;  /* J/(kg*K) */
    double D_AB  = 3.0e-5;  /* m^2/s, vapor diffusivity in air */

    /* Liquid properties (water) */
    double rho_L = 997.0;   /* kg/m^3 */
    double h_fg  = 2.26e6;  /* J/kg */
    double cp_L  = 4181.0;  /* J/(kg*K) */

    /* Interface equilibrium */
    double p_sat = 101325.0; /* Pa, saturation at droplet temperature (approx) */
    double M_v   = 0.018;    /* kg/mol water */
    double R     = 8.314;    /* J/(mol*K) */

    printf("Initial conditions:\n");
    printf("  Droplet diameter:  %.0f um\n", d0 * 1.0e6);
    printf("  Droplet temp:      %.0f K\n", T_drop);
    printf("  Gas temp:          %.0f K\n", T_gas);
    printf("\n");

    /* Mass fraction of vapor at droplet surface */
    double Y_v_surface = p_sat * M_v / (R * T_drop * rho_g);
    double Y_v_inf     = 0.0;  /* dry gas far away */

    printf("Vapor mass fraction:\n");
    printf("  At surface:        %.4f\n", Y_v_surface);
    printf("  Far field:         %.4f\n", Y_v_inf);
    printf("\n");

    /* Spalding mass transfer number */
    double B_M = (Y_v_surface - Y_v_inf) / (1.0 - Y_v_surface);
    printf("Spalding numbers:\n");
    printf("  B_M (mass):        %.4f\n", B_M);

    /* Spalding heat transfer number */
    double T_boil = 373.15;  /* K, boiling point */
    double B_T = cp_g * (T_gas - T_drop) / h_fg;
    printf("  B_T (thermal):     %.4f\n", B_T);
    printf("\n");

    /* Sherwood and Nusselt numbers for sphere */
    /* Relative velocity ~ 0 (quiescent evaporation) */
    double Re_p = 0.0;
    double Sh, Nu;
    sherwood_single_sphere(Re_p, mu_g/(rho_g*D_AB), &Sh);

    /* For pure diffusion: Sh = 2, Nu = 2 */
    printf("Transport coefficients:\n");
    printf("  Sh (pure diff):    %.1f\n", Sh);
    printf("  k_m = Sh*D/d:      %.3f m/s\n", Sh * D_AB / d0);

    double htc = Sh * k_g / d0;  /* using Sh in place of Nu for laminar case */
    printf("  htc = Nu*k/d:      %.1f W/(m^2*K)\n", 2.0 * k_g / d0);
    printf("\n");

    /* Evaporation constant K from d^2-law:
     * K = 8*k_g*ln(1+B_T)/(rho_L*cp_g) = 8*D_AB*rho_g*ln(1+B_M)/rho_L
     */
    double K_thermal = 8.0 * k_g * log(1.0 + B_T) / (rho_L * cp_g);
    double K_mass    = 8.0 * D_AB * rho_g * log(1.0 + B_M) / rho_L;

    printf("Evaporation constant (d^2-law):\n");
    printf("  K (thermal):       %.2e m^2/s\n", K_thermal);
    printf("  K (mass):          %.2e m^2/s\n", K_mass);
    printf("  K (average):       %.2e m^2/s\n", 0.5*(K_thermal+K_mass));
    printf("\n");

    /* Lifetime estimate */
    double K_avg = 0.5 * (K_thermal + K_mass);
    double t_life = (d0 * d0) / K_avg;

    printf("Droplet lifetime:\n");
    printf("  t_evap = d0^2/K:   %.4f s\n", t_life);
    printf("                      %.1f ms\n", t_life * 1000.0);
    printf("\n");

    /* Time evolution (first 5 points) */
    printf("Time evolution (d^2-law):\n");
    printf("  Time (ms)    d (um)    d^2/d0^2\n");
    printf("  ---------    ------    ---------\n");
    for (int i = 0; i <= 5; i++) {
        double t = i * t_life / 5.0;
        double d_sq = d0*d0 - K_avg * t;
        if (d_sq < 0.0) d_sq = 0.0;
        double d = sqrt(d_sq);
        printf("  %8.1f    %6.1f    %8.4f\n",
               t * 1000.0, d * 1.0e6, d_sq / (d0*d0));
    }

    printf("\nExample completed successfully.\n");
}

int main(void) {
    droplet_evaporation_example();
    return 0;
}

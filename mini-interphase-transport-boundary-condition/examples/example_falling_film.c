/**
 * example_falling_film.c
 * =======================
 * Example: CO2 absorption in a falling film absorber.
 *
 * Demonstrates L6 problem solving using interphase transport theory.
 *
 * Scenario: A vertical wetted-wall column absorbs CO2 from air into water.
 * This is a standard chemical engineering unit operation for gas scrubbing.
 *
 * Reference: Danckwerts (1970) "Gas-Liquid Reactions"
 *            Treybal (1980) "Mass-Transfer Operations"
 */

#include "interphase_types.h"
#include "interphase_transport.h"
#include "interphase_boundary_conditions.h"
#include "interphase_jump_conditions.h"
#include <stdio.h>
#include <math.h>

/**
 * Compute the absorption rate for a falling film contactor.
 */
static void falling_film_example(void) {
    printf("=== Falling Film CO2 Absorption Example ===\n\n");

    /* Geometric parameters */
    double L = 1.0;      /* film length (m) */
    double W = 0.1;      /* film width (m) */
    double Q_L = 1.0e-6; /* liquid flow rate (m^3/s), ~ 1 mL/s */

    /* Physical properties (water at 298K) */
    double rho_L = 997.0;    /* kg/m^3 */
    double mu_L  = 8.9e-4;   /* Pa*s */
    double g     = 9.81;     /* m/s^2 */

    /* Mass transfer properties (CO2 in water) */
    double D_CO2  = 1.92e-9;  /* m^2/s, CO2 diffusivity in water */
    double H_CO2  = 1.67e8;   /* Pa, Henry's constant */
    double p_CO2  = 400.0;    /* Pa, atmospheric CO2 partial pressure (~400 ppm) */
    double C_in   = 0.0;      /* mol/m^3, pure water inlet */

    /* Film thickness from Nusselt solution */
    double delta_film = pow(3.0 * mu_L * Q_L / (rho_L * g * W), 1.0/3.0);
    double u_surf = rho_L * g * delta_film * delta_film / (2.0 * mu_L);
    double t_contact = L / u_surf;
    double Re_film = 4.0 * rho_L * Q_L / (mu_L * W);

    /* Mass transfer coefficient from penetration theory */
    double k_L = 2.0 * sqrt(D_CO2 / (M_PI * t_contact));

    printf("Film hydrodynamics:\n");
    printf("  Film thickness:    %.3f mm\n", delta_film * 1000.0);
    printf("  Surface velocity:  %.3f m/s\n", u_surf);
    printf("  Contact time:      %.4f s\n", t_contact);
    printf("  Film Reynolds:     %.1f\n", Re_film);
    printf("\n");

    printf("Mass transfer:\n");
    printf("  k_L (penetration): %.2e m/s\n", k_L);
    printf("  CO2 diffusivity:   %.2e m^2/s\n", D_CO2);
    printf("  Henry constant:    %.2e Pa\n", H_CO2);
    printf("  C* (equilibrium):  %.2e mol/m^3\n", p_CO2 / H_CO2);
    printf("\n");

    /* Concentration profile solution */
    double C_star = p_CO2 / H_CO2;
    double a_specific = 1.0 / delta_film;
    double exponent = -k_L * a_specific * L / u_surf;
    double C_out = C_star + (C_in - C_star) * exp(exponent);
    double N_abs = Q_L * (C_out - C_in);

    printf("Absorption results:\n");
    printf("  Exit concentration: %.4e mol/m^3\n", C_out);
    printf("  Saturation:         %.2f %%\n", 100.0 * C_out / C_star);
    printf("  Absorption rate:    %.4e mol/s\n", N_abs);
    printf("  Absorption flux:    %.4e mol/(m^2*s)\n", N_abs / (L * W));
    printf("\n");

    /* Dimensionless numbers */
    double Sc = mu_L / (rho_L * D_CO2);
    double Sh;
    sherwood_falling_film(Re_film, Sc, L / delta_film, &Sh);

    printf("Dimensionless analysis:\n");
    printf("  Schmidt number:     %.0f\n", Sc);
    printf("  Sherwood number:    %.1f\n", Sh);
    printf("  k_L (from Sh):      %.2e m/s\n", Sh * D_CO2 / delta_film);
    printf("\n");
}

int main(void) {
    falling_film_example();
    printf("Example completed successfully.\n");
    return 0;
}

/**
 * @file test_diffusion.c
 * @brief Tests for diffusion transport models.
 */

#include "../include/cdr_diffusion.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define TOL 1e-10

static void test_effective_diffusivity(void)
{
    /* D_eff = porosity/tortuosity * D_bulk */
    double D_eff = cdr_diffusion_effective(1.0e-5, 0.4, 3.0);
    double expected = 1.0e-5 * 0.4 / 3.0;
    assert(fabs(D_eff - expected) < TOL);

    /* Edge cases */
    assert(cdr_diffusion_effective(1.0e-5, 0.0, 3.0) == 0.0);
    assert(cdr_diffusion_effective(1.0e-5, 1.1, 3.0) == 0.0); /* porosity > 1 */
    assert(cdr_diffusion_effective(1.0e-5, 0.4, 0.5) == 0.0); /* tortuosity < 1 */

    printf("  [PASS] test_effective_diffusivity\n");
}

static void test_knudsen_diffusion(void)
{
    /* D_K ~ pore_diameter * sqrt(T/M) */
    double D_K = cdr_diffusion_knudsen(1.0e-9, 300.0, 0.028);
    assert(D_K > 0.0);
    assert(D_K < 1.0e-4);

    /* Linear in pore diameter */
    double D_K2 = cdr_diffusion_knudsen(2.0e-9, 300.0, 0.028);
    assert(fabs(D_K2 / D_K - 2.0) < 0.01);

    assert(cdr_diffusion_knudsen(-1.0, 300.0, 0.028) == 0.0);
    printf("  [PASS] test_knudsen_diffusion\n");
}

static void test_bosanquet(void)
{
    double D_bos = cdr_diffusion_bosanquet(1.0e-5, 1.0e-5);
    assert(fabs(D_bos - 5.0e-6) < TOL);

    assert(cdr_diffusion_bosanquet(0.0, 1.0e-5) == 0.0);
    printf("  [PASS] test_bosanquet\n");
}

static void test_fuller(void)
{
    /* O2-N2 diffusion at STP */
    /* V_O2 ~ 16.6, V_N2 ~ 17.9 cm3/mol, M_O2=32, M_N2=28 */
    double D = cdr_diffusion_fuller(298.0, 1.0, 32.0, 28.0, 16.6, 17.9);
    /* Expected ~ 2.0e-5 m2/s */
    assert(D > 1.0e-5 && D < 3.0e-5);

    assert(cdr_diffusion_fuller(0.0, 1.0, 32.0, 28.0, 16.6, 17.9) == 0.0);
    printf("  [PASS] test_fuller\n");
}

static void test_wilke_chang(void)
{
    /* Oxygen in water at 25C */
    double D = cdr_diffusion_wilke_chang(298.0, 1.0, 18.0, 25.6, 2.6);
    /* Expected ~ 2e-9 m2/s */
    assert(D > 1.0e-10 && D < 1.0e-8);

    printf("  [PASS] test_wilke_chang\n");
}

static void test_fick_first_law_1d(void)
{
    /* J = D * (C1 - C2) / dx */
    double J = cdr_fick_first_law_1d(1.0e-5, 1.0, 0.5, 0.01);
    double expected = 1.0e-5 * (1.0 - 0.5) / 0.01;
    assert(fabs(J - expected) < TOL);

    /* Zero flux when concentrations equal */
    J = cdr_fick_first_law_1d(1.0e-5, 0.5, 0.5, 0.01);
    assert(fabs(J) < TOL);

    assert(cdr_fick_first_law_1d(1.0, 1.0, 0.0, 0.0) == 0.0);
    printf("  [PASS] test_fick_first_law_1d\n");
}

static void test_fick_first_law_3d(void)
{
    DiffusionFlux flux;
    cdr_fick_first_law_3d(1.0e-5, 100.0, 0.0, 0.0, &flux);

    assert(fabs(flux.J_diff_x + 0.001) < TOL);
    assert(fabs(flux.J_diff_y) < TOL);
    assert(fabs(flux.magnitude - 0.001) < 1e-8);

    /* NULL safety */
    cdr_fick_first_law_3d(1.0, 1.0, 0.0, 0.0, NULL);
    printf("  [PASS] test_fick_first_law_3d\n");
}

static void test_fick_second_law_1d(void)
{
    /* Set up a Gaussian-like concentration profile */
    size_t n = 101;
    ConcentrationProfile1D profile;
    profile.C = (double *)malloc(n * sizeof(double));
    profile.x = (double *)malloc(n * sizeof(double));
    profile.n_points = n;
    profile.L = 1.0;
    profile.dx = 1.0 / (n - 1);
    profile.t = 0.0;

    /* Initial: high concentration at center */
    for (size_t i = 0; i < n; i++) {
        profile.x[i] = i * profile.dx;
        double xc = profile.x[i] - 0.5;
        profile.C[i] = exp(-xc * xc / 0.01);
    }

    /* Test stability check */
    double dt = cdr_fick_second_law_1d_step(&profile, 0.1, 1.0);
    /* Should return max stable dt since r would be too large */
    assert(dt > 0.0);

    /* Test with stable dt */
    double D = 0.01;
    double dx = profile.dx;
    double stable_dt = 0.4 * dx * dx / D;
    double result = cdr_fick_second_law_1d_step(&profile, D, stable_dt);
    assert(result == stable_dt);
    assert(profile.t == stable_dt);

    /* NULL safety */
    assert(cdr_fick_second_law_1d_step(NULL, 1.0, 0.01) == 0.0);

    free(profile.C);
    free(profile.x);
    printf("  [PASS] test_fick_second_law_1d\n");
}

static void test_semi_infinite_diffusion(void)
{
    /* At x=0, C = C_surface */
    double C = cdr_diffusion_semi_infinite(1.0, 0.0, 1.0e-9, 0.0, 100.0);
    assert(fabs(C - 0.0) < TOL);

    /* Far from surface, C ? C0 */
    C = cdr_diffusion_semi_infinite(1.0, 0.0, 1.0e-9, 1.0, 1.0);
    assert(C > 0.99);

    /* Penetration depth: ~4*sqrt(D*t) */
    double delta = cdr_diffusion_penetration_depth(1.0e-9, 100.0);
    C = cdr_diffusion_semi_infinite(1.0, 0.0, 1.0e-9, delta, 100.0);
    /* Should be ~0.99 * 1.0 */
    assert(C > 0.97);

    assert(cdr_diffusion_semi_infinite(1.0, 0.0, 0.0, 0.1, 1.0) == 1.0);
    printf("  [PASS] test_semi_infinite_diffusion\n");
}

static void test_steady_state_solutions(void)
{
    /* Membrane flux */
    double J = cdr_diffusion_membrane_flux(1.0e-5, 10.0, 0.0, 0.001);
    assert(fabs(J - 0.1) < TOL);

    assert(cdr_diffusion_membrane_flux(1.0, 1.0, 0.0, 0.0) == 0.0);

    /* Cylindrical steady-state */
    double C_cyl = cdr_diffusion_cylindrical_ss(0.015, 0.01, 0.02, 10.0, 0.0);
    assert(C_cyl > 0.0 && C_cyl < 10.0);

    /* Spherical steady-state */
    double C_sph = cdr_diffusion_spherical_ss(0.015, 0.01, 0.02, 10.0, 0.0);
    assert(C_sph > 0.0 && C_sph < 10.0);

    printf("  [PASS] test_steady_state_solutions\n");
}

static void test_diffusion_timescale(void)
{
    double tau = cdr_diffusion_timescale(0.001, 1.0e-9);
    assert(fabs(tau - 1000.0) < TOL);

    assert(cdr_diffusion_timescale(0.0, 1.0) == 0.0);
    printf("  [PASS] test_diffusion_timescale\n");
}

int main(void)
{
    printf("CDR Diffusion Tests:\n");
    test_effective_diffusivity();
    test_knudsen_diffusion();
    test_bosanquet();
    test_fuller();
    test_wilke_chang();
    test_fick_first_law_1d();
    test_fick_first_law_3d();
    test_fick_second_law_1d();
    test_semi_infinite_diffusion();
    test_steady_state_solutions();
    test_diffusion_timescale();
    printf("All diffusion tests passed.\n");
    return 0;
}

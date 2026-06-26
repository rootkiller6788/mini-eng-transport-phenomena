/**
 * @file test_effectiveness.c
 * @brief Tests for catalyst effectiveness factor and diffusion limitations.
 */

#include "../include/cdr_effectiveness.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define TOL 1e-8

static void test_effectiveness_slab_asymptotics(void)
{
    /* phi -> 0: eta -> 1 */
    double eta = cdr_effectiveness_slab(0.001);
    assert(fabs(eta - 1.0) < 0.01);

    /* phi -> inf: eta -> 1/phi */
    double eta_large = cdr_effectiveness_slab(100.0);
    assert(fabs(eta_large - 0.01) < 0.001);

    /* phi = 1.0: eta = tanh(1)/1 ? 0.7616 */
    double eta_1 = cdr_effectiveness_slab(1.0);
    assert(fabs(eta_1 - 0.761594) < 1e-4);

    /* phi <= 0 */
    assert(fabs(cdr_effectiveness_slab(0.0) - 1.0) < TOL);
    assert(fabs(cdr_effectiveness_slab(-1.0) - 1.0) < TOL);

    printf("  [PASS] test_effectiveness_slab_asymptotics\n");
}

static void test_effectiveness_cylinder(void)
{
    /* phi -> 0: eta -> 1 */
    double eta = cdr_effectiveness_cylinder(0.001);
    assert(fabs(eta - 1.0) < 0.01);

    /* phi -> inf: eta -> 2/phi */
    double eta_large = cdr_effectiveness_cylinder(100.0);
    assert(fabs(eta_large - 0.02) < 0.001);

    /* Compare with slab: cylinder has slightly higher eta for same phi */
    double eta_slab = cdr_effectiveness_slab(1.0);
    double eta_cyl = cdr_effectiveness_cylinder(1.0);
    assert(eta_cyl > eta_slab);

    printf("  [PASS] test_effectiveness_cylinder\n");
}

static void test_effectiveness_sphere(void)
{
    /* phi -> 0: eta -> 1 */
    double eta = cdr_effectiveness_sphere(0.001);
    assert(fabs(eta - 1.0) < 0.01);

    /* phi -> inf: eta -> 3/phi */
    double eta_large = cdr_effectiveness_sphere(100.0);
    assert(fabs(eta_large - 0.03) < 0.001);

    /* For same phi, eta_sphere > eta_cylinder > eta_slab */
    double eta_slab = cdr_effectiveness_slab(1.0);
    double eta_cyl = cdr_effectiveness_cylinder(1.0);
    double eta_sph = cdr_effectiveness_sphere(1.0);
    assert(eta_sph > eta_cyl);
    assert(eta_cyl > eta_slab);

    printf("  [PASS] test_effectiveness_sphere\n");
}

static void test_effectiveness_general(void)
{
    /* Should match sphere for GEOM_SPHERE */
    double eta_gen = cdr_effectiveness_general(1.0, GEOM_SPHERE);
    double eta_sph = cdr_effectiveness_sphere(1.0);
    assert(fabs(eta_gen - eta_sph) < TOL);

    /* Should match slab for GEOM_SLAB */
    eta_gen = cdr_effectiveness_general(1.0, GEOM_SLAB);
    double eta_slab = cdr_effectiveness_slab(1.0);
    assert(fabs(eta_gen - eta_slab) < TOL);

    printf("  [PASS] test_effectiveness_general\n");
}

static void test_effectiveness_nth_order(void)
{
    /* n=1: should match slab */
    double eta_n1 = cdr_effectiveness_nth_order(1.0, 1.0, 1.0, 1.0e-9, 0.001);
    double phi = 0.001 * sqrt(1.0 / 1.0e-9);
    double eta_slab = cdr_effectiveness_slab(phi);
    assert(fabs(eta_n1 - eta_slab) < 1e-6);

    /* n=2: different modulus */
    double eta_n2 = cdr_effectiveness_nth_order(0.1, 2.0, 10.0, 1.0e-5, 0.001);
    assert(eta_n2 > 0.0 && eta_n2 <= 1.0);

    printf("  [PASS] test_effectiveness_nth_order\n");
}

static void test_effectiveness_mm(void)
{
    double eta = cdr_effectiveness_michaelis_menten(1.0, 0.1, 0.1, 1.0e-9, 0.001);
    assert(eta > 0.0 && eta <= 1.0);

    /* Low substrate ? first-order-like */
    double eta_low = cdr_effectiveness_michaelis_menten(1.0, 1.0, 0.01, 1.0e-9, 0.001);
    assert(eta_low > 0.0);

    printf("  [PASS] test_effectiveness_mm\n");
}

static void test_weisz_prater(void)
{
    /* No limitation: C_WP << 1 */
    double C_WP = cdr_weisz_prater_criterion(0.001, 0.001, 1.0e-5, 1.0);
    assert(C_WP < 0.1);

    /* Strong limitation: C_WP >> 1 */
    double C_WP2 = cdr_weisz_prater_criterion(100.0, 0.01, 1.0e-7, 1.0);
    assert(C_WP2 > 1.0);

    assert(cdr_weisz_prater_criterion(1.0, 1.0, 0.0, 1.0) == 0.0);
    printf("  [PASS] test_weisz_prater\n");
}

static void test_mears_criterion(void)
{
    double C_M = cdr_mears_external_criterion(0.1, 0.001, 1.0, 0.01, 1.0);
    assert(C_M > 0.0);

    /* Large k_c ? small C_M (no external limitation) */
    double C_M_fast = cdr_mears_external_criterion(0.1, 0.001, 1.0, 100.0, 1.0);
    assert(C_M_fast < 0.001);

    assert(cdr_mears_external_criterion(1.0, 1.0, 1.0, 0.0, 1.0) == 0.0);
    printf("  [PASS] test_mears_criterion\n");
}

static void test_diagnose_limitations(void)
{
    DiffusionLimitationDiagnosis diag;

    /* Strong internal diffusion limitation */
    cdr_diagnose_limitations(100.0, 0.01, 1.0e-7, 1.0, 1.0, 0.1, 1.0, &diag);
    assert(diag.is_internal_limited == 1);
    assert(diag.Weisz_Prater > 1.0);

    /* No limitation (fast diffusion, small particle) */
    cdr_diagnose_limitations(0.001, 0.0001, 1.0e-3, 1.0, 1.0, 10.0, 1.0, &diag);
    assert(diag.is_internal_limited == 0);

    printf("  [PASS] test_diagnose_limitations\n");
}

static void test_observed_rate(void)
{
    double r_obs = cdr_observed_rate(10.0, 0.5);
    assert(fabs(r_obs - 5.0) < TOL);

    printf("  [PASS] test_observed_rate\n");
}

static void test_intrinsic_activation_energy(void)
{
    /* Under strong diffusion limitation: E_obs = (E_true + E_diff)/2 */
    double E_true = cdr_intrinsic_activation_energy(50.0, 10.0);
    /* E_obs = 50, E_diff = 10 ? E_true = 2*50 - 10 = 90 kJ/mol */
    assert(fabs(E_true - 90.0) < TOL);

    /* Verify: (E_true + E_diff)/2 = (90+10)/2 = 50 = E_obs ? */
    double E_obs_back = (E_true + 10.0) / 2.0;
    assert(fabs(E_obs_back - 50.0) < TOL);

    printf("  [PASS] test_intrinsic_activation_energy\n");
}

static void test_apparent_order(void)
{
    /* n_true = 2 ? n_app = 1.5 */
    double n_app = cdr_apparent_reaction_order(2.0);
    assert(fabs(n_app - 1.5) < TOL);

    /* n_true = 1 ? n_app = 1 (first-order is unchanged!) */
    n_app = cdr_apparent_reaction_order(1.0);
    assert(fabs(n_app - 1.0) < TOL);

    /* n_true = 0 ? n_app = 0.5 */
    n_app = cdr_apparent_reaction_order(0.0);
    assert(fabs(n_app - 0.5) < TOL);

    printf("  [PASS] test_apparent_order\n");
}

static void test_pellet_profile(void)
{
    /* At r=R: concentration = C_surface */
    double C = cdr_pellet_concentration_profile(0.001, 0.001, 1.0, 1.0);
    assert(fabs(C - 1.0) < TOL);

    /* At center (r=0): lowest concentration */
    double C_center = cdr_pellet_concentration_profile(0.0, 0.001, 1.0, 5.0);
    assert(C_center > 0.0 && C_center < 1.0);

    printf("  [PASS] test_pellet_profile\n");
}

static void test_optimal_pellet_radius(void)
{
    double R_opt = cdr_optimal_pellet_radius(1.0, 1.0e-5, 0.5);
    assert(R_opt > 0.0);

    /* High target eta ? small pellet */
    double R_small = cdr_optimal_pellet_radius(1.0, 1.0e-5, 0.95);
    assert(R_small < R_opt);

    printf("  [PASS] test_optimal_pellet_radius\n");
}

static void test_washcoat_effectiveness(void)
{
    /* Thin washcoat ? high effectiveness */
    double eta = cdr_washcoat_effectiveness(1.0e-6, 100.0, 1.0e-8);
    assert(eta > 0.9);

    /* Thick washcoat ? low effectiveness */
    double eta_thick = cdr_washcoat_effectiveness(1.0e-4, 100.0, 1.0e-8);
    assert(eta_thick < 0.5);

    printf("  [PASS] test_washcoat_effectiveness\n");
}

int main(void)
{
    printf("CDR Effectiveness Tests:\n");
    test_effectiveness_slab_asymptotics();
    test_effectiveness_cylinder();
    test_effectiveness_sphere();
    test_effectiveness_general();
    test_effectiveness_nth_order();
    test_effectiveness_mm();
    test_weisz_prater();
    test_mears_criterion();
    test_diagnose_limitations();
    test_observed_rate();
    test_intrinsic_activation_energy();
    test_apparent_order();
    test_pellet_profile();
    test_optimal_pellet_radius();
    test_washcoat_effectiveness();
    printf("All effectiveness tests passed.\n");
    return 0;
}

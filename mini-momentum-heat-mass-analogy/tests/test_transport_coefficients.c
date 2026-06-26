/**
 * @file test_transport_coefficients.c
 * @brief Test suite for transport coefficient computations
 *
 * Tests: Sutherland's law, power-law viscosity, Eucken conductivity,
 * Chapman-Enskog diffusivity, Wilke-Chang correlation, mixture rules,
 * conservation laws, and air/water benchmarks.
 */

#include "transport_coefficients.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>

#define TOL 1e-6
#define TOL_PCT 0.05  /* 5% relative tolerance for engineering correlations */

static int tests_run = 0;
static int tests_failed = 0;

#define TEST(name) printf("  TEST: %s\n", name)

#define CHECK(cond, msg) do { \
    tests_run++; \
    if (!(cond)) { \
        tests_failed++; \
        printf("    FAIL: %s\n", msg); \
    } \
} while(0)

#define ASSERT_NEAR(a, b, tol) CHECK(fabs((a)-(b)) <= (tol), "value mismatch")
#define ASSERT_GT(a, b) CHECK((a) > (b), "expected >")
#define ASSERT_TRUE(cond) CHECK((cond), "expected true")

/* ─── Sutherland's Law Tests ────────────────────────────────────── */

static void test_sutherland_air_viscosity(void)
{
    TEST("Air viscosity at 293K");
    double mu = air_viscosity(293.15);
    /* Expected: ~1.82e-5 Pa·s at 293K */
    ASSERT_NEAR(mu, 1.82e-5, 0.05e-5);
}

static void test_sutherland_null_params(void)
{
    TEST("Sutherland with NULL params");
    double mu = sutherland_viscosity(300.0, NULL);
    ASSERT_NEAR(mu, -1.0, TOL);  /* error return */
}

static void test_sutherland_negative_T(void)
{
    TEST("Sutherland with negative temperature");
    SutherlandParams air = {1.716e-5, 273.15, 110.4};
    double mu = sutherland_viscosity(-100.0, &air);
    ASSERT_NEAR(mu, -1.0, TOL);  /* error return */
}

static void test_sutherland_T_dependence(void)
{
    TEST("Sutherland viscosity increases with T");
    SutherlandParams air = {1.716e-5, 273.15, 110.4};
    double mu_273 = sutherland_viscosity(273.15, &air);
    double mu_373 = sutherland_viscosity(373.15, &air);
    ASSERT_GT(mu_373, mu_273);  /* viscosity increases with temperature for gases */
}

/* ─── Power-Law Viscosity Tests ─────────────────────────────────── */

static void test_power_law_water(void)
{
    TEST("Water viscosity at 293K");
    double mu = water_viscosity(293.15);
    /* Expected: ~1.0e-3 Pa·s */
    ASSERT_NEAR(mu, 1.0e-3, 0.2e-3);
}

static void test_power_law_T_dependence(void)
{
    TEST("Water viscosity decreases with T");
    double mu_273 = water_viscosity(273.15);
    double mu_373 = water_viscosity(373.15);
    ASSERT_GT(mu_273, mu_373);  /* viscosity decreases with temperature for liquids */
}

/* ─── Eucken Conductivity Tests ─────────────────────────────────── */

static void test_eucken_air_conductivity(void)
{
    TEST("Air thermal conductivity at 300K");
    double k = air_thermal_conductivity(300.0);
    /* Expected: ~0.026 W/(m·K) */
    ASSERT_NEAR(k, 0.026, 0.005);
}

static void test_eucken_formula(void)
{
    TEST("Eucken formula with known values");
    /* For N2 at 300K: μ≈1.78e-5, cp≈1041, M=0.028 */
    double mu = 1.78e-5;
    double cp = 1041.0;
    double M = 0.028;
    double k = eucken_thermal_conductivity(mu, cp, M);
    /* Expected: k ≈ μ*(cp + 1.25*R/M) */
    double R = 8.314462618;
    double k_expected = mu * (cp + 1.25 * R / M);
    ASSERT_NEAR(k, k_expected, 1e-6);
}

/* ─── Diffusivity Tests ─────────────────────────────────────────── */

static void test_water_vapor_air_diffusivity(void)
{
    TEST("Water vapor in air diffusivity at 298K, 1atm");
    double D = water_vapor_air_diffusivity(298.0, 101325.0);
    /* Expected: ~2.6e-5 m²/s */
    ASSERT_NEAR(D, 2.6e-5, 0.5e-5);
}

static void test_fuller_diffusivity(void)
{
    TEST("Fuller diffusivity with known inputs (M in g/mol)");
    /* M in g/mol: 18 for water, 29 for air */
    double D = fuller_diffusivity(300.0, 101325.0, 18.0, 29.0, 12.7, 20.1);
    ASSERT_GT(D, 0.0);
    ASSERT_GT(5.0e-4, D);  /* should be positive and reasonable (~2.5e-5 m²/s) */
}

/* ─── Diffusivity Definitions ───────────────────────────────────── */

static void test_thermal_diffusivity(void)
{
    TEST("Thermal diffusivity α = k/(ρ·cp)");
    /* Air: k=0.026, ρ=1.18, cp=1007 → α ≈ 2.2e-5 m²/s */
    double alpha = thermal_diffusivity(0.026, 1.18, 1007.0);
    ASSERT_NEAR(alpha, 2.2e-5, 0.2e-5);
}

static void test_momentum_diffusivity(void)
{
    TEST("Momentum diffusivity ν = μ/ρ");
    /* Air: μ=1.82e-5, ρ=1.18 → ν ≈ 1.54e-5 m²/s */
    double nu = momentum_diffusivity(1.82e-5, 1.18);
    ASSERT_NEAR(nu, 1.54e-5, 0.1e-5);
}

/* ─── Conservation Laws ──────────────────────────────────────────── */

static void test_newtons_law(void)
{
    TEST("Newton's law of viscosity");
    double tau = newtons_law_shear_stress(1.0, 100.0);
    ASSERT_NEAR(tau, 100.0, TOL);  /* τ = μ·|dv/dy| */
}

static void test_fouriers_law(void)
{
    TEST("Fourier's law of conduction");
    double q = fouriers_law_heat_flux(0.6, 1000.0);
    ASSERT_NEAR(q, 600.0, TOL);  /* q = k·|dT/dy| */
}

static void test_ficks_law(void)
{
    TEST("Fick's first law of diffusion");
    double J = ficks_law_mass_flux(2.5e-5, 0.1);
    ASSERT_NEAR(J, 2.5e-6, 1e-12);  /* J = D·|dC/dy| */
}

static void test_generalized_flux(void)
{
    TEST("Generalized transport flux");
    double flux = generalized_transport_flux(0.5, 200.0);
    ASSERT_NEAR(flux, 100.0, TOL);
}

/* ─── Water Benchmark Tests ──────────────────────────────────────── */

static void test_water_conductivity(void)
{
    TEST("Water thermal conductivity at 300K");
    double k = water_thermal_conductivity(300.0);
    /* Expected: ~0.61 W/(m·K) */
    ASSERT_NEAR(k, 0.61, 0.05);
}

/* ─── Mixture Rules ──────────────────────────────────────────────── */

static void test_wilke_interaction_parameter(void)
{
    TEST("Wilke Φ_ij = 1 for identical species");
    double phi = wilke_interaction_parameter(1.0e-5, 1.0e-5, 0.029, 0.029);
    ASSERT_NEAR(phi, 1.0, 0.01);
}

static void test_wilke_mixture_viscosity(void)
{
    TEST("Wilke mixture viscosity single species");
    double x[] = {1.0};
    double mu[] = {1.8e-5};
    double M[] = {0.029};
    double mu_mix = wilke_mixture_viscosity(1, x, mu, M);
    ASSERT_NEAR(mu_mix, 1.8e-5, TOL);
}

/* ─── Transport State Validation ─────────────────────────────────── */

static void test_transport_state_validate(void)
{
    TEST("Validate transport state");
    TransportState state = {1.8e-5, 0.026, 2.5e-5, 1.18, 1007.0, 300.0, 101325.0};
    int err = validate_transport_state(&state);
    CHECK(err == 0, "expected valid state");
}

static void test_transport_state_invalid(void)
{
    TEST("Validate invalid transport state");
    TransportState state = {0.0, 0.026, 2.5e-5, 1.18, 1007.0, 300.0, 101325.0};
    int err = validate_transport_state(&state);
    CHECK(err != 0, "expected error for zero viscosity");
}

static void test_transport_state_null(void)
{
    TEST("Validate NULL transport state");
    int err = validate_transport_state(NULL);
    CHECK(err == -1, "expected -1 for NULL");
}

/* ─── Chapman-Enskog ─────────────────────────────────────────────── */

static void test_chapman_enskog_viscosity(void)
{
    TEST("Chapman-Enskog viscosity for N2 at 300K");
    /* N2: σ=3.681e-10 m, ε/kB=91.5 K, M=4.65e-26 kg/molecule */
    KineticTheoryParams n2 = {3.681e-10, 91.5, 4.65e-26, 1.0, 1.0};
    double mu = chapman_enskog_viscosity(300.0, &n2);
    /* Should be positive and reasonable (~1.8e-5) */
    ASSERT_GT(mu, 0.0);
    ASSERT_GT(1.0e-4, mu);
}

/* ─── Summary ────────────────────────────────────────────────────── */

int main(void)
{
    printf("=== Transport Coefficients Tests ===\n\n");

    /* Sutherland / Gas viscosity */
    test_sutherland_air_viscosity();
    test_sutherland_null_params();
    test_sutherland_negative_T();
    test_sutherland_T_dependence();

    /* Power-law / Liquid viscosity */
    test_power_law_water();
    test_power_law_T_dependence();

    /* Thermal conductivity */
    test_eucken_air_conductivity();
    test_eucken_formula();
    test_water_conductivity();

    /* Diffusivity */
    test_water_vapor_air_diffusivity();
    test_fuller_diffusivity();
    test_thermal_diffusivity();
    test_momentum_diffusivity();

    /* Conservation laws */
    test_newtons_law();
    test_fouriers_law();
    test_ficks_law();
    test_generalized_flux();

    /* Mixture rules */
    test_wilke_interaction_parameter();
    test_wilke_mixture_viscosity();

    /* Transport state */
    test_transport_state_validate();
    test_transport_state_invalid();
    test_transport_state_null();

    /* Chapman-Enskog */
    test_chapman_enskog_viscosity();

    printf("\n=== Results: %d/%d checks passed ===\n",
           tests_run - tests_failed, tests_run);
    return (tests_failed == 0) ? 0 : 1;
}

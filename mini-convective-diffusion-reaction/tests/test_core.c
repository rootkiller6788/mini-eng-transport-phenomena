/**
 * @file test_core.c
 * @brief Tests for CDR core dimensionless numbers and system setup.
 */

#include "../include/cdr_core.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>

#define TOL 1e-10

static void test_peclet_mass(void)
{
    double Pe = cdr_peclet_mass(1.0, 1.0, 1.0);
    assert(fabs(Pe - 1.0) < TOL);

    /* Edge cases */
    assert(cdr_peclet_mass(1.0, 1.0, 0.0) == 0.0);
    assert(cdr_peclet_mass(1.0, 0.0, 1.0) == 0.0);

    /* Pe >> 1: convection-dominated */
    Pe = cdr_peclet_mass(10.0, 1.0, 0.001);
    assert(Pe > 100.0);

    /* Pe << 1: diffusion-dominated */
    Pe = cdr_peclet_mass(0.001, 1.0, 10.0);
    assert(Pe < 0.1);

    printf("  [PASS] test_peclet_mass\n");
}

static void test_peclet_heat(void)
{
    double Pe_h = cdr_peclet_heat(1.0, 1.0, 2.0e-5);
    assert(fabs(Pe_h - 50000.0) < 100.0);

    assert(cdr_peclet_heat(1.0, 1.0, 0.0) == 0.0);
    printf("  [PASS] test_peclet_heat\n");
}

static void test_damkohler_numbers(void)
{
    /* Da_I = k * tau */
    double Da_I = cdr_damkohler_I(0.1, 10.0);
    assert(fabs(Da_I - 1.0) < TOL);

    assert(cdr_damkohler_I(-1.0, 1.0) == 0.0);

    /* Da_II = k * L^2 / D */
    double Da_II = cdr_damkohler_II(1.0, 0.01, 1.0e-5);
    assert(fabs(Da_II - 10.0) < 1.0);

    /* Da_II = phi^2 relationship */
    double phi = cdr_thiele_modulus(0.01, 1.0, 1.0e-5);
    assert(fabs(Da_II - phi * phi) < 1e-8);

    assert(cdr_damkohler_II(1.0, 1.0, 0.0) == 0.0);
    printf("  [PASS] test_damkohler_numbers\n");
}

static void test_thiele_modulus(void)
{
    /* phi = L_c * sqrt(k/D_eff) */
    double phi = cdr_thiele_modulus(0.001, 100.0, 1.0e-5);
    double expected = 0.001 * sqrt(100.0 / 1.0e-5);
    assert(fabs(phi - expected) < TOL);

    assert(cdr_thiele_modulus(1.0, 1.0, 0.0) == 0.0);
    assert(cdr_thiele_modulus(0.0, 1.0, 1.0) == 0.0);

    printf("  [PASS] test_thiele_modulus\n");
}

static void test_schmidt(void)
{
    /* Sc = mu / (rho * D) */
    /* Water at 25C: mu ~ 0.001 Pa.s, rho ~ 1000 kg/m3, D_O2 ~ 2e-9 m2/s */
    double Sc = cdr_schmidt(0.001, 1000.0, 2.0e-9);
    assert(fabs(Sc - 500.0) < 10.0);

    assert(cdr_schmidt(1.0, 0.0, 1.0) == 0.0);
    printf("  [PASS] test_schmidt\n");
}

static void test_reynolds(void)
{
    /* Re = rho * u * L / mu */
    double Re = cdr_reynolds(1000.0, 0.1, 0.01, 0.001);
    assert(fabs(Re - 1000.0) < TOL);

    /* Pipe flow: transition at Re ~ 2100 */
    Re = cdr_reynolds(1000.0, 2.1, 0.01, 0.01);
    assert(Re > 2000.0 && Re < 2200.0);

    assert(cdr_reynolds(1.0, 1.0, 1.0, 0.0) == 0.0);
    printf("  [PASS] test_reynolds\n");
}

static void test_classify_regime(void)
{
    /* Diffusion-dominated: Pe < 0.1, Da_II < 1 */
    assert(cdr_classify_regime(0.01, 0.1, 0.5) == TRANSPORT_DIFFUSION_DOMINATED);

    /* Convection-dominated: Pe > 10, Da_I > 10, Da_II < 1 */
    assert(cdr_classify_regime(100.0, 20.0, 0.5) == TRANSPORT_CONVECTION_DOMINATED);

    /* Reaction-dominated: Da_II > 10, Da_I > 10 */
    assert(cdr_classify_regime(1.0, 20.0, 20.0) == TRANSPORT_REACTION_DOMINATED);

    /* Mass-transfer limited */
    assert(cdr_classify_regime(10.0, 5.0, 0.5) == TRANSPORT_MASS_TRANSFER_LIMITED);

    /* Mixed regime (default) */
    TransportRegime mixed = cdr_classify_regime(1.0, 1.0, 1.0);
    assert(mixed == TRANSPORT_DISPERSION_REGIME || mixed == TRANSPORT_MIXED_REGIME);

    printf("  [PASS] test_classify_regime\n");
}

static void test_characteristic_length(void)
{
    /* Slab: L_c = L (half-thickness) */
    double Lc = cdr_characteristic_length(GEOM_SLAB, 0.01);
    assert(fabs(Lc - 0.01) < TOL);

    /* Cylinder: L_c = R/2 */
    Lc = cdr_characteristic_length(GEOM_CYLINDER, 0.01);
    assert(fabs(Lc - 0.005) < TOL);

    /* Sphere: L_c = R/3 */
    Lc = cdr_characteristic_length(GEOM_SPHERE, 0.01);
    assert(fabs(Lc - 0.01/3.0) < TOL);

    /* Edge case */
    assert(cdr_characteristic_length(GEOM_SPHERE, -1.0) == 0.0);
    assert(cdr_characteristic_length(GEOM_SPHERE, 0.0) == 0.0);

    printf("  [PASS] test_characteristic_length\n");
}

static void test_system_init(void)
{
    CDRSystem sys;
    cdr_system_init(&sys);

    assert(sys.mixture.temperature == 298.15);
    assert(sys.mixture.pressure == 101325.0);
    assert(sys.regime == TRANSPORT_MIXED_REGIME);
    assert(sys.phase == PHASE_HOMOGENEOUS);

    /* NULL safety */
    cdr_system_init(NULL);

    printf("  [PASS] test_system_init\n");
}

static void test_compute_dimensionless(void)
{
    CDRSystem sys;
    cdr_system_init(&sys);
    cdr_compute_dimensionless(&sys);

    /* With default values, Reynolds and heat Peclet should be positive */
    assert(sys.reynolds.Re > 0.0);
    assert(sys.peclet.heat_peclet > 0.0);

    /* NULL safety */
    cdr_compute_dimensionless(NULL);

    printf("  [PASS] test_compute_dimensionless\n");
}

int main(void)
{
    printf("CDR Core Tests:\n");
    test_peclet_mass();
    test_peclet_heat();
    test_damkohler_numbers();
    test_thiele_modulus();
    test_schmidt();
    test_reynolds();
    test_classify_regime();
    test_characteristic_length();
    test_system_init();
    test_compute_dimensionless();
    printf("All core tests passed.\n");
    return 0;
}

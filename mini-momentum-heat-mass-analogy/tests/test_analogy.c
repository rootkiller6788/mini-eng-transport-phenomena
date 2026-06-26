/**
 * @file test_analogy.c
 * @brief Test suite for momentum-heat-mass analogy relationships
 *
 * Tests: Reynolds analogy, Chilton-Colburn j-factors,
 * Prandtl-Taylor, von Karman, boundary layer analogies,
 * dimensionless groups, pipe flow analogies.
 */

#include "momentum_heat_mass_analogy.h"
#include "boundary_layer_analogy.h"
#include "dimensionless_groups.h"
#include "tube_channel_analogy.h"
#include <math.h>
#include <stdio.h>

#define TOL 1e-6

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

/* ─── Reynolds Analogy Tests ────────────────────────────────────── */

static void test_reynolds_analogy_pr1(void)
{
    TEST("Reynolds analogy St = f_F when Pr=1");
    double cf = 0.006;
    double St = stanton_number_from_friction(cf, 1.0, 0);  /* Reynolds form */
    ASSERT_NEAR(St, cf/2.0, TOL);
}

static void test_colburn_form(void)
{
    TEST("Chilton-Colburn St·Pr^(2/3) = f_F");
    double cf = 0.006;
    double Pr = 7.0;  /* water */
    double St = stanton_number_from_friction(cf, Pr, 1);
    double j_H = St * pow(Pr, 2.0/3.0);
    ASSERT_NEAR(j_H, cf/2.0, TOL);
}

static void test_nusselt_from_analogy(void)
{
    TEST("Nu = St·Re·Pr");
    double St = 0.002;
    double Re = 50000.0;
    double Pr = 0.71;
    double Nu = nusselt_from_analogy(St, Re, Pr);
    ASSERT_NEAR(Nu, 0.002 * 50000.0 * 0.71, TOL);
}

static void test_sherwood_from_analogy(void)
{
    TEST("Sh from Chilton-Colburn: Sh = (f/2)·Re·Sc^(1/3)");
    double cf = 0.006;
    double Re = 30000.0;
    double Sc = 0.6;
    double Sh = sherwood_from_analogy(cf, Re, Sc);
    double Sh_expected = (cf/2.0) * Re * pow(Sc, 1.0/3.0);
    ASSERT_NEAR(Sh, Sh_expected, TOL);
}

static void test_friction_from_heat_transfer(void)
{
    TEST("Reverse analogy: f = 4·St·Pr^(2/3)");
    double St = 0.002;
    double Pr = 7.0;
    double f_F = friction_from_heat_transfer(St, Pr);
    ASSERT_NEAR(f_F, St * pow(Pr, 2.0/3.0), TOL);
}

/* ─── Analogy State Tests ───────────────────────────────────────── */

static void test_analogy_state(void)
{
    TEST("Complete analogy state computation");
    AnalogyState state;
    compute_analogy_state(50000.0, 0.71, 0.6, 0.006, &state);
    ASSERT_NEAR(state.Cf, 0.006, TOL);
    ASSERT_GT(state.St_heat, 0.0);
    ASSERT_GT(state.St_mass, 0.0);
}

static void test_colburn_state(void)
{
    TEST("Colburn j-factor computation");
    ChiltonColburnState state;
    state.Re = 50000.0;
    compute_colburn_j_factors(50.0, 0.02, 1.18, 10.0, 1007.0, 0.71, 0.6, &state);
    ASSERT_GT(state.j_H, 0.0);
    ASSERT_GT(state.j_D, 0.0);
}

/* ─── Analogy Closeness ──────────────────────────────────────────── */

static void test_analogy_closeness_metric(void)
{
    TEST("Analogy closeness metric Pr=1, Sc=1 → 1.0");
    double metric = analogy_closeness_metric(1.0, 1.0);
    ASSERT_NEAR(metric, 1.0, 0.01);
}

static void test_analogy_closeness_far(void)
{
    TEST("Analogy closeness Pr>>1, Sc>>1 → near 0");
    double metric = analogy_closeness_metric(1000.0, 1000.0);
    ASSERT_GT(0.5, metric);  /* should be low */
}

/* ─── Reynolds Analogy Error ─────────────────────────────────────── */

static void test_reynolds_analogy_error(void)
{
    TEST("Reynolds analogy error when St = Cf/2 → 0");
    double err = reynolds_analogy_error(0.003, 0.006);
    ASSERT_NEAR(err, 0.0, TOL);
}

/* ─── Engineering Predictions ────────────────────────────────────── */

static void test_predict_h_from_friction(void)
{
    TEST("Predict h from ΔP");
    double h = predict_h_from_friction(100.0, 10.0, 0.0254, 1.18, 10.0, 1007.0, 0.71);
    ASSERT_GT(h, 0.0);
}

static void test_predict_kc_from_h(void)
{
    TEST("Predict k_c from h (Chilton-Colburn)");
    double kc = predict_kc_from_h(50.0, 1.18, 1007.0, 10.0, 0.71, 0.6);
    ASSERT_GT(kc, 0.0);
}

/* ─── Prandtl-Taylor Analogy ─────────────────────────────────────── */

static void test_prandtl_taylor(void)
{
    TEST("Prandtl-Taylor analogy for Pr=1 → St = f_F");
    double f_F = 0.003;
    double St = prandtl_taylor_analogy_St(f_F, 1.0);
    ASSERT_NEAR(St, f_F, TOL);
}

static void test_von_karman(void)
{
    TEST("von Karman analogy for Pr=1 → St = f_F");
    double f_F = 0.003;
    double St = von_karman_analogy_St(f_F, 1.0);
    ASSERT_NEAR(St, f_F, TOL);
}

/* ─── Boundary Layer Tests ──────────────────────────────────────── */

static void test_blasius_eta(void)
{
    TEST("Blasius similarity variable");
    double eta = blasius_eta(0.01, 10.0, 1.5e-5, 1.0);
    ASSERT_GT(eta, 0.0);
}

static void test_laminar_bl_thickness(void)
{
    TEST("Laminar BL thickness δ/x = 5/√Re_x");
    double Re_x = 1.0e5;
    double delta = laminar_bl_thickness(1.0, Re_x);
    double expected = 5.0 / sqrt(Re_x);
    ASSERT_NEAR(delta, expected, 1e-6);
}

static void test_laminar_cf(void)
{
    TEST("Laminar Cf = 0.664/√Re_x");
    double Re_x = 1.0e5;
    double cf = laminar_cf_local(Re_x);
    ASSERT_NEAR(cf, 0.664 / sqrt(Re_x), TOL);
}

static void test_laminar_nusselt(void)
{
    TEST("Laminar Nu_x = 0.332·√Re_x·Pr^(1/3)");
    double Re_x = 1.0e5;
    double Pr = 0.71;
    double Nu = laminar_nusselt_local(Re_x, Pr);
    double expected = 0.332 * sqrt(Re_x) * pow(Pr, 1.0/3.0);
    ASSERT_NEAR(Nu, expected, TOL);
}

static void test_laminar_sherwood(void)
{
    TEST("Laminar Sh_x = 0.332·√Re_x·Sc^(1/3)");
    double Re_x = 1.0e5;
    double Sc = 0.6;
    double Sh = laminar_sherwood_local(Re_x, Sc);
    double expected = 0.332 * sqrt(Re_x) * pow(Sc, 1.0/3.0);
    ASSERT_NEAR(Sh, expected, TOL);
}

/* ─── Dimensionless Groups ───────────────────────────────────────── */

static void test_compute_Re(void)
{
    TEST("Reynolds number computation");
    double Re = compute_Re(1.18, 10.0, 1.0, 1.8e-5);
    ASSERT_GT(Re, 1.0e5);
}

static void test_compute_Pr(void)
{
    TEST("Prandtl number computation");
    double Pr = compute_Pr(1.8e-5, 1007.0, 0.026);
    ASSERT_NEAR(Pr, 0.71, 0.05);
}

static void test_compute_Sc(void)
{
    TEST("Schmidt number computation");
    double Sc = compute_Sc(1.8e-5, 1.18, 2.5e-5);
    ASSERT_GT(Sc, 0.0);
    ASSERT_GT(2.0, Sc);  /* gases: Sc ≈ 0.2-2 */
}

static void test_flow_regime_pipe(void)
{
    TEST("Pipe flow regime: Re=1000 → laminar");
    int regime = flow_regime_pipe(1000.0);
    CHECK(regime == 0, "expected laminar");
}

static void test_flow_regime_pipe_turbulent(void)
{
    TEST("Pipe flow regime: Re=10000 → turbulent");
    int regime = flow_regime_pipe(10000.0);
    CHECK(regime == 2, "expected turbulent");
}

static void test_flow_regime_flat_plate(void)
{
    TEST("Flat plate flow regime: Re_x=1e5 → laminar");
    int regime = flow_regime_flat_plate(1.0e5);
    CHECK(regime == 0, "expected laminar");
}

static void test_j_factor_symmetry(void)
{
    TEST("j_H ≈ j_D for turbulent gas flow");
    PipeFlowAnalogy pipe;
    compute_pipe_flow_analogy(0.0254, 10.0, 5.0,
                              1.18, 1.8e-5, 1007.0, 0.026, 2.5e-5,
                              350.0, 300.0, 0.1, 0.0,
                              &pipe);
    CHECK(pipe.j_H > 0.0, "j_H not positive");
    CHECK(pipe.j_D > 0.0, "j_D not positive");
    double ratio = pipe.j_H / (pipe.j_D + 1e-30);
    CHECK(ratio > 0.1 && ratio < 10.0, "j_H/j_D ratio out of range");
}

/* ─── Tube/Channel Flow Tests ────────────────────────────────────── */

static void test_darcy_friction_laminar(void)
{
    TEST("Darcy friction for laminar: f = 64/Re");
    double f = darcy_friction_laminar(1000.0);
    ASSERT_NEAR(f, 0.064, 1e-6);
}

static void test_blasius_friction_turbulent(void)
{
    TEST("Blasius friction: f = 0.046·Re^(-1/5)");
    double f = blasius_friction_turbulent(50000.0);
    double expected = 0.046 * pow(50000.0, -0.2);
    ASSERT_NEAR(f, expected, 1e-6);
}

static void test_dittus_boelter(void)
{
    TEST("Dittus-Boelter: Nu = 0.023·Re^(4/5)·Pr^0.4");
    double Nu = dittus_boelter_nu(50000.0, 0.71, 1);
    double expected = 0.023 * pow(50000.0, 0.8) * pow(0.71, 0.4);
    ASSERT_NEAR(Nu, expected, 1e-3);
}

static void test_laminar_nu_pipe(void)
{
    TEST("Laminar pipe Nu = 3.66 (const T_w)");
    double Nu = laminar_nu_pipe(0);
    ASSERT_NEAR(Nu, 3.66, TOL);
}

static void test_predict_nu_from_friction(void)
{
    TEST("Predict Nu from friction (Colburn)");
    double f = 0.023;
    double Re = 50000.0;
    double Pr = 0.71;
    double Nu = predict_nu_from_friction(f, Re, Pr);
    ASSERT_GT(Nu, 0.0);
}

static void test_predict_h_from_pressure_drop(void)
{
    TEST("Predict h from pressure drop measurement");
    double h = predict_h_from_pressure_drop(50.0, 0.0254, 10.0, 1.18, 5.0,
                                            0.026, 0.71);
    ASSERT_GT(h, 0.0);
}

/* ─── Summary ────────────────────────────────────────────────────── */

int main(void)
{
    printf("=== Analogy Tests ===\n\n");

    /* Reynolds / Chilton-Colburn */
    test_reynolds_analogy_pr1();
    test_colburn_form();
    test_nusselt_from_analogy();
    test_sherwood_from_analogy();
    test_friction_from_heat_transfer();

    /* Analogy states */
    test_analogy_state();
    test_colburn_state();

    /* Analogy metrics */
    test_analogy_closeness_metric();
    test_analogy_closeness_far();
    test_reynolds_analogy_error();

    /* Engineering predictions */
    test_predict_h_from_friction();
    test_predict_kc_from_h();

    /* Prandtl-Taylor / von Karman */
    test_prandtl_taylor();
    test_von_karman();

    /* Boundary layers */
    test_blasius_eta();
    test_laminar_bl_thickness();
    test_laminar_cf();
    test_laminar_nusselt();
    test_laminar_sherwood();

    /* Dimensionless groups */
    test_compute_Re();
    test_compute_Pr();
    test_compute_Sc();
    test_flow_regime_pipe();
    test_flow_regime_pipe_turbulent();
    test_flow_regime_flat_plate();

    /* Tube/channel flow */
    test_darcy_friction_laminar();
    test_blasius_friction_turbulent();
    test_dittus_boelter();
    test_laminar_nu_pipe();
    test_predict_nu_from_friction();
    test_predict_h_from_pressure_drop();

    /* j-factor symmetry */
    test_j_factor_symmetry();

    printf("\n=== Results: %d/%d checks passed ===\n",
           tests_run - tests_failed, tests_run);
    return (tests_failed == 0) ? 0 : 1;
}

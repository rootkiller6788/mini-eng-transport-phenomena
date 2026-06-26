/*
 * test_dimensionless_numbers.c — Unit tests for core dimensionless group computations
 *
 * Each test verifies a fundamental physical fact encoded in the code.
 * Using standard assert() from <assert.h>, not custom macros.
 */

#include "../include/dimensionless_numbers.h"
#include "../include/reynolds_number.h"
#include "../include/buckingham_pi.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>

#define TOL 1e-10

/* Test: Re for water at 1 m/s in 0.1 m pipe at 20 degC */
void test_reynolds_water_pipe(void)
{
    /* Water at 20 degC: rho=998, mu=0.001 Pa.s */
    double Re = compute_reynolds_number(998.0, 1.0, 0.1, 0.001);
    assert(Re > 0.0);
    assert(fabs(Re - 99800.0) < 1000.0); /* Re ~ 10⁵, turbulent */
    printf("  PASS: test_reynolds_water_pipe — Re ≈ %.1f\n", Re);
}

/* Test: Re for air, low speed — creeping flow */
void test_reynolds_creeping(void)
{
    double Re = compute_reynolds_number(1.2, 0.0001, 0.01, 1.8e-5);
    assert(Re < 1.0);
    assert(Re > 0.0);
    printf("  PASS: test_reynolds_creeping — Re ≈ %.6f (creeping flow)\n", Re);
}

/* Test: Nu basic computation */
void test_nusselt_basic(void)
{
    double Nu = compute_nusselt_number(100.0, 0.05, 0.6);
    assert(Nu > 0.0);
    assert(fabs(Nu - 8.333) < 0.1);
    printf("  PASS: test_nusselt_basic — Nu ≈ %.3f\n", Nu);
}

/* Test: Pr for air */
void test_prandtl_air(void)
{
    double nu = 1.5e-5;   /* m^2/s */
    double alpha = 2.2e-5; /* m^2/s */
    double Pr = compute_prandtl_number(nu, alpha);
    assert(Pr > 0.6 && Pr < 0.8); /* air Pr ≈ 0.7 */
    assert(fabs(Pr - 0.682) < 0.05);
    printf("  PASS: test_prandtl_air — Pr ≈ %.3f\n", Pr);
}

/* Test: Sc for CO2 in air */
void test_schmidt_gas(void)
{
    double nu = 1.5e-5;
    double D = 1.6e-5; /* CO2 in air */
    double Sc = compute_schmidt_number(nu, D);
    assert(Sc > 0.5 && Sc < 2.0);
    printf("  PASS: test_schmidt_gas — Sc ≈ %.3f\n", Sc);
}

/* Test: Gr for natural convection */
void test_grashof_vertical_plate(void)
{
    /* Vertical plate at 60 degC in 20 degC air: beta=1/300, nu=1.6e-5, L=0.3 m */
    double Gr = compute_grashof_number(9.81, 1.0/300.0, 40.0, 0.3, 1.6e-5);
    assert(Gr > 1.0e7); /* Should be large for this DeltaT */
    assert(Gr < 1.0e10);
    printf("  PASS: test_grashof_vertical_plate — Gr ≈ %.2e\n", Gr);
}

/* Test: Ra for water */
void test_rayleigh(void)
{
    double Gr = 1.0e8;
    double Pr = 7.0; /* water */
    double Ra = compute_rayleigh_number(Gr, Pr);
    assert(fabs(Ra - 7.0e8) < 1.0e5);
    printf("  PASS: test_rayleigh — Ra ≈ %.2e\n", Ra);
}

/* Test: We for water droplet */
void test_weber_droplet(void)
{
    double We = compute_weber_number(1000.0, 10.0, 0.001, 0.073);
    assert(We > 0.0);
    assert(We > 1.0); /* Inertia dominates surface tension */
    printf("  PASS: test_weber_droplet — We ≈ %.3f\n", We);
}

/* Test: Ma for subsonic flow */
void test_mach_subsonic(void)
{
    double Ma = compute_mach_number(100.0, 340.0);
    assert(Ma < 1.0);
    assert(fabs(Ma - 0.294) < 0.01);
    printf("  PASS: test_mach_subsonic — Ma ≈ %.3f\n", Ma);
}

/* Test: Fr for ship model */
void test_froude_ship(void)
{
    double Fr = compute_froude_number(5.0, 9.81, 50.0);
    assert(Fr < 1.0);
    printf("  PASS: test_froude_ship — Fr ≈ %.3f\n", Fr);
}

/* Test: Bi for lumped capacitance check */
void test_biot_lumped_capacitance(void)
{
    /* Aluminum sphere in air: h=20, L_c=R/3=0.01/3, k=200 */
    double Bi = compute_biot_number(20.0, 0.01/3.0, 200.0);
    assert(Bi < 0.1); /* Lumped capacitance should be valid */
    printf("  PASS: test_biot_lumped_capacitance — Bi ≈ %.6f (Bi < 0.1)\n", Bi);
}

/* Test: Fo for transient conduction */
void test_fourier_transient(void)
{
    double Fo = compute_fourier_number(1.0e-7, 100.0, 0.01);
    assert(Fo > 0.0);
    assert(fabs(Fo - 0.1) < 0.01);
    printf("  PASS: test_fourier_transient — Fo = %.3f\n", Fo);
}

/* Test: St for heat exchanger analogy */
void test_stanton_heat(void)
{
    double St = compute_stanton_number_heat(100.0, 10000.0, 0.7);
    assert(St > 0.0);
    assert(St < 1.0);
    printf("  PASS: test_stanton_heat — St ≈ %.6f\n", St);
}

/* Test: flow regime classification */
void test_flow_regime_classification(void)
{
    FlowGeometry geom = {0};
    geom.pipe_diameter = 0.05;

    assert(classify_flow_regime(500.0, &geom) == FLOW_LAMINAR);
    assert(classify_flow_regime(3000.0, &geom) == FLOW_TRANSITIONAL);
    assert(classify_flow_regime(10000.0, &geom) == FLOW_TURBULENT_SMOOTH);
    assert(classify_flow_regime(-1.0, &geom) == FLOW_INVALID);
    printf("  PASS: test_flow_regime_classification\n");
}

/* Test: convection type classification */
void test_convection_type(void)
{
    assert(classify_convection_type(0.0, 0.0) == CONV_NONE);
    assert(classify_convection_type(1.0e8, 0.0) == CONV_NATURAL);
    assert(classify_convection_type(0.0, 10000.0) == CONV_FORCED);
    assert(classify_convection_type(1.0e8, 1.0e5) == CONV_MIXED);
    printf("  PASS: test_convection_type\n");
}

/* Test: dynamic similarity check */
void test_dynamic_similarity(void)
{
    assert(check_dynamic_similarity(100000.0, 100100.0, 0.01) == 1);
    assert(check_dynamic_similarity(100000.0, 200000.0, 0.01) == 0);
    assert(check_dynamic_similarity(-1.0, 100.0, 0.01) == 0);
    printf("  PASS: test_dynamic_similarity\n");
}

/* Test: BL thickness ratio — Pr effect */
void test_bl_thickness_ratio(void)
{
    /* Pr = 0.7 (air): delta_T/delta ≈ Pr^(-1/3) ≈ 1.13 */
    double ratio = boundary_layer_thickness_ratio(0.7);
    assert(ratio > 1.0);
    assert(ratio < 1.2);
    /* Pr = 7 (water): delta_T/delta ≈ 7^(-1/3) ≈ 0.52 */
    ratio = boundary_layer_thickness_ratio(7.0);
    assert(ratio < 1.0);
    assert(ratio > 0.4);
    printf("  PASS: test_bl_thickness_ratio\n");
}

/* Test: Re from mass flow rate */
void test_re_pipe_from_mass_flow(void)
{
    /* ṁ=0.5 kg/s, D=0.05 m, mu=0.001 Pa.s -> Re ≈ 4×0.5/(π×0.05×0.001) ≈ 12732 */
    double Re = re_pipe_from_mass_flow(0.5, 0.05, 0.001);
    assert(Re > 10000.0 && Re < 15000.0);
    printf("  PASS: test_re_pipe_from_mass_flow — Re ≈ %.1f\n", Re);
}

/* Test: BL transition location */
void test_bl_transition_location(void)
{
    double x_crit = bl_transition_location(500000.0, 1.5e-5, 30.0);
    assert(x_crit > 0.2 && x_crit < 0.3); /* ~0.25 m for air at 30 m/s */
    printf("  PASS: test_bl_transition_location — x_crit ≈ %.3f m\n", x_crit);
}

/* Test: critical Re database */
void test_critical_reynolds(void)
{
    assert(fabs(critical_reynolds(0) - 2300.0) < 0.01);
    assert(fabs(critical_reynolds(1) - 500000.0) < 0.01);
    assert(fabs(critical_reynolds(2) - 250000.0) < 0.01);
    assert(critical_reynolds(99) < 0.0); /* invalid */
    printf("  PASS: test_critical_reynolds\n");
}

/* Test: N-S non-dimensional scales */
void test_ns_scales(void)
{
    double Re_scale = ns_re_scale(1000.0, 1.0, 0.1, 0.001);
    assert(fabs(Re_scale - 100000.0) < 100.0);
    double p_scale = ns_pressure_scale(1000.0, 1.0);
    assert(fabs(p_scale - 1000.0) < 0.1);
    double t_conv = ns_time_scale(0.1, 1.0);
    assert(fabs(t_conv - 0.1) < 0.01);
    printf("  PASS: test_ns_scales\n");
}

/* Test: dimensional analysis verification */
void test_verify_dimensionless_groups(void)
{
    /* Verify Re is dimensionless */
    const char *vars[] = {"density", "velocity", "length", "viscosity"};
    double exps[] = {1.0, 1.0, 1.0, -1.0};
    assert(verify_dimensionless_group(vars, exps, 4, 1e-10) == 1);

    /* Verify Pr is dimensionless */
    const char *vars2[] = {"kinematic_viscosity", "thermal_diffusivity"};
    double exps2[] = {1.0, -1.0};
    assert(verify_dimensionless_group(vars2, exps2, 2, 1e-10) == 1);

    /* A non-dimensionless combination must fail */
    const char *vars3[] = {"velocity", "length"};
    double exps3[] = {1.0, 1.0};
    assert(verify_dimensionless_group(vars3, exps3, 2, 1e-10) == 0);
    printf("  PASS: test_verify_dimensionless_groups\n");
}

/* Test: dimension vector operations */
void test_dimension_vectors(void)
{
    DimensionVector dv_rho, dv_vel, dv_result;
    double drho[] = {1.0, -3.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double dvel[] = {0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0};

    dimension_vector_create(drho, &dv_rho);
    dimension_vector_create(dvel, &dv_vel);
    dimension_vector_divide(&dv_rho, &dv_vel, &dv_result);

    /* rho/v should have dims [1,-4,1] = mass/(L⁴/T) */
    assert(fabs(dv_result.exponents[0] - 1.0) < TOL);
    assert(fabs(dv_result.exponents[1] - (-4.0)) < TOL);
    assert(fabs(dv_result.exponents[2] - 1.0) < TOL);
    printf("  PASS: test_dimension_vectors\n");
}

/* Test: BL thickness estimates */
void test_bl_thickness_estimates(void)
{
    double x = 0.5;
    double Re_x = 500000.0;
    double delta_lam = bl_thickness_laminar(x, Re_x);
    /* delta_99 ≈ 5.0×0.5/√500000 ≈ 0.0035 m = 3.5 mm */
    assert(delta_lam > 0.002 && delta_lam < 0.005);
    printf("  PASS: test_bl_thickness_estimates — delta_99,lam ≈ %.4f m\n", delta_lam);
}

/* Test: form drag coefficient */
void test_form_drag_coefficient(void)
{
    /* Sphere: Re=100 -> Cd ≈ 0.47 */
    double Cd = form_drag_coefficient(100.0, 0);
    assert(Cd > 0.4 && Cd < 1.5);
    /* Sphere: Re=0.5 -> Cd ≈ 48 */
    Cd = form_drag_coefficient(0.5, 0);
    assert(Cd > 40.0 && Cd < 60.0);
    printf("  PASS: test_form_drag_coefficient\n");
}

int main(void)
{
    printf("=== Test Suite: Dimensionless Numbers Core ===\n");

    test_reynolds_water_pipe();
    test_reynolds_creeping();
    test_nusselt_basic();
    test_prandtl_air();
    test_schmidt_gas();
    test_grashof_vertical_plate();
    test_rayleigh();
    test_weber_droplet();
    test_mach_subsonic();
    test_froude_ship();
    test_biot_lumped_capacitance();
    test_fourier_transient();
    test_stanton_heat();
    test_flow_regime_classification();
    test_convection_type();
    test_dynamic_similarity();
    test_bl_thickness_ratio();
    test_re_pipe_from_mass_flow();
    test_bl_transition_location();
    test_critical_reynolds();
    test_ns_scales();
    test_verify_dimensionless_groups();
    test_dimension_vectors();
    test_bl_thickness_estimates();
    test_form_drag_coefficient();

    printf("\n=== All tests passed ===\n");
    return 0;
}

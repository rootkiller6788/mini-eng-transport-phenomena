/*
 * test_correlations.c — Unit tests for engineering correlations
 *
 * Tests Nu correlations, heat exchanger sizing, electronics cooling,
 * and mass transfer predictions against known benchmark cases.
 */

#include "../include/nusselt_number.h"
#include "../include/prandtl_schmidt.h"
#include "../include/transport_correlations.h"
#include "../include/reynolds_number.h"
#include "../include/buckingham_pi.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>

#define TOL 1e-10

/* Test: Laminar pipe Nu (const T) — known analytical value */
void test_nu_pipe_laminar_const_T(void)
{
    double Nu = nu_pipe_laminar_constant_T();
    assert(fabs(Nu - 3.66) < 0.01);
    printf("  PASS: test_nu_pipe_laminar_const_T — Nu = %.2f\n", Nu);
}

/* Test: Laminar pipe Nu (const q") — known analytical value */
void test_nu_pipe_laminar_const_qflux(void)
{
    double Nu = nu_pipe_laminar_constant_qflux();
    assert(fabs(Nu - 4.36) < 0.01);
    printf("  PASS: test_nu_pipe_laminar_const_qflux — Nu = %.2f\n", Nu);
}

/* Test: Dittus-Boelter turbulent pipe */
void test_nu_dittus_boelter(void)
{
    /* Use Pr=7 (water): Pr^0.4 > Pr^0.3 when Pr > 1 */
    double Nu_h = nu_pipe_dittus_boelter(20000.0, 7.0, 1);
    double Nu_c = nu_pipe_dittus_boelter(20000.0, 7.0, 0);
    assert(Nu_h > 0.0 && Nu_c > 0.0);
    assert(Nu_h > Nu_c); /* For Pr>1: heating exponent 0.4 gives larger Nu than cooling 0.3 */
    printf("  PASS: test_nu_dittus_boelter — Nu_heat=%.2f, Nu_cool=%.2f (Pr=7)\n", Nu_h, Nu_c);
}

/* Test: Gnielinski turbulent pipe — most accurate */
void test_nu_gnielinski(void)
{
    double f = friction_factor_blasius(40000.0);
    double Nu = nu_pipe_gnielinski(40000.0, 0.7, f);
    assert(Nu > 0.0);
    assert(Nu > 50.0 && Nu < 200.0); /* reasonable range */
    printf("  PASS: test_nu_gnielinski — Nu = %.2f (f=%.4f)\n", Nu, f);
}

/* Test: Flat plate laminar average Nu */
void test_nu_flat_plate_laminar(void)
{
    double Re_L = 200000.0;
    double Pr = 0.7;
    double Nu_avg = nu_flat_plate_laminar_avg(Re_L, Pr);
    /* Nu_avg = 0.664*sqrt(2e5)*cbrt(0.7) ≈ 264 */
    assert(Nu_avg > 200.0 && Nu_avg < 350.0);
    double Nu_local = nu_flat_plate_laminar_local(Re_L, Pr);
    assert(fabs(Nu_avg / Nu_local - 2.0) < 0.01); /* avg = 2 × local */
    printf("  PASS: test_nu_flat_plate_laminar — Nu_avg=%.2f, Nu_local=%.2f\n",
           Nu_avg, Nu_local);
}

/* Test: Flat plate mixed boundary layer */
void test_nu_flat_plate_mixed(void)
{
    double Re_L = 5.0e6;
    double Pr = 0.7;
    double Nu = nu_flat_plate_mixed(Re_L, Pr, 500000.0);
    assert(Nu > 1000.0);
    printf("  PASS: test_nu_flat_plate_mixed — Nu = %.2f (Re_L = 5×10⁶)\n", Nu);
}

/* Test: Cylinder cross-flow correlation */
void test_nu_cylinder(void)
{
    double Nu = nu_cylinder_crossflow(50000.0, 0.7);
    assert(Nu > 50.0 && Nu < 200.0);
    printf("  PASS: test_nu_cylinder — Nu = %.2f\n", Nu);
}

/* Test: Sphere correlation (Whitaker) */
void test_nu_sphere(void)
{
    double Nu = nu_sphere(1000.0, 0.7, 1.8e-5, 1.8e-5);
    assert(Nu > 2.0); /* must be greater than pure conduction limit */
    assert(Nu > 10.0 && Nu < 50.0);
    printf("  PASS: test_nu_sphere — Nu = %.2f\n", Nu);
}

/* Test: Nu -> 2 as Re -> 0 for sphere (pure conduction) */
void test_nu_sphere_stagnant(void)
{
    double Nu = nu_sphere(0.0, 0.7, 1.0, 1.0);
    assert(fabs(Nu - 2.0) < 0.01);
    printf("  PASS: test_nu_sphere_stagnant — Nu -> 2.0 (pure conduction)\n");
}

/* Test: Churchill-Chu vertical plate natural convection */
void test_churchill_chu(void)
{
    double Nu = nu_natural_vertical_plate_churchill_chu(1.0e8, 0.7);
    assert(Nu > 10.0 && Nu < 100.0);
    printf("  PASS: test_churchill_chu — Nu = %.2f\n", Nu);
}

/* Test: Colburn j-factor */
void test_colburn_j_factor(void)
{
    double j = colburn_j_factor_heat(100.0, 20000.0, 0.7);
    assert(j > 0.0 && j < 0.01);
    printf("  PASS: test_colburn_j_factor — j_H = %.6f\n", j);
}

/* Test: Prandtl for air vs temperature */
void test_prandtl_air_vs_T(void)
{
    double Pr_300 = prandtl_air(300.0);
    double Pr_1000 = prandtl_air(1000.0);
    assert(Pr_300 > 0.69 && Pr_300 < 0.73);
    assert(Pr_1000 > 0.70 && Pr_1000 < 0.76);
    printf("  PASS: test_prandtl_air_vs_T — Pr(300K)=%.3f, Pr(1000K)=%.3f\n",
           Pr_300, Pr_1000);
}

/* Test: Prandtl for water vs temperature */
void test_prandtl_water_vs_T(void)
{
    double Pr_20 = prandtl_water(20.0);
    double Pr_80 = prandtl_water(80.0);
    assert(Pr_20 > 5.0 && Pr_20 < 10.0);  /* ~7 at 20 degC */
    assert(Pr_80 > 1.5 && Pr_80 < 3.0);   /* ~2 at 80 degC */
    assert(Pr_20 > Pr_80); /* Pr decreases with T for water */
    printf("  PASS: test_prandtl_water_vs_T — Pr(20 degC)=%.2f, Pr(80 degC)=%.2f\n",
           Pr_20, Pr_80);
}

/* Test: Lewis number typical values */
void test_lewis_typical(void)
{
    double Le_gas    = lewis_typical(0);
    double Le_water  = lewis_typical(1);
    assert(Le_gas > 0.5 && Le_gas < 5.0);
    assert(Le_water > 10.0); /* heat diffuses much faster in liquids */
    printf("  PASS: test_lewis_typical — Le_gas=%.2f, Le_water=%.1f\n",
           Le_gas, Le_water);
}

/* Test: Reynolds analogy */
void test_reynolds_analogy(void)
{
    double f = 0.03;
    double St = reynolds_analogy_stanton(f);
    assert(fabs(St - 0.00375) < 0.001);
    double Nu = reynolds_analogy_nusselt(f, 20000.0, 0.7);
    assert(Nu > 30.0);
    printf("  PASS: test_reynolds_analogy — St=%.6f, Nu=%.2f\n", St, Nu);
}

/* Test: Chilton-Colburn Sherwood */
void test_chilton_colburn_sherwood(void)
{
    double Sh = chilton_colburn_sherwood(0.03, 20000.0, 1.0);
    assert(Sh > 50.0);
    printf("  PASS: test_chilton_colburn_sherwood — Sh = %.2f\n", Sh);
}

/* Test: Heat-Mass analogy */
void test_heat_mass_analogy(void)
{
    double Sh = heat_mass_analogy_nusselt_to_sherwood(100.0, 0.7, 1.0);
    assert(Sh > 100.0); /* Sh > Nu when Sc > Pr */
    printf("  PASS: test_heat_mass_analogy — Sh = %.2f\n", Sh);
}

/* Test: Mass transfer flat plate */
void test_mass_transfer_flat_plate(void)
{
    double Sh, h_m;
    int ret = mass_transfer_flat_plate(500000.0, 1.0, 2.0e-5, 1.0,
                                       &Sh, &h_m);
    assert(ret == 0);
    assert(Sh > 100.0);
    assert(h_m > 0.0);
    printf("  PASS: test_mass_transfer_flat_plate — Sh=%.2f, h_m=%.6f m/s\n", Sh, h_m);
}

/* Test: Mass transfer droplet evaporation */
void test_mass_transfer_droplet(void)
{
    double m_dot;
    int ret = mass_transfer_droplet_evaporation(100.0, 0.6, 0.001,
                                                  2.5e-5, 0.02, 0.005,
                                                  &m_dot);
    assert(ret == 0);
    assert(m_dot > 0.0);
    printf("  PASS: test_mass_transfer_droplet — ṁ_evap = %.3e kg/s\n", m_dot);
}

/* Test: Pipe flow design */
void test_pipe_flow_design(void)
{
    double Re, f, delta_P, W;
    int ret = pipe_flow_design(1000.0, 0.001, 2.0, 0.05, 100.0, 0.000045,
                               &Re, &f, &delta_P, &W);
    assert(ret == 0);
    assert(Re > 2300.0); /* Should be turbulent */
    assert(delta_P > 0.0);
    assert(W > 0.0);
    printf("  PASS: test_pipe_flow_design — Re=%.0f, DeltaP=%.0f Pa, W=%.1f W\n",
           Re, delta_P, W);
}

/* Test: Heat exchanger sizing */
void test_heat_exchanger_sizing(void)
{
    double A, U, LMTD;
    int ret = heat_exchanger_sizing(100000.0, 400.0, 350.0, 300.0, 330.0,
                                      2000.0, 500.0, 200.0, 0.002,
                                      &A, &U, &LMTD);
    assert(ret == 0);
    assert(A > 0.0);
    assert(U > 0.0 && U < 500.0); /* should be dominated by lower h */
    printf("  PASS: test_heat_exchanger_sizing — A=%.2f m^2, U=%.1f W/m^2K\n", A, U);
}

/* Test: Friction factor Churchill — covers all regimes */
void test_friction_factor_churchill(void)
{
    double f_lam = friction_factor_churchill(1000.0, 0.0);
    assert(fabs(f_lam - 0.064) < 0.001); /* 64/1000 = 0.064 */

    double f_turb = friction_factor_churchill(50000.0, 0.0);
    assert(f_turb > 0.01 && f_turb < 0.05);
    printf("  PASS: test_friction_factor_churchill — f_lam=%.4f, f_turb=%.4f\n",
           f_lam, f_turb);
}

/* Test: Model scaling velocity */
void test_model_scale_velocity(void)
{
    double U_m = model_scale_velocity_reynolds(10.0, 100.0, 1.5e-5,
                                                1.0, 1.0e-6);
    assert(U_m > 0.0);
    /* water model (nu smaller) -> lower velocity needed */
    assert(U_m < 100.0);
    printf("  PASS: test_model_scale_velocity — U_model = %.2f m/s\n", U_m);
}

/* Test: Kolmogorov microscale */
void test_kolmogorov_microscale(void)
{
    double eta = scaling_law_kolmogorov_microscale(1.0, 1.0e6);
    assert(eta < 1.0);
    assert(eta > 1.0e-6);
    printf("  PASS: test_kolmogorov_microscale — eta = %.3e m\n", eta);
}

/* Test: Diffusion length scale */
void test_diffusion_length(void)
{
    double L = scaling_law_diffusion_length(1.0e-7, 3600.0);
    /* For water: alpha≈1.4e-7, t=1h -> L ≈ 2.2 cm */
    assert(L > 0.01 && L < 0.05);
    printf("  PASS: test_diffusion_length — L_diff = %.4f m\n", L);
}

/* Test: Mars rover night survival */
void test_mars_rover(void)
{
    double T_final, Q_lost;
    int ret = mars_rover_night_survival(8.0, 2.0, 293.0, 50400.0,
                                          50000.0, &T_final, &Q_lost);
    assert(ret == 0);
    assert(T_final > 200.0); /* above Mars night ambient */
    assert(T_final < 320.0);
    printf("  PASS: test_mars_rover — T_final=%.1f K, Q_lost=%.0f J\n",
           T_final, Q_lost);
}

int main(void)
{
    printf("=== Test Suite: Engineering Correlations ===\n");

    test_nu_pipe_laminar_const_T();
    test_nu_pipe_laminar_const_qflux();
    test_nu_dittus_boelter();
    test_nu_gnielinski();
    test_nu_flat_plate_laminar();
    test_nu_flat_plate_mixed();
    test_nu_cylinder();
    test_nu_sphere();
    test_nu_sphere_stagnant();
    test_churchill_chu();
    test_colburn_j_factor();
    test_prandtl_air_vs_T();
    test_prandtl_water_vs_T();
    test_lewis_typical();
    test_reynolds_analogy();
    test_chilton_colburn_sherwood();
    test_heat_mass_analogy();
    test_mass_transfer_flat_plate();
    test_mass_transfer_droplet();
    test_pipe_flow_design();
    test_heat_exchanger_sizing();
    test_friction_factor_churchill();
    test_model_scale_velocity();
    test_kolmogorov_microscale();
    test_diffusion_length();
    test_mars_rover();

    printf("\n=== All tests passed ===\n");
    return 0;
}

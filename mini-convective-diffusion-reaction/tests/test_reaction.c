/**
 * @file test_reaction.c
 * @brief Tests for chemical reaction kinetics.
 */

#include "../include/cdr_reaction.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>

#define TOL 1e-10

static void test_arrhenius_rate(void)
{
    /* k(T) = A * exp(-Ea/(R*T)) */
    double k = cdr_arrhenius_rate(1.0e10, 80000.0, 500.0);
    /* 80000 J/mol at 500K: exp(-80000/(8.314*500)) = exp(-19.24) ~ 4.5e-9 */
    /* k = 1e10 * 4.5e-9 = 45 */
    assert(k > 10.0 && k < 100.0);

    /* At infinite T (large), k -> A */
    double k_hot = cdr_arrhenius_rate(1.0e10, 80000.0, 1.0e10);
    assert(fabs(k_hot - 1.0e10) / 1.0e10 < 0.01);

    assert(cdr_arrhenius_rate(1.0, 100.0, -1.0) == 0.0);
    printf("  [PASS] test_arrhenius_rate\n");
}

static void test_arrhenius_from_ref(void)
{
    double k_ref = 0.1;
    double T_ref = 300.0;
    double k1 = cdr_arrhenius_from_ref(k_ref, 50000.0, T_ref, 350.0);
    /* At higher T, k should increase */
    assert(k1 > k_ref);

    double k2 = cdr_arrhenius_from_ref(k_ref, 50000.0, T_ref, 250.0);
    /* At lower T, k should decrease */
    assert(k2 < k_ref);

    assert(cdr_arrhenius_from_ref(1.0, 100.0, -1.0, 300.0) == 0.0);
    printf("  [PASS] test_arrhenius_from_ref\n");
}

static void test_power_law_rate(void)
{
    /* r = k * C_A^alpha * C_B^beta */
    double r = cdr_power_law_rate(0.1, 2.0, 3.0, 1.0, 1.0);
    assert(fabs(r - 0.6) < TOL);

    /* First order in A only */
    r = cdr_power_law_rate(0.1, 2.0, 1.0, 1.0, 0.0);
    assert(fabs(r - 0.2) < TOL);

    /* Second order in A */
    r = cdr_power_law_rate(0.1, 2.0, 1.0, 2.0, 0.0);
    assert(fabs(r - 0.4) < TOL);

    /* Zero concentration ? zero rate */
    r = cdr_power_law_rate(0.1, 0.0, 1.0, 1.0, 0.0);
    assert(fabs(r) < TOL);

    printf("  [PASS] test_power_law_rate\n");
}

static void test_reversible_rate(void)
{
    PowerLawRate law;
    law.k_forward = 1.0;
    law.k_reverse = 0.1;
    law.alpha = 1.0;
    law.beta = 0.0;
    law.gamma = 1.0;
    law.delta = 0.0;
    law.is_reversible = 1;
    law.Keq = 10.0;

    /* At equilibrium concentrations: C_A * k_f = C_C * k_r ? C_C = 10*C_A */
    double r = cdr_reversible_rate(&law, 1.0, 1.0, 10.0, 1.0);
    assert(fabs(r) < TOL);

    /* NULL safety */
    assert(cdr_reversible_rate(NULL, 1.0, 1.0, 1.0, 1.0) == 0.0);

    printf("  [PASS] test_reversible_rate\n");
}

static void test_michaelis_menten(void)
{
    MichaelisMentenParams mm;
    mm.V_max = 1.0;
    mm.K_m = 1.0;
    mm.k_cat = 10.0;
    mm.E_total = 0.1;

    /* At S = K_m: r = V_max/2 */
    double r = cdr_michaelis_menten_rate(&mm, 1.0);
    assert(fabs(r - 0.5) < TOL);

    /* At S >> K_m: r ? V_max */
    r = cdr_michaelis_menten_rate(&mm, 1000.0);
    assert(fabs(r - 1.0) < 1e-3);

    /* At S << K_m: r ? (V_max/K_m) * S */
    r = cdr_michaelis_menten_rate(&mm, 0.01);
    assert(fabs(r - 0.01) < 0.01);

    /* Competitive inhibition */
    double r_inh = cdr_michaelis_menten_inhibited(&mm, 1.0, 1.0, 0.5);
    assert(r_inh < 0.5);  /* Should be slower than uninhibited */

    assert(cdr_michaelis_menten_rate(NULL, 1.0) == 0.0);
    printf("  [PASS] test_michaelis_menten\n");
}

static void test_langmuir_hinshelwood(void)
{
    LangmuirHinshelwoodParams lh;
    lh.k_surface = 1.0;
    lh.K_ads_A = 1.0;
    lh.K_ads_B = 1.0;
    lh.K_ads_C = 0.0;
    lh.n_sites = 1;
    lh.is_dissociative_A = 0;

    /* Unimolecular at low P: r ? k * K_A * P_A */
    double r = cdr_langmuir_hinshelwood_uni(&lh, 0.01, NULL, NULL, 0);
    assert(r > 0.0);

    /* Bimolecular */
    double r_bi = cdr_langmuir_hinshelwood_bi(&lh, 1.0, 1.0);
    assert(r_bi > 0.0);

    /* NULL safety */
    assert(cdr_langmuir_hinshelwood_uni(NULL, 1.0, NULL, NULL, 0) == 0.0);

    printf("  [PASS] test_langmuir_hinshelwood\n");
}

static void test_half_life(void)
{
    /* First-order t_1/2 = ln(2)/k */
    double t_half = cdr_half_life_first_order(0.1);
    assert(fabs(t_half - log(2.0)/0.1) < TOL);

    assert(isinf(cdr_half_life_first_order(0.0)));
    assert(isinf(cdr_half_life_first_order(-1.0)));

    /* Second-order: t_1/2 = 1/(k*C_A0) */
    double t_half2 = cdr_half_life_nth_order(0.1, 2.0, 1.0);
    assert(fabs(t_half2 - 10.0) < TOL);

    /* nth-order should match first-order when n=1 */
    double t_half_n1 = cdr_half_life_nth_order(0.1, 1.0, 1.0);
    assert(fabs(t_half_n1 - t_half) < TOL);

    printf("  [PASS] test_half_life\n");
}

static void test_conversion(void)
{
    /* First-order: X = 1 - exp(-k*t) */
    double X = cdr_conversion_first_order(0.1, 10.0);
    assert(fabs(X - 0.6321205588) < 1e-6);

    /* Time to reach 90% conversion */
    double t90 = cdr_time_to_conversion_first_order(0.1, 0.9);
    assert(fabs(t90 - 23.02585) < 0.1);

    /* Concentration evolution */
    double C = cdr_concentration_first_order(1.0, 0.1, 10.0);
    assert(fabs(C - 0.367879) < 1e-4);

    printf("  [PASS] test_conversion\n");
}

static void test_concentration_orders(void)
{
    /* Zero-order */
    double C0 = cdr_concentration_zero_order(10.0, 1.0, 5.0);
    assert(fabs(C0 - 5.0) < TOL);

    double C0_depleted = cdr_concentration_zero_order(1.0, 1.0, 2.0);
    assert(fabs(C0_depleted) < TOL);

    /* Second-order */
    double C2 = cdr_concentration_second_order(1.0, 0.1, 10.0);
    assert(fabs(C2 - 0.5) < TOL);

    printf("  [PASS] test_concentration_orders\n");
}

static void test_series_reaction(void)
{
    double C_out[3];
    cdr_series_reaction_concentrations(1.0, 0.5, 0.1, 1.0, C_out);

    /* Conservation: C_A + C_B + C_C = C_A0 */
    double sum = C_out[0] + C_out[1] + C_out[2];
    assert(fabs(sum - 1.0) < TOL);

    /* C_A should decrease */
    assert(C_out[0] < 1.0);
    assert(C_out[0] > 0.0);

    /* Test degenerate case k1 == k2 */
    double C_out2[3];
    cdr_series_reaction_concentrations(1.0, 0.1, 0.1, 1.0, C_out2);
    double sum2 = C_out2[0] + C_out2[1] + C_out2[2];
    assert(fabs(sum2 - 1.0) < TOL);

    printf("  [PASS] test_series_reaction\n");
}

static void test_activation_energy_estimation(void)
{
    double k1 = 0.01, T1 = 300.0;
    double k2 = 0.1, T2 = 350.0;
    double Ea = cdr_activation_energy_from_two_points(k1, T1, k2, T2);
    assert(Ea > 0.0);

    /* Check self-consistency */
    double k2_back = cdr_arrhenius_from_ref(k1, Ea, T1, T2);
    assert(fabs(k2_back - k2) / k2 < 0.01);

    printf("  [PASS] test_activation_energy_estimation\n");
}

static void test_selectivity(void)
{
    /* Parallel: if n1 > n2, higher C_A favors B */
    double S = cdr_selectivity_parallel(1.0, 1.0, 2.0, 1.0, 10.0);
    assert(S > 1.0);

    /* If n1 < n2, lower C_A favors B */
    double S2 = cdr_selectivity_parallel(1.0, 1.0, 1.0, 2.0, 10.0);
    assert(S2 < 1.0);

    assert(isinf(cdr_selectivity_parallel(1.0, 0.0, 1.0, 1.0, 1.0)));
    printf("  [PASS] test_selectivity\n");
}

static void test_yield_series(void)
{
    /* A -> B -> C, k1 > k2: B accumulates */
    double t_opt = cdr_optimal_time_series(1.0, 0.1);
    assert(t_opt > 0.0);

    double Y_opt = cdr_yield_intermediate_series(1.0, 0.1, 1.0, t_opt);
    assert(Y_opt > 0.0 && Y_opt < 1.0);

    printf("  [PASS] test_yield_series\n");
}

int main(void)
{
    printf("CDR Reaction Tests:\n");
    test_arrhenius_rate();
    test_arrhenius_from_ref();
    test_power_law_rate();
    test_reversible_rate();
    test_michaelis_menten();
    test_langmuir_hinshelwood();
    test_half_life();
    test_conversion();
    test_concentration_orders();
    test_series_reaction();
    test_activation_energy_estimation();
    test_selectivity();
    test_yield_series();
    printf("All reaction tests passed.\n");
    return 0;
}

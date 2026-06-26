/**
 * @file bench_cdr.c
 * @brief Performance benchmarks for CDR computations.
 *
 * Benchmarks key computational kernels: Thiele modulus, effectiveness
 * factor, RTD evaluation, and concentration profile calculations.
 */

#include "../include/cdr_core.h"
#include "../include/cdr_effectiveness.h"
#include "../include/cdr_reactor.h"
#include "../include/cdr_reaction.h"
#include <stdio.h>
#include <time.h>

#define ITERATIONS 1000000

static double get_time_sec(void)
{
    return (double)clock() / CLOCKS_PER_SEC;
}

static void bench_thiele_modulus(void)
{
    double start = get_time_sec();
    double sum = 0.0;
    for (int i = 0; i < ITERATIONS; i++) {
        double Lc = 0.001 + i * 1e-9;
        sum += cdr_thiele_modulus(Lc, 100.0, 1.0e-5);
    }
    double elapsed = get_time_sec() - start;
    printf("  Thiele modulus: %.3f us/call (sum=%.6f)\n",
           elapsed / ITERATIONS * 1e6, sum);
}

static void bench_effectiveness_sphere(void)
{
    double start = get_time_sec();
    double sum = 0.0;
    for (int i = 0; i < ITERATIONS; i++) {
        double phi = 0.01 + (double)i / ITERATIONS * 50.0;
        sum += cdr_effectiveness_sphere(phi);
    }
    double elapsed = get_time_sec() - start;
    printf("  Effectiveness (sphere): %.3f us/call (sum=%.6f)\n",
           elapsed / ITERATIONS * 1e6, sum);
}

static void bench_rtd_cstr(void)
{
    double start = get_time_sec();
    double sum = 0.0;
    double tau = 10.0;
    for (int i = 0; i < ITERATIONS; i++) {
        double t = (double)i / ITERATIONS * 50.0;
        sum += cdr_rtd_cstr_E(t, tau);
    }
    double elapsed = get_time_sec() - start;
    printf("  RTD CSTR: %.3f us/call (sum=%.6f)\n",
           elapsed / ITERATIONS * 1e6, sum);
}

static void bench_arrhenius(void)
{
    double start = get_time_sec();
    double sum = 0.0;
    for (int i = 0; i < ITERATIONS; i++) {
        double T = 300.0 + (double)i / ITERATIONS * 500.0;
        sum += cdr_arrhenius_rate(1.0e10, 80000.0, T);
    }
    double elapsed = get_time_sec() - start;
    printf("  Arrhenius: %.3f us/call (sum=%.6f)\n",
           elapsed / ITERATIONS * 1e6, sum);
}

static void bench_pfr_profile(void)
{
    double start = get_time_sec();
    double sum = 0.0;
    for (int i = 0; i < ITERATIONS; i++) {
        double z = (double)i / ITERATIONS;
        sum += cdr_pfr_profile_first_order(1.0, 0.1, 100.0, z);
    }
    double elapsed = get_time_sec() - start;
    printf("  PFR profile: %.3f us/call (sum=%.6f)\n",
           elapsed / ITERATIONS * 1e6, sum);
}

int main(void)
{
    printf("CDR Benchmarks (%d iterations each):\n\n", ITERATIONS);
    bench_thiele_modulus();
    bench_effectiveness_sphere();
    bench_rtd_cstr();
    bench_arrhenius();
    bench_pfr_profile();
    printf("\nBenchmarks complete.\n");
    return 0;
}

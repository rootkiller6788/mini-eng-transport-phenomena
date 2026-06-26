/**
 * @file boundary_layer_analogy.c
 * @brief Boundary layer analogy implementations
 *
 * Implements Blasius velocity profile, Pohlhausen thermal profile,
 * and related boundary layer thickness correlations that demonstrate
 * the analogy between momentum, thermal, and concentration boundary layers.
 *
 * Reference: Blasius (1908), Z. Math. Phys., 56, 1-37.
 *            Pohlhausen (1921), Z. Angew. Math. Mech., 1, 115-121.
 *            Schlichting & Gersten (2000), Boundary Layer Theory.
 */

#include "boundary_layer_analogy.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/* ─── Blasius Similarity Solution ───────────────────────────────── */

/**
 * The Blasius equation: f''' + (1/2)·f·f'' = 0
 * BC: f(0) = 0, f'(0) = 0, f'(∞) = 1
 *
 * We solve this as a system of first-order ODEs:
 *   y0 = f    → y0' = y1
 *   y1 = f'   → y1' = y2
 *   y2 = f''  → y2' = -0.5 · y0 · y2
 *
 * Using RK4 with step size h = 0.025 from η=0 to η=10.
 * The shooting method adjusts f''(0) until f'(∞) ≈ 1.
 * The correct value is f''(0) ≈ 0.332057...
 */

#define BLASIUS_FPP0 0.332057336215195  /* f''(0) for Blasius solution */
#define BLASIUS_ETA_MAX 10.0
#define BLASIUS_N_STEPS 2000

static void rk4_step_blasius(double *y, double h)
{
    /* y = [f, f', f''] */
    double k1[3], k2[3], k3[3], k4[3], y_temp[3];

    /* k1 */
    k1[0] = y[1];
    k1[1] = y[2];
    k1[2] = -0.5 * y[0] * y[2];

    /* k2 */
    y_temp[0] = y[0] + 0.5 * h * k1[0];
    y_temp[1] = y[1] + 0.5 * h * k1[1];
    y_temp[2] = y[2] + 0.5 * h * k1[2];
    k2[0] = y_temp[1];
    k2[1] = y_temp[2];
    k2[2] = -0.5 * y_temp[0] * y_temp[2];

    /* k3 */
    y_temp[0] = y[0] + 0.5 * h * k2[0];
    y_temp[1] = y[1] + 0.5 * h * k2[1];
    y_temp[2] = y[2] + 0.5 * h * k2[2];
    k3[0] = y_temp[1];
    k3[1] = y_temp[2];
    k3[2] = -0.5 * y_temp[0] * y_temp[2];

    /* k4 */
    y_temp[0] = y[0] + h * k3[0];
    y_temp[1] = y[1] + h * k3[1];
    y_temp[2] = y[2] + h * k3[2];
    k4[0] = y_temp[1];
    k4[1] = y_temp[2];
    k4[2] = -0.5 * y_temp[0] * y_temp[2];

    /* Update */
    y[0] += (h / 6.0) * (k1[0] + 2.0*k2[0] + 2.0*k3[0] + k4[0]);
    y[1] += (h / 6.0) * (k1[1] + 2.0*k2[1] + 2.0*k3[1] + k4[1]);
    y[2] += (h / 6.0) * (k1[2] + 2.0*k2[2] + 2.0*k3[2] + k4[2]);
}

/**
 * Generate Blasius profile at specified η values.
 *
 * We integrate the full profile once, then interpolate to the
 * requested η values using piecewise linear interpolation.
 */
void blasius_profile(const double *eta, size_t n, BlasiusProfilePoint *profile)
{
    if (!eta || !profile || n == 0) return;

    /* Generate reference solution on a fine grid */
    size_t n_ref = BLASIUS_N_STEPS;
    double *eta_ref = (double *)malloc(n_ref * sizeof(double));
    double *f_ref   = (double *)malloc(n_ref * sizeof(double));
    double *fp_ref  = (double *)malloc(n_ref * sizeof(double));
    double *fpp_ref = (double *)malloc(n_ref * sizeof(double));

    if (!eta_ref || !f_ref || !fp_ref || !fpp_ref) {
        free(eta_ref); free(f_ref); free(fp_ref); free(fpp_ref);
        return;
    }

    double h = BLASIUS_ETA_MAX / (n_ref - 1.0);
    double y[3] = {0.0, 0.0, BLASIUS_FPP0};

    for (size_t i = 0; i < n_ref; i++) {
        eta_ref[i] = i * h;
        f_ref[i] = y[0];
        fp_ref[i] = y[1];
        fpp_ref[i] = y[2];
        rk4_step_blasius(y, h);
    }

    /* Interpolate to requested eta values */
    for (size_t k = 0; k < n; k++) {
        double e = eta[k];
        if (e < 0.0) e = 0.0;
        if (e > BLASIUS_ETA_MAX) e = BLASIUS_ETA_MAX;

        /* Find bracket */
        size_t idx = (size_t)(e / h);
        if (idx >= n_ref - 1) idx = n_ref - 2;

        /* Linear interpolation */
        double t = (e - eta_ref[idx]) / h;
        profile[k].eta  = e;
        profile[k].f    = f_ref[idx]   + t * (f_ref[idx+1]   - f_ref[idx]);
        profile[k].fp   = fp_ref[idx]  + t * (fp_ref[idx+1]  - fp_ref[idx]);
        profile[k].fpp  = fpp_ref[idx] + t * (fpp_ref[idx+1] - fpp_ref[idx]);
    }

    free(eta_ref); free(f_ref); free(fp_ref); free(fpp_ref);
}

double blasius_eta(double y, double U_inf, double nu, double x)
{
    if (x <= 0.0 || nu <= 0.0) return -1.0;
    return y * sqrt(U_inf / (nu * x));
}

/* ─── Velocity Boundary Layer ──────────────────────────────────── */

double laminar_bl_thickness(double x, double Re_x)
{
    if (Re_x <= 0.0) return -1.0;
    return 5.0 * x / sqrt(Re_x);
}

double turbulent_bl_thickness(double x, double Re_x)
{
    if (Re_x <= 0.0) return -1.0;
    return 0.37 * x / pow(Re_x, 0.2);
}

double laminar_cf_local(double Re_x)
{
    if (Re_x <= 0.0) return -1.0;
    return 0.664 / sqrt(Re_x);
}

double turbulent_cf_local(double Re_x)
{
    if (Re_x <= 0.0) return -1.0;
    return 0.0592 / pow(Re_x, 0.2);
}

double laminar_momentum_thickness(double x, double Re_x)
{
    if (Re_x <= 0.0) return -1.0;
    return 0.664 * x / sqrt(Re_x);
}

double laminar_displacement_thickness(double x, double Re_x)
{
    if (Re_x <= 0.0) return -1.0;
    return 1.7208 * x / sqrt(Re_x);
}

/* ─── Thermal Boundary Layer ───────────────────────────────────── */

/**
 * Pohlhausen thermal boundary layer solution.
 *
 * For laminar flow over an isothermal flat plate:
 * Θ'' + (Pr/2)·f(η)·Θ' = 0
 * BC: Θ(0) = 0, Θ(∞) = 1
 *
 * This is solved by integration using the known Blasius f(η):
 * Θ'(η) = Θ'(0)·exp(-(Pr/2)·∫₀^η f(s)·ds)
 * Θ(η)  = Θ'(0)·∫₀^η exp(-(Pr/2)·∫₀^t f(s)·ds)·dt
 *
 * where Θ'(0) = 1/∫₀^∞ exp(-(Pr/2)·∫₀^η f(s)·ds)·dη
 *
 * Complexity: O(n_steps)
 */
void pohlhausen_thermal_profile(const double *eta, const BlasiusProfilePoint *blasius,
                                double Pr, size_t n,
                                PohlhausenProfilePoint *thermal_profile)
{
    if (!eta || !blasius || !thermal_profile || n == 0) return;

    /* For Pr = 1, Θ(η) = f'(η) exactly (Reynolds analogy) */
    if (fabs(Pr - 1.0) < 1e-6) {
        for (size_t i = 0; i < n; i++) {
            thermal_profile[i].eta = eta[i];
            thermal_profile[i].Theta = blasius[i].fp;
            thermal_profile[i].Theta_p = blasius[i].fpp;
            thermal_profile[i].Pr = Pr;
        }
        return;
    }

    /* General case: numerical integration using trapezoidal rule */
    /* Step 1: compute ∫₀^η f(s)·ds (trapezoidal on the given points) */
    double *integral_f = (double *)malloc(n * sizeof(double));
    if (!integral_f) return;

    integral_f[0] = 0.0;
    for (size_t i = 1; i < n; i++) {
        double ds = eta[i] - eta[i-1];
        integral_f[i] = integral_f[i-1]
                      + 0.5 * (blasius[i].f + blasius[i-1].f) * ds;
    }

    /* Step 2: compute exp(-(Pr/2)·∫f·ds) at each point */
    double *exp_factor = (double *)malloc(n * sizeof(double));
    if (!exp_factor) { free(integral_f); return; }

    for (size_t i = 0; i < n; i++) {
        exp_factor[i] = exp(-0.5 * Pr * integral_f[i]);
    }

    /* Step 3: compute ∫₀^η exp_factor ds for the denominator integral */
    double *integral_exp = (double *)malloc(n * sizeof(double));
    if (!integral_exp) { free(integral_f); free(exp_factor); return; }

    integral_exp[0] = 0.0;
    for (size_t i = 1; i < n; i++) {
        double ds = eta[i] - eta[i-1];
        integral_exp[i] = integral_exp[i-1]
                        + 0.5 * (exp_factor[i] + exp_factor[i-1]) * ds;
    }

    /* Θ'(0) = 1 / ∫₀^∞ exp_factor·dη */
    double Theta_p_0 = 1.0 / integral_exp[n-1];

    /* Step 4: compute Θ(η) = Θ'(0) · integral_exp[η] */
    for (size_t i = 0; i < n; i++) {
        thermal_profile[i].eta = eta[i];
        thermal_profile[i].Theta = Theta_p_0 * integral_exp[i];
        thermal_profile[i].Theta_p = Theta_p_0 * exp_factor[i];
        thermal_profile[i].Pr = Pr;
    }

    free(integral_f);
    free(exp_factor);
    free(integral_exp);
}

double thermal_bl_thickness_laminar(double delta_velocity, double Pr)
{
    if (Pr <= 0.0) return -1.0;
    return delta_velocity * 0.975 * pow(Pr, -1.0 / 3.0);
}

double laminar_nusselt_local(double Re_x, double Pr)
{
    if (Re_x <= 0.0) return -1.0;
    return 0.332 * sqrt(Re_x) * pow(Pr, 1.0 / 3.0);
}

double laminar_nusselt_average(double Re_L, double Pr)
{
    if (Re_L <= 0.0) return -1.0;
    return 0.664 * sqrt(Re_L) * pow(Pr, 1.0 / 3.0);
}

double turbulent_nusselt_local(double Re_x, double Pr)
{
    if (Re_x <= 0.0) return -1.0;
    return 0.0296 * pow(Re_x, 0.8) * pow(Pr, 1.0 / 3.0);
}

/* ─── Concentration Boundary Layer ──────────────────────────────── */

double laminar_sherwood_local(double Re_x, double Sc)
{
    if (Re_x <= 0.0) return -1.0;
    return 0.332 * sqrt(Re_x) * pow(Sc, 1.0 / 3.0);
}

double laminar_sherwood_average(double Re_L, double Sc)
{
    if (Re_L <= 0.0) return -1.0;
    return 0.664 * sqrt(Re_L) * pow(Sc, 1.0 / 3.0);
}

double turbulent_sherwood_local(double Re_x, double Sc)
{
    if (Re_x <= 0.0) return -1.0;
    return 0.0296 * pow(Re_x, 0.8) * pow(Sc, 1.0 / 3.0);
}

double concentration_bl_thickness(double delta_velocity, double Sc)
{
    if (Sc <= 0.0) return -1.0;
    return delta_velocity * 0.975 * pow(Sc, -1.0 / 3.0);
}

/* ─── L4: Analogy Verification ─────────────────────────────────── */

double profile_correlation(const BlasiusProfilePoint *velocity,
                           const PohlhausenProfilePoint *thermal,
                           size_t n)
{
    if (!velocity || !thermal || n < 2) return -1.0;

    /* Pearson r = Σ((x_i-x̄)(y_i-ȳ)) / sqrt(Σ(x_i-x̄)²·Σ(y_i-ȳ)²) */
    double sum_v = 0.0, sum_t = 0.0;
    for (size_t i = 0; i < n; i++) {
        sum_v += velocity[i].fp;
        sum_t += thermal[i].Theta;
    }
    double mean_v = sum_v / n;
    double mean_t = sum_t / n;

    double cov = 0.0, var_v = 0.0, var_t = 0.0;
    for (size_t i = 0; i < n; i++) {
        double dv = velocity[i].fp - mean_v;
        double dt = thermal[i].Theta - mean_t;
        cov  += dv * dt;
        var_v += dv * dv;
        var_t += dt * dt;
    }

    if (var_v < 1e-30 || var_t < 1e-30) return 0.0;
    return cov / sqrt(var_v * var_t);
}

/* ─── L5: Complete Boundary Layer Summary ──────────────────────── */

void compute_boundary_layer_summary(const FlatPlateState *state,
                                    BoundaryLayerSummary *summary)
{
    if (!state || !summary) return;

    summary->x = state->x;
    summary->flow = *state;

    /* Velocity BL */
    if (state->flow_regime == 0) {
        summary->bl.delta = laminar_bl_thickness(state->x, state->Re_x);
        summary->bl.delta_m = laminar_momentum_thickness(state->x, state->Re_x);
        summary->bl.delta_d = laminar_displacement_thickness(state->x, state->Re_x);
        summary->Cf_x = laminar_cf_local(state->Re_x);
    } else {
        summary->bl.delta = turbulent_bl_thickness(state->x, state->Re_x);
        summary->bl.delta_m = 0.036 * state->x
                              * pow(state->Re_x, -0.2);
        summary->bl.delta_d = summary->bl.delta * 0.125;  /* approx 1/8 */
        summary->Cf_x = turbulent_cf_local(state->Re_x);
    }

    /* Thermal BL */
    if (state->flow_regime == 0) {
        summary->bl.delta_T = thermal_bl_thickness_laminar(
                                  summary->bl.delta, state->Pr);
        summary->bl.delta_h = laminar_nusselt_local(state->Re_x, state->Pr)
                              * state->alpha / state->U_inf;
        summary->Nu_x = laminar_nusselt_local(state->Re_x, state->Pr);
    } else {
        summary->bl.delta_T = summary->bl.delta * pow(state->Pr, -1.0/3.0);
        summary->Nu_x = turbulent_nusselt_local(state->Re_x, state->Pr);
    }

    /* Concentration BL */
    if (state->flow_regime == 0) {
        summary->bl.delta_C = concentration_bl_thickness(
                                  summary->bl.delta, state->Sc);
        summary->bl.delta_c = laminar_sherwood_local(state->Re_x, state->Sc)
                              * state->D_AB / state->U_inf;
        summary->Sh_x = laminar_sherwood_local(state->Re_x, state->Sc);
    } else {
        summary->bl.delta_C = summary->bl.delta * pow(state->Sc, -1.0/3.0);
        summary->Sh_x = turbulent_sherwood_local(state->Re_x, state->Sc);
    }

    summary->bl.Re_x = state->Re_x;
    summary->bl.Pr = state->Pr;
    summary->bl.Sc = state->Sc;

    /* Colburn j-factors */
    double St_h = summary->Nu_x / (state->Re_x * state->Pr);
    double St_m = summary->Sh_x / (state->Re_x * state->Sc);
    summary->j_H_x = St_h * pow(state->Pr, 2.0 / 3.0);
    summary->j_D_x = St_m * pow(state->Sc, 2.0 / 3.0);
}

int boundary_layer_ordering(const BoundaryLayerThicknesses *bl)
{
    if (!bl) return 0;
    /* Compare relative thicknesses */
    if (bl->delta_T > bl->delta && bl->delta_C > bl->delta_T)
        return -1;  /* δ_C > δ_T > δ (Sc < Pr < 1) */
    if (bl->delta > bl->delta_T && bl->delta_T > bl->delta_C)
        return +1;  /* δ > δ_T > δ_C (1 < Pr < Sc) */
    return 0;        /* all similar */
}

void print_boundary_layer_comparison(const BoundaryLayerSummary *summary)
{
    if (!summary) return;
    const BoundaryLayerThicknesses *bl = &summary->bl;

    printf("━━ Boundary Layer Comparison at x=%.3f m ━━━━━━━━━━━━━━━━\n",
           summary->x);
    printf("  Re_x = %.1e    Pr = %.3f    Sc = %.3f\n",
           bl->Re_x, bl->Pr, bl->Sc);
    printf("  ─── Thicknesses ────────────────────────────────\n");
    printf("  Velocity BL δ_99:     %.4f mm\n",
           bl->delta * 1000.0);
    printf("  Thermal BL  δ_T:      %.4f mm  (δ_T/δ = %.3f)\n",
           bl->delta_T * 1000.0, bl->delta_T / (bl->delta + 1e-30));
    printf("  Concentration BL δ_C: %.4f mm  (δ_C/δ = %.3f)\n",
           bl->delta_C * 1000.0, bl->delta_C / (bl->delta + 1e-30));
    printf("  Displacement δ*:      %.4f mm\n",
           bl->delta_d * 1000.0);
    printf("  Momentum θ:           %.4f mm\n",
           bl->delta_m * 1000.0);
    printf("  ─── Transfer Coefficients ──────────────────────\n");
    printf("  Cf_x  = %.6f\n", summary->Cf_x);
    printf("  Nu_x  = %.3f\n", summary->Nu_x);
    printf("  Sh_x  = %.3f\n", summary->Sh_x);
    printf("  j_H   = %.6f\n", summary->j_H_x);
    printf("  j_D   = %.6f\n", summary->j_D_x);
    printf("  j_H ≈ j_D: %s\n",
           fabs(summary->j_H_x - summary->j_D_x) < 0.1 * summary->j_H_x
           ? "YES ✓" : "No");
    int order = boundary_layer_ordering(bl);
    printf("  BL ordering: ");
    if (order < 0) printf("δ_C > δ_T > δ  (Sc < Pr < 1, liquid metal)\n");
    else if (order > 0) printf("δ > δ_T > δ_C  (Pr > 1, Sc > 1, liquid)\n");
    else printf("δ ≈ δ_T ≈ δ_C  (Pr ≈ Sc ≈ 1, gas)\n");
    printf("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");
}

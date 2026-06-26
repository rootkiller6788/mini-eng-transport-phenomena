/**
 * @file momentum_heat_mass_analogy.c
 * @brief Core analogy implementations: Reynolds, Chilton-Colburn, Prandtl-Taylor
 *
 * This is the central module of the transport analogy framework.
 * It implements the mathematical relationships connecting momentum,
 * heat, and mass transfer through dimensionless analogies.
 *
 * Reference: Bird, Stewart, Lightfoot (2007), Transport Phenomena, Ch. 13.
 *            Colburn (1933), "A method of correlating forced convection
 *            heat-transfer data and a comparison with fluid friction,"
 *            Trans. AIChE, 29, 174-210.
 */

#include "momentum_heat_mass_analogy.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

/* ─── L1: Stanton Number from Friction ──────────────────────────── */

double stanton_number_from_friction(double cf, double Pr, int colburn_form)
{
    /* cf = skin friction coefficient = 2·τ_w/(ρ·U∞²) */
    /* f_F = cf/2 = Fanning friction factor */
    double f_F = cf / 2.0;

    if (colburn_form) {
        /* Chilton-Colburn: St·Pr^(2/3) = f_F */
        return f_F * pow(Pr, -2.0 / 3.0);
    } else {
        /* Reynolds: St = f_F  (only valid when Pr ≈ 1) */
        return f_F;
    }
}

/* ─── L1: Nusselt from Analogy ──────────────────────────────────── */

double nusselt_from_analogy(double St, double Re, double Pr)
{
    return St * Re * Pr;
}

/* ─── L1: Sherwood from Analogy ─────────────────────────────────── */

double sherwood_from_analogy(double cf, double Re, double Sc)
{
    double f_F = cf / 2.0;
    /* Chilton-Colburn mass transfer: Sh/(Re·Sc^(1/3)) = f_F */
    return f_F * Re * pow(Sc, 1.0 / 3.0);
}

/* ─── L1: Friction from Heat Transfer ───────────────────────────── */

double friction_from_heat_transfer(double St_heat, double Pr)
{
    /* Reverse Chilton-Colburn: f_F = St·Pr^(2/3) */
    return St_heat * pow(Pr, 2.0 / 3.0);
}

/* ─── L2: Complete Analogy State ────────────────────────────────── */

void compute_analogy_state(double Re, double Pr, double Sc, double cf,
                           AnalogyState *analogy)
{
    if (!analogy) return;
    analogy->Re = Re;
    analogy->Pr = Pr;
    analogy->Sc = Sc;
    analogy->Cf = cf;

    /* Stanton numbers */
    analogy->St_heat = stanton_number_from_friction(cf, Pr, 1);
    analogy->St_mass = stanton_number_from_friction(cf, Sc, 1);

    /* Analogy type: 0=Reynolds, 1=Chilton-Colburn */
    if (fabs(Pr - 1.0) < 0.1 && fabs(Sc - 1.0) < 0.1) {
        analogy->analogy_type = 0;
    } else {
        analogy->analogy_type = 1;
    }
}

/* ─── L2: Colburn j-Factors ─────────────────────────────────────── */

void compute_colburn_j_factors(double h, double kc, double rho, double v,
                               double cp, double Pr, double Sc,
                               ChiltonColburnState *state)
{
    if (!state) return;

    /* j_H = St_heat · Pr^(2/3) = (h/(ρ·v·cp)) · Pr^(2/3) */
    double St_h = h / (rho * v * cp);
    state->j_H = St_h * pow(Pr, 2.0 / 3.0);

    /* j_D = St_mass · Sc^(2/3) = (k_c/v) · Sc^(2/3) */
    double St_m = kc / v;
    state->j_D = St_m * pow(Sc, 2.0 / 3.0);

    state->Pr = Pr;
    state->Sc = Sc;
    state->f_F = state->j_H;  /* Analogy prediction: f_F = j_H = j_D */
    state->Cf = 2.0 * state->f_F;
    state->Nu = St_h * state->Re * Pr;   /* needs Re set externally */
    state->Sh = St_m * state->Re * Sc;
}

/* ─── L4: Transport Equation Balance ────────────────────────────── */

double transport_equation_residual(double transient, double convective,
                                   double diffuse, double source)
{
    /* ∂φ/∂t + v·∇φ = Γ·∇²φ + S_φ */
    /* residual = |LHS - RHS| = |transient + convective - diffuse - source| */
    return fabs(transient + convective - diffuse - source);
}

double analogy_closeness_metric(double Pr, double Sc)
{
    /* How close is the system to the ideal analogy condition Pr=Sc=1? */
    /* Returns a value in [0, 1] where 1 = perfect analogy */
    double d_Pr = fabs(Pr - 1.0);
    double d_Sc = fabs(Sc - 1.0);
    double metric = 1.0 / (1.0 + d_Pr + d_Sc);
    return metric;
}

double reynolds_analogy_error(double St_heat, double cf)
{
    /* Reynolds analogy predicts: St = cf/2 */
    double predicted = cf / 2.0;
    if (predicted < 1e-15) return 1.0;  /* degenerate */
    return fabs(St_heat - predicted) / predicted;
}

/* ─── L5: Engineering Methods ──────────────────────────────────── */

/**
 * Pressure drop → heat transfer coefficient prediction.
 *
 * From measured ΔP:
 *   f = 2·ΔP·D_h / (ρ·v²·L)
 *   St·Pr^(2/3) = f/2  (Colburn)
 *   h = St · ρ · v · Cp
 *
 * This is the most important practical application of the analogy:
 * friction measurements (easy, cheap) predict heat transfer (hard to measure).
 *
 * Complexity: O(1)
 */
double predict_h_from_friction(double delta_P, double L, double D_h,
                               double rho, double v, double cp, double Pr)
{
    double f = 2.0 * delta_P * D_h / (rho * v * v * L);
    double f_F = f / 4.0;
    double St = f_F * pow(Pr, -2.0 / 3.0);  /* Colburn form */
    return St * rho * v * cp;
}

/**
 * Heat transfer coefficient → mass transfer coefficient prediction.
 *
 * From the equality j_H = j_D:
 *   k_c = (h/(ρ·cp)) · (Sc/Pr)^(-2/3)  = h/(ρ·cp) · (Pr/Sc)^(2/3)
 *
 * Complexity: O(1)
 */
double predict_kc_from_h(double h, double rho, double cp, double v,
                         double Pr, double Sc)
{
    /* j_H = (h/(ρ·v·cp)) · Pr^(2/3) = j_D = (k_c/v) · Sc^(2/3) */
    double j_H = (h / (rho * v * cp)) * pow(Pr, 2.0 / 3.0);
    /* k_c = j_H · v · Sc^(-2/3) */
    return j_H * v * pow(Sc, -2.0 / 3.0);
}

/* ─── L5: Prandtl-Taylor Analogy ───────────────────────────────── */

/**
 * Prandtl-Taylor two-layer model.
 *
 * The flow is divided into:
 *   1. Laminar sublayer: molecular transport dominates
 *   2. Turbulent core: eddy transport dominates
 *
 * Result:
 *   St = (f/2) / [1 + 5·√(f/2)·(Pr - 1)]
 *
 * This extends Reynolds analogy to Pr ≠ 1 by accounting for
 * the different transport mechanisms in each layer.
 *
 * Complexity: O(1)
 */
double prandtl_taylor_analogy_St(double f_F, double Pr)
{
    double sqrt_f2 = sqrt(f_F);
    double denominator = 1.0 + 5.0 * sqrt_f2 * (Pr - 1.0);
    if (denominator <= 0.0) return f_F;  /* fallback for Pr << 1 */
    return f_F / denominator;
}

/**
 * von Karman analogy (three-layer model).
 *
 * Adds a buffer layer between the laminar sublayer and turbulent core.
 *
 * St = (f/2) / [1 + 5·√(f/2)·{(Pr-1) + ln[1 + 5·(Pr-1)/6]}]
 *
 * More accurate than Prandtl-Taylor for wider Pr ranges.
 *
 * Complexity: O(1)
 */
double von_karman_analogy_St(double f_F, double Pr)
{
    double sqrt_f2 = sqrt(f_F);
    double inner = 1.0 + 5.0 * (Pr - 1.0) / 6.0;
    if (inner <= 0.0) inner = 1e-6;
    double bracket = (Pr - 1.0) + log(inner);
    double denominator = 1.0 + 5.0 * sqrt_f2 * bracket;
    if (denominator <= 0.0) return f_F;
    return f_F / denominator;
}

/**
 * Lyon-Martinelli correlation for liquid metals (Pr << 1).
 *
 * For liquid metals (Pr ≈ 0.001-0.03), molecular conduction dominates
 * even in the turbulent core because thermal diffusivity is much larger
 * than momentum diffusivity.
 *
 * Nu = 7.0 + 0.025 · Pe^(0.8)  (circular tube, uniform heat flux)
 *
 * This is fundamentally different from the standard analogy because
 * the eddy diffusivity for heat (ε_H) can be neglected compared to α.
 *
 * Complexity: O(1)
 */
double lyon_martinelli_nu(double Re, double Pr)
{
    double Pe = Re * Pr;
    return 7.0 + 0.025 * pow(Pe, 0.8);
}

/* ─── L7: Application Utilities ────────────────────────────────── */

double heat_transfer_coefficient(double Nu, double k, double L)
{
    return Nu * k / L;
}

double mass_transfer_coefficient(double Sh, double D_AB, double L)
{
    return Sh * D_AB / L;
}

double friction_from_nusselt(double Nu, double Re, double Pr)
{
    double St = Nu / (Re * Pr);
    /* Reverse Colburn: f_F = St · Pr^(2/3) */
    return St * pow(Pr, 2.0 / 3.0);
}

char *analogy_regime_description(double Re, double Pr, double Sc,
                                 char *buf, size_t size)
{
    if (!buf || size == 0) return NULL;

    char flow[32], pr_desc[64], sc_desc[64], analogy[128];

    /* Flow regime */
    if (Re < 2300.0)
        snprintf(flow, sizeof(flow), "laminar");
    else if (Re < 10000.0)
        snprintf(flow, sizeof(flow), "transitional");
    else
        snprintf(flow, sizeof(flow), "turbulent");

    /* Prandtl meaning */
    if (Pr < 0.1)
        snprintf(pr_desc, sizeof(pr_desc),
                 "liquid metal (Pr=%.3g), thermal BL >> velocity BL", Pr);
    else if (Pr < 0.7)
        snprintf(pr_desc, sizeof(pr_desc),
                 "low-Pr gas (Pr=%.3g), thermal BL > velocity BL", Pr);
    else if (Pr < 1.5)
        snprintf(pr_desc, sizeof(pr_desc),
                 "gas (Pr≈1), thermal BL ≈ velocity BL (ideal analogy)");
    else if (Pr < 10.0)
        snprintf(pr_desc, sizeof(pr_desc),
                 "liquid (Pr=%.3g), thermal BL < velocity BL", Pr);
    else if (Pr < 100.0)
        snprintf(pr_desc, sizeof(pr_desc),
                 "oil-like (Pr=%.3g), thin thermal BL", Pr);
    else
        snprintf(pr_desc, sizeof(pr_desc),
                 "very viscous (Pr=%.3g), very thin thermal BL", Pr);

    /* Schmidt meaning */
    if (Sc < 1.0)
        snprintf(sc_desc, sizeof(sc_desc),
                 "gas diffusion, concentration BL > velocity BL");
    else if (Sc < 10.0)
        snprintf(sc_desc, sizeof(sc_desc),
                 "gas in gas, concentration BL ≈ velocity BL");
    else
        snprintf(sc_desc, sizeof(sc_desc),
                 "liquid/gas in liquid, thin concentration BL");

    /* Analogy applicability */
    double metric = analogy_closeness_metric(Pr, Sc);
    if (Re > 10000.0 && metric > 0.3)
        snprintf(analogy, sizeof(analogy),
                 "Chilton-Colburn analogy APPLICABLE (turbulent, metric=%.2f)",
                 metric);
    else if (Re > 10000.0)
        snprintf(analogy, sizeof(analogy),
                 "Chilton-Colburn MARGINAL (turbulent but Pr/Sc mismatched)");
    else
        snprintf(analogy, sizeof(analogy),
                 "Reynolds analogy NOT RELIABLE (laminar/transitional)");

    snprintf(buf, size,
             "[%s] Re=%.1f | %s | %s | %s",
             flow, Re, pr_desc, sc_desc, analogy);
    return buf;
}

void print_analogy_diagnostics(const ChiltonColburnState *ccs)
{
    if (!ccs) return;
    printf("━━ Chilton-Colburn Analogy Diagnostics ━━━━━━━━━━━━━━━━\n");
    printf("  Re = %.1f    Pr = %.3f    Sc = %.3f\n",
           ccs->Re, ccs->Pr, ccs->Sc);
    printf("  j_H (heat) = %.6f\n", ccs->j_H);
    printf("  j_D (mass) = %.6f\n", ccs->j_D);
    printf("  j_H/j_D    = %.4f  (ideal = 1.0)\n",
           ccs->j_H / (ccs->j_D + 1e-30));
    printf("  f_F = %.6f (from analogy = j_H)\n", ccs->f_F);
    printf("  Cf  = %.6f\n", ccs->Cf);
    printf("  Nu (predicted) = %.1f\n", ccs->Nu);
    printf("  Sh (predicted) = %.1f\n", ccs->Sh);
    printf("  Analogy holds: %s\n",
           fabs(ccs->j_H - ccs->j_D) < 0.1 * ccs->j_H ? "YES ✓" : "marginal");
    printf("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");
}

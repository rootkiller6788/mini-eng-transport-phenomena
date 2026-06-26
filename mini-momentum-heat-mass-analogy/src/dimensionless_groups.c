/**
 * @file dimensionless_groups.c
 * @brief Computation of all dimensionless groups in the analogy framework
 *
 * This module computes every dimensionless number connecting momentum,
 * heat, and mass transfer: Re, Pr, Sc, Le, Nu, Sh, St, Pe, Gr, Ra, Br, Ec.
 *
 * Reference: Bird, Stewart, Lightfoot (2007), Transport Phenomena.
 *            Incropera & DeWitt (2007), Fundamentals of Heat and Mass Transfer.
 *            Buckingham (1914), "On physically similar systems; illustrations
 *            of the use of dimensional equations," Phys. Rev., 4, 345-376.
 */

#include "dimensionless_groups.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

/* ─── L1: Primary Dimensionless Numbers ─────────────────────────── */

double compute_Re(double rho, double v, double L, double mu)
{
    if (mu <= 0.0) return -1.0;
    return rho * v * L / mu;
}

double compute_Pr(double mu, double cp, double k)
{
    if (k <= 0.0) return -1.0;
    return mu * cp / k;
}

double compute_Sc(double mu, double rho, double D)
{
    if (rho <= 0.0 || D <= 0.0) return -1.0;
    return mu / (rho * D);
}

double compute_Le(double alpha, double D)
{
    if (D <= 0.0) return -1.0;
    return alpha / D;
}

double compute_Nu(double h, double L, double k)
{
    if (k <= 0.0) return -1.0;
    return h * L / k;
}

double compute_Sh(double kc, double L, double D)
{
    if (D <= 0.0) return -1.0;
    return kc * L / D;
}

double compute_St_heat(double h, double rho, double v, double cp)
{
    if (rho <= 0.0 || v <= 0.0 || cp <= 0.0) return -1.0;
    return h / (rho * v * cp);
}

double compute_St_mass(double kc, double v)
{
    if (v <= 0.0) return -1.0;
    return kc / v;
}

double compute_Pe(double Re, double Pr)
{
    return Re * Pr;
}

double compute_Gr(double g, double beta, double delta_T, double L, double nu)
{
    if (nu <= 0.0 || nu <= 0.0) return -1.0;
    return g * beta * delta_T * L * L * L / (nu * nu);
}

double compute_Ra(double Gr, double Pr)
{
    return Gr * Pr;
}

double compute_Br(double mu, double v, double k, double delta_T)
{
    if (k <= 0.0 || fabs(delta_T) < 1e-15) return -1.0;
    return mu * v * v / (k * fabs(delta_T));
}

double compute_Ec(double v, double cp, double delta_T)
{
    if (cp <= 0.0 || fabs(delta_T) < 1e-15) return -1.0;
    return v * v / (cp * fabs(delta_T));
}

double compute_Bi(double h, double L_c, double k_solid)
{
    if (k_solid <= 0.0) return -1.0;
    return h * L_c / k_solid;
}

double compute_Fo(double alpha, double t, double L)
{
    if (L <= 0.0) return -1.0;
    return alpha * t / (L * L);
}

/* ─── L1: Complete Dimensionless Group Set ──────────────────────── */

void compute_dimensionless_group_set(double rho, double v, double L,
                                     double mu, double cp, double k,
                                     double D, double h, double kc,
                                     double g, double beta, double delta_T,
                                     double t, double k_solid, double L_c,
                                     double f_F,
                                     DimensionlessGroupSet *set)
{
    if (!set) return;

    set->Re = compute_Re(rho, v, L, mu);
    set->Pr = compute_Pr(mu, cp, k);
    set->Sc = compute_Sc(mu, rho, D);
    set->Le = compute_Le(k / (rho * cp), D);

    set->St_heat = compute_St_heat(h, rho, v, cp);
    set->St_mass = compute_St_mass(kc, v);

    set->Nu = set->St_heat * set->Re * set->Pr;
    set->Sh = set->St_mass * set->Re * set->Sc;

    set->Pe = compute_Pe(set->Re, set->Pr);
    set->Gr = compute_Gr(g, beta, delta_T, L, mu / rho);
    set->Ra = compute_Ra(set->Gr, set->Pr);
    set->Br = compute_Br(mu, v, k, delta_T);
    set->Ec = compute_Ec(v, cp, delta_T);
    set->Bi = compute_Bi(h, L_c, k_solid);
    set->Fo = compute_Fo(k / (rho * cp), t, L);

    set->f_F = f_F;
    set->Cf = 2.0 * f_F;
}

/* ─── L2: Flow Regime Classification ────────────────────────────── */

int flow_regime_pipe(double Re)
{
    if (Re < 0.0) return -1;
    if (Re < 2300.0) return 0;      /* laminar */
    if (Re < 4000.0) return 1;      /* transitional */
    return 2;                        /* turbulent */
}

int flow_regime_flat_plate(double Re_x)
{
    if (Re_x < 0.0) return -1;
    if (Re_x < 5.0e5) return 0;     /* laminar */
    if (Re_x < 3.0e6) return 1;     /* transitional */
    return 2;                        /* turbulent */
}

int natural_convection_regime(double Ra)
{
    if (Ra < 0.0) return -1;
    if (Ra < 1.0e3) return 0;       /* conduction dominant */
    if (Ra < 1.0e9) return 1;       /* laminar natural convection */
    return 2;                        /* turbulent natural convection */
}

int transport_dominance(double Pe, char *description, size_t desc_size)
{
    if (!description || desc_size == 0) return -1;

    if (Pe < 0.1) {
        snprintf(description, desc_size, "Diffusion-dominated (Pe=%.2e << 1)", Pe);
        return -1;
    } else if (Pe < 10.0) {
        snprintf(description, desc_size, "Mixed advection-diffusion (Pe=%.2f ≈ 1)", Pe);
        return 0;
    } else {
        snprintf(description, desc_size, "Advection-dominated (Pe=%.2e >> 1)", Pe);
        return 1;
    }
}

int lumped_capacitance_valid(double Bi)
{
    return Bi < 0.1 ? 1 : 0;
}

/* ─── L3: Typical Engineering Values ────────────────────────────── */

/**
 * Typical Prandtl numbers for common fluids.
 *
 * Air:    Pr ≈ 0.71 (mostly T-independent for ideal gas)
 * Water:  Pr ≈ 13 at 273K, 7 at 293K, 1.8 at 373K
 * Oil:    Pr ≈ 100-10000 (strongly T-dependent)
 * Hg:     Pr ≈ 0.025  (liquid metal)
 * CO2:    Pr ≈ 0.77   (gas)
 * He:     Pr ≈ 0.68   (gas)
 */
double typical_prandtl(int fluid_id, double T)
{
    switch (fluid_id) {
    case 0: /* Air */
        return 0.71;
    case 1: /* Water */
        /* Pr decreases from ~13 at 0°C to ~1.8 at 100°C */
        if (T <= 273.15) return 13.0;
        if (T >= 373.15) return 1.8;
        return 13.0 - 11.2 * (T - 273.15) / 100.0;
    case 2: /* Engine oil (SAE 30) */
        /* Pr ≈ 10000 at 0°C, ≈ 100 at 100°C */
        if (T <= 273.15) return 10000.0;
        if (T >= 373.15) return 100.0;
        return 10000.0 * pow(100.0 / 10000.0,
                             (T - 273.15) / 100.0);
    case 3: /* Mercury */
        return 0.0248;
    case 4: /* CO2 */
        return 0.77;
    case 5: /* Helium */
        return 0.68;
    default:
        return 1.0;
    }
}

/**
 * Typical Schmidt numbers for gas-in-air systems.
 *
 * H2O-air:   Sc ≈ 0.6   (water vapor ≈ same size as air molecules)
 * CO2-air:   Sc ≈ 0.96
 * O2-air:    Sc ≈ 0.74
 * NH3-air:   Sc ≈ 0.78
 * Ethanol:   Sc ≈ 1.3
 */
double typical_schmidt(int system_id, double T)
{
    /* Schmidt number for gases is weakly T-dependent */
    (void)T;  /* temperature dependence is small for gases */

    switch (system_id) {
    case 0: return 0.60;   /* H2O in air */
    case 1: return 0.96;   /* CO2 in air */
    case 2: return 0.74;   /* O2 in air */
    case 3: return 0.78;   /* NH3 in air */
    case 4: return 1.30;   /* ethanol in air */
    default: return 1.0;
    }
}

void prandtl_physical_meaning(double Pr, char *buf, size_t buf_size)
{
    if (!buf || buf_size == 0) return;

    if (Pr < 0.01) {
        snprintf(buf, buf_size,
                 "Pr=%.4f: Liquid metal — thermal diffusivity dominates. "
                 "δ_T >> δ. Heat spreads rapidly through conduction. "
                 "Reynolds analogy fails (use Lyon-Martinelli).", Pr);
    } else if (Pr < 0.5) {
        snprintf(buf, buf_size,
                 "Pr=%.3f: Low-Pr gas — thermal BL thicker than velocity BL. "
                 "δ_T > δ. Chilton-Colburn valid with Pr^(2/3) correction.", Pr);
    } else if (Pr < 1.5) {
        snprintf(buf, buf_size,
                 "Pr=%.3f: Gas (Pr≈1) — ideal Reynolds analogy condition. "
                 "δ_T ≈ δ. Momentum and heat transfer are directly linked. "
                 "St ≈ f/2 holds approximately.", Pr);
    } else if (Pr < 10.0) {
        snprintf(buf, buf_size,
                 "Pr=%.3f: Liquid (e.g., water) — velocity BL thicker than "
                 "thermal BL. δ > δ_T. Chilton-Colburn adequate.", Pr);
    } else if (Pr < 100.0) {
        snprintf(buf, buf_size,
                 "Pr=%.3f: Oil — thermal BL is thin. δ >> δ_T. "
                 "Heat transfer dominated by thin thermal sublayer. "
                 "Use Petukhov or Gnielinski correlations.", Pr);
    } else {
        snprintf(buf, buf_size,
                 "Pr=%.3f: Very viscous fluid — extremely thin thermal BL. "
                 "δ >> δ_T. Molecular conduction controls heat transfer. "
                 "Simple analogies may not be accurate.", Pr);
    }
}

/* ─── L7: Diagnostics ──────────────────────────────────────────── */

void print_dimensionless_summary(const DimensionlessGroupSet *set)
{
    if (!set) return;

    printf("━━ Dimensionless Group Summary ━━━━━━━━━━━━━━━━━━━━━━━\n");
    printf("  ── Primary Transport Numbers ──\n");
    printf("  Re  = %8.1f  (inertia/viscous)\n", set->Re);
    printf("  Pr  = %8.4f  (momentum/thermal diffusivity)\n", set->Pr);
    printf("  Sc  = %8.4f  (momentum/mass diffusivity)\n", set->Sc);
    printf("  Le  = %8.4f  (thermal/mass diffusivity = Sc/Pr)\n", set->Le);
    printf("  ── Heat Transfer Numbers ──\n");
    printf("  Nu  = %8.2f  (convective/conductive heat transfer)\n", set->Nu);
    printf("  St  = %8.6f  (heat transferred/thermal capacity)\n", set->St_heat);
    printf("  Pe  = %8.1f  (advection/diffusion heat)\n", set->Pe);
    printf("  ── Mass Transfer Numbers ──\n");
    printf("  Sh  = %8.2f  (convective/diffusive mass transfer)\n", set->Sh);
    printf("  St_m= %8.6f  (mass transferred/flow capacity)\n", set->St_mass);
    printf("  ── Natural Convection ──\n");
    printf("  Gr  = %8.2e  (buoyancy/viscous)\n", set->Gr);
    printf("  Ra  = %8.2e  (buoyancy/diffusion)\n", set->Ra);
    printf("  ── Special Numbers ──\n");
    printf("  Br  = %8.4f  (viscous heating/conduction)\n", set->Br);
    printf("  Ec  = %8.4f  (kinetic energy/enthalpy)\n", set->Ec);
    printf("  Bi  = %8.4f  (internal/external resistance)\n", set->Bi);
    printf("  Fo  = %8.4f  (dimensionless time)\n", set->Fo);
    printf("  ── Friction ──\n");
    printf("  f_F  = %8.6f  (Fanning friction factor)\n", set->f_F);
    printf("  Cf   = %8.6f  (skin friction coefficient)\n", set->Cf);
    printf("  ── Analogy Check ──\n");
    printf("  St·Pr^(2/3) / f_F = %.3f  (should = 1 for Colburn)\n",
           set->St_heat * pow(set->Pr, 2.0/3.0) / (set->f_F + 1e-30));
    printf("  St_m·Sc^(2/3) / f_F = %.3f  (should = 1 for Colburn)\n",
           set->St_mass * pow(set->Sc, 2.0/3.0) / (set->f_F + 1e-30));
    printf("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");
}

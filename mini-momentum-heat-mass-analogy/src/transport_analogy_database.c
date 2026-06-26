/**
 * @file transport_analogy_database.c
 * @brief Engineering database of fluid properties for analogy calculations
 *
 * Provides built-in property data for 12 common engineering fluids and
 * utility functions for quick engineering estimates using the analogy.
 *
 * Reference: Incropera & DeWitt (2007), Appendix A.
 *            Cengel (2014), Appendix 1.
 *            Perry's Chemical Engineers' Handbook (2019).
 *            NIST REFPROP database values.
 */

#include "transport_analogy_database.h"
#include "momentum_heat_mass_analogy.h"
#include "boundary_layer_analogy.h"
#include "tube_channel_analogy.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

/* ─── L3: Built-in Fluid Database ───────────────────────────────── */

void init_fluid_database(FluidDatabaseEntry db[N_BUILTIN_FLUIDS])
{
    if (!db) return;

    /* 0: Air */
    {
        FluidDatabaseEntry *f = &db[0];
        snprintf(f->name, sizeof(f->name), "Air");
        f->T_min = 150.0;   f->T_max = 1800.0;
        f->M = 0.02897;     f->Tc = 132.5;
        f->Pc = 3.77e6;     f->rho_ref = 1.184;
        f->T_ref = 298.15;  f->is_gas = 1;  f->fluid_id = 0;
        f->sutherland = (SutherlandParams){1.716e-5, 273.15, 110.4};
        f->k_at_ref = 0.02624;  f->cp_at_ref = 1007.0;
    }

    /* 1: Water (liquid) */
    {
        FluidDatabaseEntry *f = &db[1];
        snprintf(f->name, sizeof(f->name), "Water");
        f->T_min = 273.15;  f->T_max = 373.15;
        f->M = 0.018015;    f->Tc = 647.1;
        f->Pc = 22.06e6;    f->rho_ref = 998.0;
        f->T_ref = 293.15;  f->is_gas = 0;  f->fluid_id = 1;
        f->plaw = (PowerLawViscosityParams){1.002e-3, 293.15, -1.5};
        f->k_at_ref = 0.603;  f->cp_at_ref = 4182.0;
    }

    /* 2: Engine Oil (SAE 30) */
    {
        FluidDatabaseEntry *f = &db[2];
        snprintf(f->name, sizeof(f->name), "Engine Oil SAE 30");
        f->T_min = 273.15;  f->T_max = 423.15;
        f->M = 0.400;       f->Tc = 800.0;
        f->Pc = 1.5e6;      f->rho_ref = 888.0;
        f->T_ref = 293.15;  f->is_gas = 0;  f->fluid_id = 2;
        f->plaw = (PowerLawViscosityParams){0.486, 293.15, -3.0};
        f->k_at_ref = 0.145;  f->cp_at_ref = 1900.0;
    }

    /* 3: Mercury (liquid metal) */
    {
        FluidDatabaseEntry *f = &db[3];
        snprintf(f->name, sizeof(f->name), "Mercury");
        f->T_min = 234.0;   f->T_max = 600.0;
        f->M = 0.20059;     f->Tc = 1750.0;
        f->Pc = 172.0e6;    f->rho_ref = 13529.0;
        f->T_ref = 300.0;   f->is_gas = 0;  f->fluid_id = 3;
        f->plaw = (PowerLawViscosityParams){1.53e-3, 300.0, -0.1};
        f->k_at_ref = 8.54;  f->cp_at_ref = 139.0;
    }

    /* 4: Carbon Dioxide (gas) */
    {
        FluidDatabaseEntry *f = &db[4];
        snprintf(f->name, sizeof(f->name), "Carbon Dioxide");
        f->T_min = 220.0;   f->T_max = 1500.0;
        f->M = 0.04401;     f->Tc = 304.2;
        f->Pc = 7.38e6;     f->rho_ref = 1.773;
        f->T_ref = 300.0;   f->is_gas = 1;  f->fluid_id = 4;
        f->sutherland = (SutherlandParams){1.37e-5, 273.15, 254.0};
        f->k_at_ref = 0.0166;  f->cp_at_ref = 851.0;
    }

    /* 5: Helium */
    {
        FluidDatabaseEntry *f = &db[5];
        snprintf(f->name, sizeof(f->name), "Helium");
        f->T_min = 4.0;     f->T_max = 1500.0;
        f->M = 0.0040026;   f->Tc = 5.19;
        f->Pc = 0.227e6;    f->rho_ref = 0.1614;
        f->T_ref = 300.0;   f->is_gas = 1;  f->fluid_id = 5;
        f->sutherland = (SutherlandParams){1.86e-5, 273.15, 79.4};
        f->k_at_ref = 0.156;  f->cp_at_ref = 5193.0;
    }

    /* 6: Hydrogen */
    {
        FluidDatabaseEntry *f = &db[6];
        snprintf(f->name, sizeof(f->name), "Hydrogen");
        f->T_min = 20.0;    f->T_max = 1500.0;
        f->M = 0.002016;    f->Tc = 33.2;
        f->Pc = 1.30e6;     f->rho_ref = 0.08078;
        f->T_ref = 300.0;   f->is_gas = 1;  f->fluid_id = 6;
        f->sutherland = (SutherlandParams){8.41e-6, 273.15, 97.0};
        f->k_at_ref = 0.187;  f->cp_at_ref = 14310.0;
    }

    /* 7: Steam (H2O gas) */
    {
        FluidDatabaseEntry *f = &db[7];
        snprintf(f->name, sizeof(f->name), "Steam (H2O gas)");
        f->T_min = 380.0;   f->T_max = 1500.0;
        f->M = 0.018015;    f->Tc = 647.1;
        f->Pc = 22.06e6;    f->rho_ref = 0.5542;
        f->T_ref = 400.0;   f->is_gas = 1;  f->fluid_id = 7;
        f->sutherland = (SutherlandParams){1.32e-5, 400.0, 673.0};
        f->k_at_ref = 0.0261;  f->cp_at_ref = 2010.0;
    }

    /* 8: Ethylene Glycol (50% aq) */
    {
        FluidDatabaseEntry *f = &db[8];
        snprintf(f->name, sizeof(f->name), "Ethylene Glycol 50%%");
        f->T_min = 253.0;   f->T_max = 373.0;
        f->M = 0.06207;     f->Tc = 720.0;
        f->Pc = 7.5e6;      f->rho_ref = 1075.0;
        f->T_ref = 293.15;  f->is_gas = 0;  f->fluid_id = 8;
        f->plaw = (PowerLawViscosityParams){3.94e-3, 293.15, -2.5};
        f->k_at_ref = 0.427;  f->cp_at_ref = 3380.0;
    }

    /* 9: Glycerin */
    {
        FluidDatabaseEntry *f = &db[9];
        snprintf(f->name, sizeof(f->name), "Glycerin");
        f->T_min = 291.0;   f->T_max = 373.0;
        f->M = 0.09209;     f->Tc = 850.0;
        f->Pc = 7.5e6;      f->rho_ref = 1261.0;
        f->T_ref = 293.15;  f->is_gas = 0;  f->fluid_id = 9;
        f->plaw = (PowerLawViscosityParams){1.5, 293.15, -4.0};
        f->k_at_ref = 0.286;  f->cp_at_ref = 2430.0;
    }

    /* 10: Nitrogen */
    {
        FluidDatabaseEntry *f = &db[10];
        snprintf(f->name, sizeof(f->name), "Nitrogen");
        f->T_min = 70.0;    f->T_max = 1500.0;
        f->M = 0.028013;    f->Tc = 126.2;
        f->Pc = 3.39e6;     f->rho_ref = 1.123;
        f->T_ref = 300.0;   f->is_gas = 1;  f->fluid_id = 10;
        f->sutherland = (SutherlandParams){1.66e-5, 273.15, 107.0};
        f->k_at_ref = 0.0259;  f->cp_at_ref = 1041.0;
    }

    /* 11: Ammonia */
    {
        FluidDatabaseEntry *f = &db[11];
        snprintf(f->name, sizeof(f->name), "Ammonia");
        f->T_min = 240.0;   f->T_max = 800.0;
        f->M = 0.017031;    f->Tc = 405.4;
        f->Pc = 11.33e6;    f->rho_ref = 0.6894;
        f->T_ref = 300.0;   f->is_gas = 1;  f->fluid_id = 11;
        f->sutherland = (SutherlandParams){9.3e-6, 273.15, 527.0};
        f->k_at_ref = 0.0247;  f->cp_at_ref = 2190.0;
    }
}

/* ─── L3: Fluid Lookup ──────────────────────────────────────────── */

int fluid_analogy_lookup(const FluidDatabaseEntry db[N_BUILTIN_FLUIDS],
                          int id, double T, double P,
                          AnalogyLookup *lookup)
{
    if (!db || !lookup || id < 0 || id >= N_BUILTIN_FLUIDS) return -1;

    const FluidDatabaseEntry *f = &db[id];

    if (T < f->T_min || T > f->T_max) return -2;

    lookup->fluid = *f;
    lookup->T = T;
    lookup->P = P;

    /* Compute transport properties */
    if (f->is_gas) {
        lookup->transport.mu = sutherland_viscosity(T, &f->sutherland);
        lookup->transport.k  = eucken_thermal_conductivity(
                                  lookup->transport.mu, f->cp_at_ref, f->M);
    } else {
        lookup->transport.mu = power_law_viscosity(T, &f->plaw);
        lookup->transport.k  = f->k_at_ref * pow(T / f->T_ref, 0.3);
        /* Liquid k varies only slightly with T */
    }

    lookup->transport.rho = f->is_gas
        ? P * f->M / (8.314462618 * T)  /* ideal gas */
        : f->rho_ref;                    /* liquid: constant approximation */

    lookup->transport.cp = f->cp_at_ref;  /* first approximation */
    lookup->transport.T = T;
    lookup->transport.P = P;
    lookup->transport.D = 1.0e-5;  /* placeholder — requires binary system */

    /* Prandtl */
    lookup->Pr = lookup->transport.mu * lookup->transport.cp
                 / lookup->transport.k;

    /* Re at 1 m/s, 1 mm */
    lookup->Re_1ms_1mm = lookup->transport.rho * 1.0 * 0.001
                         / lookup->transport.mu;

    /* Sc for water vapor in this fluid (approximate) */
    lookup->Sc_water_air = 0.6;  /* default for gases */

    /* Analogy regime description */
    double metric = 1.0 / (1.0 + fabs(lookup->Pr - 1.0) + fabs(lookup->Sc_water_air - 1.0));
    lookup->analogy_applicable = (metric > 0.3) ? 1 : 0;

    if (lookup->analogy_applicable) {
        snprintf(lookup->analogy_regime, sizeof(lookup->analogy_regime),
                 "Analogy APPLICABLE (Pr=%.3f, Sc≈%.3f, metric=%.2f)",
                 lookup->Pr, lookup->Sc_water_air, metric);
    } else {
        snprintf(lookup->analogy_regime, sizeof(lookup->analogy_regime),
                 "Analogy MARGINAL (Pr=%.3f, Sc≈%.3f, metric=%.2f) — "
                 "use extended correlations",
                 lookup->Pr, lookup->Sc_water_air, metric);
    }

    return 0;
}

/* ─── L6: Engineering Benchmarks ────────────────────────────────── */

/**
 * Air flat plate benchmark.
 *
 * Standard conditions: air at 300K, 1 atm, U∞ = 10 m/s, L = 1 m.
 * Expected: Re_L ≈ 6.9×10⁵, Cf ≈ 0.00265, h ≈ 15 W/(m²·K).
 */
void air_flat_plate_benchmark(double U_inf, double L, double T_inf,
                               double T_wall, double P_atm,
                               BoundaryLayerSummary *result)
{
    if (!result) return;

    FlatPlateState state;
    state.x = L;
    state.U_inf = U_inf;
    state.T_inf = T_inf;
    state.T_w = T_wall;
    state.C_inf = 0.0;
    state.C_w = 0.0;

    /* Air properties at film temperature */
    double T_film = (T_inf + T_wall) / 2.0;
    state.rho = P_atm * 0.02897 / (8.314462618 * T_film);
    state.mu = air_viscosity(T_film);
    double k_air = air_thermal_conductivity(T_film);
    double cp_air = 1007.0;

    state.nu = state.mu / state.rho;
    state.alpha = k_air / (state.rho * cp_air);
    state.D_AB = water_vapor_air_diffusivity(T_film, P_atm);

    state.Re_x = state.rho * U_inf * L / state.mu;
    state.Pr = state.mu * cp_air / k_air;
    state.Sc = state.mu / (state.rho * state.D_AB);

    /* Determine flow regime */
    if (state.Re_x < 5.0e5)
        state.flow_regime = 0;
    else if (state.Re_x < 3.0e6)
        state.flow_regime = 1;
    else
        state.flow_regime = 2;

    compute_boundary_layer_summary(&state, result);
}

/**
 * Water pipe flow benchmark.
 *
 * Standard: D = 0.0254 m (1 in), v = 1 m/s, L = 10 m.
 * Expected: Re ≈ 3.0×10⁴, Nu ≈ 176, h ≈ 4240 W/(m²·K).
 */
void water_pipe_flow_benchmark(double D, double v, double L,
                                double T_bulk, double T_wall,
                                PipeFlowAnalogy *result)
{
    if (!result) return;

    double rho = 998.0;     /* water at ~20°C */
    double mu  = water_viscosity(T_bulk);
    double cp  = 4182.0;
    double k   = water_thermal_conductivity(T_bulk);
    double D_AB = 1.0e-9;   /* typical liquid diffusivity */

    compute_pipe_flow_analogy(D, L, v, rho, mu, cp, k, D_AB,
                              T_wall, T_bulk,
                              0.0, 0.0,  /* concentration */
                              result);
}

void compare_fluid_analogies(int fluid_id_A, int fluid_id_B,
                             double T, double P, double v, double L)
{
    FluidDatabaseEntry db[N_BUILTIN_FLUIDS];
    init_fluid_database(db);

    AnalogyLookup lookupA, lookupB;
    int retA = fluid_analogy_lookup(db, fluid_id_A, T, P, &lookupA);
    int retB = fluid_analogy_lookup(db, fluid_id_B, T, P, &lookupB);

    if (retA != 0 || retB != 0) {
        printf("  Fluid lookup failed.\n");
        return;
    }

    /* Transport properties comparison */
    double Re_A = lookupA.transport.rho * v * L / lookupA.transport.mu;
    double Re_B = lookupB.transport.rho * v * L / lookupB.transport.mu;

    printf("━━━ Fluid Analogy Comparison ━━━━━━━━━━━━━━━━━━━━━━━━━\n");
    printf("  %-20s vs %-20s at T=%.1f K\n",
           lookupA.fluid.name, lookupB.fluid.name, T);
    printf("  ─────────┬─────────────────┬─────────────────\n");
    printf("  Property │ %-15s │ %-15s\n",
           lookupA.fluid.name, lookupB.fluid.name);
    printf("  ─────────┼─────────────────┼─────────────────\n");
    printf("  μ [Pa·s] │ %15.3e │ %15.3e\n",
           lookupA.transport.mu, lookupB.transport.mu);
    printf("  k [W/mK] │ %15.4f │ %15.4f\n",
           lookupA.transport.k, lookupB.transport.k);
    printf("  Pr       │ %15.3f │ %15.3f\n",
           lookupA.Pr, lookupB.Pr);
    printf("  Re(at v,L)│ %15.1f │ %15.1f\n", Re_A, Re_B);
    printf("  Analogy  │ %-15s │ %-15s\n",
           lookupA.analogy_applicable ? "GOOD" : "MARGINAL",
           lookupB.analogy_applicable ? "GOOD" : "MARGINAL");
    printf("  ─────────┴─────────────────┴─────────────────\n");
    printf("  Key insight: Pr=%.1f fluid has δ_T/δ ≈ %.1f\n",
           lookupA.Pr, pow(lookupA.Pr, -1.0/3.0));
    printf("               Pr=%.1f fluid has δ_T/δ ≈ %.1f\n",
           lookupB.Pr, pow(lookupB.Pr, -1.0/3.0));
    printf("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");
}

/* ─── L7: Engineering Applications ──────────────────────────────── */

/**
 * Heat exchanger sizing using the analogy.
 *
 * Starting from friction data (pressure drop), use the analogy to
 * estimate h, then size the heat exchanger.
 *
 * A = Q / (h · ΔT_lm)
 *
 * Complexity: O(1)
 */
void heat_exchanger_sizing_by_analogy(double Q_dot, double delta_T_lm,
                                       double rho, double v, double cp,
                                       double k, double mu, double D_h,
                                       double Pr, double Re,
                                       double *area)
{
    if (!area) return;
    /* Use Colburn analogy to get h */
    /* f ≈ 0.046·Re^(-0.2) for turbulent flow */
    (void)k; (void)mu; (void)D_h;
    double f = (Re < 2300.0) ? 64.0 / Re : 0.046 * pow(Re, -0.2);
    double f_F = f / 4.0;
    double St = f_F * pow(Pr, -2.0 / 3.0);
    double h_pred = St * rho * v * cp;
    *area = Q_dot / (h_pred * delta_T_lm);
}

double cooling_tower_evaporation_rate(double h, double rho, double cp,
                                       double Le, double A,
                                       double delta_omega)
{
    /* Merkel's method: using Lewis relation */
    /* k_c = h/(ρ·cp·Le^(2/3)) */
    double kc = h / (rho * cp * pow(Le, 2.0 / 3.0));
    return kc * A * delta_omega;
}

double reactor_wall_mass_transfer(double delta_P, double D_tube,
                                  double L_tube, double rho, double v,
                                  double D_AB, double Sc)
{
    /* From ΔP, get f, then Sh, then k_c */
    if (L_tube <= 0.0 || rho <= 0.0 || v <= 0.0) return -1.0;
    double f_D = 2.0 * delta_P * D_tube / (rho * v * v * L_tube);
    double Re_D = rho * v * D_tube / (rho * D_AB * Sc);  /* μ = ρ·D_AB·Sc */
    /* Better: use Sh correlation */
    double Sh = predict_sh_from_friction(f_D, Re_D, Sc);
    return Sh * D_AB / D_tube;
}

/**
 * Electronics cooling analogy.
 *
 * Forced convection over a PCB/chip: predict h from flow conditions
 * using Colburn analogy with standard air properties.
 *
 * Complexity: O(1)
 */
double electronics_cooling_analogy(double v_air, double L_chip,
                                   double T_chip, double T_ambient,
                                   double P_atm, double *h)
{
    double T_film = (T_chip + T_ambient) / 2.0;
    double rho = P_atm * 0.02897 / (8.314462618 * T_film);
    double mu = air_viscosity(T_film);
    double k = air_thermal_conductivity(T_film);
    double cp = 1007.0;

    double Re_L = rho * v_air * L_chip / mu;
    double Pr = mu * cp / k;

    /* Nu from laminar or turbulent correlation */
    double Nu;
    if (Re_L < 5.0e5) {
        Nu = 0.664 * sqrt(Re_L) * pow(Pr, 1.0/3.0);
    } else {
        Nu = 0.037 * pow(Re_L, 0.8) * pow(Pr, 1.0/3.0);
    }

    double h_local = Nu * k / L_chip;
    if (h) *h = h_local;

    return h_local * (T_chip - T_ambient);  /* heat flux [W/m²] */
}

void quick_analogy_estimate(const char *fluid_name,
                            double T_C, double P_kPa,
                            double v_ms, double L_mm,
                            char *summary, size_t summary_size)
{
    if (!summary || summary_size == 0) return;

    double T = T_C + 273.15;
    double P = P_kPa * 1000.0;
    double v = v_ms;
    double L = L_mm / 1000.0;

    /* Find fluid ID by name */
    FluidDatabaseEntry db[N_BUILTIN_FLUIDS];
    init_fluid_database(db);

    int id = -1;
    for (int i = 0; i < N_BUILTIN_FLUIDS; i++) {
        if (strstr(db[i].name, fluid_name)) {
            id = i;
            break;
        }
    }

    if (id < 0) {
        snprintf(summary, summary_size,
                 "Fluid '%s' not found in database.", fluid_name);
        return;
    }

    AnalogyLookup lookup;
    int ret = fluid_analogy_lookup(db, id, T, P, &lookup);
    if (ret != 0) {
        snprintf(summary, summary_size,
                 "Error in fluid lookup (code %d).", ret);
        return;
    }

    double Re = lookup.transport.rho * v * L / lookup.transport.mu;
    double Pr = lookup.Pr;

    /* Quick friction estimate */
    double f_D;
    if (Re < 2300.0)
        f_D = 64.0 / Re;
    else
        f_D = 0.046 * pow(Re, -0.2);

    double f_F = f_D / 4.0;
    double St = f_F * pow(Pr, -2.0/3.0);
    double Nu = St * Re * Pr;
    double h = Nu * lookup.transport.k / L;
    double delta_P = f_D * (L / L) * 0.5 * lookup.transport.rho * v * v;
    /* Actually for a duct of length L: */
    delta_P = f_D * 1.0 * 0.5 * lookup.transport.rho * v * v;  /* L/D ≈ 1 */

    snprintf(summary, summary_size,
             "Fluid: %s | T=%.1f°C | Re=%.1f | Pr=%.3f | "
             "f=%.5f | Nu=%.1f | h=%.1f W/m²K | ΔP/L≈%.1f Pa/m | "
             "Analogy: %s",
             lookup.fluid.name, T_C, Re, Pr, f_D, Nu, h, delta_P,
             lookup.analogy_applicable ? "GOOD" : "MARGINAL");
}

void print_analogy_vs_experiment_comparison(void)
{
    printf("━━ Analogy vs. Experimental Correlation Comparison ━━━━━━━\n");
    printf("  Configuration        │ Analogy Nu  │ Dittus-Boelter  │ Error\n");
    printf("  ─────────────────────┼─────────────┼─────────────────┼──────\n");

    /* Test cases */
    struct {
        double Re, Pr;
        const char *label;
    } cases[] = {
        {10000.0, 0.71, "Air, Re=10k        "},
        {50000.0, 0.71, "Air, Re=50k        "},
        {20000.0, 7.0,  "Water, Re=20k      "},
        {100000.0, 7.0, "Water, Re=100k     "},
        {30000.0, 50.0, "Oil, Re=30k, Pr=50 "},
    };
    int n_cases = sizeof(cases) / sizeof(cases[0]);

    for (int i = 0; i < n_cases; i++) {
        double Re = cases[i].Re, Pr = cases[i].Pr;
        double f = 0.046 * pow(Re, -0.2);
        double Nu_analogy = predict_nu_from_friction(f, Re, Pr);
        double Nu_db = dittus_boelter_nu(Re, Pr, 1);
        double err = fabs(Nu_analogy - Nu_db) / Nu_db * 100.0;
        printf("  %s │ %8.1f    │ %8.1f       │ %5.1f%%\n",
               cases[i].label, Nu_analogy, Nu_db, err);
    }
    printf("  ─────────────────────┴─────────────┴─────────────────┴──────\n");
    printf("  Note: Analogy Nu uses f=0.046·Re^(-0.2) with Colburn form.\n");
    printf("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");
}

void print_multi_fluid_analogy_table(double T, double P)
{
    FluidDatabaseEntry db[N_BUILTIN_FLUIDS];
    init_fluid_database(db);

    printf("━━ Multi-Fluid Analogy Table at T=%.1f K, P=%.2f kPa ━━━━━━━━━━━\n",
           T, P / 1000.0);
    printf("  %-20s %6s %10s %8s %8s %10s %s\n",
           "Fluid", "Phase", "μ [Pa·s]", "k [W/mK]", "Pr", "Re@1m/s,1mm", "Analogy");
    printf("  ──────────────────── ────── ────────── ──────── ──────── ────────── ──────\n");

    for (int i = 0; i < N_BUILTIN_FLUIDS; i++) {
        AnalogyLookup lookup;
        int ret = fluid_analogy_lookup(db, i, T, P, &lookup);
        if (ret != 0) {
            printf("  %-20s [out of T range]\n", db[i].name);
            continue;
        }
        printf("  %-20s %-6s %10.3e %8.4f %8.3f %10.1f %s\n",
               lookup.fluid.name,
               lookup.fluid.is_gas ? "gas" : "liquid",
               lookup.transport.mu,
               lookup.transport.k,
               lookup.Pr,
               lookup.Re_1ms_1mm,
               lookup.analogy_applicable ? "✓ GOOD" : "⚠ marginal");
    }
    printf("  ──────────────────── ────── ────────── ──────── ──────── ────────── ──────\n");
    printf("  Legend: ✓ GOOD = Reynolds/Colburn analogy applicable (Pr≈1)\n");
    printf("          ⚠ marginal = analogy requires Pr correction (>20%% error without)\n");
    printf("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");
}

/**
 * @file example_fluid_comparison.c
 * @brief End-to-end example: Comparing analogy behavior across fluids
 *
 * Demonstrates how the momentum-heat-mass analogy applies differently
 * to gases, liquids, and liquid metals, and why the Chilton-Colburn
 * correction (Pr^(2/3)) is essential.
 *
 * Fluids compared at 300 K:
 *   - Air (Pr ≈ 0.71):  near-ideal analogy, gases
 *   - Water (Pr ≈ 7):   thermal BL thinner than velocity BL
 *   - Engine Oil (Pr ≈ 3000): very thin thermal BL, analogy marginal
 *   - Mercury (Pr ≈ 0.025): liquid metal, thermal BL much thicker
 *
 * For each fluid, we compute:
 *   1. Transport properties
 *   2. Dimensionless numbers (Re, Pr, Sc, Le)
 *   3. Analogy predictions vs. correlations
 *   4. Analogy error assessment
 *
 * Reference: Bird, Stewart, Lightfoot (2007), Chapter 13.
 */

#include "transport_coefficients.h"
#include "dimensionless_groups.h"
#include "momentum_heat_mass_analogy.h"
#include "tube_channel_analogy.h"
#include "transport_analogy_database.h"
#include <math.h>
#include <stdio.h>

int main(void)
{
    printf("╔══════════════════════════════════════════════════════════════╗\n");
    printf("║  Multi-Fluid Analogy Comparison                             ║\n");
    printf("╚══════════════════════════════════════════════════════════════╝\n\n");

    /* ── Common Conditions ────────────────────────────────────── */
    double T      = 300.0;      /* Temperature [K] */
    double P      = 101325.0;   /* Pressure [Pa] */
    double v_flow = 10.0;       /* Characteristic velocity [m/s] */
    double L_char = 1.0;        /* Characteristic length [m] */

    printf("Common Conditions: T = %.1f K, P = %.1f kPa, v = %.1f m/s, L = %.2f m\n\n",
           T, P/1000.0, v_flow, L_char);

    /* ── Fluid Database ────────────────────────────────────────── */
    FluidDatabaseEntry db[N_BUILTIN_FLUIDS];
    init_fluid_database(db);

    /* Fluids to compare */
    int fluid_ids[] = {0, 1, 2, 3};  /* Air, Water, Oil, Mercury */
    const char *fluid_names[] = {"Air", "Water", "Engine Oil SAE 30", "Mercury"};
    int n_fluids = 4;

    printf("╔══════════════════════════════════════════════════════════════╗\n");
    printf("║  Transport Properties                                       ║\n");
    printf("╚══════════════════════════════════════════════════════════════╝\n\n");

    printf("  %-20s │ %8s │ %10s │ %8s │ %8s │ %8s │ %8s\n",
           "Fluid", "μ [Pa·s]", "k [W/m·K]", "ρ [kg/m³]", "Pr", "Re", "Analogy?");
    printf("  ─────────────────────┼──────────┼────────────┼──────────┼──────────┼──────────┼────────\n");

    for (int i = 0; i < n_fluids; i++) {
        AnalogyLookup lookup;
        (void)fluid_analogy_lookup(db, fluid_ids[i], T, P, &lookup);

        double Re = lookup.transport.rho * v_flow * L_char / lookup.transport.mu;

        printf("  %-20s │ %8.2e │ %10.4f │ %8.2f │ %8.3f │ %8.1e │ %s\n",
               fluid_names[i],
               lookup.transport.mu,
               lookup.transport.k,
               lookup.transport.rho,
               lookup.Pr,
               Re,
               lookup.analogy_applicable ? "✓ YES" : "⚠ NO");
    }

    printf("\n╔══════════════════════════════════════════════════════════════╗\n");
    printf("║  Analogy Prediction Accuracy (Pipe Flow)                    ║\n");
    printf("╚══════════════════════════════════════════════════════════════╝\n\n");

    printf("  %-20s │ %8s │ %8s │ %8s │ %8s\n",
           "Fluid", "Nu (corr)", "Nu (anal)", "Error %", "j_H/f_F");
    printf("  ─────────────────────┼──────────┼──────────┼──────────┼──────────\n");

    for (int i = 0; i < n_fluids; i++) {
        AnalogyLookup lookup;
        fluid_analogy_lookup(db, fluid_ids[i], T, P, &lookup);

        double mu  = lookup.transport.mu;
        double rho = lookup.transport.rho;
        double Pr  = lookup.Pr;

        /* Use 25.4 mm pipe for comparison */
        double D_pipe = 0.0254;
        double Re_D   = rho * v_flow * D_pipe / mu;

        /* Correlation Nusselt (Dittus-Boelter or laminar) */
        double Nu_corr;
        if (Re_D > 10000.0) {
            Nu_corr = dittus_boelter_nu(Re_D, Pr, 1);
        } else if (Re_D > 2300.0) {
            Nu_corr = gnielinski_nu(Re_D, Pr, blasius_friction_turbulent(Re_D));
        } else {
            Nu_corr = laminar_nu_pipe(0);
        }

        /* Analogy Nusselt from friction */
        double f_D;
        if (Re_D < 2300.0)
            f_D = darcy_friction_laminar(Re_D);
        else
            f_D = blasius_friction_turbulent(Re_D);

        double Nu_anal = predict_nu_from_friction(f_D, Re_D, Pr);
        double err_pct = fabs(Nu_anal - Nu_corr) / Nu_corr * 100.0;

        /* j_H / f_F ratio */
        double j_H = (Nu_corr / (Re_D * Pr)) * pow(Pr, 2.0/3.0);
        double f_F = f_D / 4.0;
        double j_ratio = j_H / (f_F + 1e-30);

        printf("  %-20s │ %8.1f │ %8.1f │ %7.1f%% │ %8.3f\n",
               fluid_names[i], Nu_corr, Nu_anal, err_pct, j_ratio);
    }

    printf("\n╔══════════════════════════════════════════════════════════════╗\n");
    printf("║  Boundary Layer Thickness Comparison                        ║\n");
    printf("╚══════════════════════════════════════════════════════════════╝\n\n");

    printf("  Fluid                │ Pr      │ δ/δ_ref │ δ_T/δ   │ δ_T/δ_v\n");
    printf("  ─────────────────────┼─────────┼─────────┼─────────┼─────────\n");

    /* Use air as reference for velocity BL thickness */
    AnalogyLookup air_lookup;
    fluid_analogy_lookup(db, 0, T, P, &air_lookup);
    double Re_air = air_lookup.transport.rho * v_flow * L_char / air_lookup.transport.mu;
    double delta_air = laminar_bl_thickness(L_char, Re_air);

    for (int i = 0; i < n_fluids; i++) {
        AnalogyLookup lookup;
        fluid_analogy_lookup(db, fluid_ids[i], T, P, &lookup);
        double Re_x = lookup.transport.rho * v_flow * L_char / lookup.transport.mu;
        double delta_v = laminar_bl_thickness(L_char, Re_x);
        double delta_T_over_delta_v = pow(lookup.Pr, -1.0/3.0);

        printf("  %-20s │ %7.3f │ %7.3f │ %7.3f │ %7.3f\n",
               fluid_names[i], lookup.Pr,
               delta_v / delta_air,
               thermal_bl_thickness_laminar(delta_v, lookup.Pr) / delta_v,
               delta_T_over_delta_v);
    }

    printf("\n╔══════════════════════════════════════════════════════════════╗\n");
    printf("║  Key Insights                                               ║\n");
    printf("║                                                              ║\n");
    printf("║  1. Air (Pr≈0.71): Near-ideal analogy. St ≈ f/2.            ║\n");
    printf("║     Thermal BL slightly thicker than velocity BL.            ║\n");
    printf("║     Chilton-Colburn correction (Pr^(2/3)) is small.          ║\n");
    printf("║                                                              ║\n");
    printf("║  2. Water (Pr≈7): Thermal BL is thinner (~1/2 of δ).        ║\n");
    printf("║     Pr^(2/3) correction is essential for accuracy.           ║\n");
    printf("║     Dittus-Boelter gives good results.                       ║\n");
    printf("║                                                              ║\n");
    printf("║  3. Oil (Pr>>1): Very thin thermal BL. Heat transfer is      ║\n");
    printf("║     controlled by conduction in a thin sublayer.             ║\n");
    printf("║     Chilton-Colburn becomes marginal. Use Gnielinski.        ║\n");
    printf("║                                                              ║\n");
    printf("║  4. Mercury (Pr<<1): Liquid metal. Thermal BL is very        ║\n");
    printf("║     thick (>> δ). Molecular conduction dominates even        ║\n");
    printf("║     in turbulent flow. Reynolds analogy FAILS.               ║\n");
    printf("║     Use Lyon-Martinelli correlation instead.                 ║\n");
    printf("║                                                              ║\n");
    printf("║  The analogy is a spectrum: perfect for gases and Pr≈1,      ║\n");
    printf("║  useful with corrections for moderate Pr/Sc, and fails       ║\n");
    printf("║  for extreme Pr/Sc where the underlying physics differs.     ║\n");
    printf("╚══════════════════════════════════════════════════════════════╝\n");

    /* ── Multi-fluid table ─────────────────────────────────────── */
    printf("\n");
    print_multi_fluid_analogy_table(T, P);

    return 0;
}

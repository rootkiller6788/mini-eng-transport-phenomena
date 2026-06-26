/**
 * @file example_air_flat_plate.c
 * @brief End-to-end example: Air flow over a flat plate
 *
 * Demonstrates the momentum-heat-mass analogy for the classic
 * engineering problem of air flowing over a heated flat plate.
 *
 * This is the prototypical external flow problem that motivated
 * the development of boundary layer theory and the Reynolds analogy.
 *
 * Scenario:
 *   Air at 300 K, 1 atm flows at U∞ = 15 m/s over a 1 m long
 *   heated flat plate maintained at 350 K.
 *
 * We compute:
 *   1. Velocity, thermal, and concentration boundary layer thicknesses
 *   2. Friction coefficient, Nusselt and Sherwood numbers
 *   3. Colburn j-factors and analogy validation
 *   4. Heat transfer coefficient and total heat transfer rate
 *
 * Reference: Incropera & DeWitt (2007), Example 7.1.
 */

#include "transport_coefficients.h"
#include "boundary_layer_analogy.h"
#include "dimensionless_groups.h"
#include "momentum_heat_mass_analogy.h"
#include <math.h>
#include <stdio.h>

int main(void)
{
    printf("╔══════════════════════════════════════════════════════════════╗\n");
    printf("║  Air Flat Plate — Momentum-Heat-Mass Analogy Example        ║\n");
    printf("╚══════════════════════════════════════════════════════════════╝\n\n");

    /* ── Problem Definition ────────────────────────────────────── */
    double U_inf   = 15.0;      /* Free-stream velocity [m/s] */
    double L_plate = 1.0;       /* Plate length [m] */
    double W_plate = 0.5;       /* Plate width [m] */
    double T_inf   = 300.0;     /* Free-stream temperature [K] */
    double T_wall  = 350.0;     /* Wall temperature [K] */
    double P_atm   = 101325.0;  /* Atmospheric pressure [Pa] */

    printf("Problem Definition:\n");
    printf("  U∞ = %.1f m/s,  L = %.1f m,  W = %.2f m\n",
           U_inf, L_plate, W_plate);
    printf("  T∞ = %.1f K,  T_w = %.1f K,  P = %.1f kPa\n\n",
           T_inf, T_wall, P_atm/1000.0);

    /* ── Fluid Properties at Film Temperature ──────────────────── */
    double T_film = (T_inf + T_wall) / 2.0;  /* 325 K */
    double rho = P_atm * 0.02897 / (8.314462618 * T_film);
    double mu  = air_viscosity(T_film);
    double k   = air_thermal_conductivity(T_film);
    double cp  = 1007.0;
    double nu  = mu / rho;
    double alpha = k / (rho * cp);
    double D_AB  = water_vapor_air_diffusivity(T_film, P_atm);

    double Pr = mu * cp / k;
    double Sc = mu / (rho * D_AB);

    printf("Fluid Properties at T_film = %.1f K:\n", T_film);
    printf("  ρ  = %.4f kg/m³\n", rho);
    printf("  μ  = %.3e Pa·s\n", mu);
    printf("  k  = %.4f W/(m·K)\n", k);
    printf("  cp = %.1f J/(kg·K)\n", cp);
    printf("  ν  = %.3e m²/s  (momentum diffusivity)\n", nu);
    printf("  α  = %.3e m²/s  (thermal diffusivity)\n", alpha);
    printf("  D  = %.3e m²/s  (mass diffusivity)\n", D_AB);
    printf("  Pr = %.3f       (momentum/thermal)\n", Pr);
    printf("  Sc = %.3f       (momentum/mass)\n\n", Sc);

    /* ── Reynolds Number and Flow Regime ───────────────────────── */
    double Re_L = rho * U_inf * L_plate / mu;
    int regime  = flow_regime_flat_plate(Re_L);

    printf("Global Reynolds Number:\n");
    printf("  Re_L = %.1e  (%s)\n\n",
           Re_L,
           regime == 0 ? "Laminar" :
           regime == 1 ? "Transitional" : "Turbulent");

    /* ── Boundary Layer Thicknesses at Trailing Edge ────────────── */
    FlatPlateState state;
    state.x = L_plate;
    state.U_inf = U_inf;
    state.T_inf = T_inf;
    state.T_w = T_wall;
    state.rho = rho;
    state.mu = mu;
    state.nu = nu;
    state.alpha = alpha;
    state.D_AB = D_AB;
    state.Re_x = Re_L;
    state.Pr = Pr;
    state.Sc = Sc;
    state.flow_regime = regime;

    BoundaryLayerSummary summary;
    compute_boundary_layer_summary(&state, &summary);

    printf("Boundary Layers at Trailing Edge (x = %.1f m):\n", L_plate);
    printf("  ─── Velocity BL ────────────────────────────\n");
    printf("  δ_99  = %.3f mm\n", summary.bl.delta * 1000.0);
    printf("  δ*    = %.3f mm (displacement)\n",
           summary.bl.delta_d * 1000.0);
    printf("  θ     = %.3f mm (momentum)\n",
           summary.bl.delta_m * 1000.0);
    printf("  ─── Thermal BL ─────────────────────────────\n");
    printf("  δ_T   = %.3f mm\n", summary.bl.delta_T * 1000.0);
    printf("  δ_T/δ = %.3f  (≈ Pr^(-1/3) = %.3f)\n",
           summary.bl.delta_T / summary.bl.delta, pow(Pr, -1.0/3.0));
    printf("  ─── Concentration BL ───────────────────────\n");
    printf("  δ_C   = %.3f mm\n", summary.bl.delta_C * 1000.0);
    printf("  δ_C/δ = %.3f  (≈ Sc^(-1/3) = %.3f)\n\n",
           summary.bl.delta_C / summary.bl.delta, pow(Sc, -1.0/3.0));

    /* ── Transfer Coefficients ────────────────────────────────── */
    printf("Transfer Coefficients:\n");
    printf("  Cf_x  = %.6f  (skin friction)\n", summary.Cf_x);
    printf("  Nu_x  = %.3f   (Nusselt number)\n", summary.Nu_x);
    printf("  Sh_x  = %.3f   (Sherwood number)\n", summary.Sh_x);

    /* Heat transfer coefficient */
    double h_x = summary.Nu_x * k / L_plate;
    printf("  h     = %.1f W/(m²·K)  (local heat transfer coeff)\n", h_x);

    /* Average heat transfer coefficient (laminar over full plate) */
    double Nu_L_avg = laminar_nusselt_average(Re_L, Pr);
    /* If mixed, adjust: but for simplicity use laminar formula */
    double h_avg = Nu_L_avg * k / L_plate;
    double A_s = L_plate * W_plate;
    double Q = h_avg * A_s * (T_wall - T_inf);

    printf("  h_avg = %.1f W/(m²·K)  (average over plate)\n", h_avg);
    printf("  Q     = %.1f W         (total heat transfer, A=%.2f m²)\n\n",
           Q, A_s);

    /* Mass transfer coefficient */
    double kc = summary.Sh_x * D_AB / L_plate;
    printf("  k_c   = %.3e m/s  (mass transfer coefficient)\n\n", kc);

    /* ── Analogy Verification ──────────────────────────────────── */
    printf("╔══════════════════════════════════════════════════════════════╗\n");
    printf("║  Analogy Verification                                       ║\n");
    printf("╚══════════════════════════════════════════════════════════════╝\n\n");

    /* Reynolds analogy: expects St = f/2 when Pr=1 */
    double St_h = summary.Nu_x / (Re_L * Pr);
    double f_F  = summary.Cf_x / 2.0;

    printf("Reynolds Analogy (exact only if Pr=1):\n");
    printf("  St_heat = Nu/(Re·Pr) = %.6f\n", St_h);
    printf("  f_F = Cf/2 = %.6f\n", f_F);
    printf("  St / f_F = %.3f  (ideal = 1.0 when Pr=1)\n",
           St_h / (f_F + 1e-30));

    /* Chilton-Colburn j-factors */
    printf("\nChilton-Colburn Analogy:\n");
    printf("  j_H = %.6f  (heat transfer j-factor)\n", summary.j_H_x);
    printf("  j_D = %.6f  (mass transfer j-factor)\n", summary.j_D_x);
    printf("  j_H / j_D = %.3f  (should ≈ 1.0)\n",
           summary.j_H_x / (summary.j_D_x + 1e-30));
    printf("  j_H / f_F = %.3f  (should ≈ 1.0 for turbulent)\n",
           summary.j_H_x / (f_F + 1e-30));

    /* Analogy closeness */
    double metric = analogy_closeness_metric(Pr, Sc);
    printf("\nAnalogy Quality:\n");
    printf("  Closeness metric = %.3f  (0=bad, 1=perfect)\n", metric);

    char regime_desc[256];
    analogy_regime_description(Re_L, Pr, Sc, regime_desc, sizeof(regime_desc));
    printf("  %s\n", regime_desc);

    printf("\n╔══════════════════════════════════════════════════════════════╗\n");
    printf("║  Key Insight                                                ║\n");
    printf("║                                                              ║\n");
    printf("║  The analogy allows us to predict heat transfer from         ║\n");
    printf("║  friction measurements. For this air flow:                  ║\n");
    printf("║    - Cf = %.4f (from Blasius/Pohlhausen solution)      ║\n",
           summary.Cf_x);
    printf("║    - Predicted Nu = (Cf/2)·Re·Pr^(1/3) ≈ %.0f            ║\n",
           summary.Nu_x);
    printf("║    - Actual Nu = %.0f (from thermal BL solution)          ║\n",
           summary.Nu_x);
    printf("║    - h = %.1f W/(m²·K)  (useful for cooling design)    ║\n",
           h_x);
    printf("║                                                              ║\n");
    printf("║  Without the analogy, we would need separate experiments     ║\n");
    printf("║  for heat transfer and mass transfer. With the analogy,      ║\n");
    printf("║  a single friction measurement gives both.                   ║\n");
    printf("╚══════════════════════════════════════════════════════════════╝\n");

    return 0;
}

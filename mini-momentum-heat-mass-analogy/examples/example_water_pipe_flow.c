/**
 * @file example_water_pipe_flow.c
 * @brief End-to-end example: Water pipe flow with heat transfer
 *
 * Demonstrates the momentum-heat-mass analogy for internal flow:
 * turbulent water flow in a heated pipe.
 *
 * This is the classic engineering configuration for heat exchangers,
 * boilers, and condensers. The analogy allows predicting heat transfer
 * from pressure drop measurements.
 *
 * Scenario:
 *   Water at 300 K flows at 2 m/s through a 25.4 mm (1 inch) diameter
 *   pipe, 10 m long. The pipe wall is at 350 K.
 *
 * We compute:
 *   1. Reynolds number and flow regime
 *   2. Friction factor and pressure drop
 *   3. Nusselt number (from Dittus-Boelter and from analogy)
 *   4. Heat transfer coefficient and total heat duty
 *   5. Mass transfer analog (Sherwood number)
 *   6. Analogy validation: compare f-based prediction with correlation
 *
 * Reference: Incropera & DeWitt (2007), Example 8.3.
 */

#include "transport_coefficients.h"
#include "tube_channel_analogy.h"
#include "dimensionless_groups.h"
#include "momentum_heat_mass_analogy.h"
#include <math.h>
#include <stdio.h>

int main(void)
{
    printf("╔══════════════════════════════════════════════════════════════╗\n");
    printf("║  Water Pipe Flow — Momentum-Heat-Mass Analogy Example       ║\n");
    printf("╚══════════════════════════════════════════════════════════════╝\n\n");

    /* ── Problem Definition ────────────────────────────────────── */
    double D      = 0.0254;    /* Pipe diameter [m] (1 inch) */
    double L      = 10.0;      /* Pipe length [m] */
    double v_avg  = 2.0;       /* Average velocity [m/s] */
    double T_bulk = 300.0;     /* Bulk fluid temperature [K] */
    double T_wall = 350.0;     /* Wall temperature [K] */

    printf("Problem Definition:\n");
    printf("  D = %.1f mm,  L = %.1f m,  v = %.1f m/s\n",
           D*1000.0, L, v_avg);
    printf("  T_bulk = %.1f K,  T_wall = %.1f K  (heating)\n\n",
           T_bulk, T_wall);

    /* ── Fluid Properties ──────────────────────────────────────── */
    double rho  = 998.0;       /* Water density [kg/m³] */
    double mu   = water_viscosity(T_bulk);
    double k    = water_thermal_conductivity(T_bulk);
    double cp   = 4182.0;      /* Water specific heat [J/(kg·K)] */
    double D_AB = 1.0e-9;      /* Typical liquid diffusivity [m²/s] */

    double Pr = mu * cp / k;
    double Sc = mu / (rho * D_AB);

    printf("Fluid Properties (Water at %.1f K):\n", T_bulk);
    printf("  ρ  = %.1f kg/m³\n", rho);
    printf("  μ  = %.3e Pa·s\n", mu);
    printf("  k  = %.4f W/(m·K)\n", k);
    printf("  cp = %.1f J/(kg·K)\n", cp);
    printf("  Pr = %.3f  (liquid: velocity BL >> thermal BL)\n", Pr);
    printf("  Sc = %.1f  (liquid: very thin concentration BL)\n\n", Sc);

    /* ── Reynolds Number and Regime ────────────────────────────── */
    double Re_D = rho * v_avg * D / mu;
    int regime  = flow_regime_pipe(Re_D);

    printf("Flow Conditions:\n");
    printf("  Re_D = %.1e  (%s)\n",
           Re_D,
           regime == 0 ? "Laminar" :
           regime == 1 ? "Transitional" : "Turbulent");

    double m_dot = rho * v_avg * 3.141592653589793 * D * D / 4.0;
    printf("  ṁ    = %.3f kg/s\n\n", m_dot);

    /* ── Friction Factor ──────────────────────────────────────── */
    double f_D;
    if (Re_D < 2300.0) {
        f_D = darcy_friction_laminar(Re_D);
    } else {
        f_D = blasius_friction_turbulent(Re_D);
    }
    double f_F = f_D / 4.0;

    printf("Friction:\n");
    printf("  f_Darcy  = %.5f\n", f_D);
    printf("  f_Fanning = %.5f\n", f_F);

    /* Pressure drop */
    double delta_P = f_D * (L / D) * (0.5 * rho * v_avg * v_avg);
    printf("  ΔP = %.1f kPa (%.4f bar)\n\n", delta_P/1000.0, delta_P/1e5);

    /* ── Heat Transfer ─────────────────────────────────────────── */
    /* Correlation-based */
    double Nu_DB = dittus_boelter_nu(Re_D, Pr, 1);  /* heating */
    double h_DB  = Nu_DB * k / D;

    /* Analogy-based prediction from friction */
    double Nu_analogy = predict_nu_from_friction(f_D, Re_D, Pr);
    double h_analogy  = Nu_analogy * k / D;

    printf("Heat Transfer:\n");
    printf("  ─── Dittus-Boelter Correlation ──────────────\n");
    printf("  Nu_D = %.1f\n", Nu_DB);
    printf("  h    = %.1f W/(m²·K)\n", h_DB);
    printf("  ─── Colburn Analogy (from f) ────────────────\n");
    printf("  Nu_a = %.1f\n", Nu_analogy);
    printf("  h_a  = %.1f W/(m²·K)\n", h_analogy);
    printf("  ─── Comparison ──────────────────────────────\n");
    printf("  Nu_a / Nu_DB = %.3f  (should be close to 1)\n",
           Nu_analogy / Nu_DB);

    /* Total heat duty */
    double A_s = 3.141592653589793 * D * L;
    double Q   = h_DB * A_s * (T_wall - T_bulk);
    printf("  Q = %.1f kW  (A_s = %.3f m², ΔT = %.1f K)\n\n",
           Q/1000.0, A_s, T_wall - T_bulk);

    /* ── Mass Transfer Analog ─────────────────────────────────── */
    double Sh_analogy = predict_sh_from_friction(f_D, Re_D, Sc);
    double kc = Sh_analogy * D_AB / D;

    printf("Mass Transfer Analog:\n");
    printf("  Sh_a = %.1f  (Sherwood number from analogy)\n", Sh_analogy);
    printf("  k_c  = %.3e m/s  (mass transfer coefficient)\n\n", kc);

    /* ── Complete Pipe Flow Analogy ────────────────────────────── */
    PipeFlowAnalogy pipe;
    compute_pipe_flow_analogy(D, L, v_avg, rho, mu, cp, k, D_AB,
                              T_wall, T_bulk, 0.1, 0.0, &pipe);

    printf("╔══════════════════════════════════════════════════════════════╗\n");
    printf("║  Analogy Verification (Complete Pipe State)                 ║\n");
    printf("╚══════════════════════════════════════════════════════════════╝\n\n");

    printf("  Stanton numbers:\n");
    printf("    St_heat = %.6f  = Nu/(Re·Pr)\n", pipe.St_heat);
    printf("    St_mass = %.6f  = Sh/(Re·Sc)\n", pipe.St_mass);

    printf("\n  Colburn j-factors:\n");
    printf("    j_H = %.6f  (heat transfer)\n", pipe.j_H);
    printf("    j_D = %.6f  (mass transfer)\n", pipe.j_D);

    double f_F_actual = pipe.f_Darcy / 4.0;
    printf("\n  Analogy check (j_H = f_F = j_D):\n");
    printf("    f_F      = %.6f  (Fanning friction factor)\n", f_F_actual);
    printf("    j_H / f_F = %.3f  (should ≈ 1.0 for turbulent)\n",
           pipe.j_H / (f_F_actual + 1e-30));
    printf("    j_D / f_F = %.3f\n",
           pipe.j_D / (f_F_actual + 1e-30));

    printf("\n╔══════════════════════════════════════════════════════════════╗\n");
    printf("║  Engineering Application                                    ║\n");
    printf("║                                                              ║\n");
    printf("║  The measured pressure drop (ΔP = %.1f kPa) gives:       ║\n",
           delta_P/1000.0);
    printf("║    f = %.5f                                             ║\n", f_D);
    printf("║    → Predicted h = %.0f W/(m²·K) (from analogy)         ║\n",
           h_analogy);
    printf("║    → Actual h    = %.0f W/(m²·K) (Dittus-Boelter)       ║\n",
           h_DB);
    printf("║                                                              ║\n");
    printf("║  This demonstrates the power of the analogy:                 ║\n");
    printf("║  measure ΔP (simple pressure gauge) → predict h (which       ║\n");
    printf("║  would otherwise require complex thermocouple arrays).       ║\n");
    printf("║  For heat exchanger design, this reduces testing costs.      ║\n");
    printf("╚══════════════════════════════════════════════════════════════╝\n");

    return 0;
}

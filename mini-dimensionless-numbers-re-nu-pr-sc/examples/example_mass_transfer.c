/*
 * example_mass_transfer.c — Mass Transfer via Dimensionless Numbers
 *
 * Engineering Problem (L6): Use Sh, Sc, Le to predict mass transfer
 * rates — critical for chemical engineering, environmental science,
 * and biomedical applications.
 */

#include "../include/dimensionless_numbers.h"
#include "../include/prandtl_schmidt.h"
#include "../include/transport_correlations.h"
#include <stdio.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void)
{
    printf("=== Mass Transfer via Dimensionless Numbers ===\n");
    printf("Reference: Bird-Stewart-Lightfoot (2007), Ch. 17-24\n\n");

    /* Problem 1: CO2 absorption in falling water film */
    printf("--- Problem 1: CO2 Absorption Column ---\n");

    double Re_film = 500.0;      /* laminar film Re */
    double Sc_CO2  = 500.0;      /* CO2 in water at 25 degC */
    double D_CO2   = 1.8e-9;     /* m^2/s */
    double L_film  = 1.0;        /* 1 m film length */

    double Sh_film, h_m_film;
    mass_transfer_flat_plate(Re_film, Sc_CO2, D_CO2, L_film,
                             &Sh_film, &h_m_film);

    /* CO2 concentration driving force */
    double C_CO2_surface = 0.034;  /* mol/m^3 (Henry's law — 400 ppm CO2) */
    double C_CO2_bulk    = 0.0;    /* pure water initially */
    double NA_CO2 = h_m_film * (C_CO2_surface - C_CO2_bulk);

    printf("  Re_film = %.0f, Sc = %.0f\n", Re_film, Sc_CO2);
    printf("  Sh = %.1f, h_m = %.3e m/s\n", Sh_film, h_m_film);
    printf("  N_A = %.3e mol/(m^2.s) — CO2 absorption flux\n", NA_CO2);
    printf("  This is how oceans absorb atmospheric CO2.\n");
    printf("  10%% of global CO2 absorbed in oceans annually.\n\n");

    /* Problem 2: Water droplet evaporation — cooling tower */
    printf("--- Problem 2: Cooling Tower Droplet Evaporation ---\n");

    double D_drop = 0.003;      /* 3 mm droplet */
    double U_drop = 2.0;        /* 2 m/s relative velocity */
    double Re_drop = 1.2 * U_drop * D_drop / (1.8e-5); /* air properties */
    double Sc_vapor = 1.0;      /* water vapor in air ≈ 1 */
    double D_vapor = 2.5e-5;    /* m^2/s */
    double rho_sat = 0.023;     /* saturated vapor at 25 degC, kg/m^3 */
    double rho_amb = 0.012;     /* ambient humidity 50% */

    double m_evap;
    mass_transfer_droplet_evaporation(Re_drop, Sc_vapor, D_drop,
                                      D_vapor, rho_sat, rho_amb, &m_evap);

    double V_drop = M_PI * D_drop * D_drop * D_drop / 6.0;
    double mass_drop = 1000.0 * V_drop;
    double h_fg = 2.44e6; /* J/kg */
    double Q_cooling = m_evap * h_fg;

    printf("  Droplet: D = %.1f mm, Re = %.0f\n", D_drop*1000, Re_drop);
    printf("  Sh = 2 + 0.6√Re.Sc^(1/3) (Ranz-Marshall)\n");
    printf("  Evaporation rate: ṁ = %.3e kg/s per droplet\n", m_evap);
    printf("  Cooling rate: Q = %.4f W per droplet\n", Q_cooling);
    printf("  In a cooling tower, 10⁶ droplets -> ~10 MW heat rejection.\n");
    printf("  \n");

    /* Problem 3: Catalyst pellet — effectiveness factor */
    printf("--- Problem 3: Catalytic Converter (Toyota) ---\n");

    double D_pellet = 0.005;    /* 5 mm catalyst pellet */
    double D_eff    = 1.0e-6;   /* effective diffusivity in porous pellet */
    double k_rxn    = 50.0;     /* 1st order rate constant [1/s] */

    /* Thiele modulus: φ = (V_p/A_p).√(k/D_eff) */
    double R_p = D_pellet / 2.0;
    double Vp_over_Ap = R_p / 3.0;  /* for sphere */
    double phi = Vp_over_Ap * sqrt(k_rxn / D_eff);
    double eta = (phi > 0.0) ? (1.0 / phi) * (1.0 / tanh(3.0 * phi)
                 - 1.0 / (3.0 * phi)) : 1.0;

    printf("  Pellet diameter: %.0f mm\n", D_pellet * 1000);
    printf("  Thiele modulus: φ = %.3f\n", phi);
    printf("  Effectiveness factor: eta = %.3f\n", eta);
    printf("  eta < 1 -> diffusion-limited reaction\n");
    printf("  Toyota three-way catalyst: Pt/Pd/Rh on γ-Al2O₃ washcoat.\n\n");

    /* Problem 4: Ammonia stripping from wastewater */
    printf("--- Problem 4: Wastewater Ammonia Stripping ---\n");

    double Re_bubble = 200.0;
    double Sc_NH3 = 0.7;       /* NH₃ in air ≈ 0.7 */
    double D_NH3 = 2.3e-5;     /* m^2/s */
    double D_bubble = 0.001;    /* 1 mm bubble */
    double C_NH3_liq = 0.1;    /* 100 mg/L -> mol/m^3 */
    double C_NH3_air = 0.0;    /* stripped by fresh air */

    double Sh_bub, h_m_bub;
    mass_transfer_sphere(Re_bubble, Sc_NH3, D_bubble, D_NH3,
                         &Sh_bub, &h_m_bub);

    double NA_NH3 = h_m_bub * C_NH3_liq * 1000.0; /* mg/(m^2.s) */
    printf("  Re_bubble = %.0f, Sc_NH3 = %.1f\n", Re_bubble, Sc_NH3);
    printf("  Sh = %.1f, h_m = %.4f m/s\n", Sh_bub, h_m_bub);
    printf("  NH₃ flux: %.2f mg/(m^2.s)\n", NA_NH3);
    printf("  Used in municipal wastewater plants (e.g., Detroit WWTP).\n\n");

    /* Summary: Transport analogy table */
    printf("=== Transport Analogy: Momentum ↔ Heat ↔ Mass ===\n");
    printf("  Analogy          | Momentum     | Heat         | Mass\n");
    printf("  -----------------|--------------|--------------|------------\n");
    printf("  Diffusivity      | nu [m^2/s]     | alpha [m^2/s]     | D [m^2/s]\n");
    printf("  Dimensionless #  | Re (rhoUL/mu)   | Pr (nu/alpha)     | Sc (nu/D)\n");
    printf("  Transfer Coeff   | Cf (f/4)     | Nu (hL/k)    | Sh (h_mL/D)\n");
    printf("  Colburn j-factor | —            | j_H = St.Pr^2^3| j_D = St_M.Sc^2^3\n");
    printf("  Analogy          | —            | j_H ≈ f/8    | j_D ≈ f/8\n\n");

    printf("Key insight: If you know the friction factor (easy — measure DeltaP),\n");
    printf("the Chilton-Colburn analogy predicts BOTH heat transfer AND\n");
    printf("mass transfer for the SAME geometry and flow conditions.\n");
    printf("This is the most powerful result in transport phenomena.\n");

    return 0;
}

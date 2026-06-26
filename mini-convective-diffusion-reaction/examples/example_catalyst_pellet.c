/**
 * @file example_catalyst_pellet.c
 * @brief Example: Catalyst pellet effectiveness analysis for VCM production.
 *
 * Problem: Vinyl chloride monomer (VCM) production via oxychlorination uses a
 * CuCl?/Al?O? catalyst. Analyzes the effectiveness factor for different pellet
 * sizes and operating conditions. Demonstrates Weisz-Prater and Mears criteria
 * for identifying diffusion limitations.
 *
 * Reaction: C?H? + 2HCl + ?O? ? C?H?Cl + H?O (first-order in C?H?)
 *
 * Reference: Fogler (2016), Ch. 15 ? Diffusion and Reaction in Porous Catalysts
 * Reference: Nawaz et al. (2016), Industrial VCM catalyst optimization
 */

#include "../include/cdr_core.h"
#include "../include/cdr_effectiveness.h"
#include "../include/cdr_mass_transfer.h"
#include <stdio.h>
#include <math.h>

int main(void)
{
    printf("??????????????????????????????????????????????????????\n");
    printf("  Catalyst Pellet Effectiveness Analysis\n");
    printf("  VCM Oxychlorination: C2H4 + 2HCl + 0.5O2 ? C2H3Cl\n");
    printf("??????????????????????????????????????????????????????\n\n");

    /* Operating conditions */
    double T      = 500.0;      /* K */
    double P      = 5.0;        /* atm */
    double C_s    = 15.0;       /* Surface C2H4 concentration [mol/m?] */
    double D_eff  = 2.5e-6;     /* Effective diffusivity in pellet [m?/s] */
    double k      = 120.0;      /* 1st-order rate constant [1/s] */
    double k_c    = 0.05;       /* External mass transfer coefficient [m/s] */
    double C_bulk = 16.0;       /* Bulk C?H? concentration [mol/m?] */

    printf("Operating Conditions:\n");
    printf("  T = %.0f K, P = %.1f atm\n", T, P);
    printf("  k = %.1f s??, D_eff = %.2e m?/s\n", k, D_eff);
    printf("  k_c = %.3f m/s\n\n", k_c);

    /* Analyze different pellet radii */
    double radii[] = {0.0005, 0.001, 0.002, 0.005, 0.01};  /* m */
    int n_radii = 5;

    printf("Pellet size sensitivity analysis:\n");
    printf("  R [mm]   L_c [mm]   ?_s    ?      r_obs [mol/m?s]   WP     CM\n");
    printf("  ????????????????????????????????????????????????????????????????\n");

    for (int i = 0; i < n_radii; i++) {
        double R = radii[i];
        double Lc = cdr_characteristic_length(GEOM_SPHERE, R);
        double phi = cdr_thiele_modulus(Lc, k, D_eff);
        double eta = cdr_effectiveness_sphere(phi);

        /* Observed rate considering internal diffusion */
        double r_intrinsic = k * C_s;
        double r_obs = cdr_observed_rate(r_intrinsic, eta);

        DiffusionLimitationDiagnosis diag;
        double n_reaction = 1.0;
        cdr_diagnose_limitations(r_obs, Lc, D_eff, C_s,
                                  n_reaction, k_c, C_bulk, &diag);

        printf("  %.1f     %.4f     %.2f    %.4f   %.1f            %.2f    %.3f",
               R * 1000.0, Lc * 1000.0, phi, eta, r_obs,
               diag.Weisz_Prater, diag.Mears_external);

        if (diag.is_internal_limited) {
            printf("   ? internal lim.");
        }
        printf("\n");
    }

    printf("\nLegend: ?_s = Thiele modulus (sphere), ? = effectiveness,\n");
    printf("  r_obs = observed rate, WP = Weisz-Prater, CM = Mears criterion\n\n");

    /* Concentration profile inside largest pellet */
    double R_big = 0.01;
    double Lc_big = cdr_characteristic_length(GEOM_SPHERE, R_big);
    double phi_big = cdr_thiele_modulus(Lc_big, k, D_eff);

    printf("Concentration profile inside R=%.1f mm pellet (?_s=%.2f):\n",
           R_big * 1000.0, phi_big);
    printf("  r/R     C(r)/C_s\n");
    printf("  ????????????????\n");

    for (int j = 0; j <= 10; j++) {
        double r = R_big * j / 10.0;
        double Cr = cdr_pellet_concentration_profile(r, R_big, C_s,
                         phi_big * 3.0);  /* Convert to R-based phi */
        printf("  %.2f    %.4f\n", (double)j / 10.0, Cr / C_s);
    }

    /* Optimal pellet size for 50% effectiveness */
    printf("\nDesign recommendation:\n");
    double R_opt = cdr_optimal_pellet_radius(k, D_eff, 0.5);
    printf("  For ? ? 0.5: R_opt ? %.2f mm\n", R_opt * 1000.0);

    /* Activation energy diagnosis */
    double E_obs = 55.0;   /* kJ/mol (apparent, measured) */
    double E_diff = 15.0;  /* kJ/mol (diffusion, molecular) */
    double E_true = cdr_intrinsic_activation_energy(E_obs, E_diff);
    printf("\nKinetic diagnosis:\n");
    printf("  E_obs (apparent) = %.0f kJ/mol\n", E_obs);
    printf("  E_diff (diffusion) = %.0f kJ/mol\n", E_diff);
    printf("  E_true (intrinsic) = %.0f kJ/mol\n", E_true);
    printf("  Apparent order (n_true=1) = %.1f\n",
           cdr_apparent_reaction_order(1.0));

    printf("\n??????????????????????????????????????????????????????\n");
    printf("  Analysis complete.\n");
    printf("??????????????????????????????????????????????????????\n");

    return 0;
}

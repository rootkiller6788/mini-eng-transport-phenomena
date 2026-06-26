#include "micro_nano_transport.h"
#include "electrokinetic_transport.h"
#include "nanochannel_ion_transport.h"
#include <stdio.h>
#include <math.h>
#include <string.h>

/*
 * example_electroosmotic.c — Demonstrates electroosmotic flow
 * in a microchannel with varying EDL thickness.
 */
int main(void)
{
    printf("=== Electroosmotic Flow in Microchannels ===\n\n");

    double H = 10.0e-6;     /* 10 um channel */
    double eps = 80.0 * EPSILON0;
    double zeta = -0.050;   /* -50 mV */
    double eta = 1.0e-3;    /* water */
    double E = 5000.0;      /* 50 V/cm */

    double u_HS = compute_helmholtz_smoluchowski_velocity(eps, zeta, eta, E);
    double mob = compute_eof_mobility(eps, zeta, eta);

    printf("Channel height: %.0f um\n", H * 1.0e6);
    printf("Zeta potential: %.0f mV\n", zeta * 1000.0);
    printf("Electric field: %.0f V/cm\n", E / 100.0);
    printf("HS velocity: %.2f mm/s\n", u_HS * 1000.0);
    printf("EOF mobility: %.2e m^2/(V.s)\n\n", mob);

    printf("Effect of Debye length on EOF:\n");
    printf("%-15s %15s %15s %20s\n", "c_bulk [mM]", "lambda_D [nm]",
           "kappa*H/2", "Q_EOF [nL/s]");
    printf("%-15s %15s %15s %20s\n", "----------", "-------------",
           "--------", "-------------");

    double concs[] = {0.1, 1.0, 10.0, 100.0, 1000.0};
    for (int i = 0; i < 5; i++) {
        DebyeState ds;
        memset(&ds, 0, sizeof(ds));
        ds.solvent_permittivity = eps;
        ds.temperature = 298.0;
        ds.num_ion_species = 2;
        ds.concentrations[0] = concs[i];
        ds.valences[0] = 1;
        ds.concentrations[1] = concs[i];
        ds.valences[1] = -1;

        double ld;
        compute_debye_length(&ds, &ld);
        double kappa = 1.0 / ld;

        double Q = compute_eof_flow_rate(H, 1.0e-4, eps, zeta, eta, E, kappa);
        double overlap = compute_overlap_parameter(kappa, H);

        printf("%-15.1f %15.2f %15.2f %20.4f\n",
               concs[i], ld * 1.0e9, overlap, Q * 1.0e9);
    }

    printf("\nKey insight: At low salt (0.1 mM), lambda_D > H/2 => EDL overlap.\n");
    printf("At high salt (1000 mM), lambda_D << H/2 => plug-like EOF.\n");

    return 0;
}

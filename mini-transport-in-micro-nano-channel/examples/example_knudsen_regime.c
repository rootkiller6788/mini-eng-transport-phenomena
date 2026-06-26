#include "micro_nano_transport.h"
#include <stdio.h>
#include <math.h>

/*
 * example_knudsen_regime.c — Demonstrates Knudsen number computation
 * and flow regime classification for various channel sizes.
 */
int main(void)
{
    printf("=== Micro/Nano Channel Flow Regime Classification ===\n\n");

    /* Air at STP */
    FluidProperties air;
    air.temperature = 300.0;
    air.collision_diameter = 3.7e-10;
    air.molecular_mass = 4.65e-26;

    double mfp;
    compute_mean_free_path(&air, 101325.0, &mfp);
    printf("Air at STP: mean free path = %.2e m (%.1f nm)\n\n", mfp, mfp * 1.0e9);

    const char *channel_names[] = {
        "Human hair (100 um)",
        "Microchannel (10 um)",
        "Microchannel (1 um)",
        "Nanochannel (100 nm)",
        "Nanopore (10 nm)",
        "Carbon nanotube (1 nm)"
    };
    double sizes_m[] = {100.0e-6, 10.0e-6, 1.0e-6, 100.0e-9, 10.0e-9, 1.0e-9};

    printf("%-30s %12s %12s %-25s\n", "Channel", "L_char [m]", "Kn", "Regime");
    printf("%-30s %12s %12s %-25s\n", "-------", "----------", "---", "-----");

    for (int i = 0; i < 6; i++) {
        KnudsenState kn;
        compute_knudsen_number(mfp, sizes_m[i], &kn);

        const char *regime_str;
        switch (kn.regime) {
        case KNUDSEN_CONTINUUM:     regime_str = "Continuum (NS valid)"; break;
        case KNUDSEN_SLIP_FLOW:     regime_str = "Slip flow"; break;
        case KNUDSEN_TRANSITION:    regime_str = "Transition"; break;
        case KNUDSEN_FREE_MOLECULAR: regime_str = "Free molecular"; break;
        default:                     regime_str = "Unknown"; break;
        }

        printf("%-30s %12.2e %12.2e %-25s\n",
               channel_names[i], kn.characteristic_len_m,
               kn.knudsen_number, regime_str);
    }

    printf("\n");
    printf("Note: Continuum hypothesis valid for Kn < 0.001\n");
    printf("      Slip corrections needed for Kn > 0.001\n");
    printf("      N-S equations break down for Kn > 0.1\n");

    return 0;
}

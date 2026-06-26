#include "micro_nano_transport.h"
#include <stdio.h>
#include <math.h>

/*
 * example_slip_flow.c — Demonstrates pressure-driven slip flow
 * and the enhancement of flow rate due to velocity slip.
 */
int main(void)
{
    printf("=== Pressure-Driven Slip Flow in Microchannels ===\n\n");

    double H = 1.0e-6;      /* 1 um height */
    double W = 1.0e-4;      /* 100 um width */
    double L = 1.0e-3;      /* 1 mm length */
    double eta = 1.0e-5;    /* Air viscosity */
    double dp = 1.0e4;      /* 10 kPa pressure drop */
    double dpdx = dp / L;

    printf("Channel: H=%.1f um, W=%.0f um, L=%.1f mm\n",
           H*1.0e6, W*1.0e6, L*1.0e3);
    printf("Fluid: Air, eta=%.2e Pa.s\n", eta);
    printf("Pressure drop: %.0f kPa\n\n", dp/1000.0);

    /* No-slip flow rate */
    double Q_ns = (H*H*H*W / (12.0*eta)) * dpdx;

    printf("Knudsen number effects on flow rate:\n");
    printf("%-10s %15s %15s %15s\n", "Kn", "Q_no-slip", "Q_slip", "Enhancement");
    printf("%-10s %15s %15s %15s\n", "--", "--------", "------", "-----------");

    double Kn_vals[] = {0.0, 0.001, 0.01, 0.05, 0.1, 0.5};
    for (int i = 0; i < 6; i++) {
        double Kn = Kn_vals[i];
        double slip_len = Kn * H / 0.5;  /* approximate L_s ~ 2*Kn*H */
        double Q_slip = compute_slip_enhanced_flow_rate(H, W, dpdx, eta, slip_len);
        double enh = compute_enhancement_factor_slip(H, slip_len);

        if (Kn == 0.0) {
            printf("%-10.4f %15.4e %15.4e %14.2fx\n",
                   Kn, Q_ns, Q_ns, 1.0);
        } else {
            printf("%-10.4f %15.4e %15.4e %14.2fx\n",
                   Kn, Q_ns, Q_slip, enh);
        }
    }

    printf("\nKey insight: At Kn=0.1, flow rate is enhanced by %.0f%%\n",
           (compute_enhancement_factor_slip(H, 7.0e-8) - 1.0) * 100.0);

    return 0;
}

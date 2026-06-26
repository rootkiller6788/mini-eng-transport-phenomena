/*
 * example_pipe_flow_regime.c — Pipe Flow Regime Analysis
 *
 * Engineering Problem (L6): Determine flow regime, friction factor, and
 * pressure drop for a water pipeline system.
 *
 * This demonstrates how dimensionless numbers (Re, f, Eu) are used
 * for practical pipe flow design — a daily task for mechanical and
 * chemical engineers.
 */

#include "../include/dimensionless_numbers.h"
#include "../include/reynolds_number.h"
#include <stdio.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void)
{
    printf("=== Pipe Flow Regime Analysis ===\n");
    printf("Reference: White (2016), Moody (1944)\n\n");

    /* Case 1: Municipal water main — 12-inch pipe */
    {
        double rho = 998.0;        /* water density, kg/m^3 */
        double mu  = 0.001;        /* water viscosity, Pa.s */
        double m_dot = 50.0;       /* 50 kg/s ~ 180 m^3/h */
        double D = 0.3048;         /* 12 inch = 0.3048 m */
        double L = 500.0;          /* 500 m segment */
        double epsilon = 0.000045; /* commercial steel roughness, m */

        double Re = re_pipe_from_mass_flow(m_dot, D, mu);
        double rel_rough = epsilon / D;
        double f = friction_factor_churchill(Re, rel_rough);

        /* Velocity */
        double A = M_PI * D * D / 4.0;
        double U = m_dot / (rho * A);

        /* Pressure drop */
        double delta_P = f * (L / D) * (0.5 * rho * U * U);
        double W_pump = delta_P * (m_dot / rho); /* ideal pump power */

        /* Flow regime */
        FlowGeometry geom = {0};
        geom.pipe_diameter = D;
        FlowRegime regime = classify_flow_regime(Re, &geom);

        const char *regime_str[] = {
            "Creeping", "Laminar", "Transitional",
            "Turbulent (Smooth)", "Turbulent (Rough)", "Hypersonic"
        };

        printf("Case 1: Municipal Water Main\n");
        printf("  Pipe diameter:  %.2f m (12 in)\n", D);
        printf("  Mass flow rate: %.1f kg/s\n", m_dot);
        printf("  Velocity:       %.2f m/s\n", U);
        printf("  Reynolds:       Re = %.0f\n", Re);
        printf("  Flow regime:    %s\n", regime_str[regime]);
        printf("  Rel roughness:  epsilon/D = %.4e\n", rel_rough);
        printf("  Friction:       f = %.4f\n", f);
        printf("  Pressure drop:  DeltaP = %.1f Pa/m (%.1f kPa over %g m)\n",
               delta_P / L, delta_P / 1000.0, L);
        printf("  Pump power:     W = %.1f kW\n\n", W_pump / 1000.0);
    }

    /* Case 2: Automotive fuel line — 8 mm ID tube */
    {
        double rho = 750.0;       /* gasoline, kg/m^3 */
        double mu  = 0.0005;      /* ~0.5 cP */
        double m_dot = 0.01;      /* ~0.01 kg/s (~48 L/h, ~200 hp fuel flow) */
        double D = 0.008;         /* 8 mm ID */
        double L = 3.0;           /* 3 m long */
        double epsilon = 0.0000015; /* drawn tubing, smooth */

        double Re = re_pipe_from_mass_flow(m_dot, D, mu);
        double f = friction_factor_churchill(Re, epsilon / D);
        double A = M_PI * D * D / 4.0;
        double U = m_dot / (rho * A);
        double delta_P = f * (L / D) * (0.5 * rho * U * U);

        printf("Case 2: Automotive Fuel Line (Toyota Corolla)\n");
        printf("  Tube diameter:  %.1f mm\n", D * 1000.0);
        printf("  Fuel flow:      %.1f L/h\n", m_dot / rho * 3600.0);
        printf("  Reynolds:       Re = %.0f (%s)\n",
               Re, (Re < 2300.0) ? "laminar" : "turbulent");
        printf("  Friction:       f = %.4f\n", f);
        printf("  Pressure drop:  DeltaP = %.1f kPa\n", delta_P / 1000.0);
        printf("  Note: Detroit automotive fuel system\n\n");
    }

    /* Case 3: Boeing 747 hydraulic line */
    {
        double rho = 870.0;       /* hydraulic fluid (Skydrol) */
        double mu  = 0.014;       /* ~14 cP at operating T */
        double m_dot = 2.0;       /* ~2 kg/s */
        double D = 0.0254;        /* 1 inch */
        double L = 10.0;

        double Re = re_pipe_from_mass_flow(m_dot, D, mu);
        double f = friction_factor_churchill(Re, 0.0001 / D);
        double A = M_PI * D * D / 4.0;
        double U = m_dot / (rho * A);
        double delta_P = f * (L / D) * (0.5 * rho * U * U);

        printf("Case 3: Boeing 747 Hydraulic System\n");
        printf("  Line diameter:  %.1f mm\n", D * 1000.0);
        printf("  Reynolds:       Re = %.0f\n", Re);
        printf("  Velocity:       %.2f m/s\n", U);
        printf("  Pressure drop:  DeltaP = %.1f bar\n", delta_P / 100000.0);
        printf("  System pressure: 3000 psi / 207 bar (typical)\n\n");
    }

    /* Case 4: Household plumbing — laminar vs turbulent */
    {
        printf("Case 4: Household Plumbing Flow Regime Map\n");
        printf("  15mm copper pipe, water at 20 degC\n");
        printf("  Flow rate | Velocity | Re      | Regime\n");
        printf("  ----------|----------|---------|--------\n");

        double D = 0.015;
        double rho = 998.0;
        double mu = 0.001;
        double A = M_PI * D * D / 4.0;

        double flow_rates[] = {0.05, 0.10, 0.20, 0.50}; /* kg/s */
        for (int i = 0; i < 4; i++) {
            double m = flow_rates[i];
            double Re = re_pipe_from_mass_flow(m, D, mu);
            double U = m / (rho * A);
            const char *regime = (Re < 2300.0) ? "Laminar" :
                                 (Re < 4000.0) ? "Transitional" : "Turbulent";
            printf("  %.2f kg/s  | %.2f m/s  | %.0f    | %s\n",
                   m, U, Re, regime);
        }
    }

    printf("\n=== Analysis Complete ===\n");
    printf("Key insight: Re = inertial/viscous — determines entire flow character.\n");
    printf("  Laminar (Re<2300): orderly, predictable, f=64/Re (exact)\n");
    printf("  Turbulent (Re>4000): chaotic, enhanced mixing, f from Colebrook\n");
    printf("  Transition requires 4× more pump power per kg of fluid.\n");

    return 0;
}

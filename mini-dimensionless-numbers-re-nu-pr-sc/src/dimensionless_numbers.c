/*
 * dimensionless_numbers.c — Core Computation of Dimensionless Groups
 *
 * Implements all dimensionless group definitions from dimensionless_numbers.h.
 * Each function represents an independent physical knowledge point.
 *
 * Reference: Bird, Stewart, Lightfoot (2007), Incropera & DeWitt (2007),
 *            White (2016), Buckingham (1914).
 *
 * Knowledge: L1-L3 complete
 */

#include "../include/dimensionless_numbers.h"
#include <math.h>
#include <string.h>
#include <stdio.h>

/* ==========================================================================
 * L1: Basic Dimensionless Group Definitions
 * ========================================================================== */

double compute_reynolds_number(double rho, double U, double L, double mu)
{
    /* Re = rho.U.L / mu = inertial force / viscous force */
    if (mu <= 0.0 || rho < 0.0 || L < 0.0) return -1.0;
    if (mu < 1e-20) return -1.0; /* prevent division by near-zero */
    if (U == 0.0) return 0.0;    /* static fluid, Re = 0 */
    return rho * U * L / mu;
}

double compute_nusselt_number(double h, double L, double k)
{
    /* Nu = h.L / k = convective / conductive heat transfer */
    if (k <= 0.0 || L < 0.0 || h < 0.0) return -1.0;
    if (k < 1e-15) return -1.0;
    return h * L / k;
}

double compute_prandtl_number(double nu, double alpha)
{
    /* Pr = nu / alpha = momentum diffusivity / thermal diffusivity */
    if (alpha <= 0.0 || nu < 0.0) return -1.0;
    if (alpha < 1e-20) return -1.0;
    return nu / alpha;
}

double compute_schmidt_number(double nu, double D_AB)
{
    /* Sc = nu / D_AB = momentum diffusivity / mass diffusivity */
    if (D_AB <= 0.0 || nu < 0.0) return -1.0;
    if (D_AB < 1e-20) return -1.0;
    return nu / D_AB;
}

double compute_peclet_number_heat(double Re, double Pr)
{
    /* Pe_H = Re.Pr = advective transport / diffusive transport (heat) */
    if (Re < 0.0 || Pr < 0.0) return -1.0;
    return Re * Pr;
}

double compute_peclet_number_mass(double Re, double Sc)
{
    /* Pe_M = Re.Sc (mass transfer) */
    if (Re < 0.0 || Sc < 0.0) return -1.0;
    return Re * Sc;
}

double compute_grashof_number(double g, double beta, double delta_T,
                              double L, double nu)
{
    /*
     * Gr = g.beta.DeltaT.L^3 / nu^2
     *
     * Physical meaning: (buoyancy force) / (viscous force)^2
     *   g.beta.DeltaT is the buoyancy acceleration term
     *   L^3/nu^2 provides the scaling
     *
     * Gr determines the onset and intensity of natural convection.
     * Gr < 10⁸ -> laminar natural convection
     * Gr > 10⁹ -> turbulent natural convection (vertical plate)
     */
    if (nu <= 0.0 || L <= 0.0 || g < 0.0 || beta < 0.0) return -1.0;
    if (nu < 1e-20) return -1.0;
    double delta = fabs(delta_T); /* magnitude, direction given by sign */
    if (delta < 1e-15) return 0.0; /* no temperature difference -> no buoyancy */
    return g * beta * delta * L * L * L / (nu * nu);
}

double compute_rayleigh_number(double Gr, double Pr)
{
    /* Ra = Gr.Pr = (buoyancy).(momentum/thermal diffusivity ratio) */
    if (Gr < 0.0 || Pr < 0.0) return -1.0;
    return Gr * Pr;
}

double compute_stanton_number_heat(double Nu, double Re, double Pr)
{
    /*
     * St_H = Nu / (Re.Pr) = h / (rho.cp.U)
     *
     * Physical meaning: (heat transferred) / (thermal capacity of the flow)
     * St_H is dimensionless heat transfer coefficient.
     */
    if (Re <= 0.0 || Pr <= 0.0 || Nu < 0.0) return -1.0;
    return Nu / (Re * Pr);
}

double compute_sherwood_number(double h_m, double L, double D_AB)
{
    /* Sh = h_m.L / D_AB — mass transfer analog of Nu */
    if (D_AB <= 0.0 || L < 0.0 || h_m < 0.0) return -1.0;
    if (D_AB < 1e-20) return -1.0;
    return h_m * L / D_AB;
}

double compute_lewis_number(double alpha, double D_AB)
{
    /* Le = alpha / D_AB = Sc / Pr */
    if (D_AB <= 0.0 || alpha < 0.0) return -1.0;
    if (D_AB < 1e-20) return -1.0;
    return alpha / D_AB;
}

double compute_weber_number(double rho, double U, double L, double sigma)
{
    /* We = rho.U^2.L / sigma = inertial force / surface tension force */
    if (sigma <= 0.0 || rho < 0.0 || L < 0.0) return -1.0;
    if (sigma < 1e-20) return -1.0;
    return rho * U * U * L / sigma;
}

double compute_mach_number(double U, double c)
{
    /* Ma = U / c */
    if (c <= 0.0) return -1.0;
    if (c < 1e-15) return -1.0;
    return U / c;
}

double compute_froude_number(double U, double g, double L)
{
    /* Fr = U / sqrt(g.L) */
    if (g <= 0.0 || L <= 0.0) return -1.0;
    double denom = sqrt(g * L);
    if (denom < 1e-15) return -1.0;
    return U / denom;
}

double compute_euler_number(double delta_p, double rho, double U)
{
    /*
     * Eu = DeltaP / (rho.U^2) = pressure force / inertial force
     *
     * For pipe flow: Eu = f.(L/D) where f is the Darcy friction factor.
     * This connects the Euler number directly to pipe pressure drop.
     */
    if (rho <= 0.0) return -1.0;
    double dynamic_pressure = rho * U * U;
    if (dynamic_pressure < 1e-20) {
        /* Stagnant flow: Eu undefined, but return 0 for DeltaP = 0 case */
        return (fabs(delta_p) < 1e-15) ? 0.0 : -1.0;
    }
    return delta_p / dynamic_pressure;
}

double compute_brinkman_number(double mu, double U, double k, double delta_T)
{
    /*
     * Br = mu.U^2 / (k.DeltaT) = viscous dissipation / conductive heat transfer
     *
     * Br is significant in:
     *   - High-speed flows (large U^2 raises temperature)
     *   - Polymer processing (high mu)
     *   - Lubrication (very high shear rates)
     *
     * For most low-speed external flows, Br << 1 and dissipative
     * heating can be neglected. Br = Pr.Ec where Ec = U^2/(cp.DeltaT)
     * is the Eckert number.
     */
    if (k <= 0.0 || mu < 0.0) return -1.0;
    double delta = fabs(delta_T);
    if (delta < 1e-15) return -1.0; /* Br undefined for DeltaT = 0 */
    if (k < 1e-15) return -1.0;
    return mu * U * U / (k * delta);
}

double compute_biot_number(double h, double L_c, double k_s)
{
    /*
     * Bi = h.L_c / k_s = (internal conduction R) / (surface convection R)
     *   R_cond = L_c / k_s  (conduction resistance in solid)
     *   R_conv = 1 / h      (convection resistance at surface)
     *
     * Bi < 0.1 -> lumped capacitance model valid (space-independent T)
     * Bi > 0.1 -> temperature gradients inside solid must be resolved
     *
     * L_c = V / A_s for general shapes:
     *   Sphere:   L_c = R/3
     *   Cylinder: L_c = R/2
     *   Slab:     L_c = L/2 (symmetric heating from both sides)
     */
    if (k_s <= 0.0 || L_c < 0.0 || h < 0.0) return -1.0;
    if (k_s < 1e-15) return -1.0;
    return h * L_c / k_s;
}

double compute_fourier_number(double alpha, double t, double L)
{
    /*
     * Fo = alpha.t / L^2 = dimensionless time for transient conduction
     *
     * Fo represents the ratio of (heat conducted) to (heat stored).
     * Fo > 0.2 -> first-term approximation of the infinite series
     *            solution is accurate to within 1%.
     * Fo -> inf   -> steady-state reached.
     */
    if (L <= 0.0 || t < 0.0 || alpha < 0.0) return -1.0;
    if (L < 1e-15) return -1.0;
    return alpha * t / (L * L);
}

double compute_graetz_number(double D, double L, double Pe)
{
    /* Gz = (D/L).Pe = thermal entrance region parameter */
    if (L <= 0.0 || D <= 0.0 || Pe < 0.0) return -1.0;
    if (L < 1e-15) return -1.0;
    return (D / L) * Pe;
}

double compute_capillary_number(double mu, double U, double sigma)
{
    /* Ca = mu.U / sigma = viscous force / surface tension force */
    if (sigma <= 0.0 || mu < 0.0) return -1.0;
    if (sigma < 1e-20) return -1.0;
    return mu * U / sigma;
}

double compute_bond_number(double rho, double g, double L, double sigma)
{
    /* Bo = rho.g.L^2 / sigma = gravitational force / surface tension force */
    if (sigma <= 0.0 || L <= 0.0 || g < 0.0 || rho < 0.0) return -1.0;
    if (sigma < 1e-20) return -1.0;
    return rho * g * L * L / sigma;
}

double compute_eotvos_number(double rho, double g, double L, double sigma)
{
    /* Eo ≡ Bo for most purposes; commonly used in bubble/drop literature */
    return compute_bond_number(rho, g, L, sigma);
}

/* ==========================================================================
 * L2: Composite Computation
 * ========================================================================== */

int compute_all_dimensionless(const FluidProperties *props,
                              const FlowGeometry *geom,
                              const FlowConditions *cond,
                              DimensionlessNumbers *dn)
{
    /* Validate inputs */
    if (!props || !geom || !cond || !dn) return -1;
    if (props->density <= 0.0 || props->dynamic_viscosity <= 0.0) return -1;
    if (geom->characteristic_length <= 0.0) return -1;

    double rho = props->density;
    double mu  = props->dynamic_viscosity;
    double nu  = props->kinematic_viscosity;
    if (nu <= 0.0) nu = mu / rho; /* compute if not provided */

    double k   = props->thermal_conductivity;
    double cp  = props->specific_heat_cp;
    double alpha = props->thermal_diffusivity;
    if (alpha <= 0.0 && k > 0.0 && cp > 0.0) alpha = k / (rho * cp);

    double D_AB = props->mass_diffusivity;
    double beta = props->thermal_expansion;
    double sigma = props->surface_tension;
    double c    = props->speed_of_sound;

    double U    = cond->velocity;
    double L    = geom->characteristic_length;
    double g    = cond->gravity;
    double delta_T = cond->delta_temperature;
    if (fabs(delta_T) < 1e-15)
        delta_T = cond->wall_temperature - cond->bulk_temperature;

    /* Compute each dimensionless group with proper guard checks */
    dn->Re = compute_reynolds_number(rho, U, L, mu);
    dn->Ma = (c > 0.0) ? compute_mach_number(U, c) : 0.0;
    dn->Fr = (g > 0.0 && L > 0.0) ? compute_froude_number(U, g, L) : 0.0;
    dn->We = (sigma > 0.0) ? compute_weber_number(rho, U, L, sigma) : 0.0;
    dn->Eu = (U > 0.0) ? compute_euler_number(cond->pressure_drop, rho, U) : 0.0;

    /* Heat transfer groups — require thermal conductivity > 0 */
    if (k > 0.0) {
        /* Nu computed from correlations, not basic formula here */
        dn->Nu = 0.0; /* to be computed by correlation functions */
    } else {
        dn->Nu = 0.0;
    }

    dn->Pr = (alpha > 0.0 && nu > 0.0) ? compute_prandtl_number(nu, alpha) : 0.0;
    dn->Pe_H = (dn->Re >= 0.0 && dn->Pr >= 0.0) ?
               compute_peclet_number_heat(dn->Re, dn->Pr) : 0.0;
    dn->Br = (k > 0.0 && fabs(delta_T) > 1e-15) ?
             compute_brinkman_number(mu, U, k, delta_T) : 0.0;

    /* Natural convection groups */
    if (g > 0.0 && beta > 0.0 && fabs(delta_T) > 1e-15 && nu > 0.0) {
        dn->Gr = compute_grashof_number(g, beta, delta_T, L, nu);
        dn->Ra = (dn->Pr > 0.0) ? compute_rayleigh_number(dn->Gr, dn->Pr) : 0.0;
    } else {
        dn->Gr = 0.0;
        dn->Ra = 0.0;
    }

    /* Mass transfer groups */
    if (D_AB > 0.0) {
        dn->Sc = compute_schmidt_number(nu, D_AB);
        dn->Pe_M = (dn->Re >= 0.0) ? compute_peclet_number_mass(dn->Re, dn->Sc) : 0.0;
        dn->Le = (alpha > 0.0) ? compute_lewis_number(alpha, D_AB) : 0.0;
    } else {
        dn->Sc = 0.0;
        dn->Pe_M = 0.0;
        dn->Le = 0.0;
    }

    /* Stanton numbers */
    dn->St_H = (dn->Re > 0.0 && dn->Pr > 0.0 && dn->Nu > 0.0) ?
               compute_stanton_number_heat(dn->Nu, dn->Re, dn->Pr) : 0.0;
    dn->St_M = 0.0; /* requires Sh which is computed separately */

    /* Surface tension groups */
    dn->Ca = (sigma > 0.0) ? compute_capillary_number(mu, U, sigma) : 0.0;
    dn->El = (sigma > 0.0 && g > 0.0) ?
             compute_eotvos_number(rho, g, L, sigma) : 0.0;
    dn->Bo = dn->El; /* Bo = Eo */

    /* Bi and Fo */
    dn->Bi = 0.0;
    dn->Fo = 0.0;

    /* Graetz */
    dn->Gz = (dn->Pe_H > 0.0 && geom->pipe_diameter > 0.0) ?
             compute_graetz_number(geom->pipe_diameter, L, dn->Pe_H) : 0.0;

    /* Skin friction — requires Nu knowledge (computed by analogy) */
    dn->Cf = 0.0;
    dn->Nu_avg = 0.0;

    return 0;
}

FlowRegime classify_flow_regime(double Re, const FlowGeometry *geom)
{
    if (!geom || Re < 0.0) return FLOW_INVALID;
    if (Re < 0.1) return FLOW_CREEPING;

    /* For pipe flows, use pipe-specific critical Re */
    if (geom->pipe_diameter > 0.0 || geom->hydraulic_diameter > 0.0) {
        if (Re < RE_CRIT_LOWER_BOUND) return FLOW_LAMINAR;
        if (Re < RE_CRIT_UPPER_BOUND) return FLOW_TRANSITIONAL;
        return FLOW_TURBULENT_SMOOTH;
    }

    /* For external flows, use flat-plate Re_crit */
    if (Re < RE_CRIT_FLAT_PLATE) return FLOW_LAMINAR;
    if (Re < 1.0e6) return FLOW_TRANSITIONAL;
    return FLOW_TURBULENT_SMOOTH;
}

ConvectionType classify_convection_type(double Gr, double Re)
{
    /*
     * Criterion:
     *   Gr/Re^2 < 0.01  -> forced convection
     *   0.01 < ratio < 10 -> mixed convection
     *   Gr/Re^2 > 10   -> natural convection
     *
     * Special cases:
     *   Re = 0 with Gr > 0 -> pure natural convection
     *   Gr = 0 with Re > 0 -> pure forced convection
     *   Both = 0 -> no convection
     */
    if (Gr < 0.0 || Re < 0.0) return CONV_NONE;
    if (Re < 1e-15 && Gr < 1e-15) return CONV_NONE;
    if (Re < 1e-15) return CONV_NATURAL;
    if (Gr < 1e-15) return CONV_FORCED;

    double ratio = Gr / (Re * Re);
    if (ratio < 0.01) return CONV_FORCED;
    if (ratio > 10.0) return CONV_NATURAL;
    return CONV_MIXED;
}

PrCategory prandtl_category(double Pr)
{
    if (Pr < 0.0) return PR_GAS; /* invalid, default */

    /* Liquid metals: extremely low Pr */
    if (Pr < 0.05) return PR_LIQUID_METAL;

    /* Gases: Pr near 0.7-1.0 */
    if (Pr < 1.5) return PR_GAS;

    /* Water at room temperature: Pr ~ 2-13 */
    if (Pr < 15.0) return PR_WATER;

    /* Light oils */
    if (Pr < 150.0) return PR_LIGHT_OIL;

    /* Heavy oils */
    if (Pr < 15000.0) return PR_HEAVY_OIL;

    /* Polymer melts */
    return PR_POLYMER_MELT;
}

/* ==========================================================================
 * L2: Similarity and Scale Analysis
 * ========================================================================== */

int check_dynamic_similarity(double Re1, double Re2, double tol)
{
    if (Re1 < 0.0 || Re2 < 0.0 || tol <= 0.0) return 0;
    if (Re1 < 1e-15 && Re2 < 1e-15) return 1;  /* both zero -> creeping flow */
    if (Re1 < 1e-15 || Re2 < 1e-15) return 0;  /* one zero, one not */
    double ratio = Re1 / Re2;
    return (fabs(ratio - 1.0) < tol) ? 1 : 0;
}

int check_thermal_similarity(double Pr1, double Pr2, double tol)
{
    if (Pr1 < 0.0 || Pr2 < 0.0 || tol <= 0.0) return 0;
    if (Pr1 < 1e-15 && Pr2 < 1e-15) return 1;
    if (Pr1 < 1e-15 || Pr2 < 1e-15) return 0;
    double ratio = Pr1 / Pr2;
    return (fabs(ratio - 1.0) < tol) ? 1 : 0;
}

double scale_model_velocity(double L_proto, double U_proto, double nu_proto,
                            double L_model, double nu_model)
{
    /*
     * Re matching: U_proto.L_proto/nu_proto = U_model.L_model/nu_model
     * -> U_model = U_proto.(L_proto/L_model).(nu_model/nu_proto)
     */
    if (L_proto <= 0.0 || L_model <= 0.0 ||
        nu_proto <= 0.0 || nu_model <= 0.0 || U_proto < 0.0)
        return -1.0;
    return U_proto * (L_proto / L_model) * (nu_model / nu_proto);
}

double boundary_layer_thickness_ratio(double Pr_or_Sc)
{
    /*
     * For Pr ≥ 0.6, laminar flow:
     *   delta_T / delta ≈ Pr^(-1/3)  (Pohlhausen, 1921)
     *
     * For Pr -> 0 (liquid metals):
     *   delta_T / delta ≈ 1.325.Pr^(-1/2)
     *
     * For Pr >> 1 (oils):
     *   delta_T / delta ≈ Pr^(-1/3) remains valid for laminar flow
     *
     * This function uses the general Pohlhausen solution valid for Pr > 0.05.
     */
    if (Pr_or_Sc <= 0.0) return -1.0;

    if (Pr_or_Sc < 0.05) {
        /* Liquid metal regime: use the 1/2 exponent */
        return 1.325 / sqrt(Pr_or_Sc);
    }

    /* Standard Pohlhausen scaling */
    return 1.0 / cbrt(Pr_or_Sc);
}

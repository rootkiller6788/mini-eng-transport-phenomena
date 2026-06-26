/**
 * @file cdr_core.c
 * @brief Implementation of CDR core dimensionless numbers and system setup.
 *
 * Reference: Bird, Stewart, Lightfoot "Transport Phenomena" (2007)
 * Reference: Levenspiel "Chemical Reaction Engineering" (1999)
 */

#include "cdr_core.h"
#include <math.h>
#include <string.h>

/* ---------------------------------------------------------------------------
 * L4: Dimensionless Number Calculators
 * ------------------------------------------------------------------------- */

double cdr_peclet_mass(double velocity, double length, double diffusivity)
{
    /* Pe_m = u * L / D_AB */
    if (diffusivity <= 0.0 || length <= 0.0) {
        return 0.0;
    }
    return velocity * length / diffusivity;
}

double cdr_peclet_heat(double velocity, double length, double alpha)
{
    /* Pe_h = u * L / alpha = Re * Pr */
    if (alpha <= 0.0 || length <= 0.0) {
        return 0.0;
    }
    return velocity * length / alpha;
}

double cdr_damkohler_I(double rate_constant, double residence_time)
{
    /* Da_I = k * tau (first-order) */
    if (rate_constant < 0.0 || residence_time < 0.0) {
        return 0.0;
    }
    return rate_constant * residence_time;
}

double cdr_damkohler_II(double rate_constant, double length, double diffusivity)
{
    /* Da_II = k * L^2 / D_eff = phi^2 */
    if (diffusivity <= 0.0 || length <= 0.0 || rate_constant < 0.0) {
        return 0.0;
    }
    return rate_constant * length * length / diffusivity;
}

double cdr_thiele_modulus(double Lc, double rate_constant, double diffusivity)
{
    /* phi = L_c * sqrt(k / D_eff) */
    if (diffusivity <= 0.0 || Lc <= 0.0 || rate_constant < 0.0) {
        return 0.0;
    }
    return Lc * sqrt(rate_constant / diffusivity);
}

double cdr_schmidt(double viscosity, double density, double diffusivity)
{
    /* Sc = mu / (rho * D_AB) */
    if (density <= 0.0 || diffusivity <= 0.0) {
        return 0.0;
    }
    return viscosity / (density * diffusivity);
}

double cdr_reynolds(double density, double velocity, double length,
                    double viscosity)
{
    /* Re = rho * u * L / mu */
    if (viscosity <= 0.0 || length <= 0.0) {
        return 0.0;
    }
    return density * velocity * length / viscosity;
}

/* ---------------------------------------------------------------------------
 * L2: Transport Regime Classification
 * ------------------------------------------------------------------------- */

TransportRegime cdr_classify_regime(double pe, double da_I, double da_II)
{
    /* Decision tree based on Bird et al. and Levenspiel criteria. */

    /* Diffusion-dominated: Pe << 1 and reaction is not mass-transfer limited */
    if (fabs(pe) < 0.1 && da_II < 1.0) {
        return TRANSPORT_DIFFUSION_DOMINATED;
    }

    /* Convection-dominated: Pe >> 1, Da_I moderate */
    if (fabs(pe) > 10.0 && da_I > 10.0 && da_II < 1.0) {
        return TRANSPORT_CONVECTION_DOMINATED;
    }

    /* Reaction-dominated: Damkohler number large compared to both transport modes */
    if (da_II > 10.0 && da_I > 10.0) {
        return TRANSPORT_REACTION_DOMINATED;
    }

    /* Mass-transfer limited: Pe moderate, Da_II > 1 (external limitation) */
    if (da_I > 1.0 && da_II < 1.0 && fabs(pe) > 1.0) {
        return TRANSPORT_MASS_TRANSFER_LIMITED;
    }

    /* Dispersion regime: moderate Pe with significant axial mixing */
    if (fabs(pe) >= 0.1 && fabs(pe) <= 100.0 && da_II > 0.1) {
        return TRANSPORT_DISPERSION_REGIME;
    }

    /* Default: mixed */
    return TRANSPORT_MIXED_REGIME;
}

/* ---------------------------------------------------------------------------
 * L3: Geometry and System Setup
 * ------------------------------------------------------------------------- */

double cdr_characteristic_length(CatalystGeometry geometry, double radius)
{
    /* L_c = V_p / S_p (volume-to-external-surface-area ratio) */
    if (radius <= 0.0) {
        return 0.0;
    }

    switch (geometry) {
        case GEOM_SLAB:
            /* For a slab of half-thickness L: V_p/S_p = L */
            return radius;  /* radius stores half-thickness L */

        case GEOM_CYLINDER:
            /* For infinite cylinder radius R: V_p/S_p = pi*R^2*L/(2*pi*R*L) = R/2 */
            return radius / 2.0;

        case GEOM_SPHERE:
            /* For sphere radius R: V_p/S_p = (4/3*pi*R^3)/(4*pi*R^2) = R/3 */
            return radius / 3.0;

        case GEOM_ARBITRARY_3D:
            /* Default to spherical approximation if no shape info */
            return radius / 3.0;

        default:
            return 0.0;
    }
}

void cdr_system_init(CDRSystem *sys)
{
    if (sys == NULL) {
        return;
    }

    memset(sys, 0, sizeof(CDRSystem));

    /* Mixture defaults (air-like at STP) */
    sys->mixture.temperature = 298.15;
    sys->mixture.pressure    = 101325.0;
    sys->mixture.density     = 1.184;
    sys->mixture.viscosity   = 1.85e-5;
    sys->mixture.heat_capacity = 1005.0;
    sys->mixture.thermal_conductivity = 0.026;
    sys->mixture.n_species   = 1;

    /* Geometry defaults */
    sys->characteristic_L = 0.01;
    sys->characteristic_u = 0.1;
    sys->residence_time   = 1.0;
    sys->porosity         = 0.4;
    sys->tortuosity       = 3.0;

    /* Regime defaults */
    sys->regime = TRANSPORT_MIXED_REGIME;
    sys->phase  = PHASE_HOMOGENEOUS;
}

void cdr_compute_dimensionless(CDRSystem *sys)
{
    if (sys == NULL) {
        return;
    }

    double rho   = sys->mixture.density;
    double mu    = sys->mixture.viscosity;
    double u     = sys->characteristic_u;
    double L     = sys->characteristic_L;
    double tau   = sys->residence_time;

    /* Compute Reynolds number */
    sys->reynolds.Re = cdr_reynolds(rho, u, L, mu);

    /* If we had D_AB, compute Schmidt and Peclet */
    /* These would require species-specific diffusivity, initialized to 0 */
    sys->schmidt.Sc          = 0.0;
    sys->peclet.mass_peclet  = 0.0;
    sys->peclet.heat_peclet   = 0.0;
    sys->peclet.dispersion_pe = 0.0;
    sys->thiele.phi          = 0.0;
    sys->damkohler.Da_I      = 0.0;
    sys->damkohler.Da_II     = 0.0;

    /* Compute heat Peclet if thermal conductivity available */
    if (sys->mixture.thermal_conductivity > 0.0 && rho > 0.0 &&
        sys->mixture.heat_capacity > 0.0) {
        double alpha = sys->mixture.thermal_conductivity /
                       (rho * sys->mixture.heat_capacity);
        sys->peclet.heat_peclet = cdr_peclet_heat(u, L, alpha);
    }

    /* Regime classification (needs D_AB, so deferred) */
    (void)tau;
}


#include "micro_nano_transport.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

/*
 * knudsen_regime.c - Knudsen Number and Flow Regime Classification
 *
 * The Knudsen number Kn = lambda / L is the fundamental parameter
 * determining whether continuum transport theory applies at small scales.
 * Knowledge: L1 (definitions), L2 (regime concepts), L3 (scales)
 *
 * Reference: Karniadakis, Beskok, Aluru (2005) Microflows and Nanoflows
 */

/* Mean free path from hard-sphere kinetic theory.
 * lambda = 1 / (sqrt(2) * pi * d^2 * n), n = p/(kB*T) */
void compute_mean_free_path(FluidProperties *fluid,
                            double pressure,
                            double *mean_free_path)
{
    if (fluid == NULL || mean_free_path == NULL || pressure <= 0.0 ||
        fluid->temperature <= 0.0 || fluid->collision_diameter <= 0.0) {
        if (mean_free_path) *mean_free_path = -1.0;
        return;
    }
    double n = pressure / (KB * fluid->temperature);
    double d2 = fluid->collision_diameter * fluid->collision_diameter;
    *mean_free_path = 1.0 / (M_SQRT2 * M_PI * d2 * n);
}

/* Kn = lambda / L_char, with regime classification */
void compute_knudsen_number(double mean_free_path,
                            double characteristic_length,
                            KnudsenState *state)
{
    if (state == NULL || characteristic_length <= 0.0) return;
    state->mean_free_path_m = mean_free_path;
    state->characteristic_len_m = characteristic_length;
    state->knudsen_number = mean_free_path / characteristic_length;
    double Kn = state->knudsen_number;
    if (Kn < 1.0e-3)      state->regime = KNUDSEN_CONTINUUM;
    else if (Kn < 0.1)    state->regime = KNUDSEN_SLIP_FLOW;
    else if (Kn < 10.0)   state->regime = KNUDSEN_TRANSITION;
    else                  state->regime = KNUDSEN_FREE_MOLECULAR;
}

/* Validate continuum hypothesis using 4 criteria:
 * Kn < 0.001, N_mol >> 1, tau_col << tau_flow, lambda << L */
int validate_continuum_hypothesis(const KnudsenState *kn,
                                   const FluidProperties *fluid,
                                   double pressure)
{
    if (kn == NULL || fluid == NULL || pressure <= 0.0) return -1;
    int criteria_met = 0;
    if (kn->knudsen_number < 1.0e-3) criteria_met++;
    double n = pressure / (KB * fluid->temperature);
    double V = kn->characteristic_len_m * kn->characteristic_len_m * kn->characteristic_len_m;
    if (n * V > 1.0e6) criteria_met++;
    if (kn->knudsen_number < 0.01) criteria_met++;
    double c_bar = sqrt(8.0 * KB * fluid->temperature / (M_PI * fluid->molecular_mass));
    double tau_col = kn->mean_free_path_m / c_bar;
    double tau_flow = kn->characteristic_len_m / c_bar;
    if (tau_col < 0.01 * tau_flow) criteria_met++;
    return criteria_met;
}

/* VHS model mean free path with temperature-dependent cross-section */
double compute_mean_free_path_vhs(FluidProperties *fluid,
                                   double pressure,
                                   double omega,
                                   double ref_T)
{
    if (fluid == NULL || pressure <= 0.0 || omega < 0.3 || omega > 1.2)
        return -1.0;
    double n = pressure / (KB * fluid->temperature);
    double d2 = fluid->collision_diameter * fluid->collision_diameter;
    double lambda_hs = 1.0 / (M_SQRT2 * M_PI * d2 * n);
    double T_ratio = fluid->temperature / ref_T;
    return lambda_hs * pow(T_ratio, omega - 0.5);
}

/* Effective transport length incorporating rarefaction and roughness */
double compute_effective_transport_length(const ChannelGeometry *geom,
                                           const KnudsenState *kn,
                                           double alpha, double beta)
{
    if (geom == NULL || kn == NULL) return -1.0;
    double L = geom->hydraulic_diameter_m;
    if (L <= 0.0) return -1.0;
    double rough = geom->surface_roughness_nm * 1.0e-9 / L;
    return L * (1.0 + alpha * kn->knudsen_number + beta * rough);
}

/* Detailed 6-subregime classification with textual description */
int classify_flow_regime_detailed(const KnudsenState *kn,
                                   char *name, size_t nlen,
                                   double *inv_kn)
{
    if (kn == NULL || name == NULL || inv_kn == NULL) return -1;
    *inv_kn = 1.0 / kn->knudsen_number;
    double Kn = kn->knudsen_number;
    if (Kn < 1.0e-3) {
        snprintf(name, nlen, "Continuum (no-slip)");
        return 0;
    } else if (Kn < 1.0e-2) {
        snprintf(name, nlen, "Continuum with slip (weak rarefaction)");
        return 1;
    } else if (Kn < 0.1) {
        snprintf(name, nlen, "Slip flow (moderate rarefaction)");
        return 2;
    } else if (Kn < 1.0) {
        snprintf(name, nlen, "Transition I (Burnett/DSMC required)");
        return 3;
    } else if (Kn < 10.0) {
        snprintf(name, nlen, "Transition II (near-free-molecular)");
        return 4;
    } else {
        snprintf(name, nlen, "Free molecular flow (Knudsen diffusion)");
        return 5;
    }
}

/* Rarefaction parameter: delta = 1/Kn, hydrodynamic limit -> inf */
double compute_rarefaction_parameter(const KnudsenState *kn)
{
    if (kn == NULL || kn->knudsen_number <= 0.0) return -1.0;
    return 1.0 / kn->knudsen_number;
}

/* von Karman relation: Kn = sqrt(pi*gamma/2) * (Ma/Re) */
void compare_reynolds_knudsen(double Re, double Ma, double gamma,
                               double *kn, double *delta)
{
    if (kn == NULL || delta == NULL) return;
    *kn = sqrt(M_PI * gamma / 2.0) * (Ma / Re);
    *delta = (Re > 0.0 && Ma > 0.0) ? (1.0 / *kn) : -1.0;
}

/* Estimate flow regime from accessible measurements Re, Ma */
KnudsenRegime estimate_flow_regime_from_re_ma(double Re, double Ma, double gamma)
{
    if (Re <= 0.0 || Ma <= 0.0) return KNUDSEN_FREE_MOLECULAR;
    double Kn = sqrt(M_PI * gamma / 2.0) * (Ma / Re);
    if (Kn < 1.0e-3) return KNUDSEN_CONTINUUM;
    if (Kn < 0.1)    return KNUDSEN_SLIP_FLOW;
    if (Kn < 10.0)   return KNUDSEN_TRANSITION;
    return KNUDSEN_FREE_MOLECULAR;
}

/* Collision frequency: Z = n * sigma * c_bar */
double compute_collision_frequency(double T, double p,
                                    double m_mol, double d_col)
{
    if (T <= 0.0 || p <= 0.0 || m_mol <= 0.0 || d_col <= 0.0)
        return -1.0;
    double n = p / (KB * T);
    double sigma = M_PI * d_col * d_col;
    double c_bar = sqrt(8.0 * KB * T / (M_PI * m_mol));
    return n * sigma * c_bar;
}

/* Auto-select characteristic length from geometry */
double compute_characteristic_length_automatic(const ChannelGeometry *geom)
{
    if (geom == NULL) return -1.0;
    switch (geom->type) {
    case GEOM_PARALLEL_PLATES:    return geom->height_m;
    case GEOM_CYLINDRICAL_CAPILLARY: return 2.0 * geom->radius_m;
    case GEOM_RECTANGULAR_CHANNEL: return geom->height_m;
    case GEOM_NANOPORE:           return geom->height_m;
    default: return geom->hydraulic_diameter_m;
    }
}

/**
 * @file cdr_mass_transfer.c
 * @brief Implementation of interphase mass transfer and correlations.
 *
 * Sherwood number correlations for various geometries,
 * two-film theory, enhancement by chemical reaction,
 * and HTU/NTU methods.
 *
 * Reference: Bird, Stewart, Lightfoot "Transport Phenomena" (2007), Ch. 21-22
 * Reference: Treybal "Mass Transfer Operations" (1980)
 */

#include "cdr_mass_transfer.h"
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---------------------------------------------------------------------------
 * L3: Sherwood Number Correlations
 * ------------------------------------------------------------------------- */

double cdr_sherwood_flat_plate_laminar(double Re, double Sc)
{
    /*
     * Sh_L = 0.664 * Re^(1/2) * Sc^(1/3)   [Re < 5e5]
     *
     * From boundary layer theory, analogy to heat transfer.
     * Valid: 0.6 < Sc < 50, Re < 5e5
     */
    if (Re <= 0.0 || Sc <= 0.0) {
        return 0.0;
    }
    return 0.664 * sqrt(Re) * pow(Sc, 1.0/3.0);
}

double cdr_sherwood_flat_plate_turbulent(double Re, double Sc)
{
    /*
     * Sh_L = 0.037 * Re^(4/5) * Sc^(1/3)   [Re > 5e5]
     */
    if (Re <= 0.0 || Sc <= 0.0) {
        return 0.0;
    }
    return 0.037 * pow(Re, 0.8) * pow(Sc, 1.0/3.0);
}

double cdr_sherwood_sphere(double Re, double Sc)
{
    /*
     * Sh = 2.0 + 0.552 * Re^(1/2) * Sc^(1/3)   (Froessling, 1938)
     *
     * The 2.0 comes from pure molecular diffusion around a sphere
     * in a stagnant fluid (Re=0 limit).
     * Valid: 2 < Re < 800, 0.6 < Sc < 2.7 (gas), extended in practice.
     */
    if (Re < 0.0 || Sc <= 0.0) {
        return 2.0;  /* Diffusion-only limit */
    }
    return 2.0 + 0.552 * sqrt(Re) * pow(Sc, 1.0/3.0);
}

double cdr_sherwood_packed_bed(double Re, double Sc, double epsilon)
{
    /*
     * Thoenes-Kramers correlation for packed beds:
     *
     * Sh = 1.0 * (Re/(1-epsilon))^(1/2) * Sc^(1/3) * epsilon
     *
     * Valid: 0.25 < epsilon < 0.5, 40 < Re/(1-epsilon) < 4000, 0.6 < Sc < 1000
     */
    if (Re <= 0.0 || Sc <= 0.0 || epsilon <= 0.0 || epsilon >= 1.0) {
        return 0.0;
    }

    double Re_modified = Re / (1.0 - epsilon);
    return 1.0 * sqrt(Re_modified) * pow(Sc, 1.0/3.0) * epsilon;
}

double cdr_mass_transfer_coeff_from_sh(double Sh, double D, double L)
{
    /*
     * k_c = Sh * D / L
     */
    if (L <= 0.0 || D < 0.0 || Sh < 0.0) {
        return 0.0;
    }
    return Sh * D / L;
}

double cdr_mass_transfer_penetration(double D, double t_e)
{
    /*
     * k_L = 2 * sqrt(D / (pi * t_e))
     *
     * Higbie's penetration theory for gas-liquid mass transfer.
     * Assumes unsteady-state diffusion into a stagnant liquid element
     * with exposure time t_e.
     */
    if (D <= 0.0 || t_e <= 0.0) {
        return 0.0;
    }
    return 2.0 * sqrt(D / (M_PI * t_e));
}

double cdr_mass_transfer_falling_film(double D, double u_s, double L)
{
    /*
     * Laminar falling film, short contact time (penetration theory):
     * k_L_avg = 2 * sqrt(D * u_s / (pi * L))
     *
     * Valid when penetration depth << film thickness.
     *
     * For long contact (fully developed):
     * k_L_avg = 3.41 * D / delta  (where delta is film thickness)
     *
     * Here we use the penetration form (short contact limit).
     */
    if (D <= 0.0 || u_s <= 0.0 || L <= 0.0) {
        return 0.0;
    }
    return 2.0 * sqrt(D * u_s / (M_PI * L));
}

/* ---------------------------------------------------------------------------
 * L5: Two-Film Theory
 * ------------------------------------------------------------------------- */

void cdr_two_film_overall_gas(double k_G, double k_L, double H, double P,
                               double *K_overall)
{
    /*
     * Gas-side overall mass transfer coefficient:
     *
     * 1/K_OG = 1/k_G + H/(k_L * P)   [if k_G in mol/(m^2.s.Pa), k_L in m/s]
     *
     * Explanation:
     *   Resistance in gas film: 1/k_G
     *   Resistance in liquid film: H/(k_L * P)   [converting to gas-side units]
     *
     * H = Henry's constant [Pa.m^3/mol]
     * P = total pressure [Pa]
     */
    if (K_overall == NULL) {
        return;
    }

    if (k_G <= 0.0 || k_L <= 0.0 || H <= 0.0 || P <= 0.0) {
        *K_overall = 0.0;
        return;
    }

    *K_overall = 1.0 / (1.0 / k_G + H / (k_L * P));
}

double cdr_two_film_flux(double K_overall, double p_A_bulk,
                          double C_A_bulk, double H)
{
    /*
     * N_A = K_OG * (p_A_bulk - H * C_A_bulk)
     *
     * H*C_A_bulk is the partial pressure of A in equilibrium with the
     * bulk liquid concentration (driving force on gas-side basis).
     */
    if (p_A_bulk < 0.0 || C_A_bulk < 0.0 || H <= 0.0) {
        return 0.0;
    }

    double p_A_eq = H * C_A_bulk;  /* Equilibrium partial pressure */
    return K_overall * (p_A_bulk - p_A_eq);
}

/* ---------------------------------------------------------------------------
 * L5: Mass Transfer with Chemical Reaction
 * ------------------------------------------------------------------------- */

double cdr_enhancement_factor_homogeneous(double k, double D, double k_L)
{
    /*
     * Enhancement factor for irreversible 1st-order reaction in liquid film:
     *
     * E = Ha / tanh(Ha)
     *
     * where Ha = sqrt(k*D) / k_L is the Hatta number.
     *
     * E ~ Ha (Ha > 3): Fast reaction, reaction mainly in film
     * E ~ 1 (Ha < 0.3): Slow reaction, reaction mainly in bulk
     */
    if (k_L <= 0.0) {
        return 1.0;
    }
    if (k <= 0.0 || D <= 0.0) {
        return 1.0;  /* No reaction enhancement */
    }

    double Ha = sqrt(k * D) / k_L;

    /* For small Ha, tanh(Ha) ? Ha ? E ? 1 */
    if (Ha < 0.001) {
        return 1.0;
    }

    /* For large Ha, tanh(Ha) ? 1 ? E ? Ha */
    if (Ha > 50.0) {
        return Ha;
    }

    double tanh_Ha = tanh(Ha);
    if (fabs(tanh_Ha) < 1e-15) {
        return Ha;
    }

    return Ha / tanh_Ha;
}

double cdr_hatta_number(double k1, double D_A, double k_L)
{
    /*
     * Ha = sqrt(k1 * D_A) / k_L
     *
     * For pseudo-first-order reaction with excess B.
     *
     * Ha < 0.3: Slow reaction regime
     * 0.3 < Ha < 3: Intermediate regime
     * Ha > 3: Fast reaction regime (reaction in film)
     * Ha > 10*E_i: Instantaneous regime (reaction plane at interface)
     */
    if (k_L <= 0.0) {
        return INFINITY;
    }
    if (k1 <= 0.0 || D_A <= 0.0) {
        return 0.0;
    }

    return sqrt(k1 * D_A) / k_L;
}

const char* cdr_classify_hatta_regime(double Ha)
{
    /*
     * Classify gas-liquid reaction regime based on Hatta number.
     */
    if (Ha < 0.0) {
        return "Invalid (Ha < 0)";
    }
    if (Ha < 0.3) {
        return "Slow reaction (kinetics-controlled, bulk liquid)";
    }
    if (Ha < 3.0) {
        return "Intermediate regime (reaction in film and bulk)";
    }
    if (Ha < 50.0) {
        return "Fast reaction (reaction predominantly in liquid film)";
    }
    return "Instantaneous reaction (reaction plane at interface)";
}

/* ---------------------------------------------------------------------------
 * L5: HTU/NTU Method
 * ------------------------------------------------------------------------- */

double cdr_ntu_gas_phase(double y_in, double y_out, double y_eq)
{
    /*
     * NTU_G = (y_in - y_out) / Delta_y_lm
     *
     * Log-mean driving force:
     * Delta_y_lm = (Delta_y_in - Delta_y_out) / ln(Delta_y_in/Delta_y_out)
     *
     * where Delta_y_in = y_in - y_eq_in, Delta_y_out = y_out - y_eq_out
     *
     * For linear equilibrium: y_eq = m * x
     */
    if (y_in <= y_out) {
        return 0.0;
    }

    /* Simplified: using y_eq as equilibrium outlet (for linear case) */
    double Delta_y_in  = y_in - y_eq;
    double Delta_y_out = y_out - y_eq;

    if (Delta_y_in <= 0.0 || Delta_y_out <= 0.0) {
        return 0.0;
    }
    if (fabs(Delta_y_in - Delta_y_out) < 1e-15) {
        return (y_in - y_out) / Delta_y_in;
    }

    double Delta_y_lm = (Delta_y_in - Delta_y_out) / log(Delta_y_in / Delta_y_out);
    return (y_in - y_out) / Delta_y_lm;
}

double cdr_htu_gas_phase(double G, double k_y, double a)
{
    /*
     * HTU_G = G / (k_y * a)
     *
     * G = molar gas velocity [mol/(m^2.s)]
     * k_y = gas-side mass transfer coefficient [mol/(m^2.s)]
     * a = interfacial area per unit volume [m^2/m^3]
     *
     * Z_packed = HTU_G * NTU_G  (height of packed section)
     */
    if (k_y <= 0.0 || a <= 0.0) {
        return 0.0;
    }
    return G / (k_y * a);
}

/* ---------------------------------------------------------------------------
 * L5: Chilton-Colburn Analogy
 * ------------------------------------------------------------------------- */

/**
 * @brief Compute j_D factor from Sherwood, Reynolds, and Schmidt numbers.
 *
 * j_D = Sh / (Re * Sc^(1/3))
 *
 * @param Sh  Sherwood number [-]
 * @param Re  Reynolds number [-]
 * @param Sc  Schmidt number [-]
 * @return j_D factor [-]
 */
double cdr_j_factor_mass(double Sh, double Re, double Sc)
{
    if (Re <= 0.0 || Sc <= 0.0) {
        return 0.0;
    }
    return Sh / (Re * pow(Sc, 1.0/3.0));
}

/**
 * @brief Estimate mass transfer coefficient using Chilton-Colburn analogy.
 *
 * k_c = (f/2) * u / Sc^(2/3)
 *
 * where f is the Fanning friction factor.
 *
 * @param f   Fanning friction factor [-]
 * @param u   Free-stream velocity [m/s]
 * @param Sc  Schmidt number [-]
 * @return Mass transfer coefficient [m/s]
 */
double cdr_mass_transfer_chilton_colburn(double f, double u, double Sc)
{
    if (Sc <= 0.0) {
        return 0.0;
    }
    return (f / 2.0) * u / pow(Sc, 2.0/3.0);
}

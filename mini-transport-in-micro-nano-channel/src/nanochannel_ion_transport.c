#include "micro_nano_transport.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Ionic Transport in Nanochannels */

double compute_ionic_current_nanochannel(int n_species,
    const double *conc, const int *valence, const double *D,
    double V_applied, double L, double A, double T)
{
    if (conc == NULL || valence == NULL || D == NULL
        || n_species <= 0 || L <= 0.0 || A <= 0.0 || T <= 0.0)
        return -1.0;
    double conductivity = 0.0;
    double e2_over_kT = E_CHARGE * E_CHARGE / (KB * T);
    for (int i = 0; i < n_species; i++) {
        double zi = (double)valence[i];
        conductivity += zi * zi * D[i] * conc[i] * NA * e2_over_kT;
    }
    return conductivity * A * V_applied / L;
}
double compute_nanochannel_conductance(double c_bulk, double z,
    double D, double L, double A, double sigma_s, double H, double T)
{
    if (c_bulk < 0.0 || L <= 0.0 || A <= 0.0 || T <= 0.0 || H <= 0.0)
        return -1.0;
    double e2_over_kT = E_CHARGE * E_CHARGE / (KB * T);
    double excess_conc = fabs(sigma_s) / (E_CHARGE * H);
    double effective_conc = c_bulk + excess_conc;
    return z * z * D * effective_conc * NA * e2_over_kT * A / L;
}

double compute_donnan_potential(double sigma_s, double c_bulk,
                                 double z, double H, double T)
{
    if (c_bulk <= 0.0 || T <= 0.0 || H <= 0.0 || z == 0.0) return 0.0;
    double denominator = 2.0 * fabs(z) * E_CHARGE * c_bulk * H * NA;
    if (denominator <= 0.0) return 0.0;
    double arg = sigma_s / denominator;
    if (fabs(arg) > 50.0) arg = 50.0 * ((arg > 0) ? 1.0 : -1.0);
    return (KB * T / (fabs(z) * E_CHARGE))
           * log(fabs(arg) + sqrt(arg * arg + 1.0));
}

double compute_ion_selectivity(int n_species, const double *conc,
                                const int *valence, const double *D)
{
    if (conc == NULL || valence == NULL || D == NULL || n_species < 2)
        return 0.0;
    double t_plus = 0.0, t_minus = 0.0;
    double denominator = 0.0;
    for (int i = 0; i < n_species; i++) {
        double weight = fabs((double)valence[i]) * D[i] * conc[i];
        denominator += weight;
        if (valence[i] > 0) t_plus += weight;
        else if (valence[i] < 0) t_minus += weight;
    }
    if (denominator <= 0.0) return 0.0;
    t_plus /= denominator;
    t_minus /= denominator;
    return (t_plus - t_minus);
}

double compute_ion_rectification_ratio(double sigma_in, double sigma_out,
                                        double z, double T)
{
    if (T <= 0.0 || z == 0.0) return 1.0;
    double sigma_diff = sigma_in - sigma_out;
    double delta_V = fabs(sigma_diff) * 1.0e-9;
    double exponent = fabs(z) * E_CHARGE * delta_V / (KB * T);
    if (exponent > 10.0) exponent = 10.0;
    return exp(exponent);
}

double compute_nanopore_resistance(double L, double r, double sigma_bulk)
{
    if (L <= 0.0 || r <= 0.0 || sigma_bulk <= 0.0) return -1.0;
    double R_channel = L / (sigma_bulk * M_PI * r * r);
    double R_access = 1.0 / (2.0 * sigma_bulk * r);
    return R_channel + 2.0 * R_access;
}

double compute_nanopore_blockade_current(double I_open, double r_pore,
                                          double d_molecule)
{
    if (I_open <= 0.0 || r_pore <= 0.0 || d_molecule < 0.0) return -1.0;
    if (d_molecule >= 2.0 * r_pore) return I_open;
    double area_ratio = (d_molecule * d_molecule) / (r_pore * r_pore);
    if (area_ratio > 1.0) area_ratio = 1.0;
    return I_open * area_ratio;
}

double compute_overlap_parameter(double kappa, double H)
{
    if (kappa <= 0.0 || H <= 0.0) return -1.0;
    return kappa * H / 2.0;
}

double compute_coion_exclusion_coefficient(double psi_center,
                                            int z_co, double T)
{
    if (T <= 0.0) return 1.0;
    double exponent = -fabs((double)z_co) * E_CHARGE * psi_center / (KB * T);
    if (exponent < -50.0) return 0.0;
    return exp(exponent);
}

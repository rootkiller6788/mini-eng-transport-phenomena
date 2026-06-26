#include "micro_nano_transport.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * poisson_boltzmann.c - Poisson-Boltzmann Equation Solvers
 * Knowledge: L1 (Debye length), L4 (PB eq), L5 (solvers)
 * Reference: Hunter (1981), Russel, Saville, Schowalter (1989)
 */

void compute_debye_length(DebyeState *state, double *lambda_d)
{
    if (state == NULL || lambda_d == NULL
        || state->num_ion_species <= 0
        || state->temperature <= 0.0
        || state->solvent_permittivity <= 0.0) {
        if (lambda_d) *lambda_d = -1.0;
        return;
    }

    double kappa_sq = 0.0;
    for (int i = 0; i < state->num_ion_species && i < MAX_ION_SPECIES; i++) {
        double zi = (double)state->valences[i];
        double ci = state->concentrations[i];
        kappa_sq += zi * zi * ci * NA;
    }
    kappa_sq *= (E_CHARGE * E_CHARGE)
                / (state->solvent_permittivity * KB * state->temperature);

    double kappa = sqrt(kappa_sq);
    state->debye_parameter = kappa;
    state->debye_length_m = 1.0 / kappa;

    double I = 0.0;
    for (int i = 0; i < state->num_ion_species && i < MAX_ION_SPECIES; i++) {
        double zi = (double)state->valences[i];
        I += state->concentrations[i] * zi * zi;
    }
    state->ionic_strength = 0.5 * I;

    *lambda_d = state->debye_length_m;
}

void compute_poisson_boltzmann_1d(const double *y, int n,
                                   double psi0, double kappa,
                                   double *psi)
{
    if (y == NULL || psi == NULL || n <= 0 || kappa <= 0.0) return;

    double scaled_potential = fabs(E_CHARGE * psi0 / (KB * 293.15));
    int use_nonlinear = (scaled_potential > 0.5);

    for (int i = 0; i < n; i++) {
        if (use_nonlinear) {
            double z = 1.0;
            double arg = z * E_CHARGE * psi0 / (4.0 * KB * 293.15);
            if (fabs(arg) > 10.0) {
                psi[i] = psi0 * exp(-kappa * fabs(y[i]));
            } else {
                double gamma = tanh(arg);
                double exp_term = exp(-kappa * fabs(y[i]));
                double numerator = 1.0 + gamma * exp_term;
                double denominator = 1.0 - gamma * exp_term;
                if (denominator <= 0.0) {
                    psi[i] = psi0 * exp(-kappa * fabs(y[i]));
                } else {
                    psi[i] = (2.0 * KB * 293.15 / (z * E_CHARGE))
                             * log(numerator / denominator);
                }
            }
        } else {
            psi[i] = psi0 * exp(-kappa * fabs(y[i]));
        }
    }
}

void compute_ion_concentration_profile(const double *psi, int n,
                                        double c_bulk, int z, double T,
                                        double *c_profile)
{
    if (psi == NULL || c_profile == NULL || n <= 0
        || T <= 0.0 || c_bulk < 0.0) return;

    double beta = 1.0 / (KB * T);
    for (int i = 0; i < n; i++) {
        double exponent = -(double)z * E_CHARGE * psi[i] * beta;
        if (exponent > 50.0) exponent = 50.0;
        if (exponent < -50.0) exponent = -50.0;
        c_profile[i] = c_bulk * exp(exponent);
    }
}

double compute_surface_charge_potential(double sigma_s, double c_bulk,
                                         double z, double eps, double T)
{
    if (c_bulk <= 0.0 || eps <= 0.0 || T <= 0.0) return 0.0;

    double prefactor = sqrt(8.0 * eps * KB * T * c_bulk * NA);
    if (prefactor <= 0.0) return 0.0;

    double scaled = sigma_s / prefactor;
    double arg = scaled;
    double psi0 = (2.0 * KB * T / (fabs(z) * E_CHARGE))
                  * log(fabs(arg) + sqrt(arg * arg + 1.0));

    return (sigma_s > 0.0) ? psi0 : -psi0;
}

double compute_surface_charge_from_potential(double psi0, double c_bulk,
                                              double z, double eps, double T)
{
    if (c_bulk <= 0.0 || eps <= 0.0 || T <= 0.0) return 0.0;

    double prefactor = sqrt(8.0 * eps * KB * T * c_bulk * NA);
    double arg = fabs(z) * E_CHARGE * psi0 / (2.0 * KB * T);
    if (arg > 20.0) arg = 20.0;

    return prefactor * sinh(arg) * ((psi0 > 0) ? 1.0 : -1.0);
}

double compute_differential_capacitance(double psi0, double kappa,
                                         double eps, double z, double T)
{
    if (kappa <= 0.0 || eps <= 0.0 || T <= 0.0) return -1.0;
    double arg = fabs(z) * E_CHARGE * psi0 / (2.0 * KB * T);
    if (arg > 20.0) arg = 20.0;
    return eps * kappa * cosh(arg);
}

void compute_stern_layer_correction(double sigma_s, double d_stern,
                                     double eps_stern, double zeta_in,
                                     double *psi_0, double *psi_stern)
{
    if (psi_0 == NULL || psi_stern == NULL) return;
    if (d_stern <= 0.0 || eps_stern <= 0.0) {
        *psi_stern = 0.0;
        *psi_0 = zeta_in;
        return;
    }
    *psi_stern = sigma_s * d_stern / eps_stern;
    *psi_0 = zeta_in + *psi_stern;
}

double compute_pb_free_energy(double sigma_s, double psi0)
{
    return 0.5 * sigma_s * psi0;
}

double compute_effective_screening_length(double kappa, double sigma_s,
                                           double eps, double T)
{
    if (kappa <= 0.0 || eps <= 0.0 || T <= 0.0) return -1.0;
    double l_B = E_CHARGE * E_CHARGE / (4.0 * M_PI * eps * KB * T);
    double sigma_c = E_CHARGE / (2.0 * M_PI * l_B * l_B);
    if (sigma_c <= 0.0) return 1.0 / kappa;
    double correction = sqrt(1.0 + (sigma_s / sigma_c) * (sigma_s / sigma_c));
    return correction / kappa;
}

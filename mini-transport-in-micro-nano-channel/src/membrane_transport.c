#include "micro_nano_transport.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * membrane_transport.c - Membrane Transport Processes
 *
 * Membrane-based separation processes including reverse osmosis (RO),
 * forward osmosis (FO), and pressure-retarded osmosis (PRO) represent
 * key applications of nanofluidic transport principles.
 *
 * Knowledge: L7 (membrane applications), L5 (solution-diffusion model)
 * Reference: Baker (2012) Membrane Technology and Applications
 *            Wijmans & Baker (1995) J. Membr. Sci. 107, 1
 */

/*
 * Solution-diffusion model for RO membrane water flux:
 *   J_w = A * (delta_p - delta_pi)
 * where A is the water permeability coefficient,
 * delta_p is the hydraulic pressure difference,
 * delta_pi is the osmotic pressure difference.
 */
double compute_ro_water_flux(double A, double delta_p, double delta_pi)
{
    if (A <= 0.0) return -1.0;
    double driving_force = delta_p - delta_pi;
    if (driving_force < 0.0) return 0.0;
    return A * driving_force;
}

double compute_osmotic_pressure(double c, double T, double phi, int n_ions)
{
    if (c <= 0.0 || T <= 0.0 || n_ions <= 0) return -1.0;
    return (double)n_ions * phi * c * R_GAS * T;
}

double compute_salt_rejection(double c_permeate, double c_feed)
{
    if (c_feed <= 0.0) return -1.0;
    if (c_permeate < 0.0) c_permeate = 0.0;
    return 1.0 - c_permeate / c_feed;
}

double compute_water_recovery(double Q_permeate, double Q_feed)
{
    if (Q_feed <= 0.0) return -1.0;
    return Q_permeate / Q_feed;
}

double compute_specific_energy_consumption(double P_feed, double recovery,
    double pump_efficiency, double ERD_efficiency)
{
    if (recovery <= 0.0 || recovery > 1.0 || pump_efficiency <= 0.0)
        return -1.0;
    return P_feed / (recovery * pump_efficiency) * (1.0 - ERD_efficiency);
}

double compute_concentration_factor(double recovery, double rejection)
{
    if (recovery >= 1.0 || recovery < 0.0) return -1.0;
    double numerator = 1.0 - recovery * (1.0 - (1.0 - rejection));
    double denominator = 1.0 - recovery;
    if (denominator <= 0.0) return -1.0;
    return numerator / denominator;
}

double compute_pro_power_density(double A, double delta_pi,
    double delta_p_hyd, double turbine_efficiency)
{
    if (A <= 0.0 || turbine_efficiency <= 0.0) return -1.0;
    double J_w = A * (delta_pi - delta_p_hyd);
    if (J_w <= 0.0) return 0.0;
    return J_w * delta_p_hyd * turbine_efficiency;
}

static double Jw_guess(double A, double B, double pi_d, double pi_f,
    double K, double D);

double compute_fo_water_flux(double A, double B, double pi_draw,
    double pi_feed, double K, double D)
{
    if (A <= 0.0 || K <= 0.0) return -1.0;
    double pi_draw_eff = pi_draw * exp(-Jw_guess(A, B, pi_draw, pi_feed, K, D) * K);
    (void)pi_draw_eff;
    return A * (pi_draw - pi_feed);
}

static double Jw_guess(double A, double B, double pi_d, double pi_f,
    double K, double D)
{
    (void)B; (void)D;
    return A * (pi_d - pi_f) / (1.0 + A * K * pi_d);
}

double compute_internal_concentration_polarization(double J_w,
    double K_solute, double D_solute)
{
    if (K_solute <= 0.0 || D_solute <= 0.0) return -1.0;
    double Pe = J_w * K_solute / D_solute;
    if (Pe > 50.0) Pe = 50.0;
    return exp(Pe);
}

double compute_membrane_structural_parameter(double S, double D)
{
    if (D <= 0.0) return -1.0;
    return S / D;
}

double compute_solute_permeability(double B, double c_feed,
    double c_permeate)
{
    double delta_c = c_feed - c_permeate;
    if (delta_c <= 0.0) return -1.0;
    return B * delta_c;
}

double compute_membrane_selectivity(double P_A, double P_B,
    double l_membrane)
{
    if (P_B <= 0.0 || l_membrane <= 0.0) return -1.0;
    return P_A / P_B;
}

double compute_hagen_poiseuille_membrane(double r_pore, double porosity,
    double tortuosity, double L, double eta, double dp)
{
    if (r_pore <= 0.0 || porosity <= 0.0 || tortuosity <= 0.0
        || L <= 0.0 || eta <= 0.0) return -1.0;
    double J = (porosity * r_pore * r_pore / (8.0 * eta * tortuosity * L)) * dp;
    return J;
}

double compute_kozeny_carman_permeability(double d_particle, double porosity)
{
    if (d_particle <= 0.0 || porosity <= 0.0 || porosity >= 1.0) return -1.0;
    double phi3 = porosity * porosity * porosity;
    double one_minus_phi_sq = (1.0 - porosity) * (1.0 - porosity);
    double K = d_particle * d_particle * phi3 / (180.0 * one_minus_phi_sq);
    return K;
}

double compute_membrane_fouling_index(double J_t, double J_0, double t)
{
    if (J_0 <= 0.0 || t <= 0.0) return -1.0;
    if (J_t <= 0.0) return 1.0e10;
    return (1.0 / J_t - 1.0 / J_0) / t;
}

double compute_mass_transfer_coefficient_laminar(double D, double D_h,
    double Re, double Sc, double L)
{
    if (D <= 0.0 || D_h <= 0.0 || L <= 0.0) return -1.0;
    double Gz = D_h / L * Re * Sc;
    double Sh;
    if (Gz > 100.0) {
        Sh = 1.62 * pow(Gz, 1.0/3.0);
    } else {
        Sh = 3.66 + 0.0668 * Gz / (1.0 + 0.04 * pow(Gz, 2.0/3.0));
    }
    return Sh * D / D_h;
}

double compute_recovery_ratio_n_stages(double recovery_per_stage, int n_stages)
{
    if (recovery_per_stage <= 0.0 || recovery_per_stage >= 1.0 || n_stages <= 0)
        return -1.0;
    double overall = 1.0 - pow(1.0 - recovery_per_stage, (double)n_stages);
    return overall;
}

double compute_energy_efficiency_ro(double TDS, double T, double recovery)
{
    if (TDS <= 0.0 || T <= 0.0 || recovery <= 0.0 || recovery >= 1.0)
        return -1.0;
    double pi_feed = TDS * R_GAS * T * 2.0 / 1000.0;
    double pi_avg = pi_feed * log(1.0 / (1.0 - recovery)) / recovery;
    double E_min = pi_avg / recovery;
    return E_min;
}

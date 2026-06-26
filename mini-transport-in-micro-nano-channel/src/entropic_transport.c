#include "micro_nano_transport.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * entropic_transport.c - Entropic Barriers and Configurational Transport
 *
 * At the nanoscale, transport is often dominated by entropic effects.
 * When molecules or particles pass through constrictions comparable to
 * their size, the reduction in configurational entropy creates an energy
 * barrier that controls the transport rate.
 *
 * Knowledge: L2 (entropic traps), L8 (configurational entropy)
 * Reference: Muthukumar (2011) Polymer Translocation
 *            Han & Craighead (2000) Science 288, 1026
 */

/* Entropic barrier for a spherical particle entering a cylindrical pore.
 * delta_S = kB * ln(A_pore / A_bulk)
 * where A_pore is the accessible cross-sectional area in the pore
 * and A_bulk is the total cross-sectional area in the reservoir.
 *
 * For a particle of radius a and pore of radius R:
 *   A_accessible = pi * (R - a)^2   (if R > a, else 0)
 */
double compute_entropic_barrier(double R_pore, double a_particle,
    double T)
{
    if (R_pore <= 0.0 || a_particle < 0.0 || T <= 0.0) return -1.0;
    if (a_particle >= R_pore) return 1.0e10;
    double A_pore = M_PI * (R_pore - a_particle) * (R_pore - a_particle);
    double A_bulk = M_PI * R_pore * R_pore;
    double delta_S = KB * log(A_pore / A_bulk);
    return -T * delta_S;
}

double compute_entropic_spring_constant(double L_polymer,
    double L_persistence, double T)
{
    if (L_polymer <= 0.0 || L_persistence <= 0.0 || T <= 0.0) return -1.0;
    return 3.0 * KB * T / (2.0 * L_persistence * L_polymer);
}

double compute_entropic_trapping_time(double barrier_energy,
    double attempt_frequency, double T)
{
    if (attempt_frequency <= 0.0 || T <= 0.0) return -1.0;
    if (barrier_energy < 0.0) barrier_energy = 0.0;
    double tau = exp(barrier_energy / (KB * T)) / attempt_frequency;
    if (tau > 1.0e10) tau = 1.0e10;
    return tau;
}

double compute_ogston_sieving_coefficient(double R_pore,
    double R_hydrodynamic)
{
    if (R_pore <= 0.0 || R_hydrodynamic <= 0.0) return -1.0;
    if (R_hydrodynamic >= R_pore) return 0.0;
    double ratio = R_hydrodynamic / R_pore;
    return (1.0 - ratio) * (1.0 - ratio);
}

double compute_free_energy_translocation(double N_monomers,
    double epsilon, double T)
{
    if (N_monomers <= 0.0 || T <= 0.0) return -1.0;
    double delta_F = N_monomers * epsilon - T * KB * log(N_monomers);
    return delta_F;
}

double compute_partition_coefficient_steric(double R_pore,
    double R_solute)
{
    if (R_pore <= 0.0 || R_solute < 0.0) return -1.0;
    if (R_solute >= R_pore) return 0.0;
    double ratio = R_solute / R_pore;
    return (1.0 - ratio) * (1.0 - ratio);
}

double compute_hindered_diffusion_coefficient(double D_bulk,
    double R_pore, double R_solute)
{
    if (D_bulk <= 0.0 || R_pore <= 0.0 || R_solute < 0.0) return -1.0;
    if (R_solute >= R_pore) return 0.0;
    double lambda = R_solute / R_pore;
    double hindrance = 1.0 - 2.1044 * lambda + 2.089 * lambda * lambda * lambda
                       - 0.948 * lambda * lambda * lambda * lambda * lambda;
    if (hindrance < 0.0) hindrance = 0.0;
    return D_bulk * hindrance;
}

double compute_entropic_voltage(double c_high, double c_low,
    double T)
{
    if (c_high <= 0.0 || c_low <= 0.0 || T <= 0.0) return 0.0;
    return (KB * T / E_CHARGE) * log(c_high / c_low);
}

double compute_debye_overlap_potential(double zeta, double kappa,
    double H)
{
    if (kappa <= 0.0 || H <= 0.0) return -1.0;
    double half_H = 0.5 * H;
    return zeta / cosh(kappa * half_H);
}

double compute_bjerrum_length(double eps, double T)
{
    if (eps <= 0.0 || T <= 0.0) return -1.0;
    return E_CHARGE * E_CHARGE / (4.0 * M_PI * eps * KB * T);
}

double compute_gouy_chapman_length(double sigma_s, double l_B, int z)
{
    if (sigma_s == 0.0 || l_B <= 0.0 || z == 0) return -1.0;
    return E_CHARGE / (2.0 * M_PI * l_B * fabs(sigma_s) * fabs((double)z));
}

double compute_overscreening_parameter(double sigma_s, double c_bulk,
    double z, double eps, double T)
{
    if (c_bulk <= 0.0 || eps <= 0.0 || T <= 0.0) return -1.0;
    double l_B = compute_bjerrum_length(eps, T);
    double lambda_D = 1.0 / sqrt(2.0 * z * z * E_CHARGE * E_CHARGE
                                  * c_bulk * NA / (eps * KB * T));
    return l_B * fabs(sigma_s) * lambda_D / E_CHARGE;
}

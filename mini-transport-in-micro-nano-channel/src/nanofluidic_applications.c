/**
 * @file nanofluidic_applications.c
 * @brief Nanofluidic diode, ion selectivity, DNA translocation, lab-on-chip.
 * Ref: Schoch et al. (2008) Rev Mod Phys, Daiguji (2010) Chem Soc Rev.
 */
#include "micro_nano_transport.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef KB
#define KB 1.380649e-23
#endif
#ifndef EC
#define EC 1.602176634e-19
#endif

/* ---- Nanofluidic Diode ---- */

typedef struct {
    double channel_length, channel_width, channel_height;
    double surface_charge_density_left, surface_charge_density_right;
    double bulk_concentration, zeta_left, zeta_right;
    double applied_voltage;
} NanofluidicDiode;

double nanofluidic_diode_rectification_ratio(const NanofluidicDiode *d)
{
    if (!d || d->applied_voltage == 0.0) return 1.0;
    double zeta_L = d->zeta_left, zeta_R = d->zeta_right;
    double V_app = d->applied_voltage;

    double D = 1.0e-9; double mu = EC * D / (KB * 300.0);
    double A = d->channel_width * d->channel_height;
    double L = d->channel_length;
    double c0 = d->bulk_concentration * 6.022e23;

    double sigma_L = fabs(d->surface_charge_density_left);
    double sigma_R = fabs(d->surface_charge_density_right);

    double I_forward = A / L * (c0 * EC * mu * V_app
        + sigma_L * mu * EC * zeta_L * V_app / (KB * 300.0));
    double I_reverse = A / L * (c0 * EC * mu * (-V_app)
        + sigma_R * mu * EC * zeta_R * (-V_app) / (KB * 300.0));

    if (fabs(I_reverse) < 1e-18) return 100.0;
    return fabs(I_forward / I_reverse);
}

double nanofluidic_diode_forward_current(const NanofluidicDiode *d)
{
    if (!d) return 0.0;
    double D = 1.0e-9; double mu = EC * D / (KB * 300.0);
    double A = d->channel_width * d->channel_height, L = d->channel_length;
    double c0 = d->bulk_concentration * 6.022e23;
    double sigma_avg = 0.5 * (d->surface_charge_density_left + d->surface_charge_density_right);
    return A / L * (c0 * EC * mu * d->applied_voltage
        + sigma_avg * mu * EC * (d->zeta_left + d->zeta_right) * 0.5
          * d->applied_voltage / (KB * 300.0));
}

/* ---- Ion Selectivity ---- */

typedef struct {
    double channel_height;
    double debye_length;
    double surface_charge_density;
    double cation_mobility, anion_mobility;
    double bulk_concentration;
} IonSelectiveChannel;

double ion_selectivity_cation_transference(const IonSelectiveChannel *ch)
{
    if (!ch) return 0.5;
    double h = ch->channel_height;
    double lambda_D = ch->debye_length;
    if (lambda_D <= 0.0) return 0.5;

    double Dukhin = fabs(ch->surface_charge_density) / (EC * ch->bulk_concentration * 6.022e23 * h);
    double sigma_star = Dukhin * lambda_D / h;

    double mu_plus = ch->cation_mobility, mu_minus = ch->anion_mobility;
    double c0 = ch->bulk_concentration * 6.022e23;

    double excess_cation = sigma_star * c0;
    double I_plus = mu_plus * (c0 + excess_cation);
    double I_minus = mu_minus * c0;

    return I_plus / (I_plus + I_minus);
}

double debye_length_from_concentration(double concentration_molar, double permittivity, double T)
{
    if (concentration_molar <= 0.0) return -1.0;
    double c0 = concentration_molar * 6.022e23 * 1000.0;
    return sqrt(permittivity * 8.854e-12 * KB * T / (2.0 * c0 * EC * EC));
}

/* ---- DNA Translocation ---- */

typedef struct {
    double pore_radius, pore_length;
    double dna_charge_per_nucleotide;
    int n_nucleotides;
    double applied_voltage;
    double solution_viscosity, temperature;
} DNATranslocation;

double dna_translocation_velocity(const DNATranslocation *dna)
{
    if (!dna || dna->pore_radius <= 0.0) return 0.0;
    double Q_total = dna->dna_charge_per_nucleotide * dna->n_nucleotides * EC;
    double E_field = dna->applied_voltage / dna->pore_length;
    double F_electric = Q_total * E_field;
    double r_eff = dna->pore_radius * 0.5;
    double F_drag = 6.0 * M_PI * dna->solution_viscosity * r_eff;
    if (F_drag <= 0.0) return 0.0;
    return F_electric / F_drag;
}

double dna_translocation_time(const DNATranslocation *dna)
{
    double v = dna_translocation_velocity(dna);
    if (v <= 0.0) return -1.0;
    return dna->pore_length / v;
}

double dna_blockade_current_drop(const DNATranslocation *dna, double open_pore_current)
{
    if (!dna || dna->pore_radius <= 0.0) return 0.0;
    double A_pore = M_PI * dna->pore_radius * dna->pore_radius;
    double r_dna = 1.1e-9;
    double A_dna = M_PI * r_dna * r_dna;
    if (A_dna >= A_pore) return open_pore_current;
    return open_pore_current * A_dna / A_pore;
}

/* ---- Lab-on-Chip Transport ---- */

typedef struct {
    double channel_length, channel_width, channel_height;
    double inlet_concentration, flow_rate;
    double diffusivity, reaction_rate_constant;
} LabOnChipReactor;

double lab_on_chip_conversion(const LabOnChipReactor *r)
{
    if (!r || r->flow_rate <= 0.0) return 0.0;
    double A = r->channel_width * r->channel_height;
    double u = r->flow_rate / A;
    if (u <= 0.0) return 0.0;
    double tau = r->channel_length / u;
    double Da = r->reaction_rate_constant * tau;
    double Pe = u * r->channel_height / r->diffusivity;
    double dispersion = 1.0 + Pe * Pe / 210.0;
    double Da_eff = Da / dispersion;
    return 1.0 - exp(-Da_eff);
}

double lab_on_chip_peclet_number(const LabOnChipReactor *r)
{
    if (!r || r->diffusivity <= 0.0) return 0.0;
    double A = r->channel_width * r->channel_height;
    double u = r->flow_rate / A;
    return u * r->channel_height / r->diffusivity;
}

double lab_on_chip_pressure_drop(const LabOnChipReactor *r, double viscosity)
{
    if (!r) return 0.0;
    double h = r->channel_height, w = r->channel_width;
    double Dh = 2.0 * h * w / (h + w);
    double A = w * h;
    double u = r->flow_rate / A;
    return 12.0 * viscosity * u * r->channel_length / (h * h);
}

/* ---- Electrokinetic Energy Conversion ---- */

typedef struct {
    double channel_height, zeta_potential;
    double permittivity, viscosity, conductivity;
    double pressure_gradient;
} ElectrokineticConverter;

double streaming_current(const ElectrokineticConverter *ek)
{
    if (!ek || ek->channel_height <= 0.0) return 0.0;
    double eps = ek->permittivity * 8.854e-12;
    double h = ek->channel_height;
    return -eps * ek->zeta_potential * h * ek->pressure_gradient / ek->viscosity;
}

double streaming_potential_gradient(const ElectrokineticConverter *ek)
{
    double I_str = streaming_current(ek);
    if (fabs(ek->conductivity) < 1e-18) return 0.0;
    return -I_str / (ek->conductivity * ek->channel_height);
}

double electrokinetic_conversion_efficiency(const ElectrokineticConverter *ek)
{
    if (!ek || ek->pressure_gradient <= 0.0) return 0.0;
    double I_str = streaming_current(ek);
    double dP = ek->pressure_gradient;
    double Q = ek->channel_height * ek->channel_height * ek->channel_height * dP
               / (12.0 * ek->viscosity);
    double P_mech = Q * dP;
    double E_str = streaming_potential_gradient(ek);
    double P_elec = I_str * I_str / ek->conductivity;
    return P_mech > 0.0 ? P_elec / P_mech : 0.0;
}
/* ---- Nanofluidic Energy Harvester ---- */

typedef struct {
    double membrane_thickness, pore_diameter, porosity;
    double surface_charge_density, salt_concentration;
    double temperature, pressure_difference;
} NanofluidicEnergyHarvester;

double nanofluidic_harvester_power_density(const NanofluidicEnergyHarvester *h)
{
    if (!h || h->membrane_thickness <= 0.0) return 0.0;
    double eps = 80.0 * 8.854e-12;
    double lambda_D = debye_length_from_concentration(h->salt_concentration, 80.0, h->temperature);
    double sigma = h->surface_charge_density;
    double E_str = sigma * lambda_D * h->pressure_difference / (h->salt_concentration * 6.022e23 * 1000.0 * EC);
    double I_str = sigma * EC * h->porosity * h->pressure_difference * h->pore_diameter
                   / (h->membrane_thickness * 1.0e-3);
    return I_str * I_str / (h->porosity * h->membrane_thickness * 1e-4);
}

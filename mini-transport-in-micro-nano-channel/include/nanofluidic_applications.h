/**
 * @file nanofluidic_applications.h
 * @brief Nanofluidic diode, ion selectivity, DNA translocation, lab-on-chip.
 */
#ifndef NANOFLUIDIC_APPLICATIONS_H
#define NANOFLUIDIC_APPLICATIONS_H

typedef struct {
    double channel_length, channel_width, channel_height;
    double surface_charge_density_left, surface_charge_density_right;
    double bulk_concentration, zeta_left, zeta_right;
    double applied_voltage;
} NanofluidicDiode;

double nanofluidic_diode_rectification_ratio(const NanofluidicDiode *d);
double nanofluidic_diode_forward_current(const NanofluidicDiode *d);

typedef struct {
    double channel_height, debye_length;
    double surface_charge_density;
    double cation_mobility, anion_mobility;
    double bulk_concentration;
} IonSelectiveChannel;

double ion_selectivity_cation_transference(const IonSelectiveChannel *ch);
double debye_length_from_concentration(double c, double permittivity, double T);

typedef struct {
    double pore_radius, pore_length;
    double dna_charge_per_nucleotide;
    int n_nucleotides;
    double applied_voltage;
    double solution_viscosity, temperature;
} DNATranslocation;

double dna_translocation_velocity(const DNATranslocation *dna);
double dna_translocation_time(const DNATranslocation *dna);
double dna_blockade_current_drop(const DNATranslocation *dna, double open_pore_current);

typedef struct {
    double channel_length, channel_width, channel_height;
    double inlet_concentration, flow_rate;
    double diffusivity, reaction_rate_constant;
} LabOnChipReactor;

double lab_on_chip_conversion(const LabOnChipReactor *r);
double lab_on_chip_peclet_number(const LabOnChipReactor *r);
double lab_on_chip_pressure_drop(const LabOnChipReactor *r, double viscosity);

typedef struct {
    double channel_height, zeta_potential;
    double permittivity, viscosity, conductivity;
    double pressure_gradient;
} ElectrokineticConverter;

double streaming_current(const ElectrokineticConverter *ek);
double streaming_potential_gradient(const ElectrokineticConverter *ek);
double electrokinetic_conversion_efficiency(const ElectrokineticConverter *ek);

#endif
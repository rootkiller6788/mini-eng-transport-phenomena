#ifndef NANOCHANNEL_ION_TRANSPORT_H
#define NANOCHANNEL_ION_TRANSPORT_H

#include "micro_nano_transport.h"

/*
 * nanochannel_ion_transport.h - Ionic Transport in Nanochannels
 *
 * Specialized types and API for ion transport in nanoscale channels
 * where electric double layer overlap, surface charge dominance,
 * and ion selectivity are the governing phenomena.
 */

/* Ion transport regime classification */
typedef enum {
    ION_REGIME_BULK = 0,
    ION_REGIME_SURFACE_CONDUCTION,
    ION_REGIME_EDL_OVERLAP,
    ION_REGIME_SINGLE_FILE,
    ION_REGIME_COULOMB_BLOCKADE
} IonTransportRegime;

/* Ion species properties */
typedef struct {
    double concentration;
    double diffusion_coefficient;
    int valence;
    double hydrodynamic_radius;
    double hydrated_radius;
    char name[32];
} IonSpecies;

/* Nanopore sensing parameters */
typedef struct {
    double pore_radius;
    double pore_length;
    double membrane_thickness;
    double surface_charge_density;
    double access_resistance;
    double channel_resistance;
    double total_resistance;
    double open_pore_current;
    double blockade_current;
    double signal_to_noise_ratio;
    double bandwidth;
    char pore_material[64];
} NanoporeSensor;

/* Ion current rectification state */
typedef struct {
    double current_forward;
    double current_reverse;
    double rectification_ratio;
    double diode_quality_factor;
    double built_in_potential;
    int rectifying;
} IonRectification;

/* API for nanochannel ion transport */
double compute_ionic_current_nanochannel(int n_species, const double *conc, const int *valence, const double *D, double V_applied, double L, double A, double T);
double compute_nanochannel_conductance(double c_bulk, double z, double D, double L, double A, double sigma_s, double H, double T);
double compute_surface_charge_potential(double sigma_s, double c_bulk, double z, double eps, double T);
double compute_donnan_potential(double sigma_s, double c_bulk, double z, double H, double T);
double compute_ion_selectivity(int n_species, const double *conc, const int *valence, const double *D);
double compute_ion_rectification_ratio(double sigma_in, double sigma_out, double z, double T);
double compute_nanopore_resistance(double L, double r, double sigma_bulk);
double compute_nanopore_blockade_current(double I_open, double r_pore, double d_molecule);
double compute_overlap_parameter(double kappa, double H);
double compute_coion_exclusion_coefficient(double psi_center, int z_co, double T);
void compute_ion_concentration_profile(const double *psi, int n, double c_bulk, int z, double T, double *c_profile);

/* Ion transport advanced */
double compute_ion_pairing_constant(double z_plus, double z_minus, double eps, double T, double a_separation);
double compute_activity_coefficient_debye_huckel(double z, double I, double eps, double T, double a_ion);
double compute_conductivity_onsager(double conc, double z, double D, double eps, double eta, double T);
double compute_electroosmotic_conductivity(double sigma_s, double H, double mu, double c_bulk);

#endif

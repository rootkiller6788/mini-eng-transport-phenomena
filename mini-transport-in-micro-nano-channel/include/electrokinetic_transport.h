#ifndef ELECTROKINETIC_TRANSPORT_H
#define ELECTROKINETIC_TRANSPORT_H

#include "micro_nano_transport.h"

/*
 * electrokinetic_transport.h - Electrokinetic Phenomena in Channels
 *
 * Types and API for electroosmosis, electrophoresis, streaming potential,
 * and related electrokinetic phenomena in micro/nano channels.
 */

/* Electrokinetic phenomenon type */
typedef enum {
    EK_ELECTROOSMOSIS = 0,
    EK_ELECTROPHORESIS,
    EK_STREAMING_POTENTIAL,
    EK_STREAMING_CURRENT,
    EK_SEDIMENTATION_POTENTIAL,
    EK_DIELECTROPHORESIS,
    EK_AC_ELECTROOSMOSIS,
    EK_ICEO
} ElectrokineticType;

/* Electric double layer model */
typedef enum {
    EDL_MODEL_HELMHOLTZ = 0,
    EDL_MODEL_GOUY_CHAPMAN,
    EDL_MODEL_STERN,
    EDL_MODEL_GRAHAME,
    EDL_MODEL_BOCKRIS,
    EDL_MODEL_MEAN_SPHERICAL
} EDLModel;

/* Electrokinetic transport state */
typedef struct {
    ElectrokineticType ek_type;
    EDLModel edl_model;
    double eof_velocity;
    double eof_mobility;
    double streaming_potential;
    double streaming_current;
    double electrophoretic_velocity;
    double electrophoretic_mobility;
    double dielectrophoretic_force;
    double zeta_potential_effective;
    double surface_conductance;
    double dukhin_number;
    double peclet_number_ek;
    int edl_overlap;
} ElectrokineticState;

/* Extended EDL parameters */
typedef struct {
    double debye_length;
    double gouy_chapman_length;
    double bjerrum_length;
    double stern_layer_thickness;
    double inner_helmholtz_capacitance;
    double outer_helmholtz_capacitance;
    double diffuse_layer_capacitance;
    double total_capacitance;
    double surface_charge_regulation;
    double pzc;  /* potential of zero charge */
    double isoelectric_point;
} EDLParameters;

/* Electrokinetic API */
double compute_helmholtz_smoluchowski_velocity(double eps, double zeta, double eta, double E);
double compute_eof_mobility(double eps, double zeta, double eta);
void solve_electroosmotic_flow(double H, double eps, double zeta, double eta, double Ex, double kappa, int n, double *y, double *u);
double compute_eof_flow_rate(double H, double W, double eps, double zeta, double eta, double Ex, double kappa);
double compute_eof_with_backpressure(double H, double W, double dpdx, double eps, double zeta, double eta, double Ex, double kappa);
double compute_max_eof_backpressure(double H, double L, double eps, double zeta, double eta, double Ex);
void compute_eof_profile_full(double H, double eps, double zeta, double eta, double Ex, double kappa, double dpdx, int n, double *y, double *u);
double compute_edl_thickness_ratio(double kappa, double H);
double compute_eof_power_consumption(double sigma_bulk, double A, double L, double V);
double compute_eof_efficiency(double Q, double delta_p, double current, double voltage);

/* Streaming potential API */
void compute_streaming_potential(double eps, double zeta, double eta, double sigma_bulk, double dpdx, double *V_str);
void compute_zeta_from_streaming(double V_str, double dp, double eps, double eta, double sigma, double *zeta);
double compute_streaming_current_full(double eps, double zeta, double eta, double H, double W, double dpdx, double kappa);
double compute_surface_conductance(double sigma_s, double mu_counterion, double H);

/* Advanced EDL functions */
double compute_differential_capacitance(double psi0, double kappa, double eps, double z, double T);
void compute_stern_layer_correction(double sigma_s, double d_stern, double eps_stern, double zeta_in, double *psi_0, double *psi_stern);
double compute_effective_screening_length(double kappa, double sigma_s, double eps, double T);
double compute_overscreening_parameter(double sigma_s, double c_bulk, double z, double eps, double T);

#endif

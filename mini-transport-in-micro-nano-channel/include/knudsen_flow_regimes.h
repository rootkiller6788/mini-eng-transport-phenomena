#ifndef KNUDSEN_FLOW_REGIMES_H
#define KNUDSEN_FLOW_REGIMES_H

#include "micro_nano_transport.h"

/*
 * knudsen_flow_regimes.h - Knudsen Number Classification Utilities
 *
 * This header declares functions for determining flow regime,
 * computing rarefaction parameters, and assessing the validity
 * of the continuum hypothesis in micro/nano channels.
 */

/* Knudsen classification thresholds */
#define KN_CONTINUUM_THRESHOLD    1.0e-3
#define KN_SLIP_FLOW_THRESHOLD    0.1
#define KN_TRANSITION_THRESHOLD   10.0

/* von Karman relation structure */
typedef struct {
    double Reynolds;
    double Mach;
    double Kn_estimated;
    double delta_estimated;
    int regime_consistent;
} VonKarmanRelation;

/* Multi-regime transport model selector */
typedef enum {
    MODEL_NAVIER_STOKES = 0,
    MODEL_NAVIER_STOKES_SLIP,
    MODEL_BURNETT,
    MODEL_BGK,
    MODEL_DSMC,
    MODEL_MOLECULAR_DYNAMICS,
    MODEL_FREE_MOLECULAR
} TransportModel;

/* Transport model recommendation based on Kn */
typedef struct {
    KnudsenRegime regime;
    TransportModel recommended_model;
    double expected_error_percent;
    const char *model_description;
    int requires_kinetic_solver;
} ModelRecommendation;

/* Extended regime information with sub-regimes */
typedef struct {
    KnudsenRegime primary_regime;
    int sub_regime;
    char description[256];
    double inverse_kn;
    double rarefaction_coefficient;
    int continuum_valid;
    int slip_correction_needed;
    int dsmc_recommended;
} RegimeInfo;

/* API for regime classification */
int validate_continuum_hypothesis(const KnudsenState *kn, const FluidProperties *fluid, double pressure);
double compute_mean_free_path_vhs(FluidProperties *fluid, double pressure, double omega, double ref_T);
double compute_effective_transport_length(const ChannelGeometry *geom, const KnudsenState *kn, double alpha, double beta);
int classify_flow_regime_detailed(const KnudsenState *kn, char *name, size_t nlen, double *inv_kn);
double compute_rarefaction_parameter(const KnudsenState *kn);
void compare_reynolds_knudsen(double Re, double Ma, double gamma, double *kn, double *delta);
KnudsenRegime estimate_flow_regime_from_re_ma(double Re, double Ma, double gamma);
double compute_collision_frequency(double T, double p, double m_mol, double d_col);
double compute_characteristic_length_automatic(const ChannelGeometry *geom);

/* Model recommendation */
void recommend_transport_model(const KnudsenState *kn, ModelRecommendation *rec);
double compute_knudsen_layer_correction(double Kn, double slip_len, double H);
int is_continuum_valid(const KnudsenState *kn);
double compute_effective_viscosity_rarefied(double eta, double Kn);
double compute_thermal_creep_coefficient(double Kn, double sigma_T);

#endif

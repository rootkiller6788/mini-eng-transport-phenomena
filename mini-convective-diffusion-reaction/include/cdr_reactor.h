/**
 * @file cdr_reactor.h
 * @brief Ideal and Non-Ideal Chemical Reactor Models
 *
 * CSTR (Continuous Stirred Tank Reactor), PFR (Plug Flow Reactor),
 * Batch reactor, Dispersion model, Tanks-in-Series model,
 * and Residence Time Distribution (RTD) analysis.
 *
 * Reference: Levenspiel "Chemical Reaction Engineering" (1999), Ch. 5-15
 * Reference: Fogler "Elements of Chemical Reaction Engineering" (2016), Ch. 1-6, 13-18
 *
 * L1: CSTR, PFR, Batch reactor definitions, space time, space velocity
 * L2: Ideal reactor behavior, RTD concepts
 * L4: Design equations (mass balances for each reactor type)
 * L5: RTD analysis, reactor design algorithms
 * L6: Reactor sizing problems, conversion optimization
 */

#ifndef CDR_REACTOR_H
#define CDR_REACTOR_H

#include "cdr_core.h"
#include "cdr_reaction.h"

/* ---------------------------------------------------------------------------
 * L1: Reactor Type Definitions
 * ------------------------------------------------------------------------- */

/**
 * @brief Reactor types in chemical engineering.
 */
typedef enum {
    REACTOR_BATCH,           /**< Closed system, time-dependent */
    REACTOR_CSTR,            /**< Continuous stirred tank, steady-state mixing */
    REACTOR_PFR,             /**< Plug flow, no axial mixing, steady-state */
    REACTOR_PBR,             /**< Packed bed reactor (heterogeneous catalytic) */
    REACTOR_SEMIBATCH,       /**< One reactant charged, other fed continuously */
    REACTOR_DISPERSION,      /**< Non-ideal PFR with axial dispersion */
    REACTOR_TANKS_IN_SERIES  /**< N CSTRs in series */
} ReactorType;

/**
 * @brief Reactor operating conditions and design specifications.
 */
typedef struct {
    ReactorType  type;           /**< Reactor type */
    double       volume;         /**< Reactor volume [m^3] */
    double       volumetric_flow;/**< Volumetric flow rate Q [m^3/s] */
    double       space_time;     /**< tau = V/Q [s] */
    double       space_velocity; /**< SV = 1/tau = Q/V [1/s] */
    double       C_A0;           /**< Inlet concentration of A [mol/m^3] */
    double       C_B0;           /**< Inlet concentration of B [mol/m^3] */
    double       T;              /**< Operating temperature [K] */
    double       conversion;     /**< Target or achieved conversion X_A [0-1] */
    double       pressure_drop;  /**< Pressure drop [Pa] */
    int          is_isothermal;  /**< 1 = isothermal, 0 = adiabatic/non-isothermal */
} ReactorSpec;

/* ---------------------------------------------------------------------------
 * L2: Residence Time Distribution (RTD)
 * ------------------------------------------------------------------------- */

/**
 * @brief RTD function E(t) ? exit age distribution.
 *
 * E(t) dt = fraction of fluid leaving with residence time in [t, t+dt]
 *
 * Normalization: integral_0^inf E(t) dt = 1
 */
typedef struct {
    double *time;        /**< Time points [s] */
    double *E;           /**< E(t) values [1/s] */
    double *F;           /**< Cumulative F(t) = integral_0^t E(t') dt' */
    double  mean_t;      /**< Mean residence time t_bar = integral t*E(t) dt */
    double  variance;    /**< Variance sigma^2 */
    double  skewness;    /**< Skewness (asymmetry measure) */
    size_t  n_points;    /**< Number of RTD data points */
} ResidenceTimeDistribution;

/**
 * @brief RTD model types for reactor characterization.
 */
typedef enum {
    RTD_MODEL_CSTR,              /**< E(t) = (1/tau)*exp(-t/tau) */
    RTD_MODEL_PFR,               /**< E(t) = delta(t - tau) */
    RTD_MODEL_DISPERSION,        /**< Axial dispersion model */
    RTD_MODEL_TANKS_IN_SERIES,   /**< N equal CSTRs in series */
    RTD_MODEL_CSTR_DEAD_ZONE,    /**< CSTR with stagnant volume */
    RTD_MODEL_CSTR_BYPASS        /**< CSTR with bypass/short-circuiting */
} RTDModelType;

/* ---------------------------------------------------------------------------
 * L4: Reactor Design Equations (Mass Balances)
 * ------------------------------------------------------------------------- */

/**
 * @brief CSTR design equation for first-order irreversible reaction.
 *
 * V = Q * (C_A0 - C_A) / (k * C_A)
 * C_A = C_A0 / (1 + k * tau)
 * X_A = k * tau / (1 + k * tau)
 *
 * @param C_A0    Inlet concentration [mol/m^3]
 * @param Q       Volumetric flow rate [m^3/s]
 * @param k       First-order rate constant [1/s]
 * @param X_target Target conversion [0-1]
 * @return Required CSTR volume [m^3]
 *
 * Complexity: O(1)
 */
double cdr_cstr_volume_first_order(double C_A0, double Q, double k,
                                    double X_target);

/**
 * @brief CSTR outlet concentration (first-order, irreversible).
 *
 * @param C_A0  Inlet concentration [mol/m^3]
 * @param tau   Space time V/Q [s]
 * @param k     Rate constant [1/s]
 * @return Outlet concentration C_A [mol/m^3]
 */
double cdr_cstr_outlet_first_order(double C_A0, double tau, double k);

/**
 * @brief CSTR design equation for second-order (2A -> products).
 *
 * V = Q * (C_A0 - C_A) / (k * C_A^2)
 * Requires solving quadratic: k*tau*C_A^2 + C_A - C_A0 = 0
 *
 * @param C_A0     Inlet concentration [mol/m^3]
 * @param Q        Volumetric flow rate [m^3/s]
 * @param k        Second-order rate constant [m^3/(mol.s)]
 * @param X_target Target conversion [0-1]
 * @return Required CSTR volume [m^3]
 */
double cdr_cstr_volume_second_order(double C_A0, double Q, double k,
                                     double X_target);

/**
 * @brief PFR design equation for first-order irreversible reaction.
 *
 * V = Q * ln(1/(1-X)) / k  [constant density]
 * C_A(z) = C_A0 * exp(-k * tau * z/L)
 *
 * @param C_A0     Inlet concentration [mol/m^3]
 * @param Q        Volumetric flow rate [m^3/s]
 * @param k        First-order rate constant [1/s]
 * @param X_target Target conversion [0-1]
 * @return Required PFR volume [m^3]
 */
double cdr_pfr_volume_first_order(double C_A0, double Q, double k,
                                   double X_target);

/**
 * @brief PFR design equation for second-order (2A -> products).
 *
 * V = Q * (X / (1-X)) / (k * C_A0)  [constant density]
 *
 * @param C_A0     Inlet concentration [mol/m^3]
 * @param Q        Volumetric flow rate [m^3/s]
 * @param k        Second-order rate constant [m^3/(mol.s)]
 * @param X_target Target conversion [0-1]
 * @return Required PFR volume [m^3]
 */
double cdr_pfr_volume_second_order(double C_A0, double Q, double k,
                                    double X_target);

/**
 * @brief PFR concentration profile along reactor length.
 *
 * C_A(z) for first-order reaction in PFR.
 *
 * @param C_A0  Inlet concentration [mol/m^3]
 * @param k     Rate constant [1/s]
 * @param tau   Space time V/Q [s]
 * @param z_L   Fractional position z/L [0-1]
 * @return Concentration at position z [mol/m^3]
 */
double cdr_pfr_profile_first_order(double C_A0, double k, double tau, double z_L);

/**
 * @brief Batch reactor time to reach conversion X (first-order).
 *
 * t = -ln(1-X) / k
 *
 * @param k  Rate constant [1/s]
 * @param X  Target conversion [0-1)
 * @return Batch time [s]
 */
double cdr_batch_time_first_order(double k, double X);

/**
 * @brief Batch reactor time for second-order (2A -> products).
 *
 * t = X / (k * C_A0 * (1-X))
 *
 * @param k     Rate constant [m^3/(mol.s)]
 * @param C_A0  Initial concentration [mol/m^3]
 * @param X     Target conversion [0-1)
 * @return Batch time [s]
 */
double cdr_batch_time_second_order(double k, double C_A0, double X);

/**
 * @brief Compare CSTR vs PFR volume for same conversion (first-order).
 *
 * Ratio V_CSTR / V_PFR = X / ((1-X) * ln(1/(1-X)))
 *
 * @param X  Conversion [0-1)
 * @return Volume ratio V_CSTR / V_PFR (>= 1 for positive-order kinetics)
 */
double cdr_cstr_pfr_volume_ratio(double X);

/**
 * @brief PBR (Packed Bed Reactor) catalyst weight (first-order).
 *
 * W = Q * ln(1/(1-X)) / k'   where k' is mass-based rate constant [m^3/(kg.s)]
 *
 * @param Q        Volumetric flow rate [m^3/s]
 * @param k_prime  Mass-based rate constant [m^3/(kg_cat.s)]
 * @param X_target Target conversion [0-1)
 * @return Required catalyst mass [kg]
 */
double cdr_pbr_catalyst_mass_first_order(double Q, double k_prime,
                                          double X_target);

/* ---------------------------------------------------------------------------
 * L5: RTD Analysis Methods
 * ------------------------------------------------------------------------- */

/**
 * @brief E(t) for an ideal CSTR.
 *
 * E(t) = (1/tau) * exp(-t/tau)
 *
 * @param t    Time [s]
 * @param tau  Mean residence time [s]
 * @return E(t) value [1/s]
 */
double cdr_rtd_cstr_E(double t, double tau);

/**
 * @brief E(t) for N equal CSTRs in series.
 *
 * E(t) = (N/tau)^N * t^(N-1) * exp(-N*t/tau) / (N-1)!
 *
 * @param t    Time [s]
 * @param tau  Total mean residence time [s]
 * @param N    Number of tanks (>= 1)
 * @return E(t) value [1/s]
 */
double cdr_rtd_tanks_in_series_E(double t, double tau, int N);

/**
 * @brief E(t) for axial dispersion model (open-open boundary).
 *
 * Uses approximate analytical form.
 *
 * @param t    Time [s]
 * @param tau  Mean residence time [s]
 * @param Pe   Axial Peclet number Pe = u*L/D_ax
 * @return E(t) value [1/s]
 */
double cdr_rtd_dispersion_E(double t, double tau, double Pe);

/**
 * @brief Calculate mean residence time from experimental RTD data.
 *
 * t_bar = integral t*E(t) dt / integral E(t) dt
 *
 * Uses trapezoidal integration.
 *
 * @param time      Time points [s] (size n)
 * @param E         E(t) values [1/s] (size n)
 * @param n_points  Number of data points
 * @return Mean residence time [s]
 */
double cdr_rtd_mean_time(const double *time, const double *E, size_t n_points);

/**
 * @brief Calculate RTD variance from experimental data.
 *
 * sigma^2 = integral (t - t_bar)^2 * E(t) dt
 *
 * @param time      Time points [s]
 * @param E         E(t) values [1/s]
 * @param n_points  Number of data points
 * @param t_bar     Mean residence time (pre-computed)
 * @return Variance [s^2]
 */
double cdr_rtd_variance(const double *time, const double *E,
                         size_t n_points, double t_bar);

/**
 * @brief Estimate dispersion number from RTD variance.
 *
 * sigma_theta^2 = sigma^2 / tau^2 = 2/Pe - 2/Pe^2 * (1 - exp(-Pe))
 *
 * For open-open vessels. Uses numerical inversion.
 *
 * @param variance    RTD variance [s^2]
 * @param tau         Mean residence time [s]
 * @return Estimated Peclet number Pe [-]
 */
double cdr_rtd_peclet_from_variance(double variance, double tau);

/**
 * @brief Estimate number of tanks from RTD variance.
 *
 * N = tau^2 / sigma^2  (for tanks-in-series model)
 *
 * @param tau      Mean residence time [s]
 * @param variance RTD variance [s^2]
 * @return Estimated number of tanks (N >= 1)
 */
double cdr_rtd_N_tanks_from_variance(double tau, double variance);

/**
 * @brief Conversion in a real reactor using segregation model.
 *
 * X_bar = integral_0^inf X_batch(t) * E(t) dt
 *
 * where X_batch(t) is the conversion in a batch reactor after time t.
 *
 * Uses trapezoidal integration for first-order reaction.
 *
 * @param time      RTD time points [s]
 * @param E         E(t) values [1/s]
 * @param n_points  Number of points
 * @param k         First-order rate constant [1/s]
 * @return Predicted conversion [-]
 */
double cdr_segregation_model_conversion(const double *time, const double *E,
                                         size_t n_points, double k);

/* ---------------------------------------------------------------------------
 * L5: Dispersion Model (Non-Ideal PFR)
 * ------------------------------------------------------------------------- */

/**
 * @brief Solve 1D convection-diffusion-reaction steady-state (dispersion model).
 *
 * D_ax * d^2C/dz^2 - u * dC/dz + r(C) = 0
 *
 * For first-order reaction, analytical solution with Danckwerts BC:
 *
 * C/C0 = 4*q*exp(Pe/2) / ((1+q)^2*exp(q*Pe/2) - (1-q)^2*exp(-q*Pe/2))
 *
 * where q = sqrt(1 + 4*Da/Pe), Da = k*tau, Pe = u*L/D_ax
 *
 * @param C_A0  Inlet concentration [mol/m^3]
 * @param Pe    Axial Peclet number [-]
 * @param Da    Damkohler number Da_I = k*tau [-]
 * @return Outlet concentration [mol/m^3]
 */
double cdr_dispersion_outlet_concentration(double C_A0, double Pe, double Da);

/**
 * @brief Dispersion model conversion.
 *
 * X = 1 - C_out/C_A0
 *
 * @param Pe  Axial Peclet number [-]
 * @param Da  Damkohler-I number [-]
 * @return Conversion [0-1]
 */
double cdr_dispersion_conversion(double Pe, double Da);

/* ---------------------------------------------------------------------------
 * L5: Reactor Network Optimization
 * ------------------------------------------------------------------------- */

/**
 * @brief Optimal intermediate conversion for two CSTRs in series (first-order).
 *
 * For N CSTRs in series, optimal equal sizing gives:
 * X_final = 1 - 1/(1 + k*tau_i)^N  where tau_i = V_i/Q
 *
 * @param k           Rate constant [1/s]
 * @param total_tau   Total space time [s]
 * @param N           Number of equal CSTRs
 * @return Final conversion [-]
 */
double cdr_cstr_series_conversion(double k, double total_tau, int N);

/**
 * @brief PFR with recycle, first-order reaction.
 *
 * X_final = ...  (solving the recycle reactor equation)
 *
 * @param X_target   Target overall conversion
 * @param k          Rate constant [1/s]
 * @param tau        Space time based on fresh feed [s]
 * @param R          Recycle ratio (recycle flow / fresh feed flow)
 * @return Required volume factor relative to PFR without recycle
 */
double cdr_pfr_recycle_volume_factor(double X_target, double k,
                                      double tau, double R);

#endif /* CDR_REACTOR_H */

/**
 * @file cdr_reaction.h
 * @brief Chemical Reaction Kinetics for CDR Systems
 *
 * Reaction rate laws (zero through nth-order), Arrhenius temperature dependence,
 * Michaelis-Menten enzyme kinetics, Langmuir-Hinshelwood surface kinetics,
 * and complex reaction networks (series, parallel, reversible).
 *
 * Reference: Levenspiel "Chemical Reaction Engineering" (1999), Ch. 2-7
 * Reference: Fogler "Elements of Chemical Reaction Engineering" (2016), Ch. 3-10
 *
 * L1: Reaction order, rate constant, activation energy, Arrhenius parameters
 * L2: Rate-limiting step, rate-determining step concepts
 * L4: Arrhenius law, mass-action kinetics
 * L5: Rate law fitting, half-life analysis
 */

#ifndef CDR_REACTION_H
#define CDR_REACTION_H

#include "cdr_core.h"

/* ---------------------------------------------------------------------------
 * L1: Reaction Kinetics Structures
 * ------------------------------------------------------------------------- */

/**
 * @brief Reaction order specification.
 */
typedef enum {
    REACTION_ZERO_ORDER   = 0,
    REACTION_FIRST_ORDER  = 1,
    REACTION_SECOND_ORDER = 2,
    REACTION_THIRD_ORDER  = 3,
    REACTION_FRACTIONAL_ORDER = -1  /**< Non-integer, specified separately */
} ReactionOrder;

/**
 * @brief Arrhenius parameters for temperature-dependent rate constants.
 *
 * k(T) = A * exp(-Ea / (R * T))
 *
 * A  = pre-exponential (frequency) factor [units depend on order]
 * Ea = activation energy [J/mol]
 */
typedef struct {
    double A;            /**< Pre-exponential factor [(m^3/mol)^(n-1)/s] */
    double Ea;           /**< Activation energy [J/mol] */
    double T_ref;        /**< Reference temperature [K] for k_ref variant */
    double k_ref;        /**< Reference rate constant at T_ref */
} ArrheniusParams;

/**
 * @brief General rate law parameters for power-law kinetics.
 *
 * r = k * C_A^alpha * C_B^beta  (irreversible)
 * r = k_f * C_A^alpha * C_B^beta - k_r * C_C^gamma * C_D^delta  (reversible)
 */
typedef struct {
    double k_forward;     /**< Forward rate constant */
    double k_reverse;     /**< Reverse rate constant (0 for irreversible) */
    double alpha;         /**< Reaction order w.r.t. species A */
    double beta;          /**< Reaction order w.r.t. species B */
    double gamma;         /**< Reaction order w.r.t. species C (reverse) */
    double delta;         /**< Reaction order w.r.t. species D (reverse) */
    double Keq;           /**< Equilibrium constant = k_f / k_r */
    int    is_reversible; /**< 1 if reversible, 0 if irreversible */
} PowerLawRate;

/**
 * @brief Michaelis-Menten enzyme kinetics parameters.
 *
 * E + S <-> ES -> E + P
 *
 * r = V_max * [S] / (K_m + [S])
 *
 * V_max = k_cat * [E]_0  (maximum rate)
 * K_m   = (k_-1 + k_cat) / k_1  (Michaelis constant)
 */
typedef struct {
    double V_max;     /**< Maximum reaction rate [mol/(m^3.s)] */
    double K_m;       /**< Michaelis constant [mol/m^3] */
    double k_cat;     /**< Catalytic rate constant (turnover number) [1/s] */
    double E_total;   /**< Total enzyme concentration [mol/m^3] */
} MichaelisMentenParams;

/**
 * @brief Langmuir-Hinshelwood rate parameters for surface reactions.
 *
 * r = k * K_A * P_A * K_B * P_B / (1 + K_A*P_A + K_B*P_B)^2  (bimolecular)
 *
 * For unimolecular: r = k * K_A * P_A / (1 + K_A*P_A)
 */
typedef struct {
    double k_surface;     /**< Surface reaction rate constant [mol/(m^2.s)] */
    double K_ads_A;       /**< Adsorption equilibrium constant A [1/Pa] */
    double K_ads_B;       /**< Adsorption equilibrium constant B [1/Pa] */
    double K_ads_C;       /**< Adsorption constant for product (if inhibiting) */
    int    n_sites;       /**< Number of active site types */
    int    is_dissociative_A; /**< 1 if A dissociates on adsorption */
} LangmuirHinshelwoodParams;

/**
 * @brief Reaction network type classification.
 */
typedef enum {
    NET_SINGLE,            /**< A -> B, single irreversible */
    NET_SERIES,            /**< A -> B -> C, consecutive */
    NET_PARALLEL,          /**< A -> B, A -> C, competing pathways */
    NET_REVERSIBLE,        /**< A <-> B */
    NET_SERIES_PARALLEL,   /**< A -> B -> C, A -> D (mixed) */
    NET_AUTOCATALYTIC      /**< A + B -> 2B */
} ReactionNetworkType;

/**
 * @brief Concentration state for multiple species.
 */
typedef struct {
    double *C;          /**< Concentrations [mol/m^3], length = n_species */
    size_t  n_species;  /**< Number of species */
    double  time;       /**< Reaction time [s] */
} ConcentrationState;

/* ---------------------------------------------------------------------------
 * L4: Arrhenius Rate Law
 * ------------------------------------------------------------------------- */

/**
 * @brief Compute rate constant from Arrhenius law.
 *
 * k(T) = A * exp(-Ea / (R * T))
 *
 * @param A   Pre-exponential factor [units depend on order]
 * @param Ea  Activation energy [J/mol]
 * @param T   Temperature [K]
 * @return Rate constant k(T)
 *
 * R = 8.314 J/(mol.K), Complexity: O(1)
 */
double cdr_arrhenius_rate(double A, double Ea, double T);

/**
 * @brief Compute rate constant from Arrhenius law using reference state.
 *
 * k(T) = k_ref * exp( (Ea/R) * (1/T_ref - 1/T) )
 *
 * Numerically more stable when k_ref and T_ref are known.
 *
 * @param k_ref Rate constant at T_ref
 * @param Ea    Activation energy [J/mol]
 * @param T_ref Reference temperature [K]
 * @param T     Target temperature [K]
 * @return Rate constant at T
 */
double cdr_arrhenius_from_ref(double k_ref, double Ea,
                               double T_ref, double T);

/* ---------------------------------------------------------------------------
 * L4: Rate Laws
 * ------------------------------------------------------------------------- */

/**
 * @brief Power-law reaction rate (irreversible).
 *
 * r = k * C_A^alpha * C_B^beta
 *
 * @param k      Rate constant [units depend on order]
 * @param C_A    Concentration of A [mol/m^3]
 * @param C_B    Concentration of B [mol/m^3]
 * @param alpha  Order w.r.t. A
 * @param beta   Order w.r.t. B
 * @return Reaction rate [mol/(m^3.s)]
 */
double cdr_power_law_rate(double k, double C_A, double C_B,
                           double alpha, double beta);

/**
 * @brief Reversible power-law reaction rate.
 *
 * r = k_f * C_A^alpha * C_B^beta - k_r * C_C^gamma * C_D^delta
 *
 * @param law  PowerLawRate structure with all parameters
 * @param C_A  Concentration of reactant A [mol/m^3]
 * @param C_B  Concentration of reactant B [mol/m^3]
 * @param C_C  Concentration of product C [mol/m^3]
 * @param C_D  Concentration of product D [mol/m^3]
 * @return Net forward reaction rate [mol/(m^3.s)]
 */
double cdr_reversible_rate(const PowerLawRate *law,
                            double C_A, double C_B,
                            double C_C, double C_D);

/**
 * @brief Michaelis-Menten rate.
 *
 * r = V_max * [S] / (K_m + [S])
 *
 * @param params  Michaelis-Menten parameters
 * @param S       Substrate concentration [mol/m^3]
 * @return Reaction rate [mol/(m^3.s)]
 */
double cdr_michaelis_menten_rate(const MichaelisMentenParams *params,
                                  double S);

/**
 * @brief Michaelis-Menten with competitive inhibition.
 *
 * r = V_max * [S] / (K_m * (1 + [I]/K_i) + [S])
 *
 * @param params  MM parameters
 * @param S       Substrate concentration [mol/m^3]
 * @param I       Inhibitor concentration [mol/m^3]
 * @param K_i     Inhibition constant [mol/m^3]
 * @return Inhibited reaction rate [mol/(m^3.s)]
 */
double cdr_michaelis_menten_inhibited(const MichaelisMentenParams *params,
                                       double S, double I, double K_i);

/**
 * @brief Langmuir-Hinshelwood unimolecular rate.
 *
 * r = k * K_A * P_A / (1 + K_A * P_A + sum(K_j * P_j))
 *
 * @param params   LH parameters
 * @param P_A      Partial pressure of reactant [Pa]
 * @param P_inhib  Array of inhibitor partial pressures [Pa]
 * @param K_inhib  Array of inhibitor adsorption constants [1/Pa]
 * @param n_inhib  Number of inhibitor species
 * @return Surface reaction rate [mol/(m^2.s)]
 */
double cdr_langmuir_hinshelwood_uni(const LangmuirHinshelwoodParams *params,
                                     double P_A,
                                     const double *P_inhib,
                                     const double *K_inhib,
                                     size_t n_inhib);

/**
 * @brief Langmuir-Hinshelwood bimolecular (Langmuir-Hinshelwood) rate.
 *
 * r = k * K_A * K_B * P_A * P_B / (1 + K_A*P_A + K_B*P_B)^2
 *
 * @param params  LH parameters
 * @param P_A     Partial pressure of A [Pa]
 * @param P_B     Partial pressure of B [Pa]
 * @return Surface reaction rate [mol/(m^2.s)]
 */
double cdr_langmuir_hinshelwood_bi(const LangmuirHinshelwoodParams *params,
                                    double P_A, double P_B);

/* ---------------------------------------------------------------------------
 * L5: Kinetic Analysis Methods
 * ------------------------------------------------------------------------- */

/**
 * @brief Half-life for nth-order irreversible reaction (n != 1).
 *
 * t_1/2 = (2^(n-1) - 1) / (k * (n-1) * C_A0^(n-1))
 *
 * @param k    Rate constant [units depend on order]
 * @param n    Reaction order (n != 1)
 * @param C_A0 Initial concentration [mol/m^3]
 * @return Half-life [s], or INFINITY if invalid inputs
 */
double cdr_half_life_nth_order(double k, double n, double C_A0);

/**
 * @brief Half-life for first-order reaction.
 *
 * t_1/2 = ln(2) / k
 *
 * @param k  First-order rate constant [1/s]
 * @return Half-life [s]
 */
double cdr_half_life_first_order(double k);

/**
 * @brief Conversion after time t for first-order batch reaction.
 *
 * X = 1 - exp(-k * t)
 *
 * @param k  First-order rate constant [1/s]
 * @param t  Reaction time [s]
 * @return Conversion X [0-1]
 */
double cdr_conversion_first_order(double k, double t);

/**
 * @brief Time to reach specified conversion (first-order, batch).
 *
 * t = -ln(1 - X) / k
 *
 * @param k  Rate constant [1/s]
 * @param X  Target conversion [0-1)
 * @return Time [s]
 */
double cdr_time_to_conversion_first_order(double k, double X);

/**
 * @brief Concentration evolution for first-order reaction.
 *
 * C_A(t) = C_A0 * exp(-k*t)
 *
 * @param C_A0  Initial concentration [mol/m^3]
 * @param k     Rate constant [1/s]
 * @param t     Time [s]
 * @return Concentration C_A(t) [mol/m^3]
 */
double cdr_concentration_first_order(double C_A0, double k, double t);

/**
 * @brief Concentration evolution for second-order (2A -> products).
 *
 * 1/C_A(t) = 1/C_A0 + k*t
 *
 * @param C_A0  Initial concentration [mol/m^3]
 * @param k     Rate constant [m^3/(mol.s)]
 * @param t     Time [s]
 * @return Concentration C_A(t) [mol/m^3]
 */
double cdr_concentration_second_order(double C_A0, double k, double t);

/**
 * @brief Concentration evolution for zero-order reaction.
 *
 * C_A(t) = C_A0 - k*t  (for t <= C_A0/k)
 *
 * @param C_A0  Initial concentration [mol/m^3]
 * @param k     Rate constant [mol/(m^3.s)]
 * @param t     Time [s]
 * @return Concentration C_A(t) [mol/m^3], minimum 0
 */
double cdr_concentration_zero_order(double C_A0, double k, double t);

/**
 * @brief Concentration profiles for A -> B -> C series reaction.
 *
 * C_A(t) = C_A0 * exp(-k1*t)
 * C_B(t) = C_A0 * k1/(k2-k1) * (exp(-k1*t) - exp(-k2*t))  [k1 != k2]
 * C_C(t) = C_A0 - C_A(t) - C_B(t)
 *
 * @param C_A0   Initial concentration of A [mol/m^3]
 * @param k1     Rate constant A->B [1/s]
 * @param k2     Rate constant B->C [1/s]
 * @param t      Time [s]
 * @param C_out  Output array of size 3: [C_A, C_B, C_C]
 *
 * Complexity: O(1), special handling for k1==k2
 */
void cdr_series_reaction_concentrations(double C_A0, double k1, double k2,
                                         double t, double C_out[3]);

/**
 * @brief Estimate activation energy from two rate measurements.
 *
 * Ea = R * ln(k2/k1) / (1/T1 - 1/T2)
 *
 * @param k1  Rate constant at T1
 * @param T1  Temperature 1 [K]
 * @param k2  Rate constant at T2
 * @param T2  Temperature 2 [K]
 * @return Activation energy Ea [J/mol]
 */
double cdr_activation_energy_from_two_points(double k1, double T1,
                                              double k2, double T2);

/**
 * @brief Compute selectivity for parallel reactions.
 *
 * A -> B (desired, rate constant k1, order n1)
 * A -> C (undesired, rate constant k2, order n2)
 *
 * Instantaneous selectivity: S_B/C = (dC_B/dt) / (dC_C/dt) = (k1/k2)*C_A^(n1-n2)
 *
 * @param k1   Rate constant for desired product [units depend on order]
 * @param k2   Rate constant for undesired product
 * @param n1   Order of desired reaction
 * @param n2   Order of undesired reaction
 * @param C_A  Current concentration of A [mol/m^3]
 * @return Instantaneous selectivity S_B/C [-]
 */
double cdr_selectivity_parallel(double k1, double k2, double n1, double n2,
                                 double C_A);

/**
 * @brief Yield of intermediate B in series reaction A->B->C.
 *
 * Y_B = C_B / C_A0
 *
 * Maximum yield occurs at t_opt = ln(k2/k1) / (k2 - k1)
 * Y_B_max = (k1/k2)^(k2/(k2-k1))
 *
 * @param k1  Rate constant A->B [1/s]
 * @param k2  Rate constant B->C [1/s]
 * @param C_A0 Initial concentration [mol/m^3]
 * @param t    Time [s]
 * @return Yield of B [0-1]
 */
double cdr_yield_intermediate_series(double k1, double k2,
                                      double C_A0, double t);

/**
 * @brief Time for maximum intermediate yield in series reaction.
 *
 * t_opt = ln(k2/k1) / (k2 - k1)
 *
 * @param k1  Rate constant A->B [1/s]
 * @param k2  Rate constant B->C [1/s]
 * @return Optimal time [s] or INFINITY if k1==k2
 */
double cdr_optimal_time_series(double k1, double k2);

#endif /* CDR_REACTION_H */

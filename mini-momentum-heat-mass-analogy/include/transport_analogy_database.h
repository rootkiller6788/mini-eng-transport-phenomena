/**
 * @file transport_analogy_database.h
 * @brief Engineering database of transport property analogies
 *
 * Provides lookup tables and semi-empirical models for practical
 * engineering calculations using the analogy framework.
 *
 * Reference: Bird, Stewart, Lightfoot (2007), Transport Phenomena.
 *            Cengel (2014), Thermodynamics: An Engineering Approach.
 *            Incropera & DeWitt (2007), Fundamentals of Heat and Mass Transfer.
 *            Perry's Chemical Engineers' Handbook (2019), 9th Ed.
 *
 * The analogy is most powerful when combined with engineering databases
 * of transport properties. This module bridges the gap between the
 * theoretical analogy and practical engineering design.
 */

#ifndef TRANSPORT_ANALOGY_DATABASE_H
#define TRANSPORT_ANALOGY_DATABASE_H

#include "transport_coefficients.h"
#include "dimensionless_groups.h"
#include "boundary_layer_analogy.h"
#include "tube_channel_analogy.h"

/* ─── Fluid Property Database Entry ─────────────────────────────── */

/** Engineering fluid entry in database */
typedef struct {
    char name[32];         /**< Fluid name */
    double T_min;          /**< Minimum valid temperature [K] */
    double T_max;          /**< Maximum valid temperature [K] */
    double M;              /**< Molecular weight [kg/mol] */
    double Tc;             /**< Critical temperature [K] */
    double Pc;             /**< Critical pressure [Pa] */
    double rho_ref;        /**< Reference density [kg/m³] */
    double T_ref;          /**< Reference temperature for density [K] */
    SutherlandParams sutherland;  /**< Sutherland parameters (gases) */
    PowerLawViscosityParams plaw; /**< Power-law parameters (liquids) */
    double k_at_ref;       /**< Thermal conductivity at T_ref [W/(m·K)] */
    double cp_at_ref;      /**< Specific heat at T_ref [J/(kg·K)] */
    int is_gas;            /**< 1 if gas, 0 if liquid */
    int fluid_id;          /**< Unique identifier */
} FluidDatabaseEntry;

/** Analogy lookup result for a substance pair at given conditions */
typedef struct {
    FluidDatabaseEntry fluid;      /**< Main fluid properties */
    double T;                      /**< Temperature [K] */
    double P;                      /**< Pressure [Pa] */
    TransportState transport;      /**< Computed transport properties */
    double Re_1ms_1mm;            /**< Re at 1 m/s and 1 mm length scale */
    double Pr;                     /**< Prandtl number */
    double Sc_water_air;          /**< Sc for water vapor in this fluid */
    double St_heat;                /**< Stanton number (if Re computed) */
    char analogy_regime[128];     /**< Description of analogy regime */
    int analogy_applicable;       /**< 1 if Chilton-Colburn holds */
} AnalogyLookup;

/* ─── L3: Engineering Property Databases ────────────────────────── */

/** Number of fluids in the built-in database */
#define N_BUILTIN_FLUIDS 12

/**
 * Initializes the built-in fluid database with 12 common engineering fluids:
 *   0: Air, 1: Water, 2: Engine Oil (SAE 30), 3: Mercury,
 *   4: CO2, 5: Helium, 6: Hydrogen, 7: Steam (H2O gas),
 *   8: Ethylene Glycol, 9: Glycerin, 10: Nitrogen, 11: Ammonia
 *
 * Complexity: O(1)
 */
void init_fluid_database(FluidDatabaseEntry db[N_BUILTIN_FLUIDS]);

/**
 * Looks up a fluid by ID and computes transport properties at given (T,P).
 * If gas, uses Sutherland's law for μ, Eucken for k.
 * If liquid, uses power-law model for μ, linear interpolation for k.
 *
 * @param db      Fluid database (initialized with init_fluid_database)
 * @param id      Fluid ID (0-11)
 * @param T       Temperature [K]
 * @param P       Pressure [Pa]
 * @param lookup  Output: computed analogy lookup
 * @return 0 on success, -1 if fluid not found, -2 if T out of range
 *
 * Complexity: O(1)
 */
int fluid_analogy_lookup(const FluidDatabaseEntry db[N_BUILTIN_FLUIDS],
                          int id, double T, double P,
                          AnalogyLookup *lookup);

/* ─── L3: Common Engineering Benchmarks ─────────────────────────── */

/**
 * Computes the Reynolds analogy benchmark for air flowing over a flat plate.
 * Air at 300K, 1 atm, U∞ = 10 m/s, plate length L = 1m.
 *
 * Expected values:
 *   Re_L ≈ 6.9×10⁵
 *   Cf_L ≈ 0.00265 (turbulent)
 *   St ≈ 0.00132 (Reynolds analogy)
 *   h ≈ 15 W/(m²·K)
 *
 * Complexity: O(1)
 */
void air_flat_plate_benchmark(double U_inf, double L, double T_inf,
                               double T_wall, double P_atm,
                               BoundaryLayerSummary *result);

/**
 * Computes the water pipe flow analogy benchmark.
 * Water at 300K, D = 0.0254 m, v = 1 m/s, L = 10 m.
 *
 * Expected values:
 *   Re_D ≈ 3.0×10⁴
 *   f ≈ 0.023 (Blasius)
 *   Nu_D ≈ 176 (Dittus-Boelter, heating)
 *   h ≈ 4240 W/(m²·K)
 *
 * Complexity: O(1)
 */
void water_pipe_flow_benchmark(double D, double v, double L,
                                double T_bulk, double T_wall,
                                PipeFlowAnalogy *result);

/**
 * Compares the analogy prediction for two different fluids
 * (e.g., air vs. water) to demonstrate the universality of the
 * dimensionless analogy framework.
 *
 * Complexity: O(1)
 */
void compare_fluid_analogies(int fluid_id_A, int fluid_id_B,
                             double T, double P, double v, double L);

/* ─── L7: Engineering Application Utilities ─────────────────────── */

/**
 * Heat exchanger sizing using the analogy method.
 * Given flow conditions and required heat duty, estimates the required
 * heat transfer area using Chilton-Colburn analogy.
 *
 * This is the practical bridge: friction data → h → heat exchanger size.
 *
 * @param Q_dot        Required heat duty [W]
 * @param delta_T_lm   Log-mean temperature difference [K]
 * @param rho, v, cp   Fluid properties
 * @param k            Thermal conductivity [W/(m·K)]
 * @param mu           Dynamic viscosity [Pa·s]
 * @param D_h          Hydraulic diameter [m]
 * @param Pr, Re       Dimensionless numbers
 * @param area         Output: required heat transfer area [m²]
 *
 * Complexity: O(1)
 */
void heat_exchanger_sizing_by_analogy(double Q_dot, double delta_T_lm,
                                       double rho, double v, double cp,
                                       double k, double mu, double D_h,
                                       double Pr, double Re,
                                       double *area);

/**
 * Cooling tower analysis using mass transfer analogy.
 * Uses the analogy between heat and mass transfer to estimate
 * evaporation rate from known heat transfer characteristics.
 *
 * Merkel's method: the mass transfer coefficient is predicted from
 * the heat transfer coefficient using the Lewis relation.
 *
 * @param h           Heat transfer coefficient [W/(m²·K)]
 * @param rho, cp     Fluid properties
 * @param Le          Lewis number
 * @param A           Contact area [m²]
 * @param delta_omega Humidity ratio difference [kg_v/kg_da]
 * @return Evaporation rate [kg/s]
 *
 * Complexity: O(1)
 */
double cooling_tower_evaporation_rate(double h, double rho, double cp,
                                       double Le, double A,
                                       double delta_omega);

/**
 * Chemical reactor wall mass transfer using analogy.
 * For tubular reactors, predicts mass transfer coefficient at wall
 * (e.g., for catalytic wall reactions) from pressure drop data.
 *
 * @param delta_P     Pressure drop [Pa]
 * @param D_tube      Tube diameter [m]
 * @param L_tube      Tube length [m]
 * @param rho, v      Fluid density and velocity
 * @param D_AB        Binary diffusivity [m²/s]
 * @param Sc          Schmidt number
 * @return Predicted mass transfer coefficient [m/s]
 *
 * Complexity: O(1)
 */
double reactor_wall_mass_transfer(double delta_P, double D_tube,
                                  double L_tube, double rho, double v,
                                  double D_AB, double Sc);

/**
 * Electronics cooling analysis using the analogy.
 * For forced convection over PCBs/chips, predicts heat transfer
 * coefficient from known flow characteristics.
 *
 * Typical application: air cooling of a CPU heat sink.
 *
 * @param v_air       Air velocity [m/s]
 * @param L_chip      Chip characteristic length [m]
 * @param T_chip      Chip temperature [K]
 * @param T_ambient   Ambient temperature [K]
 * @param P_atm       Atmospheric pressure [Pa]
 * @param h           Output: heat transfer coefficient [W/(m²·K)]
 * @return Heat flux [W/m²]
 *
 * Complexity: O(1)
 */
double electronics_cooling_analogy(double v_air, double L_chip,
                                   double T_chip, double T_ambient,
                                   double P_atm, double *h);

/**
 * Spreadsheet-friendly lookup: returns a complete analogy summary
 * for quick engineering estimation.
 *
 * Example: given "air at 25°C, 1 m/s over 10 cm plate",
 * returns all relevant dimensionless numbers and predicted
 * heat/mass transfer coefficients.
 */
void quick_analogy_estimate(const char *fluid_name,
                            double T_C, double P_kPa,
                            double v_ms, double L_mm,
                            char *summary, size_t summary_size);

/**
 * Prints comparison table of analogy predictions vs. experimental
 * correlations for common engineering configurations.
 */
void print_analogy_vs_experiment_comparison(void);

/**
 * Multi-fluid analogy summary: prints a table of transport properties
 * and analogy regimes for all built-in fluids at given conditions.
 */
void print_multi_fluid_analogy_table(double T, double P);

#endif /* TRANSPORT_ANALOGY_DATABASE_H */

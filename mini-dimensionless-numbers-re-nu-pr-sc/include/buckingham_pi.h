/*
 * buckingham_pi.h — Buckingham Pi Theorem and Dimensional Analysis
 *
 * Reference: Buckingham (1914) "On Physically Similar Systems"
 *            Barenblatt (2003) "Scaling"
 *            White (2016) "Fluid Mechanics", Ch. 5
 *            Bird, Stewart, Lightfoot (2007) "Transport Phenomena", Ch. 3
 *
 * The Buckingham Pi theorem is the mathematical foundation of dimensional
 * analysis and dimensionless numbers. It states that any physically meaningful
 * equation relating n variables with k fundamental dimensions can be rewritten
 * in terms of n-k dimensionless Pi groups.
 *
 * Knowledge Coverage: L1 Definition, L2 Core Concept (Pi Theorem),
 *                     L3 Engineering Method, L4 Dimensional Homogeneity,
 *                     L5 Algorithm (Matrix Method for Pi Groups),
 *                     L7 Application (Model Scaling)
 */

#ifndef BUCKINGHAM_PI_H
#define BUCKINGHAM_PI_H

#include <stddef.h>  /* size_t */

#ifdef __cplusplus
extern "C" {
#endif

/* ==========================================================================
 * L1: Fundamental Dimensions and Unit Systems
 * ========================================================================== */

/**
 * FundamentalDimension — the seven SI base dimensions.
 */
typedef enum {
    DIM_MASS         = 0,  /* M — kilogram [kg] */
    DIM_LENGTH       = 1,  /* L — metre [m] */
    DIM_TIME         = 2,  /* T — second [s] */
    DIM_TEMPERATURE  = 3,  /* Θ — kelvin [K] */
    DIM_CURRENT      = 4,  /* I — ampere [A] */
    DIM_AMOUNT       = 5,  /* N — mole [mol] */
    DIM_LUMINOUS     = 6,  /* J — candela [cd] */
    DIM_COUNT        = 7
} FundamentalDimension;

/**
 * DerivedDimension — exponents for dimensional representation.
 *
 * Any physical quantity Q has a dimension [Q] = M^a.L^b.T^c.Θ^d.I^e.N^f.J^g
 * represented as a vector of exponents.
 */
typedef struct {
    double exponents[7];  /* M, L, T, Θ, I, N, J */
} DimensionVector;

/**
 * PhysicalQuantity — a named variable with its dimensional representation.
 *
 * Used in Pi theorem analysis: collect all relevant variables,
 * compute their dimensional matrix, and derive dimensionless groups.
 */
typedef struct {
    char name[32];            /* variable name (e.g., "velocity", "density") */
    DimensionVector dim;      /* dimensional exponents */
    double value;             /* representative value (for scaling) */
} PhysicalQuantity;

/**
 * DimensionalMatrix — the complete dimensional matrix for Pi analysis.
 *
 * n_vars rows (physical quantities), n_dims columns (fundamental dimensions).
 * Each entry D[i][j] is the exponent of dimension j in variable i.
 */
typedef struct {
    int n_vars;               /* number of physical variables */
    int n_dims;               /* number of fundamental dimensions involved */
    double **matrix;          /* [n_vars][n_dims] row-major matrix */
    PhysicalQuantity *vars;   /* array of variable metadata */
} DimensionalMatrix;

/**
 * PiGroup — a dimensionless product Pi.
 *
 * Pi = ∏_{i=1}^{n} Q_i^{e_i}  such that [Pi] = 1 (dimensionless).
 */
typedef struct {
    double *exponents;        /* exponents for each variable, length = n_vars */
    char expression[256];     /* human-readable Pi group expression */
    double value;             /* computed numerical value of this Pi group */
} PiGroup;

/* ==========================================================================
 * L2: Pi Theorem Implementation
 * ========================================================================== */

/**
 * dimension_vector_create — build a dimension vector from exponents.
 *
 * @param exponents  array of 7 doubles for M,L,T,Θ,I,N,J
 * @param[out] dv    filled dimension vector
 * @return           0 on success
 */
int dimension_vector_create(const double exponents[7], DimensionVector *dv);

/**
 * dimension_vector_is_dimensionless — check if [Q] = 1.
 *
 * @param dv  dimension vector
 * @return    1 if all exponents are zero (dimensionless), 0 otherwise
 */
int dimension_vector_is_dimensionless(const DimensionVector *dv);

/**
 * dimension_vector_multiply — multiply two dimension vectors.
 *
 * [Q1.Q2] = [Q1] + [Q2] (exponent-wise addition)
 *
 * @param a     first vector
 * @param b     second vector
 * @param[out] result = a + b
 */
void dimension_vector_multiply(const DimensionVector *a,
                               const DimensionVector *b,
                               DimensionVector *result);

/**
 * dimension_vector_power — raise dimension vector to a power.
 *
 * [Q^n] = n.[Q] (exponent-wise scalar multiplication)
 *
 * @param dv      dimension vector
 * @param power   exponent
 * @param[out] result
 */
void dimension_vector_power(const DimensionVector *dv, double power,
                            DimensionVector *result);

/**
 * dimension_vector_divide — [Q1/Q2] = [Q1] - [Q2]
 */
void dimension_vector_divide(const DimensionVector *a,
                             const DimensionVector *b,
                             DimensionVector *result);

/* ==========================================================================
 * L2-L5: Standard Dimension Vectors for Transport Phenomena
 * ========================================================================== */

/**
 * get_standard_dimension — return dimension vector for common variables.
 *
 * Supported variables:
 *   "density"         — M.L⁻^3
 *   "velocity"        — L.T⁻¹
 *   "length"          — L
 *   "time"            — T
 *   "viscosity"       — M.L⁻¹.T⁻¹
 *   "kinematic_viscosity" — L^2.T⁻¹
 *   "pressure"        — M.L⁻¹.T⁻^2
 *   "force"           — M.L.T⁻^2
 *   "energy"          — M.L^2.T⁻^2
 *   "power"           — M.L^2.T⁻^3
 *   "thermal_conductivity" — M.L.T⁻^3.Θ⁻¹
 *   "heat_capacity"   — L^2.T⁻^2.Θ⁻¹
 *   "heat_transfer_coeff" — M.T⁻^3.Θ⁻¹
 *   "mass_diffusivity" — L^2.T⁻¹
 *   "temperature"     — Θ
 *   "heat_flux"       — M.T⁻^3
 *   "mass_flux"       — M.L⁻^2.T⁻¹
 *   "surface_tension" — M.T⁻^2
 *   "gravity"         — L.T⁻^2
 *   "angular_velocity" — T⁻¹
 *   "frequency"       — T⁻¹
 *   "area"            — L^2
 *   "volume"          — L^3
 *   "volumetric_flow_rate" — L^3.T⁻¹
 *   "mass_flow_rate"  — M.T⁻¹
 *
 * @param name   variable name string
 * @param[out] dv filled dimension vector (zeroed on unrecognized name)
 * @return       0 on success, -1 if name not recognized
 */
int get_standard_dimension(const char *name, DimensionVector *dv);

/**
 * get_dimension_name — human-readable dimension string.
 *
 * @param dv  dimension vector
 * @param buf output buffer
 * @param sz  buffer size
 * @return    0 on success
 */
int get_dimension_name(const DimensionVector *dv, char *buf, size_t sz);

/**
 * get_derived_dimension — compute dimension of a product expression.
 *
 * Given an expression like v^a . L^b . rho^c, compute its overall dimension.
 *
 * @param varnames  array of variable name strings
 * @param exponents array of exponents (same length)
 * @param n_terms   number of terms
 * @param[out] dv   resulting dimension vector
 * @return          0 on success, -1 on unrecognized variable
 *
 * Example: varnames = {"velocity","length","kinematic_viscosity"}
 *          exponents = {1, 1, -1}
 *          -> computes [v.L/nu] = [Re] = dimensionless
 */
int get_derived_dimension(const char **varnames, const double *exponents,
                          int n_terms, DimensionVector *dv);

/* ==========================================================================
 * L3: Dimensionless Group Verification
 * ========================================================================== */

/**
 * verify_dimensionless_group — check if a product is truly dimensionless.
 *
 * Given: Pi = v1^e1 . v2^e2 . ... . vn^en,
 * verify that [Pi] = 1.
 *
 * @param varnames   array of variable names
 * @param exponents  array of exponents
 * @param n_vars     number of variables
 * @param tol        tolerance for floating-point zero check
 * @return           1 if dimensionless, 0 otherwise
 *
 * Example: For Re = rho.v.D/mu:
 *   varnames    = {"density","velocity","length","viscosity"}
 *   exponents   = {1, 1, 1, -1}
 *   -> returns 1 (dimensionless)
 */
int verify_dimensionless_group(const char **varnames,
                               const double *exponents,
                               int n_vars, double tol);

/**
 * verify_all_standard_groups — verify that all standard dimensionless
 *   groups are indeed dimensionless using dimensional analysis.
 *
 * @param[out] results  array of 1/0 results, must be at least 16 elements
 *                      Re, Nu, Pr, Sc, Pe_H, Pe_M, Gr, Ra, St_H, St_M,
 *                      We, Ma, Fr, Eu, Br, Bi, Fo
 * @return              number of groups verified
 *
 * This function provides a self-consistency check: each standard
 * dimensionless group is proven to be dimensionless via the exponent method.
 */
int verify_all_standard_groups(int *results);

/* ==========================================================================
 * L5: Matrix Method for Pi Theorem
 * ========================================================================== */

/**
 * pi_theorem_matrix_method — find dimensionless groups using the matrix method.
 *
 * Given n physical variables involving k fundamental dimensions, this solves
 * A.e = 0 where A is the (k×n) dimensional matrix and e are exponents.
 * The nullspace of A gives n-k independent Pi groups.
 *
 * For transport phenomena, the three key cases are:
 *   1. Fluid mechanics:    vars = {F, rho, U, L, mu, g, c}
 *      k = 3 (M, L, T) -> 7-3 = 4 Pi groups: Re, Ma, Fr, Eu
 *
 *   2. Forced convection:  vars = {h, rho, U, L, mu, k, cp}
 *      k = 4 (M, L, T, Θ) -> 7-4 = 3 Pi groups: Nu, Re, Pr
 *
 *   3. Natural convection: vars = {h, rho, L, mu, k, cp, g, beta, DeltaT}
 *      k = 4 (M, L, T, Θ) -> 9-4 = 5 Pi groups: Nu, Gr, Pr, ...
 *
 * @param vars       array of physical quantities involved
 * @param n_vars     number of variables
 * @param k_dims     number of independent fundamental dimensions
 * @param[out] pi_groups  array of (n_vars - k_dims) PiGroup structs
 * @param max_groups      max capacity of pi_groups array
 * @param[out] n_groups   actual number found
 * @return          0 on success, -1 on dimension mismatch, -2 on rank deficiency
 */
int pi_theorem_matrix_method(const PhysicalQuantity *vars, int n_vars,
                             int k_dims, PiGroup *pi_groups,
                             int max_groups, int *n_groups);

/**
 * pi_rank_dimensional_matrix — compute the rank of the dimensional matrix.
 *
 * Uses Gaussian elimination with partial pivoting.
 *
 * @param rows  number of variables
 * @param cols  number of dimensions
 * @param D     row-major matrix [rows][cols]
 * @param tol   tolerance for zero detection
 * @return      rank of the matrix, -1 on invalid inputs
 */
int pi_rank_dimensional_matrix(int rows, int cols, double **D, double tol);

/**
 * pi_nullspace — compute nullspace vectors of the dimensional matrix.
 *
 * @param rows      number of variables = n
 * @param cols      number of dimensions = k
 * @param D         dimensional matrix [rows][cols]
 * @param[out] nullspace [rows][rows-cols] nullspace vectors (column vectors)
 * @param tol       tolerance
 * @return          number of nullspace vectors (= rows - rank), -1 on error
 *
 * Each nullspace vector corresponds to one Pi group exponent set.
 */
int pi_nullspace(int rows, int cols, double **D,
                 double **nullspace, double tol);

/**
 * pi_groups_to_string — generate human-readable Pi group expressions.
 *
 * @param pi_groups  array of PiGroup
 * @param n_groups   number of groups
 * @param vars       array of variable names
 * @param n_vars     number of variables
 * @return           0 on success
 */
int pi_groups_to_string(PiGroup *pi_groups, int n_groups,
                        const PhysicalQuantity *vars, int n_vars);

/* ==========================================================================
 * L7: Model Scaling — Wind Tunnel / Water Channel Design
 * ========================================================================== */

/**
 * model_scale_velocity_reynolds — compute model test velocity for Re matching.
 *
 * Prototype: (L_p, U_p, nu_p); Model: (L_m, nu_m)
 * Target: Re_p = Re_m -> U_m = U_p.(L_p/L_m).(nu_m/nu_p)
 *
 * @param L_p  prototype length [m]
 * @param U_p  prototype velocity [m/s]
 * @param nu_p prototype kinematic viscosity [m^2/s]
 * @param L_m  model length [m]
 * @param nu_m model kinematic viscosity [m^2/s]
 * @return     required model velocity [m/s]
 */
double model_scale_velocity_reynolds(double L_p, double U_p, double nu_p,
                                     double L_m, double nu_m);

/**
 * model_scale_force — compute force scale factor from Re matching.
 *
 * F_p/F_m = (rho_p/rho_m).(L_p/L_m)^2.(U_p/U_m)^2
 * When Re_p = Re_m, U_p/U_m = (nu_p/nu_m).(L_m/L_p), so:
 * F_p/F_m = (rho_p/rho_m).(nu_p/nu_m)^2
 *
 * @param rho_p  prototype density [kg/m^3]
 * @param nu_p   prototype viscosity [m^2/s]
 * @param rho_m  model density [kg/m^3]
 * @param nu_m   model viscosity [m^2/s]
 * @return       force ratio F_prototype / F_model [-]
 */
double model_scale_force(double rho_p, double nu_p, double rho_m, double nu_m);

/**
 * model_scale_power — power ratio for dynamical similarity.
 *
 * P ∝ F.U -> P_p/P_m = (F_p/F_m).(U_p/U_m)
 *
 * @param rho_p, nu_p, L_p, U_p  prototype properties
 * @param rho_m, nu_m, L_m       model properties
 * @return                       P_p/P_m [-]
 */
double model_scale_power(double rho_p, double nu_p, double L_p, double U_p,
                         double rho_m, double nu_m, double L_m);

/**
 * model_scale_frequency — frequency/timescale for dynamic similarity.
 *
 * Strouhal number matching: St = f.L/U = const
 * f_m = f_p.(L_p/L_m).(U_m/U_p)
 *
 * @param f_p  prototype frequency [Hz]
 * @param L_p  prototype length [m]
 * @param U_p  prototype velocity [m/s]
 * @param L_m  model length [m]
 * @param U_m  model velocity [m/s]
 * @return     model frequency [Hz]
 */
double model_scale_frequency(double f_p, double L_p, double U_p,
                             double L_m, double U_m);

/* ==========================================================================
 * L5: Scaling Laws from Dimensional Analysis
 * ========================================================================== */

/**
 * scaling_law_pipe_pressure_drop — derive DeltaP ∝ f(L/D).(rhoU^2/2) from Pi theorem.
 *
 * Given variables {DeltaP, rho, U, D, L, mu, epsilon}, k=3 (M,L,T) -> 4 Pi groups.
 * The result is the Darcy-Weisbach equation: Eu = f.(L/D).
 *
 * @param f    friction factor [-]
 * @param rho  density [kg/m^3]
 * @param U    velocity [m/s]
 * @param L    pipe length [m]
 * @param D    pipe diameter [m]
 * @return     DeltaP [Pa]
 */
double scaling_law_pipe_pressure_drop(double f, double rho, double U,
                                      double L, double D);

/**
 * scaling_law_buoyant_plume — velocity and temperature in a buoyant plume.
 *
 * Dimensional analysis of {u, z, g, beta, Q, nu, alpha} gives:
 *   u ~ (g.beta.Q/(rho.cp))^(1/3).z^(-1/3)
 *
 * @param g     gravity [m/s^2]
 * @param beta  thermal expansion coefficient [1/K]
 * @param Q_dot heat input [W]
 * @param rho   density [kg/m^3]
 * @param cp    specific heat [J/(kg.K)]
 * @param z     height above source [m]
 * @return      plume centerline velocity [m/s]
 *
 * Reference: Turner (1979) "Buoyancy Effects in Fluids"
 */
double scaling_law_buoyant_plume_velocity(double g, double beta, double Q_dot,
                                          double rho, double cp, double z);

/**
 * scaling_law_kolmogorov_microscale — smallest turbulent length scale.
 *
 * eta = (nu^3/epsilon)^(1/4)  where epsilon is dissipation rate.
 *
 * Using large-scale parameters: eta/L ~ Re^(-3/4)
 *
 * @param L    integral length scale [m]
 * @param Re   Reynolds number based on L
 * @return     Kolmogorov length scale eta [m]
 */
double scaling_law_kolmogorov_microscale(double L, double Re);

/**
 * scaling_law_batchelor_microscale — thermal microscale in turbulence.
 *
 * eta_T = eta.Pr^(-1/2) for Pr ≥ 1
 * eta_T = eta.Pr^(-3/4) for Pr << 1
 *
 * @param eta  Kolmogorov microscale [m]
 * @param Pr   Prandtl number [-]
 * @return     Batchelor microscale eta_T [m]
 */
double scaling_law_batchelor_microscale(double eta, double Pr);

/**
 * scaling_law_diffusion_length — characteristic diffusion distance.
 *
 * L_diff ~ sqrt(alpha.t)  (heat) or sqrt(D.t)  (mass)
 *
 * @param diffusivity  thermal or mass diffusivity [m^2/s]
 * @param t            time [s]
 * @return             diffusion length [m]
 */
double scaling_law_diffusion_length(double diffusivity, double t);

#ifdef __cplusplus
}
#endif

#endif /* BUCKINGHAM_PI_H */

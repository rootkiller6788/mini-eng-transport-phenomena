/**
 * @file boundary_layer_analogy.h
 * @brief Analogy between velocity, thermal, and concentration boundary layers
 *
 * The boundary layer concept unifies momentum, heat, and mass transfer near
 * solid surfaces. This module provides the mathematical framework for
 * comparing and relating the three boundary layer thicknesses.
 *
 * Reference: Bird, Stewart, Lightfoot (2007), Transport Phenomena,
 *            Chapters 4, 12, 19, 22.
 *            Schlichting & Gersten (2000), Boundary Layer Theory.
 *
 * Key relationships:
 *   Velocity BL thickness:    δ ∝ x / sqrt(Re_x)
 *   Thermal BL thickness:     δ_T ≈ δ · Pr^(-1/3)  (laminar)
 *   Concentration BL thickness: δ_C ≈ δ · Sc^(-1/3) (laminar)
 *
 * The relative thicknesses depend on Pr and Sc:
 *   Pr > 1: δ_T < δ  (thermal BL thinner, e.g., oils)
 *   Pr < 1: δ_T > δ  (thermal BL thicker, e.g., liquid metals)
 *   Pr ≈ 1: δ_T ≈ δ  (gases)
 */

#ifndef BOUNDARY_LAYER_ANALOGY_H
#define BOUNDARY_LAYER_ANALOGY_H

#include "transport_coefficients.h"

/* ─── L1: Boundary Layer Definitions ───────────────────────────── */

/** Blasius boundary layer profile for laminar flat-plate flow */
typedef struct {
    double eta;          /**< Similarity variable η = y·sqrt(U∞/(ν·x)) */
    double f;            /**< Stream function f(η), where u/U∞ = f'(η) */
    double fp;           /**< f'(η) = u/U∞, dimensionless velocity */
    double fpp;          /**< f''(η), dimensionless shear (τ_w related) */
} BlasiusProfilePoint;

/** Boundary layer thickness comparison for all three types */
typedef struct {
    double delta;        /**< Velocity BL thickness (99% of U∞) [m] */
    double delta_T;      /**< Thermal BL thickness (99% of ΔT) [m] */
    double delta_C;      /**< Concentration BL thickness (99% of ΔC) [m] */
    double delta_d;      /**< Displacement thickness [m] */
    double delta_m;      /**< Momentum thickness [m] */
    double delta_h;      /**< Enthalpy thickness (thermal) [m] */
    double delta_c;      /**< Concentration thickness [m] */
    double x;            /**< Distance from leading edge [m] */
    double Re_x;         /**< Local Reynolds number at x */
    double Pr;           /**< Prandtl number */
    double Sc;           /**< Schmidt number */
} BoundaryLayerThicknesses;

/** Flat-plate boundary layer solution state */
typedef struct {
    double x;            /**< Distance from leading edge [m] */
    double U_inf;        /**< Free-stream velocity [m/s] */
    double T_inf;        /**< Free-stream temperature [K] */
    double T_w;          /**< Wall temperature [K] */
    double C_inf;        /**< Free-stream concentration [mol/m³] */
    double C_w;          /**< Wall concentration [mol/m³] */
    double rho;          /**< Density [kg/m³] */
    double mu;           /**< Dynamic viscosity [Pa·s] */
    double nu;           /**< Kinematic viscosity [m²/s] */
    double alpha;        /**< Thermal diffusivity [m²/s] */
    double D_AB;         /**< Mass diffusivity [m²/s] */
    double Re_x;         /**< Local Reynolds number */
    double Pr;           /**< Prandtl number */
    double Sc;           /**< Schmidt number */
    int flow_regime;     /**< 0=laminar, 1=transitional, 2=turbulent */
} FlatPlateState;

/** Pohlhausen solution for thermal boundary layer (laminar, constant T_w) */
typedef struct {
    double eta;          /**< Similarity variable */
    double Theta;        /**< Dimensionless temperature (T-T_w)/(T_inf-T_w) */
    double Theta_p;      /**< dΘ/dη */
    double Pr;           /**< Prandtl number */
} PohlhausenProfilePoint;

/** Analogy summary: all three boundary layers at one streamwise station */
typedef struct {
    double x;                    /**< Streamwise position [m] */
    FlatPlateState flow;         /**< Flow state */
    BoundaryLayerThicknesses bl; /**< Thicknesses */
    double Cf_x;                 /**< Local skin friction */
    double Nu_x;                 /**< Local Nusselt number */
    double Sh_x;                 /**< Local Sherwood number */
    double j_H_x;                /**< Local Colburn j_H */
    double j_D_x;                /**< Local Colburn j_D */
} BoundaryLayerSummary;

/* ─── L2: Velocity Boundary Layer Functions ────────────────────── */

/**
 * Blasius boundary layer similarity solution for flat plate.
 * Solves f''' + (1/2)·f·f'' = 0 with BC: f(0)=f'(0)=0, f'(∞)=1.
 *
 * Uses the Runge-Kutta shooting method.
 * Complexity: O(n_steps) where n_steps = number of η integration steps.
 *
 * @param eta    Array of η values to compute profile at
 * @param n      Number of points
 * @param profile Output profile array (must be pre-allocated to size n)
 */
void blasius_profile(const double *eta, size_t n, BlasiusProfilePoint *profile);

/**
 * Computes the Blasius similarity variable η at given physical location.
 * η = y · sqrt(U_inf / (ν · x))
 * Complexity: O(1)
 */
double blasius_eta(double y, double U_inf, double nu, double x);

/**
 * Computes laminar boundary layer thickness at distance x.
 * δ_99% / x = 5.0 / sqrt(Re_x)
 * Complexity: O(1)
 */
double laminar_bl_thickness(double x, double Re_x);

/**
 * Computes turbulent boundary layer thickness at distance x.
 * δ / x = 0.37 / Re_x^(1/5)   (1/7th power law, smooth flat plate)
 * Complexity: O(1)
 */
double turbulent_bl_thickness(double x, double Re_x);

/**
 * Computes local skin friction coefficient for laminar flow.
 * Cf_x = 0.664 / sqrt(Re_x)   (Blasius)
 * Complexity: O(1)
 */
double laminar_cf_local(double Re_x);

/**
 * Computes local skin friction coefficient for turbulent flow.
 * Cf_x = 0.0592 / Re_x^(1/5)   (1/7th power law, Re_x < 10⁷)
 * Complexity: O(1)
 */
double turbulent_cf_local(double Re_x);

/**
 * Computes momentum thickness from Blasius solution.
 * θ / x = 0.664 / sqrt(Re_x)
 * Complexity: O(1)
 */
double laminar_momentum_thickness(double x, double Re_x);

/**
 * Computes displacement thickness from Blasius solution.
 * δ* / x = 1.7208 / sqrt(Re_x)
 * Complexity: O(1)
 */
double laminar_displacement_thickness(double x, double Re_x);

/* ─── L2: Thermal Boundary Layer Functions ─────────────────────── */

/**
 * Pohlhausen thermal boundary layer solution for laminar flow.
 * Solves Θ'' + (Pr/2)·f·Θ' = 0 where f is the Blasius stream function.
 * BC: Θ(0)=0, Θ(∞)=1 (isothermal wall).
 *
 * For Pr = 1: Θ(η) = 1 - f'(η)  (Reynolds analogy exact)
 *
 * Complexity: O(n_steps)
 */
void pohlhausen_thermal_profile(const double *eta, const BlasiusProfilePoint *blasius,
                                double Pr, size_t n,
                                PohlhausenProfilePoint *thermal_profile);

/**
 * Thermal boundary layer thickness (99% criterion).
 * δ_T / δ ≈ Pr^(-1/3)  (for Pr ≥ 0.6, laminar)
 *
 * More precisely: δ_T/δ = 0.975·Pr^(-1/3)
 *
 * Complexity: O(1)
 */
double thermal_bl_thickness_laminar(double delta_velocity, double Pr);

/**
 * Local Nusselt number from Pohlhausen solution for laminar flat plate.
 * Nu_x = 0.332 · Re_x^(1/2) · Pr^(1/3)   (0.6 < Pr < 50)
 *
 * This is the fundamental result that connects h to the analogy.
 * Complexity: O(1)
 */
double laminar_nusselt_local(double Re_x, double Pr);

/**
 * Average Nusselt number over plate length L (laminar).
 * Nu_L = 0.664 · Re_L^(1/2) · Pr^(1/3)
 *
 * Complexity: O(1)
 */
double laminar_nusselt_average(double Re_L, double Pr);

/**
 * Local Nusselt number for turbulent flat plate flow.
 * Nu_x = 0.0296 · Re_x^(4/5) · Pr^(1/3)   (0.6 < Pr < 60)
 *
 * Complexity: O(1)
 */
double turbulent_nusselt_local(double Re_x, double Pr);

/* ─── L2: Concentration Boundary Layer Functions ───────────────── */

/**
 * Local Sherwood number for laminar flat plate.
 * Sh_x = 0.332 · Re_x^(1/2) · Sc^(1/3)   (0.6 < Sc < 3000)
 *
 * Exact analog of the thermal Nusselt number correlation.
 * Complexity: O(1)
 */
double laminar_sherwood_local(double Re_x, double Sc);

/**
 * Average Sherwood number over plate length L (laminar).
 * Sh_L = 0.664 · Re_L^(1/2) · Sc^(1/3)
 *
 * Complexity: O(1)
 */
double laminar_sherwood_average(double Re_L, double Sc);

/**
 * Local Sherwood number for turbulent flat plate flow.
 * Sh_x = 0.0296 · Re_x^(4/5) · Sc^(1/3)
 *
 * Complexity: O(1)
 */
double turbulent_sherwood_local(double Re_x, double Sc);

/**
 * Concentration boundary layer thickness.
 * δ_C / δ ≈ Sc^(-1/3)  (for Sc ≥ 0.6, laminar)
 *
 * Complexity: O(1)
 */
double concentration_bl_thickness(double delta_velocity, double Sc);

/* ─── L4: The Analogy Identity ─────────────────────────────────── */

/**
 * Demonstrates the fundamental boundary layer analogy:
 * When Pr = 1, the dimensionless temperature profile equals the
 * dimensionless velocity profile.
 *
 * Θ(η) = u(η)/U∞  when Pr = 1 and T_w = constant.
 *
 * This function computes the correlation between the two profiles.
 * Returns Pearson correlation coefficient r ∈ [0,1].
 *
 * Complexity: O(n)
 */
double profile_correlation(const BlasiusProfilePoint *velocity,
                           const PohlhausenProfilePoint *thermal,
                           size_t n);

/**
 * Computes the complete boundary layer summary at one x-location.
 * Fills all fields of BoundaryLayerSummary.
 *
 * Complexity: O(1) (uses algebraic correlations)
 */
void compute_boundary_layer_summary(const FlatPlateState *state,
                                    BoundaryLayerSummary *summary);

/**
 * Identifies which boundary layer is thicker and provides physical explanation.
 * Returns: -1 if δ_C > δ_T > δ (Sc < Pr < 1, e.g., liquid metals)
 *           0 if all similar (Pr ≈ Sc ≈ 1, e.g., gases)
 *          +1 if δ > δ_T > δ_C (1 < Pr < Sc, e.g., oils diffusing in air)
 */
int boundary_layer_ordering(const BoundaryLayerThicknesses *bl);

/**
 * Prints the boundary layer comparison table showing all relevant
 * thicknesses, profiles, and their physical interpretation.
 */
void print_boundary_layer_comparison(const BoundaryLayerSummary *summary);

#endif /* BOUNDARY_LAYER_ANALOGY_H */

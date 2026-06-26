/**
 * @file cdr_core.h
 * @brief Convection-Diffusion-Reaction (CDR) Core Definitions
 *
 * Defines the fundamental dimensionless numbers, state structures, and
 * material property containers for convective-diffusive-reactive transport.
 *
 * Reference: Bird, Stewart, Lightfoot "Transport Phenomena" (2007), Ch. 17-24
 * Reference: Levenspiel "Chemical Reaction Engineering" (1999), Ch. 1-6
 *
 * L1 Definitions: Peclet, Damkohler, Thiele, Reynolds, Schmidt, Sherwood
 * L2 Core Concepts: Transport regimes, rate-limiting steps
 * L3 Engineering Quantities: Characteristic scales, dimensionless groups
 */

#ifndef CDR_CORE_H
#define CDR_CORE_H

#include <stddef.h>

/* ---------------------------------------------------------------------------
 * L1: Fundamental Dimensionless Numbers
 * ------------------------------------------------------------------------- */

/**
 * @brief Peclet number -- ratio of convective to diffusive transport rate.
 *
 * Pe = u * L / D
 *
 * Pe << 1 : Diffusion-dominated (microfluidic, porous media)
 * Pe >> 1 : Convection-dominated (turbulent flow, large-scale reactors)
 * Pe ~ 1  : Mixed regime
 *
 * Two variants: Pe_m (mass transfer, with D_AB), Pe_h (heat transfer, with alpha)
 */
typedef struct {
    double mass_peclet;    /**< Pe_m = u*L / D_AB, species mass transport */
    double heat_peclet;    /**< Pe_h = u*L / alpha, thermal transport */
    double dispersion_pe;  /**< Pe_d = u*d_p / D_ax, axial dispersion in packed beds */
} PecletNumbers;

/**
 * @brief Damkohler number -- ratio of reaction rate to transport rate.
 *
 * Da_I   = r * V / (Q * C_0)            reaction vs. bulk convection
 * Da_II  = r * L^2 / (D * C_0)          reaction vs. diffusion
 * Da_III = (-DH_r) * r * L / (rho*Cp*u*T) heat generation vs. convection
 * Da_IV  = (-DH_r) * r * L^2 / (k*T)    heat generation vs. conduction
 */
typedef struct {
    double Da_I;    /**< Convection timescale / reaction timescale */
    double Da_II;   /**< Diffusion timescale / reaction timescale */
    double Da_III;  /**< Heat generation / convective heat removal */
    double Da_IV;   /**< Heat generation / conductive heat removal */
} DamkohlerNumbers;

/**
 * @brief Thiele modulus -- ratio of reaction rate to diffusion rate in a porous catalyst.
 *
 * phi = L_c * sqrt(k_v / D_eff)  for 1st-order reaction
 *
 * phi < 0.4 : No diffusion limitation (effectiveness ~ 1)
 * 0.4 < phi < 3.0 : Transition regime
 * phi > 3.0 : Strong diffusion limitation
 *
 * Shape-dependent characteristic length L_c = V_p / S_p
 *   Slab:  L_c = L (half-thickness)
 *   Cylinder: L_c = R/2
 *   Sphere: L_c = R/3
 */
typedef struct {
    double phi;             /**< Thiele modulus (dimensionless) */
    double phi_slab;        /**< phi for slab geometry, L_c = half-thickness */
    double phi_cylinder;    /**< phi for cylinder geometry, L_c = R/2 */
    double phi_sphere;      /**< phi for sphere geometry, L_c = R/3 */
} ThieleModulus;

/**
 * @brief Reynolds number -- ratio of inertial to viscous forces.
 */
typedef struct {
    double Re;          /**< Reynolds number */
    double Re_particle; /**< Particle Reynolds: Re_p = rho*u*d_p/mu */
    double Re_film;     /**< Film Reynolds: Re_f = 4*Gamma/mu (falling film) */
} ReynoldsNumbers;

/**
 * @brief Schmidt number -- ratio of momentum diffusivity to mass diffusivity.
 *
 * Sc = nu / D_AB = mu / (rho * D_AB)
 *
 * Gases:  Sc ~ 0.5 - 2.0
 * Liquids: Sc ~ 10^2 - 10^4
 */
typedef struct {
    double Sc;           /**< Schmidt number, characteristic of fluid pair */
    double Sc_turbulent; /**< Turbulent Schmidt (Prt analog, typically ~0.7-0.9) */
} SchmidtNumber;

/* ---------------------------------------------------------------------------
 * L2: Transport Regime Classification
 * ------------------------------------------------------------------------- */

/** @brief Identifies which transport mechanism dominates. */
typedef enum {
    TRANSPORT_DIFFUSION_DOMINATED,
    TRANSPORT_CONVECTION_DOMINATED,
    TRANSPORT_REACTION_DOMINATED,
    TRANSPORT_MIXED_REGIME,
    TRANSPORT_DISPERSION_REGIME,
    TRANSPORT_MASS_TRANSFER_LIMITED
} TransportRegime;

/**
 * @brief Phase classification for multiphase CDR systems.
 */
typedef enum {
    PHASE_HOMOGENEOUS,
    PHASE_GAS_LIQUID,
    PHASE_GAS_SOLID,
    PHASE_LIQUID_SOLID,
    PHASE_GAS_LIQUID_SOLID,
    PHASE_IMMISCIBLE_LIQUID
} PhaseSystem;

/* ---------------------------------------------------------------------------
 * L3: Material and State Structures
 * ------------------------------------------------------------------------- */

/**
 * @brief Species properties for CDR calculations.
 */
typedef struct {
    const char *name;
    double molecular_weight;
    double diffusion_coeff;
    double diffusion_coeff_eff;
    double henry_constant;
    double heat_capacity;
    int    charge;
} SpeciesProperties;

/**
 * @brief Fluid mixture state for CDR analysis.
 */
typedef struct {
    double temperature;
    double pressure;
    double density;
    double viscosity;
    double thermal_conductivity;
    double heat_capacity;
    size_t n_species;
} MixtureState;

/**
 * @brief Reactor geometry parameters.
 */
typedef enum {
    GEOM_SLAB,
    GEOM_CYLINDER,
    GEOM_SPHERE,
    GEOM_ARBITRARY_3D
} CatalystGeometry;

/**
 * @brief CDR system configuration combining all elements.
 */
typedef struct {
    MixtureState      mixture;
    PecletNumbers     peclet;
    DamkohlerNumbers  damkohler;
    ThieleModulus     thiele;
    ReynoldsNumbers   reynolds;
    SchmidtNumber     schmidt;
    TransportRegime   regime;
    PhaseSystem       phase;
    double            characteristic_L;
    double            characteristic_u;
    double            residence_time;
    double            porosity;
    double            tortuosity;
} CDRSystem;

/* ---------------------------------------------------------------------------
 * L4: Dimensionless Number Calculators
 * ------------------------------------------------------------------------- */

double cdr_peclet_mass(double velocity, double length, double diffusivity);
double cdr_peclet_heat(double velocity, double length, double alpha);
double cdr_damkohler_I(double rate_constant, double residence_time);
double cdr_damkohler_II(double rate_constant, double length, double diffusivity);
double cdr_thiele_modulus(double Lc, double rate_constant, double diffusivity);
double cdr_schmidt(double viscosity, double density, double diffusivity);
double cdr_reynolds(double density, double velocity, double length, double viscosity);
TransportRegime cdr_classify_regime(double pe, double da_I, double da_II);
double cdr_characteristic_length(CatalystGeometry geometry, double radius);
void cdr_system_init(CDRSystem *sys);
void cdr_compute_dimensionless(CDRSystem *sys);

#endif /* CDR_CORE_H */

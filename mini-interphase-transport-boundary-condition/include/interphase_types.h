/**
 * interphase_types.h
 * ==================
 * Core data structures for interphase transport and boundary conditions.
 *
 * Reference: Bird, Stewart, Lightfoot (2007) "Transport Phenomena"
 *            Chapter 2 (Shell Momentum Balances & Boundary Conditions)
 *            Chapter 19 (Equations of Change for Multicomponent Systems)
 *            Chapter 22 (Interphase Transport in Isothermal Systems)
 *
 * Knowledge coverage: L1 Definitions, L2 Core Concepts, L3 Engineering Quantities
 */

#ifndef INTERPHASE_TYPES_H
#define INTERPHASE_TYPES_H

#include <stddef.h>
#include <stdint.h>
#include "interphase_config.h"

/* ============================================================================
 * L1: Core Definitions - Interface Taxonomy
 * ============================================================================ */

/** Phase identifiers for interphase transport classification */
typedef enum {
    PHASE_GAS        = 0,  /**< Gaseous phase */
    PHASE_LIQUID     = 1,  /**< Liquid phase */
    PHASE_SOLID      = 2,  /**< Solid phase */
    PHASE_PLASMA     = 3,  /**< Plasma/ionized phase */
    PHASE_SUPERCRITICAL = 4 /**< Supercritical fluid phase */
} PhaseType;

/** Boundary condition category from mathematical physics taxonomy.
 *  After Hadamard (1902) classification of PDE boundary conditions.
 */
typedef enum {
    BC_DIRICHLET    = 0,  /**< Fixed value: u = u_b at boundary */
    BC_NEUMANN      = 1,  /**< Fixed gradient: du/dn = q at boundary */
    BC_ROBIN        = 2,  /**< Mixed/Convective: a*u + b*du/dn = c */
    BC_JUMP         = 3,  /**< Discontinuous jump across interface */
    BC_PERIODIC     = 4,  /**< Periodic boundary condition */
    BC_SYMMETRY     = 5,  /**< Symmetry/reflection boundary condition */
    BC_OUTFLOW      = 6,  /**< Outflow/convective boundary condition */
    BC_CAUCHY       = 7   /**< Both value and gradient specified */
} BCType;

/** Interface geometry classification.
 *  Critical for determining whether curvature effects dominate.
 */
typedef enum {
    IFACE_PLANAR      = 0,  /**< Flat interface (zero curvature) */
    IFACE_CYLINDRICAL = 1,  /**< Cylindrical interface */
    IFACE_SPHERICAL   = 2,  /**< Spherical interface (droplets, bubbles) */
    IFACE_ELLIPSOIDAL = 3,  /**< Ellipsoidal (deformed drops) */
    IFACE_ARBITRARY   = 4,  /**< Arbitrary curved surface */
    IFACE_CORRUGATED  = 5,  /**< Periodic wavy/rough interface */
    IFACE_FRACTAL     = 6   /**< Fractal/self-affine rough interface */
} InterfaceGeometry;

/** Phase pair classification for transport regime identification. */
typedef enum {
    PAIR_GAS_LIQUID     = 0,  /**< Absorption, stripping, evaporation */
    PAIR_LIQUID_LIQUID  = 1,  /**< Extraction, emulsification */
    PAIR_GAS_SOLID      = 2,  /**< Adsorption, drying, sublimation */
    PAIR_LIQUID_SOLID   = 3,  /**< Dissolution, crystallization, corrosion */
    PAIR_LIQUID_VAPOR   = 4,  /**< Boiling, condensation, distillation */
    PAIR_SOLID_SOLID    = 5,  /**< Diffusion bonding, sintering */
    PAIR_GAS_PLASMA     = 6,  /**< Plasma processing, sputtering */
    PAIR_IMMISCIBLE_LIQ = 7   /**< Liquid-liquid with sharp interface */
} PhasePair;

/* ============================================================================
 * L1: Thermodynamic State at Interface
 * ============================================================================ */

/** Single-component state at a point on/beside the interface.
 *  Represents intensive thermodynamic variables at one phase side.
 */
typedef struct {
    double temperature;
    double pressure;
    double density;
    double viscosity;
    double thermal_conductivity;
    double specific_heat_cp;
    double speed_of_sound;
} PhaseState;

/** Multicomponent extension for mixture state at interface. */
typedef struct {
    PhaseState bulk;
    size_t    n_species;
    double   *mass_fraction;
    double   *mole_fraction;
    double   *diffusivity;
    double    molar_mass;
} MixtureState;

/* ============================================================================
 * L1: Interface Descriptor
 * ============================================================================ */

typedef struct {
    PhasePair        pair;
    InterfaceGeometry geometry;
    PhaseType        phase_a;
    PhaseType        phase_b;
    PhaseState       state_a;
    PhaseState       state_b;
    PhaseState       state_ai;
    PhaseState       state_bi;
    double           surface_tension;
    double           surface_energy;
    double           contact_angle;
    double           interfacial_area;
    double           curvature_mean;
    double           curvature_gaussian;
    double           roughness_rms;
    double           slip_length;
    uint32_t         flags;
} InterfaceDescriptor;

#define IFACE_FLAG_WETTING       (1u << 0)
#define IFACE_FLAG_NONWETTING    (1u << 1)
#define IFACE_FLAG_PINNED        (1u << 2)
#define IFACE_FLAG_MOVING        (1u << 3)
#define IFACE_FLAG_PHASE_CHANGE  (1u << 4)
#define IFACE_FLAG_CHEMICAL_RXN  (1u << 5)
#define IFACE_FLAG_SURFACTANT    (1u << 6)
#define IFACE_FLAG_CHARGED       (1u << 7)

/* ============================================================================
 * L2: Core Concepts - Flux Vectors and Transport Coefficients
 * ============================================================================ */

typedef struct {
    double tau_xx, tau_xy, tau_xz;
    double tau_yx, tau_yy, tau_yz;
    double tau_zx, tau_zy, tau_zz;
    double pressure;
} StressTensor;

typedef struct {
    double qx, qy, qz;
    double htc;
    double thermal_resistance;
    double nusselt_number;
} HeatFlux;

typedef struct {
    double jx, jy, jz;
    double mtc;
    double sherwood_number;
    double mass_transfer_resistance;
} MassFlux;

typedef struct {
    HeatFlux heat;
    MassFlux mass;
    StressTensor momentum;
    double    mass_flux_total;
    double    energy_flux_total;
} InterphaseFlux;

/* ============================================================================
 * L2: Interfacial Resistance Models
 * ============================================================================ */

typedef struct {
    double film_thickness_a;
    double film_thickness_b;
    double k_a;
    double k_b;
    double K_overall_a;
    double K_overall_b;
    double h_a;
    double h_b;
    double U_overall;
} TwoFilmModel;

typedef struct {
    double exposure_time;
    double diffusivity_a;
    double diffusivity_b;
    double k_penetration_a;
    double k_penetration_b;
} PenetrationModel;

typedef struct {
    double renewal_rate;
    double diffusivity;
    double k_renewal;
    double mean_element_age;
} SurfaceRenewalModel;

/* ============================================================================
 * L3: Engineering Quantities - Dimensionless Groups
 * ============================================================================ */

typedef struct {
    double Reynolds;
    double Prandtl;
    double Schmidt;
    double Nusselt;
    double Sherwood;
    double Peclet_heat;
    double Peclet_mass;
    double Grashof;
    double Rayleigh;
    double Stanton_heat;
    double Stanton_mass;
    double Biot;
    double Fourier;
    double Weber;
    double Capillary;
    double Bond;
    double Marangoni;
    double Mach;
    double Knudsen;
    double Lewis;
} DimensionlessGroups;

/* ============================================================================
 * L2: Jump Condition Data
 * ============================================================================ */

typedef struct {
    double mass_jump;
    double momentum_jump[3];
    double energy_jump;
    double species_jump;
    double entropy_jump;
} JumpConditions;

/* ============================================================================
 * L2: Henry's Law and Vapor-Liquid Equilibrium
 * ============================================================================ */

typedef struct {
    double henry_coefficient;
    double temperature;
    double henry_dHdT;
    double henry_coefficient_log;
} HenryLaw;

typedef struct {
    double saturation_pressure;
    double activity_coefficient_a;
    double activity_coefficient_b;
    double partition_coefficient;
    double relative_volatility;
} VaporLiquidEquilibrium;

/* ============================================================================
 * L3: Engineering Data for Common Systems
 * ============================================================================ */

typedef struct {
    char     system_name[64];
    PhasePair pair;
    double   T_range[2];
    double   henry_const_298;
    double   diffusivity_gas;
    double   diffusivity_liquid;
    double   surface_tension_298;
    double   viscosity_ratio;
    double   density_ratio;
    double   schmidt_number_gas;
    double   schmidt_number_liquid;
} SystemProperties;

/* ============================================================================
 * L3: Transport Property Estimation Parameters
 * ============================================================================ */

typedef struct {
    double sigma_AB;
    double epsilon_over_k_AB;
    double omega_D;
    double molar_mass_A;
    double molar_mass_B;
} DiffusivityParams;

typedef struct {
    double sigma;
    double epsilon_over_k;
    double omega_v;
    double molar_mass;
    double acentric_factor;
    double dipole_moment;
    double association_factor;
} ViscosityParams;

typedef struct {
    double cp_ideal;
    double cv_ideal;
    double gamma;
    double molar_mass;
    double acentric_factor;
} ThermalConductivityParams;

/* ============================================================================
 * L4: Conservation Law Data Structures
 * ============================================================================ */

typedef struct {
    double volume;
    double area_interface;
    double volume_a;
    double volume_b;
    double normal[3];
    double velocity_interface[3];
} InterfaceControlVolume;

/* ============================================================================
 * L1: Boundary Conditions for Specific Fields
 * ============================================================================ */

typedef struct {
    BCType type;
    double u_wall[3];
    double shear_stress[3];
    double slip_coefficient;
    double transpiration;
} VelocityBC;

typedef struct {
    BCType type;
    double T_wall;
    double heat_flux;
    double h_ext;
    double T_ext;
    double emissivity;
    double T_surroundings;
} ThermalBC;

typedef struct {
    BCType type;
    double C_wall;
    double mass_flux;
    double k_ext;
    double C_ext;
    double reaction_rate_constant;
    double henry_constant;
    double partition_coefficient;
} ConcentrationBC;

typedef struct {
    VelocityBC       velocity;
    ThermalBC        thermal;
    ConcentrationBC   concentration;
    double           pressure_outlet;
    double           turbulence_intensity;
    double           turbulence_length_scale;
} BoundaryConditionSet;

/* ============================================================================
 * Function Declarations for interphase_types.c
 * ============================================================================ */

#ifdef __cplusplus
extern "C" {
#endif

/* PhaseState operations */
void phase_state_init(PhaseState *state, PhaseType phase);
double phase_state_prandtl(const PhaseState *state);
double phase_state_kinematic_viscosity(const PhaseState *state);
double phase_state_thermal_diffusivity(const PhaseState *state);

/* InterfaceDescriptor operations */
void interface_descriptor_init(InterfaceDescriptor *iface, PhasePair pair, InterfaceGeometry geometry);
void interface_set_surface_tension(InterfaceDescriptor *iface, double sigma);
void interface_set_roughness(InterfaceDescriptor *iface, double rms);
double interface_capillary_length(const InterfaceDescriptor *iface, double g);

/* Dimensionless number computation */
double compute_prandtl(double mu, double cp, double k);
double compute_schmidt(double mu, double rho, double D);
double compute_lewis(double alpha, double D);
void compute_dimensionless_groups(double rho, double U, double L,
                                   double mu, double k, double cp,
                                   double D, double beta,
                                   double delta_T, double sigma,
                                   double g, double c, double lambda,
                                   DimensionlessGroups *groups);

/* StressTensor operations */
void stress_tensor_init(StressTensor *tau, double pressure);
void stress_tensor_newtonian(StressTensor *tau, double mu, const double grad_u[9]);
double stress_tensor_magnitude(const StressTensor *tau);

/* Henry's law operations */
void henry_law_at_temperature(const HenryLaw *hl, double T, double *H_at_T);

/* VLE operations */
void raoult_vapor_composition(const VaporLiquidEquilibrium *vle,
                               double x_liquid, double total_pressure,
                               double *y_vapor);
void bubble_point_pressure(const VaporLiquidEquilibrium *vle,
                            double x_liquid, double *P_bubble);
void dew_point_composition(const VaporLiquidEquilibrium *vle,
                            double y_vapor, double total_pressure,
                            double *x_liquid);

/* Transport property estimation */
void chapman_enskog_diffusivity(const DiffusivityParams *params,
                                 double T, double p, double *D_AB);
void chung_viscosity(const ViscosityParams *params, double T, double *mu);
void eucken_thermal_conductivity(double mu, double cp, double cv,
                                  double molar_mass, double *k);

#ifdef __cplusplus
}
#endif

#endif /* INTERPHASE_TYPES_H */

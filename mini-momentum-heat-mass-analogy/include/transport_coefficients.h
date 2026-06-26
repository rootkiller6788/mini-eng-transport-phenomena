/**
 * @file transport_coefficients.h
 * @brief Transport coefficients for momentum, heat, and mass transfer
 *
 * Defines the three fundamental transport coefficients and their
 * temperature/pressure dependencies. Based on the kinetic theory of
 * gases and corresponding states principle.
 *
 * Reference: Bird, Stewart, Lightfoot (2007), Transport Phenomena, 2nd Ed.
 *            Chapters 1, 9, 16, 17.
 *
 * Key definitions:
 *   - Dynamic viscosity μ [Pa·s]: momentum diffusivity
 *   - Thermal conductivity k [W/(m·K)]: heat diffusivity
 *   - Mass diffusivity D [m²/s]: species diffusivity
 *
 * The analogy rests on the mathematical identity of the governing equations:
 *   τ_yx = -μ · dv_x/dy    (Newton's law of viscosity)
 *   q_y  = -k · dT/dy      (Fourier's law of heat conduction)
 *   J_Ay = -D_AB · dC_A/dy (Fick's first law of diffusion)
 */

#ifndef TRANSPORT_COEFFICIENTS_H
#define TRANSPORT_COEFFICIENTS_H

#include <stddef.h>

/* ─── Fundamental Transport Coefficient Structures ─────────────── */

/** Sutherland's law coefficients for gas viscosity */
typedef struct {
    double mu0;    /**< Reference viscosity at T0 [Pa·s] */
    double T0;     /**< Reference temperature [K] */
    double S;      /**< Sutherland constant [K] */
} SutherlandParams;

/** Power-law model coefficients for liquid viscosity */
typedef struct {
    double mu_ref; /**< Reference viscosity [Pa·s] */
    double T_ref;  /**< Reference temperature [K] */
    double n;      /**< Temperature exponent (dimensionless) */
} PowerLawViscosityParams;

/** Transport coefficient triplet for a substance at a given state */
typedef struct {
    double mu;     /**< Dynamic viscosity [Pa·s] */
    double k;      /**< Thermal conductivity [W/(m·K)] */
    double D;      /**< Mass diffusivity (binary) [m²/s] */
    double rho;    /**< Density [kg/m³] */
    double cp;     /**< Specific heat capacity [J/(kg·K)] */
    double T;      /**< Temperature [K] */
    double P;      /**< Pressure [Pa] */
} TransportState;

/** Kinetic theory parameters for gases */
typedef struct {
    double sigma;       /**< Collision diameter [m] (Lennard-Jones) */
    double epsilon_kB;  /**< Energy parameter / Boltzmann constant [K] */
    double M;           /**< Molecular weight [kg/mol] */
    double omega_mu;    /**< Collision integral for viscosity */
    double omega_k;     /**< Collision integral for thermal conductivity */
} KineticTheoryParams;

/** Corresponding-states parameters for estimating transport properties */
typedef struct {
    double Tc;       /**< Critical temperature [K] */
    double Pc;       /**< Critical pressure [Pa] */
    double Vc;       /**< Critical volume [m³/mol] */
    double Zc;       /**< Critical compressibility factor */
    double omega;    /**< Pitzer acentric factor */
    double mu_c;     /**< Critical viscosity (estimated) [Pa·s] */
    double k_c;      /**< Critical thermal conductivity (estimated) [W/(m·K)] */
} CorrespondingStates;

/** Mixture transport properties with composition */
typedef struct {
    size_t n_species;         /**< Number of species */
    double *mole_fractions;   /**< Mole fractions [n_species] */
    double *mu_species;       /**< Pure-component viscosities [Pa·s] */
    double *k_species;        /**< Pure-component conductivities [W/(m·K)] */
    double **D_binary;        /**< Binary diffusivity matrix [m²/s] (n×n) */
    double mu_mix;            /**< Mixture viscosity [Pa·s] */
    double k_mix;             /**< Mixture thermal conductivity [W/(m·K)] */
} MixtureTransport;

/* ─── L1: Core Definitions ─────────────────────────────────────── */

/** Viscosity definition and types.
 *  Dynamic viscosity μ [Pa·s] = shear stress / velocity gradient.
 *  Kinematic viscosity ν = μ/ρ [m²/s] = momentum diffusivity.
 */
typedef struct {
    double mu;             /**< Dynamic viscosity [Pa·s] */
    double nu;             /**< Kinematic viscosity [m²/s] */
    double rho;            /**< Density [kg/m³] */
    double T;              /**< Temperature [K] */
    int is_newtonian;      /**< 1 if Newtonian, 0 if non-Newtonian */
} ViscosityDefinition;

/** Thermal conductivity definition.
 *  k [W/(m·K)] = heat flux / temperature gradient.
 *  Thermal diffusivity α = k/(ρ·cp) [m²/s] = heat diffusivity.
 */
typedef struct {
    double k;              /**< Thermal conductivity [W/(m·K)] */
    double alpha;          /**< Thermal diffusivity [m²/s] */
    double rho;            /**< Density [kg/m³] */
    double cp;             /**< Specific heat [J/(kg·K)] */
    double T;              /**< Temperature [K] */
} ThermalConductivityDef;

/** Mass diffusivity definition.
 *  D_AB [m²/s] = mass flux / concentration gradient (Fick's law).
 *  For ideal gases, D_AB ∝ T^(3/2) / P.
 */
typedef struct {
    double D_AB;           /**< Binary mass diffusivity [m²/s] */
    double M_A;            /**< Molecular weight of species A [kg/mol] */
    double M_B;            /**< Molecular weight of species B [kg/mol] */
    double T;              /**< Temperature [K] */
    double P;              /**< Pressure [Pa] */
} MassDiffusivityDef;

/* ─── L2: Transport Property Models ────────────────────────────── */

/**
 * Computes gas viscosity using Sutherland's law.
 * Sutherland (1893): μ = μ₀ · (T/T₀)^(3/2) · (T₀+S)/(T+S)
 * Complexity: O(1)
 */
double sutherland_viscosity(double T, const SutherlandParams *params);

/**
 * Computes gas viscosity using the kinetic theory (Chapman-Enskog).
 * μ = (5/16) · sqrt(π·M·kB·T) / (π·σ²·Ω_μ)
 * Complexity: O(1)
 */
double chapman_enskog_viscosity(double T, const KineticTheoryParams *params);

/**
 * Computes liquid viscosity using the power-law model.
 * μ = μ_ref · (T/T_ref)^n
 * Complexity: O(1)
 */
double power_law_viscosity(double T, const PowerLawViscosityParams *params);

/**
 * Computes gas thermal conductivity from kinetic theory (Eucken formula).
 * k = μ · (Cp + 5/4 · R)
 * Eucken (1913): modified for polyatomic gases.
 * Complexity: O(1)
 */
double eucken_thermal_conductivity(double mu, double cp, double M);

/**
 * Computes binary gas diffusivity using Chapman-Enskog theory.
 * D_AB = (3/16)·sqrt(2π·kB³·T³/M_AB) / (P·π·σ_AB²·Ω_D)
 * where M_AB = 2/(1/M_A + 1/M_B) (reduced mass)
 * Complexity: O(1)
 */
double chapman_enskog_diffusivity(double T, double P,
                                  const KineticTheoryParams *gasA,
                                  const KineticTheoryParams *gasB);

/**
 * Computes liquid diffusivity using Wilke-Chang correlation.
 * D_AB = 7.4×10⁻⁸ · (φ·M_B)^(1/2) · T / (μ_B · V_A^(0.6))
 * where φ = association factor, V_A = molar volume at boiling point.
 * Complexity: O(1)
 */
double wilke_chang_diffusivity(double T, double mu_B, double phi_MB,
                               double V_A);

/**
 * Computes thermal diffusivity α = k/(ρ·cp) [m²/s].
 * This is the key quantity in the momentum-heat analogy.
 * Complexity: O(1)
 */
double thermal_diffusivity(double k, double rho, double cp);

/**
 * Computes momentum diffusivity (kinematic viscosity) ν = μ/ρ [m²/s].
 * This is the key quantity in the momentum-heat-mass analogy.
 * Complexity: O(1)
 */
double momentum_diffusivity(double mu, double rho);

/**
 * Computes mass diffusivity for gases using Fuller-Schettler-Giddings method.
 * D_AB = 1.0×10⁻³ · T^(1.75) · sqrt(1/M_A + 1/M_B) / [P · (Σv_A^(1/3) + Σv_B^(1/3))²]
 * where Σv = atomic diffusion volumes.
 * Complexity: O(1)
 */
double fuller_diffusivity(double T, double P,
                          double M_A, double M_B,
                          double sum_v_A, double sum_v_B);

/* ─── L3: Engineering Quantities ───────────────────────────────── */

/**
 * Returns typical viscosity of air at given temperature [Pa·s].
 * Uses Sutherland's law with standard air parameters.
 * Benchmark: μ_air(293K) ≈ 1.82×10⁻⁵ Pa·s
 */
double air_viscosity(double T);

/**
 * Returns typical viscosity of water at given temperature [Pa·s].
 * Uses empirical correlation valid for 273K < T < 373K.
 * Benchmark: μ_water(293K) ≈ 1.00×10⁻³ Pa·s
 */
double water_viscosity(double T);

/**
 * Returns typical thermal conductivity of air [W/(m·K)].
 * Benchmark: k_air(300K) ≈ 0.026 W/(m·K)
 */
double air_thermal_conductivity(double T);

/**
 * Returns typical thermal conductivity of water [W/(m·K)].
 * Benchmark: k_water(300K) ≈ 0.61 W/(m·K)
 */
double water_thermal_conductivity(double T);

/**
 * Returns binary diffusivity of water vapor in air [m²/s].
 * Benchmark: D_H2O-air(298K, 1atm) ≈ 2.6×10⁻⁵ m²/s
 */
double water_vapor_air_diffusivity(double T, double P);

/* ─── L4: Conservation Law Forms ───────────────────────────────── */

/**
 * Computes shear stress from Newton's law of viscosity.
 * τ_yx = -μ · dv_x/dy  [Pa]
 * Returns the magnitude of shear stress.
 * Complexity: O(1)
 */
double newtons_law_shear_stress(double mu, double dv_dy);

/**
 * Computes heat flux from Fourier's law.
 * q_y = -k · dT/dy  [W/m²]
 * Returns the magnitude of heat flux.
 * Complexity: O(1)
 */
double fouriers_law_heat_flux(double k, double dT_dy);

/**
 * Computes mass flux from Fick's first law.
 * J_Ay = -D_AB · dC_A/dy  [mol/(m²·s)]
 * Returns the magnitude of molar flux.
 * Complexity: O(1)
 */
double ficks_law_mass_flux(double D_AB, double dC_dy);

/**
 * Computes the generalized transport flux.
 * Flux = -(transport coefficient) × (driving force gradient)
 * Unifies all three transport laws into a single algebraic form.
 * Complexity: O(1)
 */
double generalized_transport_flux(double transport_coeff,
                                  double driving_force_gradient);

/* ─── Mixture Rules ────────────────────────────────────────────── */

/**
 * Wilke's mixing rule for gas mixture viscosity.
 * μ_mix = Σ (x_i·μ_i / Σ_j x_j·Φ_ij)
 * where Φ_ij is the interaction parameter.
 * Complexity: O(n²) where n is number of species.
 */
double wilke_mixture_viscosity(size_t n, const double *x, const double *mu,
                               const double *M);

/**
 * Mason-Saxena mixing rule for gas mixture thermal conductivity.
 * Similar form to Wilke's rule but with different interaction parameters.
 * Complexity: O(n²)
 */
double mason_saxena_mixture_conductivity(size_t n, const double *x,
                                         const double *k, const double *M);

/**
 * Computes the interaction parameter Φ_ij for Wilke's mixing rule.
 * Φ_ij = [1 + (μ_i/μ_j)^(1/2)·(M_j/M_i)^(1/4)]² / [8·(1 + M_i/M_j)]^(1/2)
 * Complexity: O(1)
 */
double wilke_interaction_parameter(double mu_i, double mu_j,
                                   double M_i, double M_j);

/**
 * Computes the complete transport state for a pure substance.
 * Fills in all fields of TransportState from basic inputs.
 * Complexity: O(1)
 */
void compute_transport_state(const PowerLawViscosityParams *visc_params,
                             double k, double D, double rho,
                             double cp, double T, double P,
                             TransportState *state);

/**
 * Validates that transport coefficients are physically plausible.
 * Returns 0 if valid, error code otherwise.
 * Checks: μ>0, k>0, D>0, ρ>0, cp>0, T>0 (Kelvin), P>0.
 */
int validate_transport_state(const TransportState *state);

/**
 * Prints a human-readable summary of transport coefficients.
 * Provides engineering intuition (e.g., "gas-like", "liquid-like").
 */
void print_transport_summary(const TransportState *state);

#endif /* TRANSPORT_COEFFICIENTS_H */

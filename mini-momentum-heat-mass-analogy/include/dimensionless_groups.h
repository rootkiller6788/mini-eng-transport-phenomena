/**
 * @file dimensionless_groups.h
 * @brief Dimensionless groups for the momentum-heat-mass analogy
 *
 * The analogy is expressed through the equality of dimensionless groups
 * that characterize the relative importance of transport mechanisms.
 *
 * Reference: Bird, Stewart, Lightfoot (2007), Transport Phenomena,
 *            Chapter 10 (dimensional analysis), Chapters 12-22.
 *            Incropera & DeWitt (2007), Fundamentals of Heat and Mass Transfer.
 *
 * Key dimensionless groups connecting the three transport modes:
 *
 *   Momentum ↔ Heat:     Pr = ν/α = μ·Cp/k
 *   Momentum ↔ Mass:     Sc = ν/D
 *   Heat ↔ Mass:         Le = α/D = Sc/Pr
 *
 *   Friction ↔ Heat:     St = Nu/(Re·Pr) = h/(ρ·v·Cp)
 *   Friction ↔ Mass:     St_m = Sh/(Re·Sc) = k_c/v
 */

#ifndef DIMENSIONLESS_GROUPS_H
#define DIMENSIONLESS_GROUPS_H

#include <stddef.h>

/* ─── L1: Primary Dimensionless Number Definitions ──────────────── */

/** Reynolds number Re = ρ·v·L/μ = v·L/ν.
 *  Ratio of inertial to viscous forces.
 *  Determines flow regime: laminar (Re<2300 pipe), transitional, turbulent.
 */
typedef struct {
    double Re;           /**< Reynolds number */
    double rho;          /**< Density [kg/m³] */
    double v;            /**< Characteristic velocity [m/s] */
    double L;            /**< Characteristic length [m] */
    double mu;           /**< Dynamic viscosity [Pa·s] */
    double nu;           /**< Kinematic viscosity [m²/s] */
} ReynoldsNumber;

/** Prandtl number Pr = ν/α = μ·Cp/k.
 *  Ratio of momentum diffusivity to thermal diffusivity.
 *  Determines relative thickness of velocity vs. thermal boundary layers.
 *  Gases: Pr ≈ 0.7, Water: Pr ≈ 7, Oils: Pr ≈ 100-10000.
 */
typedef struct {
    double Pr;           /**< Prandtl number */
    double nu;           /**< Kinematic viscosity [m²/s] */
    double alpha;        /**< Thermal diffusivity [m²/s] */
    double mu;           /**< Dynamic viscosity [Pa·s] */
    double cp;           /**< Specific heat [J/(kg·K)] */
    double k;            /**< Thermal conductivity [W/(m·K)] */
} PrandtlNumber;

/** Schmidt number Sc = ν/D = μ/(ρ·D).
 *  Ratio of momentum diffusivity to mass diffusivity.
 *  Determines relative thickness of velocity vs. concentration BLs.
 *  Gases: Sc ≈ 0.2-5, Liquids: Sc ≈ 100-10000.
 */
typedef struct {
    double Sc;           /**< Schmidt number */
    double nu;           /**< Kinematic viscosity [m²/s] */
    double D;            /**< Mass diffusivity [m²/s] */
    double mu;           /**< Dynamic viscosity [Pa·s] */
    double rho;          /**< Density [kg/m³] */
} SchmidtNumber;

/** Lewis number Le = α/D = Sc/Pr.
 *  Ratio of thermal diffusivity to mass diffusivity.
 *  Determines relative thickness of thermal vs. concentration BLs.
 *  Gases: Le ≈ 1, Liquids: Le ≈ 100.
 */
typedef struct {
    double Le;           /**< Lewis number */
    double alpha;        /**< Thermal diffusivity [m²/s] */
    double D;            /**< Mass diffusivity [m²/s] */
    double Sc;           /**< Schmidt number */
    double Pr;           /**< Prandtl number */
} LewisNumber;

/** Nusselt number Nu = h·L/k.
 *  Dimensionless temperature gradient at the wall.
 *  Ratio of convective to conductive heat transfer.
 */
typedef struct {
    double Nu;           /**< Nusselt number */
    double h;            /**< Heat transfer coefficient [W/(m²·K)] */
    double L;            /**< Characteristic length [m] */
    double k;            /**< Thermal conductivity [W/(m·K)] */
} NusseltNumber;

/** Sherwood number Sh = k_c·L/D.
 *  Dimensionless concentration gradient at the wall (mass transfer Nu).
 *  Ratio of convective to diffusive mass transfer.
 */
typedef struct {
    double Sh;           /**< Sherwood number */
    double kc;           /**< Mass transfer coefficient [m/s] */
    double L;            /**< Characteristic length [m] */
    double D;            /**< Mass diffusivity [m²/s] */
} SherwoodNumber;

/** Stanton number St = Nu/(Re·Pr) = h/(ρ·v·Cp).
 *  Ratio of heat transferred to fluid thermal capacity.
 *  Key quantity in Reynolds analogy: St = f/2.
 */
typedef struct {
    double St;           /**< Stanton number */
    double Nu;           /**< Nusselt number */
    double Re;           /**< Reynolds number */
    double Pr;           /**< Prandtl number */
    double h;            /**< Heat transfer coefficient [W/(m²·K)] */
    double rho;          /**< Density [kg/m³] */
    double v;            /**< Velocity [m/s] */
    double cp;           /**< Specific heat [J/(kg·K)] */
} StantonNumber;

/** Péclet number Pe = Re·Pr = v·L/α.
 *  Ratio of advective to diffusive heat transport.
 */
typedef struct {
    double Pe;           /**< Péclet number */
    double Re;           /**< Reynolds number */
    double Pr;           /**< Prandtl number */
    double v;            /**< Velocity [m/s] */
    double L;            /**< Characteristic length [m] */
    double alpha;        /**< Thermal diffusivity [m²/s] */
} PecletNumber;

/** Grashof number Gr = g·β·ΔT·L³/ν².
 *  Ratio of buoyancy to viscous forces (natural convection).
 */
typedef struct {
    double Gr;           /**< Grashof number */
    double g;            /**< Gravitational acceleration [m/s²] */
    double beta;         /**< Thermal expansion coefficient [1/K] */
    double delta_T;      /**< Temperature difference [K] */
    double L;            /**< Characteristic length [m] */
    double nu;           /**< Kinematic viscosity [m²/s] */
} GrashofNumber;

/** Rayleigh number Ra = Gr·Pr = g·β·ΔT·L³/(ν·α).
 *  Determines natural convection regime.
 */
typedef struct {
    double Ra;           /**< Rayleigh number */
    double Gr;           /**< Grashof number */
    double Pr;           /**< Prandtl number */
} RayleighNumber;

/** Brinkman number Br = μ·v²/(k·ΔT).
 *  Ratio of viscous dissipation to conductive heat transfer.
 *  Important in high-speed flows and polymer processing.
 */
typedef struct {
    double Br;           /**< Brinkman number */
    double mu;           /**< Dynamic viscosity [Pa·s] */
    double v;            /**< Velocity [m/s] */
    double k;            /**< Thermal conductivity [W/(m·K)] */
    double delta_T;      /**< Temperature difference [K] */
} BrinkmanNumber;

/** Eckert number Ec = v²/(cp·ΔT) = Br/Pr.
 *  Ratio of kinetic energy to enthalpy difference.
 *  Important in compressible flow heat transfer.
 */
typedef struct {
    double Ec;           /**< Eckert number */
    double v;            /**< Velocity [m/s] */
    double cp;           /**< Specific heat [J/(kg·K)] */
    double delta_T;      /**< Temperature difference [K] */
} EckertNumber;

/** Biot number Bi = h·L_c/k_solid.
 *  Ratio of internal conduction resistance to external convection resistance.
 *  Bi < 0.1: lumped capacitance valid.
 */
typedef struct {
    double Bi;           /**< Biot number */
    double h;            /**< Heat transfer coefficient [W/(m²·K)] */
    double L_c;          /**< Characteristic length of solid [m] */
    double k_solid;      /**< Thermal conductivity of solid [W/(m·K)] */
} BiotNumber;

/** Fourier number Fo = α·t/L².
 *  Dimensionless time in transient conduction.
 */
typedef struct {
    double Fo;           /**< Fourier number */
    double alpha;        /**< Thermal diffusivity [m²/s] */
    double t;            /**< Time [s] */
    double L;            /**< Characteristic length [m] */
} FourierNumber;

/** Dimensionless group set: all groups for a given flow+fluid state. */
typedef struct {
    double Re; double Pr; double Sc; double Le;
    double Nu; double Sh; double St_heat; double St_mass;
    double Pe; double Gr; double Ra; double Br; double Ec;
    double Bi; double Fo;
    double f_F;  /**< Fanning friction factor */
    double Cf;   /**< Skin friction coefficient */
} DimensionlessGroupSet;

/* ─── L1: Computing Dimensionless Numbers ──────────────────────── */

/**
 * Computes Reynolds number from primitive variables.
 * Re = ρ·v·L/μ
 * Complexity: O(1)
 */
double compute_Re(double rho, double v, double L, double mu);

/**
 * Computes Prandtl number.
 * Pr = μ·Cp/k = ν/α
 * Complexity: O(1)
 */
double compute_Pr(double mu, double cp, double k);

/**
 * Computes Schmidt number.
 * Sc = μ/(ρ·D) = ν/D
 * Complexity: O(1)
 */
double compute_Sc(double mu, double rho, double D);

/**
 * Computes Lewis number.
 * Le = α/D = k/(ρ·cp·D)
 * Complexity: O(1)
 */
double compute_Le(double alpha, double D);

/**
 * Computes Nusselt number.
 * Nu = h·L/k
 * Complexity: O(1)
 */
double compute_Nu(double h, double L, double k);

/**
 * Computes Sherwood number.
 * Sh = k_c·L/D
 * Complexity: O(1)
 */
double compute_Sh(double kc, double L, double D);

/**
 * Computes Stanton number for heat transfer.
 * St = h/(ρ·v·Cp) = Nu/(Re·Pr)
 * Complexity: O(1)
 */
double compute_St_heat(double h, double rho, double v, double cp);

/**
 * Computes Stanton number for mass transfer.
 * St_m = k_c/v = Sh/(Re·Sc)
 * Complexity: O(1)
 */
double compute_St_mass(double kc, double v);

/**
 * Computes Péclet number.
 * Pe = Re·Pr = v·L/α
 * Complexity: O(1)
 */
double compute_Pe(double Re, double Pr);

/**
 * Computes Grashof number.
 * Gr = g·β·ΔT·L³/ν²
 * Complexity: O(1)
 */
double compute_Gr(double g, double beta, double delta_T, double L, double nu);

/**
 * Computes Rayleigh number.
 * Ra = Gr·Pr
 * Complexity: O(1)
 */
double compute_Ra(double Gr, double Pr);

/**
 * Computes Brinkman number.
 * Br = μ·v²/(k·ΔT)
 * Complexity: O(1)
 */
double compute_Br(double mu, double v, double k, double delta_T);

/**
 * Computes Eckert number.
 * Ec = v²/(cp·ΔT)
 * Complexity: O(1)
 */
double compute_Ec(double v, double cp, double delta_T);

/**
 * Computes Biot number.
 * Bi = h·L_c/k_solid
 * Complexity: O(1)
 */
double compute_Bi(double h, double L_c, double k_solid);

/**
 * Computes Fourier number.
 * Fo = α·t/L²
 * Complexity: O(1)
 */
double compute_Fo(double alpha, double t, double L);

/**
 * Computes the full set of dimensionless groups for a given state.
 * Fills all fields of DimensionlessGroupSet.
 *
 * Inputs:
 *   rho, v, L, mu, cp, k, D, h, kc, g, beta, delta_T, t
 *
 * Complexity: O(1)
 */
void compute_dimensionless_group_set(double rho, double v, double L,
                                     double mu, double cp, double k,
                                     double D, double h, double kc,
                                     double g, double beta, double delta_T,
                                     double t, double k_solid, double L_c,
                                     double f_F,
                                     DimensionlessGroupSet *set);

/* ─── L2: Regime Classification ────────────────────────────────── */

/**
 * Classifies flow regime based on Reynolds number.
 * Returns: 0 = laminar, 1 = transitional, 2 = turbulent
 * Thresholds based on pipe flow criteria.
 */
int flow_regime_pipe(double Re);

/**
 * Classifies flow regime for flat plate flow.
 * Returns: 0 = laminar (Re_x < 5×10⁵),
 *          1 = transitional (5×10⁵ ≤ Re_x ≤ 3×10⁶),
 *          2 = turbulent (Re_x > 3×10⁶)
 */
int flow_regime_flat_plate(double Re_x);

/**
 * Classifies natural convection regime based on Rayleigh number.
 * Returns: 0 = conduction dominant (Ra < 10³),
 *          1 = laminar natural convection (10³ ≤ Ra < 10⁹),
 *          2 = turbulent natural convection (Ra ≥ 10⁹)
 */
int natural_convection_regime(double Ra);

/**
 * Determines the dominant transport mechanism.
 * Based on Péclet number:
 *   Pe << 1: diffusion dominant
 *   Pe ≈ 1: mixed
 *   Pe >> 1: advection dominant
 */
int transport_dominance(double Pe, char *description, size_t desc_size);

/**
 * Checks if lumped capacitance approximation is valid.
 * Bi < 0.1 → valid (temperature uniform within solid).
 * Returns 1 if valid, 0 otherwise.
 */
int lumped_capacitance_valid(double Bi);

/**
 * Returns typical Prandtl number for common engineering fluids.
 * @param fluid_id  0=air, 1=water, 2=engine oil, 3=mercury, 4=CO2, 5=helium
 * @param T         Temperature [K]
 * @return Approximate Prandtl number
 */
double typical_prandtl(int fluid_id, double T);

/**
 * Returns typical Schmidt number for common gas-in-air systems.
 * @param system_id  0=H2O-air, 1=CO2-air, 2=O2-air, 3=NH3-air, 4=ethanol-air
 * @param T          Temperature [K]
 * @return Approximate Schmidt number
 */
double typical_schmidt(int system_id, double T);

/**
 * Provides engineering intuition about what a given Pr number means
 * for boundary layer relative thicknesses.
 */
void prandtl_physical_meaning(double Pr, char *buf, size_t buf_size);

/**
 * Prints a comprehensive table of dimensionless groups
 * with physical interpretation.
 */
void print_dimensionless_summary(const DimensionlessGroupSet *set);

#endif /* DIMENSIONLESS_GROUPS_H */

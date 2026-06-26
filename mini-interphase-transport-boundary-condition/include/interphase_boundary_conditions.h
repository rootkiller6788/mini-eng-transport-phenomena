/**
 * interphase_boundary_conditions.h
 * =================================
 * Boundary condition functions for interphase transport problems.
 *
 * Reference: Bird, Stewart, Lightfoot (2007) "Transport Phenomena"
 *   Chapter 2 - Shell Momentum Balances and Velocity Distributions in Laminar Flow
 *   Chapter 10 - Shell Energy Balances and Temperature Distributions
 *   Chapter 19 - Equations of Change for Multicomponent Systems
 *   Chapter 22 - Interphase Transport in Isothermal Systems
 *
 * Slattery, J.C. (1990) "Interfacial Transport Phenomena"
 *   Chapter 1 - Kinematics and the Jump Entropy Inequality
 *   Chapter 2 - Surface Tension and Curvature Stress
 *
 * Knowledge coverage: L1 Definitions, L2 Core Concepts, L4 Conservation Laws
 */

#ifndef INTERPHASE_BOUNDARY_CONDITIONS_H
#define INTERPHASE_BOUNDARY_CONDITIONS_H

#include "interphase_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * L1: Velocity Boundary Conditions
 * ============================================================================ */

/**
 * Initialize a no-slip velocity boundary condition.
 * The no-slip condition (Stokes, 1851) states that fluid velocity at a solid
 * wall equals the wall velocity: u_fluid = u_wall at y=0.
 *
 * @param bc     Output boundary condition struct
 * @param u_wall Wall velocity components [3] (m/s), NULL for stationary wall
 * Complexity: O(1) constant time
 */
void velocity_bc_no_slip(VelocityBC *bc, const double u_wall[3]);

/**
 * Initialize a Navier slip velocity boundary condition.
 * u_slip = L_s * (du/dn) at wall, where L_s is the slip length (Navier, 1823).
 * Recovered by molecular dynamics for nanochannels (Thompson & Troian, 1997).
 *
 * @param bc           Output boundary condition struct
 * @param slip_length  Navier slip length L_s (m), must be >= 0
 * @param u_wall       Wall velocity [3] (m/s) or NULL
 * Complexity: O(1)
 *
 * Theorem: As L_s -> 0, recovers no-slip condition.
 * Theorem: As L_s -> infinity, approaches perfect slip.
 */
void velocity_bc_navier_slip(VelocityBC *bc, double slip_length,
                              const double u_wall[3]);

/**
 * Initialize a specified shear stress boundary condition (Neumann type).
 * mu * (du/dn) = tau_wall at the boundary.
 *
 * @param bc          Output boundary condition struct
 * @param shear_stress Specified wall shear stress [3] (Pa)
 * Complexity: O(1)
 */
void velocity_bc_shear_stress(VelocityBC *bc, const double shear_stress[3]);

/**
 * Initialize a transpiration (blowing/suction) velocity boundary condition.
 * u_normal = v_w (specified normal velocity), u_tangential = u_wall.
 * Relevant for film cooling, ablation, and transpiration cooling.
 *
 * @param bc            Output boundary condition struct
 * @param transpiration Normal velocity magnitude (m/s), + for blowing, - for suction
 * @param u_wall        Tangential wall velocity [3] (m/s) or NULL
 * Complexity: O(1)
 */
void velocity_bc_transpiration(VelocityBC *bc, double transpiration,
                                const double u_wall[3]);

/**
 * Initialize a symmetry/free-slip boundary condition.
 * du_t/dn = 0 (zero tangential stress), u_n = 0 (no penetration).
 * Used at symmetry planes and inviscid free surfaces.
 *
 * @param bc Output boundary condition struct
 * Complexity: O(1)
 */
void velocity_bc_symmetry(VelocityBC *bc);

/**
 * Initialize an outflow velocity boundary condition.
 * Zero streamwise gradient: d^2u/dn^2 = 0 (fully developed flow assumption).
 *
 * @param bc Output boundary condition struct
 * Complexity: O(1)
 */
void velocity_bc_outflow(VelocityBC *bc);

/* ============================================================================
 * L1: Thermal Boundary Conditions
 * ============================================================================ */

/**
 * Initialize an isothermal (Dirichlet) thermal boundary condition.
 * T = T_wall at the boundary.
 *
 * @param bc     Output boundary condition struct
 * @param T_wall Specified wall temperature (K), must be > 0
 * Complexity: O(1)
 */
void thermal_bc_isothermal(ThermalBC *bc, double T_wall);

/**
 * Initialize a specified heat flux (Neumann) thermal boundary condition.
 * -k * (dT/dn) = q_wall at the boundary. q > 0 means heat into domain.
 *
 * @param bc      Output boundary condition struct
 * @param q_flux  Specified heat flux (W/m^2), + into domain, - out of domain
 * Complexity: O(1)
 */
void thermal_bc_heat_flux(ThermalBC *bc, double q_flux);

/**
 * Initialize an adiabatic (insulated) thermal boundary condition.
 * -k * (dT/dn) = 0 at the boundary. Special case of Neumann with q=0.
 *
 * @param bc Output boundary condition struct
 * Complexity: O(1)
 */
void thermal_bc_adiabatic(ThermalBC *bc);

/**
 * Initialize a convective (Robin) thermal boundary condition.
 * -k * (dT/dn) = h_ext * (T - T_ext) at the boundary.
 * Newton's law of cooling (Newton, 1701).
 *
 * @param bc      Output boundary condition struct
 * @param h_ext   External heat transfer coefficient (W/(m^2*K)), must be >= 0
 * @param T_ext   External fluid temperature (K), must be > 0
 * Complexity: O(1)
 */
void thermal_bc_convective(ThermalBC *bc, double h_ext, double T_ext);

/**
 * Initialize a radiation thermal boundary condition.
 * -k * (dT/dn) = eps * sigma * (T^4 - T_surr^4) at the boundary.
 * Stefan-Boltzmann radiation law. Linearized for small dT.
 *
 * @param bc            Output boundary condition struct
 * @param emissivity    Surface emissivity [0, 1], dimensionless
 * @param T_surroundings Surroundings temperature (K), must be > 0
 * Complexity: O(1)
 */
void thermal_bc_radiation(ThermalBC *bc, double emissivity,
                           double T_surroundings);

/**
 * Initialize a combined convection + radiation thermal BC.
 * -k * (dT/dn) = h*(T-T_ext) + eps*sigma*(T^4-T_surr^4).
 * Used for building heat loss, spacecraft thermal control.
 *
 * @param bc            Output boundary condition struct
 * @param h_ext         Convective heat transfer coefficient (W/(m^2*K))
 * @param T_ext         External fluid temperature (K)
 * @param emissivity    Surface emissivity [0, 1]
 * @param T_surroundings Radiation surroundings temperature (K)
 * Complexity: O(1)
 */
void thermal_bc_combined(ThermalBC *bc, double h_ext, double T_ext,
                          double emissivity, double T_surroundings);

/**
 * Validate the thermal boundary condition for physical consistency.
 * Returns 0 if valid, negative error code otherwise.
 * Error codes: -1 T<0, -2 h<0, -3 emissivity out of [0,1]
 *
 * @param bc Boundary condition to validate
 * @return 0 on valid, negative error code on invalid
 * Complexity: O(1)
 */
int thermal_bc_validate(const ThermalBC *bc);

/* ============================================================================
 * L1: Concentration Boundary Conditions
 * ============================================================================ */

/**
 * Initialize a fixed concentration (Dirichlet) BC.
 * C = C_wall at the boundary.
 *
 * @param bc     Output boundary condition struct
 * @param C_wall Surface concentration (mol/m^3), must be >= 0
 * Complexity: O(1)
 */
void concentration_bc_fixed(ConcentrationBC *bc, double C_wall);

/**
 * Initialize an impermeable wall boundary condition.
 * -D * (dC/dn) = 0 (zero mass flux at wall).
 *
 * @param bc Output boundary condition struct
 * Complexity: O(1)
 */
void concentration_bc_impermeable(ConcentrationBC *bc);

/**
 * Initialize a specified mass flux (Neumann) concentration BC.
 * -D * (dC/dn) = N_wall at the boundary.
 *
 * @param bc       Output boundary condition struct
 * @param mass_flux Specified molar flux (mol/(m^2*s)), + into domain
 * Complexity: O(1)
 */
void concentration_bc_mass_flux(ConcentrationBC *bc, double mass_flux);

/**
 * Initialize a convective mass transfer (Robin) concentration BC.
 * -D * (dC/dn) = k_ext * (C - C_ext).
 *
 * @param bc    Output boundary condition struct
 * @param k_ext External mass transfer coefficient (m/s), must be >= 0
 * @param C_ext External bulk concentration (mol/m^3)
 * Complexity: O(1)
 */
void concentration_bc_convective(ConcentrationBC *bc, double k_ext,
                                  double C_ext);

/**
 * Initialize a surface reaction concentration BC.
 * -D * (dC/dn) = k_rxn * C (first-order surface reaction).
 * C represents reactant concentration; consumption rate = k_rxn * C_surface.
 *
 * @param bc       Output boundary condition struct
 * @param k_rxn    First-order surface reaction rate constant (m/s), >= 0
 * Complexity: O(1)
 */
void concentration_bc_surface_reaction(ConcentrationBC *bc, double k_rxn);

/**
 * Initialize a Henry's law concentration BC for gas-liquid interfaces.
 * C_liquid = p_gas / H, where H is Henry's coefficient (Pa).
 * At interface: p_i = H * x_i (dilute solution assumption).
 *
 * @param bc           Output boundary condition struct
 * @param henry_const  Henry's law constant H (Pa), must be > 0
 * @param p_gas        Partial pressure of species in gas phase (Pa)
 * Complexity: O(1)
 */
void concentration_bc_henry_law(ConcentrationBC *bc, double henry_const,
                                 double p_gas);

/**
 * Initialize a partition-coefficient concentration BC for liquid-liquid interfaces.
 * C_A = K_partition * C_B at equilibrium.
 *
 * @param bc              Output boundary condition struct
 * @param partition_coeff Partition coefficient K = C_A/C_B (-), must be > 0
 * Complexity: O(1)
 */
void concentration_bc_partition(ConcentrationBC *bc, double partition_coeff);

/**
 * Validate the concentration boundary condition.
 * Returns 0 if valid, negative error code otherwise.
 *
 * @param bc Boundary condition to validate
 * @return 0 on valid, negative on invalid
 * Complexity: O(1)
 */
int concentration_bc_validate(const ConcentrationBC *bc);

/* ============================================================================
 * L2: Fluid-Fluid Interface Conditions
 * ============================================================================ */

/**
 * Compute continuity of velocity at a fluid-fluid interface.
 * u_A = u_B at the interface (no slip between phases).
 * This is the kinematic condition for immiscible fluids.
 *
 * @param u_a      Velocity at interface in phase A [3] (m/s)
 * @param u_b      Velocity at interface in phase B [3] (m/s)
 * @param residual Output L2 norm of velocity difference ||u_a - u_b|| (m/s)
 * @return 1 if continuity satisfied within tolerance (1e-6), 0 otherwise
 * Complexity: O(1)
 */
int interface_velocity_continuity(const double u_a[3], const double u_b[3],
                                   double *residual);

/**
 * Compute continuity of shear stress at a fluid-fluid interface.
 * tau_A * n = tau_B * n + sigma_surface (if surface tension gradients exist).
 * For clean interface: mu_A * du_A/dn = mu_B * du_B/dn.
 *
 * @param mu_a         Dynamic viscosity phase A (Pa*s)
 * @param mu_b         Dynamic viscosity phase B (Pa*s)
 * @param du_dn_a      Velocity gradient normal to interface in A [3] (1/s)
 * @param du_dn_b      Velocity gradient normal to interface in B [3] (1/s)
 * @param surface_tension_grad Tangential gradient of surface tension (N/m^2)
 * @param stress_residual Output L2 norm of stress difference
 * @return 1 if continuity satisfied, 0 otherwise
 * Complexity: O(1)
 */
int interface_shear_stress_continuity(double mu_a, double mu_b,
                                       const double du_dn_a[3],
                                       const double du_dn_b[3],
                                       double surface_tension_grad,
                                       double *stress_residual);

/**
 * Compute continuity of temperature at a fluid-fluid interface.
 * T_A = T_B at the interface (assuming no interfacial resistance).
 *
 * @param T_a      Temperature phase A (K), must be > 0
 * @param T_b      Temperature phase B (K), must be > 0
 * @param delta_T  Output |T_a - T_b| (K)
 * @return 1 if continuity satisfied within 1e-3 K, 0 otherwise
 * Complexity: O(1)
 */
int interface_temperature_continuity(double T_a, double T_b, double *delta_T);

/**
 * Compute continuity of heat flux at a fluid-fluid interface.
 * k_A * dT_A/dn = k_B * dT_B/dn (assuming no phase change).
 * With phase change: add latent heat term m_dot * h_fg.
 *
 * @param k_a       Thermal conductivity phase A (W/(m*K))
 * @param k_b       Thermal conductivity phase B (W/(m*K))
 * @param dT_dn_a   Temperature gradient in A at interface (K/m)
 * @param dT_dn_b   Temperature gradient in B at interface (K/m)
 * @param m_dot     Mass flux due to phase change (kg/(m^2*s)), 0 if no phase change
 * @param h_fg      Latent heat of phase change (J/kg), 0 if no phase change
 * @param flux_diff Output heat flux difference (W/m^2)
 * @return 1 if flux balance satisfied within tolerance, 0 otherwise
 * Complexity: O(1)
 */
int interface_heat_flux_continuity(double k_a, double k_b,
                                    double dT_dn_a, double dT_dn_b,
                                    double m_dot, double h_fg,
                                    double *flux_diff);

/**
 * Compute Henry's law equilibrium at a gas-liquid interface.
 * p_i = H * x_i (dilute solution).
 *
 * @param henry_const Henry's law constant H (Pa), must be > 0
 * @param x_liquid    Mole fraction in liquid phase (-), [0,1]
 * @param p_interface Output equilibrium partial pressure (Pa)
 * Complexity: O(1)
 */
void interface_henry_equilibrium(double henry_const, double x_liquid,
                                  double *p_interface);

/**
 * Compute Raoult's law equilibrium for vapor-liquid interfaces.
 * p_i = x_i * gamma_i * p_i^sat (non-ideal liquid mixture).
 *
 * @param x_liquid        Mole fraction in liquid phase (-)
 * @param activity_coeff  Activity coefficient gamma_i (-)
 * @param p_sat           Pure component saturation pressure (Pa)
 * @param p_interface     Output equilibrium partial pressure (Pa)
 * Complexity: O(1)
 */
void interface_raoult_equilibrium(double x_liquid, double activity_coeff,
                                   double p_sat, double *p_interface);

/**
 * Compute concentration jump at interface due to partition coefficient.
 * K = C_A / C_B at equilibrium.
 *
 * @param K          Partition coefficient (-), must be > 0
 * @param C_phase_B  Concentration in phase B (mol/m^3)
 * @param C_phase_A  Output equilibrium concentration in phase A (mol/m^3)
 * Complexity: O(1)
 */
void interface_partition_equilibrium(double K, double C_phase_B,
                                      double *C_phase_A);

/* ============================================================================
 * L4: Conservation-Based Matching Conditions
 * ============================================================================ */

/**
 * Apply the overall mass continuity condition across interface.
 * For steady state with no phase change: rho_A * u_nA = rho_B * u_nB.
 * With phase change: rho_A*(u_nA - u_nI) = rho_B*(u_nB - u_nI).
 *
 * @param rho_A     Density phase A (kg/m^3)
 * @param rho_B     Density phase B (kg/m^3)
 * @param u_n_A     Normal velocity in phase A (m/s)
 * @param u_n_B     Normal velocity in phase B (m/s)
 * @param u_n_I     Interface normal velocity (m/s), 0 for steady interface
 * @param mass_flux Output mass flux across interface (kg/(m^2*s))
 * @return 1 if continuity satisfied, 0 if violation detected
 * Complexity: O(1)
 */
int interface_mass_continuity(double rho_A, double rho_B,
                               double u_n_A, double u_n_B, double u_n_I,
                               double *mass_flux);

/**
 * Compute the normal stress jump at a curved fluid-fluid interface
 * (Young-Laplace equation).
 * p_A - p_B = sigma * (1/R1 + 1/R2) = 2 * sigma * H
 * where sigma is surface tension, H is mean curvature.
 *
 * @param sigma            Surface tension (N/m)
 * @param curvature_mean   Mean curvature H = 0.5*(1/R1 + 1/R2) (1/m)
 * @param pressure_jump    Output p_A - p_B (Pa), > 0 means A at higher pressure
 * Complexity: O(1)
 *
 * Theorem: For a spherical droplet of radius R, p_in - p_out = 2*sigma/R.
 * Theorem: For a planar interface (H=0), p_A = p_B.
 */
void interface_young_laplace(double sigma, double curvature_mean,
                              double *pressure_jump);

/**
 * Compute the Marangoni stress due to surface tension gradient.
 * tau_surface = d_sigma/dT * grad_T (tangential to interface).
 * Drives Benard-Marangoni convection cells (Benard, 1900; Pearson, 1958).
 *
 * @param d_sigma_dT   Temperature coefficient of surface tension (N/(m*K))
 * @param dT_ds        Tangential temperature gradient (K/m)
 * @param marangoni_stress Output tangential stress (Pa), direction of grad_T
 * Complexity: O(1)
 */
void interface_marangoni_stress(double d_sigma_dT, double dT_ds,
                                 double *marangoni_stress);

#ifdef __cplusplus
}
#endif

#endif /* INTERPHASE_BOUNDARY_CONDITIONS_H */

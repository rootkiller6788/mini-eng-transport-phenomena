/**
 * interphase_jump_conditions.h
 * =============================
 * Jump conditions at phase interfaces derived from integral conservation laws.
 *
 * Reference: Slattery, J.C. (1990) "Interfacial Transport Phenomena"
 *   Chapter 1 - Kinematics, Mass, Momentum, Energy, and Entropy Balances
 *   Appendix B - Jump Balances Derivation from Integral Conservation Principles
 *
 * Bird, Stewart, Lightfoot (2007) "Transport Phenomena"
 *   Appendix A - Vector and Tensor Notation
 *   Section 22.5 - The Jump Energy Balance
 *
 * Delhaye, J.M. (1974) "Jump Conditions and Entropy Sources in Two-Phase Systems"
 *   Int. J. Multiphase Flow, 1(3), 395-409.
 *
 * Knowledge coverage: L4 Conservation Laws, L1 Definitions, L5 Engineering Methods
 */

#ifndef INTERPHASE_JUMP_CONDITIONS_H
#define INTERPHASE_JUMP_CONDITIONS_H

#include "interphase_types.h"
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * L4: Conservative Jump Balances (Slattery, 1990, Appendix B)
 *
 * General form: [[rho * psi * (v - v_I) * n + phi * n]] = surface_source
 * where [[X]] = X_B - X_A, psi is conserved quantity per unit mass,
 * phi is surface flux vector, n is unit normal from A to B.
 * ============================================================================ */

/**
 * Evaluate the mass jump condition at a phase interface.
 * [[rho * (v - v_I) * n]] = 0  (no phase change).
 * [[rho * (v - v_I) * n]] = 0 still holds with phase change,
 * because mass is conserved in the phase transformation.
 *
 * With phase change: rho_A*(v_A - v_I)*n = rho_B*(v_B - v_I)*n = m_dot.
 *
 * @param rho_A   Density in phase A (kg/m^3), must be > 0
 * @param rho_B   Density in phase B (kg/m^3), must be > 0
 * @param v_A     Velocity in phase A [3] (m/s)
 * @param v_B     Velocity in phase B [3] (m/s)
 * @param v_I     Interface velocity [3] (m/s)
 * @param n       Unit normal from A to B [3] (dimensionless)
 * @param m_dot   Output mass flux across interface (kg/(m^2*s))
 *                m_dot > 0 for phase A to phase B transport
 * @param jump    Output computed mass jump residual (kg/(m^2*s))
 *                Should be ~0 for valid input
 * Complexity: O(1)
 *
 * Theorem (Mass conservation): m_dot = rho_A*(v_A - v_I)*n = rho_B*(v_B - v_I)*n.
 * Theorem: If no phase change, v_A*n = v_B*n = v_I*n and m_dot = 0.
 * Theorem: For condensation, m_dot > 0 (vapor -> liquid).
 * Theorem: For evaporation, m_dot < 0 (liquid -> vapor).
 */
void jump_mass(double rho_A, double rho_B,
               const double v_A[3], const double v_B[3],
               const double v_I[3], const double n[3],
               double *m_dot, double *jump);

/**
 * Evaluate the linear momentum jump condition at a phase interface.
 * [[rho*v*(v - v_I)*n + p*I*n - tau*n]] = 2*H*sigma*n + grad_s(sigma)
 *
 * Where:
 *   [[X]] = X_B - X_A is the jump across interface
 *   2*H*sigma*n = capillary pressure (normal component)
 *   grad_s(sigma) = Marangoni stress (tangential component)
 *   tau = mu*(grad(v) + grad(v)^T) is the viscous stress tensor
 *
 * @param rho_A         Density phase A (kg/m^3)
 * @param rho_B         Density phase B (kg/m^3)
 * @param p_A           Pressure phase A (Pa)
 * @param p_B           Pressure phase B (Pa)
 * @param v_A           Velocity phase A [3] (m/s)
 * @param v_B           Velocity phase B [3] (m/s)
 * @param v_I           Interface velocity [3] (m/s)
 * @param n             Unit normal A->B [3]
 * @param tau_A         Viscous stress phase A as 3x3 matrix (row-major) [Pa]
 * @param tau_B         Viscous stress phase B as 3x3 matrix (row-major) [Pa]
 * @param sigma         Surface tension (N/m)
 * @param curvature_mean Mean curvature H = 0.5*(1/R1+1/R2) (1/m)
 * @param grad_s_sigma  Surface gradient of surface tension [2] (N/m^2)
 * @param t1            First tangent direction [3]
 * @param t2            Second tangent direction [3]
 * @param jump_n        Output normal momentum jump residual (Pa)
 * @param jump_t1       Output tangential momentum jump in t1 dir (Pa)
 * @param jump_t2       Output tangential momentum jump in t2 dir (Pa)
 * Complexity: O(1)
 *
 * Theorem (Normal jump): p_A - p_B = 2*sigma*H + (tau_nn_A - tau_nn_B).
 *   Reduces to Young-Laplace when viscous normal stresses negligible.
 * Theorem (Tangential jump): tau_nt_A - tau_nt_B = grad_s(sigma).
 *   Drives Marangoni flows when surface tension varies spatially.
 */
void jump_momentum(double rho_A, double rho_B,
                   double p_A, double p_B,
                   const double v_A[3], const double v_B[3],
                   const double v_I[3], const double n[3],
                   const double tau_A[9], const double tau_B[9],
                   double sigma, double curvature_mean,
                   const double grad_s_sigma[2],
                   const double t1[3], const double t2[3],
                   double *jump_n, double *jump_t1, double *jump_t2);

/**
 * Evaluate the energy jump condition at a phase interface.
 * [[rho*(U + 0.5*v^2)*(v - v_I)*n + q*n - tau*v*n]] = surface_energy_source
 *
 * Where:
 *   U = specific internal energy (J/kg)
 *   0.5*v^2 = specific kinetic energy (J/kg)
 *   q*n = conductive heat flux (W/m^2)
 *   tau*v*n = work done by viscous stresses
 *
 * For phase change (e.g., boiling/condensation):
 *   k_A*dT_A/dn - k_B*dT_B/dn = m_dot * h_fg
 *   (latent heat released/absorbed at interface)
 *
 * @param k_A         Thermal conductivity phase A (W/(m*K))
 * @param k_B         Thermal conductivity phase B (W/(m*K))
 * @param dT_dn_A     Normal temperature gradient phase A (K/m)
 * @param dT_dn_B     Normal temperature gradient phase B (K/m)
 * @param m_dot       Mass flux across interface (kg/(m^2*s))
 * @param h_fg        Latent heat of phase change (J/kg), > 0 for evaporation
 * @param jump        Output energy jump residual (W/m^2)
 *                    Should be ~0 for valid inputs satisfying conservation
 * Complexity: O(1)
 *
 * Theorem: k_A*dT/dn_A - k_B*dT/dn_B = m_dot*h_fg  (Stefan condition).
 * Theorem: For m_dot = 0, heat flux is continuous: k_A*dT/dn_A = k_B*dT/dn_B.
 * Theorem: At boiling, m_dot > 0, so net heat flux goes into phase change.
 * Theorem: At condensation, m_dot < 0, net heat flux released from phase change.
 */
void jump_energy(double k_A, double k_B,
                 double dT_dn_A, double dT_dn_B,
                 double m_dot, double h_fg,
                 double *jump);

/**
 * Evaluate the species jump condition at a phase interface.
 * [[rho*omega_i*(v - v_I)*n + j_i*n]] = r_i^surface
 *
 * Where:
 *   omega_i = mass fraction of species i
 *   j_i = diffusive mass flux of species i relative to mass-average velocity
 *         j_i = -rho*D_i*grad(omega_i) (Fick's law)
 *   r_i^surface = surface reaction rate of species i (kg/(m^2*s))
 *
 * @param rho_A           Density phase A (kg/m^3)
 * @param rho_B           Density phase B (kg/m^3)
 * @param omega_i_A       Mass fraction of species i in phase A
 * @param omega_i_B       Mass fraction of species i in phase B
 * @param normal_vel_A    Normal mass-average velocity in phase A (m/s)
 * @param normal_vel_B    Normal mass-average velocity in phase B (m/s)
 * @param v_I_n           Interface normal velocity (m/s)
 * @param D_i_A           Diffusivity of species i in phase A (m^2/s)
 * @param D_i_B           Diffusivity of species i in phase B (m^2/s)
 * @param domega_dn_A     Normal gradient of omega_i in phase A (1/m)
 * @param domega_dn_B     Normal gradient of omega_i in phase B (1/m)
 * @param surface_rxn     Surface reaction rate (kg/(m^2*s)), + for production at interface
 * @param jump            Output species jump residual (kg/(m^2*s))
 * Complexity: O(1)
 *
 * Theorem: For non-reacting interface, species flux is continuous.
 * Theorem: At catalytic surface, diffusive flux = surface reaction rate.
 * Theorem: At phase interface with equilibrium, the jump is determined by
 *          partition coefficient or Henry's law rather than flux continuity.
 */
void jump_species(double rho_A, double rho_B,
                  double omega_i_A, double omega_i_B,
                  double normal_vel_A, double normal_vel_B,
                  double v_I_n,
                  double D_i_A, double D_i_B,
                  double domega_dn_A, double domega_dn_B,
                  double surface_rxn, double *jump);

/**
 * Evaluate the entropy jump condition (second law at interface).
 * [[rho*s*(v - v_I)*n + q*n/T]] >= 0  (entropy production at interface).
 *
 * The entropy jump must be non-negative for a physically admissible
 * process (local formulation of the second law at the interface).
 *
 * @param rho_A       Density phase A (kg/m^3)
 * @param rho_B       Density phase B (kg/m^3)
 * @param s_A         Specific entropy phase A (J/(kg*K))
 * @param s_B         Specific entropy phase B (J/(kg*K))
 * @param normal_vel_A Normal velocity phase A relative to interface (m/s)
 * @param normal_vel_B Normal velocity phase B relative to interface (m/s)
 * @param q_n_A       Normal heat flux in phase A at interface (W/m^2)
 * @param q_n_B       Normal heat flux in phase B at interface (W/m^2)
 * @param T_A         Temperature at interface phase A side (K)
 * @param T_B         Temperature at interface phase B side (K)
 * @param entropy_production Output entropy production rate at interface (W/(m^2*K))
 *                           Must be >= 0 for physically valid process
 * @return 1 if second law satisfied (entropy_production >= 0), 0 if violated
 *
 * Theorem (Second Law at Interface): Interface entropy production >= 0.
 * Theorem: For T_A = T_B and no phase change, entropy production from heat
 *          flux discontinuity: q_n_A*(1/T_A - 1/T_B). With T_A = T_B, this
 *          vanishes, so q continuity follows from entropy considerations
 *          only if no other irreversibilities exist.
 */
int jump_entropy(double rho_A, double rho_B,
                 double s_A, double s_B,
                 double normal_vel_A, double normal_vel_B,
                 double q_n_A, double q_n_B,
                 double T_A, double T_B,
                 double *entropy_production);

/**
 * Evaluate the complete set of jump conditions for a phase interface
 * and return a structured report. This is the master function for
 * interface conservation verification.
 *
 * @param iface   Interface descriptor with all state data
 * @param sigma   Surface tension (N/m)
 * @param H       Mean curvature (1/m)
 * @param jumps   Output jump condition structure
 * Complexity: O(1)
 *
 * Theorem (Consistency): All five jump balances must be simultaneously
 * satisfied for a physically realizable interface configuration.
 */
void jump_all_conditions(const InterfaceDescriptor *iface,
                          double sigma, double H,
                          JumpConditions *jumps);

/* ============================================================================
 * L4: Simplified Jump Conditions for Common Cases
 * ============================================================================ */

/**
 * Stefan condition for phase change problems.
 * k_s * dT_s/dn - k_l * dT_l/dn = rho * L * v_n
 * where v_n is interface normal velocity, L is latent heat.
 * This is a simplified energy jump for pure substance phase change.
 *
 * @param k_s     Solid phase thermal conductivity (W/(m*K))
 * @param k_l     Liquid phase thermal conductivity (W/(m*K))
 * @param dT_dn_s Temperature gradient in solid at interface (K/m)
 * @param dT_dn_l Temperature gradient in liquid at interface (K/m)
 * @param rho_s   Solid density (kg/m^3)
 * @param L       Latent heat of fusion (J/kg)
 * @param v_n     Output interface normal velocity (m/s), + into liquid
 * Complexity: O(1)
 *
 * Theorem: v_n > 0 for melting (heat flows into interface).
 * Theorem: v_n < 0 for solidification (heat flows out of interface).
 * Theorem: The Stefan number Ste = cp*delta_T/L controls interface velocity.
 */
void stefan_condition(double k_s, double k_l,
                      double dT_dn_s, double dT_dn_l,
                      double rho_s, double L, double *v_n);

/**
 * Hertz-Knudsen-Schrage condition for evaporation/condensation kinetics.
 * m_dot = (2*alpha/(2-alpha)) * sqrt(M/(2*pi*R)) * (p_sat/T_l^(1/2) - p_v/T_v^(1/2)).
 *
 * @param alpha        Accommodation coefficient [0, 1]
 * @param M            Molar mass (kg/mol)
 * @param R            Universal gas constant 8.314 (J/(mol*K))
 * @param p_sat        Saturation pressure at liquid temperature (Pa)
 * @param T_l          Liquid interface temperature (K)
 * @param p_v          Vapor partial pressure near interface (Pa)
 * @param T_v          Vapor temperature near interface (K)
 * @param m_dot        Output mass flux (kg/(m^2*s)), + for net evaporation
 * Complexity: O(1)
 *
 * Theorem: At equilibrium (p_v = p_sat, T_l = T_v), m_dot = 0.
 * Theorem: alpha = 1 corresponds to theoretical maximum evaporation rate.
 * Theorem: For alpha close to 1, process is heat-transfer limited.
 * Theorem: For small alpha, kinetic resistance dominates.
 */
void hertz_knudsen_schrage(double alpha, double M, double R,
                            double p_sat, double T_l,
                            double p_v, double T_v,
                            double *m_dot);

#ifdef __cplusplus
}
#endif

#endif /* INTERPHASE_JUMP_CONDITIONS_H */

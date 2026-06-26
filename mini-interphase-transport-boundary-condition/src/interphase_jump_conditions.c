/**
 * interphase_jump_conditions.c
 * =============================
 * Implementation of jump conditions at phase interfaces.
 *
 * Reference: Slattery, J.C. (1990) "Interfacial Transport Phenomena"
 *   Appendix B - Derivation of Jump Balances from Integral Balances
 *
 * Delhaye, J.M. (1974) "Jump Conditions and Entropy Sources in
 *   Two-Phase Systems" Int. J. Multiphase Flow, 1(3), 395-409.
 *
 * Bird, Stewart, Lightfoot (2007) "Transport Phenomena"
 *   Section 22.5 - The Jump Energy Balance
 *
 * Knowledge coverage: L4 Conservation Laws, L5 Engineering Methods
 */

#include "interphase_jump_conditions.h"
#include <math.h>
#include <string.h>

/* ============================================================================
 * L4: Mass Jump Condition
 *
 * Integral mass balance for control volume containing interface:
 *   d/dt integral_CV rho dV + integral_CS rho (v - v_CS) * n dA = 0
 *
 * As CV shrinks to interface: [[rho (v - v_I) * n]] = 0
 * ============================================================================ */

void jump_mass(double rho_A, double rho_B,
               const double v_A[3], const double v_B[3],
               const double v_I[3], const double n[3],
               double *m_dot, double *jump) {
    if (!v_A || !v_B || !v_I || !n) {
        if (m_dot) *m_dot = 0.0;
        if (jump)  *jump  = 0.0;
        return;
    }

    /* Normal velocity in phase A relative to interface */
    double u_rel_A = (v_A[0] - v_I[0]) * n[0]
                   + (v_A[1] - v_I[1]) * n[1]
                   + (v_A[2] - v_I[2]) * n[2];

    /* Normal velocity in phase B relative to interface */
    double u_rel_B = (v_B[0] - v_I[0]) * n[0]
                   + (v_B[1] - v_I[1]) * n[1]
                   + (v_B[2] - v_I[2]) * n[2];

    double m_dot_A = rho_A * u_rel_A;
    double m_dot_B = rho_B * u_rel_B;

    if (m_dot) *m_dot = m_dot_A;
    if (jump)  *jump  = m_dot_A - m_dot_B;
}

/* ============================================================================
 * L4: Linear Momentum Jump Condition
 *
 * [[rho v (v - v_I) * n + p I * n - tau * n]] = 2 H sigma n + grad_s(sigma)
 *
 * This is the vector jump condition. We decompose into normal and
 * tangential components for engineering use.
 * ============================================================================ */

void jump_momentum(double rho_A, double rho_B,
                   double p_A, double p_B,
                   const double v_A[3], const double v_B[3],
                   const double v_I[3], const double n[3],
                   const double tau_A[9], const double tau_B[9],
                   double sigma, double curvature_mean,
                   const double grad_s_sigma[2],
                   const double t1[3], const double t2[3],
                   double *jump_n, double *jump_t1, double *jump_t2) {
    if (!v_A || !v_B || !v_I || !n || !tau_A || !tau_B) {
        if (jump_n)  *jump_n  = 0.0;
        if (jump_t1) *jump_t1 = 0.0;
        if (jump_t2) *jump_t2 = 0.0;
        return;
    }

    /* Compute relative normal velocities */
    double urel_A = (v_A[0]-v_I[0])*n[0] + (v_A[1]-v_I[1])*n[1]
                    + (v_A[2]-v_I[2])*n[2];
    double urel_B = (v_B[0]-v_I[0])*n[0] + (v_B[1]-v_I[1])*n[1]
                    + (v_B[2]-v_I[2])*n[2];

    /* Convective momentum flux in normal direction: rho * v * (v-v_I)*n */
    /* For phase A: rho_A * v_A * urel_A, add to momentum flux vector */
    double mom_flux_A[3] = {
        rho_A * v_A[0] * urel_A + p_A * n[0],
        rho_A * v_A[1] * urel_A + p_A * n[1],
        rho_A * v_A[2] * urel_A + p_A * n[2]
    };

    double mom_flux_B[3] = {
        rho_B * v_B[0] * urel_B + p_B * n[0],
        rho_B * v_B[1] * urel_B + p_B * n[1],
        rho_B * v_B[2] * urel_B + p_B * n[2]
    };

    /* Subtract viscous stress contribution: tau * n */
    /* tau is row-major: tau[row][col], tau*n gives vector with comp i = sum_j tau_ij * n_j */
    double tau_n_A[3] = {
        tau_A[0]*n[0] + tau_A[1]*n[1] + tau_A[2]*n[2],
        tau_A[3]*n[0] + tau_A[4]*n[1] + tau_A[5]*n[2],
        tau_A[6]*n[0] + tau_A[7]*n[1] + tau_A[8]*n[2]
    };

    double tau_n_B[3] = {
        tau_B[0]*n[0] + tau_B[1]*n[1] + tau_B[2]*n[2],
        tau_B[3]*n[0] + tau_B[4]*n[1] + tau_B[5]*n[2],
        tau_B[6]*n[0] + tau_B[7]*n[1] + tau_B[8]*n[2]
    };

    /* Total momentum flux vector on each side */
    double total_A[3] = {
        mom_flux_A[0] - tau_n_A[0],
        mom_flux_A[1] - tau_n_A[1],
        mom_flux_A[2] - tau_n_A[2]
    };

    double total_B[3] = {
        mom_flux_B[0] - tau_n_B[0],
        mom_flux_B[1] - tau_n_B[1],
        mom_flux_B[2] - tau_n_B[2]
    };

    /* Jump vector: B - A (raw difference) */
    double jump_vec[3] = {
        total_B[0] - total_A[0],
        total_B[1] - total_A[1],
        total_B[2] - total_A[2]
    };

    /* Surface force: 2*H*sigma*n + grad_s(sigma) */
    double surface_force[3] = {
        2.0 * curvature_mean * sigma * n[0],
        2.0 * curvature_mean * sigma * n[1],
        2.0 * curvature_mean * sigma * n[2]
    };

    /* Add Marangoni: grad_s(sigma) = (d_sigma/ds1)*t1 + (d_sigma/ds2)*t2 */
    if (t1 && grad_s_sigma) {
        surface_force[0] += grad_s_sigma[0] * t1[0];
        surface_force[1] += grad_s_sigma[0] * t1[1];
        surface_force[2] += grad_s_sigma[0] * t1[2];
    }
    if (t2 && grad_s_sigma) {
        surface_force[0] += grad_s_sigma[1] * t2[0];
        surface_force[1] += grad_s_sigma[1] * t2[1];
        surface_force[2] += grad_s_sigma[1] * t2[2];
    }

    /* Residual: jump_vec - surface_force should be zero */
    double residual[3] = {
        jump_vec[0] - surface_force[0],
        jump_vec[1] - surface_force[1],
        jump_vec[2] - surface_force[2]
    };

    /* Decompose into normal and tangential components */
    if (jump_n) {
        *jump_n = residual[0]*n[0] + residual[1]*n[1] + residual[2]*n[2];
    }
    if (jump_t1 && t1) {
        *jump_t1 = residual[0]*t1[0] + residual[1]*t1[1] + residual[2]*t1[2];
    }
    if (jump_t2 && t2) {
        *jump_t2 = residual[0]*t2[0] + residual[1]*t2[1] + residual[2]*t2[2];
    }
}

/* ============================================================================
 * L4: Energy Jump Condition
 *
 * [[rho*(U + 0.5*v^2)*(v - v_I)*n + q*n - tau*v*n]] = 0
 *
 * Simplified (no kinetic energy, no viscous work):
 *   k_A * dT/dn_A - k_B * dT/dn_B = m_dot * h_fg  (phase change)
 * ============================================================================ */

void jump_energy(double k_A, double k_B,
                 double dT_dn_A, double dT_dn_B,
                 double m_dot, double h_fg,
                 double *jump) {
    if (!jump) return;

    /* Heat flux on each side: q_n = -k * dT/dn (Fourier's law) */
    double q_n_A = -k_A * dT_dn_A;
    double q_n_B = -k_B * dT_dn_B;

    /* Energy jump includes latent heat of phase change */
    /* [[q*n]] + m_dot * h_fg should = 0 */
    *jump = q_n_B - q_n_A + m_dot * h_fg;

    /* Alternative convention:
     * k_A * dT_A/dn - k_B * dT_B/dn = m_dot * h_fg
     * (net conductive heat flux into interface = latent heat released)
     * jump = k_A*dT_A/dn - k_B*dT_B/dn - m_dot*h_fg should = 0
     */
}

/* ============================================================================
 * L4: Species Jump Condition
 *
 * [[rho*omega_i*(v - v_I)*n + j_i*n]] = r_i^surface
 *
 * For non-reacting interface without phase change:
 *   Flux of i in phase A = Flux of i in phase B
 * ============================================================================ */

void jump_species(double rho_A, double rho_B,
                  double omega_i_A, double omega_i_B,
                  double normal_vel_A, double normal_vel_B,
                  double v_I_n,
                  double D_i_A, double D_i_B,
                  double domega_dn_A, double domega_dn_B,
                  double surface_rxn, double *jump) {
    if (!jump) return;

    /* Relative normal velocities */
    double urel_A = normal_vel_A - v_I_n;
    double urel_B = normal_vel_B - v_I_n;

    /* Convective flux: rho * omega * u_rel */
    double conv_A = rho_A * omega_i_A * urel_A;
    double conv_B = rho_B * omega_i_B * urel_B;

    /* Diffusive flux (Fick's law): -rho * D * d(omega)/dn */
    double diff_A = -rho_A * D_i_A * domega_dn_A;
    double diff_B = -rho_B * D_i_B * domega_dn_B;

    /* Total flux on each side */
    double total_A = conv_A + diff_A;
    double total_B = conv_B + diff_B;

    /* Jump: B - A = surface reaction rate */
    *jump = total_B - total_A - surface_rxn;
}

/* ============================================================================
 * L4: Entropy Jump Condition
 *
 * [[rho*s*(v - v_I)*n + q*n/T]] >= 0
 *
 * Local second law: entropy production at interface must be non-negative.
 * ============================================================================ */

int jump_entropy(double rho_A, double rho_B,
                 double s_A, double s_B,
                 double normal_vel_A, double normal_vel_B,
                 double q_n_A, double q_n_B,
                 double T_A, double T_B,
                 double *entropy_production) {
    /* Entropy flux on side A: rho_A * s_A * u_rel_A + q_n_A / T_A */
    double entropy_flux_A = rho_A * s_A * (normal_vel_A) + q_n_A / T_A;

    /* Entropy flux on side B */
    double entropy_flux_B = rho_B * s_B * (normal_vel_B) + q_n_B / T_B;

    /* Entropy production = flux_B - flux_A (must be >= 0) */
    double sigma_s = entropy_flux_B - entropy_flux_A;

    if (entropy_production) *entropy_production = sigma_s;

    /* Physical admissibility: sigma_s >= -epsilon (numerical tolerance) */
    return (sigma_s >= -1e-10) ? 1 : 0;
}

/* ============================================================================
 * L4: Complete Jump Condition Evaluation
 * ============================================================================ */

void jump_all_conditions(const InterfaceDescriptor *iface,
                          double sigma, double H,
                          JumpConditions *jumps) {
    (void)sigma; (void)H;  /* reserved for surface tension + curvature effects */
    if (!iface || !jumps) return;
    memset(jumps, 0, sizeof(JumpConditions));

    double n[3] = {0.0, 0.0, 1.0};  /* default normal (planar interface) */
    double v_I[3] = {0.0, 0.0, 0.0};
    double v_A[3] = {0.0, 0.0, 0.0};
    double v_B[3] = {0.0, 0.0, 0.0};

    /* Mass jump */
    double m_dot;
    jump_mass(iface->state_a.density, iface->state_b.density,
              v_A, v_B, v_I, n, &m_dot, &jumps->mass_jump);

    /* Energy jump (simplified: no kinetic or viscous terms) */
    jump_energy(iface->state_a.thermal_conductivity,
                iface->state_b.thermal_conductivity,
                0.0, 0.0,  /* zero temperature gradients for now */
                m_dot, 0.0, &jumps->energy_jump);

    /* Species jump (simplified) */
    jumps->species_jump = 0.0;
    jumps->entropy_jump = 0.0;
}

/* ============================================================================
 * L4: Simplified Jump Conditions for Common Cases
 * ============================================================================ */

void stefan_condition(double k_s, double k_l,
                      double dT_dn_s, double dT_dn_l,
                      double rho_s, double L, double *v_n) {
    if (!v_n) return;
    if (rho_s <= 0.0 || L <= 0.0) {
        *v_n = 0.0;
        return;
    }

    /* Stefan condition:
     * k_s * dT_s/dn - k_l * dT_l/dn = rho_s * L * v_n
     * => v_n = (k_s * dT_s/dn - k_l * dT_l/dn) / (rho_s * L)
     */
    *v_n = (k_s * dT_dn_s - k_l * dT_dn_l) / (rho_s * L);
}

void hertz_knudsen_schrage(double alpha, double M, double R,
                            double p_sat, double T_l,
                            double p_v, double T_v,
                            double *m_dot) {
    if (!m_dot) return;
    if (alpha <= 0.0 || alpha > 1.0 || M <= 0.0 || R <= 0.0
        || T_l <= 0.0 || T_v <= 0.0) {
        *m_dot = 0.0;
        return;
    }

    /* Hertz-Knudsen-Schrage equation:
     * m_dot = (2*alpha/(2-alpha)) * sqrt(M/(2*pi*R))
     *         * (p_sat/sqrt(T_l) - p_v/sqrt(T_v))
     */
    double coeff = (2.0 * alpha) / (2.0 - alpha)
                   * sqrt(M / (2.0 * M_PI * R));

    double driving_force = p_sat / sqrt(T_l) - p_v / sqrt(T_v);

    *m_dot = coeff * driving_force;
}

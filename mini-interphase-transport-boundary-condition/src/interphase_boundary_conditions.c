/**
 * interphase_boundary_conditions.c
 * =================================
 * Implementation of boundary conditions at phase interfaces.
 *
 * Reference: Bird, Stewart, Lightfoot (2007) "Transport Phenomena"
 *   Chapter 2 - Velocity Boundary Conditions
 *   Chapter 10 - Thermal Boundary Conditions
 *   Chapter 22 - Interfacial Concentration Conditions
 *
 * Knowledge coverage: L1 Definitions, L2 Core Concepts, L4 Conservation Laws
 */

#include "interphase_boundary_conditions.h"
#include <math.h>
#include <string.h>

/* ============================================================================
 * L1: Velocity Boundary Conditions
 *
 * Reference: Stokes (1851) "On the Effect of the Internal Friction of
 *   Fluids on the Motion of Pendulums" (Cambridge Philosophical
 *   Transactions, Vol. IX) - established the no-slip condition.
 * ============================================================================ */

void velocity_bc_no_slip(VelocityBC *bc, const double u_wall[3]) {
    if (!bc) return;
    memset(bc, 0, sizeof(VelocityBC));
    bc->type = BC_DIRICHLET;
    if (u_wall) {
        bc->u_wall[0] = u_wall[0];
        bc->u_wall[1] = u_wall[1];
        bc->u_wall[2] = u_wall[2];
    } else {
        bc->u_wall[0] = 0.0;
        bc->u_wall[1] = 0.0;
        bc->u_wall[2] = 0.0;
    }
    bc->slip_coefficient = 0.0;
}

void velocity_bc_navier_slip(VelocityBC *bc, double slip_length,
                              const double u_wall[3]) {
    if (!bc || slip_length < 0.0) return;
    memset(bc, 0, sizeof(VelocityBC));
    bc->type = BC_ROBIN;  /* mixed: u - L_s * du/dn = u_wall */
    bc->slip_coefficient = slip_length;
    if (u_wall) {
        bc->u_wall[0] = u_wall[0];
        bc->u_wall[1] = u_wall[1];
        bc->u_wall[2] = u_wall[2];
    }
    /* Robin form: u - L_s * du/dn = u_wall at boundary */
}

void velocity_bc_shear_stress(VelocityBC *bc, const double shear_stress[3]) {
    if (!bc || !shear_stress) return;
    memset(bc, 0, sizeof(VelocityBC));
    bc->type = BC_NEUMANN;
    bc->shear_stress[0] = shear_stress[0];
    bc->shear_stress[1] = shear_stress[1];
    bc->shear_stress[2] = shear_stress[2];
}

void velocity_bc_transpiration(VelocityBC *bc, double transpiration,
                                const double u_wall[3]) {
    if (!bc) return;
    memset(bc, 0, sizeof(VelocityBC));
    bc->type = BC_CAUCHY;  /* both normal and tangential specified */
    bc->transpiration = transpiration;
    if (u_wall) {
        bc->u_wall[0] = u_wall[0];
        bc->u_wall[1] = u_wall[1];
        bc->u_wall[2] = u_wall[2];
    }
}

void velocity_bc_symmetry(VelocityBC *bc) {
    if (!bc) return;
    memset(bc, 0, sizeof(VelocityBC));
    bc->type = BC_SYMMETRY;
    /* Implies: u_n = 0, du_t/dn = 0 */
}

void velocity_bc_outflow(VelocityBC *bc) {
    if (!bc) return;
    memset(bc, 0, sizeof(VelocityBC));
    bc->type = BC_OUTFLOW;
    /* Implies: zero second derivative in streamwise direction */
}

/* ============================================================================
 * L1: Thermal Boundary Conditions
 * ============================================================================ */

void thermal_bc_isothermal(ThermalBC *bc, double T_wall) {
    if (!bc || T_wall <= 0.0) return;
    memset(bc, 0, sizeof(ThermalBC));
    bc->type = BC_DIRICHLET;
    bc->T_wall = T_wall;
}

void thermal_bc_heat_flux(ThermalBC *bc, double q_flux) {
    if (!bc) return;
    memset(bc, 0, sizeof(ThermalBC));
    bc->type = BC_NEUMANN;
    bc->heat_flux = q_flux;
}

void thermal_bc_adiabatic(ThermalBC *bc) {
    if (!bc) return;
    memset(bc, 0, sizeof(ThermalBC));
    bc->type = BC_NEUMANN;
    bc->heat_flux = 0.0;
}

void thermal_bc_convective(ThermalBC *bc, double h_ext, double T_ext) {
    if (!bc || h_ext < 0.0 || T_ext <= 0.0) return;
    memset(bc, 0, sizeof(ThermalBC));
    bc->type = BC_ROBIN;
    bc->h_ext = h_ext;
    bc->T_ext = T_ext;
}

void thermal_bc_radiation(ThermalBC *bc, double emissivity,
                           double T_surroundings) {
    if (!bc || emissivity < 0.0 || emissivity > 1.0
        || T_surroundings <= 0.0) return;
    memset(bc, 0, sizeof(ThermalBC));
    bc->type = BC_ROBIN;  /* nonlinear, but Robin after linearization */
    bc->emissivity = emissivity;
    bc->T_surroundings = T_surroundings;
}

void thermal_bc_combined(ThermalBC *bc, double h_ext, double T_ext,
                          double emissivity, double T_surroundings) {
    if (!bc) return;
    memset(bc, 0, sizeof(ThermalBC));
    bc->type = BC_ROBIN;  /* combined convection + radiation */
    bc->h_ext = h_ext;
    bc->T_ext = T_ext;
    bc->emissivity = emissivity;
    bc->T_surroundings = T_surroundings;
}

int thermal_bc_validate(const ThermalBC *bc) {
    if (!bc) return -4;  /* null pointer */
    if (bc->type == BC_DIRICHLET && bc->T_wall <= 0.0) return -1;
    if (bc->type == BC_ROBIN && bc->h_ext < 0.0) return -2;
    if (bc->emissivity < 0.0 || bc->emissivity > 1.0) return -3;
    return 0;
}

/* ============================================================================
 * L1: Concentration Boundary Conditions
 * ============================================================================ */

void concentration_bc_fixed(ConcentrationBC *bc, double C_wall) {
    if (!bc || C_wall < 0.0) return;
    memset(bc, 0, sizeof(ConcentrationBC));
    bc->type = BC_DIRICHLET;
    bc->C_wall = C_wall;
}

void concentration_bc_impermeable(ConcentrationBC *bc) {
    if (!bc) return;
    memset(bc, 0, sizeof(ConcentrationBC));
    bc->type = BC_NEUMANN;
    bc->mass_flux = 0.0;
}

void concentration_bc_mass_flux(ConcentrationBC *bc, double mass_flux) {
    if (!bc) return;
    memset(bc, 0, sizeof(ConcentrationBC));
    bc->type = BC_NEUMANN;
    bc->mass_flux = mass_flux;
}

void concentration_bc_convective(ConcentrationBC *bc, double k_ext,
                                  double C_ext) {
    if (!bc || k_ext < 0.0) return;
    memset(bc, 0, sizeof(ConcentrationBC));
    bc->type = BC_ROBIN;
    bc->k_ext = k_ext;
    bc->C_ext = C_ext;
}

void concentration_bc_surface_reaction(ConcentrationBC *bc, double k_rxn) {
    if (!bc || k_rxn < 0.0) return;
    memset(bc, 0, sizeof(ConcentrationBC));
    bc->type = BC_ROBIN;
    bc->reaction_rate_constant = k_rxn;
}

void concentration_bc_henry_law(ConcentrationBC *bc, double henry_const,
                                 double p_gas) {
    if (!bc || henry_const <= 0.0) return;
    memset(bc, 0, sizeof(ConcentrationBC));
    bc->type = BC_DIRICHLET;  /* concentration fixed by Henry's law */
    bc->henry_constant = henry_const;
    /* C_liquid = p_gas / H */
    if (henry_const > 0.0) {
        bc->C_wall = p_gas / henry_const;
    }
}

void concentration_bc_partition(ConcentrationBC *bc, double partition_coeff) {
    if (!bc || partition_coeff <= 0.0) return;
    memset(bc, 0, sizeof(ConcentrationBC));
    bc->type = BC_JUMP;  /* concentration jump at interface */
    bc->partition_coefficient = partition_coeff;
}

int concentration_bc_validate(const ConcentrationBC *bc) {
    if (!bc) return -4;
    if (bc->type == BC_DIRICHLET && bc->C_wall < 0.0) return -1;
    if (bc->type == BC_ROBIN && bc->k_ext < 0.0) return -2;
    if (bc->reaction_rate_constant < 0.0) return -3;
    return 0;
}

/* ============================================================================
 * L2: Fluid-Fluid Interface Conditions
 *
 * Reference: Slattery (1990) Chapter 1 - Interfacial Kinematics
 *            Batchelor (1967) "An Introduction to Fluid Dynamics"
 * ============================================================================ */

int interface_velocity_continuity(const double u_a[3], const double u_b[3],
                                   double *residual) {
    if (!u_a || !u_b) {
        if (residual) *residual = -1.0;
        return 0;
    }
    double dx = u_a[0] - u_b[0];
    double dy = u_a[1] - u_b[1];
    double dz = u_a[2] - u_b[2];
    double res = sqrt(dx*dx + dy*dy + dz*dz);
    if (residual) *residual = res;
    return (res < 1e-6) ? 1 : 0;
}

int interface_shear_stress_continuity(double mu_a, double mu_b,
                                       const double du_dn_a[3],
                                       const double du_dn_b[3],
                                       double surface_tension_grad,
                                       double *stress_residual) {
    if (!du_dn_a || !du_dn_b) {
        if (stress_residual) *stress_residual = -1.0;
        return 0;
    }
    /* Tangential stress in A: mu_a * du_t/dn (t = tangential component) */
    /* Tangential stress in B: mu_b * du_t/dn */
    /* Jump: mu_A*du_t/dn - mu_B*du_t/dn = grad_s(sigma) */

    /* For simplicity, assume gradients are tangential components only */
    double shear_a = mu_a * sqrt(du_dn_a[0]*du_dn_a[0]
                                  + du_dn_a[1]*du_dn_a[1]
                                  + du_dn_a[2]*du_dn_a[2]);
    double shear_b = mu_b * sqrt(du_dn_b[0]*du_dn_b[0]
                                  + du_dn_b[1]*du_dn_b[1]
                                  + du_dn_b[2]*du_dn_b[2]);

    double jump = shear_a - shear_b - surface_tension_grad;
    if (stress_residual) *stress_residual = fabs(jump);
    return (fabs(jump) < 1e-6) ? 1 : 0;
}

int interface_temperature_continuity(double T_a, double T_b, double *delta_T) {
    double dT = fabs(T_a - T_b);
    if (delta_T) *delta_T = dT;
    return (dT < 1e-3) ? 1 : 0;
}

int interface_heat_flux_continuity(double k_a, double k_b,
                                    double dT_dn_a, double dT_dn_b,
                                    double m_dot, double h_fg,
                                    double *flux_diff) {
    /* Continuous: k_A * dT_A/dn - k_B * dT_B/dn = m_dot * h_fg */
    double q_a = k_a * dT_dn_a;
    double q_b = k_b * dT_dn_b;
    double diff = q_a - q_b - m_dot * h_fg;
    if (flux_diff) *flux_diff = fabs(diff);
    return (fabs(diff) < 1e-3) ? 1 : 0;
}

void interface_henry_equilibrium(double henry_const, double x_liquid,
                                  double *p_interface) {
    if (!p_interface) return;
    if (henry_const <= 0.0 || x_liquid < 0.0 || x_liquid > 1.0) {
        *p_interface = 0.0;
        return;
    }
    *p_interface = henry_const * x_liquid;
}

void interface_raoult_equilibrium(double x_liquid, double activity_coeff,
                                   double p_sat, double *p_interface) {
    if (!p_interface) return;
    if (x_liquid < 0.0 || x_liquid > 1.0
        || p_sat <= 0.0 || activity_coeff < 0.0) {
        *p_interface = 0.0;
        return;
    }
    *p_interface = x_liquid * activity_coeff * p_sat;
}

void interface_partition_equilibrium(double K, double C_phase_B,
                                      double *C_phase_A) {
    if (!C_phase_A) return;
    if (K <= 0.0 || C_phase_B < 0.0) {
        *C_phase_A = 0.0;
        return;
    }
    *C_phase_A = K * C_phase_B;
}

/* ============================================================================
 * L4: Conservation-Based Matching Conditions
 * ============================================================================ */

int interface_mass_continuity(double rho_A, double rho_B,
                               double u_n_A, double u_n_B, double u_n_I,
                               double *mass_flux) {
    if (rho_A <= 0.0 || rho_B <= 0.0) {
        if (mass_flux) *mass_flux = 0.0;
        return 0;
    }

    double m_dot_A = rho_A * (u_n_A - u_n_I);
    double m_dot_B = rho_B * (u_n_B - u_n_I);

    if (mass_flux) *mass_flux = m_dot_A;

    double diff = fabs(m_dot_A - m_dot_B);
    double scale = fmax(fabs(m_dot_A), fabs(m_dot_B));
    if (scale < 1e-15) scale = 1.0;

    return (diff / scale < 1e-9) ? 1 : 0;
}

void interface_young_laplace(double sigma, double curvature_mean,
                              double *pressure_jump) {
    if (!pressure_jump) return;
    if (sigma < 0.0) { *pressure_jump = 0.0; return; }
    /* p_A - p_B = 2 * sigma * H (H = mean curvature) */
    *pressure_jump = 2.0 * sigma * curvature_mean;
}

void interface_marangoni_stress(double d_sigma_dT, double dT_ds,
                                 double *marangoni_stress) {
    if (!marangoni_stress) return;
    /* Marangoni stress magnitude: tau_s = (d_sigma/dT) * (dT/ds) */
    *marangoni_stress = d_sigma_dT * dT_ds;
    /* Convention: positive stress drives flow from hot to cold
     * (for typical fluids where d_sigma/dT < 0) */
}

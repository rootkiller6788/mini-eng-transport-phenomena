#ifndef CFD_TURBULENCE_H
#define CFD_TURBULENCE_H
#include "cfd_types.h"
#ifdef __cplusplus
extern "C" {
#endif
void cfd_kepsilon_init(CfdKEpsilonModel *model);
void cfd_komega_sst_init(CfdKOmegaSSTModel *model);
double cfd_eddy_viscosity_kepsilon(const CfdKEpsilonModel *model, double rho, double k, double epsilon);
double cfd_eddy_viscosity_komega(const CfdKOmegaSSTModel *model, double rho, double k, double omega, double strain_rate);
double cfd_production_k(const CfdVectorField2D *vel, const CfdGrid2D *grid, double mu_t, int i, int j);
double cfd_wall_function_uplus(double y_plus, double kappa, double B);
double cfd_friction_velocity(double tau_wall, double rho);
double cfd_yplus(double y, double u_tau, double nu);
double cfd_turbulent_bl_thickness(double x, double Re_x);
double cfd_turbulence_intensity(double k, double U_mean);
double cfd_turbulence_length_scale(double k, double epsilon, double C_mu);
#ifdef __cplusplus
}
#endif
#endif

/*
 * Turbulence Modeling Theory
 * ==========================
 * Reynolds decomposition: u = U + u, p = P + p
 * RANS equations: dU_i/dt + U_j*dU_i/dx_j = -dP/dx_i + nu*d^2U_i/dx_j^2 - d(u_i*u_j)/dx_j
 * Reynolds stress: tau_ij = -rho * mean{u_i * u_j}
 *
 * Boussinesq hypothesis (1877):
 *   tau_ij = 2*mu_t*S_ij - (2/3)*rho*k*delta_ij
 *   S_ij = 0.5*(dU_i/dx_j + dU_j/dx_i)
 *
 * k-epsilon model (Launder & Spalding 1974):
 *   Two-equation model: transport equations for k and epsilon.
 *   mu_t = C_mu * rho * k^2 / epsilon
 *   Standard for industrial CFD (broad validation).
 *   Weakness: poor near-wall behavior, requires wall functions.
 *
 * k-omega SST (Menter 1994):
 *   Blends k-omega (near wall) with k-epsilon (far field).
 *   Better separation prediction, no wall functions needed for y+<1.
 *
 * Wall functions (Launder & Spalding 1974):
 *   Log-law: u+ = (1/kappa)*ln(y+) + B for y+ > 30
 *   Viscous sublayer: u+ = y+ for y+ < 5
 *   Buffer layer: 5 < y+ < 30 (blended)
 *
 * Kolmogorov scales:
 *   eta = (nu^3/epsilon)^(1/4), tau_eta = (nu/epsilon)^(1/2)
 *   u_eta = (nu*epsilon)^(1/4)
 *   Re_eta = u_eta*eta/nu = 1 (by definition)
 */

/* LES subgrid models */
double cfd_les_smagorinsky_viscosity(double rho, double Delta, double S, double Cs);
double cfd_les_wale_viscosity(double rho, double Delta, double S, double Cw);
double cfd_les_dynamic_cs(double Lij, double Mij);
double cfd_les_sigma_viscosity(double rho, double Delta, double sigma3, double C);
double cfd_les_vreman_viscosity(double rho, double Delta, double B_beta, double C);
/* Transition modeling */
double cfd_transition_intermittency(double Re_theta, double Re_theta_c);
double cfd_abughannam_shaw_Re_theta(double Tu, double lambda_theta);
double cfd_michel_transition_Re_theta(double Re_x);
/* DES */
double cfd_des_length_scale(double d_wall, double Delta, double C_DES);
double cfd_ddes_shielding(double r_d);
/* RSM */
double cfd_rotta_return_isotropy(double C_R, double epsilon, double k, double ui_uj_bar);
double cfd_lumley_flatness(double II_b, double III_b);
/* Realizability */
int cfd_check_realizability(double uu, double vv, double ww, double uv, double uw, double vw);
double cfd_durbin_realizability_limit(double rho, double k, double S_mag);
/* Init and spectrum */
void cfd_init_uniform_turbulence(double U_inf, double TI, double mu_ratio, double L, double rho, double mu, double *k, double *epsilon, double *omega, double *mu_t);
double cfd_von_karman_spectrum(double wavenumber, double epsilon, double eta, double A);
/* Mixing length */
double cfd_mixing_length_viscosity(double rho, double lm, double dudy);
double cfd_sa_production(double C_b1, double S_tilde, double nu_tilde);
double cfd_sa_destruction(double C_w1, double f_w, double nu_tilde, double d);
double cfd_sst_blending_f1(double y, double k, double omega, double rho, double mu, double sigma_omega2, double CD_kw);
double cfd_sst_cross_diffusion(double rho, double sigma_omega2, double omega, double dk_dx, double domega_dx, double dk_dy, double domega_dy);
void cfd_inlet_turbulence(double U_inf, double TI, double length_scale, double C_mu, double *k, double *epsilon);
void cfd_lumley_invariants(double b11, double b22, double b33, double b12, double b13, double b23, double *II_b, double *III_b);
double cfd_taylor_microscale(double nu, double k, double epsilon);
double cfd_kolmogorov_scale(double nu, double epsilon);
double cfd_integral_length_scale(double k, double epsilon);
double cfd_log_law_velocity(double y_plus);
double cfd_colebrook_friction_factor(double Re, double roughness, double D, int max_iter, double tol);
double cfd_darcy_to_skin_friction(double f_D);
double cfd_roughness_wall_shift(double k_s_plus, double C_s);
double cfd_baldwin_lomax_inner(double rho, double y, double omega_mag, double y_plus, double kappa, double A_plus);
double cfd_van_driest_damping(double y, double y_plus, double kappa, double A_plus);
double cfd_wall_shear_stress(double u_tau, double rho);
double cfd_van_driest_transformation(double u, double nu, double u_tau);
double cfd_clauser_chart_method(double Cf, double Re_x);

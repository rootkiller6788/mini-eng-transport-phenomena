#include "cfd_turbulence.h"
#include "cfd_discretization.h"
#include <math.h>
#include <stdlib.h>

/* ===== k-epsilon Model ===== */

void cfd_kepsilon_init(CfdKEpsilonModel *model)
{
    if (!model) return;
    model->C_mu = 0.09;
    model->C1_eps = 1.44;
    model->C2_eps = 1.92;
    model->sigma_k = 1.0;
    model->sigma_eps = 1.3;
    model->k_min = 1e-14;
    model->eps_min = 1e-20;
}

double cfd_eddy_viscosity_kepsilon(const CfdKEpsilonModel *model,
                                    double rho, double k, double epsilon)
{
    if (!model) return 0.0;
    double k_clip = CFD_MAX(k, model->k_min);
    double eps_clip = CFD_MAX(epsilon, model->eps_min);
    return model->C_mu * rho * k_clip * k_clip / eps_clip;
}

double cfd_production_k(const CfdVectorField2D *vel, const CfdGrid2D *grid,
                        double mu_t, int i, int j)
{
    if (!vel || !grid) return 0.0;
    int nx = vel->nx, ny = vel->ny;
    if (i < 0 || i >= nx || j < 0 || j >= ny) return 0.0;

    double dudx, dudy, dvdx, dvdy;
    if (i > 0 && i < nx - 1) {
        dudx = (vel->u[CFD_IDX2(i+1,j,nx)] - vel->u[CFD_IDX2(i-1,j,nx)]) / (2.0*grid->dx);
        dvdx = (vel->v[CFD_IDX2(i+1,j,nx)] - vel->v[CFD_IDX2(i-1,j,nx)]) / (2.0*grid->dx);
    } else { dudx = 0.0; dvdx = 0.0; }
    if (j > 0 && j < ny - 1) {
        dudy = (vel->u[CFD_IDX2(i,j+1,nx)] - vel->u[CFD_IDX2(i,j-1,nx)]) / (2.0*grid->dy);
        dvdy = (vel->v[CFD_IDX2(i,j+1,nx)] - vel->v[CFD_IDX2(i,j-1,nx)]) / (2.0*grid->dy);
    } else { dudy = 0.0; dvdy = 0.0; }

    /* P_k = mu_t * (2*S_ij*S_ij) */
    double S11 = dudx;
    double S22 = dvdy;
    double S12 = 0.5 * (dudy + dvdx);
    return mu_t * (2.0*S11*S11 + 2.0*S22*S22 + 4.0*S12*S12);
}

/* ===== k-omega SST Model ===== */

void cfd_komega_sst_init(CfdKOmegaSSTModel *model)
{
    if (!model) return;
    model->beta_star = 0.09;
    model->beta1 = 0.075;
    model->beta2 = 0.0828;
    model->sigma_k1 = 0.85;
    model->sigma_k2 = 1.0;
    model->sigma_omega1 = 0.5;
    model->sigma_omega2 = 0.856;
    model->a1 = 0.31;
    model->k_min = 1e-14;
    model->omega_min = 1e-20;
}

double cfd_eddy_viscosity_komega(const CfdKOmegaSSTModel *model,
                                   double rho, double k, double omega,
                                   double strain_rate)
{
    if (!model) return 0.0;
    double k_clip = CFD_MAX(k, model->k_min);
    double w_clip = CFD_MAX(omega, model->omega_min);
    /* SST: nu_t = a1*k / max(a1*omega, SF2) */
    double arg = CFD_MAX(model->a1 * w_clip, strain_rate);
    return rho * model->a1 * k_clip / arg;
}

/* ===== Wall Functions ===== */

double cfd_wall_function_uplus(double y_plus, double kappa, double B)
{
    /* Log-law: u+ = (1/kappa)*ln(y+) + B  for y+ > 30 */
    if (y_plus <= 0.0) return 0.0;
    if (y_plus < 11.225) return y_plus; /* Viscous sublayer: u+ = y+ */
    return log(y_plus) / kappa + B;
}

double cfd_friction_velocity(double tau_wall, double rho)
{
    /* u_tau = sqrt(tau_wall / rho) */
    if (rho <= 0.0) return 0.0;
    return sqrt(fabs(tau_wall) / rho);
}

double cfd_yplus(double y, double u_tau, double nu)
{
    /* y+ = y * u_tau / nu */
    if (nu <= 0.0) return 0.0;
    return y * u_tau / nu;
}

double cfd_turbulent_bl_thickness(double x, double Re_x)
{
    /* Turbulent BL thickness (1/7th power law):
     * delta_turb = 0.37 * x / Re_x^(1/5) */
    if (x <= 0.0 || Re_x <= 0.0) return 0.0;
    return 0.37 * x / pow(Re_x, 0.2);
}

/* ===== Turbulence Statistics ===== */

double cfd_turbulence_intensity(double k, double U_mean)
{
    /* TI = sqrt(2k/3) / |U_mean| */
    if (U_mean == 0.0) return 0.0;
    return sqrt(2.0 * CFD_MAX(k, 0.0) / 3.0) / fabs(U_mean);
}

double cfd_turbulence_length_scale(double k, double epsilon, double C_mu)
{
    /* Integral length scale: L = C_mu^(3/4) * k^(3/2) / epsilon */
    if (epsilon <= 0.0 || k <= 0.0) return 0.0;
    double Cmu34 = pow(C_mu, 0.75);
    return Cmu34 * pow(k, 1.5) / epsilon;
}

/* ===== Algebraic Turbulence Models ===== */

double cfd_mixing_length_viscosity(double rho, double lm, double dudy)
{
    if (rho <= 0.0 || lm < 0.0) return 0.0;
    return rho * lm * lm * fabs(dudy);
}

double cfd_van_driest_damping(double y, double y_plus, double kappa, double A_plus)
{
    if (y_plus <= 0.0) return 0.0;
    return kappa * y * (1.0 - exp(-y_plus / A_plus));
}

double cfd_baldwin_lomax_inner(double rho, double y, double omega_mag,
                                double y_plus, double kappa, double A_plus)
{
    double lm = cfd_van_driest_damping(y, y_plus, kappa, A_plus);
    return rho * lm * lm * fabs(omega_mag);
}

double cfd_sa_production(double C_b1, double S_tilde, double nu_tilde)
{
    return C_b1 * S_tilde * nu_tilde;
}

double cfd_sa_destruction(double C_w1, double f_w, double nu_tilde, double d)
{
    if (d <= 0.0) return 0.0;
    double ratio = nu_tilde / (d * d);
    return C_w1 * f_w * ratio * nu_tilde;
}

double cfd_sst_blending_f1(double y, double k, double omega, double rho,
                            double mu, double sigma_omega2, double CD_kw)
{
    if (y <= 0.0 || omega <= 0.0) return 1.0;
    double arg1 = fmin(fmax(sqrt(k)/(0.09*omega*y), 500.0*mu/(rho*omega*y*y)),
                         4.0*rho*sigma_omega2*k/(CD_kw*y*y));
    return tanh(arg1 * arg1 * arg1 * arg1);
}

double cfd_sst_cross_diffusion(double rho, double sigma_omega2, double omega,
                                 double dk_dx, double domega_dx,
                                 double dk_dy, double domega_dy)
{
    double denom = CFD_MAX(omega, 1e-15);
    return 2.0 * rho * sigma_omega2 / denom * (dk_dx * domega_dx + dk_dy * domega_dy);
}

void cfd_inlet_turbulence(double U_inf, double TI, double length_scale,
                           double C_mu, double *k, double *epsilon)
{
    *k = 1.5 * (TI * U_inf) * (TI * U_inf);
    double Cmu34 = pow(C_mu, 0.75);
    *epsilon = Cmu34 * pow(*k, 1.5) / CFD_MAX(length_scale, 1e-15);
}

void cfd_lumley_invariants(double b11, double b22, double b33,
                            double b12, double b13, double b23,
                            double *II_b, double *III_b)
{
    double b21=b12, b31=b13, b32=b23;
    *II_b = b11*b11 + b22*b22 + b33*b33 + 2.0*(b12*b12 + b13*b13 + b23*b23);
    double bsq11 = b11*b11 + b12*b12 + b13*b13;
    double bsq21 = b21*b11 + b22*b21 + b23*b31;
    double bsq31 = b31*b11 + b32*b21 + b33*b31;
    *III_b = b11*bsq11 + b21*bsq21 + b31*bsq31;
}

double cfd_taylor_microscale(double nu, double k, double epsilon)
{
    if (epsilon <= 0.0 || k <= 0.0) return 0.0;
    return sqrt(10.0 * nu * k / epsilon);
}

double cfd_kolmogorov_scale(double nu, double epsilon)
{
    if (epsilon <= 0.0) return 0.0;
    return pow(nu * nu * nu / epsilon, 0.25);
}

double cfd_integral_length_scale(double k, double epsilon)
{
    if (epsilon <= 0.0 || k <= 0.0) return 0.0;
    return pow(k, 1.5) / epsilon;
}

double cfd_log_law_velocity(double y_plus)
{
    static const double kappa = 0.41;
    static const double B = 5.2;
    return cfd_wall_function_uplus(y_plus, kappa, B);
}

double cfd_colebrook_friction_factor(double Re, double roughness, double D,
                                       int max_iter, double tol)
{
    if (Re <= 2300.0) return 64.0 / Re;
    double f = 0.02;
    for (int iter = 0; iter < max_iter; iter++) {
        double lhs = 1.0 / sqrt(f);
        double rhs = -2.0 * log10(roughness/(3.7*D) + 2.51/(Re*sqrt(f)));
        double f_new = 1.0 / (rhs * rhs);
        if (fabs(f_new - f) < tol) return f_new;
        f = f_new;
    }
    return f;
}

double cfd_darcy_to_skin_friction(double f_D)
{
    return f_D / 4.0;
}

double cfd_roughness_wall_shift(double k_s_plus, double C_s)
{
    if (k_s_plus <= 0.0) return 0.0;
    return log(1.0 + C_s * k_s_plus) / 0.41;
}

/* ===== LES Subgrid Models ===== */
double cfd_les_smagorinsky_viscosity(double rho, double Delta, double S, double Cs) {
    if (Delta <= 0.0) return 0.0;
    return rho * Cs * Cs * Delta * Delta * fabs(S);
}
double cfd_les_wale_viscosity(double rho, double Delta, double S, double Cw) {
    if (Delta <= 0.0) return 0.0;
    double Sd = S * S; return rho * Cw * Cw * Delta * Delta * sqrt(Sd);
}
double cfd_les_dynamic_cs(double Lij, double Mij) {
    if (fabs(Mij) < 1e-15) return 0.1;
    return -0.5 * Lij / Mij;
}
double cfd_les_sigma_viscosity(double rho, double Delta, double sigma3, double C) {
    return rho * C * C * Delta * Delta * sigma3;
}
double cfd_les_vreman_viscosity(double rho, double Delta, double B_beta, double C) {
    return rho * C * C * sqrt(B_beta);
}

/* ===== Transition Modeling ===== */
double cfd_transition_intermittency(double Re_theta, double Re_theta_c) {
    if (Re_theta <= 0.0 || Re_theta_c <= 0.0) return 0.0;
    double ratio = (Re_theta - Re_theta_c) / Re_theta_c;
    if (ratio <= 0.0) return 0.0; if (ratio >= 1.0) return 1.0;
    return ratio * ratio * (3.0 - 2.0 * ratio);
}
double cfd_abughannam_shaw_Re_theta(double Tu, double lambda_theta) {
    double Tu_pct = Tu * 100.0; if (Tu_pct <= 0.1) Tu_pct = 0.1;
    double Re_theta_t = 163.0 + exp(6.91 - Tu_pct);
    return Re_theta_t * (1.0 + 0.5 * lambda_theta);
}
double cfd_michel_transition_Re_theta(double Re_x) {
    return 1.174 * pow(Re_x, 0.46) * (1.0 + 22400.0 / Re_x);
}

/* ===== DES and DDES ===== */
double cfd_des_length_scale(double d_wall, double Delta, double C_DES) {
    return CFD_MIN(d_wall, C_DES * Delta);
}
double cfd_ddes_shielding(double r_d) {
    return 1.0 - tanh(8.0 * r_d * r_d * r_d);
}

/* ===== RSM Basics ===== */
double cfd_rotta_return_isotropy(double C_R, double epsilon, double k, double ui_uj_bar) {
    if (k <= 0.0) return 0.0; return -C_R * epsilon / k * (ui_uj_bar - 2.0 * k / 3.0);
}
double cfd_lumley_flatness(double II_b, double III_b) {
    return 1.0 - 1.125 * (II_b - III_b);
}

/* ===== Realizability ===== */
int cfd_check_realizability(double uu, double vv, double ww, double uv, double uw, double vw) {
    if (uu < 0.0 || vv < 0.0 || ww < 0.0) return 0;
    if (uv*uv > uu*vv || uw*uw > uu*ww || vw*vw > vv*ww) return 0;
    return 1;
}
double cfd_durbin_realizability_limit(double rho, double k, double S_mag) {
    if (S_mag <= 0.0) return 1e10; return 0.31 * rho * k / S_mag;
}

/* ===== Turbulence Init and Spectrum ===== */
void cfd_init_uniform_turbulence(double U_inf, double TI, double mu_ratio,
                                  double L, double rho, double mu,
                                  double *k, double *epsilon, double *omega,
                                  double *mu_t) {
    *k = 1.5 * (TI * U_inf) * (TI * U_inf);
    *mu_t = mu_ratio * mu;
    CfdKEpsilonModel ke; cfd_kepsilon_init(&ke);
    *epsilon = ke.C_mu * rho * (*k) * (*k) / (*mu_t);
    *omega = *epsilon / (*k);
}
double cfd_von_karman_spectrum(double wavenumber, double epsilon, double eta, double A) {
    if (wavenumber <= 0.0) return 0.0;
    double k_eta = wavenumber * eta;
    return A * pow(epsilon, 2.0/3.0) * pow(wavenumber, -5.0/3.0) * exp(-1.5*A*k_eta*k_eta*k_eta*k_eta);
}

/* ===== Additional Wall Function Variants ===== */
double cfd_wall_shear_stress(double u_tau, double rho) { return rho * u_tau * u_tau; }
double cfd_van_driest_transformation(double u, double nu, double u_tau) {
    if (u_tau <= 0.0) return 0.0; return u / u_tau;
}
double cfd_clauser_chart_method(double Cf, double Re_x) {
    double Cf_ref = 0.027 / pow(Re_x, 1.0/7.0);
    return Cf / Cf_ref;
}

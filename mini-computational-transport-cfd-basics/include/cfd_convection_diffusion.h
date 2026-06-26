#ifndef CFD_CONVECTION_DIFFUSION_H
#define CFD_CONVECTION_DIFFUSION_H
#include "cfd_types.h"
#ifdef __cplusplus
extern "C" {
#endif
int cfd_solve_convdiff_1d(const CfdGrid1D *grid, const CfdTransportEquation *eqn, double velocity, double phi_left, double phi_right, CfdConvectiveScheme scheme, double *phi);
int cfd_solve_convdiff_1d_explicit(const CfdGrid1D *grid, const CfdTransportEquation *eqn, double velocity, double phi_left, double phi_right, CfdConvectiveScheme space_scheme, double dt, int n_steps, double *phi);
int cfd_solve_convdiff_1d_implicit(const CfdGrid1D *grid, const CfdTransportEquation *eqn, double velocity, double phi_left, double phi_right, CfdConvectiveScheme space_scheme, double dt, int n_steps, double *phi);
double cfd_convdiff_1d_exact(double x, double L, double Pe, double phi_0, double phi_L);
int cfd_solve_convdiff_2d(const CfdGrid2D *grid, const CfdTransportEquation *eqn, const CfdVectorField2D *vel, const CfdBCSet2D *bc, CfdConvectiveScheme scheme, CfdScalarField2D *phi);
int cfd_solve_convdiff_2d_explicit(const CfdGrid2D *grid, const CfdTransportEquation *eqn, const CfdVectorField2D *vel, const CfdBCSet2D *bc, CfdConvectiveScheme scheme, double dt, int n_steps, CfdScalarField2D *phi);
int cfd_solve_poisson_2d(const CfdGrid2D *grid, const CfdScalarField2D *source, const CfdBCSet2D *bc, CfdScalarField2D *phi);
int cfd_simple_iteration(const CfdGrid2D *grid, CfdVectorField2D *vel, const CfdFluidProperty *fp, CfdSimpleState *simple);
int cfd_simple_solve(const CfdGrid2D *grid, CfdVectorField2D *vel, const CfdFluidProperty *fp, CfdSimpleState *simple);
int cfd_solve_lid_driven_cavity(const CfdGrid2D *grid, CfdVectorField2D *vel, const CfdFluidProperty *fp, double U_lid, CfdSimpleState *state);
#ifdef __cplusplus
}
#endif
#endif

/*
 * Convection-Diffusion Theory
 * ===========================
 * 1D steady: u*d(phi)/dx = Gamma*d^2(phi)/dx^2
 * Exact solution: phi(x) = phi_0 + (phi_L-phi_0)*(exp(Pe*x/L)-1)/(exp(Pe)-1)
 * where Pe = u*L/Gamma is the global Peclet number.
 *
 * Boundary layer behavior (Pe >> 1):
 *   - Upstream: phi ~ phi_0, exponential growth near outlet.
 *   - Downstream boundary layer thickness ~ L/Pe.
 *
 * Godunov theorem (1959):
 *   Linear monotone schemes for advection are at most 1st-order.
 *   High-order schemes must be non-linear (TVD, WENO, etc.)
 *
 * SIMPLE Algorithm (Patankar & Spalding 1972):
 *   1. Guess pressure field p*.
 *   2. Solve momentum equations for u*, v*.
 *   3. Solve pressure correction equation for p.
 *   4. Correct velocities: u = u* + du, v = v* + dv.
 *   5. Correct pressure: p = p* + alpha_p * p.
 *   6. Repeat until convergence.
 *
 * Rhie-Chow interpolation (1983):
 *   Prevents checkerboard pressure on collocated grids.
 *   Face velocity = average(neighbor velocities) + pressure smoothing term.
 */

/* Additional API: 2D Gauss-Seidel convection-diffusion solver */
int cfd_solve_convdiff_2d_gs(const CfdGrid2D *grid, const CfdTransportEquation *eqn,
                             const CfdVectorField2D *vel, const CfdBCSet2D *bc,
                             CfdScalarField2D *phi, int max_iter, double tol);
/* Additional API: post-processing utilities */
double cfd_streamfunction_2d(const CfdVectorField2D *vel, int i, int j, const CfdGrid2D *grid);
double cfd_vorticity_2d(const CfdVectorField2D *vel, int i, int j, const CfdGrid2D *grid);
double cfd_pressure_drop_1d(const double *p, int nx, double dx);
double cfd_nusselt_from_field(const CfdScalarField2D *T, const CfdGrid2D *grid, double T_wall, double T_bulk, double L, double k);
double cfd_chilton_colburn_j_factor(double Nu, double Re, double Pr);
double cfd_sherwood_from_j_factor(double j_D, double Re, double Sc);
double cfd_wall_shear_from_velocity(const CfdVectorField2D *vel, const CfdGrid2D *grid, double mu, int wall_face);
double cfd_energy_balance_2d(const CfdVectorField2D *vel, const CfdGrid2D *grid, const CfdBCSet2D *bc, double k);
/* Additional API: L2 error and Richardson */
double cfd_l2_error_1d(const double *num, const CfdGrid1D *grid, double (*exact)(double, void *), void *params);
double cfd_richardson_extrapolation(double f_h, double f_h2, double h_ratio);
double cfd_order_of_accuracy(double err_h, double err_h2, double h_ratio);
double cfd_grid_convergence_index(double err_fine, double err_coarse, double h_ratio, double p, double Fs);
double cfd_mass_flux_balance_2d(const CfdVectorField2D *vel, const CfdGrid2D *grid);
double cfd_kinetic_energy_2d(const CfdVectorField2D *vel);
double cfd_enstrophy_2d(const CfdVectorField2D *vel, const CfdGrid2D *grid);

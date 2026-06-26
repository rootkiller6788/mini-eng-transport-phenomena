#ifndef CFD_DISCRETIZATION_H
#define CFD_DISCRETIZATION_H
#include "cfd_types.h"
#ifdef __cplusplus
extern "C" {
#endif
double cfd_central_diff_1d(double u_west, double u_east, double dx);
double cfd_upwind_diff_1d(double u_w, double u_p, double dx, double velocity);
double cfd_second_deriv_1d(double u_w, double u_p, double u_e, double dx);
double cfd_grad_x_2d(const CfdScalarField2D *field, int i, int j, const CfdGrid2D *grid);
double cfd_grad_y_2d(const CfdScalarField2D *field, int i, int j, const CfdGrid2D *grid);
double cfd_grad_mag_2d(const CfdScalarField2D *field, int i, int j, const CfdGrid2D *grid);
double cfd_laplacian_2d(const CfdScalarField2D *field, int i, int j, const CfdGrid2D *grid);
double cfd_divergence_2d(const CfdVectorField2D *vel, int i, int j, const CfdGrid2D *grid);
void cfd_powerlaw_1d(double Pe, double F, double D, double *aW, double *aE, double *aP);
void cfd_hybrid_1d(double Pe, double F, double D, double *aW, double *aE, double *aP);
void cfd_assemble_convdiff_1d(const CfdGrid1D *grid, const CfdTransportEquation *eqn, double velocity, double phi_left, double phi_right, CfdConvectiveScheme scheme, CfdMatrix2D *mat);
void cfd_assemble_poisson_2d(const CfdGrid2D *grid, const CfdScalarField2D *source, const CfdBCSet2D *bc, CfdMatrix2D *mat);
void cfd_assemble_u_momentum_2d(const CfdGrid2D *grid, const CfdVectorField2D *vel, const CfdFluidProperty *fp, CfdConvectiveScheme scheme, CfdMatrix2D *u_mat);
void cfd_assemble_v_momentum_2d(const CfdGrid2D *grid, const CfdVectorField2D *vel, const CfdFluidProperty *fp, CfdConvectiveScheme scheme, CfdMatrix2D *v_mat);
double cfd_stable_timestep_1d(double velocity, double diffusivity, double dx, double CFL_target);
double cfd_cell_peclet(double velocity, double dx, double diffusivity);
int cfd_check_maximum_principle(const CfdMatrix2D *mat);
#ifdef __cplusplus
}
#endif
#endif

/*
 * Discretization Theory Summary
 * =============================
 *
 * Finite Difference Method (FDM):
 *   Approximates derivatives using Taylor series expansions.
 *   du/dx = (u_{i+1} - u_{i-1})/(2*dx) - (dx^2/6)*u_triple_prime + O(dx^4)
 *   d^2u/dx^2 = (u_{i+1} - 2*u_i + u_{i-1})/dx^2 - (dx^2/12)*u_quad_prime + O(dx^4)
 *
 * Finite Volume Method (FVM):
 *   Integrates conservation laws over control volumes.
 *   Integral form: d/dt(int_V(rho*phi)dV) + int_S(rho*phi*u*n)dS
 *                 = int_S(Gamma*grad(phi)*n)dS + int_V(S)dV
 *   Naturally conservative: fluxes cancel at interior faces.
 *
 * Convection Schemes:
 *   - Central (CDS): 2nd-order, oscillatory for Pe > 2.
 *   - Upwind (UDS): 1st-order, monotone but diffusive.
 *   - Hybrid: UDS for |Pe| > 2, CDS otherwise (Spalding 1972).
 *   - Power-law: Near-exact for 1D (Patankar 1980).
 *   - QUICK: 3rd-order upwind-biased (Leonard 1979).
 *   - TVD: Total Variation Diminishing (Harten 1983).
 *
 * Lax Equivalence Theorem (Lax & Richtmyer 1956):
 *   For a well-posed linear IVP with a consistent scheme:
 *   Stability <=> Convergence
 *
 * Discrete Maximum Principle:
 *   If A is an M-matrix (a_ii > 0, a_ij <= 0 for i != j,
 *   diagonal dominant), then Au = f has no interior extrema.
 *   This guarantees physical (non-oscillatory) solutions.
 */

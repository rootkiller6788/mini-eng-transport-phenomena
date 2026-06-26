/**
 * interphase_numerical.h
 * =======================
 * Numerical methods for interface tracking, reconstruction, and
 * boundary condition implementation in interphase transport.
 *
 * Reference: Tryggvason, Scardovelli, Zaleski (2011) "Direct Numerical
 *   Simulations of Gas-Liquid Multiphase Flows"
 *
 * Osher & Fedkiw (2003) "Level Set Methods and Dynamic Implicit Surfaces"
 * Hirt & Nichols (1981) "Volume of Fluid (VOF) Method"
 * Peskin (2002) "The Immersed Boundary Method" (Acta Numerica)
 *
 * Knowledge coverage: L5 Engineering Methods, L8 Advanced Methods
 */

#ifndef INTERPHASE_NUMERICAL_H
#define INTERPHASE_NUMERICAL_H

#include "interphase_types.h"
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * L5: Interface Representation Methods
 * ============================================================================ */

/**
 * Level-set signed distance function at a point.
 * phi(x) = signed distance to interface, phi<0 in phase A, phi>0 in phase B.
 * Interface is the zero level set: phi(x,t) = 0.
 *
 * @param x              Point coordinates [dim] to evaluate
 * @param dim            Spatial dimension (2 or 3)
 * @param interface_center Interface reference point [dim]
 * @param normal         Interface normal direction [dim] (unit)
 * @param offset         Distance of interface from coordinate origin along normal
 * @return Signed distance, negative in phase A, positive in phase B
 *
 * Theorem: |grad(phi)| = 1 for a true signed distance function.
 * Theorem: phi is reinitialized each timestep to maintain signed distance property.
 */
double level_set_signed_distance(const double *x, int dim,
                                  const double *interface_center,
                                  const double *normal,
                                  double offset);

/**
 * Regularized Heaviside function for smooth transition across interface.
 * H_eps(phi) = 0 for phi < -eps, 1 for phi > +eps,
 *              = 0.5*(1 + phi/eps + sin(pi*phi/eps)/pi) for |phi| <= eps.
 * Used to smoothly interpolate properties across the interface.
 *
 * @param phi Signed distance
 * @param eps Half-thickness of transition zone (typically 1.5*dx)
 * @return Regularized Heaviside value in [0, 1]
 * Complexity: O(1)
 *
 * Theorem: H_eps -> H (sharp step) as eps -> 0.
 * Theorem: dH_eps/dphi = delta_eps (regularized delta function).
 */
double regularized_heaviside(double phi, double eps);

/**
 * Regularized delta function (derivative of regularized Heaviside).
 * delta_eps(phi) = 0 for |phi| > eps,
 *                  = 0.5*(1 + cos(pi*phi/eps))/eps for |phi| <= eps.
 *
 * @param phi Signed distance
 * @param eps Half-thickness of transition zone
 * @return Regularized delta value, integrated over phi gives 1
 * Complexity: O(1)
 */
double regularized_delta(double phi, double eps);

/**
 * Smooth property interpolation across interface using Heaviside.
 * beta(phi) = beta_A + (beta_B - beta_A) * H_eps(phi).
 *
 * @param phi      Signed distance
 * @param eps      Transition half-thickness
 * @param beta_A   Property value in phase A
 * @param beta_B   Property value in phase B
 * @return Interpolated property value
 * Complexity: O(1)
 */
double level_set_interpolate(double phi, double eps,
                              double beta_A, double beta_B);

/**
 * Compute interface normal from level-set function using central differences.
 * n = grad(phi) / |grad(phi)|.
 *
 * @param phi      Level-set field at grid point (i,j,k) and neighbors
 * @param nx       Number of grid points in x
 * @param ny       Number of grid points in y
 * @param nz       Number of grid points in z (1 for 2D)
 * @param i        x-index of central point (0-indexed)
 * @param j        y-index of central point
 * @param k        z-index of central point (0 for 2D)
 * @param dx       Grid spacing [dim]
 * @param normal   Output unit normal vector [3] (z component = 0 for 2D)
 * Complexity: O(1)
 */
void level_set_normal(const double *phi, int nx, int ny, int nz,
                      int i, int j, int k,
                      const double *dx, double normal[3]);

/**
 * Compute mean curvature from level-set function.
 * kappa = div(grad(phi)/|grad(phi)|) = (phi_xx*phi_y^2 - 2*phi_x*phi_y*phi_xy + phi_yy*phi_x^2) / (phi_x^2 + phi_y^2)^(3/2).
 *
 * @param phi       Level-set field at grid point and neighbors
 * @param nx        Grid points in x
 * @param ny        Grid points in y
 * @param i         x-index of central point
 * @param j         y-index of central point
 * @param dx        Grid spacing [2] (dx, dy)
 * @param curvature Output mean curvature (1/m), + for convex toward phase A
 * Complexity: O(1)
 */
void level_set_curvature_2d(const double *phi, int nx, int ny,
                             int i, int j, const double dx[2],
                             double *curvature);

/* ============================================================================
 * L5: Volume of Fluid (VOF) Method
 * ============================================================================ */

/**
 * Compute VOF fraction from interface geometry using PLIC
 * (Piecewise Linear Interface Calculation).
 * Given a known interface normal and a line offset, compute the
 * fractional volume of phase B in a cell.
 *
 * @param normal     Interface normal [2] pointing toward phase B (unit)
 * @param alpha      Line constant: n_x*x + n_y*y = alpha defines interface
 * @param dx         Cell dimensions [2]
 * @param x0         Cell lower-left corner [2]
 * @param F          Output VOF fraction in [0, 1]
 * Complexity: O(1)
 *
 * Theorem: F = 0 -> cell is full of phase A.
 * Theorem: F = 1 -> cell is full of phase B.
 * Theorem: F = 0.5 for interface through cell center.
 */
void plic_vof_fraction_2d(const double normal[2], double alpha,
                           const double dx[2], const double x0[2],
                           double *F);

/**
 * Reconstruct interface line from VOF fraction using inverse PLIC.
 * Given F in a cell and interface normal, find alpha such that
 * PLIC(n, alpha, cell) = F.
 *
 * @param normal     Interface normal [2] (unit)
 * @param F          VOF fraction in cell, [0, 1]
 * @param dx         Cell dimensions [2]
 * @param alpha      Output line constant
 * Complexity: O(log(1/eps)) using bisection
 *
 * Theorem: The mapping F -> alpha is monotonic for fixed normal direction.
 * Theorem: Reconstruction is unique for each (normal, F) pair.
 */
void plic_inverse_vof_2d(const double normal[2], double F,
                          const double dx[2], double *alpha);

/**
 * Compute VOF flux across a cell face for advection.
 * Flux = area of phase B crossing the face per unit time.
 * Uses geometric flux computation based on the PLIC reconstruction.
 *
 * @param normal      Interface normal [2]
 * @param alpha       Interface line constant
 * @param dx          Cell size [2]
 * @param dt          Timestep (s)
 * @param u_face      Velocity at face [2] (m/s)
 * @param face_normal Face normal direction index (0=x, 1=y)
 * @param flux        Output VOF flux (volume fraction * area)
 * Complexity: O(1)
 */
void vof_flux_2d(const double normal[2], double alpha,
                  const double dx[2], double dt,
                  const double u_face[2], int face_normal,
                  double *flux);

/**
 * Compute interface area within a cell from VOF reconstruction.
 * Uses PLIC to find line segment length in 2D (or polygon area in 3D).
 *
 * @param normal      Interface normal [2]
 * @param alpha       Interface line constant
 * @param dx          Cell size [2]
 * @param area        Output interface length in 2D (area in 3D)
 * Complexity: O(1)
 */
void plic_interface_area_2d(const double normal[2], double alpha,
                             const double dx[2], double *area);

/* ============================================================================
 * L8: Ghost Fluid Method (GFM) for Sharp Interface BCs
 * ============================================================================ */

/**
 * Populate ghost cells for a Dirichlet boundary condition at an interface
 * using the Ghost Fluid Method (Fedkiw et al., 1999).
 * Ghost value is extrapolated across the interface to enforce
 * u_interface = specified_value.
 *
 * @param phi_real    Real node signed distance
 * @param phi_ghost   Ghost node signed distance
 * @param value_real  Real node value
 * @param value_interface Specified interface value
 * @param value_ghost Output ghost node value
 * Complexity: O(1)
 *
 * Theorem: GFM gives O(dx) accuracy for smooth interfaces.
 * Theorem: GFM avoids smearing the discontinuity across the interface.
 */
void gfm_dirichlet(double phi_real, double phi_ghost,
                    double value_real, double value_interface,
                    double *value_ghost);

/**
 * Populate ghost cells for a Neumann boundary condition at an interface.
 * Used to enforce specified gradient (e.g., zero flux for insulation).
 *
 * @param phi_real         Real node signed distance
 * @param phi_ghost        Ghost node signed distance
 * @param value_real       Real node value
 * @param gradient_interface Specified normal gradient at interface
 * @param value_ghost      Output ghost node value
 * Complexity: O(1)
 */
void gfm_neumann(double phi_real, double phi_ghost,
                  double value_real, double gradient_interface,
                  double *value_ghost);

/**
 * Populate ghost cells for a Robin (mixed) boundary condition at an interface.
 * Enforces: a*u + b*du/dn = c at the interface.
 *
 * @param phi_real    Real node signed distance
 * @param phi_ghost   Ghost node signed distance
 * @param value_real  Real node value
 * @param a_coeff     Robin coefficient a
 * @param b_coeff     Robin coefficient b
 * @param c_coeff     Robin RHS c
 * @param value_ghost Output ghost node value
 * Complexity: O(1)
 */
void gfm_robin(double phi_real, double phi_ghost,
                double value_real,
                double a_coeff, double b_coeff, double c_coeff,
                double *value_ghost);

/* ============================================================================
 * L8: Immersed Boundary Method (IBM)
 * ============================================================================ */

/**
 * Compute the regularized delta function for IBM force spreading.
 * Roma et al. (1999) 3-point delta function:
 * delta_2(r) = (1/3)*(1 + sqrt(-3*r^2+1)) for |r| <= 0.5,
 *            = (1/6)*(5 - 3*|r| - sqrt(-3*(1-|r|)^2+1)) for 0.5 < |r| <= 1.5,
 *            = 0 for |r| > 1.5.
 *
 * @param r  Normalized distance r = (x - X)/dx
 * @return Regularized delta function value
 * Complexity: O(1)
 *
 * Theorem: Integral of delta_2 over R is 1 (conservation).
 * Theorem: delta_2 has continuous first derivative.
 * Theorem: Support on [-1.5, 1.5] minimizes spurious oscillations.
 */
double ibm_delta_roma(double r);

/**
 * Spread a Lagrangian interface force to Eulerian grid using IBM.
 * F_euler(x) = sum_k F_lagrangian(X_k) * delta_h(x - X_k) * ds_k.
 *
 * @param x_grid        Eulerian grid points [nx*ny*2] (x,y coords)
 * @param nx            Number of x grid points
 * @param ny            Number of y grid points
 * @param dx            Grid spacing [2]
 * @param X_lagrangian  Lagrangian marker positions [n_markers*2]
 * @param F_lagrangian  Lagrangian forces [n_markers*2]
 * @param n_markers     Number of markers
 * @param ds            Marker spacing (arc length per marker)
 * @param F_euler       Output Eulerian force field [nx*ny*2]
 *                      (must be pre-allocated and zero-initialized)
 * Complexity: O(nx*ny*n_markers)
 *
 * Note: This is the fundamental IBM operation. The cost scales with
 *       the product of grid points and markers. For efficiency in
 *       production codes, narrow-band spreading is used.
 */
void ibm_spread_force(const double *x_grid, int nx, int ny,
                      const double dx[2],
                      const double *X_lagrangian,
                      const double *F_lagrangian,
                      int n_markers, double ds,
                      double *F_euler);

/**
 * Interpolate Eulerian velocity to Lagrangian markers using IBM.
 * U_lagrangian(X_k) = sum_ij u_euler(x_ij) * delta_h(x_ij - X_k) * dx*dy.
 *
 * @param u_euler       Eulerian velocity field [nx*ny*2]
 * @param nx            Number of x grid points
 * @param ny            Number of y grid points
 * @param dx            Grid spacing [2]
 * @param X_lagrangian  Lagrangian marker positions [n_markers*2]
 * @param n_markers     Number of markers
 * @param U_lagrangian  Output Lagrangian velocities [n_markers*2]
 * Complexity: O(nx*ny*n_markers)
 */
void ibm_interpolate_velocity(const double *u_euler, int nx, int ny,
                               const double dx[2],
                               const double *X_lagrangian,
                               int n_markers,
                               double *U_lagrangian);

/* ============================================================================
 * L5: Interface Advection Schemes
 * ============================================================================ */

/**
 * Advect the level-set function using first-order upwind scheme.
 * d(phi)/dt + u * grad(phi) = 0.
 * Upwind bias based on velocity sign.
 *
 * @param phi_n     Level-set field at timestep n [nx*ny]
 * @param phi_np1   Output level-set field at timestep n+1 [nx*ny]
 * @param u         x-velocity field [nx*ny]
 * @param v         y-velocity field [nx*ny]
 * @param nx        Grid points in x
 * @param ny        Grid points in y
 * @param dx        Grid spacing (assumed uniform) [2]
 * @param dt        Timestep (s)
 * Complexity: O(nx*ny)
 *
 * CFL condition: dt <= CFL * min(dx, dy) / max(|u|, |v|).
 * Stability requires CFL <= 1 for explicit upwind.
 */
void level_set_advect_upwind(const double *phi_n, double *phi_np1,
                              const double *u, const double *v,
                              int nx, int ny,
                              const double dx[2], double dt);

/**
 * Reinitialize the level-set function to signed distance.
 * Solves d(phi)/d(tau) + sign(phi0)*(|grad(phi)| - 1) = 0
 * to steady state using Godunov scheme (Sussman et al., 1994).
 *
 * @param phi        Level-set field (updated in place) [nx*ny]
 * @param nx         Grid points in x
 * @param ny         Grid points in y
 * @param dx         Grid spacing [2]
 * @param dtau       Pseudo-timestep (typically 0.5*min(dx,dy))
 * @param n_iter     Number of pseudo-time iterations
 * Complexity: O(n_iter * nx * ny)
 *
 * Theorem: Steady state satisfies |grad(phi)| = 1.
 * Theorem: Interface position (phi=0) is preserved.
 */
void level_set_reinitialize(double *phi, int nx, int ny,
                             const double dx[2],
                             double dtau, int n_iter);

/* ============================================================================
 * L5: Surface Tension Force Computation
 * ============================================================================ */

/**
 * Compute surface tension force from level-set function using CSF
 * (Continuum Surface Force) model (Brackbill et al., 1992).
 * F_sigma = sigma * kappa * delta(phi) * n
 * where kappa = curvature, n = normal, delta = smoothed delta.
 *
 * @param phi     Level-set field [nx*ny]
 * @param nx      Grid points in x
 * @param ny      Grid points in y
 * @param dx      Grid spacing [2]
 * @param sigma   Surface tension coefficient (N/m)
 * @param eps     Smoothing half-width for delta function
 * @param Fx      Output x-component of surface tension force [nx*ny]
 * @param Fy      Output y-component of surface tension force [nx*ny]
 * Complexity: O(nx*ny)
 *
 * Theorem: Integral of F_sigma over domain = sigma * interface_length (2D).
 * Theorem: CSF converges to sharp interface force as eps -> 0.
 * Theorem: Spurious currents scale with sigma/(rho*dx) for given eps.
 */
void csf_surface_tension_2d(const double *phi, int nx, int ny,
                             const double dx[2], double sigma,
                             double eps, double *Fx, double *Fy);

#ifdef __cplusplus
}
#endif

#endif /* INTERPHASE_NUMERICAL_H */

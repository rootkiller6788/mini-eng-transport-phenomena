/**
 * interphase_numerical.c
 * =======================
 * Numerical methods for interface tracking, reconstruction, and
 * boundary condition enforcement.
 *
 * Reference: Tryggvason, Scardovelli, Zaleski (2011) "Direct Numerical
 *   Simulations of Gas-Liquid Multiphase Flows" (Cambridge)
 *
 * Osher & Fedkiw (2003) "Level Set Methods and Dynamic Implicit Surfaces"
 *   (Springer Applied Mathematical Sciences, Vol. 153)
 *
 * Peskin (2002) "The Immersed Boundary Method" Acta Numerica, 11, 479-517.
 * Brackbill, Kothe, Zemach (1992) "A Continuum Method for Modeling Surface
 *   Tension" J. Comput. Phys., 100(2), 335-354.
 *
 * Knowledge coverage: L5 Engineering Methods, L8 Advanced Methods
 */

#include "interphase_numerical.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

/* ============================================================================
 * L5: Level-Set Signed Distance Function
 * ============================================================================ */

double level_set_signed_distance(const double *x, int dim,
                                  const double *interface_center,
                                  const double *normal,
                                  double offset) {
    (void)offset;  /* plane offset = n*center, computed internally */
    if (!x || !interface_center || !normal || dim < 1 || dim > 3)
        return 0.0;

    /* Signed distance from point x to plane n*x = offset */
    /* d = n*x - offset, where offset = n*interface_center */
    double dot = 0.0;
    double computed_offset = 0.0;
    for (int d = 0; d < dim; d++) {
        dot += normal[d] * x[d];
        computed_offset += normal[d] * interface_center[d];
    }

    /* signed distance = dot - computed_offset (plane through interface_center) */
    return dot - computed_offset;
}

double regularized_heaviside(double phi, double eps) {
    if (eps <= 0.0) return (phi >= 0.0) ? 1.0 : 0.0;

    if (phi < -eps) return 0.0;
    if (phi >  eps) return 1.0;

    /* Smooth transition: H(phi) = 0.5*(1 + phi/eps + sin(pi*phi/eps)/pi) */
    return 0.5 * (1.0 + phi / eps + sin(M_PI * phi / eps) / M_PI);
}

double regularized_delta(double phi, double eps) {
    if (eps <= 0.0) return 0.0;

    if (fabs(phi) > eps) return 0.0;

    /* delta(phi) = 0.5 * (1 + cos(pi*phi/eps)) / eps */
    return 0.5 * (1.0 + cos(M_PI * phi / eps)) / eps;
}

double level_set_interpolate(double phi, double eps,
                              double beta_A, double beta_B) {
    double H = regularized_heaviside(phi, eps);
    return beta_A + (beta_B - beta_A) * H;
}

void level_set_normal(const double *phi, int nx, int ny, int nz,
                      int i, int j, int k,
                      const double *dx, double normal[3]) {
    if (!phi || !dx || !normal || nx < 2 || ny < 2) {
        if (normal) { normal[0] = 0.0; normal[1] = 0.0; normal[2] = 1.0; }
        return;
    }

    int idx = k * nx * ny + j * nx + i;
    int im1 = (i > 0)      ? idx - 1     : idx;
    int ip1 = (i < nx-1)   ? idx + 1     : idx;
    int jm1 = (j > 0)      ? idx - nx    : idx;
    int jp1 = (j < ny-1)   ? idx + nx    : idx;

    /* Central differences for grad(phi) */
    double dphi_dx = (phi[ip1] - phi[im1]) / (2.0 * dx[0]);
    double dphi_dy = (phi[jp1] - phi[jm1]) / (2.0 * dx[1]);
    double dphi_dz = 0.0;

    if (nz > 1) {
        int km1 = (k > 0)    ? idx - nx*ny : idx;
        int kp1 = (k < nz-1) ? idx + nx*ny : idx;
        dphi_dz = (phi[kp1] - phi[km1]) / (2.0 * dx[2]);
    }

    /* Normalize */
    double mag = sqrt(dphi_dx*dphi_dx + dphi_dy*dphi_dy + dphi_dz*dphi_dz);
    if (mag < 1e-15) {
        normal[0] = 0.0; normal[1] = 0.0; normal[2] = 1.0;
        return;
    }

    normal[0] = dphi_dx / mag;
    normal[1] = dphi_dy / mag;
    normal[2] = dphi_dz / mag;
}

void level_set_curvature_2d(const double *phi, int nx, int ny,
                             int i, int j, const double dx[2],
                             double *curvature) {
    if (!phi || !dx || !curvature || i < 1 || i >= nx-1
        || j < 1 || j >= ny-1) {
        if (curvature) *curvature = 0.0;
        return;
    }

    int idx = j * nx + i;

    /* First derivatives (central) */
    double phi_x = (phi[idx+1] - phi[idx-1]) / (2.0 * dx[0]);
    double phi_y = (phi[idx+nx] - phi[idx-nx]) / (2.0 * dx[1]);

    /* Second derivatives */
    double phi_xx = (phi[idx+1] - 2.0*phi[idx] + phi[idx-1]) / (dx[0]*dx[0]);
    double phi_yy = (phi[idx+nx] - 2.0*phi[idx] + phi[idx-nx]) / (dx[1]*dx[1]);
    double phi_xy = (phi[idx+1+nx] - phi[idx-1+nx]
                     - phi[idx+1-nx] + phi[idx-1-nx])
                    / (4.0 * dx[0] * dx[1]);

    /* Mean curvature in 2D:
     * kappa = (phi_xx*phi_y^2 - 2*phi_x*phi_y*phi_xy + phi_yy*phi_x^2) /
     *         (phi_x^2 + phi_y^2)^(3/2)
     */
    double denom = pow(phi_x*phi_x + phi_y*phi_y, 1.5);
    if (denom < 1e-15) {
        *curvature = 0.0;
        return;
    }

    *curvature = (phi_xx * phi_y*phi_y - 2.0*phi_x*phi_y*phi_xy
                  + phi_yy * phi_x*phi_x) / denom;
}

/* ============================================================================
 * L5: Volume of Fluid (VOF) - PLIC Method
 * ============================================================================ */

void plic_vof_fraction_2d(const double normal[2], double alpha,
                           const double dx[2], const double x0[2],
                           double *F) {
    if (!normal || !dx || !F) {
        if (F) *F = 0.0;
        return;
    }

    double nx = normal[0];
    double ny = normal[1];

    /* Normalize */
    double nmag = sqrt(nx*nx + ny*ny);
    if (nmag < 1e-15) { *F = 0.5; return; }
    nx /= nmag;
    ny /= nmag;

    if (fabs(nx) < 1e-10 && fabs(ny) < 1e-10) { *F = 0.5; return; }

    /* Compute line offset at cell corners */
    double corner_vals[4];

    /* Bottom-left */
    corner_vals[0] = nx * x0[0] + ny * x0[1];
    /* Bottom-right */
    corner_vals[1] = nx * (x0[0] + dx[0]) + ny * x0[1];
    /* Top-right */
    corner_vals[2] = nx * (x0[0] + dx[0]) + ny * (x0[1] + dx[1]);
    /* Top-left */
    corner_vals[3] = nx * x0[0] + ny * (x0[1] + dx[1]);

    /* Count how many corners have value > alpha (in phase B) */
    int count_B = 0;
    for (int c = 0; c < 4; c++) {
        if (corner_vals[c] > alpha) count_B++;
    }

    /* For simple cases */
    if (count_B == 0) { *F = 0.0; return; }
    if (count_B == 4) { *F = 1.0; return; }

    /* For the general case, compute area fraction geometrically.
     * Here we use a simplified version based on the PLIC formula.
     * The general solution involves polygon clipping.
     * We approximate using the normalized alpha. */
    double cell_center_val = nx * (x0[0] + 0.5*dx[0])
                              + ny * (x0[1] + 0.5*dx[1]);
    double alpha_norm = (alpha - cell_center_val)
                        / sqrt(dx[0]*dx[0] + dx[1]*dx[1]);
    /* Simple linear model: F ~ 0.5 - alpha_norm */
    *F = 0.5 - alpha_norm;
    if (*F < 0.0) *F = 0.0;
    if (*F > 1.0) *F = 1.0;
}

void plic_inverse_vof_2d(const double normal[2], double F,
                          const double dx[2], double *alpha) {
    if (!normal || !dx || !alpha) {
        if (alpha) *alpha = 0.0;
        return;
    }

    /* Solve PLIC(F) = alpha_target using bisection */
    double alpha_lo = -10.0, alpha_hi = 10.0;
    double x0[2] = {0.0, 0.0};
    double F_lo, F_hi, F_mid;

    plic_vof_fraction_2d(normal, alpha_lo, dx, x0, &F_lo);
    plic_vof_fraction_2d(normal, alpha_hi, dx, x0, &F_hi);

    /* Bisection for 200 iterations (precision ~ 2^(-200)) */
    for (int iter = 0; iter < 200; iter++) {
        double alpha_mid = 0.5 * (alpha_lo + alpha_hi);
        plic_vof_fraction_2d(normal, alpha_mid, dx, x0, &F_mid);

        if (F_mid < F) {
            alpha_lo = alpha_mid;
        } else {
            alpha_hi = alpha_mid;
        }

        if (fabs(alpha_hi - alpha_lo) < 1e-15) break;
    }

    *alpha = 0.5 * (alpha_lo + alpha_hi);
}

void vof_flux_2d(const double normal[2], double alpha,
                  const double dx[2], double dt,
                  const double u_face[2], int face_normal,
                  double *flux) {
    if (!normal || !dx || !u_face || !flux) {
        if (flux) *flux = 0.0;
        return;
    }

    /* Compute the VOF flux across a cell face.
     * Flux = (u_face * dt) * (face length) * F_face * face_normal_sign
     * where F_face is the VOF fraction of the flux polygon.
     */
    double u_n = u_face[face_normal];
    double face_length = dx[1 - face_normal];
    double swept_volume = fabs(u_n) * dt * face_length;

    /* Approximate: flux = F * swept_volume (upwinded F would be better) */
    double F_cell;
    double x0[2] = {0.0, 0.0};
    plic_vof_fraction_2d(normal, alpha, dx, x0, &F_cell);

    *flux = F_cell * swept_volume;
}

void plic_interface_area_2d(const double normal[2], double alpha,
                             const double dx[2], double *area) {
    if (!normal || !dx || !area) {
        if (area) *area = 0.0;
        return;
    }

    double nx = normal[0], ny = normal[1];
    double nmag = sqrt(nx*nx + ny*ny);
    if (nmag < 1e-15) { *area = 0.0; return; }
    nx /= nmag; ny /= nmag;

    /* Interface length intersecting the cell.
     * For a line nx*x + ny*y = alpha in rectangle [0,dx[0]] x [0,dx[1]].
     * Simplified: area ~ sqrt(dx^2 + dy^2) * F*(1-F) (empirical max at F=0.5)
     */
    double F;
    double x0[2] = {0.0, 0.0};
    plic_vof_fraction_2d(normal, alpha, dx, x0, &F);

    /* Typical maximum chord length ~ sqrt(dx_min^2 + dx_max^2) */
    double diag = sqrt(dx[0]*dx[0] + dx[1]*dx[1]);
    *area = diag * 4.0 * F * (1.0 - F);
}

/* ============================================================================
 * L8: Ghost Fluid Method
 *
 * Reference: Fedkiw, Aslam, Merriman, Osher (1999)
 *   "A Non-oscillatory Eulerian Approach to Interfaces in Multimaterial
 *    Flows (the Ghost Fluid Method)" J. Comput. Phys., 152(2), 457-492.
 * ============================================================================ */

void gfm_dirichlet(double phi_real, double phi_ghost,
                    double value_real, double value_interface,
                    double *value_ghost) {
    if (!value_ghost) return;

    /* Dirichlet GFM: u_I = specified
     * Ghost value = 2*u_I - u_real (linear extrapolation)
     * More precisely: u_ghost = u_I + (u_I - u_real) * (phi_ghost/phi_real)
     * For simple case with |phi_ghost| ~ |phi_real| (uniform grid):
     *   u_ghost = 2*u_I - u_real
     */
    if (fabs(phi_real) < 1e-15) {
        *value_ghost = value_interface;
        return;
    }

    /* Linear extrapolation: u(phi)=u_I + m*phi, m=(u_real-u_I)/phi_real */
    double theta = phi_ghost / phi_real;
    *value_ghost = value_interface + (value_real - value_interface) * theta;
}

void gfm_neumann(double phi_real, double phi_ghost,
                  double value_real, double gradient_interface,
                  double *value_ghost) {
    if (!value_ghost) return;

    /* Neumann GFM: du/dn = specified at interface
     * u_ghost = u_real + (phi_ghost - phi_real) * gradient_interface
     * Actually: u_ghost = u_I + phi_ghost * (du/dn)_I
     *           u_I = u_real - phi_real * (du/dn)_I
     * => u_ghost = u_real - phi_real*(du/dn)_I + phi_ghost*(du/dn)_I
     *            = u_real + (phi_ghost - phi_real) * (du/dn)_I
     */
    double du_dn = gradient_interface;
    *value_ghost = value_real + (phi_ghost - phi_real) * du_dn;
}

void gfm_robin(double phi_real, double phi_ghost,
                double value_real,
                double a_coeff, double b_coeff, double c_coeff,
                double *value_ghost) {
    if (!value_ghost) return;

    /* Robin GFM: a*u + b*du/dn = c at interface
     * Interface value u_I satisfies: a*u_I + b*(u_I - u_R)/(phi_I - phi_R) = c
     * where phi_I = 0, phi_R = phi_real
     * => a*u_I + b*(u_I - u_R)/(-phi_R) = c
     * => a*u_I - b*(u_I - u_R)/phi_R = c
     * => u_I*(a - b/phi_R) = c - b*u_R/phi_R
     * => u_I = (c - b*u_R/phi_R) / (a - b/phi_R)
     *
     * Then: u_ghost = u_I + (phi_ghost/phi_real)*(u_I - u_real)
     */

    if (fabs(b_coeff) < 1e-15) {
        /* Degenerate to Dirichlet: a*u = c => u_I = c/a */
        double u_I = (fabs(a_coeff) > 1e-15) ? (c_coeff / a_coeff) : 0.0;
        gfm_dirichlet(phi_real, phi_ghost, value_real, u_I, value_ghost);
        return;
    }

    if (fabs(phi_real) < 1e-15) {
        /* Real node at interface */
        double u_I = value_real;
        double du_dn = (c_coeff - a_coeff * u_I) / b_coeff;
        *value_ghost = u_I + phi_ghost * du_dn;
        return;
    }

    double inv_phi = 1.0 / phi_real;
    double denom = a_coeff - b_coeff * inv_phi;

    if (fabs(denom) < 1e-15) {
        *value_ghost = value_real;
        return;
    }

    double u_I = (c_coeff - b_coeff * value_real * inv_phi) / denom;

    /* Extrapolate to ghost */
    *value_ghost = u_I + (phi_ghost / phi_real) * (u_I - value_real);
}

/* ============================================================================
 * L8: Immersed Boundary Method
 *
 * Reference: Roma, Peskin, Berger (1999)
 *   "An Adaptive Version of the Immersed Boundary Method"
 *   J. Comput. Phys., 153(2), 509-534.
 * ============================================================================ */

double ibm_delta_roma(double r) {
    double abs_r = fabs(r);

    if (abs_r <= 0.5) {
        double val = 1.0 - 3.0 * abs_r * abs_r;
        if (val < 0.0) val = 0.0;
        return (1.0 / 3.0) * (1.0 + sqrt(val));
    } else if (abs_r <= 1.5) {
        double val = -3.0 * (1.0 - abs_r) * (1.0 - abs_r) + 1.0;
        if (val < 0.0) val = 0.0;
        return (1.0 / 6.0) * (5.0 - 3.0 * abs_r - sqrt(val));
    } else {
        return 0.0;
    }
}

void ibm_spread_force(const double *x_grid, int nx, int ny,
                      const double dx[2],
                      const double *X_lagrangian,
                      const double *F_lagrangian,
                      int n_markers, double ds,
                      double *F_euler) {
    if (!x_grid || !X_lagrangian || !F_lagrangian || !F_euler
        || n_markers <= 0 || nx <= 0 || ny <= 0) {
        return;
    }

    /* Zero output */
    memset(F_euler, 0, sizeof(double) * nx * ny * 2);

    /* Spread: for each Lagrangian marker, distribute force to nearby
     * Eulerian grid points using the regularized delta function. */
    for (int m = 0; m < n_markers; m++) {
        double X = X_lagrangian[2*m];
        double Y = X_lagrangian[2*m+1];
        double Fx = F_lagrangian[2*m];
        double Fy = F_lagrangian[2*m+1];

        /* Find the bounding Eulerian cell indices */
        int i0 = (int) floor((X - x_grid[0]) / dx[0]);
        int j0 = (int) floor((Y - x_grid[1]) / dx[1]);

        /* Spread to 4x4 stencil around (i0, j0) */
        for (int di = -1; di <= 2; di++) {
            for (int dj = -1; dj <= 2; dj++) {
                int i = i0 + di;
                int j = j0 + dj;
                if (i < 0 || i >= nx || j < 0 || j >= ny) continue;

                double x_node = x_grid[2*(j*nx + i)];
                double y_node = x_grid[2*(j*nx + i) + 1];

                double rx = (x_node - X) / dx[0];
                double ry = (y_node - Y) / dx[1];

                double delta_h = ibm_delta_roma(rx) * ibm_delta_roma(ry);

                int idx = (j * nx + i);
                F_euler[2*idx]     += Fx * delta_h * ds;
                F_euler[2*idx + 1] += Fy * delta_h * ds;
            }
        }
    }
}

void ibm_interpolate_velocity(const double *u_euler, int nx, int ny,
                               const double dx[2],
                               const double *X_lagrangian,
                               int n_markers,
                               double *U_lagrangian) {
    if (!u_euler || !X_lagrangian || !U_lagrangian
        || n_markers <= 0 || nx <= 0 || ny <= 0) {
        return;
    }

    memset(U_lagrangian, 0, sizeof(double) * n_markers * 2);

    /* x_grid is implicit: i*dx, j*dy */
    double x0 = 0.0, y0 = 0.0;  /* origin */

    /* Interpolate: for each Lagrangian marker, gather velocity from
     * nearby Eulerian grid points using the same delta function. */
    for (int m = 0; m < n_markers; m++) {
        double X = X_lagrangian[2*m];
        double Y = X_lagrangian[2*m+1];

        int i0 = (int) floor((X - x0) / dx[0]);
        int j0 = (int) floor((Y - y0) / dx[1]);

        double Ux = 0.0, Uy = 0.0;
        double weight_sum = 0.0;

        for (int di = -1; di <= 2; di++) {
            for (int dj = -1; dj <= 2; dj++) {
                int i = i0 + di;
                int j = j0 + dj;
                if (i < 0 || i >= nx || j < 0 || j >= ny) continue;

                double x_node = x0 + (i + 0.5) * dx[0];
                double y_node = y0 + (j + 0.5) * dx[1];

                double rx = (x_node - X) / dx[0];
                double ry = (y_node - Y) / dx[1];

                double delta_h = ibm_delta_roma(rx) * ibm_delta_roma(ry);

                int idx = (j * nx + i);
                Ux += u_euler[2*idx]     * delta_h;
                Uy += u_euler[2*idx + 1] * delta_h;
                weight_sum += delta_h;
            }
        }

        /* Normalize: discrete sum approximates integral */
        if (weight_sum > 0.0) {
            U_lagrangian[2*m]     = Ux / weight_sum;
            U_lagrangian[2*m + 1] = Uy / weight_sum;
        }
    }
}

/* ============================================================================
 * L5: Level-Set Advection and Reinitialization
 * ============================================================================ */

void level_set_advect_upwind(const double *phi_n, double *phi_np1,
                              const double *u, const double *v,
                              int nx, int ny,
                              const double dx[2], double dt) {
    if (!phi_n || !phi_np1 || !u || !v || nx < 2 || ny < 2) return;

    double dxi = 1.0 / dx[0];
    double dyi = 1.0 / dx[1];

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int idx = j * nx + i;

            /* First-order upwind for d(phi)/dt + u*d(phi)/dx + v*d(phi)/dy = 0 */

            /* x-derivative */
            double dphi_dx;
            if (u[idx] > 0.0) {
                /* Upwind from left */
                int im1 = (i > 0) ? (idx - 1) : idx;
                dphi_dx = (phi_n[idx] - phi_n[im1]) * dxi;
            } else {
                /* Upwind from right */
                int ip1 = (i < nx-1) ? (idx + 1) : idx;
                dphi_dx = (phi_n[ip1] - phi_n[idx]) * dxi;
            }

            /* y-derivative */
            double dphi_dy;
            if (v[idx] > 0.0) {
                /* Upwind from bottom */
                int jm1 = (j > 0) ? (idx - nx) : idx;
                dphi_dy = (phi_n[idx] - phi_n[jm1]) * dyi;
            } else {
                /* Upwind from top */
                int jp1 = (j < ny-1) ? (idx + nx) : idx;
                dphi_dy = (phi_n[jp1] - phi_n[idx]) * dyi;
            }

            phi_np1[idx] = phi_n[idx]
                           - dt * (u[idx] * dphi_dx + v[idx] * dphi_dy);
        }
    }
}

void level_set_reinitialize(double *phi, int nx, int ny,
                             const double dx[2],
                             double dtau, int n_iter) {
    if (!phi || nx < 2 || ny < 2) return;

    double *phi_temp = (double*) malloc(nx * ny * sizeof(double));
    if (!phi_temp) return;

    for (int iter = 0; iter < n_iter; iter++) {
        memcpy(phi_temp, phi, nx * ny * sizeof(double));

        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                int idx = j * nx + i;
                double phi0 = phi[idx];

                /* sign function (smoothed) */
                double s = phi0 / sqrt(phi0*phi0 + dx[0]*dx[0]);

                /* Godunov scheme for |grad(phi)| */
                double a, b, c, d;

                /* x-direction */
                int im1 = (i > 0) ? idx-1 : idx;
                int ip1 = (i < nx-1) ? idx+1 : idx;
                double Dx_minus = (phi_temp[idx] - phi_temp[im1]) / dx[0];
                double Dx_plus  = (phi_temp[ip1] - phi_temp[idx]) / dx[0];

                if (s > 0) {
                    a = fmax(Dx_minus, 0.0); a = a*a;
                    b = fmin(Dx_plus, 0.0);  b = b*b;
                } else {
                    a = fmin(Dx_minus, 0.0); a = a*a;
                    b = fmax(Dx_plus, 0.0);  b = b*b;
                }

                /* y-direction */
                int jm1 = (j > 0) ? idx-nx : idx;
                int jp1 = (j < ny-1) ? idx+nx : idx;
                double Dy_minus = (phi_temp[idx] - phi_temp[jm1]) / dx[1];
                double Dy_plus  = (phi_temp[jp1] - phi_temp[idx]) / dx[1];

                if (s > 0) {
                    c = fmax(Dy_minus, 0.0); c = c*c;
                    d = fmin(Dy_plus, 0.0);  d = d*d;
                } else {
                    c = fmin(Dy_minus, 0.0); c = c*c;
                    d = fmax(Dy_plus, 0.0);  d = d*d;
                }

                double grad_mag = sqrt(fmax(a,b) + fmax(c,d));

                /* Update: phi^(n+1) = phi^n - dtau * s * (grad_mag - 1) */
                phi[idx] = phi_temp[idx]
                           - dtau * s * (grad_mag - 1.0);
            }
        }
    }

    free(phi_temp);
}

/* ============================================================================
 * L5: Continuum Surface Force (CSF) for Surface Tension
 *
 * Reference: Brackbill, Kothe, Zemach (1992)
 *   "A Continuum Method for Modeling Surface Tension"
 *   J. Comput. Phys., 100(2), 335-354.
 * ============================================================================ */

void csf_surface_tension_2d(const double *phi, int nx, int ny,
                             const double dx[2], double sigma,
                             double eps, double *Fx, double *Fy) {
    if (!phi || !Fx || !Fy || nx < 3 || ny < 3) return;

    /* Zero output */
    memset(Fx, 0, nx * ny * sizeof(double));
    memset(Fy, 0, nx * ny * sizeof(double));

    for (int j = 1; j < ny - 1; j++) {
        for (int i = 1; i < nx - 1; i++) {
            int idx = j * nx + i;

            /* Curvature */
            double kappa;
            level_set_curvature_2d(phi, nx, ny, i, j, dx, &kappa);

            /* Normal */
            double normal[3];
            level_set_normal(phi, nx, ny, 1, i, j, 0, dx, normal);

            /* Smoothed delta function */
            double delta_phi = regularized_delta(phi[idx], eps);

            /* F = sigma * kappa * delta(phi) * n */
            Fx[idx] = sigma * kappa * delta_phi * normal[0];
            Fy[idx] = sigma * kappa * delta_phi * normal[1];
        }
    }
}

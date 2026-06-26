#include "cfd_convection_diffusion.h"
#include "cfd_discretization.h"
#include "cfd_solvers.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ===== 1D Steady Convection-Diffusion Solver ===== */

int cfd_solve_convdiff_1d(const CfdGrid1D *grid, const CfdTransportEquation *eqn,
                          double velocity, double phi_left, double phi_right,
                          CfdConvectiveScheme scheme, double *phi)
{
    if (!grid || !eqn || !phi) return -1;
    int nx = grid->nx;
    CfdMatrix2D *mat = cfd_matrix2d_create(nx, 1);
    if (!mat) return -1;

    cfd_assemble_convdiff_1d(grid, eqn, velocity, phi_left, phi_right, scheme, mat);

    CfdSolverState state;
    cfd_solver_state_init(&state, CFD_SOLVER_TDMA, 1e-10, 1);
    int ret = cfd_solve_tdma(mat, phi, nx);
    cfd_matrix2d_destroy(mat);
    return ret;
}

/* ===== 1D Exact Solution ===== */

double cfd_convdiff_1d_exact(double x, double L, double Pe,
                              double phi_0, double phi_L)
{
    if (L <= 0.0) return phi_0;
    if (fabs(Pe) < 1e-10) {
        return phi_0 + (phi_L - phi_0) * x / L;
    }
    double exp_Pe = exp(Pe);
    return phi_0 + (phi_L - phi_0) * (exp(Pe * x / L) - 1.0) / (exp_Pe - 1.0);
}

/* ===== 1D Unsteady Convection-Diffusion (Explicit Euler) ===== */

int cfd_solve_convdiff_1d_explicit(const CfdGrid1D *grid,
                                    const CfdTransportEquation *eqn,
                                    double velocity,
                                    double phi_left, double phi_right,
                                    CfdConvectiveScheme space_scheme,
                                    double dt, int n_steps, double *phi)
{
    if (!grid || !eqn || !phi || dt <= 0.0) return -1;
    int nx = grid->nx;
    double dx = grid->dx;
    double Gamma = eqn->diffusivity;
    double u = velocity;

    double *phi_new = malloc((size_t)nx * sizeof(double));
    if (!phi_new) return -1;

    for (int step = 0; step < n_steps; step++) {
        for (int i = 0; i < nx; i++) {
            double conv = 0.0, diff = 0.0;
            /* Convection: upwind */
            if (u > 0.0) {
                double phi_w = (i > 0) ? phi[i-1] : phi_left;
                conv = -u * (phi[i] - phi_w) / dx;
            } else if (u < 0.0) {
                double phi_e = (i < nx-1) ? phi[i+1] : phi_right;
                conv = -u * (phi_e - phi[i]) / dx;
            }
            /* Diffusion: central */
            double phi_w = (i > 0) ? phi[i-1] : phi_left;
            double phi_e = (i < nx-1) ? phi[i+1] : phi_right;
            diff = Gamma * (phi_e - 2.0*phi[i] + phi_w) / (dx * dx);

            phi_new[i] = phi[i] + dt * (conv + diff + eqn->source_constant);
        }
        memcpy(phi, phi_new, (size_t)nx * sizeof(double));
    }
    free(phi_new);
    return n_steps;
}

/* ===== 1D Unsteady Convection-Diffusion (Implicit Euler) ===== */

int cfd_solve_convdiff_1d_implicit(const CfdGrid1D *grid,
                                    const CfdTransportEquation *eqn,
                                    double velocity,
                                    double phi_left, double phi_right,
                                    CfdConvectiveScheme space_scheme,
                                    double dt, int n_steps, double *phi)
{
    if (!grid || !eqn || !phi || dt <= 0.0) return -1;
    int nx = grid->nx;
    double dx = grid->dx;
    double Gamma = eqn->diffusivity;
    double rho = eqn->density;
    double u = velocity;
    double F = rho * u;
    double D = Gamma / dx;

    for (int step = 0; step < n_steps; step++) {
        CfdMatrix2D *mat = cfd_matrix2d_create(nx, 1);
        if (!mat) return -1;

        for (int i = 0; i < nx; i++) {
            double Pe = (D > 1e-15) ? F * dx / Gamma : 0.0;
            double aWc = 0, aEc = 0, aPc = 0;
            if (space_scheme == CFD_SCHEME_POWERLAW)
                cfd_powerlaw_1d(Pe, F, D, &aWc, &aEc, &aPc);
            else {
                aWc = D + CFD_MAX(F, 0.0);
                aEc = D + CFD_MAX(-F, 0.0);
                aPc = aWc + aEc;
            }

            /* Implicit Euler: add dt/dx contributions */
            mat->aW[i] = -dt * aWc;
            mat->aE[i] = -dt * aEc;
            if (i > 0) mat->aW[i] = -dt * aWc;
            if (i < nx-1) mat->aE[i] = -dt * aEc;
            else mat->aE[i] = 0.0;
            if (i == 0) mat->aW[i] = 0.0;

            mat->aP[i] = 1.0 + dt * aPc;
            mat->b[i] = phi[i] + dt * eqn->source_constant;
        }

        /* Dirichlet BC */
        mat->aP[0] += dt * (2.0*D + CFD_MAX(F,0.0));
        mat->b[0] += dt * (2.0*D + CFD_MAX(F,0.0)) * phi_left;
        int last = nx - 1;
        mat->aP[last] += dt * (2.0*D + CFD_MAX(-F,0.0));
        mat->b[last] += dt * (2.0*D + CFD_MAX(-F,0.0)) * phi_right;

        cfd_solve_tdma(mat, phi, nx);
        cfd_matrix2d_destroy(mat);
    }
    return n_steps;
}


/* ===== 2D Steady Convection-Diffusion ===== */

int cfd_solve_convdiff_2d(const CfdGrid2D *grid, const CfdTransportEquation *eqn,
                          const CfdVectorField2D *vel, const CfdBCSet2D *bc,
                          CfdConvectiveScheme scheme, CfdScalarField2D *phi)
{
    if (!grid || !eqn || !vel || !bc || !phi) return -1;
    int nx = grid->nx, ny = grid->ny;
    CfdMatrix2D *mat = cfd_matrix2d_create(nx, ny);
    if (!mat) return -1;

    /* Build scalar transport matrix similar to u-momentum */
    cfd_matrix2d_zero(mat);
    double rho = eqn->density, Gamma = eqn->diffusivity;
    double dx = grid->dx, dy = grid->dy;

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int idx = CFD_IDX2(i, j, nx);
            double DW = Gamma * dy / dx, DE = Gamma * dy / dx;
            double DS = Gamma * dx / dy, DN = Gamma * dx / dy;
            double FW = 0, FE = 0, FS = 0, FN = 0;

            if (i > 0) FW = rho * vel->u[CFD_IDX2(i-1,j,nx)] * dy;
            if (i < nx-1) FE = rho * vel->u[CFD_IDX2(i+1,j,nx)] * dy;
            if (j > 0) FS = rho * vel->v[CFD_IDX2(i,j-1,nx)] * dx;
            if (j < ny-1) FN = rho * vel->v[CFD_IDX2(i,j+1,nx)] * dx;

            double aW=0, aE=0, aS=0, aN=0;
            if (i > 0)      aW = DW + CFD_MAX(FW, 0.0);
            if (i < nx - 1) aE = DE + CFD_MAX(-FE, 0.0);
            if (j > 0)      aS = DS + CFD_MAX(FS, 0.0);
            if (j < ny - 1) aN = DN + CFD_MAX(-FN, 0.0);
            double aP = aW + aE + aS + aN + (FE - FW + FN - FS);

            mat->aW[idx] = aW; mat->aE[idx] = aE;
            mat->aS[idx] = aS; mat->aN[idx] = aN;
            mat->aP[idx] = aP;
            mat->b[idx] = eqn->source_constant * dx * dy;
        }
    }

    /* Note: BCs are implicitly handled through the face coefficient
     * computation at domain boundaries. For Dirichlet BCs, the
     * boundary values are incorporated into the source term b[idx].
     * For more general BC types, explicit boundary row modification
     * would be added here. */

    CfdSolverState state;
    cfd_solver_state_init(&state, CFD_SOLVER_GAUSS_SEIDEL, 1e-6, 5000);
    int ret = cfd_solve_gauss_seidel(mat, phi->data, &state);
    cfd_matrix2d_destroy(mat);
    return ret;
}

/* ===== 2D Unsteady Convection-Diffusion (Explicit) ===== */

int cfd_solve_convdiff_2d_explicit(const CfdGrid2D *grid,
                                    const CfdTransportEquation *eqn,
                                    const CfdVectorField2D *vel,
                                    const CfdBCSet2D *bc,
                                    CfdConvectiveScheme scheme,
                                    double dt, int n_steps,
                                    CfdScalarField2D *phi)
{
    if (!grid || !eqn || !vel || !bc || !phi) return -1;
    int nx = grid->nx, ny = grid->ny, nc = nx * ny;
    double dx = grid->dx, dy = grid->dy;
    double Gamma = eqn->diffusivity;

    double *phi_new = malloc((size_t)nc * sizeof(double));
    if (!phi_new) return -1;

    for (int step = 0; step < n_steps; step++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                int idx = CFD_IDX2(i, j, nx);
                double phi_P = phi->data[idx];
                double conv_x = 0, conv_y = 0, diff_x = 0, diff_y = 0;

                /* Convection x */
                double u_P = vel->u[idx];
                if (u_P > 0 && i > 0)
                    conv_x = -u_P * (phi_P - phi->data[CFD_IDX2(i-1,j,nx)]) / dx;
                else if (u_P < 0 && i < nx-1)
                    conv_x = -u_P * (phi->data[CFD_IDX2(i+1,j,nx)] - phi_P) / dx;

                /* Convection y */
                double v_P = vel->v[idx];
                if (v_P > 0 && j > 0)
                    conv_y = -v_P * (phi_P - phi->data[CFD_IDX2(i,j-1,nx)]) / dy;
                else if (v_P < 0 && j < ny-1)
                    conv_y = -v_P * (phi->data[CFD_IDX2(i,j+1,nx)] - phi_P) / dy;

                /* Diffusion */
                double phi_w = (i>0)?phi->data[CFD_IDX2(i-1,j,nx)]:bc->west.value;
                double phi_e = (i<nx-1)?phi->data[CFD_IDX2(i+1,j,nx)]:bc->east.value;
                double phi_s = (j>0)?phi->data[CFD_IDX2(i,j-1,nx)]:bc->south.value;
                double phi_n = (j<ny-1)?phi->data[CFD_IDX2(i,j+1,nx)]:bc->north.value;
                diff_x = Gamma * (phi_e - 2.0*phi_P + phi_w) / (dx*dx);
                diff_y = Gamma * (phi_n - 2.0*phi_P + phi_s) / (dy*dy);

                phi_new[idx] = phi_P + dt * (conv_x + conv_y + diff_x + diff_y + eqn->source_constant);
            }
        }
        memcpy(phi->data, phi_new, (size_t)nc * sizeof(double));
    }
    free(phi_new);
    return n_steps;
}

/* ===== 2D Poisson Solver ===== */

int cfd_solve_poisson_2d(const CfdGrid2D *grid, const CfdScalarField2D *source,
                          const CfdBCSet2D *bc, CfdScalarField2D *phi)
{
    if (!grid || !bc || !phi) return -1;
    int nx = grid->nx, ny = grid->ny;
    CfdMatrix2D *mat = cfd_matrix2d_create(nx, ny);
    if (!mat) return -1;

    cfd_assemble_poisson_2d(grid, source, bc, mat);

    CfdSolverState state;
    cfd_solver_state_init(&state, CFD_SOLVER_GAUSS_SEIDEL, 1e-6, 10000);
    int ret = cfd_solve_gauss_seidel(mat, phi->data, &state);

    cfd_matrix2d_destroy(mat);
    return ret;
}


/* ===== SIMPLE Algorithm (Patankar & Spalding 1972) ===== */

int cfd_simple_iteration(const CfdGrid2D *grid, CfdVectorField2D *vel,
                          const CfdFluidProperty *fp, CfdSimpleState *simple)
{
    if (!grid || !vel || !fp || !simple) return -1;
    int nx = grid->nx, ny = grid->ny, nc = nx * ny;
    double *p = vel->p;
    double *u = vel->u;
    double *v = vel->v;
    double rho = fp->density;
    double dx = grid->dx, dy = grid->dy;

    /* Step 1: Solve u-momentum */
    CfdMatrix2D *u_mat = cfd_matrix2d_create(nx, ny);
    CfdMatrix2D *v_mat = cfd_matrix2d_create(nx, ny);
    if (!u_mat || !v_mat) { cfd_matrix2d_destroy(u_mat); cfd_matrix2d_destroy(v_mat); return -1; }

    cfd_assemble_u_momentum_2d(grid, vel, fp, CFD_SCHEME_POWERLAW, u_mat);
    cfd_assemble_v_momentum_2d(grid, vel, fp, CFD_SCHEME_POWERLAW, v_mat);

    /* Add pressure gradient to source terms */
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int idx = CFD_IDX2(i, j, nx);
            double dpdx = 0, dpdy = 0;
            if (i > 0 && i < nx-1) dpdx = (p[CFD_IDX2(i+1,j,nx)] - p[CFD_IDX2(i-1,j,nx)]) / (2.0*dx);
            if (j > 0 && j < ny-1) dpdy = (p[CFD_IDX2(i,j+1,nx)] - p[CFD_IDX2(i,j-1,nx)]) / (2.0*dy);
            u_mat->b[idx] -= dpdx * dx * dy;
            v_mat->b[idx] -= dpdy * dx * dy;
        }
    }

    /* Under-relax and solve momentum */
    CfdSolverState u_state;
    cfd_solver_state_init(&u_state, CFD_SOLVER_GAUSS_SEIDEL, simple->momentum_tol, simple->max_inner_iter);
    u_state.under_relaxation = simple->urf_u;
    cfd_solve_gauss_seidel(u_mat, u, &u_state);

    CfdSolverState v_state;
    cfd_solver_state_init(&v_state, CFD_SOLVER_GAUSS_SEIDEL, simple->momentum_tol, simple->max_inner_iter);
    v_state.under_relaxation = simple->urf_v;
    cfd_solve_gauss_seidel(v_mat, v, &v_state);

    /* Step 2: Pressure correction equation */
    CfdMatrix2D *p_mat = cfd_matrix2d_create(nx, ny);
    if (!p_mat) { cfd_matrix2d_destroy(u_mat); cfd_matrix2d_destroy(v_mat); return -1; }

    cfd_matrix2d_zero(p_mat);
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int idx = CFD_IDX2(i, j, nx);
            double aW = (i > 0)      ? rho * dy * dy : 0.0;
            double aE = (i < nx - 1) ? rho * dy * dy : 0.0;
            double aS = (j > 0)      ? rho * dx * dx : 0.0;
            double aN = (j < ny - 1) ? rho * dx * dx : 0.0;
            p_mat->aW[idx] = -aW;
            p_mat->aE[idx] = -aE;
            p_mat->aS[idx] = -aS;
            p_mat->aN[idx] = -aN;
            p_mat->aP[idx] = aW + aE + aS + aN;

            /* Source: mass imbalance (continuity) */
            double divU = cfd_divergence_2d(vel, i, j, grid);
            p_mat->b[idx] = -rho * divU * dx * dy;
        }
    }

    double *p_corr = calloc((size_t)nc, sizeof(double));
    if (!p_corr) { cfd_matrix2d_destroy(p_mat); cfd_matrix2d_destroy(u_mat); cfd_matrix2d_destroy(v_mat); return -1; }

    CfdSolverState p_state;
    cfd_solver_state_init(&p_state, CFD_SOLVER_GAUSS_SEIDEL, simple->pressure_tol, simple->max_inner_iter);
    cfd_solve_gauss_seidel(p_mat, p_corr, &p_state);

    /* Step 3: Correct velocities and pressure */
    double mass_res = 0.0;
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int idx = CFD_IDX2(i, j, nx);
            double dpc_dx = 0, dpc_dy = 0;
            if (i > 0 && i < nx-1) dpc_dx = (p_corr[CFD_IDX2(i+1,j,nx)] - p_corr[CFD_IDX2(i-1,j,nx)]) / (2.0*dx);
            if (j > 0 && j < ny-1) dpc_dy = (p_corr[CFD_IDX2(i,j+1,nx)] - p_corr[CFD_IDX2(i,j-1,nx)]) / (2.0*dy);

            u[idx] -= simple->urf_u * dx * dpc_dx / rho;
            v[idx] -= simple->urf_v * dy * dpc_dy / rho;
            p[idx] += simple->urf_p * p_corr[idx];

            double div = cfd_divergence_2d(vel, i, j, grid);
            mass_res += fabs(div) * dx * dy;
        }
    }

    simple->mass_residual = mass_res;
    simple->n_outer++;

    free(p_corr);
    cfd_matrix2d_destroy(p_mat);
    cfd_matrix2d_destroy(u_mat);
    cfd_matrix2d_destroy(v_mat);
    return 0;
}

int cfd_simple_solve(const CfdGrid2D *grid, CfdVectorField2D *vel,
                      const CfdFluidProperty *fp, CfdSimpleState *simple)
{
    if (!grid || !vel || !fp || !simple) return -1;
    simple->n_outer = 0;
    for (int outer = 0; outer < simple->max_outer_iter; outer++) {
        int ret = cfd_simple_iteration(grid, vel, fp, simple);
        if (ret < 0) return ret;
        if (simple->mass_residual < simple->pressure_tol) break;
    }
    return simple->n_outer;
}

/* ===== Lid-Driven Cavity ===== */

int cfd_solve_lid_driven_cavity(const CfdGrid2D *grid, CfdVectorField2D *vel,
                                 const CfdFluidProperty *fp, double U_lid,
                                 CfdSimpleState *state)
{
    if (!grid || !vel || !fp || !state) return -1;
    int nx = grid->nx, ny = grid->ny, nc = nx * ny;

    /* Initialize: zero velocity, zero pressure */
    for (int i = 0; i < nc; i++) { vel->u[i] = 0.0; vel->v[i] = 0.0; vel->p[i] = 0.0; }

    /* Set lid velocity on top row */
    for (int i = 0; i < nx; i++) {
        vel->u[CFD_IDX2(i, ny-1, nx)] = U_lid;
    }

    /* Configure SIMPLE */
    state->max_outer_iter = CFD_PRESSURE_CORR_ITER * 100;
    state->max_inner_iter = 50;
    state->momentum_tol = 1e-4;
    state->pressure_tol = 1e-6;
    state->urf_u = CFD_URF_U;
    state->urf_v = CFD_URF_U;
    state->urf_p = CFD_URF_P;
    state->n_outer = 0;

    return cfd_simple_solve(grid, vel, fp, state);
}

/* ===== Additional Utilities ===== */

/* Compute L2 error between numerical and exact solution for 1D. */
double cfd_l2_error_1d(const double *num, const CfdGrid1D *grid,
                        double (*exact)(double, void *), void *params)
{
    if (!num || !grid || !exact) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < grid->nx; i++) {
        double diff = num[i] - exact(grid->xc[i], params);
        sum += diff * diff;
    }
    return sqrt(sum / (double)grid->nx);
}

/* Richardson extrapolation for order of accuracy estimation. */
double cfd_richardson_extrapolation(double f_h, double f_h2, double h_ratio)
{
    /* f_exact = f_h + (f_h - f_h2) / (h_ratio^p - 1) for order p */
    if (h_ratio <= 1.0) return f_h;
    double p = 2.0; /* assumed second order */
    return f_h + (f_h - f_h2) / (pow(h_ratio, p) - 1.0);
}

/* Estimate order of accuracy from two grid solutions. */
double cfd_order_of_accuracy(double err_h, double err_h2, double h_ratio)
{
    if (err_h <= 0.0 || err_h2 <= 0.0 || h_ratio <= 1.0) return 0.0;
    return log(err_h2 / err_h) / log(h_ratio);
}

/* Compute the grid convergence index (GCI) (Roache 1994). */
double cfd_grid_convergence_index(double err_fine, double err_coarse,
                                    double h_ratio, double p, double Fs)
{
    /* GCI = Fs * |err| / (h_ratio^p - 1) */
    if (h_ratio <= 1.0) return fabs(err_fine);
    return Fs * fabs(err_coarse - err_fine) / (pow(h_ratio, p) - 1.0);
}

/* Compute mass flux balance for 2D incompressible flow. */
double cfd_mass_flux_balance_2d(const CfdVectorField2D *vel,
                                  const CfdGrid2D *grid)
{
    if (!vel || !grid) return 0.0;
    int nx = grid->nx, ny = grid->ny;
    double total_div = 0.0;
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            double div = cfd_divergence_2d(vel, i, j, grid);
            total_div += fabs(div) * grid->dx * grid->dy;
        }
    }
    return total_div;
}

/* Compute kinetic energy (per unit mass) of 2D flow field. */
double cfd_kinetic_energy_2d(const CfdVectorField2D *vel)
{
    if (!vel) return 0.0;
    double ke = 0.0;
    for (int i = 0; i < vel->ncells; i++)
        ke += 0.5 * (vel->u[i]*vel->u[i] + vel->v[i]*vel->v[i]);
    return ke / (double)vel->ncells;
}

/* Compute enstrophy (vorticity squared) for 2D flow. */
double cfd_enstrophy_2d(const CfdVectorField2D *vel, const CfdGrid2D *grid)
{
    if (!vel || !grid) return 0.0;
    int nx = grid->nx, ny = grid->ny;
    double ens = 0.0;
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int idx = CFD_IDX2(i, j, nx);
            double dudy = 0.0, dvdx = 0.0;
            if (j > 0 && j < ny-1)
                dudy = (vel->u[CFD_IDX2(i,j+1,nx)] - vel->u[CFD_IDX2(i,j-1,nx)]) / (2.0*grid->dy);
            if (i > 0 && i < nx-1)
                dvdx = (vel->v[CFD_IDX2(i+1,j,nx)] - vel->v[CFD_IDX2(i-1,j,nx)]) / (2.0*grid->dx);
            double vort_z = dvdx - dudy;
            ens += vort_z * vort_z;
        }
    }
    return ens / (double)(nx * ny);
}

/* ===== Additional 2D Steady Solver ===== */
int cfd_solve_convdiff_2d_gs(const CfdGrid2D *grid, const CfdTransportEquation *eqn,
                             const CfdVectorField2D *vel, const CfdBCSet2D *bc,
                             CfdScalarField2D *phi, int max_iter, double tol)
{
    if (!grid || !eqn || !vel || !bc || !phi) return -1;
    int nx = grid->nx, ny = grid->ny, nc = nx * ny;
    double dx = grid->dx, dy = grid->dy;
    double Gamma = eqn->diffusivity, rho = eqn->density;

    for (int iter = 0; iter < max_iter; iter++) {
        double max_res = 0.0;
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                int idx = CFD_IDX2(i, j, nx);
                double phi_P = phi->data[idx];
                double aW = 0, aE = 0, aS = 0, aN = 0;
                double Dw = Gamma / dx, De = Gamma / dx;
                double Ds = Gamma / dy, Dn = Gamma / dy;
                if (i > 0) { double u_w = vel->u[CFD_IDX2(i-1,j,nx)]; aW = De + CFD_MAX(rho*u_w, 0.0); }
                if (i < nx-1) { double u_e = vel->u[CFD_IDX2(i+1,j,nx)]; aE = Dw + CFD_MAX(-rho*u_e, 0.0); }
                if (j > 0) { double v_s = vel->v[CFD_IDX2(i,j-1,nx)]; aS = Dn + CFD_MAX(rho*v_s, 0.0); }
                if (j < ny-1) { double v_n = vel->v[CFD_IDX2(i,j+1,nx)]; aN = Ds + CFD_MAX(-rho*v_n, 0.0); }
                double aP = aW + aE + aS + aN;
                double sum_nb = 0.0;
                if (i > 0) sum_nb += aW * phi->data[CFD_IDX2(i-1,j,nx)];
                if (i < nx-1) sum_nb += aE * phi->data[CFD_IDX2(i+1,j,nx)];
                if (j > 0) sum_nb += aS * phi->data[CFD_IDX2(i,j-1,nx)];
                if (j < ny-1) sum_nb += aN * phi->data[CFD_IDX2(i,j+1,nx)];
                double phi_new = (eqn->source_constant + sum_nb) / CFD_MAX(aP, 1e-15);
                double res = fabs(phi_new - phi_P);
                if (res > max_res) max_res = res;
                phi->data[idx] = phi_new;
            }
        }
        if (max_res < tol) return iter + 1;
    }
    return max_iter;
}

/* ===== Post-Processing Utilities ===== */
double cfd_streamfunction_2d(const CfdVectorField2D *vel, int i, int j, const CfdGrid2D *grid) {
    /* Compute streamfunction via integration: psi = -int(u*dy) + int(v*dx) */
    if (!vel || !grid || i < 0 || j < 0) return 0.0;
    double psi = 0.0;
    for (int jj = 0; jj <= j; jj++)
        psi += vel->u[CFD_IDX2(i, CFD_MIN(jj,grid->ny-1), grid->nx)] * grid->dy;
    for (int ii = 0; ii <= i; ii++)
        psi -= vel->v[CFD_IDX2(CFD_MIN(ii,grid->nx-1), 0, grid->nx)] * grid->dx;
    return psi;
}

double cfd_vorticity_2d(const CfdVectorField2D *vel, int i, int j, const CfdGrid2D *grid) {
    if (!vel || !grid) return 0.0;
    int nx = grid->nx, ny = grid->ny;
    if (i <= 0 || i >= nx-1 || j <= 0 || j >= ny-1) return 0.0;
    double dvdx = (vel->v[CFD_IDX2(i+1,j,nx)] - vel->v[CFD_IDX2(i-1,j,nx)]) / (2.0*grid->dx);
    double dudy = (vel->u[CFD_IDX2(i,j+1,nx)] - vel->u[CFD_IDX2(i,j-1,nx)]) / (2.0*grid->dy);
    return dvdx - dudy;
}

/* ===== Pressure Drop Computation ===== */
double cfd_pressure_drop_1d(const double *p, int nx, double dx) {
    if (!p || nx < 2) return 0.0;
    return p[0] - p[nx-1];
}

/* ===== Nusselt Number from Temperature Field ===== */
double cfd_nusselt_from_field(const CfdScalarField2D *T, const CfdGrid2D *grid,
                                double T_wall, double T_bulk, double L, double k) {
    if (!T || !grid || k <= 0.0) return 0.0;
    int nx = grid->nx;
    double dT_dn = (T->data[CFD_IDX2(0,0,nx)] - T_wall) / (0.5*grid->dx);
    double q_wall = -k * dT_dn;
    double h = q_wall / (T_wall - T_bulk);
    return h * L / k;
}

/* ===== Mass Transfer Analogy (Chilton-Colburn) ===== */
double cfd_chilton_colburn_j_factor(double Nu, double Re, double Pr) {
    if (Re <= 0.0 || Pr <= 0.0) return 0.0;
    return Nu / (Re * Pr) * pow(Pr, 2.0/3.0);
}
double cfd_sherwood_from_j_factor(double j_D, double Re, double Sc) {
    return j_D * Re * pow(Sc, 1.0/3.0);
}

/* ===== Wall Shear Stress from Velocity Gradient ===== */
double cfd_wall_shear_from_velocity(const CfdVectorField2D *vel,
                                      const CfdGrid2D *grid, double mu,
                                      int wall_face) {
    if (!vel || !grid || mu <= 0.0) return 0.0;
    int nx = grid->nx;
    double dudy, tau;
    switch (wall_face) {
        case 0: /* west wall */
            dudy = vel->u[CFD_IDX2(0,0,nx)] / (0.5*grid->dx); tau = mu*dudy; break;
        case 1: /* south wall */
            dudy = vel->u[CFD_IDX2(0,0,nx)] / (0.5*grid->dy); tau = mu*dudy; break;
        case 2: /* north wall (lid) */
            dudy = vel->u[CFD_IDX2(0,grid->ny-1,nx)] / (0.5*grid->dy); tau = mu*dudy; break;
        default: return 0.0;
    }
    return tau;
}

/* ===== Energy Balance Check ===== */
double cfd_energy_balance_2d(const CfdVectorField2D *vel, const CfdGrid2D *grid,
                               const CfdBCSet2D *bc, double k) {
    if (!vel || !grid || !bc) return 0.0;
    int nx = grid->nx, ny = grid->ny;
    double energy_in = 0.0, energy_out = 0.0;
    /* West boundary flux */
    for (int j = 0; j < ny; j++) energy_in += fabs(vel->u[CFD_IDX2(0,j,nx)]) * grid->dy;
    /* East boundary flux */
    for (int j = 0; j < ny; j++) energy_out += fabs(vel->u[CFD_IDX2(nx-1,j,nx)]) * grid->dy;
    return fabs(energy_in - energy_out) / CFD_MAX(energy_in, 1e-15);
}

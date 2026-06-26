#include "cfd_discretization.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

double cfd_central_diff_1d(double uw, double ue, double dx) {
    if (dx == 0.0) return 0.0;
    return (ue - uw) / (2.0 * dx);
}

double cfd_upwind_diff_1d(double uw, double up, double dx, double vel) {
    if (dx == 0.0) return 0.0;
    if (vel > 0.0) return (up - uw) / dx;
    if (vel < 0.0) return (up - uw) / dx;
    return 0.0;
}

double cfd_second_deriv_1d(double uw, double up, double ue, double dx) {
    if (dx == 0.0) return 0.0;
    return (ue - 2.0 * up + uw) / (dx * dx);
}

/* 2D gradient operators */
double cfd_grad_x_2d(const CfdScalarField2D *field, int i, int j,
                     const CfdGrid2D *grid)
{
    if (!field || !grid || !field->data) return 0.0;
    int nx = field->nx, ny = field->ny;
    if (i < 0 || i >= nx || j < 0 || j >= ny) return 0.0;
    double dx = grid->dx;
    if (i == 0) {
        int ie = CFD_IDX2(1, j, nx);
        int ic = CFD_IDX2(0, j, nx);
        return (field->data[ie] - field->data[ic]) / dx;
    }
    if (i == nx - 1) {
        int ic = CFD_IDX2(i, j, nx);
        int iw = CFD_IDX2(i - 1, j, nx);
        return (field->data[ic] - field->data[iw]) / dx;
    }
    int ie = CFD_IDX2(i + 1, j, nx);
    int iw = CFD_IDX2(i - 1, j, nx);
    return (field->data[ie] - field->data[iw]) / (2.0 * dx);
}

double cfd_grad_y_2d(const CfdScalarField2D *field, int i, int j,
                     const CfdGrid2D *grid)
{
    if (!field || !grid || !field->data) return 0.0;
    int nx = field->nx, ny = field->ny;
    if (i < 0 || i >= nx || j < 0 || j >= ny) return 0.0;
    double dy = grid->dy;
    if (j == 0) {
        int inr = CFD_IDX2(i, 1, nx);
        int ic = CFD_IDX2(i, 0, nx);
        return (field->data[inr] - field->data[ic]) / dy;
    }
    if (j == ny - 1) {
        int ic = CFD_IDX2(i, j, nx);
        int js = CFD_IDX2(i, j - 1, nx);
        return (field->data[ic] - field->data[js]) / dy;
    }
    int inr = CFD_IDX2(i, j + 1, nx);
    int js = CFD_IDX2(i, j - 1, nx);
    return (field->data[inr] - field->data[js]) / (2.0 * dy);
}

double cfd_grad_mag_2d(const CfdScalarField2D *field, int i, int j,
                       const CfdGrid2D *grid)
{
    double gx = cfd_grad_x_2d(field, i, j, grid);
    double gy = cfd_grad_y_2d(field, i, j, grid);
    return sqrt(gx * gx + gy * gy);
}


/* 2D Laplacian (5-point stencil) */
double cfd_laplacian_2d(const CfdScalarField2D *field, int i, int j,
                        const CfdGrid2D *grid)
{
    if (!field || !grid || !field->data) return 0.0;
    int nx = field->nx, ny = field->ny;
    if (i < 0 || i >= nx || j < 0 || j >= ny) return 0.0;
    double dx = grid->dx, dy = grid->dy;
    int ic = CFD_IDX2(i, j, nx);
    double d2x, d2y;

    if (i == 0)
        d2x = (field->data[CFD_IDX2(1,j,nx)] - field->data[ic]) / (dx*dx);
    else if (i == nx - 1)
        d2x = (field->data[CFD_IDX2(i-1,j,nx)] - field->data[ic]) / (dx*dx);
    else
        d2x = (field->data[CFD_IDX2(i+1,j,nx)] - 2.0*field->data[ic]
               + field->data[CFD_IDX2(i-1,j,nx)]) / (dx*dx);

    if (j == 0)
        d2y = (field->data[CFD_IDX2(i,1,nx)] - field->data[ic]) / (dy*dy);
    else if (j == ny - 1)
        d2y = (field->data[CFD_IDX2(i,j-1,nx)] - field->data[ic]) / (dy*dy);
    else
        d2y = (field->data[CFD_IDX2(i,j+1,nx)] - 2.0*field->data[ic]
               + field->data[CFD_IDX2(i,j-1,nx)]) / (dy*dy);
    return d2x + d2y;
}

/* 2D Velocity Divergence */
double cfd_divergence_2d(const CfdVectorField2D *vel, int i, int j,
                         const CfdGrid2D *grid)
{
    if (!vel || !grid || !vel->u || !vel->v) return 0.0;
    int nx = vel->nx, ny = vel->ny;
    if (i < 0 || i >= nx || j < 0 || j >= ny) return 0.0;
    double dx = grid->dx, dy = grid->dy;
    double dudx, dvdy;

    if (i == 0)
        dudx = (vel->u[CFD_IDX2(1,j,nx)] - vel->u[CFD_IDX2(0,j,nx)]) / dx;
    else if (i == nx - 1)
        dudx = (vel->u[CFD_IDX2(i,j,nx)] - vel->u[CFD_IDX2(i-1,j,nx)]) / dx;
    else
        dudx = (vel->u[CFD_IDX2(i+1,j,nx)] - vel->u[CFD_IDX2(i-1,j,nx)]) / (2.0*dx);

    if (j == 0)
        dvdy = (vel->v[CFD_IDX2(i,1,nx)] - vel->v[CFD_IDX2(i,0,nx)]) / dy;
    else if (j == ny - 1)
        dvdy = (vel->v[CFD_IDX2(i,j,nx)] - vel->v[CFD_IDX2(i,j-1,nx)]) / dy;
    else
        dvdy = (vel->v[CFD_IDX2(i,j+1,nx)] - vel->v[CFD_IDX2(i,j-1,nx)]) / (2.0*dy);
    return dudx + dvdy;
}


/* ===== Convection-Diffusion FVM Face Coefficients ===== */

void cfd_powerlaw_1d(double Pe, double F, double D,
                     double *aW, double *aE, double *aP)
{
    /* Power-law scheme (Patankar 1980): A(|Pe|) = max(0, (1-0.1|Pe|)^5) */
    if (!aW || !aE || !aP) return;
    double absPe = fabs(Pe);
    double A_val;
    if (absPe < 1e-10) {
        A_val = 1.0;
    } else {
        double term = 1.0 - 0.1 * absPe;
        if (term < 0.0) term = 0.0;
        double t2 = term * term;
        double t4 = t2 * t2;
        A_val = t4 * term;
    }
    *aW = D * A_val + CFD_MAX(F, 0.0);
    *aE = D * A_val + CFD_MAX(-F, 0.0);
    *aP = *aW + *aE;
}

void cfd_hybrid_1d(double Pe, double F, double D,
                   double *aW, double *aE, double *aP)
{
    /* Hybrid scheme (Spalding 1972): central for |Pe|<=2, upwind otherwise */
    if (!aW || !aE || !aP) return;
    double absPe = fabs(Pe);
    double A_val = (absPe <= 2.0) ? 1.0 : 0.0;
    *aW = D * A_val + CFD_MAX(F, 0.0);
    *aE = D * A_val + CFD_MAX(-F, 0.0);
    *aP = *aW + *aE;
}

/* ===== 1D Steady Convection-Diffusion Assembly ===== */

void cfd_assemble_convdiff_1d(const CfdGrid1D *grid,
                               const CfdTransportEquation *eqn,
                               double velocity,
                               double phi_left, double phi_right,
                               CfdConvectiveScheme scheme,
                               CfdMatrix2D *mat)
{
    if (!grid || !eqn || !mat) return;
    int nx = grid->nx;
    double dx = grid->dx;
    double rho = eqn->density;
    double Gamma = eqn->diffusivity;
    double F_conv = rho * velocity;
    double D_diff = Gamma / dx;
    double Sc = eqn->source_constant;
    double Sp = eqn->source_linear;
    int i;

    for (i = 0; i < nx; i++) {
        double Pe_cell = (Gamma > 1e-15) ? F_conv * dx / Gamma : 0.0;
        double aWc = 0.0, aEc = 0.0, aPc = 0.0;

        if (scheme == CFD_SCHEME_POWERLAW) {
            cfd_powerlaw_1d(Pe_cell, F_conv, D_diff, &aWc, &aEc, &aPc);
        } else if (scheme == CFD_SCHEME_HYBRID) {
            cfd_hybrid_1d(Pe_cell, F_conv, D_diff, &aWc, &aEc, &aPc);
        } else {
            aWc = D_diff + CFD_MAX(F_conv, 0.0);
            aEc = D_diff + CFD_MAX(-F_conv, 0.0);
            aPc = aWc + aEc;
        }

        if (i > 0)      mat->aW[i] = aWc;
        if (i < nx - 1) mat->aE[i] = aEc;
        mat->aP[i] = aPc - Sp * dx;
        mat->b[i]  = Sc * dx;
    }

    /* Dirichlet BC at left (x=0) */
    mat->aP[0] += 2.0 * D_diff + CFD_MAX(F_conv, 0.0);
    mat->b[0]  += (2.0 * D_diff + CFD_MAX(F_conv, 0.0)) * phi_left;
    mat->aW[0] = 0.0;

    /* Dirichlet BC at right (x=L) */
    int last = nx - 1;
    mat->aP[last] += 2.0 * D_diff + CFD_MAX(-F_conv, 0.0);
    mat->b[last]  += (2.0 * D_diff + CFD_MAX(-F_conv, 0.0)) * phi_right;
    mat->aE[last] = 0.0;
}


/* ===== 2D Poisson Equation Assembly (FVM) ===== */

void cfd_assemble_poisson_2d(const CfdGrid2D *grid,
                              const CfdScalarField2D *source,
                              const CfdBCSet2D *bc,
                              CfdMatrix2D *mat)
{
    if (!grid || !mat) return;
    int nx = grid->nx, ny = grid->ny;
    double dx = grid->dx, dy = grid->dy;
    double dxi2 = 1.0 / (dx * dx);
    double dyi2 = 1.0 / (dy * dy);
    int i, j;

    cfd_matrix2d_zero(mat);

    for (j = 0; j < ny; j++) {
        for (i = 0; i < nx; i++) {
            int idx = CFD_IDX2(i, j, nx);
            double src_val = (source && source->data) ? source->data[idx] : 0.0;

            double aW = (i > 0)      ? dxi2 : 0.0;
            double aE = (i < nx - 1) ? dxi2 : 0.0;
            double aS = (j > 0)      ? dyi2 : 0.0;
            double aN = (j < ny - 1) ? dyi2 : 0.0;
            double aP = aW + aE + aS + aN;

            if (bc) {
                if (i == 0)      { aW = 0.0; aP += dxi2; mat->b[idx] += dxi2 * bc->west.value; }
                if (i == nx - 1) { aE = 0.0; aP += dxi2; mat->b[idx] += dxi2 * bc->east.value; }
                if (j == 0)      { aS = 0.0; aP += dyi2; mat->b[idx] += dyi2 * bc->south.value; }
                if (j == ny - 1) { aN = 0.0; aP += dyi2; mat->b[idx] += dyi2 * bc->north.value; }
            }

            mat->aW[idx] = -aW;
            mat->aE[idx] = -aE;
            mat->aS[idx] = -aS;
            mat->aN[idx] = -aN;
            mat->aP[idx] = aP;
            mat->b[idx] += src_val;
        }
    }
}

/* ===== Momentum Equation Assembly (2D FVM) ===== */

void cfd_assemble_u_momentum_2d(const CfdGrid2D *grid,
                                 const CfdVectorField2D *vel,
                                 const CfdFluidProperty *fp,
                                 CfdConvectiveScheme scheme,
                                 CfdMatrix2D *u_mat)
{
    if (!grid || !vel || !fp || !u_mat) return;
    int nx = grid->nx, ny = grid->ny;
    double dx = grid->dx, dy = grid->dy;
    double mu = fp->viscosity, rho = fp->density;
    int i, j;

    cfd_matrix2d_zero(u_mat);

    for (j = 0; j < ny; j++) {
        for (i = 0; i < nx; i++) {
            int idx = CFD_IDX2(i, j, nx);
            double u_P = vel->u[idx];
            double DW = mu * dy / dx, DE = mu * dy / dx;
            double DS = mu * dx / dy, DN = mu * dx / dy;
            double FW = 0, FE = 0, FS = 0, FN = 0;

            if (i > 0) { double u_w=0.5*(vel->u[CFD_IDX2(i-1,j,nx)]+u_P); FW=rho*u_w*dy; }
            if (i < nx-1) { double u_e=0.5*(u_P+vel->u[CFD_IDX2(i+1,j,nx)]); FE=rho*u_e*dy; }
            if (j > 0) { double v_s=0.5*(vel->v[CFD_IDX2(i,j-1,nx)]+vel->v[idx]); FS=rho*v_s*dx; }
            if (j < ny-1) { double v_n=0.5*(vel->v[idx]+vel->v[CFD_IDX2(i,j+1,nx)]); FN=rho*v_n*dx; }

            double aW=0, aE=0, aS=0, aN=0;
            if (i > 0)      aW = DW + CFD_MAX(FW, 0.0);
            if (i < nx - 1) aE = DE + CFD_MAX(-FE, 0.0);
            if (j > 0)      aS = DS + CFD_MAX(FS, 0.0);
            if (j < ny - 1) aN = DN + CFD_MAX(-FN, 0.0);
            double aP = aW + aE + aS + aN + (FE - FW + FN - FS);

            u_mat->aW[idx] = aW; u_mat->aE[idx] = aE;
            u_mat->aS[idx] = aS; u_mat->aN[idx] = aN;
            u_mat->aP[idx] = aP; u_mat->b[idx] = 0.0;
        }
    }
}

void cfd_assemble_v_momentum_2d(const CfdGrid2D *grid,
                                 const CfdVectorField2D *vel,
                                 const CfdFluidProperty *fp,
                                 CfdConvectiveScheme scheme,
                                 CfdMatrix2D *v_mat)
{
    if (!grid || !vel || !fp || !v_mat) return;
    int nx = grid->nx, ny = grid->ny;
    double dx = grid->dx, dy = grid->dy;
    double mu = fp->viscosity, rho = fp->density;
    int i, j;

    cfd_matrix2d_zero(v_mat);

    for (j = 0; j < ny; j++) {
        for (i = 0; i < nx; i++) {
            int idx = CFD_IDX2(i, j, nx);
            double v_P = vel->v[idx];
            double DW = mu * dy / dx, DE = mu * dy / dx;
            double DS = mu * dx / dy, DN = mu * dx / dy;
            double FW = 0, FE = 0, FS = 0, FN = 0;

            if (i > 0) { double u_w=0.5*(vel->u[CFD_IDX2(i-1,j,nx)]+vel->u[idx]); FW=rho*u_w*dy; }
            if (i < nx-1) { double u_e=0.5*(vel->u[idx]+vel->u[CFD_IDX2(i+1,j,nx)]); FE=rho*u_e*dy; }
            if (j > 0) { double v_s=0.5*(vel->v[CFD_IDX2(i,j-1,nx)]+v_P); FS=rho*v_s*dx; }
            if (j < ny-1) { double v_n=0.5*(v_P+vel->v[CFD_IDX2(i,j+1,nx)]); FN=rho*v_n*dx; }

            double aW=0, aE=0, aS=0, aN=0;
            if (i > 0)      aW = DW + CFD_MAX(FW, 0.0);
            if (i < nx - 1) aE = DE + CFD_MAX(-FE, 0.0);
            if (j > 0)      aS = DS + CFD_MAX(FS, 0.0);
            if (j < ny - 1) aN = DN + CFD_MAX(-FN, 0.0);
            double aP = aW + aE + aS + aN + (FE - FW + FN - FS);

            v_mat->aW[idx] = aW; v_mat->aE[idx] = aE;
            v_mat->aS[idx] = aS; v_mat->aN[idx] = aN;
            v_mat->aP[idx] = aP; v_mat->b[idx] = 0.0;
        }
    }
}


/* ===== Stability Analysis ===== */

double cfd_stable_timestep_1d(double velocity, double diffusivity,
                               double dx, double CFL_target)
{
    if (dx <= 0.0) return 0.0;
    double dt_conv = 1e100, dt_diff = 1e100;
    if (fabs(velocity) > 1e-15)
        dt_conv = CFL_target * dx / fabs(velocity);
    if (diffusivity > 1e-15)
        dt_diff = dx * dx / (2.0 * diffusivity);
    return (dt_conv < dt_diff) ? dt_conv : dt_diff;
}

double cfd_cell_peclet(double velocity, double dx, double diffusivity)
{
    if (diffusivity <= 1e-15) return 1e100;
    return fabs(velocity) * dx / diffusivity;
}

/* ===== Discrete Maximum Principle Check ===== */

int cfd_check_maximum_principle(const CfdMatrix2D *mat)
{
    if (!mat) return 0;
    int n = mat->ncells, i;
    for (i = 0; i < n; i++) {
        if (mat->aP[i] <= 0.0) return 0;
        if (mat->aW[i] < -1e-15 || mat->aE[i] < -1e-15) return 0;
        if (mat->aS[i] < -1e-15 || mat->aN[i] < -1e-15) return 0;
        double sum_nb = mat->aW[i] + mat->aE[i] + mat->aS[i] + mat->aN[i];
        if (mat->aP[i] < sum_nb - 1e-12) return 0;
    }
    return 1;
}

/* ===== Additional Gradient Operators ===== */
void cfd_grad_2d(const CfdScalarField2D *field, int i, int j,
                 const CfdGrid2D *grid, double grad[2]) {
    grad[0] = cfd_grad_x_2d(field, i, j, grid);
    grad[1] = cfd_grad_y_2d(field, i, j, grid);
}
double cfd_divergence_scalar_flux_2d(const CfdScalarField2D *flux_x,
                                      const CfdScalarField2D *flux_y,
                                      int i, int j, const CfdGrid2D *grid) {
    double dfx_dx = cfd_grad_x_2d(flux_x, i, j, grid);
    double dfy_dy = cfd_grad_y_2d(flux_y, i, j, grid);
    return dfx_dx + dfy_dy;
}

/* ===== Upwind Convection Flux ===== */
double cfd_upwind_flux_1d(double phi_w, double phi_p, double phi_e, double u) {
    if (u > 0.0) return u * phi_p;
    else         return u * phi_e;
}

/* ===== QUICK Scheme (Leonard 1979) ===== */
double cfd_quick_face_value(double phi_ww, double phi_w, double phi_p,
                             double phi_e, double u) {
    if (u > 0.0)
        return 0.375 * phi_p + 0.75 * phi_w - 0.125 * phi_ww;
    else
        return 0.375 * phi_p + 0.75 * phi_e - 0.125 * phi_e;
}

/* ===== TVD Limiter Functions ===== */
double cfd_limiter_minmod(double r) {
    return CFD_MAX(0.0, CFD_MIN(1.0, r));
}
double cfd_limiter_van_leer(double r) {
    if (r <= 0.0) return 0.0;
    return 2.0 * r / (1.0 + r);
}
double cfd_limiter_superbee(double r) {
    return CFD_MAX(0.0, CFD_MAX(CFD_MIN(2.0*r, 1.0), CFD_MIN(r, 2.0)));
}
double cfd_limiter_mc(double r) {
    return CFD_MAX(0.0, CFD_MIN(CFD_MIN(2.0*r, 0.5*(1.0+r)), 2.0));
}

/* ===== Gradient Ratio for TVD ===== */
double cfd_gradient_ratio(double phi_u, double phi_c, double phi_d) {
    double num = phi_c - phi_u;
    double den = phi_d - phi_c;
    if (fabs(den) < 1e-15) return (num > 0.0) ? 1e10 : -1e10;
    return num / den;
}

/* ===== WENO5 Weights (Jiang & Shu 1996) ===== */
void cfd_weno5_weights(double v1, double v2, double v3, double v4, double v5,
                        double w[3]) {
    double eps = 1e-6;
    double beta1 = 13.0/12.0*(v1-2.0*v2+v3)*(v1-2.0*v2+v3) + 0.25*(v1-4.0*v2+3.0*v3)*(v1-4.0*v2+3.0*v3);
    double beta2 = 13.0/12.0*(v2-2.0*v3+v4)*(v2-2.0*v3+v4) + 0.25*(v2-v4)*(v2-v4);
    double beta3 = 13.0/12.0*(v3-2.0*v4+v5)*(v3-2.0*v4+v5) + 0.25*(3.0*v3-4.0*v4+v5)*(3.0*v3-4.0*v4+v5);
    double d1=0.1, d2=0.6, d3=0.3;
    double alpha1 = d1/((eps+beta1)*(eps+beta1));
    double alpha2 = d2/((eps+beta2)*(eps+beta2));
    double alpha3 = d3/((eps+beta3)*(eps+beta3));
    double alpha_sum = alpha1+alpha2+alpha3;
    w[0] = alpha1/alpha_sum; w[1] = alpha2/alpha_sum; w[2] = alpha3/alpha_sum;
}

/* ===== WENO5 Reconstruction ===== */
double cfd_weno5_reconstruct(double v1, double v2, double v3, double v4, double v5) {
    double w[3]; cfd_weno5_weights(v1, v2, v3, v4, v5, w);
    double p1 = (2.0*v1 - 7.0*v2 + 11.0*v3)/6.0;
    double p2 = (-v2 + 5.0*v3 + 2.0*v4)/6.0;
    double p3 = (2.0*v3 + 5.0*v4 - v5)/6.0;
    return w[0]*p1 + w[1]*p2 + w[2]*p3;
}

/* ===== Gauss Quadrature Points (2-point) ===== */
void cfd_gauss_quadrature_2pt(double a, double b, double *x, double *w) {
    double mid = 0.5*(a+b); double half = 0.5*(b-a);
    double xi = 0.5773502691896257; /* 1/sqrt(3) */
    x[0] = mid - half*xi; x[1] = mid + half*xi;
    w[0] = half; w[1] = half;
}

/* ===== Interpolation Functions ===== */
double cfd_linear_interp(double x, double x0, double x1, double f0, double f1) {
    if (fabs(x1-x0) < 1e-15) return f0;
    return f0 + (f1-f0)*(x-x0)/(x1-x0);
}
double cfd_bilinear_interp(double x, double y,
                            double x0, double x1, double y0, double y1,
                            double f00, double f10, double f01, double f11) {
    double tx = (fabs(x1-x0)>1e-15) ? (x-x0)/(x1-x0) : 0.0;
    double ty = (fabs(y1-y0)>1e-15) ? (y-y0)/(y1-y0) : 0.0;
    return (1.0-tx)*(1.0-ty)*f00 + tx*(1.0-ty)*f10 + (1.0-tx)*ty*f01 + tx*ty*f11;
}

/* ===== Numerical Integration ===== */
double cfd_trapezoidal_integrate(const double *f, double dx, int n) {
    if (n < 2) return 0.0;
    double sum = 0.5*(f[0]+f[n-1]);
    for (int i = 1; i < n-1; i++) sum += f[i];
    return sum * dx;
}
double cfd_simpson_integrate(const double *f, double dx, int n) {
    if (n < 3 || n%2 == 0) return cfd_trapezoidal_integrate(f, dx, n);
    double sum = f[0] + f[n-1];
    for (int i = 1; i < n-1; i++) sum += f[i] * ((i%2==0) ? 2.0 : 4.0);
    return sum * dx / 3.0;
}

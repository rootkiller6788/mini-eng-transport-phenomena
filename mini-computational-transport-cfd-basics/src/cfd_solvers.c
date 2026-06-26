#include "cfd_solvers.h"
#include "cfd_discretization.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ===== Residual Computation ===== */

void cfd_compute_residual(const CfdMatrix2D *A, const double *x, double *r)
{
    if (!A || !x || !r) return;
    int nx = A->nx, ny = A->ny, i, j;
    for (j = 0; j < ny; j++) {
        for (i = 0; i < nx; i++) {
            int idx = CFD_IDX2(i, j, nx);
            double Ax = A->aP[idx] * x[idx];
            if (i > 0)      Ax -= A->aW[idx] * x[CFD_IDX2(i-1, j, nx)];
            if (i < nx - 1) Ax -= A->aE[idx] * x[CFD_IDX2(i+1, j, nx)];
            if (j > 0)      Ax -= A->aS[idx] * x[CFD_IDX2(i, j-1, nx)];
            if (j < ny - 1) Ax -= A->aN[idx] * x[CFD_IDX2(i, j+1, nx)];
            r[idx] = A->b[idx] - Ax;
        }
    }
}

double cfd_residual_l2(const double *r, int n)
{
    if (!r || n < 1) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < n; i++) sum += r[i] * r[i];
    return sqrt(sum);
}

double cfd_residual_linf(const double *r, int n)
{
    if (!r || n < 1) return 0.0;
    double maxv = 0.0;
    for (int i = 0; i < n; i++) {
        double av = fabs(r[i]);
        if (av > maxv) maxv = av;
    }
    return maxv;
}

/* ===== Jacobi Iteration ===== */

int cfd_solve_jacobi(const CfdMatrix2D *A, double *x, CfdSolverState *state)
{
    if (!A || !x || !state) return -1;
    int nc = A->ncells, nx = A->nx, ny = A->ny;
    double *x_new = (double *)malloc((size_t)nc * sizeof(double));
    if (!x_new) return -1;

    double omega = state->under_relaxation;
    double init_r = 0.0;
    for (int iter = 0; iter < state->max_iter; iter++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                int idx = CFD_IDX2(i, j, nx);
                double sum = A->b[idx];
                if (i > 0)      sum += A->aW[idx] * x[CFD_IDX2(i-1, j, nx)];
                if (i < nx - 1) sum += A->aE[idx] * x[CFD_IDX2(i+1, j, nx)];
                if (j > 0)      sum += A->aS[idx] * x[CFD_IDX2(i, j-1, nx)];
                if (j < ny - 1) sum += A->aN[idx] * x[CFD_IDX2(i, j+1, nx)];
                double diag = A->aP[idx];
                if (fabs(diag) < 1e-15) diag = 1e-15;
                x_new[idx] = (1.0 - omega) * x[idx] + omega * sum / diag;
            }
        }
        memcpy(x, x_new, (size_t)nc * sizeof(double));

        /* Residual check every 10 iterations */
        if (iter % 10 == 0 || iter == state->max_iter - 1) {
            double *r = (double *)malloc((size_t)nc * sizeof(double));
            if (r) {
                cfd_compute_residual(A, x, r);
                double rnorm = cfd_residual_l2(r, nc);
                if (iter == 0) init_r = rnorm;
                state->final_residual = rnorm;
                free(r);
                if (init_r > 0.0 && rnorm / init_r < state->tolerance) {
                    state->n_iter = iter + 1;
                    free(x_new);
                    return iter + 1;
                }
            }
        }

        if (iter == state->max_iter - 1) {
            state->n_iter = state->max_iter;
            state->diverged = 1;
            free(x_new);
            return -1;
        }
    }
    free(x_new);
    return -1;
}

/* ===== Gauss-Seidel Iteration ===== */

int cfd_solve_gauss_seidel(const CfdMatrix2D *A, double *x, CfdSolverState *state)
{
    if (!A || !x || !state) return -1;
    int nc = A->ncells, nx = A->nx, ny = A->ny;
    double omega = state->under_relaxation;
    double init_r = 0.0;

    for (int iter = 0; iter < state->max_iter; iter++) {
        double max_delta = 0.0;
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                int idx = CFD_IDX2(i, j, nx);
                double sum = A->b[idx];
                if (i > 0)      sum += A->aW[idx] * x[CFD_IDX2(i-1, j, nx)];
                if (i < nx - 1) sum += A->aE[idx] * x[CFD_IDX2(i+1, j, nx)];
                if (j > 0)      sum += A->aS[idx] * x[CFD_IDX2(i, j-1, nx)];
                if (j < ny - 1) sum += A->aN[idx] * x[CFD_IDX2(i, j+1, nx)];
                double diag = A->aP[idx];
                if (fabs(diag) < 1e-15) diag = 1e-15;
                double x_new = (1.0 - omega) * x[idx] + omega * sum / diag;
                double delta = fabs(x_new - x[idx]);
                if (delta > max_delta) max_delta = delta;
                x[idx] = x_new;
            }
        }

        if (iter % 10 == 0 || iter == state->max_iter - 1) {
            double *r = (double *)malloc((size_t)nc * sizeof(double));
            if (r) {
                cfd_compute_residual(A, x, r);
                double rnorm = cfd_residual_l2(r, nc);
                if (iter == 0) init_r = rnorm;
                state->final_residual = rnorm;
                free(r);
                if (init_r > 0.0 && rnorm / init_r < state->tolerance) {
                    state->n_iter = iter + 1;
                    return iter + 1;
                }
            }
        }
    }
    state->n_iter = state->max_iter;
    state->diverged = 1;
    return -1;
}


/* ===== SOR Iteration ===== */

int cfd_solve_sor(const CfdMatrix2D *A, double *x, double omega,
                  CfdSolverState *state)
{
    if (!A || !x || !state) return -1;
    int nc = A->ncells, nx = A->nx, ny = A->ny;
    double init_r = 0.0;
    for (int iter = 0; iter < state->max_iter; iter++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                int idx = CFD_IDX2(i, j, nx);
                double sum = A->b[idx];
                if (i > 0)      sum += A->aW[idx] * x[CFD_IDX2(i-1, j, nx)];
                if (i < nx - 1) sum += A->aE[idx] * x[CFD_IDX2(i+1, j, nx)];
                if (j > 0)      sum += A->aS[idx] * x[CFD_IDX2(i, j-1, nx)];
                if (j < ny - 1) sum += A->aN[idx] * x[CFD_IDX2(i, j+1, nx)];
                double diag = A->aP[idx];
                if (fabs(diag) < 1e-15) diag = 1e-15;
                x[idx] = (1.0 - omega) * x[idx] + omega * sum / diag;
            }
        }
        if (iter % 10 == 0 || iter == state->max_iter - 1) {
            double *r = malloc((size_t)nc * sizeof(double));
            if (r) {
                cfd_compute_residual(A, x, r);
                double rn = cfd_residual_l2(r, nc);
                if (iter == 0) init_r = rn;
                state->final_residual = rn; free(r);
                if (init_r > 0.0 && rn / init_r < state->tolerance) {
                    state->n_iter = iter + 1; return iter + 1;
                }
            }
        }
    }
    state->n_iter = state->max_iter; state->diverged = 1; return -1;
}

/* ===== Thomas TDMA ===== */

int cfd_solve_tdma(const CfdMatrix2D *A, double *x, int nx)
{
    if (!A || !x || nx < 1) return -1;
    double *c = malloc((size_t)nx * sizeof(double));
    double *d = malloc((size_t)nx * sizeof(double));
    if (!c || !d) { free(c); free(d); return -1; }
    c[0] = A->aE[0] / A->aP[0];
    d[0] = A->b[0] / A->aP[0];
    for (int i = 1; i < nx; i++) {
        double denom = A->aP[i] - A->aW[i] * c[i-1];
        if (fabs(denom) < 1e-15) denom = 1e-15;
        c[i] = (i < nx - 1) ? A->aE[i] / denom : 0.0;
        d[i] = (A->b[i] - A->aW[i] * d[i-1]) / denom;
    }
    x[nx - 1] = d[nx - 1];
    for (int i = nx - 2; i >= 0; i--) x[i] = d[i] - c[i] * x[i+1];
    free(c); free(d); return 0;
}


/* ===== Conjugate Gradient ===== */

int cfd_solve_cg(const CfdMatrix2D *A, double *x, CfdSolverState *state)
{
    if (!A || !x || !state) return -1;
    int nc = A->ncells;
    double *r = malloc((size_t)nc * sizeof(double));
    double *p = malloc((size_t)nc * sizeof(double));
    double *Ap = malloc((size_t)nc * sizeof(double));
    if (!r || !p || !Ap) { free(r); free(p); free(Ap); return -1; }

    cfd_compute_residual(A, x, r);
    double rsold = 0.0;
    for (int i = 0; i < nc; i++) { p[i] = r[i]; rsold += r[i] * r[i]; }
    double init_rs = rsold;

    for (int iter = 0; iter < state->max_iter; iter++) {
        cfd_compute_residual(A, p, Ap);
        double pAp = 0.0;
        for (int i = 0; i < nc; i++) pAp += p[i] * Ap[i];
        if (fabs(pAp) < 1e-20) { free(r); free(p); free(Ap); return -1; }
        double alpha = rsold / pAp;
        double rsnew = 0.0;
        for (int i = 0; i < nc; i++) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
            rsnew += r[i] * r[i];
        }
        if (sqrt(rsnew / init_rs) < state->tolerance) {
            state->n_iter = iter + 1; state->final_residual = sqrt(rsnew);
            free(r); free(p); free(Ap); return iter + 1;
        }
        double beta = rsnew / rsold;
        for (int i = 0; i < nc; i++) p[i] = r[i] + beta * p[i];
        rsold = rsnew;
    }
    state->n_iter = state->max_iter; state->diverged = 1;
    free(r); free(p); free(Ap); return -1;
}

/* ===== Preconditioned CG (Diagonal) ===== */

int cfd_solve_pcg_diag(const CfdMatrix2D *A, double *x, CfdSolverState *state)
{
    if (!A || !x || !state) return -1;
    int nc = A->ncells;
    double *r = malloc((size_t)nc * sizeof(double));
    double *z = malloc((size_t)nc * sizeof(double));
    double *p = malloc((size_t)nc * sizeof(double));
    double *Ap = malloc((size_t)nc * sizeof(double));
    if (!r || !z || !p || !Ap) { free(r); free(z); free(p); free(Ap); return -1; }

    cfd_compute_residual(A, x, r);
    double rz_old = 0.0;
    for (int i = 0; i < nc; i++) {
        double diag = A->aP[i];
        z[i] = r[i] / (fabs(diag) > 1e-15 ? diag : 1e-15);
        p[i] = z[i];
        rz_old += r[i] * z[i];
    }
    double init_rn = cfd_residual_l2(r, nc);

    for (int iter = 0; iter < state->max_iter; iter++) {
        cfd_compute_residual(A, p, Ap);
        double pAp = 0.0;
        for (int i = 0; i < nc; i++) pAp += p[i] * Ap[i];
        if (fabs(pAp) < 1e-20) { free(r); free(z); free(p); free(Ap); return -1; }
        double alpha = rz_old / pAp;
        double rz_new = 0.0;
        for (int i = 0; i < nc; i++) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
            double diag = A->aP[i];
            z[i] = r[i] / (fabs(diag) > 1e-15 ? diag : 1e-15);
            rz_new += r[i] * z[i];
        }
        if (sqrt(rz_new) / init_rn < state->tolerance) {
            state->n_iter = iter + 1; state->final_residual = sqrt(rz_new);
            free(r); free(z); free(p); free(Ap); return iter + 1;
        }
        double beta = rz_new / rz_old;
        for (int i = 0; i < nc; i++) p[i] = z[i] + beta * p[i];
        rz_old = rz_new;
    }
    state->n_iter = state->max_iter; state->diverged = 1;
    free(r); free(z); free(p); free(Ap); return -1;
}


/* ===== BiCGSTAB ===== */

int cfd_solve_bicgstab(const CfdMatrix2D *A, double *x, CfdSolverState *state)
{
    if (!A || !x || !state) return -1;
    int nc = A->ncells;
    double *r = malloc((size_t)nc*sizeof(double));
    double *rh = malloc((size_t)nc*sizeof(double));
    double *p = malloc((size_t)nc*sizeof(double));
    double *v = malloc((size_t)nc*sizeof(double));
    double *s = malloc((size_t)nc*sizeof(double));
    double *t = malloc((size_t)nc*sizeof(double));
    if (!r||!rh||!p||!v||!s||!t) { free(r);free(rh);free(p);free(v);free(s);free(t); return -1; }

    cfd_compute_residual(A, x, r);
    double init_rn = cfd_residual_l2(r, nc);
    for (int i=0;i<nc;i++) { rh[i]=r[i]; p[i]=0.0; v[i]=0.0; }
    double rho=1.0, alpha=1.0, omega=1.0;

    for (int iter=0; iter<state->max_iter; iter++) {
        double rn=0.0;
        for (int i=0;i<nc;i++) rn+=rh[i]*r[i];
        if (fabs(rn)<1e-30) { free(r);free(rh);free(p);free(v);free(s);free(t); return -1; }
        double beta=(rn/rho)*(alpha/omega);
        for (int i=0;i<nc;i++) p[i]=r[i]+beta*(p[i]-omega*v[i]);
        cfd_compute_residual(A, p, v);
        double rhv=0.0;
        for (int i=0;i<nc;i++) rhv+=rh[i]*v[i];
        if (fabs(rhv)<1e-30) { free(r);free(rh);free(p);free(v);free(s);free(t); return -1; }
        alpha=rn/rhv;
        for (int i=0;i<nc;i++) s[i]=r[i]-alpha*v[i];
        cfd_compute_residual(A, s, t);
        double ts=0.0, tt=0.0;
        for (int i=0;i<nc;i++) { ts+=t[i]*s[i]; tt+=t[i]*t[i]; }
        omega=(fabs(tt)>1e-30)?ts/tt:1.0;
        for (int i=0;i<nc;i++) { x[i]+=alpha*p[i]+omega*s[i]; r[i]=s[i]-omega*t[i]; }
        double rn2=cfd_residual_l2(r,nc);
        if(rn2/init_rn<state->tolerance){
            state->n_iter=iter+1;state->final_residual=rn2;
            free(r);free(rh);free(p);free(v);free(s);free(t);return iter+1;
        }
        rho=rn;
    }
    state->n_iter=state->max_iter;state->diverged=1;
    free(r);free(rh);free(p);free(v);free(s);free(t);
    return -1;
}


/* ===== LU Decomposition ===== */

int cfd_solve_lu_banded(double *A_band, double *b, double *x, int n, int bw)
{
    if (!A_band || !b || !x || n < 1 || n > 100) return -1;
    double *LU = malloc((size_t)n * n * sizeof(double));
    if (!LU) return -1;
    for (int i=0;i<n;i++) {
        for (int j=0;j<n;j++) {
            if (abs(i-j)<=bw) LU[i*n+j]=A_band[i*(2*bw+1)+(j-i+bw)];
            else LU[i*n+j]=0.0;
        }
    }
    for (int k=0;k<n;k++) {
        if (fabs(LU[k*n+k])<1e-15) { free(LU); return -1; }
        for (int i=k+1;i<n;i++) {
            double factor=LU[i*n+k]/LU[k*n+k];
            LU[i*n+k]=factor;
            for (int j=k+1;j<n;j++) LU[i*n+j]-=factor*LU[k*n+j];
        }
    }
    for (int i=0;i<n;i++) { x[i]=b[i]; for (int j=0;j<i;j++) x[i]-=LU[i*n+j]*x[j]; }
    for (int i=n-1;i>=0;i--) { for (int j=i+1;j<n;j++) x[i]-=LU[i*n+j]*x[j]; x[i]/=LU[i*n+i]; }
    free(LU);
    return 0;
}


/* ===== Multigrid Helpers ===== */

static void mg_restrict_2d(const double *fine, double *coarse, int nxf, int nyf)
{
    int nxc=nxf/2, nyc=nyf/2;
    for (int j=0;j<nyc;j++) for (int i=0;i<nxc;i++) {
        int i2=2*i, j2=2*j;
        coarse[CFD_IDX2(i,j,nxc)]=0.25*(
            fine[CFD_IDX2(i2,j2,nxf)]+fine[CFD_IDX2(i2+1,j2,nxf)]+
            fine[CFD_IDX2(i2,j2+1,nxf)]+fine[CFD_IDX2(i2+1,j2+1,nxf)]);
    }
}

static void mg_prolong_2d(const double *coarse, double *fine, int nxc, int nyc)
{
    int nxf=2*nxc;
    for (int j=0;j<nyc;j++) for (int i=0;i<nxc;i++) {
        double cv=coarse[CFD_IDX2(i,j,nxc)];
        int i2=2*i, j2=2*j;
        fine[CFD_IDX2(i2,j2,nxf)]+=cv;
        fine[CFD_IDX2(i2+1,j2,nxf)]+=cv;
        fine[CFD_IDX2(i2,j2+1,nxf)]+=cv;
        fine[CFD_IDX2(i2+1,j2+1,nxf)]+=cv;
    }
}

int cfd_solve_multigrid_vcycle(const CfdGrid2D *grid, const CfdBCSet2D *bc,
                                double *phi, int n_levels, int n_pre, int n_post)
{
    if (!grid||!phi||n_levels<2) return -1;
    int nx=grid->nx, ny=grid->ny, nc=nx*ny;

    CfdMatrix2D *A=cfd_matrix2d_create(nx,ny);
    if (!A) return -1;
    cfd_assemble_poisson_2d(grid,NULL,bc,A);

    /* Pre-smoothing: n_pre GS sweeps */
    for (int s=0;s<n_pre;s++) {
        for (int j=0;j<ny;j++) for (int i=0;i<nx;i++) {
            int idx=CFD_IDX2(i,j,nx);
            double sum=A->b[idx];
            if(i>0)sum+=A->aW[idx]*phi[CFD_IDX2(i-1,j,nx)];
            if(i<nx-1)sum+=A->aE[idx]*phi[CFD_IDX2(i+1,j,nx)];
            if(j>0)sum+=A->aS[idx]*phi[CFD_IDX2(i,j-1,nx)];
            if(j<ny-1)sum+=A->aN[idx]*phi[CFD_IDX2(i,j+1,nx)];
            phi[idx]=sum/A->aP[idx];
        }
    }

    double *res=malloc((size_t)nc*sizeof(double));
    if(!res){cfd_matrix2d_destroy(A);return -1;}
    cfd_compute_residual(A,phi,res);

    int nxc=nx/2, nyc=ny/2, ncc=nxc*nyc;
    double *res_c=calloc((size_t)ncc,sizeof(double));
    double *corr_c=calloc((size_t)ncc,sizeof(double));
    mg_restrict_2d(res,res_c,nx,ny);

    /* Coarse grid solve (20 GS sweeps) */
    double hx2=4.0*grid->dx*grid->dx;
    double hy2=4.0*grid->dy*grid->dy;
    for (int s=0;s<20;s++) {
        for (int j=0;j<nyc;j++) for (int i=0;i<nxc;i++) {
            int idx=CFD_IDX2(i,j,nxc);
            double sum=res_c[idx];
            if(i>0)sum+=corr_c[CFD_IDX2(i-1,j,nxc)]/hx2;
            if(i<nxc-1)sum+=corr_c[CFD_IDX2(i+1,j,nxc)]/hx2;
            if(j>0)sum+=corr_c[CFD_IDX2(i,j-1,nxc)]/hy2;
            if(j<nyc-1)sum+=corr_c[CFD_IDX2(i,j+1,nxc)]/hy2;
            corr_c[idx]=sum/(2.0/hx2+2.0/hy2);
        }
    }

    double *corr_f=calloc((size_t)nc,sizeof(double));
    mg_prolong_2d(corr_c,corr_f,nxc,nyc);
    for (int i=0;i<nc;i++) phi[i]+=corr_f[i];

    /* Post-smoothing */
    for (int s=0;s<n_post;s++) {
        for (int j=0;j<ny;j++) for (int i=0;i<nx;i++) {
            int idx=CFD_IDX2(i,j,nx);
            double sum=A->b[idx];
            if(i>0)sum+=A->aW[idx]*phi[CFD_IDX2(i-1,j,nx)];
            if(i<nx-1)sum+=A->aE[idx]*phi[CFD_IDX2(i+1,j,nx)];
            if(j>0)sum+=A->aS[idx]*phi[CFD_IDX2(i,j-1,nx)];
            if(j<ny-1)sum+=A->aN[idx]*phi[CFD_IDX2(i,j+1,nx)];
            phi[idx]=sum/A->aP[idx];
        }
    }

    free(res);free(res_c);free(corr_c);free(corr_f);
    cfd_matrix2d_destroy(A);
    return 0;
}

int cfd_solve_full_multigrid(const CfdGrid2D *grid, const CfdBCSet2D *bc,
                              double *phi, int n_levels)
{
    if (!grid||!phi||n_levels<2) return -1;
    return cfd_solve_multigrid_vcycle(grid,bc,phi,n_levels,2,2);
}

/* ===== Additional Solver Utilities ===== */
int cfd_solve_richardson(const CfdMatrix2D *A, double *x, CfdSolverState *state) {
    /* Richardson iteration: x^(k+1) = x^k + omega*(b - A*x^k) */
    if (!A || !x || !state) return -1;
    int nc = A->ncells; double omega = state->under_relaxation;
    double *r = malloc((size_t)nc * sizeof(double));
    if (!r) return -1;
    for (int iter = 0; iter < state->max_iter; iter++) {
        cfd_compute_residual(A, x, r);
        for (int i = 0; i < nc; i++) x[i] += omega * r[i];
        double rn = cfd_residual_l2(r, nc);
        if (iter == 0) state->initial_residual = rn;
        if (rn / state->initial_residual < state->tolerance) {
            state->n_iter = iter + 1; state->final_residual = rn;
            free(r); return iter + 1;
        }
    }
    state->n_iter = state->max_iter; state->diverged = 1; free(r); return -1;
}

/* Chebyshev acceleration: optimal polynomial acceleration for SPD systems. */
int cfd_solve_chebyshev(const CfdMatrix2D *A, double *x, CfdSolverState *state,
                         double lambda_min, double lambda_max) {
    if (!A || !x || !state) return -1;
    int nc = A->ncells;
    double *r = malloc((size_t)nc * sizeof(double));
    double *p = malloc((size_t)nc * sizeof(double));
    if (!r || !p) { free(r); free(p); return -1; }
    double d = (lambda_max + lambda_min) / 2.0;
    double c = (lambda_max - lambda_min) / 2.0;
    cfd_compute_residual(A, x, r);
    for (int i = 0; i < nc; i++) p[i] = r[i] / d;
    double init_rn = cfd_residual_l2(r, nc);
    for (int iter = 0; iter < state->max_iter; iter++) {
        for (int i = 0; i < nc; i++) x[i] += p[i];
        cfd_compute_residual(A, x, r);
        double rn = cfd_residual_l2(r, nc);
        if (rn / init_rn < state->tolerance) { free(r); free(p); return iter + 1; }
        double theta = (c / (2.0 * d)) * (c / (2.0 * d));
        for (int i = 0; i < nc; i++) p[i] = r[i] / d + theta * p[i];
    }
    free(r); free(p); return -1;
}

/* Estimate spectral radius via power iteration for convergence analysis. */
double cfd_estimate_spectral_radius(const CfdMatrix2D *A, int n_power_iters) {
    if (!A || n_power_iters < 1) return 0.0;
    int nc = A->ncells;
    double *v = malloc((size_t)nc * sizeof(double));
    double *Av = malloc((size_t)nc * sizeof(double));
    if (!v || !Av) { free(v); free(Av); return 0.0; }
    for (int i = 0; i < nc; i++) v[i] = 1.0;
    double lambda = 0.0;
    for (int iter = 0; iter < n_power_iters; iter++) {
        cfd_compute_residual(A, v, Av);
        double norm = 0.0;
        for (int i = 0; i < nc; i++) { norm += Av[i] * Av[i]; }
        norm = sqrt(norm);
        if (norm > 1e-15) {
            for (int i = 0; i < nc; i++) v[i] = Av[i] / norm;
            lambda = norm;
        }
    }
    free(v); free(Av); return lambda;
}

/* Compute diagonal preconditioner matrix M^{-1}. */
void cfd_preconditioner_diag(const CfdMatrix2D *A, double *M_inv) {
    if (!A || !M_inv) return;
    for (int i = 0; i < A->ncells; i++)
        M_inv[i] = (fabs(A->aP[i]) > 1e-15) ? 1.0 / A->aP[i] : 0.0;
}

/* Apply SSOR preconditioner sweep (forward then backward). */
void cfd_preconditioner_ssor(const CfdMatrix2D *A, const double *r, double *z, double omega) {
    if (!A || !r || !z) return;
    int nc = A->ncells, nx = A->nx, ny = A->ny;
    for (int i = 0; i < nc; i++) z[i] = r[i];
    for (int j = 0; j < ny; j++) for (int i = 0; i < nx; i++) {
        int idx = CFD_IDX2(i, j, nx);
        double sum = 0.0;
        if (i > 0) sum += A->aW[idx] * z[CFD_IDX2(i-1, j, nx)];
        if (j > 0) sum += A->aS[idx] * z[CFD_IDX2(i, j-1, nx)];
        z[idx] = omega * (z[idx] - sum) / A->aP[idx] + (1.0-omega) * z[idx];
    }
    for (int j = ny-1; j >= 0; j--) for (int i = nx-1; i >= 0; i--) {
        int idx = CFD_IDX2(i, j, nx);
        double sum = 0.0;
        if (i < nx-1) sum += A->aE[idx] * z[CFD_IDX2(i+1, j, nx)];
        if (j < ny-1) sum += A->aN[idx] * z[CFD_IDX2(i, j+1, nx)];
        z[idx] = omega * (z[idx] - sum) / A->aP[idx] + (1.0-omega) * z[idx];
    }
}

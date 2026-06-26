#ifndef CFD_SOLVERS_H
#define CFD_SOLVERS_H
#include "cfd_types.h"
#ifdef __cplusplus
extern "C" {
#endif
int cfd_solve_jacobi(const CfdMatrix2D *A, double *x, CfdSolverState *state);
int cfd_solve_gauss_seidel(const CfdMatrix2D *A, double *x, CfdSolverState *state);
int cfd_solve_sor(const CfdMatrix2D *A, double *x, double omega, CfdSolverState *state);
int cfd_solve_tdma(const CfdMatrix2D *A, double *x, int nx);
int cfd_solve_cg(const CfdMatrix2D *A, double *x, CfdSolverState *state);
int cfd_solve_pcg_diag(const CfdMatrix2D *A, double *x, CfdSolverState *state);
int cfd_solve_bicgstab(const CfdMatrix2D *A, double *x, CfdSolverState *state);
int cfd_solve_lu_banded(double *A_band, double *b, double *x, int n, int bandwidth);
int cfd_solve_multigrid_vcycle(const CfdGrid2D *grid, const CfdBCSet2D *bc, double *phi, int n_levels, int n_pre, int n_post);
int cfd_solve_full_multigrid(const CfdGrid2D *grid, const CfdBCSet2D *bc, double *phi, int n_levels);
void cfd_compute_residual(const CfdMatrix2D *A, const double *x, double *r);
double cfd_residual_l2(const double *r, int n);
double cfd_residual_linf(const double *r, int n);
#ifdef __cplusplus
}
#endif
#endif

/*
 * Linear Solver Theory
 * =====================
 *
 * Jacobi: x_i^(k+1) = (b_i - sum_{j!=i} a_ij*x_j^(k)) / a_ii
 *   Convergence rate: rho(D^{-1}*(L+U)). Grid-dependent.
 *   For 2D Poisson (N x N): rho = cos(pi/(N+1)) ~ 1 - O(1/N^2)
 *
 * Gauss-Seidel: Uses latest values as they become available.
 *   Convergence: ~2x faster than Jacobi (rho_GS = rho_J^2).
 *
 * SOR: Extrapolated Gauss-Seidel with omega in (0,2).
 *   Optimal omega for Poisson: omega_opt = 2/(1+sqrt(1-rho_J^2))
 *
 * Thomas TDMA: O(N) direct solver for tridiagonal systems.
 *   Forward sweep: c_i = c_i/(b_i - a_i*c_{i-1})
 *                 d_i = (d_i - a_i*d_{i-1})/(b_i - a_i*c_{i-1})
 *   Back substitution: x_N = d_N, x_i = d_i - c_i*x_{i+1}
 *
 * Conjugate Gradient: Optimal Krylov method for SPD systems.
 *   Convergence: ||e_k||_A <= 2*((sqrt(kappa)-1)/(sqrt(kappa)+1))^k * ||e_0||_A
 *   Preconditioned CG: M^{-1}*A*x = M^{-1}*b for faster convergence.
 *
 * BiCGSTAB: For non-symmetric systems. Van der Vorst (1992).
 *   Two-term recurrence with stabilization step.
 *
 * Multigrid: Optimal O(N) solver for elliptic problems.
 *   V-cycle: restriction -> coarse solve -> prolongation.
 *   Smoothing: Gauss-Seidel eliminates high-frequency errors.
 *   Coarse grid: Low-frequency errors become high-frequency.
 */

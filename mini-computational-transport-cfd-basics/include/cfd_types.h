#ifndef CFD_TYPES_H
#define CFD_TYPES_H
#include "cfd_config.h"
#include <stddef.h>
#include <stdint.h>

typedef enum { CFD_GRID_1D=1, CFD_GRID_2D=2, CFD_GRID_3D=3 } CfdGridDim;
typedef enum { CFD_TOPOLOGY_STRUCTURED=0, CFD_TOPOLOGY_BLOCK_STRUCT=1, CFD_TOPOLOGY_UNSTRUCTURED=2, CFD_TOPOLOGY_OVERSET=3 } CfdGridTopology;
typedef enum { CFD_CELL_LINE=0, CFD_CELL_QUAD=1, CFD_CELL_TRI=2, CFD_CELL_HEX=3, CFD_CELL_TET=4, CFD_CELL_WEDGE=5, CFD_CELL_PYRAMID=6 } CfdCellShape;
typedef enum { CFD_LOC_CENTER=0, CFD_LOC_FACE_X=1, CFD_LOC_FACE_Y=2, CFD_LOC_FACE_Z=3, CFD_LOC_NODE=4 } CfdStorageLocation;
typedef enum { CFD_FIELD_VELOCITY_X=0, CFD_FIELD_VELOCITY_Y=1, CFD_FIELD_VELOCITY_Z=2, CFD_FIELD_PRESSURE=3, CFD_FIELD_TEMPERATURE=4, CFD_FIELD_CONCENTRATION=5, CFD_FIELD_TKE=6, CFD_FIELD_EPSILON=7, CFD_FIELD_OMEGA=8, CFD_FIELD_VISCOSITY=9, CFD_FIELD_STREAMFUNC=10, CFD_FIELD_VORTICITY=11 } CfdFieldType;
typedef enum { CFD_BC_INLET_VELOCITY=0, CFD_BC_INLET_PRESSURE=1, CFD_BC_OUTLET_PRESSURE=2, CFD_BC_OUTLET_OUTFLOW=3, CFD_BC_WALL_NO_SLIP=4, CFD_BC_WALL_FREE_SLIP=5, CFD_BC_WALL_MOVING=6, CFD_BC_WALL_TEMP=7, CFD_BC_WALL_HEAT_FLUX=8, CFD_BC_WALL_ADIABATIC=9, CFD_BC_SYMMETRY=10, CFD_BC_PERIODIC=11, CFD_BC_FAR_FIELD=12 } CfdBCType;

typedef struct { int nx,ny,nz,ncells,nnodes; double Lx,Ly,Lz,dx,dy,dz; double *x,*y,*z,*xc,*yc,*zc,*cell_vol; } CfdGrid3D;
typedef struct { int nx,ny,ncells,nnodes; double Lx,Ly,dx,dy; double *x,*y,*xc,*yc,*cell_vol,*face_area_x,*face_area_y; CfdGridTopology topology; CfdCellShape cell_shape; } CfdGrid2D;
typedef struct { int nx,nnodes; double L,dx; double *x,*xc; } CfdGrid1D;

typedef struct { int nx,ny,ncells; double *data; CfdFieldType ftype; CfdStorageLocation loc; char name[32]; } CfdScalarField2D;
typedef struct { int ncells; double *data; CfdFieldType ftype; CfdStorageLocation loc; char name[32]; } CfdScalarField1D;
typedef struct { int nx,ny,ncells; double *u,*v,*p; CfdStorageLocation u_loc,v_loc,p_loc; char name[32]; } CfdVectorField2D;

typedef struct { CfdBCType type; double value,gradient,htc,T_ref,wall_velocity; } CfdBCSingle;
typedef struct { CfdBCSingle west,east,south,north; } CfdBCSet2D;
typedef struct { CfdBCSingle left,right; } CfdBCSet1D;

typedef struct { int nx,ny,ncells; double *aW,*aS,*aP,*aE,*aN,*b; } CfdMatrix2D;

typedef struct { int iter; double residual,delta_phi,residual_ratio; } CfdSolverStep;
typedef struct { CfdSolverType solver_type; double tolerance,initial_residual,final_residual,final_delta,under_relaxation; int max_iter,n_iter,diverged; CfdSolverStep *history; int history_capacity,history_count; } CfdSolverState;

typedef struct { double density,diffusivity,source_constant,source_linear,cfl_number; CfdConvectiveScheme conv_scheme; } CfdTransportEquation;
typedef struct { double density,viscosity,kinematic_visc,thermal_cond,specific_heat_cp,prandtl,thermal_expansion,mass_diffusivity; } CfdFluidProperty;

typedef struct { double Reynolds,cell_Reynolds,Peclet,cell_Peclet,Prandtl,Schmidt,Nusselt,Sherwood,Grashof,Rayleigh,Courant,Fourier_diff,Stanton,Mach,Weber,Froude; } CfdDimensionless;
typedef struct { double L1_norm,L2_norm,L_inf_norm,relative_L2,mass_imbalance; } CfdResidualNorms;
typedef struct { double abs_tolerance,rel_tolerance,div_tolerance,stall_tolerance; int min_iter,stall_window; } CfdConvergenceCriterion;
typedef struct { int max_outer_iter,max_inner_iter,n_outer; double momentum_tol,pressure_tol,urf_u,urf_v,urf_p,mass_residual,momentum_residual; } CfdSimpleState;
typedef struct { double C_mu,C1_eps,C2_eps,sigma_k,sigma_eps,k_min,eps_min; } CfdKEpsilonModel;
typedef struct { double beta_star,beta1,beta2,sigma_k1,sigma_k2,sigma_omega1,sigma_omega2,a1,k_min,omega_min; } CfdKOmegaSSTModel;

#ifdef __cplusplus
extern "C" {
#endif
CfdGrid1D *cfd_grid1d_create(int nx, double L);
void cfd_grid1d_destroy(CfdGrid1D *grid);
CfdGrid2D *cfd_grid2d_create(int nx, int ny, double Lx, double Ly);
void cfd_grid2d_destroy(CfdGrid2D *grid);
int cfd_grid2d_allocate_all(CfdGrid2D *grid);
CfdGrid3D *cfd_grid3d_create(int nx, int ny, int nz, double Lx, double Ly, double Lz);
void cfd_grid3d_destroy(CfdGrid3D *grid);
CfdScalarField1D *cfd_scalar1d_create(int ncells, CfdFieldType ftype);
void cfd_scalar1d_destroy(CfdScalarField1D *field);
int cfd_scalar1d_fill(CfdScalarField1D *field, double value);
CfdScalarField2D *cfd_scalar2d_create(int nx, int ny, CfdFieldType ftype);
void cfd_scalar2d_destroy(CfdScalarField2D *field);
int cfd_scalar2d_copy(const CfdScalarField2D *src, CfdScalarField2D *dst);
CfdVectorField2D *cfd_vector2d_create(int nx, int ny);
void cfd_vector2d_destroy(CfdVectorField2D *field);
CfdMatrix2D *cfd_matrix2d_create(int nx, int ny);
void cfd_matrix2d_destroy(CfdMatrix2D *mat);
void cfd_matrix2d_zero(CfdMatrix2D *mat);
void cfd_solver_state_init(CfdSolverState *state, CfdSolverType stype, double tol, int max_iter);
int cfd_solver_converged(const CfdSolverState *state);
int cfd_bc1d_default(CfdBCSet1D *bc);
int cfd_bc2d_default(CfdBCSet2D *bc);
void cfd_fill_dimensionless_full(const CfdFluidProperty *fp, double U_ref, double L_ref, double dx_cell, double dt_step, CfdDimensionless *dless);
int cfd_fluid_air(CfdFluidProperty *fp);
int cfd_fluid_water(CfdFluidProperty *fp);
int cfd_fluid_oil(CfdFluidProperty *fp);
int cfd_fluid_mercury(CfdFluidProperty *fp);
void cfd_residual_norms_init(CfdResidualNorms *r);
void cfd_compute_residual_norms(const double *residual, int n, const double *ref, CfdResidualNorms *norms);
#ifdef __cplusplus
}
#endif
#endif

/*
 * ============================================================
 * L1: Core Definitions - Extended Documentation
 * ============================================================
 *
 * Grid Terminology (after Thompson, Warsi & Mastin 1985):
 *   - Structured grid: i,j,k indexing with implicit connectivity.
 *   - Unstructured grid: Explicit node/cell connectivity lists.
 *   - Block-structured: Domain decomposed into structured blocks.
 *   - Overset/Chimera: Overlapping grids with hole-cutting.
 *
 * Cell Shapes:
 *   1D: Line (2 nodes)
 *   2D: Triangle (3), Quadrilateral (4)
 *   3D: Tetrahedron (4), Hexahedron (8), Wedge/Prism (6), Pyramid (5)
 *
 * Field Storage Locations:
 *   - Collocated (cell-centered): All variables at cell centers.
 *     Requires Rhie-Chow interpolation for pressure-velocity coupling.
 *   - Staggered (MAC, Harlow & Welch 1965): Velocity at faces,
 *     pressure at cell centers. Natural pressure-velocity coupling.
 *
 * Boundary Condition Types (after Hirsch 2007):
 *   - Inlet: Velocity or total pressure specified.
 *   - Outlet: Static pressure or fully-developed (zero gradient).
 *   - Wall: No-slip (u=0), free-slip (zero shear), moving wall.
 *   - Symmetry: Zero normal velocity, zero normal gradients.
 *   - Periodic: phi(x+L) = phi(x).
 *   - Far-field: Characteristic-based non-reflecting BCs.
 */

/*
 * L2: Coefficient Matrix Conventions
 * ==================================
 * The penta-diagonal system for 2D structured FVM:
 *   A_P * phi_P = A_W * phi_W + A_E * phi_E + A_S * phi_S + A_N * phi_N + b
 *
 * Where A_W, A_E, A_S, A_N are neighbor coefficients (non-negative for DMP).
 * A_P = A_W + A_E + A_S + A_N - S_P * Vol (diagonal dominance for stability).
 * b = S_C * Vol + boundary contributions.
 *
 * For 1D: only A_W, A_P, A_E are non-zero.
 * For tri-diagonal: A_S = A_N = 0.
 */

/*
 * L3: Dimensionless Groups in CFD
 * ===============================
 * Re    = rho*U*L/mu         : Inertia / Viscous forces
 * Pe    = Re*Pr               : Convection / Diffusion (thermal)
 * Pr    = nu/alpha            : Momentum / Thermal diffusivity
 * Sc    = nu/D                : Momentum / Mass diffusivity
 * Nu    = h*L/k               : Convection / Conduction
 * Gr    = g*beta*dT*L^3/nu^2  : Buoyancy / Viscous forces
 * Ra    = Gr*Pr               : Buoyancy-driven convection strength
 * CFL   = U*dt/dx             : Numerical stability parameter
 * Fo    = alpha*dt/dx^2       : Diffusion stability parameter
 * Ma    = U/c                 : Compressibility (Ma<0.3=incompressible)
 * We    = rho*U^2*L/sigma     : Inertia / Surface tension
 * Fr    = U/sqrt(g*L)         : Inertia / Gravity
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cfd_types.h"
#include "cfd_discretization.h"
#include "cfd_solvers.h"
#include "cfd_convection_diffusion.h"

int main(void) {
    printf("=== Lid-Driven Cavity (SIMPLE) ===
");
    printf("Classic benchmark: square cavity with moving top wall

");

    int nx = 16, ny = 16;
    double L = 1.0, U_lid = 1.0;

    CfdGrid2D *grid = cfd_grid2d_create(nx, ny, L, L);
    cfd_grid2d_allocate_all(grid);
    CfdVectorField2D *vel = cfd_vector2d_create(nx, ny);
    CfdFluidProperty fp;
    cfd_fluid_water(&fp);

    /* Adjust for benchmark: use Re=100 */
    double nu = 0.01;
    fp.viscosity = fp.density * nu;
    fp.kinematic_visc = nu;

    double Re = U_lid * L / nu;
    printf("Reynolds number Re = %.0f
", Re);
    printf("Grid: %dx%d, cells: %d
", nx, ny, nx*ny);

    CfdSimpleState simple;
    simple.max_outer_iter = 200;
    simple.max_inner_iter = 50;
    simple.momentum_tol = 1e-4;
    simple.pressure_tol = 1e-5;
    simple.urf_u = 0.5;
    simple.urf_v = 0.5;
    simple.urf_p = 0.3;

    int iter = cfd_solve_lid_driven_cavity(grid, vel, &fp, U_lid, &simple);
    printf("SIMPLE iterations: %d
", iter);

    /* Print u-velocity along vertical centerline */
    printf("
u-velocity at x=0.5:
");
    int i_mid = nx / 2;
    for (int j = 0; j < ny; j++) {
        int idx = CFD_IDX2(i_mid, j, nx);
        double y = grid->yc[idx];
        printf("  y=%.4f  u=%.6f  v=%.6f
", y, vel->u[idx], vel->v[idx]);
    }

    cfd_vector2d_destroy(vel);
    cfd_grid2d_destroy(grid);
    printf("
Done.
");
    return 0;
}

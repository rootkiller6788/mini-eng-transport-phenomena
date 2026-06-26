#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cfd_types.h"
#include "cfd_discretization.h"
#include "cfd_solvers.h"
#include "cfd_convection_diffusion.h"

int main(void) {
    printf("=== 2D Poisson Equation ===
");
    printf("Equation: -del^2(phi) = 0 (Laplace)
");
    printf("Domain: [0,1]x[0,1]
");
    printf("BC: phi=0 at x=0,y=0; phi=1 at x=1; phi=0 at y=1

");

    int nx = 20, ny = 20;
    CfdGrid2D *grid = cfd_grid2d_create(nx, ny, 1.0, 1.0);
    cfd_grid2d_allocate_all(grid);

    CfdScalarField2D *phi = cfd_scalar2d_create(nx, ny, CFD_FIELD_TEMPERATURE);
    CfdBCSet2D bc;
    cfd_bc2d_default(&bc);
    bc.west.value = 0.0;
    bc.east.value = 1.0;
    bc.south.value = 0.0;
    bc.north.value = 0.0;

    cfd_solve_poisson_2d(grid, NULL, &bc, phi);

    printf("Solution at mid-plane (y=0.5):
");
    int j_mid = ny / 2;
    for (int i = 0; i < nx; i += 2) {
        int idx = CFD_IDX2(i, j_mid, nx);
        double x = grid->xc[idx];
        printf("  x=%.3f  phi=%.6f
", x, phi->data[idx]);
    }

    cfd_scalar2d_destroy(phi);
    cfd_grid2d_destroy(grid);
    printf("
Done.
");
    return 0;
}

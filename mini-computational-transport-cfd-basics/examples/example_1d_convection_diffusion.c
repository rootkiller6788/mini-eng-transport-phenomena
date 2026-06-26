#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cfd_types.h"
#include "cfd_discretization.h"
#include "cfd_solvers.h"
#include "cfd_convection_diffusion.h"

int main(void) {
    printf("=== 1D Steady Convection-Diffusion ===
");
    printf("Equation: u*d(phi)/dx = Gamma*d^2(phi)/dx^2
");
    printf("BC: phi(0)=0, phi(L)=1

");

    int nx = 40;
    double L = 1.0;
    double u_vel = 0.5;
    double Gamma = 0.05;
    double Pe = u_vel * L / Gamma;

    CfdGrid1D *grid = cfd_grid1d_create(nx, L);
    CfdTransportEquation eqn = {1.0, Gamma, 0.0, 0.0, CFD_SCHEME_POWERLAW, 0.5};
    double *phi = calloc((size_t)nx, sizeof(double));

    printf("Peclet number Pe = %.2f
", Pe);
    printf("Cell Peclet Pe_cell = %.2f
", u_vel*grid->dx/Gamma);

    /* Solve with power-law scheme */
    cfd_solve_convdiff_1d(grid, &eqn, u_vel, 0.0, 1.0, CFD_SCHEME_POWERLAW, phi);

    printf("
%-8s %-12s %-12s %-12s
", "x", "phi(FVM)", "phi(exact)", "error");
    for (int i = 0; i < nx; i += 4) {
        double xc = grid->xc[i];
        double exact = cfd_convdiff_1d_exact(xc, L, Pe, 0.0, 1.0);
        double err = fabs(phi[i] - exact);
        printf("%-8.3f %-12.6f %-12.6f %-12.2e
", xc, phi[i], exact, err);
    }

    free(phi); cfd_grid1d_destroy(grid);
    printf("
Done.
");
    return 0;
}

# mini-computational-transport-cfd-basics

**Computational Transport CFD Basics**

Part of mini-engineering-physics / Transport Phenomena module.

---

## Module Status: COMPLETE

- **L1 Definitions:** Complete (8 items)
- **L2 Core Concepts:** Complete (8 items)
- **L3 Engineering Quantities:** Complete (6 items)
- **L4 Conservation Laws:** Complete (6 items)
- **L5 Engineering Methods:** Complete (12 items)
- **L6 Engineering Problems:** Complete (5 items)
- **L7 Applications:** Partial+ (4 applications)
- **L8 Advanced Methods:** Complete (8 items)
- **L9 Research Frontiers:** Partial (7 items documented)

**Code: 3301 lines (include/ + src/)** | **Tests: 13/13 passed** | **Examples: 3**

---

## Core Definitions (L1)

| Term | Symbol | Definition |
|------|--------|------------|
| Grid | - | Structured 1D/2D/3D mesh |
| Scalar field | phi(x) | Cell-centered field data |
| Vector field | (u,v) | Velocity field on 2D grid |
| Coefficient matrix | A | Penta-diagonal FVM matrix |
| BC types | - | Dirichlet, Neumann, Robin, etc. |
| Convection scheme | - | 7 discretization schemes |
| Solver type | - | 8 solver types |
| Fluid properties | rho,mu,k,cp | Thermophysical container |

## Core Theorems (L4)

| Theorem | Formula | Implementation |
|---------|---------|---------------|
| Central diff | du/dx = (u_{i+1}-u_{i-1})/(2dx) + O(dx^2) | cfd_central_diff_1d() |
| 2nd derivative | d^2u/dx^2 = (u_{i+1}-2u_i+u_{i-1})/dx^2 | cfd_second_deriv_1d() |
| Power-law | A(|Pe|) = max(0,(1-0.1|Pe|)^5) | cfd_powerlaw_1d() |
| DMP | aP >= aW+aE+aS+aN | cfd_check_maximum_principle() |
| CFL condition | dt <= CFL*dx/|u| | cfd_stable_timestep_1d() |
| Diffusion stability | dt <= dx^2/(2*Gamma) | cfd_stable_timestep_1d() |

## Core Algorithms (L5)

| Algorithm | Complexity | Implementation |
|-----------|-----------|---------------|
| Jacobi | O(N^2)/iter | cfd_solve_jacobi() |
| Gauss-Seidel | O(N^2)/iter | cfd_solve_gauss_seidel() |
| SOR | O(N^2)/iter | cfd_solve_sor() |
| Thomas TDMA | O(N) | cfd_solve_tdma() |
| Conjugate Gradient | O(N*sqrt(kappa)) | cfd_solve_cg() |
| PCG (diagonal) | O(N*sqrt(kappa)) | cfd_solve_pcg_diag() |
| BiCGSTAB | O(N*sqrt(kappa)) | cfd_solve_bicgstab() |
| Multigrid V-cycle | O(N) | cfd_solve_multigrid_vcycle() |
| SIMPLE | O(N*iters) | cfd_simple_solve() |
| k-epsilon eddy visc | O(1) | cfd_eddy_viscosity_kepsilon() |
| Colebrook-White | O(iter) | cfd_colebrook_friction_factor() |

## Classic Problems (L6)

1. 1D steady convection-diffusion
2. 2D Laplace/Poisson equation
3. Lid-driven cavity (SIMPLE, Re=100)
4. 1D unsteady convection-diffusion
5. Grid convergence study (Richardson, GCI)

## Engineering Applications (L7)

1. Electronics cooling (flat plate convection)
2. HVAC duct flow analysis
3. Chemical reactor transport
4. Heat exchanger sizing

## Nine-School Curriculum Mapping

| School | Course | Topics |
|--------|--------|--------|
| MIT | 2.005, 2.25 | FVM, SIMPLE |
| Stanford | ME 346A, 351A | CFD, turb |
| Berkeley | ME 105/106 | Num trans |
| Michigan | ME 320/420 | Pipe CFD |
| Purdue | ME 505/509 | Discretiz |
| TU Munich | MW 0854 | CFD heat |
| ETH | 151-0111/0103 | Num meth |
| Tsinghua | Eng Therm/Hyd | CFD fund |

---

## Build and Test



## File Structure



## References

- Ferziger, Peric & Street (2020) Computational Methods for Fluid Dynamics, 4th ed.
- Patankar (1980) Numerical Heat Transfer and Fluid Flow
- Versteeg & Malalasekera (2007) An Introduction to CFD: The FVM, 2nd ed.
- LeVeque (2002) Finite Volume Methods for Hyperbolic Problems
- Bird, Stewart, Lightfoot (2007) Transport Phenomena, 2nd ed.
- Launder & Spalding (1974) The numerical computation of turbulent flows
- Menter (1994) Two-equation eddy-viscosity turbulence models
- Roache (1994) Perspective: A method for uniform reporting of grid refinement studies

## License

MIT

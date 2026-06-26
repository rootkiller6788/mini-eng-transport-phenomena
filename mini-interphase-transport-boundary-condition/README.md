# mini-interphase-transport-boundary-condition

**Interphase Transport and Boundary Conditions**

Part of `mini-engineering-physics` / Transport Phenomena module.

---

## Module Status: COMPLETE ✅

- **L1 Definitions:** Complete (17 items)
- **L2 Core Concepts:** Complete (17 items)
- **L3 Engineering Quantities:** Complete (10 items)
- **L4 Conservation Laws:** Complete (10 items)
- **L5 Engineering Methods:** Complete (28 items)
- **L6 Engineering Problems:** Complete (8 items)
- **L7 Applications:** Complete (3 applications)
- **L8 Advanced Methods:** Complete (6 items)
- **L9 Research Frontiers:** Partial (3 items documented)

**Code: 4911 lines (include/ + src/)** | **Tests: 31/31 passed** | **Examples: 3**

---

## Core Definitions (L1)

- **Phase types**: Gas, Liquid, Solid, Plasma, Supercritical fluid
- **Boundary condition taxonomy**: Dirichlet, Neumann, Robin, Jump, Periodic, Symmetry, Outflow, Cauchy (after Hadamard, 1902)
- **Interface geometry**: Planar, Cylindrical, Spherical, Ellipsoidal, Arbitrary, Corrugated, Fractal
- **Phase pairs**: Gas-Liquid, Liquid-Liquid, Gas-Solid, Liquid-Solid, Liquid-Vapor, Solid-Solid, Gas-Plasma, Immiscible Liquids
- **Interface descriptor**: Complete thermodynamic state on both sides, surface tension, curvature, roughness, slip length
- **Flux vectors**: Stress tensor (9-component), Heat flux, Mass flux, Combined interphase flux

## Core Theorems (L4)

| Theorem | Formula | Implementation |
|---------|---------|---------------|
| Mass jump condition | [[ρ(v - v_I)·n]] = 0 | `jump_mass()` |
| Momentum jump condition | [[ρv(v-v_I)·n + pI·n - τ·n]] = 2Hσn + ∇_s(σ) | `jump_momentum()` |
| Energy jump condition | k_A(∂T/∂n)_A - k_B(∂T/∂n)_B = ṁ·h_fg | `jump_energy()` |
| Stefan condition | k_s(∂T_s/∂n) - k_l(∂T_l/∂n) = ρ_s·L·v_n | `stefan_condition()` |
| Young-Laplace | Δp = 2σH | `interface_young_laplace()` |
| Hertz-Knudsen-Schrage | ṁ = (2α/(2-α))·√(M/(2πR))·(p_sat/√T_l - p_v/√T_v) | `hertz_knudsen_schrage()` |
| Entropy 2nd Law at interface | [[ρs(v-v_I)·n + q·n/T]] ≥ 0 | `jump_entropy()` |
| Henry's Law | p_i = H·x_i | `interface_henry_equilibrium()` |
| Raoult's Law | p_i = x_i·γ_i·p_i^sat | `interface_raoult_equilibrium()` |

## Core Algorithms (L5)

| Algorithm | Source | Implementation |
|-----------|--------|---------------|
| Two-film mass transfer | Whitman (1923) | `two_film_mass_transfer()` |
| Penetration theory | Higbie (1935) | `penetration_mass_transfer()` |
| Surface renewal | Danckwerts (1951) | `surface_renewal_mass_transfer()` |
| Chapman-Enskog diffusion | Chapman & Enskog (1917) | `chapman_enskog_diffusivity()` |
| Dittus-Boelter correlation | Dittus & Boelter (1930) | `nusselt_pipe_turbulent()` |
| Gnielinski correlation | Gnielinski (1976) | `nusselt_pipe_turbulent()` |
| Churchill-Chu correlation | Churchill & Chu (1975) | `nusselt_natural_convection_vertical()` |
| Nusselt film condensation | Nusselt (1916) | `nusselt_film_condensation()` |
| Rohsenow pool boiling | Rohsenow (1952) | `rohsenow_pool_boiling()` |
| Froessling sphere mass transfer | Froessling (1938) | `sherwood_single_sphere()` |
| PLIC VOF reconstruction | Youngs (1982) | `plic_vof_fraction_2d()` |
| Level-set reinitialization | Sussman et al. (1994) | `level_set_reinitialize()` |
| CSF surface tension | Brackbill et al. (1992) | `csf_surface_tension_2d()` |
| Ghost Fluid Method | Fedkiw et al. (1999) | `gfm_dirichlet()`, `gfm_neumann()`, `gfm_robin()` |
| Immersed Boundary Method | Peskin (1977), Roma et al. (1999) | `ibm_spread_force()`, `ibm_interpolate_velocity()` |

## Classic Problems (L6)

1. **Falling film gas absorption** — CO2 absorption in vertical wetted-wall column
2. **Evaporative cooling** — Wet-bulb temperature via simultaneous heat/mass transfer
3. **Distillation tray efficiency** — Murphree vapor efficiency from mass transfer
4. **Membrane contactor sizing** — Hollow-fiber module for gas separation
5. **Condensation on horizontal tubes** — Power plant condenser design
6. **Droplet evaporation** — d^2-law (Godsave/Spalding, 1953)
7. **Bubble column mass transfer** — Interfacial area in gas-liquid reactors
8. **Dimensionless group analysis** — Complete Buckingham Pi for interphase transport

## Engineering Applications (L7)

1. **Wastewater aeration design** (Metcalf & Eddy / ASCE standard) — `wastewater_aeration_design()`
2. **PWR nuclear steam generator sizing** (Westinghouse AP1000 reference) — `pwr_steam_generator_sizing()`
3. **PEM fuel cell water management** (Toyota Mirai, Springer model) — `pem_fuel_cell_water_flux()`

## Advanced Methods (L8)

1. **Ghost Fluid Method** — Sharp interface BC enforcement for Dirichlet, Neumann, and Robin conditions
2. **Immersed Boundary Method** — Roma et al. (1999) 3-point regularized delta function with force spreading and velocity interpolation
3. **Level-Set Method** — Signed distance function, advection (upwind), reinitialization (Sussman)
4. **Volume of Fluid** — PLIC reconstruction and geometric flux computation
5. **CSF Surface Tension** — Continuum surface force model (Brackbill et al., 1992)
6. **Interface property interpolation** — Smooth Heaviside/delta regularization across diffuse interface

## Nine-School Curriculum Mapping

| School | Key Course | Module Section |
|--------|-----------|---------------|
| MIT | 2.005 Thermal-Fluids | Ch. 2, 10, 19, 22 |
| Stanford | ME 346A Heat Transfer | L5 condensation/boiling, L6 droplet |
| Berkeley | ME 105/106 Thermo/Fluids | L2 VLE/Henry, L2 BCs |
| Michigan | ME 320/420 Fluids/Heat | L1 stress tensor, L4 momentum jump |
| Purdue | ME 505/509 Heat/Fluids | L5 correlations (Gnielinski, Churchill) |
| TU Munich | MW 0854 Heat Transfer | L3 Chapman-Enskog, Chung, Eucken |
| ETH | 151-0111/0103 Heat/Fluids | L4 jump conditions, L5 numerical |
| Tsinghua | 工程热力学/传热学 | L5 boundary layer, L5 phase change |

Reference: Bird, Stewart, Lightfoot (2007) "Transport Phenomena" (2nd ed.)

---

## Build & Test

```bash
make test          # Build and run all tests (31/31)
make examples      # Build all examples (3)
make run-examples  # Run all examples
make clean         # Remove build artifacts
make lines         # Line count statistics
```

## File Structure

```
mini-interphase-transport-boundary-condition/
├── README.md
├── Makefile
├── include/
│   ├── interphase_config.h            # Math constants, universal constants
│   ├── interphase_types.h             # Core data structures (L1-L3)
│   ├── interphase_boundary_conditions.h # BC API (L1-L2, L4)
│   ├── interphase_transport.h         # Transport coefficients (L2-L3, L5)
│   ├── interphase_jump_conditions.h   # Jump conditions (L4-L5)
│   └── interphase_numerical.h         # Numerical methods (L5, L8)
├── src/
│   ├── interphase_types.c             # 666 lines
│   ├── interphase_boundary_conditions.c # 357 lines
│   ├── interphase_transport.c         # 382 lines
│   ├── interphase_jump_conditions.c   # 344 lines
│   ├── interphase_numerical.c         # 656 lines
│   └── interphase_applications.c      # 645 lines
├── tests/
│   └── test_interphase.c              # 31 tests, all passing
├── examples/
│   ├── example_falling_film.c         # CO2 absorption
│   ├── example_condenser.c            # Power plant condenser
│   └── example_droplet_evaporation.c  # d^2-law evaporation
├── docs/
│   ├── knowledge-graph.md             # L1-L9 knowledge coverage
│   ├── coverage-report.md             # Detailed coverage assessment
│   ├── gap-report.md                  # Gaps and priorities
│   ├── course-alignment.md            # 9-school curriculum mapping
│   └── course-tree.md                 # Prerequisite dependency tree
├── demos/                             # (reserved for visualization)
└── benches/                           # (reserved for benchmarks)
```

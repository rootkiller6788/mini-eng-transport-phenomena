# Mini Momentum-Heat-Mass Analogy

**Transport Phenomena — The Analogy Between Momentum, Heat, and Mass Transfer**

> "The governing equations for momentum, heat, and mass transfer are
>  mathematically identical. This is the central insight of transport
>  phenomena, and it enables engineers to predict heat and mass transfer
>  from simple friction measurements."
>  — Bird, Stewart, Lightfoot (1960)

## Module Status: COMPLETE ✅

- **L1 Definitions**: Complete (20 definitions, C structs + Lean structures)
- **L2 Core Concepts**: Complete (18 concepts implemented)
- **L3 Engineering Quantities**: Complete (12-fluid database, all benchmarks)
- **L4 Conservation Laws**: Complete (13 laws/theorems, C + Lean dual verification)
- **L5 Engineering Methods**: Complete (15 methods)
- **L6 Engineering Problems**: Complete (3 examples + 7 benchmark problems)
- **L7 Applications**: Complete (5 application functions)
- **L8 Advanced Methods**: Complete (5 advanced topics)
- **L9 Research Frontiers**: Partial (documented)

**Code Metrics:**
- `include/` + `src/` total: **5,219 lines** (minimum 3,000)
- Headers: 6 files, 1,871 lines
- C Sources: 6 files, 2,632 lines
- Lean Formalization: 1 file, 716 lines
- Tests: 2 files, 45+ test cases
- Examples: 3 end-to-end >30 lines each

## Core Definitions

| Term | Symbol | Definition | Units |
|------|--------|-----------|-------|
| Dynamic viscosity | μ | τ_yx = -μ·dv_x/dy | Pa·s |
| Kinematic viscosity | ν | ν = μ/ρ | m²/s |
| Thermal conductivity | k | q_y = -k·dT/dy | W/(m·K) |
| Thermal diffusivity | α | α = k/(ρ·Cp) | m²/s |
| Mass diffusivity | D_AB | J_Ay = -D_AB·dC_A/dy | m²/s |
| Reynolds number | Re | ρ·v·L/μ | — |
| Prandtl number | Pr | ν/α = μ·Cp/k | — |
| Schmidt number | Sc | ν/D_AB | — |
| Nusselt number | Nu | h·L/k | — |
| Sherwood number | Sh | k_c·L/D_AB | — |
| Stanton number | St | Nu/(Re·Pr) = h/(ρ·v·Cp) | — |
| Skin friction coeff. | Cf | 2·τ_w/(ρ·U∞²) | — |

## Core Theorems

### 1. Reynolds Analogy (1874)
```
When Pr = Sc = 1 (ν = α = D):
    St = f_F = Cf/2
    Nu = (f/2)·Re·Pr
    Sh = (f/2)·Re·Sc
```
**Condition**: Turbulent flow, Pr ≈ Sc ≈ 1 (gases).

### 2. Chilton-Colburn Analogy (1933)
```
For 0.6 < Pr < 60, 0.6 < Sc < 3000, turbulent flow:
    j_H = St·Pr^(2/3) = f_F
    j_D = St_m·Sc^(2/3) = f_F
    ∴ j_H = j_D = f_F
```
**Engineering significance**: Measure friction → predict heat and mass transfer.

### 3. Blasius-Pohlhausen Laminar Analogy
```
Laminar flat plate:
    Cf_x = 0.664/√Re_x
    Nu_x = 0.332·√Re_x·Pr^(1/3)
    Sh_x = 0.332·√Re_x·Sc^(1/3)
    
    j_H = Nu_x/(Re_x·Pr)·Pr^(2/3) = 0.332/√Re_x
    Cf_x/2 = 0.332/√Re_x
    ∴ j_H = Cf_x/2 (exact for laminar flow!)
```

### 4. Boundary Layer Thickness Ratios
```
δ_T/δ ≈ Pr^(-1/3)    δ_C/δ ≈ Sc^(-1/3)    δ_T/δ_C ≈ Le^(1/3)

Pr > 1: δ_T < δ   (thermal BL thinner — water, oils)
Pr < 1: δ_T > δ   (thermal BL thicker — liquid metals)
Pr = 1: δ_T = δ   (Reynolds analogy exact — gases)
```

### 5. Pipe Flow Analogy
```
Laminar:    f_D = 64/Re,     Nu_D = 3.66,  Sh_D = 3.66
Turbulent:  f_D = 0.046·Re^(-1/5)
            Nu_D = 0.023·Re^(4/5)·Pr^n  (Dittus-Boelter)
            j_H = Nu/(Re·Pr^(1/3)) = f_D/8
```

## Core Algorithms

| Algorithm | Complexity | Description |
|-----------|-----------|-------------|
| `sutherland_viscosity()` | O(1) | Gas viscosity via Sutherland's law |
| `chapman_enskog_viscosity()` | O(1) | Kinetic theory viscosity (Lennard-Jones) |
| `eucken_thermal_conductivity()` | O(1) | Gas conductivity via Eucken formula |
| `fuller_diffusivity()` | O(1) | Gas diffusivity (Fuller-Schettler-Giddings) |
| `wilke_chang_diffusivity()` | O(1) | Liquid diffusivity (Wilke-Chang) |
| `blasius_profile()` | O(n) | RK4 solution of Blasius BL equation |
| `pohlhausen_thermal_profile()` | O(n) | Thermal BL profile (Pohlhausen) |
| `colebrook_friction()` | O(iter) | Newton's method for Colebrook-White |
| `predict_h_from_friction()` | O(1) | Analogy: ΔP → f → h |
| `compute_pipe_flow_analogy()` | O(1) | Complete pipe flow state |
| `wilke_mixture_viscosity()` | O(n²) | Gas mixture viscosity (Wilke) |
| `mason_saxena_mixture_conductivity()` | O(n²) | Gas mixture conductivity |

## Classic Problems

| Problem | Example | Key Result |
|---------|---------|-----------|
| Air over heated flat plate | `example_air_flat_plate.c` | Cf→Nu→h prediction from analogy |
| Water pipe flow with heating | `example_water_pipe_flow.c` | ΔP→h verified against Dittus-Boelter |
| Multi-fluid analogy comparison | `example_fluid_comparison.c` | Air/Water/Oil/Hg analogy accuracy |
| Electronics cooling | `electronics_cooling_analogy()` | Chip h from air velocity |
| Heat exchanger sizing | `heat_exchanger_sizing_by_analogy()` | Area from friction data |
| Reactor wall mass transfer | `reactor_wall_mass_transfer()` | k_c from ΔP |
| Cooling tower design | `cooling_tower_evaporation_rate()` | Evaporation from heat transfer |

## Nine-School Curriculum Mapping

| School | Course | Topics Covered |
|--------|--------|---------------|
| MIT | 2.25 Fluids, 2.005 Thermal-Fluids | BL theory, pipe flow, analogies |
| Stanford | ME 351A, ME 346A | Internal/external convection |
| Berkeley | ME 106, ME 105 | Convection, analogies |
| Michigan | ME 320, ME 420 | Pipe flow, boundary layers |
| Purdue | ME 505, ME 509 | Advanced convection, analogies |
| TU Munich | MW 0854 | Convection, heat transfer |
| ETH | 151-0111, 151-0103 | BL theory, analogies |
| Tsinghua | 传热学, 流体力学 | Convection, mass transfer |

## Project Structure

```
mini-momentum-heat-mass-analogy/
├── Makefile                          # make test
├── README.md                         # This file
├── include/
│   ├── transport_coefficients.h      # μ, k, D definitions and models
│   ├── momentum_heat_mass_analogy.h  # Core analogy relationships
│   ├── boundary_layer_analogy.h      # Velocity/thermal/concentration BL
│   ├── dimensionless_groups.h        # Re, Pr, Sc, Nu, Sh, St, etc.
│   ├── tube_channel_analogy.h        # Pipe/duct flow analogies
│   └── transport_analogy_database.h  # 12-fluid engineering database
├── src/
│   ├── transport_coefficients.c      # 454 lines — property computations
│   ├── momentum_heat_mass_analogy.c  # 345 lines — core analogy
│   ├── boundary_layer_analogy.c      # 444 lines — BL solutions
│   ├── dimensionless_groups.c        # 325 lines — dimensionless numbers
│   ├── tube_channel_analogy.c        # 495 lines — pipe flow
│   ├── transport_analogy_database.c  # 569 lines — fluid database
│   └── analogy_theorems.lean         # 716 lines — formal verification
├── tests/
│   ├── test_transport_coefficients.c # 21 tests
│   └── test_analogy.c               # 24 tests
├── examples/
│   ├── example_air_flat_plate.c      # >100 lines, end-to-end
│   ├── example_water_pipe_flow.c     # >100 lines, end-to-end
│   └── example_fluid_comparison.c    # >150 lines, end-to-end
├── docs/
│   ├── knowledge-graph.md            # L1-L9 complete coverage table
│   ├── coverage-report.md            # Assessment by level
│   ├── gap-report.md                 # Filled gaps + remaining
│   ├── course-alignment.md           # Nine-school curriculum mapping
│   └── course-tree.md                # Prerequisite dependency tree
├── demos/                            # Visualization/demos directory
└── benches/                          # Performance benchmarks directory
```

## Building and Testing

```bash
# Build everything
make all

# Run tests
make test

# Build examples
make examples

# Check Lean formalization
make lean

# Clean
make clean
```

## Quick Start

```c
#include "transport_coefficients.h"
#include "momentum_heat_mass_analogy.h"

// Get air properties at 300K
double mu = air_viscosity(300.0);        // 1.85e-5 Pa·s
double k  = air_thermal_conductivity(300.0); // 0.026 W/(m·K)
double Pr = mu * 1007.0 / k;             // 0.71

// Predict h from friction measurement
double h = predict_h_from_friction(
    delta_P, length, diameter, rho, velocity, cp, Pr);

printf("Predicted heat transfer coefficient: %.1f W/(m²·K)\n", h);
```

## Key Physical Insights

1. **The analogy works because**: All three transport mechanisms share the same mathematical form — flux ∝ gradient. The dimensionless governing equations become identical when Pr = Sc = 1.

2. **Turbulence enables the analogy**: In turbulent flow, eddy mixing transports momentum, heat, and mass with nearly equal efficiency. This is why the analogy works best for turbulent flows.

3. **Pr and Sc determine BL thickness**: Pr > 1 means the thermal BL is thinner than the velocity BL (oils). Pr < 1 means thicker (liquid metals). This affects where in the flow heat transfer resistance lies.

4. **Engineering power**: The most important practical application is predicting heat transfer coefficients from pressure drop measurements — saving the cost and complexity of thermal instrumentation.

5. **Limits**: The analogy fails for laminar flow (no eddy mixing), extreme Pr/Sc (different transport mechanisms dominate), and when property variations are large (Sieder-Tate correction needed).

## References

- Bird, R.B., Stewart, W.E., Lightfoot, E.N. (2007). *Transport Phenomena*, 2nd Ed. Wiley.
- Colburn, A.P. (1933). "A method of correlating forced convection heat-transfer data and a comparison with fluid friction." *Trans. AIChE*, 29, 174-210.
- Incropera, F.P., DeWitt, D.P. (2007). *Fundamentals of Heat and Mass Transfer*, 6th Ed. Wiley.
- Schlichting, H., Gersten, K. (2000). *Boundary Layer Theory*, 8th Ed. Springer.
- White, F.M. (2016). *Fluid Mechanics*, 8th Ed. McGraw-Hill.

## License

MIT

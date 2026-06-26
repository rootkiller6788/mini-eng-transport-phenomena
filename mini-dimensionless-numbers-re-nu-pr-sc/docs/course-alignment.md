# Course Alignment — mini-dimensionless-numbers-re-nu-pr-sc

## Nine-School Curriculum Mapping

### MIT — 2.005 Thermal-Fluids Engineering / 2.25 Advanced Fluid Mechanics

| Chapter | Topic | This Module |
|---------|-------|-------------|
| Dimensional Analysis | Buckingham Pi, similarity | `buckingham_pi.h/c` — full Pi theorem implementation |
| Pipe Flow | Re, friction factor, Moody diagram | `reynolds_number.h/c` — 7 friction models, pipe design |
| External Flow | Boundary layer, drag coefficient | `nusselt_number.h/c` — flat plate, cylinder, sphere |
| Convection | Nu correlations, forced/natural | 18+ Nu correlations across all geometries |
| Heat Exchangers | LMTD method, ε-NTU | `heat_exchanger_sizing()` |

### Stanford — ME 346A Heat Transfer / ME 351A Fluid Mechanics

| Chapter | Topic | This Module |
|---------|-------|-------------|
| Conduction | Bi, Fo, lumped capacitance | `compute_biot_number()`, `compute_fourier_number()` |
| Forced Convection | Flat plate, pipe, cylinder | All Nu correlations |
| Natural Convection | Ra, Churchill-Chu | `nu_natural_vertical_plate_churchill_chu()` |
| Mass Transfer | Sc, Sh, Lewis number | `prandtl_schmidt.h/c` — full analogy |
| Silicon Valley Apps | Electronics cooling | `electronics_cooling_finned_array()` |

### Berkeley — ME 105 Thermodynamics / ME 106 Fluid Mechanics / ME 107 Heat Transfer

| Chapter | Topic | This Module |
|---------|-------|-------------|
| Fluid Properties | ρ, μ, ν, α, D, β, σ | `FluidProperties` struct |
| Dimensional Analysis | Pi theorem, non-dimensionalization | Full matrix method |
| Convection Correlations | Dittus-Boelter, Sieder-Tate | Both implemented |
| Heat Exchangers | Sizing, LMTD | `heat_exchanger_sizing()` |

### Michigan — ME 320 Fluid Mechanics / ME 420 Heat Transfer / AERO 533 Combustion

| Chapter | Topic | This Module |
|---------|-------|-------------|
| Turbulence | Kolmogorov scale, Reynolds stress | `scaling_law_kolmogorov_microscale()` |
| Mass Transfer | Droplet evaporation, combustion | `mass_transfer_droplet_evaporation()` |
| Automotive Apps | Fuel line, engine cooling | `example_pipe_flow_regime.c` (Detroit), `dc_motor_natural_cooling()` |

### Purdue — ME 505 Heat Transfer / ME 509 Fluid Mechanics / ME 597 Plasma

| Chapter | Topic | This Module |
|---------|-------|-------------|
| Advanced Convection | Mixed convection, superposition | `mixed_convection_nusselt()` |
| Natural Convection | Enclosure, Ra-critical | `nu_natural_enclosure()`, Rayleigh-Bénard onset |
| Correlations | Gnielinski, Churchill-Chu | Both implemented as most accurate options |

### TU Munich — MW 0798 Thermodynamics / MW 0854 Heat Transfer

| Chapter | Topic | This Module |
|---------|-------|-------------|
| Ähnlichkeitstheorie | Similarity theory | `check_dynamic_similarity()` |
| Wärmeübergang | Heat transfer correlations | Full Nu correlation suite |
| Stoffübergang | Mass transfer | Sh, Sc, Le — complete |

### ETH — 151-0103 Fluid Dynamics / 151-0111 Heat Transfer

| Chapter | Topic | This Module |
|---------|-------|-------------|
| Scaling Laws | Kolmogorov, Batchelor | Both microscales implemented |
| Boundary Layers | Blasius, displacement thickness | Full BL suite |
| Convection | Churchill-Chu, Gnielinski | Best-available correlations |

### Tsinghua — 工程热力学 / 传热学 / 流体力学

| Chapter | Topic | This Module |
|---------|-------|-------------|
| 量纲分析 | Dimensional analysis | Pi theorem with Chinese-relevant examples |
| 管内流动 | Pipe flow, 摩擦系数 | 7 friction factor methods |
| 对流传热 | Convective heat transfer | Full Nu correlation suite |
| 自然对流 | Natural convection | Ra, Gr, Churchill-Chu |
| 换热器 | Heat exchanger design | LMTD method |

## Key Textbook Alignment

| Textbook | Coverage in This Module |
|----------|------------------------|
| Bird, Stewart, Lightfoot (2007) "Transport Phenomena" | Full: Re, Pr, Sc, Nu, Sh, Le, Chilton-Colburn, BL theory, analogies |
| Incropera & DeWitt (2007) "Fundamentals of Heat and Mass Transfer" | Full: all Nu correlations, natural convection, HX sizing |
| White (2016) "Fluid Mechanics" | Full: Re, friction factors, Moody, BL, dimensional analysis |
| Buckingham (1914) "On Physically Similar Systems" | Full: Pi theorem matrix method |
| Bejan (2013) "Convection Heat Transfer" | Full: scaling laws, mixed convection |

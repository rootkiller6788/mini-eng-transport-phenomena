# Knowledge Graph — mini-dimensionless-numbers-re-nu-pr-sc

## L1: Definitions

| # | Concept | Definition | Code |
|---|---------|------------|------|
| 1 | Reynolds Number | Re = ρUL/μ, ratio of inertial to viscous forces | `DimensionlessNumbers.Re`, `compute_reynolds_number()` |
| 2 | Nusselt Number | Nu = hL/k, convective/conductive heat transfer | `DimensionlessNumbers.Nu`, `compute_nusselt_number()` |
| 3 | Prandtl Number | Pr = ν/α = cp·μ/k | `DimensionlessNumbers.Pr`, `compute_prandtl_number()` |
| 4 | Schmidt Number | Sc = ν/D, momentum/mass diffusivity | `DimensionlessNumbers.Sc`, `compute_schmidt_number()` |
| 5 | Péclet Number | Pe = Re·Pr (heat) or Re·Sc (mass) | `DimensionlessNumbers.Pe_H`, `compute_peclet_number_heat()` |
| 6 | Grashof Number | Gr = gβΔT L³/ν², buoyancy/viscous | `DimensionlessNumbers.Gr`, `compute_grashof_number()` |
| 7 | Rayleigh Number | Ra = Gr·Pr, driving force for natural convection | `DimensionlessNumbers.Ra`, `compute_rayleigh_number()` |
| 8 | Stanton Number | St = Nu/(Re·Pr) = h/(ρcpU) | `DimensionlessNumbers.St_H`, `compute_stanton_number_heat()` |
| 9 | Weber Number | We = ρU²L/σ, inertia/surface tension | `DimensionlessNumbers.We`, `compute_weber_number()` |
| 10 | Mach Number | Ma = U/c | `DimensionlessNumbers.Ma`, `compute_mach_number()` |
| 11 | Froude Number | Fr = U/√(gL), inertia/gravity | `DimensionlessNumbers.Fr`, `compute_froude_number()` |
| 12 | Euler Number | Eu = ΔP/(ρU²) | `DimensionlessNumbers.Eu`, `compute_euler_number()` |
| 13 | Brinkman Number | Br = μU²/(kΔT), viscous heating/heat conduction | `DimensionlessNumbers.Br`, `compute_brinkman_number()` |
| 14 | Biot Number | Bi = hLc/k_s, internal conduction/surface convection | `DimensionlessNumbers.Bi`, `compute_biot_number()` |
| 15 | Fourier Number | Fo = αt/L², dimensionless time | `DimensionlessNumbers.Fo`, `compute_fourier_number()` |
| 16 | Sherwood Number | Sh = h_m L/D, convective/diffusive mass transfer | `DimensionlessNumbers.Sh`, `compute_sherwood_number()` |
| 17 | Lewis Number | Le = α/D = Sc/Pr | `DimensionlessNumbers.Le`, `compute_lewis_number()` |
| 18 | Graetz Number | Gz = Re·Pr·D/L, thermal entrance region | `DimensionlessNumbers.Gz`, `compute_graetz_number()` |
| 19 | Capillary Number | Ca = μU/σ | `DimensionlessNumbers.Ca`, `compute_capillary_number()` |
| 20 | Bond Number | Bo = ρgL²/σ (= Eötvös) | `DimensionlessNumbers.Bo`, `compute_bond_number()` |

## L2: Core Concepts

| # | Concept | Implementation |
|---|---------|---------------|
| 1 | Dynamic Similarity (Reynolds similarity) | `check_dynamic_similarity()` |
| 2 | Thermal Similarity | `check_thermal_similarity()` |
| 3 | Buckingham Pi Theorem | `pi_theorem_matrix_method()`, `pi_rank_dimensional_matrix()` |
| 4 | Dimensional Homogeneity | `verify_dimensionless_group()` |
| 5 | Flow Regime Classification | `classify_flow_regime()`, `flow_regime_pipe()` |
| 6 | Convection Type (forced/natural/mixed) | `classify_convection_type()`, `convection_regime_map()` |
| 7 | Boundary Layer Theory (δ, δ*, θ) | `bl_thickness_laminar()`, `bl_displacement_thickness_laminar()`, `bl_momentum_thickness_laminar()` |
| 8 | Transport Analogy (Reynolds, Chilton-Colburn) | `colburn_analogy_verify()`, `reynolds_analogy_stanton()` |
| 9 | Heat-Mass Transfer Analogy | `heat_mass_analogy_nusselt_to_sherwood()` |
| 10 | Non-dimensional N-S Equation | `ns_re_scale()`, `energy_eq_dimensionless_form()` |

## L3: Engineering Quantities

| # | Quantity | Implementation |
|---|----------|---------------|
| 1 | Critical Re Database (10 geometries) | `critical_reynolds()`, `critical_re_description()` |
| 2 | Pr Category Classification | `prandtl_category()` |
| 3 | Typical Pr ranges (air, water, oil) | `prandtl_air()`, `prandtl_water()`, `prandtl_engine_oil()` |
| 4 | Typical Sc ranges (gas, liquid) | `schmidt_gas_typical()`, `schmidt_liquid_typical()` |
| 5 | Typical Le values | `lewis_typical()` |
| 6 | Gas Diffusivity (Chapman-Enskog) | `mass_diffusivity_gas_estimate()` |
| 7 | Liquid Diffusivity (Wilke-Chang) | `mass_diffusivity_liquid_wilke_chang()` |
| 8 | Dimensional Form Database (23 variables) | `get_standard_dimension()` |

## L4: Conservation Laws

| # | Law | Implementation |
|---|-----|---------------|
| 1 | Navier-Stokes Non-dimensionalization → Re | `ns_re_scale()`, `ns_pressure_scale()`, `ns_time_scale()` |
| 2 | Energy Equation → Pe, Br | `energy_eq_dimensionless_form()` |
| 3 | Species Equation → Sc, Sh | `schmidt_basic()`, `compute_sherwood_number()` |
| 4 | Dimensional Homogeneity (Pi Theorem) | `verify_dimensionless_group()`, `verify_all_standard_groups()` |
| 5 | Lean Formalized: Re, Nu, Pr, Sc, Flow Regime, Colburn | `analogy_theorems.lean` |

## L5: Engineering Methods

| # | Method | Implementation |
|---|--------|---------------|
| 1 | Pi Theorem Matrix Method | `pi_theorem_matrix_method()`, `pi_nullspace()` |
| 2 | Laminar Friction Factor (exact, 64/Re) | `friction_factor_laminar()` |
| 3 | Blasius Friction Factor (turbulent) | `friction_factor_blasius()` |
| 4 | Colebrook Equation (implicit, iterative) | `friction_factor_colebrook()` |
| 5 | Haaland (explicit approximation) | `friction_factor_haaland()` |
| 6 | Swamee-Jain (explicit approximation) | `friction_factor_swamee_jain()` |
| 7 | Churchill (all-regime single equation) | `friction_factor_churchill()` |
| 8 | Nu Correlations — Internal Flow (7 types) | Pipe: Dittus-Boelter, Sieder-Tate, Gnielinski, Laminar; Annulus; Non-circular ducts |
| 9 | Nu Correlations — External Flow (6 types) | Flat plate (laminar/turbulent/mixed), Cylinder, Sphere, Tube bank |
| 10 | Nu Correlations — Natural Convection (5 types) | Vertical plate, Horizontal plate, Cylinder, Sphere, Enclosure |
| 11 | Mass Transfer Correlations (3 types) | Flat plate, Pipe, Sphere/droplet |
| 12 | Colburn/Chilton-Colburn Analogy | `colburn_analogy_verify()`, `chilton_colburn_sherwood()` |
| 13 | Boundary Layer Scaling Laws | `thermal_to_velocity_bl_ratio()`, `scaling_law_kolmogorov_microscale()` |
| 14 | Model Scaling (Re/Fr/St matching) | `model_scale_velocity_reynolds()`, `model_scale_force()`, `model_scale_frequency()` |

## L6: Engineering Problems

| # | Problem | Implementation |
|---|---------|---------------|
| 1 | Pipe Flow Design (ΔP, pump power) | `pipe_flow_design()` |
| 2 | Flat Plate Heat Transfer (solar panel, wing de-icing) | `example_flat_plate_heat_transfer.c` |
| 3 | Natural Convection Cooling (radiator, motor, window) | `example_natural_convection.c` |
| 4 | Mass Transfer (CO₂ absorption, droplet evaporation) | `example_mass_transfer.c` |
| 5 | Heat Exchanger Sizing (LMTD method) | `heat_exchanger_sizing()` |
| 6 | Pipe Flow Regime Determination | `example_pipe_flow_regime.c` |
| 7 | Electronics Cooling (finned array) | `electronics_cooling_finned_array()` |

## L7: Applications

| # | Application | Implementation |
|---|-------------|---------------|
| 1 | DC Motor Natural Cooling (Toyota Prius, Tesla) | `dc_motor_natural_cooling()` |
| 2 | DC Motor Forced Cooling (servo drives) | `dc_motor_forced_cooling()` |
| 3 | Air-Cooled Condenser (Boeing 747 ECS power) | `air_cooled_condenser()` |
| 4 | Mars Rover Night Survival (NASA Perseverance) | `mars_rover_night_survival()` |
| 5 | iPhone Thermal Management (Apple) | `example_flat_plate_heat_transfer.c` |
| 6 | Model Scaling — Wind Tunnel Design | `model_scale_velocity_reynolds()`, `model_scale_force()` |

## L8: Advanced Methods

| # | Method | Implementation |
|---|--------|---------------|
| 1 | Kolmogorov Microscale (turbulence) | `scaling_law_kolmogorov_microscale()` |
| 2 | Batchelor Microscale (scalar turbulence) | `scaling_law_batchelor_microscale()` |
| 3 | Buoyant Plume Scaling (Morton-Taylor-Turner) | `scaling_law_buoyant_plume_velocity()` |
| 4 | Diffusion Length Scale | `scaling_law_diffusion_length()` |
| 5 | Mixed Convection Correlation (superposition) | `mixed_convection_nusselt()` |

## L9: Research Frontiers

| # | Topic | Status |
|---|-------|--------|
| 1 | Nano-scale Transport (Knudsen number effects) | Documented, not implemented |
| 2 | Non-Fourier Heat Conduction (Cattaneo-Vernotte) | Documented |
| 3 | Micro-channel Flow (slip boundary conditions) | Documented |
| 4 | Turbulent Prandtl Number (Reynolds stress closure) | Documented |
| 5 | Marangoni Convection (thermocapillary, Bénard-Marangoni) | Documented |

# Knowledge Graph — Mini Momentum-Heat-Mass Analogy

## L1: Definitions (Complete ✅)

| # | Definition | C Implementation | Lean Formalization |
|---|-----------|-----------------|-------------------|
| 1 | Dynamic viscosity μ [Pa·s] | `transport_coefficients.h::ViscosityDefinition` | `TransportTriplet` |
| 2 | Kinematic viscosity ν = μ/ρ [m²/s] | `momentum_diffusivity()` | — |
| 3 | Thermal conductivity k [W/(m·K)] | `transport_coefficients.h::ThermalConductivityDef` | `TransportTriplet` |
| 4 | Thermal diffusivity α = k/(ρ·cp) [m²/s] | `thermal_diffusivity()` | — |
| 5 | Mass diffusivity D_AB [m²/s] | `transport_coefficients.h::MassDiffusivityDef` | `TransportTriplet` |
| 6 | Reynolds number Re = ρ·v·L/μ | `compute_Re()` | `DimensionlessGroups` |
| 7 | Prandtl number Pr = ν/α = μ·Cp/k | `compute_Pr()` | `DimensionlessGroups` |
| 8 | Schmidt number Sc = ν/D | `compute_Sc()` | `DimensionlessGroups` |
| 9 | Lewis number Le = α/D = Sc/Pr | `compute_Le()` | `DimensionlessGroups` |
| 10 | Nusselt number Nu = h·L/k | `compute_Nu()` | `StantonHeat` |
| 11 | Sherwood number Sh = k_c·L/D | `compute_Sh()` | `StantonMass` |
| 12 | Stanton number St = Nu/(Re·Pr) | `compute_St_heat()` | `StantonHeat` |
| 13 | Skin friction coefficient Cf = 2·τ_w/(ρ·U∞²) | `laminar_cf_local()` | `ReynoldsAnalogyState` |
| 14 | Fanning friction factor f_F | `fanning_friction_laminar()` | `PipeFlowState` |
| 15 | Darcy friction factor f_D = 4·f_F | `darcy_friction_laminar()` | `PipeFlowState` |
| 16 | Colburn j-factor j_H = St·Pr^(2/3) | `compute_colburn_j_factors()` | `ColburnState` |
| 17 | Boundary layer thickness δ (99%) | `laminar_bl_thickness()` | `BoundaryLayerTriplet` |
| 18 | Grashof number Gr | `compute_Gr()` | — |
| 19 | Rayleigh number Ra = Gr·Pr | `compute_Ra()` | — |
| 20 | Péclet number Pe = Re·Pr | `compute_Pe()` | `peclet_from_re_pr` |

## L2: Core Concepts (Complete ✅)

| # | Concept | Implementation |
|---|---------|---------------|
| 1 | Newton's law of viscosity: τ = -μ·dv/dy | `newtons_law_shear_stress()` |
| 2 | Fourier's law: q = -k·dT/dy | `fouriers_law_heat_flux()` |
| 3 | Fick's first law: J = -D·dC/dy | `ficks_law_mass_flux()` |
| 4 | Reynolds analogy: St = Cf/2 (Pr=1) | `stanton_number_from_friction()` |
| 5 | Chilton-Colburn analogy: j_H = j_D = f/2 | `compute_colburn_j_factors()` |
| 6 | Prandtl-Taylor two-layer model | `prandtl_taylor_analogy_St()` |
| 7 | von Karman three-layer model | `von_karman_analogy_St()` |
| 8 | Blasius flat-plate BL solution | `blasius_profile()` |
| 9 | Pohlhausen thermal BL solution | `pohlhausen_thermal_profile()` |
| 10 | Hagen-Poiseuille laminar pipe flow | `darcy_friction_laminar()` |
| 11 | Dittus-Boelter turbulent correlation | `dittus_boelter_nu()` |
| 12 | Colburn j-factor analogy for pipe flow | `colburn_j_factor_pipe()` |
| 13 | Boundary layer relative thicknesses | `boundary_layer_ordering()` |
| 14 | Flow regime classification | `flow_regime_pipe()`, `flow_regime_flat_plate()` |
| 15 | Transport equation unity | `generalized_transport_flux()` |
| 16 | Eddy diffusivity concept | Analogy regime description |
| 17 | Entry length effects | `developing_nusselt()`, `DevelopingFlowState` |
| 18 | Non-circular duct corrections | `setup_duct_geometry()` |

## L3: Engineering Quantities (Complete ✅)

| # | Quantity | Implementation |
|---|----------|---------------|
| 1 | Air viscosity μ(300K) ≈ 1.8×10⁻⁵ Pa·s | `air_viscosity()` |
| 2 | Water viscosity μ(293K) ≈ 1.0×10⁻³ Pa·s | `water_viscosity()` |
| 3 | Air k(300K) ≈ 0.026 W/(m·K) | `air_thermal_conductivity()` |
| 4 | Water k(300K) ≈ 0.61 W/(m·K) | `water_thermal_conductivity()` |
| 5 | D_H2O-air(298K) ≈ 2.6×10⁻⁵ m²/s | `water_vapor_air_diffusivity()` |
| 6 | Pr_air ≈ 0.71, Pr_water ≈ 7, Pr_oil ≈ 100-10000 | `typical_prandtl()` |
| 7 | Sc_gases ≈ 0.2-5, Sc_liquids ≈ 100-10000 | `typical_schmidt()` |
| 8 | 12-fluid built-in database with property models | `init_fluid_database()` |
| 9 | Re for common engineering scales | `Re_1ms_1mm` in `AnalogyLookup` |
| 10 | Critical Re: pipe=2300, flat plate=5×10⁵ | Flow regime functions |

## L4: Conservation Laws (Complete ✅)

| # | Law | Implementation |
|---|-----|---------------|
| 1 | Newton's law of viscosity | `newtons_law_shear_stress()` |
| 2 | Fourier's law of heat conduction | `fouriers_law_heat_flux()` |
| 3 | Fick's first law of diffusion | `ficks_law_mass_flux()` |
| 4 | Generalized transport flux | `generalized_transport_flux()` |
| 5 | Hagen-Poiseuille equation (laminar pipe) | `darcy_friction_laminar()` |
| 6 | Blasius BL equation (f''' + ½·f·f'' = 0) | `rk4_step_blasius()` |
| 7 | Pohlhausen energy equation | `pohlhausen_thermal_profile()` |
| 8 | Transport equation balance | `transport_equation_residual()` |
| 9 | Reynolds analogy identity (Lean) | `reynolds_analogy_stanton_equality` |
| 10 | Colburn symmetry (Lean) | `colburn_j_factor_equality` |
| 11 | Unified diffusivity theorem (Lean) | `unified_diffusivity_when_equal` |
| 12 | Richardson's annular effect | `pipe_flow_positivity` |
| 13 | Transport coefficient positivity (Lean) | `transport_coefficients_positive` |

## L5: Engineering Methods (Complete ✅)

| # | Method | Implementation |
|---|--------|---------------|
| 1 | Predict h from ΔP (friction → heat transfer) | `predict_h_from_friction()` |
| 2 | Predict k_c from h (heat → mass transfer) | `predict_kc_from_h()` |
| 3 | Predict Nu from friction factor | `predict_nu_from_friction()` |
| 4 | Predict Sh from friction factor | `predict_sh_from_friction()` |
| 5 | Predict h from pressure drop (pipe flow) | `predict_h_from_pressure_drop()` |
| 6 | Predict k_c from h (Chilton-Colburn) | `predict_kc_from_heat_transfer_coefficient()` |
| 7 | Heat exchanger sizing by analogy | `heat_exchanger_sizing_by_analogy()` |
| 8 | Cooling tower evaporation rate (Merkel) | `cooling_tower_evaporation_rate()` |
| 9 | Reactor wall mass transfer from ΔP | `reactor_wall_mass_transfer()` |
| 10 | Sieder-Tate viscosity correction | `sieder_tate_laminar_nu()` |
| 11 | Colebrook-White iterative solution | `colebrook_friction()` |
| 12 | Gnielinski correlation | `gnielinski_nu()` |
| 13 | Petukhov correlation | `petukhov_nu()` |
| 14 | Hausen developing flow correlation | `hausen_laminar_nu()` |
| 15 | Lyon-Martinelli liquid metal correlation | `lyon_martinelli_nu()` |

## L6: Engineering Problems (Complete ✅)

| # | Problem | Implementation |
|---|---------|---------------|
| 1 | Air flat plate: h prediction from Cf | `example_air_flat_plate.c` |
| 2 | Water pipe flow: h from ΔP | `example_water_pipe_flow.c` |
| 3 | Multi-fluid analogy comparison | `example_fluid_comparison.c` |
| 4 | Electronics cooling: chip h from v_air | `electronics_cooling_analogy()` |
| 5 | Heat exchanger sizing | `heat_exchanger_sizing_by_analogy()` |
| 6 | Non-circular duct analysis | `setup_duct_geometry()` |
| 7 | Analogy vs. experiment comparison | `print_analogy_vs_experiment_comparison()` |
| 8 | Developing flow analysis | `compute_developing_flow()` |
| 9 | Boundary layer comparison (air) | `air_flat_plate_benchmark()` |
| 10 | Pipe flow benchmark (water) | `water_pipe_flow_benchmark()` |

## L7: Engineering Applications (Partial+ ✅ — 5 applications)

| # | Application | Implementation |
|---|------------|---------------|
| 1 | Heat exchanger sizing by analogy | `heat_exchanger_sizing_by_analogy()` |
| 2 | Cooling tower evaporation rate | `cooling_tower_evaporation_rate()` |
| 3 | Chemical reactor wall mass transfer | `reactor_wall_mass_transfer()` |
| 4 | Electronics cooling (PCB/chip) | `electronics_cooling_analogy()` |
| 5 | Quick engineering estimation | `quick_analogy_estimate()` |

## L8: Advanced Methods (Partial+ ✅ — 3 topics)

| # | Topic | Implementation |
|---|-------|---------------|
| 1 | Prandtl-Taylor multi-layer model | `prandtl_taylor_analogy_St()` |
| 2 | von Karman three-layer analogy | `von_karman_analogy_St()` |
| 3 | Lyon-Martinelli liquid metal model | `lyon_martinelli_nu()` |
| 4 | Colebrook-White iterative solver | `colebrook_friction()` |
| 5 | Pohlhausen thermal BL numerical solution | `pohlhausen_thermal_profile()` |

## L9: Research Frontiers (Partial ✅ — documented)

| # | Topic | Documentation |
|---|-------|--------------|
| 1 | Nano-scale transport (non-Fourier effects) | Documented in course-tree.md |
| 2 | Turbulent Prandtl number models | Documented |
| 3 | Multi-component diffusion analogies | Wilke/Mason-Saxena mixture rules |
| 4 | Non-Newtonian momentum-heat analogies | Mentioned in ViscosityDefinition |
| 5 | Machine learning for transport correlations | Future work |

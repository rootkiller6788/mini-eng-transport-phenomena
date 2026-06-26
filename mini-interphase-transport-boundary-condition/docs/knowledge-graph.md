# Knowledge Graph — Interphase Transport & Boundary Conditions

## L1: Definitions
| ID | Concept | C Implementation | Status |
|----|---------|-----------------|--------|
| L1.1 | Phase types (gas, liquid, solid, plasma, supercritical) | `PhaseType` enum | Complete |
| L1.2 | Boundary condition taxonomy (Dirichlet, Neumann, Robin, Jump, Periodic, Symmetry, Outflow, Cauchy) | `BCType` enum | Complete |
| L1.3 | Interface geometry (planar, cylindrical, spherical, ellipsoidal, arbitrary, corrugated, fractal) | `InterfaceGeometry` enum | Complete |
| L1.4 | Phase pair classification (gas-liquid, liquid-liquid, gas-solid, liquid-solid, liquid-vapor, solid-solid, gas-plasma, immiscible liquid) | `PhasePair` enum | Complete |
| L1.5 | Thermodynamic state at interface (T, P, rho, mu, k, cp, c) | `PhaseState` struct | Complete |
| L1.6 | Multicomponent mixture state | `MixtureState` struct | Complete |
| L1.7 | Interface descriptor (pair, geometry, states, sigma, curvature, roughness, slip) | `InterfaceDescriptor` struct | Complete |
| L1.8 | Stress tensor | `StressTensor` struct | Complete |
| L1.9 | Heat flux vector | `HeatFlux` struct | Complete |
| L1.10 | Mass flux vector | `MassFlux` struct | Complete |
| L1.11 | Combined interphase flux | `InterphaseFlux` struct | Complete |
| L1.12 | Velocity BC specification | `VelocityBC` struct | Complete |
| L1.13 | Thermal BC specification | `ThermalBC` struct | Complete |
| L1.14 | Concentration BC specification | `ConcentrationBC` struct | Complete |
| L1.15 | Combined BC set | `BoundaryConditionSet` struct | Complete |
| L1.16 | Henry's law coefficients | `HenryLaw` struct | Complete |
| L1.17 | Vapor-liquid equilibrium data | `VaporLiquidEquilibrium` struct | Complete |

## L2: Core Concepts
| ID | Concept | C Implementation | Status |
|----|---------|-----------------|--------|
| L2.1 | Two-film theory (Whitman, 1923) | `two_film_mass_transfer()` | Complete |
| L2.2 | Penetration theory (Higbie, 1935) | `penetration_mass_transfer()` | Complete |
| L2.3 | Surface renewal theory (Danckwerts, 1951) | `surface_renewal_mass_transfer()` | Complete |
| L2.4 | Interfacial resistance models | `TwoFilmModel`, `PenetrationModel`, `SurfaceRenewalModel` | Complete |
| L2.5 | Henry's law equilibrium at gas-liquid interface | `interface_henry_equilibrium()` | Complete |
| L2.6 | Raoult's law for vapor-liquid interface | `interface_raoult_equilibrium()` | Complete |
| L2.7 | Partition coefficient at liquid-liquid interface | `interface_partition_equilibrium()` | Complete |
| L2.8 | Velocity continuity at fluid-fluid interface | `interface_velocity_continuity()` | Complete |
| L2.9 | Shear stress continuity | `interface_shear_stress_continuity()` | Complete |
| L2.10 | Temperature continuity | `interface_temperature_continuity()` | Complete |
| L2.11 | Heat flux continuity with phase change | `interface_heat_flux_continuity()` | Complete |
| L2.12 | Chemical reaction enhancement of mass transfer | `enhancement_factor_instantaneous()`, `hatta_number()` | Complete |
| L2.13 | Hatta number for reaction-diffusion | `hatta_number()` | Complete |
| L2.14 | No-slip condition (Stokes, 1851) | `velocity_bc_no_slip()` | Complete |
| L2.15 | Navier slip condition | `velocity_bc_navier_slip()` | Complete |
| L2.16 | Marangoni stress | `interface_marangoni_stress()` | Complete |
| L2.17 | Young-Laplace pressure jump | `interface_young_laplace()` | Complete |

## L3: Engineering Quantities
| ID | Quantity | C Implementation | Status |
|----|----------|-----------------|--------|
| L3.1 | Dimensionless groups collection (Re, Pr, Sc, Nu, Sh, Pe, Gr, Ra, St, Bi, Fo, We, Ca, Bo, Ma, Ma, Kn, Le) | `DimensionlessGroups` struct, `compute_dimensionless_groups()` | Complete |
| L3.2 | Prandtl number | `compute_prandtl()` | Complete |
| L3.3 | Schmidt number | `compute_schmidt()` | Complete |
| L3.4 | Lewis number | `compute_lewis()` | Complete |
| L3.5 | Capillary length | `interface_capillary_length()` | Complete |
| L3.6 | Chapman-Enskog diffusivity | `chapman_enskog_diffusivity()` | Complete |
| L3.7 | Chung viscosity estimation | `chung_viscosity()` | Complete |
| L3.8 | Eucken thermal conductivity | `eucken_thermal_conductivity()` | Complete |
| L3.9 | Interfacial area (packed column) | `interfacial_area_packed_column()` | Complete |
| L3.10 | System property data (engineering handbooks) | `SystemProperties` struct | Complete |

## L4: Conservation Laws
| ID | Law | C Implementation | Status |
|----|-----|-----------------|--------|
| L4.1 | Mass jump condition at interface | `jump_mass()` | Complete |
| L4.2 | Linear momentum jump condition | `jump_momentum()` | Complete |
| L4.3 | Energy jump condition | `jump_energy()` | Complete |
| L4.4 | Species jump condition | `jump_species()` | Complete |
| L4.5 | Entropy jump condition (2nd Law) | `jump_entropy()` | Complete |
| L4.6 | Complete jump condition evaluation | `jump_all_conditions()` | Complete |
| L4.7 | Stefan condition (phase change) | `stefan_condition()` | Complete |
| L4.8 | Hertz-Knudsen-Schrage (evaporation kinetics) | `hertz_knudsen_schrage()` | Complete |
| L4.9 | Interface mass continuity | `interface_mass_continuity()` | Complete |
| L4.10 | Integral balance control volume | `InterfaceControlVolume` struct | Complete |

## L5: Engineering Methods
| ID | Method | C Implementation | Status |
|----|--------|-----------------|--------|
| L5.1 | Two-film overall mass transfer coefficient | `two_film_mass_transfer()` | Complete |
| L5.2 | Penetration theory k_L computation | `penetration_mass_transfer()` | Complete |
| L5.3 | Surface renewal k_L computation | `surface_renewal_mass_transfer()` | Complete |
| L5.4 | Mass transfer theory comparison | `compare_mass_transfer_theories()` | Complete |
| L5.5 | Nusselt number - laminar flat plate | `nusselt_flat_plate_laminar()` | Complete |
| L5.6 | Nusselt number - turbulent flat plate | `nusselt_flat_plate_turbulent()` | Complete |
| L5.7 | Nusselt number - laminar pipe | `nusselt_pipe_laminar()` | Complete |
| L5.8 | Nusselt number - turbulent pipe (Dittus-Boelter, Gnielinski) | `nusselt_pipe_turbulent()` | Complete |
| L5.9 | Nusselt number - natural convection (Churchill-Chu) | `nusselt_natural_convection_vertical()` | Complete |
| L5.10 | Nusselt number - film condensation | `nusselt_film_condensation()` | Complete |
| L5.11 | Rohsenow pool boiling correlation | `rohsenow_pool_boiling()` | Complete |
| L5.12 | Sherwood - falling film | `sherwood_falling_film()` | Complete |
| L5.13 | Sherwood - single sphere (Froessling) | `sherwood_single_sphere()` | Complete |
| L5.14 | Sherwood - packed bed (Wakao-Funazkri) | `sherwood_packed_bed()` | Complete |
| L5.15 | Sherwood - bubble column (Akita-Yoshida) | `sherwood_bubble_column()` | Complete |
| L5.16 | Level-set signed distance function | `level_set_signed_distance()` | Complete |
| L5.17 | Regularized Heaviside and delta functions | `regularized_heaviside()`, `regularized_delta()` | Complete |
| L5.18 | Level-set interface normal computation | `level_set_normal()` | Complete |
| L5.19 | Level-set curvature (2D) | `level_set_curvature_2d()` | Complete |
| L5.20 | PLIC VOF fraction computation | `plic_vof_fraction_2d()` | Complete |
| L5.21 | Inverse PLIC reconstruction | `plic_inverse_vof_2d()` | Complete |
| L5.22 | VOF flux computation | `vof_flux_2d()` | Complete |
| L5.23 | PLIC interface area | `plic_interface_area_2d()` | Complete |
| L5.24 | Level-set advection (upwind) | `level_set_advect_upwind()` | Complete |
| L5.25 | Level-set reinitialization (Sussman) | `level_set_reinitialize()` | Complete |
| L5.26 | CSF surface tension force (Brackbill) | `csf_surface_tension_2d()` | Complete |
| L5.27 | Murphree tray efficiency | `murphree_tray_efficiency()` | Complete |
| L5.28 | Onda interfacial area correlation | `interfacial_area_packed_column()` | Complete |

## L6: Engineering Problems
| ID | Problem | Implementation | Status |
|----|---------|---------------|--------|
| L6.1 | Falling film gas absorption | `falling_film_absorption()` + example | Complete |
| L6.2 | Evaporative cooling / wet-bulb temperature | `evaporative_cooling()` | Complete |
| L6.3 | Distillation tray efficiency | `murphree_tray_efficiency()` | Complete |
| L6.4 | Membrane contactor sizing | `membrane_contactor_sizing()` | Complete |
| L6.5 | Condensation on horizontal tubes | `condensation_horizontal_tube()` + example | Complete |
| L6.6 | Droplet evaporation (d^2-law) | example_droplet_evaporation.c | Complete |
| L6.7 | Bubble column interfacial area | `sherwood_bubble_column()` | Complete |
| L6.8 | Power plant condenser design | example_condenser.c | Complete |

## L7: Applications
| ID | Application | Implementation | Status |
|----|------------|---------------|--------|
| L7.1 | Wastewater aeration design | `wastewater_aeration_design()` | Complete |
| L7.2 | PWR steam generator sizing | `pwr_steam_generator_sizing()` | Complete |
| L7.3 | PEM fuel cell water management | `pem_fuel_cell_water_flux()` | Complete |

## L8: Advanced Methods
| ID | Method | Implementation | Status |
|----|--------|---------------|--------|
| L8.1 | Ghost Fluid Method (GFM) - Dirichlet | `gfm_dirichlet()` | Complete |
| L8.2 | Ghost Fluid Method - Neumann | `gfm_neumann()` | Complete |
| L8.3 | Ghost Fluid Method - Robin | `gfm_robin()` | Complete |
| L8.4 | Immersed Boundary Method - Roma delta | `ibm_delta_roma()` | Complete |
| L8.5 | IBM force spreading | `ibm_spread_force()` | Complete |
| L8.6 | IBM velocity interpolation | `ibm_interpolate_velocity()` | Complete |

## L9: Research Frontiers
| ID | Topic | Documentation | Status |
|----|-------|--------------|--------|
| L9.1 | Nano-scale interface phenomena | Documented in course-tree.md | Partial |
| L9.2 | Non-equilibrium interface thermodynamics | Jump condition framework supports extension | Partial |
| L9.3 | Molecular dynamics of interfaces | GFM/IBM methods provide foundation | Partial |

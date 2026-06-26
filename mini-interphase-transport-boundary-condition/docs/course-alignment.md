# Course Alignment — Interphase Transport & Boundary Conditions

Reference textbook: **Bird, Stewart, Lightfoot (2007) "Transport Phenomena" (2nd ed.)**

## MIT — 2.005 Thermal-Fluids Engineering I & II
| Course Topic | Module Coverage | Section |
|-------------|----------------|---------|
| Boundary conditions for velocity, temperature, concentration | `VelocityBC`, `ThermalBC`, `ConcentrationBC` types | Ch. 2, 10, 19 |
| Interfacial mass transfer | Two-film, penetration, surface renewal theories | Ch. 22 |
| Heat exchanger analysis | LMTD method (condenser example) | Ch. 14 |
| Natural convection correlations | Churchill-Chu Nusselt correlation | Ch. 10 |

## Stanford — ME 346A Heat Transfer
| Course Topic | Module Coverage | Section |
|-------------|----------------|---------|
| Laminar and turbulent convection correlations | Flat plate and pipe Nusselt correlations | L5.5-L5.8 |
| Film condensation (Nusselt theory) | `nusselt_film_condensation()` | L5.10 |
| Pool boiling (Rohsenow) | `rohsenow_pool_boiling()` | L5.11 |
| Droplet evaporation (d^2-law) | `example_droplet_evaporation.c` | L6.6 |

## Berkeley — ME 105 Thermodynamics / ME 106 Fluids
| Course Topic | Module Coverage | Section |
|-------------|----------------|---------|
| Vapor-liquid equilibrium (Raoult's law) | `VaporLiquidEquilibrium`, `interface_raoult_equilibrium()` | L2.6 |
| Henry's law for dilute solutions | `HenryLaw`, `henry_law_at_temperature()` | L2.5 |
| No-slip and slip boundary conditions | `velocity_bc_no_slip()`, `velocity_bc_navier_slip()` | L2.14-15 |

## Michigan — ME 320 Fluid Mechanics / ME 420 Heat Transfer
| Course Topic | Module Coverage | Section |
|-------------|----------------|---------|
| Stress tensor in Newtonian fluids | `StressTensor`, `stress_tensor_newtonian()` | L1.8 |
| Momentum jump conditions at interfaces | `jump_momentum()` | L4.2 |
| Marangoni effects | `interface_marangoni_stress()` | L2.16 |
| Automotive cooling applications | Condenser design example | L6.5 |

## Purdue — ME 505 Heat Transfer / ME 509 Fluids
| Course Topic | Module Coverage | Section |
|-------------|----------------|---------|
| Gnielinski correlation for turbulent pipe flow | `nusselt_pipe_turbulent()` | L5.8 |
| Dittus-Boelter correlation | `nusselt_pipe_turbulent()` | L5.8 |
| Churchill-Chu natural convection | `nusselt_natural_convection_vertical()` | L5.9 |
| Froessling sphere mass transfer | `sherwood_single_sphere()` | L5.13 |

## TU Munich — MW 0854 Heat Transfer
| Course Topic | Module Coverage | Section |
|-------------|----------------|---------|
| Chapman-Enskog gas kinetic theory | `chapman_enskog_diffusivity()` | L3.6 |
| Chung viscosity method | `chung_viscosity()` | L3.7 |
| Eucken thermal conductivity | `eucken_thermal_conductivity()` | L3.8 |
| Stefan problem (phase change) | `stefan_condition()` | L4.7 |

## ETH — 151-0111 Heat Transfer / 151-0103 Fluid Dynamics
| Course Topic | Module Coverage | Section |
|-------------|----------------|---------|
| Jump conditions from integral balances | `jump_mass()`, `jump_momentum()`, `jump_energy()` | L4.1-4.3 |
| Young-Laplace equation | `interface_young_laplace()` | L2.17 |
| Entropy production at interfaces | `jump_entropy()` | L4.5 |
| Level-set and VOF methods | Full numerical module | L5.16-5.26 |

## Tsinghua — 工程热力学 / 传热学 / 流体力学
| Course Topic | Module Coverage | Section |
|-------------|----------------|---------|
| 边界层理论 (Boundary layer theory) | Flat plate Nusselt/Sherwood correlations | L5.5-5.6 |
| 对流传热 (Convective heat transfer) | Pipe flow correlations | L5.7-5.8 |
| 相变换热 (Phase change heat transfer) | Condensation, boiling, Stefan | L5.10-5.11 |
| 传质理论 (Mass transfer theory) | Two-film, penetration, surface renewal | L5.1-5.4 |
| 无因次数群 (Dimensionless groups) | Complete Buckingham Pi analysis | L3.1 |

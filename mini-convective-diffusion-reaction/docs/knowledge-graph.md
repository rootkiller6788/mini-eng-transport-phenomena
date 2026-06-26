# Knowledge Graph ！ mini-convective-diffusion-reaction

## L1: Definitions

| # | Definition | C Implementation | Lean Definition |
|---|-----------|-----------------|-----------------|
| 1 | Peclet number (Pe = uL/D) | PecletNumbers struct, cdr_peclet_mass() | PecletNumber type, peclet.mk |
| 2 | Damkohler number (Da_I..Da_IV) | DamkohlerNumbers struct | DamkohlerI/II types |
| 3 | Thiele modulus | ThieleModulus struct | ThieleModulus type |
| 4 | Reynolds number (Re = rho*u*L/mu) | ReynoldsNumbers struct | ReynoldsNumber type |
| 5 | Schmidt number (Sc = nu/D) | SchmidtNumber struct | SchmidtNumber type |
| 6 | Sherwood number (Sh = kc*L/D) | SherwoodNumbers struct | SherwoodNumber type |
| 7 | Mass transfer coefficient | MassTransferCoefficients struct | ！ |
| 8 | Effectiveness factor (eta) | EffectivenessFactors struct | effectivenessFactorSlab |
| 9 | Hatta number (Ha) | cdr_hatta_number() | ！ |
| 10 | Space time (tau = V/Q) | ReactorSpec.space_time | ！ |
| 11 | Conversion (X = 1 - C_A/C_A0) | cdr_conversion_first_order() | conversion function |
| 12 | Half-life (t1/2) | cdr_half_life_first_order() | halfLifeFirstOrder function |
| 13 | Selectivity (S_B/C) | cdr_selectivity_parallel() | ！ |
| 14 | Yield (Y_B = C_B/C_A0) | cdr_yield_intermediate_series() | ！ |
| 15 | Arrhenius parameters (A, Ea) | ArrheniusParams struct | arrheniusRate function |
| 16 | Penetration depth | cdr_diffusion_penetration_depth() | ！ |
| 17 | Colburn j-factor (j_D) | ColburnJFactors struct | ！ |
| 18 | Weisz-Prater parameter | cdr_weisz_prater_criterion() | ！ |

## L2: Core Concepts

| # | Concept | Implementation |
|---|---------|---------------|
| 1 | Transport regime classification | TransportRegime enum, cdr_classify_regime() |
| 2 | Diffusion regimes (molecular, Knudsen, surface) | DiffusionRegime enum |
| 3 | Reaction network types | ReactionNetworkType enum |
| 4 | RTD concepts (E(t), F(t), mean, variance) | ResidenceTimeDistribution struct |
| 5 | Two-film theory | cdr_two_film_overall_gas(), cdr_two_film_flux() |
| 6 | Catalyst geometry effects | CatalystGeometry enum |
| 7 | Diffusion-reaction coupling (Thiele approach) | cdr_effectiveness_general() |
| 8 | Micromixing (segregation vs. max mixedness) | cdr_segregation_model_conversion() |
| 9 | Mass transfer enhancement by reaction | cdr_enhancement_factor_homogeneous() |
| 10 | Dispersion model (non-ideal PFR) | cdr_dispersion_conversion() |

## L3: Engineering Quantities

| # | Quantity | Typical Range | Implementation |
|---|----------|--------------|---------------|
| 1 | Gas diffusion coefficient | 1e-5 to 1e-4 m2/s | cdr_diffusion_fuller() |
| 2 | Liquid diffusion coefficient | 1e-10 to 1e-9 m2/s | cdr_diffusion_wilke_chang() |
| 3 | Knudsen diffusivity | 1e-8 to 1e-6 m2/s | cdr_diffusion_knudsen() |
| 4 | Effective diffusivity in catalyst | 1e-7 to 1e-5 m2/s | cdr_diffusion_effective() |
| 5 | Mass transfer coefficient (gas) | 1e-3 to 1e-1 m/s | cdr_mass_transfer_coeff_from_sh() |
| 6 | Mass transfer coefficient (liquid) | 1e-5 to 1e-3 m/s | cdr_mass_transfer_penetration() |
| 7 | Schmidt number (gases) | 0.5 to 2.0 | cdr_schmidt() |
| 8 | Schmidt number (liquids) | 100 to 10000 | cdr_schmidt() |
| 9 | Activation energy (typical) | 40 to 200 kJ/mol | cdr_activation_energy_from_two_points() |
| 10 | Thiele modulus (practical range) | 0.1 to 50 | cdr_thiele_modulus() |

## L4: Conservation Laws

| # | Law | C Implementation | Lean Theorem |
|---|-----|-----------------|-------------|
| 1 | Species mass balance (CDR equation) | cdr_dispersion_outlet_concentration() | fickSecondLawFTCS |
| 2 | Ficks first law | cdr_fick_first_law_1d/3d() | FickFirstLaw |
| 3 | Ficks second law | cdr_fick_second_law_1d_step() | fickSecondLawFTCS |
| 4 | Arrhenius law | cdr_arrhenius_rate() | arrheniusRate |
| 5 | Mass-action kinetics | cdr_power_law_rate() | firstOrderRate |
| 6 | Michaelis-Menten rate law | cdr_michaelis_menten_rate() | ！ |
| 7 | CSTR design equation | cdr_cstr_volume_first_order() | cstrOutletConcentration |
| 8 | PFR design equation | cdr_pfr_volume_first_order() | pfrOutletConcentration |
| 9 | Effectiveness factor PDE solution | cdr_effectiveness_slab/cylinder/sphere() | effectivenessFactorSlab |
| 10 | Two-film resistance-in-series | cdr_two_film_overall_gas() | overallMassTransferCoeff |

## L5: Engineering Methods

| # | Method | Implementation |
|---|--------|---------------|
| 1 | Fuller-Schettler-Giddings (gas diffusivity) | cdr_diffusion_fuller() |
| 2 | Wilke-Chang (liquid diffusivity) | cdr_diffusion_wilke_chang() |
| 3 | Bosanquet formula (transition diffusion) | cdr_diffusion_bosanquet() |
| 4 | CSTR design (1st & 2nd order) | cdr_cstr_volume_first/second_order() |
| 5 | PFR design | cdr_pfr_volume_first/second_order() |
| 6 | RTD parameter estimation | cdr_rtd_mean_time(), cdr_rtd_peclet_from_variance() |
| 7 | Segregation model | cdr_segregation_model_conversion() |
| 8 | Weisz-Prater criterion | cdr_weisz_prater_criterion() |
| 9 | Mears criterion | cdr_mears_external_criterion() |
| 10 | HTU/NTU method | cdr_ntu_gas_phase(), cdr_htu_gas_phase() |
| 11 | Chilton-Colburn analogy | cdr_j_factor_mass() |
| 12 | Dispersion model with Danckwerts BC | cdr_dispersion_outlet_concentration() |
| 13 | Tanks-in-series model | cdr_rtd_tanks_in_series_E() |
| 14 | Penetration theory (gas-liquid kL) | cdr_mass_transfer_penetration() |
| 15 | Kinetic parameter estimation (Ea) | cdr_activation_energy_from_two_points() |
| 16 | Optimal pellet size | cdr_optimal_pellet_radius() |

## L6: Engineering Problems

| # | Problem | Example/Solution |
|---|---------|-----------------|
| 1 | Reactor sizing (CSTR vs PFR) | examples/example_cstr_design.c |
| 2 | Catalyst pellet effectiveness | examples/example_catalyst_pellet.c |
| 3 | Axial dispersion in tubular reactor | examples/example_dispersion.c |
| 4 | RTD-based reactor diagnosis | examples/example_rtd_analysis.c |
| 5 | Diffusion limitation diagnosis | cdr_diagnose_limitations() |
| 6 | Recycle reactor optimization | cdr_pfr_recycle_volume_factor() |
| 7 | Series reaction intermediate maximization | cdr_optimal_time_series() |

## L7: Applications

| # | Application | Keywords |
|---|------------|----------|
| 1 | Industrial VCM catalyst optimization | NASA, chemical, industrial |
| 2 | Wastewater treatment UV disinfection | EPA, Detroit |
| 3 | Pharmaceutical process validation (RTD) | FDA, validation |
| 4 | Monolith catalytic converter design | Automotive, Toyota |
| 5 | Gas absorption column design | Chemical, packed column |
| 6 | Fermentation kinetics (Michaelis-Menten) | Bioprocess, fermentation |

## L8: Advanced Methods

| # | Method | Implementation |
|---|--------|---------------|
| 1 | Multi-component diffusion (Stefan-Maxwell) | MulticomponentDiffusion struct |
| 2 | Generalized Thiele modulus (nth-order) | cdr_effectiveness_nth_order() |
| 3 | Michaelis-Menten with pore diffusion | cdr_effectiveness_michaelis_menten() |
| 4 | Dispersion model with reaction | cdr_dispersion_outlet_concentration() |
| 5 | Bessel function effectiveness (cylinder) | cdr_effectiveness_cylinder() |
| 6 | Langmuir-Hinshelwood kinetics | cdr_langmuir_hinshelwood_uni/bi() |

## L9: Research Frontiers

| # | Topic | Status |
|---|-------|--------|
| 1 | Micro-reactor transport (Pe < 1) | Documented |
| 2 | Nano-catalysis (configurational diffusion) | Documented |
| 3 | Single-molecule reaction kinetics | Documented |
| 4 | Non-Fickian diffusion | Documented |
| 5 | Population balance models | Documented |

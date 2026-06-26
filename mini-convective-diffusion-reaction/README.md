# mini-convective-diffusion-reaction

## Module Status: COMPLETE ✅

- L1-L6: Complete
- L7: Complete (6 applications)
- L8: Partial (6 advanced topics with implementations)
- L9: Partial (5 frontier topics documented)

**Total Score: 16/18 → COMPLETE**

## Overview

Convection-Diffusion-Reaction transport phenomena module. Implements the core of chemical reaction engineering and transport theory based on Bird-Stewart-Lightfoot (2007) and Levenspiel (1999).

## Nine-Layer Knowledge Coverage

| Level | Name | Status | Items |
|-------|------|--------|-------|
| **L1** | Definitions | COMPLETE | 18 definitions (Peclet, Damkohler, Thiele, Schmidt, Sherwood, conversion, selectivity, yield, effectiveness, Arrhenius...) |
| **L2** | Core Concepts | COMPLETE | 10 concepts (transport regimes, diffusion mechanisms, RTD, two-film theory, catalyst geometry effects...) |
| **L3** | Engineering Quantities | COMPLETE | 10 quantity ranges (D_AB gas/liquid, k_c gas/liquid, Sc, Ea, Thiele modulus...) |
| **L4** | Conservation Laws | COMPLETE | 10 laws (CDR equation, Fick 1st/2nd, Arrhenius, mass-action, MM, CSTR/PFR design, effectiveness PDE, two-film) |
| **L5** | Engineering Methods | COMPLETE | 16 methods (Fuller, Wilke-Chang, CSTR/PFR design, RTD estimation, Weisz-Prater, Mears, HTU/NTU, Chilton-Colburn...) |
| **L6** | Engineering Problems | COMPLETE | 7 problems + 4 end-to-end examples |
| **L7** | Applications | COMPLETE | 6 applications (VCM catalyst, EPA disinfection, FDA validation, Toyota catalytic converter, gas absorption, fermentation) |
| **L8** | Advanced Methods | PARTIAL | 6 topics (multi-component diffusion, generalized Thiele, MM+pore diffusion, Bessel effectiveness, LH kinetics, dispersion with reaction) |
| **L9** | Research Frontiers | PARTIAL | 5 topics documented (micro-reactors, nano-catalysis, single-molecule kinetics, non-Fickian diffusion, population balances) |

## Core Definitions (L1)

| # | Definition | Formula | Function |
|---|-----------|---------|----------|
| 1 | Peclet number (mass) | Pe_m = u·L/D_AB | cdr_peclet_mass() |
| 2 | Damköhler-I | Da_I = k·τ | cdr_damkohler_I() |
| 3 | Damköhler-II | Da_II = k·L²/D | cdr_damkohler_II() |
| 4 | Thiele modulus | φ = L_c·√(k/D_eff) | cdr_thiele_modulus() |
| 5 | Reynolds number | Re = ρ·u·L/μ | cdr_reynolds() |
| 6 | Schmidt number | Sc = μ/(ρ·D) | cdr_schmidt() |
| 7 | Sherwood number | Sh = k_c·L/D | cdr_sherwood_sphere() |
| 8 | Effectiveness factor | η = tanh(φ)/φ (slab) | cdr_effectiveness_slab() |
| 9 | Hatta number | Ha = √(k·D)/k_L | cdr_hatta_number() |
| 10 | Conversion (1st-order) | X = 1 − exp(−k·t) | cdr_conversion_first_order() |
| 11 | Half-life (1st-order) | t₁/₂ = ln(2)/k | cdr_half_life_first_order() |
| 12 | Arrhenius rate | k = A·exp(−Ea/RT) | cdr_arrhenius_rate() |
| 13 | Fick's 1st law | J = −D·dC/dx | cdr_fick_first_law_1d() |
| 14 | Fick's 2nd law | ∂C/∂t = D·∂²C/∂x² | cdr_fick_second_law_1d_step() |
| 15 | CSTR design (1st-order) | V = Q·X/(k·(1−X)) | cdr_cstr_volume_first_order() |
| 16 | PFR design (1st-order) | V = Q·ln(1/(1−X))/k | cdr_pfr_volume_first_order() |
| 17 | Weisz-Prater criterion | C_WP = r_obs·L_c²/(D_eff·C_s) | cdr_weisz_prater_criterion() |
| 18 | Selectivity (parallel) | S_B/C = (k₁/k₂)·C_A^(n₁−n₂) | cdr_selectivity_parallel() |

## Core Theorems (L4, Lean 4)

| # | Theorem | Lean Proof |
|---|---------|-----------|
| 1 | Fick FTCS preserves uniform concentration | fick_ftcs_preserves_uniform |
| 2 | Arrhenius: zero Ea → k = A | arrhenius_zero_activation |
| 3 | Effectiveness: φ=0 → η=1 | effectiveness_no_limitation |
| 4 | PFR ≠ CSTR outlet concentration | pfr_more_efficient_structural |
| 5 | Conversion with no reaction = 0 | conversion_no_reaction |
| 6 | V_CSTR/V_PFR > 1 at X=0.5 | volume_ratio_above_one |
| 7 | Overall K_OG < individual k_G | overall_coeff_less_than_individual |
| 8 | Half-life positive for k > 0 | half_life_positive |
| 9 | Gas-side resistance dominates → K_OG = k_G | overall_approx_gas_side |
| 10 | All transport regimes distinct | regimes_distinct |

## Core Algorithms (L5)

| # | Algorithm | Function | Complexity |
|---|----------|----------|------------|
| 1 | Fuller-Schettler-Giddings gas D estimation | cdr_diffusion_fuller() | O(1) |
| 2 | Wilke-Chang liquid D estimation | cdr_diffusion_wilke_chang() | O(1) |
| 3 | CSTR volume (1st & 2nd order) | cdr_cstr_volume_* | O(1) |
| 4 | PFR volume (1st & 2nd order) | cdr_pfr_volume_* | O(1) |
| 5 | RTD mean & variance (trapezoidal) | cdr_rtd_mean_time/variance | O(n) |
| 6 | Peclet number from RTD variance | cdr_rtd_peclet_from_variance | O(1) |
| 7 | Segregation model (trapezoidal) | cdr_segregation_model_conversion | O(n) |
| 8 | Dispersion model with Danckwerts BC | cdr_dispersion_outlet_concentration | O(1) |
| 9 | Effectiveness factor (Bessel-based) | cdr_effectiveness_cylinder | O(1) |
| 10 | Optimal pellet radius (bisection) | cdr_optimal_pellet_radius | O(log(1/ε)) |
| 11 | Ea from two rate measurements | cdr_activation_energy_from_two_points | O(1) |

## Classic Problems (L6)

| # | Problem | Example |
|---|---------|---------|
| 1 | CSTR vs PFR reactor sizing | examples/example_cstr_design.c |
| 2 | Catalyst pellet effectiveness optimization | examples/example_catalyst_pellet.c |
| 3 | Axial dispersion in tubular reactors | examples/example_dispersion.c |
| 4 | RTD-based reactor diagnosis | examples/example_rtd_analysis.c |

## Nine-School Course Mapping

| School | Key Course | Module Coverage |
|--------|-----------|----------------|
| **MIT** | 2.005 Thermal-Fluids, 10.37 Kinetics & Reactors | CDR equation, CSTR/PFR/RTD |
| **Stanford** | ME 346A Heat/Mass Transfer, ME 451 Transport | Sherwood, Chilton-Colburn, multi-component |
| **Berkeley** | ME 105/106, CHM ENG 150A/160 | Transport fundamentals, reactor design |
| **Michigan** | ME 320/420, AERO 533 | Fluid mechanics, mass transfer, combustion |
| **Purdue** | ME 505/509/597 | Heat/mass transfer, dispersion, plasma |
| **TU Munich** | MW 0798/0854, CH 4005 | Thermodynamics, mass transfer, reactor eng. |
| **ETH Zurich** | 151-0103/0111, 529-0634 | Fluid dynamics, heat transfer, reaction eng. |
| **Tsinghua** | Engineering Thermodynamics, Heat Transfer, Reaction Engineering | All core topics |

## Building and Testing

`ash
make          # Build library, tests, and examples
make test     # Build and run all tests
make count    # Count lines of code
make clean    # Remove build artifacts
`

## File Structure

`
mini-convective-diffusion-reaction/
├── Makefile
├── README.md                        ← This file
├── include/
│   ├── cdr_core.h                   # Dimensionless numbers, CDRSystem
│   ├── cdr_diffusion.h              # Fick's laws, D_eff models
│   ├── cdr_reaction.h               # Kinetics, rate laws
│   ├── cdr_reactor.h                # CSTR, PFR, RTD, dispersion
│   ├── cdr_mass_transfer.h          # Sh correlations, two-film
│   └── cdr_effectiveness.h          # eta(phi), diagnostics
├── src/
│   ├── cdr_core.c
│   ├── cdr_diffusion.c
│   ├── cdr_reaction.c
│   ├── cdr_reactor.c
│   ├── cdr_mass_transfer.c
│   ├── cdr_effectiveness.c
│   └── cdr.lean                     # Lean 4 formalization
├── tests/
│   ├── test_core.c
│   ├── test_diffusion.c
│   ├── test_reaction.c
│   ├── test_reactor.c
│   └── test_effectiveness.c
├── examples/
│   ├── example_cstr_design.c
│   ├── example_catalyst_pellet.c
│   ├── example_dispersion.c
│   └── example_rtd_analysis.c
├── demos/
│   └── demo_conversion.c
├── benches/
│   └── bench_cdr.c
└── docs/
    ├── knowledge-graph.md
    ├── coverage-report.md
    ├── gap-report.md
    ├── course-alignment.md
    └── course-tree.md
`

## References

1. Bird, R.B., Stewart, W.E., Lightfoot, E.N. (2007). *Transport Phenomena*, 2nd ed. Wiley.
2. Levenspiel, O. (1999). *Chemical Reaction Engineering*, 3rd ed. Wiley.
3. Fogler, H.S. (2016). *Elements of Chemical Reaction Engineering*, 5th ed. Pearson.
4. Cussler, E.L. (2009). *Diffusion: Mass Transfer in Fluid Systems*, 3rd ed. Cambridge.
5. Treybal, R.E. (1980). *Mass Transfer Operations*, 3rd ed. McGraw-Hill.

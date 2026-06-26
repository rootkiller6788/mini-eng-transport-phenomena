# Course Tree ¡ª mini-convective-diffusion-reaction

## Prerequisite Chain

This module depends on:
- mini-eng-thermodynamics (0): Equilibrium, energy balances
- mini-eng-fluid-mechanics (1): Fluid flow, Reynolds, boundary layers
- mini-eng-heat-mass-transfer (2): Diffusion, heat/mass analogy

And feeds into:
- mini-eng-reaction-engineering (9): Advanced reactor networks, multiphase reactors

## Internal Dependency Tree

cdr_core.h is the root dependency:
- cdr_diffusion.h: Fick laws, effective diffusivity models
- cdr_reaction.h: Arrhenius, rate laws, kinetic analysis
- cdr_reactor.h: CSTR, PFR, RTD, dispersion (depends on cdr_reaction.h)
- cdr_mass_transfer.h: Sherwood correlations, two-film theory (depends on cdr_diffusion.h)
- cdr_effectiveness.h: Catalyst effectiveness, Weisz-Prater, Mears criteria

## Knowledge Prerequisites by Level

| Level | Prerequisites | This Module |
|-------|--------------|-------------|
| L1 | Calculus, ODEs | Dimensionless numbers, structural definitions |
| L2 | Transport concepts | Regime concepts, rate-limiting steps |
| L3 | Physical chemistry | Engineering quantity ranges |
| L4 | PDEs, reaction kinetics | CDR mass balance, Fick+Arrhenius formalization |
| L5 | Numerical methods | Design algorithms, estimation methods |
| L6 | Reactor engineering | Classic design/analysis problems |
| L7 | Industrial practice | Application case studies |
| L8 | Advanced math (Bessel) | Generalized Thiele, multi-component formulations |
| L9 | Research literature | Micro/nano transport frontiers (documented) |

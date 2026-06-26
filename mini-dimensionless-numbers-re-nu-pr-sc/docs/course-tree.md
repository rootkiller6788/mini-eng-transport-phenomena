# Course Tree — mini-dimensionless-numbers-re-nu-pr-sc

## Prerequisite Dependencies

This module requires knowledge of:

```
mini-eng-fluid-mechanics (Re, flow regimes, N-S equations)
    │
    ├── mini-eng-heat-mass-transfer (conduction, convection, Nu, Pr)
    │       │
    │       └── mini-dimensionless-numbers-re-nu-pr-sc  ← THIS MODULE
    │               │
    │               ├── mini-momentum-heat-mass-analogy
    │               ├── mini-convective-diffusion-reaction
    │               ├── mini-interphase-transport-boundary-condition
    │               └── mini-transport-in-micro-nano-channel
    │
    └── mini-eng-thermodynamics (energy balance, fluid properties)
```

## Internal Dependency Tree

```
dimensionless_numbers.h/c (L1-L2 core)
    │
    ├── reynolds_number.h/c (L1-L5: Re-specific)
    │       ├── Flow regimes, critical Re
    │       ├── Friction factors (laminar → Churchill)
    │       └── Boundary layer thickness
    │
    ├── nusselt_number.h/c (L1-L6: Nu-specific)
    │       ├── Internal flow (pipe, duct, annulus)
    │       ├── External flow (flat plate, cylinder, sphere)
    │       ├── Natural convection (plate, enclosure)
    │       └── Energy equation non-dimensionalization
    │
    ├── prandtl_schmidt.h/c (L1-L5: Pr, Sc, Le)
    │       ├── Temperature-dependent properties
    │       ├── Diffusivity models
    │       ├── BL thickness ratios
    │       └── Transport analogies
    │
    ├── buckingham_pi.h/c (L1-L5, L7: Pi theorem)
    │       ├── Dimension vectors
    │       ├── Group verification
    │       ├── Matrix method
    │       ├── Model scaling
    │       └── Scaling laws
    │
    └── transport_correlations.h/c (L5-L7: applications)
            ├── Forced/natural convection engine
            ├── Mass transfer engine
            ├── Heat exchanger sizing
            ├── Pipe flow design
            ├── Electronics cooling
            ├── DC motor cooling
            ├── Air-cooled condenser
            └── Mars rover thermal

analogy_theorems.lean (L4: formal proofs)
    ├── ReynoldsNumber, NusseltNumber, PrandtlNumber, SchmidtNumber (structures)
    ├── FlowRegime classification (inductive, exhaustive proof)
    ├── ColburnFactors verification
    ├── Boundary layer theory (structures)
    ├── NSRegime (Euler/Stokes/Navier-Stokes)
    ├── LewisNumber, WeberNumber
    ├── Rayleigh-Bénard onset
    └── Darcy friction laminar (bounded proof)
```

## Learning Path

1. **Start with L1**: Read `dimensionless_numbers.h` — understand what Re, Nu, Pr, Sc mean
2. **Proceed to L2**: `buckingham_pi.c` — understand why these groups exist (Pi theorem)
3. **Deepen with L3**: `prandtl_schmidt.c` — see how properties vary with temperature
4. **Understand L4**: `reynolds_number.c` (N-S) + `nusselt_number.c` (energy eq) — see where these numbers come from
5. **Apply L5**: `transport_correlations.c` — use correlations for real problems
6. **Solve L6**: Run the examples — pipe flow, flat plate, natural convection, mass transfer
7. **Connect L7**: See how Boeing, Toyota, Tesla, NASA, Apple use these numbers daily
8. **Verify L4**: Read `analogy_theorems.lean` — formalized theorems

## Key Insight

All transport phenomena — momentum, heat, and mass — are governed by
the same mathematical structure. The dimensionless numbers capture
the essential physics:

| Transport Mode | Diffusivity | Eqn | BC at wall | Dimensionless |
|---------------|-------------|-----|-----------|---------------|
| Momentum | ν [m²/s] | N-S | u=0 (no-slip) | Re = UL/ν, Cf = 2τ_w/(ρU²) |
| Heat | α [m²/s] | Energy | T=T_w or q" | Nu = hL/k, Pr = ν/α |
| Mass | D [m²/s] | Species | C=C_w or N" | Sh = h_mL/D, Sc = ν/D |

The analogy: **if you know the velocity field, you can predict both
heat and mass transfer**.

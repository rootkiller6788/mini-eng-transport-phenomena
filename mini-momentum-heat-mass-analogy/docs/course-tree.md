# Course Tree — Mini Momentum-Heat-Mass Analogy

## Prerequisite Dependency Tree

```
mini-momentum-heat-mass-analogy
├── P0: Transport coefficients (μ, k, D)
│   ├── Kinetic theory of gases (Chapman-Enskog)
│   ├── Sutherland's law for gas viscosity
│   ├── Power-law model for liquid viscosity
│   ├── Eucken formula for gas conductivity
│   └── Fuller/Wilke-Chang diffusivity correlations
│
├── P1: Conservation laws
│   ├── Newton's law of viscosity: τ = -μ·dv/dy
│   ├── Fourier's law of conduction: q = -k·dT/dy
│   ├── Fick's first law of diffusion: J = -D·dC/dy
│   └── Generalized transport flux: Flux = -Γ·∇φ
│
├── P2: Dimensionless numbers
│   ├── Reynolds number (Re = ρvL/μ)
│   ├── Prandtl number (Pr = ν/α)
│   ├── Schmidt number (Sc = ν/D)
│   ├── Lewis number (Le = α/D = Sc/Pr)
│   ├── Nusselt number (Nu = hL/k)
│   ├── Sherwood number (Sh = k_cL/D)
│   └── Stanton number (St = Nu/(Re·Pr))
│
├── P3: Reynolds Analogy (Pr = Sc = 1)
│   ├── St = f/2
│   ├── Eddy diffusivity concept
│   └── Limitations (laminar flow, Pr ≠ 1)
│
├── P4: Chilton-Colburn Analogy
│   ├── j_H = St_h·Pr^(2/3) = f/2
│   ├── j_D = St_m·Sc^(2/3) = f/2
│   ├── Validity: 0.6 < Pr < 60, turbulent
│   └── Practical: measure f → predict h, k_c
│
├── P5: Boundary layer theory
│   ├── Blasius velocity BL solution
│   ├── Pohlhausen thermal BL solution
│   ├── Concentration BL (analogous)
│   ├── BL thickness ratios: δ_T/δ ≈ Pr^(-1/3)
│   └── Entrance region effects
│
├── P6: Internal flow analogies
│   ├── Laminar: f=64/Re, Nu=3.66 (const T_w)
│   ├── Turbulent: Blasius + Dittus-Boelter
│   ├── Gnielinski/Petukhov improvements
│   ├── Non-circular ducts
│   └── Developing flow (Hausen)
│
├── P7: Advanced analogy models
│   ├── Prandtl-Taylor (two-layer)
│   ├── von Karman (three-layer)
│   └── Lyon-Martinelli (liquid metals, Pr<<1)
│
└── P8: Engineering applications
    ├── Heat exchanger sizing
    ├── Cooling tower (Merkel)
    ├── Electronics cooling
    ├── Reactor wall mass transfer
    └── Quick engineering estimation
```

## Dependencies from other mini-modules

- **mini-dimensionless-numbers-re-nu-pr-sc**: Re, Pr, Sc, Nu, Sh definitions
- **mini-interphase-transport-boundary-condition**: Boundary conditions for transport
- **mini-convective-diffusion-reaction**: Convective transport equations
- **mini-transport-in-micro-nano-channel**: Micro/nano-scale transport (L9)
- **mini-eng-heat-mass-transfer** (module 2): Heat exchanger concepts
- **mini-eng-fluid-mechanics** (module 1): BL theory, pipe flow

## Self-contained knowledge: This module can be studied independently

The module includes all prerequisite definitions (L1), builds up the
analogy framework (L2-L4), provides computational methods (L5),
solves canonical problems (L6), and demonstrates applications (L7-L8).

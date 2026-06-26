# Course Tree — Interphase Transport & Boundary Conditions

## Prerequisites Dependency Tree

```
mini-interphase-transport-boundary-condition
├── Thermodynamics
│   ├── First law (energy conservation) → energy jump conditions
│   ├── Second law (entropy) → entropy jump condition
│   ├── Phase equilibrium (Raoult, Henry) → interfacial concentration BCs
│   └── Equation of state → density, property estimation
├── Fluid Mechanics
│   ├── Navier-Stokes equations → momentum jump conditions
│   ├── Boundary layer theory → convective correlations
│   ├── No-slip condition → velocity BCs
│   └── Surface tension → Young-Laplace, Marangoni
├── Heat Transfer
│   ├── Fourier's law → thermal BCs
│   ├── Convection correlations → Nusselt number
│   ├── Phase change heat transfer → Stefan condition
│   └── Radiation → radiation BC
├── Mass Transfer
│   ├── Fick's law → concentration BCs
│   ├── Diffusion coefficients → Chapman-Enskog, property estimation
│   ├── Mass transfer coefficients → Sherwood correlations
│   └── Reaction-diffusion → Hatta number, enhancement factor
├── Numerical Methods
│   ├── Finite difference/volume → Level-set advection
│   ├── Interface tracking → VOF/PLIC, Level-set
│   ├── Immersed boundary method → IBM force spreading
│   └── Ghost fluid method → sharp interface BCs
└── Engineering Mathematics
    ├── Vector calculus → stress tensor, flux vectors
    ├── Buckingham Pi theorem → dimensionless groups
    └── Integral conservation laws → jump condition derivation
```

## Internal Module Dependencies

```
interphase_types.h (foundation)
├── interphase_config.h (constants)
├── interphase_boundary_conditions.h
│   └── interphase_types.h
├── interphase_transport.h
│   └── interphase_types.h
├── interphase_jump_conditions.h
│   └── interphase_types.h
└── interphase_numerical.h
    └── interphase_types.h
```

## Downstream Dependencies (modules that depend on this)

- mini-computational-transport-cfd-basics (needs interface tracking methods)
- mini-convective-diffusion-reaction (needs mass transfer boundary conditions)
- mini-momentum-heat-mass-analogy (needs transport coefficient analogies)
- mini-transport-in-micro-nano-channel (needs slip boundary conditions, nano-scale interfaces)
- mini-boiling-condensation-two-phase (needs phase change jump conditions)
- mini-multiphysics-coupling (needs interface coupling framework)

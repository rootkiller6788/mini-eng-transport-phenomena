# Mini Engineering Transport Phenomena

A collection of **from-scratch, zero-dependency C implementations** of university-level transport phenomena — the unified framework of momentum, heat, and mass transfer. Each module maps to MIT (and other top-tier university) courses, bridging theory and practice by translating textbook equations into runnable C code.

## Sub-Modules

| Sub-Module | Topics | Key Courses |
|--------|--------|-------------|
| [mini-computational-transport-cfd-basics](mini-computational-transport-cfd-basics/) | FDM/FVM discretization, upwind/power-law/hybrid schemes, convection-diffusion solvers, SIMPLE algorithm, lid-driven cavity | MIT 2.29, MIT 16.90 |
| [mini-convective-diffusion-reaction](mini-convective-diffusion-reaction/) | Peclet/Damköhler/Thiele numbers, Fick's law, Stefan-Maxwell diffusion, catalyst effectiveness factor, Weisz-Prater criterion | MIT 10.50, MIT 10.37 |
| [mini-dimensionless-numbers-re-nu-pr-sc](mini-dimensionless-numbers-re-nu-pr-sc/) | Buckingham Pi theorem, Reynolds/Prandtl/Schmidt/Nusselt numbers, convective heat transfer correlations for canonical geometries | MIT 10.50, MIT 2.51 |
| [mini-interphase-transport-boundary-condition](mini-interphase-transport-boundary-condition/) | No-slip/slip BCs, interphase jump conditions (mass/momentum/energy), heat/mass transfer boundary conditions, surface tension | MIT 10.50, Stanford ME 351 |
| [mini-momentum-heat-mass-analogy](mini-momentum-heat-mass-analogy/) | Reynolds analogy (St = f/2), Chilton-Colburn analogy, velocity/thermal/concentration boundary layers, dimensionless group coupling | MIT 10.50, MIT 2.51 |
| [mini-transport-in-micro-nano-channel](mini-transport-in-micro-nano-channel/) | Knudsen flow regimes (continuum to free molecular), slip models, electroosmosis, electrophoresis, streaming potential, EDL models | MIT 2.25, Stanford ME 457 |

## Design Philosophy

- **Zero external dependencies** — pure C (C99/C11), only `libc` and `libm`
- **Self-contained modules** — each directory has its own `Makefile`, `include/`, `src/`, `examples/`, `demos/`, `tests/`
- **Theory-to-code mapping** — every module includes `docs/` with course-alignment notes
- **Practical demos** — SIMPLE solver, effectiveness factor calculator, Nu correlation explorer, jump condition verifier, analogy plots, Knudsen regime classifier

## Building

Each module is standalone. Navigate to a module directory and run:

```bash
cd mini-computational-transport-cfd-basics
make all    # build everything
make test   # run tests
```

Requires **GCC** and **GNU Make**.

## Project Structure

```
mini-eng-transport-phenomena/
├── mini-computational-transport-cfd-basics/       # FDM/FVM, convection-diffusion solvers, SIMPLE
├── mini-convective-diffusion-reaction/            # CDR transport, Fick's law, catalyst effectiveness
├── mini-dimensionless-numbers-re-nu-pr-sc/        # Buckingham Pi, Re/Nu/Pr/Sc, Nu correlations
├── mini-interphase-transport-boundary-condition/  # Boundary conditions, interphase jump balances
├── mini-momentum-heat-mass-analogy/               # Reynolds/Chilton-Colburn analogies, BL theory
└── mini-transport-in-micro-nano-channel/          # Knudsen regimes, electrokinetics, EDL models
```

## License

MIT

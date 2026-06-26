# mini-transport-in-micro-nano-channel

Transport phenomena in micro and nano scale channels — where continuum assumptions break down and surface effects dominate bulk behavior.

## Module Status: COMPLETE ✅

- L1-L6: Complete
- L7: Partial+ (6 applications: nanopore sensing, EOF pumping, membrane transport, water splitting, streaming potential, ion separation)
- L8: Partial (7/7 advanced topics implemented: Burnett equations, BGK model, unified mass flow, electroconvection, hindered diffusion, effective screening, entropic springs)
- L9: Partial (3 frontiers documented in Lean structures: non-Fourier conduction, quantum confinement, topological transport)

## Nine-Layer Knowledge Coverage

| Level | Name | Status | Count |
|-------|------|--------|-------|
| L1 | Definitions | **Complete** | 11 definitions |
| L2 | Core Concepts | **Complete** | 10 concepts |
| L3 | Engineering Quantities | **Complete** | 9 quantities |
| L4 | Conservation Laws | **Complete** | 9 laws |
| L5 | Engineering Methods | **Complete** | 11 methods |
| L6 | Engineering Problems | **Complete** | 9 problems |
| L7 | Applications | **Partial+** | 6 applications |
| L8 | Advanced Methods | **Partial+** | 7 topics |
| L9 | Frontiers | **Partial** | 3 documented |

## Core Definitions (L1)

| Definition | Symbol | Formula |
|-----------|--------|---------|
| Knudsen number | Kn | λ / L_char |
| Mean free path (hard-sphere) | λ | 1/(√2·π·d²·n) |
| Slip length (Maxwell) | L_s | ((2-σ)/σ)·λ |
| Temperature jump length | ζ_T | ((2-α_T)/α_T)·(2γ/((γ+1)Pr))·λ |
| Debye length | κ⁻¹ | sqrt(ε·kB·T/(2·e²·z²·c_bulk·NA)) |
| Zeta potential | ζ | ψ at shear plane |
| Hydraulic diameter | D_h | 4A/P_wetted |
| Bjerrum length | l_B | e²/(4π·ε·kB·T) |
| Gouy-Chapman length | λ_GC | e/(2π·l_B·|σ_s|) |
| Helmholtz-Smoluchowski velocity | u_HS | -(εζ/η)·E |
| EOF mobility | μ_eo | -εζ/η |

## Core Theorems (L4)

1. **Poiseuille no-slip theorem**: u(H/2) = 0 at walls for no-slip flow
2. **Poiseuille profile**: u(y) = (dpdx/(2η))·(H²/4 - y²)
3. **Maxwell slip**: u_slip = L_s·(∂u/∂n)|_wall
4. **Debye-Huckel**: ψ(y) = ψ₀·e^{-κy} for low potentials
5. **Grahame equation**: σ_s = √(8εkBTc₀NA)·sinh(zeψ₀/(2kBT))
6. **Nernst-Planck**: J = -D·(∇c + (ze/(kBT))·c·∇φ)
7. **Goldman-Hodgkin-Katz**: Ionic current under constant field
8. **Enhancement factor**: Q_slip/Q_no-slip = 1 + 6·L_s/H
9. **Poiseuille number slip**: Po = 96/(1 + 6·L_s/H)

## Core Algorithms (L5)

- **solve_pressure_driven_slip_flow**: Analytical solution of NS + Maxwell slip
- **solve_second_order_slip_flow**: Deissler 2nd-order slip
- **solve_electroosmotic_flow**: EOF with finite EDL
- **compute_poisson_boltzmann_1d**: Gouy-Chapman + Debye-Huckel solvers
- **compute_nernst_planck_steady_1d**: GHK constant-field ion flux
- **compute_rarefied_mass_flow_rate**: Beskok-Karniadakis unified model
- **compute_slip_correction_mass_flow_gas**: Compressible gas slip flow
- **compute_eof_profile_full**: Combined EOF + pressure flow
- **compute_concentration_polarization**: Film model CP profile
- **compute_ghk_current_voltage**: Multi-ion GHK current equation
- **compute_transition_regime_mass_flow**: Interpolation across all Kn

## Classic Problems (L6)

1. Pressure-driven flow in a 1 μm channel with slip
2. Electroosmotic flow with varying Debye length
3. Ionic current through a nanopore (blockade sensing)
4. Streaming potential measurement for zeta determination
5. Concentration polarization at RO membrane
6. Ion current rectification in conical nanopore
7. Knudsen minimum in rarefied gas flow
8. Hindered diffusion of macromolecules in nanopores
9. Entropic trapping of polymers in nanofluidic constrictions

## Nine-School Course Mapping

| School | Relevant Course | Topics |
|--------|----------------|--------|
| MIT | 2.25 Advanced Fluid Mechanics | Knudsen regime, slip flow |
| Stanford | ME 346A Heat Transfer | Microscale transport |
| Berkeley | ME 106 Fluid Mechanics | Low-Re flow, microfluidics |
| Michigan | ME 320 Fluid Mechanics | Microchannel flow |
| Purdue | ME 505 Heat Transfer | Micro/nanoscale transport |
| TU Munich | MW 0854 Heat/Mass Transfer | Microstructures |
| ETH | 151-0111 Heat Transfer | Nanoscale transport |
| Tsinghua | 传热学 | Micro/nano channel |
| Cambridge | 4M12 Microfluidics | Lab-on-a-chip transport |

## Build

```bash
make          # Compile all source files
make test     # Run test suite
make examples # Run examples
make clean    # Clean build artifacts
```

## File Structure

```
├── Makefile
├── README.md
├── include/
│   ├── micro_nano_transport.h          # Core types and API
│   ├── knudsen_flow_regimes.h          # Kn regime classification
│   ├── electrokinetic_transport.h      # Electrokinetic API
│   └── nanochannel_ion_transport.h     # Ion transport API
├── src/
│   ├── knudsen_regime.c               # Kn number and flow regimes
│   ├── slip_flow_solver.c             # Pressure-driven slip flow
│   ├── electroosmotic_flow.c          # EOF solutions
│   ├── poisson_boltzmann.c            # PB equation solvers
│   ├── nanochannel_ion_transport.c    # Nanochannel ion transport
│   ├── rarefied_gas_dynamics.c        # Rarefied gas dynamics
│   ├── transport_coefficients.c       # Transport property computation
│   ├── concentration_polarization.c   # CP modeling
│   ├── nernst_planck_solver.c         # Nernst-Planck solver
│   ├── channel_geometry.c             # Geometry and hydraulics
│   ├── streaming_potential.c          # Streaming potential/current
│   ├── entropic_transport.c           # Entropic barriers (L8)
│   ├── membrane_transport.c           # Membrane processes (L7)
│   └── micro_nano_transport.lean      # Lean 4 formalization
├── tests/
│   └── test_core.c                    # Core API test suite
├── examples/
│   ├── example_knudsen_regime.c       # Kn regime demonstration
│   ├── example_slip_flow.c            # Slip flow demonstration
│   └── example_electroosmotic.c       # EOF demonstration
├── docs/
│   ├── knowledge-graph.md            # Nine-layer knowledge map
│   ├── coverage-report.md            # Coverage assessment
│   ├── gap-report.md                 # Missing items
│   ├── course-alignment.md           # Course mapping
│   └── course-tree.md                # Dependency tree
├── demos/
└── benches/
```

## Known Limitations

- Lean 4 formalization uses Rat (rational numbers) and linear approximations; a full real-number Lean formalization would require mathlib4
- DSMC and full Boltzmann solvers are computationally intensive and left for future work
- Quantum-confined transport is documented but not numerically implemented

---

*Reference: Bird, Stewart, Lightfoot (2007) Transport Phenomena, 2nd ed.*
*Reference: Karniadakis, Beskok, Aluru (2005) Microflows and Nanoflows*
*Reference: Schoch, Han, Renaud (2008) Rev. Mod. Phys. 80, 839*

# Knowledge Graph: Transport in Micro/Nano Channels

## L1: Definitions — Complete
- Knudsen number: Kn = λ / L_char (mean free path / characteristic length)
- Mean free path (hard-sphere): λ = 1/(√2·π·d²·n)
- Mean free path (VHS model): λ_VHS = λ_HS · (T/T_ref)^(ω-0.5)
- Knudsen number regime thresholds: continuum (<0.001), slip (0.001-0.1), transition (0.1-10), free molecular (>10)
- Slip length (Maxwell): L_s = ((2-σ)/σ)·λ
- Temperature jump length: ζ_T = ((2-α_T)/α_T)·(2γ/((γ+1)Pr))·λ
- Debye length: κ⁻¹ = sqrt(ε·kB·T/(2·e²·z²·c_bulk·NA))
- Zeta potential: ψ at shear plane, controls electrokinetics
- Hydraulic diameter: D_h = 4A/P_wetted
- Bjerrum length: l_B = e²/(4π·ε·kB·T)
- Gouy-Chapman length: λ_GC = e/(2π·l_B·|σ_s|)

## L2: Core Concepts — Complete
- Flow regime classification (4 regimes + 6 sub-regimes)
- Continuum hypothesis breakdown at microscale
- Velocity slip and temperature jump
- Rarefaction effects on transport
- Electroosmotic flow (EOF)
- Electric double layer (EDL) structure
- Debye-Huckel linearization
- Concentration polarization (CP)
- Entropic barriers in nanofluidics
- Ion selectivity in nanochannels
- Donnan equilibrium and exclusion

## L3: Engineering Quantities — Complete
- Slip coefficients (Maxwell first-order, second-order, Burnett)
- Accommodation coefficients (tangential momentum σ, thermal α_T)
- Hydraulic resistance scaling: R_hyd ~ H⁻³
- Poiseuille number with slip: Po = 96/(1+6L_s/H)
- EOF mobility: μ_eo = -εζ/η [m²/(V·s)]
- Debye length for 1 mM NaCl: ~9.6 nm
- Collision frequency: Z = n·σ·c_bar
- Viscosity from kinetic theory: η ~ √(mT)/d²
- Womersley number: Wo = (H/2)·√(ωρ/η)

## L4: Conservation Laws — Complete
- Navier-Stokes with slip boundary: 0 = -∇p + η∇²u, u|_wall = L_s·∂u/∂n
- Poisson-Boltzmann: ∇²ψ = -ρ_e/ε = (2ze c₀/ε)·sinh(zeψ/kBT)
- Debye-Huckel linearized PB: ∇²ψ = κ²ψ
- Nernst-Planck flux: J_i = -D_i[∇c_i + (z_i e/kBT)·c_i·∇φ]
- Goldman-Hodgkin-Katz flux (constant field)
- Streaming current/potential from nonequilibrium thermodynamics
- Grahame equation: σ_s vs ψ₀
- Donnan equilibrium partitioning
- Boltzmann ion distribution: c_i(y) = c_i_bulk·exp(-z_i e ψ(y)/kBT)

## L5: Engineering Methods — Complete
- Analytical solution: pressure-driven slip flow (parallel plates)
- Analytical solution: second-order slip flow (Deissler)
- Analytical solution: electroosmotic flow (thin EDL + full)
- Poisson-Boltzmann 1D solver (Gouy-Chapman + Debye-Huckel)
- Nernst-Planck steady-state 1D solver (GHK)
- Grahame equation (forward + inverse)
- Transport coefficient computation from kinetic theory
- Sutherland viscosity law, power-law viscosity
- Hydraulic geometry computation (6 cross-section types)
- Mass flow rate with rarefaction (unified model)
- CP layer thickness from film theory

## L6: Engineering Problems — Complete
- Pressure-driven flow in microchannel with slip
- Electroosmotic flow with EDL (thin + overlapping)
- Ionic current through nanochannel
- Nanopore resistance and blockade current
- Limiting current density in electrodialysis
- Knudsen minimum (paradox) location
- Ion current rectification in asymmetric nanopores
- Streaming potential measurement
- Membrane water flux and salt rejection (RO/FO/PRO)

## L7: Applications — Partial+
- Nanopore sensing (DNA sequencing): blockade current model
- Electroosmotic pumping: power consumption, efficiency
- Membrane desalination: RO, FO, PRO models
- Water splitting at overlimiting current
- Streaming potential for zeta measurement
- Nanofluidic ion separation membranes

## L8: Advanced Methods — Partial+
- Burnett equations (second-order slip coefficients)
- BGK relaxation model for Boltzmann equation
- Unified mass flow across all Kn regimes
- Electroconvective instability (Rubinstein-Zaltzman)
- Hindered diffusion in nanopores
- Effective screening length beyond mean-field
- Entropic spring constant for polymers

## L9: Research Frontiers — Partial
- Non-Fourier heat conduction (Cattaneo-Vernotte): structure defined
- Quantum-confined transport in sub-1nm channels: structure defined
- Topological nanofluidic transport: structure defined

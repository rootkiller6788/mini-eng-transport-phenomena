/*
 * dimensionless_numbers.h — Core Dimensionless Numbers in Transport Phenomena
 *
 * Reference: Bird, Stewart, Lightfoot (2007) "Transport Phenomena", 2nd ed.
 *            Incropera & DeWitt (2007) "Fundamentals of Heat and Mass Transfer"
 *            White (2016) "Fluid Mechanics"
 *
 * This header defines all major dimensionless groups used in transport phenomena:
 * fluid mechanics, heat transfer, mass transfer, and their coupling.
 *
 * Theorem (Buckingham Pi, 1914):
 *   Given n physical variables involving k fundamental dimensions,
 *   the system can be described by n-k dimensionless Pi groups.
 *
 * Knowledge Coverage: L1 Definitions, L2 Core Concepts, L3 Engineering Quantities
 */

#ifndef DIMENSIONLESS_NUMBERS_H
#define DIMENSIONLESS_NUMBERS_H

#ifdef __cplusplus
extern "C" {
#endif

/* ==========================================================================
 * L1: Core Dimensionless Number Definitions (Struct-based)
 * ========================================================================== */

/**
 * FluidProperties — material properties needed for dimensionless number computation.
 *
 * Fields:
 *   density         — rho [kg/m^3], fluid mass per unit volume
 *   dynamic_viscosity — mu [Pa.s = kg/(m.s)], resistance to shear flow
 *   kinematic_viscosity — nu = mu/rho [m^2/s], momentum diffusivity
 *   thermal_conductivity — k [W/(m.K)], Fourier's law coefficient
 *   specific_heat_cp  — cp [J/(kg.K)], heat capacity at constant pressure
 *   thermal_diffusivity — alpha = k/(rho.cp) [m^2/s], thermal diffusivity
 *   mass_diffusivity   — D_AB [m^2/s], binary mass diffusivity of A in B
 *   thermal_expansion   — beta [1/K], volumetric thermal expansion coefficient
 *   surface_tension     — sigma [N/m], liquid-gas interface energy per area
 *   speed_of_sound      — c [m/s], acoustic wave speed in the fluid
 */
typedef struct {
    double density;
    double dynamic_viscosity;
    double kinematic_viscosity;
    double thermal_conductivity;
    double specific_heat_cp;
    double thermal_diffusivity;
    double mass_diffusivity;
    double thermal_expansion;
    double surface_tension;
    double speed_of_sound;
} FluidProperties;

/**
 * FlowGeometry — characteristic geometric parameters for a flow configuration.
 *
 * Fields:
 *   characteristic_length — L [m], typically hydraulic diameter for pipes,
 *                            plate length for flat plates, etc.
 *   hydraulic_diameter    — D_h = 4A/P [m], for non-circular ducts
 *   plate_length          — [m], for flat-plate boundary layer problems
 *   pipe_diameter         — [m], for internal pipe flow
 *   cross_sectional_area  — A [m^2]
 *   wetted_perimeter      — P [m]
 *   wall_roughness        — epsilon [m], for Moody diagram applications
 */
typedef struct {
    double characteristic_length;
    double hydraulic_diameter;
    double plate_length;
    double pipe_diameter;
    double cross_sectional_area;
    double wetted_perimeter;
    double wall_roughness;
} FlowGeometry;

/**
 * FlowConditions — operating conditions for the flow.
 *
 * Fields:
 *   velocity         — U [m/s], free-stream or bulk velocity
 *   wall_temperature — T_w [K], temperature at the solid surface
 *   bulk_temperature — T_b [K], bulk fluid temperature (far from wall)
 *   delta_temperature — DeltaT = T_w - T_b [K]
 *   pressure_drop    — DeltaP [Pa]
 *   gravity          — g [m/s^2], gravitational acceleration
 *   heat_flux        — q" [W/m^2], wall heat flux
 *   mass_flux        — m" [kg/(m^2.s)], wall mass flux
 *   solute_conc_wall — C_w [mol/m^3], concentration at wall
 *   solute_conc_bulk — C_b [mol/m^3], bulk concentration
 */
typedef struct {
    double velocity;
    double wall_temperature;
    double bulk_temperature;
    double delta_temperature;
    double pressure_drop;
    double gravity;
    double heat_flux;
    double mass_flux;
    double solute_conc_wall;
    double solute_conc_bulk;
} FlowConditions;

/**
 * DimensionlessNumbers — computed dimensionless groups for a given problem.
 *
 * Each group represents a physically meaningful ratio of forces, rates,
 * or length scales. These groups enable similarity analysis and empirical
 * correlation across different fluids, scales, and operating conditions.
 *
 * Group definitions (L1):
 *   Re  = rhoUL/mu     — Reynolds number: inertial / viscous force ratio
 *   Nu  = hL/k      — Nusselt number: convective / conductive heat transfer ratio
 *   Pr  = nu/alpha       — Prandtl number: momentum / thermal diffusivity ratio
 *   Sc  = nu/D       — Schmidt number: momentum / mass diffusivity ratio
 *   Pe_H = Re.Pr    — Péclet number (heat): advective / diffusive heat transport
 *   Pe_M = Re.Sc    — Péclet number (mass): advective / diffusive mass transport
 *   Gr  = gbetaDeltaT L^3/nu^2 — Grashof number: buoyancy / viscous force ratio
 *   Ra  = Gr.Pr     — Rayleigh number: driving force for natural convection
 *   St_H = Nu/(Re.Pr) — Stanton number (heat): heat transfer / thermal capacity
 *   St_M = Sh/(Re.Sc) — Stanton number (mass): mass transfer / advective capacity
 *   We  = rhoU^2L/sigma    — Weber number: inertia / surface tension ratio
 *   Ma  = U/c       — Mach number: flow velocity / speed of sound
 *   Fr  = U/sqrt(gL) — Froude number: inertia / gravity ratio
 *   Eu  = DeltaP/(rhoU^2)  — Euler number: pressure drop / dynamic pressure
 *   Br  = muU^2/(kDeltaT) — Brinkman number: viscous dissipation / conductive heat
 *   Bi  = hL/k_s    — Biot number: internal conduction / surface convection resistance
 *   Fo  = alphat/L^2     — Fourier number: dimensionless time for heat diffusion
 *   Sh  = h_m L/D   — Sherwood number: convective / diffusive mass transfer
 *   Le  = alpha/D       — Lewis number: thermal / mass diffusivity ratio
 *   Gz  = Re.Pr.D/L — Graetz number: thermal entrance region parameter
 *   Bo  = rhogL^2/sigma    — Bond number: gravity / surface tension
 *   Ca  = muU/sigma      — Capillary number: viscous / surface tension
 *   El  = sigma/(rhogL^2)  — Eötvös number (= 1/Bo), also known as Bond number inverse
 */
typedef struct {
    /* Flow regime indicators */
    double Re;       /* Reynolds number */
    double Ma;       /* Mach number */
    double Fr;       /* Froude number */
    double We;       /* Weber number */
    double Eu;       /* Euler number */

    /* Heat transfer groups */
    double Nu;       /* Nusselt number */
    double Pr;       /* Prandtl number */
    double Pe_H;     /* Péclet number (heat) */
    double St_H;     /* Stanton number (heat) */
    double Br;       /* Brinkman number */
    double Bi;       /* Biot number */
    double Fo;       /* Fourier number */
    double Gz;       /* Graetz number */

    /* Natural convection */
    double Gr;       /* Grashof number */
    double Ra;       /* Rayleigh number */
    double Bo;       /* Bond number */

    /* Mass transfer groups */
    double Sc;       /* Schmidt number */
    double Sh;       /* Sherwood number */
    double Pe_M;     /* Péclet number (mass) */
    double St_M;     /* Stanton number (mass) */
    double Le;       /* Lewis number */

    /* Interfacial groups */
    double Ca;       /* Capillary number */
    double El;       /* Eötvös number */

    /* Derived */
    double Nu_avg;   /* Average Nusselt number (integrated) */
    double Cf;       /* Skin friction coefficient = 2tau_w/(rhoU^2) */
    double C_H;      /* Stanton number alternate = St_H */
    double C_M;      /* Mass transfer Stanton number */
} DimensionlessNumbers;

/* ==========================================================================
 * L2: Flow Regime and Similarity Enums
 * ========================================================================== */

/**
 * FlowRegime — classification of flow behaviour based on Re.
 */
typedef enum {
    FLOW_CREEPING       = 0,  /* Re < 0.1, Stokes flow, viscous-dominated */
    FLOW_LAMINAR        = 1,  /* Re < 2300 (pipes), smooth layered flow */
    FLOW_TRANSITIONAL   = 2,  /* 2300 < Re < 4000, intermittently turbulent */
    FLOW_TURBULENT_SMOOTH = 3, /* Re > 4000, fully turbulent smooth wall */
    FLOW_TURBULENT_ROUGH  = 4, /* Turbulent with significant roughness effects */
    FLOW_HYPERSONIC     = 5,  /* Ma > 5, compressibility dominates */
    FLOW_INVALID        = -1  /* Negative Re or NaN inputs */
} FlowRegime;

/**
 * ConvectionType — driving mechanism for fluid motion.
 */
typedef enum {
    CONV_NONE           = 0,
    CONV_FORCED         = 1,  /* Gr/Re^2 << 1: externally driven flow */
    CONV_NATURAL        = 2,  /* Gr/Re^2 >> 1: buoyancy-driven flow */
    CONV_MIXED          = 3   /* Gr/Re^2 ~ O(1): both mechanisms significant */
} ConvectionType;

/**
 * ThermalBC — thermal boundary condition type for heat transfer analysis.
 */
typedef enum {
    BC_UNIFORM_TEMPERATURE = 0, /* Constant wall temperature (T_w = const) */
    BC_UNIFORM_HEAT_FLUX   = 1, /* Constant wall heat flux (q"_w = const) */
    BC_CONVECTIVE          = 2, /* Robin condition: -k(dT/dn) = h(T_w - T_inf) */
    BC_ADIABATIC           = 3  /* Insulated wall, q"_w = 0 */
} ThermalBC;

/* ==========================================================================
 * L3: Engineering Quantities — Typical Value Ranges
 * ========================================================================== */

/**
 * TypicalPr — empirical Prandtl number ranges for common fluids.
 * Source: Incropera & DeWitt (2007), Appendix A.
 *
 * Values at ~300K, 1 atm unless noted.
 *   LIQUID_METALS     ~0.004-0.03  (sodium, mercury, NaK)
 *   GASES             ~0.7         (air, N2, O2, H2 ~0.7; steam ~1.0)
 *   WATER             ~2-7         (varies strongly with T: 13 at 0 degC, 2 at 100 degC)
 *   LIGHT_OILS        ~10-100      (engine oils at moderate T)
 *   HEAVY_OILS        ~100-10000   (glycerin, heavy lubricants)
 *   POLYMER_MELTS     ~10⁴-10⁶     (very high viscosity)
 */
typedef enum {
    PR_LIQUID_METAL  = 0,
    PR_GAS           = 1,
    PR_WATER         = 2,
    PR_LIGHT_OIL     = 3,
    PR_HEAVY_OIL     = 4,
    PR_POLYMER_MELT  = 5
} PrCategory;

/**
 * CriticalRe — critical Reynolds numbers for common geometries.
 * Source: White (2016), Incropera & DeWitt (2007).
 *
 *   PIPE_FLOW      Recrit ~ 2300 (lower bound; can reach 10⁴ in smooth pipes)
 *   FLAT_PLATE     Recrit ~ 5×10⁵ (transition from laminar to turbulent BL)
 *   SPHERE         Recrit ~ 2.5×10⁵ (drag crisis)
 *   CYLINDER       Recrit ~ 2×10⁵ (drag crisis for circular cylinder)
 *   TAYLOR_COUETTE Ta_crit ~ 1708 (narrow gap, inner cylinder rotating)
 */
#define RE_CRIT_PIPE        2300.0
#define RE_CRIT_FLAT_PLATE  500000.0
#define RE_CRIT_SPHERE      250000.0
#define RE_CRIT_CYLINDER    200000.0
#define RE_CRIT_LOWER_BOUND 1800.0   /* Absolute lower bound for any pipe transition */
#define RE_CRIT_UPPER_BOUND 4000.0   /* Upper transitional bound */

/* ==========================================================================
 * L2: Dimensionless Group Computation API
 * ========================================================================== */

/**
 * compute_reynolds_number — Re = rhoUL/mu = UL/nu
 *
 * @param rho   density [kg/m^3]
 * @param U     velocity [m/s]
 * @param L     characteristic length [m]
 * @param mu    dynamic viscosity [Pa.s]
 * @return      Reynolds number [-], negative on invalid input
 *
 * Theorem: Re is the ratio of inertial to viscous forces.
 *   Re = (rhoU^2/L^2 in characteristic form) -> rhoUL/mu
 *   Re << 1  -> Stokes flow (viscous forces dominate)
 *   Re >> 1  -> Inertia-dominated (boundary layer theory valid)
 */
double compute_reynolds_number(double rho, double U, double L, double mu);

/**
 * compute_nusselt_number — Nu = hL/k
 *
 * @param h    convective heat transfer coefficient [W/(m^2.K)]
 * @param L    characteristic length [m]
 * @param k    fluid thermal conductivity [W/(m.K)]
 * @return     Nusselt number [-]
 *
 * Theorem: Nu represents the ratio of convective to conductive heat transfer.
 *   Nu = 1 -> pure conduction; Nu >> 1 -> strong convection.
 */
double compute_nusselt_number(double h, double L, double k);

/**
 * compute_prandtl_number — Pr = nu/alpha = cp.mu/k
 *
 * @param nu   kinematic viscosity [m^2/s]
 * @param alpha thermal diffusivity [m^2/s]
 * @return     Prandtl number [-]
 *
 * Theorem: Pr relates momentum diffusivity to thermal diffusivity.
 *   Pr << 1 -> thermal boundary layer >> velocity boundary layer (liquid metals)
 *   Pr ~ 1 -> similar thickness (gases)
 *   Pr >> 1 -> velocity BL >> thermal BL (oils)
 *   delta_T/delta_H ≈ Pr^(-1/3) for laminar flow over flat plate
 */
double compute_prandtl_number(double nu, double alpha);

/**
 * compute_schmidt_number — Sc = nu/D_AB
 *
 * @param nu   kinematic viscosity [m^2/s]
 * @param D_AB binary mass diffusivity [m^2/s]
 * @return     Schmidt number [-]
 */
double compute_schmidt_number(double nu, double D_AB);

/**
 * compute_peclet_number — Pe = Re.Pr (heat) or Re.Sc (mass)
 */
double compute_peclet_number_heat(double Re, double Pr);
double compute_peclet_number_mass(double Re, double Sc);

/**
 * compute_grashof_number — Gr = g.beta.DeltaT.L^3 / nu^2
 *
 * @param g     gravitational acceleration [m/s^2]
 * @param beta  volumetric thermal expansion coefficient [1/K]
 * @param delta_T temperature difference |T_w - T_inf| [K]
 * @param L     characteristic length [m]
 * @param nu    kinematic viscosity [m^2/s]
 * @return      Grashof number [-]
 *
 * Theorem: Gr is the ratio of buoyancy to viscous forces.
 *   Gr > 10⁹ -> turbulent natural convection
 */
double compute_grashof_number(double g, double beta, double delta_T,
                              double L, double nu);

/**
 * compute_rayleigh_number — Ra = Gr.Pr
 *
 * @return Rayleigh number [-]
 *
 * Theorem: Ra determines the onset of natural convection.
 *   Ra > Ra_crit (~1708 for horizontal layer heated from below) -> convection begins.
 *   Ra < 10⁸ -> laminar; Ra > 10¹⁰ -> turbulent natural convection.
 */
double compute_rayleigh_number(double Gr, double Pr);

/**
 * compute_stanton_number — St = Nu/(Re.Pr) = h/(rho.cp.U)
 *
 * @return Stanton number [-]
 *
 * Theorem: St is the ratio of heat transferred to thermal capacity of the fluid.
 *   Useful in heat exchanger analysis (Colburn analogy: St.Pr^(2/3) = f/8).
 */
double compute_stanton_number_heat(double Nu, double Re, double Pr);

/**
 * compute_sherwood_number — Sh = h_m.L/D_AB
 *
 * @param h_m  mass transfer coefficient [m/s]
 * @param L    characteristic length [m]
 * @param D_AB binary mass diffusivity [m^2/s]
 * @return     Sherwood number [-]
 *
 * Theorem: Sh is the mass-transfer analog of Nu.
 *   Chilton-Colburn analogy: Nu/Sh = (Pr/Sc)^(1/3) for similar geometry.
 */
double compute_sherwood_number(double h_m, double L, double D_AB);

/**
 * compute_lewis_number — Le = alpha/D_AB = Sc/Pr
 *
 * @return Lewis number [-]
 *
 * Theorem: Le compares thermal to mass diffusion rates.
 *   Le ~ 1 -> heat and mass transfer at similar rates (most gases)
 *   Le >> 1 -> heat diffuses faster than mass (liquids)
 *   Le << 1 -> mass diffuses faster than heat (rare)
 */
double compute_lewis_number(double alpha, double D_AB);

/**
 * compute_weber_number — We = rhoU^2L/sigma
 */
double compute_weber_number(double rho, double U, double L, double sigma);

/**
 * compute_mach_number — Ma = U/c
 */
double compute_mach_number(double U, double c);

/**
 * compute_froude_number — Fr = U/sqrt(gL)
 */
double compute_froude_number(double U, double g, double L);

/**
 * compute_euler_number — Eu = DeltaP/(rhoU^2)
 *   Eu is twice the friction factor for internal flows: Eu = f.(L/D)
 */
double compute_euler_number(double delta_p, double rho, double U);

/**
 * compute_brinkman_number — Br = muU^2/(k.DeltaT)
 *
 * Theorem: Br represents the relative importance of viscous dissipation
 *   to conductive heat transfer. Negligible for low-speed flows.
 *   Significant in high-speed lubrication and polymer processing.
 */
double compute_brinkman_number(double mu, double U, double k, double delta_T);

/**
 * compute_biot_number — Bi = h.L_c/k_s
 *
 * @param h   convective heat transfer coefficient [W/(m^2.K)]
 * @param L_c characteristic conduction length (V/A_s) [m]
 * @param k_s solid thermal conductivity [W/(m.K)]
 * @return    Biot number [-]
 *
 * Theorem: Bi is the ratio of internal conduction resistance to
 *   surface convection resistance.
 *   Bi < 0.1 -> lumped capacitance method valid (T uniform in solid)
 *   Bi > 0.1 -> spatial temperature distribution must be resolved
 */
double compute_biot_number(double h, double L_c, double k_s);

/**
 * compute_fourier_number — Fo = alpha.t/L^2
 *
 * @param alpha thermal diffusivity [m^2/s]
 * @param t     time [s]
 * @param L     characteristic length [m]
 * @return      Fourier number [-]
 *
 * Theorem: Fo is dimensionless time for transient conduction.
 *   Fo > 0.2 -> first-term approximation valid for series solutions
 */
double compute_fourier_number(double alpha, double t, double L);

/**
 * compute_graetz_number — Gz = Re.Pr.D/L = (D/L).Pe_H
 *
 * @param D   pipe diameter [m]
 * @param L   pipe length [m]
 * @param Pe  Péclet number [-]
 * @return    Graetz number [-]
 *
 * Theorem: Gz characterizes thermal entrance region in ducts.
 *   Gz⁻¹ = x* = (x/D)/(Re.Pr) is the dimensionless axial coordinate.
 *   Gz > 1000 -> thermally developing; Gz < 10 -> thermally fully developed.
 */
double compute_graetz_number(double D, double L, double Pe);

/**
 * compute_capillary_number — Ca = muU/sigma
 */
double compute_capillary_number(double mu, double U, double sigma);

/**
 * compute_bond_number — Bo = rhogL^2/sigma
 */
double compute_bond_number(double rho, double g, double L, double sigma);

/**
 * compute_eotvos_number — Eo = rhogL^2/sigma (= Bo, but more commonly used in drops/bubbles)
 */
double compute_eotvos_number(double rho, double g, double L, double sigma);

/* ==========================================================================
 * L2: Composite Computation Functions
 * ========================================================================== */

/**
 * compute_all_dimensionless — fill complete DimensionlessNumbers struct.
 *
 * @param props   fluid material properties
 * @param geom    flow geometry parameters
 * @param cond    operating conditions
 * @param[out] dn filled DimensionlessNumbers struct
 * @return        0 on success, -1 on invalid inputs
 *
 * Uses all the individual computation functions above.  Additionally computes:
 *   - skin friction coefficient Cf = 2.tau_w/(rhoU^2)
 *   - Nu_avg correlations (laminar & turbulent flat plate)
 *   - Convection type classification based on Gr/Re^2 ratio
 *
 * Algorithm complexity: O(1), ~30 double operations
 */
int compute_all_dimensionless(const FluidProperties *props,
                              const FlowGeometry *geom,
                              const FlowConditions *cond,
                              DimensionlessNumbers *dn);

/**
 * classify_flow_regime — determine flow type from Reynolds number and geometry.
 *
 * @param Re   Reynolds number
 * @param geom flow geometry (for critical Re selection)
 * @return     FlowRegime enum value
 */
FlowRegime classify_flow_regime(double Re, const FlowGeometry *geom);

/**
 * classify_convection_type — determine if natural, forced, or mixed convection.
 *
 * @param Gr  Grashof number
 * @param Re  Reynolds number
 * @return    ConvectionType enum value
 *
 * Criterion:
 *   Gr/Re^2 < 0.01 -> forced convection dominates
 *   Gr/Re^2 > 10   -> natural convection dominates
 *   0.01 < Gr/Re^2 < 10 -> mixed convection
 */
ConvectionType classify_convection_type(double Gr, double Re);

/**
 * prandtl_category — classify Pr into typical engineering categories.
 *
 * @param Pr  Prandtl number [-]
 * @return    PrCategory enum value
 */
PrCategory prandtl_category(double Pr);

/* ==========================================================================
 * L2: Similarity and Scale Analysis
 * ========================================================================== */

/**
 * check_dynamic_similarity — verify that two flows satisfy Reynolds similarity.
 *
 * @param Re1  Reynolds number of flow 1 (prototype)
 * @param Re2  Reynolds number of flow 2 (model)
 * @param tol  tolerance for equality check (e.g., 0.01 for 1%)
 * @return     1 if dynamically similar, 0 otherwise
 *
 * Theorem (Reynolds Similarity): Two geometrically similar flows with equal Re
 *   have dynamically similar velocity fields and identical dimensionless
 *   pressure distributions.
 */
int check_dynamic_similarity(double Re1, double Re2, double tol);

/**
 * check_thermal_similarity — verify similarity of heat transfer conditions.
 *
 * @param Pr1  Prandtl number of flow 1
 * @param Pr2  Prandtl number of flow 2
 * @param tol  tolerance
 * @return     1 if thermally similar, 0 otherwise
 */
int check_thermal_similarity(double Pr1, double Pr2, double tol);

/**
 * scale_model_velocity — compute model velocity to achieve Re similarity.
 *
 * Given a prototype with (L1, U1, nu1) and model scale L2, nu2,
 * find U2 such that Re1 = Re2.
 *
 * @param L_proto  prototype length [m]
 * @param U_proto  prototype velocity [m/s]
 * @param nu_proto prototype kinematic viscosity [m^2/s]
 * @param L_model  model length [m]
 * @param nu_model model kinematic viscosity [m^2/s]
 * @return         model velocity [m/s], -1 on invalid inputs
 */
double scale_model_velocity(double L_proto, double U_proto, double nu_proto,
                            double L_model, double nu_model);

/**
 * boundary_layer_thickness_ratio — estimate delta_T/delta from Pr or delta_C/delta from Sc.
 *
 * @param Pr_or_Sc  Prandtl number (for thermal) or Schmidt number (for mass)
 * @return          delta_T/delta ≈ Pr^(-1/3) (laminar) or delta_C/delta ≈ Sc^(-1/3)
 *
 * Theorem (Laminar BL scaling, e.g., Schlichting):
 *   For Pr ≥ 0.6 and laminar flow:  delta_T/delta ≈ Pr^(-1/3)
 *   General correlation: delta_T/delta = Pr^(-1/3) for Pr > 0.5
 *   For Pr -> 0 (liquid metals): delta_T/delta ≈ 1.325.Pr^(-1/2)
 */
double boundary_layer_thickness_ratio(double Pr_or_Sc);

#ifdef __cplusplus
}
#endif

#endif /* DIMENSIONLESS_NUMBERS_H */

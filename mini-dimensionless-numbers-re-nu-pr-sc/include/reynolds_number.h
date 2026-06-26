/*
 * reynolds_number.h — Reynolds Number Deep Dive
 *
 * Reference: White (2016) "Fluid Mechanics", Ch. 5-6
 *            Bird, Stewart, Lightfoot (2007) "Transport Phenomena", Ch. 6
 *            Schlichting & Gersten (2017) "Boundary-Layer Theory"
 *
 * The Reynolds number Re = rhoUL/mu is the most important dimensionless group
 * in fluid mechanics. It characterises the ratio of inertial to viscous forces
 * and determines whether flow is laminar or turbulent.
 *
 * Knowledge Coverage: L1 Definitions, L2 Flow Regime Theory, L3 Critical Re,
 *                     L4 Momentum Conservation (Navier-Stokes scaling),
 *                     L5 Engineering Methods (Re-based friction factor)
 */

#ifndef REYNOLDS_NUMBER_H
#define REYNOLDS_NUMBER_H

#include <stddef.h>  /* size_t */

#ifdef __cplusplus
extern "C" {
#endif

/* ==========================================================================
 * L1: Re Variants for Different Geometries
 * ========================================================================== */

/**
 * Pipe Reynolds Number: Re_D = 4.ṁ/(π.D.mu)
 *
 * @param m_dot  mass flow rate [kg/s]
 * @param D      pipe diameter [m]
 * @param mu     dynamic viscosity [Pa.s]
 * @return       pipe Reynolds number [-]
 *
 * Theorem: For fully developed internal flow in a circular pipe,
 *   Re_D < 2300 -> Laminar
 *   2300 < Re_D < 4000 -> Transitional
 *   Re_D > 4000 -> Turbulent
 *
 * Alternative form: Re_D = rho.U_avg.D/mu where U_avg = 4ṁ/(rhoπD^2)
 */
double re_pipe_from_mass_flow(double m_dot, double D, double mu);

/**
 * Flat Plate Reynolds Number: Re_x = rhoUx/mu
 *
 * @param rho  density [kg/m^3]
 * @param U    free-stream velocity [m/s]
 * @param x    distance from leading edge [m]
 * @param mu   dynamic viscosity [Pa.s]
 * @return     local Re_x [-]
 *
 * Theorem: Re_x determines boundary layer state.
 *   Re_x < 5×10⁵ -> Laminar BL
 *   Re_x > 5×10⁵ -> Transition to turbulent (may vary with free-stream turbulence)
 */
double re_flat_plate_local(double rho, double U, double x, double mu);

/**
 * Sphere/Cylinder Reynolds Number: Re = rhoUD/mu
 *
 * @param rho  density [kg/m^3]
 * @param U    free-stream velocity [m/s]
 * @param D    sphere/cylinder diameter [m]
 * @param mu   dynamic viscosity [Pa.s]
 * @return     Re [-]
 */
double re_sphere(double rho, double U, double D, double mu);
double re_cylinder(double rho, double U, double D, double mu);

/**
 * Hydraulic Diameter Reynolds Number: Re_Dh = rhoU.D_h/mu
 *
 * @param rho  density [kg/m^3]
 * @param U    bulk velocity [m/s]
 * @param D_h  hydraulic diameter = 4A/P [m]
 * @param mu   dynamic viscosity [Pa.s]
 * @return     Re based on hydraulic diameter [-]
 *
 * Theorem: The hydraulic diameter concept extends pipe flow results
 *   to non-circular ducts. D_h = 4.(cross-sectional area)/(wetted perimeter).
 *   For circular pipe: D_h = D.
 *   For rectangular duct a×b: D_h = 2ab/(a+b).
 */
double re_hydraulic_diameter(double rho, double U, double D_h, double mu);

/**
 * Rotational Reynolds Number: Re_Omega = rhoOmegaD^2/mu
 *
 * @param rho    density [kg/m^3]
 * @param omega  angular velocity [rad/s]
 * @param D      impeller/rotor diameter [m]
 * @param mu     dynamic viscosity [Pa.s]
 * @return       rotational Re [-]
 *
 * Used in stirred tanks, rotating machinery, Taylor-Couette flow.
 */
double re_rotational(double rho, double omega, double D, double mu);

/* ==========================================================================
 * L2: Flow Regime Determination and Analysis
 * ========================================================================== */

/**
 * FlowRegimeExtended — detailed pipe flow regime classification.
 */
typedef enum {
    RE_PIPE_LAMINAR            = 0,  /* Re < 2300 */
    RE_PIPE_LOWER_TRANSITIONAL = 1,  /* 2000 < Re < 2700 */
    RE_PIPE_UPPER_TRANSITIONAL = 2,  /* 2700 < Re < 4000 */
    RE_PIPE_TURBULENT          = 3,  /* 4000 < Re < 10⁵ */
    RE_PIPE_FULLY_ROUGH        = 4,  /* Turbulent + roughness dominates */
    RE_PIPE_INVALID            = -1
} FlowRegimePipe;

/**
 * flow_regime_pipe — detailed pipe flow regime classification.
 *
 * @param Re            Reynolds number [-]
 * @param roughness_e   wall roughness height [m]
 * @param D             pipe diameter [m]
 * @return              FlowRegimePipe enum
 *
 * Accounts for: Re magnitude, relative roughness epsilon/D.
 * epsilon/D < 10⁻⁶ -> hydraulically smooth.
 */
FlowRegimePipe flow_regime_pipe(double Re, double roughness_e, double D);

/**
 * flow_regime_external — external flow regime (flat plate, sphere, cylinder).
 *
 * @param Re        Reynolds number [-]
 * @param recrit    critical Re for this geometry
 * @return          0 = laminar, 1 = transitional, 2 = turbulent, -1 = invalid
 */
int flow_regime_external(double Re, double recrit);

/**
 * is_laminar — convenience check: returns 1 if flow is laminar.
 */
int is_laminar(double Re, double Re_crit);

/**
 * is_turbulent — convenience check: returns 1 if fully turbulent.
 */
int is_turbulent(double Re, double Re_crit);

/* ==========================================================================
 * L3: Critical Reynolds Number Database
 * ========================================================================== */

/**
 * critical_reynolds — returns the critical Re for a given configuration.
 *
 * Configurations:
 *   0 = circular pipe, smooth
 *   1 = flat plate (zero pressure gradient)
 *   2 = sphere
 *   3 = circular cylinder (cross-flow)
 *   4 = triangular duct
 *   5 = square duct
 *   6 = concentric annulus
 *   7 = open channel (hydraulic)
 *   8 = packed bed (Ergun range: Re_p ~ 1-10)
 *   9 = Taylor-Couette (inner cylinder rotating)
 *
 * @param config_id  configuration identifier (0-9)
 * @return            critical Re, -1 for invalid config_id
 */
double critical_reynolds(int config_id);

/**
 * critical_re_description — human-readable description of the configuration.
 *
 * @param config_id configuration identifier (0-9)
 * @param buf       output buffer (caller-allocated)
 * @param bufsz     buffer size
 * @return          0 on success, -1 on invalid config_id
 */
int critical_re_description(int config_id, char *buf, size_t bufsz);

/* ==========================================================================
 * L4: Non-dimensional Navier-Stokes Equation
 * ========================================================================== */

/**
 * ns_dimensionless_scale — return the scale factors for N-S non-dimensionalization.
 *
 * Given characteristic (rho0, U0, L0, mu0), the dimensionless N-S equation is:
 *   Re (du* / dt* + u* grad* u_star) = -grad* p* + laplacian* u* + ...
 * where * denotes dimensionless variables.
 *
 * This function computes Re from the scale parameters.
 *
 * @param rho0  reference density [kg/m^3]
 * @param U0    reference velocity [m/s]
 * @param L0    reference length [m]
 * @param mu0   reference viscosity [Pa.s]
 * @return      Re scale, -1 on invalid inputs
 */
double ns_re_scale(double rho0, double U0, double L0, double mu0);

/**
 * ns_pressure_scale — compute the pressure scale factor.
 *
 * In non-dimensional form: p* = p/(rhoU^2) or p* = p.L/(muU)
 *   depending on the regime.
 *
 * @param rho   density [kg/m^3]
 * @param U     velocity [m/s]
 * @return      pressure scale rhoU^2 [Pa]
 *
 * This is the dynamic pressure — fundamental to the Euler number Eu = Deltap/(rhoU^2).
 */
double ns_pressure_scale(double rho, double U);

/**
 * ns_time_scale — compute convection time scale: t_conv = L/U.
 *
 * @param L   characteristic length [m]
 * @param U   characteristic velocity [m/s]
 * @return    convective time scale [s]
 */
double ns_time_scale(double L, double U);

/**
 * ns_diffusion_time_scale — compute diffusion time scale: t_diff = L^2/nu.
 *
 * @param L   characteristic length [m]
 * @param nu  kinematic viscosity [m^2/s]
 * @return    diffusion time scale [s]
 */
double ns_diffusion_time_scale(double L, double nu);

/* ==========================================================================
 * L5: Friction Factor Methods based on Re
 * ========================================================================== */

/**
 * friction_factor_laminar — Darcy friction factor for laminar pipe flow.
 *
 * @param Re  Reynolds number (Re < 2300 expected)
 * @return    f_D = 64/Re (exact solution of Hagen-Poiseuille flow)
 *
 * Theorem: For fully developed laminar flow in a circular pipe, the
 *   Darcy friction factor is exactly f = 64/Re. This follows from
 *   the parabolic velocity profile u(r) = 2U_avg[1-(r/R)^2].
 *   DeltaP = f.(L/D).(rhoU^2/2)
 */
double friction_factor_laminar(double Re);

/**
 * friction_factor_blasius — turbulent friction factor (Blasius, 1913).
 *
 * @param Re  Reynolds number (4000 < Re < 10⁵ expected)
 * @return    f = 0.316.Re^(-1/4)
 *
 * Range: 4000 < Re < 10⁵, smooth pipes.
 * Error: ±10% for 4000 < Re < 10⁵.
 */
double friction_factor_blasius(double Re);

/**
 * friction_factor_colebrook — implicit Colebrook-White equation.
 *
 * 1/sqrt(f) = -2.0.log10[(epsilon/D)/3.7 + 2.51/(Re.sqrt(f))]
 *
 * @param Re           Reynolds number (Re > 4000)
 * @param rel_roughness relative roughness epsilon/D [-]
 * @param max_iter     maximum iterations for convergence (e.g., 50)
 * @param tol          convergence tolerance (e.g., 1e-10)
 * @return             Darcy friction factor f, -1 on failure
 *
 * Valid for the full turbulent range. Converges rapidly (usually < 10 iterations).
 * Uses the fixed-point iteration method with f = 0.02 as initial guess.
 *
 * Reference: Colebrook & White (1937), Moody (1944).
 */
double friction_factor_colebrook(double Re, double rel_roughness,
                                 int max_iter, double tol);

/**
 * friction_factor_haaland — explicit approximation (Haaland, 1983).
 *
 * 1/sqrt(f) ≈ -1.8.log10[(epsilon/(3.7D))^1.11 + 6.9/Re]
 *
 * @param Re           Reynolds number
 * @param rel_roughness relative roughness epsilon/D
 * @return             Darcy friction factor f
 *
 * Error: ±1.5% compared to Colebrook for 4000 < Re < 10⁸, epsilon/D < 0.05.
 * Most commonly used explicit formula in engineering practice.
 */
double friction_factor_haaland(double Re, double rel_roughness);

/**
 * friction_factor_swamee_jain — another explicit approximation (1976).
 *
 * f = 0.25/[log10(epsilon/(3.7D) + 5.74/Re^0.9)]^2
 *
 * @param Re           Reynolds number (5000 < Re < 10⁸)
 * @param rel_roughness relative roughness (10⁻⁶ < epsilon/D < 0.05)
 * @return             Darcy friction factor f
 *
 * Error: ±1% of Colebrook.
 */
double friction_factor_swamee_jain(double Re, double rel_roughness);

/**
 * friction_factor_churchill — Churchill equation (1977), covers all regimes.
 *
 * Covers laminar, transitional, and turbulent flows in a single,
 * continuous equation. Valid for all Re and epsilon/D.
 *
 * f = 8.[(8/Re)^12 + 1/(A+B)^1.5]^(1/12)
 * where A = [2.457.ln(1/((7/Re)^0.9 + 0.27.epsilon/D))]^16
 *       B = (37530/Re)^16
 *
 * @param Re           Reynolds number
 * @param rel_roughness relative roughness epsilon/D
 * @return             Darcy friction factor f
 *
 * This is the most comprehensive single-equation friction factor model.
 * Reference: Churchill (1977), Chemical Engineering.
 */
double friction_factor_churchill(double Re, double rel_roughness);

/**
 * form_drag_coefficient — Cd for sphere/cylinder as function of Re.
 *
 * @param Re       Reynolds number based on diameter
 * @param geometry 0 = sphere, 1 = cylinder (cross-flow)
 * @return         drag coefficient Cd [-]
 *
 * Empirical fits to experimental data (Clift-Gauvin for spheres).
 *   Re < 1 (Stokes):  Cd = 24/Re
 *   1 < Re < 1000:    Cd ≈ 24/Re.(1 + 0.15.Re^0.687)
 *   10^3 < Re < 2.5×10⁵: Cd ≈ 0.44 (Newton's law)
 *   Re > 2.5×10⁵:     Cd ≈ 0.2 (drag crisis)
 */
double form_drag_coefficient(double Re, int geometry);

/* ==========================================================================
 * L5: Boundary Layer Thickness Estimation
 * ========================================================================== */

/**
 * bl_thickness_laminar — 99% boundary layer thickness (Blasius solution).
 *
 * @param x   distance from leading edge [m]
 * @param Re_x local Reynolds number at x
 * @return    delta_99 ≈ 5.0.x/sqrt(Re_x)
 */
double bl_thickness_laminar(double x, double Re_x);

/**
 * bl_thickness_turbulent — 99% BL thickness (1/7th power law).
 *
 * @param x   distance from leading edge [m]
 * @param Re_x local Reynolds number at x
 * @return    delta_99 ≈ 0.37.x/Re_x^(1/5)
 *
 * Valid for 5×10⁵ < Re_x < 10⁷.
 */
double bl_thickness_turbulent(double x, double Re_x);

/**
 * bl_displacement_thickness_laminar — delta* for laminar BL.
 *
 * @param x   distance from leading edge [m]
 * @param Re_x local Reynolds number
 * @return    delta* ≈ 1.721.x/sqrt(Re_x)
 */
double bl_displacement_thickness_laminar(double x, double Re_x);

/**
 * bl_momentum_thickness_laminar — theta for laminar BL.
 *
 * @param x   distance from leading edge [m]
 * @param Re_x local Reynolds number
 * @return    theta ≈ 0.664.x/sqrt(Re_x)
 */
double bl_momentum_thickness_laminar(double x, double Re_x);

/**
 * bl_transition_location — estimate x_crit where transition occurs.
 *
 * @param Re_crit critical Re for transition (typically 5×10⁵)
 * @param nu      kinematic viscosity [m^2/s]
 * @param U       free-stream velocity [m/s]
 * @return        x_crit [m]
 */
double bl_transition_location(double Re_crit, double nu, double U);

#ifdef __cplusplus
}
#endif

#endif /* REYNOLDS_NUMBER_H */

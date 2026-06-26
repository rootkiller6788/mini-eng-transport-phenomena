/*
 * reynolds_number.c — Reynolds Number: Flow Regimes, Friction, Boundary Layers
 *
 * Reference: White (2016), Bird-Stewart-Lightfoot (2007),
 *            Schlichting & Gersten (2017), Moody (1944)
 *
 * Knowledge: L1-L5 complete
 */

#include "../include/reynolds_number.h"
#include "../include/dimensionless_numbers.h"
#include <math.h>
#include <string.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ==========================================================================
 * L1: Re variants for different geometries
 * ========================================================================== */

double re_pipe_from_mass_flow(double m_dot, double D, double mu)
{
    /*
     * Re_D = 4.ṁ / (π.D.mu)
     *
     * Derivation:
     *   Re_D = rho.U_avg.D/mu, U_avg = ṁ/(rho.A) = ṁ/(rho.π.D^2/4)
     *   -> Re_D = rho.(4ṁ/(rhoπD^2)).D/mu = 4ṁ/(π.D.mu)
     *
     * This form is convenient when mass flow rate is known
     * but density is not directly available.
     */
    if (mu <= 0.0 || D <= 0.0 || m_dot < 0.0) return -1.0;
    if (mu < 1e-20 || D < 1e-15) return -1.0;
    if (m_dot == 0.0) return 0.0;
    return 4.0 * m_dot / (M_PI * D * mu);
}

double re_flat_plate_local(double rho, double U, double x, double mu)
{
    /* Re_x = rho.U.x/mu */
    if (mu <= 0.0 || rho < 0.0 || x < 0.0) return -1.0;
    if (mu < 1e-20) return -1.0;
    return rho * U * x / mu;
}

double re_sphere(double rho, double U, double D, double mu)
{
    if (mu <= 0.0 || rho < 0.0 || D < 0.0) return -1.0;
    if (mu < 1e-20) return -1.0;
    return rho * U * D / mu;
}

double re_cylinder(double rho, double U, double D, double mu)
{
    return re_sphere(rho, U, D, mu); /* same formula; different critical Re */
}

double re_hydraulic_diameter(double rho, double U, double D_h, double mu)
{
    /*
     * Hydraulic diameter: D_h = 4A/P
     * For non-circular ducts, D_h extends pipe flow results.
     * Examples:
     *   Square duct a×a:    D_h = a
     *   Rectangular a×b:    D_h = 2ab/(a+b)
     *   Concentric annulus: D_h = D_o - D_i
     *   Open channel:       D_h = 4.(flow area)/(wetted perimeter)
     */
    if (mu <= 0.0 || rho < 0.0 || D_h < 0.0) return -1.0;
    if (mu < 1e-20) return -1.0;
    return rho * U * D_h / mu;
}

double re_rotational(double rho, double omega, double D, double mu)
{
    /*
     * Re_Omega = rho.Omega.D^2/mu
     *
     * Used in stirred tanks, turbomachinery, Taylor-Couette flow.
     * The characteristic velocity is Omega.D (tip speed).
     */
    if (mu <= 0.0 || rho < 0.0 || D < 0.0 || omega < 0.0) return -1.0;
    if (mu < 1e-20) return -1.0;
    return rho * omega * D * D / mu;
}

/* ==========================================================================
 * L2: Flow regime determination
 * ========================================================================== */

FlowRegimePipe flow_regime_pipe(double Re, double roughness_e, double D)
{
    if (Re < 0.0 || roughness_e < 0.0 || D <= 0.0) return RE_PIPE_INVALID;
    if (Re < 1e-10) return RE_PIPE_LAMINAR; /* stationary fluid */

    double rel_roughness = roughness_e / D;

    if (Re < 2000.0) {
        return RE_PIPE_LAMINAR;
    } else if (Re < 2700.0) {
        return RE_PIPE_LOWER_TRANSITIONAL;
    } else if (Re < 4000.0) {
        return RE_PIPE_UPPER_TRANSITIONAL;
    } else {
        /* Turbulent regime: check if fully rough */
        /* Fully rough when Re > 1000/(epsilon/D), approximately */
        if (rel_roughness > 0.0 && Re * rel_roughness > 70.0) {
            return RE_PIPE_FULLY_ROUGH;
        }
        return RE_PIPE_TURBULENT;
    }
}

int flow_regime_external(double Re, double recrit)
{
    if (Re < 0.0 || recrit <= 0.0) return -1;
    if (Re < recrit * 0.9) return 0;         /* laminar */
    if (Re < recrit * 1.1) return 1;         /* transitional */
    return 2;                                 /* turbulent */
}

int is_laminar(double Re, double Re_crit)
{
    if (Re < 0.0 || Re_crit <= 0.0) return -1;
    return (Re < Re_crit) ? 1 : 0;
}

int is_turbulent(double Re, double Re_crit)
{
    if (Re < 0.0 || Re_crit <= 0.0) return -1;
    return (Re > Re_crit) ? 1 : 0;
}

/* ==========================================================================
 * L3: Critical Reynolds number database
 * ========================================================================== */

double critical_reynolds(int config_id)
{
    switch (config_id) {
        case 0: return RE_CRIT_PIPE;             /* circular pipe: 2300 */
        case 1: return RE_CRIT_FLAT_PLATE;        /* flat plate: 5×10⁵ */
        case 2: return RE_CRIT_SPHERE;            /* sphere: 2.5×10⁵ */
        case 3: return RE_CRIT_CYLINDER;          /* cylinder: 2×10⁵ */
        case 4: return 2300.0;                    /* triangular duct ~ pipe */
        case 5: return 2300.0;                    /* square duct ~ pipe */
        case 6: return 2300.0;                    /* annulus ~ pipe */
        case 7: return 500.0;                     /* open channel: Re_h ~ 500 */
        case 8: return 10.0;                      /* packed bed: Re_p ~ 1-10 */
        case 9: return 1708.0;                    /* Taylor-Couette: Ta_crit ~ 1708 */
        default: return -1.0;
    }
}

int critical_re_description(int config_id, char *buf, size_t bufsz)
{
    if (!buf || bufsz == 0) return -1;

    const char *desc;
    switch (config_id) {
        case 0: desc = "Circular Pipe (smooth)"; break;
        case 1: desc = "Flat Plate (zero pressure gradient)"; break;
        case 2: desc = "Sphere (drag crisis)"; break;
        case 3: desc = "Circular Cylinder (cross-flow)"; break;
        case 4: desc = "Triangular Duct"; break;
        case 5: desc = "Square Duct"; break;
        case 6: desc = "Concentric Annulus"; break;
        case 7: desc = "Open Channel (hydraulic)"; break;
        case 8: desc = "Packed Bed (Ergun regime)"; break;
        case 9: desc = "Taylor-Couette (inner rotating)"; break;
        default: return -1;
    }

    size_t len = strlen(desc);
    if (len >= bufsz) len = bufsz - 1;
    memcpy(buf, desc, len);
    buf[len] = '\0';
    return 0;
}

/* ==========================================================================
 * L4: Non-dimensional Navier-Stokes
 * ========================================================================== */

double ns_re_scale(double rho0, double U0, double L0, double mu0)
{
    /*
     * From N-S: rho(du/dt + u.gradu) = -gradp + mugrad^2u + f_body
     *
     * Non-dimensional variables: u* = u/U0, x* = x/L0, t* = t.U0/L0
     *   p* = p/(rhoU0^2)   (pressure scaling)
     *
     * -> du* / dt* + u_star.grad* u_star = -grad* p_star + (1/Re).grad*^2 u_star + ...
     * where Re = rho0.U0.L0/mu0
     *
     * When Re -> 0: viscous term dominates (Stokes flow)
     * When Re -> inf: inertial terms dominate (Euler equations)
     */
    if (mu0 <= 0.0 || rho0 <= 0.0 || L0 <= 0.0) return -1.0;
    return rho0 * U0 * L0 / mu0;
}

double ns_pressure_scale(double rho, double U)
{
    /* Dynamic pressure: q = ½rhoU^2 -> pressure scale is rhoU^2 */
    if (rho <= 0.0) return -1.0;
    return rho * U * U;
}

double ns_time_scale(double L, double U)
{
    /* Convective time scale: t_conv = L/U */
    if (U <= 0.0 || L <= 0.0) return -1.0;
    return L / U;
}

double ns_diffusion_time_scale(double L, double nu)
{
    /* Diffusive time scale: t_diff = L^2/nu */
    if (nu <= 0.0 || L <= 0.0) return -1.0;
    return L * L / nu;
}

/* ==========================================================================
 * L5: Friction factor methods based on Re
 * ========================================================================== */

double friction_factor_laminar(double Re)
{
    /*
     * f = 64/Re  (Darcy friction factor for circular pipe)
     *
     * Derivation: Hagen-Poiseuille flow
     *   u(r) = (DeltaP/(4muL)).(R^2 - r^2)
     *   U_avg = DeltaP.R^2/(8muL)
     *   DeltaP = (8muLU_avg)/R^2 = (32muLU_avg)/D^2
     *   f = DeltaP/((L/D).(rhoU^2/2)) = 64.mu/(rho.U.D) = 64/Re
     *
     * This is an exact solution of the Navier-Stokes equations.
     */
    if (Re <= 0.0) return -1.0;
    return 64.0 / Re;
}

double friction_factor_blasius(double Re)
{
    /*
     * Blasius (1913) correlation for smooth pipes, turbulent:
     *   f = 0.316.Re^(-1/4)
     *
     * Valid for: 4000 < Re < 10⁵
     * Power-law fit to experimental data (Blasius, 1913).
     */
    if (Re < 4000.0 || Re > 100000.0) return -1.0;
    return 0.316 * pow(Re, -0.25);
}

double friction_factor_colebrook(double Re, double rel_roughness,
                                 int max_iter, double tol)
{
    /*
     * Colebrook-White equation (1937):
     *   1/√f = -2.0.log10[(epsilon/D)/3.7 + 2.51/(Re.√f)]
     *
     * This is an implicit equation. We solve by fixed-point iteration:
     *   f_{n+1} = 1 / [(-2.0.log10(epsilon/(3.7D) + 2.51/(Re.√f_n)))]^2
     *
     * Initial guess: f_0 = 0.02
     *
     * The equation combines the smooth pipe law (Prandtl, 1935)
     * and the fully rough law (von Kármán, 1930).
     */
    if (Re <= 0.0 || rel_roughness < 0.0 || max_iter < 1 || tol <= 0.0)
        return -1.0;

    if (Re < 2300.0) {
        return friction_factor_laminar(Re);
    }

    /* Initial guess from Blasius (smooth) or a reasonable value */
    double f = (Re < 100000.0) ? 0.316 * pow(Re, -0.25) : 0.02;
    double f_old;

    for (int i = 0; i < max_iter; i++) {
        f_old = f;
        double sqrt_f_inv = -2.0 * log10(rel_roughness / 3.7
                                         + 2.51 / (Re * sqrt(f)));
        f = 1.0 / (sqrt_f_inv * sqrt_f_inv);

        if (fabs(f - f_old) < tol) break;
    }

    return f;
}

double friction_factor_haaland(double Re, double rel_roughness)
{
    /*
     * Haaland (1983) explicit approximation:
     *   1/√f ≈ -1.8.log10[(epsilon/(3.7D))^1.11 + 6.9/Re]
     *
     * Error: ±1.5% of Colebrook for 4000 < Re < 10⁸.
     */
    if (Re <= 0.0 || rel_roughness < 0.0) return -1.0;

    if (Re < 2300.0) return friction_factor_laminar(Re);

    double term = pow(rel_roughness / 3.7, 1.11) + 6.9 / Re;
    double sqrt_f_inv = -1.8 * log10(term);
    return 1.0 / (sqrt_f_inv * sqrt_f_inv);
}

double friction_factor_swamee_jain(double Re, double rel_roughness)
{
    /*
     * Swamee-Jain (1976):
     *   f = 0.25 / [log10(epsilon/(3.7D) + 5.74/Re^0.9)]^2
     *
     * Valid for: 5000 < Re < 10⁸, 10⁻⁶ < epsilon/D < 0.05.
     */
    if (Re <= 0.0 || rel_roughness < 0.0) return -1.0;

    if (Re < 2300.0) return friction_factor_laminar(Re);

    double denom = log10(rel_roughness / 3.7 + 5.74 / pow(Re, 0.9));
    return 0.25 / (denom * denom);
}

double friction_factor_churchill(double Re, double rel_roughness)
{
    /*
     * Churchill (1977) — covers all flow regimes.
     *
     * f = 8.[(8/Re)^12 + 1/(A+B)^1.5]^(1/12)
     * where:
     *   A = [2.457.ln(1/((7/Re)^0.9 + 0.27.epsilon/D))]^16
     *   B = (37530/Re)^16
     *
     * This single equation covers laminar, transitional, and turbulent
     * flows, and matches both the Hagen-Poiseuille limit (f->64/Re) and
     * the Colebrook equation in the turbulent range.
     */
    if (Re <= 0.0 || rel_roughness < 0.0) return -1.0;

    double A = pow(2.457 * log(1.0 / (pow(7.0 / Re, 0.9)
                    + 0.27 * rel_roughness)), 16);
    double B = pow(37530.0 / Re, 16);

    double term = pow(8.0 / Re, 12.0) + 1.0 / pow(A + B, 1.5);
    double f = 8.0 * pow(term, 1.0 / 12.0);

    return f;
}

double form_drag_coefficient(double Re, int geometry)
{
    /*
     * Drag coefficient Cd for sphere (geometry=0) or cylinder (geometry=1).
     *
     * For sphere (Clift-Gauvin correlation):
     *   Re < 1:        Cd = 24/Re
     *   1 < Re < 1000: Cd = (24/Re).(1 + 0.15.Re^0.687)
     *   10^3 < Re < 2.5×10⁵: Cd ≈ 0.44
     *   Re > 2.5×10⁵:  Cd ≈ 0.2
     *
     * For cylinder (cross-flow):
     *   Re < 1:        Cd ≈ 10/Re
     *   1 < Re < 10^3:  Cd ≈ varying
     *   10^3 < Re < 2×10⁵: Cd ≈ 1.0
     *   Re > 2×10⁵:    Cd ≈ 0.3
     */
    if (Re <= 0.0) return -1.0;

    if (geometry == 0) {
        /* Sphere */
        if (Re < 1.0) return 24.0 / Re;
        if (Re < 1000.0) return (24.0 / Re) * (1.0 + 0.15 * pow(Re, 0.687));
        if (Re < 250000.0) return 0.44;
        return 0.2;
    } else if (geometry == 1) {
        /* Cylinder cross-flow */
        if (Re < 1.0) return 10.0 / Re;
        if (Re < 1000.0) {
            /* Approximate fit */
            double logRe = log10(Re);
            double logCd = -0.6 * logRe + 1.0;
            return pow(10.0, logCd);
        }
        if (Re < 200000.0) return 1.0;
        return 0.3;
    }
    return -1.0;
}

/* ==========================================================================
 * L5: Boundary layer thickness estimation
 * ========================================================================== */

double bl_thickness_laminar(double x, double Re_x)
{
    /*
     * delta_99 = 5.0.x / √(Re_x)
     *
     * This is the distance from the wall where u/U_inf = 0.99.
     * Derived from the Blasius similarity solution for laminar
     * flat-plate boundary layer (1908).
     */
    if (Re_x <= 0.0 || x <= 0.0) return -1.0;
    return 5.0 * x / sqrt(Re_x);
}

double bl_thickness_turbulent(double x, double Re_x)
{
    /*
     * delta = 0.37.x / Re_x^(1/5)
     *
     * Based on the 1/7th power-law velocity profile.
     * Valid for 5×10⁵ < Re_x < 10⁷.
     */
    if (Re_x <= 0.0 || x <= 0.0) return -1.0;
    return 0.37 * x / pow(Re_x, 0.2);
}

double bl_displacement_thickness_laminar(double x, double Re_x)
{
    /*
     * delta* = ∫0^inf (1 - u/U_inf) dy
     *
     * For laminar Blasius: delta* = 1.721.x / √(Re_x)
     */
    if (Re_x <= 0.0 || x <= 0.0) return -1.0;
    return 1.721 * x / sqrt(Re_x);
}

double bl_momentum_thickness_laminar(double x, double Re_x)
{
    /*
     * theta = ∫0^inf (u/U_inf).(1 - u/U_inf) dy
     *
     * For laminar Blasius: theta = 0.664.x / √(Re_x)
     *
     * Note: This is the same 0.664 that appears in the
     * flat-plate skin friction: Cf = 0.664/√(Re_x) = 2.theta/x.
     */
    if (Re_x <= 0.0 || x <= 0.0) return -1.0;
    return 0.664 * x / sqrt(Re_x);
}

double bl_transition_location(double Re_crit, double nu, double U)
{
    /*
     * x_crit = Re_crit.nu / U
     *
     * This is the distance from the leading edge where the
     * boundary layer transitions from laminar to turbulent.
     *
     * Example: For air (nu ≈ 1.5×10⁻⁵ m^2/s) at U = 30 m/s,
     *   x_crit ≈ 5×10⁵ × 1.5×10⁻⁵ / 30 ≈ 0.25 m.
     */
    if (nu <= 0.0 || U <= 0.0 || Re_crit <= 0.0) return -1.0;
    return Re_crit * nu / U;
}

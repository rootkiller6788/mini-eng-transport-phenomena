/*
 * nusselt_number.c — Nusselt Number and Convective Heat Transfer Correlations
 *
 * Reference: Incropera & DeWitt (2007), Bejan (2013),
 *            Kays, Crawford & Weigand (2005)
 *
 * Knowledge: L1-L6 complete
 */

#include "../include/nusselt_number.h"
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ==========================================================================
 * L1: Core Nu definition
 * ========================================================================== */

double nu_basic(double h, double L, double k)
{
    if (k <= 0.0 || L < 0.0 || h < 0.0) return -1.0;
    if (k < 1e-15) return -1.0;
    return h * L / k;
}

double h_from_nusselt(double Nu, double L, double k)
{
    if (L <= 0.0 || k <= 0.0 || Nu < 0.0) return -1.0;
    return Nu * k / L;
}

double nu_from_heat_flux(double q_flux, double L, double k, double delta_T)
{
    /*
     * Nu = q".L / (k.DeltaT)
     *
     * At the wall: q" = -k.(dT/dy)|_w = h.DeltaT
     * -> Nu = h.L/k = q".L/(k.DeltaT)
     */
    if (k <= 0.0 || L < 0.0) return -1.0;
    double delta = fabs(delta_T);
    if (delta < 1e-15) return -1.0;
    if (k < 1e-15) return -1.0;
    return q_flux * L / (k * delta);
}

/* ==========================================================================
 * L3-L5: Internal Flow Correlations
 * ========================================================================== */

double nu_pipe_laminar_constant_T(void)
{
    /* Fully developed, laminar, constant T_w: Nu_D = 3.66 */
    return 3.66;
}

double nu_pipe_laminar_constant_qflux(void)
{
    /* Fully developed, laminar, constant q"_w: Nu_D = 4.36 */
    return 4.36;
}

double nu_pipe_laminar_developing(double Re, double Pr, double D, double L,
                                  double mu, double mu_w)
{
    /*
     * Sieder-Tate (1936) for laminar, thermally developing:
     *   Nu_D = 1.86.(Re.Pr.D/L)^(1/3).(mu/mu_w)^0.14
     *
     * Applicability: Nu_D.(mu_w/mu)^0.14 > 2.
     *
     * The term (mu/mu_w)^0.14 corrects for property variation:
     *   Heating liquid:  mu_w < mu -> (mu/mu_w) > 1 -> Nu increased
     *   Cooling liquid:  mu_w > mu -> (mu/mu_w) < 1 -> Nu decreased
     */
    if (Re <= 0.0 || Pr <= 0.0 || D <= 0.0 || L <= 0.0 || mu <= 0.0 || mu_w <= 0.0)
        return -1.0;

    double Gz_inv_term = Re * Pr * D / L;
    if (Gz_inv_term <= 0.0) return -1.0;

    double viscosity_correction = mu_w > 0.0 ? pow(mu / mu_w, 0.14) : 1.0;
    return 1.86 * pow(Gz_inv_term, 1.0 / 3.0) * viscosity_correction;
}

double nu_pipe_dittus_boelter(double Re, double Pr, int heating)
{
    /*
     * Dittus-Boelter (1930):
     *   Nu_D = 0.023.Re^(4/5).Pr^n
     *   n = 0.4 for heating (fluid being heated by wall)
     *   n = 0.3 for cooling (fluid being cooled by wall)
     *
     * Valid: Re > 10000, 0.7 < Pr < 160, L/D > 10.
     *
     * This is the most widely used turbulent heat transfer correlation
     * despite its simplicity and ±25% error margin.
     */
    if (Re < 10000.0 || Pr <= 0.0) return -1.0;

    double n = heating ? 0.4 : 0.3;
    return 0.023 * pow(Re, 0.8) * pow(Pr, n);
}

double nu_pipe_sieder_tate(double Re, double Pr, double mu, double mu_w)
{
    /*
     * Sieder-Tate (1936) for turbulent flow:
     *   Nu_D = 0.027.Re^(4/5).Pr^(1/3).(mu/mu_w)^0.14
     *
     * Compared to Dittus-Boelter:
     *   - Pr exponent fixed at 1/3 (from laminar theory)
     *   - Explicit viscosity ratio correction
     *   - Slightly higher coefficient (0.027 vs 0.023)
     *
     * Valid: Re > 10000, 0.7 < Pr < 16700, L/D > 10.
     */
    if (Re < 10000.0 || Pr <= 0.0 || mu <= 0.0 || mu_w <= 0.0) return -1.0;

    double viscosity_correction = pow(mu / mu_w, 0.14);
    return 0.027 * pow(Re, 0.8) * pow(Pr, 1.0 / 3.0) * viscosity_correction;
}

double nu_pipe_gnielinski(double Re, double Pr, double f)
{
    /*
     * Gnielinski (1976):
     *   Nu_D = (f/8).(Re-1000).Pr / [1 + 12.7.√(f/8).(Pr^(2/3) - 1)]
     *
     * This is the most accurate single-phase turbulent heat transfer
     * correlation for circular pipes. Recommended by most modern
     * heat transfer textbooks.
     *
     * Valid: 3000 < Re < 5×10⁶, 0.5 < Pr < 2000.
     *
     * Note: Uses the Darcy friction factor f (not the Fanning Cf = f/4).
     * Error: ±10% for most fluids.
     */
    if (Re < 3000.0 || Pr <= 0.0 || f <= 0.0) return -1.0;

    double f_over_8 = f / 8.0;
    double numerator = f_over_8 * (Re - 1000.0) * Pr;
    double denominator = 1.0 + 12.7 * sqrt(f_over_8) * (pow(Pr, 2.0 / 3.0) - 1.0);

    if (denominator < 1e-15) return -1.0;
    return numerator / denominator;
}

double nu_annulus(double Re, double Pr, double Di, double Do, int inner_heated)
{
    /*
     * Annular duct: Hydraulic diameter D_h = Do - Di.
     *
     * For turbulent flow, Dittus-Boelter is used with corrections:
     *   Inner wall heating: Nu_ii/Nu_DB = 0.86.(Do/Di)^0.16
     *   Outer wall heating: Nu_oo/Nu_DB = 1.0 - 0.14.(Di/Do)^0.6
     *
     * These corrections account for the asymmetric velocity and
     * temperature profiles in annular flow.
     */
    if (Re < 10000.0 || Pr <= 0.0 || Di <= 0.0 || Do <= Di) return -1.0;

    double Re_Dh = Re; /* Assumes Re already based on hydraulic diameter Dh = Do - Di */
    double Nu_db = 0.023 * pow(Re_Dh, 0.8) * pow(Pr, 0.4);

    if (inner_heated) {
        return Nu_db * 0.86 * pow(Do / Di, 0.16);
    } else {
        return Nu_db * (1.0 - 0.14 * pow(Di / Do, 0.6));
    }
}

double nu_noncircular_duct_laminar(int geometry_id)
{
    /*
     * Fully developed laminar, constant T_w for non-circular ducts.
     *
     * Values from analytical solutions of the energy equation:
     */
    switch (geometry_id) {
        case 0: return 3.66;   /* circular */
        case 1: return 2.98;   /* square */
        case 2: return 2.47;   /* triangular (equilateral) */
        case 3: return 7.54;   /* parallel plates (b >> a) */
        case 4: return 3.39;   /* rectangular 2:1 aspect ratio */
        case 5: return 4.44;   /* rectangular 4:1 aspect ratio */
        default: return -1.0;
    }
}

/* ==========================================================================
 * L3-L5: External Flow Correlations
 * ========================================================================== */

double nu_flat_plate_laminar_local(double Re_x, double Pr)
{
    /*
     * Blasius/Pohlhausen similarity solution:
     *   Nu_x = 0.332.Re_x^(1/2).Pr^(1/3)
     *
     * Valid for Pr > 0.6, Re_x < 5×10⁵.
     *
     * Derivation via similarity: u* = f'(eta), T* = g(eta)
     * where eta = y.√(U/(nux)). The ODE solution yields Nu_x.
     */
    if (Re_x <= 0.0 || Pr <= 0.0) return -1.0;
    return 0.332 * sqrt(Re_x) * cbrt(Pr);
}

double nu_flat_plate_laminar_avg(double Re_L, double Pr)
{
    /*
     * Average Nusselt number: Nu_L(avg) = 2.Nu_L(local)
     *   Nu_L = 0.664.Re_L^(1/2).Pr^(1/3)
     *
     * Proof: Nu_avg = (1/L).∫0^L Nu_x dx
     *   = (1/L).∫0^L const.√(U/(nux)) dx
     *   = (const/L).2√(UL/nu) = 2.Nu_L(local)
     */
    if (Re_L <= 0.0 || Pr <= 0.0) return -1.0;
    return 0.664 * sqrt(Re_L) * cbrt(Pr);
}

double nu_flat_plate_turbulent_local(double Re_x, double Pr)
{
    /*
     * Turbulent flat plate (1/7th power law + Reynolds analogy):
     *   Nu_x = 0.0296.Re_x^(4/5).Pr^(1/3)
     *
     * Valid: 5×10⁵ < Re_x < 10⁷, 0.6 < Pr < 60.
     */
    if (Re_x <= 0.0 || Pr <= 0.0) return -1.0;
    return 0.0296 * pow(Re_x, 0.8) * cbrt(Pr);
}

double nu_flat_plate_turbulent_avg(double Re_L, double Pr)
{
    /*
     * Averaging the turbulent local correlation:
     *   Nu_L = 0.037.Re_L^(4/5).Pr^(1/3)
     */
    if (Re_L <= 0.0 || Pr <= 0.0) return -1.0;
    return 0.037 * pow(Re_L, 0.8) * cbrt(Pr);
}

double nu_flat_plate_mixed(double Re_L, double Pr, double Re_crit)
{
    /*
     * Mixed laminar-turbulent boundary layer:
     *   Nu_L = [0.037.Re_L^(4/5) - 871].Pr^(1/3)
     *
     * Derived by integrating:
     *   ∫0^{x_crit} Nu_x,lam dx + ∫_{x_crit}^L Nu_x,turb dx
     *
     * The constant 871 = 0.037.Re_crit^(4/5) - 0.664.Re_crit^(1/2)
     * evaluated at Re_crit = 5×10⁵.
     *
     * For Re_crit ≠ 5×10⁵, we use the general formula:
     *   Nu_L = [0.037.Re_L^(4/5) - A].Pr^(1/3)
     *   where A = 0.037.Re_crit^(4/5) - 0.664.Re_crit^(1/2)
     *
     * Valid: 0.6 < Pr < 60, Re_crit < Re_L < 10⁸.
     */
    if (Re_L <= 0.0 || Pr <= 0.0 || Re_crit <= 0.0) return -1.0;
    if (Re_L <= Re_crit) {
        /* Entirely laminar */
        return nu_flat_plate_laminar_avg(Re_L, Pr);
    }

    double A_term = 0.037 * pow(Re_crit, 0.8) - 0.664 * sqrt(Re_crit);
    double Nu = (0.037 * pow(Re_L, 0.8) - A_term) * cbrt(Pr);
    return (Nu > 0.0) ? Nu : nu_flat_plate_laminar_avg(Re_L, Pr);
}

double nu_flat_plate_constant_heatflux_laminar(double Re_x, double Pr)
{
    /*
     * Constant wall heat flux, laminar:
     *   Nu_x = 0.453.Re_x^(1/2).Pr^(1/3)
     *
     * For uniform heat flux, the surface temperature adjusts to
     * match the specified flux, resulting in a larger Nusselt number
     * compared to uniform temperature (0.453 vs 0.332).
     *
     * Valid for Pr > 0.6.
     */
    if (Re_x <= 0.0 || Pr <= 0.0) return -1.0;
    return 0.453 * sqrt(Re_x) * cbrt(Pr);
}

double nu_cylinder_crossflow(double Re_D, double Pr)
{
    /*
     * Zhukauskas (1972) correlation for cylinder in cross-flow:
     *   Nu_D = C.Re_D^m.Pr^n.(Pr/Pr_w)^0.25
     *   n = 0.37 for Pr ≤ 10, n = 0.36 for Pr > 10
     *
     * Piecewise C, m based on Re_D:
     */
    if (Re_D <= 0.0 || Pr <= 0.0) return -1.0;

    double C, m, n;
    n = (Pr <= 10.0) ? 0.37 : 0.36;

    if (Re_D < 4.0)        { C = 0.989; m = 0.330; }
    else if (Re_D < 40.0)  { C = 0.911; m = 0.385; }
    else if (Re_D < 4000.0){ C = 0.683; m = 0.466; }
    else if (Re_D < 40000.0){C = 0.193; m = 0.618; }
    else                    { C = 0.027; m = 0.805; }

    /* Using Pr_w ≈ Pr as first approximation */
    return C * pow(Re_D, m) * pow(Pr, n) * 1.0;
}

double nu_sphere(double Re, double Pr, double mu, double mu_w)
{
    /*
     * Whitaker (1972):
     *   Nu_D = 2 + (0.4.Re^(1/2) + 0.06.Re^(2/3)).Pr^0.4.(mu/mu_w)^0.25
     *
     * The constant 2 comes from the pure conduction limit:
     *   For a sphere in an infinite stagnant fluid (Re -> 0):
     *   Q = 2π.k.D.DeltaT -> Nu = hD/k = 2
     *
     * Valid: 3.5 < Re < 7.6×10⁴, 0.7 < Pr < 380.
     */
    if (Re < 0.0 || Pr <= 0.0 || mu <= 0.0 || mu_w <= 0.0) return -1.0;

    double Re_term = 0.4 * sqrt(Re) + 0.06 * pow(Re, 2.0 / 3.0);
    double viscosity_correction = pow(mu / mu_w, 0.25);

    return 2.0 + Re_term * pow(Pr, 0.4) * viscosity_correction;
}

double nu_tube_bank(double Re_max, double Pr, double Pr_w,
                    int arrangement, int n_rows)
{
    /*
     * Tube bank correlation (Zukauskas, 1987):
     *   Nu_D = C.Re_max^m.Pr^0.36.(Pr/Pr_w)^0.25
     *
     * Arrangement: 0 = aligned, 1 = staggered.
     * Correction for fewer than 20 tube rows.
     */
    if (Re_max <= 0.0 || Pr <= 0.0 || Pr_w <= 0.0 || n_rows < 1) return -1.0;

    double C, m;

    if (arrangement == 0) { /* Aligned */
        if (Re_max < 100.0)        { C = 0.80; m = 0.40; }
        else if (Re_max < 1000.0)  { C = 0.51; m = 0.50; }
        else if (Re_max < 200000.0){ C = 0.27; m = 0.63; }
        else                        { C = 0.021; m = 0.84; }
    } else { /* Staggered */
        if (Re_max < 100.0)        { C = 0.90; m = 0.40; }
        else if (Re_max < 1000.0)  { C = 0.51; m = 0.50; }
        else if (Re_max < 200000.0){ C = 0.35 * 0.95; m = 0.60; }
        else                        { C = 0.022; m = 0.84; }
    }

    double Nu = C * pow(Re_max, m) * pow(Pr, 0.36) * pow(Pr / Pr_w, 0.25);

    /* Correction for number of tube rows */
    if (n_rows < 20) {
        double row_correction;
        if (n_rows == 1) row_correction = 0.70;
        else if (n_rows == 2) row_correction = 0.80;
        else if (n_rows == 3) row_correction = 0.87;
        else if (n_rows == 4) row_correction = 0.90;
        else if (n_rows == 5) row_correction = 0.92;
        else if (n_rows == 6) row_correction = 0.94;
        else if (n_rows == 8) row_correction = 0.96;
        else if (n_rows == 10) row_correction = 0.98;
        else if (n_rows == 12) row_correction = 0.99;
        else row_correction = 0.995; /* 14-19 rows */

        Nu *= row_correction;
    }

    return Nu;
}

/* ==========================================================================
 * L3-L5: Natural Convection Correlations
 * ========================================================================== */

double nu_natural_convection(double Ra, int geometry_id, int flow_regime)
{
    /*
     * General form: Nu = C.Ra^m
     *
     * Vertical plate (geometry_id = 0):
     *   Laminar  (10⁴ < Ra < 10⁹):    C=0.59, m=1/4
     *   Turbulent (10⁹ < Ra < 10¹^3):  C=0.10, m=1/3
     */
    if (Ra <= 0.0) return -1.0;

    double C, m;
    switch (geometry_id) {
        case 0: /* vertical plate */
            if (flow_regime == 0) { C = 0.59; m = 0.25; }
            else                  { C = 0.10; m = 1.0 / 3.0; }
            break;
        case 1: /* horizontal plate, heated up */
            if (flow_regime == 0) { C = 0.54; m = 0.25; }
            else                  { C = 0.15; m = 1.0 / 3.0; }
            break;
        case 2: /* horizontal plate, heated down */
            if (flow_regime == 0) { C = 0.27; m = 0.25; }
            else                 { C = 0.07; m = 1.0 / 3.0; }
            break;
        case 3: /* horizontal cylinder */
            C = 0.53; m = 0.25;
            break;
        case 4: /* sphere */
            C = 0.589; m = 0.25;
            break;
        case 5: /* vertical cylinder (L >> D) */
            C = 0.59; m = 0.25;
            break;
        default: return -1.0;
    }

    return C * pow(Ra, m);
}

double nu_natural_vertical_plate_churchill_chu(double Ra, double Pr)
{
    /*
     * Churchill & Chu (1975):
     *   Nu_L = {0.825 + 0.387.Ra^(1/6)/[1 + (0.492/Pr)^(9/16)]^(8/27)}^2
     *
     * This single correlation covers both laminar and turbulent
     * natural convection regimes for a vertical plate. It is
     * the most accurate correlation for this geometry.
     *
     * Valid for all Ra, 0 < Pr < inf.
     */
    if (Ra <= 0.0 || Pr <= 0.0) return -1.0;

    double denominator = pow(1.0 + pow(0.492 / Pr, 9.0 / 16.0), 8.0 / 27.0);
    double term = 0.825 + 0.387 * pow(Ra, 1.0 / 6.0) / denominator;
    return term * term;
}

double nu_natural_enclosure(double Ra, double aspect, double angle)
{
    /*
     * Enclosure natural convection.
     *
     * For vertical enclosure (heated from side), angle ≈ 90 deg,
     * the Berkovsky-Polevikov correlation:
     *   Nu = 0.22.(Pr.Ra/(0.2+Pr))^0.28.(H/L)^(-1/4)
     *   Valid: 2 < H/L < 10, Pr < 10⁵, Ra < 10¹⁰.
     *
     * For horizontal enclosure (heated from below), angle ≈ 0 deg,
     * onset of convection at Ra > 1708.
     *
     * Nu = max(1.0, Nu_correlated)
     * Pure conduction gives Nu = 1.
     */
    if (Ra <= 0.0 || aspect <= 0.0) return -1.0;

    /* Simplified version for vertical enclosure */
    double angle_rad = angle * M_PI / 180.0;
    double Pr_typical = 0.7; /* air — for this simplified version */

    if (fabs(angle - 90.0) < 5.0 || fabs(sin(angle_rad)) > 0.99) {
        /* Vertical or near-vertical */
        double Nu = 0.22 * pow(Pr_typical * Ra / (0.2 + Pr_typical), 0.28)
                      * pow(aspect, -0.25);
        return (Nu > 1.0) ? Nu : 1.0;
    } else if (fabs(angle) < 5.0 || fabs(cos(angle_rad)) > 0.99) {
        /* Horizontal */
        if (Ra < 1708.0) {
            return 1.0; /* Pure conduction, no convection */
        }
        return 0.059 * pow(Ra, 1.0 / 3.0);
    } else {
        /* Inclined — use vertical as approximation */
        double Nu = 0.22 * pow(Pr_typical * Ra / (0.2 + Pr_typical), 0.28)
                      * pow(aspect, -0.25);
        return (Nu > 1.0) ? Nu : 1.0;
    }
}

/* ==========================================================================
 * L4: Energy equation non-dimensionalization
 * ========================================================================== */

int energy_eq_dimensionless_form(double Re, double Pr,
                                 double mu, double U, double k, double delta_T,
                                 double *Pe_out, double *Br_out)
{
    /*
     * Energy equation (dimensional):
     *   rho.cp.DT/Dt = k.grad^2T + mu.Phi + q""
     *
     * Non-dimensionalize with:
     *   x* = x/L, u* = u/U, t* = t.U/L, T* = (T-T0)/DeltaT
     *
     * Result:
     *   Re.Pr.DT_star/Dt_star = grad*^2T* + Br.Phi* + q*""
     *
     *   Pe_H = Re.Pr = rho.cp.U.L/k (advective/diffusive heat transport)
     *   Br   = mu.U^2/(k.DeltaT) (viscous dissipation/heat conduction)
     *
     * Note: sometimes written using Eckert number Ec = U^2/(cp.DeltaT),
     * so Br = Pr.Ec.
     */
    if (Re < 0.0 || Pr < 0.0 || !Pe_out || !Br_out) return -1;

    *Pe_out = Re * Pr;

    double delta = fabs(delta_T);
    if (delta < 1e-15 || k < 1e-15 || mu < 0.0) {
        *Br_out = 0.0;
        return 0;
    }

    *Br_out = mu * U * U / (k * delta);
    return 0;
}

double colburn_j_factor_heat(double Nu, double Re, double Pr)
{
    /*
     * Colburn j-factor for heat transfer:
     *   j_H = St_H × Pr^(2/3) = (Nu/(Re.Pr)) × Pr^(2/3) = Nu/(Re.Pr^(1/3))
     *
     * Chilton-Colburn analogy:
     *   j_H = j_D = f/8
     *
     * This relates momentum, heat, and mass transfer in a single
     * dimensionless parameter. It's the engineering foundation of
     * the heat-mass-momentum transfer analogy.
     *
     * Valid for: 0.6 < Pr < 60, turbulent flow.
     */
    if (Re <= 0.0 || Pr <= 0.0 || Nu < 0.0) return -1.0;
    return Nu / (Re * pow(Pr, 1.0 / 3.0));
}

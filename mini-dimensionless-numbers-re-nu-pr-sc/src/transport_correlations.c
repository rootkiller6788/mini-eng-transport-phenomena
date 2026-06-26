/*
 * transport_correlations.c — Engineering Transport Correlations
 *
 * Reference: Incropera & DeWitt (2007), Cengel (2007),
 *            Bejan (2013), Bird-Stewart-Lightfoot (2007)
 *
 * Knowledge: L5-L7 complete
 */

#include "../include/transport_correlations.h"
#include "../include/dimensionless_numbers.h"
#include "../include/nusselt_number.h"
#include "../include/reynolds_number.h"
#include "../include/prandtl_schmidt.h"
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ==========================================================================
 * L5: Forced Convection Correlation Engine
 * ========================================================================== */

int forced_convection_correlation(const CorrelationInput *input,
                                  int geometry, int bc_type,
                                  CorrelationResult *result)
{
    if (!input || !result) return -1;
    if (input->rho <= 0.0 || input->mu <= 0.0 || input->k <= 0.0 ||
        input->L <= 0.0) return -1;

    /* Compute Re and Pr */
    double Re = compute_reynolds_number(input->rho, input->U, input->L, input->mu);
    double Pr = input->Pr;
    if (Pr <= 0.0) {
        double alpha = input->k / (input->rho * input->cp);
        if (alpha <= 0.0) return -1;
        Pr = compute_prandtl_number(input->mu / input->rho, alpha);
    }

    if (Re < 0.0 || Pr < 0.0) return -1;

    double delta_T = input->T_w - input->T_b;
    if (fabs(delta_T) < 1e-15) return -1;

    double Nu = 0.0;

    switch (geometry) {
        case 0: /* Flat plate */
            if (Re < RE_CRIT_FLAT_PLATE) {
                /* Laminar */
                if (bc_type == 0)
                    Nu = nu_flat_plate_laminar_avg(Re, Pr);
                else
                    Nu = nu_flat_plate_constant_heatflux_laminar(Re, Pr);
            } else if (Re < 5.0e6) {
                /* Mixed laminar-turbulent */
                Nu = nu_flat_plate_mixed(Re, Pr, RE_CRIT_FLAT_PLATE);
            } else {
                Nu = nu_flat_plate_turbulent_avg(Re, Pr);
            }
            result->error_estimate = 15.0; /* ±15% */
            break;

        case 1: /* Circular pipe */
            if (Re < 2300.0) {
                if (bc_type == 0)
                    Nu = nu_pipe_laminar_constant_T();
                else
                    Nu = nu_pipe_laminar_constant_qflux();
            } else if (Re > 10000.0) {
                Nu = nu_pipe_dittus_boelter(Re, Pr,
                        (delta_T > 0.0) ? 1 : 0);
            } else {
                /* Transitional — interpolate */
                double Nu_lam = (bc_type == 0) ?
                    nu_pipe_laminar_constant_T() :
                    nu_pipe_laminar_constant_qflux();
                double Nu_turb = nu_pipe_dittus_boelter(Re, Pr,
                        (delta_T > 0.0) ? 1 : 0);
                double frac = (Re - 2300.0) / (10000.0 - 2300.0);
                Nu = Nu_lam + frac * (Nu_turb - Nu_lam);
            }
            result->error_estimate = 25.0; /* Dittus-Boelter ±25% */
            break;

        case 2: /* Cylinder crossflow */
            Nu = nu_cylinder_crossflow(Re, Pr);
            result->error_estimate = 20.0;
            break;

        case 3: /* Sphere */
            Nu = nu_sphere(Re, Pr, input->mu, input->mu);
            result->error_estimate = 25.0;
            break;

        case 4: /* Tube bank aligned */
        case 5: /* Tube bank staggered */
            Nu = nu_tube_bank(Re, Pr, Pr, (geometry == 5) ? 1 : 0, 20);
            result->error_estimate = 25.0;
            break;

        default:
            return -1;
    }

    if (Nu <= 0.0) return -1;

    result->Nu = Nu;
    result->h = h_from_nusselt(Nu, input->L, input->k);
    result->q_flux = result->h * fabs(delta_T);
    result->total_power = result->q_flux * input->A_s;

    return 0;
}

/* ==========================================================================
 * L5: Natural Convection
 * ========================================================================== */

int natural_convection_vertical_plate(double k, double L, double Ra, double Pr,
                                      double *Nu, double *h)
{
    if (k <= 0.0 || L <= 0.0 || Ra <= 0.0 || Pr <= 0.0 || !Nu || !h)
        return -1;

    *Nu = nu_natural_vertical_plate_churchill_chu(Ra, Pr);
    *h = h_from_nusselt(*Nu, L, k);

    return 0;
}

int natural_convection_horizontal_plate(double k, double L, double Ra,
                                        int heated_up, double *Nu, double *h)
{
    if (k <= 0.0 || L <= 0.0 || Ra <= 0.0 || !Nu || !h) return -1;

    int regime = (Ra > 1.0e7) ? 1 : 0;
    (void)heated_up; /* used via geom_id below */
    int geom_id = heated_up ? 1 : 2;

    *Nu = nu_natural_convection(Ra, geom_id, regime);
    *h = h_from_nusselt(*Nu, L, k);

    return 0;
}

int natural_convection_cylinder(double k, double D, double Ra, double Pr,
                                double *Nu, double *h)
{
    if (k <= 0.0 || D <= 0.0 || Ra <= 0.0 || Pr <= 0.0 || !Nu || !h)
        return -1;
    /* Use Churchill-Chu adapted for cylinder */
    *Nu = (0.60 + 0.387 * pow(Ra, 1.0/6.0)
           / pow(1.0 + pow(0.559/Pr, 9.0/16.0), 8.0/27.0));
    *Nu = (*Nu) * (*Nu);
    *h = h_from_nusselt(*Nu, D, k);
    return 0;
}

int natural_convection_sphere(double k, double D, double Ra, double Pr,
                              double *Nu, double *h)
{
    if (k <= 0.0 || D <= 0.0 || Ra <= 0.0 || Pr <= 0.0 || !Nu || !h)
        return -1;
    *Nu = 2.0 + 0.589 * pow(Ra, 0.25) / pow(1.0 + pow(0.469/Pr, 9.0/16.0), 4.0/9.0);
    *h = h_from_nusselt(*Nu, D, k);
    return 0;
}

int natural_convection_enclosure(double k, double L, double H,
                                 double T_h, double T_c,
                                 double Ra, double Pr, double angle,
                                 double *Nu, double *h_eff, double *q_flux)
{
    if (k <= 0.0 || L <= 0.0 || H <= 0.0 || !Nu || !h_eff || !q_flux)
        return -1;
    if (T_h <= T_c) return -1;

    (void)Pr;  /* Pr reserved for future Pr-dependent enclosure correlation */
    double aspect = H / L;
    *Nu = nu_natural_enclosure(Ra, aspect, angle);
    *h_eff = h_from_nusselt(*Nu, L, k);
    *q_flux = (*h_eff) * (T_h - T_c);

    return 0;
}

/* ==========================================================================
 * L5-L6: Mass Transfer
 * ========================================================================== */

int mass_transfer_flat_plate(double Re_L, double Sc, double D_AB, double L,
                             double *Sh, double *h_m)
{
    if (Re_L <= 0.0 || Sc <= 0.0 || D_AB <= 0.0 || L <= 0.0 || !Sh || !h_m)
        return -1;

    *Sh = 0.664 * sqrt(Re_L) * cbrt(Sc);
    *h_m = (*Sh) * D_AB / L;

    return 0;
}

int mass_transfer_pipe(double Re, double Sc, double D, double L, double D_AB,
                       double *Sh, double *h_m)
{
    if (Re <= 0.0 || Sc <= 0.0 || D <= 0.0 || L <= 0.0 || D_AB <= 0.0 ||
        !Sh || !h_m) return -1;

    if (Re < 2300.0) {
        /* Laminar, developing: Lévêque solution */
        *Sh = 1.62 * cbrt(Re * Sc * D / L);
    } else {
        /* Turbulent */
        *Sh = 0.023 * pow(Re, 0.83) * pow(Sc, 0.44);
    }

    *h_m = (*Sh) * D_AB / D;
    return 0;
}

int mass_transfer_sphere(double Re, double Sc, double D, double D_AB,
                         double *Sh, double *h_m)
{
    /* Frössling (1938): Sh = 2 + 0.552.Re^(1/2).Sc^(1/3) */
    if (Re < 0.0 || Sc <= 0.0 || D <= 0.0 || D_AB <= 0.0 || !Sh || !h_m)
        return -1;

    *Sh = 2.0 + 0.552 * sqrt(Re) * cbrt(Sc);
    *h_m = (*Sh) * D_AB / D;

    return 0;
}

int mass_transfer_droplet_evaporation(double Re, double Sc, double D_droplet,
                                      double D_AB, double rho_s, double rho_inf,
                                      double *m_dot_evap)
{
    /*
     * Evaporating droplet — Ranz-Marshall (1952):
     *   Sh = 2 + 0.6.Re^(1/2).Sc^(1/3)
     *   ṁ_evap = Sh.π.D.D_AB.(rho_s - rho_inf)
     *
     * The Sh->2 limit corresponds to pure diffusion from a sphere
     * into a stagnant medium. The Re term adds convective enhancement.
     *
     * Applications: spray drying, fuel injection, atmospheric water droplets,
     *   COVID aerosol transmission modeling.
     */
    if (Re < 0.0 || Sc <= 0.0 || D_droplet <= 0.0 || D_AB <= 0.0 || !m_dot_evap)
        return -1;

    double Sh = 2.0 + 0.6 * sqrt(Re) * cbrt(Sc);
    double delta_rho = rho_s - rho_inf;
    if (delta_rho < 0.0) delta_rho = fabs(delta_rho);

    *m_dot_evap = Sh * M_PI * D_droplet * D_AB * delta_rho;

    return 0;
}

/* ==========================================================================
 * L6: Heat Exchanger Sizing
 * ========================================================================== */

int heat_exchanger_sizing(double Q_dot,
                          double T_h_in, double T_h_out,
                          double T_c_in, double T_c_out,
                          double h_hot, double h_cold,
                          double k_wall, double t_wall,
                          double *A_required, double *U_overall,
                          double *delta_T_lm)
{
    /*
     * Heat exchanger design using LMTD method.
     *
     * DeltaT_lm = (DeltaT1 - DeltaT2) / ln(DeltaT1/DeltaT2)
     *   where DeltaT1 = T_h_in - T_c_out
     *         DeltaT2 = T_h_out - T_c_in
     *
     * Overall heat transfer coefficient:
     *   1/U = 1/h_hot + t_wall/k_wall + 1/h_cold
     *
     * Required area:
     *   A = Q / (U.DeltaT_lm)
     */
    if (Q_dot <= 0.0 || h_hot <= 0.0 || h_cold <= 0.0 ||
        !A_required || !U_overall || !delta_T_lm) return -1;

    double dT1 = T_h_in - T_c_out;
    double dT2 = T_h_out - T_c_in;

    if (dT1 <= 0.0 || dT2 <= 0.0) return -1; /* temperature cross */
    if (fabs(dT1 / dT2 - 1.0) < 1e-10) {
        *delta_T_lm = dT1; /* L'Hôpital limit */
    } else {
        *delta_T_lm = (dT1 - dT2) / log(dT1 / dT2);
    }

    /* Overall U */
    double R_total = 1.0 / h_hot + t_wall / k_wall + 1.0 / h_cold;
    *U_overall = 1.0 / R_total;

    /* Required area */
    *A_required = Q_dot / ((*U_overall) * (*delta_T_lm));

    return 0;
}

/* ==========================================================================
 * L6: Pipe Flow Design
 * ========================================================================== */

int pipe_flow_design(double rho, double mu, double m_dot,
                     double D, double L, double epsilon,
                     double *Re, double *f, double *delta_P,
                     double *W_dot_pump)
{
    if (rho <= 0.0 || mu <= 0.0 || m_dot < 0.0 || D <= 0.0 || L < 0.0 ||
        epsilon < 0.0 || !Re || !f || !delta_P || !W_dot_pump) return -1;

    /* Re from mass flow */
    *Re = re_pipe_from_mass_flow(m_dot, D, mu);
    if (*Re <= 0.0) return -1;

    /* Friction factor using Churchill (all-regime) */
    double rel_rough = epsilon / D;
    *f = friction_factor_churchill(*Re, rel_rough);

    /* Velocity */
    double A_cross = M_PI * D * D / 4.0;
    double U = m_dot / (rho * A_cross);

    /* Pressure drop: Darcy-Weisbach */
    *delta_P = (*f) * (L / D) * (0.5 * rho * U * U);

    /* Pumping power */
    double V_dot = m_dot / rho; /* volumetric flow rate */
    *W_dot_pump = (*delta_P) * V_dot; /* ideal minimum power */

    return 0;
}

/* ==========================================================================
 * L6: Electronics Cooling — Finned Array
 * ========================================================================== */

int electronics_cooling_finned_array(double k_fluid, double nu, double alpha,
                                     double beta, double T_s, double T_amb,
                                     double L_fin, double H_fin, double t_fin,
                                     double S, int n_fins, double k_fin,
                                     double W_base,
                                     double *Q_total, double *h_bar,
                                     double *eta_fin)
{
    /*
     * Natural convection from vertical fin array (Elenbaas, 1942).
     *
     * This is the classic electronics cooling problem:
     * heat sink with parallel vertical fins, cooled by natural convection.
     *
     * The Elenbaas correlation for isothermal parallel plates:
     *   Nu_S = (1/24).Ra_S.(S/L).{1 - exp[-35/(Ra_S.S/L)]}^(3/4)
     * where Ra_S = g.beta.(T_s-T_amb).S^3/(nu.alpha)
     *
     * Fin efficiency:
     *   m = √(2h/(k_fin.t_fin))
     *   eta_fin = tanh(m.H_fin) / (m.H_fin)
     *
     * Total heat transfer:
     *   Q = h.[n_fins.2.H_fin.L_fin.eta_fin + (W_base - n_fins.t_fin).L_fin]
     *        . (T_s - T_amb)
     *
     * L7 Application: Toyota ECU, iPhone thermal management,
     * Tesla power electronics cooling.
     */
    double g = 9.81;
    double delta_T = T_s - T_amb;
    if (delta_T <= 0.0 || k_fluid <= 0.0 || nu <= 0.0 || alpha <= 0.0 ||
        beta <= 0.0 || L_fin <= 0.0 || H_fin <= 0.0 || t_fin <= 0.0 ||
        S <= 0.0 || n_fins < 1 || k_fin <= 0.0 || W_base <= 0.0 ||
        !Q_total || !h_bar || !eta_fin) return -1;

    /* Rayleigh number based on fin spacing */
    double Ra_S = g * beta * delta_T * S * S * S / (nu * alpha);

    /* Elenbaas correlation */
    double arg = 35.0 / (Ra_S * S / L_fin);
    if (arg > 50.0) arg = 50.0; /* prevent underflow */
    double Nu_S = (1.0 / 24.0) * Ra_S * (S / L_fin) *
                  pow(1.0 - exp(-arg), 0.75);

    *h_bar = Nu_S * k_fluid / S;

    /* Fin efficiency */
    double m = sqrt(2.0 * (*h_bar) / (k_fin * t_fin));
    double mH = m * H_fin;
    if (mH > 50.0) mH = 50.0; /* prevent overflow */
    *eta_fin = tanh(mH) / mH;

    /* Total heat transfer */
    double A_fins = n_fins * 2.0 * H_fin * L_fin * (*eta_fin);
    double A_base = (W_base - n_fins * t_fin) * L_fin;
    double A_total = A_fins + A_base;

    *Q_total = (*h_bar) * A_total * delta_T;

    return 0;
}

/* ==========================================================================
 * L7: Air-Cooled Condenser
 * ========================================================================== */

int air_cooled_condenser(double T_steam, double T_air_in,
                         double m_dot_steam, double h_fg,
                         double D_tube, double L_tube,
                         int n_tubes, int n_rows,
                         double U_air, double S_T,
                         double S_L __attribute__((unused)),
                         double *Q_condenser, double *T_air_out,
                         double *A_required)
{
    /*
     * Air-cooled condenser — steam condensing inside horizontal tubes,
     * air flowing across tube bank.
     *
     * Uses: power plant dry cooling, Boeing 747 ECS condensers,
     * industrial refrigeration.
     *
     * Assumptions:
     *   - Steam side: condensing at constant T (h_steam ~ 5000-15000 W/m^2K)
     *   - Air side: cross-flow over tube bank (Zukauskas correlation)
     *   - In-line tube arrangement
     */
    if (T_steam <= T_air_in || m_dot_steam <= 0.0 || h_fg <= 0.0 ||
        D_tube <= 0.0 || L_tube <= 0.0 || n_tubes < 1 || n_rows < 1 ||
        U_air <= 0.0 || !Q_condenser || !T_air_out || !A_required) return -1;

    /* Air properties at film temperature ~350K (approximate) */
    (void)(T_steam + T_air_in); /* Film temp ~ (T_steam + T_air_in)/2 */
    double k_air = 0.028;     /* W/(m.K) */
    double nu_air = 2.0e-5;   /* m^2/s */
    double Pr_air = 0.70;
    double rho_air = 1.0;     /* kg/m^3 at ~350K, 1 atm */
    double cp_air = 1007.0;   /* J/(kg.K) */

    /* Air-side Re based on max velocity in minimum flow area */
    double A_min = L_tube * (S_T - D_tube) * (n_tubes / n_rows);
    if (A_min <= 0.0) A_min = L_tube * D_tube * 0.5;
    double U_max = U_air * S_T / (S_T - D_tube);
    double Re_max = U_max * D_tube / nu_air;

    /* Tube bank Nu */
    double Nu_air = nu_tube_bank(Re_max, Pr_air, Pr_air, 0, n_rows);
    double h_air = Nu_air * k_air / D_tube;

    /* Overall U (steam side resistance negligible) */
    double U_overall = h_air; /* Simplified single-phase air dominates */
    /* A_total = n_tubes * pi * D_tube * L_tube; reserved for area check */

    /* LMTD */
    /* For condenser: T_steam = constant, air heats from T_in to T_out */
    /* Need iterative solution. Use simplified energy balance. */
    *Q_condenser = m_dot_steam * h_fg; /* latent heat */

    /* Air-side energy balance: Q = m_dot_air * cp * (T_out - T_in) */
    double m_dot_air = rho_air * U_air * n_rows * S_T * L_tube;
    *T_air_out = T_air_in + (*Q_condenser) / (m_dot_air * cp_air);

    /* Check if physically possible */
    if (*T_air_out >= T_steam) {
        *T_air_out = T_steam - 5.0; /* approach temperature */
    }

    /* LMTD */
    double dT1 = T_steam - T_air_in;
    double dT2 = T_steam - (*T_air_out);
    double LMTD;
    if (dT1 <= 0.0 || dT2 <= 0.0) {
        LMTD = 10.0; /* fallback */
    } else if (fabs(dT1 - dT2) < 1e-10) {
        LMTD = dT1;
    } else {
        LMTD = (dT1 - dT2) / log(dT1 / dT2);
    }

    *A_required = (*Q_condenser) / (U_overall * LMTD);

    return 0;
}

/* ==========================================================================
 * L7: DC Motor Cooling
 * ========================================================================== */

int dc_motor_natural_cooling(double T_amb, double P_loss,
                             double D_motor, double L_motor,
                             double k_air, double nu_air, double alpha_air,
                             double *T_surface, double *Ra_D, double *Nu_D)
{
    /*
     * Natural convection from horizontal cylindrical motor housing.
     *
     * The motor housing acts as a horizontal cylinder with natural
     * convection on its surface, dissipating I^2R losses.
     *
     * Churchill-Chu for horizontal cylinder:
     *   Nu_D = (0.60 + 0.387.Ra_D^(1/6)/[1+(0.559/Pr)^(9/16)]^(8/27))^2
     *
     * Applications: Detroit auto assembly line servo motors,
     * Toyota hybrid motor-generators, industrial robotics.
     */
    if (P_loss <= 0.0 || D_motor <= 0.0 || L_motor <= 0.0 ||
        k_air <= 0.0 || nu_air <= 0.0 || alpha_air <= 0.0 ||
        !T_surface || !Ra_D || !Nu_D) return -1;

    double g = 9.81;
    double beta = 1.0 / T_amb; /* ideal gas approximation */
    double Pr = nu_air / alpha_air;

    /* Iterative solution: T_surface appears in Ra */
    /* Simplified: assume initial T guess, iterate once */
    double T_guess = T_amb + 40.0; /* typical motor temp rise ~40K */
    double A_surface = M_PI * D_motor * L_motor;

    for (int iter = 0; iter < 5; iter++) {
        double delta_T = T_guess - T_amb;
        if (delta_T <= 0.0) delta_T = 1.0;

        *Ra_D = g * beta * delta_T * D_motor * D_motor * D_motor
                / (nu_air * alpha_air);

        /* Churchill-Chu for horizontal cylinder */
        double denom = pow(1.0 + pow(0.559 / Pr, 9.0 / 16.0), 8.0 / 27.0);
        *Nu_D = 0.60 + 0.387 * pow(*Ra_D, 1.0 / 6.0) / denom;
        *Nu_D = (*Nu_D) * (*Nu_D);

        double h = (*Nu_D) * k_air / D_motor;
        T_guess = T_amb + P_loss / (h * A_surface);
    }

    *T_surface = T_guess;
    return 0;
}

int dc_motor_forced_cooling(double T_amb, double P_loss,
                            double D_motor, double L_motor,
                            double U_air,
                            double *T_surface, double *Re_D, double *Nu_D)
{
    /*
     * Forced air cooling over cylindrical motor housing.
     *
     * Uses: Tesla Model 3 oil-cooled + air-cooled motor,
     * industrial servo drives, wind turbine generators.
     */
    if (P_loss <= 0.0 || D_motor <= 0.0 || L_motor <= 0.0 ||
        U_air <= 0.0 || !T_surface || !Re_D || !Nu_D) return -1;

    /* Air properties at ~300K */
    double rho_air = 1.18;
    double mu_air = 1.85e-5;
    double k_air = 0.026;
    double Pr_air = 0.71;

    double A_surface = M_PI * D_motor * L_motor;

    *Re_D = rho_air * U_air * D_motor / mu_air;
    *Nu_D = nu_cylinder_crossflow(*Re_D, Pr_air);

    double h = (*Nu_D) * k_air / D_motor;
    *T_surface = T_amb + P_loss / (h * A_surface);

    return 0;
}

/* ==========================================================================
 * L7: Mars Rover Night Survival
 * ========================================================================== */

int mars_rover_night_survival(double P_heat, double A_surface,
                              double T_init, double t_night,
                              double C_thermal,
                              double *T_final, double *Q_lost)
{
    /*
     * Mars rover thermal survival model.
     *
     * Mars conditions:
     *   T_amb ~ -80 degC (193 K) at night (equatorial)
     *   P_atm ~ 600 Pa (CO2, 0.6% of Earth)
     *   g_mars = 3.72 m/s^2
     *
     * The key dimensionless numbers for this problem:
     *   - Ra: extremely small due to thin atmosphere
     *     -> Nu ≈ 1 (pure conduction to stagnant gas)
     *   - Bi << 1: internal conduction >> surface convection
     *     -> lumped capacitance model valid
     *
     * Lumped capacitance model:
     *   C_thermal . dT/dt = P_heat - h.A.(T - T_amb)
     *
     * With constant h and T_amb:
     *   T(t) = T_amb + P/(h.A) + [T_init - T_amb - P/(h.A)].exp(-h.A.t/C)
     *
     * Application: NASA Mars 2020 (Perseverance) rover uses
     * Radioisotope Heater Units (RHUs) producing ~1 W each,
     * with a total of ~5-10 to keep electronics above -40 degC.
     *
     * Note: The thin Martian atmosphere gives h ≈ 0.5-2 W/(m^2.K),
     * compared to ~5-25 W/(m^2.K) for natural convection on Earth.
     */
    if (!T_final || !Q_lost) return -1;
    if (A_surface <= 0.0 || C_thermal <= 0.0 || t_night < 0.0) return -1;

    /* Mars natural convection h (very low due to thin CO2) */
    /* For Bi << 1, h ~ 1-2 W/(m^2.K) on Mars */
    double h_mars = 1.5; /* W/(m^2.K) — empirical estimate */

    double T_amb = 193.0; /* Mars night ambient, K */
    double tau = C_thermal / (h_mars * A_surface); /* thermal time constant [s] */

    if (tau < 1e-10) {
        *T_final = T_amb + P_heat / (h_mars * A_surface);
        *Q_lost = h_mars * A_surface * (T_init - T_amb) * t_night;
        return 0;
    }

    /* Steady-state temperature with heating */
    double T_ss = T_amb + P_heat / (h_mars * A_surface);

    /* Exponential decay/growth toward T_ss */
    double exp_term = exp(-t_night / tau);
    *T_final = T_ss + (T_init - T_ss) * exp_term;

    /* Total heat lost = integrated (h.A.(T(t)-T_amb) - P) dt */
    /* = C_thermal . (T_init - T_final) + P.t_night  [from energy balance] */
    *Q_lost = C_thermal * (T_init - *T_final) + P_heat * t_night;

    return 0;
}

/* ==========================================================================
 * L6: Mixed Convection
 * ========================================================================== */

double mixed_convection_nusselt(double Nu_forced, double Nu_natural,
                                int aiding)
{
    /*
     * Combined forced + natural convection:
     *   Nu_mixed^n = Nu_forced^n ± Nu_natural^n
     *
     * n = 3 for vertical surfaces (Churchill, 1977)
     * n = 4 for horizontal surfaces
     * n = 7/2 for spheres
     *
     * + sign: aiding flow (buoyancy helps forced flow — upward heated,
     *         or downward cooled)
     * - sign: opposing flow (buoyancy opposes forced flow)
     */
    if (Nu_forced < 0.0 || Nu_natural < 0.0) return -1.0;

    double n = 3.0; /* vertical surfaces */
    double sign = aiding ? 1.0 : -1.0;

    double combined = pow(Nu_forced, n) + sign * pow(Nu_natural, n);
    if (combined < 0.0) combined = pow(Nu_forced, n); /* can't be negative */

    return pow(combined, 1.0 / n);
}

int convection_regime_map(double Gr, double Re)
{
    /*
     *   Gr/Re^2 < 0.01 -> forced convection
     *   0.01 < Gr/Re^2 < 10 -> mixed
     *   Gr/Re^2 > 10 -> natural convection
     */
    if (Gr < 0.0 || Re < 0.0) return -1;
    if (Re < 1e-15) return (Gr > 0) ? 3 : 0;
    if (Gr < 1e-15) return 0;

    double ratio = Gr / (Re * Re);

    if (ratio < 0.01) return 0; /* forced */
    if (ratio < 10.0) {
        /* Mixed — direction cannot be determined without DeltaT.U sign */
        return 1; /* assume aiding as default */
    }
    return 3; /* natural */
}

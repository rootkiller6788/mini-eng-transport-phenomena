/**
 * interphase_applications.c
 * ==========================
 * Engineering applications of interphase transport theory.
 *
 * Reference: Bird, Stewart, Lightfoot (2007) "Transport Phenomena"
 * Treybal, R.E. (1980) "Mass-Transfer Operations" (3rd ed.)
 * Incropera & DeWitt (2007) "Fundamentals of Heat and Mass Transfer"
 *
 * Knowledge coverage: L6 Engineering Problems, L7 Applications
 */

#include "interphase_types.h"
#include "interphase_boundary_conditions.h"
#include "interphase_transport.h"
#include "interphase_jump_conditions.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

/* ============================================================================
 * L6: Falling Film Gas Absorption
 *
 * Problem: CO2 absorption from air into a falling water film.
 * This is a classic interphase mass transfer problem in chemical
 * engineering, used in packed columns for gas scrubbing.
 *
 * Reference: Danckwerts (1970) "Gas-Liquid Reactions"
 *   Chapter 5 - Absorption with Chemical Reaction
 * ============================================================================ */

/**
 * Compute CO2 absorption rate in a falling film absorber.
 *
 * Given: Liquid flow rate, film geometry, inlet concentrations,
 *        Henry's constant for CO2 in water, mass transfer coefficients.
 *
 * Returns: Absorption flux (mol/(m^2*s)) and exit concentration.
 *
 * @param Q_L          Liquid volumetric flow rate (m^3/s)
 * @param L            Film length (m)
 * @param W            Film width (m)
 * @param p_CO2        Partial pressure of CO2 in gas (Pa)
 * @param H_CO2        Henry's constant for CO2-H2O (Pa), ~1.67e8 at 298K
 * @param C_in         Inlet dissolved CO2 concentration (mol/m^3)
 * @param k_L          Liquid-side mass transfer coeff. (m/s)
 * @param D_CO2        CO2 diffusivity in water (m^2/s), ~1.92e-9 at 298K
 * @param rho_L        Liquid density (kg/m^3)
 * @param mu_L         Liquid viscosity (Pa*s)
 * @param N_abs        Output total absorption rate (mol/s)
 * @param C_out        Output exit concentration (mol/m^3)
 * Complexity: O(1)
 *
 * Theory:
 *   Film thickness: delta = (3*mu*Q/(rho*g*W))^(1/3) (Nusselt film)
 *   Surface velocity: u_s = rho*g*delta^2/(2*mu)
 *   Contact time: t_c = L / u_s
 *   Penetration theory: k_L = 2*sqrt(D/(pi*t_c))
 *   Concentration profile: C_out = C* + (C_in - C*) * exp(-k_L*a*L/u_s)
 *     where C* = p_CO2/H_CO2 (equilibrium interfacial concentration)
 */
void falling_film_absorption(double Q_L, double L, double W,
                              double p_CO2, double H_CO2,
                              double C_in, double k_L, double D_CO2,
                              double rho_L, double mu_L,
                              double *N_abs, double *C_out) {
    if (!N_abs || !C_out) return;

    double g = 9.81;

    /* Film thickness from Nusselt solution */
    double delta_film = pow(3.0 * mu_L * Q_L / (rho_L * g * W), 1.0/3.0);

    /* Surface velocity (maximum film velocity) */
    double u_surf = rho_L * g * delta_film * delta_film / (2.0 * mu_L);

    /* Exposure time = L/u_surf determines mass transfer regime */
    (void)(L / u_surf);  /* t_contact reference */

    /* Equilibrium interfacial concentration from Henry's law */
    double C_star = p_CO2 / H_CO2;

    /* Specific interfacial area (planar film: a = W / cross_section_area) */
    double a = 1.0 / delta_film;  /* 1/m, interface area per liquid volume */

    /* Mass balance for plug flow with mass transfer:
     * u_s * dC/dx = k_L * a * (C* - C)
     * Solution: C(x) = C* + (C_in - C*) * exp(-k_L*a*x/u_s)
     */
    double exponent = -k_L * a * L / u_surf;
    if (exponent < -50.0) exponent = -50.0;  /* prevent underflow */
    *C_out = C_star + (C_in - C_star) * exp(exponent);

    /* Total absorption rate */
    *N_abs = Q_L * (*C_out - C_in);

    /* N_abs > 0 means absorption (gas -> liquid) */
}

/* ============================================================================
 * L6: Evaporative Cooling at a Liquid-Vapor Interface
 *
 * Problem: Water evaporates from a wet surface into a flowing air stream.
 * This is the principle behind cooling towers, sweat cooling, and
 * wet-bulb thermometers.
 *
 * The interface temperature drops due to latent heat removal, reaching
 * the wet-bulb temperature when steady state is achieved.
 *
 * Reference: ASHRAE Handbook (2017) "Fundamentals" Chapter 1 -
 *   Psychrometrics
 * ============================================================================ */

/**
 * Compute evaporative cooling rate and wet-bulb temperature.
 *
 * Energy balance at interface:
 *   h * (T_air - T_interface) = k_m * (C_interface - C_bulk) * h_fg
 * (convective heat transfer = latent heat of evaporation)
 *
 * @param T_air        Dry-bulb air temperature (K)
 * @param RH           Relative humidity (0 to 1)
 * @param p_atm        Atmospheric pressure (Pa)
 * @param h            Heat transfer coefficient (W/(m^2*K))
 * @param k_m          Mass transfer coefficient (m/s)
 * @param h_fg         Latent heat of vaporization at T_interface (J/kg)
 * @param cp_air       Specific heat of air (J/(kg*K))
 * @param Le           Lewis number alpha/D (-)
 * @param T_wet_bulb   Output wet-bulb temperature (K)
 * @param evap_flux    Output evaporation mass flux (kg/(m^2*s))
 * Complexity: O(n_iter) ~ O(20)
 *
 * Theorem (Psychrometric ratio): h/k_m = rho*cp*Le^(2/3)
 *   via Chilton-Colburn analogy.
 * Theorem: T_wet_bulb <= T_air, with equality only at RH = 1.
 * Theorem: Maximum evaporative cooling T_air - T_wet_bulb occurs at RH = 0.
 */
void evaporative_cooling(double T_air, double RH, double p_atm,
                          double h, double k_m, double h_fg,
                          double cp_air, double Le,
                          double *T_wet_bulb, double *evap_flux) {
    (void)cp_air;  /* used via h/k_m ratio per Chilton-Colburn analogy */
    (void)Le;      /* determines relative rates of heat/mass transfer */
    if (!T_wet_bulb || !evap_flux) return;
    if (T_air <= 0.0 || p_atm <= 0.0) {
        *T_wet_bulb = T_air;
        *evap_flux = 0.0;
        return;
    }

    /* Saturation vapor pressure at T using Antoine equation for water */
    /* log10(P_sat [Pa]) = 10.227 - 1730.63/(T - 39.72) [T in K] */
    double p_sat_air = pow(10.0, 10.227 - 1730.63 / (T_air - 39.72));

    /* Humidity ratio in bulk air */
    double p_v_air = RH * p_sat_air;
    double omega_air = 0.622 * p_v_air / (p_atm - p_v_air);

    /* Iterate to find wet-bulb temperature (simple fixed-point iteration) */
    double T_wb = T_air - 5.0;  /* initial guess */
    double rho_air = p_atm / (287.0 * T_air);

    for (int iter = 0; iter < 50; iter++) {
        /* Saturation pressure at T_wb */
        double p_sat_wb = pow(10.0, 10.227 - 1730.63 / (T_wb - 39.72));
        double omega_sat_wb = 0.622 * p_sat_wb / (p_atm - p_sat_wb);

        /* Energy balance: h*(T_air - T_wb) = k_m*rho*(omega_sat - omega)*h_fg */
        double lhs = h * (T_air - T_wb);
        double rhs = k_m * rho_air * (omega_sat_wb - omega_air) * h_fg;

        double dT = 0.01 * (T_air - T_wb);
        if (fabs(lhs - rhs) < 1e-3 * fmax(lhs, 1.0)) break;

        if (lhs > rhs) {
            T_wb -= dT;  /* too much cooling predicted, reduce T_wb */
        } else {
            T_wb += dT;
        }
        if (T_wb < 273.15) T_wb = 273.15;
        if (T_wb >= T_air) T_wb = T_air - 0.1;
    }

    *T_wet_bulb = T_wb;

    /* Evaporation flux */
    double p_sat_final = pow(10.0, 10.227 - 1730.63 / (T_wb - 39.72));
    double omega_sat_final = 0.622 * p_sat_final / (p_atm - p_sat_final);
    *evap_flux = k_m * rho_air * (omega_sat_final - omega_air);
    if (*evap_flux < 0.0) *evap_flux = 0.0;
}

/* ============================================================================
 * L6: Distillation Tray Efficiency
 *
 * Problem: Compute Murphree vapor efficiency for a distillation tray.
 *
 * The Murphree efficiency (Murphree, 1925) characterizes deviation
 * from equilibrium on a distillation tray:
 *   E_MV = (y_n - y_{n+1}) / (y*_n - y_{n+1})
 * where y*_n is vapor composition in equilibrium with liquid leaving tray n.
 *
 * Point efficiency from mass transfer:
 *   E_OG = 1 - exp(-N_OG)
 *   N_OG = number of overall gas-phase transfer units
 *
 * Reference: Kister (1992) "Distillation Design" (McGraw-Hill)
 *   Chapter 6 - Tray Efficiency
 * ============================================================================ */

/**
 * Compute Murphree vapor tray efficiency.
 *
 * @param k_g          Gas-phase mass transfer coeff. (kmol/(m^2*s))
 * @param k_l          Liquid-phase mass transfer coeff. (kmol/(m^2*s))
 * @param m            Equilibrium line slope dy_star/dx (dimensionless)
 * @param a            Specific interfacial area (m^2/m^3)
 * @param h_f          Froth height on tray (m)
 * @param G            Gas molar velocity (kmol/(m^2*s))
 * @param L            Liquid molar velocity (kmol/(m^2*s))
 * @param lambda_stripping Stripping factor m*G/L (-)
 * @param E_MV         Output Murphree vapor efficiency (0 to 1)
 * Complexity: O(1)
 *
 * Theorem: E_MV -> 1 as N_OG -> infinity (equilibrium stage).
 * Theorem: E_MV decreases with increasing G (shorter contact time).
 * Theorem: For lambda = 1, E_MV = E_OG.
 * Theorem: For lambda != 1, E_MV = (exp(lambda*E_OG) - 1) / lambda
 *          (Lewis, 1936 - fully mixed liquid).
 */
void murphree_tray_efficiency(double k_g, double k_l, double m,
                               double a, double h_f,
                               double G, double L,
                               double *E_MV) {
    if (!E_MV || G <= 0.0) {
        if (E_MV) *E_MV = 0.0;
        return;
    }

    /* Overall gas-phase mass transfer coefficient */
    /* 1/K_OG = 1/k_g + m/k_l */
    double inv_K_OG = 1.0 / k_g + m / k_l;
    double K_OG = 1.0 / inv_K_OG;

    /* Number of overall gas-phase transfer units */
    double N_OG = K_OG * a * h_f / G;

    /* Point efficiency */
    double E_OG = 1.0 - exp(-N_OG);

    /* Murphree vapor efficiency for fully mixed liquid */
    double lambda = (G > 0.0) ? (m * G / L) : 1.0;

    if (fabs(lambda - 1.0) < 1e-10) {
        *E_MV = E_OG / (1.0 + E_OG);  /* Lewis case 1: lambda=1 */
    } else {
        /* E_MV = (exp(lambda*E_OG) - 1) / lambda */
        double exp_term = exp(lambda * E_OG);
        *E_MV = (exp_term - 1.0) / lambda;
    }

    if (*E_MV > 1.0) *E_MV = 1.0;
    if (*E_MV < 0.0) *E_MV = 0.0;
}

/* ============================================================================
 * L6: Membrane Separator Design
 *
 * Problem: Sizing a hollow-fiber membrane contactor for gas separation.
 *
 * Membrane contactors use microporous hollow fibers to provide high
 * interfacial area per unit volume (typically 1000-10000 m^2/m^3),
 * orders of magnitude higher than conventional packed columns (100-300).
 *
 * Reference: Gabelman & Hwang (1999) "Hollow Fiber Membrane Contactors"
 *   J. Membrane Science, 159(1-2), 61-106.
 * ============================================================================ */

/**
 * Size a hollow-fiber membrane contactor for gas absorption.
 *
 * @param Q_L          Liquid flow rate (m^3/s)
 * @param C_in         Inlet concentration (mol/m^3)
 * @param C_out_target Target outlet concentration (mol/m^3)
 * @param C_star       Equilibrium interfacial concentration (mol/m^3)
 * @param K_overall    Overall mass transfer coefficient (m/s)
 * @param d_fiber      Fiber outer diameter (m)
 * @param packing_frac Packing fraction (fiber area / module area) (-)
 * @param L_module     Output required module length (m)
 * @param A_total      Output total membrane area (m^2)
 * Complexity: O(1)
 *
 * Design equation: L = (Q_L / (K_overall * a * A_cs)) * ln((C*-C_in)/(C*-C_out))
 * where a = 4*epsilon/d_fiber (specific area)
 */
void membrane_contactor_sizing(double Q_L, double C_in, double C_out_target,
                                double C_star, double K_overall,
                                double d_fiber, double packing_frac,
                                double *L_module, double *A_total) {
    if (!L_module || !A_total) return;
    if (Q_L <= 0.0 || K_overall <= 0.0 || d_fiber <= 0.0) {
        if (L_module) *L_module = 0.0;
        if (A_total)  *A_total  = 0.0;
        return;
    }

    /* Specific interfacial area: a = 4*epsilon/d_fiber (cylindrical fibers) */
    double a_specific = 4.0 * packing_frac / d_fiber;

    /* Cross-sectional area: arbitrary, choose for reasonable velocity */
    /* Target liquid velocity ~ 0.01-0.1 m/s */
    double v_target = 0.05;  /* m/s */
    double A_cs = Q_L / v_target;
    (void)sqrt(4.0 * A_cs / M_PI);  /* D_module for reference */

    /* Number of transfer units required */
    double driving_force_ratio = (C_star - C_in) / (C_star - C_out_target);
    if (driving_force_ratio <= 0.0) {
        *L_module = 0.0;
        *A_total  = 0.0;
        return;
    }
    double NTU = log(driving_force_ratio);

    /* Module length from design equation:
     * NTU = K_overall * a_specific * L / v
     * => L = NTU * v / (K_overall * a_specific)
     */
    double v_superficial = Q_L / A_cs;
    *L_module = NTU * v_superficial / (K_overall * a_specific);

    /* Total membrane area */
    *A_total = a_specific * A_cs * (*L_module);
}

/* ============================================================================
 * L6: Condensation Heat Transfer on a Horizontal Tube
 *
 * Problem: Steam condenses on the outside of a horizontal cooling tube.
 * This is the dominant heat transfer mechanism in shell-and-tube
 * condensers used in power plants (Rankine cycle) and refrigeration.
 *
 * Reference: Nusselt (1916) "Die Oberflaechenkondensation des Wasserdampfes"
 *   VDI Zeitschrift, 60, 541-546, 569-575.
 * ============================================================================ */

/**
 * Compute condensation heat transfer for a horizontal tube bank.
 *
 * Nusselt's film condensation theory for a single horizontal tube:
 *   h_avg = 0.725 * [g*rho_L*(rho_L-rho_v)*k_L^3*h_fg/(mu_L*D*delta_T)]^(1/4)
 *
 * For N tubes in a vertical column, h reduces by N^(-1/4) due to
 * condensate inundation (Kern, 1958).
 *
 * @param D_tube      Tube outer diameter (m)
 * @param T_sat       Saturation temperature (K)
 * @param T_wall      Tube wall temperature (K)
 * @param rho_L       Liquid density (kg/m^3)
 * @param rho_v       Vapor density (kg/m^3)
 * @param k_L         Liquid thermal conductivity (W/(m*K))
 * @param h_fg        Latent heat (J/kg)
 * @param mu_L        Liquid viscosity (Pa*s)
 * @param N_tubes     Number of tubes in vertical column
 * @param h_avg       Output average heat transfer coefficient (W/(m^2*K))
 * @param Q_cond      Output condensation heat rate per tube (W/m)
 * Complexity: O(1)
 *
 * Theorem: h ~ delta_T^(-1/4) (decreases with larger subcooling).
 * Theorem: h ~ D^(-1/4) (smaller tubes give higher h).
 * Theorem: h_N = h_1 * N^(-1/4) (inundation effect).
 */
void condensation_horizontal_tube(double D_tube, double T_sat, double T_wall,
                                   double rho_L, double rho_v,
                                   double k_L, double h_fg, double mu_L,
                                   int N_tubes,
                                   double *h_avg, double *Q_cond) {
    if (!h_avg || !Q_cond) return;

    double delta_T = T_sat - T_wall;
    if (delta_T <= 0.0 || D_tube <= 0.0 || mu_L <= 0.0 || k_L <= 0.0) {
        *h_avg = 0.0;
        *Q_cond = 0.0;
        return;
    }

    double g = 9.81;
    double drho = rho_L - rho_v;
    if (drho <= 0.0) { *h_avg = 0.0; *Q_cond = 0.0; return; }

    /* Single tube Nusselt coefficient */
    double h_1 = 0.725 * pow(g * drho * k_L*k_L*k_L * h_fg
                              / (mu_L * D_tube * delta_T), 0.25);

    /* Inundation correction for N tubes */
    double N_factor = pow((double)N_tubes, -0.25);
    *h_avg = h_1 * N_factor;

    /* Heat rate per unit tube length */
    *Q_cond = *h_avg * M_PI * D_tube * delta_T;
}

/* ============================================================================
 * L7: Bubble Column Reactor Interfacial Area
 *
 * Application: Wastewater treatment using aerobic bacteria in bubble
 * columns. Oxygen transfer from air bubbles to liquid phase is the
 * rate-limiting step. Design requires accurate prediction of
 * interfacial area and mass transfer rate.
 *
 * Reference: Tchobanoglous, Burton, Stensel (2003) "Wastewater
 *   Engineering: Treatment and Reuse" (4th ed., Metcalf & Eddy)
 *   Chapter 5 - Biological Unit Processes
 * ============================================================================ */

/**
 * Design oxygen transfer in a bubble column for wastewater treatment.
 *
 * @param Q_wastewater Flow rate (m^3/day)
 * @param BOD_in      Inlet BOD (mg/L), typically 200-300 for municipal
 * @param BOD_out     Target effluent BOD (mg/L), typically < 30
 * @param T           Operating temperature (K)
 * @param alpha_factor Oxygen transfer correction for wastewater (~0.6-0.9)
 * @param beta_factor  Salinity/solids correction (~0.95-1.0)
 * @param V_tank       Output required tank volume (m^3)
 * @param Q_air        Output required air flow rate (m^3/s)
 * @param P_blower     Output estimated blower power (kW)
 * Complexity: O(1)
 *
 * Design approach (ASCE standard):
 *   SOTR = alpha * F * beta * (C*_s - C_L) / C*_s_20 * SOTR_std * theta^(T-20)
 *   where SOTR = standard oxygen transfer rate
 */
void wastewater_aeration_design(double Q_wastewater, double BOD_in,
                                 double BOD_out, double T,
                                 double alpha_factor, double beta_factor,
                                 double *V_tank, double *Q_air,
                                 double *P_blower) {
    if (!V_tank || !Q_air || !P_blower) return;

    /* Typical design parameters for activated sludge */
    double theta = 1.024;  /* temperature correction factor */
    double C_star_20 = 9.09;  /* mg/L, DO saturation at 20C, 1 atm */
    double C_L_target = 2.0;  /* mg/L, operating DO level */
    double F = 0.9;  /* fouling factor for fine-bubble diffusers */

    /* Temperature correction for DO saturation
     * C*(T) = 14.652 - 0.41022*T + 0.007991*T^2 - 0.000077774*T^3
     * (APHA standard, T in Celsius)
     */
    double T_celsius = T - 273.15;
    double C_star_T = 14.652 - 0.41022 * T_celsius
                      + 0.007991 * T_celsius * T_celsius
                      - 0.000077774 * T_celsius * T_celsius * T_celsius;

    /* Actual oxygen transfer rate needed:
     * O2 required ~ 1.5 kg O2 / kg BOD removed (typical for extended aeration)
     */
    double BOD_removed = BOD_in - BOD_out;  /* mg/L */
    double O2_required = 1.5 * Q_wastewater * BOD_removed / 1000.0;  /* kg/day */

    /* Standard oxygen transfer rate (SOTR) accounting for process conditions */
    double temp_corr = pow(theta, T_celsius - 20.0);
    double driving_force_corr = (beta_factor * C_star_T - C_L_target)
                                / C_star_20;
    double SOTR = O2_required
                  / (alpha_factor * F * driving_force_corr * temp_corr);

    /* Tank volume: assume typical loading ~ 0.5 kg BOD/(m^3*day) */
    double vol_loading = 0.5;  /* kg BOD/(m^3*day) */
    *V_tank = (Q_wastewater * BOD_in / 1000.0) / vol_loading;

    /* Air flow rate:
     * O2 in air = 0.232 kg O2/kg air, air density ~1.2 kg/m^3
     * SOTE ~ 25% for fine-bubble diffusers at 4m depth
     */
    double SOTE = 0.25;  /* standard oxygen transfer efficiency */
    double O2_per_m3_air = 1.2 * 0.232;  /* kg O2 / m^3 air */
    *Q_air = SOTR / (24.0 * 3600.0 * O2_per_m3_air * SOTE);  /* m^3/s */

    /* Blower power:
     * P = Q * delta_P / eta
     * delta_P ~ rho*g*h ~ 1000*9.81*4 ~ 40 kPa (4m submergence)
     */
    double delta_P = 40000.0;  /* Pa */
    double eta_blower = 0.70;
    *P_blower = *Q_air * delta_P / eta_blower / 1000.0;  /* kW */
}

/* ============================================================================
 * L7: Condensation in Nuclear Power Plant Steam Generator
 *
 * Application: Steam generator in a PWR nuclear power plant.
 * Primary loop (radioactive, ~15.5 MPa, ~315C) transfers heat to
 * secondary loop through Inconel tubes. Secondary side boils at
 * ~6.9 MPa, ~285C. The steam generator is a critical component
 * connecting the nuclear island to the turbine.
 *
 * Reference: Todreas & Kazimi (2012) "Nuclear Systems" (2nd ed.)
 *   Volume I - Thermal Hydraulic Fundamentals
 *   Chapter 10 - Single Heated Channel: Steady State
 *
 * Key data from real plants:
 * - Westinghouse AP1000: 2 steam generators, each ~ 700 MWe
 * - U-tube design, ~10,000 tubes per SG, 19mm OD, ~20m tall
 * - Recirculating type: ~ 3-5 circulation ratio
 * ============================================================================ */

/**
 * Simplified steam generator thermal sizing for a PWR.
 *
 * @param P_thermal    Reactor thermal power (MW_th)
 * @param T_primary_in Primary coolant inlet temperature (K), ~ 588K (315C)
 * @param T_primary_out Primary coolant outlet temperature (K), ~ 560K (287C)
 * @param T_sat_secondary Secondary side saturation temperature (K), ~ 558K (285C)
 * @param cp_primary   Primary coolant specific heat (J/(kg*K)), ~5000
 * @param U_overall    Overall HTC (W/(m^2*K)), ~5000-8000 for PWR SG
 * @param D_tube       Tube outer diameter (m), typically 0.019
 * @param L_tube       Average tube length (m), typically 20
 * @param N_tubes      Output number of tubes required
 * @param A_ht         Output total heat transfer area (m^2)
 * @param m_dot_primary Output primary mass flow rate (kg/s)
 * Complexity: O(1)
 */
void pwr_steam_generator_sizing(double P_thermal,
                                 double T_primary_in, double T_primary_out,
                                 double T_sat_secondary,
                                 double cp_primary, double U_overall,
                                 double D_tube, double L_tube,
                                 double *N_tubes, double *A_ht,
                                 double *m_dot_primary) {
    if (!N_tubes || !A_ht || !m_dot_primary) return;
    if (U_overall <= 0.0 || cp_primary <= 0.0 || D_tube <= 0.0
        || L_tube <= 0.0) {
        *N_tubes = 0; *A_ht = 0.0; *m_dot_primary = 0.0;
        return;
    }

    double P_W = P_thermal * 1.0e6;  /* convert MW_th to W */

    /* Primary mass flow rate (energy balance) */
    *m_dot_primary = P_W / (cp_primary * (T_primary_in - T_primary_out));

    /* LMTD for counterflow with secondary at constant T_sat */
    double dT1 = T_primary_in - T_sat_secondary;
    double dT2 = T_primary_out - T_sat_secondary;
    if (dT1 <= 0.0 || dT2 <= 0.0) {
        *N_tubes = 0; *A_ht = 0.0;
        return;
    }
    double LMTD = (dT1 - dT2) / log(dT1 / dT2);

    /* Required heat transfer area */
    *A_ht = P_W / (U_overall * LMTD);

    /* Number of tubes */
    double A_per_tube = M_PI * D_tube * L_tube;
    *N_tubes = (int) ceil(*A_ht / A_per_tube);
}

/* ============================================================================
 * L7: Interfacial Transport in Fuel Cells
 *
 * Application: Proton Exchange Membrane (PEM) fuel cell water management.
 * Water transport across the membrane and gas diffusion layers involves
 * multi-phase transport with electro-osmotic drag and back diffusion.
 *
 * Reference: Springer, Zawodzinski, Gottesfeld (1991)
 *   "Polymer Electrolyte Fuel Cell Model"
 *   J. Electrochemical Society, 138(8), 2334-2342.
 *
 * Key data (from Toyota Mirai, Honda Clarity, Hyundai Nexo):
 * - Operating temperature: 60-80C
 * - Membrane: Nafion, ~25-50 um thick
 * - Water content lambda = 2-22 H2O/SO3H
 * - Electro-osmotic drag coefficient n_d ~ 0.5-1.5 H2O/H+
 * - Back diffusion opposes drag at high current density
 * ============================================================================ */

/**
 * Compute water balance in a PEM fuel cell membrane.
 *
 * Net water flux: J_w = n_d * i/F - D_lambda * dc/dz (Springer model)
 * where:
 *   n_d = electro-osmotic drag coefficient (mol H2O/mol H+)
 *   i = current density (A/cm^2)
 *   F = Faraday constant 96485 C/mol
 *   D_lambda = water diffusion coefficient in membrane (cm^2/s)
 *   dc/dz = water concentration gradient (mol/cm^4)
 *
 * @param i_current      Current density (A/cm^2)
 * @param T_cell         Cell temperature (K)
 * @param lambda_anode   Water content at anode side (mol H2O/mol SO3H)
 * @param lambda_cathode Water content at cathode side
 * @param t_membrane     Membrane thickness (cm)
 * @param rho_dry        Dry membrane density (g/cm^3), Nafion ~2.0
 * @param EW             Equivalent weight (g/mol SO3H), Nafion ~1100
 * @param J_net          Output net water flux (mol/(cm^2*s)), + toward cathode
 * @param drag_flux      Output electro-osmotic drag component
 * @param diff_flux      Output back diffusion component
 * Complexity: O(1)
 *
 * Theorem: At equilibrium (no current), drag = 0 and diffusion = 0.
 * Theorem: At high current, drag dominates and cathode floods.
 * Theorem: At low current / high T, back diffusion dominates and anode dries.
 */
void pem_fuel_cell_water_flux(double i_current, double T_cell,
                               double lambda_anode, double lambda_cathode,
                               double t_membrane, double rho_dry,
                               double EW,
                               double *J_net, double *drag_flux,
                               double *diff_flux) {
    if (!J_net || !drag_flux || !diff_flux) return;

    double F = 96485.0;  /* Faraday constant C/mol */
    if (t_membrane <= 0.0 || EW <= 0.0) {
        *J_net = 0.0; *drag_flux = 0.0; *diff_flux = 0.0;
        return;
    }

    /* Electro-osmotic drag coefficient (Springer correlation) */
    double n_d = 2.5 * lambda_anode / 22.0;  /* simplified linear model */
    if (n_d < 0.0) n_d = 0.0;

    /* Drag flux: n_d * i / F [mol/(cm^2*s)] */
    *drag_flux = n_d * i_current / F;

    /* Water diffusion coefficient (Springer correlation) */
    double D_lambda = 1.0e-6 * exp(2416.0 * (1.0/303.0 - 1.0/T_cell))
                      * (0.0001 * lambda_anode * lambda_anode
                         + 0.001 * lambda_anode + 0.01);  /* cm^2/s */

    /* Concentration gradient: dc/dz = rho_dry/EW * d(lambda)/dz */
    double c_factor = rho_dry / EW;  /* mol SO3H / cm^3 */
    double dlambda_dz = (lambda_cathode - lambda_anode) / t_membrane;
    double dc_dz = c_factor * dlambda_dz;

    /* Back diffusion flux: -D_lambda * dc/dz [mol/(cm^2*s)] */
    *diff_flux = -D_lambda * dc_dz;

    /* Net flux */
    *J_net = *drag_flux + *diff_flux;
}

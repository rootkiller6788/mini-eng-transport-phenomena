/*
 * nusselt_number.h — Nusselt Number and Convective Heat Transfer Correlations
 *
 * Reference: Incropera & DeWitt (2007) "Fundamentals of Heat and Mass Transfer"
 *            Bejan (2013) "Convection Heat Transfer"
 *            Kays, Crawford, Weigand (2005) "Convective Heat and Mass Transfer"
 *
 * The Nusselt number Nu = hL/k is the dimensionless temperature gradient
 * at the wall. It represents the ratio of convective to conductive heat transfer.
 *
 * Knowledge Coverage: L1 Definitions, L2 Heat Transfer Concepts,
 *                     L3 Correlations by Geometry, L4 Energy Conservation,
 *                     L5 Engineering Methods (empirical Nu correlations),
 *                     L6 Canonical Problems (pipe, flat plate, cylinder, sphere)
 */

#ifndef NUSSELT_NUMBER_H
#define NUSSELT_NUMBER_H

#ifdef __cplusplus
extern "C" {
#endif

/* ==========================================================================
 * L1: Core Nu Definition
 * ========================================================================== */

/**
 * nu_basic — compute Nu from heat transfer coefficient.
 *
 * Nu = h.L/k
 *
 * @param h    convective heat transfer coefficient [W/(m^2.K)]
 * @param L    characteristic length [m]
 * @param k    fluid thermal conductivity [W/(m.K)]
 * @return     Nusselt number [-]
 */
double nu_basic(double h, double L, double k);

/**
 * h_from_nusselt — recover heat transfer coefficient from Nu.
 *
 * h = Nu.k/L
 */
double h_from_nusselt(double Nu, double L, double k);

/**
 * nu_conductivity — log-mean temperature interpretation.
 *
 * Nu = (q".L)/(k.DeltaT)
 *
 * @param q_flux   wall heat flux [W/m^2]
 * @param L        characteristic length [m]
 * @param k        thermal conductivity [W/(m.K)]
 * @param delta_T  wall-to-bulk temperature difference [K]
 * @return         Nusselt number [-]
 */
double nu_from_heat_flux(double q_flux, double L, double k, double delta_T);

/* ==========================================================================
 * L3-L5: Nu Correlations — Internal Flow (Pipes & Ducts)
 * ========================================================================== */

/**
 * nu_pipe_laminar_constant_T — fully developed laminar pipe flow, const T wall.
 *
 * @return  Nu_D = 3.66
 *
 * Theorem: For hydrodynamically and thermally fully developed laminar flow
 *   in a circular tube with constant wall temperature, the Nusselt number
 *   is a constant independent of Re and Pr.
 *
 *   Nu_D = 3.66 (T_w = constant), Nu_D = 4.36 (q"_w = constant).
 *
 *   Valid for: Re_D < 2300, Gz⁻¹ = (x/D)/(Re.Pr) > 0.05.
 */
double nu_pipe_laminar_constant_T(void);

/**
 * nu_pipe_laminar_constant_qflux — fully developed laminar, const heat flux.
 *
 * @return  Nu_D = 4.36
 */
double nu_pipe_laminar_constant_qflux(void);

/**
 * nu_pipe_laminar_developing — thermal entrance region, laminar.
 *
 * @param Gz_inv  inverse Graetz number = (x/D)/(Re.Pr)
 * @return        Nu_D [-]
 *
 * For thermally developing flow with fully developed velocity profile
 * and constant wall temperature (Sieder-Tate, 1936):
 *
 *   Nu_D = 1.86.(Re.Pr.D/L)^(1/3).(mu/mu_w)^0.14
 *
 * Valid for: 0.48 < Pr < 16700, (Re.Pr.D/L)^(1/3).(mu/mu_w)^0.14 > 2.
 */
double nu_pipe_laminar_developing(double Re, double Pr, double D, double L,
                                  double mu, double mu_w);

/**
 * nu_pipe_dittus_boelter — turbulent pipe flow, Dittus-Boelter (1930).
 *
 * Nu_D = 0.023.Re_D^0.8.Pr^n
 *   n = 0.4 for heating (T_w > T_b)
 *   n = 0.3 for cooling (T_w < T_b)
 *
 * @param Re      Reynolds number (Re_D > 10000)
 * @param Pr      Prandtl number (0.7 < Pr < 160)
 * @param heating 1 = heating, 0 = cooling
 * @return        Nusselt number [-]
 *
 * Valid for: Re > 10⁴, 0.7 < Pr < 160, L/D > 10.
 * Simple and widely used but ±25% error.
 */
double nu_pipe_dittus_boelter(double Re, double Pr, int heating);

/**
 * nu_pipe_sieder_tate — Sieder-Tate correlation (1936), includes viscosity ratio.
 *
 * Nu_D = 0.027.Re_D^0.8.Pr^(1/3).(mu/mu_w)^0.14
 *
 * @param Re    Reynolds number (Re_D > 10000)
 * @param Pr    Prandtl number (0.7 < Pr < 16700)
 * @param mu    bulk dynamic viscosity [Pa.s]
 * @param mu_w  wall-temperature dynamic viscosity [Pa.s]
 * @return      Nusselt number [-]
 *
 * The (mu/mu_w)^0.14 factor accounts for property variations across the
 * thermal boundary layer. Use mu_w evaluated at wall temperature.
 *
 * Valid for: Re > 10⁴, 0.7 < Pr < 16700, L/D > 10.
 */
double nu_pipe_sieder_tate(double Re, double Pr, double mu, double mu_w);

/**
 * nu_pipe_gnielinski — Gnielinski correlation (1976), most accurate.
 *
 * Nu_D = (f/8).(Re-1000).Pr / [1 + 12.7.(f/8)^0.5.(Pr^(2/3) - 1)]
 *
 * @param Re    Reynolds number (3000 < Re < 5×10⁶)
 * @param Pr    Prandtl number (0.5 < Pr < 2000)
 * @param f     Darcy friction factor for smooth pipe
 * @return      Nusselt number [-]
 *
 * Valid for: 3000 < Re < 5×10⁶, 0.5 < Pr < 2000.
 * Error: ±10% for most fluids. This is the recommended correlation
 *   in modern heat transfer textbooks.
 *
 * Reference: Gnielinski (1976), Int. Chem. Eng.
 */
double nu_pipe_gnielinski(double Re, double Pr, double f);

/**
 * nu_annulus — concentric annular duct.
 *
 * Nu_Dh = Nu_ii for inner wall heating, Nu_oo for outer wall.
 *
 * @param Re          Reynolds number (based on D_h)
 * @param Pr          Prandtl number
 * @param Di          inner tube diameter [m]
 * @param Do          outer tube diameter [m]
 * @param inner_heated 1 if inner wall heated, 0 for outer
 * @return            Nusselt number [-]
 *
 * D_h = Do - Di for annulus.
 * Uses modified Dittus-Boelter with diameter ratio correction.
 */
double nu_annulus(double Re, double Pr, double Di, double Do, int inner_heated);

/**
 * nu_noncircular_duct — laminar Nu for non-circular ducts.
 *
 * Fully developed laminar, constant T_w:
 *   Triangular: Nu = 2.47
 *   Square:     Nu = 2.98
 *   Hexagonal:  Nu = 3.66 (approaches circle)
 *   Parallel plates (b >> a): Nu = 7.54
 *
 * @param geometry_id  0 = circular (3.66), 1 = square, 2 = triangular,
 *                     3 = parallel plates, 4 = rectangular 2:1, 5 = rectangular 4:1
 * @return             Nusselt number [-]
 */
double nu_noncircular_duct_laminar(int geometry_id);

/* ==========================================================================
 * L3-L5: Nu Correlations — External Flow (Flat Plate, Cylinder, Sphere)
 * ========================================================================== */

/**
 * nu_flat_plate_laminar_local — local Nu for laminar BL over isothermal plate.
 *
 * Nu_x = 0.332.Re_x^(1/2).Pr^(1/3)
 *
 * @param Re_x  local Reynolds number at x
 * @param Pr    Prandtl number (Pr > 0.6)
 * @return      local Nu_x [-]
 *
 * Valid for: Re_x < 5×10⁵, Pr > 0.6, no viscous dissipation.
 * This is the exact Blasius/Pohlhausen similarity solution.
 */
double nu_flat_plate_laminar_local(double Re_x, double Pr);

/**
 * nu_flat_plate_laminar_avg — average Nu over plate length L.
 *
 * Nu_L = 0.664.Re_L^(1/2).Pr^(1/3)  (for Pr > 0.6)
 *
 * @param Re_L  Reynolds number based on plate length L
 * @param Pr    Prandtl number
 * @return      average Nu_L [-]
 *
 * Theorem: Nu_L(avg) = 2.Nu_L(local). This follows from integrating
 *   the local Nu_x over the plate length.
 */
double nu_flat_plate_laminar_avg(double Re_L, double Pr);

/**
 * nu_flat_plate_turbulent_local — turbulent BL, local Nu.
 *
 * Nu_x = 0.0296.Re_x^(4/5).Pr^(1/3)
 *
 * @param Re_x  local Reynolds number
 * @param Pr    Prandtl number (0.6 < Pr < 60)
 * @return      local Nu_x [-]
 *
 * Valid for: 5×10⁵ < Re_x < 10⁷.
 * Based on the 1/7th power law velocity profile and Reynolds analogy.
 */
double nu_flat_plate_turbulent_local(double Re_x, double Pr);

/**
 * nu_flat_plate_turbulent_avg — average Nu for fully turbulent plate.
 *
 * Nu_L = 0.037.Re_L^(4/5).Pr^(1/3)
 *
 * @param Re_L  Reynolds number based on entire plate length
 * @param Pr    Prandtl number
 * @return      average Nu_L [-]
 */
double nu_flat_plate_turbulent_avg(double Re_L, double Pr);

/**
 * nu_flat_plate_mixed — avg Nu for plate with laminar + turbulent regions.
 *
 * Nu_L = [0.037.Re_L^(4/5) - 871].Pr^(1/3)
 *
 * @param Re_L    Reynolds number based on plate length
 * @param Pr      Prandtl number
 * @param Re_crit critical Re for transition (typ. 5×10⁵)
 * @return        average Nu_L [-]
 *
 * This is the most commonly used flat-plate heat transfer correlation
 * in engineering design. It accounts for the laminar leading edge region
 * followed by turbulent BL development.
 *
 * Valid for: 0.6 < Pr < 60, Re_crit < Re_L < 10⁸.
 */
double nu_flat_plate_mixed(double Re_L, double Pr, double Re_crit);

/**
 * nu_flat_plate_constant_heatflux_laminar — const q", laminar.
 *
 * Nu_x = 0.453.Re_x^(1/2).Pr^(1/3)
 *
 * @param Re_x  local Re
 * @param Pr    Prandtl number
 * @return      local Nu_x [-]
 *
 * The 0.453 factor is larger than 0.332 because the temperature
 * difference DeltaT is smaller for uniform heat flux (surface temperature
 * adjusts to the flux).
 */
double nu_flat_plate_constant_heatflux_laminar(double Re_x, double Pr);

/**
 * nu_cylinder_crossflow — Hilpert/Zhukauskas for cylinder in cross-flow.
 *
 * Nu_D = C.Re_D^m.Pr^n.(Pr/Pr_w)^0.25
 *
 * @param Re_D  Reynolds number based on cylinder diameter
 * @param Pr    Prandtl number
 * @return      average Nu_D [-]
 *
 * Coefficient lookup based on Re_D ranges:
 *   Re:   0.4-4    4-40     40-4000    4000-40000    40000-400000
 *   C:    0.989    0.911    0.683      0.193         0.027
 *   m:    0.330    0.385    0.466      0.618         0.805
 *
 * Reference: Zhukauskas (1972), Hilpert (1933).
 */
double nu_cylinder_crossflow(double Re_D, double Pr);

/**
 * nu_sphere — Whitaker correlation for sphere.
 *
 * Nu_D = 2 + (0.4.Re^(1/2) + 0.06.Re^(2/3)).Pr^0.4.(mu/mu_w)^0.25
 *
 * @param Re    Reynolds number based on diameter
 * @param Pr    Prandtl number
 * @param mu    bulk dynamic viscosity
 * @param mu_w  wall-temperature dynamic viscosity
 * @return      average Nu_D [-]
 *
 * Valid for: 3.5 < Re < 7.6×10⁴, 0.7 < Pr < 380.
 * Note: Nu -> 2.0 as Re -> 0 (pure conduction from sphere to infinite medium).
 *
 * Reference: Whitaker (1972), AIChE Journal.
 */
double nu_sphere(double Re, double Pr, double mu, double mu_w);

/**
 * nu_tube_bank — Zukauskas correlation for flow across tube banks.
 *
 * Nu_D = C.Re_Dmax^m.Pr^0.36.(Pr/Pr_w)^0.25
 *
 * @param Re_max  Re based on max velocity in minimum flow area
 * @param Pr      Prandtl number
 * @param Pr_w    Pr at wall temperature (use Pr if unknown)
 * @param arrangement 0 = aligned, 1 = staggered
 * @param n_rows  number of tube rows (use n_rows for correction)
 * @return        average Nu_D [-]
 *
 * Correction factor C2 applied for fewer than 20 rows.
 * Reference: Zukauskas (1987), Handbook of Single-Phase Convective Heat Transfer.
 */
double nu_tube_bank(double Re_max, double Pr, double Pr_w,
                    int arrangement, int n_rows);

/* ==========================================================================
 * L3-L5: Nu Correlations — Natural Convection
 * ========================================================================== */

/**
 * nu_natural_convection — general correlation: Nu_L = C.Ra_L^n.
 *
 * @param Ra          Rayleigh number [-]
 * @param geometry_id 0 = vertical plate, 1 = horizontal plate (heated up),
 *                    2 = horizontal plate (heated down), 3 = horizontal cylinder,
 *                    4 = sphere, 5 = vertical cylinder (L >> D)
 * @param flow_regime 0 = laminar, 1 = turbulent
 * @return            Nu_L [-]
 *
 * Vertical plate:
 *   Laminar (10⁴ < Ra < 10⁹):    Nu = 0.59.Ra^(1/4)
 *   Turbulent (10⁹ < Ra < 10¹^3): Nu = 0.10.Ra^(1/3)
 *
 * Horizontal plate, heated upper surface (or cooled lower):
 *   Laminar (10⁴ < Ra < 10⁷):    Nu = 0.54.Ra^(1/4)
 *   Turbulent (10⁷ < Ra < 10¹¹): Nu = 0.15.Ra^(1/3)
 *
 * Horizontal cylinder:
 *   Laminar: Nu = 0.36 + 0.518.Ra^(1/4)/[1+(0.56/Pr)^(9/16)]^(4/9)
 *   (Churchill and Chu correlation is more accurate — see below)
 */
double nu_natural_convection(double Ra, int geometry_id, int flow_regime);

/**
 * nu_natural_vertical_plate_churchill_chu — Churchill & Chu (1975).
 *
 * Nu_L = {0.825 + 0.387.Ra^(1/6)/[1+(0.492/Pr)^(9/16)]^(8/27)}^2
 *
 * @param Ra  Rayleigh number based on plate height
 * @param Pr  Prandtl number
 * @return    Nu_L [-]
 *
 * Valid for the full range of Ra (laminar + turbulent). This is the
 * most accurate single correlation for vertical plates.
 *
 * Reference: Churchill & Chu (1975), Int. J. Heat Mass Transfer.
 */
double nu_natural_vertical_plate_churchill_chu(double Ra, double Pr);

/**
 * nu_natural_enclosure — heat transfer across rectangular enclosure.
 *
 * Nu = C.Ra^n.(H/L)^m
 *
 * @param Ra     Rayleigh number based on gap width L
 * @param aspect H/L (height/gap ratio)
 * @param angle  tilt angle from horizontal [degrees] (0 deg = horizontal)
 * @return       Nu [-]
 *
 * For vertical enclosure (heated from side):
 *   Nu = max(1, 0.22.(Pr.Ra/(0.2+Pr))^0.28.(H/L)^(-1/4))
 *   Valid for: 2 < H/L < 10, Pr < 10⁵, Ra < 10¹⁰.
 */
double nu_natural_enclosure(double Ra, double aspect, double angle);

/* ==========================================================================
 * L4: Energy Equation Non-Dimensionalization
 * ========================================================================== */

/**
 * energy_eq_dimensionless_form — returns the Peclet and Prandtl numbers
 *   that emerge from non-dimensionalizing the energy equation.
 *
 * The energy equation: rhocp.DT/Dt = kgrad^2T + muPhi + q""
 *
 * Using dimensionless variables: T* = (T-T0)/DeltaT, x* = x/L, u* = u/U,
 * the non-dimensional form is:
 *   Pe_H.(dT_star/dt_star + u*.grad*T*) = grad*^2T* + Br.Phi* + q*""
 * where Pe_H = Re.Pr, Br = muU^2/(kDeltaT).
 *
 * @param Re           Reynolds number
 * @param Pr           Prandtl number
 * @param mu           dynamic viscosity [Pa.s]
 * @param U            velocity [m/s]
 * @param k            thermal conductivity [W/(m.K)]
 * @param delta_T      temperature difference [K]
 * @param[out] Pe_out  Peclet number Pe_H
 * @param[out] Br_out  Brinkman number Br
 * @return             0 on success
 */
int energy_eq_dimensionless_form(double Re, double Pr,
                                 double mu, double U, double k, double delta_T,
                                 double *Pe_out, double *Br_out);

/**
 * colburn_j_factor — Colburn j-factor for heat transfer.
 *
 * j_H = St.Pr^(2/3) = (Nu/(Re.Pr)).Pr^(2/3)
 *
 * @param Nu  Nusselt number
 * @param Re  Reynolds number
 * @param Pr  Prandtl number
 * @return    Colburn j_H factor [-]
 *
 * Theorem (Colburn Analogy, 1933):
 *   j_H = f/8 where f is the Darcy friction factor.
 *   This relates momentum transfer (friction) to heat transfer (Nu).
 *   Valid for: 0.6 < Pr < 60, turbulent flow.
 */
double colburn_j_factor_heat(double Nu, double Re, double Pr);

#ifdef __cplusplus
}
#endif

#endif /* NUSSELT_NUMBER_H */

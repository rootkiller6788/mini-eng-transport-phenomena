/*
 * transport_correlations.h — Engineering Correlations based on Dimensionless Numbers
 *
 * Reference: Incropera & DeWitt (2007) "Fundamentals of Heat and Mass Transfer"
 *            Bird, Stewart, Lightfoot (2007) "Transport Phenomena"
 *            Bejan (2013) "Convection Heat Transfer"
 *            Cengel (2007) "Heat and Mass Transfer: A Practical Approach"
 *
 * This module provides engineering correlations that use dimensionless
 * numbers (Re, Nu, Pr, Sc, Gr, Ra, etc.) to solve practical problems in
 * convective heat transfer, mass transfer, and fluid dynamics.
 *
 * Knowledge Coverage: L5 Engineering Methods (correlations),
 *                     L6 Engineering Problems (design calculations),
 *                     L7 Applications (heat exchanger sizing, electronics cooling)
 */

#ifndef TRANSPORT_CORRELATIONS_H
#define TRANSPORT_CORRELATIONS_H

#ifdef __cplusplus
extern "C" {
#endif

/* ==========================================================================
 * L5: Forced Convection Correlation Suite
 * ========================================================================== */

/**
 * CorrelationResult — output of a correlation computation.
 *
 * Includes the primary Nusselt number and derived quantities for
 * engineering design.
 */
typedef struct {
    double Nu;            /* Nusselt number [-], primary result */
    double h;             /* heat transfer coefficient [W/(m^2.K)] */
    double q_flux;        /* heat flux [W/m^2] */
    double total_power;   /* total heat transfer rate [W] */
    double error_estimate; /* estimated error [%] */
} CorrelationResult;

/**
 * CorrelationInput — comprehensive input for convective heat transfer analysis.
 */
typedef struct {
    /* Fluid properties at film temperature T_f = (T_w + T_b)/2 */
    double rho;           /* density [kg/m^3] */
    double mu;            /* dynamic viscosity [Pa.s] */
    double k;             /* thermal conductivity [W/(m.K)] */
    double cp;            /* specific heat [J/(kg.K)] */
    double Pr;            /* Prandtl number [-] */
    double beta;          /* thermal expansion [1/K] */

    /* Flow conditions */
    double U;             /* free-stream or bulk velocity [m/s] */
    double T_w;           /* wall temperature [K] */
    double T_b;           /* bulk fluid temperature [K] */

    /* Geometry */
    double L;             /* characteristic length [m] (plate length, pipe diam) */
    double A_s;           /* surface area for heat transfer [m^2] */
} CorrelationInput;

/**
 * forced_convection_correlation — general forced convection engine.
 *
 * Selects and applies the appropriate Nu correlation based on geometry
 * and flow regime, then computes h, q", and Q.
 *
 * @param input    correlation input parameters
 * @param geometry 0 = flat plate, 1 = circular pipe, 2 = cylinder crossflow,
 *                 3 = sphere, 4 = tube bank aligned, 5 = tube bank staggered
 * @param bc_type  0 = uniform T_w, 1 = uniform q"_w
 * @param[out] result populated CorrelationResult
 * @return         0 on success, -1 on invalid inputs
 */
int forced_convection_correlation(const CorrelationInput *input,
                                  int geometry, int bc_type,
                                  CorrelationResult *result);

/* ==========================================================================
 * L5: Natural Convection Correlation Suite
 * ========================================================================== */

/**
 * natural_convection_vertical_plate — compute h for vertical plate.
 *
 * @param k       thermal conductivity [W/(m.K)]
 * @param L       plate height [m]
 * @param Ra      Rayleigh number based on L
 * @param Pr      Prandtl number
 * @param[out] Nu computed Nusselt number
 * @param[out] h  heat transfer coefficient [W/(m^2.K)]
 * @return        0 on success, -1 on invalid inputs
 *
 * Uses Churchill-Chu correlation valid for all Ra.
 */
int natural_convection_vertical_plate(double k, double L, double Ra, double Pr,
                                      double *Nu, double *h);

/**
 * natural_convection_horizontal_plate — heated surface facing up.
 *
 * @param k       thermal conductivity [W/(m.K)]
 * @param L       characteristic length = A_s/P [m]
 * @param Ra      Rayleigh number based on L
 * @param heated_up 1 if heated surface facing up (or cooled down),
 *                  0 if heated down
 * @param[out] Nu  computed Nusselt number
 * @param[out] h   heat transfer coefficient [W/(m^2.K)]
 * @return         0 on success
 */
int natural_convection_horizontal_plate(double k, double L, double Ra,
                                        int heated_up, double *Nu, double *h);

/**
 * natural_convection_cylinder — horizontal cylinder.
 *
 * @param k       fluid thermal conductivity [W/(m.K)]
 * @param D       cylinder diameter [m]
 * @param Ra      Rayleigh number based on D
 * @param Pr      Prandtl number
 * @param[out] Nu  Nusselt number
 * @param[out] h   heat transfer coefficient [W/(m^2.K)]
 * @return         0 on success
 */
int natural_convection_cylinder(double k, double D, double Ra, double Pr,
                                double *Nu, double *h);

/**
 * natural_convection_sphere — sphere in natural convection.
 *
 * @param k   thermal conductivity [W/(m.K)]
 * @param D   sphere diameter [m]
 * @param Ra  Rayleigh number based on D
 * @param Pr  Prandtl number
 * @param[out] Nu
 * @param[out] h
 * @return    0 on success
 */
int natural_convection_sphere(double k, double D, double Ra, double Pr,
                              double *Nu, double *h);

/**
 * natural_convection_enclosure — rectangular enclosure.
 *
 * @param k      fluid thermal conductivity [W/(m.K)]
 * @param L      gap width [m]
 * @param H      enclosure height [m]
 * @param T_h    hot wall temperature [K]
 * @param T_c    cold wall temperature [K]
 * @param Ra     Rayleigh number based on L
 * @param Pr     Prandtl number
 * @param angle  tilt from horizontal [degrees]
 * @param[out] Nu  Nusselt number
 * @param[out] h_eff  effective heat transfer coefficient [W/(m^2.K)]
 * @param[out] q_flux  heat flux [W/m^2]
 * @return        0 on success
 */
int natural_convection_enclosure(double k, double L, double H,
                                 double T_h, double T_c,
                                 double Ra, double Pr, double angle,
                                 double *Nu, double *h_eff, double *q_flux);

/* ==========================================================================
 * L5-L6: Mass Transfer Correlation Suite
 * ========================================================================== */

/**
 * mass_transfer_flat_plate — Sh for laminar flat plate.
 *
 * Sh_L = 0.664.Re_L^(1/2).Sc^(1/3)
 *
 * @param Re_L  Reynolds number based on plate length
 * @param Sc    Schmidt number
 * @param D_AB  mass diffusivity [m^2/s]
 * @param L     plate length [m]
 * @param[out] Sh   Sherwood number
 * @param[out] h_m  mass transfer coefficient [m/s]
 * @param[out] N_A  molar flux [mol/(m^2.s)] (requires C_s, C_inf separately)
 * @return      0 on success
 */
int mass_transfer_flat_plate(double Re_L, double Sc, double D_AB, double L,
                             double *Sh, double *h_m);

/**
 * mass_transfer_pipe — mass transfer in pipe flow.
 *
 * Sh = 0.023.Re^0.83.Sc^0.44 (turbulent, Gaseous)
 * Sh = 1.62.(Re.Sc.D/L)^(1/3) (laminar, developing)
 *
 * @param Re    Reynolds number
 * @param Sc    Schmidt number
 * @param D     pipe diameter [m]
 * @param L     pipe length [m]
 * @param D_AB  mass diffusivity [m^2/s]
 * @param[out] Sh
 * @param[out] h_m
 * @return      0 on success
 */
int mass_transfer_pipe(double Re, double Sc, double D, double L, double D_AB,
                       double *Sh, double *h_m);

/**
 * mass_transfer_sphere — mass transfer from a sphere.
 *
 * Frössling (1938): Sh = 2 + 0.552.Re^(1/2).Sc^(1/3)
 *
 * @param Re    Reynolds number based on diameter
 * @param Sc    Schmidt number
 * @param D     sphere diameter [m]
 * @param D_AB  mass diffusivity [m^2/s]
 * @param[out] Sh
 * @param[out] h_m
 * @return      0 on success
 *
 * Note: Sh -> 2 as Re -> 0 (pure diffusion from sphere).
 */
int mass_transfer_sphere(double Re, double Sc, double D, double D_AB,
                         double *Sh, double *h_m);

/**
 * mass_transfer_droplet_evaporation — evaporating droplet mass transfer.
 *
 * Ranz-Marshall correlation: Sh = 2 + 0.6.Re^(1/2).Sc^(1/3)
 *
 * @param Re     Reynolds number of droplet
 * @param Sc     Schmidt number (vapor in surrounding gas)
 * @param D_droplet droplet diameter [m]
 * @param D_AB   mass diffusivity [m^2/s]
 * @param rho_s  saturated vapor density at droplet surface [kg/m^3]
 * @param rho_inf ambient vapor density [kg/m^3]
 * @param[out] m_dot_evap mass evaporation rate [kg/s]
 * @return        0 on success
 *
 * Reference: Ranz & Marshall (1952), Chemical Engineering Progress.
 */
int mass_transfer_droplet_evaporation(double Re, double Sc, double D_droplet,
                                      double D_AB, double rho_s, double rho_inf,
                                      double *m_dot_evap);

/* ==========================================================================
 * L6: Engineering Design — Heat Exchanger Sizing
 * ========================================================================== */

/**
 * heat_exchanger_sizing — determine required area for given duty.
 *
 * Using the epsilon-NTU method simplified: Q = UA.DeltaT_lm
 *
 * @param Q_dot      required heat transfer rate [W]
 * @param T_h_in     hot fluid inlet temperature [K]
 * @param T_h_out    hot fluid outlet temperature [K]
 * @param T_c_in     cold fluid inlet temperature [K]
 * @param T_c_out    cold fluid outlet temperature [K]
 * @param h_hot      hot-side heat transfer coefficient [W/(m^2.K)]
 * @param h_cold     cold-side heat transfer coefficient [W/(m^2.K)]
 * @param k_wall     wall thermal conductivity [W/(m.K)]
 * @param t_wall     wall thickness [m]
 * @param[out] A_required    required area [m^2]
 * @param[out] U_overall     overall heat transfer coefficient [W/(m^2.K)]
 * @param[out] delta_T_lm    log-mean temperature difference [K]
 * @return          0 on success, -1 on invalid inputs
 */
int heat_exchanger_sizing(double Q_dot,
                          double T_h_in, double T_h_out,
                          double T_c_in, double T_c_out,
                          double h_hot, double h_cold,
                          double k_wall, double t_wall,
                          double *A_required, double *U_overall,
                          double *delta_T_lm);

/* ==========================================================================
 * L6: Pipe Flow Design
 * ========================================================================== */

/**
 * pipe_flow_design — compute pressure drop and pumping power.
 *
 * @param rho       density [kg/m^3]
 * @param mu        dynamic viscosity [Pa.s]
 * @param m_dot     mass flow rate [kg/s]
 * @param D         pipe diameter [m]
 * @param L         pipe length [m]
 * @param epsilon   wall roughness [m]
 * @param[out] Re        computed Reynolds number
 * @param[out] f         friction factor (Churchill)
 * @param[out] delta_P   pressure drop [Pa]
 * @param[out] W_dot_pump required pumping power [W]
 * @return          0 on success, -1 on invalid inputs
 */
int pipe_flow_design(double rho, double mu, double m_dot,
                     double D, double L, double epsilon,
                     double *Re, double *f, double *delta_P,
                     double *W_dot_pump);

/* ==========================================================================
 * L6: Electronics Cooling Problem
 * ========================================================================== */

/**
 * electronics_cooling_finned_array — natural convection from finned heat sink.
 *
 * Computes total heat dissipated by a vertical finned array
 * under natural convection.
 *
 * @param k_fluid     fluid thermal conductivity [W/(m.K)]
 * @param nu          kinematic viscosity [m^2/s]
 * @param alpha       thermal diffusivity [m^2/s]
 * @param beta        thermal expansion coefficient [1/K]
 * @param T_s         surface temperature [K]
 * @param T_amb       ambient temperature [K]
 * @param L_fin       fin length (vertical) [m]
 * @param H_fin       fin height (from base) [m]
 * @param t_fin       fin thickness [m]
 * @param S           fin spacing (gap) [m]
 * @param n_fins      number of fins
 * @param k_fin       fin material conductivity [W/(m.K)]
 * @param W_base      base width [m]
 * @param[out] Q_total    total heat dissipated [W]
 * @param[out] h_bar      average heat transfer coefficient [W/(m^2.K)]
 * @param[out] eta_fin    fin efficiency [-]
 * @return           0 on success
 *
 * Uses the Elenbaas correlation for isothermal vertical parallel plates
 * and fin efficiency analysis.
 */
int electronics_cooling_finned_array(double k_fluid, double nu, double alpha,
                                     double beta, double T_s, double T_amb,
                                     double L_fin, double H_fin, double t_fin,
                                     double S, int n_fins, double k_fin,
                                     double W_base,
                                     double *Q_total, double *h_bar,
                                     double *eta_fin);

/* ==========================================================================
 * L7: Application — Air-Cooled Condenser
 * ========================================================================== */

/**
 * air_cooled_condenser — condensation in horizontal tubes, air-cooled.
 *
 * @param T_steam     steam temperature [K]
 * @param T_air_in    inlet air temperature [K]
 * @param m_dot_steam steam mass flow rate [kg/s]
 * @param h_fg        latent heat of vaporization [J/kg]
 * @param D_tube      tube outer diameter [m]
 * @param L_tube      tube length (per pass) [m]
 * @param n_tubes     number of tubes
 * @param n_rows      number of tube rows
 * @param U_air       air face velocity [m/s]
 * @param S_T         transverse pitch [m]
 * @param S_L         longitudinal pitch [m]
 * @param[out] Q_condenser  condenser duty [W]
 * @param[out] T_air_out    outlet air temperature [K]
 * @param[out] A_required   required area [m^2]
 * @return           0 on success
 *
 * Models: Boeing 747 air conditioning pack condenser (~100 kW per unit).
 * NASA space suit sublimator condenser (~800 W metabolic heat).
 */
int air_cooled_condenser(double T_steam, double T_air_in,
                         double m_dot_steam, double h_fg,
                         double D_tube, double L_tube,
                         int n_tubes, int n_rows,
                         double U_air, double S_T, double S_L,
                         double *Q_condenser, double *T_air_out,
                         double *A_required);

/* ==========================================================================
 * L7: Application — DC Motor Cooling
 * ========================================================================== */

/**
 * dc_motor_natural_cooling — natural convection cooling of a DC motor housing.
 *
 * Models a cylindrical motor housing with internal heat generation
 * (I^2R losses) cooled by natural convection on the outer surface.
 *
 * @param T_amb     ambient temperature [K]
 * @param P_loss    electrical power loss [W]
 * @param D_motor   motor housing diameter [m]
 * @param L_motor   motor housing length [m]
 * @param k_air     air thermal conductivity [W/(m.K)]
 * @param nu_air    air kinematic viscosity [m^2/s]
 * @param alpha_air air thermal diffusivity [m^2/s]
 * @param[out] T_surface  computed surface temperature [K]
 * @param[out] Ra_D       Rayleigh number based on diameter
 * @param[out] Nu_D       Nusselt number
 * @return         0 on success
 *
 * Application: Toyota Prius traction motor, Tesla Model 3 drive unit,
 * industrial servo motors (e.g., Detroit automation lines).
 */
int dc_motor_natural_cooling(double T_amb, double P_loss,
                             double D_motor, double L_motor,
                             double k_air, double nu_air, double alpha_air,
                             double *T_surface, double *Ra_D, double *Nu_D);

/**
 * dc_motor_forced_cooling — forced air cooling of a DC motor.
 *
 * @param T_amb     ambient temperature [K]
 * @param P_loss    electrical power loss [W]
 * @param D_motor   motor diameter [m]
 * @param L_motor   motor length [m]
 * @param U_air     cooling air velocity [m/s]
 * @param[out] T_surface
 * @param[out] Re_D
 * @param[out] Nu_D
 * @return         0 on success
 */
int dc_motor_forced_cooling(double T_amb, double P_loss,
                            double D_motor, double L_motor,
                            double U_air,
                            double *T_surface, double *Re_D, double *Nu_D);

/* ==========================================================================
 * L7: Application — Mars Rover Thermal
 * ========================================================================== */

/**
 * mars_rover_night_survival — estimate rover internal temperature overnight.
 *
 * Mars conditions: T_amb ~ -80 degC to 20 degC, P_atm ~ 600 Pa (CO2),
 * g_mars = 3.72 m/s^2, Ra is very small due to thin atmosphere.
 *
 * @param P_heat    internal heater power [W] (typ. RHU ~ 1 W each)
 * @param A_surface rover surface area [m^2]
 * @param T_init    initial internal temperature [K]
 * @param t_night   night duration [s] (typ. 14 h = 50400 s)
 * @param C_thermal rover thermal capacitance [J/K]
 * @param[out] T_final     temperature after night [K]
 * @param[out] Q_lost      total heat lost during night [J]
 * @return         0 on success
 *
 * Application: NASA Mars 2020 Perseverance rover thermal design.
 * Uses Bi << 1 assumption (thin atmosphere -> low h -> lumped capacitance valid).
 */
int mars_rover_night_survival(double P_heat, double A_surface,
                              double T_init, double t_night,
                              double C_thermal,
                              double *T_final, double *Q_lost);

/* ==========================================================================
 * L6: Combined Free + Forced Convection
 * ========================================================================== */

/**
 * mixed_convection_assessment — assess relative importance and compute Nu.
 *
 * Uses the correlation: Nu^n_mixed = Nu^n_forced ± Nu^n_natural
 *   n = 3 for most configurations
 *   + for aiding flow, - for opposing flow
 *
 * @param Nu_forced  Nusselt number from forced convection only
 * @param Nu_natural Nusselt number from natural convection only
 * @param aiding     1 = buoyancy-aided (flow directions same),
 *                   0 = buoyancy-opposed
 * @return           mixed convection Nusselt number [-]
 */
double mixed_convection_nusselt(double Nu_forced, double Nu_natural,
                                int aiding);

/**
 * convection_regime_map — determine the convection regime.
 *
 * Returns:
 *   0 = pure forced convection (Gr/Re^2 < 0.01)
 *   1 = mixed convection aiding (0.01 < Gr/Re^2 < 10 and DeltaT.U > 0)
 *   2 = mixed convection opposing (0.01 < Gr/Re^2 < 10 and DeltaT.U < 0)
 *   3 = pure natural convection (Gr/Re^2 > 10)
 *
 * @param Gr  Grashof number
 * @param Re  Reynolds number
 * @return     regime code (0-3), -1 on invalid inputs
 */
int convection_regime_map(double Gr, double Re);

#ifdef __cplusplus
}
#endif

#endif /* TRANSPORT_CORRELATIONS_H */

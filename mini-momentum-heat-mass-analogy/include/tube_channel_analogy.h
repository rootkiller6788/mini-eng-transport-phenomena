/**
 * @file tube_channel_analogy.h
 * @brief Analogy applied to internal flows: pipes, ducts, and channels
 *
 * Internal flows are the most important engineering application of the
 * momentum-heat-mass analogy. This module provides correlations for
 * friction factor, heat transfer, and mass transfer in tubes and channels
 * based on the analogy framework.
 *
 * Reference: Bird, Stewart, Lightfoot (2007), Transport Phenomena,
 *            Chapter 6 (pipe flow), Chapters 13-14 (heat/mass in tubes).
 *            Incropera & DeWitt (2007), Chapter 8 (internal flow).
 *
 * Key results for fully-developed laminar flow in circular tubes:
 *   Nu_D = 3.66  (constant wall temperature)
 *   Nu_D = 4.36  (constant wall heat flux)
 *   Sh_D = 3.66  (constant wall concentration, analog of const T_w)
 *
 * For turbulent flow (Dittus-Boelter):
 *   Nu_D = 0.023 · Re_D^(4/5) · Pr^n  (n=0.4 heating, n=0.3 cooling)
 *   f = 0.046 · Re_D^(-1/5)           (Blasius, smooth pipe)
 *
 * The analogy yields: St·Pr^(2/3) = f/8
 */

#ifndef TUBE_CHANNEL_ANALOGY_H
#define TUBE_CHANNEL_ANALOGY_H

#include "dimensionless_groups.h"

/* ─── L2: Internal Flow Analogy Structures ──────────────────────── */

/** Pipe flow analogy state: complete description of internal flow */
typedef struct {
    double D;            /**< Pipe diameter [m] */
    double L;            /**< Pipe length [m] */
    double A_c;          /**< Cross-sectional area [m²] */
    double A_s;          /**< Surface area [m²] */
    double v_avg;        /**< Average velocity [m/s] */
    double m_dot;        /**< Mass flow rate [kg/s] */
    double Re_D;         /**< Reynolds number based on diameter */
    double Pr;           /**< Prandtl number */
    double Sc;           /**< Schmidt number */
    double f_Darcy;      /**< Darcy friction factor */
    double f_Fanning;    /**< Fanning friction factor (= f_D/4) */
    double Cf;           /**< Skin friction coefficient */
    double Nu_D;         /**< Nusselt number (heat transfer) */
    double Sh_D;         /**< Sherwood number (mass transfer) */
    double h;            /**< Heat transfer coefficient [W/(m²·K)] */
    double kc;           /**< Mass transfer coefficient [m/s] */
    double delta_P;      /**< Pressure drop [Pa] */
    double Q;            /**< Heat transfer rate [W] */
    double m_A;          /**< Mass transfer rate of species A [mol/s] */
    double j_H;          /**< Colburn j-factor (heat) */
    double j_D;          /**< Colburn j-factor (mass) */
    double St_heat;      /**< Stanton number (heat) */
    double St_mass;      /**< Stanton number (mass) */
    int flow_regime;     /**< 0=laminar, 1=transitional, 2=turbulent */
    int thermal_bc;      /**< 0=const T_w, 1=const q_w'', 2=const T with developing */
} PipeFlowAnalogy;

/** Developing flow state (thermal/hydrodynamic entrance regions) */
typedef struct {
    double x;                /**< Distance from inlet [m] */
    double D;                /**< Tube diameter [m] */
    double Re_D;             /**< Reynolds number */
    double Pr;               /**< Prandtl number */
    double x_plus_hydro;     /**< Dimensionless hydrodynamic entry length: x/(D·Re) */
    double x_plus_thermal;   /**< Dimensionless thermal entry length: x/(D·Re·Pr) */
    double L_hydro;          /**< Hydrodynamic entry length [m]: L_h ≈ 0.05·D·Re */
    double L_thermal;        /**< Thermal entry length [m]: L_th ≈ 0.05·D·Re·Pr */
    int hydro_developed;     /**< 1 if hydrodynamically fully developed */
    int thermal_developed;   /**< 1 if thermally fully developed */
} DevelopingFlowState;

/** Non-circular duct parameters for analogy correlations */
typedef struct {
    double D_h;           /**< Hydraulic diameter D_h = 4·A_c/P_w [m] */
    double A_c;           /**< Cross-sectional area [m²] */
    double P_w;           /**< Wetted perimeter [m] */
    double aspect_ratio;  /**< Aspect ratio (width/height for rectangular) */
    int geometry;         /**< 0=circular, 1=rectangular, 2=triangular,
                                3=annular, 4=parallel plates */
    double Nu_fd;         /**< Fully-developed Nusselt number (laminar) */
    double f_Re;          /**< f·Re product for laminar flow (constant) */
} DuctGeometry;

/** Friction-heat analogy validity metrics for internal flows */
typedef struct {
    double Re_D;          /**< Reynolds number */
    double Pr;            /**< Prandtl number */
    double Sc;            /**< Schmidt number */
    double f_measured;    /**< Measured friction factor */
    double Nu_measured;   /**< Measured Nusselt number */
    double Sh_measured;   /**< Measured Sherwood number */
    double Nu_analogy;    /**< Nusselt from analogy with f */
    double Sh_analogy;    /**< Sherwood from analogy with f */
    double f_analogy;     /**< Friction factor from analogy with Nu */
    double err_Nu;        /**< Relative error in Nu prediction */
    double err_Sh;        /**< Relative error in Sh prediction */
    double err_f;         /**< Relative error in f prediction */
} AnalogyValidation;

/* ─── L2: Laminar Flow Correlations ─────────────────────────────── */

/**
 * Darcy friction factor for laminar pipe flow.
 * f = 64/Re_D   (Hagen-Poiseuille)
 * Complexity: O(1)
 */
double darcy_friction_laminar(double Re_D);

/**
 * Fanning friction factor for laminar pipe flow.
 * f_F = f_D/4 = 16/Re_D
 * Complexity: O(1)
 */
double fanning_friction_laminar(double Re_D);

/**
 * Nusselt number for fully-developed laminar flow in circular tube.
 * Nu_D = 3.66  (constant wall temperature)
 * Nu_D = 4.36  (constant wall heat flux)
 * Complexity: O(1)
 */
double laminar_nu_pipe(double thermal_bc_type);

/**
 * Sherwood number for fully-developed laminar flow (analogous to Nu=3.66).
 * Sh_D = 3.66  (constant wall concentration)
 * Complexity: O(1)
 */
double laminar_sh_pipe(void);

/**
 * Nusselt number for laminar developing flow (thermal entrance).
 * Hausen correlation: Nu_D = 3.66 + 0.0668·(D/L)·Re·Pr / [1 + 0.04·((D/L)·Re·Pr)^(2/3)]
 * Valid for: thermal entrance region, constant T_w.
 * Complexity: O(1)
 */
double hausen_laminar_nu(double Re_D, double Pr, double D, double L);

/**
 * Sieder-Tate correlation for laminar flow (accounts for viscosity variation).
 * Nu_D = 1.86 · (Re·Pr·D/L)^(1/3) · (μ_b/μ_w)^(0.14)
 * Valid: Re·Pr·D/L > 10, 0.48 < Pr < 16700
 * Complexity: O(1)
 */
double sieder_tate_laminar_nu(double Re_D, double Pr, double D, double L,
                               double mu_bulk, double mu_wall);

/* ─── L2: Turbulent Flow Correlations ───────────────────────────── */

/**
 * Blasius friction factor for turbulent pipe flow (smooth).
 * f = 0.0791 · Re_D^(-1/4)   (4×10³ < Re < 10⁵)
 * f = 0.046  · Re_D^(-1/5)   (more common form)
 * Complexity: O(1)
 */
double blasius_friction_turbulent(double Re_D);

/**
 * Colebrook-White equation for turbulent friction factor.
 * 1/sqrt(f) = -2.0·log10(ε/(3.7·D) + 2.51/(Re·sqrt(f)))
 * Solved iteratively with Newton's method.
 * Complexity: O(iterations) ≈ O(10)
 */
double colebrook_friction(double Re_D, double roughness, double D);

/**
 * Dittus-Boelter correlation for turbulent pipe flow.
 * Nu_D = 0.023 · Re_D^(4/5) · Pr^n
 * n = 0.4 for heating (T_w > T_b), n = 0.3 for cooling
 * Valid: Re > 10⁴, 0.6 < Pr < 160, L/D > 10
 * Complexity: O(1)
 */
double dittus_boelter_nu(double Re_D, double Pr, int heating);

/**
 * Gnielinski correlation (more accurate than Dittus-Boelter).
 * Nu = (f/8)·(Re-1000)·Pr / [1 + 12.7·sqrt(f/8)·(Pr^(2/3)-1)]
 * Valid: 3000 < Re < 5×10⁶, 0.5 < Pr < 2000
 * Complexity: O(1)
 */
double gnielinski_nu(double Re_D, double Pr, double f);

/**
 * Petukhov correlation (high accuracy for wide Pr range).
 * Nu = (f/8)·Re·Pr / [1.07 + 12.7·sqrt(f/8)·(Pr^(2/3)-1)]
 * Complexity: O(1)
 */
double petukhov_nu(double Re_D, double Pr, double f);

/**
 * Chilton-Colburn analogy for turbulent pipe flow.
 * j_H = Nu/(Re·Pr^(1/3)) = f/8  (heat transfer j-factor)
 * j_D = Sh/(Re·Sc^(1/3)) = f/8  (mass transfer j-factor)
 *
 * @return j-factor value (= f/8 if analogy holds perfectly)
 *
 * Complexity: O(1)
 */
double colburn_j_factor_pipe(double f_F);

/* ─── L5: Analogy-Based Prediction Methods ──────────────────────── */

/**
 * Predicts Nusselt number from friction factor using the analogy.
 * Nu = (f/8) · Re · Pr^(1/3)   (Colburn form)
 *
 * This is the core engineering method: measure ΔP → get f → predict h.
 * Complexity: O(1)
 */
double predict_nu_from_friction(double f_D, double Re, double Pr);

/**
 * Predicts Sherwood number from friction factor (mass transfer analog).
 * Sh = (f/8) · Re · Sc^(1/3)
 * Complexity: O(1)
 */
double predict_sh_from_friction(double f_D, double Re, double Sc);

/**
 * Predicts heat transfer coefficient from pressure drop measurement.
 * h = k/D · (f/8) · Re · Pr^(1/3)
 *
 * This is the most practical form: measure ΔP across pipe section,
 * compute f = 2·ΔP·D/(ρ·v²·L), then predict h without any heat
 * transfer measurement.
 *
 * Complexity: O(1)
 */
double predict_h_from_pressure_drop(double delta_P, double D, double L,
                                    double rho, double v, double k,
                                    double Pr);

/**
 * Predicts mass transfer coefficient from heat transfer measurement.
 * Using the analogy: if we know h (from heat transfer experiment),
 * we can predict k_c (mass transfer coefficient) without doing a
 * mass transfer experiment.
 *
 * j_H = j_D → k_c = h/(ρ·cp) · (Sc/Pr)^(2/3)
 *
 * Complexity: O(1)
 */
double predict_kc_from_heat_transfer_coefficient(double h, double rho,
                                                  double cp, double Pr,
                                                  double Sc);

/**
 * Computes the complete pipe flow analogy state.
 * Given geometry, flow conditions, and fluid properties,
 * fills in all fields of PipeFlowAnalogy.
 *
 * Complexity: O(1) (all algebraic correlations)
 */
void compute_pipe_flow_analogy(double D, double L, double v_avg,
                               double rho, double mu, double cp,
                               double k, double D_AB,
                               double T_wall, double T_bulk,
                               double C_wall, double C_bulk,
                               PipeFlowAnalogy *pipe);

/**
 * Computes the analogy validation metrics.
 * Compares measured data with analogy predictions.
 *
 * Complexity: O(1)
 */
void validate_analogy(double Re_D, double Pr, double Sc,
                      double f_measured, double Nu_measured,
                      double Sh_measured,
                      AnalogyValidation *valid);

/* ─── L6: Non-Circular Ducts ───────────────────────────────────── */

/**
 * Sets up duct geometry parameters for non-circular cross-sections.
 * Provides f·Re product and fully-developed Nu for laminar flow.
 *
 * Laminar f·Re values:
 *   Circular: 64
 *   Parallel plates: 96
 *   Square: 56.91
 *   Equilateral triangle: 53.33
 *   Rectangular (aspect ratio): varies
 *
 * Complexity: O(1)
 */
void setup_duct_geometry(int geometry_type, double aspect_ratio,
                         double A_c, DuctGeometry *duct);

/**
 * Computes developing flow state.
 * Determines if flow is hydrodynamically/thermally developed at position x.
 *
 * Complexity: O(1)
 */
void compute_developing_flow(double x, double D, double Re_D, double Pr,
                             DevelopingFlowState *dev);

/**
 * Computes Nusselt number for developing laminar flow in entry region
 * with combined hydrodynamic and thermal development.
 *
 * Complexity: O(1)
 */
double developing_nusselt(double Re_D, double Pr, double D, double x);

/**
 * Prints complete pipe flow analogy report with all three transport
 * modes compared side-by-side.
 */
void print_pipe_flow_report(const PipeFlowAnalogy *pipe);

#endif /* TUBE_CHANNEL_ANALOGY_H */

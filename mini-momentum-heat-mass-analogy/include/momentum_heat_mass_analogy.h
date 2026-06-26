/**
 * @file momentum_heat_mass_analogy.h
 * @brief Core analogy structures connecting momentum, heat, and mass transfer
 *
 * The fundamental insight of transport phenomena: momentum, heat, and mass
 * transfer are governed by mathematically identical equations when cast in
 * dimensionless form. This allows using friction factor data to predict
 * heat and mass transfer coefficients.
 *
 * Reference: Bird, Stewart, Lightfoot (2007), Transport Phenomena,
 *            Chapters 1-2 (governing equations), Chapter 13 (analogies).
 *
 * Governing equations (steady, incompressible, constant properties):
 *
 *   Momentum:  ρ(v·∇)v = μ∇²v - ∇p + ρg
 *   Heat:      ρCp(v·∇)T = k∇²T + Φ
 *   Mass:      (v·∇)C_A = D_AB·∇²C_A + r_A
 *
 * When ν = α = D (i.e., Pr = Sc = Le = 1), the dimensionless forms are
 * identical, and the Reynolds analogy holds:
 *   St = f/2   or   Nu = (f/2)·Re·Pr
 */

#ifndef MOMENTUM_HEAT_MASS_ANALOGY_H
#define MOMENTUM_HEAT_MASS_ANALOGY_H

#include "transport_coefficients.h"

/* ─── L1: Core Analogy Definitions ─────────────────────────────── */

/**
 * Reynolds Analogy: relates skin friction to heat transfer.
 * Definition: St = (h)/(ρ·v·Cp) = f/2  (when Pr = 1)
 * Physical meaning: the same eddy mechanism transports momentum and heat.
 *
 * Stanton number St = Nu/(Re·Pr) = h/(ρ·v·Cp)
 * Friction factor f = 2·τ_w/(ρ·v²)
 */
typedef struct {
    double Cf;           /**< Skin friction coefficient Cf = 2·τ_w/(ρ·U∞²) */
    double St_heat;      /**< Stanton number for heat transfer */
    double St_mass;      /**< Stanton number for mass transfer */
    double Re;           /**< Reynolds number */
    double Pr;           /**< Prandtl number */
    double Sc;           /**< Schmidt number */
    int analogy_type;    /**< 0=Reynolds, 1=Chilton-Colburn, 2=Prandtl-Taylor */
} AnalogyState;

/**
 * Chilton-Colburn Analogy (j-factor analogy).
 * j_H = (h/(ρ·v·Cp))·Pr^(2/3) = f/2  (heat transfer j-factor)
 * j_D = (k_c/v)·Sc^(2/3) = f/2        (mass transfer j-factor)
 *
 * The Colburn j-factors are the most widely used engineering form.
 * Valid for: 0.6 < Pr < 60, 0.6 < Sc < 3000, turbulent flow.
 */
typedef struct {
    double j_H;          /**< Colburn j-factor for heat transfer */
    double j_D;          /**< Colburn j-factor for mass transfer */
    double f_F;          /**< Fanning friction factor f_F = f_Darcy/4 */
    double Cf;           /**< Skin friction coefficient Cf = 2·f_F */
    double Pr;           /**< Prandtl number */
    double Sc;           /**< Schmidt number */
    double Re;           /**< Reynolds number */
    double Nu;           /**< Nusselt number from analogy */
    double Sh;           /**< Sherwood number from analogy */
} ChiltonColburnState;

/**
 * Governing equation structure for unified transport.
 * All three transport equations can be written as:
 *   ∂φ/∂t + v·∇φ = Γ·∇²φ + S_φ
 * where φ = v (momentum), T (heat), C_A (mass),
 *       Γ = ν (momentum), α (heat), D_AB (mass).
 */
typedef struct {
    double convective_term;   /**< v·∇φ term */
    double diffusive_term;    /**< Γ·∇²φ term */
    double source_term;       /**< S_φ source/sink */
    double transient_term;    /**< ∂φ/∂t term */
    double Gamma;             /**< Transport diffusivity [m²/s] */
    double residual;          /**< Equation residual (should be ~0) */
} TransportEquationBalance;

/* ─── L2: Core Concept Functions ───────────────────────────────── */

/**
 * Computes the Stanton number for heat transfer from Reynolds analogy.
 * St = (f/2)·Pr^(-2/3)  (Chilton-Colburn form)
 * St = f/2               (Reynolds form, Pr = 1)
 *
 * Complexity: O(1)
 * @param cf    Skin friction coefficient Cf = 2·τ_w/(ρ·U∞²)
 * @param Pr    Prandtl number
 * @param colburn_form  If non-zero, use Colburn Pr^(-2/3) correction
 * @return Stanton number for heat transfer
 */
double stanton_number_from_friction(double cf, double Pr, int colburn_form);

/**
 * Computes Nusselt number from the analogy.
 * Nu = St·Re·Pr
 *
 * This is the preferred form for engineering calculations because
 * Nu directly gives the heat transfer coefficient:
 * h = Nu·k/L
 *
 * Complexity: O(1)
 */
double nusselt_from_analogy(double St, double Re, double Pr);

/**
 * Computes Sherwood number from the analogy (mass transfer analog of Nu).
 * Sh = St_mass·Re·Sc
 * Sh = (f/2)·Re·Sc^(1/3)  (Chilton-Colburn form)
 *
 * Complexity: O(1)
 */
double sherwood_from_analogy(double cf, double Re, double Sc);

/**
 * Computes the Fanning friction factor from Reynolds analogy.
 * f_F = 2·j_H = 2·St·Pr^(2/3)
 *
 * Complexity: O(1)
 */
double friction_from_heat_transfer(double St_heat, double Pr);

/**
 * Computes the complete analogy state (both Reynolds and Chilton-Colburn).
 * Given Re, Pr, Sc, and cf (skin friction), fills in all derived quantities.
 *
 * Complexity: O(1)
 */
void compute_analogy_state(double Re, double Pr, double Sc, double cf,
                           AnalogyState *analogy);

/**
 * Computes the Colburn j-factors from experimental data.
 * j_H = (h/(ρ·v·Cp))·Pr^(2/3)
 * j_D = (k_c/v)·Sc^(2/3)
 *
 * @param h      Heat transfer coefficient [W/(m²·K)]
 * @param kc     Mass transfer coefficient [m/s]
 * @param rho    Density [kg/m³]
 * @param v      Velocity [m/s]
 * @param cp     Specific heat [J/(kg·K)]
 * @param Pr     Prandtl number
 * @param Sc     Schmidt number
 * @param state  Output: filled Colburn state
 *
 * Complexity: O(1)
 */
void compute_colburn_j_factors(double h, double kc, double rho, double v,
                               double cp, double Pr, double Sc,
                               ChiltonColburnState *state);

/* ─── L4: Conservation Law Forms ───────────────────────────────── */

/**
 * Computes the dimensionless transport equation residual.
 * Verifies that the governing equation balance is zero within tolerance.
 *
 * ∂φ/∂t + v·∇φ - Γ·∇²φ - S_φ ≈ 0
 *
 * Used to validate numerical solutions against the analogy predictions.
 *
 * Complexity: O(1)
 * @return residual value (should approach 0 for a converged solution)
 */
double transport_equation_residual(double transient, double convective,
                                   double diffuse, double source);

/**
 * Checks the mathematical analogy condition: ν ≈ α ≈ D
 * (i.e., Pr ≈ 1, Sc ≈ 1, Le ≈ 1).
 *
 * Returns a "closeness" metric between 0 (completely different) and
 * 1 (exact analogy holds).
 *
 * Complexity: O(1)
 */
double analogy_closeness_metric(double Pr, double Sc);

/**
 * Verifies the Reynolds analogy by comparing both sides.
 * Left side:  St_heat (or Nu/(Re·Pr))
 * Right side: f/2  (or Cf/2)
 *
 * Returns relative error |LHS - RHS| / |RHS|.
 * For laminar flow with Pr≠1, this error will be significant.
 *
 * Complexity: O(1)
 */
double reynolds_analogy_error(double St_heat, double cf);

/* ─── L5: Engineering Methods ──────────────────────────────────── */

/**
 * Predicts heat transfer coefficient from pressure drop data
 * using the Reynolds analogy.
 *
 * Given: measured pressure drop ΔP, geometry, fluid properties
 * Returns: predicted heat transfer coefficient h [W/(m²·K)]
 *
 * This is a powerful engineering method: measure friction (easy),
 * predict heat transfer (hard to measure directly).
 *
 * Complexity: O(1)
 */
double predict_h_from_friction(double delta_P, double L, double D_h,
                               double rho, double v, double cp, double Pr);

/**
 * Predicts mass transfer coefficient from heat transfer data
 * using the Chilton-Colburn analogy (j_H = j_D).
 *
 * Given: measured heat transfer coefficient h
 * Returns: predicted mass transfer coefficient k_c [m/s]
 *
 * Complexity: O(1)
 */
double predict_kc_from_h(double h, double rho, double cp, double v,
                         double Pr, double Sc);

/**
 * Computes the Prandtl-Taylor analogy (two-layer model).
 * Accounts for laminar sublayer + turbulent core.
 *
 * St = (f/2) / [1 + 5·sqrt(f/2)·(Pr - 1 + ln(1 + 5·(Pr-1)/6))]
 *
 * More accurate than simple Reynolds analogy for wide Pr range.
 *
 * Complexity: O(1)
 */
double prandtl_taylor_analogy_St(double f_F, double Pr);

/**
 * von Karman analogy (three-layer model).
 * Accounts for laminar sublayer + buffer layer + turbulent core.
 *
 * St = (f/2) / [1 + 5·sqrt(f/2)·{(Pr-1) + ln[1 + 5·(Pr-1)/6]}]
 *
 * Complexity: O(1)
 */
double von_karman_analogy_St(double f_F, double Pr);

/**
 * Martinelli analogy for liquid metals (Pr << 1).
 * For very low Pr fluids, molecular conduction dominates over eddy conduction.
 *
 * Nu = 7.0 + 0.025·Pe^(0.8)  (Lyon-Martinelli, circular tube, uniform heat flux)
 * where Pe = Re·Pr
 *
 * Complexity: O(1)
 */
double lyon_martinelli_nu(double Re, double Pr);

/* ─── L7: Application Utilities ────────────────────────────────── */

/**
 * Computes heat transfer coefficient from Nusselt number.
 * h = Nu·k / L   [W/(m²·K)]
 *
 * Essential for heat exchanger sizing calculations.
 * Complexity: O(1)
 */
double heat_transfer_coefficient(double Nu, double k, double L);

/**
 * Computes mass transfer coefficient from Sherwood number.
 * k_c = Sh·D_AB / L   [m/s]
 *
 * Essential for absorber/stripper design.
 * Complexity: O(1)
 */
double mass_transfer_coefficient(double Sh, double D_AB, double L);

/**
 * Analogous friction factor computation: applies the analogy "in reverse".
 * Given Nu (from heat transfer), estimates friction factor.
 * f_F ≈ 2·(Nu/(Re·Pr))·Pr^(2/3)  (Colburn form)
 *
 * Complexity: O(1)
 */
double friction_from_nusselt(double Nu, double Re, double Pr);

/**
 * Checks if the analogy conditions are satisfied for a given flow.
 * Returns a summary string describing the analogy regime.
 *
 * Analogy conditions:
 *   - Pr ≈ 1 for Reynolds analogy to hold exactly
 *   - 0.5 < Pr < 50 for Chilton-Colburn validity
 *   - Turbulent flow (Re > Re_crit) for eddy-based analogies
 *
 * @param buf  Output buffer (at least 256 bytes)
 * @param size Buffer size
 * @return pointer to buf
 */
char *analogy_regime_description(double Re, double Pr, double Sc,
                                 char *buf, size_t size);

/**
 * Prints the full analogy diagnostic table showing all three transport
 * mechanisms side by side with their dimensionless numbers.
 */
void print_analogy_diagnostics(const ChiltonColburnState *ccs);

#endif /* MOMENTUM_HEAT_MASS_ANALOGY_H */

/**
 * @file cdr_diffusion.h
 * @brief Diffusion Transport Models
 *
 * Fick's law of diffusion, effective diffusivity in porous media,
 * multi-component diffusion (Stefan-Maxwell), and transient diffusion.
 *
 * Reference: Bird, Stewart, Lightfoot "Transport Phenomena" (2007), Ch. 17-18
 * Reference: Cussler "Diffusion: Mass Transfer in Fluid Systems" (2009)
 *
 * L1: Fick's first/second law, diffusivity definitions
 * L2: Diffusion regimes (molecular, Knudsen, surface)
 * L4: Species mass conservation with diffusion
 * L5: Diffusion coefficient estimation methods
 */

#ifndef CDR_DIFFUSION_H
#define CDR_DIFFUSION_H

#include "cdr_core.h"

/* ---------------------------------------------------------------------------
 * L1: Fick's Law Structures
 * ------------------------------------------------------------------------- */

/**
 * @brief 1D concentration profile for transient diffusion analysis.
 * Stores the spatial discretization of concentration field C(x,t).
 */
typedef struct {
    double *C;          /**< Concentration array [mol/m^3], length = n_points */
    double *x;          /**< Spatial positions [m], length = n_points */
    double  dx;         /**< Spatial step size [m] */
    double  L;          /**< Domain length [m] */
    double  t;          /**< Current time [s] */
    size_t  n_points;   /**< Number of spatial points */
} ConcentrationProfile1D;

/**
 * @brief Diffusion flux result from Fick's law.
 */
typedef struct {
    double J_diff;           /**< Diffusive flux [mol/(m^2.s)] */
    double J_diff_x;         /**< x-component [mol/(m^2.s)] */
    double J_diff_y;         /**< y-component [mol/(m^2.s)] */
    double J_diff_z;         /**< z-component [mol/(m^2.s)] */
    double magnitude;        /**< Flux magnitude [mol/(m^2.s)] */
} DiffusionFlux;

/**
 * @brief Diffusion coefficient regimes for porous media.
 */
typedef enum {
    DIFF_MOLECULAR,        /**< Bulk molecular diffusion (pore >> mean free path) */
    DIFF_KNUDSEN,          /**< Knudsen diffusion (pore ~ mean free path) */
    DIFF_SURFACE,          /**< Surface diffusion along pore walls */
    DIFF_CONFIGURATIONAL   /**< Configurational diffusion (zeolites, pore ~ molecule) */
} DiffusionRegime;

/**
 * @brief Multi-component diffusion matrix (Stefan-Maxwell formulation).
 *
 * For N species, D_ij forms an (N x N) matrix where:
 *   J_i = -sum_j D_ij * grad(C_j)
 *
 * The diagonal terms D_ii represent self-diffusion.
 */
typedef struct {
    double **D;          /**< Diffusion coefficient matrix [m^2/s], size n x n */
    size_t   n_species;  /**< Number of species */
    double  *x_mole;     /**< Mole fractions, length n_species */
    double  *C_total;    /**< Total concentration [mol/m^3] */
} MulticomponentDiffusion;

/* ---------------------------------------------------------------------------
 * L2: Effective Diffusivity Models
 * ------------------------------------------------------------------------- */

/**
 * @brief Compute effective diffusivity in porous media.
 *
 * D_eff = (epsilon / tau) * D_bulk
 *
 * where epsilon = porosity, tau = tortuosity.
 *
 * This is the parallel-pore model (simplest).
 * @param D_bulk     Bulk diffusion coefficient [m^2/s]
 * @param porosity   Void fraction epsilon [0-1]
 * @param tortuosity Tortuosity factor tau [>= 1]
 * @return Effective diffusivity D_eff [m^2/s]
 */
double cdr_diffusion_effective(double D_bulk, double porosity, double tortuosity);

/**
 * @brief Compute Knudsen diffusion coefficient for gas in narrow pores.
 *
 * D_K = (d_pore / 3) * sqrt(8*R*T / (pi * M))
 *
 * where d_pore = pore diameter, M = molecular weight.
 * Applicable when pore diameter < mean free path.
 *
 * @param pore_diameter    Pore diameter [m]
 * @param temperature      Temperature [K]
 * @param molecular_weight Molecular weight [kg/mol]
 * @return Knudsen diffusivity D_K [m^2/s]
 */
double cdr_diffusion_knudsen(double pore_diameter, double temperature,
                             double molecular_weight);

/**
 * @brief Combined diffusivity in transition regime (Bosanquet formula).
 *
 * 1/D_eff = 1/D_bulk + 1/D_K
 *
 * @param D_bulk  Bulk (molecular) diffusivity [m^2/s]
 * @param D_K     Knudsen diffusivity [m^2/s]
 * @return Combined effective diffusivity [m^2/s]
 */
double cdr_diffusion_bosanquet(double D_bulk, double D_K);

/**
 * @brief Estimate binary gas diffusion coefficient (Fuller-Schettler-Giddings).
 *
 * D_AB = 1.0e-7 * T^1.75 * sqrt(1/M_A + 1/M_B) /
 *        (P * (V_A^(1/3) + V_B^(1/3))^2)   [m^2/s]
 *
 * @param T            Temperature [K]
 * @param P            Pressure [atm]
 * @param M_A, M_B     Molecular weights [g/mol]
 * @param V_A, V_B     Diffusion volumes [cm^3/mol] (tabulated)
 * @return Binary diffusion coefficient D_AB [m^2/s]
 */
double cdr_diffusion_fuller(double T, double P,
                            double M_A, double M_B,
                            double V_A, double V_B);

/**
 * @brief Estimate liquid diffusion coefficient (Wilke-Chang).
 *
 * D_AB = 7.4e-12 * (phi_B * M_B)^0.5 * T / (mu_B * V_A^0.6)   [m^2/s]
 *
 * @param T            Temperature [K]
 * @param mu_B         Solvent viscosity [cP]
 * @param M_B          Solvent molecular weight [g/mol]
 * @param V_A          Solute molar volume at boiling point [cm^3/mol]
 * @param phi_B        Solvent association factor (2.6 water, 1.9 methanol, 1.0 benzene)
 * @return Liquid diffusion coefficient D_AB [m^2/s]
 */
double cdr_diffusion_wilke_chang(double T, double mu_B, double M_B,
                                  double V_A, double phi_B);

/* ---------------------------------------------------------------------------
 * L4: Fick's Law Implementations
 * ------------------------------------------------------------------------- */

/**
 * @brief Fick's first law: compute diffusive flux from concentration gradient.
 *
 * J = -D * dC/dx  (1D)
 *
 * @param D   Diffusion coefficient [m^2/s]
 * @param C1  Concentration at point 1 [mol/m^3]
 * @param C2  Concentration at point 2 [mol/m^3]
 * @param dx  Distance between points [m]
 * @return Diffusive flux J [mol/(m^2.s)], positive if C2 > C1
 */
double cdr_fick_first_law_1d(double D, double C1, double C2, double dx);

/**
 * @brief Fick's first law in 3D: compute diffusive flux vector.
 *
 * J = -D * grad(C)
 *
 * @param D     Diffusion coefficient [m^2/s]
 * @param dCdx  Concentration gradient in x [mol/m^4]
 * @param dCdy  Concentration gradient in y [mol/m^4]
 * @param dCdz  Concentration gradient in z [mol/m^4]
 * @param flux  Output: populated DiffusionFlux struct
 */
void cdr_fick_first_law_3d(double D, double dCdx, double dCdy, double dCdz,
                            DiffusionFlux *flux);

/**
 * @brief Fick's second law: 1D transient diffusion, explicit Euler step.
 *
 * dC/dt = D * d^2C/dx^2
 *
 * Forward-time centered-space (FTCS) scheme.
 * Stability: D*dt/dx^2 <= 0.5
 *
 * @param profile  Current concentration profile (updated in-place)
 * @param D        Diffusion coefficient [m^2/s]
 * @param dt       Time step [s]
 * @return Maximum stable dt, or 0 if profile is NULL
 */
double cdr_fick_second_law_1d_step(ConcentrationProfile1D *profile,
                                    double D, double dt);

/**
 * @brief Analytical solution for 1D semi-infinite diffusion.
 *
 * C(x,t) = C_s + (C_0 - C_s) * erf(x / (2*sqrt(D*t)))
 *
 * @param C0  Initial uniform concentration [mol/m^3]
 * @param Cs  Surface concentration (boundary) [mol/m^3]
 * @param D   Diffusion coefficient [m^2/s]
 * @param x   Position from surface [m]
 * @param t   Time [s]
 * @return Concentration C(x,t) [mol/m^3]
 */
double cdr_diffusion_semi_infinite(double C0, double Cs, double D,
                                    double x, double t);

/* ---------------------------------------------------------------------------
 * L5: Penetration depth calculation
 * ------------------------------------------------------------------------- */

/**
 * @brief Diffusion penetration depth (99% criterion).
 *
 * delta = 4 * sqrt(D * t)
 *
 * This is the distance at which C(x,t) reaches ~99% of (C0 - Cs).
 *
 * @param D  Diffusion coefficient [m^2/s]
 * @param t  Time [s]
 * @return Penetration depth [m]
 */
double cdr_diffusion_penetration_depth(double D, double t);

/**
 * @brief Characteristic diffusion time over length L.
 *
 * tau_diff = L^2 / D
 *
 * @param L  Characteristic length [m]
 * @param D  Diffusion coefficient [m^2/s]
 * @return Diffusion timescale [s]
 */
double cdr_diffusion_timescale(double L, double D);

/**
 * @brief Steady-state diffusion through a flat membrane.
 *
 * J_ss = D * (C_high - C_low) / L
 *
 * @param D       Diffusion coefficient [m^2/s]
 * @param C_high  Concentration on high side [mol/m^3]
 * @param C_low   Concentration on low side [mol/m^3]
 * @param L       Membrane thickness [m]
 * @return Steady-state flux [mol/(m^2.s)]
 */
double cdr_diffusion_membrane_flux(double D, double C_high, double C_low, double L);

/**
 * @brief Steady-state radial diffusion in a cylinder (1D).
 *
 * C(r) = C_inner + (C_outer - C_inner) * ln(r/R_inner) / ln(R_outer/R_inner)
 *
 * @param r           Radial position [m]
 * @param R_inner     Inner radius [m]
 * @param R_outer     Outer radius [m]
 * @param C_inner     Concentration at inner wall [mol/m^3]
 * @param C_outer     Concentration at outer wall [mol/m^3]
 * @return Concentration C(r) [mol/m^3]
 */
double cdr_diffusion_cylindrical_ss(double r, double R_inner, double R_outer,
                                     double C_inner, double C_outer);

/**
 * @brief Steady-state radial diffusion in a sphere (1D).
 *
 * C(r) = C_inner + (C_outer - C_inner) *
 *        (1/r - 1/R_inner) / (1/R_outer - 1/R_inner) * (-1)
 *
 * @param r           Radial position [m]
 * @param R_inner     Inner radius [m]
 * @param R_outer     Outer radius [m]
 * @param C_inner     Concentration at inner wall [mol/m^3]
 * @param C_outer     Concentration at outer wall [mol/m^3]
 * @return Concentration C(r) [mol/m^3]
 */
double cdr_diffusion_spherical_ss(double r, double R_inner, double R_outer,
                                   double C_inner, double C_outer);

#endif /* CDR_DIFFUSION_H */

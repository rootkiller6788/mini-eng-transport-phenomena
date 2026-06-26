/**
 * @file cdr_mass_transfer.h
 * @brief Interphase Mass Transfer and Convection-Enhanced Transport
 *
 * Mass transfer coefficients, Sherwood number correlations,
 * boundary layer theory for mass transfer, two-film theory,
 * and mass transfer with chemical reaction.
 *
 * Reference: Bird, Stewart, Lightfoot "Transport Phenomena" (2007), Ch. 21-22
 * Reference: Incropera & DeWitt "Fundamentals of Heat and Mass Transfer" (2007)
 * Reference: Treybal "Mass Transfer Operations" (1980)
 *
 * L1: Mass transfer coefficient, Sherwood number, mass transfer driving force
 * L2: Boundary layer concepts for mass transfer, Chilton-Colburn analogy
 * L3: Correlations for k_c and k_G in various geometries
 * L5: Two-film theory, HTU/NTU method for mass transfer
 */

#ifndef CDR_MASS_TRANSFER_H
#define CDR_MASS_TRANSFER_H

#include "cdr_core.h"
#include "cdr_diffusion.h"

/* ---------------------------------------------------------------------------
 * L1: Mass Transfer Coefficient Definitions
 * ------------------------------------------------------------------------- */

/**
 * @brief Mass transfer coefficient with different driving forces.
 *
 * N_A = k_c * (C_A1 - C_A2)       [liquid phase, concentration driving force]
 * N_A = k_G * (p_A1 - p_A2)       [gas phase, partial pressure driving force]
 * N_A = k_y * (y_A1 - y_A2)       [gas phase, mole fraction driving force]
 * N_A = k_x * (x_A1 - x_A2)       [liquid phase, mole fraction driving force]
 *
 * k_c: [m/s], k_G: [mol/(m^2.s.Pa)], k_y: [mol/(m^2.s)], k_x: [mol/(m^2.s)]
 */
typedef struct {
    double k_c;       /**< Liquid-side mass transfer coefficient [m/s] */
    double k_G;       /**< Gas-side mass transfer coefficient [mol/(m^2.s.Pa)] */
    double k_y;       /**< Gas-side mole-fraction coefficient [mol/(m^2.s)] */
    double k_x;       /**< Liquid-side mole-fraction coefficient [mol/(m^2.s)] */
    double k_overall; /**< Overall mass transfer coefficient [m/s] or [mol/(m^2.s.Pa)] */
} MassTransferCoefficients;

/**
 * @brief Driving force specification for mass transfer.
 */
typedef enum {
    DRIVING_CONCENTRATION,   /**< Delta C [mol/m^3] */
    DRIVING_PARTIAL_PRESSURE,/**< Delta p [Pa] */
    DRIVING_MOLE_FRACTION,   /**< Delta y or Delta x [-] */
    DRIVING_CHEMICAL_POTENTIAL /**< Delta mu [J/mol] (rigorous) */
} DrivingForceType;

/**
 * @brief Mass transfer flux computation result.
 */
typedef struct {
    double N_A;           /**< Molar flux of A [mol/(m^2.s)] */
    double N_A_diffusion; /**< Diffusive contribution [mol/(m^2.s)] */
    double N_A_convection;/**< Convective contribution [mol/(m^2.s)] */
    double enhancement;   /**< Enhancement factor due to reaction (> 1) */
} MassTransferFlux;

/* ---------------------------------------------------------------------------
 * L2: Sherwood Number and Correlations
 * ------------------------------------------------------------------------- */

/**
 * @brief Sherwood number ? ratio of convective to diffusive mass transfer.
 *
 * Sh = k_c * L / D_AB
 *
 * Analogous to Nusselt number in heat transfer.
 */
typedef struct {
    double Sh;          /**< Sherwood number [-] */
    double Sh_laminar;  /**< Sh for laminar flow */
    double Sh_turbulent;/**< Sh for turbulent flow */
} SherwoodNumbers;

/**
 * @brief Chilton-Colburn j-factor analogy.
 *
 * j_D = (k_c / u) * Sc^(2/3) = f/2  (mass transfer)
 * j_H = (h / (rho*Cp*u)) * Pr^(2/3) = f/2  (heat transfer)
 *
 * j_D = j_H = f/2  (Reynolds analogy for Sc = Pr = 1)
 */
typedef struct {
    double j_D;         /**< Mass transfer j-factor [-] */
    double j_H;         /**< Heat transfer j-factor [-] */
    double friction_f;  /**< Fanning friction factor [-] */
} ColburnJFactors;

/* ---------------------------------------------------------------------------
 * L3: Mass Transfer Coefficient Correlations
 * ------------------------------------------------------------------------- */

/**
 * @brief Sherwood number for laminar flow over a flat plate.
 *
 * Sh_x = 0.332 * Re_x^(1/2) * Sc^(1/3)   [local]
 * Sh_L = 0.664 * Re_L^(1/2) * Sc^(1/3)   [average, Re_L < 5e5]
 *
 * @param Re   Reynolds number (based on plate length) [-]
 * @param Sc   Schmidt number [-]
 * @return Average Sherwood number [-]
 */
double cdr_sherwood_flat_plate_laminar(double Re, double Sc);

/**
 * @brief Sherwood number for turbulent flow over a flat plate.
 *
 * Sh_L = 0.037 * Re_L^(4/5) * Sc^(1/3)   [average, Re_L > 5e5]
 *
 * @param Re   Reynolds number [-]
 * @param Sc   Schmidt number [-]
 * @return Average Sherwood number [-]
 */
double cdr_sherwood_flat_plate_turbulent(double Re, double Sc);

/**
 * @brief Sherwood number for flow past a single sphere (Froessling).
 *
 * Sh = 2 + 0.552 * Re^(1/2) * Sc^(1/3)
 *
 * The "2" corresponds to pure diffusion from a sphere into stagnant medium.
 *
 * @param Re   Particle Reynolds number Re_p [-]
 * @param Sc   Schmidt number [-]
 * @return Sherwood number [-]
 */
double cdr_sherwood_sphere(double Re, double Sc);

/**
 * @brief Sherwood number for flow through packed bed (Thoenes-Kramers).
 *
 * Sh' = 1.0 * Re'^(1/2) * Sc^(1/3)
 *
 * where Sh' = Sh * epsilon / (gamma*(1-epsilon)), Re' = Re / (gamma*(1-epsilon))
 * gamma = particle shape factor (1 for sphere)
 *
 * @param Re       Particle Reynolds number [-]
 * @param Sc       Schmidt number [-]
 * @param epsilon  Bed porosity [-]
 * @return Sherwood number (conventional definition) [-]
 */
double cdr_sherwood_packed_bed(double Re, double Sc, double epsilon);

/**
 * @brief Mass transfer coefficient from Sherwood number.
 *
 * k_c = Sh * D_AB / L
 *
 * @param Sh  Sherwood number [-]
 * @param D   Diffusion coefficient [m^2/s]
 * @param L   Characteristic length [m]
 * @return Mass transfer coefficient k_c [m/s]
 */
double cdr_mass_transfer_coeff_from_sh(double Sh, double D, double L);

/**
 * @brief Mass transfer coefficient for gas-liquid interface (penetration theory).
 *
 * k_L = 2 * sqrt(D / (pi * t_e))
 *
 * where t_e = exposure time.
 *
 * @param D    Diffusion coefficient [m^2/s]
 * @param t_e  Exposure (contact) time [s]
 * @return Liquid-side mass transfer coefficient [m/s]
 */
double cdr_mass_transfer_penetration(double D, double t_e);

/**
 * @brief Mass transfer coefficient for falling film (laminar).
 *
 * k_L_avg = (6 * D * u_s / (pi * L))^(1/2)  [short contact]
 * k_L_avg = 3.41 * D / delta_f                [long contact]
 *
 * @param D       Diffusion coefficient [m^2/s]
 * @param u_s     Surface velocity [m/s]
 * @param L       Film length [m]
 * @return Average mass transfer coefficient [m/s]
 */
double cdr_mass_transfer_falling_film(double D, double u_s, double L);

/* ---------------------------------------------------------------------------
 * L5: Two-Film Theory
 * ------------------------------------------------------------------------- */

/**
 * @brief Overall mass transfer coefficient from two-film theory.
 *
 * 1/K_OG = 1/k_y + m'/k_x  (gas-side overall)
 * 1/K_OL = 1/(m''*k_y) + 1/k_x  (liquid-side overall)
 *
 * where m', m'' are slopes of equilibrium line.
 *
 * @param k_G    Gas-side coefficient [mol/(m^2.s.Pa)]
 * @param k_L    Liquid-side coefficient (k_c) [m/s]
 * @param H      Henry's constant [Pa.m^3/mol]
 * @param P      Total pressure [Pa]
 * @param K_overall Output: overall coefficient [mol/(m^2.s.Pa)] (gas-side basis)
 */
void cdr_two_film_overall_gas(double k_G, double k_L, double H, double P,
                               double *K_overall);

/**
 * @brief Mass flux across gas-liquid interface (two-film theory).
 *
 * N_A = K_OG * (p_A - H*C_A)  [gas-side driving force]
 *
 * @param K_overall  Overall mass transfer coefficient [mol/(m^2.s.Pa)]
 * @param p_A_bulk   Bulk gas partial pressure [Pa]
 * @param C_A_bulk   Bulk liquid concentration [mol/m^3]
 * @param H          Henry's constant [Pa.m^3/mol]
 * @return Interfacial molar flux N_A [mol/(m^2.s)]
 */
double cdr_two_film_flux(double K_overall, double p_A_bulk,
                          double C_A_bulk, double H);

/**
 * @brief Mass transfer with first-order chemical reaction (enhancement factor).
 *
 * E = Ha / tanh(Ha)  for irreversible first-order reaction
 *
 * where Ha = sqrt(k*D) / k_L  (Hatta number)
 *
 * @param k    First-order rate constant [1/s]
 * @param D    Diffusion coefficient [m^2/s]
 * @param k_L  Mass transfer coefficient without reaction [m/s]
 * @return Enhancement factor E (>= 1)
 */
double cdr_enhancement_factor_homogeneous(double k, double D, double k_L);

/**
 * @brief Hatta number for gas-liquid reaction.
 *
 * Ha = sqrt(2/(m+1) * k_{m,n} * D_A * C_Ai^(m-1) * C_B_bulk^n) / k_L
 *
 * For pseudo-first-order (C_B_bulk >> C_Ai): Ha = sqrt(k1 * D_A) / k_L
 *
 * @param k1     Pseudo-first-order rate constant [1/s]
 * @param D_A    Diffusivity of A [m^2/s]
 * @param k_L    Liquid-side mass transfer coefficient [m/s]
 * @return Hatta number Ha [-]
 */
double cdr_hatta_number(double k1, double D_A, double k_L);

/**
 * @brief Regime classification for gas-liquid reaction based on Hatta number.
 *
 * Ha < 0.3   : Slow reaction (kinetics control, bulk liquid)
 * 0.3 < Ha < 3 : Intermediate
 * Ha > 3      : Fast reaction (in film, enhancement)
 * Ha > 10*E_i : Instantaneous reaction (at interface)
 *
 * @param Ha  Hatta number [-]
 * @return String describing the regime
 */
const char* cdr_classify_hatta_regime(double Ha);

/**
 * @brief Number of Transfer Units (NTU) for mass transfer.
 *
 * NTU_G = integral dy / (y_i - y)
 * HTU_G = G / (k_y * a)   [m]
 *
 * Height of column = HTU * NTU
 *
 * @param y_in    Inlet mole fraction [-]
 * @param y_out   Outlet mole fraction [-]
 * @param y_eq    Equilibrium mole fraction at outlet [-]
 * @return NTU_G for linear operating/equilibrium lines
 */
double cdr_ntu_gas_phase(double y_in, double y_out, double y_eq);

/**
 * @brief Height of Transfer Unit (HTU) for gas phase.
 *
 * @param G    Gas molar velocity [mol/(m^2.s)]
 * @param k_y  Gas-side mass transfer coefficient [mol/(m^2.s)]
 * @param a    Interfacial area per volume [m^2/m^3]
 * @return HTU_G [m]
 */
double cdr_htu_gas_phase(double G, double k_y, double a);

#endif /* CDR_MASS_TRANSFER_H */

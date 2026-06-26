/**
 * interphase_transport.c
 * =======================
 * Transport coefficient calculations and flux models across interfaces.
 *
 * Reference: Bird, Stewart, Lightfoot (2007) "Transport Phenomena"
 *   Chapter 14 - Interphase Transport in Nonisothermal Systems
 *   Chapter 22 - Interphase Transport in Isothermal Systems
 *
 * Cussler, E.L. (2009) "Diffusion: Mass Transfer in Fluid Systems" (3rd ed.)
 * Treybal, R.E. (1980) "Mass-Transfer Operations" (3rd ed.)
 *
 * Knowledge coverage: L2 Core Concepts, L3 Engineering Quantities, L5 Methods
 */

#include "interphase_transport.h"
#include <math.h>
#include <string.h>

/* ============================================================================
 * L5: Mass Transfer Coefficient Models
 * ============================================================================ */

void two_film_mass_transfer(double k_a, double k_b, double m_partition,
                             TwoFilmModel *film) {
    if (!film || k_a <= 0.0 || k_b <= 0.0 || m_partition <= 0.0) {
        if (film) memset(film, 0, sizeof(TwoFilmModel));
        return;
    }

    film->k_a = k_a;
    film->k_b = k_b;

    /* Overall coefficient on A-phase basis (gas-side for gas-liquid) */
    /* 1/K_A = 1/k_A + m/k_B */
    double inv_K_a = 1.0 / k_a + m_partition / k_b;
    film->K_overall_a = 1.0 / inv_K_a;

    /* Overall coefficient on B-phase basis (liquid-side) */
    /* 1/K_B = 1/(m*k_A) + 1/k_B */
    double inv_K_b = 1.0 / (m_partition * k_a) + 1.0 / k_b;
    film->K_overall_b = 1.0 / inv_K_b;

    /* Film thicknesses from k = D/delta (approximate) */
    /* Without D, we just store k values */
    film->film_thickness_a = 0.0;
    film->film_thickness_b = 0.0;
    film->h_a = 0.0;
    film->h_b = 0.0;
    film->U_overall = 0.0;
}

void penetration_mass_transfer(double diffusivity, double exposure_time,
                                PenetrationModel *model) {
    if (!model || diffusivity <= 0.0 || exposure_time <= 0.0) {
        if (model) memset(model, 0, sizeof(PenetrationModel));
        return;
    }

    model->diffusivity_a   = diffusivity;
    model->diffusivity_b   = diffusivity;
    model->exposure_time   = exposure_time;

    /* k_avg = 2 * sqrt(D / (pi * t_e)) */
    double sqrt_factor = sqrt(diffusivity / (M_PI * exposure_time));
    model->k_penetration_a = 2.0 * sqrt_factor;
    model->k_penetration_b = 2.0 * sqrt_factor;
}

void surface_renewal_mass_transfer(double diffusivity, double renewal_rate,
                                    SurfaceRenewalModel *model) {
    if (!model || diffusivity <= 0.0 || renewal_rate <= 0.0) {
        if (model) memset(model, 0, sizeof(SurfaceRenewalModel));
        return;
    }

    model->diffusivity      = diffusivity;
    model->renewal_rate     = renewal_rate;
    model->k_renewal        = sqrt(diffusivity * renewal_rate);
    model->mean_element_age = 1.0 / renewal_rate;
}

void compare_mass_transfer_theories(double film_k, double penetration_k,
                                     double renewal_k, double ratios[3]) {
    if (!ratios) return;
    ratios[0] = (film_k > 0.0) ? (penetration_k / film_k) : 0.0;
    ratios[1] = (film_k > 0.0) ? (renewal_k / film_k) : 0.0;
    ratios[2] = (penetration_k > 0.0) ? (renewal_k / penetration_k) : 0.0;
}

/* ============================================================================
 * L5: Heat Transfer Correlations for Interface Problems
 * ============================================================================ */

void nusselt_flat_plate_laminar(double Re_x, double Pr,
                                 double *Nu_x, double *Nu_avg) {
    if (Re_x <= 0.0 || Pr <= 0.0) {
        if (Nu_x)   *Nu_x   = 0.0;
        if (Nu_avg) *Nu_avg = 0.0;
        return;
    }

    /* Pohlhausen solution (1921): Nu_x = 0.332 * Re_x^(1/2) * Pr^(1/3) */
    if (Nu_x) {
        *Nu_x = 0.332 * sqrt(Re_x) * pow(Pr, 1.0/3.0);
    }

    /* Average over length L: Nu_L = 0.664 * Re_L^(1/2) * Pr^(1/3) */
    if (Nu_avg) {
        *Nu_avg = 0.664 * sqrt(Re_x) * pow(Pr, 1.0/3.0);
    }
}

void nusselt_flat_plate_turbulent(double Re_x, double Pr,
                                   double *Nu_x, double *Nu_avg) {
    if (Re_x <= 0.0 || Pr <= 0.0) {
        if (Nu_x)   *Nu_x   = 0.0;
        if (Nu_avg) *Nu_avg = 0.0;
        return;
    }

    /* Colburn analogy: Nu_x = 0.0296 * Re_x^(4/5) * Pr^(1/3) */
    if (Nu_x) {
        *Nu_x = 0.0296 * pow(Re_x, 0.8) * pow(Pr, 1.0/3.0);
    }

    /* Average: Nu_L = 0.037 * Re_L^(4/5) * Pr^(1/3) */
    if (Nu_avg) {
        *Nu_avg = 0.037 * pow(Re_x, 0.8) * pow(Pr, 1.0/3.0);
    }
}

void nusselt_pipe_laminar(double Re_D, double Pr, double L_over_D,
                           double mu_ratio, int is_const_T, double *Nu) {
    if (!Nu || Re_D <= 0.0 || Pr <= 0.0 || L_over_D <= 0.0) {
        if (Nu) *Nu = 0.0;
        return;
    }

    double Gz = Re_D * Pr / L_over_D;  /* Graetz number */

    if (Gz > 1000.0) {
        /* Thermal entry region: Hausen correlation */
        double Nu_entry = 3.66;
        if (is_const_T) {
            Nu_entry = 3.66 + (0.0668 * Gz)
                       / (1.0 + 0.04 * pow(Gz, 2.0/3.0));
        } else {
            Nu_entry = 4.36 + (0.023 * Gz)
                       / (1.0 + 0.0012 * Gz);
        }
        *Nu = Nu_entry * pow(mu_ratio, 0.14);
    } else {
        /* Fully developed */
        *Nu = is_const_T ? 3.66 : 4.36;
    }
}

void nusselt_pipe_turbulent(double Re_D, double Pr, int heating,
                             double *Nu_dittus, double *Nu_gnielinski) {
    if (Re_D < 2300.0 || Pr <= 0.5) {
        if (Nu_dittus)    *Nu_dittus    = 0.0;
        if (Nu_gnielinski) *Nu_gnielinski = 0.0;
        return;
    }

    /* Dittus-Boelter (1930) */
    double n = heating ? 0.4 : 0.3;
    if (Nu_dittus) {
        *Nu_dittus = 0.023 * pow(Re_D, 0.8) * pow(Pr, n);
    }

    /* Gnielinski (1976) - valid for 0.5 <= Pr <= 2000, 2300 < Re < 5e6 */
    if (Nu_gnielinski) {
        /* Blasius friction factor: f = 0.079 * Re^(-0.25) */
        double f = 0.079 * pow(Re_D, -0.25);
        double denom = 1.0 + 12.7 * sqrt(f/8.0) * (pow(Pr, 2.0/3.0) - 1.0);
        if (denom > 0.0) {
            *Nu_gnielinski = (f/8.0) * (Re_D - 1000.0) * Pr / denom;
        } else {
            *Nu_gnielinski = 0.0;
        }
    }
}

void nusselt_natural_convection_vertical(double Ra, double Pr, double *Nu) {
    if (!Nu || Ra <= 0.0 || Pr <= 0.0) {
        if (Nu) *Nu = 0.0;
        return;
    }

    /* Churchill-Chu correlation - universal for vertical plate */
    double denom = pow(1.0 + pow(0.492 / Pr, 9.0/16.0), 8.0/27.0);
    double term  = 0.387 * pow(Ra, 1.0/6.0) / denom;
    double sum   = 0.825 + term;
    *Nu = sum * sum;
}

void nusselt_film_condensation(double g, double rho_L, double rho_v,
                                double h_fg, double L, double mu_L,
                                double k_L, double delta_T, double *Nu) {
    if (!Nu || delta_T <= 0.0 || mu_L <= 0.0 || k_L <= 0.0) {
        if (Nu) *Nu = 0.0;
        return;
    }

    /* Nusselt (1916) laminar film condensation on vertical plate */
    double drho = rho_L - rho_v;
    if (drho <= 0.0) { *Nu = 0.0; return; }

    (void)0.0;  /* Ja = cp_L*delta_T/h_fg, reserved for correction */
    /* Nu_L = 0.943 * (g*drho*h_fg*L^3/(mu_L*k_L*delta_T))^(1/4) */
    double arg = g * drho * h_fg * L*L*L
                 / (mu_L * k_L * delta_T);
    if (arg <= 0.0) { *Nu = 0.0; return; }

    *Nu = 0.943 * pow(arg, 0.25);
}

void rohsenow_pool_boiling(double q_flux, double rho_L, double rho_v,
                            double h_fg, double cp_L, double mu_L,
                            double Pr_L, double sigma, double C_sf,
                            double n_exp, double *delta_T) {
    if (!delta_T || q_flux <= 0.0 || mu_L <= 0.0
        || h_fg <= 0.0 || sigma <= 0.0) {
        if (delta_T) *delta_T = 0.0;
        return;
    }

    double drho = rho_L - rho_v;
    if (drho <= 0.0) { *delta_T = 0.0; return; }

    double g = 9.81;
    /* Rohsenow correlation rearranged to solve for delta_T:
     * delta_T = (q/(mu_L*h_fg))^(1/3) * sqrt(sigma/(g*drho)) * (C_sf*h_fg*Pr_L^n/cp_L)
     */
    double term1 = pow(q_flux / (mu_L * h_fg), 1.0/3.0);
    double term2 = sqrt(sigma / (g * drho));
    double term3 = C_sf * h_fg * pow(Pr_L, n_exp) / cp_L;

    *delta_T = term1 * term2 * term3;
}

/* ============================================================================
 * L5: Sherwood Number Correlations
 * ============================================================================ */

void sherwood_falling_film(double Re_film, double Sc, double L_over_d,
                            double *Sh) {
    if (!Sh || Re_film <= 0.0 || Sc <= 0.0 || L_over_d <= 0.0) {
        if (Sh) *Sh = 0.0;
        return;
    }

    /* Short contact time (< 0.1s): penetration theory */
    /* Sh = 2/sqrt(pi) * sqrt(Re * Sc * d/L) */
    double Gz_mass = Re_film * Sc / L_over_d;  /* mass Graetz number */

    if (Gz_mass > 100.0) {
        /* Short contact time (entry region) */
        *Sh = (2.0 / sqrt(M_PI)) * sqrt(Gz_mass);
    } else if (Gz_mass < 10.0) {
        /* Long contact time (fully developed) */
        *Sh = 3.41;
    } else {
        /* Interpolate */
        double Sh_short = (2.0 / sqrt(M_PI)) * sqrt(Gz_mass);
        double Sh_long  = 3.41;
        double w = (Gz_mass - 10.0) / 90.0;  /* weight 0..1 for short contact */
        *Sh = w * Sh_short + (1.0 - w) * Sh_long;
    }
}

void sherwood_single_sphere(double Re_sphere, double Sc, double *Sh) {
    if (!Sh || Re_sphere < 0.0 || Sc <= 0.0) {
        if (Sh) *Sh = 0.0;
        return;
    }

    /* Froessling (1938): Sh = 2.0 + 0.552 * Re^(1/2) * Sc^(1/3) */
    *Sh = 2.0 + 0.552 * sqrt(Re_sphere) * pow(Sc, 1.0/3.0);

    /* Limiting case: Re -> 0 yields Sh -> 2 (pure diffusion) */
}

void sherwood_packed_bed(double Re_p, double Sc, double epsilon, double *Sh) {
    if (!Sh || Re_p <= 0.0 || Sc <= 0.0 || epsilon <= 0.0 || epsilon >= 1.0) {
        if (Sh) *Sh = 0.0;
        return;
    }

    /* Wakao & Funazkri (1978) */
    *Sh = 2.0 + 1.1 * pow(Re_p, 0.6) * pow(Sc, 1.0/3.0);

    /* Void fraction effect on fluid-particle mass transfer:
     * j_D = (Sh/(Re*Sc^(1/3))), related by Chilton-Colburn analogy */
}

void sherwood_bubble_column(double d_b, double g, double rho_L,
                             double sigma, double nu_L, double D_L,
                             double *Sh) {
    if (!Sh || d_b <= 0.0 || nu_L <= 0.0 || D_L <= 0.0) {
        if (Sh) *Sh = 0.0;
        return;
    }

    double Sc = nu_L / D_L;
    double Bo = g * d_b * d_b * rho_L / sigma;  /* Bond number */
    double Ga = g * d_b * d_b * d_b / (nu_L * nu_L);  /* Galileo number */

    if (Bo <= 0.0 || Ga <= 0.0) { *Sh = 0.0; return; }

    /* Akita & Yoshida (1974) */
    *Sh = 0.5 * sqrt(Sc) * pow(Bo, 3.0/8.0) * pow(Ga, 0.25);
}

/* ============================================================================
 * L2: Interfacial Resistance and Enhancement
 * ============================================================================ */

void enhancement_factor_instantaneous(double D_A, double D_B,
                                       double C_B_bulk, double C_Ai,
                                       double nu, double *E_inf) {
    if (!E_inf) return;
    if (D_A <= 0.0 || C_Ai <= 0.0) {
        *E_inf = 1.0;  /* no enhancement possible */
        return;
    }

    /* Film theory: E_inf = 1 + (D_B/D_A)^(n) * (C_B_bulk/(nu*C_Ai))
     * For equal diffusivities and nu=1:
     * E_inf = 1 + D_B*C_B_bulk/(D_A*C_Ai)
     */
    double ratio = D_B * C_B_bulk / (nu * D_A * C_Ai);
    *E_inf = 1.0 + ratio;
    if (*E_inf < 1.0) *E_inf = 1.0;  /* enhancement >= 1 */
}

void hatta_number(double k_rxn, double D_A, double k_L, double *Ha) {
    if (!Ha || k_L <= 0.0) {
        if (Ha) *Ha = 0.0;
        return;
    }
    /* Ha = sqrt(k_rxn * D_A) / k_L */
    *Ha = sqrt(k_rxn * D_A) / k_L;
}

void overall_coefficient_with_reaction(double K_overall, double E,
                                        double *K_with_rxn) {
    if (!K_with_rxn) return;
    *K_with_rxn = K_overall * E;
}

/* ============================================================================
 * L3: Interfacial Area Estimation
 * ============================================================================ */

void interfacial_area_packed_column(double a_t, double sigma_c,
                                     double sigma_L, double Re_L,
                                     double Fr_L, double We_L, double *a_eff) {
    if (!a_eff || a_t <= 0.0) {
        if (a_eff) *a_eff = 0.0;
        return;
    }
    if (sigma_L <= 0.0 || Re_L <= 0.0) {
        *a_eff = 0.0;
        return;
    }

    /* Onda et al. (1968) correlation */
    double sigma_ratio = sigma_c / sigma_L;
    double exponent = -1.45 * pow(sigma_ratio, 0.75)
                      * pow(Re_L, 0.1)
                      * pow(fabs(Fr_L) + 1e-10, -0.05)
                      * pow(We_L + 1e-10, 0.2);

    /* Clamp exponent to avoid underflow */
    if (exponent < -50.0) exponent = -50.0;
    if (exponent > 50.0) exponent = 50.0;

    *a_eff = a_t * (1.0 - exp(exponent));
}

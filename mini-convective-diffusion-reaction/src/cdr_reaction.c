/**
 * @file cdr_reaction.c
 * @brief Implementation of chemical reaction kinetics.
 *
 * Arrhenius law, power-law rate expressions, Michaelis-Menten,
 * Langmuir-Hinshelwood, and kinetic analysis (half-life, conversion,
 * selectivity, yield).
 *
 * Reference: Levenspiel "Chemical Reaction Engineering" (1999), Ch. 2-7
 * Reference: Fogler "Elements of Chemical Reaction Engineering" (2016)
 */

#include "cdr_reaction.h"
#include <math.h>

#ifndef R_GAS
#define R_GAS 8.314462618  /* J/(mol.K) */
#endif

/* ---------------------------------------------------------------------------
 * L4: Arrhenius Rate Law
 * ------------------------------------------------------------------------- */

double cdr_arrhenius_rate(double A, double Ea, double T)
{
    /* k(T) = A * exp(-Ea / (R * T)) */
    if (T <= 0.0 || A < 0.0) {
        return 0.0;
    }
    if (Ea < 0.0) {
        Ea = 0.0;  /* Negative Ea (rare, e.g., radical recombination) allowed */
    }
    return A * exp(-Ea / (R_GAS * T));
}

double cdr_arrhenius_from_ref(double k_ref, double Ea,
                               double T_ref, double T)
{
    /*
     * k(T) = k_ref * exp( (Ea/R) * (1/T_ref - 1/T) )
     *
     * More numerically stable than recomputing A.
     */
    if (T_ref <= 0.0 || T <= 0.0 || k_ref < 0.0 || Ea < 0.0) {
        return 0.0;
    }

    double exponent = (Ea / R_GAS) * (1.0 / T_ref - 1.0 / T);
    return k_ref * exp(exponent);
}

/* ---------------------------------------------------------------------------
 * L4: Rate Laws
 * ------------------------------------------------------------------------- */

double cdr_power_law_rate(double k, double C_A, double C_B,
                           double alpha, double beta)
{
    /*
     * r = k * C_A^alpha * C_B^beta
     *
     * Handles special case of concentration <= 0.
     */
    if (k < 0.0) {
        return 0.0;
    }
    if (C_A <= 0.0 && alpha > 0.0) {
        return 0.0;
    }
    if (C_B <= 0.0 && beta > 0.0) {
        return 0.0;
    }

    double rate = k;
    if (alpha != 0.0) {
        rate *= pow(C_A, alpha);
    }
    if (beta != 0.0) {
        rate *= pow(C_B, beta);
    }
    return rate;
}

double cdr_reversible_rate(const PowerLawRate *law,
                            double C_A, double C_B,
                            double C_C, double C_D)
{
    /*
     * r = k_f * C_A^alpha * C_B^beta - k_r * C_C^gamma * C_D^delta
     */
    if (law == NULL) {
        return 0.0;
    }

    double r_forward = 0.0;
    double r_reverse = 0.0;

    if (law->k_forward > 0.0) {
        r_forward = law->k_forward;
        if (law->alpha != 0.0 && C_A > 0.0) {
            r_forward *= pow(C_A, law->alpha);
        }
        if (law->beta != 0.0 && C_B > 0.0) {
            r_forward *= pow(C_B, law->beta);
        }
    }

    if (law->is_reversible && law->k_reverse > 0.0) {
        r_reverse = law->k_reverse;
        if (law->gamma != 0.0 && C_C > 0.0) {
            r_reverse *= pow(C_C, law->gamma);
        }
        if (law->delta != 0.0 && C_D > 0.0) {
            r_reverse *= pow(C_D, law->delta);
        }
    }

    return r_forward - r_reverse;
}

double cdr_michaelis_menten_rate(const MichaelisMentenParams *params,
                                  double S)
{
    /*
     * r = V_max * [S] / (K_m + [S])
     *
     * For S >> K_m: r ~ V_max (zero-order)
     * For S << K_m: r ~ (V_max/K_m) * S (first-order)
     */
    if (params == NULL || S < 0.0) {
        return 0.0;
    }
    if (params->K_m <= 0.0) {
        return 0.0;
    }

    return params->V_max * S / (params->K_m + S);
}

double cdr_michaelis_menten_inhibited(const MichaelisMentenParams *params,
                                       double S, double I, double K_i)
{
    /*
     * Competitive inhibition:
     * r = V_max * [S] / (K_m * (1 + [I]/K_i) + [S])
     */
    if (params == NULL || S < 0.0 || I < 0.0 || K_i <= 0.0) {
        return 0.0;
    }
    if (params->K_m <= 0.0) {
        return 0.0;
    }

    double apparent_Km = params->K_m * (1.0 + I / K_i);
    return params->V_max * S / (apparent_Km + S);
}

double cdr_langmuir_hinshelwood_uni(const LangmuirHinshelwoodParams *params,
                                     double P_A,
                                     const double *P_inhib,
                                     const double *K_inhib,
                                     size_t n_inhib)
{
    /*
     * r = k * K_A * P_A / (1 + K_A*P_A + sum(K_j * P_j))
     */
    if (params == NULL || P_A < 0.0) {
        return 0.0;
    }

    double numerator = params->k_surface * params->K_ads_A * P_A;

    double denominator = 1.0 + params->K_ads_A * P_A;

    /* Add inhibition terms */
    if (P_inhib != NULL && K_inhib != NULL) {
        for (size_t j = 0; j < n_inhib; j++) {
            if (K_inhib[j] > 0.0 && P_inhib[j] > 0.0) {
                denominator += K_inhib[j] * P_inhib[j];
            }
        }
    }

    if (denominator <= 0.0) {
        return 0.0;
    }

    return numerator / denominator;
}

double cdr_langmuir_hinshelwood_bi(const LangmuirHinshelwoodParams *params,
                                    double P_A, double P_B)
{
    /*
     * r = k * K_A * K_B * P_A * P_B / (1 + K_A*P_A + K_B*P_B)^2
     *
     * For bimolecular surface reaction with competitive adsorption.
     */
    if (params == NULL || P_A < 0.0 || P_B < 0.0) {
        return 0.0;
    }

    double numerator = params->k_surface * params->K_ads_A *
                       params->K_ads_B * P_A * P_B;

    double denominator = 1.0 + params->K_ads_A * P_A + params->K_ads_B * P_B;
    denominator = denominator * denominator;

    if (denominator <= 0.0) {
        return 0.0;
    }

    return numerator / denominator;
}

/* ---------------------------------------------------------------------------
 * L5: Kinetic Analysis Methods
 * ------------------------------------------------------------------------- */

double cdr_half_life_first_order(double k)
{
    /* t_1/2 = ln(2) / k */
    if (k <= 0.0) {
        return INFINITY;
    }
    return log(2.0) / k;
}

double cdr_half_life_nth_order(double k, double n, double C_A0)
{
    /*
     * t_1/2 = (2^(n-1) - 1) / (k * (n-1) * C_A0^(n-1))   for n != 1
     */
    if (k <= 0.0 || C_A0 <= 0.0) {
        return INFINITY;
    }

    if (fabs(n - 1.0) < 1e-10) {
        /* Degenerate to first-order */
        return cdr_half_life_first_order(k);
    }

    double numerator = pow(2.0, n - 1.0) - 1.0;
    double denominator = k * (n - 1.0) * pow(C_A0, n - 1.0);

    if (fabs(denominator) < 1e-15) {
        return INFINITY;
    }

    return numerator / denominator;
}

double cdr_conversion_first_order(double k, double t)
{
    /* X = 1 - exp(-k*t) */
    if (k <= 0.0 || t < 0.0) {
        return 0.0;
    }
    return 1.0 - exp(-k * t);
}

double cdr_time_to_conversion_first_order(double k, double X)
{
    /* t = -ln(1 - X) / k */
    if (k <= 0.0 || X >= 1.0 || X < 0.0) {
        return INFINITY;
    }
    if (X <= 0.0) {
        return 0.0;
    }
    return -log(1.0 - X) / k;
}

double cdr_concentration_first_order(double C_A0, double k, double t)
{
    /* C_A(t) = C_A0 * exp(-k*t) */
    if (C_A0 < 0.0 || k < 0.0 || t < 0.0) {
        return (C_A0 < 0.0) ? 0.0 : C_A0;
    }
    return C_A0 * exp(-k * t);
}

double cdr_concentration_second_order(double C_A0, double k, double t)
{
    /* 1/C_A = 1/C_A0 + k*t */
    if (C_A0 <= 0.0 || k <= 0.0 || t < 0.0) {
        return (C_A0 <= 0.0) ? 0.0 : C_A0;
    }
    return 1.0 / (1.0 / C_A0 + k * t);
}

double cdr_concentration_zero_order(double C_A0, double k, double t)
{
    /* C_A = C_A0 - k*t (until depleted) */
    if (C_A0 < 0.0 || k < 0.0 || t < 0.0) {
        return (C_A0 < 0.0) ? 0.0 : C_A0;
    }
    double C = C_A0 - k * t;
    return (C > 0.0) ? C : 0.0;
}

void cdr_series_reaction_concentrations(double C_A0, double k1, double k2,
                                         double t, double C_out[3])
{
    /*
     * A --k1--> B --k2--> C
     *
     * C_A = C_A0 * exp(-k1*t)
     * C_B = C_A0 * k1/(k2-k1) * (exp(-k1*t) - exp(-k2*t))   [k1 != k2]
     * C_B = C_A0 * k1 * t * exp(-k1*t)                       [k1 == k2]
     * C_C = C_A0 - C_A - C_B
     */
    if (C_out == NULL || C_A0 <= 0.0) {
        return;
    }

    if (t < 0.0) {
        t = 0.0;
    }

    double exp_k1t = exp(-k1 * t);
    C_out[0] = C_A0 * exp_k1t;

    if (fabs(k1 - k2) < 1e-12) {
        /* Degenerate case: k1 == k2 */
        C_out[1] = C_A0 * k1 * t * exp_k1t;
    } else {
        double exp_k2t = exp(-k2 * t);
        C_out[1] = C_A0 * k1 / (k2 - k1) * (exp_k1t - exp_k2t);
    }

    C_out[2] = C_A0 - C_out[0] - C_out[1];
}

double cdr_activation_energy_from_two_points(double k1, double T1,
                                              double k2, double T2)
{
    /*
     * Ea = R * ln(k2/k1) / (1/T1 - 1/T2)
     */
    if (k1 <= 0.0 || k2 <= 0.0 || T1 <= 0.0 || T2 <= 0.0) {
        return 0.0;
    }
    if (fabs(1.0/T1 - 1.0/T2) < 1e-15) {
        return 0.0;  /* Same temperature, infinite data needed */
    }
    return R_GAS * log(k2 / k1) / (1.0 / T1 - 1.0 / T2);
}

double cdr_selectivity_parallel(double k1, double k2, double n1, double n2,
                                 double C_A)
{
    /*
     * A -> B (k1, n1)  desired
     * A -> C (k2, n2)  undesired
     *
     * S_B/C = (k1/k2) * C_A^(n1 - n2)
     */
    if (k2 <= 0.0) {
        return INFINITY;  /* No undesired product formation */
    }
    if (k1 <= 0.0) {
        return 0.0;  /* No desired product formation */
    }
    if (C_A <= 0.0 && (n1 - n2) < 0.0) {
        return 0.0;
    }

    double power = n1 - n2;
    return (k1 / k2) * pow(C_A, power);
}

double cdr_yield_intermediate_series(double k1, double k2,
                                      double C_A0, double t)
{
    /*
     * Y_B = C_B / C_A0
     *
     * C_B = C_A0 * k1/(k2-k1) * (exp(-k1*t) - exp(-k2*t))  [k1 != k2]
     */
    if (C_A0 <= 0.0 || k1 <= 0.0 || k2 <= 0.0 || t < 0.0) {
        return 0.0;
    }

    if (fabs(k1 - k2) < 1e-12) {
        return k1 * t * exp(-k1 * t);
    }

    double exp_k1t = exp(-k1 * t);
    double exp_k2t = exp(-k2 * t);

    return k1 / (k2 - k1) * (exp_k1t - exp_k2t);
}

double cdr_optimal_time_series(double k1, double k2)
{
    /*
     * t_opt = ln(k2/k1) / (k2 - k1)
     *
     * This maximizes the intermediate B concentration.
     */
    if (k1 <= 0.0 || k2 <= 0.0) {
        return INFINITY;
    }
    if (fabs(k1 - k2) < 1e-12) {
        /* k1 == k2: t_opt = 1/k1 */
        return 1.0 / k1;
    }

    return log(k2 / k1) / (k2 - k1);
}

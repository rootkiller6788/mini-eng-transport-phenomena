/**
 * @file cdr_effectiveness.c
 * @brief Implementation of catalyst effectiveness factor and pore diffusion.
 *
 * Effectiveness factors for slab, cylinder, sphere geometries,
 * generalized Thiele modulus approach, diffusion limitation diagnostics
 * (Weisz-Prater, Mears criteria), and washcoat modeling.
 *
 * Reference: Bird, Stewart, Lightfoot "Transport Phenomena" (2007)
 * Reference: Fogler "Elements of Chemical Reaction Engineering" (2016), Ch. 12, 15
 */

#include "cdr_effectiveness.h"
#include <math.h>
#include <string.h>

/* ---------------------------------------------------------------------------
 * Helper: modified Bessel functions I_0 and I_1 approximations
 * ------------------------------------------------------------------------- */

/**
 * @brief Modified Bessel function I_0(x) polynomial approximation.
 *
 * Abramowitz & Stegun 9.8.1 (|x| <= 3.75) and 9.8.2 (x > 3.75).
 * Absolute error < 1.6e-7.
 */
static double bessel_I0(double x)
{
    double ax = fabs(x);
    double result;

    if (ax < 3.75) {
        double t = x / 3.75;
        double t2 = t * t;
        result = 1.0 + 3.5156229 * t2
                    + 3.0899424 * t2 * t2
                    + 1.2067492 * t2 * t2 * t2
                    + 0.2659732 * t2 * t2 * t2 * t2
                    + 0.0360768 * t2 * t2 * t2 * t2 * t2
                    + 0.0045813 * t2 * t2 * t2 * t2 * t2 * t2;
    } else {
        double t = 3.75 / ax;
        result = (exp(ax) / sqrt(ax)) *
                 (0.39894228 + 0.01328592 * t
                             + 0.00225319 * t * t
                             - 0.00157565 * t * t * t
                             + 0.00916281 * t * t * t * t
                             - 0.02057706 * t * t * t * t * t
                             + 0.02635537 * t * t * t * t * t * t
                             - 0.01647633 * t * t * t * t * t * t * t
                             + 0.00392377 * t * t * t * t * t * t * t * t);
    }
    return result;
}

/**
 * @brief Modified Bessel function I_1(x) polynomial approximation.
 *
 * Abramowitz & Stegun 9.8.3 (|x| <= 3.75) and 9.8.4 (x > 3.75).
 * Absolute error < 1.6e-7.
 */
static double bessel_I1(double x)
{
    double ax = fabs(x);
    double result;

    if (ax < 3.75) {
        double t = x / 3.75;
        double t2 = t * t;
        result = x * (0.5 + 0.87890594 * t2
                         + 0.51498869 * t2 * t2
                         + 0.15084934 * t2 * t2 * t2
                         + 0.02658733 * t2 * t2 * t2 * t2
                         + 0.00301532 * t2 * t2 * t2 * t2 * t2
                         + 0.00032411 * t2 * t2 * t2 * t2 * t2 * t2);
    } else {
        double t = 3.75 / ax;
        double sign = (x > 0.0) ? 1.0 : -1.0;
        result = sign * (exp(ax) / sqrt(ax)) *
                 (0.39894228 - 0.03988024 * t
                             - 0.00362018 * t * t
                             + 0.00163801 * t * t * t
                             - 0.01031555 * t * t * t * t
                             + 0.02282967 * t * t * t * t * t
                             - 0.02895312 * t * t * t * t * t * t
                             + 0.01787654 * t * t * t * t * t * t * t
                             - 0.00420059 * t * t * t * t * t * t * t * t);
    }
    return result;
}

/* ---------------------------------------------------------------------------
 * L4: Effectiveness Factor for Standard Geometries (1st order)
 * ------------------------------------------------------------------------- */

double cdr_effectiveness_slab(double phi)
{
    /*
     * eta = tanh(phi) / phi
     *
     * phi = L * sqrt(k/D_eff) where L = half-thickness.
     *
     * Limits:
     *   phi -> 0: eta -> 1
     *   phi -> inf: eta -> 1/phi
     */
    if (phi <= 0.0) {
        return 1.0;  /* No diffusion limitation */
    }

    /* For very small phi, use series: tanh(phi) ~ phi - phi^3/3 */
    if (phi < 0.01) {
        return 1.0 - phi * phi / 3.0;
    }

    /* For very large phi, tanh(phi) ~ 1 */
    if (phi > 50.0) {
        return 1.0 / phi;
    }

    return tanh(phi) / phi;
}

double cdr_effectiveness_cylinder(double phi)
{
    /*
     * eta = 2 * I_1(phi) / (phi * I_0(phi))
     *
     * where L_c = R/2.
     *
     * Limits:
     *   phi -> 0: eta -> 1
     *   phi -> inf: eta -> 2/phi
     */
    if (phi <= 0.0) {
        return 1.0;
    }

    /* For small phi, use series expansion: */
    if (phi < 0.01) {
        double phi2 = phi * phi;
        return 1.0 - phi2 / 8.0 + phi2 * phi2 / 192.0;
    }

    /* For large phi, I_1/I_0 -> 1 */
    if (phi > 50.0) {
        return 2.0 / phi;
    }

    double I0 = bessel_I0(phi);
    double I1 = bessel_I1(phi);

    if (I0 <= 0.0) {
        return 2.0 / phi;  /* Fallback for large phi */
    }

    return 2.0 * I1 / (phi * I0);
}

double cdr_effectiveness_sphere(double phi)
{
    /*
     * eta = (3/phi) * (1/tanh(phi) - 1/phi)
     *
     * where L_c = R/3.
     *
     * Limits:
     *   phi -> 0: eta -> 1
     *   phi -> inf: eta -> 3/phi
     */
    if (phi <= 0.0) {
        return 1.0;
    }

    /* For small phi, use Taylor expansion: */
    if (phi < 0.01) {
        double phi2 = phi * phi;
        return 1.0 - phi2 / 15.0 + 2.0 * phi2 * phi2 / 315.0;
    }

    /* For large phi: */
    if (phi > 50.0) {
        return 3.0 / phi;
    }

    double tanh_phi = tanh(phi);
    return (3.0 / phi) * (1.0 / tanh_phi - 1.0 / phi);
}

/* ---------------------------------------------------------------------------
 * L4: Generalized Effectiveness Factor
 * ------------------------------------------------------------------------- */

double cdr_effectiveness_general(double phi, CatalystGeometry geometry)
{
    /*
     * Generalized method: compute eta based on geometry.
     * Uses the same phi definition (based on L_c = V_p/S_p).
     */
    switch (geometry) {
        case GEOM_SLAB:
            return cdr_effectiveness_slab(phi);
        case GEOM_CYLINDER:
            return cdr_effectiveness_cylinder(phi);
        case GEOM_SPHERE:
            return cdr_effectiveness_sphere(phi);
        case GEOM_ARBITRARY_3D:
        default:
            /* For arbitrary geometry, use sphere as conservative estimate */
            return cdr_effectiveness_sphere(phi);
    }
}

double cdr_effectiveness_nth_order(double k_n, double n, double C_s,
                                    double D_eff, double L_c)
{
    /*
     * Generalized Thiele modulus for nth-order:
     *
     * phi_n = L_c * sqrt( (n+1) * k_n * C_s^(n-1) / (2 * D_eff) )
     *
     * This comes from linearizing the rate expression around C_s.
     *
     * Then eta is approximated using the slab formula:
     * eta ? tanh(phi_n) / phi_n
     *
     * Note: This is exact for n=1 and a reasonable approximation for
     * other orders.
     */
    if (k_n <= 0.0 || C_s <= 0.0 || D_eff <= 0.0 || L_c <= 0.0) {
        return 1.0;
    }
    if (n < 0.0) {
        return 1.0;  /* Negative order (inhibited) ? approximate */
    }
    if (fabs(n - 1.0) < 1e-10) {
        /* First-order: exact formula */
        double phi = cdr_thiele_modulus(L_c, k_n, D_eff);
        return cdr_effectiveness_slab(phi);
    }

    /* Generalized modulus (Bischoff, 1965) */
    double factor = (n + 1.0) / 2.0;
    double eff_rate = k_n * pow(C_s, n - 1.0);
    double phi_n = L_c * sqrt(factor * eff_rate / D_eff);

    /* Use slab formula as a practical approximation */
    return cdr_effectiveness_slab(phi_n);
}

double cdr_effectiveness_michaelis_menten(double V_max, double K_m,
                                           double S_s, double D_eff,
                                           double L_c)
{
    /*
     * Approximate effectiveness for MM kinetics using the general modulus
     * approach (Atkinson & Davies, 1974).
     *
     * phi_MM = L_c * sqrt(V_max / (D_eff * K_m)) / sqrt(1 + S_s/K_m)
     *
     * eta ? tanh(phi_MM) / phi_MM  (slab approximation)
     */
    if (V_max <= 0.0 || K_m <= 0.0 || S_s < 0.0 || D_eff <= 0.0 || L_c <= 0.0) {
        return 1.0;
    }

    double k_eff_zero = V_max / K_m;  /* Effective 1st-order rate constant at low S */
    double sat_factor  = 1.0 + S_s / K_m;

    double phi_MM = L_c * sqrt(k_eff_zero / D_eff) / sqrt(sat_factor);

    return cdr_effectiveness_slab(phi_MM);
}

/* ---------------------------------------------------------------------------
 * L5: Diffusion Limitation Diagnostics
 * ------------------------------------------------------------------------- */

double cdr_weisz_prater_criterion(double r_obs, double L_c,
                                   double D_eff, double C_s)
{
    /*
     * Weisz-Prater parameter:
     *
     * C_WP = r_obs * L_c^2 / (D_eff * C_s)
     *
     * Interpretation:
     *   C_WP << 1: No internal concentration gradient (eta ~ 1)
     *   C_WP >> 1: Significant internal diffusion limitation
     *
     * Note: C_WP = eta * phi^2 for first-order reaction.
     */
    if (D_eff <= 0.0 || C_s <= 0.0 || L_c <= 0.0 || r_obs < 0.0) {
        return 0.0;
    }

    return r_obs * L_c * L_c / (D_eff * C_s);
}

double cdr_mears_external_criterion(double r_obs, double L_c, double n,
                                     double k_c, double C_bulk)
{
    /*
     * Mears criterion for external mass transfer:
     *
     * C_M = r_obs * L_c * n / (k_c * C_bulk)
     *
     * If C_M < 0.15: External mass transfer resistance < 5% error
     * If C_M > 0.15: External limitations are significant
     *
     * r_obs per unit catalyst volume, L_c = V_p/S_p
     */
    if (k_c <= 0.0 || C_bulk <= 0.0 || r_obs < 0.0) {
        return 0.0;
    }

    return r_obs * L_c * n / (k_c * C_bulk);
}

void cdr_diagnose_limitations(double r_obs, double L_c, double D_eff,
                               double C_s, double n, double k_c,
                               double C_bulk, DiffusionLimitationDiagnosis *diag)
{
    /*
     * Full diagnostic combining Weisz-Prater and Mears criteria.
     */
    if (diag == NULL) {
        return;
    }

    memset(diag, 0, sizeof(DiffusionLimitationDiagnosis));

    /* Weisz-Prater for internal limitation */
    diag->Weisz_Prater = cdr_weisz_prater_criterion(r_obs, L_c, D_eff, C_s);

    /* C_WP > 1 indicates significant internal gradient */
    diag->is_internal_limited = (diag->Weisz_Prater > 1.0) ? 1 : 0;

    /* Mears for external limitation */
    diag->Mears_external = cdr_mears_external_criterion(r_obs, L_c, n,
                                                         k_c, C_bulk);
    /* C_M > 0.15 indicates significant external resistance */
    diag->is_external_limited = (diag->Mears_external > 0.15) ? 1 : 0;

    /* Carberry number: Ca = r_obs / (k_c * a * C_bulk) */
    if (k_c > 0.0 && C_bulk > 0.0 && L_c > 0.0) {
        /* For sphere: a = 3/R ? a * L_c ? 3 * (L_c being V_p/S_p = R/3) */
        double a_approx = 1.0 / L_c;  /* a = S_p/V_p = 1/L_c */
        diag->Carberry = r_obs / (k_c * a_approx * C_bulk);
    } else {
        diag->Carberry = 0.0;
    }
}

/* ---------------------------------------------------------------------------
 * L5: Observed Rate and Kinetics Recovery
 * ------------------------------------------------------------------------- */

double cdr_observed_rate(double r_intrinsic, double eta)
{
    /* r_obs = eta * r_intrinsic */
    if (r_intrinsic < 0.0 || eta < 0.0 || eta > 1.0) {
        return 0.0;
    }
    return eta * r_intrinsic;
}

double cdr_intrinsic_activation_energy(double E_obs, double E_diff)
{
    /*
     * For strong pore diffusion limitation in a sphere (1st order):
     *
     * E_obs = (E_true + E_diff) / 2
     *
     * Therefore: E_true = 2 * E_obs - E_diff
     *
     * This recovers the true (intrinsic) activation energy from the
     * apparent one measured under diffusion-limited conditions.
     *
     * E_diff for Knudsen diffusion: E_diff ~ 0 (temperature dependence
     * only through sqrt(T), so E_diff ~ 0.5*R*T ~ 2-4 kJ/mol).
     *
     * E_diff for molecular diffusion: ~10-20 kJ/mol.
     */
    if (E_obs < 0.0 || E_diff < 0.0) {
        return 0.0;
    }

    double E_true = 2.0 * E_obs - E_diff;
    return (E_true > 0.0) ? E_true : 0.0;
}

double cdr_apparent_reaction_order(double n_true)
{
    /*
     * Under strong internal diffusion limitation:
     *
     * n_app = (n_true + 1) / 2
     *
     * This means:
     *   1st-order true (n=1) ? n_app = 1 (order preserved!)
     *   2nd-order true (n=2) ? n_app = 1.5
     *   0th-order true (n=0) ? n_app = 0.5
     *
     * The diffusion disguises the true kinetics.
     */
    return (n_true + 1.0) / 2.0;
}

double cdr_pellet_concentration_profile(double r, double R, double C_s,
                                         double phi)
{
    /*
     * Spherical catalyst pellet, 1st-order irreversible reaction.
     *
     * C(r) = C_s * (R/r) * sinh(phi * r / R) / sinh(phi)
     *
     * where phi = R * sqrt(k/D_eff) based on R (V_p/S_p = R/3 used elsewhere).
     *
     * Note: Here phi is based on R (not L_c). The standard Thiele modulus
     * for a sphere is phi_s = R * sqrt(k/D_eff) / 3.
     *
     * If using the standard phi_s (based on R/3), then the profile is:
     *   C(r) = C_s * (R/r) * sinh(3*phi_s * r/R) / sinh(3*phi_s)
     *
     * This function expects phi based on R: phi_R = R*sqrt(k/D_eff) = 3*phi_s.
     */
    if (R <= 0.0 || r < 0.0 || r > R || C_s <= 0.0 || phi < 0.0) {
        return 0.0;
    }
    if (r <= 0.0) {
        /* Center: C(0) = C_s * phi / sinh(phi) (limit as r->0) */
        if (phi < 1e-10) {
            return C_s;
        }
        return C_s * phi / sinh(phi);
    }
    if (r >= R) {
        return C_s;
    }
    if (phi < 1e-10) {
        return C_s;  /* No reaction, uniform concentration */
    }

    double z = phi * r / R;
    return C_s * (R / r) * sinh(z) / sinh(phi);
}

double cdr_optimal_pellet_radius(double k, double D_eff, double target_eta)
{
    /*
     * Find optimal pellet radius for first-order reaction in a sphere.
     *
     * eta = (3/phi) * (1/tanh(phi) - 1/phi)
     *
     * For a target eta, we solve for phi, then:
     *   R_opt = 3 * phi * sqrt(D_eff/k)
     *
     * Since the equation is transcendental, we use a precomputed lookup
     * approach: phi is approximately 1/eta for eta << 1.
     *
     * For the sphere, eta(phi) is monotonic decreasing.
     * We use bisection to find phi for the target eta.
     */
    if (k <= 0.0 || D_eff <= 0.0 || target_eta <= 0.0 || target_eta > 1.0) {
        return 0.0;
    }

    /* For eta ~ 1, phi ~ 0 ? use very small pellet */
    if (target_eta > 0.999) {
        return 0.001 * sqrt(D_eff / k);
    }

    /* Bisection to find phi */
    double phi_low  = 0.0;
    double phi_high = 100.0;  /* Sufficient for essentially any target */
    double eta_mid, phi_mid;

    for (int iter = 0; iter < 100; iter++) {
        phi_mid = (phi_low + phi_high) / 2.0;
        eta_mid = cdr_effectiveness_sphere(phi_mid);

        if (fabs(eta_mid - target_eta) < 1e-8) {
            break;
        }

        if (eta_mid > target_eta) {
            phi_low = phi_mid;
        } else {
            phi_high = phi_mid;
        }
    }

    /* phi is based on L_c = R/3 for sphere */
    /* phi = L_c * sqrt(k/D_eff) = R/3 * sqrt(k/D_eff) */
    /* R = 3 * phi * sqrt(D_eff/k) */
    return 3.0 * phi_mid * sqrt(D_eff / k);
}

double cdr_washcoat_effectiveness(double L_w, double k, double D_eff_w)
{
    /*
     * Washcoat is a thin layer on a monolith wall.
     * Since it's thin relative to curvature, slab geometry applies:
     *
     * eta = tanh(phi_w) / phi_w
     *
     * where phi_w = L_w * sqrt(k/D_eff,w)
     *
     * L_w = washcoat thickness [m]
     */
    if (L_w <= 0.0) {
        return 1.0;
    }

    double phi_w = L_w * sqrt(k / D_eff_w);
    return cdr_effectiveness_slab(phi_w);
}

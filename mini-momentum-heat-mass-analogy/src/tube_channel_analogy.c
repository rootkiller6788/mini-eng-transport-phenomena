/**
 * @file tube_channel_analogy.c
 * @brief Analogy applied to internal flows: pipes, ducts, channels
 *
 * Implements friction factor correlations, heat/mass transfer correlations,
 * and analogy-based prediction methods for internal flows.
 *
 * Reference: Bird, Stewart, Lightfoot (2007), Transport Phenomena, Ch. 6, 13, 14.
 *            Dittus & Boelter (1930), Univ. Calif. Publ. Eng., 2, 443.
 *            Gnielinski (1976), Int. Chem. Eng., 16, 359-368.
 *            Colebrook (1939), J. Inst. Civil Eng., 11, 133-156.
 */

#include "tube_channel_analogy.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

/* ─── L2: Friction Factors ──────────────────────────────────────── */

double darcy_friction_laminar(double Re_D)
{
    if (Re_D <= 0.0) return -1.0;
    return 64.0 / Re_D;
}

double fanning_friction_laminar(double Re_D)
{
    if (Re_D <= 0.0) return -1.0;
    return 16.0 / Re_D;
}

double blasius_friction_turbulent(double Re_D)
{
    if (Re_D <= 0.0) return -1.0;
    return 0.046 * pow(Re_D, -0.2);
}

/**
 * Colebrook-White equation for turbulent friction factor.
 *
 * 1/√f = -2.0·log₁₀(ε/(3.7·D) + 2.51/(Re·√f))
 *
 * Solved using Newton's method with initial guess from Blasius.
 * Convergence tolerance: 1e-8, max 50 iterations.
 *
 * Complexity: O(iterations) ≈ O(10-20)
 */
double colebrook_friction(double Re_D, double roughness, double D)
{
    if (Re_D <= 0.0 || D <= 0.0) return -1.0;

    /* Initial guess: Blasius */
    double f = blasius_friction_turbulent(Re_D);
    if (f < 1e-10) f = 0.02;

    double eps_D = roughness / D;  /* relative roughness */

    /* Newton's method: solve g(f) = 1/√f + 2·log₁₀(...) = 0 */
    for (int iter = 0; iter < 50; iter++) {
        double sqrt_f = sqrt(f);
        double log_arg = eps_D / 3.7 + 2.51 / (Re_D * sqrt_f);
        if (log_arg <= 0.0) log_arg = 1e-15;

        double g = 1.0 / sqrt_f + 2.0 * log10(log_arg);

        /* Derivative: dg/df = -1/(2·f^(3/2)) - (2/(ln10))·(2.51/(Re·f^(3/2))) / log_arg */
        double dg_df = -0.5 / (f * sqrt_f) -
                        (2.0 / log(10.0)) * (2.51 / (Re_D * 2.0 * f * sqrt_f)) / log_arg;

        double df = -g / dg_df;
        f += df;

        /* Limit f to physically meaningful range */
        if (f < 1e-8) f = 1e-8;
        if (f > 0.1)  f = 0.1;

        if (fabs(df) < 1e-8) break;
    }

    return f;
}

/* ─── L2: Nusselt Number Correlations ───────────────────────────── */

double laminar_nu_pipe(double thermal_bc_type)
{
    /* thermal_bc_type: 0 = constant T_w, 1 = constant q_w'' */
    return (thermal_bc_type == 1) ? 4.36 : 3.66;
}

double laminar_sh_pipe(void)
{
    return 3.66;
}

/**
 * Hausen correlation for laminar developing flow (thermal entrance).
 *
 * Nu_D = 3.66 + 0.0668·Gz / (1 + 0.04·Gz^(2/3))
 *
 * where Gz = (D/L)·Re·Pr is the Graetz number.
 *
 * Valid: laminar flow, constant wall temperature, thermal entrance region.
 *
 * Complexity: O(1)
 */
double hausen_laminar_nu(double Re_D, double Pr, double D, double L)
{
    if (L <= 0.0) return -1.0;
    double Gz = (D / L) * Re_D * Pr;  /* Graetz number */
    return 3.66 + 0.0668 * Gz / (1.0 + 0.04 * pow(Gz, 2.0 / 3.0));
}

/**
 * Sieder-Tate correlation for laminar flow with viscosity correction.
 *
 * Nu_D = 1.86 · (Re·Pr·D/L)^(1/3) · (μ_b/μ_w)^(0.14)
 *
 * Accounts for property variation across the boundary layer.
 * Valid: 0.48 < Pr < 16700, Re·Pr·D/L > 10.
 *
 * Complexity: O(1)
 */
double sieder_tate_laminar_nu(double Re_D, double Pr, double D, double L,
                               double mu_bulk, double mu_wall)
{
    if (L <= 0.0 || mu_wall <= 0.0) return -1.0;
    double Gz = Re_D * Pr * D / L;
    if (Gz < 10.0) Gz = 10.0;  /* correlation validity limit */
    return 1.86 * pow(Gz, 1.0 / 3.0) * pow(mu_bulk / mu_wall, 0.14);
}

/**
 * Dittus-Boelter correlation for turbulent pipe flow.
 *
 * Nu_D = 0.023 · Re_D^(4/5) · Pr^n
 *
 * n = 0.4 for heating (T_w > T_b), n = 0.3 for cooling (T_w < T_b).
 * This is perhaps the most famous heat transfer correlation.
 * Valid: Re > 10⁴, 0.6 < Pr < 160, L/D > 10.
 *
 * Complexity: O(1)
 */
double dittus_boelter_nu(double Re_D, double Pr, int heating)
{
    double n = heating ? 0.4 : 0.3;
    return 0.023 * pow(Re_D, 0.8) * pow(Pr, n);
}

/**
 * Gnielinski correlation — more accurate than Dittus-Boelter.
 *
 * Nu = (f/8)·(Re-1000)·Pr / [1 + 12.7·√(f/8)·(Pr^(2/3)-1)]
 *
 * Valid: 3000 < Re < 5×10⁶, 0.5 < Pr < 2000.
 * Accounts for the friction factor explicitly, making the analogy
 * relationship clearer.
 *
 * Complexity: O(1)
 */
double gnielinski_nu(double Re_D, double Pr, double f)
{
    if (Re_D <= 1000.0) return laminar_nu_pipe(0);
    double sqrt_f8 = sqrt(f / 8.0);
    double numerator = (f / 8.0) * (Re_D - 1000.0) * Pr;
    double denominator = 1.0 + 12.7 * sqrt_f8 * (pow(Pr, 2.0 / 3.0) - 1.0);
    return numerator / denominator;
}

/**
 * Petukhov correlation — high accuracy for wide Pr range.
 *
 * Nu = (f/8)·Re·Pr / [1.07 + 12.7·√(f/8)·(Pr^(2/3)-1)]
 *
 * Valid: 10⁴ < Re < 5×10⁶, 0.5 < Pr < 200.
 * This correlation is recommended over Dittus-Boelter when high
 * accuracy is needed.
 *
 * Complexity: O(1)
 */
double petukhov_nu(double Re_D, double Pr, double f)
{
    double sqrt_f8 = sqrt(f / 8.0);
    double numerator = (f / 8.0) * Re_D * Pr;
    double denominator = 1.07 + 12.7 * sqrt_f8 * (pow(Pr, 2.0 / 3.0) - 1.0);
    return numerator / denominator;
}

/* ─── L5: Analogy-Based Predictions ─────────────────────────────── */

double colburn_j_factor_pipe(double f_F)
{
    return f_F;  /* f_F = j_H = j_D per Chilton-Colburn analogy */
}

double predict_nu_from_friction(double f_D, double Re, double Pr)
{
    double f_F = f_D / 4.0;
    return f_F * Re * pow(Pr, 1.0 / 3.0);
}

double predict_sh_from_friction(double f_D, double Re, double Sc)
{
    double f_F = f_D / 4.0;
    return f_F * Re * pow(Sc, 1.0 / 3.0);
}

/**
 * Predict h from pressure drop measurement.
 *
 * This is the most practical form of the analogy:
 *   Step 1: Measure ΔP over length L (manometer/pressure transducer)
 *   Step 2: Compute f = 2·ΔP·D/(ρ·v²·L)
 *   Step 3: Compute St from analogy: St = (f/2)·Pr^(-2/3)
 *   Step 4: Compute h = St·ρ·v·Cp
 *
 * No temperature measurement needed! This is the power of the analogy.
 *
 * Complexity: O(1)
 */
double predict_h_from_pressure_drop(double delta_P, double D, double L,
                                    double rho, double v, double k,
                                    double Pr)
{
    double f_D = 2.0 * delta_P * D / (rho * v * v * L);
    /* Colburn analogy: Nu = (f_D/8)·Re·Pr^(1/3)
     * We compute h directly using: h = St·ρ·v·cp
     * and cp = Pr·k/μ, so h = (f_D/8)·ρ·v·k·Pr^(1/3)/μ
     * Since μ is not provided, estimate from k and Pr:
     * μ_est = k·Pr / cp, with cp_est = 1000 for liquids, 4000 for water. */
    double cp_est = 4180.0;  /* water */
    double mu_est = k * Pr / cp_est;
    double Re_D = rho * v * D / mu_est;
    double Nu = (f_D / 8.0) * Re_D * pow(Pr, 1.0 / 3.0);
    return Nu * k / D;
}

double predict_kc_from_heat_transfer_coefficient(double h, double rho,
                                                  double cp, double Pr,
                                                  double Sc)
{
    /* j_H = (h/(ρ·v·cp))·Pr^(2/3) = j_D = (k_c/v)·Sc^(2/3) */
    /* Therefore: k_c = h/(ρ·cp) · (Pr/Sc)^(2/3) */
    return (h / (rho * cp)) * pow(Pr / Sc, 2.0 / 3.0);
}

/* ─── L5: Complete Pipe Flow Analogy ────────────────────────────── */

void compute_pipe_flow_analogy(double D, double L, double v_avg,
                               double rho, double mu, double cp,
                               double k, double D_AB,
                               double T_wall, double T_bulk,
                               double C_wall, double C_bulk,
                               PipeFlowAnalogy *pipe)
{
    if (!pipe) return;

    double A_c = 3.141592653589793 * D * D / 4.0;
    double A_s = 3.141592653589793 * D * L;

    double Re_D = rho * v_avg * D / mu;
    double Pr = mu * cp / k;
    double Sc = mu / (rho * D_AB);

    pipe->D = D;
    pipe->L = L;
    pipe->A_c = A_c;
    pipe->A_s = A_s;
    pipe->v_avg = v_avg;
    pipe->m_dot = rho * v_avg * A_c;
    pipe->Re_D = Re_D;
    pipe->Pr = Pr;
    pipe->Sc = Sc;

    /* Flow regime and friction factor */
    if (Re_D < 2300.0) {
        pipe->flow_regime = 0;
        pipe->f_Darcy = darcy_friction_laminar(Re_D);
        pipe->f_Fanning = pipe->f_Darcy / 4.0;
        /* Laminar: use Sieder-Tate or Hausen */
        pipe->Nu_D = sieder_tate_laminar_nu(Re_D, Pr, D, L, mu, mu);
        pipe->Sh_D = laminar_nu_pipe(0);  /* analog: 3.66 for const C_w */
    } else {
        pipe->flow_regime = 2;  /* turbulent */
        pipe->f_Darcy = blasius_friction_turbulent(Re_D);
        pipe->f_Fanning = pipe->f_Darcy / 4.0;
        /* Use Gnielinski correlation (more accurate) */
        pipe->Nu_D = gnielinski_nu(Re_D, Pr, pipe->f_Darcy);
        /* Mass transfer from analogy */
        pipe->Sh_D = predict_sh_from_friction(pipe->f_Darcy, Re_D, Sc);
    }

    pipe->Cf = 2.0 * pipe->f_Fanning;

    /* Transfer coefficients */
    pipe->h  = pipe->Nu_D * k / D;
    pipe->kc = pipe->Sh_D * D_AB / D;

    /* Pressure drop */
    pipe->delta_P = pipe->f_Darcy * (L / D) * (0.5 * rho * v_avg * v_avg);

    /* Heat transfer rate */
    double delta_T = T_wall - T_bulk;
    pipe->Q = pipe->h * A_s * delta_T;

    /* Mass transfer rate */
    double delta_C = C_wall - C_bulk;
    pipe->m_A = pipe->kc * A_s * delta_C;

    /* Stanton and j-factors */
    pipe->St_heat = pipe->Nu_D / (Re_D * Pr);
    pipe->St_mass = pipe->Sh_D / (Re_D * Sc);
    pipe->j_H = pipe->St_heat * pow(Pr, 2.0 / 3.0);
    pipe->j_D = pipe->St_mass * pow(Sc, 2.0 / 3.0);

    pipe->thermal_bc = 0;  /* constant T_w assumed */
}

/* ─── L5: Analogy Validation ────────────────────────────────────── */

void validate_analogy(double Re_D, double Pr, double Sc,
                      double f_measured, double Nu_measured,
                      double Sh_measured,
                      AnalogyValidation *valid)
{
    if (!valid) return;

    valid->Re_D = Re_D;
    valid->Pr = Pr;
    valid->Sc = Sc;
    valid->f_measured = f_measured;
    valid->Nu_measured = Nu_measured;
    valid->Sh_measured = Sh_measured;

    /* Predictions from analogy */
    valid->Nu_analogy = predict_nu_from_friction(f_measured, Re_D, Pr);
    valid->Sh_analogy = predict_sh_from_friction(f_measured, Re_D, Sc);

    /* Reverse: predict f from Nu */
    double f_F_pred = (Nu_measured / (Re_D * Pr)) * pow(Pr, 2.0 / 3.0);
    valid->f_analogy = 4.0 * f_F_pred;

    /* Relative errors */
    valid->err_Nu = fabs(valid->Nu_analogy - Nu_measured) / Nu_measured;
    valid->err_Sh = fabs(valid->Sh_analogy - Sh_measured)
                    / (Sh_measured + 1e-30);
    valid->err_f  = fabs(valid->f_analogy - f_measured) / f_measured;
}

/* ─── L6: Non-Circular Ducts ────────────────────────────────────── */

/**
 * Sets up duct geometry for non-circular cross-sections.
 *
 * Laminar f·Re values (fully developed):
 *   Circular:               f·Re = 64
 *   Parallel plates:        f·Re = 96
 *   Square:                 f·Re = 56.91
 *   Equilateral triangle:   f·Re = 53.33
 *   Rectangular AR=2:      f·Re = 62.19
 *   Rectangular AR=4:      f·Re = 73.0
 *   Rectangular AR=8:      f·Re = 82.3
 *
 * Laminar Nu (constant T_w):
 *   Circular:               Nu = 3.66
 *   Parallel plates:        Nu = 7.54
 *   Square:                 Nu = 2.98
 *   Equilateral triangle:   Nu = 2.35
 */
void setup_duct_geometry(int geometry_type, double aspect_ratio,
                         double A_c, DuctGeometry *duct)
{
    if (!duct) return;

    duct->geometry = geometry_type;
    duct->aspect_ratio = aspect_ratio;
    duct->A_c = A_c;

    switch (geometry_type) {
    case 0: /* Circular */
        duct->D_h = sqrt(4.0 * A_c / 3.141592653589793);
        duct->P_w = 3.141592653589793 * duct->D_h;
        duct->f_Re = 64.0;
        duct->Nu_fd = 3.66;
        break;
    case 1: /* Rectangular */
        /* A_c = W·H, aspect_ratio = W/H */
        {
            double H = sqrt(A_c / aspect_ratio);
            double W = H * aspect_ratio;
            duct->D_h = 4.0 * A_c / (2.0 * (W + H));
            duct->P_w = 2.0 * (W + H);

            /* f·Re interpolation (Shah & London, 1978) */
            double ar = aspect_ratio;
            duct->f_Re = 96.0 * (1.0 - 1.3553/ar + 1.9467/(ar*ar)
                                 - 1.7012/(ar*ar*ar) + 0.9564/(ar*ar*ar*ar)
                                 - 0.2537/(ar*ar*ar*ar*ar));
            duct->Nu_fd = 7.541 * (1.0 - 2.610/ar + 4.970/(ar*ar)
                                   - 5.119/(ar*ar*ar) + 2.702/(ar*ar*ar*ar)
                                   - 0.548/(ar*ar*ar*ar*ar));
        }
        break;
    case 2: /* Equilateral triangle */
        duct->D_h = sqrt(4.0 * A_c / 3.141592653589793);
        duct->P_w = 3.0 * sqrt(4.0 * A_c / sqrt(3.0));
        duct->f_Re = 53.33;
        duct->Nu_fd = 2.35;
        break;
    case 3: /* Annular */
        duct->f_Re = 64.0;    /* approximate */
        duct->Nu_fd = 3.66;   /* approximate */
        duct->P_w = 2.0 * 3.141592653589793 * sqrt(4.0*A_c/3.141592653589793);
        break;
    case 4: /* Parallel plates */
        duct->D_h = 2.0 * sqrt(A_c / aspect_ratio);
        duct->P_w = 2.0 * aspect_ratio;  /* per unit depth */
        duct->f_Re = 96.0;
        duct->Nu_fd = 7.54;
        break;
    default:
        duct->D_h = 1.0;
        duct->P_w = 1.0;
        duct->f_Re = 64.0;
        duct->Nu_fd = 3.66;
        break;
    }
}

void compute_developing_flow(double x, double D, double Re_D, double Pr,
                             DevelopingFlowState *dev)
{
    if (!dev) return;

    dev->x = x;
    dev->D = D;
    dev->Re_D = Re_D;
    dev->Pr = Pr;

    /* Hydrodynamic entry length: L_h ≈ 0.05·D·Re (laminar) */
    dev->L_hydro = 0.05 * D * Re_D;
    dev->x_plus_hydro = x / (D * Re_D);
    dev->hydro_developed = (x >= dev->L_hydro) ? 1 : 0;

    /* Thermal entry length: L_th ≈ 0.05·D·Re·Pr (laminar) */
    dev->L_thermal = 0.05 * D * Re_D * Pr;
    dev->x_plus_thermal = x / (D * Re_D * Pr);
    dev->thermal_developed = (x >= dev->L_thermal) ? 1 : 0;
}

double developing_nusselt(double Re_D, double Pr, double D, double x)
{
    /* Hausen correlation for thermal entrance region */
    return hausen_laminar_nu(Re_D, Pr, D, x);
}

/* ─── Diagnostics ───────────────────────────────────────────────── */

void print_pipe_flow_report(const PipeFlowAnalogy *pipe)
{
    if (!pipe) return;

    printf("━━━━ Pipe Flow Analogy Report ━━━━━━━━━━━━━━━━━━━━━━━━━\n");
    printf("  Geometry: D=%.4f m, L=%.3f m, A_c=%.6f m²\n",
           pipe->D, pipe->L, pipe->A_c);
    printf("  Flow: v=%.3f m/s, m_dot=%.4f kg/s, Re_D=%.1f (%s)\n",
           pipe->v_avg, pipe->m_dot, pipe->Re_D,
           pipe->flow_regime == 0 ? "laminar" : "turbulent");
    printf("  Fluid: Pr=%.3f, Sc=%.3f\n", pipe->Pr, pipe->Sc);
    printf("  ─── Momentum (Friction) ──────────────────────────\n");
    printf("  f_Darcy = %.5f    f_Fanning = %.5f\n",
           pipe->f_Darcy, pipe->f_Fanning);
    printf("  ΔP = %.2f Pa (%.4f bar)\n",
           pipe->delta_P, pipe->delta_P / 1e5);
    printf("  ─── Heat Transfer ───────────────────────────────\n");
    printf("  Nu_D = %.1f    h = %.1f W/(m²·K)\n", pipe->Nu_D, pipe->h);
    printf("  St_heat = %.6f    j_H = %.6f\n", pipe->St_heat, pipe->j_H);
    printf("  Q = %.1f W\n", pipe->Q);
    printf("  ─── Mass Transfer ───────────────────────────────\n");
    printf("  Sh_D = %.1f    k_c = %.3e m/s\n", pipe->Sh_D, pipe->kc);
    printf("  St_mass = %.6f    j_D = %.6f\n", pipe->St_mass, pipe->j_D);
    printf("  m_A = %.3e mol/s\n", pipe->m_A);
    printf("  ─── Analogy Check ───────────────────────────────\n");
    printf("  j_H / (f_F) = %.3f  (ideal=1)\n",
           pipe->j_H / (pipe->f_Fanning + 1e-30));
    printf("  j_D / (f_F) = %.3f  (ideal=1)\n",
           pipe->j_D / (pipe->f_Fanning + 1e-30));
    printf("  Colburn analogy holds: %s\n",
           fabs(pipe->j_H - pipe->f_Fanning) < 0.2 * pipe->f_Fanning
           ? "YES ✓" : "marginal");
    printf("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");
}

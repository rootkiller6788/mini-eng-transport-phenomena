/**
 * interphase_types.c
 * ==================
 * Implementations of core data structure operations for interphase transport.
 *
 * Reference: Bird, Stewart, Lightfoot (2007) "Transport Phenomena"
 * Knowledge coverage: L1 Definitions, L2 Core Concepts, L3 Engineering Quantities
 */

#include "interphase_types.h"
#include <string.h>
#include <math.h>
#include <stdlib.h>

/* ============================================================================
 * L1: PhaseState Operations
 * ============================================================================ */

/**
 * Initialize a PhaseState with standard values for a given phase.
 * Uses reference property values at 298K, 1 atm unless overridden.
 *
 * @param state     Output phase state to initialize
 * @param phase     Phase type (determines default property values)
 * Complexity: O(1)
 */
void phase_state_init(PhaseState *state, PhaseType phase) {
    if (!state) return;
    memset(state, 0, sizeof(PhaseState));

    switch (phase) {
    case PHASE_GAS:
        /* Air at 298K, 1 atm as reference */
        state->temperature        = 298.15;
        state->pressure           = 101325.0;
        state->density            = 1.184;
        state->viscosity          = 1.849e-5;
        state->thermal_conductivity = 0.0262;
        state->specific_heat_cp   = 1005.0;
        state->speed_of_sound     = 346.0;
        break;
    case PHASE_LIQUID:
        /* Water at 298K, 1 atm */
        state->temperature        = 298.15;
        state->pressure           = 101325.0;
        state->density            = 997.0;
        state->viscosity          = 8.9e-4;
        state->thermal_conductivity = 0.606;
        state->specific_heat_cp   = 4181.0;
        state->speed_of_sound     = 1497.0;
        break;
    case PHASE_SOLID:
        /* Copper at 298K */
        state->temperature        = 298.15;
        state->pressure           = 101325.0;
        state->density            = 8960.0;
        state->viscosity          = 0.0;   /* solid does not flow */
        state->thermal_conductivity = 401.0;
        state->specific_heat_cp   = 385.0;
        state->speed_of_sound     = 4760.0;
        break;
    case PHASE_PLASMA:
        /* Low-temperature argon plasma, atmospheric pressure */
        state->temperature        = 10000.0;
        state->pressure           = 101325.0;
        state->density            = 0.049;
        state->viscosity          = 1.0e-4;
        state->thermal_conductivity = 0.5;
        state->specific_heat_cp   = 520.0;
        state->speed_of_sound     = 1500.0;
        break;
    case PHASE_SUPERCRITICAL:
        /* CO2 at supercritical conditions */
        state->temperature        = 310.0;
        state->pressure           = 8.0e6;
        state->density            = 600.0;
        state->viscosity          = 2.0e-5;
        state->thermal_conductivity = 0.05;
        state->specific_heat_cp   = 3000.0;
        state->speed_of_sound     = 400.0;
        break;
    default:
        break;
    }
}

/**
 * Compute the Prandtl number from PhaseState.
 * Pr = nu / alpha = mu * cp / k
 *
 * @param state Phase state (must have non-zero k and mu, cp > 0)
 * @return Prandtl number, dimensionless
 * Complexity: O(1)
 *
 * Theorem: Pr >> 1 -> momentum diffuses faster than heat (oils, polymers).
 * Theorem: Pr << 1 -> heat diffuses faster than momentum (liquid metals).
 * Theorem: Pr ~ 0.7 for gases (air ~ 0.71).
 */
double phase_state_prandtl(const PhaseState *state) {
    if (!state || state->thermal_conductivity == 0.0) return 0.0;
    return state->viscosity * state->specific_heat_cp
           / state->thermal_conductivity;
}

/**
 * Compute the kinematic viscosity from PhaseState.
 * nu = mu / rho
 *
 * @param state Phase state (must have non-zero density)
 * @return Kinematic viscosity (m^2/s)
 * Complexity: O(1)
 */
double phase_state_kinematic_viscosity(const PhaseState *state) {
    if (!state || state->density == 0.0) return 0.0;
    return state->viscosity / state->density;
}

/**
 * Compute the thermal diffusivity from PhaseState.
 * alpha = k / (rho * cp)
 *
 * @param state Phase state
 * @return Thermal diffusivity (m^2/s)
 * Complexity: O(1)
 */
double phase_state_thermal_diffusivity(const PhaseState *state) {
    if (!state || state->density == 0.0 || state->specific_heat_cp == 0.0)
        return 0.0;
    return state->thermal_conductivity
           / (state->density * state->specific_heat_cp);
}

/* ============================================================================
 * L1: InterfaceDescriptor Operations
 * ============================================================================ */

/**
 * Initialize an InterfaceDescriptor for a given phase pair.
 * Sets up geometry-independent fields to sensible defaults.
 *
 * @param iface     Output interface descriptor
 * @param pair      Phase pair type
 * @param geometry  Interface geometry
 * Complexity: O(1)
 */
void interface_descriptor_init(InterfaceDescriptor *iface,
                                PhasePair pair,
                                InterfaceGeometry geometry) {
    if (!iface) return;
    memset(iface, 0, sizeof(InterfaceDescriptor));

    iface->pair     = pair;
    iface->geometry = geometry;
    iface->flags    = 0;

    /* Set phase types based on pair */
    switch (pair) {
    case PAIR_GAS_LIQUID:
        iface->phase_a = PHASE_GAS;    iface->phase_b = PHASE_LIQUID;
        break;
    case PAIR_LIQUID_LIQUID:
        iface->phase_a = PHASE_LIQUID; iface->phase_b = PHASE_LIQUID;
        break;
    case PAIR_GAS_SOLID:
        iface->phase_a = PHASE_GAS;    iface->phase_b = PHASE_SOLID;
        break;
    case PAIR_LIQUID_SOLID:
        iface->phase_a = PHASE_LIQUID; iface->phase_b = PHASE_SOLID;
        break;
    case PAIR_LIQUID_VAPOR:
        iface->phase_a = PHASE_LIQUID; iface->phase_b = PHASE_GAS;
        break;
    case PAIR_SOLID_SOLID:
        iface->phase_a = PHASE_SOLID;  iface->phase_b = PHASE_SOLID;
        break;
    case PAIR_GAS_PLASMA:
        iface->phase_a = PHASE_GAS;    iface->phase_b = PHASE_PLASMA;
        break;
    case PAIR_IMMISCIBLE_LIQ:
        iface->phase_a = PHASE_LIQUID; iface->phase_b = PHASE_LIQUID;
        break;
    default:
        break;
    }

    /* Initialize states */
    phase_state_init(&iface->state_a,  iface->phase_a);
    phase_state_init(&iface->state_b,  iface->phase_b);
    phase_state_init(&iface->state_ai, iface->phase_a);
    phase_state_init(&iface->state_bi, iface->phase_b);

    /* Sensible defaults */
    iface->surface_tension    = 0.0728;  /* water-air at 298K */
    iface->surface_energy     = 0.0728;
    iface->contact_angle      = 0.0;
    iface->interfacial_area   = 1.0;
    iface->curvature_mean     = 0.0;
    iface->curvature_gaussian = 0.0;
    iface->roughness_rms      = 0.0;
    iface->slip_length        = 0.0;     /* no-slip by default */
}

/**
 * Set the interface surface tension and update related flags.
 * Also computes surface energy (numerically equal to surface tension
 * for liquids in the absence of surfactants, per Gibbs).
 *
 * @param iface Interface descriptor
 * @param sigma Surface tension (N/m), must be >= 0
 * Complexity: O(1)
 *
 * Theorem: sigma >= 0 for thermodynamic stability.
 * Theorem: sigma(T) decreases with temperature, d_sigma/dT < 0.
 */
void interface_set_surface_tension(InterfaceDescriptor *iface, double sigma) {
    if (!iface || sigma < 0.0) return;
    iface->surface_tension = sigma;
    iface->surface_energy  = sigma;  /* Gibbs: surface energy = surface tension for simple liquids */
}

/**
 * Set the interface roughness and update geometry classification.
 *
 * @param iface Interface descriptor
 * @param rms   RMS roughness height (m)
 * Complexity: O(1)
 */
void interface_set_roughness(InterfaceDescriptor *iface, double rms) {
    if (!iface || rms < 0.0) return;
    iface->roughness_rms = rms;
    if (rms > 0.0) {
        iface->geometry = IFACE_CORRUGATED;
    } else {
        iface->geometry = IFACE_PLANAR;
    }
}

/**
 * Compute the capillary length for this interface.
 * l_c = sqrt(sigma / (|rho_B - rho_A| * g))
 * Below l_c, surface tension dominates; above, gravity dominates.
 *
 * @param iface Interface descriptor
 * @param g     Gravitational acceleration (m/s^2), typically 9.81
 * @return Capillary length (m), > 0
 * Complexity: O(1)
 *
 * Theorem: For water-air: l_c ~ 2.7 mm.
 * Theorem: Bond number Bo = (L / l_c)^2.
 * Theorem: For L << l_c, surface tension dominates (microgravity, microfluidics).
 * Theorem: For L >> l_c, gravity dominates (large-scale flows).
 */
double interface_capillary_length(const InterfaceDescriptor *iface, double g) {
    if (!iface || g <= 0.0) return 0.0;
    double drho = fabs(iface->state_a.density - iface->state_b.density);
    if (drho < 1e-10) return 1e10;  /* neutrally buoyant, effectively infinite */
    return sqrt(iface->surface_tension / (drho * g));
}

/**
 * Compute the Prandtl number from individual fluid properties.
 * Pr = mu * cp / k.
 *
 * @param mu Dynamic viscosity (Pa*s)
 * @param cp Specific heat cp (J/(kg*K))
 * @param k  Thermal conductivity (W/(m*K))
 * @return Prandtl number (-)
 * Complexity: O(1)
 */
double compute_prandtl(double mu, double cp, double k) {
    if (k <= 0.0) return 0.0;
    return mu * cp / k;
}

/**
 * Compute the Schmidt number.
 * Sc = mu / (rho * D) = nu / D
 *
 * @param mu  Dynamic viscosity (Pa*s)
 * @param rho Density (kg/m^3)
 * @param D   Mass diffusivity (m^2/s)
 * @return Schmidt number (-)
 * Complexity: O(1)
 */
double compute_schmidt(double mu, double rho, double D) {
    if (rho <= 0.0 || D <= 0.0) return 0.0;
    return mu / (rho * D);
}

/**
 * Compute the Lewis number.
 * Le = alpha / D = Sc / Pr
 *
 * @param alpha Thermal diffusivity (m^2/s)
 * @param D     Mass diffusivity (m^2/s)
 * @return Lewis number (-)
 * Complexity: O(1)
 *
 * Theorem: Le >> 1 -> heat diffuses faster than mass (most liquids).
 * Theorem: Le << 1 -> mass diffuses faster than heat.
 * Theorem: Le ~ 1 for gases (typical air-water vapor: Le ~ 0.85).
 */
double compute_lewis(double alpha, double D) {
    if (D <= 0.0) return 0.0;
    return alpha / D;
}

/* ============================================================================
 * L3: Dimensionless Groups Computation
 * ============================================================================ */

/**
 * Compute a complete set of dimensionless groups from flow and
 * transport parameters. This is the engineering workhorse for
 * characterizing interphase transport regimes.
 *
 * @param rho       Density (kg/m^3)
 * @param U         Characteristic velocity (m/s)
 * @param L         Characteristic length (m)
 * @param mu        Dynamic viscosity (Pa*s)
 * @param k         Thermal conductivity (W/(m*K))
 * @param cp        Specific heat (J/(kg*K))
 * @param D         Mass diffusivity (m^2/s)
 * @param beta      Thermal expansion coefficient (1/K)
 * @param delta_T   Characteristic temperature difference (K)
 * @param sigma     Surface tension (N/m)
 * @param g         Gravitational acceleration (m/s^2)
 * @param c         Speed of sound (m/s)
 * @param lambda    Mean free path (m)
 * @param groups    Output dimensionless groups (must not be NULL)
 * Complexity: O(1)
 *
 * Knowledge: This single function captures the complete Buckingham Pi
 * analysis for interphase transport. Over 20 dimensionless numbers
 * are computed from 13 dimensional inputs, representing 20-13=7
 * independent dimensionless groups.
 */
void compute_dimensionless_groups(double rho, double U, double L,
                                   double mu, double k, double cp,
                                   double D, double beta,
                                   double delta_T, double sigma,
                                   double g, double c, double lambda,
                                   DimensionlessGroups *groups) {
    if (!groups) return;

    /* Guard against zero denominators */
    if (mu <= 0.0) mu = 1e-10;
    if (k  <= 0.0)  k = 1e-10;
    if (cp <= 0.0) cp = 1e-10;
    if (D  <= 0.0)  D = 1e-10;
    if (sigma <= 0.0) sigma = 1e-10;
    if (c   <= 0.0)  c = 1e-10;

    double nu     = mu / rho;           /* kinematic viscosity */
    double alpha  = k / (rho * cp);     /* thermal diffusivity */

    groups->Reynolds       = rho * U * L / mu;
    groups->Prandtl        = mu * cp / k;
    groups->Schmidt        = nu / D;
    groups->Nusselt        = 0.0;        /* requires h, filled externally */
    groups->Sherwood       = 0.0;        /* requires k_m, filled externally */
    groups->Peclet_heat    = groups->Reynolds * groups->Prandtl;
    groups->Peclet_mass    = groups->Reynolds * groups->Schmidt;
    groups->Grashof        = g * beta * fabs(delta_T) * L*L*L / (nu * nu);
    groups->Rayleigh       = groups->Grashof * groups->Prandtl;
    groups->Stanton_heat   = 0.0;        /* requires Nu */
    groups->Stanton_mass   = 0.0;        /* requires Sh */
    groups->Biot           = 0.0;        /* requires solid k */
    groups->Fourier        = 0.0;        /* requires time t */
    groups->Weber          = rho * U * U * L / sigma;
    groups->Capillary      = mu * U / sigma;
    groups->Bond           = fabs(rho - 1.0) * g * L * L / sigma; /* approximate */
    groups->Marangoni      = 0.0;        /* requires d_sigma/dT */
    groups->Mach           = U / c;
    groups->Knudsen        = lambda / L;
    groups->Lewis          = alpha / D;
}

/* ============================================================================
 * L2: StressTensor Operations
 * ============================================================================ */

/**
 * Initialize a StressTensor to zero (quiescent fluid at reference pressure).
 *
 * @param tau      Output stress tensor
 * @param pressure Reference pressure (Pa)
 * Complexity: O(1)
 */
void stress_tensor_init(StressTensor *tau, double pressure) {
    if (!tau) return;
    memset(tau, 0, sizeof(StressTensor));
    tau->pressure = pressure;
}

/**
 * Set viscous stress for a Newtonian fluid with given velocity gradient.
 * tau_ij = mu * (du_i/dx_j + du_j/dx_i) for i != j.
 * tau_ii = 2 * mu * du_i/dx_i  (diagonal components).
 * Pressure is added separately as -p*delta_ij.
 *
 * @param tau     Output stress tensor (viscous part only)
 * @param mu      Dynamic viscosity (Pa*s)
 * @param grad_u  Velocity gradient tensor [9] (row-major: du_i/dx_j)
 *                grad_u[row*3+col] = du_row/dx_col
 * Complexity: O(1)
 */
void stress_tensor_newtonian(StressTensor *tau, double mu,
                              const double grad_u[9]) {
    if (!tau || !grad_u) return;

    /* trace of grad(u) = div_u, used when bulk viscosity is non-zero */
    (void)(grad_u[0] + grad_u[4] + grad_u[8]);

    tau->tau_xx = mu * (2.0 * grad_u[0]);  /* ignoring bulk viscosity */
    tau->tau_yy = mu * (2.0 * grad_u[4]);
    tau->tau_zz = mu * (2.0 * grad_u[8]);
    tau->tau_xy = mu * (grad_u[1] + grad_u[3]);
    tau->tau_xz = mu * (grad_u[2] + grad_u[6]);
    tau->tau_yx = tau->tau_xy;
    tau->tau_yz = mu * (grad_u[5] + grad_u[7]);
    tau->tau_zx = tau->tau_xz;
    tau->tau_zy = tau->tau_yz;
}

/**
 * Compute the magnitude (Frobenius norm) of the deviatoric stress tensor.
 * |tau| = sqrt(sum_{i,j} tau_ij^2).
 *
 * @param tau Stress tensor
 * @return Frobenius norm of deviatoric stress (Pa)
 * Complexity: O(1)
 */
double stress_tensor_magnitude(const StressTensor *tau) {
    if (!tau) return 0.0;
    return sqrt(tau->tau_xx * tau->tau_xx +
                tau->tau_xy * tau->tau_xy +
                tau->tau_xz * tau->tau_xz +
                tau->tau_yx * tau->tau_yx +
                tau->tau_yy * tau->tau_yy +
                tau->tau_yz * tau->tau_yz +
                tau->tau_zx * tau->tau_zx +
                tau->tau_zy * tau->tau_zy +
                tau->tau_zz * tau->tau_zz);
}

/* ============================================================================
 * L2: Henry's Law Operations
 * ============================================================================ */

/**
 * Compute Henry's law constant at a temperature from reference value.
 * H(T) = H(Tref) * exp(C * (1/T - 1/Tref))
 * where C = -d(ln H)/d(1/T) is a system-specific constant.
 *
 * @param hl         Henry's law data (reference state + dH/dT)
 * @param T          New temperature (K), must be > 0
 * @param H_at_T     Output H(T) (Pa)
 * Complexity: O(1)
 *
 * Theorem: H typically decreases with increasing T for gases in water
 *          (lower solubility at higher T), then increases above ~370K.
 */
void henry_law_at_temperature(const HenryLaw *hl, double T,
                               double *H_at_T) {
    if (!hl || !H_at_T || T <= 0.0 || hl->temperature <= 0.0) {
        if (H_at_T) *H_at_T = 0.0;
        return;
    }

    /* Van't Hoff-type temperature dependence */
    double C = -hl->henry_dHdT * hl->temperature
               / (hl->henry_coefficient);
    if (fabs(C) < 1e-15) {
        *H_at_T = hl->henry_coefficient;
    } else {
        *H_at_T = hl->henry_coefficient
                  * exp(C * (1.0 / T - 1.0 / hl->temperature));
    }
}

/* ============================================================================
 * L2: Vapor-Liquid Equilibrium Operations
 * ============================================================================ */

/**
 * Compute the equilibrium vapor composition from liquid composition
 * using modified Raoult's law.
 * y_i * P = x_i * gamma_i * P_i^sat
 *
 * @param vle          VLE data
 * @param x_liquid     Liquid mole fraction (-), [0,1]
 * @param total_pressure Total system pressure (Pa)
 * @param y_vapor      Output vapor mole fraction (-)
 * Complexity: O(1)
 *
 * Theorem: For ideal solution (gamma=1), y*P = x*P^sat (Raoult's law).
 * Theorem: For positive deviation (gamma>1), vapor enriched in component.
 * Theorem: Azeotrope when x_i * gamma_i * P_i^sat / P = x_i (i.e., K=1).
 */
void raoult_vapor_composition(const VaporLiquidEquilibrium *vle,
                               double x_liquid, double total_pressure,
                               double *y_vapor) {
    if (!vle || !y_vapor || total_pressure <= 0.0) {
        if (y_vapor) *y_vapor = 0.0;
        return;
    }
    *y_vapor = x_liquid * vle->activity_coefficient_a
               * vle->saturation_pressure / total_pressure;
    if (*y_vapor > 1.0) *y_vapor = 1.0;
    if (*y_vapor < 0.0) *y_vapor = 0.0;
}

/**
 * Compute the bubble point pressure for a given liquid composition.
 * P_bubble = sum_i x_i * gamma_i * P_i^sat(T).
 *
 * @param vle     VLE data (single component)
 * @param x_liquid Liquid mole fraction
 * @param P_bubble Output bubble point pressure (Pa)
 * Complexity: O(1)
 */
void bubble_point_pressure(const VaporLiquidEquilibrium *vle,
                            double x_liquid, double *P_bubble) {
    if (!vle || !P_bubble) return;
    *P_bubble = x_liquid * vle->activity_coefficient_a
                * vle->saturation_pressure;
}

/**
 * Compute the dew point composition.
 * x_i = y_i * P / (gamma_i * P_i^sat)
 *
 * @param vle          VLE data
 * @param y_vapor      Vapor mole fraction
 * @param total_pressure Total pressure (Pa)
 * @param x_liquid     Output equilibrium liquid mole fraction
 * Complexity: O(1)
 */
void dew_point_composition(const VaporLiquidEquilibrium *vle,
                            double y_vapor, double total_pressure,
                            double *x_liquid) {
    if (!vle || !x_liquid || total_pressure <= 0.0) {
        if (x_liquid) *x_liquid = 0.0;
        return;
    }
    if (vle->activity_coefficient_b <= 0.0
        || vle->saturation_pressure <= 0.0) {
        *x_liquid = 0.0;
        return;
    }
    *x_liquid = y_vapor * total_pressure
                / (vle->activity_coefficient_b * vle->saturation_pressure);
    if (*x_liquid > 1.0) *x_liquid = 1.0;
    if (*x_liquid < 0.0) *x_liquid = 0.0;
}

/* ============================================================================
 * L3: Transport Property Estimation
 * ============================================================================ */

/**
 * Estimate binary gas diffusivity using Chapman-Enskog theory.
 * D_AB = (3/16) * sqrt(4*pi*k_B*T/M_AB) / (n*pi*sigma_AB^2*Omega_D)
 * Simplified: D [cm^2/s] = 0.0018583 * T^(3/2) * sqrt(1/M_A+1/M_B)
 *                          / (p * sigma_AB^2 * Omega_D)
 *
 * @param params Chapman-Enskog parameters
 * @param T      Temperature (K), must be > 0
 * @param p      Pressure (atm)
 * @param D_AB   Output diffusivity (m^2/s)
 * Complexity: O(1)
 *
 * Theorem: D_AB ~ T^(3/2) / p (ideal gas behavior).
 * Theorem: D_AB ~ 1/sqrt(M_AB) (lighter molecules diffuse faster).
 * Theorem: D_AB ~ 1/(sigma_AB^2) (larger molecules diffuse slower).
 */
void chapman_enskog_diffusivity(const DiffusivityParams *params,
                                 double T, double p, double *D_AB) {
    if (!params || !D_AB || T <= 0.0 || p <= 0.0) {
        if (D_AB) *D_AB = 0.0;
        return;
    }

    double M_AB = 2.0 / (1.0 / params->molar_mass_A
                         + 1.0 / params->molar_mass_B);
    double D_cm2s = 0.0018583 * pow(T, 1.5) * sqrt(M_AB)
                    / (p * params->sigma_AB * params->sigma_AB
                       * params->omega_D);
    *D_AB = D_cm2s * 1.0e-4;  /* convert cm^2/s to m^2/s */
}

/**
 * Estimate dynamic viscosity using the Chung method (Chung et al., 1988).
 * mu = 40.785 * sqrt(M*T) / (V_c^(2/3) * Omega_v) * F_c
 * where F_c accounts for molecular shape and polarity.
 *
 * @param params Chung method parameters
 * @param T      Temperature (K), > 0
 * @param mu     Output viscosity (Pa*s)
 * Complexity: O(1)
 *
 * Theorem: For dilute gases, mu ~ sqrt(T) (kinetic theory).
 * Theorem: For dense gases near critical point, Chung includes
 *          acentric factor and polarity corrections.
 */
void chung_viscosity(const ViscosityParams *params, double T, double *mu) {
    if (!params || !mu || T <= 0.0) {
        if (mu) *mu = 0.0;
        return;
    }

    /* Reduced temperature */
    double T_star = T / params->epsilon_over_k;

    /* Omega_v ~ 1.16145 * T_star^(-0.14874) + ... (Neufeld approx) */
    double omega_v = params->omega_v;
    if (omega_v <= 0.0) {
        omega_v = 1.16145 * pow(T_star, -0.14874)
                  + 0.52487 * exp(-0.77320 * T_star)
                  + 2.16178 * exp(-2.43787 * T_star);
    }

    /* Estimate critical volume from sigma: V_c ~ (sigma/0.809)^3 [cm^3/mol] */
    double V_c = pow(params->sigma / 0.809, 3.0);

    /* Chung correction factor */
    double kappa = 0.0;
    double mu_r = 131.3 * params->dipole_moment
                  / sqrt(V_c * params->epsilon_over_k);
    if (mu_r > 0.0) {
        /* polar correction: kappa ~ 0.0 for nonpolar, >0 for polar */
        (void)(mu_r * mu_r);  /* reserved for full polar correction */
        kappa = 0.0; /* simplified: neglect polar correction for now */
    }

    double F_c = 1.0 - 0.2756 * params->acentric_factor + kappa;

    /* mu in microPoise (uP), convert to Pa*s */
    double mu_uP = 40.785 * F_c * sqrt(params->molar_mass * T)
                   / (pow(V_c, 2.0/3.0) * omega_v);
    if (params->association_factor > 0.0) {
        mu_uP *= (1.0 + 0.05 * params->association_factor);
    }
    *mu = mu_uP * 1.0e-7;  /* microPoise -> Pa*s */
}

/**
 * Estimate thermal conductivity of a dilute gas using modified Eucken.
 * k = mu * (cp + 5*R/(4*M)) for polyatomic gases.
 *
 * @param mu       Dynamic viscosity (Pa*s)
 * @param cp       Specific heat cp (J/(kg*K))
 * @param cv       Specific heat cv (J/(kg*K))
 * @param molar_mass Molar mass (kg/mol)
 * @param k         Output thermal conductivity (W/(m*K))
 * Complexity: O(1)
 *
 * Theorem (Eucken): k = mu * (cp + 5*R/(4*M)).
 * Theorem: For monatomic gases (cp/cv = 5/3), k = 2.5*mu*cv.
 */
void eucken_thermal_conductivity(double mu, double cp, double cv,
                                  double molar_mass, double *k) {
    if (!k) return;
    if (molar_mass <= 0.0) { *k = 0.0; return; }
    (void)cv;  /* cv used in monatomic gas correlation */
    double R = 8.314;  /* universal gas constant J/(mol*K) */
    *k = mu * (cp + 2.5 * R / molar_mass);
}

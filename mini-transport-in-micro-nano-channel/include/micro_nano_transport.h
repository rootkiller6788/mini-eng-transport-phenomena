#ifndef MICRO_NANO_TRANSPORT_H
#define MICRO_NANO_TRANSPORT_H
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>

/* Physical constants */
#define KB 1.380649e-23
#define E_CHARGE 1.602176634e-19
#define NA 6.02214076e23
#define R_GAS 8.314462618
#define EPSILON0 8.854187817e-12
#define MAX_ION_SPECIES 8
#define MAX_CHANNEL_SEGMENTS 32
#define MAX_GRID_POINTS 10000

/* L1: Knudsen number regime */
typedef enum {
    KNUDSEN_CONTINUUM = 0,
    KNUDSEN_SLIP_FLOW,
    KNUDSEN_TRANSITION,
    KNUDSEN_FREE_MOLECULAR
} KnudsenRegime;

typedef struct {
    double mean_free_path_m;
    double knudsen_number;
    double characteristic_len_m;
    KnudsenRegime regime;
} KnudsenState;

/* L1: Slip boundary models */
typedef enum {
    SLIP_MODEL_NONE = 0,
    SLIP_MODEL_MAXWELL_FIRST,
    SLIP_MODEL_MAXWELL_SECOND,
    SLIP_MODEL_LINEAR_NAVIER,
    SLIP_MODEL_POWER_LAW,
    SLIP_MODEL_MOLECULAR
} SlipModelType;

typedef struct {
    SlipModelType model;
    double slip_length_m;
    double tangential_momentum_acc;
    double wall_shear_rate;
    double slip_velocity;
} SlipBoundary;

/* L1: Temperature jump */
typedef struct {
    double jump_length_m;
    double thermal_accommodation;
    double wall_temperature;
    double gas_temperature_at_wall;
    double temperature_jump;
} TemperatureJump;

/* L1: Debye state */
typedef struct {
    double debye_length_m;
    double debye_parameter;
    double solvent_permittivity;
    double temperature;
    double ionic_strength;
    int num_ion_species;
    double concentrations[MAX_ION_SPECIES];
    int valences[MAX_ION_SPECIES];
} DebyeState;

/* L1: Zeta potential */
typedef struct {
    double zeta_potential_V;
    double surface_charge_density;
    double stern_layer_capacitance;
    double slipping_plane_distance_nm;
} ZetaPotential;

/* L2: Channel geometry */
typedef enum {
    GEOM_PARALLEL_PLATES,
    GEOM_CYLINDRICAL_CAPILLARY,
    GEOM_RECTANGULAR_CHANNEL,
    GEOM_TRIANGULAR_GROOVE,
    GEOM_TRAPEZOIDAL_CHANNEL,
    GEOM_NANOPORE
} ChannelGeometryType;

typedef struct {
    ChannelGeometryType type;
    double height_m;
    double width_m;
    double length_m;
    double radius_m;
    double hydraulic_diameter_m;
    double cross_sectional_area_m2;
    double wetted_perimeter_m;
    double aspect_ratio;
    double surface_roughness_nm;
} ChannelGeometry;

/* L2: Fluid properties */
typedef struct {
    double density;
    double dynamic_viscosity;
    double kinematic_viscosity;
    double thermal_conductivity;
    double specific_heat_cp;
    double thermal_diffusivity;
    double mass_diffusivity;
    double prandtl_number;
    double schmidt_number;
    double relative_permittivity;
    double temperature;
    double molecular_mass;
    double collision_diameter;
} FluidProperties;

/* L4: Velocity profile result */
typedef struct {
    double *y_coords;
    double *velocity;
    double *shear_rate;
    int n_points;
    double y_max;
    double slip_velocity_lo;
    double slip_velocity_hi;
    double max_velocity;
    double avg_velocity;
    double flow_rate;
    double friction_factor;
    double poiseuille_number;
} VelocityProfile;

/* L4: Ionic profile result */
typedef struct {
    double *y_coords;
    double *potential_V;
    double *charge_density;
    double **ion_concentrations;
    double *total_ion_flux;
    int n_species;
    int n_points;
    double debye_length;
    double zeta_potential;
    double surface_potential;
    double gouy_chapman_length;
} IonicProfile;

/* L5: Solver config */
typedef struct {
    double pressure_gradient;
    double applied_voltage;
    double electric_field;
    int n_grid;
    double tolerance;
    int max_iterations;
    int include_slip;
    int include_electroosmosis;
    int include_pressure_gradient;
    int debye_huckel_linear;
} SolverConfig;

/* L6: Full transport problem */
typedef struct {
    ChannelGeometry geometry;
    FluidProperties fluid;
    KnudsenState knudsen;
    SlipBoundary slip;
    TemperatureJump temp_jump;
    DebyeState debye;
    ZetaPotential zeta;
    SolverConfig config;
} TransportProblem;

/* API declarations */
void compute_mean_free_path(FluidProperties *fluid, double pressure, double *mfp);
void compute_knudsen_number(double mfp, double L_char, KnudsenState *state);
void compute_hydraulic_diameter(ChannelGeometry *geom);
void compute_debye_length(DebyeState *state, double *lambda_d);
void compute_kinetic_fluid_props(double T, double p, double m_mol, double d_col, FluidProperties *out);
void transport_problem_init(TransportProblem *problem);
void velocity_profile_free(VelocityProfile *profile);
void ionic_profile_free(IonicProfile *profile);

double compute_maxwell_slip_length(double sigma, double mean_free_path);
double compute_temperature_jump_length(double alpha_T, double mean_free_path, double gamma, double Pr);
double compute_helmholtz_smoluchowski_velocity(double eps, double zeta, double eta, double E);
double compute_ionic_current_nanochannel(int n_species, const double *conc, const int *valence,
    const double *D, double V_applied, double L, double A, double T);
void compute_poisson_boltzmann_1d(const double *y, int n, double psi0, double kappa, double *psi);
void compute_ion_concentration_profile(const double *psi, int n, double c_bulk, int z, double T, double *c);
void solve_pressure_driven_slip_flow(double H, double dpdx, double eta, double slip_len, int n, double *y, double *u);
void solve_electroosmotic_flow(double H, double eps, double zeta, double eta, double Ex, double kappa, int n, double *y, double *u);
void compute_streaming_potential(double eps, double zeta, double eta, double sigma_bulk, double dpdx, double *V_str);
double compute_slip_enhanced_flow_rate(double H, double W, double dpdx, double eta, double slip_len);
double compute_poiseuille_number_slip(double H, double slip_len);
double compute_rarefied_mass_flow_rate(double H, double dpdx, double eta, double rho, double Kn);
void compute_transition_regime_mass_flow(double H, double dpdx, double eta, double rho, double Kn, double *mdot);
double compute_enhancement_factor_slip(double H, double slip_len);
double compute_nanochannel_conductance(double c_bulk, double z, double D, double L, double A, double sigma_s, double H, double T);
double compute_surface_charge_potential(double sigma_s, double c_bulk, double z, double eps, double T);
void compute_concentration_polarization(double c_bulk, double c_wall, double x_max, double delta_CP, double *x, double *c, int n);
double compute_bjerrum_length(double eps, double T);
double compute_gouy_chapman_length(double sigma_s, double l_B, int z);
double compute_overscreening_parameter(double sigma_s, double c_bulk, double z, double eps, double T);
void compute_nernst_planck_steady_1d(const double *c, const double *phi, double D, double z, double T, double L, int n, double *flux);
void compute_diffusion_potential(double c_hi, double c_lo, double z, double D_plus, double D_minus, double T, double *V_diff);
double compute_aspect_ratio_effect(double aspect_ratio, double *Po);
double compute_surface_to_volume_ratio(const ChannelGeometry *geom);
double compute_hydraulic_resistance(double H, double W, double L, double eta);
double compute_hydraulic_inductance(double rho, double L, const ChannelGeometry *geom);
double compute_womersley_number(double H, double f, double rho, double eta);

#endif

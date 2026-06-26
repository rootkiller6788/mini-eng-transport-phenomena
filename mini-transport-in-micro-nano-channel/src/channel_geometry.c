#include "micro_nano_transport.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Channel Geometry and Hydraulic Computations */

void compute_hydraulic_diameter(ChannelGeometry *geom)
{
    if (geom == NULL) return;
    switch (geom->type) {
    case GEOM_PARALLEL_PLATES:
        geom->cross_sectional_area_m2 = geom->height_m * geom->width_m;
        geom->wetted_perimeter_m = 2.0 * geom->width_m;
        geom->hydraulic_diameter_m = 2.0 * geom->height_m;
        break;
    case GEOM_CYLINDRICAL_CAPILLARY:
        geom->cross_sectional_area_m2 = M_PI * geom->radius_m * geom->radius_m;
        geom->wetted_perimeter_m = 2.0 * M_PI * geom->radius_m;
        geom->hydraulic_diameter_m = 2.0 * geom->radius_m;
        break;
    case GEOM_RECTANGULAR_CHANNEL:
        geom->cross_sectional_area_m2 = geom->height_m * geom->width_m;
        geom->wetted_perimeter_m = 2.0 * (geom->height_m + geom->width_m);
        geom->hydraulic_diameter_m = 4.0 * geom->cross_sectional_area_m2
                                     / geom->wetted_perimeter_m;
        geom->aspect_ratio = geom->width_m / geom->height_m;
        break;
    case GEOM_TRIANGULAR_GROOVE:
        geom->cross_sectional_area_m2 = 0.5 * geom->height_m * geom->width_m;
        geom->wetted_perimeter_m = geom->width_m
            + 2.0 * sqrt(geom->height_m * geom->height_m
                         + 0.25 * geom->width_m * geom->width_m);
        geom->hydraulic_diameter_m = 4.0 * geom->cross_sectional_area_m2
                                     / geom->wetted_perimeter_m;
        break;
    case GEOM_TRAPEZOIDAL_CHANNEL:
        {
            double top_w = geom->width_m;
            double bot_w = 0.7 * geom->width_m;
            double side = sqrt(geom->height_m * geom->height_m
                + 0.25 * (top_w - bot_w) * (top_w - bot_w));
            geom->cross_sectional_area_m2 = 0.5 * (top_w + bot_w)
                                            * geom->height_m;
            geom->wetted_perimeter_m = bot_w + 2.0 * side;
            geom->hydraulic_diameter_m = 4.0 * geom->cross_sectional_area_m2
                                         / geom->wetted_perimeter_m;
        }
        break;
    case GEOM_NANOPORE:
        geom->cross_sectional_area_m2 = M_PI * geom->height_m
                                        * geom->height_m / 4.0;
        geom->wetted_perimeter_m = M_PI * geom->height_m;
        geom->hydraulic_diameter_m = geom->height_m;
        break;
    default:
        geom->hydraulic_diameter_m = geom->height_m;
        geom->cross_sectional_area_m2 = geom->height_m * geom->width_m;
        geom->wetted_perimeter_m = 2.0 * (geom->height_m + geom->width_m);
        break;
    }
}

double compute_aspect_ratio_effect(double aspect_ratio, double *Po)
{
    if (aspect_ratio < 1.0) aspect_ratio = 1.0 / aspect_ratio;
    if (aspect_ratio < 1.0) return -1.0;
    double inv_AR = 1.0 / aspect_ratio;
    double f_Re = 96.0 * (1.0 - 1.3553 * inv_AR
                          + 1.9467 * inv_AR * inv_AR
                          - 1.7012 * inv_AR * inv_AR * inv_AR
                          + 0.9564 * inv_AR * inv_AR * inv_AR * inv_AR
                          - 0.2537 * inv_AR * inv_AR * inv_AR * inv_AR * inv_AR);
    if (Po != NULL) *Po = f_Re;
    return f_Re;
}

double compute_surface_to_volume_ratio(const ChannelGeometry *geom)
{
    if (geom == NULL) return -1.0;
    if (geom->cross_sectional_area_m2 <= 0.0) return -1.0;
    return geom->wetted_perimeter_m / geom->cross_sectional_area_m2;
}

double compute_hydraulic_resistance(double H, double W, double L, double eta)
{
    if (H <= 0.0 || W <= 0.0 || L <= 0.0 || eta <= 0.0) return -1.0;
    return 12.0 * eta * L / (H * H * H * W);
}

double compute_hydraulic_inductance(double rho, double L,
    const ChannelGeometry *geom)
{
    if (geom == NULL || rho <= 0.0 || L <= 0.0) return -1.0;
    if (geom->cross_sectional_area_m2 <= 0.0) return -1.0;
    return rho * L / geom->cross_sectional_area_m2;
}

double compute_womersley_number(double H, double f, double rho, double eta)
{
    if (H <= 0.0 || f <= 0.0 || rho <= 0.0 || eta <= 0.0) return -1.0;
    return 0.5 * H * sqrt(2.0 * M_PI * f * rho / eta);
}

void transport_problem_init(TransportProblem *problem)
{
    if (problem == NULL) return;
    memset(problem, 0, sizeof(TransportProblem));
    problem->geometry.type = GEOM_PARALLEL_PLATES;
    problem->geometry.height_m = 1.0e-6;
    problem->geometry.width_m = 1.0e-4;
    problem->geometry.length_m = 1.0e-3;
    compute_hydraulic_diameter(&problem->geometry);
    problem->fluid.relative_permittivity = 80.0;
    problem->fluid.temperature = 293.15;
    problem->config.n_grid = 101;
    problem->config.tolerance = 1.0e-8;
    problem->config.max_iterations = 1000;
    problem->config.include_slip = 0;
    problem->config.include_electroosmosis = 0;
    problem->config.include_pressure_gradient = 1;
    problem->config.debye_huckel_linear = 1;
}

void velocity_profile_free(VelocityProfile *profile)
{
    if (profile == NULL) return;
    free(profile->y_coords);
    free(profile->velocity);
    free(profile->shear_rate);
    memset(profile, 0, sizeof(VelocityProfile));
}

void ionic_profile_free(IonicProfile *profile)
{
    if (profile == NULL) return;
    free(profile->y_coords);
    free(profile->potential_V);
    free(profile->charge_density);
    if (profile->ion_concentrations != NULL) {
        for (int i = 0; i < profile->n_species; i++) {
            free(profile->ion_concentrations[i]);
        }
        free(profile->ion_concentrations);
    }
    free(profile->total_ion_flux);
    memset(profile, 0, sizeof(IonicProfile));
}

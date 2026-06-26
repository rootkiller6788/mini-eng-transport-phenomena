/**
 * interphase_config.h
 * ====================
 * Common configuration and math constants for interphase transport module.
 *
 * Provides platform-independent definitions of mathematical constants
 * not guaranteed by the C11 standard.
 */

#ifndef INTERPHASE_CONFIG_H
#define INTERPHASE_CONFIG_H

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif

#ifndef M_E
#define M_E 2.71828182845904523536
#endif

/* Universal constants in SI units */
#define R_GAS_UNIVERSAL  8.314462618   /* J/(mol*K) */
#define G_GRAVITY        9.80665       /* m/s^2 (standard gravity) */
#define SIGMA_SB         5.670374419e-8 /* Stefan-Boltzmann W/(m^2*K^4) */
#define F_FARADAY        96485.33212   /* C/mol */
#define KB_BOLTZMANN     1.380649e-23  /* J/K */
#define NA_AVOGADRO      6.02214076e23 /* 1/mol */

/* Engineering tolerances */
#define TOL_TIGHT 1e-15
#define TOL_MEDIUM 1e-9
#define TOL_LOOSE 1e-6

#endif /* INTERPHASE_CONFIG_H */

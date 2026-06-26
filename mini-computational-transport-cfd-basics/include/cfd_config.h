/**
 * cfd_config.h
 * =============
 * Configuration, math, and physical constants for computational transport CFD.
 *
 * Reference: Ferziger, Peric & Street (2020) "Computational Methods for
 *            Fluid Dynamics" (4th ed.), Springer.
 *            LeVeque (2002) "Finite Volume Methods for Hyperbolic Problems"
 *
 * Knowledge coverage: L1 Definitions, L2 Numerical constants
 */

#ifndef CFD_CONFIG_H
#define CFD_CONFIG_H

#include <math.h>
#include <stddef.h>

/* Solver tolerances */
#ifndef CFD_EPS
#define CFD_EPS               1e-8
#endif

#ifndef CFD_EPS_MACHINE
#define CFD_EPS_MACHINE       1e-14
#endif

#ifndef CFD_MAX_ITER
#define CFD_MAX_ITER          100000
#endif

/* Under-relaxation factors */
#ifndef CFD_URF_U
#define CFD_URF_U             0.7
#endif
#ifndef CFD_URF_P
#define CFD_URF_P             0.3
#endif
#ifndef CFD_URF_SCALAR
#define CFD_URF_SCALAR        0.8
#endif

/* Stability limits */
#ifndef CFD_CFL_MAX
#define CFD_CFL_MAX           1.0
#endif

#ifndef CFD_PRESSURE_CORR_ITER
#define CFD_PRESSURE_CORR_ITER  10
#endif

#ifndef CFD_MG_MAX_LEVELS
#define CFD_MG_MAX_LEVELS     6
#endif

#ifndef CFD_SOR_OMEGA
#define CFD_SOR_OMEGA         1.5
#endif

/* ============================================================================
 * L2: Numerical Scheme Selectors
 * ============================================================================ */

/** Convective discretization scheme taxonomy */
typedef enum {
    CFD_SCHEME_CENTRAL  = 0,  /**< 2nd-order central differencing */
    CFD_SCHEME_UPWIND1  = 1,  /**< 1st-order upwind (UDS) */
    CFD_SCHEME_HYBRID   = 2,  /**< Hybrid scheme (Spalding 1972) */
    CFD_SCHEME_POWERLAW = 3,  /**< Power-law scheme (Patankar 1980) */
    CFD_SCHEME_QUICK    = 4,  /**< QUICK (Leonard 1979) 3rd-order */
    CFD_SCHEME_TVD      = 5,  /**< TVD scheme (Sweby 1984) */
    CFD_SCHEME_WENO5    = 6   /**< WENO5 (Jiang & Shu 1996) */
} CfdConvectiveScheme;

/** Temporal discretization scheme taxonomy */
typedef enum {
    CFD_TIME_EXPLICIT_EULER = 0,  /**< 1st-order forward Euler */
    CFD_TIME_IMPLICIT_EULER = 1,  /**< 1st-order backward Euler */
    CFD_TIME_CRANK_NICOLSON = 2,  /**< 2nd-order Crank-Nicolson */
    CFD_TIME_ADAMS_BASHFORTH = 3, /**< 2nd-order Adams-Bashforth */
    CFD_TIME_RK4            = 4   /**< 4th-order Runge-Kutta */
} CfdTimeScheme;

/** Linear solver type taxonomy */
typedef enum {
    CFD_SOLVER_JACOBI       = 0,  /**< Jacobi iteration */
    CFD_SOLVER_GAUSS_SEIDEL = 1,  /**< Gauss-Seidel iteration */
    CFD_SOLVER_SOR          = 2,  /**< Successive Over-Relaxation */
    CFD_SOLVER_TDMA         = 3,  /**< Thomas TDMA (tri-diagonal) */
    CFD_SOLVER_CG           = 4,  /**< Conjugate Gradient (SPD) */
    CFD_SOLVER_BICGSTAB     = 5,  /**< BiCGSTAB (non-symmetric) */
    CFD_SOLVER_MULTIGRID    = 6,  /**< Multi-grid V-cycle */
    CFD_SOLVER_DIRECT_LU    = 7   /**< LU decomposition */
} CfdSolverType;

/* ============================================================================
 * Helper Math Macros
 * ============================================================================ */

#define CFD_SQ(x)        ((x)*(x))
#define CFD_MIN(a,b)     ((a) < (b) ? (a) : (b))
#define CFD_MAX(a,b)     ((a) > (b) ? (a) : (b))
#define CFD_CLAMP(x,lo,hi) ((x) < (lo) ? (lo) : ((x) > (hi) ? (hi) : (x)))
#define CFD_SIGN(x)      ((x) >= 0.0 ? 1.0 : -1.0)
#define CFD_ABS(x)       ((x) >= 0.0 ? (x) : -(x))
#define CFD_LERP(a,b,t)  ((a) + (t)*((b)-(a)))
#define CFD_IDX2(i,j,nx)  ((i) + (nx)*(j))
#define CFD_IDX3(i,j,k,nx,ny)  ((i) + (nx)*(j) + (nx)*(ny)*(k))
#define CFD_SAFE_DIV(num,den)  ((den) != 0.0 ? (num)/(den) : 0.0)

/* ============================================================================
 * L3: Fluid Property Constants (Engineering Reference at 20C, 1 atm)
 * ============================================================================ */

/* Water at 20C */
#define CFD_WATER_DENSITY        998.2    /* kg/m3 */
#define CFD_WATER_VISCOSITY      1.002e-3 /* Pa s */
#define CFD_WATER_CONDUCTIVITY   0.598    /* W/(m K) */
#define CFD_WATER_PRANDTL        7.01

/* Air at 20C */
#define CFD_AIR_DENSITY          1.204
#define CFD_AIR_VISCOSITY        1.825e-5
#define CFD_AIR_CONDUCTIVITY     0.0257
#define CFD_AIR_PRANDTL          0.71

/* Engine oil at 20C */
#define CFD_OIL_DENSITY          880.0
#define CFD_OIL_VISCOSITY        0.065
#define CFD_OIL_CONDUCTIVITY     0.144
#define CFD_OIL_PRANDTL          650.0

/* Mercury (liquid metal) at 20C */
#define CFD_HG_DENSITY           13534.0
#define CFD_HG_VISCOSITY         1.55e-3
#define CFD_HG_CONDUCTIVITY      8.34
#define CFD_HG_PRANDTL           0.026

#endif /* CFD_CONFIG_H */

/* Additional physical constants for CFD applications */
#define CFD_GRAVITY         9.80665   /* m/s^2 */
#define CFD_GAS_CONSTANT    287.058   /* J/(kg*K) for dry air */
#define CFD_STEFAN_BOLTZMANN 5.670374419e-8 /* W/(m^2*K^4) */
#define CFD_REF_TEMPERATURE 293.15    /* K */
#define CFD_REF_PRESSURE    101325.0  /* Pa */
#define CFD_MOLAR_MASS_AIR  0.0289644 /* kg/mol */
#define CFD_BOLTZMANN       1.380649e-23 /* J/K */
#define CFD_AVOGADRO        6.02214076e23 /* 1/mol */

/* Grid quality metrics thresholds */
#define CFD_MAX_ASPECT_RATIO    100.0
#define CFD_MIN_ORTHOGONALITY   0.1
#define CFD_MAX_SKEWNESS        0.85
#define CFD_MAX_NONORTHOGONAL   70.0  /* degrees */

/* Time integration parameters */
#define CFD_RK4_A21  0.5
#define CFD_RK4_A32  0.5
#define CFD_RK4_A43  1.0
#define CFD_RK4_B1   0.16666666666666666
#define CFD_RK4_B2   0.33333333333333333
#define CFD_RK4_B3   0.33333333333333333
#define CFD_RK4_B4   0.16666666666666666

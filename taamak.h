/**
 * Copyright (c) 2025 Eyal Levenberg (eylev@dtu.dk)
 *                    Fang Zeyuan    (@veslrs)
 */

#ifndef INCLUDE_TAAMAK_H
#define INCLUDE_TAAMAK_H

#ifndef TMKDEF
#ifdef TAAMAK_STATIC
#define TMKDEF static
#else
#ifdef __cplusplus
#define TMKDEF extern "C"
#else
#define TMKDEF extern
#endif
#endif
#endif

/**
 * api
 * TODO: descriptions
 */

typedef struct {double x, y, z;} tmk_vec3;
typedef struct {tmk_vec3 center; double half_width; double intensity;} tmk_load;
typedef struct {double youngs_modulus, poissons_ratio;} tmk_hs;

/**
 * TODO: description
 * taamak.h API revolves around an objec of type `tmk_model` (a customised 
 * struct type). Via this object, all of the parameters of the model are 
 * specified (number of layers, layer compositions, load, material properties, 
 * slope angle, etcetera), and then one finally passes this object to `tmk_solve` 
 * in order to solve the model and evaluate stress/strain at given points.
 */
typedef struct tmk_model_s tmk_model;

/**
 * Returns a newly allocated `tmk_model` object (or NULL if there was an error, 
 * e.g. out of memory), given the numbers of layers of the pavement system.
 */
TMKDEF void tmk_init(tmk_model *mdl, const int n_layers);
TMKDEF void tmk_set_load(tmk_model *mdl, 
                         const tmk_vec3 center,
                         const double half_width, 
                         const double intensity);
TMKDEF void tmk_set_composition(tmk_model *mdl, 
                                const double *top_depths);
TMKDEF void tmk_set_material_properties(tmk_model *mdl,
                                        const tmk_hs *hs);
TMKDEF void tmk_set_evaluation_points(tmk_model *mdl, 
                                      const int npts, 
                                      const tmk_vec3 *pts);
TMKDEF void tmk_set_slope_angle(tmk_model *mdl, const double degree);
TMKDEF void tmk_set_measurement(tmk_model *mdl, 
                                const int npts, 
                                tmk_vec3 *pts, 
                                double *vals);
TMKDEF void tmk_solve(tmk_model *mdl);
TMKDEF void tmk_reset(tmk_model *mdl);

#endif // INCLUDE_TAAMAK_H

#ifdef TAAMAK_IMPLEMENTATION

#ifdef _WIN32
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#ifndef _CRT_NONSTDC_NO_DEPRECATE
#define _CRT_NONSTDC_NO_DEPRECATE
#endif
#endif

#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "cblas.h"
#include "lapacke.h"

#define TMK_GEOM_EPSILON 0.001

typedef enum {
        TMK_INITIAL = 0,
        TMK_SUCCESS,
        TMK_ALLOCATION_FAILURE,
        TMK_LAPACK_FAILURE,
        TMK_NO_COMPOSITION,
        TMK_NO_LOAD,
} tmk_result;

struct tmk_model_s {
        int         n_layers;
        bool        info_depths;
        double     *top_depths;
        bool        info_mat;
        tmk_hs     *halfspace;
        bool        info_load;
        tmk_load    load;
        bool        info_eval_pts;
        int         npts_eval;
        tmk_vec3   *pts_eval;
        double     *vals_eval;
        bool        info_measure_edge;
        int         npts_measure_edge;
        tmk_vec3   *pts_measure_edge;
        double     *vals_measure_edge;
        bool        info_measure_center;
        int         npts_measure_center;
        tmk_vec3   *pts_measure_center;
        double     *vals_measure_center;
        bool        info_angle;
        double      slope_angle;
        int         info_lapack;
        tmk_result  result;
        // double     *sol;
        // tmk_vec3   *force_boundary;
};

// static const int tmk__const_m_side = 50;
// static const int tmk__const_n_side = 40;

#ifndef TMK_CONST_M_SIDE
#define TMK_CONST_M_SIDE 50
#endif

#ifndef TMK_CONST_N_SIDE
#define TMK_CONST_N_SIDE 50
#endif

/**
 * fundamental solutions
 * TODO: descriptions
 */

static inline double tmk__coefficients_radius(double a, double b, double c)
{
        return sqrt(a * a + b * b + c * c);
}

static inline double tmk__coefficients_concentrated_sxx_x(double x, double z, double c, 
                                                          double r1, double r2, double poissons_ratio)
{
        return x / 8.0 / M_PI / (poissons_ratio - 1.0) * (
                (2.0 * poissons_ratio - 1.0) / pow(r1, 3.0) 
                + (1.0 - 2.0 * poissons_ratio) * (5.0 - 4.0 * poissons_ratio) / pow(r2, 3.0) 
                - 3.0 * x * x / pow(r1, 5.0) 
                - 3.0 * (3.0 - 4.0 * poissons_ratio) * x * x / pow(r2, 5.0) 
                - 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / r2 / (r2 + z + c) / (r2 + z + c) * (
                        3.0 - x * x * (3.0 * r2 + z + c) / r2 / r2 / (r2 + z + c)
                ) 
                + 6.0 * c / pow(r2, 5.0) * (
                        3.0 * c - (3.0 - 2.0 * poissons_ratio) * (z + c) 
                        + 5.0 * x * x * z / r2 / r2
                )
        );
}

static inline double tmk__coefficients_concentrated_syy_x(double x, double y, double z, double c, 
                                                          double r1, double r2, double poissons_ratio)
{
        return x / 8.0 / M_PI / (poissons_ratio - 1) * (
                (1.0 - 2.0 * poissons_ratio) / pow(r1, 3.0) 
                + (1.0 - 2.0 * poissons_ratio) * (3.0 - 4.0 * poissons_ratio) / pow(r2, 3.0) 
                - 3.0 * y * y / pow(r1, 5.0) 
                - 3.0 * (3.0 - 4.0 * poissons_ratio) * y * y / pow(r2, 5.0) 
                - 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / r2 / (r2 + z + c) / (r2 + z + c) * (
                        1.0 - y * y * (3.0 * r2 + z + c) / r2 / r2 / (r2 + z + c)
                ) 
                + 6.0 * c / pow(r2, 5.0) * (
                        c - (1 - 2.0 * poissons_ratio) * (z + c) 
                        + 5.0 * y * y * z / r2 / r2
                )
        );
}

static inline double tmk__coefficients_concentrated_szz_x(double x, double z, double c, 
                                                          double r1, double r2, double poissons_ratio)
{
        return x / 8.0 / M_PI / (poissons_ratio - 1.0) * (
                (1.0 - 2.0 * poissons_ratio) / pow(r1, 3.0) 
                - (1.0 - 2.0 * poissons_ratio) / pow(r2, 3.0) 
                - 3.0 * (z - c) * (z - c) / pow(r1, 5.0) 
                - 3.0 * (3.0 - 4.0 * poissons_ratio) * (z + c) * (z + c) / pow(r2, 5.0) 
                + 6.0 * c / pow(r2, 5.0) * (
                        c + (1.0 - 2.0 * poissons_ratio) * (z + c) 
                        + 5.0 * z * (z + c) * (z + c) / r2 / r2
                )
        );
}

static inline double tmk__coefficients_concentrated_syz_x(double x, double y, double z, double c, 
                                                          double r1, double r2, double poissons_ratio)
{
        return x * y / 8.0 / M_PI / (poissons_ratio - 1.0) * (
                3.0 * (c - z) / pow(r1, 5.0) 
                - 3.0 * (3.0 - 4.0 * poissons_ratio) * (z + c) / pow(r2, 5.0) 
                + 6.0 * c / pow(r2, 5.0) * (
                        1.0 - 2.0 * poissons_ratio 
                        + 5.0 * z * (z + c) / r2 / r2
                )
        );
}

static inline double tmk__coefficients_concentrated_sxz_x(double x, double z, double c, 
                                                          double r1, double r2, double poissons_ratio)
{
        return 1.0 / 8.0 / M_PI / (poissons_ratio - 1.0) * (
                (2.0 * poissons_ratio - 1.0) * (z - c) / pow(r1, 3.0) 
                + (1.0 - 2.0 * poissons_ratio) * (z - c) / pow(r2, 3.0) 
                - 3.0 * x * x * (z - c) / pow(r1, 5.0) 
                - 3.0 * (3.0 - 4.0 * poissons_ratio) * x * x * (z + c) / pow(r2, 5.0) 
                - 6.0 * c / pow(r2, 5.0) * (
                        z * (z + c) 
                        - (1.0 - 2.0 * poissons_ratio) * x * x 
                        - 5.0 * x * x * z * (z + c) / r2 / r2
                )
        );
}

static inline double tmk__coefficients_concentrated_sxy_x(double x, double y, double z, double c, 
                                                          double r1, double r2, double poissons_ratio)
{
        return y / 8.0 / M_PI / (poissons_ratio - 1.0) * (
                (2.0 * poissons_ratio - 1.0) / pow(r1, 3.0) 
                + (1.0 - 2.0 * poissons_ratio) / pow(r2, 3.0) 
                - 3.0 * x * x / pow(r1, 5.0) 
                - 3.0 * (3.0 - 4.0 * poissons_ratio) * x * x / pow(r2, 5.0) 
                - 6.0 * c * z / pow(r2, 5.0) * (1.0 - 5.0 * x * x / r2 / r2) 
                - 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / r2 / (r2 + z + c) / (r2 + z + c) * (
                        1.0 - x * x * (3.0 * r2 + z + c) / r2 / r2 / (r2 + z + c)
                )
        );
}

static inline double tmk__coefficients_concentrated_sxx_y(double x, double y, double z, double c, 
                                                          double r1, double r2, double poissons_ratio)
{
        return y / 8.0 / M_PI / (poissons_ratio - 1.0) * (
                (1.0 - 2.0 * poissons_ratio) / pow(r1, 3.0) 
                + (1.0 - 2.0 * poissons_ratio) * (3.0 - 4.0 * poissons_ratio) / pow(r2, 3.0) 
                - 3.0 * x * x / pow(r1, 5.0) 
                - 3.0 * (3.0 - 4.0 * poissons_ratio) * x * x / pow(r2, 5.0) 
                - 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / r2 / (r2 + z + c) / (r2 + z + c) * (
                        1.0 - x * x * (3.0 * r2 + z + c) / r2 / r2 / (r2 + z + c)
                ) 
                + 6.0 * c / pow(r2, 5.0) * (
                        c - (1.0 - 2.0 * poissons_ratio) * (z + c) 
                        + 5.0 * x * x * z / r2 / r2
                )
        );
}

static inline double tmk__coefficients_concentrated_syy_y(double y, double z, double c, 
                                                          double r1, double r2, double poissons_ratio)
{
        return y / 8.0 / M_PI / (poissons_ratio - 1.0) * (
                (2.0 * poissons_ratio - 1.0) / pow(r1, 3.0) 
                + (1.0 - 2.0 * poissons_ratio) * (5.0 - 4.0 * poissons_ratio) / pow(r2, 3.0) 
                - 3.0 * y * y / pow(r1, 5.0) 
                - 3.0 * (3.0 - 4.0 * poissons_ratio) * y * y / pow(r2, 5.0) 
                - 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / r2 / (r2 + z + c) / (r2 + z + c) * (
                        3.0 - y * y * (3.0 * r2 + z + c) / r2 / r2 / (r2 + z + c)
                ) 
                + 6.0 * c / pow(r2, 5.0) * (
                        3.0 * c - (3.0 - 2.0 * poissons_ratio) * (z + c) 
                        + 5.0 * y * y * z / r2 / r2
                )
        );
}

static inline double tmk__coefficients_concentrated_szz_y(double y, double z, double c, 
                                                          double r1, double r2, double poissons_ratio)
{
        return y / 8.0 / M_PI / (poissons_ratio - 1.0) * (
                (1.0 - 2.0 * poissons_ratio) / pow(r1, 3.0) 
                - (1.0 - 2.0 * poissons_ratio) / pow(r2, 3.0) 
                - 3.0 * (z - c) * (z - c) / pow(r1, 5.0) 
                - 3.0 * (3.0 - 4.0 * poissons_ratio) * (z + c) * (z + c) / pow(r2, 5.0) 
                + 6.0 * c / pow(r2, 5.0) * (
                        c + (1.0 - 2.0 * poissons_ratio) * (z + c) 
                        + 5.0 * z * (z + c) * (z + c) / r2 / r2
                )
        );
}

static inline double tmk__coefficients_concentrated_syz_y(double y, double z, double c, 
                                                          double r1, double r2, double poissons_ratio)
{
        return 1.0 / 8.0 / M_PI / (poissons_ratio - 1.0) * (
                (2.0 * poissons_ratio - 1.0) * (z - c) / pow(r1, 3.0) 
                + (1.0 - 2.0 * poissons_ratio) * (z - c) / pow(r2, 3.0) 
                - 3.0 * (3.0 - 4.0 * poissons_ratio) * y * y * (z + c) / pow(r2, 5.0) 
                - 3.0 * y * y * (z - c) / pow(r1, 5.0) 
                - 6.0 * c / pow(r2, 5.0) * (
                        z * (z + c) 
                        - (1.0 - 2.0 * poissons_ratio) * y * y 
                        - 5.0 * y * y * z * (z + c) / r2 / r2
                )
        );
}

static inline double tmk__coefficients_concentrated_sxz_y(double x, double y, double z, double c, 
                                                          double r1, double r2, double poissons_ratio)
{
        return y * x / 8.0 / M_PI / (poissons_ratio - 1.0) * (
                3.0 * (c - z) / pow(r1, 5.0) 
                - 3.0 * (3.0 - 4.0 * poissons_ratio) * (z + c) / pow(r2, 5.0) 
                + 6.0 * c / pow(r2, 5.0) * (
                        1.0 - 2.0 * poissons_ratio 
                        + 5.0 * z * (z + c) / r2 / r2
                )
        );
}

static inline double tmk__coefficients_concentrated_sxy_y(double x, double y, double z, double c, 
                                                          double r1, double r2, double poissons_ratio)
{
        return x / 8.0 / M_PI / (poissons_ratio - 1.0) * (
                (2.0 * poissons_ratio - 1.0) / pow(r1, 3.0) 
                + (1.0 - 2.0 * poissons_ratio) / pow(r2, 3.0) 
                - 3.0 * y * y / pow(r1, 5.0) 
                - 3.0 * (3.0 - 4.0 * poissons_ratio) * y * y / pow(r2, 5.0) 
                - 6.0 * c * z / pow(r2, 5.0) * (1.0 - 5.0 * y * y / r2 / r2) 
                - 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / r2 / (r2 + z + c) / (r2 + z + c) * (
                        1.0 - y * y * (3.0 * r2 + z + c) / r2 / r2 / (r2 + z + c)
                )
        );
}

static inline double tmk__coefficients_concentrated_sxx_z(double x, double z, double c, 
                                                          double r1, double r2, double poissons_ratio)
{
        return 1.0 / 8.0 / M_PI / (poissons_ratio - 1.0) * (
                (1.0 - 2.0 * poissons_ratio) * (z - c) / pow(r1, 3.0) 
                - 3.0 * x * x * (z - c) / pow(r1, 5.0) 
                - 30.0 * c * x * x * z * (z + c) / pow(r2, 7.0) 
                + (1.0 - 2.0 * poissons_ratio) * (3.0 * (z - c) - 4.0 * poissons_ratio * (z + c)) / pow(r2, 3.0) 
                - (
                        3.0 * x * x * (3.0 - 4.0 * poissons_ratio) * (z - c) 
                        - 6.0 * c * (z + c) * ((1.0 - 2.0 * poissons_ratio) * z 
                        - 2.0 * poissons_ratio * c)
                ) / pow(r2, 5.0) 
                - 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / r2 / (r2 + z + c) * (
                        1.0 - x * x / r2 / (r2 + z + c) 
                        - x * x / r2 / r2
                )
        );
}

static inline double tmk__coefficients_concentrated_syy_z(double y, double z, double c, 
                                                          double r1, double r2, double poissons_ratio)
{
        return 1.0 / 8.0 / M_PI / (poissons_ratio - 1.0) * (
                (1.0 - 2.0 * poissons_ratio) * (z - c) / pow(r1, 3.0) 
                - 3.0 * y * y * (z - c) / pow(r1, 5.0) 
                - 30.0 * c * y * y * z * (z + c) / pow(r2, 7.0) 
                + (1.0 - 2.0 * poissons_ratio) * (3.0 * (z - c) - 4.0 * poissons_ratio * (z + c)) / pow(r2, 3.0) 
                - (
                        3.0 * (3.0 - 4.0 * poissons_ratio) * y * y * (z - c) 
                        - 6.0 * c * (z + c) * ((1.0 - 2.0 * poissons_ratio) * z 
                        - 2.0 * poissons_ratio * c)
                ) / pow(r2, 5.0) 
                - 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / r2 / (r2 + z + c) * (
                        1.0 - y * y / r2 / (r2 + z + c) 
                        - y * y / r2 / r2
                )
        );
}

static inline double tmk__coefficients_concentrated_szz_z(double z, double c, 
                                                          double r1, double r2, double poissons_ratio)
{
        return 1.0 / 8.0 / M_PI / (poissons_ratio - 1.0) * (
                (2.0 * poissons_ratio - 1.0) * (z - c) / pow(r1, 3.0) 
                + (1.0 - 2.0 * poissons_ratio) * (z - c) / pow(r2, 3.0) 
                - 3.0 * pow(z - c, 3.0) / pow(r1, 5.0) 
                - 30.0 * c * z * pow(z + c, 3.0) / pow(r2, 7.0) 
                - (
                        3.0 * (3.0 - 4.0 * poissons_ratio) * z * (z + c) * (z + c) 
                        - 3.0 * c * (z + c) * (5.0 * z - c)
                ) / pow(r2, 5.0)
        );
}

static inline double tmk__coefficients_concentrated_syz_z(double y, double z, double c, 
                                                          double r1, double r2, double poissons_ratio)
{
        return y / 8.0 / M_PI / (poissons_ratio - 1.0) * (
                (2.0 * poissons_ratio - 1.0) / pow(r1, 3.0) 
                + (1.0 - 2.0 * poissons_ratio) / pow(r2, 3.0) 
                - (
                        3.0 * (3.0 - 4.0 * poissons_ratio) * z * (z + c) 
                        - 3.0 * c * (3.0 * z + c)
                ) / pow(r2, 5.0) 
                - 3.0 * (z - c) * (z - c) / pow(r1, 5.0) 
                - 30.0 * c * z * (z + c) * (z + c) / pow(r2, 7.0)
        );
}

static inline double tmk__coefficients_concentrated_sxz_z(double x, double z, double c, 
                                                          double r1, double r2, double poissons_ratio)
{
        return x / 8.0 / M_PI / (poissons_ratio - 1.0) * (
                (2.0 * poissons_ratio - 1.0) / pow(r1, 3.0) 
                + (1.0 - 2.0 * poissons_ratio) / pow(r2, 3.0) 
                - (
                        3.0 * (3.0 - 4.0 * poissons_ratio) * z * (z + c) 
                        - 3.0 * c * (3.0 * z + c)
                ) / pow(r2, 5.0) 
                - 3.0 * (z - c) * (z - c) / pow(r1, 5.0) 
                - 30.0 * c * z * (z + c) * (z + c) / pow(r2, 7.0)
        );
}

static inline double tmk__coefficients_concentrated_sxy_z(double x, double y, double z, double c, 
                                                          double r1, double r2, double poissons_ratio)
{
        return x * y / 8.0 / M_PI / (poissons_ratio - 1.0) * (
                3.0 * (c - z) / pow(r1, 5.0) 
                - 3.0 * (3.0 - 4.0 * poissons_ratio) * (z - c) / pow(r2, 5.0) 
                - 30.0 * c * z * (z + c) / pow(r2, 7.0) 
                + 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / r2 / r2 / (r2 + z + c) * (
                        1.0 / (r2 + z + c) + 1.0 / r2
                )
        );
}

static inline double tmk__coefficients_concentrated_ux_x(double x, double z, double c, 
                                                         double r1, double r2, 
                                                         double youngs_modulus, double poissons_ratio)
{
        return (1.0 + poissons_ratio) / 8.0 / M_PI / youngs_modulus / (1.0 - poissons_ratio) * (
                (3.0 - 4.0 * poissons_ratio) / r1 
                + 1.0 / r2 
                + (3.0 - 4.0 * poissons_ratio) * x * x / pow(r2, 3.0) 
                + 2.0 * c * z / pow(r2, 3.0) * (1.0 - 3.0 * x * x / r2 / r2) 
                + x * x / pow(r1, 3.0) + 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / (r2 + z + c) * (1.0 - x * x / r2 / (r2 + z + c))
        );
}

static inline double tmk__coefficients_concentrated_uy_x(double x, double y, double z, double c, 
                                                         double r1, double r2, 
                                                         double youngs_modulus, double poissons_ratio)
{
        return (1.0 + poissons_ratio) * x * y / 8.0 / M_PI / youngs_modulus / (1.0 - poissons_ratio) * (
                1.0 / pow(r1, 3.0) 
                + (3.0 - 4.0 * poissons_ratio) / pow(r2, 3.0) 
                - 6.0 * c * z / pow(r2, 5.0) 
                - 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / r2 / (r2 + z + c) / (r2 + z + c)
        );
}

static inline double tmk__coefficients_concentrated_uz_x(double x, double z, double c, 
                                                         double r1, double r2, 
                                                         double youngs_modulus, double poissons_ratio)
{
        return (1.0 + poissons_ratio) * x / 8.0 / M_PI / youngs_modulus / (1.0 - poissons_ratio) * (
                (z - c) / pow(r1, 3.0) 
                + (3.0 - 4.0 * poissons_ratio) * (z - c) / pow(r2, 3.0) 
                - 6.0 * c * z * (z + c) / pow(r2, 5.0) 
                + 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / r2 / (r2 + z + c)
        );
}

static inline double tmk__coefficients_concentrated_ux_y(double x, double y, double z, double c, 
                                                         double r1, double r2, 
                                                         double youngs_modulus, double poissons_ratio)
{
        return (1.0 + poissons_ratio) * x * y / 8.0 / M_PI / youngs_modulus / (1.0 - poissons_ratio) * (
                1.0 / pow(r1, 3.0) 
                + (3.0 - 4.0 * poissons_ratio) / pow(r2, 3.0) 
                - 6.0 * c * z / pow(r2, 5.0) 
                - 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / r2 / (r2 + z + c) / (r2 + z + c)
        );
}

static inline double tmk__coefficients_concentrated_uy_y(double y, double z, double c, 
                                                         double r1, double r2, 
                                                         double youngs_modulus, double poissons_ratio)
{
        return (1.0 + poissons_ratio) / 8.0 / M_PI / youngs_modulus / (1.0 - poissons_ratio) * (
                (3.0 - 4.0 * poissons_ratio) / r1 
                + 1.0 / r2 + y * y / pow(r1, 3.0) 
                + (3.0 - 4.0 * poissons_ratio) * y * y / pow(r2, 3.0) 
                + 2.0 * c * z / pow(r2, 3.0) * (1.0 - 3.0 * y * y / r2 / r2) 
                + 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / (r2 + z + c) * (1.0 - y * y / r2 / (r2 + z + c))
        );
}

static inline double tmk__coefficients_concentrated_uz_y(double y, double z, double c, 
                                                         double r1, double r2, 
                                                         double youngs_modulus, double poissons_ratio)
{
        return (1.0 + poissons_ratio) * y / 8.0 / M_PI / youngs_modulus / (1.0 - poissons_ratio) * (
                (z - c) / pow(r1, 3.0) 
                + (3.0 - 4.0 * poissons_ratio) * (z - c) / pow(r2, 3.0) 
                - 6.0 * c * z * (z + c) / pow(r2, 5.0) 
                + 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / r2 / (r2 + z + c)
        );
}

static inline double tmk__coefficients_concentrated_ux_z(double x, double z, double c, 
                                                         double r1, double r2, 
                                                         double youngs_modulus, double poissons_ratio)
{
        return (1.0 + poissons_ratio) * x / 8.0 / M_PI / youngs_modulus / (1.0 - poissons_ratio) * (
                (z - c) / pow(r1, 3) 
                + (3.0 - 4.0 * poissons_ratio) * (z - c) / pow(r2, 3.0) 
                - 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / r2 / (r2 + z + c) 
                + 6.0 * c * z * (z + c) / pow(r2, 5.0)
        );
}

static inline double tmk__coefficients_concentrated_uy_z(double y, double z, double c, 
                                                         double r1, double r2, 
                                                         double youngs_modulus, double poissons_ratio)
{
        return (1.0 + poissons_ratio) * y / 8.0 / M_PI / youngs_modulus / (1.0 - poissons_ratio) * (
                (z - c) / pow(r1, 3.0) 
                + (3.0 - 4.0 * poissons_ratio) * (z - c) / pow(r2, 3.0) 
                - 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / r2 / (r2 + z + c) 
                + 6.0 * c * z * (z + c) / pow(r2, 5.0)
        );
}

static inline double tmk__coefficients_concentrated_uz_z(double z, double c, 
                                                         double r1, double r2, 
                                                         double youngs_modulus, double poissons_ratio)
{
        return (1.0 + poissons_ratio) / 8.0 / M_PI / youngs_modulus / (1.0 - poissons_ratio) * (
                (3.0 - 4.0 * poissons_ratio) / r1 
                + (8.0 * (1.0 - poissons_ratio) * (1.0 - poissons_ratio) 
                - (3.0 - 4.0 * poissons_ratio)) / r2 + (z - c) * (z - c) / pow(r1, 3.0) 
                + ((3.0 - 4.0 * poissons_ratio) * (z + c) * (z + c) - 2.0 * c * z) / pow(r2, 3.0) 
                + 6.0 * c * z * (z + c) * (z + c) / pow(r2, 5.0)
        );
}

static void tmk__coefficients_concentrated_stress_x(double *coefficients, 
                                                    tmk_vec3 *evaluation_point, 
                                                    tmk_vec3 *force_point, 
                                                    const tmk_hs *halfspace)
{
        double x = evaluation_point->x - force_point->x;
        double y = evaluation_point->y - force_point->y;
        double z = evaluation_point->z;

        double c = force_point->z;

        double r1 = tmk__coefficients_radius(x, y, z - c);
        double r2 = tmk__coefficients_radius(x, y, z + c);

        double poissons_ratio = halfspace->poissons_ratio;

        coefficients[0] = tmk__coefficients_concentrated_sxx_x(x, z, c, r1, r2, poissons_ratio);
        coefficients[1] = tmk__coefficients_concentrated_syy_x(x, y, z, c, r1, r2, poissons_ratio);
        coefficients[2] = tmk__coefficients_concentrated_szz_x(x, z, c, r1, r2, poissons_ratio);
        coefficients[3] = tmk__coefficients_concentrated_syz_x(x, y, z, c, r1, r2, poissons_ratio);
        coefficients[4] = tmk__coefficients_concentrated_sxz_x(x, z, c, r1, r2, poissons_ratio);
        coefficients[5] = tmk__coefficients_concentrated_sxy_x(x, y, z, c, r1, r2, poissons_ratio);
}

static void tmk__coefficients_concentrated_stress_y(double *coefficients, 
                                                    tmk_vec3 *evaluation_point, 
                                                    tmk_vec3 *force_point, 
                                                    const tmk_hs *halfspace)
{
        double x = evaluation_point->x - force_point->x;
        double y = evaluation_point->y - force_point->y;
        double z = evaluation_point->z;

        double c = force_point->z;

        double r1 = tmk__coefficients_radius(x, y, z - c);
        double r2 = tmk__coefficients_radius(x, y, z + c);

        double poissons_ratio = halfspace->poissons_ratio;

        coefficients[0] = tmk__coefficients_concentrated_sxx_y(x, y, z, c, r1, r2, poissons_ratio);
        coefficients[1] = tmk__coefficients_concentrated_syy_y(y, z, c, r1, r2, poissons_ratio);
        coefficients[2] = tmk__coefficients_concentrated_szz_y(y, z, c, r1, r2, poissons_ratio);
        coefficients[3] = tmk__coefficients_concentrated_syz_y(y, z, c, r1, r2, poissons_ratio);
        coefficients[4] = tmk__coefficients_concentrated_sxz_y(x, y, z, c, r1, r2, poissons_ratio);
        coefficients[5] = tmk__coefficients_concentrated_sxy_y(x, y, z, c, r1, r2, poissons_ratio);
}

static void tmk__coefficients_concentrated_stress_z(double *coefficients, 
                                                    tmk_vec3 *evaluation_point, 
                                                    tmk_vec3 *force_point, 
                                                    const tmk_hs *halfspace)
{
        double x = evaluation_point->x - force_point->x;
        double y = evaluation_point->y - force_point->y;
        double z = evaluation_point->z;

        double c = force_point->z;

        double r1 = tmk__coefficients_radius(x, y, z - c);
        double r2 = tmk__coefficients_radius(x, y, z + c);

        double poissons_ratio = halfspace->poissons_ratio;

        coefficients[0] = tmk__coefficients_concentrated_sxx_z(x, z, c, r1, r2, poissons_ratio);
        coefficients[1] = tmk__coefficients_concentrated_syy_z(y, z, c, r1, r2, poissons_ratio);
        coefficients[2] = tmk__coefficients_concentrated_szz_z(z, c, r1, r2, poissons_ratio);
        coefficients[3] = tmk__coefficients_concentrated_syz_z(y, z, c, r1, r2, poissons_ratio);
        coefficients[4] = tmk__coefficients_concentrated_sxz_z(x, z, c, r1, r2, poissons_ratio);
        coefficients[5] = tmk__coefficients_concentrated_sxy_z(x, y, z, c, r1, r2, poissons_ratio);
}

static void tmk__coefficients_concentrated_displacement_x(double *coefficients, 
                                                          tmk_vec3 *evaluation_point, 
                                                          tmk_vec3 *force_point, 
                                                          tmk_hs *halfspace)
{
        double x = evaluation_point->x - force_point->x;
        double y = evaluation_point->y - force_point->y;
        double z = evaluation_point->z;

        double c = force_point->z;

        double r1 = tmk__coefficients_radius(x, y, z - c);
        double r2 = tmk__coefficients_radius(x, y, z + c);

        double youngs_modulus = halfspace->youngs_modulus;
        double poissons_ratio = halfspace->poissons_ratio;

        coefficients[0] = tmk__coefficients_concentrated_ux_x(x, z, c, r1, r2, youngs_modulus, poissons_ratio);
        coefficients[1] = tmk__coefficients_concentrated_uy_x(x, y, z, c, r1, r2, youngs_modulus, poissons_ratio);
        coefficients[2] = tmk__coefficients_concentrated_uz_x(x, z, c, r1, r2, youngs_modulus, poissons_ratio);
}

static void tmk__coefficients_concentrated_displacement_y(double *coefficients, 
                                                          tmk_vec3 *evaluation_point, 
                                                          tmk_vec3 *force_point, 
                                                          tmk_hs *halfspace)
{
        double x = evaluation_point->x - force_point->x;
        double y = evaluation_point->y - force_point->y;
        double z = evaluation_point->z;

        double c = force_point->z;

        double r1 = tmk__coefficients_radius(x, y, z - c);
        double r2 = tmk__coefficients_radius(x, y, z + c);

        double youngs_modulus = halfspace->youngs_modulus;
        double poissons_ratio = halfspace->poissons_ratio;

        coefficients[0] = tmk__coefficients_concentrated_ux_y(x, y, z, c, r1, r2, youngs_modulus, poissons_ratio);
        coefficients[1] = tmk__coefficients_concentrated_uy_y(y, z, c, r1, r2, youngs_modulus, poissons_ratio);
        coefficients[2] = tmk__coefficients_concentrated_uz_y(y, z, c, r1, r2, youngs_modulus, poissons_ratio);
}

static void tmk__coefficients_concentrated_displacement_z(double *coefficients, 
                                                          tmk_vec3 *evaluation_point, 
                                                          tmk_vec3 *force_point, 
                                                          tmk_hs *halfspace)
{
        double x = evaluation_point->x - force_point->x;
        double y = evaluation_point->y - force_point->y;
        double z = evaluation_point->z;

        double c = force_point->z;

        double r1 = tmk__coefficients_radius(x, y, z - c);
        double r2 = tmk__coefficients_radius(x, y, z + c);

        double youngs_modulus = halfspace->youngs_modulus;
        double poissons_ratio = halfspace->poissons_ratio;

        coefficients[0] = tmk__coefficients_concentrated_ux_z(x, z, c, r1, r2, youngs_modulus, poissons_ratio);
        coefficients[1] = tmk__coefficients_concentrated_uy_z(y, z, c, r1, r2, youngs_modulus, poissons_ratio);
        coefficients[2] = tmk__coefficients_concentrated_uz_z(z, c, r1, r2, youngs_modulus, poissons_ratio);
}

/**
 * Half-space solution for a squared stress patch
 * TODO: descriptions
 */

static inline double tmk__coefficients_distributed_a_xx(double x, double y, double z, 
                                                        double poissons_ratio)
{
        return (
                2.0 * poissons_ratio * (atan2(x, y) 
                - atan2(z * x, y * tmk__coefficients_radius(x, y, z))) 
                + atan2(y, x) - atan2(z * y, x * tmk__coefficients_radius(x, y, z)) 
                - x * y * z / (tmk__coefficients_radius(x, y, z) * (x * x + z * z))
        );
}

static inline double tmk__coefficients_distributed_a_yy(double x, double y, double z, 
                                                        double poissons_ratio)
{
        return (
                2.0 * poissons_ratio * (atan2(y, x) - atan2(z * y, x * tmk__coefficients_radius(x, y, z))) 
                + atan2(x, y) 
                - atan2(z * x, y * tmk__coefficients_radius(x, y, z)) 
                - x * y * z / (tmk__coefficients_radius(x, y, z) * (y * y + z * z))
        );
}

static inline double tmk__coefficients_distributed_a_zz(double x, double y, double z)
{
        return (
                atan2(y, x) + atan2(x, y) 
                - atan2(z * y, x * tmk__coefficients_radius(x, y, z)) 
                - atan2(z * x, y * tmk__coefficients_radius(x, y, z)) 
                + x * y * z / (tmk__coefficients_radius(x, y, z) * (y * y + z * z)) 
                + x * y * z / (tmk__coefficients_radius(x, y, z) * (x * x + z * z))
        );
}

static inline double tmk__coefficients_distributed_a_yz(double x, double y, double z)
{
        return -1.0 * z * z * x / tmk__coefficients_radius(x, y, z) / (y * y + z * z);
}

static inline double tmk__coefficients_distributed_a_xz(double x, double y, double z)
{
        return -1.0 * z * z * y / tmk__coefficients_radius(x, y, z) / (x * x + z * z);
}

static inline double tmk__coefficients_distributed_a_xy(double x, double y, double z, 
                                                        double poissons_ratio)
{
        return (
                (1.0 - 2.0 * poissons_ratio) * log(tmk__coefficients_radius(x, y, z) + z) 
                + z / tmk__coefficients_radius(x, y, z)
        );
}

static inline double tmk__coefficients_distributed_b_x(double x, double y, double z, 
                                                       double poissons_ratio)
{
        return (
                (2.0 * poissons_ratio - 1.0) * (
                        y * log(tmk__coefficients_radius(x, y, z) + z) 
                        + z * log(tmk__coefficients_radius(x, y, z) + y) 
                        - 2.0 * x * atan2(x, tmk__coefficients_radius(x, y, z) + y + z)
                ) 
                - z * log(tmk__coefficients_radius(x, y, z) + y)
        );
}

static inline double tmk__coefficients_distributed_b_y(double x, double y, double z, 
                                                       double poissons_ratio)
{
        return (
                (2.0 * poissons_ratio - 1.0) * (
                        x * log(tmk__coefficients_radius(x, y, z) + z) 
                        + z * log(tmk__coefficients_radius(x, y, z) + x) 
                        - 2.0 * y * atan2(y, tmk__coefficients_radius(x, y, z) + x + z)
                ) 
                - z * log(tmk__coefficients_radius(x, y, z) + x)
        );
}

static inline double tmk__coefficients_distributed_b_z(double x, double y, double z, 
                                                       double poissons_ratio)
{
        return (
                2.0 * (1.0 - poissons_ratio) * (
                        y * log(tmk__coefficients_radius(x, y, z) + x) 
                        + x * log(tmk__coefficients_radius(x, y, z) + y) 
                        + 2.0 * z * atan2((sqrt(x * x + z * z) - x) * (tmk__coefficients_radius(x, y, z) - sqrt(x * x + z * z)), z * y)
                ) 
                + z * atan2(x * y, tmk__coefficients_radius(x, y, z) * z)
        );
}

static void tmk__coefficients_distributed_stress(double *coefficients, 
                                                 tmk_vec3 *evaluation_point, 
                                                 const tmk_load *load, 
                                                 const tmk_hs *halfspace)
{
        double x = evaluation_point->x - load->center.x;
        double y = evaluation_point->y - load->center.y;
        double z = evaluation_point->z - load->center.z;

        double b = load->half_width;

        double poissons_ratio = halfspace->poissons_ratio;

        coefficients[0] = (
                tmk__coefficients_distributed_a_xx(x + b, y + b, z, poissons_ratio) 
                + tmk__coefficients_distributed_a_xx(x - b, y - b, z, poissons_ratio) 
                - tmk__coefficients_distributed_a_xx(x - b, y + b, z, poissons_ratio) 
                - tmk__coefficients_distributed_a_xx(x + b, y - b, z, poissons_ratio)
        ) / 2.0 / M_PI;
        coefficients[1] = (
                tmk__coefficients_distributed_a_yy(x + b, y + b, z, poissons_ratio) 
                + tmk__coefficients_distributed_a_yy(x - b, y - b, z, poissons_ratio) 
                - tmk__coefficients_distributed_a_yy(x - b, y + b, z, poissons_ratio) 
                - tmk__coefficients_distributed_a_yy(x + b, y - b, z, poissons_ratio)
        ) / 2.0 / M_PI;
        coefficients[2] = (
                tmk__coefficients_distributed_a_zz(x + b, y + b, z) 
                + tmk__coefficients_distributed_a_zz(x - b, y - b, z) 
                - tmk__coefficients_distributed_a_zz(x - b, y + b, z) 
                - tmk__coefficients_distributed_a_zz(x + b, y - b, z)
        ) / 2.0 / M_PI;
        coefficients[3] = (
                tmk__coefficients_distributed_a_yz(x + b, y + b, z) 
                + tmk__coefficients_distributed_a_yz(x - b, y - b, z) 
                - tmk__coefficients_distributed_a_yz(x - b, y + b, z) 
                - tmk__coefficients_distributed_a_yz(x + b, y - b, z)
        ) / 2.0 / M_PI;
        coefficients[4] = (
                tmk__coefficients_distributed_a_xz(x + b, y + b, z) 
                + tmk__coefficients_distributed_a_xz(x - b, y - b, z) 
                - tmk__coefficients_distributed_a_xz(x - b, y + b, z) 
                - tmk__coefficients_distributed_a_xz(x + b, y - b, z)
        ) / 2.0 / M_PI;
        coefficients[5] = (
                tmk__coefficients_distributed_a_xy(x + b, y + b, z, poissons_ratio) 
                + tmk__coefficients_distributed_a_xy(x - b, y - b, z, poissons_ratio) 
                - tmk__coefficients_distributed_a_xy(x - b, y + b, z, poissons_ratio) 
                - tmk__coefficients_distributed_a_xy(x + b, y - b, z, poissons_ratio)
        ) / 2.0 / M_PI;
}

static void tmk__coefficients_distributed_displacement(double *coefficients, 
                                                       tmk_vec3 *evaluation_point, 
                                                       tmk_load *load, 
                                                       tmk_hs *halfspace)
{
        double x = evaluation_point->x - load->center.x;
        double y = evaluation_point->y - load->center.y;
        double z = evaluation_point->z - load->center.z;

        double b = load->half_width;

        double youngs_modulus = halfspace->youngs_modulus;
        double poissons_ratio = halfspace->poissons_ratio;

        coefficients[0] = (1.0 + poissons_ratio) / 2.0 / M_PI / youngs_modulus * (
                tmk__coefficients_distributed_b_x(x + b, y + b, z, poissons_ratio) 
                + tmk__coefficients_distributed_b_x(x - b, y - b, z, poissons_ratio) 
                - tmk__coefficients_distributed_b_x(x - b, y + b, z, poissons_ratio) 
                - tmk__coefficients_distributed_b_x(x + b, y - b, z, poissons_ratio)
        );
        coefficients[1] = (1.0 + poissons_ratio) / 2.0 / M_PI / youngs_modulus * (
                tmk__coefficients_distributed_b_y(x + b, y + b, z, poissons_ratio) 
                + tmk__coefficients_distributed_b_y(x - b, y - b, z, poissons_ratio) 
                - tmk__coefficients_distributed_b_y(x - b, y + b, z, poissons_ratio) 
                - tmk__coefficients_distributed_b_y(x + b, y - b, z, poissons_ratio)
        );
        coefficients[2] = (1.0 + poissons_ratio) / 2.0 / M_PI / youngs_modulus * (
                tmk__coefficients_distributed_b_z(x + b, y + b, z, poissons_ratio) 
                + tmk__coefficients_distributed_b_z(x - b, y - b, z, poissons_ratio) 
                - tmk__coefficients_distributed_b_z(x - b, y + b, z, poissons_ratio) 
                - tmk__coefficients_distributed_b_z(x + b, y - b, z, poissons_ratio)
        );
}

// static double tmk__coefficients_distributed_uz(tmk_vec3 *evaluation_point, 
//                                                tmk_load *load, 
//                                                tmk_hs *halfspace)
// {
//         double x = evaluation_point->x - load->center.x;
//         double y = evaluation_point->y - load->center.y;
//         double z = evaluation_point->z - load->center.z;

//         double b = load->half_width;

//         double youngs_modulus = halfspace->youngs_modulus;
//         double poissons_ratio = halfspace->poissons_ratio;

//         double uz = (1.0 + poissons_ratio) / 2.0 / M_PI / youngs_modulus * (
//                 tmk__coefficients_distributed_b_z(x + b, y + b, z, poissons_ratio) 
//                 + tmk__coefficients_distributed_b_z(x - b, y - b, z, poissons_ratio) 
//                 - tmk__coefficients_distributed_b_z(x - b, y + b, z, poissons_ratio) 
//                 - tmk__coefficients_distributed_b_z(x + b, y - b, z, poissons_ratio)
//         );

//         return uz;
// }

/**
 * Tensor transformations
 * TODO: descriptions
 */

 static void tmk__coefficients_tensor_transformation_y(double *tensor, double degree)
{
        double radian = degree * M_PI / 180.0;

        double xx = tensor[0];
        double zz = tensor[2];
        double yz = tensor[3];
        double xz = tensor[4];
        double xy = tensor[5];

        tensor[2] = zz * cos(radian) * cos(radian) 
                    + xx * sin(radian) * sin(radian) 
                    + xz * sin(2.0 * radian);
        tensor[3] = yz * cos(radian) + xy * sin(radian);
        tensor[4] = xz * cos(2.0 * radian) + (xx - zz) * cos(radian) * sin(radian);
}

// static void tmk__coefficients_tensor_transformation_y_full(double *tensor, double degree)
// {
//         double radian = degree * M_PI / 180.0;

//         double xx = tensor[0];
//         double zz = tensor[2];
//         double yz = tensor[3];
//         double xz = tensor[4];
//         double xy = tensor[5];

//         tensor[0] = xx * cos(radian) * cos(radian) 
//                     - 2.0 * xz * cos(radian) * sin(radian) 
//                     + zz * sin(radian) * sin(radian);
//         tensor[2] = zz * cos(radian) * cos(radian) 
//                     + xx * sin(radian) * sin(radian) 
//                     + xz * sin(2.0 * radian);
//         tensor[3] = yz * cos(radian) + xy * sin(radian);
//         tensor[4] = xz * cos(2.0 * radian) + (xx - zz) * cos(radian) * sin(radian);
//         tensor[5] = xy * cos(radian) - yz * sin(radian);
// }

// static void tmk__coefficients_tensor_transformation_z_full(double *tensor, double degree)
// {
//         double radian = degree * M_PI / 180.0;

//         double xx = tensor[0];
//         double yy = tensor[1];
//         double yz = tensor[3];
//         double xz = tensor[4];
//         double xy = tensor[5];

//         tensor[0] = xx * cos(radian) * cos(radian) 
//                     + yy * sin(radian) * sin(radian) 
//                     + xy * sin(2.0 * radian);
//         tensor[1] = yy * cos(radian) * cos(radian) 
//                     - 2.0 * xy * cos(radian) * sin(radian) 
//                     + xx * sin(radian) * sin(radian);
//         tensor[3] = yz * cos(radian) - xz * sin(radian);
//         tensor[4] = xz * cos(radian) + yz * sin(radian);
//         tensor[5] = xy * cos(2.0 * radian) + (yy - xx) * cos(radian) * sin(radian);
// }

TMKDEF void tmk_init(tmk_model *mdl, const int n_layers)
{
        mdl->n_layers = n_layers;

        mdl->info_depths = false;
        mdl->info_load = false;
        mdl->info_mat = false;
        mdl->info_eval_pts = false;
        mdl->info_measure_edge = false;
        mdl->info_measure_center = false;
        mdl->info_angle = false;
        mdl->info_lapack = -999;        // random large negative number to initialise lapack info status
        mdl->result = TMK_INITIAL;        
}

TMKDEF void tmk_set_load(tmk_model *mdl, const tmk_vec3 center,
                         const double half_width, const double intensity)
{
        mdl->load = (tmk_load){
                .center = center,
                .half_width = half_width,
                .intensity = intensity,
        };
        mdl->info_load = true;
}

TMKDEF void tmk_set_composition(tmk_model *mdl, 
                                const double *top_depths)
{
        int len = sizeof(*top_depths) / sizeof(top_depths);
        if (len != mdl->n_layers) 
                mdl->info_depths = false;
        else {
                mdl->top_depths = (double *)top_depths;
                mdl->info_depths = true;
        }
}

TMKDEF void tmk_set_material_properties(tmk_model *mdl,
                                        const tmk_hs *hs)
{
        mdl->halfspace = (tmk_hs *)hs;
        mdl->info_mat = true;
}

TMKDEF void tmk_set_evaluation_points(tmk_model *mdl, const int npts, 
                                      const tmk_vec3 *pts)
{
        mdl->npts_eval = npts;
        mdl->pts_eval = (tmk_vec3 *)pts;
        mdl->info_eval_pts = true;
}

TMKDEF void tmk_set_slope_angle(tmk_model *mdl, const double degree)
{
        mdl->slope_angle = degree;
        mdl->info_angle = true;
}

// TMKDEF void tmk_set_measurement(tmk_model *mdl, 
//                                 const int npts, 
//                                 tmk_vec3 *pts, 
//                                 double *vals)
// {
//         // TODO: implementation
// }

/**
 * Function to solve a single layer half-space system.
 * 
 * Input:
 *   pre-allocated array for solutions, `double *sol`,
 *   pre-allocated array for points on force boundary, `tmk_vec3 force_boundary`,
 *   pointer to half-space parameters, `tmk_hs *halfspace`, 
 *     (note that half-space parameters are stored in an array in the model type
 *      `tmk_model` due to the compatibility concern for multi-layer systems, thus 
 *      in the single layer case, developers should pass the first element of 
 *      such array into this function.)
 *   pointer to load parameters, `tmk_load *load`,
 *   slope angle in [deg], `double alpha`,
 *   number of points on each side of the simulated boundary grid, `const int mside`,
 *   number of points on each side of the force boundary grid, `const int nside`,
 *   pointer to general result info parameter, `tmk_result *result`,
 *   pointer to lapack result info parameter, `int *info_lapack`
 * 
 * Output:
 *   the function modify `*sol` and `*force_boundary` in place, to store solutions
 *   and force boundary geometry for further use.
 */
static void tmk__solve_single(double *sol, tmk_vec3 *force_boundary, 
                              const tmk_hs *halfspace, const tmk_load *load, 
                              const double alpha, const int mside, const int nside,
                              tmk_result *result, int *info_lapack);

/**
 * Function to evaluate stresses and displacements at given list of points in a 
 * single layer pavement system. 
 * 
 * Inputs:
 *   list of points, `tmk_vec3 *pts`
 *   number of evaluation points, `const int npts`
 *   pointer to load config, `tmk_load *load`
 *   pointer to halfspace config, `tmk_hs *halfspace`
 *   list of solutions from solver, `double *sol`
 *   list of force boundary points initiated in solver, `tmk_vec3 *force_boundary`
 *   number of points on force boundary, `const int n_force`
 * 
 * Output:
 *   print evaluation results in stdout
 */
static void tmk__evaluate_single(tmk_vec3 *pts, const int npts, 
                                 tmk_load *load, tmk_hs *halfspace,
                                 double *sol, tmk_vec3 *force_boundary, const int n_forces);

TMKDEF void tmk_solve(tmk_model *mdl)
{
        // TODO: implementation
        if (mdl->info_load == false) {
                mdl->result = TMK_NO_LOAD;
                return;
        }
        const int n = TMK_CONST_N_SIDE * TMK_CONST_N_SIDE;
        double *sol = calloc(3 * n, sizeof(*sol));
        tmk_vec3 *force_boundary = calloc(n, sizeof(*force_boundary));
        tmk__solve_single(sol, force_boundary, 
                          &mdl->halfspace[0], &mdl->load, 
                          mdl->slope_angle, TMK_CONST_M_SIDE, TMK_CONST_N_SIDE,
                          &mdl->result, &mdl->info_lapack);
        if (mdl->info_eval_pts == true) 
                tmk__evaluate_single(mdl->pts_eval, mdl->npts_eval, 
                                     &mdl->load, &mdl->halfspace[0], 
                                     sol, force_boundary, n);
        free(sol);
        free(force_boundary);
        // printf("Call `tmk_log` to check results.\n");
}

static void tmk__log_evaluation_single(tmk_model *mdl)
{
        int npts = mdl->npts_eval;
        for (int i = 0; i < npts; ++i) {
                for (int j = 0; j < 9; ++j) {
                        int idx = i * 9 + j;
                        printf("%g\t", mdl->vals_eval[idx]);
                }
                printf("\n");
        }
}

TMKDEF void tmk_log(tmk_model *mdl)
{
        if (mdl == NULL) {
                perror("Fail to initialise model...");
                return;
        }

        if (mdl->result == TMK_SUCCESS) tmk__log_evaluation_single(mdl);
        else perror("Solver error...");

        // TODO: implementation
        /**
         * ---
         * LOG
         * ---
         * Initial Model: 
         * 
         *                      |---- LOAD ----|
         *     ----------------------------------------------- z = 0.0[mm]
         *     alpha:   /    E  = 500.0[MPa] / E = unknown
         *   75.0[deg] /     nu = 0.30[-] / nu = unknown
         *            / 
         * 
         * Load:
         *      center   : (x, y, z)
         *      width    : 2b[mm]
         *      intensity: q[MPa]
         * 
         * Calculations: 
         *      Evaluation of stresses and displacements at given points
         * 
         * Results:
         * 
         *      
         */
}

TMKDEF void tmk_reset(tmk_model *mdl)
{
        if (mdl == NULL) return;
        // free(mdl);
}

/**
 * Function to add epsilon to points close to the corner of the loading patch, 
 * due to the singularities. The function takes pointer of a single point, 
 * `tmk_vec3 *point`, and the pointer of load parameters `tmk_load *load`, then 
 * modify the given point in place.
 */
static void tmk__util_correct_singularity(tmk_vec3 *point, const tmk_load *load)
{
        if (fabs(fabs(point->x - load->center.x) - load->half_width) < 1.0e-9)
                point->x += TMK_GEOM_EPSILON;
        
        if (fabs(fabs(point->y - load->center.y) - load->half_width) < 1.0e-9)
                point->y += TMK_GEOM_EPSILON;
}

static void tmk__util_fill_lhs(double *lhs, const int m, const int n, 
                               double *s33_px, double *s33_py, double *s33_pz, 
                               double *s23_px, double *s23_py, double *s23_pz, 
                               double *s13_px, double *s13_py, double *s13_pz)
{
        double *blocks[3][3] = {
                {s33_px, s33_py, s33_pz},
                {s23_px, s23_py, s23_pz},
                {s13_px, s13_py, s13_pz},
        };

        for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                        double *src = blocks[i][j];

                        for (int k = 0; k < n; ++k) {
                                double *src_col = src + k * m;
                                double *dst = lhs + (j * n + k) * (3 * m) + i * m;
                                memcpy(dst, src_col, m * sizeof(double));
                        }
                }
        }
}

static void tmk__geometry_grid_uniform(tmk_vec3 *grid, 
                                       const int npts_side,
                                       const double h, 
                                       const double degree,
                                       const tmk_load *load)
{
        double radian = degree * M_PI / 180.0;
        double s = sin(radian);
        double hx = h / s;

        for (int i = 0; i < npts_side; ++i) {
                double x = hx * (double)(i + 1 - npts_side);

                for (int j = 0; j < npts_side; ++j) {
                        double y = h * (double)(2 * j + 1 - npts_side) / 2.0;

                        int idx = i * npts_side + j;

                        grid[idx].x = x * cos(-1.0 * radian);
                        grid[idx].y = y;
                        grid[idx].z = x * sin (-1.0 * radian);
                        tmk__util_correct_singularity(&grid[idx], load);
                }
        }
}

static void tmk__geometry_grid_offset(tmk_vec3 *grid, 
                                      const int npts, 
                                      const double beta, 
                                      const double x0)
{
        double offset = beta * x0;

        for (int i = 0; i < npts; ++i) {
                grid[i].x -= offset;
        }
}

/**
 * Function to set configure geometrical properties of simulated boundary and
 * force boundary. The formulae to derive these parameters are given in Levenberg's
 * work and Fang's MSc thesis (see reference).
 * 
 * Input: 
 *   load center in x-direction, `const double x0`,
 *   pointer to pre-defined step size of simulated boundary grid, `double *hm`,
 *   pointer to pre-defined step size of force boundary grid, `double *hn`,
 *   pointer to pre-defined offset coefficient for force boundary, `double *hm`
 * 
 * Output:
 *   the function modify `hm`, `hn` and `beta` in place.
 */
static void tmk__geometry_grid_properties(const double x0, 
                                          double *hm, double *hn, double *beta);

static void tmk__solve_single(double *sol, tmk_vec3 *force_boundary, 
                              const tmk_hs *halfspace, const tmk_load *load, 
                              const double alpha, const int mside, const int nside,
                              tmk_result *result, int *info_lapack)
{
        printf("Initialising model...\n");

        double x0 = load->center.x;
        double q = load->intensity;
        
        double hm, hn, beta;
        tmk__geometry_grid_properties(x0, &hm, &hn, &beta);
        
        const int m = mside * mside;
        const int n = nside * nside;
        
        tmk_vec3 *simulated_boundary = calloc(m, sizeof(*simulated_boundary));

        double *s33_q = calloc(m, sizeof(*s33_q));
        double *s23_q = calloc(m, sizeof(*s23_q));
        double *s13_q = calloc(m, sizeof(*s13_q));
        
        double *s33_px = calloc(m * n, sizeof(double));
        double *s23_px = calloc(m * n, sizeof(double));
        double *s13_px = calloc(m * n, sizeof(double));
        
        double *s33_py = calloc(m * n, sizeof(double));
        double *s23_py = calloc(m * n, sizeof(double));
        double *s13_py = calloc(m * n, sizeof(double));
        
        double *s33_pz = calloc(m * n, sizeof(double));
        double *s23_pz = calloc(m * n, sizeof(double));
        double *s13_pz = calloc(m * n, sizeof(double));
        
        double *lhs = calloc(9 * m * n, sizeof(double));
        double *rhs = calloc(3 * m, sizeof(*rhs));
        
        if (simulated_boundary == NULL) {
                *result = TMK_ALLOCATION_FAILURE;
                goto cleanup;
        }
        if (s33_q == NULL || s23_q == NULL || s13_q == NULL) {
                *result = TMK_ALLOCATION_FAILURE;
                goto cleanup;
        }
        if (s33_px == NULL || s23_px == NULL || s13_px == NULL) {
                *result = TMK_ALLOCATION_FAILURE;
                goto cleanup;
        }
        if (s33_py == NULL || s23_py == NULL || s13_py == NULL) {
                *result = TMK_ALLOCATION_FAILURE;
                goto cleanup;
        }
        if (s33_pz == NULL || s23_pz == NULL || s13_pz == NULL) {
                *result = TMK_ALLOCATION_FAILURE;
                goto cleanup;
        }
        if (lhs == NULL || rhs == NULL) {
                *result = TMK_ALLOCATION_FAILURE;
                goto cleanup;
        }

        tmk__geometry_grid_uniform(simulated_boundary, mside, hm, alpha, load);
        tmk__geometry_grid_uniform(force_boundary, nside, hn, alpha, load);
        tmk__geometry_grid_offset(force_boundary, n, beta, x0);

        for (int j = 0; j < n; ++j) {
                for (int i = 0; i < m; ++i) {
                        int idx = j * m + i;
                        double coeffs[6] = {0};

                        tmk__coefficients_concentrated_stress_x(coeffs, 
                                                                &simulated_boundary[i], 
                                                                &force_boundary[j], 
                                                                halfspace);
                        tmk__coefficients_tensor_transformation_y(coeffs, alpha);

                        s33_px[idx] = coeffs[2];
                        s23_px[idx] = coeffs[3];
                        s13_px[idx] = coeffs[4];

                        tmk__coefficients_concentrated_stress_y(coeffs, 
                                                                &simulated_boundary[i], 
                                                                &force_boundary[j], 
                                                                halfspace);
                        tmk__coefficients_tensor_transformation_y(coeffs, alpha);

                        s33_py[idx] = coeffs[2];
                        s23_py[idx] = coeffs[3];
                        s13_py[idx] = coeffs[4];

                        tmk__coefficients_concentrated_stress_z(coeffs, 
                                                                &simulated_boundary[i], 
                                                                &force_boundary[j], 
                                                                halfspace);
                        tmk__coefficients_tensor_transformation_y(coeffs, alpha);

                        s33_pz[idx] = coeffs[2];
                        s23_pz[idx] = coeffs[3];
                        s13_pz[idx] = coeffs[4];
                }
        }

        tmk__util_fill_lhs(lhs, m, n, 
                           s33_px, s33_py, s33_pz, 
                           s23_px, s23_py, s23_pz, 
                           s13_px, s13_py, s13_pz);

        for (int i = 0; i < m; ++i) {
                double coeffs[6] = {0};

                tmk__coefficients_distributed_stress(coeffs, &simulated_boundary[i], load, halfspace);
                tmk__coefficients_tensor_transformation_y(coeffs, alpha);
                s33_q[i] = coeffs[2];
                s23_q[i] = coeffs[3];
                s13_q[i] = coeffs[4];
        }

        memcpy(rhs, s33_q, m * sizeof(double));
        memcpy(rhs + m, s23_q, m * sizeof(double));
        memcpy(rhs + 2 * m, s13_q, m * sizeof(double));

        cblas_dscal(3 * m, -1.0 * q, rhs, 1);

        printf("Solving...\n");
        *info_lapack = LAPACKE_dgels(LAPACK_COL_MAJOR, 'N', 
                                     3 * m, 3 * n, 1, 
                                     lhs, 3 * m, rhs, 3 * m);

        if (*info_lapack == 0) {
                *result = TMK_SUCCESS;
                for (int i = 0; i < 3 * n; ++i) {
                        sol[i] = rhs[i];
                }
                goto cleanup;
        } else {
                *result = TMK_LAPACK_FAILURE;
                goto cleanup;
        }

        printf("Solved.\n");

cleanup: 
        free(simulated_boundary);
        free(s33_q);
        free(s23_q);
        free(s13_q);
        free(s33_px);
        free(s23_px);
        free(s13_px);
        free(s33_py);
        free(s23_py);
        free(s13_py);
        free(s33_pz);
        free(s23_pz);
        free(s13_pz);
        free(lhs);
        free(rhs);
}

static void tmk__geometry_grid_properties(const double x0, 
                                          double *hm, double *hn, double *beta)
{
    *hm = 2.0 * 0.08 * x0;
    *hn = 2.0 * 0.1 * x0;
    *beta = 1.25 * *hn / x0;
}

static void tmk__evaluate_single(tmk_vec3 *pts, const int npts, 
                                 tmk_load *load, tmk_hs *halfspace,
                                 double *sol, tmk_vec3 *force_boundary, const int n_forces)
{
        printf("Evaluating at given points...\n");

        for (int p = 0; p < npts; ++p) {
                tmk__util_correct_singularity(&pts[p], load);
                
                double sxx = 0, syy = 0, szz = 0, syz = 0, sxz = 0, sxy = 0;
                double ux = 0, uy = 0, uz = 0;
                double cs[6] = {0};
                double cu[3] = {0};

                // by load
                tmk__coefficients_distributed_stress(cs, 
                                                     &pts[p], 
                                                     load, 
                                                     halfspace);
                sxx += cs[0] * load->intensity;
                syy += cs[1] * load->intensity;
                szz += cs[2] * load->intensity;
                syz += cs[3] * load->intensity;
                sxz += cs[4] * load->intensity;
                sxy += cs[5] * load->intensity;

                tmk__coefficients_distributed_displacement(cu, 
                                                           &pts[p], 
                                                           load, 
                                                           halfspace);
                ux += cu[0] * load->intensity;
                uy += cu[1] * load->intensity;
                uz += cu[2] * load->intensity;

                // by force boundary
                for (int j = 0; j < n_forces; ++j) {
                        // px
                        tmk__coefficients_concentrated_stress_x(cs, 
                                                                &pts[p], 
                                                                &force_boundary[j], 
                                                                halfspace);
                        sxx += cs[0] * sol[0 * n_forces + j];
                        syy += cs[1] * sol[0 * n_forces + j];
                        szz += cs[2] * sol[0 * n_forces + j];
                        syz += cs[3] * sol[0 * n_forces + j];
                        sxz += cs[4] * sol[0 * n_forces + j];
                        sxy += cs[5] * sol[0 * n_forces + j];
                        tmk__coefficients_concentrated_displacement_x(cu, 
                                                                      &pts[p], 
                                                                      &force_boundary[j], 
                                                                      halfspace);
                        ux += cu[0] * sol[0 * n_forces + j];
                        uy += cu[1] * sol[0 * n_forces + j];
                        uz += cu[2] * sol[0 * n_forces + j];

                        // py
                        tmk__coefficients_concentrated_stress_y(cs, 
                                                                &pts[p], 
                                                                &force_boundary[j], 
                                                                halfspace);
                        sxx += cs[0] * sol[1 * n_forces + j];
                        syy += cs[1] * sol[1 * n_forces + j];
                        szz += cs[2] * sol[1 * n_forces + j];
                        syz += cs[3] * sol[1 * n_forces + j];
                        sxz += cs[4] * sol[1 * n_forces + j];
                        sxy += cs[5] * sol[1 * n_forces + j];
                        tmk__coefficients_concentrated_displacement_y(cu, 
                                                                      &pts[p], 
                                                                      &force_boundary[j], 
                                                                      halfspace);
                        ux += cu[0] * sol[1 * n_forces + j];
                        uy += cu[1] * sol[1 * n_forces + j];
                        uz += cu[2] * sol[1 * n_forces + j];

                        // pz
                        tmk__coefficients_concentrated_stress_z(cs, 
                                                                &pts[p], 
                                                                &force_boundary[j], 
                                                                halfspace);
                        sxx += cs[0] * sol[2 * n_forces + j];
                        syy += cs[1] * sol[2 * n_forces + j];
                        szz += cs[2] * sol[2 * n_forces + j];
                        syz += cs[3] * sol[2 * n_forces + j];
                        sxz += cs[4] * sol[2 * n_forces + j];
                        sxy += cs[5] * sol[2 * n_forces + j];

                        tmk__coefficients_concentrated_displacement_z(cu, 
                                                                      &pts[p], 
                                                                      &force_boundary[j], 
                                                                      halfspace);
                        ux += cu[0] * sol[2 * n_forces + j];
                        uy += cu[1] * sol[2 * n_forces + j];
                        uz += cu[2] * sol[2 * n_forces + j];
                }

                printf("Point %d: (%.0lf, %.0lf, %.0lf)\n", p + 1, pts[p].x, pts[p].y, pts[p].z);
                printf("Stresses:\n");
                printf("  xx = %.3e [MPa]\n", sxx);
                printf("  yy = %.3e [MPa]\n", syy);
                printf("  zz = %.3e [MPa]\n", szz);
                printf("  yz = %.3e [MPa]\n", syz);
                printf("  xz = %.3e [MPa]\n", sxz);
                printf("  xy = %.3e [MPa]\n\n", sxy);
                printf("Displacement:\n");
                printf("  x = %.3e [mm]\n", ux);
                printf("  y = %.3e [mm]\n", uy);
                printf("  z = %.3e [mm]\n\n", uz);

                // val[0 + p * 9] = sxx;
                // val[1 + p * 9] = syy;
                // val[2 + p * 9] = szz;
                // val[3 + p * 9] = syz;
                // val[4 + p * 9] = sxz;
                // val[5 + p * 9] = sxy;
                // val[6 + p * 9] = ux;
                // val[7 + p * 9] = uy;
                // val[8 + p * 9] = uz;
        }

        // mdl->vals_eval = (double *)val;
        printf("Evaluated.\n");
}

#endif // TAAMAK_IMPLEMENTATION

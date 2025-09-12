/*
        taamak.h - v1.0.0-rc - https://github.com/veslrs/taamak.h

        This library provides functionalities for modelling linear elastic asphalt
        pavement system with an edge, subject to a distribited load normal to the 
        pavement surface. 

        # Quick Example
        

        ```c
        TODO: example
        ```

        TODO: API
*/

#ifndef TAAMAK_H
#define TAAMAK_H

#include <stdio.h>

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

typedef struct {double x, y, z;} tmk_vec3;
typedef struct {tmk_vec3 center; double half_width; double intensity;} tmk_load;
typedef struct {double youngs_modulus, poissons_ratio;} tmk_hs;

typedef struct tmk_model_s tmk_model;

TMKDEF void tmk_init(tmk_model *mdl, const int n_layers);
TMKDEF void tmk_set_load(tmk_model *mdl, const tmk_vec3 center, const double half_width, const double intensity);
TMKDEF void tmk_set_composition(tmk_model *mdl, const double *top_depths);
TMKDEF void tmk_set_material_properties(tmk_model *mdl, const tmk_hs *hs);
TMKDEF void tmk_set_evaluation_points(tmk_model *mdl, const int npts, tmk_vec3 *pts);
TMKDEF void tmk_set_slope_angle(tmk_model *mdl, const double degree);
TMKDEF void tmk_set_measurement(tmk_model *mdl, const int npts, tmk_vec3 *pts, double *vals);
TMKDEF void tmk_solve(tmk_model *mdl);
TMKDEF int  tmk_checkout(tmk_model *mdl);

TMKDEF void tmk_linspace(tmk_vec3 *pts, const tmk_vec3 *a, const tmk_vec3 *b, const int npts);

#endif // TAAMAK_H

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
#include <nlopt.h>

#include "cblas.h"
#include "lapacke.h"

#define TMK__VERSION "version 1.0.0-rc"

typedef enum {
        TMK_INITIAL = 0,
        TMK_SUCCESS,
        TMK_ALLOCATION_FAILURE,
        TMK_LAPACK_FAILURE,
        TMK_NO_LOAD,
        TMK_NO_COMPOSITION,
        TMK_NO_MATERIAL,
        TMK_ESA_EXCEED_ITER_LIM,
} tmk_result;

struct tmk_model_s {
        int         n_layers;
        bool        info_comp;
        double     *top_depths;
        bool        info_mat;
        tmk_hs     *halfspace;
        bool        info_load;
        tmk_load    load;
        bool        info_eval_pts;
        int         npts_eval;
        tmk_vec3   *pts_eval;
        bool        info_measure;
        int         npts_measure;
        tmk_vec3   *pts_measure;
        double     *vals_measure;
        bool        info_angle;
        double      slope_angle;
        int         info_lapack;
        tmk_result  result;
        // double     *vals_eval;
        // bool        info_measure_edge;
        // int         npts_measure_edge;
        // tmk_vec3   *pts_measure_edge;
        // double     *vals_measure_edge;
        // bool        info_measure_center;
        // int         npts_measure_center;
        // tmk_vec3   *pts_measure_center;
        // double     *vals_measure_center;
};

// Re-definable macros
#ifndef TMK_GEOM_M
#define TMK_GEOM_M 50
#endif

#ifndef TMK_GEOM_N
#define TMK_GEOM_N 40
#endif

#ifndef TMK_ELT_FACTOR
#define TMK_ELT_FACTOR 0.9
#endif

#ifndef TMK_ELT_DEG_LIM
#define TMK_ELT_DEG_LIM 10.0
#endif

#ifndef TMK_ESA_TOL
#define TMK_ESA_TOL 1e-3
#endif

#ifndef TMK_ESA_ITER
#define TMK_ESA_ITER 20
#endif

#ifndef TMK_ESA_DEG_LIM
#define TMK_ESA_DEG_LIM 10.0
#endif

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
                - 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / r2 / (r2 + z + c) / (r2 + z + c) * (3.0 - x * x * (3.0 * r2 + z + c) / r2 / r2 / (r2 + z + c)) 
                + 6.0 * c / pow(r2, 5.0) * (3.0 * c - (3.0 - 2.0 * poissons_ratio) * (z + c) + 5.0 * x * x * z / r2 / r2)
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
                - 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / r2 / (r2 + z + c) / (r2 + z + c) * (1.0 - y * y * (3.0 * r2 + z + c) / r2 / r2 / (r2 + z + c)) 
                + 6.0 * c / pow(r2, 5.0) * (c - (1 - 2.0 * poissons_ratio) * (z + c) + 5.0 * y * y * z / r2 / r2)
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
                + 6.0 * c / pow(r2, 5.0) * (c + (1.0 - 2.0 * poissons_ratio) * (z + c) + 5.0 * z * (z + c) * (z + c) / r2 / r2)
        );
}

static inline double tmk__coefficients_concentrated_syz_x(double x, double y, double z, double c, 
                                                          double r1, double r2, double poissons_ratio)
{
        return x * y / 8.0 / M_PI / (poissons_ratio - 1.0) * (
                3.0 * (c - z) / pow(r1, 5.0) 
                - 3.0 * (3.0 - 4.0 * poissons_ratio) * (z + c) / pow(r2, 5.0) 
                + 6.0 * c / pow(r2, 5.0) * (1.0 - 2.0 * poissons_ratio + 5.0 * z * (z + c) / r2 / r2)
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
                - 6.0 * c / pow(r2, 5.0) * (z * (z + c) - (1.0 - 2.0 * poissons_ratio) * x * x - 5.0 * x * x * z * (z + c) / r2 / r2)
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
                - 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / r2 / (r2 + z + c) / (r2 + z + c) * (1.0 - x * x * (3.0 * r2 + z + c) / r2 / r2 / (r2 + z + c))
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
                - 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / r2 / (r2 + z + c) / (r2 + z + c) * (1.0 - x * x * (3.0 * r2 + z + c) / r2 / r2 / (r2 + z + c)) 
                + 6.0 * c / pow(r2, 5.0) * (c - (1.0 - 2.0 * poissons_ratio) * (z + c) 
                + 5.0 * x * x * z / r2 / r2)
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
                - 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / r2 / (r2 + z + c) / (r2 + z + c) * (3.0 - y * y * (3.0 * r2 + z + c) / r2 / r2 / (r2 + z + c)) 
                + 6.0 * c / pow(r2, 5.0) * (3.0 * c - (3.0 - 2.0 * poissons_ratio) * (z + c) + 5.0 * y * y * z / r2 / r2)
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
                + 6.0 * c / pow(r2, 5.0) * (c + (1.0 - 2.0 * poissons_ratio) * (z + c) 
                + 5.0 * z * (z + c) * (z + c) / r2 / r2)
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
                - 6.0 * c / pow(r2, 5.0) * (z * (z + c) - (1.0 - 2.0 * poissons_ratio) * y * y 
                - 5.0 * y * y * z * (z + c) / r2 / r2)
        );
}

static inline double tmk__coefficients_concentrated_sxz_y(double x, double y, double z, double c, 
                                                          double r1, double r2, double poissons_ratio)
{
        return y * x / 8.0 / M_PI / (poissons_ratio - 1.0) * (
                3.0 * (c - z) / pow(r1, 5.0) 
                - 3.0 * (3.0 - 4.0 * poissons_ratio) * (z + c) / pow(r2, 5.0) 
                + 6.0 * c / pow(r2, 5.0) * (1.0 - 2.0 * poissons_ratio + 5.0 * z * (z + c) / r2 / r2)
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
                - 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / r2 / (r2 + z + c) / (r2 + z + c) * (1.0 - y * y * (3.0 * r2 + z + c) / r2 / r2 / (r2 + z + c))
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
                - (3.0 * x * x * (3.0 - 4.0 * poissons_ratio) * (z - c) - 6.0 * c * (z + c) * ((1.0 - 2.0 * poissons_ratio) * z - 2.0 * poissons_ratio * c)) / pow(r2, 5.0) 
                - 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / r2 / (r2 + z + c) * (1.0 - x * x / r2 / (r2 + z + c) - x * x / r2 / r2)
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
                - (3.0 * (3.0 - 4.0 * poissons_ratio) * y * y * (z - c) - 6.0 * c * (z + c) * ((1.0 - 2.0 * poissons_ratio) * z - 2.0 * poissons_ratio * c)) / pow(r2, 5.0) 
                - 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / r2 / (r2 + z + c) * (1.0 - y * y / r2 / (r2 + z + c) - y * y / r2 / r2)
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
                - (3.0 * (3.0 - 4.0 * poissons_ratio) * z * (z + c) * (z + c) - 3.0 * c * (z + c) * (5.0 * z - c)) / pow(r2, 5.0)
        );
}

static inline double tmk__coefficients_concentrated_syz_z(double y, double z, double c, 
                                                          double r1, double r2, double poissons_ratio)
{
        return y / 8.0 / M_PI / (poissons_ratio - 1.0) * (
                (2.0 * poissons_ratio - 1.0) / pow(r1, 3.0) 
                + (1.0 - 2.0 * poissons_ratio) / pow(r2, 3.0) 
                - (3.0 * (3.0 - 4.0 * poissons_ratio) * z * (z + c) - 3.0 * c * (3.0 * z + c)) / pow(r2, 5.0) 
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
                - (3.0 * (3.0 - 4.0 * poissons_ratio) * z * (z + c) - 3.0 * c * (3.0 * z + c)) / pow(r2, 5.0) 
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
                + 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / r2 / r2 / (r2 + z + c) * (1.0 / (r2 + z + c) + 1.0 / r2)
        );
}

static inline double tmk__coefficients_concentrated_ux_x(double x, double z, double c, 
                                                         double r1, double r2, double youngs_modulus, double poissons_ratio)
{
        return (1.0 + poissons_ratio) / 8.0 / M_PI / youngs_modulus / (1.0 - poissons_ratio) * (
                (3.0 - 4.0 * poissons_ratio) / r1 + 1.0 / r2 
                + (3.0 - 4.0 * poissons_ratio) * x * x / pow(r2, 3.0) 
                + 2.0 * c * z / pow(r2, 3.0) * (1.0 - 3.0 * x * x / r2 / r2) 
                + x * x / pow(r1, 3.0) 
                + 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / (r2 + z + c) * (1.0 - x * x / r2 / (r2 + z + c))
        );
}

static inline double tmk__coefficients_concentrated_uy_x(double x, double y, double z, double c, 
                                                         double r1, double r2, double youngs_modulus, double poissons_ratio)
{
        return (1.0 + poissons_ratio) * x * y / 8.0 / M_PI / youngs_modulus / (1.0 - poissons_ratio) * (
                1.0 / pow(r1, 3.0) 
                + (3.0 - 4.0 * poissons_ratio) / pow(r2, 3.0) 
                - 6.0 * c * z / pow(r2, 5.0) 
                - 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / r2 / (r2 + z + c) / (r2 + z + c)
        );
}

static inline double tmk__coefficients_concentrated_uz_x(double x, double z, double c, 
                                                         double r1, double r2, double youngs_modulus, double poissons_ratio)
{
        return (1.0 + poissons_ratio) * x / 8.0 / M_PI / youngs_modulus / (1.0 - poissons_ratio) * (
                (z - c) / pow(r1, 3.0) 
                + (3.0 - 4.0 * poissons_ratio) * (z - c) / pow(r2, 3.0) 
                - 6.0 * c * z * (z + c) / pow(r2, 5.0) 
                + 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / r2 / (r2 + z + c)
        );
}

static inline double tmk__coefficients_concentrated_ux_y(double x, double y, double z, double c, 
                                                         double r1, double r2, double youngs_modulus, double poissons_ratio)
{
        return (1.0 + poissons_ratio) * x * y / 8.0 / M_PI / youngs_modulus / (1.0 - poissons_ratio) * (
                1.0 / pow(r1, 3.0) 
                + (3.0 - 4.0 * poissons_ratio) / pow(r2, 3.0) 
                - 6.0 * c * z / pow(r2, 5.0) 
                - 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / r2 / (r2 + z + c) / (r2 + z + c)
        );
}

static inline double tmk__coefficients_concentrated_uy_y(double y, double z, double c, 
                                                         double r1, double r2, double youngs_modulus, double poissons_ratio)
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
                                                         double r1, double r2, double youngs_modulus, double poissons_ratio)
{
        return (1.0 + poissons_ratio) * y / 8.0 / M_PI / youngs_modulus / (1.0 - poissons_ratio) * (
                (z - c) / pow(r1, 3.0) 
                + (3.0 - 4.0 * poissons_ratio) * (z - c) / pow(r2, 3.0) 
                - 6.0 * c * z * (z + c) / pow(r2, 5.0) 
                + 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / r2 / (r2 + z + c)
        );
}

static inline double tmk__coefficients_concentrated_ux_z(double x, double z, double c, 
                                                         double r1, double r2, double youngs_modulus, double poissons_ratio)
{
        return (1.0 + poissons_ratio) * x / 8.0 / M_PI / youngs_modulus / (1.0 - poissons_ratio) * (
                (z - c) / pow(r1, 3) 
                + (3.0 - 4.0 * poissons_ratio) * (z - c) / pow(r2, 3.0) 
                - 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / r2 / (r2 + z + c) 
                + 6.0 * c * z * (z + c) / pow(r2, 5.0)
        );
}

static inline double tmk__coefficients_concentrated_uy_z(double y, double z, double c, 
                                                         double r1, double r2, double youngs_modulus, double poissons_ratio)
{
        return (1.0 + poissons_ratio) * y / 8.0 / M_PI / youngs_modulus / (1.0 - poissons_ratio) * (
                (z - c) / pow(r1, 3.0) 
                + (3.0 - 4.0 * poissons_ratio) * (z - c) / pow(r2, 3.0) 
                - 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / r2 / (r2 + z + c) 
                + 6.0 * c * z * (z + c) / pow(r2, 5.0)
        );
}

static inline double tmk__coefficients_concentrated_uz_z(double z, double c, 
                                                         double r1, double r2, double youngs_modulus, double poissons_ratio)
{
        return (1.0 + poissons_ratio) / 8.0 / M_PI / youngs_modulus / (1.0 - poissons_ratio) * (
                (3.0 - 4.0 * poissons_ratio) / r1 
                + (8.0 * (1.0 - poissons_ratio) * (1.0 - poissons_ratio) - (3.0 - 4.0 * poissons_ratio)) / r2 
                + (z - c) * (z - c) / pow(r1, 3.0) 
                + ((3.0 - 4.0 * poissons_ratio) * (z + c) * (z + c) - 2.0 * c * z) / pow(r2, 3.0) 
                + 6.0 * c * z * (z + c) * (z + c) / pow(r2, 5.0)
        );
}

static void tmk__coefficients_concentrated_stress_x(double *coefficients, 
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

        double poissons_ratio = halfspace->poissons_ratio;

        // xx
        coefficients[0] = tmk__coefficients_concentrated_sxx_x(x, z, c, r1, r2, poissons_ratio);
        // yy
        coefficients[1] = tmk__coefficients_concentrated_syy_x(x, y, z, c, r1, r2, poissons_ratio);
        // zz
        coefficients[2] = tmk__coefficients_concentrated_szz_x(x, z, c, r1, r2, poissons_ratio);
        // yz
        coefficients[3] = tmk__coefficients_concentrated_syz_x(x, y, z, c, r1, r2, poissons_ratio);
        // xz
        coefficients[4] = tmk__coefficients_concentrated_sxz_x(x, z, c, r1, r2, poissons_ratio);
        // xy
        coefficients[5] = tmk__coefficients_concentrated_sxy_x(x, y, z, c, r1, r2, poissons_ratio);
}

static void tmk__coefficients_concentrated_stress_y(double *coefficients, 
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

        double poissons_ratio = halfspace->poissons_ratio;

        // xx
        coefficients[0] = tmk__coefficients_concentrated_sxx_y(x, y, z, c, r1, r2, poissons_ratio);
        // yy
        coefficients[1] = tmk__coefficients_concentrated_syy_y(y, z, c, r1, r2, poissons_ratio);
        // zz
        coefficients[2] = tmk__coefficients_concentrated_szz_y(y, z, c, r1, r2, poissons_ratio);
        // yz
        coefficients[3] = tmk__coefficients_concentrated_syz_y(y, z, c, r1, r2, poissons_ratio);
        // xz
        coefficients[4] = tmk__coefficients_concentrated_sxz_y(x, y, z, c, r1, r2, poissons_ratio);
        // xy
        coefficients[5] = tmk__coefficients_concentrated_sxy_y(x, y, z, c, r1, r2, poissons_ratio);
}

static void tmk__coefficients_concentrated_stress_z(double *coefficients, 
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

        double poissons_ratio = halfspace->poissons_ratio;

        // xx
        coefficients[0] = tmk__coefficients_concentrated_sxx_z(x, z, c, r1, r2, poissons_ratio);
        // yy
        coefficients[1] = tmk__coefficients_concentrated_syy_z(y, z, c, r1, r2, poissons_ratio);
        // zz
        coefficients[2] = tmk__coefficients_concentrated_szz_z(z, c, r1, r2, poissons_ratio);
        // yz
        coefficients[3] = tmk__coefficients_concentrated_syz_z(y, z, c, r1, r2, poissons_ratio);
        // xz
        coefficients[4] = tmk__coefficients_concentrated_sxz_z(x, z, c, r1, r2, poissons_ratio);
        // xy
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

        // x
        coefficients[0] = tmk__coefficients_concentrated_ux_x(x, z, c, r1, r2, youngs_modulus, poissons_ratio);
        // y
        coefficients[1] = tmk__coefficients_concentrated_uy_x(x, y, z, c, r1, r2, youngs_modulus, poissons_ratio);
        // z
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

        // x
        coefficients[0] = tmk__coefficients_concentrated_ux_y(x, y, z, c, r1, r2, youngs_modulus, poissons_ratio);
        // y
        coefficients[1] = tmk__coefficients_concentrated_uy_y(y, z, c, r1, r2, youngs_modulus, poissons_ratio);
        // z
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

        // x
        coefficients[0] = tmk__coefficients_concentrated_ux_z(x, z, c, r1, r2, youngs_modulus, poissons_ratio);
        // y
        coefficients[1] = tmk__coefficients_concentrated_uy_z(y, z, c, r1, r2, youngs_modulus, poissons_ratio);
        // z
        coefficients[2] = tmk__coefficients_concentrated_uz_z(z, c, r1, r2, youngs_modulus, poissons_ratio);
}

static inline double tmk__coefficients_distributed_a_xx(double x, double y, double z, double poissons_ratio)
{
        return 2.0 * poissons_ratio * (atan2(x, y) - atan2(z * x, y * tmk__coefficients_radius(x, y, z))) 
               + atan2(y, x) 
               - atan2(z * y, x * tmk__coefficients_radius(x, y, z)) 
               - x * y * z / (tmk__coefficients_radius(x, y, z) * (x * x + z * z));
}

static inline double tmk__coefficients_distributed_a_yy(double x, double y, double z, double poissons_ratio)
{
        return 2.0 * poissons_ratio * (atan2(y, x) - atan2(z * y, x * tmk__coefficients_radius(x, y, z))) 
               + atan2(x, y) 
               - atan2(z * x, y * tmk__coefficients_radius(x, y, z)) 
               - x * y * z / (tmk__coefficients_radius(x, y, z) * (y * y + z * z));
}

static inline double tmk__coefficients_distributed_a_zz(double x, double y, double z)
{
        return atan2(y, x) 
               + atan2(x, y) 
               - atan2(z * y, x * tmk__coefficients_radius(x, y, z)) 
               - atan2(z * x, y * tmk__coefficients_radius(x, y, z)) 
               + x * y * z / (tmk__coefficients_radius(x, y, z) * (y * y + z * z)) 
               + x * y * z / (tmk__coefficients_radius(x, y, z) * (x * x + z * z));
}

static inline double tmk__coefficients_distributed_a_yz(double x, double y, double z)
{
        return -1.0 * z * z * x / tmk__coefficients_radius(x, y, z) / (y * y + z * z);
}

static inline double tmk__coefficients_distributed_a_xz(double x, double y, double z)
{
        return -1.0 * z * z * y / tmk__coefficients_radius(x, y, z) / (x * x + z * z);
}

static inline double tmk__coefficients_distributed_a_xy(double x, double y, double z, double poissons_ratio)
{
        return (1.0 - 2.0 * poissons_ratio) * log(tmk__coefficients_radius(x, y, z) + z) 
               + z / tmk__coefficients_radius(x, y, z);
}

static inline double tmk__coefficients_distributed_b_x(double x, double y, double z, double poissons_ratio)
{
        return (2.0 * poissons_ratio - 1.0) * (
                        y * log(tmk__coefficients_radius(x, y, z) + z) 
                        + z * log(tmk__coefficients_radius(x, y, z) + y) 
                        - 2.0 * x * atan2(x, tmk__coefficients_radius(x, y, z) + y + z)
                ) - z * log(tmk__coefficients_radius(x, y, z) + y);
}

static inline double tmk__coefficients_distributed_b_y(double x, double y, double z, double poissons_ratio)
{
        return (2.0 * poissons_ratio - 1.0) * (
                        x * log(tmk__coefficients_radius(x, y, z) + z) 
                        + z * log(tmk__coefficients_radius(x, y, z) + x) 
                        - 2.0 * y * atan2(y, tmk__coefficients_radius(x, y, z) + x + z)
                ) - z * log(tmk__coefficients_radius(x, y, z) + x);
}

static inline double tmk__coefficients_distributed_b_z(double x, double y, double z, double poissons_ratio)
{
        return 2.0 * (1.0 - poissons_ratio) * (
                        y * log(tmk__coefficients_radius(x, y, z) + x) 
                        + x * log(tmk__coefficients_radius(x, y, z) + y) 
                        + 2.0 * z * atan2((sqrt(x * x + z * z) - x) * (tmk__coefficients_radius(x, y, z) 
                        - sqrt(x * x + z * z)), z * y)
                ) + z * atan2(x * y, tmk__coefficients_radius(x, y, z) * z);
}

static void tmk__coefficients_distributed_stress(double *coefficients, tmk_vec3 *evaluation_point, tmk_load *load, tmk_hs *halfspace)
{
        double x = evaluation_point->x - load->center.x;
        double y = evaluation_point->y - load->center.y;
        double z = evaluation_point->z - load->center.z;

        double b = load->half_width;

        double poissons_ratio = halfspace->poissons_ratio;

        // xx
        coefficients[0] = (tmk__coefficients_distributed_a_xx(x + b, y + b, z, poissons_ratio) + tmk__coefficients_distributed_a_xx(x - b, y - b, z, poissons_ratio) - tmk__coefficients_distributed_a_xx(x - b, y + b, z, poissons_ratio) - tmk__coefficients_distributed_a_xx(x + b, y - b, z, poissons_ratio)) / 2.0 / M_PI;
        // yy
        coefficients[1] = (tmk__coefficients_distributed_a_yy(x + b, y + b, z, poissons_ratio) + tmk__coefficients_distributed_a_yy(x - b, y - b, z, poissons_ratio) - tmk__coefficients_distributed_a_yy(x - b, y + b, z, poissons_ratio) - tmk__coefficients_distributed_a_yy(x + b, y - b, z, poissons_ratio)) / 2.0 / M_PI;
        // zz
        coefficients[2] = (tmk__coefficients_distributed_a_zz(x + b, y + b, z) + tmk__coefficients_distributed_a_zz(x - b, y - b, z) - tmk__coefficients_distributed_a_zz(x - b, y + b, z) - tmk__coefficients_distributed_a_zz(x + b, y - b, z)) / 2.0 / M_PI;
        // yz
        coefficients[3] = (tmk__coefficients_distributed_a_yz(x + b, y + b, z) + tmk__coefficients_distributed_a_yz(x - b, y - b, z) - tmk__coefficients_distributed_a_yz(x - b, y + b, z) - tmk__coefficients_distributed_a_yz(x + b, y - b, z)) / 2.0 / M_PI;
        // xz
        coefficients[4] = (tmk__coefficients_distributed_a_xz(x + b, y + b, z) + tmk__coefficients_distributed_a_xz(x - b, y - b, z) - tmk__coefficients_distributed_a_xz(x - b, y + b, z) - tmk__coefficients_distributed_a_xz(x + b, y - b, z)) / 2.0 / M_PI;
        // xy
        coefficients[5] = (tmk__coefficients_distributed_a_xy(x + b, y + b, z, poissons_ratio) + tmk__coefficients_distributed_a_xy(x - b, y - b, z, poissons_ratio) - tmk__coefficients_distributed_a_xy(x - b, y + b, z, poissons_ratio) - tmk__coefficients_distributed_a_xy(x + b, y - b, z, poissons_ratio)) / 2.0 / M_PI;
}

static void tmk__coefficients_distributed_displacement(double *coefficients, tmk_vec3 *evaluation_point, tmk_load *load, tmk_hs *halfspace)
{
        double x = evaluation_point->x - load->center.x;
        double y = evaluation_point->y - load->center.y;
        double z = evaluation_point->z - load->center.z;

        double b = load->half_width;

        double youngs_modulus = halfspace->youngs_modulus;
        double poissons_ratio = halfspace->poissons_ratio;

        // x
        coefficients[0] = (1.0 + poissons_ratio) / 2.0 / M_PI / youngs_modulus * (tmk__coefficients_distributed_b_x(x + b, y + b, z, poissons_ratio) + tmk__coefficients_distributed_b_x(x - b, y - b, z, poissons_ratio) - tmk__coefficients_distributed_b_x(x - b, y + b, z, poissons_ratio) - tmk__coefficients_distributed_b_x(x + b, y - b, z, poissons_ratio));
        // y
        coefficients[1] = (1.0 + poissons_ratio) / 2.0 / M_PI / youngs_modulus * (tmk__coefficients_distributed_b_y(x + b, y + b, z, poissons_ratio) + tmk__coefficients_distributed_b_y(x - b, y - b, z, poissons_ratio) - tmk__coefficients_distributed_b_y(x - b, y + b, z, poissons_ratio) - tmk__coefficients_distributed_b_y(x + b, y - b, z, poissons_ratio));
        // z
        coefficients[2] = (1.0 + poissons_ratio) / 2.0 / M_PI / youngs_modulus * (tmk__coefficients_distributed_b_z(x + b, y + b, z, poissons_ratio) + tmk__coefficients_distributed_b_z(x - b, y - b, z, poissons_ratio) - tmk__coefficients_distributed_b_z(x - b, y + b, z, poissons_ratio) - tmk__coefficients_distributed_b_z(x + b, y - b, z, poissons_ratio));
}

static double tmk__coefficients_distributed_uz(tmk_vec3 *evaluation_point, tmk_load *load, tmk_hs *halfspace)
{
        double x = evaluation_point->x - load->center.x;
        double y = evaluation_point->y - load->center.y;
        double z = evaluation_point->z - load->center.z;

        double b = load->half_width;

        double youngs_modulus = halfspace->youngs_modulus;
        double poissons_ratio = halfspace->poissons_ratio;

        double uz = (1.0 + poissons_ratio) / 2.0 / M_PI / youngs_modulus * (tmk__coefficients_distributed_b_z(x + b, y + b, z, poissons_ratio) + tmk__coefficients_distributed_b_z(x - b, y - b, z, poissons_ratio) - tmk__coefficients_distributed_b_z(x - b, y + b, z, poissons_ratio) - tmk__coefficients_distributed_b_z(x + b, y - b, z, poissons_ratio));

        return uz;
}

static void tmk__coefficients_tensor_transformation_y(double *tensor, double degree)
{
        double radian = degree * M_PI / 180.0;

        double xx = tensor[0];
        double zz = tensor[2];
        double yz = tensor[3];
        double xz = tensor[4];
        double xy = tensor[5];

        tensor[2] = zz * cos(radian) * cos(radian) + xx * sin(radian) * sin(radian) + xz * sin(2.0 * radian);
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

//         tensor[0] = xx * cos(radian) * cos(radian) - 2.0 * xz * cos(radian) * sin(radian) + zz * sin(radian) * sin(radian);
//         tensor[2] = zz * cos(radian) * cos(radian) + xx * sin(radian) * sin(radian) + xz * sin(2.0 * radian);
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

//         tensor[0] = xx * cos(radian) * cos(radian) + yy * sin(radian) * sin(radian) + xy * sin(2.0 * radian);
//         tensor[1] = yy * cos(radian) * cos(radian) - 2.0 * xy * cos(radian) * sin(radian) + xx * sin(radian) * sin(radian);
//         tensor[3] = yz * cos(radian) - xz * sin(radian);
//         tensor[4] = xz * cos(radian) + yz * sin(radian);
//         tensor[5] = xy * cos(2.0 * radian) + (yy - xx) * cos(radian) * sin(radian);
// }

TMKDEF void tmk_init(tmk_model *mdl, const int n_layers)
{
        mdl->n_layers = n_layers;
        mdl->info_comp = false;
        mdl->info_load = false;
        mdl->info_mat = false;
        mdl->info_eval_pts = false;
        mdl->info_measure = false;
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

TMKDEF void tmk_set_composition(tmk_model *mdl, const double *top_depths)
{
        mdl->top_depths = (double *)top_depths;
        mdl->info_comp = true;
}

TMKDEF void tmk_set_material_properties(tmk_model *mdl,
                                        const tmk_hs *hs)
{
        mdl->halfspace = (tmk_hs *)hs;
        mdl->info_mat = true;
}

#define TMK_GEOM_EPSILON 0.001

static void tmk__util_correct_singularity(tmk_vec3 *point, const tmk_load *load)
{
        if (fabs(fabs(point->x - load->center.x) - load->half_width) < 1.0e-9)
                point->x += TMK_GEOM_EPSILON;
        
        if (fabs(fabs(point->y - load->center.y) - load->half_width) < 1.0e-9)
                point->y += TMK_GEOM_EPSILON;
}

TMKDEF void tmk_linspace(tmk_vec3 *pts, const tmk_vec3 *a, const tmk_vec3 *b, const int npts)
{
        double dx = b->x - a->x;
        double dy = b->y - a->y;
        double dz = b->z - a->z;

        double n_step = (double)npts - 1.0;

        double hx = dx / n_step;
        double hy = dy / n_step;
        double hz = dz / n_step;

        for (int i = 0; i < npts; ++i) {
                pts[i].x = (a->x + (double)i * hx);
                pts[i].y = (a->y + (double)i * hy);
                pts[i].z = (a->z + (double)i * hz);
        }
}

TMKDEF void tmk_set_evaluation_points(tmk_model *mdl, const int npts, tmk_vec3 *pts)
{
        mdl->npts_eval = npts;
        for (int i = 0; i < npts; ++i) {
                tmk__util_correct_singularity(&pts[i], &mdl->load);
        }
        mdl->pts_eval = (tmk_vec3 *)pts;
        mdl->info_eval_pts = true;
}

TMKDEF void tmk_set_slope_angle(tmk_model *mdl, const double degree)
{
        mdl->slope_angle = degree;
        mdl->info_angle = true;
}

TMKDEF void tmk_set_measurement(tmk_model *mdl, const int npts, tmk_vec3 *pts, double *vals)
{
        mdl->npts_measure = npts;
        for (int i = 0; i < npts; ++i) {
                tmk__util_correct_singularity(&pts[i], &mdl->load);
        }
        mdl->pts_measure = (tmk_vec3 *)pts;
        mdl->vals_measure = (double *)vals;
        mdl->info_measure = true;
}

static void tmk__evaluate_elt(tmk_model *mdl);
static void tmk__solve_esa(tmk_model *mdl);

static void tmk__geometry_grid_uniform(tmk_vec3 *grid, size_t n_per_side, double h, double degree, tmk_load *load);
static void tmk__geometry_grid_offset(tmk_vec3 *grid, size_t npts, double beta, double x0);

static void tmk__memory_lhs_fill(double *lhs, size_t m, size_t n, double *s33_px, double *s33_py, double *s33_pz, double *s23_px, double *s23_py, double *s23_pz, double *s13_px, double *s13_py, double *s13_pz);

#define TMK__FMT_INDENT "       "

static tmk_result tmk__solve_single(double *sol, tmk_hs *halfspace, tmk_load *load, 
                                    double degree, size_t mside, size_t nside, tmk_vec3 *force_boundary,
                                    tmk_result *result, int *info_lapack);

TMKDEF void tmk_solve(tmk_model *mdl)
{
        printf("[INFO] taamak.h %s\n", (const char *)TMK__VERSION);
        if (mdl->info_load == false) {
                mdl->result = TMK_NO_LOAD;
                return;
        }
        printf("[INFO] Problem Statement:\n");
        if (mdl->n_layers > 1) {
                if (mdl->info_comp == false) {
                        mdl->result = TMK_NO_COMPOSITION;
                        return;
                }
                if (mdl->info_angle == true) {
                        printf(TMK__FMT_INDENT "Evaluation of vertical deflection on the surface of a %d-layer\n" \
                               TMK__FMT_INDENT "linear elastic pavement system. The pavement structure is subject to\n" \
                               TMK__FMT_INDENT "load normal to the surface, distributed in a square patch. The load is\n" \
                               TMK__FMT_INDENT "considered close to the edge of the pavement.\n", mdl->n_layers);
                        // TODO: fancy loggings
                        // printf("[INFO] Model Detail:\n" \
                        //        TMK__FMT_INDENT "Layer compositions and material properties:\n" \
                        //        TMK__FMT_INDENT "    Layer #\tYoung's modulus [MPa]\tPoisson's ratio [-]\tThickness[mm]\n");
                        tmk__evaluate_elt(mdl);
                } else {
                        printf(TMK__FMT_INDENT "Solve effective slope angle for a %d-layer linear elastic pavement system.\n" \
                               TMK__FMT_INDENT "The pavement structure is subject to load normal to the surface,\n" \
                               TMK__FMT_INDENT "distributed in a square patch. The load is considered close to the edge\n" \
                               TMK__FMT_INDENT "of the pavement.\n", mdl->n_layers);
                        // TODO: fancy loggings
                        // printf("[INFO] Model Detail:\n" \
                        //        TMK__FMT_INDENT "Layer compositions and material properties:\n" \
                        //        TMK__FMT_INDENT "    Layer #\tYoung's modulus [MPa]\tPoisson's ratio [-]\tThickness[mm]\n");
                        tmk__solve_esa(mdl);
                }
        } else {
                printf(TMK__FMT_INDENT "Evaluation of vertical deflection on the surface of a homogeneous\n" \
                       TMK__FMT_INDENT "linear elastic pavement system. The pavement structure is subject to\n" \
                       TMK__FMT_INDENT "load normal to the surface, distributed in a square patch. The load is\n" \
                       TMK__FMT_INDENT "considered close to the edge of the pavement.\n");
                // TODO: fancy loggings
                // printf("[INFO] Model Detail:\n" \
                //        TMK__FMT_INDENT "Layer compositions and material properties:\n" \
                //        TMK__FMT_INDENT "    Layer #\tYoung's modulus [MPa]\tPoisson's ratio [-]\tThickness[mm]\n");
                size_t n = TMK_GEOM_N * TMK_GEOM_N;
                double *sol = calloc(3 * n, sizeof(double));
                tmk_vec3 *force_boundary = calloc(n, sizeof(tmk_vec3));

                if (sol == NULL || force_boundary == NULL) {
                        mdl->result = TMK_ALLOCATION_FAILURE;
                        free(sol);
                        free(force_boundary);
                }

                tmk__solve_single(sol, &mdl->halfspace[0], &mdl->load, 
                                  mdl->slope_angle, TMK_GEOM_M, TMK_GEOM_N, 
                                  force_boundary, &mdl->result, &mdl->info_lapack);

                free(sol);
                free(force_boundary);
        }
}

static double tmk__elt_equivalent_thickness(double current_thickness, tmk_hs *current, tmk_hs *target);
static void tmk__elt_reset_z(tmk_vec3 *a, size_t npts);
static tmk_result tmk__elt_intermediate_evaluation(tmk_model *mdl, double *uzs, double *uzs_0, 
                                                   double *heqs, double *thicknesses, 
                                                   double *sol, tmk_vec3 *force_boundary,
                                                   size_t mside, size_t nside);

static void tmk__evaluate_elt(tmk_model *mdl)
{
        printf("[INFO] Initialising solver...\n");
        int n_finite = mdl->n_layers - 1;
        size_t mside = (size_t)TMK_GEOM_M;
        size_t nside = (size_t)TMK_GEOM_N;
        size_t n = nside * nside;

        double *uzs_0 = calloc(mdl->npts_eval, sizeof(double));
        double *thicknesses = calloc((size_t)(n_finite), sizeof(double));
        double *heqs = calloc(mdl->n_layers, sizeof(double));
        double *uzs = calloc(mdl->npts_eval, sizeof(double));
        double *sol = calloc(3 * n, sizeof(double));
        tmk_vec3 *force_boundary = calloc(n, sizeof(tmk_vec3));

        if (uzs_0 == NULL || thicknesses == NULL || heqs == NULL || uzs == NULL) {
                mdl->result = TMK_ALLOCATION_FAILURE;
                goto cleanup;
        }

        for (size_t i = 0; i < (size_t)n_finite; ++i) {
                thicknesses[i] = mdl->top_depths[i + 1] - mdl->top_depths[i];
        }

        for (int i = n_finite; i >= 0; --i) {
                for (int j = 0; j < i; ++j) {
                        heqs[i] += tmk__elt_equivalent_thickness(thicknesses[j], &mdl->halfspace[j], &mdl->halfspace[i]);
                }

                if (i == n_finite) {
                        for (int p = 0; p < mdl->npts_eval; ++p) {
                                mdl->pts_eval[p].z = heqs[i];
                                uzs_0[p] += tmk__coefficients_distributed_uz(&mdl->pts_eval[p], &mdl->load, &mdl->halfspace[i]) * mdl->load.intensity;
                        }
                } else {
                        for (int p = 0; p < mdl->npts_eval; ++p) {
                                mdl->pts_eval[p].z = heqs[i];

                                double u_top = tmk__coefficients_distributed_uz(&mdl->pts_eval[p], &mdl->load, &mdl->halfspace[i]) * mdl->load.intensity;
                                mdl->pts_eval[p].z += thicknesses[i];
                                double u_btm = tmk__coefficients_distributed_uz(&mdl->pts_eval[p], &mdl->load, &mdl->halfspace[i]) * mdl->load.intensity;

                                uzs_0[p] += u_top - u_btm;
                        }
                }
        }

        tmk__elt_reset_z(mdl->pts_eval, mdl->npts_eval);

        if (mdl->slope_angle >= (double)TMK_ELT_DEG_LIM) {
                if (tmk__elt_intermediate_evaluation(mdl, uzs, uzs_0, heqs, thicknesses, sol, force_boundary, mside, nside) != TMK_SUCCESS) {
                        goto cleanup;
                }
        } else {
                for (int i = 0; i < mdl->npts_eval; ++i) {
                        uzs[i] = uzs_0[i];
                }
        }

        printf("x[mm]\ty[mm]\tdeflection[micron]\n");
        for (int p = 0; p < mdl->npts_eval; ++p) {
                printf("%.1lf\t%.1lf\t%.1lf\n", mdl->pts_eval[p].x, mdl->pts_eval[p].y, uzs[p] * 1000);
        }

        // dbg
        // FILE *fp = fopen("path", "w");
        // if (fp == NULL) {
        //         fprintf(stderr, "[ERROR] Fail to open file %s\n", export_path);
        //         goto cleanup;
        // }
        // for (size_t p = 0; p < npts; ++p) {
        //         fprintf(fp, "%.6lf,%.6lf,%.6lf\n", eval_pts[p].x, eval_pts[p].y, uzs[p]);
        // }
        // printf("[INFO] Evaluation data saved.\n");
        // fclose(fp);

        goto cleanup;

cleanup: 
        free(uzs_0);
        free(thicknesses);
        free(heqs);
        free(uzs);
        free(sol);
        free(force_boundary);
}

static tmk_result tmk__esa_intermediate_evaluation(tmk_model *mdl, double degree, double *uzs, double *uzs_0,
                                                   double *heqs, double *thicknesses,
                                                   double *sol, tmk_vec3 *force_boundary,
                                                   size_t mside, size_t nside);

static double tmk__esa_dist_signed(double *interest, double *target, size_t n);

static void tmk__solve_esa(tmk_model *mdl)
{
        size_t npts = mdl->npts_measure;
        int n_finite = mdl->n_layers - 1;
        size_t mside = TMK_GEOM_M;
        size_t nside = TMK_GEOM_N;
        size_t n = nside * nside;

        double *uzs_0 = calloc(npts, sizeof(double));
        double *thicknesses = calloc((size_t)(n_finite), sizeof(double));
        double *heqs = calloc(mdl->n_layers, sizeof(double));
        double *sol = calloc(3 * n, sizeof(double));
        double *uzs_90 = malloc(npts * sizeof(double));
        double *uzs_m = malloc(npts * sizeof(double));
        tmk_vec3 *force_boundary = calloc(n, sizeof(tmk_vec3));
        if (uzs_0 == NULL || thicknesses == NULL || heqs == NULL 
                || sol == NULL || force_boundary == NULL
                || uzs_90 == NULL || uzs_m == NULL) {
                mdl->result = TMK_ALLOCATION_FAILURE;
                goto cleanup;
        }

        for (size_t i = 0; i < (size_t)n_finite; ++i) {
                thicknesses[i] = mdl->top_depths[i + 1] - mdl->top_depths[i];
        }

        printf("[INFO] Trial angle: 0 [deg].\n");

        for (int i = n_finite; i >= 0; --i) {
                for (int j = 0; j < i; ++j) {
                        heqs[i] += tmk__elt_equivalent_thickness(thicknesses[j], &mdl->halfspace[j], &mdl->halfspace[i]);
                }

                if (i == n_finite) {
                        for (size_t p = 0; p < npts; ++p) {
                                mdl->pts_measure[p].z = heqs[i];
                                uzs_0[p] += tmk__coefficients_distributed_uz(&mdl->pts_measure[p], &mdl->load, &mdl->halfspace[i]) * mdl->load.intensity;
                        }
                } else {
                        for (size_t p = 0; p < npts; ++p) {
                                mdl->pts_measure[p].z = heqs[i];

                                double u_top = tmk__coefficients_distributed_uz(&mdl->pts_measure[p], &mdl->load, &mdl->halfspace[i]) * mdl->load.intensity;
                                mdl->pts_measure[p].z += thicknesses[i];
                                double u_btm = tmk__coefficients_distributed_uz(&mdl->pts_measure[p], &mdl->load, &mdl->halfspace[i]) * mdl->load.intensity;
                                
                                uzs_0[p] += u_top - u_btm;
                        }
                }
        }

        tmk__elt_reset_z(mdl->pts_measure, npts);

        double diff_0 = tmk__esa_dist_signed(mdl->vals_measure, uzs_0, npts);
        printf("[INFO]   Signed difference with 0.0 [deg]: %.6lf\n", diff_0);
        if (fabs(diff_0) < (double)TMK_ESA_TOL) {
                printf("[INFO] Half-space solution accepted. Effective slope angle = 0.0 [deg]\n\n");
                goto cleanup;
        }

        double degree_90 = 90.0;
        printf("[INFO] Trial angle: 90 [deg]. ");

        if (tmk__esa_intermediate_evaluation(mdl, degree_90, uzs_90, uzs_0, heqs, thicknesses, sol, force_boundary, mside, nside) != TMK_SUCCESS)
                goto cleanup;

        double diff_90 = tmk__esa_dist_signed(mdl->vals_measure, uzs_90, npts);
        printf("[INFO]   Signed difference with 90.0 [deg]: %.6lf\n", diff_90);
        if (fabs(diff_90) < (double)TMK_ESA_TOL) {
                printf("[INFO] Equivalent slope angle = 90.0 [deg]\n\n");
                goto cleanup;
        }

        double degree_l = 0.0;
        double degree_m = 45.0;
        double degree_r = 90.0;

        int iter = 0;
        int is_success = 0;

        while (degree_m > (double)TMK_ESA_DEG_LIM && iter <= (int)TMK_ESA_ITER) {
                iter++;
                printf("[INFO] Iteration %d: trial algle = %.1lf [deg]. ", iter, degree_m);
                if (tmk__esa_intermediate_evaluation(mdl, degree_m, uzs_m, uzs_0, heqs, thicknesses, sol, force_boundary, mside, nside) != TMK_SUCCESS)
                        goto cleanup;

                // compare half-space, 45-, and 90-degree evaluations
                double diff_m = tmk__esa_dist_signed(mdl->vals_measure, uzs_m, npts);

                printf("[INFO]   Signed difference with %.1lf [deg]: %.6lf\n", degree_m, diff_m);

                if (fabs(diff_m) < (double)TMK_ESA_TOL) {
                        printf("[INFO] Equivalent slope angle = %.1lf [deg]\n\n", degree_m);
                        mdl->slope_angle = degree_m;
                        is_success = 1;
                        break;
                } else if (diff_m < 0.0) {
                        degree_r = degree_m;
                } else {
                        degree_l = degree_m;
                }

                degree_m = degree_l + (degree_r - degree_l) / 2.0;
        }

        if (is_success == 0) {
                if (degree_m <= (double)TMK_ESA_DEG_LIM) {
                        mdl->result = TMK_ESA_EXCEED_ITER_LIM;
                        goto cleanup;
                } else {
                        printf("[INFO] Equivalent slope angle < %.1lf [deg]\n", (double)TMK_ESA_DEG_LIM);
                        mdl->result = TMK_SUCCESS;
                        goto cleanup;
                }
        }

        goto cleanup;

cleanup:
        free(uzs_0);
        free(uzs_90);
        free(uzs_m);
        free(heqs);
        free(sol);
        free(force_boundary);
        free(thicknesses);
}

static void tmk__utils_uz_increment(double *uz, tmk_vec3 *evaluation_point,
                                    tmk_vec3 *force_boundary, tmk_hs *halfspace,
                                    double *sol, size_t npts_force);

static tmk_result tmk__esa_intermediate_evaluation(tmk_model *mdl, double degree, double *uzs, double *uzs_0,
                                                   double *heqs, double *thicknesses,
                                                   double *sol, tmk_vec3 *force_boundary,
                                                   size_t mside, size_t nside)
{
        int n_finite = mdl->n_layers - 1;
        size_t npts = mdl->npts_measure;
        size_t n = nside * nside;
        memcpy(uzs, uzs_0, npts * sizeof(double));
        
        printf("Evaluating");
        fflush(stdout);
        for (int i = n_finite; i >= 0; --i) {
                printf("...");
                fflush(stdout);
                if (i == n_finite) {
                        if (tmk__solve_single(sol, &mdl->halfspace[i], &mdl->load, degree, mside, nside, force_boundary, &mdl->result, &mdl->info_lapack) != TMK_SUCCESS) {
                                printf("\n");
                                return mdl->result;
                        }

                        for (size_t p = 0; p < npts; ++p) {
                                mdl->pts_measure[p].z = heqs[i];
                                tmk__utils_uz_increment(&uzs[p], &mdl->pts_measure[p], force_boundary, &mdl->halfspace[i], sol, n);
                        }
                } else {
                        if (tmk__solve_single(sol, &mdl->halfspace[i], &mdl->load, degree, mside, nside, force_boundary, &mdl->result, &mdl->info_lapack) != TMK_SUCCESS) {
                                printf("\n");
                                return mdl->result;
                        }

                        for (size_t p = 0; p < npts; ++p) {
                                mdl->pts_measure[p].z = heqs[i];

                                double u_top = 0.0;
                                double u_btm = 0.0;

                                tmk__utils_uz_increment(&u_top, &mdl->pts_measure[p], force_boundary, &mdl->halfspace[i], sol, n);
                                mdl->pts_measure[p].z += thicknesses[i];
                                tmk__utils_uz_increment(&u_btm, &mdl->pts_measure[p], force_boundary, &mdl->halfspace[i], sol, n);

                                uzs[p] += u_top - u_btm;
                        }
                }
        }

        tmk__elt_reset_z(mdl->pts_measure, npts);
        printf("\n");
        return mdl->result;
}

TMKDEF int tmk_checkout(tmk_model *mdl)
{
        switch (mdl->result) {
                case TMK_INITIAL: 
                        printf("[INFO] tmk_model initialised. Check out other APIs to set up pavement model.\n");
                        return 0;
                case TMK_NO_LOAD: 
                        fprintf(stderr, "[ERROR] Loading information is not found. Did you call `tmk_set_load`?\n");
                        return 1;
                case TMK_ESA_EXCEED_ITER_LIM:
                        fprintf(stderr, "[ERROR] Failed to find equivalent slope angle within %d iterations", (int)TMK_ESA_ITER);
                        return 1;
                default: return 0;
        }
}

void evaluation_evaluate(tmk_vec3 *pts,
         size_t npts,
         double *sol,
         tmk_vec3 *force_boundary,
         size_t n_forces,
         tmk_load *load,
         tmk_hs *halfspace,
         FILE *export)
{
        if (export != NULL) {
                fprintf(export, "x,y,z,sxx,syy,szz,syz,sxz,sxy,ux,uy,uz\n");
        }

        for (size_t p = 0; p < npts; ++p) {
                double sxx = 0, syy = 0, szz = 0, syz = 0, sxz = 0, sxy = 0;
                double ux = 0, uy = 0, uz = 0;
                double cs[6] = {0};
                double cu[3] = {0};

                // by load
                tmk__coefficients_distributed_stress(cs, &pts[p], load, halfspace);
                sxx += cs[0] * load->intensity;
                syy += cs[1] * load->intensity;
                szz += cs[2] * load->intensity;
                syz += cs[3] * load->intensity;
                sxz += cs[4] * load->intensity;
                sxy += cs[5] * load->intensity;

                tmk__coefficients_distributed_displacement(cu, &pts[p], load, halfspace);
                ux += cu[0] * load->intensity;
                uy += cu[1] * load->intensity;
                uz += cu[2] * load->intensity;

                // by force boundary
                for (size_t j = 0; j < n_forces; ++j) {
                        // px
                        tmk__coefficients_concentrated_stress_x(cs, &pts[p], &force_boundary[j], halfspace);
                        sxx += cs[0] * sol[0 * n_forces + j];
                        syy += cs[1] * sol[0 * n_forces + j];
                        szz += cs[2] * sol[0 * n_forces + j];
                        syz += cs[3] * sol[0 * n_forces + j];
                        sxz += cs[4] * sol[0 * n_forces + j];
                        sxy += cs[5] * sol[0 * n_forces + j];
                        tmk__coefficients_concentrated_displacement_x(cu, &pts[p], &force_boundary[j], halfspace);
                        ux += cu[0] * sol[0 * n_forces + j];
                        uy += cu[1] * sol[0 * n_forces + j];
                        uz += cu[2] * sol[0 * n_forces + j];

                        // py
                        tmk__coefficients_concentrated_stress_y(cs, &pts[p], &force_boundary[j], halfspace);
                        sxx += cs[0] * sol[1 * n_forces + j];
                        syy += cs[1] * sol[1 * n_forces + j];
                        szz += cs[2] * sol[1 * n_forces + j];
                        syz += cs[3] * sol[1 * n_forces + j];
                        sxz += cs[4] * sol[1 * n_forces + j];
                        sxy += cs[5] * sol[1 * n_forces + j];
                        tmk__coefficients_concentrated_displacement_y(cu, &pts[p], &force_boundary[j], halfspace);
                        ux += cu[0] * sol[1 * n_forces + j];
                        uy += cu[1] * sol[1 * n_forces + j];
                        uz += cu[2] * sol[1 * n_forces + j];

                        // pz
                        tmk__coefficients_concentrated_stress_z(cs, &pts[p], &force_boundary[j], halfspace);
                        sxx += cs[0] * sol[2 * n_forces + j];
                        syy += cs[1] * sol[2 * n_forces + j];
                        szz += cs[2] * sol[2 * n_forces + j];
                        syz += cs[3] * sol[2 * n_forces + j];
                        sxz += cs[4] * sol[2 * n_forces + j];
                        sxy += cs[5] * sol[2 * n_forces + j];

                        tmk__coefficients_concentrated_displacement_z(cu, &pts[p], &force_boundary[j], halfspace);
                        ux += cu[0] * sol[2 * n_forces + j];
                        uy += cu[1] * sol[2 * n_forces + j];
                        uz += cu[2] * sol[2 * n_forces + j];
                }

                // results
                if (export == NULL) {
                        printf("Point %zu: (%.0lf, %.0lf, %.0lf)\n", p + 1, pts[p].x, pts[p].y, pts[p].z);
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
                } else {
                        fprintf(export, "%.0lf,%.0lf,%.0lf,\
                                         %.3e,%.3e,%.3e,\
                                         %.3e,%.3e,%.3e,\
                                         %.3e,%.3e,%.3e\n", 
                                         pts[p].x, pts[p].y, pts[p].z, 
                                         sxx, syy, szz, 
                                         syz, sxz, sxy, 
                                         ux, uy, uz);
                }
        }
}

static void tmk__geometry_grid_uniform(tmk_vec3 *grid, size_t n_per_side, double h, double degree, tmk_load *load)
{
        double radian = degree * M_PI / 180.0;
        double hx;
        double s = sin(radian);
        if (s == 0.0)
        {
        fprintf(stderr, "%s", "0 in denominator!\n");
        exit(EXIT_FAILURE);
        }

        hx = h / s;

        for (size_t i = 0; i < n_per_side; i++)
        {
        double x = hx * (double)((int)i + 1 - (int)n_per_side);

        for (size_t j = 0; j < n_per_side; j++)
        {
        double y = h * (double)(2 * (int)j + 1 - (int)n_per_side) / 2.0;

        size_t idx = i * n_per_side + j;

        grid[idx].x = x * cos(-1.0 * radian);
        if (fabs(fabs(grid[idx].x - load->center.x) - load->half_width) < 1.0e-9)
        {
        grid[idx].x += 0.001;
        }

        grid[idx].y = y;
        if (fabs(fabs(grid[idx].y - load->center.y) - load->half_width) < 1.0e-9)
        {
        grid[idx].y += 0.001;
        }

        grid[idx].z = x * sin(-1.0 * radian);
        }
        }
}

static void tmk__geometry_grid_offset(tmk_vec3 *grid, size_t npts, double beta, double x0)
{
        double offset = beta * x0;

        for (size_t i = 0; i < npts; i++)
        {
        grid[i].x -= offset;
        }
}

static double tmk__elt_equivalent_thickness(double current_thickness, tmk_hs *current, tmk_hs *target)
{
        double ratio = current->youngs_modulus * (1.0 - pow(target->poissons_ratio, 2.0)) / target->youngs_modulus / (1.0 - pow(current->poissons_ratio, 2.0));
        double heq = current_thickness * cbrt(ratio) * (double)TMK_ELT_FACTOR;

        return heq;
}

static void tmk__elt_reset_z(tmk_vec3 *a, size_t npts)
{
        for (size_t p = 0; p < npts; ++p) {
                a[p].z = 0.0;
        }
}

double tmk__esa_dist_signed(double *interest,
           double *target,
           size_t n)
{
        double a = 0.0;

        for (size_t i = 0; i < n; ++i)
        {
        a += interest[i] - target[i];
        }

        return a;
}

void tmk__utils_uz_increment(double *uz,
          tmk_vec3 *evaluation_point,
          tmk_vec3 *force_boundary,
          tmk_hs *halfspace,
          double *sol,
          size_t npts_force)
{
        double coeffs[3] = {0.0};

        for (size_t i = 0; i < npts_force; i++)
        {
        tmk__coefficients_concentrated_displacement_x(coeffs, evaluation_point, &force_boundary[i], halfspace);
        *uz += coeffs[2] * sol[i];

        tmk__coefficients_concentrated_displacement_y(coeffs, evaluation_point, &force_boundary[i], halfspace);
        *uz += coeffs[2] * sol[npts_force + i];

        tmk__coefficients_concentrated_displacement_z(coeffs, evaluation_point, &force_boundary[i], halfspace);
        *uz += coeffs[2] * sol[2 * npts_force + i];
        }
}

static tmk_result tmk__elt_intermediate_evaluation(tmk_model *mdl, double *uzs, double *uzs_0, 
                                                   double *heqs, double *thicknesses, 
                                                   double *sol, tmk_vec3 *force_boundary,
                                                   size_t mside, size_t nside)
{
        int n_finite = mdl->n_layers - 1;
        size_t n = nside * nside;
        memcpy(uzs, uzs_0, mdl->npts_eval * sizeof(double));

        printf("[INFO] Evaluating");
        fflush(stdout);
        for (int i = n_finite; i >= 0; --i) {
                printf("...");
                fflush(stdout);
                if (i == n_finite) {
                        if (tmk__solve_single(sol, &mdl->halfspace[i], &mdl->load, 
                                              mdl->slope_angle, mside, nside, force_boundary, 
                                              &mdl->result, &mdl->info_lapack) != TMK_SUCCESS) {
                                printf("\n");
                                return mdl->result;
                        }

                        for (int p = 0; p < mdl->npts_eval; ++p) {
                                mdl->pts_eval[p].z = heqs[i];
                                tmk__utils_uz_increment(&uzs[p], &mdl->pts_eval[p], force_boundary, &mdl->halfspace[i], sol, n);
                        }
                } else {
                        if (tmk__solve_single(sol, &mdl->halfspace[i], &mdl->load, 
                                              mdl->slope_angle, mside, nside, force_boundary,
                                              &mdl->result, &mdl->info_lapack) != TMK_SUCCESS) {
                                printf("\n");
                                return mdl->result;
                        }

                        for (int p = 0; p < mdl->npts_eval; ++p) {
                                mdl->pts_eval[p].z = heqs[i];

                                double u_top = 0.0;
                                double u_btm = 0.0;

                                tmk__utils_uz_increment(&u_top, &mdl->pts_eval[p], force_boundary, &mdl->halfspace[i], sol, n);
                                mdl->pts_eval[p].z += thicknesses[i];
                                tmk__utils_uz_increment(&u_btm, &mdl->pts_eval[p], force_boundary, &mdl->halfspace[i], sol, n);

                                uzs[p] += u_top - u_btm;
                        }
                }
        }

        tmk__elt_reset_z(mdl->pts_eval, mdl->npts_eval);
                                
        printf("\n[INFO] Evaluated.\n");
        return mdl->result;
}

void boundary_properties(double *hm, double *hn, double *beta, double x0)
{
        *hm = 2.0 * 0.08 * x0;
        *hn = 2.0 * 0.1 * x0;
        *beta = 1.25 * *hn / x0;
}

static tmk_result tmk__solve_single(double *sol, tmk_hs *halfspace, tmk_load *load, 
                                    double degree, size_t mside, size_t nside, tmk_vec3 *force_boundary,
                                    tmk_result *result, int *info_lapack)
{
        double x0 = load->center.x;
        double q = load->intensity;
        
        double hm, hn, beta;
        // tmk__geometry_grid_properties(x0, &hm, &hn, &beta);
        boundary_properties(&hm, &hn, &beta, x0);
        
        int m = mside * mside;
        int n = nside * nside;
        
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

        tmk__geometry_grid_uniform(simulated_boundary, mside, hm, degree, load);
        tmk__geometry_grid_uniform(force_boundary, nside, hn, degree, load);
        tmk__geometry_grid_offset(force_boundary, n, beta, x0);

        for (int j = 0; j < n; ++j) {
                for (int i = 0; i < m; ++i) {
                        size_t idx = j * m + i;
                        double coeffs[6] = {0.0};

                        tmk__coefficients_concentrated_stress_x(coeffs, &simulated_boundary[i], &force_boundary[j], halfspace);
                        tmk__coefficients_tensor_transformation_y(coeffs, degree);

                        s33_px[idx] = coeffs[2];
                        s23_px[idx] = coeffs[3];
                        s13_px[idx] = coeffs[4];

                        tmk__coefficients_concentrated_stress_y(coeffs, &simulated_boundary[i], &force_boundary[j], halfspace);
                        tmk__coefficients_tensor_transformation_y(coeffs, degree);

                        s33_py[idx] = coeffs[2];
                        s23_py[idx] = coeffs[3];
                        s13_py[idx] = coeffs[4];

                        tmk__coefficients_concentrated_stress_z(coeffs, &simulated_boundary[i], &force_boundary[j], halfspace);
                        tmk__coefficients_tensor_transformation_y(coeffs, degree);

                        s33_pz[idx] = coeffs[2];
                        s23_pz[idx] = coeffs[3];
                        s13_pz[idx] = coeffs[4];
                }
        }

        tmk__memory_lhs_fill(lhs, m, n, s33_px, s33_py, s33_pz, s23_px, s23_py, s23_pz, s13_px, s13_py, s13_pz);

        for (int i = 0; i < m; ++i) {
                double coeffs[6] = {0.0};

                tmk__coefficients_distributed_stress(coeffs, &simulated_boundary[i], load, halfspace);
                tmk__coefficients_tensor_transformation_y(coeffs, degree);
                s33_q[i] = coeffs[2];
                s23_q[i] = coeffs[3];
                s13_q[i] = coeffs[4];
        }

        memcpy(rhs, s33_q, m * sizeof(double));
        memcpy(rhs + m, s23_q, m * sizeof(double));
        memcpy(rhs + 2 * m, s13_q, m * sizeof(double));

        cblas_dscal(3 * (int)m, -1.0 * q, rhs, 1);

        *info_lapack = LAPACKE_dgels(LAPACK_COL_MAJOR, 'N', 3 * (int)m, 3 * (int)n, 1, lhs, 3 * (int)m, rhs, 3 * (int)m);

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

        return *result;
}

static void tmk__memory_lhs_fill(double *lhs, size_t m, size_t n, double *s33_px, double *s33_py, double *s33_pz, double *s23_px, double *s23_py, double *s23_pz, double *s13_px, double *s13_py, double *s13_pz)
{

        double *blocks[3][3] = {
        {s33_px, s33_py, s33_pz},
        {s23_px, s23_py, s23_pz},
        {s13_px, s13_py, s13_pz}};

        for (size_t i = 0; i < 3; i++)
        {
        for (size_t j = 0; j < 3; j++)
        {
        double *src = blocks[i][j];
        for (size_t k = 0; k < n; k++)
        {
        double *src_col = src + k * m;
        double *dst = lhs + (j * n + k) * (3 * m) + i * m;
        memcpy(dst, src_col, m * sizeof(double));
        }
        }
        }
}

// static void tmk__boundary_dynamic(double *len_side_final, double *len_btm_final, tmk_load *load, tmk_hs *halfspace, double degree, double threshold, double h_default)
// {
//         int iter = 0;

//         double radian = degree * M_PI / 180.0;

//         *len_side_final = 0.0;
//         *len_btm_final = 0.0;

//         for (size_t k = 2; k < 5; k++)
//         {
//         int side_not_exceed = 0;
//         int btm_not_exceed = 0;

//         double len_side;
//         double len_btm;

//         double inc_side;
//         double inc_btm;

//         iter = 0;

//         len_side = 2.0 * load->half_width + 2.0 * load->center.x;
//         len_btm = len_side / 2.0;

//         inc_side = len_side / 2.0;
//         inc_btm = len_btm / 2.0;

//         printf("%-8sForward calculation of the simulated boundary area for stress %zu3 ...\n", "[INFO]", (5 - k));

//         do
//         {
//         size_t nstep_side = (size_t)(len_side / h_default);
//         size_t nstep_btm = (size_t)(len_btm / h_default);

//         double h_side = len_side / (double)nstep_side;
//         double h_btm = len_btm / (double)nstep_btm;

//         size_t n_side = nstep_side;
//         size_t n_btm = nstep_btm + 1;

//         iter += 1;

//         if (!side_not_exceed)
//         {
//         side_not_exceed = 1;

//         for (size_t i = 0; i < n_side; i++)
//         {
//         double stress[6] = {0.0};

//         tmk_vec3 side_pt = {0.0, 0.0, 0.0};

//         side_pt.x = -1.0 * ((double)i * h_side) * cos(-1.0 * radian);
//         side_pt.y = len_btm;
//         side_pt.z = -1.0 * ((double)i * h_side) * sin(-1.0 * radian);

//         tmk__coefficients_distributed_stress(stress, &side_pt, load, halfspace);
//         cblas_dscal(6, load->intensity, stress, 1);
//         tmk__coefficients_tensor_transformation_y(stress, degree);

//         if (fabs(stress[k] / load->intensity) > threshold)
//         {
//         // printf("%zu/%zu, (%.2f, %.2f, %.2f), %.3e\n", i, n_side, side_pt.x, side_pt.y, side_pt.z, fabs(stress[k] / load->intensity));
//         side_not_exceed = 0;
//         break;
//         }
//         }
//         }

//         if (!btm_not_exceed)
//         {
//         btm_not_exceed = 1;

//         for (size_t i = 0; i < n_btm; i++)
//         {
//         double stress[6] = {0.0};

//         tmk_vec3 btm_pt = {0.0, 0.0, 0.0};

//         btm_pt.x = -1.0 * len_side * cos(-1.0 * radian);
//         btm_pt.y = (double)i * h_btm;
//         btm_pt.z = -1.0 * len_side * sin(-1.0 * radian);

//         tmk__coefficients_distributed_stress(stress, &btm_pt, load, halfspace);
//         cblas_dscal(6, load->intensity, stress, 1);
//         tmk__coefficients_tensor_transformation_y(stress, degree);

//         if (fabs(stress[k] / load->intensity) > threshold)
//         {
//         // printf("%zu/%zu, (%.2f, %.2f, %.2f), %.3e\n", i, n_btm, btm_pt.x, btm_pt.y, btm_pt.z, fabs(stress[k] / load->intensity));
//         btm_not_exceed = 0;
//         break;
//         }
//         }
//         }

//         printf("[INFO] Iteration %d, ", iter);
//         printf("x: %d, y: %d\n", side_not_exceed, btm_not_exceed);

//         if (!side_not_exceed)
//         {
//         len_btm += inc_btm;
//         }

//         if (!btm_not_exceed)
//         {
//         len_side += inc_side;
//         }

//         } while (!(side_not_exceed && btm_not_exceed));

//         if (iter == 1)
//         {
//         printf("%-8sBackward calculation of the simulated boundary area for stress %zu3 ...\n", "[INFO]", (5 - k));

//         // len_btm -= h_default;
//         // len_side -= h_default;

//         do
//         {
//         if (side_not_exceed)
//         {
//         if (!(len_btm < h_default))
//         {
//         len_btm -= h_default;
//         }
//         }

//         if (btm_not_exceed)
//         {
//         if (!(len_side < h_default))
//         {
//         len_side -= h_default;
//         }
//         }

//         size_t nstep_side = (size_t)(len_side / h_default);
//         size_t nstep_btm = (size_t)(len_btm / h_default);

//         double h_side = len_side / (double)nstep_side;
//         double h_btm = len_btm / (double)nstep_btm;

//         size_t n_side = nstep_side;
//         size_t n_btm = nstep_btm + 1;

//         iter += 1;

//         if (side_not_exceed)
//         {
//         for (size_t i = 0; i < n_side; i++)
//         {
//         double stress[6] = {0.0};

//         tmk_vec3 side_pt = {0.0, 0.0, 0.0};

//         side_pt.x = -1.0 * ((double)i * h_side) * cos(-1.0 * radian);
//         side_pt.y = len_btm;
//         side_pt.z = -1.0 * ((double)i * h_side) * sin(-1.0 * radian);

//         tmk__coefficients_distributed_stress(stress, &side_pt, load, halfspace);
//         cblas_dscal(6, load->intensity, stress, 1);
//         tmk__coefficients_tensor_transformation_y(stress, degree);

//         if (fabs(stress[k] / load->intensity) > threshold)
//         {
//         // printf("%zu/%zu, (%.2f, %.2f, %.2f), %.3e\n", i, n_side, side_pt.x, side_pt.y, side_pt.z, fabs(stress[k] / load->intensity));
//         len_btm += h_default;
//         side_not_exceed = 0;
//         break;
//         }
//         }
//         }

//         if (btm_not_exceed)
//         {
//         for (size_t i = 0; i < n_btm; i++)
//         {
//         double stress[6] = {0.0};

//         tmk_vec3 btm_pt = {0.0, 0.0, 0.0};

//         btm_pt.x = -1.0 * len_side * cos(-1.0 * radian);
//         btm_pt.y = (double)i * h_btm;
//         btm_pt.z = -1.0 * len_side * sin(-1.0 * radian);

//         tmk__coefficients_distributed_stress(stress, &btm_pt, load, halfspace);
//         cblas_dscal(6, load->intensity, stress, 1);
//         tmk__coefficients_tensor_transformation_y(stress, degree);

//         if (fabs(stress[k] / load->intensity) > threshold)
//         {
//         // printf("%zu/%zu, (%.2f, %.2f, %.2f), %.3e\n", i, n_side, btm_pt.x, btm_pt.y, btm_pt.z, fabs(stress[k] / load->intensity));
//         len_side += h_default;
//         btm_not_exceed = 0;
//         break;
//         }
//         }
//         }

//         printf("%-8sIteration %d, ", "[INFO]", iter);
//         printf("x: %d, y: %d\n", side_not_exceed, btm_not_exceed);

//         } while ((side_not_exceed && !(len_btm < h_default)) || (btm_not_exceed && !(len_side < h_default)));

//         printf("%-8sBackward calculation, return %d, %d, %d, %d\n", "[INFO]", side_not_exceed, btm_not_exceed, len_side < h_default, len_btm < h_default);
//         }
//         else
//         {
//         printf("%-8sBackward calculation of the simulated boundary area for stress %zu3 ...\n", "[INFO]", (5 - k));

//         do
//         {
//         if (side_not_exceed)
//         {
//         len_btm -= h_default;
//         }

//         if (btm_not_exceed)
//         {
//         len_side -= h_default;
//         }

//         size_t nstep_side = (size_t)(len_side / h_default);
//         size_t nstep_btm = (size_t)(len_btm / h_default);

//         double h_side = len_side / (double)nstep_side;
//         double h_btm = len_btm / (double)nstep_btm;

//         size_t n_side = nstep_side;
//         size_t n_btm = nstep_btm + 1;

//         iter += 1;

//         if (side_not_exceed)
//         {
//         for (size_t i = 0; i < n_side; i++)
//         {
//         double stress[6] = {0.0};

//         tmk_vec3 side_pt = {0.0, 0.0, 0.0};

//         side_pt.x = -1.0 * ((double)i * h_side) * cos(-1.0 * radian);
//         side_pt.y = len_btm;
//         side_pt.z = -1.0 * ((double)i * h_side) * sin(-1.0 * radian);

//         tmk__coefficients_distributed_stress(stress, &side_pt, load, halfspace);
//         cblas_dscal(6, load->intensity, stress, 1);
//         tmk__coefficients_tensor_transformation_y(stress, degree);

//         if (fabs(stress[k] / load->intensity) > threshold)
//         {
//         side_not_exceed = 0;
//         len_btm += h_default;
//         break;
//         }
//         }
//         }

//         if (btm_not_exceed)
//         {
//         for (size_t i = 0; i < n_btm; i++)
//         {
//         double stress[6] = {0.0};

//         tmk_vec3 btm_pt = {0.0, 0.0, 0.0};

//         btm_pt.x = -1.0 * len_side * cos(-1.0 * radian);
//         btm_pt.y = (double)i * h_btm;
//         btm_pt.z = -1.0 * len_side * sin(-1.0 * radian);

//         tmk__coefficients_distributed_stress(stress, &btm_pt, load, halfspace);
//         cblas_dscal(6, load->intensity, stress, 1);
//         tmk__coefficients_tensor_transformation_y(stress, degree);

//         if (fabs(stress[k] / load->intensity) > threshold)
//         {
//         btm_not_exceed = 0;
//         len_side += h_default;
//         break;
//         }
//         }
//         }

//         printf("%-8sIteration %d, ", "[INFO]", iter);
//         printf("x: %d, y: %d\n", side_not_exceed, btm_not_exceed);

//         } while (side_not_exceed || btm_not_exceed);

//         printf("%-8sBackward calculation, return %d, %d\n", "[INFO]", side_not_exceed, btm_not_exceed);
//         }

//         if (len_side > *len_side_final)
//         {
//         *len_side_final = len_side;
//         }

//         if (len_btm > *len_btm_final)
//         {
//         *len_btm_final = len_btm;
//         }

//         printf("\n%-8sTest for stress %zu3 finished after iteration %d\n", "[INFO]", (5 - k), iter);
//         printf("%-8sFor stress %zu3, required domain of calculation is: x = %.1f [mm], y = %.1f [mm], viewed perpendicular to the area\n\n", "[INFO]", (5 - k), len_side, len_btm);
//         }
//         printf("%-8sFinal result: x = %.1f [mm], y = %.1f [mm]\n", "[INFO]", *len_side_final, *len_btm_final);
// }

#endif // TAAMAK_IMPLEMENTATION

/*
        Revision History:

        1.0.0-rc (2025-09-09) first release candidate
*/

/*
        Version Conventions:
        We follow https://semver.org/, so that the version number wll be in the 
        format MAJOR.MINOR.PATCH. 
        - MAJOR version when incompatible API changes are made
        - MINOR version when functionality in a backward compatible manner is added
        - PATCH version when backward compatible bug fixes are done

        Naming Conventions:
        - All publicly available interfaces and macros should be prefixed with 
        `tmk_` or `TMK_`, depending on the cases.
        - Internal functions should be prefixed with `tmk__` (double underscore)
*/

/*
        ----------------------------------------------------------------------------
        This software is under MIT License, with additional notice for dependencies.
        ----------------------------------------------------------------------------
        MIT License
        Copyright (c) 2025 Eyal Levenberg, Fang Zeyuan

        Permission is hereby granted, free of charge, to any person obtaining a 
        copy of this software and associated documentation files (the “Software”), 
        to deal in the Software without restriction, including without limitation 
        the rights to use, copy, modify, merge, publish, distribute, sublicense, 
        and/or sell copies of the Software, and to permit persons to whom the 
        Software is furnished to do so, subject to the following conditions:
        The above copyright notice and this permission notice shall be included 
        in all copies or substantial portions of the Software.
        THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
        IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
        FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
        AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
        LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
        OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS 
        IN THE SOFTWARE.
        ----------------------------------------------------------------------------
        Additional notice for dependencies:

        This project, "taamak.h", depends on OpenBLAS (licensed under a 
        BSD-style license). Include OpenBLAS’s license text and notices where 
        required when distributing combined works.
*/

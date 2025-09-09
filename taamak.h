/*
        taamak.h - v1.0.0 - https://github.com/veslrs/taamak.h

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

#define TMK__VERSION "version 1.0.0"

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

void boundary_dynamic(double *len_side_final, double *len_btm_final, tmk_load *load, tmk_hs *halfspace, double degree, double threshold, double h_default);

// Evaluation
void evaluation_evaluate(tmk_vec3 *pts,
        size_t npts,
        double *sol,
        tmk_vec3 *force_boundary,
        size_t n_forces,
        tmk_load *load,
        tmk_hs *halfspace,
        FILE *export);
void evaluation_evaluate_extrapolate(tmk_vec3 *pts,
        size_t npts,
        double *sol,
        tmk_vec3 *force_boundary,
        size_t n_forces,
        tmk_load *load,
        tmk_hs *halfspace,
        double *uz);
void evaluation_evaluate_hs(tmk_vec3 *pts,
           size_t npts,
           tmk_load *load,
           tmk_hs *halfspace,
           FILE *export);
void evaluation_free_surface(tmk_vec3 *free_surface, tmk_vec3 *force_boundary, double width, size_t n_per_side, double degree, tmk_load *load, tmk_hs *halfspace, size_t npts_force, double *sol);
void evaluation_point(tmk_vec3 *evaluation_point, tmk_vec3 *force_boundary, tmk_load *load, tmk_hs *halfspace, double *sol, size_t npts_force);
void evaluation_line(tmk_vec3 *evaluation_line, tmk_vec3 *force_boundary, tmk_vec3 *lower, tmk_vec3 *upper, size_t n_eval, tmk_load *load, tmk_hs *halfspace, double *uz, double *sol, size_t npts_force);
// TODO: evaluation_area()

void geometry_uniform_grid(tmk_vec3 *grid, size_t n_per_side, double h, double degree, tmk_load *load);
void geometry_uniform_test(tmk_vec3 *grid, double width, size_t n_per_side, double degree, tmk_load *load);
void geometry_uniform_line(tmk_vec3 *array, tmk_vec3 *a, tmk_vec3 *b, size_t npts, tmk_load *load);

void geometry_utils_offset(tmk_vec3 *grid, size_t npts, double beta, double x0);
void geometry_utils_to_view(tmk_vec3 *grid, size_t n_per_side, double degree);

void utils_memory_fill_lhs(double *lhs, size_t m, size_t n, double *s33_px, double *s33_py, double *s33_pz, double *s23_px, double *s23_py, double *s23_pz, double *s13_px, double *s13_py, double *s13_pz);
/**
 * [met] Method of Equivalent Thickness
 */

#define MET_TOL 1e-3
#define MET_EPSILON 1.0e-9
#define MET_ITER 20

typedef struct
{
        tmk_load load;
        size_t npts;
        tmk_vec3 *pts;
        double *vals;
} met_test_t;

typedef struct
{
        int n_layers;
        tmk_hs *hs;
        double *top_depths;
        met_test_t test;
} met_t;

int met_init(met_t *met,
        int n_layers,
        tmk_hs *hs,
        double *top_depths,
        tmk_load *load,
        const char *test_results);

int met_read_csv(met_t *met, const char *filepath);

void met_cleanup(met_t *met);

// int met_solver_evaluation(double *uzs,
//          met_t *met,
//          double *heqs,
//          double *sol,
//          double degree,
//          size_t mside,
//          size_t nside,
//          tmk_vec3 *force_boundary,
//          double *thicknesses,
//          double *uzs_0);

double met_dist_signed(double *interest,
          double *target,
          size_t n);

void met_uz_increment(double *uz,
         tmk_vec3 *evaluation_point,
         tmk_vec3 *force_boundary,
         tmk_hs *halfspace,
         double *sol,
         size_t npts_force);

// int met_esa_solve(met_t *met, double *angle);
        
// int met_evaluate_extrapolate(met_t *met,
//         double angle,
//         tmk_vec3 *eval_pts,
//         size_t npts,
//         size_t nside,
//         const char *export_path);

double met_backcalc_objective(unsigned n,
         const double *x,
         double *grad,
         void *my_func_data);

double met_backcalc_constraint(unsigned n,
          const double *x,
          double *grad,
          void *data);

void met_backcalc_hs_to_double(double *x,
          tmk_hs *hs,
          size_t n_layers);

void met_backcalc_double_to_hs(const double *x,
          tmk_hs *hs,
          size_t n_layers);
          
int met_backcalc_solve(met_t *met, double tol);

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

#define TMK__FMT_INDENT "       "

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
                        printf("[INFO] Model Detail:\n" \
                               TMK__FMT_INDENT "Layer compositions and material properties:\n" \
                               TMK__FMT_INDENT "    Layer #\tYoung's modulus [MPa]\tPoisson's ratio [-]\tThickness[mm]\n");
                        tmk__evaluate_elt(mdl);
                } else {
                        printf(TMK__FMT_INDENT "Solve effective slope angle for a %d-layer linear elastic pavement system.\n" \
                               TMK__FMT_INDENT "The pavement structure is subject to load normal to the surface,\n" \
                               TMK__FMT_INDENT "distributed in a square patch. The load is considered close to the edge\n" \
                               TMK__FMT_INDENT "of the pavement.\n", mdl->n_layers);
                        printf("[INFO] Model Detail:\n" \
                               TMK__FMT_INDENT "Layer compositions and material properties:\n" \
                               TMK__FMT_INDENT "    Layer #\tYoung's modulus [MPa]\tPoisson's ratio [-]\tThickness[mm]\n");
                        tmk__solve_esa(mdl);
                }
        }
}

#ifndef TMK_ELT_DEG_LIM
#define TMK_ELT_DEG_LIM 10.0
#endif
#ifndef TMK_GEOM_M
#define TMK_GEOM_M 50
#endif
#ifndef TMK_GEOM_N
#define TMK_GEOM_N 40
#endif

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

        for (int p = 0; p < mdl->npts_eval; ++p) {
                printf("%.2f,%.2lf,%.2lf\n", mdl->pts_eval[p].x, mdl->pts_eval[p].y, uzs[p] * 1000);
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

        double diff_0 = met_dist_signed(mdl->vals_measure, uzs_0, npts);
        printf("[INFO]   Signed difference with 0.0 [deg]: %.6lf\n", diff_0);
        if (fabs(diff_0) < (double)MET_TOL) {
                printf("[INFO] Half-space solution accepted. Effective slope angle = 0.0 [deg]\n\n");
                goto cleanup;
        }

        double degree_90 = 90.0;
        printf("[INFO] Trial angle: 90 [deg]. ");

        if (tmk__esa_intermediate_evaluation(mdl, degree_90, uzs_90, uzs_0, heqs, thicknesses, sol, force_boundary, mside, nside) != TMK_SUCCESS)
                goto cleanup;

        double diff_90 = met_dist_signed(mdl->vals_measure, uzs_90, npts);
        printf("[INFO]   Signed difference with 90.0 [deg]: %.6lf\n", diff_90);
        if (fabs(diff_90) < (double)MET_TOL) {
                printf("[INFO] Equivalent slope angle = 90.0 [deg]\n\n");
                goto cleanup;
        }

        double degree_l = 0.0;
        double degree_m = 45.0;
        double degree_r = 90.0;

        int iter = 0;
        int is_success = 0;

        while (degree_m > (double)TMK_ELT_DEG_LIM && iter <= (int)MET_ITER) {
                iter++;
                printf("[INFO] Iteration %d: trial algle = %.1lf [deg]. ", iter, degree_m);
                if (tmk__esa_intermediate_evaluation(mdl, degree_m, uzs_m, uzs_0, heqs, thicknesses, sol, force_boundary, mside, nside) != TMK_SUCCESS)
                        goto cleanup;

                // compare half-space, 45-, and 90-degree evaluations
                double diff_m = met_dist_signed(mdl->vals_measure, uzs_m, npts);

                printf("[INFO]   Signed difference with %.1lf [deg]: %.6lf\n", degree_m, diff_m);

                if (fabs(diff_m) < (double)MET_TOL) {
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
                if (degree_m <= (double)TMK_ELT_DEG_LIM) {
                        mdl->result = TMK_ESA_EXCEED_ITER_LIM;
                        goto cleanup;
                } else {
                        printf("[INFO] Equivalent slope angle < %.1lf [deg]\n", (double)TMK_ELT_DEG_LIM);
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

static tmk_result tmk__solve_single(double *sol, tmk_hs *halfspace, tmk_load *load, 
                                    double degree, size_t mside, size_t nside, tmk_vec3 *force_boundary,
                                    tmk_result *result, int *info_lapack);

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
                                met_uz_increment(&uzs[p], &mdl->pts_measure[p], force_boundary, &mdl->halfspace[i], sol, n);
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

                                met_uz_increment(&u_top, &mdl->pts_measure[p], force_boundary, &mdl->halfspace[i], sol, n);
                                mdl->pts_measure[p].z += thicknesses[i];
                                met_uz_increment(&u_btm, &mdl->pts_measure[p], force_boundary, &mdl->halfspace[i], sol, n);

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
                        fprintf(stderr, "[ERROR] Failed to find equivalent slope angle within %d iterations", (int)MET_ITER);
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

void evaluation_evaluate_extrapolate(tmk_vec3 *pts,
         size_t npts,
         double *sol,
         tmk_vec3 *force_boundary,
         size_t n_forces,
         tmk_load *load,
         tmk_hs *halfspace,
         double *uz)
{
        for (size_t p = 0; p < npts; ++p)
        {
        *uz = 0.0;
        double cu[3] = {0};

        tmk__coefficients_distributed_displacement(cu, &pts[p], load, halfspace);
        *uz += cu[2] * load->intensity;

        // by force boundary
        for (size_t j = 0; j < n_forces; ++j)
        {
        // px
        tmk__coefficients_concentrated_displacement_x(cu, &pts[p], &force_boundary[j], halfspace);
        *uz += cu[2] * sol[0 * n_forces + j];

        // py
        tmk__coefficients_concentrated_displacement_y(cu, &pts[p], &force_boundary[j], halfspace);
        *uz += cu[2] * sol[1 * n_forces + j];

        // pz
        tmk__coefficients_concentrated_displacement_z(cu, &pts[p], &force_boundary[j], halfspace);
        *uz += cu[2] * sol[2 * n_forces + j];
        }
        }
}

void evaluation_evaluate_hs(tmk_vec3 *pts,
        size_t npts,
        tmk_load *load,
        tmk_hs *halfspace,
        FILE *export)
{
        if (export != NULL)
        {
        fprintf(export, "x,y,z,sxx,syy,szz,syz,sxz,sxy,ux,uy,uz\n");
        }

        for (size_t p = 0; p < npts; ++p)
        {
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

        // results
        if (export == NULL)
        {
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
        }
        else
        {
        fprintf(export, "%.0lf,%.0lf,%.0lf,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e\n", pts[p].x, pts[p].y, pts[p].z, sxx, syy, szz, syz, sxz, sxy, ux, uy, uz);
        }
        }
}

void evaluation_free_surface(tmk_vec3 *free_surface, tmk_vec3 *force_boundary, double width, size_t n_per_side, double degree, tmk_load *load, tmk_hs *halfspace, size_t npts_force, double *sol)
{
        size_t n = n_per_side * n_per_side;

        double *stress_33 = calloc(n, sizeof(double));
        double *stress_23 = calloc(n, sizeof(double));
        double *stress_13 = calloc(n, sizeof(double));

        FILE *zz;
        FILE *yz;
        FILE *xz;

        geometry_uniform_test(free_surface, width, n_per_side, degree, load);

        for (size_t i = 0; i < n; i++)
        {
        double c[6] = {0.0};

        tmk__coefficients_distributed_stress(c, &free_surface[i], load, halfspace);
        cblas_dscal(6, load->intensity, c, 1);
        tmk__coefficients_tensor_transformation_y(c, degree);
        stress_33[i] += c[2];
        stress_23[i] += c[3];
        stress_13[i] += c[4];
        }

        for (size_t i = 0; i < npts_force; i++)
        {
        for (size_t j = 0; j < n; j++)
        {
        double c[6] = {0.0};

        tmk__coefficients_concentrated_stress_x(c, &free_surface[j], &force_boundary[i], halfspace);
        cblas_dscal(6, sol[i], c, 1);
        tmk__coefficients_tensor_transformation_y(c, degree);
        stress_33[j] += c[2];
        stress_23[j] += c[3];
        stress_13[j] += c[4];

        tmk__coefficients_concentrated_stress_y(c, &free_surface[j], &force_boundary[i], halfspace);
        cblas_dscal(6, sol[npts_force + i], c, 1);
        tmk__coefficients_tensor_transformation_y(c, degree);
        stress_33[j] += c[2];
        stress_23[j] += c[3];
        stress_13[j] += c[4];

        tmk__coefficients_concentrated_stress_z(c, &free_surface[j], &force_boundary[i], halfspace);
        cblas_dscal(6, sol[2 * npts_force + i], c, 1);
        tmk__coefficients_tensor_transformation_y(c, degree);
        stress_33[j] += c[2];
        stress_23[j] += c[3];
        stress_13[j] += c[4];
        }
        }

        geometry_utils_to_view(free_surface, n_per_side, degree);

        zz = fopen("../free_surf/33.dat", "w");
        yz = fopen("../free_surf/23.dat", "w");
        xz = fopen("../free_surf/13.dat", "w");

        if (zz == NULL || yz == NULL || xz == NULL)
        {
        fprintf(stderr, "[ERROR] Fail to open files\n");
        exit(EXIT_FAILURE);
        }

        for (size_t i = 0; i < n_per_side; i++)
        {
        for (size_t j = 0; j < n_per_side; j++)
        {
        size_t idx = i * n_per_side + j;

        fprintf(zz, "%.2f\t%.2f\t%.6f\n", free_surface[idx].x, free_surface[idx].y, stress_33[idx]);
        fprintf(yz, "%.2f\t%.2f\t%.6f\n", free_surface[idx].x, free_surface[idx].y, stress_23[idx]);
        fprintf(xz, "%.2f\t%.2f\t%.6f\n", free_surface[idx].x, free_surface[idx].y, stress_13[idx]);
        }

        fprintf(zz, "\n");
        fprintf(yz, "\n");
        fprintf(xz, "\n");
        }

        fclose(zz);
        fclose(yz);
        fclose(xz);

        free(stress_33);
        free(stress_23);
        free(stress_13);

        printf("[STATUS] Stresses on the slope, evaluated at a uniform grid with step size of %.1f [mm], in a square area with width %.1f [mm], return SUCCESS\n", width / ((double)n_per_side - 1.0), width);
        printf("[INFO] Data exported to ./free_surf\n");
}

void evaluation_line(tmk_vec3 *evaluation_line, tmk_vec3 *force_boundary, tmk_vec3 *lower, tmk_vec3 *upper, size_t n_eval, tmk_load *load, tmk_hs *halfspace, double *uz, double *sol, size_t npts_force)
{
        FILE *fp;

        geometry_uniform_line(evaluation_line, lower, upper, n_eval, load);

        for (size_t i = 0; i < n_eval; i++)
        {
        double coeffs[3] = {0.0};
        tmk__coefficients_distributed_displacement(coeffs, &evaluation_line[i], load, halfspace);
        uz[i] += coeffs[2] * load->intensity;

        for (size_t j = 0; j < npts_force; j++)
        {
        tmk__coefficients_concentrated_displacement_x(coeffs, &evaluation_line[i], &force_boundary[j], halfspace);
        uz[i] += coeffs[2] * sol[j];

        tmk__coefficients_concentrated_displacement_y(coeffs, &evaluation_line[i], &force_boundary[j], halfspace);
        uz[i] += coeffs[2] * sol[npts_force + j];

        tmk__coefficients_concentrated_displacement_z(coeffs, &evaluation_line[i], &force_boundary[j], halfspace);
        uz[i] += coeffs[2] * sol[2 * npts_force + j];
        }
        }

        fp = fopen("../python/data/uz.csv", "w");

        if (fp == NULL)
        {
        fprintf(stderr, "[ERROR] Fail to open file\n");
        exit(EXIT_FAILURE);
        }

        for (size_t i = 0; i < n_eval; i++)
        {
        fprintf(fp, "%.2f,%.2f,%.2f,%.6f\n", evaluation_line[i].x, evaluation_line[i].y, evaluation_line[i].z, uz[i]);
        }

        fclose(fp);

        printf("[STATUS] Displacement in z-direction, evaluated at a line of %zu points, SUCCESS\n", n_eval);
        printf("[INFO] Data exported to uz.csv\n");
}

void evaluation_point(tmk_vec3 *evaluation_point, tmk_vec3 *force_boundary, tmk_load *load, tmk_hs *halfspace, double *sol, size_t npts_force)
{
        double coeffs[3] = {0.0};
        double uz = 0.0;

        tmk__coefficients_distributed_displacement(coeffs, evaluation_point, load, halfspace);
        uz += coeffs[2] * load->intensity;

        for (size_t i = 0; i < npts_force; i++)
        {
        tmk__coefficients_concentrated_displacement_x(coeffs, evaluation_point, &force_boundary[i], halfspace);
        uz += coeffs[2] * sol[i];

        tmk__coefficients_concentrated_displacement_y(coeffs, evaluation_point, &force_boundary[i], halfspace);
        uz += coeffs[2] * sol[npts_force + i];

        tmk__coefficients_concentrated_displacement_z(coeffs, evaluation_point, &force_boundary[i], halfspace);
        uz += coeffs[2] * sol[2 * npts_force + i];
        }

        printf("[STATUS] Displacement in z-direction, evaluated at (%.1f, %.1f, %.1f), SUCCESS\n", evaluation_point->x, evaluation_point->y, evaluation_point->z);
        printf("[INFO] uz = %.3e [mm]\n", uz);
}

void geometry_uniform_grid(tmk_vec3 *grid, size_t n_per_side, double h, double degree, tmk_load *load)
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

void geometry_uniform_test(tmk_vec3 *grid, double width, size_t n_per_side, double degree, tmk_load *load)
{
        double radian = degree * M_PI / 180.0;
        double h = width / ((double)n_per_side - 1.0);

        for (size_t i = 0; i < n_per_side; i++)
        {
        double x = h * (double)((int)i + 1 - (int)n_per_side);

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

// void geometry_uniform_iter(tmk_vec3 *grid, double lx, double ly, size_t nx, size_t ny, double degree)
// {
//     double radian = degree * M_PI / 180.0;
//     double hx = lx / ((double)nx - 1.0);
//     double hy = ly / ((double)ny - 1.0);

//     for (size_t i = 0; i < nx; i++)
//     {
//         double x = hx * (double)((int)i + 1 - (int)nx);

//         for (size_t j = 0; j < ny; j++)
//         {
//             double y = hy * (double)(2 * (int)j + 1 - (int)ny) / 2.0;

//             size_t idx = i * ny + j;

//             grid[idx].x = x * cos(-1.0 * radian);
//             grid[idx].y = y;
//             grid[idx].z = x * sin(-1.0 * radian);
//         }
//     }
// }

void geometry_uniform_line(tmk_vec3 *array, tmk_vec3 *a, tmk_vec3 *b, size_t npts, tmk_load *load)
{
        double dx = b->x - a->x;
        double dy = b->y - a->y;
        double dz = b->z - a->z;

        double n_step = (double)npts - 1.0;

        double hx = dx / n_step;
        double hy = dy / n_step;
        double hz = dz / n_step;

        for (size_t i = 0; i < npts; i++)
        {
        array[i].x += (a->x + (double)i * hx);
        if (fabs(fabs(array[i].x - load->center.x) - load->half_width) < 1.0e-9)
        {
        array[i].x += 0.001;
        }
        array[i].y += (a->y + (double)i * hy);
        if (fabs(fabs(array[i].y - load->center.y) - load->half_width) < 1.0e-9)
        {
        array[i].y += 0.001;
        }
        array[i].z += (a->z + (double)i * hz);
        }
}

// void geometry_uniform_dynamic(tmk_vec3 *grid, double lx, double ly, double h, double degree)
// {
//     double radian = degree * M_PI / 180.0;

//     size_t nx = (size_t)(lx / h) + 1;
//     size_t ny = (size_t)(ly / h) + 1;

//     double hx = lx / (double)(nx - 1);
//     double hy = ly / (double)(ny - 1);

//     for (size_t i = 0; i < nx; i++)
//     {
//         double x = hx * (double)((int)i + 1 - (int)nx);

//         for (size_t j = 0; j < ny; j++)
//         {
//             double y = hy * (double)(2 * (int)j + 1 - (int)ny) / 2.0;

//             size_t idx = i * ny + j;

//             grid[idx].x = x * cos(-1.0 * radian);
//             grid[idx].y = y;
//             grid[idx].z = x * sin(-1.0 * radian);
//         }
//     }
// }

void geometry_utils_offset(tmk_vec3 *grid, size_t npts, double beta, double x0)
{
        double offset = beta * x0;

        for (size_t i = 0; i < npts; i++)
        {
        grid[i].x -= offset;
        }
}

void geometry_utils_to_view(tmk_vec3 *grid, size_t n_per_side, double degree)
{
        double radian = degree * M_PI / 180.0;

        for (size_t i = 0; i < n_per_side; i++)
        {
        for (size_t j = 0; j < n_per_side; j++)
        {
        size_t idx = i * n_per_side + j;

        grid[idx].x = grid[idx].x * cos(radian) - grid[idx].z * sin(radian);
        }
        }
}

double met_backcalc_objective(unsigned n,
          const double *x,
          double *grad,
          void *my_func_data)
{
        met_t *met = (met_t *)my_func_data;
        met_backcalc_double_to_hs(x, met->hs, met->n_layers);

        (void)n;
        (void)grad;

        if (met->n_layers < 2)
        {
        fprintf(stderr, "[ERROR] At least 2 layers are needed for the calculation\n");
        }

        size_t npts = met->test.npts;

        double *uzs = calloc(npts, sizeof(double));
        if (uzs == NULL)
        {
        fprintf(stderr, "[ERROR] Fail to allocate memory\n");
        }

        int n_finite = met->n_layers - 1;

        double *thicknesses = calloc((size_t)(n_finite), sizeof(double));
        if (thicknesses == NULL)
        {
        fprintf(stderr, "[ERROR] Failed to allocate memory for thicknesses\n");
        free(uzs);
        }

        for (size_t i = 0; i < (size_t)n_finite; ++i)
        {
        thicknesses[i] = met->top_depths[i + 1] - met->top_depths[i];
        }

        double *heqs = calloc(met->n_layers, sizeof(double));
        if (heqs == NULL)
        {
        fprintf(stderr, "[ERROR] Failed to allocate memory for equivalent thicknesses\n");
        free(uzs);
        free(thicknesses);
        }

        for (int i = n_finite; i >= 0; --i)
        {
        for (int j = 0; j < i; ++j)
        {
        heqs[i] += tmk__elt_equivalent_thickness(thicknesses[j], &met->hs[j], &met->hs[i]);
        }

        if (i == n_finite)
        {
        for (size_t p = 0; p < npts; ++p)
        {
        met->test.pts[p].z = heqs[i];
        uzs[p] += tmk__coefficients_distributed_uz(&met->test.pts[p], &met->test.load, &met->hs[i]) * met->test.load.intensity;
        }
        }
        else
        {
        for (size_t p = 0; p < npts; ++p)
        {
        met->test.pts[p].z = heqs[i];

        double u_top = tmk__coefficients_distributed_uz(&met->test.pts[p], &met->test.load, &met->hs[i]) * met->test.load.intensity;
        met->test.pts[p].z += thicknesses[i];
        double u_btm = tmk__coefficients_distributed_uz(&met->test.pts[p], &met->test.load, &met->hs[i]) * met->test.load.intensity;

        uzs[p] += u_top - u_btm;
        }
        }
        }

        tmk__elt_reset_z(met->test.pts, npts);

        double abs_diff = fabs(met_dist_signed(met->test.vals, uzs, npts));

        printf("[INFO] Solving, diff = %.6lf\n", abs_diff);

        return abs_diff;
}

double met_backcalc_constraint(unsigned n,
           const double *x,
           double *grad,
           void *data)
{
        int i = *((int *)data);

        (void)n;
        (void)grad;

        return (2 * x[i + 1] - x[i]);
}

void met_backcalc_hs_to_double(double *x,
           tmk_hs *hs,
           size_t n_layers)
{
        for (size_t i = 0; i < n_layers; ++i)
        {
        x[i] = hs[i].youngs_modulus;
        x[n_layers + i] = hs[i].poissons_ratio;
        }
}

void met_backcalc_double_to_hs(const double *x,
           tmk_hs *hs,
           size_t n_layers)
{
        for (size_t i = 0; i < n_layers; ++i)
        {
        hs[i].youngs_modulus = x[i];
        hs[i].poissons_ratio = x[n_layers + i];
        }
}

int met_backcalc_solve(met_t *met, double tol)
{
        size_t nx = met->n_layers * 2;

        // Pre-allocate memory for x
        double *x = calloc(nx, sizeof(double));
        if (x == NULL)
        {
        fprintf(stderr, "[ERROR] Fail to allocate memory for x\n");
        return 1;
        }

        // convert tmk_hs *hs to double *hs: [E_1, ..., E_n, nu_1, ..., nu_n]
        met_backcalc_hs_to_double(x, met->hs, met->n_layers);

        // Set lower bounds:
        //   - Young's moduli   > 0
        //   - Poisson's ratios > 0
        double *lb = calloc(nx, sizeof(double));
        if (lb == NULL)
        {
        fprintf(stderr, "[ERROR] Fail to allocate memory for lb\n");
        free(x);
        return 1;
        }

        // Set upper bounds:
        //   - Young's moduli   < inf
        //   - Poisson's ratios < 1/2
        double *ub = calloc(nx, sizeof(double));
        if (ub == NULL)
        {
        fprintf(stderr, "[ERROR] Fail to allocate memory for ub\n");
        free(x);
        free(lb);
        return 1;
        }

        for (int i = 0; i < met->n_layers; ++i)
        {
        ub[i] = HUGE_VAL;
        ub[met->n_layers + i] = 0.5;
        }

        // Initiate nlopt solver
        nlopt_opt opt;

        // Set algorithm and dimensionality
        //   - algorithm: local derivative-free
        //                (https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/#local-derivative-free-optimization)
        //   - dimension: 2 * #layers
        opt = nlopt_create(NLOPT_LN_COBYLA, nx);

        nlopt_set_lower_bounds(opt, lb);
        nlopt_set_upper_bounds(opt, ub);
        nlopt_set_min_objective(opt, met_backcalc_objective, met);

        // Constraint data to pass into met_backcalc_constraint():
        //   - The relation of neighbouring Young's moduli: E_i >= 2 * E_{i+1},
        //     where i indicates the ith layer. Thus for a system of n layers,
        //     (n - 1) pairs of relations should be checked, thus here we pass
        //     an list of indices [0, 1, ..., (n_layers - 2)]

        for (int i = 0; i < met->n_layers - 1; ++i)
        {
        int data = i;
        nlopt_add_inequality_constraint(opt, met_backcalc_constraint, (void *)&data, 1e-8);
        }

        nlopt_set_xtol_rel(opt, tol);

        double minf;
        if (nlopt_optimize(opt, x, &minf) < 0)
        {
        fprintf(stderr, "[ERROR] nlopt failed\n");
        return 1;
        }
        else
        {
        printf("\n[INFO] Backcalculation successful, final diff = %0.10g\n\n", minf);
        }

        nlopt_destroy(opt);

        printf("Backcalculation Results:\n");
        for (int i = 0; i < met->n_layers; ++i)
        {
        printf("Layer %d:\n", i + 1);
        printf("  Young's modulus = %.0lf [MPa]\n", x[i]);
        printf("  Poisson's ratio = %.2lf [-]\n", x[met->n_layers + i]);
        }
        printf("\n");

        free(x);
        free(lb);
        free(ub);

        return 0;
}

// int met_esa_solve(met_t *met, double *angle)
// {
//         // return 0: success
//         // return 1: failure
//         if (met->n_layers < 2)
//         {
//         fprintf(stderr, "[ERROR] At least 2 layers are needed for the calculation\n");

//         // free(met->hs);
//         // free(met->top_depths);
//         // free(met->test.pts);
//         // free(met->test.vals);
//         // met_cleanup(met);
//         return 1;
//         }

//         size_t npts = met->test.npts;

//         double *uzs_0 = calloc(npts, sizeof(double));
//         if (uzs_0 == NULL)
//         {
//         fprintf(stderr, "[ERROR] Fail to allocate memory\n");

//         // free(met->hs);
//         // free(met->top_depths);
//         // free(met->test.pts);
//         // free(met->test.vals);
//         // met_cleanup(met);
//         return 1;
//         }

//         int n_finite = met->n_layers - 1;

//         double *thicknesses = calloc((size_t)(n_finite), sizeof(double));
//         if (thicknesses == NULL)
//         {
//         fprintf(stderr, "[ERROR] Failed to allocate memory for thicknesses\n");

//         // free(met->hs);
//         // free(met->top_depths);
//         // free(met->test.pts);
//         // free(met->test.vals);
//         // met_cleanup(met);
//         free(uzs_0);

//         return 1;
//         }

//         for (size_t i = 0; i < (size_t)n_finite; ++i)
//         {
//         thicknesses[i] = met->top_depths[i + 1] - met->top_depths[i];
//         }

//         double *heqs = calloc(met->n_layers, sizeof(double));
//         if (heqs == NULL)
//         {
//         fprintf(stderr, "[ERROR] Failed to allocate memory for equivalent thicknesses\n");

//         // free(met->hs);
//         // free(met->top_depths);
//         // free(met->test.pts);
//         // free(met->test.vals);

//         // met_cleanup(met);
//         free(uzs_0);
//         free(thicknesses);

//         return 1;
//         }

//         printf("[INFO] Calculating half-space solution...\n");

//         for (int i = n_finite; i >= 0; --i)
//         {
//         for (int j = 0; j < i; ++j)
//         {
//         heqs[i] += tmk__elt_equivalent_thickness(thicknesses[j], &met->hs[j], &met->hs[i]);
//         }

//         if (i == n_finite)
//         {
//         for (size_t p = 0; p < npts; ++p)
//         {
//         met->test.pts[p].z = heqs[i];
//         uzs_0[p] += tmk__coefficients_distributed_uz(&met->test.pts[p], &met->test.load, &met->hs[i]) * met->test.load.intensity;
//         }
//         }
//         else
//         {
//         for (size_t p = 0; p < npts; ++p)
//         {
//         met->test.pts[p].z = heqs[i];

//         double u_top = tmk__coefficients_distributed_uz(&met->test.pts[p], &met->test.load, &met->hs[i]) * met->test.load.intensity;
//         met->test.pts[p].z += thicknesses[i];
//         double u_btm = tmk__coefficients_distributed_uz(&met->test.pts[p], &met->test.load, &met->hs[i]) * met->test.load.intensity;

//         uzs_0[p] += u_top - u_btm;
//         }
//         }
//         }

//         tmk__elt_reset_z(met->test.pts, npts);

//         double diff_0 = met_dist_signed(met->test.vals, uzs_0, npts);
//         printf("[INFO] Signed difference with 0.0 [deg]: %.6lf\n", diff_0);
//         if (fabs(diff_0) < (double)MET_TOL)
//         {
//         printf("[INFO] Half-space solution accepted. Equivalent slope angle = 0.0 [deg]\n\n");

//         // free(met->hs);
//         // free(met->top_depths);
//         // free(met->test.pts);
//         // free(met->test.vals);

//         // met_cleanup(met);
//         free(uzs_0);
//         free(thicknesses);
//         free(heqs);

//         return 1;
//         }

//         double degree_90 = 90.0;

//         size_t mside = 50;
//         size_t nside = 40;
//         size_t n = nside * nside;

//         double *sol = calloc(3 * n, sizeof(double));
//         tmk_vec3 *force_boundary = calloc(n, sizeof(tmk_vec3));
//         if (sol == NULL || force_boundary == NULL)
//         {
//         fprintf(stderr, "[ERROR] Failed to allocate memory\n");

//         // free(met->hs);
//         // free(met->top_depths);
//         // free(met->test.pts);
//         // free(met->test.vals);

//         // met_cleanup(met);
//         free(uzs_0);
//         free(heqs);
//         free(thicknesses);

//         return 1;
//         }

//         double *uzs_90 = malloc(npts * sizeof(double));
//         double *uzs_m = malloc(npts * sizeof(double));
//         if (uzs_90 == NULL || uzs_m == NULL)
//         {
//         fprintf(stderr, "[ERROR] Failed to allocate memory\n");

//         // free(met->hs);
//         // free(met->top_depths);
//         // free(met->test.pts);
//         // free(met->test.vals);

//         // met_cleanup(met);
//         free(uzs_0);
//         free(heqs);
//         free(sol);
//         free(force_boundary);
//         free(thicknesses);

//         return 1;
//         }

//         printf("[INFO] Calculating 90-degree solution...\n");

//         if (met_solver_evaluation(uzs_90, met, heqs, sol, degree_90, mside, nside, force_boundary, thicknesses, uzs_0) != 0)
//         {
//         // free(met->hs);
//         // free(met->top_depths);
//         // free(met->test.pts);
//         // free(met->test.vals);

//         // met_cleanup(met);
//         free(uzs_0);
//         free(uzs_90);
//         free(uzs_m);
//         free(heqs);
//         free(sol);
//         free(force_boundary);
//         free(thicknesses);

//         return 1;
//         }

//         double diff_90 = met_dist_signed(met->test.vals, uzs_90, npts);
//         printf("[INFO] Signed difference with 90.0 [deg]: %.6lf\n", diff_90);
//         if (fabs(diff_90) < (double)MET_TOL)
//         {
//         printf("[INFO] Equivalent slope angle = 90.0 [deg]\n\n");

//         // free(met->hs);
//         // free(met->top_depths);
//         // free(met->test.pts);
//         // free(met->test.vals);

//         // met_cleanup(met);
//         free(uzs_0);
//         free(uzs_90);
//         free(uzs_m);
//         free(heqs);
//         free(sol);
//         free(force_boundary);
//         free(thicknesses);

//         return 1;
//         }

//         double degree_l = 0.0;
//         double degree_m = 45.0;
//         double degree_r = 90.0;

//         int iter = 0;
//         int is_success = 0;

//         while (degree_m > (double)TMK_ELT_DEG_LIM && iter <= (int)MET_ITER)
//         {
//         iter++;

//         printf("[INFO] Iteration %d: slope algle = %.1lf [deg]\n", iter, degree_m);

//         if (met_solver_evaluation(uzs_m, met, heqs, sol, degree_m, mside, nside, force_boundary, thicknesses, uzs_0) != 0)
//         {
//         // free(met->hs);
//         // free(met->top_depths);
//         // free(met->test.pts);
//         // free(met->test.vals);

//         // met_cleanup(met);
//         free(uzs_0);
//         free(uzs_90);
//         free(uzs_m);
//         free(heqs);
//         free(sol);
//         free(force_boundary);
//         free(thicknesses);

//         return 1;
//         }

//         // compare half-space, 45-, and 90-degree evaluations
//         double diff_m = met_dist_signed(met->test.vals, uzs_m, npts);

//         printf("[INFO] Signed difference with %.1lf [deg]: %.6lf\n", degree_m, diff_m);

//         if (fabs(diff_m) < (double)MET_TOL)
//         {
//         printf("[INFO] Equivalent slope angle = %.1lf [deg]\n\n", degree_m);
//         *angle = degree_m;
//         is_success = 1;
//         break;
//         }
//         else if (diff_m < 0.0)
//         {
//         degree_r = degree_m;
//         }
//         else
//         {
//         degree_l = degree_m;
//         }

//         degree_m = degree_l + (degree_r - degree_l) / 2.0;
//         }

//         if (is_success == 0)
//         {
//         if (degree_m <= (double)TMK_ELT_DEG_LIM)
//         {
//         fprintf(stderr, "[ERROR] Failed to find equivalent slope angle within %d iterations", (int)MET_ITER);

//         // free(met->hs);
//         // free(met->top_depths);
//         // free(met->test.pts);
//         // free(met->test.vals);

//         // met_cleanup(met);
//         free(uzs_0);
//         free(uzs_90);
//         free(uzs_m);
//         free(heqs);
//         free(sol);
//         free(force_boundary);
//         free(thicknesses);

//         return 1;
//         }
//         else
//         {
//         printf("[INFO] Equivalent slope angle < %.1lf [deg]\n", (double)TMK_ELT_DEG_LIM);
//         }
//         }

//         // free(met->hs);
//         // free(met->top_depths);
//         // free(met->test.pts);
//         // free(met->test.vals);

//         // met_cleanup(met);
//         free(uzs_0);
//         free(uzs_90);
//         free(uzs_m);
//         free(heqs);
//         free(sol);
//         free(force_boundary);
//         free(thicknesses);

//         return 0;
// }

int met_init(met_t *met,
         int n_layers,
         tmk_hs *hs,
         double *top_depths,
         tmk_load *load,
         const char *test_results)
{
        met->n_layers = n_layers;
        met->hs = hs;
        met->top_depths = top_depths;

        int res;
        if (test_results == NULL)
        {
        met_test_t test = {0};
        test.load = *load;
        met->test = test;
        res = 0;
        }
        else
        {
        met->test.load = *load;
        res = met_read_csv(met, test_results);
        }

        return res;
}

int met_read_csv(met_t *met, const char *filepath)
{
        // return 0: OK
        // return -1: Fail to open test redult file
        // return -2: Empty file / Invalid data structure
        // Return -3: Fail to allocate memory
        FILE *fp = fopen(filepath, "r");
        if (fp == NULL)
        {
        fprintf(stderr, "Fail to open file: %s\n", filepath);
        return -1;
        }

        double x, y, val;
        size_t count = 0;
        while (fscanf(fp, "%lf,%lf,%lf", &x, &y, &val) == 3)
        {
        count++;
        }

        if (count == 0)
        {
        fprintf(stderr, "[ERROR] Empty file or invalid data layout: %s\n", filepath);
        fclose(fp);
        return -2;
        }

        tmk_vec3 *pts = malloc(count * sizeof(tmk_vec3));
        double *vals = malloc(count * sizeof(double));
        if (pts == NULL || vals == NULL)
        {
        fprintf(stderr, "[ERROR] Fail to allocate memory for %zu entries\n", count);
        // met_cleanup(met);
        fclose(fp);
        return -3;
        }

        met->test.pts = pts;
        met->test.vals = vals;
        met->test.npts = count;

        printf("[INFO] Number of test points: %zu\n", count);

        rewind(fp);
        size_t idx = 0;
        while (fscanf(fp, "%lf,%lf,%lf", &x, &y, &val) == 3)
        {
        if (fabs(fabs(x - met->test.load.center.x) - met->test.load.half_width) < (double)MET_EPSILON)
        {
        x += 0.001;
        }

        if (fabs(fabs(y - met->test.load.center.y) - met->test.load.half_width) < (double)MET_EPSILON)
        {
        y += 0.001;
        }

        met->test.pts[idx].x = x;
        met->test.pts[idx].y = y;
        met->test.pts[idx].z = 0.0;

        met->test.vals[idx] = val;

        // printf("%lf,%lf,%lf\n", x, y, val);

        idx++;
        }

        fclose(fp);
        return 0;
}

// int met_read_csv(met_t *met, const char *filepath)
// {
//     size_t capacity = 0;

//     FILE *fp = fopen(filepath, "r");
//     if (fp == NULL)
//     {
//         fprintf(stderr, "%-8sFail to open file: %s\n", "[ERROR]", filepath);
//         return 1;
//     }

//     double x, y, val;
//     while (fscanf(fp, "%lf,%lf,%lf", &x, &y, &val) == 3)
//     {
//         if (met->test.npts == capacity)
//         {
//             size_t new_capacity = capacity ? capacity * 2 : 128;

//             tmk_vec3 *pts = realloc(met->test.pts, new_capacity * sizeof(tmk_vec3));
//             double *vals = realloc(met->test.vals, new_capacity * sizeof(double));
//             if (pts == NULL || vals == NULL)
//             {
//                 fprintf(stderr, "%-8sFail to reallocate memory\n", "[ERROR]");

//                 free(met->hs);
//                 free(met->top_depths);
//                 free(met->test.pts);
//                 free(met->test.vals);

//                 return 1;
//             }

//             met->test.pts = pts;
//             met->test.vals = vals;
//         }

//         if (fabs(fabs(x - met->test.load.center.x) - met->test.load.half_width) < (double)MET_EPSILON)
//         {
//             x += 0.001;
//         }

//         if (fabs(fabs(y - met->test.load.center.y) - met->test.load.half_width) < (double)MET_EPSILON)
//         {
//             y += 0.001;
//         }

//         met->test.pts[met->test.npts].x = x;
//         met->test.pts[met->test.npts].y = y;
//         met->test.pts[met->test.npts].z = 0.0;

//         met->test.vals[met->test.npts] = val;

//         met->test.npts++;
//     }

//     fclose(fp);

//     return 0;
// }

void met_cleanup(met_t *met)
{
        if (met->test.pts != NULL)
        {
        free(met->test.pts);
        }

        if (met->test.vals != NULL)
        {
        free(met->test.vals);
        }

        met = NULL;
}

#ifndef TMK_ELT_FACTOR
#define TMK_ELT_FACTOR 0.9
#endif

static double tmk__elt_equivalent_thickness(double current_thickness, tmk_hs *current, tmk_hs *target)
{
        double ratio = current->youngs_modulus * (1.0 - pow(target->poissons_ratio, 2.0)) / target->youngs_modulus / (1.0 - pow(current->poissons_ratio, 2.0));
        double heq = current_thickness * cbrt(ratio) * (double)TMK_ELT_FACTOR;

        return heq;
}

// int met_solver_evaluation(double *uzs,
//           met_t *met,
//           double *heqs,
//           double *sol,
//           double degree,
//           size_t mside,
//           size_t nside,
//           tmk_vec3 *force_boundary,
//           double *thicknesses,
//           double *uzs_0)
// {
//         int n_finite = met->n_layers - 1;
//         size_t npts = met->test.npts;
//         size_t n = nside * nside;

//         memcpy(uzs, uzs_0, npts * sizeof(double));

//         for (int i = n_finite; i >= 0; --i)
//         {
//         if (i == n_finite)
//         {
//         if (solver_solve(sol, &met->hs[i], &met->test.load, degree, mside, nside, force_boundary) != 0)
//         {
//         return 1;
//         }

//         for (size_t p = 0; p < npts; ++p)
//         {
//         met->test.pts[p].z = heqs[i];
//         met_uz_increment(&uzs[p], &met->test.pts[p], force_boundary, &met->hs[i], sol, n);
//         }
//         }
//         else
//         {
//         if (solver_solve(sol, &met->hs[i], &met->test.load, degree, mside, nside, force_boundary) != 0)
//         {
//         return 1;
//         }

//         for (size_t p = 0; p < npts; ++p)
//         {
//         met->test.pts[p].z = heqs[i];

//         double u_top = 0.0;
//         double u_btm = 0.0;

//         met_uz_increment(&u_top, &met->test.pts[p], force_boundary, &met->hs[i], sol, n);
//         met->test.pts[p].z += thicknesses[i];
//         met_uz_increment(&u_btm, &met->test.pts[p], force_boundary, &met->hs[i], sol, n);

//         uzs[p] += u_top - u_btm;
//         }
//         }
//         }

//         tmk__elt_reset_z(met->test.pts, npts);

//         return 0;
// }

static void tmk__elt_reset_z(tmk_vec3 *a, size_t npts)
{
        for (size_t p = 0; p < npts; ++p) {
                a[p].z = 0.0;
        }
}

double met_dist_signed(double *interest,
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

void met_uz_increment(double *uz,
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

// int met_evaluate_extrapolate(met_t *met,
//          double angle,
//          tmk_vec3 *eval_pts,
//          size_t npts,
//          size_t nside,
//          const char *export_path)
// {
//         if (met->n_layers < 2)
//         {
//         fprintf(stderr, "[ERROR] At least 2 layers are needed for the calculation\n");

//         // free(met->hs);
//         // free(met->top_depths);

//         // met_cleanup(met);
//         return 1;
//         }

//         double *uzs_0 = calloc(npts, sizeof(double));
//         if (uzs_0 == NULL)
//         {
//         fprintf(stderr, "[ERROR] Fail to allocate memory\n");

//         // free(met->hs);
//         // free(met->top_depths);

//         // met_cleanup(met);
//         return 1;
//         }

//         int n_finite = met->n_layers - 1;

//         double *thicknesses = calloc((size_t)(n_finite), sizeof(double));
//         if (thicknesses == NULL)
//         {
//         fprintf(stderr, "[ERROR] Failed to allocate memory for thicknesses\n");

//         // free(met->hs);
//         // free(met->top_depths);

//         // met_cleanup(met);
//         free(uzs_0);

//         return 1;
//         }

//         for (size_t i = 0; i < (size_t)n_finite; ++i)
//         {
//         thicknesses[i] = met->top_depths[i + 1] - met->top_depths[i];
//         }

//         double *heqs = calloc(met->n_layers, sizeof(double));
//         if (heqs == NULL)
//         {
//         fprintf(stderr, "[ERROR] Failed to allocate memory for equivalent thicknesses\n");

//         // free(met->hs);
//         // free(met->top_depths);

//         // met_cleanup(met);
//         free(uzs_0);
//         free(thicknesses);

//         return 1;
//         }

//         printf("[INFO] Generating evaluation data with slope angle %.1lf [deg] ...\n", angle);

//         for (int i = n_finite; i >= 0; --i)
//         {
//         for (int j = 0; j < i; ++j)
//         {
//         heqs[i] += tmk__elt_equivalent_thickness(thicknesses[j], &met->hs[j], &met->hs[i]);
//         }

//         if (i == n_finite)
//         {
//         for (size_t p = 0; p < npts; ++p)
//         {
//         eval_pts[p].z = heqs[i];
//         uzs_0[p] += tmk__coefficients_distributed_uz(&eval_pts[p], &met->test.load, &met->hs[i]) * met->test.load.intensity;
//         }
//         }
//         else
//         {
//         for (size_t p = 0; p < npts; ++p)
//         {
//         eval_pts[p].z = heqs[i];

//         double u_top = tmk__coefficients_distributed_uz(&eval_pts[p], &met->test.load, &met->hs[i]) * met->test.load.intensity;
//         eval_pts[p].z += thicknesses[i];
//         double u_btm = tmk__coefficients_distributed_uz(&eval_pts[p], &met->test.load, &met->hs[i]) * met->test.load.intensity;

//         uzs_0[p] += u_top - u_btm;
//         }
//         }
//         }

//         tmk__elt_reset_z(eval_pts, npts);

//         double *uzs = calloc(npts, sizeof(double));
//         if (uzs == NULL)
//         {
//         fprintf(stderr, "[ERROR] Failed to allocate memory\n");

//         // free(met->hs);
//         // free(met->top_depths);

//         // met_cleanup(met);
//         free(uzs_0);
//         free(thicknesses);
//         free(heqs);
//         free(eval_pts);

//         return 1;
//         }

//         if (angle >= 10.0)
//         {
//         size_t mside = nside * 5 / 4;
//         // size_t nside = 40;
//         size_t n = nside * nside;

//         double *sol = calloc(3 * n, sizeof(double));
//         tmk_vec3 *force_boundary = calloc(n, sizeof(tmk_vec3));
//         if (sol == NULL || force_boundary == NULL)
//         {
//         fprintf(stderr, "[ERROR] Failed to allocate memory\n");

//         // free(met->hs);
//         // free(met->top_depths);

//         // met_cleanup(met);
//         free(uzs_0);
//         free(thicknesses);
//         free(heqs);
//         free(eval_pts);
//         free(uzs);

//         return 1;
//         }

//         if (tmk__elt_intermediate_evaluation(mdl, uzs, uzs_0, heqs, thicknesses, sol, force_boundary, mside, nside) != TMK_SUCCESS)
//         {
//         // free(met->hs);
//         // free(met->top_depths);

//         // met_cleanup(met);
//         free(uzs_0);
//         free(thicknesses);
//         free(heqs);
//         free(eval_pts);
//         free(uzs);
//         free(sol);
//         free(force_boundary);

//         return 1;
//         }

//         free(sol);
//         free(force_boundary);
//         }
//         else
//         {
//         for (size_t i = 0; i < npts; ++i)
//         {
//         uzs[i] = uzs_0[i];
//         }
//         }

//         if (export_path == NULL)
//         {
//         for (size_t p = 0; p < npts; ++p)
//         {
//         printf("n = %zu, 1/n = %.3lf\n", nside * nside, 1.0 / (double)(nside * nside));
//         printf("%.1lf,%.1lf,%.3e\n", eval_pts[p].x, eval_pts[p].y, uzs[p]);
//         }
//         }
//         else
//         {
//         FILE *fp = fopen(export_path, "w");
//         if (fp == NULL)
//         {
//         fprintf(stderr, "[ERROR] Fail to open file %s\n", export_path);

//         // free(met->hs);
//         // free(met->top_depths);

//         // met_cleanup(met);
//         free(uzs_0);
//         free(thicknesses);
//         free(heqs);
//         free(eval_pts);
//         free(uzs);

//         free(uzs);

//         return 1;
//         }

//         for (size_t p = 0; p < npts; ++p)
//         {
//         fprintf(fp, "%.6lf,%.6lf,%.6lf\n", eval_pts[p].x, eval_pts[p].y, uzs[p]);
//         }

//         // free(met->hs);
//         // free(met->top_depths);

//         // met_cleanup(met);
//         free(uzs_0);
//         free(thicknesses);
//         free(heqs);
//         free(eval_pts);
//         free(uzs);

//         fclose(fp);

//         printf("[INFO] Evaluation data saved as %s\n", export_path);
//         }

//         return 0;
// }

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
                                met_uz_increment(&uzs[p], &mdl->pts_eval[p], force_boundary, &mdl->halfspace[i], sol, n);
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

                                met_uz_increment(&u_top, &mdl->pts_eval[p], force_boundary, &mdl->halfspace[i], sol, n);
                                mdl->pts_eval[p].z += thicknesses[i];
                                met_uz_increment(&u_btm, &mdl->pts_eval[p], force_boundary, &mdl->halfspace[i], sol, n);

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

        geometry_uniform_grid(simulated_boundary, mside, hm, degree, load);
        geometry_uniform_grid(force_boundary, nside, hn, degree, load);
        geometry_utils_offset(force_boundary, n, beta, x0);

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

        utils_memory_fill_lhs(lhs, m, n, s33_px, s33_py, s33_pz, s23_px, s23_py, s23_pz, s13_px, s13_py, s13_pz);

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

void utils_memory_fill_lhs(double *lhs, size_t m, size_t n, double *s33_px, double *s33_py, double *s33_pz, double *s23_px, double *s23_py, double *s23_pz, double *s13_px, double *s13_py, double *s13_pz)
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

void boundary_dynamic(double *len_side_final, double *len_btm_final, tmk_load *load, tmk_hs *halfspace, double degree, double threshold, double h_default)
{
        int iter = 0;

        double radian = degree * M_PI / 180.0;

        *len_side_final = 0.0;
        *len_btm_final = 0.0;

        for (size_t k = 2; k < 5; k++)
        {
        int side_not_exceed = 0;
        int btm_not_exceed = 0;

        double len_side;
        double len_btm;

        double inc_side;
        double inc_btm;

        iter = 0;

        len_side = 2.0 * load->half_width + 2.0 * load->center.x;
        len_btm = len_side / 2.0;

        inc_side = len_side / 2.0;
        inc_btm = len_btm / 2.0;

        printf("%-8sForward calculation of the simulated boundary area for stress %zu3 ...\n", "[INFO]", (5 - k));

        do
        {
        size_t nstep_side = (size_t)(len_side / h_default);
        size_t nstep_btm = (size_t)(len_btm / h_default);

        double h_side = len_side / (double)nstep_side;
        double h_btm = len_btm / (double)nstep_btm;

        size_t n_side = nstep_side;
        size_t n_btm = nstep_btm + 1;

        iter += 1;

        if (!side_not_exceed)
        {
        side_not_exceed = 1;

        for (size_t i = 0; i < n_side; i++)
        {
        double stress[6] = {0.0};

        tmk_vec3 side_pt = {0.0, 0.0, 0.0};

        side_pt.x = -1.0 * ((double)i * h_side) * cos(-1.0 * radian);
        side_pt.y = len_btm;
        side_pt.z = -1.0 * ((double)i * h_side) * sin(-1.0 * radian);

        tmk__coefficients_distributed_stress(stress, &side_pt, load, halfspace);
        cblas_dscal(6, load->intensity, stress, 1);
        tmk__coefficients_tensor_transformation_y(stress, degree);

        if (fabs(stress[k] / load->intensity) > threshold)
        {
        // printf("%zu/%zu, (%.2f, %.2f, %.2f), %.3e\n", i, n_side, side_pt.x, side_pt.y, side_pt.z, fabs(stress[k] / load->intensity));
        side_not_exceed = 0;
        break;
        }
        }
        }

        if (!btm_not_exceed)
        {
        btm_not_exceed = 1;

        for (size_t i = 0; i < n_btm; i++)
        {
        double stress[6] = {0.0};

        tmk_vec3 btm_pt = {0.0, 0.0, 0.0};

        btm_pt.x = -1.0 * len_side * cos(-1.0 * radian);
        btm_pt.y = (double)i * h_btm;
        btm_pt.z = -1.0 * len_side * sin(-1.0 * radian);

        tmk__coefficients_distributed_stress(stress, &btm_pt, load, halfspace);
        cblas_dscal(6, load->intensity, stress, 1);
        tmk__coefficients_tensor_transformation_y(stress, degree);

        if (fabs(stress[k] / load->intensity) > threshold)
        {
        // printf("%zu/%zu, (%.2f, %.2f, %.2f), %.3e\n", i, n_btm, btm_pt.x, btm_pt.y, btm_pt.z, fabs(stress[k] / load->intensity));
        btm_not_exceed = 0;
        break;
        }
        }
        }

        printf("[INFO] Iteration %d, ", iter);
        printf("x: %d, y: %d\n", side_not_exceed, btm_not_exceed);

        if (!side_not_exceed)
        {
        len_btm += inc_btm;
        }

        if (!btm_not_exceed)
        {
        len_side += inc_side;
        }

        } while (!(side_not_exceed && btm_not_exceed));

        if (iter == 1)
        {
        printf("%-8sBackward calculation of the simulated boundary area for stress %zu3 ...\n", "[INFO]", (5 - k));

        // len_btm -= h_default;
        // len_side -= h_default;

        do
        {
        if (side_not_exceed)
        {
        if (!(len_btm < h_default))
        {
        len_btm -= h_default;
        }
        }

        if (btm_not_exceed)
        {
        if (!(len_side < h_default))
        {
        len_side -= h_default;
        }
        }

        size_t nstep_side = (size_t)(len_side / h_default);
        size_t nstep_btm = (size_t)(len_btm / h_default);

        double h_side = len_side / (double)nstep_side;
        double h_btm = len_btm / (double)nstep_btm;

        size_t n_side = nstep_side;
        size_t n_btm = nstep_btm + 1;

        iter += 1;

        if (side_not_exceed)
        {
        for (size_t i = 0; i < n_side; i++)
        {
        double stress[6] = {0.0};

        tmk_vec3 side_pt = {0.0, 0.0, 0.0};

        side_pt.x = -1.0 * ((double)i * h_side) * cos(-1.0 * radian);
        side_pt.y = len_btm;
        side_pt.z = -1.0 * ((double)i * h_side) * sin(-1.0 * radian);

        tmk__coefficients_distributed_stress(stress, &side_pt, load, halfspace);
        cblas_dscal(6, load->intensity, stress, 1);
        tmk__coefficients_tensor_transformation_y(stress, degree);

        if (fabs(stress[k] / load->intensity) > threshold)
        {
        // printf("%zu/%zu, (%.2f, %.2f, %.2f), %.3e\n", i, n_side, side_pt.x, side_pt.y, side_pt.z, fabs(stress[k] / load->intensity));
        len_btm += h_default;
        side_not_exceed = 0;
        break;
        }
        }
        }

        if (btm_not_exceed)
        {
        for (size_t i = 0; i < n_btm; i++)
        {
        double stress[6] = {0.0};

        tmk_vec3 btm_pt = {0.0, 0.0, 0.0};

        btm_pt.x = -1.0 * len_side * cos(-1.0 * radian);
        btm_pt.y = (double)i * h_btm;
        btm_pt.z = -1.0 * len_side * sin(-1.0 * radian);

        tmk__coefficients_distributed_stress(stress, &btm_pt, load, halfspace);
        cblas_dscal(6, load->intensity, stress, 1);
        tmk__coefficients_tensor_transformation_y(stress, degree);

        if (fabs(stress[k] / load->intensity) > threshold)
        {
        // printf("%zu/%zu, (%.2f, %.2f, %.2f), %.3e\n", i, n_side, btm_pt.x, btm_pt.y, btm_pt.z, fabs(stress[k] / load->intensity));
        len_side += h_default;
        btm_not_exceed = 0;
        break;
        }
        }
        }

        printf("%-8sIteration %d, ", "[INFO]", iter);
        printf("x: %d, y: %d\n", side_not_exceed, btm_not_exceed);

        } while ((side_not_exceed && !(len_btm < h_default)) || (btm_not_exceed && !(len_side < h_default)));

        printf("%-8sBackward calculation, return %d, %d, %d, %d\n", "[INFO]", side_not_exceed, btm_not_exceed, len_side < h_default, len_btm < h_default);
        }
        else
        {
        printf("%-8sBackward calculation of the simulated boundary area for stress %zu3 ...\n", "[INFO]", (5 - k));

        do
        {
        if (side_not_exceed)
        {
        len_btm -= h_default;
        }

        if (btm_not_exceed)
        {
        len_side -= h_default;
        }

        size_t nstep_side = (size_t)(len_side / h_default);
        size_t nstep_btm = (size_t)(len_btm / h_default);

        double h_side = len_side / (double)nstep_side;
        double h_btm = len_btm / (double)nstep_btm;

        size_t n_side = nstep_side;
        size_t n_btm = nstep_btm + 1;

        iter += 1;

        if (side_not_exceed)
        {
        for (size_t i = 0; i < n_side; i++)
        {
        double stress[6] = {0.0};

        tmk_vec3 side_pt = {0.0, 0.0, 0.0};

        side_pt.x = -1.0 * ((double)i * h_side) * cos(-1.0 * radian);
        side_pt.y = len_btm;
        side_pt.z = -1.0 * ((double)i * h_side) * sin(-1.0 * radian);

        tmk__coefficients_distributed_stress(stress, &side_pt, load, halfspace);
        cblas_dscal(6, load->intensity, stress, 1);
        tmk__coefficients_tensor_transformation_y(stress, degree);

        if (fabs(stress[k] / load->intensity) > threshold)
        {
        side_not_exceed = 0;
        len_btm += h_default;
        break;
        }
        }
        }

        if (btm_not_exceed)
        {
        for (size_t i = 0; i < n_btm; i++)
        {
        double stress[6] = {0.0};

        tmk_vec3 btm_pt = {0.0, 0.0, 0.0};

        btm_pt.x = -1.0 * len_side * cos(-1.0 * radian);
        btm_pt.y = (double)i * h_btm;
        btm_pt.z = -1.0 * len_side * sin(-1.0 * radian);

        tmk__coefficients_distributed_stress(stress, &btm_pt, load, halfspace);
        cblas_dscal(6, load->intensity, stress, 1);
        tmk__coefficients_tensor_transformation_y(stress, degree);

        if (fabs(stress[k] / load->intensity) > threshold)
        {
        btm_not_exceed = 0;
        len_side += h_default;
        break;
        }
        }
        }

        printf("%-8sIteration %d, ", "[INFO]", iter);
        printf("x: %d, y: %d\n", side_not_exceed, btm_not_exceed);

        } while (side_not_exceed || btm_not_exceed);

        printf("%-8sBackward calculation, return %d, %d\n", "[INFO]", side_not_exceed, btm_not_exceed);
        }

        if (len_side > *len_side_final)
        {
        *len_side_final = len_side;
        }

        if (len_btm > *len_btm_final)
        {
        *len_btm_final = len_btm;
        }

        printf("\n%-8sTest for stress %zu3 finished after iteration %d\n", "[INFO]", (5 - k), iter);
        printf("%-8sFor stress %zu3, required domain of calculation is: x = %.1f [mm], y = %.1f [mm], viewed perpendicular to the area\n\n", "[INFO]", (5 - k), len_side, len_btm);
        }
        printf("%-8sFinal result: x = %.1f [mm], y = %.1f [mm]\n", "[INFO]", *len_side_final, *len_btm_final);
}

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

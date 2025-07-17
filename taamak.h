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
 */

typedef struct {double x, y, z;} tmk_vec3;
typedef struct {tmk_vec3 center; double half_width; double intensity;} tmk_load;
typedef struct {double youngs_modulus, poissons_ratio;} tmk_hs;

struct tmk_model_s;
typedef struct tmk_model_s *tmk_model;

TMKDEF tmk_model tmk_create(const int n_layers);
TMKDEF void tmk_set_load(tmk_model m, 
                         const tmk_vec3 center,
                         const double half_width, 
                         const double intensity);
TMKDEF void tmk_set_halfspace(tmk_model m, 
                              const int n_layers, 
                              const tmk_hs *hs);
TMKDEF void tmk_set_evaluation_points(tmk_model m, 
                                      const int npts, 
                                      const tmk_vec3 *pts);
TMKDEF void tmk_set_slope_angle(tmk_model m, const double degree);
TMKDEF void tmk_solve(tmk_model m);
TMKDEF void tmk_destroy(tmk_model *m);

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

#include <math.h>
#include <float.h>
#include <stdbool.h>

#define TMK_GEOM_EPSILON 0.001

struct tmk_model_s {
        int n_layers;
        tmk_hs *halfspace;
        tmk_load load;
        int npts;
        tmk_vec3 *pts;
        double slope_angle;
};

static bool tmk__util_assert_eq(double a, double b)
{
        double diff = fabs(a - b);
        double rel_eps = DBL_EPSILON;
        double abs_eps = DBL_MIN * 2.0;

        if (diff <= abs_eps) return true;

        double norm = fmin(fabs(a) + fabs(b), DBL_MAX);
        return (diff <= norm * rel_eps);
}

/**
 * fundamental solutions
 * TODO: descriptions
 */

static inline double tmk__coefficients_concentrated_radius(double a, double b, double c)
{
        return sqrt(a * a + b * b + c * c);
}

static inline double tmk__coefficients_concentrated_sxx_x(double x, double z, double c, 
                                                          double r1, double r2, double poissons_ratio)
{
        return x / 8.0 / M_PI / (poissons_ratio - 1.0) * ((2.0 * poissons_ratio - 1.0) / pow(r1, 3.0) + (1.0 - 2.0 * poissons_ratio) * (5.0 - 4.0 * poissons_ratio) / pow(r2, 3.0) - 3.0 * x * x / pow(r1, 5.0) - 3.0 * (3.0 - 4.0 * poissons_ratio) * x * x / pow(r2, 5.0) - 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / r2 / (r2 + z + c) / (r2 + z + c) * (3.0 - x * x * (3.0 * r2 + z + c) / r2 / r2 / (r2 + z + c)) + 6.0 * c / pow(r2, 5.0) * (3.0 * c - (3.0 - 2.0 * poissons_ratio) * (z + c) + 5.0 * x * x * z / r2 / r2));
}

static inline double tmk__coefficients_concentrated_syy_x(double x, double y, double z, double c, 
                                                          double r1, double r2, double poissons_ratio)
{
        return x / 8.0 / M_PI / (poissons_ratio - 1) * ((1.0 - 2.0 * poissons_ratio) / pow(r1, 3.0) + (1.0 - 2.0 * poissons_ratio) * (3.0 - 4.0 * poissons_ratio) / pow(r2, 3.0) - 3.0 * y * y / pow(r1, 5.0) - 3.0 * (3.0 - 4.0 * poissons_ratio) * y * y / pow(r2, 5.0) - 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / r2 / (r2 + z + c) / (r2 + z + c) * (1.0 - y * y * (3.0 * r2 + z + c) / r2 / r2 / (r2 + z + c)) + 6.0 * c / pow(r2, 5.0) * (c - (1 - 2.0 * poissons_ratio) * (z + c) + 5.0 * y * y * z / r2 / r2));
}

static inline double tmk__coefficients_concentrated_szz_x(double x, double z, double c, 
                                                          double r1, double r2, double poissons_ratio)
{
        return x / 8.0 / M_PI / (poissons_ratio - 1.0) * ((1.0 - 2.0 * poissons_ratio) / pow(r1, 3.0) - (1.0 - 2.0 * poissons_ratio) / pow(r2, 3.0) - 3.0 * (z - c) * (z - c) / pow(r1, 5.0) - 3.0 * (3.0 - 4.0 * poissons_ratio) * (z + c) * (z + c) / pow(r2, 5.0) + 6.0 * c / pow(r2, 5.0) * (c + (1.0 - 2.0 * poissons_ratio) * (z + c) + 5.0 * z * (z + c) * (z + c) / r2 / r2));
}

static inline double tmk__coefficients_concentrated_syz_x(double x, double y, double z, double c, 
                                                          double r1, double r2, double poissons_ratio)
{
        return x * y / 8.0 / M_PI / (poissons_ratio - 1.0) * (3.0 * (c - z) / pow(r1, 5.0) - 3.0 * (3.0 - 4.0 * poissons_ratio) * (z + c) / pow(r2, 5.0) + 6.0 * c / pow(r2, 5.0) * (1.0 - 2.0 * poissons_ratio + 5.0 * z * (z + c) / r2 / r2));
}

static inline double tmk__coefficients_concentrated_sxz_x(double x, double z, double c, 
                                                          double r1, double r2, double poissons_ratio)
{
        return 1.0 / 8.0 / M_PI / (poissons_ratio - 1.0) * ((2.0 * poissons_ratio - 1.0) * (z - c) / pow(r1, 3.0) + (1.0 - 2.0 * poissons_ratio) * (z - c) / pow(r2, 3.0) - 3.0 * x * x * (z - c) / pow(r1, 5.0) - 3.0 * (3.0 - 4.0 * poissons_ratio) * x * x * (z + c) / pow(r2, 5.0) - 6.0 * c / pow(r2, 5.0) * (z * (z + c) - (1.0 - 2.0 * poissons_ratio) * x * x - 5.0 * x * x * z * (z + c) / r2 / r2));
}

static inline double tmk__coefficients_concentrated_sxy_x(double x, double y, double z, double c, 
                                                          double r1, double r2, double poissons_ratio)
{
        return y / 8.0 / M_PI / (poissons_ratio - 1.0) * ((2.0 * poissons_ratio - 1.0) / pow(r1, 3.0) + (1.0 - 2.0 * poissons_ratio) / pow(r2, 3.0) - 3.0 * x * x / pow(r1, 5.0) - 3.0 * (3.0 - 4.0 * poissons_ratio) * x * x / pow(r2, 5.0) - 6.0 * c * z / pow(r2, 5.0) * (1.0 - 5.0 * x * x / r2 / r2) - 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / r2 / (r2 + z + c) / (r2 + z + c) * (1.0 - x * x * (3.0 * r2 + z + c) / r2 / r2 / (r2 + z + c)));
}

static inline double tmk__coefficients_concentrated_sxx_y(double x, double y, double z, double c, 
                                                          double r1, double r2, double poissons_ratio)
{
        return y / 8.0 / M_PI / (poissons_ratio - 1.0) * ((1.0 - 2.0 * poissons_ratio) / pow(r1, 3.0) + (1.0 - 2.0 * poissons_ratio) * (3.0 - 4.0 * poissons_ratio) / pow(r2, 3.0) - 3.0 * x * x / pow(r1, 5.0) - 3.0 * (3.0 - 4.0 * poissons_ratio) * x * x / pow(r2, 5.0) - 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / r2 / (r2 + z + c) / (r2 + z + c) * (1.0 - x * x * (3.0 * r2 + z + c) / r2 / r2 / (r2 + z + c)) + 6.0 * c / pow(r2, 5.0) * (c - (1.0 - 2.0 * poissons_ratio) * (z + c) + 5.0 * x * x * z / r2 / r2));
}

static inline double tmk__coefficients_concentrated_syy_y(double y, double z, double c, 
                                                          double r1, double r2, double poissons_ratio)
{
        return y / 8.0 / M_PI / (poissons_ratio - 1.0) * ((2.0 * poissons_ratio - 1.0) / pow(r1, 3.0) + (1.0 - 2.0 * poissons_ratio) * (5.0 - 4.0 * poissons_ratio) / pow(r2, 3.0) - 3.0 * y * y / pow(r1, 5.0) - 3.0 * (3.0 - 4.0 * poissons_ratio) * y * y / pow(r2, 5.0) - 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / r2 / (r2 + z + c) / (r2 + z + c) * (3.0 - y * y * (3.0 * r2 + z + c) / r2 / r2 / (r2 + z + c)) + 6.0 * c / pow(r2, 5.0) * (3.0 * c - (3.0 - 2.0 * poissons_ratio) * (z + c) + 5.0 * y * y * z / r2 / r2));
}

static inline double tmk__coefficients_concentrated_szz_y(double y, double z, double c, 
                                                          double r1, double r2, double poissons_ratio)
{
        return y / 8.0 / M_PI / (poissons_ratio - 1.0) * ((1.0 - 2.0 * poissons_ratio) / pow(r1, 3.0) - (1.0 - 2.0 * poissons_ratio) / pow(r2, 3.0) - 3.0 * (z - c) * (z - c) / pow(r1, 5.0) - 3.0 * (3.0 - 4.0 * poissons_ratio) * (z + c) * (z + c) / pow(r2, 5.0) + 6.0 * c / pow(r2, 5.0) * (c + (1.0 - 2.0 * poissons_ratio) * (z + c) + 5.0 * z * (z + c) * (z + c) / r2 / r2));
}

static inline double tmk__coefficients_concentrated_syz_y(double y, double z, double c, 
                                                          double r1, double r2, double poissons_ratio)
{
        return 1.0 / 8.0 / M_PI / (poissons_ratio - 1.0) * ((2.0 * poissons_ratio - 1.0) * (z - c) / pow(r1, 3.0) + (1.0 - 2.0 * poissons_ratio) * (z - c) / pow(r2, 3.0) - 3.0 * (3.0 - 4.0 * poissons_ratio) * y * y * (z + c) / pow(r2, 5.0) - 3.0 * y * y * (z - c) / pow(r1, 5.0) - 6.0 * c / pow(r2, 5.0) * (z * (z + c) - (1.0 - 2.0 * poissons_ratio) * y * y - 5.0 * y * y * z * (z + c) / r2 / r2));
}

static inline double tmk__coefficients_concentrated_sxz_y(double x, double y, double z, double c, 
                                                          double r1, double r2, double poissons_ratio)
{
        return y * x / 8.0 / M_PI / (poissons_ratio - 1.0) * (3.0 * (c - z) / pow(r1, 5.0) - 3.0 * (3.0 - 4.0 * poissons_ratio) * (z + c) / pow(r2, 5.0) + 6.0 * c / pow(r2, 5.0) * (1.0 - 2.0 * poissons_ratio + 5.0 * z * (z + c) / r2 / r2));
}

static inline double tmk__coefficients_concentrated_sxy_y(double x, double y, double z, double c, 
                                                          double r1, double r2, double poissons_ratio)
{
        return x / 8.0 / M_PI / (poissons_ratio - 1.0) * ((2.0 * poissons_ratio - 1.0) / pow(r1, 3.0) + (1.0 - 2.0 * poissons_ratio) / pow(r2, 3.0) - 3.0 * y * y / pow(r1, 5.0) - 3.0 * (3.0 - 4.0 * poissons_ratio) * y * y / pow(r2, 5.0) - 6.0 * c * z / pow(r2, 5.0) * (1.0 - 5.0 * y * y / r2 / r2) - 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / r2 / (r2 + z + c) / (r2 + z + c) * (1.0 - y * y * (3.0 * r2 + z + c) / r2 / r2 / (r2 + z + c)));
}

static inline double tmk__coefficients_concentrated_sxx_z(double x, double z, double c, 
                                                          double r1, double r2, double poissons_ratio)
{
        return 1.0 / 8.0 / M_PI / (poissons_ratio - 1.0) * ((1.0 - 2.0 * poissons_ratio) * (z - c) / pow(r1, 3.0) - 3.0 * x * x * (z - c) / pow(r1, 5.0) - 30.0 * c * x * x * z * (z + c) / pow(r2, 7.0) + (1.0 - 2.0 * poissons_ratio) * (3.0 * (z - c) - 4.0 * poissons_ratio * (z + c)) / pow(r2, 3.0) - (3.0 * x * x * (3.0 - 4.0 * poissons_ratio) * (z - c) - 6.0 * c * (z + c) * ((1.0 - 2.0 * poissons_ratio) * z - 2.0 * poissons_ratio * c)) / pow(r2, 5.0) - 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / r2 / (r2 + z + c) * (1.0 - x * x / r2 / (r2 + z + c) - x * x / r2 / r2));
}

static inline double tmk__coefficients_concentrated_syy_z(double y, double z, double c, 
                                                          double r1, double r2, double poissons_ratio)
{
        return 1.0 / 8.0 / M_PI / (poissons_ratio - 1.0) * ((1.0 - 2.0 * poissons_ratio) * (z - c) / pow(r1, 3.0) - 3.0 * y * y * (z - c) / pow(r1, 5.0) - 30.0 * c * y * y * z * (z + c) / pow(r2, 7.0) + (1.0 - 2.0 * poissons_ratio) * (3.0 * (z - c) - 4.0 * poissons_ratio * (z + c)) / pow(r2, 3.0) - (3.0 * (3.0 - 4.0 * poissons_ratio) * y * y * (z - c) - 6.0 * c * (z + c) * ((1.0 - 2.0 * poissons_ratio) * z - 2.0 * poissons_ratio * c)) / pow(r2, 5.0) - 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / r2 / (r2 + z + c) * (1.0 - y * y / r2 / (r2 + z + c) - y * y / r2 / r2));
}

static inline double tmk__coefficients_concentrated_szz_z(double z, double c, 
                                                          double r1, double r2, double poissons_ratio)
{
        return 1.0 / 8.0 / M_PI / (poissons_ratio - 1.0) * ((2.0 * poissons_ratio - 1.0) * (z - c) / pow(r1, 3.0) + (1.0 - 2.0 * poissons_ratio) * (z - c) / pow(r2, 3.0) - 3.0 * pow(z - c, 3.0) / pow(r1, 5.0) - 30.0 * c * z * pow(z + c, 3.0) / pow(r2, 7.0) - (3.0 * (3.0 - 4.0 * poissons_ratio) * z * (z + c) * (z + c) - 3.0 * c * (z + c) * (5.0 * z - c)) / pow(r2, 5.0));
}

static inline double tmk__coefficients_concentrated_syz_z(double y, double z, double c, 
                                                          double r1, double r2, double poissons_ratio)
{
        return y / 8.0 / M_PI / (poissons_ratio - 1.0) * ((2.0 * poissons_ratio - 1.0) / pow(r1, 3.0) + (1.0 - 2.0 * poissons_ratio) / pow(r2, 3.0) - (3.0 * (3.0 - 4.0 * poissons_ratio) * z * (z + c) - 3.0 * c * (3.0 * z + c)) / pow(r2, 5.0) - 3.0 * (z - c) * (z - c) / pow(r1, 5.0) - 30.0 * c * z * (z + c) * (z + c) / pow(r2, 7.0));
}

static inline double tmk__coefficients_concentrated_sxz_z(double x, double z, double c, 
                                                          double r1, double r2, double poissons_ratio)
{
        return x / 8.0 / M_PI / (poissons_ratio - 1.0) * ((2.0 * poissons_ratio - 1.0) / pow(r1, 3.0) + (1.0 - 2.0 * poissons_ratio) / pow(r2, 3.0) - (3.0 * (3.0 - 4.0 * poissons_ratio) * z * (z + c) - 3.0 * c * (3.0 * z + c)) / pow(r2, 5.0) - 3.0 * (z - c) * (z - c) / pow(r1, 5.0) - 30.0 * c * z * (z + c) * (z + c) / pow(r2, 7.0));
}

static inline double tmk__coefficients_concentrated_sxy_z(double x, double y, double z, double c, 
                                                          double r1, double r2, double poissons_ratio)
{
        return x * y / 8.0 / M_PI / (poissons_ratio - 1.0) * (3.0 * (c - z) / pow(r1, 5.0) - 3.0 * (3.0 - 4.0 * poissons_ratio) * (z - c) / pow(r2, 5.0) - 30.0 * c * z * (z + c) / pow(r2, 7.0) + 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / r2 / r2 / (r2 + z + c) * (1.0 / (r2 + z + c) + 1.0 / r2));
}

static inline double tmk__coefficients_concentrated_ux_x(double x, double z, double c, 
                                                         double r1, double r2, double youngs_modulus, double poissons_ratio)
{
        return (1.0 + poissons_ratio) / 8.0 / M_PI / youngs_modulus / (1.0 - poissons_ratio) * ((3.0 - 4.0 * poissons_ratio) / r1 + 1.0 / r2 + (3.0 - 4.0 * poissons_ratio) * x * x / pow(r2, 3.0) + 2.0 * c * z / pow(r2, 3.0) * (1.0 - 3.0 * x * x / r2 / r2) + x * x / pow(r1, 3.0) + 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / (r2 + z + c) * (1.0 - x * x / r2 / (r2 + z + c)));
}

static inline double tmk__coefficients_concentrated_uy_x(double x, double y, double z, double c, 
                                                         double r1, double r2, double youngs_modulus, double poissons_ratio)
{
        return (1.0 + poissons_ratio) * x * y / 8.0 / M_PI / youngs_modulus / (1.0 - poissons_ratio) * (1.0 / pow(r1, 3.0) + (3.0 - 4.0 * poissons_ratio) / pow(r2, 3.0) - 6.0 * c * z / pow(r2, 5.0) - 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / r2 / (r2 + z + c) / (r2 + z + c));
}

static inline double tmk__coefficients_concentrated_uz_x(double x, double z, double c, 
                                                         double r1, double r2, double youngs_modulus, double poissons_ratio)
{
        return (1.0 + poissons_ratio) * x / 8.0 / M_PI / youngs_modulus / (1.0 - poissons_ratio) * ((z - c) / pow(r1, 3.0) + (3.0 - 4.0 * poissons_ratio) * (z - c) / pow(r2, 3.0) - 6.0 * c * z * (z + c) / pow(r2, 5.0) + 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / r2 / (r2 + z + c));
}

static inline double tmk__coefficients_concentrated_ux_y(double x, double y, double z, double c, 
                                                         double r1, double r2, double youngs_modulus, double poissons_ratio)
{
        return (1.0 + poissons_ratio) * x * y / 8.0 / M_PI / youngs_modulus / (1.0 - poissons_ratio) * (1.0 / pow(r1, 3.0) + (3.0 - 4.0 * poissons_ratio) / pow(r2, 3.0) - 6.0 * c * z / pow(r2, 5.0) - 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / r2 / (r2 + z + c) / (r2 + z + c));
}

static inline double tmk__coefficients_concentrated_uy_y(double y, double z, double c, 
                                                         double r1, double r2, double youngs_modulus, double poissons_ratio)
{
        return (1.0 + poissons_ratio) / 8.0 / M_PI / youngs_modulus / (1.0 - poissons_ratio) * ((3.0 - 4.0 * poissons_ratio) / r1 + 1.0 / r2 + y * y / pow(r1, 3.0) + (3.0 - 4.0 * poissons_ratio) * y * y / pow(r2, 3.0) + 2.0 * c * z / pow(r2, 3.0) * (1.0 - 3.0 * y * y / r2 / r2) + 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / (r2 + z + c) * (1.0 - y * y / r2 / (r2 + z + c)));
}

static inline double tmk__coefficients_concentrated_uz_y(double y, double z, double c, 
                                                         double r1, double r2, double youngs_modulus, double poissons_ratio)
{
        return (1.0 + poissons_ratio) * y / 8.0 / M_PI / youngs_modulus / (1.0 - poissons_ratio) * ((z - c) / pow(r1, 3.0) + (3.0 - 4.0 * poissons_ratio) * (z - c) / pow(r2, 3.0) - 6.0 * c * z * (z + c) / pow(r2, 5.0) + 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / r2 / (r2 + z + c));
}

static inline double tmk__coefficients_concentrated_ux_z(double x, double z, double c, 
                                                         double r1, double r2, double youngs_modulus, double poissons_ratio)
{
        return (1.0 + poissons_ratio) * x / 8.0 / M_PI / youngs_modulus / (1.0 - poissons_ratio) * ((z - c) / pow(r1, 3) + (3.0 - 4.0 * poissons_ratio) * (z - c) / pow(r2, 3.0) - 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / r2 / (r2 + z + c) + 6.0 * c * z * (z + c) / pow(r2, 5.0));
}

static inline double tmk__coefficients_concentrated_uy_z(double y, double z, double c, 
                                                         double r1, double r2, double youngs_modulus, double poissons_ratio)
{
        return (1.0 + poissons_ratio) * y / 8.0 / M_PI / youngs_modulus / (1.0 - poissons_ratio) * ((z - c) / pow(r1, 3.0) + (3.0 - 4.0 * poissons_ratio) * (z - c) / pow(r2, 3.0) - 4.0 * (1.0 - poissons_ratio) * (1.0 - 2.0 * poissons_ratio) / r2 / (r2 + z + c) + 6.0 * c * z * (z + c) / pow(r2, 5.0));
}

static inline double tmk__coefficients_concentrated_uz_z(double z, double c, 
                                                         double r1, double r2, double youngs_modulus, double poissons_ratio)
{
        return (1.0 + poissons_ratio) / 8.0 / M_PI / youngs_modulus / (1.0 - poissons_ratio) * ((3.0 - 4.0 * poissons_ratio) / r1 + (8.0 * (1.0 - poissons_ratio) * (1.0 - poissons_ratio) - (3.0 - 4.0 * poissons_ratio)) / r2 + (z - c) * (z - c) / pow(r1, 3.0) + ((3.0 - 4.0 * poissons_ratio) * (z + c) * (z + c) - 2.0 * c * z) / pow(r2, 3.0) + 6.0 * c * z * (z + c) * (z + c) / pow(r2, 5.0));
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

        double r1 = tmk__coefficients_concentrated_radius(x, y, z - c);
        double r2 = tmk__coefficients_concentrated_radius(x, y, z + c);

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
                                                    tmk_hs *halfspace)
{
        double x = evaluation_point->x - force_point->x;
        double y = evaluation_point->y - force_point->y;
        double z = evaluation_point->z;

        double c = force_point->z;

        double r1 = tmk__coefficients_concentrated_radius(x, y, z - c);
        double r2 = tmk__coefficients_concentrated_radius(x, y, z + c);

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
                                                    tmk_hs *halfspace)
{
        double x = evaluation_point->x - force_point->x;
        double y = evaluation_point->y - force_point->y;
        double z = evaluation_point->z;

        double c = force_point->z;

        double r1 = tmk__coefficients_concentrated_radius(x, y, z - c);
        double r2 = tmk__coefficients_concentrated_radius(x, y, z + c);

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

        double r1 = tmk__coefficients_concentrated_radius(x, y, z - c);
        double r2 = tmk__coefficients_concentrated_radius(x, y, z + c);

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

        double r1 = tmk__coefficients_concentrated_radius(x, y, z - c);
        double r2 = tmk__coefficients_concentrated_radius(x, y, z + c);

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

        double r1 = tmk__coefficients_concentrated_radius(x, y, z - c);
        double r2 = tmk__coefficients_concentrated_radius(x, y, z + c);

        double youngs_modulus = halfspace->youngs_modulus;
        double poissons_ratio = halfspace->poissons_ratio;

        coefficients[0] = tmk__coefficients_concentrated_ux_z(x, z, c, r1, r2, youngs_modulus, poissons_ratio);
        coefficients[1] = tmk__coefficients_concentrated_uy_z(y, z, c, r1, r2, youngs_modulus, poissons_ratio);
        coefficients[2] = tmk__coefficients_concentrated_uz_z(z, c, r1, r2, youngs_modulus, poissons_ratio);
}

static void tmk__geometry_grid_uniform(tmk_vec3 *grid, 
                                       const int n_side,
                                       const double h, 
                                       const double degree,
                                       const tmk_load *load)
{
        double radian = degree * M_PI / 180.0;
        double s = tmk__util_assert_eq(radian, 0.0) ? sin(radian + DBL_EPSILON)
                                                    : sin(radian);
        double hx = h / s;

        for (int i = 0; i < n_side; ++i) {
                double x = hx * (double)(i + 1 - n_side);

                for (int j = 0; j < n_side; ++j) {
                        double y = h * (double)(2 * j + 1 - n_side) / 2.0;

                        int idx = i * n_side + j;

                        grid[idx].x = x * cos(-1.0 * radian);
                        if (tmk__util_assert_eq(grid[idx].x - load->center.x, 
                                                 load->half_width))
                                grid[idx].x += TMK_GEOM_EPSILON;

                        grid[idx].y = y;
                        if (tmk__util_assert_eq(grid[idx].y - load->center.y, 
                                                 load->half_width))
                                grid[idx].y += TMK_GEOM_EPSILON;

                        grid[idx].z = x * sin (-1.0 * radian);
                }
        }
}

TMKDEF tmk_model tmk_create(const int n_layers)
{
        tmk_model m;
        m->n_layers = n_layers;

        return m;
}

TMKDEF void tmk_set_load(tmk_model m, const tmk_vec3 center,
                         const double half_width, const double intensity)
{
        m->load = (tmk_load){
                .center = center,
                .half_width = half_width,
                .intensity = intensity
        };
}

TMKDEF void tmk_set_halfspace(tmk_model m, const int n_layers, const tmk_hs *hs)
{
        m->n_layers = n_layers;
        m->halfspace = (tmk_hs *)hs;
}

TMKDEF void tmk_set_evaluation_points(tmk_model m, const int npts, 
                                      const tmk_vec3 *pts)
{
        m->npts = npts;
        m->pts = (tmk_vec3 *)pts;
}

TMKDEF void tmk_set_slope_angle(tmk_model m, const double degree)
{
        m->slope_angle = degree;
}

TMKDEF void tmk_solve(tmk_model m)
{

}

TMKDEF void tmk_destroy(tmk_model *m)
{
        (void *)m;
}

#endif // TAAMAK_IMPLEMENTATION

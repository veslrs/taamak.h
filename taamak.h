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

// TMKDEF int add(const int a, const int b);

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

typedef struct {double x, y, z;} vec3_t;
typedef struct {vec3_t center; double half_width; double intensity} load_t;

static bool tmk__util_assert_eqd(double a, double b)
{
        double diff = fabs(a - b);
        double rel_eps = DBL_EPSILON;
        double abs_eps = DBL_MIN * 2.0;

        if (diff <= abs_eps) return true;

        double norm = fmin(fabs(a) + fabs(b), DBL_MAX);
        return (diff <= norm * rel_eps);
}

static void tmk__geometry_grid_uniform (vec3_t *grid, 
                                        const int n_side,
                                        const double h, 
                                        const double degree,
                                        const load_t *load)
{
        double radian = degree * M_PI / 180.0;
        double s = tmk__util_assert_eqd(radian, 0.0) ? sin (radian + DBL_EPSILON)
                                                     : sin (radian);
        double hx = h / s;

        for (int i = 0; i < n_side; ++i) {
                double x = hx * (double)(i + 1 - n_side);

                for (int j = 0; j < n_side; ++j) {
                        double y = h * (double)(2 * j + 1 - n_side) / 2.0;

                        int idx = i * n_side + j;

                        grid[idx].x = x * cos (-1.0 * radian);
                        if (tmk__util_assert_eqd(grid[idx].x - load->center.x, 
                                                 load->half_width))
                                grid[idx].x += 0.001;

                        grid[idx].y = y;
                        if (tmk__util_assert_eqd(grid[idx].y - load->center.y, 
                                                 load->half_width))
                                grid[idx].y += 0.001;

                        grid[idx].z = x * sin (-1.0 * radian);
                }
        }
}

#endif

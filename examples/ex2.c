/**
 * Copyright (c) 2025 Eyal Levenberg (eylev@dtu.dk)
 *                    Fang Zeyuan    (@veslrs)
 * Example #2: 
 *     Layered system, with known slope angle, layer compositions and material 
 *     properties.
 * 
 * TODO: descriptions
 */

#define TAAMAK_IMPLEMENTATION
#include "taamak.h"

int main(void)
{
        tmk_model m;
        m = tmk_create(2);

        // set loadings
        tmk_vec3 center = {300.0, 0.0, 0.0};
        double b = 150.0;
        double q = 0.5;
        tmk_set_load(m, center, b, q);

        // set layer properties
        tmk_hs hs[] = {
                // layer #1
                {
                        .youngs_modulus = 500.0, 
                        .poissons_ratio = 0.3,
                        .top_depth = 0.0,
                },
                // layer #2
                {
                        .youngs_modulus = 50.0,
                        .poissons_ratio = 0.4,
                        .top_depth = 600.0,
                },
        };
        tmk_set_halfspace(m, hs);

        // set evaluation points
        int npts = 4;
        tmk_vec3 pts[] = {
                {0.0, 0.0, 0.0},
                {300.0, 0.0, 0.0},
                {-150.0, 150.0, 450.0},
                {450.0, -450.0, 150.0},
        };
        tmk_set_evaluation_points(m, npts, pts);

        // set slope angle
        double alpha = 75.0;
        tmk_set_slope_angle(m, alpha);

        // solve
        tmk_solve(m);

        // results
        tmk_log(m);

        // cleanup
        tmk_destroy(m);

        // // solve
        // size_t mside = 21;
        // size_t nside = 17;
        // size_t n = nside * nside;
        // double *sol = calloc(3 * n, sizeof(double));
        // vec3_t *force_boundary = calloc(n, sizeof(vec3_t));
        // assert(sol != NULL);
        // assert(force_boundary != NULL);
        // solver_solve(sol, &halfspace, &load, degree, mside, nside, force_boundary);

        // // evaluation
        // FILE *export = fopen("./data/simple.csv", "w");
        // assert(export != NULL);
        // evaluation_evaluate(pts, npts, sol, force_boundary, n, &load, &halfspace, NULL);

        // // cleanup
        // free(sol);
        // free(force_boundary);
        // fclose(export);

        return 0;
}

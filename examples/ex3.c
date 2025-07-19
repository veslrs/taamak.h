/**
 * Copyright (c) 2025 Eyal Levenberg (eylev@dtu.dk)
 *                    Fang Zeyuan    (@veslrs)
 * Example #3: 
 *     Layered system, with known layer compositions and material properties. 
 *     Solve effective slope angle. 
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
        tmk_vec3 center = {150.0, 0.0, 0.0};
        double b = 133.0;
        double q = 0.707;
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
                        .youngs_modulus = 80.0,
                        .poissons_ratio = 0.4,
                        .top_depth = 750.0,
                },
        };
        tmk_set_halfspace(m, hs);

        // set FWD measurement points
        int npts = 10;
        tmk_vec3 pts[] = {
                {150.0, 0.0, 0.0},
                {250.0, 0.0, 0.0},
                {350.0, 0.0, 0.0},
                {450.0, 0.0, 0.0},
                {600.0, 0.0, 0.0},
                {750.0, 0.0, 0.0},
                {1050.0, 0.0, 0.0},
                {1350.0, 0.0, 0.0},
                {1650.0, 0.0, 0.0},
                {1950.0, 0.0, 0.0},
        };
        double vals[] = {
                0.582, 0.527, 0.347, 0.288, 0.245,
                0.217, 0.175, 0.142, 0.118, 0.099,
        };
        tmk_set_measurement(m, npts, pts, vals);

        // solve effective slope angle
        tmk_solve_slope_angle(m);

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

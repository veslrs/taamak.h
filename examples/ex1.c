/**
 * Copyright (c) 2025 Eyal Levenberg (eylev@dtu.dk)
 *                    Fang Zeyuan    (@veslrs)
 * TODO: descriptions
 */

#define TAAMAK_IMPLEMENTATION
#include "taamak.h"

int main(void)
{
        tmk_model m;
        m = tmk_create(1);

        // set loadings
        tmk_vec3 center = {150.0, 0.0, 0.0};
        double b = 133.0;
        double q = 0.707;
        tmk_set_load(m, center, b, q);

        // set layer properties
        int n_layers = 1;
        tmk_hs hs[] = {{.youngs_modulus = 500.0, .poissons_ratio = 0.3}};
        tmk_set_halfspace(m, n_layers, hs);

        // set evaluation points
        int npts = 1;
        tmk_vec3 pts[] = {{10.0, -20.0, 50.0}};
        tmk_set_evaluation_points(m, npts, pts);

        // set slope angle
        double alpha = 35.0;
        tmk_set_slope_angle(m, alpha);

        // solve
        tmk_solve(m);

        //cleanup
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

/**
 * Copyright (c) 2025 Eyal Levenberg (eylev@dtu.dk)
 *                    Fang Zeyuan    (@veslrs)
 * Example #1: 
 *     single layer, known slope angle and material properties
 * 
 * TODO: descriptions
 */

#define TAAMAK_IMPLEMENTATION
#include "taamak.h"

int main(void)
{
        tmk_model mdl;
        mdl = tmk_create(1);

        // set loadings
        tmk_vec3 center = {300.0, 0.0, 0.0};
        double b = 150.0;
        double q = 0.5;
        tmk_set_load(mdl, center, b, q);
        
        // set layer properties
        tmk_hs hs[] = {{.youngs_modulus = 500.0, .poissons_ratio = 0.3}};
        tmk_set_material_properties(mdl, hs);

        // set evaluation points
        int npts = 4;
        tmk_vec3 pts[] = {
                {0.0, 0.0, 0.0},
                {300.0, 0.0, 0.0},
                {-150.0, 150.0, 450.0},
                {450.0, -450.0, 150.0},
        };
        tmk_set_evaluation_points(mdl, npts, pts);
        
        // set slope angle
        double alpha = 75.0;
        tmk_set_slope_angle(mdl, alpha);

        // solve
        tmk_solve(mdl);

        // // // results
        // tmk_log(mdl);

        // cleanup
        tmk_destroy(mdl);

        return 0;
}

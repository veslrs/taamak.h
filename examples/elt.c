/*
        Copyright (c) 2025 Eyal Levenberg (eylev@dtu.dk)
                           Fang Zeyuan    (@veslrs)
        Example: Three-layered system, with known slope angle, layer compositions and 
        material properties. Evaluate surface vertical deflection.
*/
#define TAAMAK_IMPLEMENTATION
#include "../taamak.h"

int main(void)
{
        tmk_model mdl;
        tmk_init(&mdl, 3);

        // set loadings
        tmk_vec3 center = {150.0, 0.0, 0.0};
        double b = 133.0;
        double q = 0.707;
        double x0 = 150.0;
        tmk_set_load(&mdl, center, b, q);

        // set layer compositions
        tmk_set_composition(&mdl, (double[]){0.0, 200.0, 700.0});

        // set layer properties
        tmk_hs hs[] = {
                // layer #1
                {
                        .youngs_modulus = 4000.0,
                        .poissons_ratio = 0.3,
                },
                // layer #2
                {
                        .youngs_modulus = 250.0,
                        .poissons_ratio = 0.35,
                },
                // layer #3
                {
                        .youngs_modulus = 70.0,
                        .poissons_ratio = 0.4,
                },
        };
        tmk_set_material_properties(&mdl, hs);

        // set evaluation points

        // parallel
        // tmk_vec3 p1 = {150.0, 0.0, 0.0};
        // tmk_vec3 p2 = {150.0, 1800.0, 0.0};
        
        // transverse
        // tmk_vec3 p1 = {150.0, 0.0, 0.0};
        // tmk_vec3 p2 = {1950.0, 0.0, 0.0};
        tmk_vec3 pts[] = {
                {x0         , 0.0, 0.0},
                {x0 + 200.0 , 0.0, 0.0},
                {x0 + 300.0 , 0.0, 0.0},
                {x0 + 450.0 , 0.0, 0.0},
                {x0 + 600.0 , 0.0, 0.0},
                {x0 + 900.0 , 0.0, 0.0},
                {x0 + 1200.0, 0.0, 0.0},
                {x0 + 1500.0, 0.0, 0.0},
                {x0 + 1800.0, 0.0, 0.0},
        };

        // tmk_vec3 pts[300] = {0};
        // tmk_linspace(pts, &p1, &p2, 300);
        // tmk_vec3 pts[] = {
        //         {0.0, 0.0, 0.0},
        //         {300.0, 0.0, 0.0},
        //         {600.0, 150.0, 0.0},
        //         {1200.0, 300.0, 0.0},
        // };
        tmk_set_evaluation_points(&mdl, 9, pts);

        // set slope angle
        tmk_set_slope_angle(&mdl, 0.0);

        // solve
        tmk_solve(&mdl);

        // done!
        tmk_checkout(&mdl);
}

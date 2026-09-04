/*
    Copyright (c) 2025 Eyal Levenberg (eylev@dtu.dk)
                       Fang Zeyuan    (@veslrs)
    Example: Single layer system, with known slope angle and 
    material properties. Evaluate surface vertical deflection.
*/
#define TAAMAK_IMPLEMENTATION
#include "../taamak.h"

int main(void)
{
    tmk_model mdl;
    tmk_init(&mdl, 1);

    // set loadings
    tmk_vec3 center = {150.0, 0.0, 0.0};
    double b = 133.0;
    double q = 0.707;
    double x0 = 150.0;
    tmk_set_load(&mdl, center, b, q);

    // set layer compositions
    // tmk_set_composition(&mdl, (double[]){0.0});

    // set layer properties
    tmk_hs hs[] = {
        {
            .youngs_modulus = 4000.0,
            .poissons_ratio = 0.3,
        },
    };
    tmk_set_material_properties(&mdl, hs);

    // set evaluation points
    tmk_vec3 pts[] = {
        {x0  , 0.0, 0.0},
    };
    tmk_set_evaluation_points(&mdl, 1, pts);

    // set slope angle
    tmk_set_slope_angle(&mdl, 90.0);

    // solve
    tmk_solve(&mdl);

    // done!
    tmk_checkout(&mdl);
}

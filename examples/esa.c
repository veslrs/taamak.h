/*
        Copyright (c) 2025 Eyal Levenberg (eylev@dtu.dk)
                           Fang Zeyuan    (@veslrs)
        Example: Layered system, with known layer compositions and material properties. 
        Solve effective slope angle. 
*/
#define TAAMAK_IMPLEMENTATION
#include "../taamak.h"

int main(void)
{
        tmk_model mdl;
        tmk_init(&mdl, 2);

        // set loadings
        tmk_vec3 center = {150.0, 0.0, 0.0};
        double b = 133.0;
        double q = 0.707;
        tmk_set_load(&mdl, center, b, q);

        // set layer compositions
        tmk_set_composition(&mdl, (double[]){0.0, 750.0});

        // set layer properties
        tmk_hs hs[] = {
                // layer #1
                {
                        .youngs_modulus = 500.0, 
                        .poissons_ratio = 0.3,
                },
                // layer #2
                {
                        .youngs_modulus = 80.0,
                        .poissons_ratio = 0.4,
                },
        };
        tmk_set_material_properties(&mdl, hs);

        // set FWD measurement points
        int npts = 10;
        tmk_vec3 pts[] = {
                {150.0,  0.0, 0.0},
                {250.0,  0.0, 0.0},
                {350.0,  0.0, 0.0},
                {450.0,  0.0, 0.0},
                {600.0,  0.0, 0.0},
                {750.0,  0.0, 0.0},
                {1050.0, 0.0, 0.0},
                {1350.0, 0.0, 0.0},
                {1650.0, 0.0, 0.0},
                {1950.0, 0.0, 0.0},
        };
        double vals[] = {
                0.659, 0.593, 0.404, 0.338, 0.287, 
                0.253, 0.201, 0.163, 0.134, 0.113,
        };
        tmk_set_measurement(&mdl, npts, pts, vals);

        // solve effective slope angle
        tmk_solve(&mdl);

        // done!
        tmk_checkout(&mdl);
}

#include "hmpp_fun.h"



#ifdef USING_LOCAL_DIRECTIVES
#pragma hmpp updateCentroidVH codelet, target=CUDA args[*].transfer=atcall
#endif
void _update_centroids_of_velocities_and_height(
        int N_c,
        int N_b,
        double w_C[N_c], // stage_centroid_values
        double uh_C[N_c],// xmomentum_centroid_values
        double vh_C[N_c],// ymomentum_centroid_values
        double h_C[N_c], // height_centroid_values
        double z_C[N_c], // elevation_centroid_values
        double u_C[N_c], // xvelocity_centroid_values
        double v_C[N_c], // yvelocity_centroid_values

        double w_B[N_b], // stage_boundary_values
        double uh_B[N_b],// xmomentum_boundary_values
        double vh_B[N_b],// ymomentum_boundary_values
        double h_B[N_b], // height_boundary_values
        double z_B[N_b], // elevation_boundary_values
        double u_B[N_b], // xvelocity_boundary_values
        double v_B[N_b] // yvelocity_boundary_values
        )
{
    int k;

    for(k=0; k<N_c; k++)
    {
        double H0 = 1.0e-8;
        double factor;

        if ( k < N_c )
        {
            h_C[k] = w_C[k] - z_C[k];
            if (h_C[k] < 0)
                h_C[k] = 0;


            factor = h_C[k] / (h_C[k]*h_C[k] + H0);
            u_C[k] = uh_C[k]*factor;
            v_C[k] = vh_C[k]*factor;
        }


        if (k < N_b)
        {
            h_B[k] = w_B[k] - z_B[k];
            if (h_B[k] < 0)
                h_B[k] = 0;


            factor = h_B[k] / (h_B[k]*h_B[k] + H0);
            u_B[k] = uh_B[k]*factor;
            v_B[k] = vh_B[k]*factor;
        }
    }
}

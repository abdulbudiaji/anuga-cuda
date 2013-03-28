__global__ void update_centroids_of_velocities_and_height(
        int N_c,
        int N_b,
        double * w_C, // stage_centroid_values
        double * uh_C,// xmomentum_centroid_values
        double * vh_C,// ymomentum_centroid_values
        double * h_C, // height_centroid_values
        double * z_C, // elevation_centroid_values
        double * u_C, // xvelocity_centroid_values
        double * v_C, // yvelocity_centroid_values

        double * w_B, // stage_boundary_values
        double * uh_B,// xmomentum_boundary_values
        double * vh_B,// ymomentum_boundary_values
        double * h_B, // height_boundary_values
        double * z_B, // elevation_boundary_values
        double * u_B, // xvelocity_boundary_values
        double * v_B // yvelocity_boundary_values
        )
{
    const int k = 
            threadIdx.x+threadIdx.y*blockDim.x+
            (blockIdx.x+blockIdx.y*gridDim.x)*blockDim.x*blockDim.y; 

    double H0 = 1.0e-8;
    double factor;

    if ( k >=N_c )
        return;

    h_C[k] = w_C[k] - z_C[k];
    if (h_C[k] < 0)
        h_C[k] = 0;


    factor = h_C[k] / (h_C[k]*h_C[k] + H0);
    u_C[k] = uh_C[k]*factor;
    v_C[k] = vh_C[k]*factor;



    if (k >= N_b)
        return;

    h_B[k] = w_B[k] - z_B[k];
    if (h_B[k] < 0)
        h_B[k] = 0;


    factor = h_B[k] / (h_B[k]*h_B[k] + H0);
    u_B[k] = uh_B[k]*factor;
    v_B[k] = vh_B[k]*factor;
}

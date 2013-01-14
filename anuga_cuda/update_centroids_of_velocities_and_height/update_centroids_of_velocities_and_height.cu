__global__ void update_centroids_of_velocities_and_height(
        int N,
        double * w_C,
        double * z_C,
        double * uh_C,
        double * vh_C,
        double * u_C,
        double * v_C,
        double * h_C,

        double * w_B,
        double * z_B,
        double * uh_B,
        double * vh_B,
        double * u_B,
        double * v_B,
        double * h_B)
{
    const int k = 
            threadIdx.x+threadIdx.y*blockDim.x+
            (blockIdx.x+blockIdx.y*gridDim.x)*blockDim.x*blockDim.y; 

    double H0 = 1.0e-8;
    double factor;

    if (k >=N)
        return;

    h_C[k] = w_C[k] - z_C[k];

    if (h_C[k] < 0)
        h_C[k] = 0;

    h_B[k] = w_B[k] - z_B[k];
    if (h_B[k] < 0)
        h_B[k] = 0;

    factor = h_C[k] / (h_C[k]*h_C[k] + H0);
    u_C[k] = uh_C[k]*factor;
    v_C[k] = vh_C[k]*factor;


    factor = h_B[k] / (h_B[k]*h_B[k] + H0);
    u_B[k] = uh_B[k]*factor;
    v_B[k] = vh_B[k]*factor;

}

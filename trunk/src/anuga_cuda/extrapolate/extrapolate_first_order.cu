//#define REARRANGED_DOMAIN

__global__ void extrapolate_first_order(
        int N,
        double * centroid_values,
        double * edge_values,
        double * vertex_values)
{
    const int k = 
            threadIdx.x+threadIdx.y*blockDim.x+
            (blockIdx.x+blockIdx.y*gridDim.x)*blockDim.x*blockDim.y;

    if (k >= N)
        return;
#ifndef REARRANGED_DOMAIN
    int k3 = k*3;
    edge_values[k3] = centroid_values[k];
    edge_values[k3+ 1] = centroid_values[k];
    edge_values[k3+ 2] = centroid_values[k];

    vertex_values[k3] = centroid_values[k];
    vertex_values[k3+ 1] = centroid_values[k];
    vertex_values[k3+ 2] = centroid_values[k];
#else
    edge_values[k] = centroid_values[k];
    edge_values[k + N] = centroid_values[k];
    edge_values[k + 2*N] = centroid_values[k];

    vertex_values[k] = centroid_values[k];
    vertex_values[k + N] = centroid_values[k];
    vertex_values[k + 2*N] = centroid_values[k];
#endif
}

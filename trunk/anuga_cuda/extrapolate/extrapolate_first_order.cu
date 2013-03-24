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

    edge_values[k*3] = centroid_values[k];
    edge_values[k*3+ 1] = centroid_values[k];
    edge_values[k*3+ 2] = centroid_values[k];

    vertex_values[k*3] = centroid_values[k];
    vertex_values[k*3+ 1] = centroid_values[k];
    vertex_values[k*3+ 2] = centroid_values[k];
}

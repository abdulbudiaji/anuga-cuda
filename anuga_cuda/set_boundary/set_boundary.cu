__global__ void set_boundary_values_from_edges(
        int N,
        int * vol_id,
        int * edge_id,
        double * boundary_values,
        double * edge_values)
{
    const int k = 
            threadIdx.x+threadIdx.y*blockDim.x+
            (blockIdx.x+blockIdx.y*gridDim.x)*blockDim.x*blockDim.y;
    
    if ( k >= N )
        return;

    int id = 3* vol_id[k] + edge_id[k];

    boundary_values[id] = edge_values[id];
}

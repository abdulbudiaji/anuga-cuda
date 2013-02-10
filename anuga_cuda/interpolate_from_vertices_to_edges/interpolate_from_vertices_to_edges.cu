__global__ void _interpolate_from_vertices_to_edges(
        int N,
        double* vertex_values,
        double* edge_values) 
{
    const int k = 
        threadIdx.x+threadIdx.y*blockDim.x+
        (blockIdx.x+blockIdx.y*gridDim.x)*blockDim.x*blockDim.y;

    if ( k >= N )
        return;


    int k3= k*3;
    double q0, q1, q2;

    //for (k=0; k<N; k++) {

    q0 = vertex_values[k3 + 0];
    q1 = vertex_values[k3 + 1];
    q2 = vertex_values[k3 + 2];

    edge_values[k3 + 0] = 0.5*(q1+q2);
    edge_values[k3 + 1] = 0.5*(q0+q2);
    edge_values[k3 + 2] = 0.5*(q0+q1);
    //}
}

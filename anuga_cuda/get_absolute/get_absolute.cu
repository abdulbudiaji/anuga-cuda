__global__ void get_absolute(
    int N,
    double xllcorner,
    double yllcorner,
    double * points)
{
    const int k = 
            threadIdx.x+threadIdx.y*blockDim.x+
            (blockIdx.x+blockIdx.y*gridDim.x)*blockDim.x*blockDim.y;
    if (k >= N )
        return;
    points[k*2] += xllcorner;
    points[k*2 + 1] += yllcorner;
}

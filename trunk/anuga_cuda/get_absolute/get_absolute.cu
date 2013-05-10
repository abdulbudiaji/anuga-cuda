//#define REARRANGED_DOMAIN

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
#ifndef REARRANGED_DOMAIN
    int k2 = k*2;
    points[k2] += xllcorner;
    points[k2 + 1] += yllcorner;
#else
    points[k] += xllcorner;
    points[k + N] += yllcorner;
#endif
}

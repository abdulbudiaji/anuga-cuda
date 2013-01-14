__global__ void _saxpy_centroid_values(
        int N,
        double a,
        double b,
        double * centroid_values,
        double * centroid_backup_values)
{
    const int k = 
            threadIdx.x+threadIdx.y*blockDim.x+
            (blockIdx.x+blockIdx.y*gridDim.x)*blockDim.x*blockDim.y;
    centroid_values[k]= a*centroid_values[k] + b*centroid_backup_values[k];
}

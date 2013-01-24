// from quantity_ext.c
__global__ void update(
        int N,
        double timestep,
        double * centroid_values,
        double * explicit_update,
        double * semi_implicit_update)
{
    const int k = 
            threadIdx.x+threadIdx.y*blockDim.x+
            (blockIdx.x+blockIdx.y*gridDim.x)*blockDim.x*blockDim.y;

    double denominator, x;

    if (k >= N)
        return;

    // Divide semi_implicit update by conserved quantity
    //for (k=0; k<N; k++) {
    x = centroid_values[k];
    if (x == 0.0) {
        semi_implicit_update[k] = 0.0;
    } else {
        semi_implicit_update[k] /= x;
    }
    //}


    // Explicit updates
    //for (k=0; k<N; k++) {
    centroid_values[k] += timestep*explicit_update[k];
    //}



    // Semi implicit updates
    //for (k=0; k<N; k++) {
    denominator = 1.0 - timestep*semi_implicit_update[k];
    if (denominator <= 0.0) {
        return;
    } else {
        //Update conserved_quantities from semi implicit updates
        centroid_values[k] /= denominator;
    }
    //}



    // Reset semi_implicit_update here ready for next time step
    //memset(semi_implicit_update, 0, N*sizeof(double));
}


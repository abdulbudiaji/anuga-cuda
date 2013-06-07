__global__ void _protect_sw(
        int N,
        double minimum_allowed_height,
        double maximum_allowed_speed,
        double epsilon,
        double* wc,
        double* zc,
        double* xmomc,
        double* ymomc) 
{
    const int k = threadIdx.x+threadIdx.y*blockDim.x+
        (blockIdx.x+blockIdx.y*gridDim.x)*blockDim.x*blockDim.y;
    double hc;
    double u, v, reduced_speed;

    if (k >= N)
        return;

    // Protect against initesimal and negative heights
    if (maximum_allowed_speed < epsilon) {
        //for (k = 0; k < N; k++) {
        hc = wc[k] - zc[k];
        if (hc < minimum_allowed_height) {
            // Set momentum to zero and ensure h is non negative
            xmomc[k] = 0.0;
            ymomc[k] = 0.0;
            if (hc <= 0.0) wc[k] = zc[k];
        }

    } else {
        // Protect against initesimal and negative heights
        //for (k = 0; k < N; k++) {
        hc = wc[k] - zc[k];
        if (hc < minimum_allowed_height) {
            //New code: Adjust momentum to guarantee speeds are physical
            //          ensure h is non negative
            if (hc <= 0.0) {
                wc[k] = zc[k];
                xmomc[k] = 0.0;
                ymomc[k] = 0.0;
            } else {
                //Reduce excessive speeds derived from division by small hc
                //FIXME (Ole): This may be unnecessary with new slope limiters
                //in effect.

                u = xmomc[k] / hc;
                if (fabs(u) > maximum_allowed_speed) {
                    reduced_speed = maximum_allowed_speed * u / fabs(u);
                    xmomc[k] = reduced_speed * hc;
                }

                v = ymomc[k] / hc;
                if (fabs(v) > maximum_allowed_speed) {
                    reduced_speed = maximum_allowed_speed * v / fabs(v);
                    ymomc[k] = reduced_speed * hc;
                }
            }
        }
    }
}


// Protect against the water elevation falling below the triangle bed
__global__ void  _protect_swb2(
        int N,
        double minimum_allowed_height,
        double maximum_allowed_speed,
        double epsilon,
        double* wc,
        double* wv,
        double* zc,
        double* zv,
        double* xmomc,
        double* ymomc,
        double* areas) 
{
    const long k = threadIdx.x+threadIdx.y*blockDim.x+
        (blockIdx.x+blockIdx.y*gridDim.x)*blockDim.x*blockDim.y;
    
    double hc, bmin, bmax;
    double mass_error=0.; 
    // This acts like minimum_allowed height, but scales with the vertical
    // distance between the bed_centroid_value and the max bed_edge_value of
    // every triangle.
    double minimum_relative_height=0.1; 
    int mass_added=0;
    
    if (k >= N)
        return;

    // Protect against inifintesimal and negative heights  
    //if (maximum_allowed_speed < epsilon) {
    //for (k=0; k<N; k++) {
    hc = wc[k] - zc[k];
    // Definine the maximum bed edge value on triangle k.
    bmax = 0.5*max((zv[3*k]+zv[3*k+1]),max((zv[3*k+1]+zv[3*k+2]),(zv[3*k+2]+zv[3*k])));

    if (hc < max(minimum_relative_height*(bmax-zc[k]), minimum_allowed_height)) {

        // Set momentum to zero and ensure h is non negative
        // NOTE: THIS IS IMPORTANT -- WE ARE SETTING MOMENTUM TO ZERO
        //if(hc<=epsilon){
        xmomc[k] = 0.0;
        ymomc[k] = 0.0;
        //}

        if (hc <= 0.0){
            // Definine the minimum bed edge value on triangle k.
            // Setting = minimum edge value can lead to mass conservation problems
            //bmin = 0.5*min((zv[3*k]+zv[3*k+1]),min((zv[3*k+1]+zv[3*k+2]),(zv[3*k+2]+zv[3*k])));
            //bmin =0.5*bmin + 0.5*min(zv[3*k],min(zv[3*k+1],zv[3*k+2]));
            // Setting = minimum vertex value seems better, but might tend to be less smooth 
            bmin =min(zv[3*k],min(zv[3*k+1],zv[3*k+2])) -minimum_allowed_height;
            //bmin=zc[k]-minimum_allowed_height;
            // Minimum allowed stage = bmin
            // WARNING: ADDING MASS if wc[k]<bmin
            if(wc[k]<bmin){
                mass_error+=(bmin-wc[k])*areas[k];
                mass_added=1; //Flag to warn of added mass                
                //printf("Adding mass to dry cell %d %f %f %f %f %f \n", k, zv[3*k], zv[3*k+1], zv[3*k+2], wc[k]- bmin, mass_error);

                wc[k] = max(wc[k], bmin); 


                // Set vertex values as well. Seems that this shouldn't be
                // needed. However from memory this is important at the first
                // time step, for 'dry' areas where the designated stage is
                // less than the bed centroid value
                wv[3*k] = max(wv[3*k], bmin);
                wv[3*k+1] = max(wv[3*k+1], bmin);
                wv[3*k+2] = max(wv[3*k+2], bmin);
            }
        }
    }
    //}

    //return mass_error;
}

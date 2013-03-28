#define REARRANGED_DOMAIN

__global__ void _manning_friction_flat(
        int N,
        double g, 
        double eps, 

        double* w, 
        double* zv,
        double* uh,
        double* vh,
        double* eta, 
        double* xmom, 
        double* ymom) 
{
    const int k = 
        threadIdx.x+threadIdx.y*blockDim.x+
        (blockIdx.x+blockIdx.y*gridDim.x)*blockDim.x*blockDim.y;
    if ( k >= N )
        return;
#ifndef REARRANGED_DOMAIN
    int k3;
#endif
    double S, h, z, z0, z1, z2;

    //for (k = 0; k < N; k++) {
    if (eta[k] > eps) {
#ifndef REARRANGED_DOMAIN
        k3 = 3 * k;
        z0 = zv[k3 + 0];
        z1 = zv[k3 + 1];
        z2 = zv[k3 + 2];
#else
        z0 = zv[k];
        z1 = zv[k + N];
        z2 = zv[k + 2*N];
#endif

        z = (z0 + z1 + z2) / 3.0;
        h = w[k] - z;
        if (h >= eps) {
            S = -g *eta[k] * eta[k] * sqrt((uh[k] * uh[k] + vh[k] * vh[k]));
            S /= pow(h, 7.0 / 3); //Expensive (on Ole's home computer)
            //seems to save about 15% over manning_friction
            //S /= exp((7.0/3.0)*log(h));      
            //FIXME: Could use a Taylor expansion
            //S /= h*h*(1 + h/3.0 - h*h/9.0); 

            //Update momentum
            xmom[k] += S * uh[k];
            ymom[k] += S * vh[k];
        }
    }
    //}
}

__global__ void _manning_friction_sloped(
        int N,
        double g, 
        double eps, 

        double* x, 
        double* w, 
        double* zv,
        double* uh, 
        double* vh,
        double* eta, 
        double* xmom_update, 
        double* ymom_update) 
{
    const int k = 
        threadIdx.x+threadIdx.y*blockDim.x+
        (blockIdx.x+blockIdx.y*gridDim.x)*blockDim.x*blockDim.y;
    if ( k >= N )
        return;
#ifndef REARRANGED_DOMAIN
    int k3, k6;
#endif
    double S, h, z, z0, z1, z2, zs, zx, zy;
    double x0, y0, x1, y1, x2, y2;
    double det;

    //for (k = 0; k < N; k++) {
    if (eta[k] > eps) {
#ifndef REARRANGED_DOMAIN
        k3 = 3 * k;
        z0 = zv[k3 + 0];
        z1 = zv[k3 + 1];
        z2 = zv[k3 + 2];

        // Compute bed slope
        k6 = 6 * k; // base index

        x0 = x[k6 + 0];
        y0 = x[k6 + 1];
        x1 = x[k6 + 2];
        y1 = x[k6 + 3];
        x2 = x[k6 + 4];
        y2 = x[k6 + 5];
#else
        z0 = zv[k];
        z1 = zv[k + N];
        z2 = zv[k + 2*N];

        // Compute bed slope
        x0 = x[k];
        y0 = x[k + N];
        x1 = x[k + 2*N];
        y1 = x[k + 3*N];
        x2 = x[k + 4*N];
        y2 = x[k + 5*N];
#endif

        //_gradient(x0, y0, x1, y1, x2, y2, z0, z1, z2, &zx, &zy);
        det = (y2-y0)*(x1-x0) - (y1-y0)*(x2-x0);

        zx = (y2-y0)*(z1-z0) - (y1-y0)*(z2-z0);
        zx /= det;

        zy = (x1-x0)*(z2-z0) - (x2-x0)*(z1-z0);
        zy /= det;


        zs = sqrt(1.0 + zx * zx + zy * zy);
        z = (z0 + z1 + z2) / 3.0;
        h = w[k] - z;
        if (h >= eps) {
            S =-g *eta[k] *eta[k] *zs *sqrt((uh[k] *uh[k] + vh[k] * vh[k]));
            S /= pow(h, 7.0 / 3); //Expensive (on Ole's home computer)
            //S /= exp((7.0/3.0)*log(h));      
            //seems to save about 15% over manning_friction
            //S /= h*h*(1 + h/3.0 - h*h/9.0); 
            //FIXME: Could use a Taylor expansion

            //Update momentum
            xmom_update[k] += S * uh[k];
            ymom_update[k] += S * vh[k];
        }
    }
    //}
}

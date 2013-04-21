
#include "hmpp_fun.h"
#ifdef USING_SEPARATE_KERNELS


void  limit_vertices_by_all_neighbours(
        int N, 
        double beta,
        double* centroid_values,
        double* vertex_values,
        double* edge_values,
        long*   neighbours,
        double* x_gradient,
        double* y_gradient) 
{
    const int k = 
        threadIdx.x+threadIdx.y*blockDim.x+
        (blockIdx.x+blockIdx.y*gridDim.x)*blockDim.x*blockDim.y;

    if (k >= N)
        return;

    int i, k3;
    long n;
    double qmin, qmax, qn, qc;
    double dq, dqa[3], phi, r;

    //for (k=0; k<N; k++){
    k3 = 3*k;

    qc = centroid_values[k];
    qmin = qc;
    qmax = qc;

    for (i=0; i<3; i++) {
        n = neighbours[k3+i];
        if (n >= 0) {
            qn = centroid_values[n]; //Neighbour's centroid value

            qmin = min(qmin, qn);
            qmax = max(qmax, qn);
        }
    }

    phi = 1.0;
    for (i=0; i<3; i++) {    
        r = 1.0;

        dq = vertex_values[k3+i] - qc;    //Delta between vertex and centroid values
        dqa[i] = dq;                      //Save dq for use in updating vertex values

        if (dq > 0.0) r = (qmax - qc)/dq;
        if (dq < 0.0) r = (qmin - qc)/dq;      


        phi = min( min(r*beta, 1.0), phi);    
    }

    //Update gradient, vertex and edge values using phi limiter
    x_gradient[k] = x_gradient[k]*phi;
    y_gradient[k] = y_gradient[k]*phi;

    vertex_values[k3+0] = qc + phi*dqa[0];
    vertex_values[k3+1] = qc + phi*dqa[1];
    vertex_values[k3+2] = qc + phi*dqa[2];

    edge_values[k3+0] = 0.5*(vertex_values[k3+1] + vertex_values[k3+2]);
    edge_values[k3+1] = 0.5*(vertex_values[k3+2] + vertex_values[k3+0]);
    edge_values[k3+2] = 0.5*(vertex_values[k3+0] + vertex_values[k3+1]);

    //}

    //return 0;
}



void limit_edges_by_all_neighbours(
        int N, 
        double beta,
        double* centroid_values,
        double* vertex_values,
        double* edge_values,
        long*   neighbours,
        double* x_gradient,
        double* y_gradient) 
{
    const int k = 
        threadIdx.x+threadIdx.y*blockDim.x+
        (blockIdx.x+blockIdx.y*gridDim.x)*blockDim.x*blockDim.y;

    if (k >= N)
        return;

    int i, k2, k3, k6;
    long n;
    double qmin, qmax, qn, qc, sign;
    double dq, dqa[3], phi, r;

    //for (k=0; k<N; k++){
    k6 = 6*k;
    k3 = 3*k;
    k2 = 2*k;

    qc = centroid_values[k];
    qmin = qc;
    qmax = qc;

    for (i=0; i<3; i++) {
        n = neighbours[k3+i];
        if (n >= 0) {
            qn = centroid_values[n]; //Neighbour's centroid value

            qmin = min(qmin, qn);
            qmax = max(qmax, qn);
        }
    }

    sign = 0.0;
    if (qmin > 0.0) {
        sign = 1.0;
    } else if (qmax < 0) {
        sign = -1.0;
    }

    phi = 1.0;
    for (i=0; i<3; i++) {  
        dq = edge_values[k3+i] - qc;      //Delta between edge and centroid values
        dqa[i] = dq;                      //Save dq for use in updating vertex values  


        // Just limit non boundary edges so that we can reconstruct a linear function
        // FIXME Problem with stability on edges 
        //if (neighbours[k3+i] >= 0) {
        r = 1.0;

        if (dq > 0.0) r = (qmax - qc)/dq;
        if (dq < 0.0) r = (qmin - qc)/dq;      

        phi = min( min(r*beta, 1.0), phi);
        //  }

        //
        /* if (neighbours[k3+i] < 0) { */
        /*    r = 1.0; */

        /*    if (dq > 0.0 && (sign == -1.0 || sign == 0.0 )) r = (0.0 - qc)/dq; */
        /*  if (dq < 0.0 && (sign ==  1.0 || sign == 0.0 )) r = (0.0 - qc)/dq; */

        /*    phi = min( min(r*beta, 1.0), phi); */
        /*  } */

    }

    //Update gradient, vertex and edge values using phi limiter
    x_gradient[k] = x_gradient[k]*phi;
    y_gradient[k] = y_gradient[k]*phi;

    edge_values[k3+0] = qc + phi*dqa[0];
    edge_values[k3+1] = qc + phi*dqa[1];
    edge_values[k3+2] = qc + phi*dqa[2];

    vertex_values[k3+0] = edge_values[k3+1] + edge_values[k3+2] - edge_values[k3+0];
    vertex_values[k3+1] = edge_values[k3+2] + edge_values[k3+0] - edge_values[k3+1];
    vertex_values[k3+2] = edge_values[k3+0] + edge_values[k3+1] - edge_values[k3+2];
    //}

    //return 0;
}



void extrapolate_from_gradient(
        int N,
        double* centroids,
        double* centroid_values,
        double* vertex_coordinates,
        double* vertex_values,
        double* edge_values,
        double* a,
        double* b) 
{
    const int k = 
            threadIdx.x+threadIdx.y*blockDim.x+
            (blockIdx.x+blockIdx.y*gridDim.x)*blockDim.x*blockDim.y;


    if (k >= N)
        return;

    int k2, k3, k6;
    double x, y, x0, y0, x1, y1, x2, y2;

    //for (k=0; k<N; k++){
    k6 = 6*k;
    k3 = 3*k;
    k2 = 2*k;

    // Centroid coordinates
    x = centroids[k2]; y = centroids[k2+1];

    // vertex coordinates
    // x0, y0, x1, y1, x2, y2 = X[k,:]
    x0 = vertex_coordinates[k6 + 0];
    y0 = vertex_coordinates[k6 + 1];
    x1 = vertex_coordinates[k6 + 2];
    y1 = vertex_coordinates[k6 + 3];
    x2 = vertex_coordinates[k6 + 4];
    y2 = vertex_coordinates[k6 + 5];

    // Extrapolate to Vertices
    vertex_values[k3+0] = centroid_values[k] + a[k]*(x0-x) + b[k]*(y0-y);
    vertex_values[k3+1] = centroid_values[k] + a[k]*(x1-x) + b[k]*(y1-y);
    vertex_values[k3+2] = centroid_values[k] + a[k]*(x2-x) + b[k]*(y2-y);

    // Extrapolate to Edges (midpoints)
    edge_values[k3+0] = 0.5*(vertex_values[k3 + 1]+vertex_values[k3 + 2]);
    edge_values[k3+1] = 0.5*(vertex_values[k3 + 2]+vertex_values[k3 + 0]);
    edge_values[k3+2] = 0.5*(vertex_values[k3 + 0]+vertex_values[k3 + 1]);

    //}
    return 0;
}



void compute_gradients(
        int N,
        double* centroids,
        double* centroid_values,
        long* number_of_boundaries,
        long* surrogate_neighbours,
        double* a,
        double* b)
{
    const int k = 
            threadIdx.x+threadIdx.y*blockDim.x+
            (blockIdx.x+blockIdx.y*gridDim.x)*blockDim.x*blockDim.y;

    if (k >= N)
        return;

    int i, k0, k1, k2;
    double x0, x1, x2, y0, y1, y2, q0, q1, q2;
    double det;


    //    for (k=0; k<N; k++) {

    if (number_of_boundaries[k] < 2) {
        // Two or three true neighbours

        // Get indices of neighbours (or self when used as surrogate)
        // k0, k1, k2 = surrogate_neighbours[k,:]

        k0 = surrogate_neighbours[3*k + 0];
        k1 = surrogate_neighbours[3*k + 1];
        k2 = surrogate_neighbours[3*k + 2];


        if (k0 == k1 || k1 == k2) return -1;

        // Get data
        q0 = centroid_values[k0];
        q1 = centroid_values[k1];
        q2 = centroid_values[k2];

        x0 = centroids[k0*2]; y0 = centroids[k0*2+1];
        x1 = centroids[k1*2]; y1 = centroids[k1*2+1];
        x2 = centroids[k2*2]; y2 = centroids[k2*2+1];

        // Gradient
        //_gradient(x0, y0, x1, y1, x2, y2, q0, q1, q2, &a[k], &b[k]);

        det = (y2-y0)*(x1-x0) - (y1-y0)*(x2-x0);

        a[k] = (y2-y0)*(q1-q0) - (y1-y0)*(q2-q0);
        a[k] /= det;

        b[k] = (x1-x0)*(q2-q0) - (x2-x0)*(q1-q0);
        b[k] /= det;

    } else if (number_of_boundaries[k] == 2) {
        // One true neighbour

        // Get index of the one neighbour
        i=0; k0 = k;
        while (i<3 && k0==k) {
            k0 = surrogate_neighbours[3*k + i];
            i++;
        }
        if (k0 == k) return -1;

        k1 = k; //self

        // Get data
        q0 = centroid_values[k0];
        q1 = centroid_values[k1];

        x0 = centroids[k0*2]; y0 = centroids[k0*2+1];
        x1 = centroids[k1*2]; y1 = centroids[k1*2+1];

        // Two point gradient
        //_gradient2(x0, y0, x1, y1, q0, q1, &a[k], &b[k]);

        det = (x1-x0) * (x1-x0) + (y1-y0)*(y1-y0);
        a[k] = (x1-x0)*(q1-q0) / det;
        b[k] = (y1-y0)*(q1-q0) / det;

    }
    //    else:
    //        #No true neighbours -
    //        #Fall back to first order scheme
    //    }
    return 0;
}


#else

#pragma hmpp lmtVByNeigh codelet, target=CUDA args[*].transfer=atcall
void _limit_vertices_by_all_neighbours(
        int N, 
        int N3,
        double beta,
        double centroid_values[N],
        double vertex_values[N3],
        double edge_values[N3],
        long   neighbours[N3],
        double x_gradient[N],
        double y_gradient[N]) 
{
    int i, k;
    long n;
    double qmin, qmax, qn, qc;
    double dq, phi, r;
    double dqa0, dqa1, dqa2;

    for (k=0; k<N; k++){
        qc = centroid_values[k];
        qmin = qc;
        qmax = qc;

        //for (i=0; i<3; i++) {
#ifndef REARRANGED_DOMAIN
            n = neighbours[k*3+i];
#else
            n = neighbours[k+i*N];
#endif
            if (n >= 0) {
                qn = centroid_values[n]; //Neighbour's centroid value

                qmin = fmin(qmin, qn);
                qmax = fmax(qmax, qn);
            }
#ifndef REARRANGED_DOMAIN
            n = neighbours[k*3+i];
#else
            n = neighbours[k+i*N];
#endif
            if (n >= 0) {
                qn = centroid_values[n]; //Neighbour's centroid value

                qmin = fmin(qmin, qn);
                qmax = fmax(qmax, qn);
            }
#ifndef REARRANGED_DOMAIN
            n = neighbours[k*3+i];
#else
            n = neighbours[k+i*N];
#endif
            if (n >= 0) {
                qn = centroid_values[n]; //Neighbour's centroid value

                qmin = fmin(qmin, qn);
                qmax = fmax(qmax, qn);
            }
        //} //for (i=0; i<3; i++) {

        phi = 1.0;
        //for (i=0; i<3; i++) {    
            r = 1.0;
#ifndef REARRANGED_DOMAIN
            dq = vertex_values[k*3] - qc; //Delta between vertex and centroid values
#else
            dq = vertex_values[k] -qc;
#endif
            dqa0 = dq;                   //Save dq for use in updating vertex values

            if (dq > 0.0) r = (qmax - qc)/dq;
            if (dq < 0.0) r = (qmin - qc)/dq;      


            phi = fmin( fmin(r*beta, 1.0), phi);    
            r = 1.0;
#ifndef REARRANGED_DOMAIN
            dq = vertex_values[k*3+1] - qc; //Delta between vertex and centroid values
#else
            dq = vertex_values[k+N] -qc;
#endif
            dqa1 = dq;                   //Save dq for use in updating vertex values

            if (dq > 0.0) r = (qmax - qc)/dq;
            if (dq < 0.0) r = (qmin - qc)/dq;      


            phi = fmin( fmin(r*beta, 1.0), phi);    
            r = 1.0;
#ifndef REARRANGED_DOMAIN
            dq = vertex_values[k*3+2] - qc; //Delta between vertex and centroid values
#else
            dq = vertex_values[k+2*N] -qc;
#endif
            dqa2 = dq;                   //Save dq for use in updating vertex values

            if (dq > 0.0) r = (qmax - qc)/dq;
            if (dq < 0.0) r = (qmin - qc)/dq;      


            phi = fmin( fmin(r*beta, 1.0), phi);    
        //} //for (i=0; i<3; i++) {    

        //Update gradient, vertex and edge values using phi limiter
        x_gradient[k] = x_gradient[k]*phi;
        y_gradient[k] = y_gradient[k]*phi;

#ifndef REARRANGED_DOMAIN
        vertex_values[k*3+0] = qc + phi*dqa0;
        vertex_values[k*3+1] = qc + phi*dqa1;
        vertex_values[k*3+2] = qc + phi*dqa2;

        edge_values[k*3+0] = 0.5*(vertex_values[k*3+1] + vertex_values[k*3+2]);
        edge_values[k*3+1] = 0.5*(vertex_values[k*3+2] + vertex_values[k*3+0]);
        edge_values[k*3+2] = 0.5*(vertex_values[k*3+0] + vertex_values[k*3+1]);
#else
        vertex_values[k] = qc + phi*dqa0;
        vertex_values[k+N] = qc + phi*dqa1;
        vertex_values[k+2*N] = qc + phi*dqa2;

        edge_values[k] = 0.5*(vertex_values[k+N] + vertex_values[k+2*N]);
        edge_values[k+N] = 0.5*(vertex_values[k+2*N] + vertex_values[k]);
        edge_values[k+2*N] = 0.5*(vertex_values[k] + vertex_values[k+N]);
#endif
    }
}



#pragma hmpp lmtEByNeigh codelet, target=CUDA args[*].transfer=atcall
void _limit_edges_by_all_neighbours(
        int N, 
        int N3,
        double beta,
        double centroid_values[N],
        double vertex_values[N3],
        double edge_values[N3],
        long   neighbours[N3],
        double x_gradient[N],
        double y_gradient[N]) 
{
    int i, k;
    long n;
    double qmin, qmax, qn, qc;
    double dq, phi, r;
    double dqa0, dqa1, dqa2;

    for (k=0; k<N; k++){
        qc = centroid_values[k];
        qmin = qc;
        qmax = qc;

        //for (i=0; i<3; i++) {
#ifndef REARRANGED_DOMAIN
            n = neighbours[k*3];
#else
            n = neighbours[k];
#endif
            if (n >= 0) {
                qn = centroid_values[n]; //Neighbour's centroid value

                qmin = fmin(qmin, qn);
                qmax = fmax(qmax, qn);
            }
#ifndef REARRANGED_DOMAIN
            n = neighbours[k*3+1];
#else
            n = neighbours[k+N];
#endif
            if (n >= 0) {
                qn = centroid_values[n]; //Neighbour's centroid value

                qmin = fmin(qmin, qn);
                qmax = fmax(qmax, qn);
            }
#ifndef REARRANGED_DOMAIN
            n = neighbours[k*3+2];
#else
            n = neighbours[k+2*N];
#endif
            if (n >= 0) {
                qn = centroid_values[n]; //Neighbour's centroid value

                qmin = fmin(qmin, qn);
                qmax = fmax(qmax, qn);
            }
        //} //for (i=0; i<3; i++) {
        /*
        sign = 0.0;
        if (qmin > 0.0) {
            sign = 1.0;
        } else if (qmax < 0) {
            sign = -1.0;
        }
        */
        phi = 1.0;
        //for (i=0; i<3; i++) {  
#ifndef REARRANGED_DOMAIN
            dq = edge_values[k*3] - qc; //Delta between edge and centroid values
#else
            dq = edge_values[k] -qc;
#endif
            dqa0 = dq;                 //Save dq for use in updating vertex values  

            r = 1.0;

            if (dq > 0.0) r = (qmax - qc)/dq;
            if (dq < 0.0) r = (qmin - qc)/dq;      

            phi = fmin( fmin(r*beta, 1.0), phi);
#ifndef REARRANGED_DOMAIN
            dq = edge_values[k*3+1] - qc; //Delta between edge and centroid values
#else
            dq = edge_values[k+N] -qc;
#endif
            dqa1 = dq;                 //Save dq for use in updating vertex values  

            r = 1.0;

            if (dq > 0.0) r = (qmax - qc)/dq;
            if (dq < 0.0) r = (qmin - qc)/dq;      

            phi = fmin( fmin(r*beta, 1.0), phi);
#ifndef REARRANGED_DOMAIN
            dq = edge_values[k*3+2] - qc; //Delta between edge and centroid values
#else
            dq = edge_values[k+2*N] -qc;
#endif
            dqa2 = dq;                 //Save dq for use in updating vertex values  

            r = 1.0;

            if (dq > 0.0) r = (qmax - qc)/dq;
            if (dq < 0.0) r = (qmin - qc)/dq;      

            phi = fmin( fmin(r*beta, 1.0), phi);
        //} //for (i=0; i<3; i++) {  



        //Update gradient, vertex and edge values using phi limiter
        x_gradient[k] = x_gradient[k]*phi;
        y_gradient[k] = y_gradient[k]*phi;

#ifndef REARRANGED_DOMAIN
        edge_values[k*3+0] = qc + phi*dqa0;
        edge_values[k*3+1] = qc + phi*dqa1;
        edge_values[k*3+2] = qc + phi*dqa2;

        vertex_values[k*3+0] =edge_values[k*3+1] +edge_values[k*3+2] - edge_values[k*3+0];
        vertex_values[k*3+1] =edge_values[k*3+2] +edge_values[k*3+0] - edge_values[k*3+1];
        vertex_values[k*3+2] =edge_values[k*3+0] +edge_values[k*3+1] - edge_values[k*3+2];
#else
        edge_values[k] = qc + phi*dqa0;
        edge_values[k+N] = qc + phi*dqa1;
        edge_values[k+2*N] = qc + phi*dqa2;

        vertex_values[k] =edge_values[k+N] +edge_values[k+2*N] - edge_values[k];
        vertex_values[k+N] =edge_values[k+2*N] +edge_values[k] - edge_values[k+N];
        vertex_values[k+2*N] =edge_values[k] +edge_values[k+N] - edge_values[k+2*N];
#endif
    }
}



#pragma hmpp extraFromGradient codelet, target=CUDA args[*].transfer=atcall
void _extrapolate_from_gradient(
        int N,
        int N2,
        int N3,
        int N6,
        double centroids[N2],
        double centroid_values[N],
        double vertex_coordinates[N6],
        double vertex_values[N3],
        double edge_values[N3],
        double a[N],
        double b[N]) 
{
    int k;
    double x, y, x0, y0, x1, y1, x2, y2;

    for (k=0; k<N; k++){

        // Centroid coordinates
#ifndef REARRANGED_DOMAIN
        x = centroids[k*2]; y = centroids[k*2+1];
#else
        x = centroids[k]; y = centroids[k + N];
#endif

        // vertex coordinates
        // x0, y0, x1, y1, x2, y2 = X[k,:]
#ifndef REARRANGED_DOMAIN
        x0 = vertex_coordinates[k*6 + 0];
        y0 = vertex_coordinates[k*6 + 1];
        x1 = vertex_coordinates[k*6 + 2];
        y1 = vertex_coordinates[k*6 + 3];
        x2 = vertex_coordinates[k*6 + 4];
        y2 = vertex_coordinates[k*6 + 5];
#else
        x0 = vertex_coordinates[k];
        y0 = vertex_coordinates[k + N];
        x1 = vertex_coordinates[k + 2*N];
        y1 = vertex_coordinates[k + 3*N];
        x2 = vertex_coordinates[k + 4*N];
        y2 = vertex_coordinates[k + 5*N];
#endif

#ifndef REARRANGED_DOMAIN
        // Extrapolate to Vertices
        vertex_values[k*3+0] = centroid_values[k] + a[k]*(x0-x) + b[k]*(y0-y);
        vertex_values[k*3+1] = centroid_values[k] + a[k]*(x1-x) + b[k]*(y1-y);
        vertex_values[k*3+2] = centroid_values[k] + a[k]*(x2-x) + b[k]*(y2-y);

        // Extrapolate to Edges (midpoints)
        edge_values[k*3+0] = 0.5*(vertex_values[k*3 + 1]+vertex_values[k*3 + 2]);
        edge_values[k*3+1] = 0.5*(vertex_values[k*3 + 2]+vertex_values[k*3 + 0]);
        edge_values[k*3+2] = 0.5*(vertex_values[k*3 + 0]+vertex_values[k*3 + 1]);
#else
        // Extrapolate to Vertices
        vertex_values[k] = centroid_values[k] + a[k]*(x0-x) + b[k]*(y0-y);
        vertex_values[k+N] = centroid_values[k] + a[k]*(x1-x) + b[k]*(y1-y);
        vertex_values[k+2*N] = centroid_values[k] + a[k]*(x2-x) + b[k]*(y2-y);

        // Extrapolate to Edges (midpoints)
        edge_values[k] = 0.5*(vertex_values[k + N]+vertex_values[k + 2*N]);
        edge_values[k+N] = 0.5*(vertex_values[k + 2*N]+vertex_values[k]);
        edge_values[k+2*N] = 0.5*(vertex_values[k]+vertex_values[k + N]);
#endif

    }
}


#pragma hmpp cptGradients codelet, target=CUDA args[*].transfer=atcall
void _compute_gradients(
        int N,
        int N2,
        int N3,
        double centroids[N2],
        double centroid_values[N],
        long number_of_boundaries[N],
        long surrogate_neighbours[N3],
        double a[N],
        double b[N])
{
    int k, k0, k1, k2;
    double x0, x1, x2, y0, y1, y2, q0, q1, q2;
    double det;


    for (k=0; k<N; k++) {

        if (number_of_boundaries[k] < 2) {
#ifndef REARRANGED_DOMAIN
            k0 = surrogate_neighbours[k*3 + 0];
            k1 = surrogate_neighbours[k*3 + 1];
            k2 = surrogate_neighbours[k*3 + 2];
#else
            k0 = surrogate_neighbours[k];
            k1 = surrogate_neighbours[k + N];
            k2 = surrogate_neighbours[k + 2*N];
#endif


            // Get data
            q0 = centroid_values[k0];
            q1 = centroid_values[k1];
            q2 = centroid_values[k2];
#ifndef REARRANGED_DOMAIN
            x0 = centroids[k0*2]; y0 = centroids[k0*2+1];
            x1 = centroids[k1*2]; y1 = centroids[k1*2+1];
            x2 = centroids[k2*2]; y2 = centroids[k2*2+1];
#else
            x0 = centroids[k0]; y0 = centroids[k0+N];
            x1 = centroids[k1]; y1 = centroids[k1+N];
            x2 = centroids[k2]; y2 = centroids[k2+N];
#endif
            // Gradient
            //_gradient(x0, y0, x1, y1, x2, y2, q0, q1, q2, &a[k], &b[k]);

            det = (y2-y0)*(x1-x0) - (y1-y0)*(x2-x0);

            a[k] = (y2-y0)*(q1-q0) - (y1-y0)*(q2-q0);
            a[k] /= det;

            b[k] = (x1-x0)*(q2-q0) - (x2-x0)*(q1-q0);
            b[k] /= det;
        } else if (number_of_boundaries[k] == 2) {
            // One true neighbour
            // Get index of the one neighbour
            //i=0; k0 = k;
            //while (i<3 && k0==k) {
#ifndef REARRANGED_DOMAIN
            k0 = surrogate_neighbours[3*k + 2];
            if ( k0 == k )
                k0 = surrogate_neighbours[3*k + 1];
            if ( k0 == k )
                k0 = surrogate_neighbours[3*k];
#else
            k0 = surrogate_neighbours[k + 2*N];
            if ( k0 == k )
                k0 = surrogate_neighbours[k + N];
            if ( k0 == k )
                k0 = surrogate_neighbours[k];
#endif

            // This is the code where compiler report : Loop not normalized: Cannot find
            //    induction variable
            //if (k0 == k) 
            ////    return -1;
            //    break;

            k1 = k; //self

            // Get data
            q0 = centroid_values[k0];
            q1 = centroid_values[k1];

#ifndef REARRANGED_DOMAIN
            x0 = centroids[k0*2]; y0 = centroids[k0*2+1];
            x1 = centroids[k1*2]; y1 = centroids[k1*2+1];
#else
            x0 = centroids[k0]; y0 = centroids[k0+N];
            x1 = centroids[k1]; y1 = centroids[k1+N];
#endif

            // Two point gradient
            //_gradient2(x0, y0, x1, y1, q0, q1, &a[k], &b[k]);

            det = (x1-x0) * (x1-x0) + (y1-y0)*(y1-y0);
            a[k] = (x1-x0)*(q1-q0) / det;
            b[k] = (y1-y0)*(q1-q0) / det;
        }

    }
}



void extrapolate_second_order_and_limit_by_vertex(
        int N,
        int N2,
        int N3,
        int N6,
        double beta,

        double domain_centroid_coordinates[N2],
        double domain_vertex_coordinates[N6],
        long domain_number_of_boundaries[N],
        long domain_surrogate_neighbours[N3],
        long domain_neighbours[N3],

        double quantity_centroid_values[N],
        double quantity_vertex_values[N3],
        double quantity_edge_values[N3],
        double quantity_x_gradient[N],
        double quantity_y_gradient[N]
        )
{
    _compute_gradients(
            N,
            N2,
            N3,
            domain_centroid_coordinates,
            quantity_centroid_values,
            domain_number_of_boundaries,
            domain_surrogate_neighbours,
            quantity_x_gradient,
            quantity_y_gradient);

    _extrapolate_from_gradient(
            N,
            N2,
            N3,
            N6,
            domain_centroid_coordinates,
            quantity_centroid_values,
            domain_vertex_coordinates,
            quantity_vertex_values,
            quantity_edge_values,
            quantity_x_gradient,
            quantity_y_gradient);

    _limit_vertices_by_all_neighbours(
            N,
            N3,
            beta,
            quantity_centroid_values,
            quantity_vertex_values,
            quantity_edge_values,
            domain_neighbours,
            quantity_x_gradient,
            quantity_y_gradient);
}



void extrapolate_second_order_and_limit_by_edge(
        int N,
        int N2,
        int N3,
        int N6,
        double beta,

        double domain_centroid_coordinates[N2],
        double domain_vertex_coordinates[N6],
        long domain_number_of_boundaries[N],
        long domain_surrogate_neighbours[N3],
        long domain_neighbours[N3],

        double quantity_centroid_values[N],
        double quantity_vertex_values[N3],
        double quantity_edge_values[N3],
        double quantity_x_gradient[N],
        double quantity_y_gradient[N]
        )
{
    _compute_gradients(
            N,
            N2,
            N3,
            domain_centroid_coordinates,
            quantity_centroid_values,
            domain_number_of_boundaries,
            domain_surrogate_neighbours,
            quantity_x_gradient,
            quantity_y_gradient);

    _extrapolate_from_gradient(
            N,
            N2,
            N3,
            N6,
            domain_centroid_coordinates,
            quantity_centroid_values,
            domain_vertex_coordinates,
            quantity_vertex_values,
            quantity_edge_values,
            quantity_x_gradient,
            quantity_y_gradient);

    _limit_edges_by_all_neighbours(
            N,
            N3,
            beta,
            quantity_centroid_values,
            quantity_vertex_values,
            quantity_edge_values,
            domain_neighbours,
            quantity_x_gradient,
            quantity_y_gradient);
}

#endif

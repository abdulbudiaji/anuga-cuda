#ifndef USING_HOST_MACRO
#define BLOCK_SIZE  8
#endif


__device__ int limit_gradient(double *dqv, double qmin, double qmax, double beta_w) 
{
    int i;
    double r = 1000.0, r0 = 1.0, phi = 1.0;
    double TINY = 1.0e-100; // to avoid machine accuracy problems.

    for (i = 0; i < 3; i++) {
        if (dqv[i]<-TINY)
            r0 = qmin / dqv[i];

        if (dqv[i] > TINY)
            r0 = qmax / dqv[i];

        r = min(r0, r);
    }

    phi = min(r*beta_w, 1.0);
    //for (i=0;i<3;i++)
    dqv[0] = dqv[0] * phi;
    dqv[1] = dqv[1] * phi;
    dqv[2] = dqv[2] * phi;

    return 0;
}   

__device__ int find_qmin_and_qmax(double dq0, double dq1, double dq2,
        double *qmin, double *qmax) 
{
    if (dq0 >= 0.0) {
        if (dq1 >= dq2) {
            if (dq1 >= 0.0)
                *qmax = dq0 + dq1;
            else
                *qmax = dq0;

            *qmin = dq0 + dq2;
            if (*qmin >= 0.0) *qmin = 0.0;
        } else {// dq1<dq2
            if (dq2 > 0)
                *qmax = dq0 + dq2;
            else
                *qmax = dq0;

            *qmin = dq0 + dq1;
            if (*qmin >= 0.0) *qmin = 0.0;
        }
    } else {//dq0<0
        if (dq1 <= dq2) {
            if (dq1 < 0.0)
                *qmin = dq0 + dq1;
            else
                *qmin = dq0;

            *qmax = dq0 + dq2;
            if (*qmax <= 0.0) *qmax = 0.0;
        } else {// dq1>dq2
            if (dq2 < 0.0)
                *qmin = dq0 + dq2;
            else
                *qmin = dq0;

            *qmax = dq0 + dq1;
            if (*qmax <= 0.0) *qmax = 0.0;
        }
    }
    return 0;
}


__global__ void extrapolate_velocity_second_order_true(
            int N,
            double minimum_allowed_height,
            double * stage_centroid_values,
            double * bed_centroid_values,
            double * xmom_centroid_values,
            double * xmom_centroid_store,
            double * ymom_centroid_values,
            double * ymom_centroid_store
            )
{
    const int k = threadIdx.x + threadIdx.y*blockDim.x +
        (blockIdx.x+blockIdx.y*gridDim.x)*blockDim.x*blockDim.y;
    
    double dk;

    if (k >= N )
        return;

    dk = max(stage_centroid_values[k] -bed_centroid_values[k], minimum_allowed_height);
    xmom_centroid_store[k] = xmom_centroid_values[k];
    xmom_centroid_values[k] = xmom_centroid_values[k] / dk;

    ymom_centroid_store[k] = ymom_centroid_values[k];
    ymom_centroid_values[k] = ymom_centroid_values[k] / dk;
}


__global__ void extrapolate_velocity_second_order_true_after(
            int N,
            double * xmom_centroid_values,
            double * xmom_centroid_store,
            double * ymom_centroid_values,
            double * ymom_centroid_store
            )
{
    const int k = threadIdx.x + threadIdx.y*blockDim.x +
        (blockIdx.x+blockIdx.y*gridDim.x)*blockDim.x*blockDim.y;
    

    if (k >= N )
        return;

    xmom_centroid_values[k] = xmom_centroid_store[k];

    ymom_centroid_values[k] = ymom_centroid_store[k];
}

    

__global__ void extrapolate_second_order_sw_true (
        int N,
        double epsilon,
        double minimum_allowed_height,
        double beta_w,
        double beta_w_dry,
        double beta_uh,
        double beta_uh_dry,
        double beta_vh,
        double beta_vh_dry,
        int optimise_dry_cells,

        long* surrogate_neighbours,
        long* number_of_boundaries,
        double* centroid_coordinates,

        double* stage_centroid_values,
        double* bed_centroid_values,
        double* xmom_centroid_values,
        double* ymom_centroid_values,

        double* vertex_coordinates,
        
        double* stage_vertex_values,
        double* bed_vertex_values,
        double* xmom_vertex_values,
        double* ymom_vertex_values,
        
        double* stage_centroid_store,
        double* xmom_centroid_store,
        double* ymom_centroid_store)
{
    const int k = threadIdx.x + threadIdx.y*blockDim.x +
        (blockIdx.x+blockIdx.y*gridDim.x)*blockDim.x*blockDim.y;

    if (k >= N )
        return;

    int k3 = 3*k,
        k6 = 6*k;

    double a, b;
    int k0, k1, k2, coord_index, i;
    double x, y, x0, y0, x1, y1, x2, y2, xv0, yv0, xv1, yv1, xv2, yv2;
    double dx1, dx2, dy1, dy2, dxv0, dxv1, dxv2, dyv0, dyv1, dyv2, dq0, dq1, dq2, area2, inv_area2;
    double dqv[3], qmin, qmax, hmin, hmax;
    double hc, h0, h1, h2, beta_tmp, hfactor;

    double dk, dv0, dv1, dv2;

    __shared__ double sh_xmom_vertex_values[ BLOCK_SIZE * 3];

    // extrapolate_velocity_second_order == 1 

    //dk = max(stage_centroid_values[k]-bed_centroid_values[k],minimum_allowed_height);
    //xmom_centroid_store[k] = xmom_centroid_values[k];
    //xmom_centroid_values[k] = xmom_centroid_values[k] / dk;

    //ymom_centroid_store[k] = ymom_centroid_values[k];
    //ymom_centroid_values[k] = ymom_centroid_values[k] / dk;



    // Begin extrapolation routine
    do {
        if ( number_of_boundaries[k] == 3 )
        {
            stage_vertex_values[k3] = stage_centroid_values[k];
            stage_vertex_values[k3 +1] = stage_centroid_values[k];
            stage_vertex_values[k3 +2] = stage_centroid_values[k];

            //xmom_vertex_values[k3] = xmom_centroid_values[k];
            //xmom_vertex_values[k3 + 1] = xmom_centroid_values[k];
            //xmom_vertex_values[k3 + 2] = xmom_centroid_values[k];
            sh_xmom_vertex_values[threadIdx.x] = xmom_centroid_values[k];
            sh_xmom_vertex_values[threadIdx.x + blockDim.x] = sh_xmom_vertex_values[threadIdx.x];
            sh_xmom_vertex_values[threadIdx.x + blockDim.x*2]=sh_xmom_vertex_values[threadIdx.x];

            ymom_vertex_values[k3] = ymom_centroid_values[k];
            ymom_vertex_values[k3 + 1] = ymom_centroid_values[k];
            ymom_vertex_values[k3 + 2] = ymom_centroid_values[k];

            // continue;
            break;
        }
        else
        {
            xv0 = vertex_coordinates[k6];
            yv0 = vertex_coordinates[k6 + 1];
            xv1 = vertex_coordinates[k6 + 2];
            yv1 = vertex_coordinates[k6 + 3];
            xv2 = vertex_coordinates[k6 + 4];
            yv2 = vertex_coordinates[k6 + 5];

            coord_index = 2 * k;
            x = centroid_coordinates[coord_index];
            y = centroid_coordinates[coord_index + 1];

            dxv0 = xv0 - x;
            dxv1 = xv1 - x;
            dxv2 = xv2 - x;
            dyv0 = yv0 - y;
            dyv1 = yv1 - y;
            dyv2 = yv2 - y;
        }

        if ( number_of_boundaries[k] <= 1 )
        {
            k0 = surrogate_neighbours[k3];
            k1 = surrogate_neighbours[k3 + 1];
            k2 = surrogate_neighbours[k3 + 2];

            // Get the auxiliary triangle's vertex coordinates
            // (really the centroids of neighbouring triangles)
            coord_index = 2 * k0;
            x0 = centroid_coordinates[coord_index];
            y0 = centroid_coordinates[coord_index + 1];

            coord_index = 2 * k1;
            x1 = centroid_coordinates[coord_index];
            y1 = centroid_coordinates[coord_index + 1];

            coord_index = 2 * k2;
            x2 = centroid_coordinates[coord_index];
            y2 = centroid_coordinates[coord_index + 1];

            // Store x- and y- differentials for the vertices
            // of the auxiliary triangle
            dx1 = x1 - x0;
            dx2 = x2 - x0;
            dy1 = y1 - y0;
            dy2 = y2 - y0;

            // Calculate 2*area of the auxiliary triangle
            // The triangle is guaranteed to be counter-clockwise
            area2 = dy2 * dx1 - dy1*dx2;

            if (area2 <= 0) {
                stage_vertex_values[k3] = stage_centroid_values[k];
                stage_vertex_values[k3 + 1] = stage_centroid_values[k];
                stage_vertex_values[k3 + 2] = stage_centroid_values[k];

                //xmom_vertex_values[k3] = xmom_centroid_values[k];
                //xmom_vertex_values[k3 + 1] = xmom_centroid_values[k];
                //xmom_vertex_values[k3 + 2] = xmom_centroid_values[k];
                sh_xmom_vertex_values[threadIdx.x] = xmom_centroid_values[k];
                sh_xmom_vertex_values[threadIdx.x + blockDim.x] = sh_xmom_vertex_values[threadIdx.x];
                sh_xmom_vertex_values[threadIdx.x + blockDim.x*2]=sh_xmom_vertex_values[threadIdx.x];
                
                ymom_vertex_values[k3] = ymom_centroid_values[k];
                ymom_vertex_values[k3 + 1] = ymom_centroid_values[k];
                ymom_vertex_values[k3 + 2] = ymom_centroid_values[k];

                //continue;
                break;
            }

            // Calculate heights of neighbouring cells
            hc = stage_centroid_values[k] - bed_centroid_values[k];
            h0 = stage_centroid_values[k0] - bed_centroid_values[k0];
            h1 = stage_centroid_values[k1] - bed_centroid_values[k1];
            h2 = stage_centroid_values[k2] - bed_centroid_values[k2];
            hmin = min(min(h0, min(h1, h2)), hc);
            //hfactor = hc/(hc + 1.0);

            hfactor = 0.0;
            if (hmin > 0.001) {
                hfactor = (hmin - 0.001) / (hmin + 0.004);
            }

            if (optimise_dry_cells) {
                // Check if linear reconstruction is necessary for triangle k
                // This check will exclude dry cells.

                hmax = max(h0, max(h1, h2));
                if (hmax < epsilon) {
                    // continue;
                    break;
                }
            }

            //-----------------------------------
            // stage
            //-----------------------------------

            // Calculate the difference between vertex 0 of the auxiliary
            // triangle and the centroid of triangle k
            dq0 = stage_centroid_values[k0] - stage_centroid_values[k];

            // Calculate differentials between the vertices
            // of the auxiliary triangle (centroids of neighbouring triangles)
            dq1 = stage_centroid_values[k1] - stage_centroid_values[k0];
            dq2 = stage_centroid_values[k2] - stage_centroid_values[k0];

            inv_area2 = 1.0 / area2;
            // Calculate the gradient of stage on the auxiliary triangle
            a = dy2 * dq1 - dy1*dq2;
            a *= inv_area2;
            b = dx1 * dq2 - dx2*dq1;
            b *= inv_area2;

            // Calculate provisional jumps in stage from the centroid
            // of triangle k to its vertices, to be limited
            dqv[0] = a * dxv0 + b*dyv0;
            dqv[1] = a * dxv1 + b*dyv1;
            dqv[2] = a * dxv2 + b*dyv2;

            // Now we want to find min and max of the centroid and the
            // vertices of the auxiliary triangle and compute jumps
            // from the centroid to the min and max
            find_qmin_and_qmax(dq0, dq1, dq2, &qmin, &qmax);

            // Playing with dry wet interface
            //hmin = qmin;
            //beta_tmp = beta_tmp_dry;
            //if (hmin>minimum_allowed_height)
            beta_tmp = beta_w_dry + (beta_w - beta_w_dry) * hfactor;


            limit_gradient(dqv, qmin, qmax, beta_tmp);

            //for (i=0;i<3;i++)
            stage_vertex_values[k3 + 0] = stage_centroid_values[k] + dqv[0];
            stage_vertex_values[k3 + 1] = stage_centroid_values[k] + dqv[1];
            stage_vertex_values[k3 + 2] = stage_centroid_values[k] + dqv[2];


            //-----------------------------------
            // xmomentum
            //-----------------------------------

            // Calculate the difference between vertex 0 of the auxiliary
            // triangle and the centroid of triangle k
            dq0 = xmom_centroid_values[k0] - xmom_centroid_values[k];
            stage_centroid_store[k] = xmom_centroid_values[k0];

            // Calculate differentials between the vertices
            // of the auxiliary triangle
            dq1 = xmom_centroid_values[k1] - xmom_centroid_values[k0];
            dq2 = xmom_centroid_values[k2] - xmom_centroid_values[k0];

            // Calculate the gradient of xmom on the auxiliary triangle
            a = dy2 * dq1 - dy1*dq2;
            a *= inv_area2;
            b = dx1 * dq2 - dx2*dq1;
            b *= inv_area2;

            // Calculate provisional jumps in stage from the centroid
            // of triangle k to its vertices, to be limited
            dqv[0] = a * dxv0 + b*dyv0;
            dqv[1] = a * dxv1 + b*dyv1;
            dqv[2] = a * dxv2 + b*dyv2;

            // Now we want to find min and max of the centroid and the
            // vertices of the auxiliary triangle and compute jumps
            // from the centroid to the min and max
            find_qmin_and_qmax(dq0, dq1, dq2, &qmin, &qmax);
            //beta_tmp = beta_uh;
            //if (hmin<minimum_allowed_height)
            //beta_tmp = beta_uh_dry;
            beta_tmp = beta_uh_dry + (beta_uh - beta_uh_dry) * hfactor;

            // Limit the gradient
            limit_gradient(dqv, qmin, qmax, beta_tmp);

            //for (i = 0; i < 3; i++) {
            //    xmom_vertex_values[k3 + i] = xmom_centroid_values[k] + dqv[i];
            //}

            sh_xmom_vertex_values[threadIdx.x] = xmom_centroid_values[k];
            sh_xmom_vertex_values[threadIdx.x + blockDim.x] = sh_xmom_vertex_values[threadIdx.x];
            sh_xmom_vertex_values[threadIdx.x + blockDim.x*2]=sh_xmom_vertex_values[threadIdx.x];
            
            sh_xmom_vertex_values[threadIdx.x]  += dqv[0];
            sh_xmom_vertex_values[threadIdx.x + blockDim.x]  += dqv[1];
            sh_xmom_vertex_values[threadIdx.x + blockDim.x*2] += dqv[2];


            //-----------------------------------
            // ymomentum
            //-----------------------------------

            // Calculate the difference between vertex 0 of the auxiliary
            // triangle and the centroid of triangle k
            dq0 = ymom_centroid_values[k0] - ymom_centroid_values[k];

            // Calculate differentials between the vertices
            // of the auxiliary triangle
            dq1 = ymom_centroid_values[k1] - ymom_centroid_values[k0];
            dq2 = ymom_centroid_values[k2] - ymom_centroid_values[k0];

            // Calculate the gradient of xmom on the auxiliary triangle
            a = dy2 * dq1 - dy1*dq2;
            a *= inv_area2;
            b = dx1 * dq2 - dx2*dq1;
            b *= inv_area2;

            // Calculate provisional jumps in stage from the centroid
            // of triangle k to its vertices, to be limited
            dqv[0] = a * dxv0 + b*dyv0;
            dqv[1] = a * dxv1 + b*dyv1;
            dqv[2] = a * dxv2 + b*dyv2;

            // Now we want to find min and max of the centroid and the
            // vertices of the auxiliary triangle and compute jumps
            // from the centroid to the min and max
            find_qmin_and_qmax(dq0, dq1, dq2, &qmin, &qmax);

            //beta_tmp = beta_vh;
            //
            //if (hmin<minimum_allowed_height)
            //beta_tmp = beta_vh_dry;
            beta_tmp = beta_vh_dry + (beta_vh - beta_vh_dry) * hfactor;

            // Limit the gradient
            limit_gradient(dqv, qmin, qmax, beta_tmp);

            for (i = 0; i < 3; i++) {
                ymom_vertex_values[k3 + i] = ymom_centroid_values[k] + dqv[i];
            }
        }// End number_of_boundaries <=1
        else
        {
            //==============================================
            // Number of boundaries == 2
            //==============================================

            // One internal neighbour and gradient is in direction of the neighbour's centroid

            // Find the only internal neighbour (k1?)
            for (k2 = k3; k2 < k3 + 3; k2++) {
                // Find internal neighbour of triangle k
                // k2 indexes the edges of triangle k

                if (surrogate_neighbours[k2] != k) {
                    break;
                }
            }



            k1 = surrogate_neighbours[k2];

            // The coordinates of the triangle are already (x,y).
            // Get centroid of the neighbour (x1,y1)
            coord_index = 2 * k1;
            x1 = centroid_coordinates[coord_index];
            y1 = centroid_coordinates[coord_index + 1];

            // Compute x- and y- distances between the centroid of
            // triangle k and that of its neighbour
            dx1 = x1 - x;
            dy1 = y1 - y;

            // Set area2 as the square of the distance
            area2 = dx1 * dx1 + dy1*dy1;

            // Set dx2=(x1-x0)/((x1-x0)^2+(y1-y0)^2)
            // and dy2=(y1-y0)/((x1-x0)^2+(y1-y0)^2) which
            // respectively correspond to the x- and y- gradients
            // of the conserved quantities
            dx2 = 1.0 / area2;
            dy2 = dx2*dy1;
            dx2 *= dx1;


            //-----------------------------------
            // stage
            //-----------------------------------

            // Compute differentials
            dq1 = stage_centroid_values[k1] - stage_centroid_values[k];

            // Calculate the gradient between the centroid of triangle k
            // and that of its neighbour
            a = dq1*dx2;
            b = dq1*dy2;

            // Calculate provisional vertex jumps, to be limited
            dqv[0] = a * dxv0 + b*dyv0;
            dqv[1] = a * dxv1 + b*dyv1;
            dqv[2] = a * dxv2 + b*dyv2;

            // Now limit the jumps
            if (dq1 >= 0.0) {
                qmin = 0.0;
                qmax = dq1;
            } else {
                qmin = dq1;
                qmax = 0.0;
            }

            // Limit the gradient
            limit_gradient(dqv, qmin, qmax, beta_w);

            //for (i=0; i < 3; i++)
            //{
            stage_vertex_values[k3] = stage_centroid_values[k] + dqv[0];
            stage_vertex_values[k3 + 1] = stage_centroid_values[k] + dqv[1];
            stage_vertex_values[k3 + 2] = stage_centroid_values[k] + dqv[2];
            //}

            //-----------------------------------
            // xmomentum
            //-----------------------------------

            // Compute differentials
            dq1 = xmom_centroid_values[k1] - xmom_centroid_values[k];

            // Calculate the gradient between the centroid of triangle k
            // and that of its neighbour
            a = dq1*dx2;
            b = dq1*dy2;

            // Calculate provisional vertex jumps, to be limited
            dqv[0] = a * dxv0 + b*dyv0;
            dqv[1] = a * dxv1 + b*dyv1;
            dqv[2] = a * dxv2 + b*dyv2;

            // Now limit the jumps
            if (dq1 >= 0.0) {
                qmin = 0.0;
                qmax = dq1;
            } else {
                qmin = dq1;
                qmax = 0.0;
            }

            // Limit the gradient
            limit_gradient(dqv, qmin, qmax, beta_w);

            //for (i=0;i<3;i++)
            //xmom_vertex_values[k3] = xmom_centroid_values[k] + dqv[0];
            //xmom_vertex_values[k3 + 1] = xmom_centroid_values[k] + dqv[1];
            //xmom_vertex_values[k3 + 2] = xmom_centroid_values[k] + dqv[2];


            //for (i = 0; i < 3; i++) {
            //    xmom_vertex_values[k3 + i] = xmom_centroid_values[k] + dqv[i];
            //}


            sh_xmom_vertex_values[threadIdx.x] = xmom_centroid_values[k];
            sh_xmom_vertex_values[threadIdx.x + blockDim.x] = sh_xmom_vertex_values[threadIdx.x];
            sh_xmom_vertex_values[threadIdx.x + blockDim.x*2]=sh_xmom_vertex_values[threadIdx.x];
            
            sh_xmom_vertex_values[threadIdx.x]  += dqv[0];
            sh_xmom_vertex_values[threadIdx.x + blockDim.x]  += dqv[1];
            sh_xmom_vertex_values[threadIdx.x + blockDim.x*2] += dqv[2];



            //-----------------------------------
            // ymomentum
            //-----------------------------------

            // Compute differentials
            dq1 = ymom_centroid_values[k1] - ymom_centroid_values[k];

            // Calculate the gradient between the centroid of triangle k
            // and that of its neighbour
            a = dq1*dx2;
            b = dq1*dy2;

            // Calculate provisional vertex jumps, to be limited
            dqv[0] = a * dxv0 + b*dyv0;
            dqv[1] = a * dxv1 + b*dyv1;
            dqv[2] = a * dxv2 + b*dyv2;

            // Now limit the jumps
            if (dq1 >= 0.0) {
                qmin = 0.0;
                qmax = dq1;
            }
            else {
                qmin = dq1;
                qmax = 0.0;
            }

            // Limit the gradient
            limit_gradient(dqv, qmin, qmax, beta_w);

            //for (i=0;i<3;i++)
            //ymom_vertex_values[k3] = ymom_centroid_values[k] + dqv[0];
            //ymom_vertex_values[k3 + 1] = ymom_centroid_values[k] + dqv[1];
            //ymom_vertex_values[k3 + 2] = ymom_centroid_values[k] + dqv[2];

            for (i = 0; i < 3; i++) {
                ymom_vertex_values[k3 + i] = ymom_centroid_values[k] + dqv[i];
            }
            //ymom_vertex_values[k3] = ymom_centroid_values[k] + dqv[0];
            //ymom_vertex_values[k3 + 1] = ymom_centroid_values[k] + dqv[1];
            //ymom_vertex_values[k3 + 2] = ymom_centroid_values[k] + dqv[2];
        } // else [number_of_boundaries==2]
    } while(0);




    //if (extrapolate_velocity_second_order == 1) 
    // Convert back from velocity to momentum
    //for (k = 0; k < number_of_elements; k++) 
    k3 = 3 * k;
    //dv0 = max(stage_vertex_values[k3]-bed_vertex_values[k3],minimum_allowed_height);
    //dv1 = max(stage_vertex_values[k3+1]-bed_vertex_values[k3+1],minimum_allowed_height);
    //dv2 = max(stage_vertex_values[k3+2]-bed_vertex_values[k3+2],minimum_allowed_height);
    dv0 = max(stage_vertex_values[k3] - bed_vertex_values[k3], 0.);
    dv1 = max(stage_vertex_values[k3 + 1] - bed_vertex_values[k3 + 1], 0.);
    dv2 = max(stage_vertex_values[k3 + 2] - bed_vertex_values[k3 + 2], 0.);

    //Correct centroid and vertex values
    //xmom_centroid_values[k] = xmom_centroid_store[k];
    //xmom_vertex_values[k3] = xmom_vertex_values[k3] * dv0;
    //xmom_vertex_values[k3 + 1] = xmom_vertex_values[k3 + 1] * dv1;
    //xmom_vertex_values[k3 + 2] = xmom_vertex_values[k3 + 2] * dv2;
    xmom_vertex_values[k3] = sh_xmom_vertex_values[threadIdx.x] * dv0;
    xmom_vertex_values[k3+1]=sh_xmom_vertex_values[threadIdx.x + BLOCK_SIZE]*dv1;
    xmom_vertex_values[k3+2]=sh_xmom_vertex_values[threadIdx.x + BLOCK_SIZE*2]*dv1;

    //ymom_centroid_values[k] = ymom_centroid_store[k];
    ymom_vertex_values[k3] = ymom_vertex_values[k3] * dv0;
    ymom_vertex_values[k3 + 1] = ymom_vertex_values[k3 + 1] * dv1;
    ymom_vertex_values[k3 + 2] = ymom_vertex_values[k3 + 2] * dv2;

}



__global__ void extrapolate_second_order_sw_false (
        int N,
        double epsilon,
        double minimum_allowed_height,
        double beta_w,
        double beta_w_dry,
        double beta_uh,
        double beta_uh_dry,
        double beta_vh,
        double beta_vh_dry,
        double optimise_dry_cells,
        long* surrogate_neighbours,
        long* number_of_boundaries,
        double* centroid_coordinates,
        double* stage_centroid_values,
        double* bed_centroid_values,
        double* xmom_centroid_values,
        double* ymom_centroid_values,
        double* vertex_coordinates,
        double* stage_vertex_values,
        double* bed_vertex_values,
        double* xmom_vertex_values,
        double* ymom_vertex_values)
{
    const long k = threadIdx.x + threadIdx.y*blockDim.x +
                (blockIdx.x+blockIdx.y*gridDim.x)*blockDim.x*blockDim.y;

    if ( k >= N )
        return;

    int k3 = 3*k,
        k6 = 6*k;

    double a, b;
    int k0, k1, k2, coord_index, i;
    double x, y, x0, y0, x1, y1, x2, y2, xv0, yv0, xv1, yv1, xv2, yv2;
    double dx1, dx2, dy1, dy2, dxv0, dxv1, dxv2, dyv0, dyv1, dyv2, dq0, dq1, dq2, area2, inv_area2;
    double dqv[3], qmin, qmax, hmin, hmax;
    double hc, h0, h1, h2, beta_tmp, hfactor;

    // Begin extrapolation routine

    if ( number_of_boundaries[k] == 3 )
    {
        stage_vertex_values[k3] = stage_centroid_values[k];
        stage_vertex_values[k3 +1] = stage_centroid_values[k];
        stage_vertex_values[k3 +2] = stage_centroid_values[k];

        xmom_vertex_values[k3] = xmom_centroid_values[k];
        xmom_vertex_values[k3 + 1] = xmom_centroid_values[k];
        xmom_vertex_values[k3 + 2] = xmom_centroid_values[k];

        ymom_vertex_values[k3] = ymom_centroid_values[k];
        ymom_vertex_values[k3 + 1] = ymom_centroid_values[k];
        ymom_vertex_values[k3 + 2] = ymom_centroid_values[k];

        // continue;
        return;
    }
    else
    {
        xv0 = vertex_coordinates[k6];
        yv0 = vertex_coordinates[k6 + 1];
        xv1 = vertex_coordinates[k6 + 2];
        yv1 = vertex_coordinates[k6 + 3];
        xv2 = vertex_coordinates[k6 + 4];
        yv2 = vertex_coordinates[k6 + 5];

        coord_index = 2 * k;
        x = centroid_coordinates[coord_index];
        y = centroid_coordinates[coord_index + 1];

        dxv0 = xv0 - x;
        dxv1 = xv1 - x;
        dxv2 = xv2 - x;
        dyv0 = yv0 - y;
        dyv1 = yv1 - y;
        dyv2 = yv2 - y;
    }

    if ( number_of_boundaries[k] <= 1 )
    {
        k0 = surrogate_neighbours[k3];
        k1 = surrogate_neighbours[k3 + 1];
        k2 = surrogate_neighbours[k3 + 2];

        // Get the auxiliary triangle's vertex coordinates
        // (really the centroids of neighbouring triangles)
        coord_index = 2 * k0;
        x0 = centroid_coordinates[coord_index];
        y0 = centroid_coordinates[coord_index + 1];

        coord_index = 2 * k1;
        x1 = centroid_coordinates[coord_index];
        y1 = centroid_coordinates[coord_index + 1];

        coord_index = 2 * k2;
        x2 = centroid_coordinates[coord_index];
        y2 = centroid_coordinates[coord_index + 1];

        // Store x- and y- differentials for the vertices
        // of the auxiliary triangle
        dx1 = x1 - x0;
        dx2 = x2 - x0;
        dy1 = y1 - y0;
        dy2 = y2 - y0;

        // Calculate 2*area of the auxiliary triangle
        // The triangle is guaranteed to be counter-clockwise
        area2 = dy2 * dx1 - dy1*dx2;

        if (area2 <= 0) {
            stage_vertex_values[k3] = stage_centroid_values[k];
            stage_vertex_values[k3 + 1] = stage_centroid_values[k];
            stage_vertex_values[k3 + 2] = stage_centroid_values[k];
            xmom_vertex_values[k3] = xmom_centroid_values[k];
            xmom_vertex_values[k3 + 1] = xmom_centroid_values[k];
            xmom_vertex_values[k3 + 2] = xmom_centroid_values[k];
            ymom_vertex_values[k3] = ymom_centroid_values[k];
            ymom_vertex_values[k3 + 1] = ymom_centroid_values[k];
            ymom_vertex_values[k3 + 2] = ymom_centroid_values[k];

            //continue;
            return;
        }

        // Calculate heights of neighbouring cells
        hc = stage_centroid_values[k] - bed_centroid_values[k];
        h0 = stage_centroid_values[k0] - bed_centroid_values[k0];
        h1 = stage_centroid_values[k1] - bed_centroid_values[k1];
        h2 = stage_centroid_values[k2] - bed_centroid_values[k2];
        hmin = min(min(h0, min(h1, h2)), hc);
        //hfactor = hc/(hc + 1.0);

        hfactor = 0.0;
        if (hmin > 0.001) {
            hfactor = (hmin - 0.001) / (hmin + 0.004);
        }

        if (optimise_dry_cells) {
            // Check if linear reconstruction is necessary for triangle k
            // This check will exclude dry cells.

            hmax = max(h0, max(h1, h2));
            if (hmax < epsilon) {
                // continue;
                return;
            }
        }

        //-----------------------------------
        // stage
        //-----------------------------------

        // Calculate the difference between vertex 0 of the auxiliary
        // triangle and the centroid of triangle k
        dq0 = stage_centroid_values[k0] - stage_centroid_values[k];

        // Calculate differentials between the vertices
        // of the auxiliary triangle (centroids of neighbouring triangles)
        dq1 = stage_centroid_values[k1] - stage_centroid_values[k0];
        dq2 = stage_centroid_values[k2] - stage_centroid_values[k0];

        inv_area2 = 1.0 / area2;
        // Calculate the gradient of stage on the auxiliary triangle
        a = dy2 * dq1 - dy1*dq2;
        a *= inv_area2;
        b = dx1 * dq2 - dx2*dq1;
        b *= inv_area2;

        // Calculate provisional jumps in stage from the centroid
        // of triangle k to its vertices, to be limited
        dqv[0] = a * dxv0 + b*dyv0;
        dqv[1] = a * dxv1 + b*dyv1;
        dqv[2] = a * dxv2 + b*dyv2;

        // Now we want to find min and max of the centroid and the
        // vertices of the auxiliary triangle and compute jumps
        // from the centroid to the min and max
        find_qmin_and_qmax(dq0, dq1, dq2, &qmin, &qmax);

        // Playing with dry wet interface
        //hmin = qmin;
        //beta_tmp = beta_w_dry;
        //if (hmin>minimum_allowed_height)
        beta_tmp = beta_w_dry + (beta_w - beta_w_dry) * hfactor;


        limit_gradient(dqv, qmin, qmax, beta_tmp);

        //for (i=0;i<3;i++)
        stage_vertex_values[k3 + 0] = stage_centroid_values[k] + dqv[0];
        stage_vertex_values[k3 + 1] = stage_centroid_values[k] + dqv[1];
        stage_vertex_values[k3 + 2] = stage_centroid_values[k] + dqv[2];


        //-----------------------------------
        // xmomentum
        //-----------------------------------

        // Calculate the difference between vertex 0 of the auxiliary
        // triangle and the centroid of triangle k
        dq0 = xmom_centroid_values[k0] - xmom_centroid_values[k];

        // Calculate differentials between the vertices
        // of the auxiliary triangle
        dq1 = xmom_centroid_values[k1] - xmom_centroid_values[k0];
        dq2 = xmom_centroid_values[k2] - xmom_centroid_values[k0];

        // Calculate the gradient of xmom on the auxiliary triangle
        a = dy2 * dq1 - dy1*dq2;
        a *= inv_area2;
        b = dx1 * dq2 - dx2*dq1;
        b *= inv_area2;

        // Calculate provisional jumps in stage from the centroid
        // of triangle k to its vertices, to be limited
        dqv[0] = a * dxv0 + b*dyv0;
        dqv[1] = a * dxv1 + b*dyv1;
        dqv[2] = a * dxv2 + b*dyv2;

        // Now we want to find min and max of the centroid and the
        // vertices of the auxiliary triangle and compute jumps
        // from the centroid to the min and max
        find_qmin_and_qmax(dq0, dq1, dq2, &qmin, &qmax);
        //beta_tmp = beta_uh;
        //if (hmin<minimum_allowed_height)
        //beta_tmp = beta_uh_dry;
        beta_tmp = beta_uh_dry + (beta_uh - beta_uh_dry) * hfactor;

        // Limit the gradient
        limit_gradient(dqv, qmin, qmax, beta_tmp);

        for (i = 0; i < 3; i++) {
            xmom_vertex_values[k3 + i] = xmom_centroid_values[k] + dqv[i];
        }

        //-----------------------------------
        // ymomentum
        //-----------------------------------

        // Calculate the difference between vertex 0 of the auxiliary
        // triangle and the centroid of triangle k
        dq0 = ymom_centroid_values[k0] - ymom_centroid_values[k];

        // Calculate differentials between the vertices
        // of the auxiliary triangle
        dq1 = ymom_centroid_values[k1] - ymom_centroid_values[k0];
        dq2 = ymom_centroid_values[k2] - ymom_centroid_values[k0];

        // Calculate the gradient of xmom on the auxiliary triangle
        a = dy2 * dq1 - dy1*dq2;
        a *= inv_area2;
        b = dx1 * dq2 - dx2*dq1;
        b *= inv_area2;

        // Calculate provisional jumps in stage from the centroid
        // of triangle k to its vertices, to be limited
        dqv[0] = a * dxv0 + b*dyv0;
        dqv[1] = a * dxv1 + b*dyv1;
        dqv[2] = a * dxv2 + b*dyv2;

        // Now we want to find min and max of the centroid and the
        // vertices of the auxiliary triangle and compute jumps
        // from the centroid to the min and max
        find_qmin_and_qmax(dq0, dq1, dq2, &qmin, &qmax);

        //beta_tmp = beta_vh;
        //
        //if (hmin<minimum_allowed_height)
        //beta_tmp = beta_vh_dry;
        beta_tmp = beta_vh_dry + (beta_vh - beta_vh_dry) * hfactor;

        // Limit the gradient
        limit_gradient(dqv, qmin, qmax, beta_tmp);

        for (i = 0; i < 3; i++) {
            ymom_vertex_values[k3 + i] = ymom_centroid_values[k] + dqv[i];
        }
    }// End number_of_boundaries <=1
    else
    {
        //==============================================
        // Number of boundaries == 2
        //==============================================

        // One internal neighbour and gradient is in direction of the neighbour's centroid

        // Find the only internal neighbour (k1?)
        for (k2 = k3; k2 < k3 + 3; k2++) {
            // Find internal neighbour of triangle k
            // k2 indexes the edges of triangle k

            if (surrogate_neighbours[k2] != k) {
                break;
            }
        }



        k1 = surrogate_neighbours[k2];

        // The coordinates of the triangle are already (x,y).
        // Get centroid of the neighbour (x1,y1)
        coord_index = 2 * k1;
        x1 = centroid_coordinates[coord_index];
        y1 = centroid_coordinates[coord_index + 1];

        // Compute x- and y- distances between the centroid of
        // triangle k and that of its neighbour
        dx1 = x1 - x;
        dy1 = y1 - y;

        // Set area2 as the square of the distance
        area2 = dx1 * dx1 + dy1*dy1;

        // Set dx2=(x1-x0)/((x1-x0)^2+(y1-y0)^2)
        // and dy2=(y1-y0)/((x1-x0)^2+(y1-y0)^2) which
        // respectively correspond to the x- and y- gradients
        // of the conserved quantities
        dx2 = 1.0 / area2;
        dy2 = dx2*dy1;
        dx2 *= dx1;


        //-----------------------------------
        // stage
        //-----------------------------------

        // Compute differentials
        dq1 = stage_centroid_values[k1] - stage_centroid_values[k];

        // Calculate the gradient between the centroid of triangle k
        // and that of its neighbour
        a = dq1*dx2;
        b = dq1*dy2;

        // Calculate provisional vertex jumps, to be limited
        dqv[0] = a * dxv0 + b*dyv0;
        dqv[1] = a * dxv1 + b*dyv1;
        dqv[2] = a * dxv2 + b*dyv2;

        // Now limit the jumps
        if (dq1 >= 0.0) {
            qmin = 0.0;
            qmax = dq1;
        } else {
            qmin = dq1;
            qmax = 0.0;
        }

        // Limit the gradient
        limit_gradient(dqv, qmin, qmax, beta_w);

        //for (i=0; i < 3; i++)
        //{
        stage_vertex_values[k3] = stage_centroid_values[k] + dqv[0];
        stage_vertex_values[k3 + 1] = stage_centroid_values[k] + dqv[1];
        stage_vertex_values[k3 + 2] = stage_centroid_values[k] + dqv[2];
        //}

        //-----------------------------------
        // xmomentum
        //-----------------------------------

        // Compute differentials
        dq1 = xmom_centroid_values[k1] - xmom_centroid_values[k];

        // Calculate the gradient between the centroid of triangle k
        // and that of its neighbour
        a = dq1*dx2;
        b = dq1*dy2;

        // Calculate provisional vertex jumps, to be limited
        dqv[0] = a * dxv0 + b*dyv0;
        dqv[1] = a * dxv1 + b*dyv1;
        dqv[2] = a * dxv2 + b*dyv2;

        // Now limit the jumps
        if (dq1 >= 0.0) {
            qmin = 0.0;
            qmax = dq1;
        } else {
            qmin = dq1;
            qmax = 0.0;
        }

        // Limit the gradient
        limit_gradient(dqv, qmin, qmax, beta_w);

        //for (i=0;i<3;i++)
        //xmom_vertex_values[k3] = xmom_centroid_values[k] + dqv[0];
        //xmom_vertex_values[k3 + 1] = xmom_centroid_values[k] + dqv[1];
        //xmom_vertex_values[k3 + 2] = xmom_centroid_values[k] + dqv[2];

        for (i = 0; i < 3; i++) {
            xmom_vertex_values[k3 + i] = xmom_centroid_values[k] + dqv[i];
        }

        //-----------------------------------
        // ymomentum
        //-----------------------------------

        // Compute differentials
        dq1 = ymom_centroid_values[k1] - ymom_centroid_values[k];

        // Calculate the gradient between the centroid of triangle k
        // and that of its neighbour
        a = dq1*dx2;
        b = dq1*dy2;

        // Calculate provisional vertex jumps, to be limited
        dqv[0] = a * dxv0 + b*dyv0;
        dqv[1] = a * dxv1 + b*dyv1;
        dqv[2] = a * dxv2 + b*dyv2;

        // Now limit the jumps
        if (dq1 >= 0.0) {
            qmin = 0.0;
            qmax = dq1;
        }
        else {
            qmin = dq1;
            qmax = 0.0;
        }

        // Limit the gradient
        limit_gradient(dqv, qmin, qmax, beta_w);

        //for (i=0;i<3;i++)
        //ymom_vertex_values[k3] = ymom_centroid_values[k] + dqv[0];
        //ymom_vertex_values[k3 + 1] = ymom_centroid_values[k] + dqv[1];
        //ymom_vertex_values[k3 + 2] = ymom_centroid_values[k] + dqv[2];

        for (i = 0; i < 3; i++) {
            ymom_vertex_values[k3 + i] = ymom_centroid_values[k] + dqv[i];
        }
        //ymom_vertex_values[k3] = ymom_centroid_values[k] + dqv[0];
        //ymom_vertex_values[k3 + 1] = ymom_centroid_values[k] + dqv[1];
        //ymom_vertex_values[k3 + 2] = ymom_centroid_values[k] + dqv[2];
    } // else [number_of_boundaries==2]

}



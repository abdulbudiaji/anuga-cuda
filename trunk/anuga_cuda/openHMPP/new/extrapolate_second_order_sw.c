//#define REARRANGED_DOMAIN

#include "hmpp_fun.h"

#ifndef USING_HOST_MACRO
#define BLOCK_SIZE  8
#endif



int limit_gradient(
        double *dqv0, 
        double *dqv1, 
        double *dqv2, 
        double qmin, 
        double qmax, 
        double beta_w) 
{
#ifdef USING_ORIGINAL
    int i;
#endif
    double r = 1000.0, r0 = 1.0, phi = 1.0;
    double TINY = 1.0e-100; // to avoid machine accuracy problems.

#ifdef USING_ORIGINAL
    for (i = 0; i < 3; i++) {
        if (dqv[i]<-TINY)
            r0 = qmin / dqv[i];

        if (dqv[i] > TINY)
            r0 = qmax / dqv[i];

        r = fmin(r0, r);
    } //for (i = 0; i < 3; i++) {
#else
    // i = 0
    if (*dqv0<-TINY)
        r0 = qmin / *dqv0;

    if (*dqv0 > TINY)
        r0 = qmax / *dqv0;

    r = fmin(r0, r);
    // i = 1
    if (*dqv1<-TINY)
        r0 = qmin / *dqv1;

    if (*dqv1 > TINY)
        r0 = qmax / *dqv1;

    r = fmin(r0, r);
    // i = 2
    if (*dqv2<-TINY)
        r0 = qmin / *dqv2;

    if (*dqv2 > TINY)
        r0 = qmax / *dqv2;

    r = fmin(r0, r);
#endif 

    phi = fmin(r*beta_w, 1.0);
    //for (i=0;i<3;i++)
    *dqv0 = *dqv0 * phi;
    *dqv1 = *dqv1 * phi;
    *dqv2 = *dqv2 * phi;

    return 0;
}   



int find_qmin_and_qmax(
        double dq0, 
        double dq1, 
        double dq2,
        double *qmin, 
        double *qmax) 
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



#ifdef USING_LOCAL_DIRECTIVES
#pragma hmpp extraSndVelocity codelet, target=CUDA args[*].transfer=atcall
#endif
void extrapolate_second_order_velocity_true(
            int N,
            double minimum_allowed_height,
            double stage_centroid_values[N],
            double bed_centroid_values[N],
            double xmom_centroid_values[N],
            double xmom_centroid_store[N],
            double ymom_centroid_values[N],
            double ymom_centroid_store[N]
            )
{
    int k; 
    double dk;

    for (k=0; k < N; k++)
    {
        dk = fmax(stage_centroid_values[k] -bed_centroid_values[k], minimum_allowed_height);
        //xmom_centroid_store[k] = xmom_centroid_values[k];
        //xmom_centroid_values[k] = xmom_centroid_values[k] / dk;
        xmom_centroid_store[k] = xmom_centroid_values[k] / dk;

        //ymom_centroid_store[k] = ymom_centroid_values[k];
        //ymom_centroid_values[k] = ymom_centroid_values[k] / dk;
        ymom_centroid_store[k] = ymom_centroid_values[k] / dk;
    }
}



#ifdef USING_LOCAL_DIRECTIVES
#pragma hmpp extraSndOrderSWT codelet, target=CUDA args[*].transfer=atcall
#endif
void extrapolate_second_order_sw_true (
        int N,
        int N3,
        int N6,
        double epsilon,
        double minimum_allowed_height,
        double beta_w,
        double beta_w_dry,
        double beta_uh,
        double beta_uh_dry,
        double beta_vh,
        double beta_vh_dry,
        int optimise_dry_cells,

        long surrogate_neighbours[N3],
        long number_of_boundaries[N],
        double centroid_coordinates[N3],

        double stage_centroid_values[N],
        double bed_centroid_values[N],
        double xmom_centroid_values[N],
        double ymom_centroid_values[N],

        double vertex_coordinates[N6],
        
        double stage_vertex_values[N3],
        double bed_vertex_values[N3],
        double xmom_vertex_values[N3],
        double ymom_vertex_values[N3])
{
    int k;
#ifndef REARRANGED_DOMAIN
    int k3,
        k6,
        coord_index;
#endif

    double a, b;
    int k0, k1, k2;
    double x, y, x0, y0, x1, y1, x2, y2, xv0, yv0, xv1, yv1, xv2, yv2;
    double dx1, dx2, dy1, dy2, dxv0, dxv1, dxv2, dyv0, dyv1, dyv2, dq0, dq1, dq2, area2, inv_area2;
    //double dqv[3];
    double dqv0, dqv1, dqv2;
    double hc, h0, h1, h2, beta_tmp, hfactor;

    double dv0, dv1, dv2;


    // This procedure is done in function extrapolate_velocity_second_order_true
    //int dk;
    //if (extrapolate_velocity_second_order == 1 )
    //{
    //    for (k=0; k<N; k++)
    //    {
    //        dk = max(stage_centroid_values[k]-bed_centroid_values[k],minimum_allowed_height);
    //        xmom_centroid_store[k] = xmom_centroid_values[k];
    //        xmom_centroid_values[k] = xmom_centroid_values[k] / dk;
    //        ymom_centroid_store[k] = ymom_centroid_values[k];
    //        ymom_centroid_values[k] = ymom_centroid_values[k] / dk;
    //    }
    //}



    // Begin extrapolation routine
    for (k=0; k<N; k++) 
    {
        double qmin, qmax, hmin;
        k3 = k*3;
        k6 = k*6;
        if ( number_of_boundaries[k] == 3 )
        {
#ifndef REARRANGED_DOMAIN
            stage_vertex_values[k*3] = stage_centroid_values[k];
            stage_vertex_values[k*3 +1] = stage_centroid_values[k];
            stage_vertex_values[k*3 +2] = stage_centroid_values[k];

            xmom_vertex_values[k*3] = xmom_centroid_values[k];
            xmom_vertex_values[k*3 + 1] = xmom_centroid_values[k];
            xmom_vertex_values[k*3 + 2] = xmom_centroid_values[k];

            ymom_vertex_values[k*3] = ymom_centroid_values[k];
            ymom_vertex_values[k*3 + 1] = ymom_centroid_values[k];
            ymom_vertex_values[k*3 + 2] = ymom_centroid_values[k];
#else
            stage_vertex_values[k] = stage_centroid_values[k];
            stage_vertex_values[k +N] = stage_centroid_values[k];
            stage_vertex_values[k +2*N] = stage_centroid_values[k];

            xmom_vertex_values[k] = xmom_centroid_values[k];
            xmom_vertex_values[k + N] = xmom_centroid_values[k];
            xmom_vertex_values[k + 2*N] = xmom_centroid_values[k];

            ymom_vertex_values[k] = ymom_centroid_values[k];
            ymom_vertex_values[k + N] = ymom_centroid_values[k];
            ymom_vertex_values[k + 2*N] = ymom_centroid_values[k];
#endif

            // continue;
            //break;
        } //if ( number_of_boundaries[k] == 3 )
        else
        {
#ifndef REARRANGED_DOMAIN
            xv0 = vertex_coordinates[k6];
            yv0 = vertex_coordinates[k6 + 1];
            xv1 = vertex_coordinates[k6 + 2];
            yv1 = vertex_coordinates[k6 + 3];
            xv2 = vertex_coordinates[k6 + 4];
            yv2 = vertex_coordinates[k6 + 5];

            coord_index = 2 * k;
            x = centroid_coordinates[coord_index];
            y = centroid_coordinates[coord_index + 1];
#else
            xv0 = vertex_coordinates[k];
            yv0 = vertex_coordinates[k + N];
            xv1 = vertex_coordinates[k + 2*N];
            yv1 = vertex_coordinates[k + 3*N];
            xv2 = vertex_coordinates[k + 4*N];
            yv2 = vertex_coordinates[k + 5*N];

            x = centroid_coordinates[k];
            y = centroid_coordinates[k + N];
#endif

            dxv0 = xv0 - x;
            dxv1 = xv1 - x;
            dxv2 = xv2 - x;
            dyv0 = yv0 - y;
            dyv1 = yv1 - y;
            dyv2 = yv2 - y;
        //}

            if ( number_of_boundaries[k] <= 1 )
            {
#ifndef REARRANGED_DOMAIN
                k0 = surrogate_neighbours[k*3];
                k1 = surrogate_neighbours[k*3 + 1];
                k2 = surrogate_neighbours[k*3 + 2];

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
#else
                k0 = surrogate_neighbours[k];
                k1 = surrogate_neighbours[k + N];
                k2 = surrogate_neighbours[k + 2*N];

                // Get the auxiliary triangle's vertex coordinates
                // (really the centroids of neighbouring triangles)
                x0 = centroid_coordinates[k0];
                y0 = centroid_coordinates[k0 + N];

                x1 = centroid_coordinates[k1];
                y1 = centroid_coordinates[k1 + N];

                x2 = centroid_coordinates[k2];
                y2 = centroid_coordinates[k2 + N];
#endif

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
#ifndef REARRANGED_DOMAIN
                    stage_vertex_values[k*3] = stage_centroid_values[k];
                    stage_vertex_values[k*3 + 1] = stage_centroid_values[k];
                    stage_vertex_values[k*3 + 2] = stage_centroid_values[k];

                    xmom_vertex_values[k*3] = xmom_centroid_values[k];
                    xmom_vertex_values[k*3 + 1] = xmom_centroid_values[k];
                    xmom_vertex_values[k*3 + 2] = xmom_centroid_values[k];

                    ymom_vertex_values[k*3] = ymom_centroid_values[k];
                    ymom_vertex_values[k*3 + 1] = ymom_centroid_values[k];
                    ymom_vertex_values[k*3 + 2] = ymom_centroid_values[k];
#else
                    stage_vertex_values[k] = stage_centroid_values[k];
                    stage_vertex_values[k + N] = stage_centroid_values[k];
                    stage_vertex_values[k + 2*N] = stage_centroid_values[k];

                    xmom_vertex_values[k] = xmom_centroid_values[k];
                    xmom_vertex_values[k + N] = xmom_centroid_values[k];
                    xmom_vertex_values[k + 2*N] = xmom_centroid_values[k];

                    ymom_vertex_values[k] = ymom_centroid_values[k];
                    ymom_vertex_values[k + N] = ymom_centroid_values[k];
                    ymom_vertex_values[k + 2*N] = ymom_centroid_values[k];
#endif

                    //continue;
                    //break;
                } //if (area2 <= 0) {
                else {
                    // Calculate heights of neighbouring cells
                    hc = stage_centroid_values[k] - bed_centroid_values[k];
                    h0 = stage_centroid_values[k0] - bed_centroid_values[k0];
                    h1 = stage_centroid_values[k1] - bed_centroid_values[k1];
                    h2 = stage_centroid_values[k2] - bed_centroid_values[k2];
                    hmin = fmin( fmin(h0, fmin(h1, h2)), hc);
                    //hfactor = hc/(hc + 1.0);

                    hfactor = 0.0;
                    if (hmin > 0.001) {
                        hfactor = (hmin - 0.001) / (hmin + 0.004);
                    }

                    // To elimate the break statement
                    //if (optimise_dry_cells) {
                    //    // Check if linear reconstruction is necessary for triangle k
                    //    // This check will exclude dry cells.

                    //    hmax = fmax(h0, fmax(h1, h2));
                    //    if (hmax < epsilon) {
                    //        // continue;
                    //        break;
                    //    }
                    //}
                    if(!optimise_dry_cells ||fmax(h0,fmax(h1,h2))>=epsilon)
                    {
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
                        dqv0 = a * dxv0 + b*dyv0;
                        dqv1 = a * dxv1 + b*dyv1;
                        dqv2 = a * dxv2 + b*dyv2;

                        // Now we want to find min and max of the centroid and the
                        // vertices of the auxiliary triangle and compute jumps
                        // from the centroid to the min and max
                        find_qmin_and_qmax(dq0, dq1, dq2, &qmin, &qmax);

                        // Playing with dry wet interface
                        //hmin = qmin;
                        //beta_tmp = beta_tmp_dry;
                        //if (hmin>minimum_allowed_height)
                        beta_tmp = beta_w_dry + (beta_w - beta_w_dry) * hfactor;


                        limit_gradient(&dqv0, &dqv1, &dqv2, qmin, qmax, beta_tmp);

                        //for (i=0;i<3;i++)
#ifndef REARRANGED_DOMAIN
                        stage_vertex_values[k*3 + 0] = stage_centroid_values[k] + dqv0;
                        stage_vertex_values[k*3 + 1] = stage_centroid_values[k] + dqv1;
                        stage_vertex_values[k*3 + 2] = stage_centroid_values[k] + dqv2;
#else
                        stage_vertex_values[k] = stage_centroid_values[k] + dqv0;
                        stage_vertex_values[k + N] = stage_centroid_values[k] + dqv1;
                        stage_vertex_values[k + 2*N] = stage_centroid_values[k] + dqv2;
#endif


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
                        dqv0 = a * dxv0 + b*dyv0;
                        dqv1 = a * dxv1 + b*dyv1;
                        dqv2 = a * dxv2 + b*dyv2;

                        // Now we want to find min and max of the centroid and the
                        // vertices of the auxiliary triangle and compute jumps
                        // from the centroid to the min and max
                        find_qmin_and_qmax(dq0, dq1, dq2, &qmin, &qmax);
                        //beta_tmp = beta_uh;
                        //if (hmin<minimum_allowed_height)
                        //beta_tmp = beta_uh_dry;
                        beta_tmp = beta_uh_dry + (beta_uh - beta_uh_dry) * hfactor;

                        // Limit the gradient
                        limit_gradient(&dqv0, &dqv1, &dqv2, qmin, qmax, beta_tmp);

                        //for (i = 0; i < 3; i++) {
#ifndef REARRANGED_DOMAIN
                        xmom_vertex_values[k*3] = xmom_centroid_values[k] + dqv0;
                        xmom_vertex_values[k*3 + 1] = xmom_centroid_values[k] + dqv1;
                        xmom_vertex_values[k*3 + 2] = xmom_centroid_values[k] + dqv2;
#else
                        xmom_vertex_values[k] = xmom_centroid_values[k] + dqv0;
                        xmom_vertex_values[k + N] = xmom_centroid_values[k] + dqv1;
                        xmom_vertex_values[k + 2*N] = xmom_centroid_values[k] + dqv2;
#endif
                        //}

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
                        dqv0 = a * dxv0 + b*dyv0;
                        dqv1 = a * dxv1 + b*dyv1;
                        dqv2 = a * dxv2 + b*dyv2;

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
                        limit_gradient(&dqv0, &dqv1, &dqv2, qmin, qmax, beta_tmp);

                        //for (i = 0; i < 3; i++) {
#ifndef REARRANGED_DOMAIN
                        ymom_vertex_values[k*3 ] = ymom_centroid_values[k] + dqv0;
                        ymom_vertex_values[k*3 + 1] = ymom_centroid_values[k] + dqv1;
                        ymom_vertex_values[k*3 + 2] = ymom_centroid_values[k] + dqv2;
#else
                        ymom_vertex_values[k] = ymom_centroid_values[k] + dqv0;
                        ymom_vertex_values[k + N] = ymom_centroid_values[k] + dqv1;
                        ymom_vertex_values[k + 2*N] = ymom_centroid_values[k] + dqv2;
#endif
                        //}
                    }//if(!optimise_dry_cells ||fmax(h0,fmax(h1,h2))>=epsilon)
                } //if (area2 <= 0) {
            }// End number_of_boundaries <=1
            else
            {
                //==============================================
                // Number of boundaries == 2
                //==============================================

                // One internal neighbour and gradient is in direction of the neighbour's centroid

                // Find the only internal neighbour (k1?)
#ifndef REARRANGED_DOMAIN
                //for (k2 = k3; k2 < k3 + 3; k2++) 
                //{
                // Find internal neighbour of triangle k
                // k2 indexes the edges of triangle k

                //    if (surrogate_neighbours[k2] != k) {
                //        break;
                //    }
                //}
                k2 = k3;
                if (surrogate_neighbours[k2] == k) {
                    k2 += 1;
                    if (surrogate_neighbours[k2] == k) 
                        k2 += 1;
                }
#else
                //for (k2 = k; k2 < N; k2+=N) 
                k2 = k3;
                if (surrogate_neighbours[k2] == k) {
                    k2 += N;
                    if (surrogate_neighbours[k2] == k) 
                        k2 += N;
                }
#endif

                k1 = surrogate_neighbours[k2];

                // The coordinates of the triangle are already (x,y).
                // Get centroid of the neighbour (x1,y1)
#ifndef REARRANGED_DOMAIN
                coord_index = 2 * k1;
                x1 = centroid_coordinates[coord_index];
                y1 = centroid_coordinates[coord_index + 1];
#else
                x1 = centroid_coordinates[k1];
                y1 = centroid_coordinates[k1 + N];
#endif

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
                dqv0 = a * dxv0 + b*dyv0;
                dqv1 = a * dxv1 + b*dyv1;
                dqv2 = a * dxv2 + b*dyv2;

                // Now limit the jumps
                if (dq1 >= 0.0) {
                    qmin = 0.0;
                    qmax = dq1;
                } else {
                    qmin = dq1;
                    qmax = 0.0;
                }

                // Limit the gradient
                limit_gradient(&dqv0, &dqv1, &dqv2, qmin, qmax, beta_w);

                //for (i=0; i < 3; i++)
                //{
#ifndef REARRANGED_DOMAIN
                stage_vertex_values[k*3] = stage_centroid_values[k] + dqv0;
                stage_vertex_values[k*3 + 1] = stage_centroid_values[k] + dqv1;
                stage_vertex_values[k*3 + 2] = stage_centroid_values[k] + dqv2;
#else
                stage_vertex_values[k] = stage_centroid_values[k] + dqv0;
                stage_vertex_values[k + N] = stage_centroid_values[k] + dqv1;
                stage_vertex_values[k + 2*N] = stage_centroid_values[k] + dqv2;
#endif
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
                dqv0 = a * dxv0 + b*dyv0;
                dqv1 = a * dxv1 + b*dyv1;
                dqv2 = a * dxv2 + b*dyv2;

                // Now limit the jumps
                if (dq1 >= 0.0) {
                    qmin = 0.0;
                    qmax = dq1;
                } else {
                    qmin = dq1;
                    qmax = 0.0;
                }

                // Limit the gradient
                limit_gradient(&dqv0, &dqv1, &dqv2, qmin, qmax, beta_w);

                //for (i=0;i<3;i++)
                //    xmom_vertex_values[k3 + i] = xmom_centroid_values[k] + dqv[i];
#ifndef REARRANGED_DOMAIN
                xmom_vertex_values[k*3] = xmom_centroid_values[k] + dqv0;
                xmom_vertex_values[k*3 + 1] = xmom_centroid_values[k] + dqv1;
                xmom_vertex_values[k*3 + 2] = xmom_centroid_values[k] + dqv2;
#else
                xmom_vertex_values[k] = xmom_centroid_values[k] + dqv0;
                xmom_vertex_values[k + N] = xmom_centroid_values[k] + dqv1;
                xmom_vertex_values[k + 2*N] = xmom_centroid_values[k] + dqv2;
#endif



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
                dqv0 = a * dxv0 + b*dyv0;
                dqv1 = a * dxv1 + b*dyv1;
                dqv2 = a * dxv2 + b*dyv2;

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
                limit_gradient(&dqv0, &dqv1, &dqv2, qmin, qmax, beta_w);

                //for (i=0;i<3;i++)
#ifndef REARRANGED_DOMAIN
                ymom_vertex_values[k*3] = ymom_centroid_values[k] + dqv0;
                ymom_vertex_values[k*3 + 1] = ymom_centroid_values[k] + dqv1;
                ymom_vertex_values[k*3 + 2] = ymom_centroid_values[k] + dqv2;
#else
                ymom_vertex_values[k] = ymom_centroid_values[k] + dqv0;
                ymom_vertex_values[k + N] = ymom_centroid_values[k] + dqv1;
                ymom_vertex_values[k + 2*N] = ymom_centroid_values[k] + dqv2;
#endif
            } // else [number_of_boundaries==2]
        }// else if ( number_of_boundaries[k] == 3 )
    }




    // Convert back from velocity to momentum
    //if (extrapolate_velocity_second_order == 1) 
    //{
    for (k = 0; k < N; k++) 
    {
#ifndef REARRANGED_DOMAIN
        dv0 = fmax(stage_vertex_values[k*3] - bed_vertex_values[k*3], 0.);
        dv1 = fmax(stage_vertex_values[k*3+1] - bed_vertex_values[k*3+1], 0.);
        dv2 = fmax(stage_vertex_values[k*3+2] - bed_vertex_values[k*3+2], 0.);

        //Correct centroid and vertex values
        //xmom_centroid_values[k] = xmom_centroid_store[k];
        xmom_vertex_values[k*3] = xmom_vertex_values[k*3] * dv0;
        xmom_vertex_values[k*3+1] = xmom_vertex_values[k*3+1] * dv1;
        xmom_vertex_values[k*3+2] = xmom_vertex_values[k*3+2] * dv2;

        //ymom_centroid_values[k] = ymom_centroid_store[k];
        ymom_vertex_values[k*3] = ymom_vertex_values[k*3] * dv0;
        ymom_vertex_values[k*3+1] = ymom_vertex_values[k*3+1] * dv1;
        ymom_vertex_values[k*3+2] = ymom_vertex_values[k*3+2] * dv2;
#else
        dv0 = max(stage_vertex_values[k] - bed_vertex_values[k], 0.);
        dv1 = max(stage_vertex_values[k+N] - bed_vertex_values[k+N], 0.);
        dv2 = max(stage_vertex_values[k+2*N] - bed_vertex_values[k+2*N], 0.);

        //Correct centroid and vertex values
        //xmom_centroid_values[k] = xmom_centroid_store[k];
        xmom_vertex_values[k] = xmom_vertex_values[k] * dv0;
        xmom_vertex_values[k+N] = xmom_vertex_values[k+N] * dv1;
        xmom_vertex_values[k+2*N] = xmom_vertex_values[k+2*N] * dv2;

        //ymom_centroid_values[k] = ymom_centroid_store[k];
        ymom_vertex_values[k] = ymom_vertex_values[k] * dv0;
        ymom_vertex_values[k+N] = ymom_vertex_values[k+N] * dv1;
        ymom_vertex_values[k+2*N] = ymom_vertex_values[k+2*N] * dv2;
#endif
    }
    //}
}



#ifdef USING_LOCAL_DIRECTIVES
#pragma hmpp extraSndOrderSWF codelet, target=CUDA args[*].transfer=atcall
#endif
void extrapolate_second_order_sw_false (
        int N,
        int N3,
        int N6,
        double epsilon,
        double minimum_allowed_height,
        double beta_w,
        double beta_w_dry,
        double beta_uh,
        double beta_uh_dry,
        double beta_vh,
        double beta_vh_dry,
        int optimise_dry_cells,

        long surrogate_neighbours[N3],
        long number_of_boundaries[N],
        double centroid_coordinates[N3],

        double stage_centroid_values[N],
        double bed_centroid_values[N],
        double xmom_centroid_values[N],
        double ymom_centroid_values[N],

        double vertex_coordinates[N6],
        
        double stage_vertex_values[N3],
        double bed_vertex_values[N3],
        double xmom_vertex_values[N3],
        double ymom_vertex_values[N3])
{
    int k;
#ifndef REARRANGED_DOMAIN
    int k3,
        k6,
        coord_index;
#endif

    double a, b;
    int k0, k1, k2;
    double x, y, x0, y0, x1, y1, x2, y2, xv0, yv0, xv1, yv1, xv2, yv2;
    double dx1, dx2, dy1, dy2, dxv0, dxv1, dxv2, dyv0, dyv1, dyv2, dq0, dq1, dq2, area2, inv_area2;
    //double dqv[3];
    double dqv0, dqv1, dqv2;
    double hc, h0, h1, h2, beta_tmp, hfactor;





    // Begin extrapolation routine
    for (k=0; k<N; k++) 
    {
        double qmin, qmax, hmin;
        k3 = k*3;
        k6 = k*6;
        if ( number_of_boundaries[k] == 3 )
        {
#ifndef REARRANGED_DOMAIN
            stage_vertex_values[k*3] = stage_centroid_values[k];
            stage_vertex_values[k*3 +1] = stage_centroid_values[k];
            stage_vertex_values[k*3 +2] = stage_centroid_values[k];

            xmom_vertex_values[k*3] = xmom_centroid_values[k];
            xmom_vertex_values[k*3 + 1] = xmom_centroid_values[k];
            xmom_vertex_values[k*3 + 2] = xmom_centroid_values[k];

            ymom_vertex_values[k*3] = ymom_centroid_values[k];
            ymom_vertex_values[k*3 + 1] = ymom_centroid_values[k];
            ymom_vertex_values[k*3 + 2] = ymom_centroid_values[k];
#else
            stage_vertex_values[k] = stage_centroid_values[k];
            stage_vertex_values[k +N] = stage_centroid_values[k];
            stage_vertex_values[k +2*N] = stage_centroid_values[k];

            xmom_vertex_values[k] = xmom_centroid_values[k];
            xmom_vertex_values[k + N] = xmom_centroid_values[k];
            xmom_vertex_values[k + 2*N] = xmom_centroid_values[k];

            ymom_vertex_values[k] = ymom_centroid_values[k];
            ymom_vertex_values[k + N] = ymom_centroid_values[k];
            ymom_vertex_values[k + 2*N] = ymom_centroid_values[k];
#endif

            // continue;
            //break;
        } //if ( number_of_boundaries[k] == 3 )
        else
        {
#ifndef REARRANGED_DOMAIN
            xv0 = vertex_coordinates[k6];
            yv0 = vertex_coordinates[k6 + 1];
            xv1 = vertex_coordinates[k6 + 2];
            yv1 = vertex_coordinates[k6 + 3];
            xv2 = vertex_coordinates[k6 + 4];
            yv2 = vertex_coordinates[k6 + 5];

            coord_index = 2 * k;
            x = centroid_coordinates[coord_index];
            y = centroid_coordinates[coord_index + 1];
#else
            xv0 = vertex_coordinates[k];
            yv0 = vertex_coordinates[k + N];
            xv1 = vertex_coordinates[k + 2*N];
            yv1 = vertex_coordinates[k + 3*N];
            xv2 = vertex_coordinates[k + 4*N];
            yv2 = vertex_coordinates[k + 5*N];

            x = centroid_coordinates[k];
            y = centroid_coordinates[k + N];
#endif

            dxv0 = xv0 - x;
            dxv1 = xv1 - x;
            dxv2 = xv2 - x;
            dyv0 = yv0 - y;
            dyv1 = yv1 - y;
            dyv2 = yv2 - y;
        //}

            if ( number_of_boundaries[k] <= 1 )
            {
#ifndef REARRANGED_DOMAIN
                k0 = surrogate_neighbours[k*3];
                k1 = surrogate_neighbours[k*3 + 1];
                k2 = surrogate_neighbours[k*3 + 2];

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
#else
                k0 = surrogate_neighbours[k];
                k1 = surrogate_neighbours[k + N];
                k2 = surrogate_neighbours[k + 2*N];

                // Get the auxiliary triangle's vertex coordinates
                // (really the centroids of neighbouring triangles)
                x0 = centroid_coordinates[k0];
                y0 = centroid_coordinates[k0 + N];

                x1 = centroid_coordinates[k1];
                y1 = centroid_coordinates[k1 + N];

                x2 = centroid_coordinates[k2];
                y2 = centroid_coordinates[k2 + N];
#endif

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
#ifndef REARRANGED_DOMAIN
                    stage_vertex_values[k*3] = stage_centroid_values[k];
                    stage_vertex_values[k*3 + 1] = stage_centroid_values[k];
                    stage_vertex_values[k*3 + 2] = stage_centroid_values[k];

                    xmom_vertex_values[k*3] = xmom_centroid_values[k];
                    xmom_vertex_values[k*3 + 1] = xmom_centroid_values[k];
                    xmom_vertex_values[k*3 + 2] = xmom_centroid_values[k];

                    ymom_vertex_values[k*3] = ymom_centroid_values[k];
                    ymom_vertex_values[k*3 + 1] = ymom_centroid_values[k];
                    ymom_vertex_values[k*3 + 2] = ymom_centroid_values[k];
#else
                    stage_vertex_values[k] = stage_centroid_values[k];
                    stage_vertex_values[k + N] = stage_centroid_values[k];
                    stage_vertex_values[k + 2*N] = stage_centroid_values[k];

                    xmom_vertex_values[k] = xmom_centroid_values[k];
                    xmom_vertex_values[k + N] = xmom_centroid_values[k];
                    xmom_vertex_values[k + 2*N] = xmom_centroid_values[k];

                    ymom_vertex_values[k] = ymom_centroid_values[k];
                    ymom_vertex_values[k + N] = ymom_centroid_values[k];
                    ymom_vertex_values[k + 2*N] = ymom_centroid_values[k];
#endif

                    //continue;
                    //break;
                } //if (area2 <= 0) {
                else {
                    // Calculate heights of neighbouring cells
                    hc = stage_centroid_values[k] - bed_centroid_values[k];
                    h0 = stage_centroid_values[k0] - bed_centroid_values[k0];
                    h1 = stage_centroid_values[k1] - bed_centroid_values[k1];
                    h2 = stage_centroid_values[k2] - bed_centroid_values[k2];
                    hmin = fmin( fmin(h0, fmin(h1, h2)), hc);
                    //hfactor = hc/(hc + 1.0);

                    hfactor = 0.0;
                    if (hmin > 0.001) {
                        hfactor = (hmin - 0.001) / (hmin + 0.004);
                    }

                    // To elimate the break statement
                    //if (optimise_dry_cells) {
                    //    // Check if linear reconstruction is necessary for triangle k
                    //    // This check will exclude dry cells.

                    //    hmax = fmax(h0, fmax(h1, h2));
                    //    if (hmax < epsilon) {
                    //        // continue;
                    //        break;
                    //    }
                    //}
                    if(!optimise_dry_cells ||fmax(h0,fmax(h1,h2))>=epsilon)
                    {
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
                        dqv0 = a * dxv0 + b*dyv0;
                        dqv1 = a * dxv1 + b*dyv1;
                        dqv2 = a * dxv2 + b*dyv2;

                        // Now we want to find min and max of the centroid and the
                        // vertices of the auxiliary triangle and compute jumps
                        // from the centroid to the min and max
                        find_qmin_and_qmax(dq0, dq1, dq2, &qmin, &qmax);

                        // Playing with dry wet interface
                        //hmin = qmin;
                        //beta_tmp = beta_tmp_dry;
                        //if (hmin>minimum_allowed_height)
                        beta_tmp = beta_w_dry + (beta_w - beta_w_dry) * hfactor;


                        limit_gradient(&dqv0, &dqv1, &dqv2, qmin, qmax, beta_tmp);

                        //for (i=0;i<3;i++)
#ifndef REARRANGED_DOMAIN
                        stage_vertex_values[k*3 + 0] = stage_centroid_values[k] + dqv0;
                        stage_vertex_values[k*3 + 1] = stage_centroid_values[k] + dqv1;
                        stage_vertex_values[k*3 + 2] = stage_centroid_values[k] + dqv2;
#else
                        stage_vertex_values[k] = stage_centroid_values[k] + dqv0;
                        stage_vertex_values[k + N] = stage_centroid_values[k] + dqv1;
                        stage_vertex_values[k + 2*N] = stage_centroid_values[k] + dqv2;
#endif


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
                        dqv0 = a * dxv0 + b*dyv0;
                        dqv1 = a * dxv1 + b*dyv1;
                        dqv2 = a * dxv2 + b*dyv2;

                        // Now we want to find min and max of the centroid and the
                        // vertices of the auxiliary triangle and compute jumps
                        // from the centroid to the min and max
                        find_qmin_and_qmax(dq0, dq1, dq2, &qmin, &qmax);
                        //beta_tmp = beta_uh;
                        //if (hmin<minimum_allowed_height)
                        //beta_tmp = beta_uh_dry;
                        beta_tmp = beta_uh_dry + (beta_uh - beta_uh_dry) * hfactor;

                        // Limit the gradient
                        limit_gradient(&dqv0, &dqv1, &dqv2, qmin, qmax, beta_tmp);

                        //for (i = 0; i < 3; i++) {
#ifndef REARRANGED_DOMAIN
                        xmom_vertex_values[k*3] = xmom_centroid_values[k] + dqv0;
                        xmom_vertex_values[k*3 + 1] = xmom_centroid_values[k] + dqv1;
                        xmom_vertex_values[k*3 + 2] = xmom_centroid_values[k] + dqv2;
#else
                        xmom_vertex_values[k] = xmom_centroid_values[k] + dqv0;
                        xmom_vertex_values[k + N] = xmom_centroid_values[k] + dqv1;
                        xmom_vertex_values[k + 2*N] = xmom_centroid_values[k] + dqv2;
#endif
                        //}

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
                        dqv0 = a * dxv0 + b*dyv0;
                        dqv1 = a * dxv1 + b*dyv1;
                        dqv2 = a * dxv2 + b*dyv2;

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
                        limit_gradient(&dqv0, &dqv1, &dqv2, qmin, qmax, beta_tmp);

                        //for (i = 0; i < 3; i++) {
#ifndef REARRANGED_DOMAIN
                        ymom_vertex_values[k*3 ] = ymom_centroid_values[k] + dqv0;
                        ymom_vertex_values[k*3 + 1] = ymom_centroid_values[k] + dqv1;
                        ymom_vertex_values[k*3 + 2] = ymom_centroid_values[k] + dqv2;
#else
                        ymom_vertex_values[k] = ymom_centroid_values[k] + dqv0;
                        ymom_vertex_values[k + N] = ymom_centroid_values[k] + dqv1;
                        ymom_vertex_values[k + 2*N] = ymom_centroid_values[k] + dqv2;
#endif
                        //}
                    }//if(!optimise_dry_cells ||fmax(h0,fmax(h1,h2))>=epsilon)
                } //if (area2 <= 0) {
            }// End number_of_boundaries <=1
            else
            {
                //==============================================
                // Number of boundaries == 2
                //==============================================

                // One internal neighbour and gradient is in direction of the neighbour's centroid

                // Find the only internal neighbour (k1?)
#ifndef REARRANGED_DOMAIN
                //for (k2 = k3; k2 < k3 + 3; k2++) 
                //{
                // Find internal neighbour of triangle k
                // k2 indexes the edges of triangle k

                //    if (surrogate_neighbours[k2] != k) {
                //        break;
                //    }
                //}
                k2 = k3;
                if (surrogate_neighbours[k2] == k) {
                    k2 += 1;
                    if (surrogate_neighbours[k2] == k) 
                        k2 += 1;
                }
#else
                //for (k2 = k; k2 < N; k2+=N) 
                k2 = k3;
                if (surrogate_neighbours[k2] == k) {
                    k2 += N;
                    if (surrogate_neighbours[k2] == k) 
                        k2 += N;
                }
#endif

                k1 = surrogate_neighbours[k2];

                // The coordinates of the triangle are already (x,y).
                // Get centroid of the neighbour (x1,y1)
#ifndef REARRANGED_DOMAIN
                coord_index = 2 * k1;
                x1 = centroid_coordinates[coord_index];
                y1 = centroid_coordinates[coord_index + 1];
#else
                x1 = centroid_coordinates[k1];
                y1 = centroid_coordinates[k1 + N];
#endif

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
                dqv0 = a * dxv0 + b*dyv0;
                dqv1 = a * dxv1 + b*dyv1;
                dqv2 = a * dxv2 + b*dyv2;

                // Now limit the jumps
                if (dq1 >= 0.0) {
                    qmin = 0.0;
                    qmax = dq1;
                } else {
                    qmin = dq1;
                    qmax = 0.0;
                }

                // Limit the gradient
                limit_gradient(&dqv0, &dqv1, &dqv2, qmin, qmax, beta_w);

                //for (i=0; i < 3; i++)
                //{
#ifndef REARRANGED_DOMAIN
                stage_vertex_values[k*3] = stage_centroid_values[k] + dqv0;
                stage_vertex_values[k*3 + 1] = stage_centroid_values[k] + dqv1;
                stage_vertex_values[k*3 + 2] = stage_centroid_values[k] + dqv2;
#else
                stage_vertex_values[k] = stage_centroid_values[k] + dqv0;
                stage_vertex_values[k + N] = stage_centroid_values[k] + dqv1;
                stage_vertex_values[k + 2*N] = stage_centroid_values[k] + dqv2;
#endif
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
                dqv0 = a * dxv0 + b*dyv0;
                dqv1 = a * dxv1 + b*dyv1;
                dqv2 = a * dxv2 + b*dyv2;

                // Now limit the jumps
                if (dq1 >= 0.0) {
                    qmin = 0.0;
                    qmax = dq1;
                } else {
                    qmin = dq1;
                    qmax = 0.0;
                }

                // Limit the gradient
                limit_gradient(&dqv0, &dqv1, &dqv2, qmin, qmax, beta_w);

                //for (i=0;i<3;i++)
                //    xmom_vertex_values[k3 + i] = xmom_centroid_values[k] + dqv[i];
#ifndef REARRANGED_DOMAIN
                xmom_vertex_values[k*3] = xmom_centroid_values[k] + dqv0;
                xmom_vertex_values[k*3 + 1] = xmom_centroid_values[k] + dqv1;
                xmom_vertex_values[k*3 + 2] = xmom_centroid_values[k] + dqv2;
#else
                xmom_vertex_values[k] = xmom_centroid_values[k] + dqv0;
                xmom_vertex_values[k + N] = xmom_centroid_values[k] + dqv1;
                xmom_vertex_values[k + 2*N] = xmom_centroid_values[k] + dqv2;
#endif



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
                dqv0 = a * dxv0 + b*dyv0;
                dqv1 = a * dxv1 + b*dyv1;
                dqv2 = a * dxv2 + b*dyv2;

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
                limit_gradient(&dqv0, &dqv1, &dqv2, qmin, qmax, beta_w);

                //for (i=0;i<3;i++)
#ifndef REARRANGED_DOMAIN
                ymom_vertex_values[k*3] = ymom_centroid_values[k] + dqv0;
                ymom_vertex_values[k*3 + 1] = ymom_centroid_values[k] + dqv1;
                ymom_vertex_values[k*3 + 2] = ymom_centroid_values[k] + dqv2;
#else
                ymom_vertex_values[k] = ymom_centroid_values[k] + dqv0;
                ymom_vertex_values[k + N] = ymom_centroid_values[k] + dqv1;
                ymom_vertex_values[k + 2*N] = ymom_centroid_values[k] + dqv2;
#endif
            } // else [number_of_boundaries==2]
        }// else if ( number_of_boundaries[k] == 3 )
    }
}

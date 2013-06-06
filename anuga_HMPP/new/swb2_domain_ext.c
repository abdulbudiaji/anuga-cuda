#include "hmpp_fun.h"

// Protect against the water elevation falling below the triangle bed
#ifdef USING_LOCAL_DIRECTIVES
#pragma hmpp protectSWB2 codelet, target=CUDA args[*].transfer=atcall
#endif
void protect_swb2(
        long N,
        long N3,
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

    int k;
    double hc, bmin, bmax;
    //double mass_error=0.; 
    // This acts like minimum_allowed height, but scales with the vertical
    // distance between the bed_centroid_value and the max bed_edge_value of
    // every triangle.
    double minimum_relative_height=0.1; 
    //int mass_added=0;

    // Protect against inifintesimal and negative heights  
    //if (maximum_allowed_speed < epsilon) {
    for (k=0; k<N; k++) {
        hc = wc[k] - zc[k];
        // Definine the maximum bed edge value on triangle k.
        bmax = 0.5*fmax((zv[3*k]+zv[3*k+1]),fmax((zv[3*k+1]+zv[3*k+2]),(zv[3*k+2]+zv[3*k])));

        if (hc < fmax(minimum_relative_height*(bmax-zc[k]), minimum_allowed_height)) {

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
                bmin = fmin(zv[3*k], fmin(zv[3*k+1],zv[3*k+2])) -minimum_allowed_height;
                //bmin=zc[k]-minimum_allowed_height;
                // Minimum allowed stage = bmin
                // WARNING: ADDING MASS if wc[k]<bmin
                if(wc[k]<bmin){
                    //mass_error+=(bmin-wc[k])*areas[k];
                    //mass_added=1; //Flag to warn of added mass                
                    //printf("Adding mass to dry cell %d %f %f %f %f %f \n", k, zv[3*k], zv[3*k+1], zv[3*k+2], wc[k]- bmin, mass_error);

                    wc[k] = fmax(wc[k], bmin); 


                    // Set vertex values as well. Seems that this shouldn't be
                    // needed. However from memory this is important at the first
                    // time step, for 'dry' areas where the designated stage is
                    // less than the bed centroid value
                    wv[3*k] = fmax(wv[3*k], bmin);
                    wv[3*k+1] = fmax(wv[3*k+1], bmin);
                    wv[3*k+2] = fmax(wv[3*k+2], bmin);
                }
            }
        }
    }

    /*
       if(mass_added==1){
       printf("Cumulative mass protection: %f m^3 \n", mass_error);
       }
     */
    //return mass_error;
}


int _find_qmin_and_qmax(double dq0, double dq1, double dq2, 
        double *qmin, double *qmax){
    // Considering the centroid of an FV triangle and the vertices of its 
    // auxiliary triangle, find 
    // qmin=min(q)-qc and qmax=max(q)-qc, 
    // where min(q) and max(q) are respectively min and max over the
    // four values (at the centroid of the FV triangle and the auxiliary 
    // triangle vertices),
    // and qc is the centroid
    // dq0=q(vertex0)-q(centroid of FV triangle)
    // dq1=q(vertex1)-q(vertex0)
    // dq2=q(vertex2)-q(vertex0)

    // This is a simple implementation 
    *qmax = fmax(fmax(dq0, fmax(dq0+dq1, dq0+dq2)), 0.0) ;
    *qmin = fmin(fmin(dq0, fmin(dq0+dq1, dq0+dq2)), 0.0) ;

    return 0;
}



int _limit_gradient(double *dqv, double qmin, double qmax, double beta_w)
{
    // Given provisional jumps dqv from the FV triangle centroid to its 
    // vertices and jumps qmin (qmax) between the centroid of the FV 
    // triangle and the minimum (maximum) of the values at the centroid of 
    // the FV triangle and the auxiliary triangle vertices,
    // calculate a multiplicative factor phi by which the provisional 
    // vertex jumps are to be limited

    double r=1000.0, r0=1.0, phi=1.0;
    static double TINY = 1.0e-100; // to avoid machine accuracy problems.
    // FIXME: Perhaps use the epsilon used elsewhere.

    // Any provisional jump with magnitude < TINY does not contribute to 
    // the limiting process.

#ifdef USING_ORIGINAL
    int i;
    for (i=0;i<3;i++){
        if (dqv[i]<-TINY)
            r0=qmin/dqv[i];

        if (dqv[i]>TINY)
            r0=qmax/dqv[i];

        r = fmin(r0,r);
    }
#else
    if (dqv[0]<-TINY)
        r0=qmin/dqv[0];

    if (dqv[0]>TINY)
        r0=qmax/dqv[0];

    r = fmin(r0,r);

    if (dqv[1]<-TINY)
        r0=qmin/dqv[1];

    if (dqv[1]>TINY)
        r0=qmax/dqv[1];

    r = fmin(r0,r);

    if (dqv[2]<-TINY)
        r0=qmin/dqv[2];

    if (dqv[2]>TINY)
        r0=qmax/dqv[2];

    r = fmin(r0,r);

#endif
    phi = fmin(r*beta_w,1.0);
    //phi=1.;
    dqv[0]=dqv[0]*phi;
    dqv[1]=dqv[1]*phi;
    dqv[2]=dqv[2]*phi;

    return 0;
}



// Computational routine
#ifdef USING_LOCAL_DIRECTIVES
#pragma hmpp extraSndOrderEdge codelet, target=CUDA args[*].transfer=atcall
#endif
void extrapolate_second_order_edge_sw(
        int number_of_elements,
        int optimise_dry_cells, 
        int extrapolate_velocity_second_order, 

        double epsilon,
        double minimum_allowed_height,
        double beta_w,
        double beta_w_dry,
        double beta_uh,
        double beta_uh_dry,
        double beta_vh,
        double beta_vh_dry,

        long* surrogate_neighbours,
        long* number_of_boundaries,

        double* centroid_coordinates,

        double* stage_centroid_values,
        double* elevation_centroid_values,
        double* xmom_centroid_values,
        double* ymom_centroid_values,

        double* edge_coordinates,

        double* stage_edge_values,
        double* elevation_edge_values,
        double* xmom_edge_values,
        double* ymom_edge_values,

        double* stage_vertex_values,
        double* xmom_vertex_values,
        double* ymom_vertex_values,
        double* elevation_vertex_values,

        double* stage_centroid_store,
        double* xmom_centroid_store,
        double* ymom_centroid_store,
        double* min_elevation_edgevalue,
        double* max_elevation_edgevalue,
        int* count_wet_neighbours)
{

    // Local variables
    double a, b; // Gradient vector used to calculate edge values from centroids
    int k, k0, k1, k2, k3, k6, coord_index, i, ktmp;
    double x, y, x0, y0, x1, y1, x2, y2, xv0, yv0, xv1, yv1, xv2, yv2; // Vertices of the auxiliary triangle
    double dx1, dx2, dy1, dy2, dxv0, dxv1, dxv2, dyv0, dyv1, dyv2, dq0, dq1, dq2, area2, inv_area2;
    double dqv[3], qmin, qmax, hmin, bedmax, stagemin;
    double hc, h0, h1, h2, beta_tmp, hfactor;
    double dk, de[3];



    if(extrapolate_velocity_second_order==1){
        // Replace momentum centroid with velocity centroid to allow velocity
        // extrapolation This will be changed back at the end of the routine
        for (k=0; k<number_of_elements; k++){

            dk = fmax(stage_centroid_values[k]-elevation_centroid_values[k],minimum_allowed_height);
            xmom_centroid_store[k] = xmom_centroid_values[k];
            xmom_centroid_values[k] = xmom_centroid_values[k]/dk;

            ymom_centroid_store[k] = ymom_centroid_values[k];
            ymom_centroid_values[k] = ymom_centroid_values[k]/dk;

            min_elevation_edgevalue[k] = fmin(elevation_edge_values[3*k], 
                    fmin(elevation_edge_values[3*k+1],
                        elevation_edge_values[3*k+2]));
            max_elevation_edgevalue[k] = fmax(elevation_edge_values[3*k], 
                    fmax(elevation_edge_values[3*k+1],
                        elevation_edge_values[3*k+2]));
        }

    }

    // Count how many 'fully submerged' neighbours the cell has
    for(k=0; k<number_of_elements;k++){ 
        count_wet_neighbours[k]=0;
        for (i=0; i<3; i++){
            ktmp = surrogate_neighbours[3*k+i];              
            if(stage_centroid_values[ktmp] > max_elevation_edgevalue[ktmp]){
                count_wet_neighbours[k]+=1;
            }
        }
    }

    // Begin extrapolation routine
    for (k = 0; k < number_of_elements; k++) 
    {
        k3=k*3;
        k6=k*6;

        if (number_of_boundaries[k]==3)
            //if (0==0)
        {
            // No neighbours, set gradient on the triangle to zero

            stage_edge_values[k3]   = stage_centroid_values[k];
            stage_edge_values[k3+1] = stage_centroid_values[k];
            stage_edge_values[k3+2] = stage_centroid_values[k];
            xmom_edge_values[k3]    = xmom_centroid_values[k];
            xmom_edge_values[k3+1]  = xmom_centroid_values[k];
            xmom_edge_values[k3+2]  = xmom_centroid_values[k];
            ymom_edge_values[k3]    = ymom_centroid_values[k];
            ymom_edge_values[k3+1]  = ymom_centroid_values[k];
            ymom_edge_values[k3+2]  = ymom_centroid_values[k];

            continue;
        }
        else 
        {
            // Triangle k has one or more neighbours. 
            // Get centroid and edge coordinates of the triangle

            // Get the edge coordinates
            xv0 = edge_coordinates[k6];   
            yv0 = edge_coordinates[k6+1];
            xv1 = edge_coordinates[k6+2]; 
            yv1 = edge_coordinates[k6+3];
            xv2 = edge_coordinates[k6+4]; 
            yv2 = edge_coordinates[k6+5];

            // Get the centroid coordinates
            coord_index = 2*k;
            x = centroid_coordinates[coord_index];
            y = centroid_coordinates[coord_index+1];

            // Store x- and y- differentials for the edges of 
            // triangle k relative to the centroid
            dxv0 = xv0 - x; 
            dxv1 = xv1 - x; 
            dxv2 = xv2 - x;
            dyv0 = yv0 - y; 
            dyv1 = yv1 - y; 
            dyv2 = yv2 - y;
            // Compute the minimum distance from the centroid to an edge
            //demin=min(dxv0*dxv0 +dyv0*dyv0, min(dxv1*dxv1+dyv1*dyv1, dxv2*dxv2+dyv2*dyv2));
            //demin=sqrt(demin);
        }



        if (number_of_boundaries[k]<=1)
        {
            //==============================================
            // Number of boundaries <= 1
            //==============================================    


            // If no boundaries, auxiliary triangle is formed 
            // from the centroids of the three neighbours
            // If one boundary, auxiliary triangle is formed 
            // from this centroid and its two neighbours

            k0 = surrogate_neighbours[k3];
            k1 = surrogate_neighbours[k3 + 1];
            k2 = surrogate_neighbours[k3 + 2];

            // Take note if the max neighbour bed elevation is greater than the min
            // neighbour stage -- suggests a 'steep' bed relative to the flow
            bedmax = fmax(elevation_centroid_values[k], 
                    fmax(elevation_centroid_values[k0],
                        fmax(elevation_centroid_values[k1], 
                            elevation_centroid_values[k2])));
            //bedmax = elevation_centroid_values[k];
            stagemin = fmin( fmax(stage_centroid_values[k], 
                        elevation_centroid_values[k]), 
                    fmin( fmax(stage_centroid_values[k0], 
                            elevation_centroid_values[k0]),
                        fmin( fmax(stage_centroid_values[k1], 
                                elevation_centroid_values[k1]),
                            fmax(stage_centroid_values[k2],
                                elevation_centroid_values[k2]))));
            if(stagemin < bedmax){
                // This will cause first order extrapolation
                k2 = k;
                k0 = k;
                k1 = k;
            } 

            // Get the auxiliary triangle's vertex coordinates 
            // (really the centroids of neighbouring triangles)
            coord_index = 2*k0;
            x0 = centroid_coordinates[coord_index];
            y0 = centroid_coordinates[coord_index+1];

            coord_index = 2*k1;
            x1 = centroid_coordinates[coord_index];
            y1 = centroid_coordinates[coord_index+1];

            coord_index = 2*k2;
            x2 = centroid_coordinates[coord_index];
            y2 = centroid_coordinates[coord_index+1];

            // Store x- and y- differentials for the vertices 
            // of the auxiliary triangle
            dx1 = x1 - x0; 
            dx2 = x2 - x0;
            dy1 = y1 - y0; 
            dy2 = y2 - y0;

            // Calculate 2*area of the auxiliary triangle
            // The triangle is guaranteed to be counter-clockwise      
            area2 = dy2*dx1 - dy1*dx2; 

            // Treat triangles with zero or 1 wet neighbours. 
            if ((area2 <= 0)) //|(count_wet_neighbours[k]==0))
            {
                //printf("Error negative triangle area \n");
                //report_python_error(AT, "Negative triangle area");
                //return -1;
                stage_edge_values[k3]   = stage_centroid_values[k];
                stage_edge_values[k3+1] = stage_centroid_values[k];
                stage_edge_values[k3+2] = stage_centroid_values[k];
                // First order momentum / velocity extrapolation
                xmom_edge_values[k3]    = xmom_centroid_values[k];
                xmom_edge_values[k3+1]  = xmom_centroid_values[k];
                xmom_edge_values[k3+2]  = xmom_centroid_values[k];
                ymom_edge_values[k3]    = ymom_centroid_values[k];
                ymom_edge_values[k3+1]  = ymom_centroid_values[k];
                ymom_edge_values[k3+2]  = ymom_centroid_values[k];

                continue;
            }  

            // Calculate heights of neighbouring cells
            hc = stage_centroid_values[k]  - elevation_centroid_values[k];
            h0 = stage_centroid_values[k0] - elevation_centroid_values[k0];
            h1 = stage_centroid_values[k1] - elevation_centroid_values[k1];
            h2 = stage_centroid_values[k2] - elevation_centroid_values[k2];
            hmin = fmin( fmin(h0, fmin(h1, h2)), hc);

            hfactor = 0.0;
            //if (hmin > 0.001) 
            if (hmin > 0.) 
                //if (hc>0.0)
            {
                hfactor = 1.0 ;//hmin/(hmin + 0.004);
                //hfactor=hmin/(hmin + 0.004);
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

            inv_area2 = 1.0/area2;
            // Calculate the gradient of stage on the auxiliary triangle
            a = dy2*dq1 - dy1*dq2;
            a *= inv_area2;
            b = dx1*dq2 - dx2*dq1;
            b *= inv_area2;

            // Calculate provisional jumps in stage from the centroid 
            // of triangle k to its vertices, to be limited
            dqv[0] = a*dxv0 + b*dyv0;
            dqv[1] = a*dxv1 + b*dyv1;
            dqv[2] = a*dxv2 + b*dyv2;

            // Now we want to find min and max of the centroid and the 
            // vertices of the auxiliary triangle and compute jumps 
            // from the centroid to the min and max
            _find_qmin_and_qmax(dq0, dq1, dq2, &qmin, &qmax);

            beta_tmp = beta_w_dry + (beta_w - beta_w_dry) * hfactor;


            // Limit the gradient
            _limit_gradient(dqv, qmin, qmax, beta_tmp);
            stage_edge_values[k3+0] = stage_centroid_values[k] + dqv[0];
            stage_edge_values[k3+1] = stage_centroid_values[k] + dqv[1];
            stage_edge_values[k3+2] = stage_centroid_values[k] + dqv[2];


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
            a = dy2*dq1 - dy1*dq2;
            a *= inv_area2;
            b = dx1*dq2 - dx2*dq1;
            b *= inv_area2;

            // Calculate provisional jumps in stage from the centroid 
            // of triangle k to its vertices, to be limited      
            dqv[0] = a*dxv0+b*dyv0;
            dqv[1] = a*dxv1+b*dyv1;
            dqv[2] = a*dxv2+b*dyv2;

            // Now we want to find min and max of the centroid and the 
            // vertices of the auxiliary triangle and compute jumps 
            // from the centroid to the min and max
            //
            _find_qmin_and_qmax(dq0, dq1, dq2, &qmin, &qmax);

            beta_tmp = beta_uh_dry + (beta_uh - beta_uh_dry) * hfactor;

            _limit_gradient(dqv, qmin, qmax, beta_tmp);


            for (i=0; i < 3; i++)
            {
                xmom_edge_values[k3+i] = xmom_centroid_values[k] + dqv[i];
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
            a = dy2*dq1 - dy1*dq2;
            a *= inv_area2;
            b = dx1*dq2 - dx2*dq1;
            b *= inv_area2;

            // Calculate provisional jumps in stage from the centroid 
            // of triangle k to its vertices, to be limited
            dqv[0] = a*dxv0 + b*dyv0;
            dqv[1] = a*dxv1 + b*dyv1;
            dqv[2] = a*dxv2 + b*dyv2;

            // Now we want to find min and max of the centroid and the 
            // vertices of the auxiliary triangle and compute jumps 
            // from the centroid to the min and max
            //
            _find_qmin_and_qmax(dq0, dq1, dq2, &qmin, &qmax);

            beta_tmp = beta_vh_dry + (beta_vh - beta_vh_dry) * hfactor;   

            _limit_gradient(dqv, qmin, qmax, beta_tmp);

            for (i=0;i<3;i++)
            {
                ymom_edge_values[k3 + i] = ymom_centroid_values[k] + dqv[i];
            }

        } // End number_of_boundaries <=1 
        else
        {

            //==============================================
            // Number of boundaries == 2
            //==============================================        

            // One internal neighbour and gradient is in direction of the neighbour's centroid

            // Find the only internal neighbour (k1?)
            for (k2 = k3; k2 < k3 + 3; k2++)
            {
                // Find internal neighbour of triangle k      
                // k2 indexes the edges of triangle k   

                if (surrogate_neighbours[k2] != k)
                {
                    break;
                }
            }

            if ((k2 == k3 + 3)) 
            {
                // If we didn't find an internal neighbour
                //report_python_error(AT, "Internal neighbour not found");
                //return -1;
            }

            k1 = surrogate_neighbours[k2];

            // The coordinates of the triangle are already (x,y). 
            // Get centroid of the neighbour (x1,y1)
            coord_index = 2*k1;
            x1 = centroid_coordinates[coord_index];
            y1 = centroid_coordinates[coord_index + 1];

            // Compute x- and y- distances between the centroid of 
            // triangle k and that of its neighbour
            dx1 = x1 - x; 
            dy1 = y1 - y;

            // Set area2 as the square of the distance
            area2 = dx1*dx1 + dy1*dy1;

            // Set dx2=(x1-x0)/((x1-x0)^2+(y1-y0)^2) 
            // and dy2=(y1-y0)/((x1-x0)^2+(y1-y0)^2) which
            // respectively correspond to the x- and y- gradients 
            // of the conserved quantities
            dx2 = 1.0/area2;
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

            // Calculate provisional edge jumps, to be limited
            dqv[0] = a*dxv0 + b*dyv0;
            dqv[1] = a*dxv1 + b*dyv1;
            dqv[2] = a*dxv2 + b*dyv2;

            // Now limit the jumps
            if (dq1>=0.0)
            {
                qmin=0.0;
                qmax=dq1;
            }
            else
            {
                qmin = dq1;
                qmax = 0.0;
            }

            // Limit the gradient
            _limit_gradient(dqv, qmin, qmax, beta_w);

            //for (i=0; i < 3; i++)
            //{
            stage_edge_values[k3] = stage_centroid_values[k] + dqv[0];
            stage_edge_values[k3 + 1] = stage_centroid_values[k] + dqv[1];
            stage_edge_values[k3 + 2] = stage_centroid_values[k] + dqv[2];
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

            // Calculate provisional edge jumps, to be limited
            dqv[0] = a*dxv0+b*dyv0;
            dqv[1] = a*dxv1+b*dyv1;
            dqv[2] = a*dxv2+b*dyv2;

            // Now limit the jumps
            if (dq1 >= 0.0)
            {
                qmin = 0.0;
                qmax = dq1;
            }
            else
            {
                qmin = dq1;
                qmax = 0.0;
            }

            // Limit the gradient      
            _limit_gradient(dqv, qmin, qmax, beta_w);

            //for (i=0;i<3;i++)
            //xmom_edge_values[k3] = xmom_centroid_values[k] + dqv[0];
            //xmom_edge_values[k3 + 1] = xmom_centroid_values[k] + dqv[1];
            //xmom_edge_values[k3 + 2] = xmom_centroid_values[k] + dqv[2];

            for (i = 0; i < 3;i++)
            {
                xmom_edge_values[k3 + i] = xmom_centroid_values[k] + dqv[i];
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

            // Calculate provisional edge jumps, to be limited
            dqv[0] = a*dxv0 + b*dyv0;
            dqv[1] = a*dxv1 + b*dyv1;
            dqv[2] = a*dxv2 + b*dyv2;

            // Now limit the jumps
            if (dq1>=0.0)
            {
                qmin = 0.0;
                qmax = dq1;
            } 
            else
            {
                qmin = dq1;
                qmax = 0.0;
            }

            // Limit the gradient
            _limit_gradient(dqv, qmin, qmax, beta_w);

            for (i=0;i<3;i++)
            {
                ymom_edge_values[k3 + i] = ymom_centroid_values[k] + dqv[i];
            }
        } // else [number_of_boundaries==2]
    } // for k=0 to number_of_elements-1

    // Compute vertex values of quantities
    for (k=0; k<number_of_elements; k++){
        k3=3*k;

        // Compute stage vertex values 
        stage_vertex_values[k3] = stage_edge_values[k3+1] + stage_edge_values[k3+2] -stage_edge_values[k3] ; 
        stage_vertex_values[k3+1] =  stage_edge_values[k3] + stage_edge_values[k3+2]-stage_edge_values[k3+1]; 
        stage_vertex_values[k3+2] =  stage_edge_values[k3] + stage_edge_values[k3+1]-stage_edge_values[k3+2]; 

        // Compute xmom vertex values 
        xmom_vertex_values[k3] = xmom_edge_values[k3+1] + xmom_edge_values[k3+2] -xmom_edge_values[k3] ; 
        xmom_vertex_values[k3+1] =  xmom_edge_values[k3] + xmom_edge_values[k3+2]-xmom_edge_values[k3+1]; 
        xmom_vertex_values[k3+2] =  xmom_edge_values[k3] + xmom_edge_values[k3+1]-xmom_edge_values[k3+2]; 

        // Compute ymom vertex values 
        ymom_vertex_values[k3] = ymom_edge_values[k3+1] + ymom_edge_values[k3+2] -ymom_edge_values[k3] ; 
        ymom_vertex_values[k3+1] =  ymom_edge_values[k3] + ymom_edge_values[k3+2]-ymom_edge_values[k3+1]; 
        ymom_vertex_values[k3+2] =  ymom_edge_values[k3] + ymom_edge_values[k3+1]-ymom_edge_values[k3+2]; 

        // If needed, convert from velocity to momenta
        if(extrapolate_velocity_second_order==1){
            //Convert velocity back to momenta at centroids
            xmom_centroid_values[k] = xmom_centroid_store[k];
            ymom_centroid_values[k] = ymom_centroid_store[k];

            // Re-compute momenta at edges
            for (i=0; i<3; i++){
                de[i] = fmax(stage_edge_values[k3+i]-elevation_edge_values[k3+i],0.0);
                xmom_edge_values[k3+i]=xmom_edge_values[k3+i]*de[i];
                ymom_edge_values[k3+i]=ymom_edge_values[k3+i]*de[i];
            }

            // Re-compute momenta at vertices
            for (i=0; i<3; i++){
                de[i] = fmax(stage_vertex_values[k3+i]-elevation_vertex_values[k3+i],0.0);
                xmom_vertex_values[k3+i]=xmom_vertex_values[k3+i]*de[i];
                ymom_vertex_values[k3+i]=ymom_vertex_values[k3+i]*de[i];
            }

        }


    } 
}           

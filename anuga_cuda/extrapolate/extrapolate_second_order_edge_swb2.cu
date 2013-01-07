// Computational routine
__global__ void  _extrapolate_second_order_edge_sw(
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
        double* elevation_vertex_values,
        double* xmom_vertex_values,
        double* ymom_vertex_values,

        double * stage_centroid_store,
        double * xmom_centroid_store,
        double * ymom_centroid_store,
        double * min_elevation_edgevalue,
        double * max_elevation_edgevalue,
        int * count_wet_neighbours
        ) 
{
    const int k = 
            threadIdx.x+threadIdx.y*blockDim.x+
            (blockIdx.x+blockIdx.y*gridDim.x)*blockDim.x*blockDim.y;

    // Local variables
    double a, b; // Gradient vector used to calculate edge values from centroids
    int k0, k1, k2, k3, k6, coord_index, i, ii, ktmp;
    double x, y, x0, y0, x1, y1, x2, y2, xv0, yv0, xv1, yv1, xv2, yv2; // Vertices of the auxiliary triangle
    double dx1, dx2, dy1, dy2, dxv0, dxv1, dxv2, dyv0, dyv1, dyv2, dq0, dq1, dq2, area2, inv_area2;
    double dqv[3], qmin, qmax, hmin, hmax, bedmax, stagemin;
    double hc, h0, h1, h2, beta_tmp, hfactor, xtmp, ytmp;
    double dk, dv0, dv1, dv2, de[3], demin, dcmax, r0scale, vect_norm, l1, l2;

    //double *xmom_centroid_store, *ymom_centroid_store, *stage_centroid_store, *min_elevation_edgevalue, *max_elevation_edgevalue;
    //int *count_wet_neighbours;

    // Use malloc to avoid putting these variables on the stack, which can cause
    // segfaults in large model runs
    //xmom_centroid_store = malloc(number_of_elements*sizeof(double));
    //ymom_centroid_store = malloc(number_of_elements*sizeof(double));
    //stage_centroid_store = malloc(number_of_elements*sizeof(double));
    //min_elevation_edgevalue = malloc(number_of_elements*sizeof(double));
    //max_elevation_edgevalue = malloc(number_of_elements*sizeof(double));
    //count_wet_neighbours = malloc(number_of_elements*sizeof(int));

    if(extrapolate_velocity_second_order==1){
        // Replace momentum centroid with velocity centroid to allow velocity
        // extrapolation This will be changed back at the end of the routine
        //for (k=0; k<number_of_elements; k++){

            dk = max(stage_centroid_values[k]-elevation_centroid_values[k],minimum_allowed_height);
            xmom_centroid_store[k] = xmom_centroid_values[k];
            xmom_centroid_values[k] = xmom_centroid_values[k]/dk;

            ymom_centroid_store[k] = ymom_centroid_values[k];
            ymom_centroid_values[k] = ymom_centroid_values[k]/dk;

            min_elevation_edgevalue[k] = min(elevation_edge_values[3*k], 
                    min(elevation_edge_values[3*k+1],
                        elevation_edge_values[3*k+2]));
            max_elevation_edgevalue[k] = max(elevation_edge_values[3*k], 
                    max(elevation_edge_values[3*k+1],
                        elevation_edge_values[3*k+2]));
        //} // for k=0 to number_of_elements-1
    }

    // Count how many 'fully submerged' neighbours the cell has
    //for(k=0; k<number_of_elements;k++){ 
        count_wet_neighbours[k]=0;
        for (i=0; i<3; i++){
            ktmp = surrogate_neighbours[3*k+i];              
            if(stage_centroid_values[ktmp] > max_elevation_edgevalue[ktmp]){
                count_wet_neighbours[k]+=1;
            }
        }
    //} // for k=0 to number_of_elements-1

    // Begin extrapolation routine
    //for (k = 0; k < number_of_elements; k++) {
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
            bedmax = max(elevation_centroid_values[k], 
                    max(elevation_centroid_values[k0],
                        max(elevation_centroid_values[k1], elevation_centroid_values[k2])));
            //bedmax = elevation_centroid_values[k];
            stagemin = min(max(stage_centroid_values[k], elevation_centroid_values[k]), 
                    min(max(stage_centroid_values[k0], elevation_centroid_values[k0]),
                        min(max(stage_centroid_values[k1], elevation_centroid_values[k1]),
                            max(stage_centroid_values[k2], elevation_centroid_values[k2]))));
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
            hmin = min(min(h0, min(h1, h2)), hc);

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
            find_qmin_and_qmax(dq0, dq1, dq2, &qmin, &qmax);

            beta_tmp = beta_w_dry + (beta_w - beta_w_dry) * hfactor;


            // Limit the gradient
            limit_gradient(dqv, qmin, qmax, beta_tmp);
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
            find_qmin_and_qmax(dq0, dq1, dq2, &qmin, &qmax);

            beta_tmp = beta_uh_dry + (beta_uh - beta_uh_dry) * hfactor;

            limit_gradient(dqv, qmin, qmax, beta_tmp);


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
            find_qmin_and_qmax(dq0, dq1, dq2, &qmin, &qmax);

            beta_tmp = beta_vh_dry + (beta_vh - beta_vh_dry) * hfactor;   

            limit_gradient(dqv, qmin, qmax, beta_tmp);

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
                report_python_error(AT, "Internal neighbour not found");
                return -1;
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
            limit_gradient(dqv, qmin, qmax, beta_w);

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
            limit_gradient(dqv, qmin, qmax, beta_w);

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
            // and that of its neighbours      a = dq1*dx2;
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
            limit_gradient(dqv, qmin, qmax, beta_w);

            for (i=0;i<3;i++)
            {
                ymom_edge_values[k3 + i] = ymom_centroid_values[k] + dqv[i];
            }
        } // else [number_of_boundaries==2]
    //} // for k=0 to number_of_elements-1

    // Compute vertex values of quantities
    //for (k=0; k<number_of_elements; k++){
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
                de[i] = max(stage_edge_values[k3+i]-elevation_edge_values[k3+i],0.0);
                xmom_edge_values[k3+i]=xmom_edge_values[k3+i]*de[i];
                ymom_edge_values[k3+i]=ymom_edge_values[k3+i]*de[i];
            }

            // Re-compute momenta at vertices
            for (i=0; i<3; i++){
                de[i] = max(stage_vertex_values[k3+i]-elevation_vertex_values[k3+i],0.0);
                xmom_vertex_values[k3+i]=xmom_vertex_values[k3+i]*de[i];
                ymom_vertex_values[k3+i]=ymom_vertex_values[k3+i]*de[i];
            }

        }
    //} // for k=0 to number_of_elements-1

    return 0;
}  

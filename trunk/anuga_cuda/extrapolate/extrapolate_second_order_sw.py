#!/usr/bin/env python

def extrapolate_second_order_sw_cuda_TRUE_second_order(domain=None):
    N = domain.number_of_elements
    W = 16
    
    elements = numpy.random.randn(9)
    elements = elements.astype(numpy.float64)
    elements[0] = domain.epsilon
    elements[1] = domain.minimum_allowed_height
    elements[2] = domain.beta_w
    elements[3] = domain.beta_w_dry
    elements[4] = domain.beta_uh
    elements[5] = domain.beta_uh_dry
    elements[6] = domain.beta_vh
    elements[7] = domain.beta_vh_dry
    elements[8] = domain.optimise_dry_cells

    mod = SourceModule("""
        __global__ void extrapolate_second_order_sw(
            	double* elements,
		        long* surrogate_neighbours,
        		long* number_of_boundaries,
		        double* centroid_coordinates,
        		double* stage_centroid_values,
				double* bed_centroid_values,
		        double* xmom_centroid_values,
        		double* ymom_centroid_values,
		        double* elevation_centroid_values,
        		double* vertex_coordinates,
		        double* stage_vertex_values,
        		double* xmom_vertex_values,
		        double* ymom_vertex_values,
        		double* elevation_vertex_values)
		{
            int k = threadIdx.x + threadIdx.y;
            int k3 = 3*k,
                k6 = 6*k;

            double a, b;
            double x, y, x0, y0, x1, y1, x2, y2, xv0, yv0, xv1, yv1, xv2, yv2;
            double dx1, dx2, dy1, dy2, 
            		dxv0, dxv1, dxv2, 
            		dyv0, dyv1, dyv2, 
            		dq0, dq1, dq2, 
            		area2, inv_area2;
            double dqv[3], qmin, qmax, hmin, hmax;
            double hc, h0, h1, h2, beta_tmp, hfactor;
            
            
            // extrapolate_velocity_second_order == 1 
                
            dk = max(stage_centroid_values[k], - bed_centroid_values[k], elements[1]);
            xmom_centroid_store[k] = xmom_centroid_values[k];
            xmom_centroid_values[k] = xmom_centroid_values / dk;

            ymom_centroid_store[k] = ymom_centroid_values[k];
            ymom_centroid_values[k] = ymom_centroid_values[k] / dk;
            
            
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

                if (elements[8]) {
                    // Check if linear reconstruction is necessary for triangle k
                    // This check will exclude dry cells.

                    hmax = max(h0, max(h1, h2));
                    if (hmax < elements[0]) {
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
                //beta_tmp = elements[3];
                //if (hmin>elements[1])
                beta_tmp = elements[3] + (elements[2] - elements[3]) * hfactor;

                //printf("min_alled_height = %f\n",elements[1]);
                //printf("hmin = %f\n",hmin);
                //printf("elements[2] = %f\n",elements[2]);
                //printf("beta_tmp = %f\n",beta_tmp);
                // Limit the gradient
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
                //beta_tmp = elements[4];
                //if (hmin<elements[1])
                //beta_tmp = elements[5];
                beta_tmp = elements[5] + (elements[4] - elements[5]) * hfactor;

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

                //beta_tmp = elements[6];
                //
                //if (hmin<elements[1])
                //beta_tmp = elements[7];
                beta_tmp = elements[7] + (elements[6] - elements[7]) * hfactor;
    
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

                if ((k2 == k3 + 3)) {
                    // If we didn't find an internal neighbour
                    report_python_error(AT, "Internal neighbour not found");
                    return -1;
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
                limit_gradient(dqv, qmin, qmax, elements[2]);

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
                limit_gradient(dqv, qmin, qmax, elements[2]);

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
                limit_gradient(dqv, qmin, qmax, elements[2]);

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
            
            //if (extrapolate_velocity_second_order == 1) 
            // Convert back from velocity to momentum
            //for (k = 0; k < number_of_elements; k++) 
            k3 = 3 * k;
            //dv0 = max(stage_vertex_values[k3]-bed_vertex_values[k3],elements[1]);
            //dv1 = max(stage_vertex_values[k3+1]-bed_vertex_values[k3+1],elements[1]);
            //dv2 = max(stage_vertex_values[k3+2]-bed_vertex_values[k3+2],elements[1]);
            dv0 = max(stage_vertex_values[k3] - bed_vertex_values[k3], 0.);
            dv1 = max(stage_vertex_values[k3 + 1] - bed_vertex_values[k3 + 1], 0.);
            dv2 = max(stage_vertex_values[k3 + 2] - bed_vertex_values[k3 + 2], 0.);

            //Correct centroid and vertex values
            xmom_centroid_values[k] = xmom_centroid_store[k];
            xmom_vertex_values[k3] = xmom_vertex_values[k3] * dv0;
            xmom_vertex_values[k3 + 1] = xmom_vertex_values[k3 + 1] * dv1;
            xmom_vertex_values[k3 + 2] = xmom_vertex_values[k3 + 2] * dv2;

            ymom_centroid_values[k] = ymom_centroid_store[k];
            ymom_vertex_values[k3] = ymom_vertex_values[k3] * dv0;
            ymom_vertex_values[k3 + 1] = ymom_vertex_values[k3 + 1] * dv1;
            ymom_vertex_values[k3 + 2] = ymom_vertex_values[k3 + 2] * dv2;
                    
    
        }
    """)

    extro_func = mod.get_function("extrapolate_second_order_sw")
    extro_func( \
    		cudaInOut( elements ), \
    		cudaInOut( domain.surrogate_neighbours ), \
    		cudaInOut( domain.number_of_boundaries ), \
    		cudaInOut( domain.centroid_coordinates ), \
    		cudaInOut( domain.quantities['stage'].centroid_values ), \
    		cudaInOut( domain.quantities['xmomentum'].centroid_values ), \
    		cudaInOut( domain.quantities['ymomentum'].centroid_values ), \
    		cudaInOut( domain.quantities['elevation'].centroid_values ), \
    		cudaInOut( domain.vertex_coordinates ), \
    		cudaInOut( domain.quantities['stage'].vertex_values ), \
    		cudaInOut( domain.quantities['xmomentum'].vertex_values ), \
    		cudaInOut( domain.quantities['ymomentum'].vertex_values ), \
    		cudaInOut( domain.quantities['elevation'].vertex_values ), \
            block = ( W, W, 1),\
            grid = ( (N + W*W -1 ) / (W*W), 1) 
    		)	

def extrapolate_second_order_sw_cuda_FALSE_second_order(domain=None):
    N = domain.number_of_elements
    W = 16
    
    elements = numpy.random.randn(9)
    elements = elements.astype(numpy.float64)
    elements[0] = domain.epsilon
    elements[1] = domain.minimum_allowed_height
    elements[2] = domain.beta_w
    elements[3] = domain.beta_w_dry
    elements[4] = domain.beta_uh
    elements[5] = domain.beta_uh_dry
    elements[6] = domain.beta_vh
    elements[7] = domain.beta_vh_dry
    elements[8] = domain.optimise_dry_cells
	
    

    mod = SourceModule("""
        __global__ void extrapolate_second_order_sw(
            	double *elements,
		        long* surrogate_neighbours,
        		long* number_of_boundaries,
		        double* centroid_coordinates,
        		double* stage_centroid_values,
		        double* xmom_centroid_values,
        		double* ymom_centroid_values,
		        double* elevation_centroid_values,
        		double* vertex_coordinates,
		        double* stage_vertex_values,
        		double* xmom_vertex_values,
		        double* ymom_vertex_values,
        		double* elevation_vertex_values)
        {
            int k = threadIdx.x + threadIdx.y;
                k3 = 3*k,
                k6 = 6*k;

            double a, b;
            double x, y, x0, y0, x1, y1, x2, y2, xv0, yv0, xv1, yv1, xv2, yv2;
            double dx1, dx2, dy1, dy2, 
            		dxv0, dxv1, dxv2, 
            		dyv0, dyv1, dyv2, 
            		dq0, dq1, dq2, 
            		area2, inv_area2;
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

                if (elements[8]) {
                    // Check if linear reconstruction is necessary for triangle k
                    // This check will exclude dry cells.

                    hmax = max(h0, max(h1, h2));
                    if (hmax < elements[0]) {
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
                //beta_tmp = elements[3];
                //if (hmin>elements[1])
                beta_tmp = elements[3] + (elements[2] - elements[3]) * hfactor;

                //printf("min_alled_height = %f\n",elements[1]);
                //printf("hmin = %f\n",hmin);
                //printf("elements[2] = %f\n",elements[2]);
                //printf("beta_tmp = %f\n",beta_tmp);
                // Limit the gradient
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
                //beta_tmp = elements[4];
                //if (hmin<elements[1])
                //beta_tmp = elements[5];
                beta_tmp = elements[5] + (elements[4] - elements[5]) * hfactor;

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

                //beta_tmp = elements[6];
                //
                //if (hmin<elements[1])
                //beta_tmp = elements[7];
                beta_tmp = elements[7] + (elements[6] - elements[7]) * hfactor;
    
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

                if ((k2 == k3 + 3)) {
                    // If we didn't find an internal neighbour
                    report_python_error(AT, "Internal neighbour not found");
                    return -1;
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
                limit_gradient(dqv, qmin, qmax, elements[2]);

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
                limit_gradient(dqv, qmin, qmax, elements[2]);

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
                limit_gradient(dqv, qmin, qmax, elements[2]);

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
    """)


    extro_func = mod.get_function("extrapolate_second_order_sw")
    extro_func( \
    		cudaInOut( elements ), \
    		cudaInOut( domain.surrogate_neighbours ), \
    		cudaInOut( domain.number_of_boundaries ), \
    		cudaInOut( domain.centroid_coordinates ), \
    		cudaInOut( domain.quantities['stage'].centroid_values ), \
    		cudaInOut( domain.quantities['xmomentum'].centroid_values ), \
    		cudaInOut( domain.quantities['ymomentum'].centroid_values ), \
    		cudaInOut( domain.quantities['elevation'].centroid_values ), \
    		cudaInOut( domain.vertex_coordinates ), \
    		cudaInOut( domain.quantities['stage'].vertex_values ), \
    		cudaInOut( domain.quantities['xmomentum'].vertex_values ), \
    		cudaInOut( domain.quantities['ymomentum'].vertex_values ), \
    		cudaInOut( domain.quantities['elevation'].vertex_values ), \
            block = ( W, W, 1),\
            grid = ( (N + W*W -1 ) / (W*W), 1) 
    		)


if __name__ == '__main__':
    from anuga_cuda.merimbula_data.generate_domain import domain_create

    domain = domain_create()

    import pycuda.driver as cuda
    import pycuda.autoinit
    from pycuda.compiler import SourceModule
    import numpy

    #if ( domain.extrapolate_velocity_second_order == 1):
    print "-----extrapolate velocity second order == 1-------"
    extrapolate_second_order_sw_cuda_TRUE_second_order(domain)
    #else:
    print "-----extrapolate velocity second order !! 1-------"
    extrapolate_second_order_sw_cuda_FALSE_second_order(domain)
    
    
    
    
    
    
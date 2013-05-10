#!/usr/bin/env python


"""
    This var is used for testing the sequential functions
"""
compute_fluxes_central_call = 1
compute_fluxes_central_structure_call = 1


def compute_fluxes_central_cuda(domain):
    print "now in the compute_fluxes_central function"
    N = domain.quantities['stage'].edge_values.shape[0]
    W = 16

    global compute_fluxes_central_call
    compute_fluxes_central_call += 1


    mod = SourceModule("""
    #define Dtimestep 0
    #define Depsilon 1
    #define DH0 2
    #define Dg 3
    #define Doptimise_dry_cells 4
    #define Dcall 5

    __device__ int _rotate(double *q, double n1, double n2) {
        /*Rotate the momentum component q (q[1], q[2])
          from x,y coordinates to coordinates based on normal vector (n1, n2).
    
          Result is returned in array 3x1 r
          To rotate in opposite direction, call rotate with (q, n1, -n2)
    
          Contents of q are changed by this function */
    
    
        double q1, q2;
    
        // Shorthands
        q1 = q[1]; // uh momentum
        q2 = q[2]; // vh momentum
    
        // Rotate
        q[1] = n1 * q1 + n2*q2;
        q[2] = -n2 * q1 + n1*q2;
    
        return 0;
    }
    __device__ double _compute_speed(double *uh,
            double *h,
            double epsilon,
            double h0,
            double limiting_threshold) {
    
        double u;
    
        if (*h < limiting_threshold) {
            // Apply limiting of speeds according to the ANUGA manual
            if (*h < epsilon) {
                *h = 0.0; // Could have been negative
                u = 0.0;
            } else {
                u = *uh / (*h + h0 / *h);
            }
    
    
            // Adjust momentum to be consistent with speed
            *uh = u * *h;
        } else {
            // We are in deep water - no need for limiting
            u = *uh / *h;
        }
    
        return u;
    }
    
    
    __device__ int _flux_function_central(double *q_left, double *q_right,
            double z_left, double z_right,
            double n1, double n2,
            double epsilon,
            double h0,
            double limiting_threshold,
            double g,
            double *edgeflux, double *max_speed) 
            {
    
        /*Compute fluxes between volumes for the shallow water wave equation
          cast in terms of the 'stage', w = h+z using
          the 'central scheme' as described in
    
          Kurganov, Noelle, Petrova. 'Semidiscrete Central-Upwind Schemes For
          Hyperbolic Conservation Laws and Hamilton-Jacobi Equations'.
          Siam J. Sci. Comput. Vol. 23, No. 3, pp. 707-740.
    
          The implemented formula is given in equation (3.15) on page 714
         */
    
        int i;
    
        double w_left, h_left, uh_left, vh_left, u_left;
        double w_right, h_right, uh_right, vh_right, u_right;
        double s_min, s_max, soundspeed_left, soundspeed_right;
        double denom, inverse_denominator, z;
    
        // Workspace (allocate once, use many)
        double q_left_rotated[3], q_right_rotated[3], flux_right[3], flux_left[3];
    
        // Copy conserved quantities to protect from modification
        q_left_rotated[0] = q_left[0];
        q_right_rotated[0] = q_right[0];
        q_left_rotated[1] = q_left[1];
        q_right_rotated[1] = q_right[1];
        q_left_rotated[2] = q_left[2];
        q_right_rotated[2] = q_right[2];
    
        // Align x- and y-momentum with x-axis
        _rotate(q_left_rotated, n1, n2);
        _rotate(q_right_rotated, n1, n2);
    
    
        if (fabs(z_left - z_right) > 1.0e-10) {
            //report_python_error(AT, "Discontinuous Elevation");
            return 0.0;
        }
        z = 0.5 * (z_left + z_right); // Average elevation values.
        // Even though this will nominally allow
        // for discontinuities in the elevation data,
        // there is currently no numerical support for
        // this so results may be strange near
        // jumps in the bed.
    
        // Compute speeds in x-direction
        w_left = q_left_rotated[0];
        h_left = w_left - z;
        uh_left = q_left_rotated[1];
        u_left = _compute_speed(&uh_left, &h_left,
                epsilon, h0, limiting_threshold);
    
        w_right = q_right_rotated[0];
        h_right = w_right - z;
        uh_right = q_right_rotated[1];
        u_right = _compute_speed(&uh_right, &h_right,
                epsilon, h0, limiting_threshold);
    
        // Momentum in y-direction
        vh_left = q_left_rotated[2];
        vh_right = q_right_rotated[2];
    
        // Limit y-momentum if necessary
        // Leaving this out, improves speed significantly (Ole 27/5/2009)
        // All validation tests pass, so do we really need it anymore?
        _compute_speed(&vh_left, &h_left,
                epsilon, h0, limiting_threshold);
        _compute_speed(&vh_right, &h_right,
                epsilon, h0, limiting_threshold);
    
        // Maximal and minimal wave speeds
        soundspeed_left = sqrt(g * h_left);
        soundspeed_right = sqrt(g * h_right);
    
        // Code to use fast square root optimisation if desired.
        // Timings on AMD 64 for the Okushiri profile gave the following timings
        //
        // SQRT           Total    Flux
        //=============================
        //
        // Ref            405s     152s
        // Fast (dbl)     453s     173s
        // Fast (sng)     437s     171s
        //
        // Consequently, there is currently (14/5/2009) no reason to use this
        // approximation.
    
        //soundspeed_left  = fast_squareroot_approximation(g*h_left);
        //soundspeed_right = fast_squareroot_approximation(g*h_right);
    
        s_max = max(u_left + soundspeed_left, u_right + soundspeed_right);
        if (s_max < 0.0) {
            s_max = 0.0;
        }
    
        s_min = min(u_left - soundspeed_left, u_right - soundspeed_right);
        if (s_min > 0.0) {
            s_min = 0.0;
        }
    
        // Flux formulas
        flux_left[0] = u_left*h_left;
        flux_left[1] = u_left * uh_left + 0.5 * g * h_left*h_left;
        flux_left[2] = u_left*vh_left;
    
        flux_right[0] = u_right*h_right;
        flux_right[1] = u_right * uh_right + 0.5 * g * h_right*h_right;
        flux_right[2] = u_right*vh_right;
    
        // Flux computation
        denom = s_max - s_min;
        if (denom < epsilon) { // FIXME (Ole): Try using h0 here
            memset(edgeflux, 0, 3 * sizeof (double));
            *max_speed = 0.0;
        }
        else {
            inverse_denominator = 1.0 / denom;
            for (i = 0; i < 3; i++) {
                edgeflux[i] = s_max * flux_left[i] - s_min * flux_right[i];
                edgeflux[i] += s_max * s_min * (q_right_rotated[i] - q_left_rotated[i]);
                edgeflux[i] *= inverse_denominator;
            }
    
            // Maximal wavespeed
            *max_speed = max(fabs(s_max), fabs(s_min));
    
            // Rotate back
            _rotate(edgeflux, n1, -n2);
        }
    
        return 0;
    }
    
    
    
    /*
     * elements[0] = _timestep
     * elements[1] = _epsilon
     * elements[2] = _H0
     * elements[3] = _g
     * elements[4] = _optimise_dry_cells
     * elements[5] = _compute_fluxes_central_call
    */
    
    static long compute_fluxes_central_call = 1;
    
    __global__ void compute_fluxes_central(
            double * elements,
            double * timestep,
            long* neighbours,
            long* neighbour_edges,
            double* normals,
            double* edgelengths,
            double* radii,
            double* areas,
            long* tri_full_flag,
            double* stage_edge_values,
            double* xmom_edge_values,
            double* ymom_edge_values,
            double* bed_edge_values,
            double* stage_boundary_values,
            double* xmom_boundary_values,
            double* ymom_boundary_values,
            double* stage_explicit_update,
            double* xmom_explicit_update,
            double* ymom_explicit_update,
            long* already_computed_flux,
            double* max_speed_array) 
    {
        int k = threadIdx.x + threadIdx.y + blockIdx.x * blockDim.x + blockIdx.y *blockDim.y;
    
                // Local variables
                double max_speed, length, inv_area, zl, zr;
                double h0 = elements[DH0]*elements[DH0]; // This ensures a good balance when h approaches H0.
    
                double limiting_threshold = 10 * elements[DH0]; // Avoid applying limiter below this
                // threshold for performance reasons.
                // See ANUGA manual under flux limiting
                int i, m, n;
                int ki, nm = 0, ki2; // Index shorthands
    
    
                // Workspace (making them static actually made function slightly slower (Ole))
                double ql[3], qr[3], edgeflux[3]; // Work array for summing up fluxes
    
                //static long call = 1; // Static local variable flagging already computed flux
    
                // Start computation
                //call++; // Flag 'id' of flux calculation for this timestep
    
                // For all triangles
    
                    // Loop through neighbours and compute edge flux for each
                    for (i = 0; i < 3; i++) {
                        ki = k * 3 + i; // Linear index to edge i of triangle k
    
                        if (already_computed_flux[ki] == elements[Dcall]) {
                            // We've already computed the flux across this edge
                            continue;
                        }
    
                        // Get left hand side values from triangle k, edge i
                        ql[0] = stage_edge_values[ki];
                        ql[1] = xmom_edge_values[ki];
                        ql[2] = ymom_edge_values[ki];
                        zl = bed_edge_values[ki];
    
                        // Get right hand side values either from neighbouring triangle
                        // or from boundary array (Quantities at neighbour on nearest face).
                        n = neighbours[ki];
                        if (n < 0) {
                            // Neighbour is a boundary condition
                            m = -n - 1; // Convert negative flag to boundary index
    
                            qr[0] = stage_boundary_values[m];
                            qr[1] = xmom_boundary_values[m];
                            qr[2] = ymom_boundary_values[m];
                            zr = zl; // Extend bed elevation to boundary
                        } else {
                            // Neighbour is a real triangle
                            m = neighbour_edges[ki];
                            nm = n * 3 + m; // Linear index (triangle n, edge m)
    
                            qr[0] = stage_edge_values[nm];
                            qr[1] = xmom_edge_values[nm];
                            qr[2] = ymom_edge_values[nm];
                            zr = bed_edge_values[nm];
                        }
    
                        // Now we have values for this edge - both from left and right side.
    
                        if (elements[Doptimise_dry_cells]) {
                            // Check if flux calculation is necessary across this edge
                            // This check will exclude dry cells.
                            // This will also optimise cases where zl != zr as
                            // long as both are dry
    
                            if (fabs(ql[0] - zl) < elements[Depsilon] &&
                                    fabs(qr[0] - zr) < elements[Depsilon]) {
                                // Cell boundary is dry
    
                                already_computed_flux[ki] = elements[Dcall]; // #k Done
                                if (n >= 0) {
                                    already_computed_flux[nm] = elements[Dcall]; // #n Done
                                }
    
                                max_speed = 0.0;
                                continue;
                            }
                        }
    
    
                        if (fabs(zl - zr) > 1.0e-10) {
                            //report_python_error(AT, "Discontinuous Elevation");
                            //return 0.0;
                            return;
                        }
    
                        // Outward pointing normal vector (domain.normals[k, 2*i:2*i+2])
                        ki2 = 2 * ki; //k*6 + i*2
    
                        // Edge flux computation (triangle k, edge i)
                        _flux_function_central(ql, qr, zl, zr,
                                normals[ki2], normals[ki2 + 1],
                                elements[Depsilon], h0, limiting_threshold, elements[Dg],
                                edgeflux, &max_speed);
    
    
                        // Multiply edgeflux by edgelength
                        length = edgelengths[ki];
                        edgeflux[0] *= length;
                        edgeflux[1] *= length;
                        edgeflux[2] *= length;
    
    
                        // Update triangle k with flux from edge i
                        stage_explicit_update[k] -= edgeflux[0];
                        xmom_explicit_update[k] -= edgeflux[1];
                        ymom_explicit_update[k] -= edgeflux[2];
    
                        already_computed_flux[ki] = elements[Dcall]; // #k Done
    
    
                        // Update neighbour n with same flux but reversed sign
                        if (n >= 0) {
                            stage_explicit_update[n] += edgeflux[0];
                            xmom_explicit_update[n] += edgeflux[1];
                            ymom_explicit_update[n] += edgeflux[2];
    
                            already_computed_flux[nm] = elements[Dcall]; // #n Done
                        }
    
                        // Update timestep based on edge i and possibly neighbour n
                        if (tri_full_flag[k] == 1) {
                            if (max_speed > elements[Depsilon]) {
                                // Apply CFL condition for triangles joining this edge (triangle k and triangle n)
    
                                // CFL for triangle k
                                timestep[k] = min(timestep[k], radii[k] / max_speed);
    
                                if (n >= 0) {
                                    // Apply CFL condition for neigbour n (which is on the ith edge of triangle k)
                                    timestep[k] = min(timestep[k], radii[n] / max_speed);
                                }
    
                                // Ted Rigby's suggested less conservative version
                                //if (n>=0) {
                                //  elements[0] = min(elements[0], (radii[k]+radii[n])/max_speed);
                                //} else {
                                //  elements[0] = min(elements[0], radii[k]/max_speed);
                                // }
                            }
                        }
    
                    } // End edge i (and neighbour n)
    
    
                    // Normalise triangle k by area and store for when all conserved
                    // quantities get updated
                    inv_area = 1.0 / areas[k];
                    stage_explicit_update[k] *= inv_area;
                    xmom_explicit_update[k] *= inv_area;
                    ymom_explicit_update[k] *= inv_area;
    
    
                    // Keep track of maximal speeds
                    max_speed_array[k] = max_speed;
    
    
                //return elements[0];
    }
          """)
    elements = numpy.random.randn(6)
    elements = elements.astype(numpy.float64)
    elements[0] = domain.timestep
    elements[1] = domain.epsilon
    elements[2] = domain.H0
    elements[3] = domain.g
    elements[4] = domain.optimise_dry_cells
    elements[5] = compute_fluxes_central_call

    timestep_array = numpy.zeros( (domain.number_of_elements,1) , dtype=numpy.float64)
    
    for x in timestep_array:
        x = domain.timestep

    for x in domain.quantities['stage'].explicit_update:
        x = 0.0
    for x in domain.quantities['xmomentum'].explicit_update:
        x = 0.0
    for x in domain.quantities['ymomentum'].explicit_update:
        x = 0.0

    compute_fluxes_central_function = mod.get_function('compute_fluxes_central')
    compute_fluxes_central_function( \
            cuda.In( elements ), \
            cuda.InOut( timestep_array ), \
            cuda.In( domain.neighbours ), \
            cuda.In( domain.neighbour_edges ),\
            cuda.In( domain.normals ), \
            cuda.In( domain.edgelengths ), \
            cuda.In( domain.radii ), \
            cuda.In( domain.areas ), \
            cuda.In( domain.tri_full_flag ),\
            cuda.In( domain.quantities['stage'].edge_values ), \
            cuda.In( domain.quantities['xmomentum'].edge_values ), \
            cuda.In( domain.quantities['ymomentum'].edge_values ), \
            cuda.In( domain.quantities['elevation'].edge_values ), \
            cuda.In( domain.quantities['stage'].boundary_values ), \
            cuda.In( domain.quantities['xmomentum'].boundary_values ), \
            cuda.In( domain.quantities['ymomentum'].boundary_values ), \
            cuda.InOut( domain.quantities['stage'].explicit_update ), \
            cuda.InOut( domain.quantities['xmomentum'].explicit_update ), \
            cuda.InOut( domain.quantities['ymomentum'].explicit_update ), \
            cuda.InOut( domain.already_computed_flux ), \
            cuda.InOut( domain.max_speed), \
            block = ( W, W, 1),\
            grid = ( (N + W*W -1 ) / (W*W), 1) )



def compute_fluxes_central_structure_cuda(
	    domain, 
    	parallelFlag = 1, 
	    name = "compute_fluxes_central_structure_CUDA"
        ):

    import pycuda.driver as cuda
    import pycuda.autoinit
    from pycuda.compiler import SourceModule
    import numpy

    N = domain.number_of_elements

    global compute_fluxes_central_structure_call
    compute_fluxes_central_structure_call += 1
    

    mod = SourceModule("""
    
    #define Dtimestep 0
    #define Depsilon 1
    #define DH0 2
    #define Dg 3
    #define Doptimise_dry_cells 4
    #define Dcall 5
    #define Dnumber_of_elements
    


    /**************************/
    /* CUDA double atomic add */
    /**************************/

    __device__ double cuda_atomicAdd(double*address,double val){    
        double assumed,old=*address;    
        
        do{        
            assumed=old;        
            old= __longlong_as_double(atomicCAS((unsigned long long int*)address,
                    __double_as_longlong(assumed),
                    __double_as_longlong(val+assumed)));    
        }while(assumed!=old);    
        return old;
    }



    __device__ int _rotate(double *q, double n1, double n2) {
        /*Rotate the momentum component q (q[1], q[2])
          from x,y coordinates to coordinates based on normal vector (n1, n2).
    
          Result is returned in array 3x1 r
          To rotate in opposite direction, call rotate with (q, n1, -n2)
    
          Contents of q are changed by this function */
   
    
        double q1, q2;
    
        // Shorthands
        q1 = q[1]; // uh momentum
        q2 = q[2]; // vh momentum
    
        // Rotate
        q[1] = n1 * q1 + n2*q2;
        q[2] = -n2 * q1 + n1*q2;
    
        return 0;
    }



    __device__ double _compute_speed(double *uh,
            double *h,
            double epsilon,
            double h0,
            double limiting_threshold) {
    
        double u;
    
        if (*h < limiting_threshold) {
            // Apply limiting of speeds according to the ANUGA manual
            if (*h < epsilon) {
                *h = 0.0; // Could have been negative
                u = 0.0;
            } else {
                u = *uh / (*h + h0 / *h);
            }
    
    
            // Adjust momentum to be consistent with speed
            *uh = u * *h;
        } else {
            // We are in deep water - no need for limiting
            u = *uh / *h;
        }
    
        return u;
    }
    


    


    __device__ int _flux_function_central(double *q_left, double *q_right,
            double z_left, double z_right,
            double n1, double n2,
            double epsilon,
            double h0,
            double limiting_threshold,
            double g,
            double *edgeflux, double *max_speed) 
            {
    
        /*Compute fluxes between volumes for the shallow water wave equation
          cast in terms of the 'stage', w = h+z using
          the 'central scheme' as described in
    
          Kurganov, Noelle, Petrova. 'Semidiscrete Central-Upwind Schemes For
          Hyperbolic Conservation Laws and Hamilton-Jacobi Equations'.
          Siam J. Sci. Comput. Vol. 23, No. 3, pp. 707-740.
    
          The implemented formula is given in equation (3.15) on page 714
         */
    
        int i;
    
        double w_left, h_left, uh_left, vh_left, u_left;
        double w_right, h_right, uh_right, vh_right, u_right;
        double s_min, s_max, soundspeed_left, soundspeed_right;
        double denom, inverse_denominator, z;
        double temp;
    
        // Workspace (allocate once, use many)
        double q_left_rotated[3], q_right_rotated[3], flux_right[3], flux_left[3];
    
        // Copy conserved quantities to protect from modification
        q_left_rotated[0] = q_left[0];
        q_right_rotated[0] = q_right[0];
        q_left_rotated[1] = q_left[1];
        q_right_rotated[1] = q_right[1];
        q_left_rotated[2] = q_left[2];
        q_right_rotated[2] = q_right[2];
    
        // Align x- and y-momentum with x-axis
        //_rotate(q_left_rotated, n1, n2);
        q_left_rotated[1] = n1*q_left[1] + n2*q_left[2];
        q_left_rotated[2] = -n2*q_left[1] + n1*q_left[2];

        //_rotate(q_right_rotated, n1, n2);
        q_right_rotated[1] = n1*q_right[1] + n2*q_right[2];
        q_right_rotated[2] = -n2*q_right[1] + n1*q_right[2];
        
    
        if (fabs(z_left - z_right) > 1.0e-10) {
            //report_python_error(AT, "Discontinuous Elevation");
            return 0.0;
        }
        z = 0.5 * (z_left + z_right); // Average elevation values.
    
        // Compute speeds in x-direction
        w_left = q_left_rotated[0];
        h_left = w_left - z;
        uh_left = q_left_rotated[1];
        //u_left = _compute_speed(&uh_left, &h_left,
        //        epsilon, h0, limiting_threshold);
        
        if (h_left < limiting_threshold) {   
            
            if (h_left < epsilon) {
                h_left = 0.0;  // Could have been negative
                u_left = 0.0;
            } else {
                u_left = uh_left/(h_left + h0/ h_left);    
            }

            uh_left = u_left * h_left;
        } else {
            u_left = uh_left/ h_left;
        }
    
        w_right = q_right_rotated[0];
        h_right = w_right - z;
        uh_right = q_right_rotated[1];
        //u_right = _compute_speed(&uh_right, &h_right,
        //        epsilon, h0, limiting_threshold);

        if (h_right < limiting_threshold) {   
            
            if (h_right < epsilon) {
                h_right = 0.0;  // Could have been negative
                u_right = 0.0;
            } else {
                u_right = uh_right/(h_right + h0/ h_right);    
            }

            uh_right = u_right * h_right;
        } else {
            u_right = uh_right/ h_right;
        }



    
        // Momentum in y-direction
        vh_left = q_left_rotated[2];
        vh_right = q_right_rotated[2];
    
        // Limit y-momentum if necessary
        // Leaving this out, improves speed significantly (Ole 27/5/2009)
        // All validation tests pass, so do we really need it anymore?
        //_compute_speed(&vh_left, &h_left,
        //        epsilon, h0, limiting_threshold);
        //_compute_speed(&vh_right, &h_right,
        //        epsilon, h0, limiting_threshold);
    
        // Maximal and minimal wave speeds
        soundspeed_left = sqrt(g * h_left);
        soundspeed_right = sqrt(g * h_right);
    
        // Code to use fast square root optimisation if desired.
        // Timings on AMD 64 for the Okushiri profile gave the following timings
        //
        // SQRT           Total    Flux
        //=============================
        //
        // Ref            405s     152s
        // Fast (dbl)     453s     173s
        // Fast (sng)     437s     171s
        //
        // Consequently, there is currently (14/5/2009) no reason to use this
        // approximation.
    
        //soundspeed_left  = fast_squareroot_approximation(g*h_left);
        //soundspeed_right = fast_squareroot_approximation(g*h_right);
    
        s_max = max(u_left + soundspeed_left, u_right + soundspeed_right);
        if (s_max < 0.0) {
            s_max = 0.0;
        }
    
        s_min = min(u_left - soundspeed_left, u_right - soundspeed_right);
        if (s_min > 0.0) {
            s_min = 0.0;
        }
    
        // Flux formulas
        flux_left[0] = u_left*h_left;
        flux_left[1] = u_left * uh_left + 0.5 * g * h_left*h_left;
        flux_left[2] = u_left*vh_left;
    
        flux_right[0] = u_right*h_right;
        flux_right[1] = u_right * uh_right + 0.5 * g * h_right*h_right;
        flux_right[2] = u_right*vh_right;
    
        // Flux computation
        denom = s_max - s_min;
        if (denom < epsilon) { // FIXME (Ole): Try using h0 here
            memset(edgeflux, 0, 3 * sizeof (double));
            *max_speed = 0.0;
        }
        else {
            inverse_denominator = 1.0 / denom;
            for (i = 0; i < 3; i++) {
                edgeflux[i] = s_max * flux_left[i] - s_min * flux_right[i];
                edgeflux[i] += s_max * s_min * (q_right_rotated[i] - q_left_rotated[i]);
                edgeflux[i] *= inverse_denominator;
            }
    
            // Maximal wavespeed
            *max_speed = max(fabs(s_max), fabs(s_min));
    
            // Rotate back
            //_rotate(edgeflux, n1, -n2);
            temp = edgeflux[1];
            edgeflux[1] = n1* temp -n2*edgeflux[2];
            edgeflux[2] = n2*temp + n1*edgeflux[2];
        }
    
        return 0;
    }





    /***********************************/
    /* The CUDA compute fluex function */
    /***********************************/


    __global__ void compute_fluxes_central_structure_CUDA(
            long N,
            double * elements,
            double * timestep,
            long * neighbours,
            long * neighbour_edges,
            double * normals,
            double * edgelengths,
            double * radii,
            double * areas,
            long * tri_full_flag,
            double * stage_edge_values,
            double * xmom_edge_values,
            double * ymom_edge_values,
            double * bed_edge_values,
            double * stage_boundary_values,
            double * xmom_boundary_values,
            double * ymom_boundary_values,
            double * stage_explicit_update,
            double * xmom_explicit_update,
            double * ymom_explicit_update, 
            double * max_speed_array)
    {
        const int k = threadIdx.x + (blockIdx.x )*blockDim.x;
    
    
        double max_speed, length, inv_area, zl, zr;
    
        double h0 = elements[DH0] * elements[DH0]; // This ensures a good balance when h approaches H0.
    
        double limiting_threshold = 10 * elements[DH0]; // Avoid applying limiter below this
        
        int i, m, n, j;
        int ki, nm = 0, ki2; // Index shorthands
        
        double ql[3], qr[3], edgeflux[3];

        if (k >=N )
            return;

#ifdef UNSORTED_DOMAIN
        int b[3]={0,1,2};
        spe_bubble_sort( b, neighbours+k*3, k);
        for (j = 0; j < 3; j++) {
            i = b[j];
#else
        for ( i = 0; i < 3; i++) {
#endif
            ki = k * 3 + i; // Linear index to edge i of triangle k
    
            n = neighbours[ki];
    
            ql[0] = stage_edge_values[ki];
            ql[1] = xmom_edge_values[ki];
            ql[2] = ymom_edge_values[ki];
            zl = bed_edge_values[ki];
           
            if (n < 0) {
                m = -n - 1; // Convert negative flag to boundary index
    
                qr[0] = stage_boundary_values[m];
                qr[1] = xmom_boundary_values[m];
                qr[2] = ymom_boundary_values[m];
                zr = zl; // Extend bed elevation to boundary
            } else {
                m = neighbour_edges[ki];
                nm = n * 3 + m; // Linear index (triangle n, edge m)
    
                qr[0] = stage_edge_values[nm];
                qr[1] = xmom_edge_values[nm];
                qr[2] = ymom_edge_values[nm];
                zr = bed_edge_values[nm];
            }
    
            if (elements[Doptimise_dry_cells]) {
                if (fabs(ql[0] - zl) < elements[Depsilon] &&
                        fabs(qr[0] - zr) < elements[Depsilon]) {
    
                    max_speed = 0.0;
                    continue;
                }
            }
            
            ki2 = 2 * ki; //k*6 + i*2
    
    
            _flux_function_central(ql, qr, zl, zr,
                    normals[ki2], normals[ki2 + 1],
                    elements[Depsilon], h0, limiting_threshold, elements[Dg],
                    edgeflux, &max_speed);

            length = edgelengths[ki];
            edgeflux[0] *= length;
            edgeflux[1] *= length;
            edgeflux[2] *= length;
    
    
            
            stage_explicit_update[k] -= edgeflux[0];
            xmom_explicit_update[k] -= edgeflux[1];
            ymom_explicit_update[k] -= edgeflux[2];
        
            //cuda_atomicAdd( stage_explicit_update + k, -edgeflux[0] );
            //cuda_atomicAdd( xmom_explicit_update + k, -edgeflux[1] );
            //cuda_atomicAdd( ymom_explicit_update + k, -edgeflux[2] );
    
            
            if (tri_full_flag[k] == 1) {
                if (max_speed > elements[Depsilon]) {
                    timestep[k] = min(timestep[k], radii[k] / max_speed);
    
                    if (n >= 0) {
                        timestep[k] = min(timestep[k], radii[n] / max_speed);
                    }
                }
            }
    
        } // End edge i (and neighbour n)
    
    
        
        inv_area = 1.0 / areas[k];
        stage_explicit_update[k] *= inv_area;
        xmom_explicit_update[k] *= inv_area;
        ymom_explicit_update[k] *= inv_area;
    
        max_speed_array[k] =  max_speed;
    }
    
    
    
    /*****************************************/
    /*          Sequential C version         */
    /*****************************************/
    __global__ void compute_fluxes_central_structure_serial(
            double * timestep,
            long * neighbours,
            long * neighbour_edges,
            double * normals,
            double * edgelengths,
            double * radii,
            double * areas,
            double * tri_full_flag,
            double * stage_edge_values,
            double * xmom_edge_values,
            double * ymom_edge_values,
            double * bed_edge_values,
            double * stage_boundary_values,
            double * xmom_boundary_values,
            double * ymom_boundary_values,
            double * stage_explicit_update,
            double * xmom_explicit_update,
            double * ymom_explicit_update, 
            long * already_computed_flux,
            double * max_speed_array,
            double * elements)
    {
        int k;
    
        // Local variables
        double max_speed, length, inv_area, zl, zr;
    
        double h0 = elements[DH0] * elements[DH0]; // This ensures a good balance when h approaches H0.
    
        double limiting_threshold = 10 * elements[DH0]; // Avoid applying limiter below this
        // threshold for performance reasons.
        // See ANUGA manual under flux limiting
        int i, m, n;
        int ki, nm = 0, ki2; // Index shorthands
    
    
        // Workspace (making them static actually made function slightly slower (Ole))
        double ql[3], qr[3], edgeflux[3]; // Work array for summing up fluxes
    
        //static long call = 1; // Static local variable flagging already computed flux
    
        // Start computation
        //call++; // Flag 'id' of flux calculation for this timestep
        for(k=0;k<43200;k++)
        {
        // Loop through neighbours and compute edge flux for each
        for (i = 0; i < 3; i++) {
            ki = k * 3 + i; // Linear index to edge i of triangle k
    
            //if (already_computed_flux[ki] == elements[Dcall]) {
                // We've already computed the flux across this edge
            //    continue;
            //}
    
            // Get left hand side values from triangle k, edge i
            ql[0] = stage_edge_values[ki];
            ql[1] = xmom_edge_values[ki];
            ql[2] = ymom_edge_values[ki];
            zl = bed_edge_values[ki];
    
            // Get right hand side values either from neighbouring triangle
            // or from boundary array (Quantities at neighbour on nearest face).
            n = neighbours[ki];
            if (n < 0) {
                // Neighbour is a boundary condition
                m = -n - 1; // Convert negative flag to boundary index
    
                qr[0] = stage_boundary_values[m];
                qr[1] = xmom_boundary_values[m];
                qr[2] = ymom_boundary_values[m];
                zr = zl; // Extend bed elevation to boundary
            } else {
                // Neighbour is a real triangle
                m = neighbour_edges[ki];
                nm = n * 3 + m; // Linear index (triangle n, edge m)
    
                qr[0] = stage_edge_values[nm];
                qr[1] = xmom_edge_values[nm];
                qr[2] = ymom_edge_values[nm];
                zr = bed_edge_values[nm];
            }
    
            // Now we have values for this edge - both from left and right side.
    
            if (elements[Doptimise_dry_cells]) {
                // Check if flux calculation is necessary across this edge
                // This check will exclude dry cells.
                // This will also optimise cases where zl != zr as
                // long as both are dry
    
                if (fabs(ql[0] - zl) < elements[Depsilon] &&
                        fabs(qr[0] - zr) < elements[Depsilon]) {
                    // Cell boundary is dry
    
                    already_computed_flux[ki] = elements[Dcall]; // #k Done
                    if (n >= 0) {
                        already_computed_flux[nm] = elements[Dcall]; // #n Done
                    }
    
                    max_speed = 0.0;
                    continue;
    
                }
            }
    
    
            /*
               if (fabs(zl - zr) > 1.0e-10) {
               report_python_error(AT, "Discontinuous Elevation");
               return 0.0;
               }
             */
    
            // Outward pointing normal vector (domain.normals[k, 2*i:2*i+2])
            ki2 = 2 * ki; //k*6 + i*2
    
            // Edge flux computation (triangle k, edge i)
            _flux_function_central(ql, qr, zl, zr,
                    normals[ki2], normals[ki2 + 1],
                    elements[Depsilon], h0, limiting_threshold, elements[Dg],
                    edgeflux, &max_speed);
    
    
            // Multiply edgeflux by edgelength
            length = edgelengths[ki];
            edgeflux[0] *= length;
            edgeflux[1] *= length;
            edgeflux[2] *= length;
    
    
            // Update triangle k with flux from edge i
            stage_explicit_update[k] -= edgeflux[0];
            xmom_explicit_update[k] -= edgeflux[1];
            ymom_explicit_update[k] -= edgeflux[2];
    
            already_computed_flux[ki] = elements[Dcall]; // #k Done
    
    
            // Update neighbour n with same flux but reversed sign
            //if (n >= 0) {
            //    stage_explicit_update[n] += edgeflux[0];
            //    xmom_explicit_update[n] += edgeflux[1];
            //    ymom_explicit_update[n] += edgeflux[2];
    
            //    already_computed_flux[nm] = elements[Dcall]; // #n Done
            //}
    
            // Update timestep based on edge i and possibly neighbour n
            if (tri_full_flag[k] == 1) {
                if (max_speed > elements[Depsilon]) {
                    // Apply CFL condition for triangles joining this edge (triangle k and triangle n)
    
                    // CFL for triangle k
                    elements[Dtimestep] = min(elements[Dtimestep], radii[k] / max_speed);
    
                    if (n >= 0) {
                        // Apply CFL condition for neigbour n (which is on the ith edge of triangle k)
                        elements[Dtimestep] = min(elements[Dtimestep], radii[n] / max_speed);
                    }
    
                    // Ted Rigbys suggested less conservative version
                    //if (n>=0) {
                    //  timestep = min(timestep, (radii[k]+radii[n])/max_speed);
                    //} else {
                    //  timestep = min(timestep, radii[k]/max_speed);
                    // }
                }
            }
        } // End edge i (and neighbour n)
        // Normalise triangle k by area and store for when all conserved
        // quantities get updated
        inv_area = 1.0 / areas[k];
        stage_explicit_update[k] *= inv_area;
        xmom_explicit_update[k] *= inv_area;
        ymom_explicit_update[k] *= inv_area;
    
    
        // Keep track of maximal speeds
        max_speed_array[k] = max_speed;
        }
        //return timestep;
    }
   



    /*****************************************/
    /* Rearrange the domain variable order   */
    /* so as to achieve memory corelasing    */
    /*****************************************/

    __global__ void compute_fluxes_central_structure_MeCo(
            double * timestep,
            long * neighbours,
            long * neighbour_edges,
            double * normals,
            double * edgelengths,
            double * radii,
            double * areas,
            long * tri_full_flag,
            double * stage_edge_values,
            double * xmom_edge_values,
            double * ymom_edge_values,
            double * bed_edge_values,
            double * stage_boundary_values,
            double * xmom_boundary_values,
            double * ymom_boundary_values,
            double * stage_explicit_update,
            double * xmom_explicit_update,
            double * ymom_explicit_update, 
            double * max_speed_array,
            double * elements
			)
    {
        const int k = threadIdx.x + (blockIdx.x )*blockDim.x;
    
    
        double max_speed, length, inv_area, zl, zr;
    
        double h0 = elements[DH0] * elements[DH0]; // This ensures a good balance when h approaches H0.
    
        double limiting_threshold = 10 * elements[DH0]; // Avoid applying limiter below this
        
        int i, m, n;
        int ki, nm = 0, ki2; // Index shorthands
        
        double ql[3], qr[3], edgeflux[3];

        //int b[3]={0,1,2};
        //spe_bubble_sort( b, neighbours+k*3, k);
        //for (j = 0; j < 3; j++) {
        //    i = b[j];

		 __shared__ double sh_data[32*9];
		
        for ( i = 0; i < 3; i++) {

            ki = k * 3 + i; // Linear index to edge i of triangle k
    
            n = neighbours[ki];
    
            ql[0] = stage_edge_values[ki];
            ql[1] = xmom_edge_values[ki];
            ql[2] = ymom_edge_values[ki];
            zl = bed_edge_values[ki];
           
            if (n < 0) {
                m = -n - 1; // Convert negative flag to boundary index
    
                qr[0] = stage_boundary_values[m];
                qr[1] = xmom_boundary_values[m];
                qr[2] = ymom_boundary_values[m];
                zr = zl; // Extend bed elevation to boundary
            } else {
                m = neighbour_edges[ki];
                nm = n * 3 + m; // Linear index (triangle n, edge m)
    
                qr[0] = stage_edge_values[nm];
                qr[1] = xmom_edge_values[nm];
                qr[2] = ymom_edge_values[nm];
                zr = bed_edge_values[nm];
            }
    
            if (elements[Doptimise_dry_cells]) {
                if (fabs(ql[0] - zl) < elements[Depsilon] &&
                        fabs(qr[0] - zr) < elements[Depsilon]) {
    
                    max_speed = 0.0;
                    continue;
                }
            }
            
            ki2 = 2 * ki; //k*6 + i*2
    
    
            _flux_function_central(ql, qr, zl, zr,
                    normals[ki2], normals[ki2 + 1],
                    elements[Depsilon], h0, limiting_threshold, elements[Dg],
                    edgeflux, &max_speed);

            length = edgelengths[ki];
            edgeflux[0] *= length;
            edgeflux[1] *= length;
            edgeflux[2] *= length;
    
    
            
            //stage_explicit_update[k] -= edgeflux[0];
            //xmom_explicit_update[k] -= edgeflux[1];
            //ymom_explicit_update[k] -= edgeflux[2];
			//temp_array[k*3+i] = -edgeflux[1];

			sh_data[threadIdx.x + i*blockDim.x]= -edgeflux[0];
			sh_data[threadIdx.x + (i+3)*blockDim.x]= -edgeflux[1];
			sh_data[threadIdx.x + (i+6)*blockDim.x]= -edgeflux[2];
			__syncthreads();
        
            //cuda_atomicAdd( stage_explicit_update + k, -edgeflux[0] );
            //cuda_atomicAdd( xmom_explicit_update + k, -edgeflux[1] );
            //cuda_atomicAdd( ymom_explicit_update + k, -edgeflux[2] );
    
            
            if (tri_full_flag[k] == 1) {
                if (max_speed > elements[Depsilon]) {
                    timestep[k] = min(timestep[k], radii[k] / max_speed);
    
                    if (n >= 0) {
                        timestep[k] = min(timestep[k], radii[n] / max_speed);
                    }
                }
            }
    
        } // End edge i (and neighbour n)
    
    
        
        inv_area = 1.0 / areas[k];
        //stage_explicit_update[k] *= inv_area;
        //xmom_explicit_update[k] *= inv_area;
        //ymom_explicit_update[k] *= inv_area;
		stage_explicit_update[k] = (
				sh_data[threadIdx.x] + 
				sh_data[threadIdx.x+blockDim.x] +
				sh_data[threadIdx.x+blockDim.x*2]) * inv_area;

		xmom_explicit_update[k] = (
				sh_data[threadIdx.x+3*blockDim.x] + 
				sh_data[threadIdx.x+4*blockDim.x] +
				sh_data[threadIdx.x+5*blockDim.x]) * inv_area;

		ymom_explicit_update[k] = (
				sh_data[threadIdx.x+6*blockDim.x] + 
				sh_data[threadIdx.x+7*blockDim.x] +
				sh_data[threadIdx.x+8*blockDim.x]) * inv_area;


		__syncthreads();
    
        max_speed_array[k] =  max_speed;
    }
 
    
    
    """)
    
    elements = numpy.random.randn(6)
    elements = elements.astype(numpy.float64)
    elements[0] = domain.evolve_max_timestep
    elements[1] = domain.epsilon
    elements[2] = domain.H0
    elements[3] = domain.g
    elements[4] = domain.optimise_dry_cells
    elements[5] = compute_fluxes_central_structure_call
    

    timestep_array = numpy.zeros( 
            (domain.number_of_elements) , 
            dtype=numpy.float64)

    # Necessary data initialization
    # FIXME: can be done by CUDA threads
    for i in range(N):
        timestep_array[i] = domain.evolve_max_timestep
        #domain.quantities['stage'].explicit_update[i] = 0.0
        #domain.quantities['xmomentum'].explicit_update[i]=0.0
        #domain.quantities['ymomentum'].explicit_update[i] =0.0
    

    # FIXME: add hardware dependency 
    #if (N % 320 == 0):
    #    W1 = 32
    #    W2 = 10
    #elif(N % 256 == 0):
    #    W1 = 16
    #    W2 = 16
    #elif (N % 32 ==0):
    #    W1 = 32
    #    W2 = 1
    #else:
    #    raise Exception('N can not be splited')

    W1 = 32
    W2 = 1
    W3 = 1
    
    if parallelFlag == 2:
        compute_fluxes_central_function = mod.get_function(name)
        compute_fluxes_central_function( \
                cuda.InOut( elements ), \
                cuda.InOut( timestep_array ), \
                cuda.In( domain.neighbours ), \
                cuda.In( domain.neighbour_edges ),\
                cuda.In( domain.normals ), \
                cuda.In( domain.edgelengths ), \
                cuda.In( domain.radii ), \
                cuda.In( domain.areas ), \
                cuda.In( domain.tri_full_flag ),\
                cuda.In( domain.quantities['stage'].edge_values ), \
                cuda.In( domain.quantities['xmomentum'].edge_values ), \
                cuda.In( domain.quantities['ymomentum'].edge_values ), \
                cuda.In( domain.quantities['elevation'].edge_values ), \
                cuda.In( domain.quantities['stage'].boundary_values ), \
                cuda.In( domain.quantities['xmomentum'].boundary_values ), \
                cuda.In( domain.quantities['ymomentum'].boundary_values ), \
                cuda.InOut( domain.quantities['stage'].explicit_update ), \
                cuda.InOut( domain.quantities['xmomentum'].explicit_update),\
                cuda.InOut( domain.quantities['ymomentum'].explicit_update),\
                cuda.InOut( domain.already_computed_flux ), \
                cuda.InOut( domain.max_speed), \
                block = ( 1, 1, 1),\
                grid = ( 1 , 1) )

        print elements[0]
    else:
        # Below commented code are manual memory-copy procedures 
        """
            compute_fluxes_central_function = mod.get_function(name)

		    from anuga_cuda.merimbula_data.utility import mem_all_cpy

            # facilities vertors
            timestep_array_gpu = mem_all_cpy(timestep_array )
            
            elements_gpu = mem_all_cpy(elements)

            # neighbours vertors
            neighbours_gpu = mem_all_cpy(domain.neighbours)
            neighbour_edges_gpu = mem_all_cpy(domain.neighbour_edges)

            # 
            normals_gpu = mem_all_cpy(domain.normals)
            edgelengths_gpu = mem_all_cpy(domain.edgelengths)
            radii_gpu = mem_all_cpy(domain.radii)
            areas_gpu = mem_all_cpy(domain.areas)
            tri_full_flag_gpu = mem_all_cpy(domain.tri_full_flag)
            vertex_coordinates_gpu = mem_all_cpy(domain.vertex_coordinates)

            # edge_values vertors
            stage_edge_values_gpu = mem_all_cpy(
		    		domain.quantities['stage'].edge_values)
            xmom_edge_values_gpu = mem_all_cpy(
		    		domain.quantities['xmomentum'].edge_values)
            ymom_edge_values_gpu = mem_all_cpy(
		    		domain.quantities['ymomentum'].edge_values)
            bed_edge_values_gpu = mem_all_cpy(
		    		domain.quantities['elevation'].edge_values)

            # boundary_values vertors
            stage_boundary_values_gpu = mem_all_cpy(
		    		domain.quantities['stage'].boundary_values)
            xmom_boundary_values_gpu = mem_all_cpy(
		    		domain.quantities['xmomentum'].boundary_values)
            ymom_boundary_values_gpu = mem_all_cpy(
		    		domain.quantities['ymomentum'].boundary_values)

            # vertex_values
            stage_vertex_values_gpu = mem_all_cpy(
		    		domain.quantities['stage'].vertex_values)

            # centroid_values
            stage_centroid_values_gpu = mem_all_cpy(
		    		domain.quantities['stage'].centroid_values)
            bed_centroid_values_gpu = mem_all_cpy(
		    		domain.quantities['elevation'].centroid_values)

            # explicit_update vertors
            stage_explicit_update_gpu = mem_all_cpy(
		    		domain.quantities['stage'].explicit_update)
            xmom_explicit_update_gpu = mem_all_cpy(
		    		domain.quantities['xmomentum'].explicit_update)
            ymom_explicit_update_gpu = mem_all_cpy(
		    		domain.quantities['ymomentum'].explicit_update)

            # max_speed vertors
            max_speed_gpu = mem_all_cpy(domain.max_speed)
            

            compute_fluxes_central_function(
                elements_gpu, 
                timestep_array_gpu, 
                neighbours_gpu, 
                neighbour_edges_gpu,
                normals_gpu, 
                edgelengths_gpu, 
                radii_gpu, 
                areas_gpu, 
                tri_full_flag_gpu,
                stage_edge_values_gpu, 
                xmom_edge_values_gpu, 
                ymom_edge_values_gpu, 
                bed_edge_values_gpu, 
                stage_boundary_values_gpu, 
                xmom_boundary_values_gpu, 
                ymom_boundary_values_gpu, 
                stage_explicit_update_gpu, 
                xmom_explicit_update_gpu, 
                ymom_explicit_update_gpu, 
                max_speed_gpu, 
                block = ( W1, 1, 1),
                grid = ( (N + W1 - 1)/W1, 1) 
                )
            
            cuda.memcpy_dtoh(
                    domain.quantities['stage'].explicit_update, 
                    stage_explicit_update_gpu)
            cuda.memcpy_dtoh(
                    domain.quantities['xmomentum'].explicit_update, 
                    xmom_explicit_update_gpu)
            cuda.memcpy_dtoh(
                    domain.quantities['ymomentum'].explicit_update, 
                    ymom_explicit_update_gpu)
            
            cuda.memcpy_dtoh(timestep_array, timestep_array_gpu)
            cuda.memcpy_dtoh(domain.max_speed, max_speed_gpu)
        """
        
        

        W1 =32
        W2 = 1
        W3 =1

        # Source the cuda function
        #compute_fluxes_central_function = mod.get_function(name)
        #from anuga_cuda.compute_fluxes.py_compute_fluxes import get_func_info
        #get_func_info(compute_fluxes_central_function)
        from anuga_cuda import kernel_path as kp
        compute_fluxes_mod = SourceModule(
            open( kp["compute_fluxes_dir"]+"compute_fluxes.cu").read(),
            )
        compute_fluxes_central_function = compute_fluxes_mod.get_function(
            #"compute_fluxes_central_structure_cuda_single")
            "compute_fluxes_central_structure_CUDA")
        if 1:
            compute_fluxes_central_function( 
                    numpy.uint(domain.number_of_elements),
                    numpy.float64(domain.g),
                    numpy.float64(domain.epsilon),
                    numpy.float64(domain.H0 * domain.H0),
                    numpy.float64(domain.H0*10),
                    numpy.uint(domain.optimise_dry_cells),

                    cuda.InOut( timestep_array ), 
                    cuda.In( domain.neighbours ), 
                    cuda.In( domain.neighbour_edges ),
                    cuda.In( domain.normals ), 
                    cuda.In( domain.edgelengths ), 
                    cuda.In( domain.radii ), 
                    cuda.In( domain.areas ), 
                    cuda.In( domain.tri_full_flag ),
                    cuda.In( domain.quantities['stage'].edge_values ), 
                    cuda.In( domain.quantities['xmomentum'].edge_values ), 
                    cuda.In( domain.quantities['ymomentum'].edge_values ), 
                    cuda.In( domain.quantities['elevation'].edge_values ), 
                    cuda.In( domain.quantities['stage'].boundary_values ), 
                    cuda.In( domain.quantities['xmomentum'].boundary_values ), 
                    cuda.In( domain.quantities['ymomentum'].boundary_values ), 
                    cuda.InOut( domain.quantities['stage'].explicit_update ), 
                    cuda.InOut( domain.quantities['xmomentum'].explicit_update ), 
                    cuda.InOut( domain.quantities['ymomentum'].explicit_update ),  
                    cuda.InOut( domain.max_speed), 
                    block = ( W1, W2, W3),
                    grid = ((N+W1*W2*W3-1)/(W1*W2*W3), 1) )
        else:
            compute_fluxes_central_function( 
                    numpy.uint(domain.number_of_elements),
                    cuda.InOut( elements ), 
                    cuda.InOut( timestep_array ), 
                    cuda.In( domain.neighbours ), 
                    cuda.In( domain.neighbour_edges ),
                    cuda.In( domain.normals ), 
                    cuda.In( domain.edgelengths ), 
                    cuda.In( domain.radii ), 
                    cuda.In( domain.areas ), 
                    cuda.In( domain.tri_full_flag ),
                    cuda.In( domain.quantities['stage'].edge_values ), 
                    cuda.In( domain.quantities['xmomentum'].edge_values ), 
                    cuda.In( domain.quantities['ymomentum'].edge_values ), 
                    cuda.In( domain.quantities['elevation'].edge_values ), 
                    cuda.In( domain.quantities['stage'].boundary_values ), 
                    cuda.In( domain.quantities['xmomentum'].boundary_values ), 
                    cuda.In( domain.quantities['ymomentum'].boundary_values ), 
                    cuda.InOut( domain.quantities['stage'].explicit_update ), 
                    cuda.InOut( domain.quantities['xmomentum'].explicit_update ), 
                    cuda.InOut( domain.quantities['ymomentum'].explicit_update ),  
                    cuda.InOut( domain.max_speed), 
                    block = ( W1, W2, W3),
                    grid = ((N+W1*W2*W3-1)/(W1*W2*W3), 1) )
		
        
        # Get the smallest timestep
        b = numpy.argsort(timestep_array)
        domain.flux_timestep = timestep_array[b[0]] 
    return domain.flux_timestep

        



""" Below 5 functions are the python version of compute fluxes """

def _rotate(q,  n1,  n2) :
    q1 = q[1]
    q2 = q[2]

    q[1] = n1 * q1 + n2*q2
    q[2] = -n2 * q1 + n1*q2

    return 0



def _compute_speed(uh, h, epsilon, h0, limiting_threshold):
    if h < limiting_threshold:
        if h < epsilon :
            h = 0.0
            u = 0.0
        else :
            u = uh / (h + h0 / h)

        uh = u * h
    else :
        u = uh / h

    return u


    
def _flux_function_central(
        q_left, q_right, 
        z_left, z_right, 
        n1, n2, 
        epsilon, h0, 
        limiting_threshold, g, edgeflux, 
        max_speed):

    import numpy
    import math
    q_left_rotated = numpy.zeros(3, dtype=numpy.float64)
    q_right_rotated = numpy.zeros(3, dtype=numpy.float64)
    flux_right = numpy.zeros(3, dtype=numpy.float64)
    flux_left = numpy.zeros(3, dtype=numpy.float64)

    q_left_rotated[0] = q_left[0]
    q_right_rotated[0] = q_right[0]
    q_left_rotated[1] = q_left[1]
    q_right_rotated[1] = q_right[1]
    q_left_rotated[2] = q_left[2]
    q_right_rotated[2] = q_right[2]

    #_rotate(q_left_rotated, n1, n2)
    q1 = q_left_rotated[1]
    q2 = q_left_rotated[2]
    q_left_rotated[1] = n1 * q1 + n2*q2
    q_left_rotated[2] = -n2 * q1 + n1*q2
    #_rotate(q_right_rotated, n1, n2)
    q1 = q_right_rotated[1]
    q2 = q_right_rotated[2]
    q_right_rotated[1] = n1 * q1 + n2*q2
    q_right_rotated[2] = -n2 * q1 + n1*q2


    z = 0.5 * (z_left + z_right)

    w_left = q_left_rotated[0]
    h_left = w_left - z
    uh_left = q_left_rotated[1]
    u_left = _compute_speed(uh_left, h_left,
            epsilon, h0, limiting_threshold)

    w_right = q_right_rotated[0]
    h_right = w_right - z
    uh_right = q_right_rotated[1]
    u_right = _compute_speed(uh_right, h_right,
            epsilon, h0, limiting_threshold)
    #if h_right < limiting_threshold:
    #    if h_right < epsilon :
    #        h_right = 0.0
    #        u_right = 0.0
    #    else :
    #        u_right = uh_right / (h_right + h0 / h_right)

    #    uh_right = u_right * h_right
    #else :
    #    u_right = uh_right / h_right

    
    vh_left = q_left_rotated[2]
    vh_right = q_right_rotated[2]

    _compute_speed(vh_left, h_left,
            epsilon, h0, limiting_threshold)
    _compute_speed(vh_right, h_right,
            epsilon, h0, limiting_threshold)

    
    soundspeed_left = math.sqrt(g * h_left)
    soundspeed_right = math.sqrt(g * h_right)

    s_max = max(u_left + soundspeed_left, u_right + soundspeed_right)
    if s_max < 0.0 :
        s_max = 0.0

    s_min = min(u_left - soundspeed_left, u_right - soundspeed_right)
    if s_min > 0.0 :
        s_min = 0.0

    
    flux_left[0] = u_left*h_left
    flux_left[1] = u_left * uh_left + 0.5 * g * h_left*h_left
    flux_left[2] = u_left*vh_left

    flux_right[0] = u_right*h_right
    flux_right[1] = u_right * uh_right + 0.5 * g * h_right*h_right
    flux_right[2] = u_right*vh_right

    
    denom = s_max - s_min
    if denom < epsilon :
        edgeflux[0] = 0.0
        edgeflux[1] = 0.0
        edgeflux[2] = 0.0
        max_speed = 0.0
    else:
        inverse_denominator = 1.0 / denom
        for i in range(3):
            edgeflux[i] = s_max * flux_left[i] - s_min * flux_right[i]
            edgeflux[i] += s_max * s_min * (q_right_rotated[i] - q_left_rotated[i])
            edgeflux[i] *= inverse_denominator
        
        max_speed = max(abs(s_max), abs(s_min))

        #_rotate(edgeflux, n1, -n2)
        q1 = edgeflux[1]
        q2 = edgeflux[2]

        edgeflux[1] = n1 * q1 + -n2*q2
        edgeflux[2] = n2 * q1 + n1*q2
    return max_speed
    

# This function helps to reorder the compute sequence so as to make it 
# as same as original C function
def spe_bubble_sort( _list , neighbours, k):
    """
        In the case of the original C function, all the triangles is calculated
        in their index order. 

        Here we put the valid neighbour (not the boundary) with smaller index 
        number than the current one's ahead
    """
    if neighbours[_list[2]]>=0 and neighbours[ _list[2] ] < k and \
            (neighbours[ _list[1] ]<0 or neighbours[ _list[2] ] <neighbours [_list[1] ]):
        _list[2], _list[1] = _list[1], _list[2]

    if neighbours[_list[1]]>=0 and neighbours[ _list[1] ] < k and \
            (neighbours[_list[0]]<0 or neighbours[_list[1]]<neighbours[_list[0]]):
        _list[1], _list[0] = _list[0], _list[1]

    if neighbours[_list[2]]>=0 and neighbours[ _list[2] ] < k and \
            (neighbours[_list[1]]<0 or neighbours[_list[2]]<neighbours[_list[1]]):
        _list[2], _list[1] = _list[1], _list[2]



def compute_fluxes_central_structure_single(domain, k =-1, i =0, call = 2):
    import numpy
    
    ql = numpy.zeros(3, dtype=numpy.float64)
    qr = numpy.zeros(3, dtype=numpy.float64)
    edgeflux = numpy.zeros(3, dtype=numpy.float64)
    counter = 0

    if k < 0:
        for k in range(domain.number_of_elements):
            max_speed_total = 0
        ### Used to generated single triangle values
        #if k:
            
            domain.quantities['stage'].explicit_update[k] = 0.0
            domain.quantities['xmomentum'].explicit_update[k] = 0.0
            domain.quantities['ymomentum'].explicit_update[k] = 0.0
            
            ### To change the calculation order as same as the C
            #b = [0,1,2]
            #spe_bubble_sort(b, domain.neighbours[k], k)
            #for i in b:
            max_speed = 0
            for i in range(3):

            ### The original C Function
            #for i in range(3):

            ### Used to generate single triangle values
            #if i>=0:

                ### The original C sequential Function
                #if domain.already_computed_flux[k][i] == call:
                #    continue

                n = domain.neighbours[k][i]


                ql[0] = domain.quantities['stage'].edge_values[k][i]
                ql[1] = domain.quantities['xmomentum'].edge_values[k][i]
                ql[2] = domain.quantities['ymomentum'].edge_values[k][i]
                zl = domain.quantities['elevation'].edge_values[k][i]

                if n < 0:
                    m = -n -1
                    qr[0] = domain.quantities['stage'].boundary_values[m]
                    qr[1] = domain.quantities['xmomentum'].boundary_values[m]
                    qr[2] = domain.quantities['ymomentum'].boundary_values[m]
                    zr = zl
                else:
                    m = domain.neighbour_edges[k][i]

                    qr[0] = domain.quantities['stage'].edge_values[n][m]
                    qr[1] = domain.quantities['xmomentum'].edge_values[n][m]
                    qr[2] = domain.quantities['ymomentum'].edge_values[n][m]
                    zr = domain.quantities['elevation'].edge_values[n][m]

                if domain.optimise_dry_cells:
                    if abs(ql[0] -zl) < domain.epsilon and \
                        abs(qr[0] - zr) < domain.epsilon:
                        max_speed = 0
                        continue


                max_speed = _flux_function_central(ql, qr, zl, zr, 
                            domain.normals[k][2*i], domain.normals[k][2*i + 1],
                            domain.epsilon, pow(domain.H0,2), 
                            10*domain.H0, domain.g, edgeflux, max_speed)

                edgeflux[0] *= domain.edgelengths[k][i]
                edgeflux[1] *= domain.edgelengths[k][i]
                edgeflux[2] *= domain.edgelengths[k][i]
                
                ### Original C Function
                domain.quantities['stage'].explicit_update[k] -= edgeflux[0]
                domain.quantities['xmomentum'].explicit_update[k] -= edgeflux[1]
                domain.quantities['ymomentum'].explicit_update[k] -= edgeflux[2]

                ### The original C sequential Function
                #if n >= 0:
                #    domain.quantities['stage'].explicit_update[n] += edgeflux[0]
                #    domain.quantities['xmomentum'].explicit_update[n] += edgeflux[1]
                #    domain.quantities['ymomentum'].explicit_update[n] += edgeflux[2]
                #    domain.already_computed_flux[n][m] = call
                if n>=0 and n < k :
                    max_speed = 0
                else:
                    max_speed_total = max(max_speed_total, max_speed) 
            
            inv_area = 1.0 / domain.areas[k]
            domain.quantities['stage'].explicit_update[k] *= inv_area
            domain.quantities['xmomentum'].explicit_update[k] *= inv_area
            domain.quantities['ymomentum'].explicit_update[k] *= inv_area
            
            domain.max_speed[k] = max_speed_total
    else:
        max_speed_total = 0
        
        domain.quantities['stage'].explicit_update[k] = 0.0
        domain.quantities['xmomentum'].explicit_update[k] = 0.0
        domain.quantities['ymomentum'].explicit_update[k] = 0.0
        
        ### To change the calculation order as same as the C
        #b = [0,1,2]
        #spe_bubble_sort(b, domain.neighbours[k], k)
        #for i in b:
        max_speed = 0
        for i in range(3):

        ### The original C Function
        #for i in range(3):

        ### Used to generate single triangle values
        #if i>=0:

            ### The original C sequential Function
            #if domain.already_computed_flux[k][i] == call:
            #    continue

            n = domain.neighbours[k][i]


            ql[0] = domain.quantities['stage'].edge_values[k][i]
            ql[1] = domain.quantities['xmomentum'].edge_values[k][i]
            ql[2] = domain.quantities['ymomentum'].edge_values[k][i]
            zl = domain.quantities['elevation'].edge_values[k][i]

            if n < 0:
                m = -n -1
                qr[0] = domain.quantities['stage'].boundary_values[m]
                qr[1] = domain.quantities['xmomentum'].boundary_values[m]
                qr[2] = domain.quantities['ymomentum'].boundary_values[m]
                zr = zl
            else:
                m = domain.neighbour_edges[k][i]

                qr[0] = domain.quantities['stage'].edge_values[n][m]
                qr[1] = domain.quantities['xmomentum'].edge_values[n][m]
                qr[2] = domain.quantities['ymomentum'].edge_values[n][m]
                zr = domain.quantities['elevation'].edge_values[n][m]

            if domain.optimise_dry_cells:
                if abs(ql[0] -zl) < domain.epsilon and \
                    abs(qr[0] - zr) < domain.epsilon:
                    max_speed = 0
                    continue


            max_speed = _flux_function_central(ql, qr, zl, zr, 
                        domain.normals[k][2*i], domain.normals[k][2*i + 1],
                        domain.epsilon, pow(domain.H0,2), 
                        10*domain.H0, domain.g, edgeflux, max_speed)

            edgeflux[0] *= domain.edgelengths[k][i]
            edgeflux[1] *= domain.edgelengths[k][i]
            edgeflux[2] *= domain.edgelengths[k][i]
            
            ### Original C Function
            domain.quantities['stage'].explicit_update[k] -= edgeflux[0]
            domain.quantities['xmomentum'].explicit_update[k] -= edgeflux[1]
            domain.quantities['ymomentum'].explicit_update[k] -= edgeflux[2]

            ### The original C sequential Function
            #if n >= 0:
            #    domain.quantities['stage'].explicit_update[n] += edgeflux[0]
            #    domain.quantities['xmomentum'].explicit_update[n] += edgeflux[1]
            #    domain.quantities['ymomentum'].explicit_update[n] += edgeflux[2]
            #    domain.already_computed_flux[n][m] = call
            print k, max_speed
            max_speed_total = max(max_speed_total, max_speed) 
        
        inv_area = 1.0 / domain.areas[k]
        domain.quantities['stage'].explicit_update[k] *= inv_area
        domain.quantities['xmomentum'].explicit_update[k] *= inv_area
        domain.quantities['ymomentum'].explicit_update[k] *= inv_area
        
        domain.max_speed[k] = max_speed_total
        return max_speed_total




if __name__ == '__main__':
    import pycuda.driver as drv
    import pycuda.autoinit
    ctx = pycuda.autoinit.context

    from pycuda.compiler import SourceModule

    import numpy 
    from time import time
    
    from anuga_cuda import *


    testing_gpu_domain = True
    testing_python_version = False
    using_rearranged_domain = False

    # This will reorder edges in order to let the one bordering on
    # triangle with smaller index number compute first
    domain2 = generate_merimbula_domain(gpu = True)
    #domain2 = generate_channel3_domain()
    #sort_domain(domain2)


    domain1 = generate_merimbula_domain()
    #domain1 = generate_channel3_domain()
    #for t in domain1.evolve(yieldstep = 50, finaltime = 50):
    #    pass
    
    host_macro = ""
    if using_rearranged_domain: 
        domain2 = rearrange_domain(domain2)
        sort_domain(domain1)
        host_macro="#define REARRANGED_DOMAIN\n"


    print " Number of elements is: %d" % domain1.number_of_elements
    """
    CUDA Function
    """
    print "CUDA Function"
    if  testing_gpu_domain:
        from anuga_cuda import kernel_path as kp
        compute_fluxes_mod = SourceModule(
                host_macro+ open( kp["compute_fluxes_dir"]+"compute_fluxes.cu").read(),
                arch = 'compute_20',
                code = 'sm_20'
                #options = ['-use_fast_math', '--prec-div=false', '--compiler-options', '-O2']
                )
        compute_fluxes_central_function = compute_fluxes_mod.get_function(
                #"compute_fluxes_central_structure_cuda_single")
                "compute_fluxes_central_structure_CUDA")

        domain2.allocate_device_array()

        # Testing the mem coalscing
        #N = domain2.number_of_elements
        #for i in range(N):
        #    domain2.neighbours[i][:] = [i+1,i+2,i+3]
        #domain2.neighbours[N-3][:] = [N-2, N-1, N-4]
        #domain2.neighbours[N-2][:] = [N-1, N-3, N-4]
        #domain2.neighbours[N-1][:] = [N-2, N-3, N-4]

        domain2.asynchronous_transfer()
        ctx.synchronize()

        compute_fluxes_central_function.set_cache_config(drv.func_cache.PREFER_L1)
        #compute_fluxes_central_function.set_cache_config(drv.func_cache.PREFER_NONE)
        #compute_fluxes_central_function.set_cache_config(drv.func_cache.PREFER_EQUAL)
        #compute_fluxes_central_function.set_cache_config(drv.func_cache.PREFER_SHARED)
        import sys
        W1 = 0
        W2 = 0
        for i in range( len(sys.argv)):
            if sys.argv[i] == "-b":
                W1 = int(sys.argv[i+1])
            elif sys.argv[i] == "-b1":
                W2 = int(sys.argv[i+1])

        if not W1:
            W1 = compute_fluxes_central_function.max_threads_per_block
        if not W2:
            W2 = 1
        print W1
        W3 = 1

        get_kernel_function_info(compute_fluxes_central_function,W1,W2, W3)

        #domain2.compute_fluxes_func(
        compute_fluxes_central_function(
                numpy.uint32(domain2.number_of_elements),
                numpy.float64(domain2.evolve_max_timestep),
                numpy.float64(domain2.g),
                numpy.float64(domain2.epsilon),
                numpy.float64(domain2.H0 * domain2.H0),
                numpy.float64(domain2.H0 * 10),
                numpy.uint32(domain2.optimise_dry_cells),

                domain2.timestep_array_gpu,
                domain2.neighbours_gpu,
                domain2.neighbour_edges_gpu,
                domain2.normals_gpu,
                domain2.edgelengths_gpu,
                domain2.radii_gpu,
                domain2.areas_gpu,
                domain2.tri_full_flag_gpu,
                domain2.quantities['stage'].edge_values_gpu,
                domain2.quantities['xmomentum'].edge_values_gpu,
                domain2.quantities['ymomentum'].edge_values_gpu,
                domain2.quantities['elevation'].edge_values_gpu,
                domain2.quantities['stage'].boundary_values_gpu,
                domain2.quantities['xmomentum'].boundary_values_gpu,
                domain2.quantities['ymomentum'].boundary_values_gpu,
                domain2.quantities['stage'].explicit_update_gpu,
                domain2.quantities['xmomentum'].explicit_update_gpu,
                domain2.quantities['ymomentum'].explicit_update_gpu,
                domain2.max_speed_gpu,
                block = (W1, W2, W3),
                grid = ((domain2.number_of_elements+W1*W2*W3-1)/(W1*W2*W3), 1)
                )

        #ctx.synchronize()
        drv.memcpy_dtoh(domain2.quantities['stage'].explicit_update,
                        domain2.quantities['stage'].explicit_update_gpu)
        drv.memcpy_dtoh(domain2.quantities['xmomentum'].explicit_update,
                        domain2.quantities['xmomentum'].explicit_update_gpu)
        drv.memcpy_dtoh(domain2.quantities['ymomentum'].explicit_update,
                        domain2.quantities['ymomentum'].explicit_update_gpu)
        
        drv.memcpy_dtoh(domain2.max_speed,
                        domain2.max_speed_gpu)
        
        drv.memcpy_dtoh(domain2.timestep_array,
                        domain2.timestep_array_gpu)
        
        #b = numpy.argsort(domain2.timestep_array)
        #domain2.flux_timestep = domain2.timestep_array[b[0]] 
        domain2.flux_timestep = numpy.min(domain2.timestep_array)
    else:
        if domain2.compute_fluxes_method == 'original':
            pass
            #compute_fluxes_central_structure_cuda(domain2)
            #gravity_c(domain2)
        elif domain2.compute_fluxes_method == 'wb_1':
            pass
        elif domain2.compute_fluxes_method == 'wb_2':
            #compute_fluxes_central_structure_cuda(domain2)
            compute_fluxes_central_structure_cuda(
                    domain2, 
                    parallelFlag=1, 
	        		name= "compute_fluxes_central_structure_CUDA")
        	#compute_fluxes_central_structure_cuda(
            #       domain2, parallelFlag=1, \
	    	#		name= "_flux_function_central_2")
            #gravity_wb_c(domain2)
        elif domain2.compute_fluxes_method == 'wb_3':
            pass
        else:
            raise Exception('unknown compute_fluxes_method')



    """
    ANUGA C Function
    """
    print "ANUGA C Function"


    if domain1.compute_fluxes_method == 'original':
        from shallow_water_ext import compute_fluxes_ext_central_structure
        from shallow_water_ext import gravity as gravity_c
        
        print "original"
        domain1.flux_timestep = compute_fluxes_ext_central_structure(domain1)
        #gravity_c(domain1)
        
    elif domain1.compute_fluxes_method == 'wb_1':
        from shallow_water_ext import compute_fluxes_ext_wb
        from shallow_water_ext import gravity as gravity_c

        print "wb_1"
        domain1.flux_timestep = compute_fluxes_ext_wb(domain1)
        #gravity_c(domain1)

    elif domain1.compute_fluxes_method == 'wb_2':
        from shallow_water_ext import compute_fluxes_ext_central_structure
        from shallow_water_ext import gravity_wb as gravity_wb_c

        #print "wb_2"
        t1 = time.time()
        domain1.flux_timestep = compute_fluxes_ext_central_structure(domain1)
        #gravity_wb_c(domain1)
        print "\n C compute flux duration is %lf\n" % (time.time()-t1)
    elif domain1.compute_fluxes_method == 'wb_3':
        from shallow_water_ext import compute_fluxes_ext_wb_3
        from shallow_water_ext import gravity_wb as gravity_wb_c

        print "wb_3"
        domain1.flux_timestep = compute_fluxes_ext_wb_3(domain1)
        #gravity_wb_c(domain1)

    else:
        raise Exception('unknown compute_fluxes_method')
    


    print numpy.allclose( domain2.quantities['stage'].explicit_update,
            domain1.quantities['stage'].explicit_update)
    print numpy.allclose( domain2.quantities['xmomentum'].explicit_update,
            domain1.quantities['xmomentum'].explicit_update)
    print numpy.allclose( domain2.quantities['ymomentum'].explicit_update,
            domain1.quantities['ymomentum'].explicit_update)
    # Check correctness
    print "\n~~~~~~~~~~~~~ domain 2 ~~~~~~~~~~~~"
    print "******* flux_timestep : %lf %lf %d" % \
            (domain1.flux_timestep, domain2.flux_timestep, \
                domain1.flux_timestep == domain2.flux_timestep)

    print "******* max_speed"
    counter = 0
    for i in range(domain1.number_of_elements):
        if ( domain1.max_speed[i] != domain2.max_speed[i]):
            counter += 1
            if counter < 10:
                print i, domain1.max_speed[i], domain2.max_speed[i]
    print "---------> # of differences: %d" % counter


    #print "******* stage_explicit_update"
    ##print domain1.quantities['stage'].explicit_update
    ##print domain2.quantities['stage'].explicit_update
    #counter = 0
    #for i in range(domain1.number_of_elements):
    #    if approx_cmp( domain1.quantities['stage'].explicit_update[i] ,
    #            domain2.quantities['stage'].explicit_update[i]):
    #        counter += 1
    #print "---------> # of differences: %d" % counter


    #print "******* xmom_explicit_update"
    ##print domain1.quantities['xmomentum'].explicit_update
    ##print domain2.quantities['xmomentum'].explicit_update
    #counter = 0
    #for i in range(domain1.number_of_elements):
    #    if approx_cmp( domain1.quantities['xmomentum'].explicit_update[i] ,
    #            domain2.quantities['xmomentum'].explicit_update[i]):
    #        counter += 1
    #        if counter < 10:
    #            print i, domain1.quantities['xmomentum'].explicit_update[i], \
    #                    domain2.quantities['xmomentum'].explicit_update[i]
    #print "---------> # of differences: %d" % counter

    #
    #print "******* ymom_explicit_update"
    ##print domain1.quantities['ymomentum'].explicit_update
    ##print domain2.quantities['ymomentum'].explicit_update
    #counter = 0
    #for i in range(domain1.number_of_elements):
    #    if approx_cmp( domain1.quantities['ymomentum'].explicit_update[i],
    #            domain2.quantities['ymomentum'].explicit_update[i]):
    #        counter += 1
    #print "---------> # of differences: %d" % counter




    """
        Following use the ANUGA original sequential C function 
        without already_computed_flux barrie and run in GPU 
        to check the result of the result generated by parallel version
    """


    #print "\n~~~~~~~~~~~~~~~~~~~ domain 3 ~~~~~~~~~~~~~~~~~"
    #domain3 = generate_merimbula_domain()

    #compute_fluxes_central_structure_cuda(domain3, parallelFlag=1, \
	#		name= "_flux_function_central_2")
    ##compute_fluxes_central_structure_cuda(domain3, parallelFlag=1, \
	##		name= "compute_fluxes_central_structure_MeCo")
    #counter_stage = 0
    #counter_xmom = 0
    #counter_ymom = 0
    #for i in range(domain3.number_of_elements):
    #    if approx_cmp(domain1.quantities['stage'].explicit_update[i],
    #            domain3.quantities['stage'].explicit_update[i]):
    #        counter_stage += 1
    #    if approx_cmp(domain1.quantities['xmomentum'].explicit_update[i],
    #           domain3.quantities['xmomentum'].explicit_update[i]):
    #        counter_xmom += 1
    #    if approx_cmp( domain1.quantities['ymomentum'].explicit_update[i] ,
    #            domain3.quantities['ymomentum'].explicit_update[i]):
    #        counter_ymom += 1

    #print "*******  of diff %d, %d, %d" % \
	#		(counter_stage, counter_xmom, counter_ymom)
    


    """
        Following use the python serial without already_computed_flux
    """


    if testing_python_version:
        print "\n~~~~~~~~~~~~~ domain 4 ~~~~~~~~~~~~"
        domain4 = generate_merimbula_domain()

        sort_domain(domain4)

        compute_fluxes_central_structure_single(domain4)


        counter_stage = 0
        counter_xmom = 0
        counter_ymom = 0
        counter_max = 0
        for i in range(domain1.number_of_elements):
            if domain1.quantities['stage'].explicit_update[i] != \
                    domain4.quantities['stage'].explicit_update[i]:
                counter_stage += 1

            if domain1.quantities['xmomentum'].explicit_update[i] != \
                    domain4.quantities['xmomentum'].explicit_update[i] :
                counter_xmom += 1

            if domain1.quantities['ymomentum'].explicit_update[i] != \
                    domain4.quantities['ymomentum'].explicit_update[i]:
                counter_ymom += 1

            if domain1.max_speed[i] != domain4.max_speed[i]:
                counter_max += 1
                if i < 5:
                    print i, domain1.max_speed[i], domain4.max_speed[i]
        print "-------  of diff with domain1 :%d, %d, %d, %d\n" % \
                (counter_stage, counter_xmom, counter_ymom, counter_max)
        


        counter_stage = 0
        counter_xmom = 0
        counter_ymom = 0
        counter_max = 0
        for i in range(domain2.number_of_elements):
            if domain2.quantities['stage'].explicit_update[i] != \
                    domain4.quantities['stage'].explicit_update[i]:
                counter_stage += 1

            if domain2.quantities['xmomentum'].explicit_update[i] != \
                    domain4.quantities['xmomentum'].explicit_update[i]:
                counter_xmom += 1

            if domain2.quantities['ymomentum'].explicit_update[i] != \
                    domain4.quantities['ymomentum'].explicit_update[i]:
                counter_ymom += 1

            if domain2.max_speed[i] != domain4.max_speed[i]:
                counter_max += 1
                if i < 5:
                    print i, domain2.max_speed[i], domain4.max_speed[i]
        print "-------  of diff with domain2 :%d, %d, %d, %d\n" % \
                (counter_stage, counter_xmom, counter_ymom, counter_max)
	    

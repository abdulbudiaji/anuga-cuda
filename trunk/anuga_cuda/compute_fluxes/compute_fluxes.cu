//#define UNSORTED_DOMAIN
#define USING_MULTI_FUNCTION

#ifdef UNSORTED_DOMAIN
__device__ void spe_bubble_sort(int* _list , long* neighbours, int k)
{
    int temp;
    if ( neighbours[_list[2]]>=0 and neighbours[ _list[2] ] < k and (neighbours[_list[1]]<0 or neighbours[_list[2]]<neighbours[_list[1]]) )
    {
        temp = _list[2];
        _list[2] = _list[1];
        _list[1] = temp;
    }
    if ( neighbours[_list[1]]>=0 and neighbours[ _list[1] ] < k and (neighbours[_list[0]]<0 or neighbours[_list[1]]<neighbours[_list[0]]) )
    {
        temp = _list[1];
        _list[1] = _list[0];
        _list[0] = temp;
    }
    if ( neighbours[_list[2]]>=0 and neighbours[ _list[2] ] < k and (neighbours[_list[1]]<0 or neighbours[_list[2]]<neighbours[_list[1]]) )
    {
        temp = _list[2];
        _list[2] = _list[1];
        _list[1] = temp;
    }
}
#endif


#ifdef USING_MULTI_FUNCTION

#define Dtimestep 0
#define Depsilon 1
#define DH0 2
#define Dg 3
#define Doptimise_dry_cells 4
#define Dcall 5
#define Dnumber_of_elements


__device__ int _rotate(double *q, double n1, double n2) 
{
    /*Rotate the momentum component q (q[1], q[2])
      from x,y coordinates to coordinates based on normal vector (n1, n2).

      Result is returned in array 3x1 r
      To rotate in opposite direction, call rotate with (q, n1, -n2)

      Contents of q are changed by this function */


    double q1, q2;

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
        double limiting_threshold) 
{
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







/*****************************************/
/* The CUDA compute fluex function       */
/*****************************************/


__global__ void compute_fluxes_central_structure_CUDA(
        long N,
        //double * elements,
        double  g,
        double epsilon,
        double h0,
        double limiting_threshold,
        int optimise_dry_cells,

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
    const long k = 
            threadIdx.x+threadIdx.y*blockDim.x+
            (blockIdx.x+blockIdx.y*gridDim.x)*blockDim.x*blockDim.y;


    double max_speed, length, inv_area, zl, zr;

    //double h0 = elements[DH0] * elements[DH0]; // This ensures a good balance when h approaches H0.

    //double limiting_threshold = 10 * elements[DH0]; // Avoid applying limiter below this

    int i, m, n;
    int ki, nm = 0, ki2; // Index shorthands

    double ql[3], qr[3], edgeflux[3];

    if (k >= N)
        return;

#ifdef UNSORTED_DOMAIN
    int b[3]={0,1,2}, j;
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

        //if (elements[Doptimise_dry_cells]) {
        if (optimise_dry_cells){
            //if (fabs(ql[0] - zl) < elements[Depsilon] &&
            //        fabs(qr[0] - zr) < elements[Depsilon]) {
            if ( fabs(ql[0] - zl) < epsilon &&
                    fabs(qr[0] - zr) < epsilon) {
                max_speed = 0.0;
                continue;
            }
        }

        ki2 = 2 * ki; //k*6 + i*2


        _flux_function_central(ql, qr, zl, zr,
                normals[ki2], normals[ki2 + 1],
                //elements[Depsilon], h0, limiting_threshold, elements[Dg],
                epsilon, h0, limiting_threshold, g,
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
            //if (max_speed > elements[Depsilon]) {
            if ( max_speed > epsilon) {
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
/*          Sequential version           */
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
    const long k = 
            threadIdx.x+threadIdx.y*blockDim.x+
            (blockIdx.x+blockIdx.y*gridDim.x)*blockDim.x*blockDim.y;


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

#endif

/*****************************************/
/* Rearrange the domain variable order   */
/* so as to achieve memory corelasing    */
/* also combination all subfunctions     */
/*****************************************/

#define B blockDim.x*blockDim.y
#define T threadIdx.x+threadIdx.y*blockDim.x

//__global__ void _flux_function_central_2(
__global__ void compute_fluxes_central_structure_cuda_single(
        long N,
        double g,
        double epsilon,
        double h0,
        double limiting_threshold,
        int optimise_dry_cells,
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
    int j;

    double w_left, h_left, uh_left, vh_left, u_left;
    double w_right, h_right, uh_right, vh_right, u_right;
    double s_min, s_max, soundspeed_left, soundspeed_right;
    double denom, inverse_denominator, z;
    double temp;

    // Workspace (allocate once, use many)
    double q_left_rotated[3], q_right_rotated[3], flux_right[3], flux_left[3];


    //const int k = threadIdx.x + blockIdx.x * blockDim.x;
    const long k = 
            threadIdx.x+threadIdx.y*blockDim.x+
            (blockIdx.x+blockIdx.y*gridDim.x)*blockDim.x*blockDim.y;

    double q_left[3], q_right[3];
    double z_left,  z_right;
    double n1,  n2;
    double edgeflux[3],  max_speed;
    double length, inv_area;

    int i, m, n;
    //int B = blockDim.x*blockDim.y;
    //int T = threadIdx.x+threadIdx.y*blockDim.x;
    int ki, nm;

    __shared__ double sh_data[32*9];

    if ( k >= N)
        return;
#ifdef UNSORTED_DOMAIN
    int b[3]={0,1,2}, l;
    spe_bubble_sort( b, neighbours+k*3, k);
    for (l = 0; l < 3; l++) {
        i = b[l];
#else
    for (i=0; i<3; i++) {
#endif

#ifdef REARRANGED_DOMAIN
        ki = k + i*N;
#else
        ki = k * 3 + i;
#endif
        n = neighbours[ki];

        q_left[0] = stage_edge_values[ki];
        q_left[1] = xmom_edge_values[ki];
        q_left[2] = ymom_edge_values[ki];
        z_left = bed_edge_values[ki];

        if (n<0) {
            m= -n -1;

            q_right[0] = stage_boundary_values[m];
            q_right[1] = xmom_boundary_values[m];
            q_right[2] = ymom_boundary_values[m];
            z_right = z_left;
        } else {
            m = neighbour_edges[ki];
#ifdef REARRANGED_DOMAIN
            nm = n + m*N;
#else
            nm = n *3 + m;
#endif

            q_right[0] = stage_edge_values[nm];
            q_right[1] = xmom_edge_values[nm];
            q_right[2] = ymom_edge_values[nm];
            z_right = bed_edge_values[nm];
        }

        if(optimise_dry_cells) {
            if(fabs(q_left[0] - z_left) < epsilon &&
                    fabs(q_right[0] - z_right) < epsilon ) {
                max_speed = 0.0;
                continue;
            }
        }


#ifdef REARRANGED_DOMAIN
        n1 = normals[k + 2*i*N];
        n2 = normals[k + (2*i+1)*N];
#else
        n1 = normals[2*ki];
        n2 = normals[2*ki + 1];
#endif


        /////////////////////////////////////////////////////////


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
            //return 0.0;
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

        soundspeed_left = sqrt(g * h_left);
        soundspeed_right = sqrt(g * h_right);

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
            max_speed = 0.0;
        }
        else {
            inverse_denominator = 1.0 / denom;
            for (j = 0; j < 3; j++) {
                edgeflux[j] = s_max * flux_left[j] - s_min * flux_right[j];
                edgeflux[j] += s_max * s_min * (q_right_rotated[j] - q_left_rotated[j]);
                edgeflux[j] *= inverse_denominator;
            }

            // Maximal wavespeed
            max_speed = max(fabs(s_max), fabs(s_min));

            // Rotate back
            //_rotate(edgeflux, n1, -n2);
            temp = edgeflux[1];
            edgeflux[1] = n1* temp -n2*edgeflux[2];
            edgeflux[2] = n2*temp + n1*edgeflux[2];
        }

        //////////////////////////////////////////////////////

        length = edgelengths[ki];
        edgeflux[0] *= length;
        edgeflux[1] *= length;
        edgeflux[2] *= length;

        sh_data[T + i*B]= -edgeflux[0];
        sh_data[T + (i+3)*B]= -edgeflux[1];
        sh_data[T + (i+6)*B]= -edgeflux[2];
        __syncthreads();


        //stage_explicit_update[k] -= edgeflux[0];
        //xmom_explicit_update[k] -= edgeflux[1];
        //ymom_explicit_update[k] -= edgeflux[2];

        if (tri_full_flag[k] == 1) {
            if (max_speed > epsilon) {
                timestep[k] = min(timestep[k], radii[k]/ max_speed);

                if (n>=0){
                    timestep[k] = min(timestep[k], radii[n]/max_speed);
                }
            }
        }

    }


    inv_area = 1.0 / areas[k];
    //stage_explicit_update[k] *= inv_area;
    //xmom_explicit_update[k] *= inv_area;
    //ymom_explicit_update[k] *= inv_area;

    stage_explicit_update[k] = (sh_data[T] + sh_data[T+B]+sh_data[T+B*2]) * inv_area;

    xmom_explicit_update[k] = (sh_data[T+B*3] + sh_data[T+B*4]+sh_data[T+B*5]) * inv_area;

    ymom_explicit_update[k] = (sh_data[T+B*6] + sh_data[T+B*7]+sh_data[T+B*8]) * inv_area;

    max_speed_array[k] = max_speed;
}

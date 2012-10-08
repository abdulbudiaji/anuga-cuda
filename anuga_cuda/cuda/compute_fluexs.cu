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
    static double q_left_rotated[3], q_right_rotated[3], flux_right[3], flux_left[3];

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
                report_python_error(AT, "Discontinuous Elevation");
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

static long compute_fluxes_central_call = 1;

/*
 * elements[0] = _timestep
 * elements[1] = _epsilon
 * elements[2] = _H0
 * elements[3] = _g
 * elements[4] = _optimise_dry_cells
*/
double _compute_fluxes_central(int number_of_elements,
        double* elements,
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
            // Local variables
            double max_speed, length, inv_area, zl, zr;
            double h0 = elements[2]*elements[2]; // This ensures a good balance when h approaches H0.

            double limiting_threshold = 10 * elements[2]; // Avoid applying limiter below this
            // threshold for performance reasons.
            // See ANUGA manual under flux limiting
            int k, i, m, n;
            int ki, nm = 0, ki2; // Index shorthands


            // Workspace (making them static actually made function slightly slower (Ole))
            double ql[3], qr[3], edgeflux[3]; // Work array for summing up fluxes

            //static long call = 1; // Static local variable flagging already computed flux

            // Start computation
            call++; // Flag 'id' of flux calculation for this elements[0]

            // Set explicit_update to zero for all conserved_quantities.
            // This assumes compute_fluxes called before forcing terms
            memset((char*) stage_explicit_update, 0, number_of_elements * sizeof (double));
            memset((char*) xmom_explicit_update, 0, number_of_elements * sizeof (double));
            memset((char*) ymom_explicit_update, 0, number_of_elements * sizeof (double));

            // For all triangles
            for (k = 0; k < number_of_elements; k++) 
            {
                // Loop through neighbours and compute edge flux for each
                for (i = 0; i < 3; i++) {
                    ki = k * 3 + i; // Linear index to edge i of triangle k

                    if (already_computed_flux[ki] == call) {
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

                    if (elements[4]) {
                        // Check if flux calculation is necessary across this edge
                        // This check will exclude dry cells.
                        // This will also optimise cases where zl != zr as
                        // long as both are dry

                        if (fabs(ql[0] - zl) < elements[1] &&
                                fabs(qr[0] - zr) < elements[1]) {
                            // Cell boundary is dry

                            already_computed_flux[ki] = call; // #k Done
                            if (n >= 0) {
                                already_computed_flux[nm] = call; // #n Done
                            }

                            max_speed = 0.0;
                            continue;
                        }
                    }


                    if (fabs(zl - zr) > 1.0e-10) {
                        report_python_error(AT, "Discontinuous Elevation");
                        return 0.0;
                    }

                    // Outward pointing normal vector (domain.normals[k, 2*i:2*i+2])
                    ki2 = 2 * ki; //k*6 + i*2

                    // Edge flux computation (triangle k, edge i)
                    _flux_function_central(ql, qr, zl, zr,
                            normals[ki2], normals[ki2 + 1],
                            elements[1], h0, limiting_threshold, elements[3],
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

                    already_computed_flux[ki] = call; // #k Done


                    // Update neighbour n with same flux but reversed sign
                    if (n >= 0) {
                        stage_explicit_update[n] += edgeflux[0];
                        xmom_explicit_update[n] += edgeflux[1];
                        ymom_explicit_update[n] += edgeflux[2];

                        already_computed_flux[nm] = call; // #n Done
                    }

                    // Update elements[0] based on edge i and possibly neighbour n
                    if (tri_full_flag[k] == 1) {
                        if (max_speed > elements[1]) {
                            // Apply CFL condition for triangles joining this edge (triangle k and triangle n)

                            // CFL for triangle k
                            elements[0] = min(elements[0], radii[k] / max_speed);

                            if (n >= 0) {
                                // Apply CFL condition for neigbour n (which is on the ith edge of triangle k)
                                elements[0] = min(elements[0], radii[n] / max_speed);
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

            } // End triangle k

            return elements[0];
}

static long compute_fluxes_central_structure_call = 1;

double _compute_fluxes_central_structure(struct domain *D) {

    // Local variables
    double max_speed, length, inv_area, zl, zr;
    double timestep = 1.0e30;
    double h0 = D->H0 * D->H0; // This ensures a good balance when h approaches H0.

    double limiting_threshold = 10 * D->H0; // Avoid applying limiter below this
    // threshold for performance reasons.
    // See ANUGA manual under flux limiting
    int k, i, m, n;
    int ki, nm = 0, ki2; // Index shorthands


    // Workspace (making them static actually made function slightly slower (Ole))
    double ql[3], qr[3], edgeflux[3]; // Work array for summing up fluxes

    //static long call = 1; // Static local variable flagging already computed flux

    // Start computation
    call++; // Flag 'id' of flux calculation for this timestep


    timestep = D->evolve_max_timestep;

    // Set explicit_update to zero for all conserved_quantities.
    // This assumes compute_fluxes called before forcing terms
    memset((char*) D->stage_explicit_update, 0, D->number_of_elements * sizeof (double));
    memset((char*) D->xmom_explicit_update, 0, D->number_of_elements * sizeof (double));
    memset((char*) D->ymom_explicit_update, 0, D->number_of_elements * sizeof (double));

    // For all triangles
    for (k = 0; k < D->number_of_elements; k++) {
        // Loop through neighbours and compute edge flux for each
        for (i = 0; i < 3; i++) {
            ki = k * 3 + i; // Linear index to edge i of triangle k

            if (D->already_computed_flux[ki] == call) {
                // We've already computed the flux across this edge
                continue;
            }

            // Get left hand side values from triangle k, edge i
            ql[0] = D->stage_edge_values[ki];
            ql[1] = D->xmom_edge_values[ki];
            ql[2] = D->ymom_edge_values[ki];
            zl = D->bed_edge_values[ki];

            // Get right hand side values either from neighbouring triangle
            // or from boundary array (Quantities at neighbour on nearest face).
            n = D->neighbours[ki];
            if (n < 0) {
                // Neighbour is a boundary condition
                m = -n - 1; // Convert negative flag to boundary index

                qr[0] = D->stage_boundary_values[m];
                qr[1] = D->xmom_boundary_values[m];
                qr[2] = D->ymom_boundary_values[m];
                zr = zl; // Extend bed elevation to boundary
            } else {
                // Neighbour is a real triangle
                m = D->neighbour_edges[ki];
                nm = n * 3 + m; // Linear index (triangle n, edge m)

                qr[0] = D->stage_edge_values[nm];
                qr[1] = D->xmom_edge_values[nm];
                qr[2] = D->ymom_edge_values[nm];
                zr = D->bed_edge_values[nm];
            }

            // Now we have values for this edge - both from left and right side.

            if (D->optimise_dry_cells) {
                // Check if flux calculation is necessary across this edge
                // This check will exclude dry cells.
                // This will also optimise cases where zl != zr as
                // long as both are dry

                if (fabs(ql[0] - zl) < D->epsilon &&
                        fabs(qr[0] - zr) < D->epsilon) {
                    // Cell boundary is dry

                    D->already_computed_flux[ki] = call; // #k Done
                    if (n >= 0) {
                        D->already_computed_flux[nm] = call; // #n Done
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
                    D->normals[ki2], D->normals[ki2 + 1],
                    D->epsilon, h0, limiting_threshold, D->g,
                    edgeflux, &max_speed);


            // Multiply edgeflux by edgelength
            length = D->edgelengths[ki];
            edgeflux[0] *= length;
            edgeflux[1] *= length;
            edgeflux[2] *= length;


            // Update triangle k with flux from edge i
            D->stage_explicit_update[k] -= edgeflux[0];
            D->xmom_explicit_update[k] -= edgeflux[1];
            D->ymom_explicit_update[k] -= edgeflux[2];

            D->already_computed_flux[ki] = call; // #k Done


            // Update neighbour n with same flux but reversed sign
            if (n >= 0) {
                D->stage_explicit_update[n] += edgeflux[0];
                D->xmom_explicit_update[n] += edgeflux[1];
                D->ymom_explicit_update[n] += edgeflux[2];

                D->already_computed_flux[nm] = call; // #n Done
            }

            // Update timestep based on edge i and possibly neighbour n
            if (D->tri_full_flag[k] == 1) {
                if (max_speed > D->epsilon) {
                    // Apply CFL condition for triangles joining this edge (triangle k and triangle n)

                    // CFL for triangle k
                    timestep = min(timestep, D->radii[k] / max_speed);

                    if (n >= 0) {
                        // Apply CFL condition for neigbour n (which is on the ith edge of triangle k)
                        timestep = min(timestep, D->radii[n] / max_speed);
                    }

                    // Ted Rigby's suggested less conservative version
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
        inv_area = 1.0 / D->areas[k];
        D->stage_explicit_update[k] *= inv_area;
        D->xmom_explicit_update[k] *= inv_area;
        D->ymom_explicit_update[k] *= inv_area;


        // Keep track of maximal speeds
        D->max_speed[k] = max_speed;

    } // End triangle k


    return timestep;
}

double _compute_fluxes_central_wb(struct domain *D) {

    // Local variables
    double max_speed, length, inv_area;
    double timestep = 1.0e30;
    double h0 = D->H0 * D->H0; // This ensures a good balance when h approaches H0.

    double limiting_threshold = 10 * D->H0; // Avoid applying limiter below this
    // threshold for performance reasons.
    // See ANUGA manual under flux limiting
    int k, i, m, n, nn;
    int k3, k3i, k3i1, k3i2, k2i; // Index short hands

    int n3m = 0, n3m1, n3m2; // Index short hand for neightbours

    double hl, hl1, hl2;
    double hr, hr1, hr2;
    double zl, zl1, zl2;
    double zr;
    double wr;

    // Workspace (making them static actually made function slightly slower (Ole))
    double ql[3], qr[3], edgeflux[3]; // Work array for summing up fluxes

    static long call = 1; // Static local variable flagging already computed flux

    // Start computation
    call++; // Flag 'id' of flux calculation for this timestep


    timestep = D->evolve_max_timestep;

    // Set explicit_update to zero for all conserved_quantities.
    // This assumes compute_fluxes called before forcing terms
    memset((char*) D->stage_explicit_update, 0, D->number_of_elements * sizeof (double));
    memset((char*) D->xmom_explicit_update, 0, D->number_of_elements * sizeof (double));
    memset((char*) D->ymom_explicit_update, 0, D->number_of_elements * sizeof (double));

    // For all triangles
    for (k = 0; k < D->number_of_elements; k++) {
        // Loop through neighbours and compute edge flux for each
        for (i = 0; i < 3; i++) {
            k3 = 3 * k;
            k3i = k3 + i; // Linear index to edge i of triangle k
            k3i1 = k3 + (i + 1) % 3;
            k3i2 = k3 + (i + 2) % 3;

            if (D->already_computed_flux[k3i] == call) {
                // We've already computed the flux across this edge
                continue;
            }


            // Get the inside values at the vertices from the triangle k, edge i
            ql[0] = D->stage_edge_values[k3i];
            ql[1] = D->xmom_edge_values[k3i];
            ql[2] = D->ymom_edge_values[k3i];

            zl = D->bed_edge_values[k3i];

            zl1 = D->bed_vertex_values[k3i1];
            zl2 = D->bed_vertex_values[k3i2];


            hl = D->stage_edge_values[k3i] - zl;

            hl1 = D->stage_vertex_values[k3i1] - zl1;
            hl2 = D->stage_vertex_values[k3i2] - zl2;


            // Get right hand side values either from neighbouring triangle
            // or from boundary array (Quantities at neighbour on nearest face).
            n = D->neighbours[k3i];
            if (n < 0) {
                // Neighbour is a boundary condition
                nn = -n - 1; // Convert negative flag to boundary index
                m = 0;

                qr[0] = D->stage_boundary_values[nn];
                qr[1] = D->xmom_boundary_values[nn];
                qr[2] = D->ymom_boundary_values[nn];

                zr = zl; // Extend bed elevation to boundary and assume flat stage

                wr = D->stage_boundary_values[nn];

                hr = wr - zr;
                hr1 = wr - zl2;
                hr2 = wr - zl1;


            } else {
                // Neighbour is a real triangle
                m = D->neighbour_edges[k3i];
                n3m = 3 * n + m; // Linear index (triangle n, edge m)
                n3m1 = 3 * n + (m + 1) % 3;
                n3m2 = 3 * n + (m + 2) % 3;

                qr[0] = D->stage_edge_values[n3m];
                qr[1] = D->xmom_edge_values[n3m];
                qr[2] = D->ymom_edge_values[n3m];

                zr = D->bed_edge_values[n3m];
                hr = D->stage_edge_values[n3m] - zr;

                hr1 = D->stage_vertex_values[n3m2] - D->bed_vertex_values[n3m2];
                hr2 = D->stage_vertex_values[n3m1] - D->bed_vertex_values[n3m1];

            }

            // Now we have values for this edge - both from left and right side.

            if (D->optimise_dry_cells) {
                // Check if flux calculation is necessary across this edge
                // This check will exclude dry cells.
                // This will also optimise cases where zl != zr as
                // long as both are dry

                if (fabs(ql[0] - zl) < D->epsilon &&
                        fabs(qr[0] - zr) < D->epsilon) {
                    // Cell boundary is dry

                    D->already_computed_flux[k3i] = call; // #k Done
                    if (n >= 0) {
                        D->already_computed_flux[n3m] = call; // #n Done
                    }

                    max_speed = 0.0;
                    continue;
                }
            }


            if (fabs(zl - zr) > 1.0e-10) {
                report_python_error(AT, "Discontinuous Elevation");
                return 0.0;
            }

            // Outward pointing normal vector (domain.normals[k, 2*i:2*i+2])
            k2i = 2 * k3i; //k*6 + i*2



            /*
               _flux_function_central(ql, qr, zl, zr,
               D->normals[k2i], D->normals[k2i + 1],
               D->epsilon, h0, limiting_threshold, D->g,
               edgeflux, &max_speed);
             */

            // Edge flux computation (triangle k, edge i)


            //printf("========== wb_1 ==============\n");

            //printf("E_left %i %i \n",k , i);
            //printf("E_right %i %i \n",n, m);

            _flux_function_central_wb(ql, qr,
                    zl, hl, hl1, hl2,
                    zr, hr, hr1, hr2,
                    D->normals[k2i], D->normals[k2i + 1],
                    D->epsilon, h0, limiting_threshold, D->g,
                    edgeflux, &max_speed);



            // Multiply edgeflux by edgelength
            length = D->edgelengths[k3i];
            edgeflux[0] *= length;
            edgeflux[1] *= length;
            edgeflux[2] *= length;


            // Update triangle k with flux from edge i
            D->stage_explicit_update[k] -= edgeflux[0];
            D->xmom_explicit_update[k] -= edgeflux[1];
            D->ymom_explicit_update[k] -= edgeflux[2];

            D->already_computed_flux[k3i] = call; // #k Done


            // Update neighbour n with same flux but reversed sign
            if (n >= 0) {
                D->stage_explicit_update[n] += edgeflux[0];
                D->xmom_explicit_update[n] += edgeflux[1];
                D->ymom_explicit_update[n] += edgeflux[2];

                D->already_computed_flux[n3m] = call; // #n Done
            }

            // Update timestep based on edge i and possibly neighbour n
            if (D->tri_full_flag[k] == 1) {
                if (max_speed > D->epsilon) {
                    // Apply CFL condition for triangles joining this edge (triangle k and triangle n)

                    // CFL for triangle k
                    timestep = min(timestep, D->radii[k] / max_speed);

                    if (n >= 0) {
                        // Apply CFL condition for neigbour n (which is on the ith edge of triangle k)
                        timestep = min(timestep, D->radii[n] / max_speed);
                    }

                    // Ted Rigby's suggested less conservative version
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
        inv_area = 1.0 / D->areas[k];
        D->stage_explicit_update[k] *= inv_area;
        D->xmom_explicit_update[k] *= inv_area;
        D->ymom_explicit_update[k] *= inv_area;


        // Keep track of maximal speeds
        D->max_speed[k] = max_speed;

    } // End triangle k


    return timestep;
}

__global__ void compute_fluxes_central_wb_3(struct domain *D) {

    // Local variables
    double max_speed, length, inv_area;
    double timestep = 1.0e30;
    double h0 = D->H0 * D->H0; // This ensures a good balance when h approaches H0.

    double limiting_threshold = 10 * D->H0; // Avoid applying limiter below this
    // threshold for performance reasons.
    // See ANUGA manual under flux limiting
    int k, i, m, n;
    int k3, k3i, k2i; // Index short hands
    int n3m = 0; // Index short hand for neightbours

    struct edge E_left;
    struct edge E_right;

    double edgeflux[3]; // Work array for summing up fluxes

    static long call = 1; // Static local variable flagging already computed flux

    // Start computation
    call++; // Flag 'id' of flux calculation for this timestep


    timestep = D->evolve_max_timestep;

    // Set explicit_update to zero for all conserved_quantities.
    // This assumes compute_fluxes called before forcing terms
    memset((char*) D->stage_explicit_update, 0, D->number_of_elements * sizeof (double));
    memset((char*) D->xmom_explicit_update,  0, D->number_of_elements * sizeof (double));
    memset((char*) D->ymom_explicit_update,  0, D->number_of_elements * sizeof (double));

    // For all triangles
    for (k = 0; k < D->number_of_elements; k++) {
        // Loop through neighbours and compute edge flux for each
        for (i = 0; i < 3; i++) {
            k3 = 3 * k;
            k3i = k3 + i; // Linear index to edge i of triangle k

            if (D->already_computed_flux[k3i] == call) {
                // We've already computed the flux across this edge
                continue;
            }

            // Get the "left" values at the edge vertices and midpoint
            // from the triangle k, edge i
            get_edge_data(&E_left, D, k, i);


            // Get right hand side values either from neighbouring triangle
            // or from boundary array (Quantities at neighbour on nearest face).
            n = D->neighbours[k3i];
            if (n < 0) {
                // Neighbour is a boundary condition
                m = -n - 1; // Convert negative flag to boundary index

                E_right.cell_id = n;
                E_right.edge_id = 0;

                // Midpoint Values provided by the boundary conditions
                E_right.w = D->stage_boundary_values[m];
                E_right.uh = D->xmom_boundary_values[m];
                E_right.vh = D->ymom_boundary_values[m];

                // Set bed and height as equal to neighbour
                E_right.z = E_left.z;
                E_right.h = E_right.w - E_right.z;

                // vertex values same as midpoint values
                E_right.w1  = E_right.w;
                E_right.h1  = E_right.h;

                E_right.uh1 = E_left.uh2;
                E_right.vh1 = E_left.vh2;
                E_right.z1  = E_left.z;

                E_right.w2  = E_right.w;
                E_right.h2  = E_right.h;

                E_right.uh2 = E_right.uh;
                E_right.vh2 = E_right.vh;
                E_right.z2  = E_left.z;

            } else {
                // Neighbour is a real triangle
                m = D->neighbour_edges[k3i];
                n3m  = 3*n+m;

                get_edge_data(&E_right, D, n, m);

            }



            // Now we have values for this edge - both from left and right side.

            if (D->optimise_dry_cells) {
                // Check if flux calculation is necessary across this edge
                // This check will exclude dry cells.
                // This will also optimise cases where zl != zr as
                // long as both are dry

                if (fabs(E_left.h) < D->epsilon &&
                        fabs(E_right.h) < D->epsilon) {
                    // Cell boundary is dry

                    D->already_computed_flux[k3i] = call; // #k Done
                    if (n >= 0) {
                        D->already_computed_flux[n3m] = call; // #n Done
                    }

                    max_speed = 0.0;
                    continue;
                }
            }


            if (fabs(E_left.z - E_right.z) > 1.0e-10) {
                report_python_error(AT, "Discontinuous Elevation");
                return 0.0;
            }

            // Outward pointing normal vector (domain.normals[k, 2*i:2*i+2])
            k2i = 2 * k3i; //k*6 + i*2



            // Edge flux computation (triangle k, edge i)
            _flux_function_central_wb_3(&E_left, &E_right,
                    D->normals[k2i], D->normals[k2i + 1],
                    D->epsilon, h0, limiting_threshold, D->g,
                    edgeflux, &max_speed);



            // Multiply edgeflux by edgelength
            length = D->edgelengths[k3i];
            edgeflux[0] *= length;
            edgeflux[1] *= length;
            edgeflux[2] *= length;


            // Update triangle k with flux from edge i
            D->stage_explicit_update[k] -= edgeflux[0];
            D->xmom_explicit_update[k]  -= edgeflux[1];
            D->ymom_explicit_update[k]  -= edgeflux[2];

            D->already_computed_flux[k3i] = call; // #k Done

            //printf("k i n m %i %i %i %i \n",k,i,n,m);

            // Update neighbour n with same flux but reversed sign
            if (n >= 0) {

                D->stage_explicit_update[n] += edgeflux[0];
                D->xmom_explicit_update[n]  += edgeflux[1];
                D->ymom_explicit_update[n]  += edgeflux[2];


                D->already_computed_flux[n3m] = call; // #n Done
            }

            // Update timestep based on edge i and possibly neighbour n
            if (D->tri_full_flag[k] == 1) {
                if (max_speed > D->epsilon) {
                    // Apply CFL condition for triangles joining this edge (triangle k and triangle n)

                    // CFL for triangle k
                    timestep = min(timestep, D->radii[k] / max_speed);

                    if (n >= 0) {
                        // Apply CFL condition for neigbour n (which is on the ith edge of triangle k)
                        timestep = min(timestep, D->radii[n] / max_speed);
                    }

                    // Ted Rigby's suggested less conservative version
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
        inv_area = 1.0 / D->areas[k];
        D->stage_explicit_update[k] *= inv_area;
        D->xmom_explicit_update[k] *= inv_area;
        D->ymom_explicit_update[k] *= inv_area;


        // Keep track of maximal speeds
        D->max_speed[k] = max_speed;

    } // End triangle k


    //return timestep;
}

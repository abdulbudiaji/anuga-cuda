

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
    int j, k;

    double w_left, h_left, uh_left, vh_left, u_left;
    double w_right, h_right, uh_right, vh_right, u_right;
    double s_min, s_max, soundspeed_left, soundspeed_right;
    double denom, inverse_denominator, z;
    double temp;

    // Workspace (allocate once, use many)
    double q_left_rotated[3], q_right_rotated[3], flux_right[3], flux_left[3];



    double q_left[3], q_right[3];
    double z_left,  z_right;
    double n1,  n2;
    double edgeflux[3],  max_speed, max_speed_total;
    double length, inv_area;

    int i, m, n;
    int ki, nm;

//#ifdef UNSORTED_DOMAIN
//    int b[3]={0,1,2}, l;
//    spe_bubble_sort( b, neighbours+k*3, k);
//    for (l = 0; l < 3; l++) {
//        i = b[l];
//#else
//    for (i=0; i<3; i++) {
//#endif
//
//#ifdef REARRANGED_DOMAIN
//        ki = k + i*N;
//#else
//        ki = k * 3 + i;
//#endif
    for (k=0; k < N; k++)
    {
        max_speed_total = 0;
        for (i=0; i< 3; i++)
        {
            ki= k*3 + i;


            q_left[0] = stage_edge_values[ki];
            q_left[1] = xmom_edge_values[ki];
            q_left[2] = ymom_edge_values[ki];
            z_left = bed_edge_values[ki];

            n = neighbours[ki];
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


            stage_explicit_update[k] -= edgeflux[0];
            xmom_explicit_update[k] -= edgeflux[1];
            ymom_explicit_update[k] -= edgeflux[2];

            if (tri_full_flag[k] == 1) {
                if (max_speed > epsilon) {
                    timestep[k] = min(timestep[k], radii[k]/ max_speed);

                    if (n>=0){
                        timestep[k] = min(timestep[k], radii[n]/max_speed);
                    }
                }
            }
            max_speed_total = max(max_speed_total, max_speed);
        }



        inv_area = 1.0 / areas[k];
        stage_explicit_update[k] *= inv_area;
        xmom_explicit_update[k] *= inv_area;
        ymom_explicit_update[k] *= inv_area;

        max_speed_array[k] = max_speed_total;
    }
}

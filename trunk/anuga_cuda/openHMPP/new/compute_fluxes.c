
// OpenHMPP ANUGA compute_fluxes
//
// Zhe Weng 2013


//#include<stdio.h>
//#include<stdlib.h>
//#include<string.h>
//#include<math.h>
//#include<assert.h>

#include "hmpp_fun.h"

#ifndef TOLERANCE
#define TOLERANCE 1e-8
#endif


//double max(double a, double b)
//{ return (a >= b) ? a : b; }
//
//double min(double a, double b)
//{ return (a <= b) ? a : b; }



double sqrt(double x)
{
    double mx = 32000;
    double mn = 0;
    while(mx - mn > 1e-9)
    {
        double md = (mx + mn) / 2;
        if(md * md > x)
            mx = md;
        else
            mn = md;
    }
    return mx;
}



/*
void spe_bubble_sort(int* _list , long* neighbours, int k)
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
*/

int _rotate(double *q, double n1, double n2) 
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


double _compute_speed(double *uh,
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


int _flux_function_central(
#ifdef USING_ORIGINAL_CUDA
        double *q_left, double *q_right,
        double z_left, double z_right,
        double n1, double n2,
        double epsilon,
        double h0,
        double limiting_threshold,
        double g,
        double *edgeflux, 
        double *max_speed) 
#else
        double q_left0, 
        double q_left1, 
        double q_left2, 
        double q_right0,
        double q_right1,
        double q_right2,
        double z_left, 
        double z_right,
        double n1, 
        double n2,
        double epsilon,
        double h0,
        double limiting_threshold,
        double g,
        double *edgeflux0, 
        double *edgeflux1, 
        double *edgeflux2, 
        double *max_speed) 
#endif
{
#ifdef USING_ORIGINAL_CUDA
    int i;
#endif

    double w_left, h_left, uh_left, vh_left, u_left;
    double w_right, h_right, uh_right, vh_right, u_right;
    double s_min, s_max, soundspeed_left, soundspeed_right;
    double denom, inverse_denominator, z;
    double temp;

    // Workspace (allocate once, use many)
    double q_left_rotated[3], q_right_rotated[3], flux_right[3], flux_left[3];

    // Copy conserved quantities to protect from modification
#ifdef USING_ORIGINAL_CUDA
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
#else
    q_left_rotated[0] = q_left0;
    q_right_rotated[0] = q_right0;
    q_left_rotated[1] = q_left1;
    q_right_rotated[1] = q_right1;
    q_left_rotated[2] = q_left2;
    q_right_rotated[2] = q_right2;

    // Align x- and y-momentum with x-axis
    //_rotate(q_left_rotated, n1, n2);
    q_left_rotated[1] = n1*q_left1 + n2*q_left2;
    q_left_rotated[2] = -n2*q_left1 + n1*q_left2;

    //_rotate(q_right_rotated, n1, n2);
    q_right_rotated[1] = n1*q_right1 + n2*q_right2;
    q_right_rotated[2] = -n2*q_right1 + n1*q_right2;
#endif


    if (fabs(z_left - z_right) > 1.0e-10) {
        //report_python_error(AT, "Discontinuous Elevation");
        return 0.0;
    }
    z = 0.5 * (z_left + z_right); // Average elevation values.

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

    s_max = fmax(u_left + soundspeed_left, u_right + soundspeed_right);
    if (s_max < 0.0) {
        s_max = 0.0;
    }

    s_min = fmin(u_left - soundspeed_left, u_right - soundspeed_right);
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
        //memset(edgeflux, 0, 3 * sizeof (double));
        *max_speed = 0.0;
    }
    else {
        inverse_denominator = 1.0 / denom;
#ifdef USING_ORIGINAL_CUDA
        for (i = 0; i < 3; i++) {
            edgeflux[i] = s_max * flux_left[i] - s_min * flux_right[i];
            edgeflux[i] += s_max *s_min *(q_right_rotated[i] -q_left_rotated[i]);
            edgeflux[i] *= inverse_denominator;
        }
#else
        *edgeflux0 = s_max * flux_left[0] - s_min * flux_right[0];
        *edgeflux0 += s_max *s_min *(q_right_rotated[0] -q_left_rotated[0]);
        *edgeflux0 *= inverse_denominator;

        *edgeflux1 = s_max * flux_left[1] - s_min * flux_right[1];
        *edgeflux1 += s_max *s_min *(q_right_rotated[1] -q_left_rotated[1]);
        *edgeflux1 *= inverse_denominator;

        *edgeflux2 = s_max * flux_left[2] - s_min * flux_right[2];
        *edgeflux2 += s_max *s_min *(q_right_rotated[2] -q_left_rotated[2]);
        *edgeflux2 *= inverse_denominator;
#endif

        // Maximal wavespeed
        *max_speed = fmax(fabs(s_max), fabs(s_min));

        // Rotate back
        //_rotate(edgeflux, n1, -n2);
#ifdef USING_ORIGINAL_CUDA
        temp = edgeflux[1];
        edgeflux[1] = n1* temp -n2*edgeflux[2];
        edgeflux[2] = n2*temp + n1*edgeflux[2];
#else
        temp = *edgeflux1;
        *edgeflux1 = n1* temp -n2* *edgeflux2;
        *edgeflux2 = n2*temp + n1* *edgeflux2;
#endif
    }

    return 0;
}







/*****************************************/
/* The CUDA compute fluex function       */
/*****************************************/


//#pragma hmpp cf_central codelet, target=CUDA args[*].transfer=atcall
void compute_fluxes_central_structure_CUDA(
        int N,
        int N3,
        int N6,
        int N2,

        double timestep[N],
        long neighbours[N3],
        long neighbour_edges[N3],
        double normals[N6],
        double edgelengths[N3],
        double radii[N],
        double areas[N],
        long tri_full_flag[N],
        double stage_edge_values[N3],
        double xmom_edge_values[N3],
        double ymom_edge_values[N3],
        double bed_edge_values[N3],
        double stage_boundary_values[N2],
        double xmom_boundary_values[N2],
        double ymom_boundary_values[N2],
        double stage_explicit_update[N],
        double xmom_explicit_update[N],
        double ymom_explicit_update[N],
        double max_speed_array[N],

        double evolve_max_timestep,
        double g,
        double epsilon,
        double h0,
        double limiting_threshold,
        int optimise_dry_cells)
{
    int k;
    for(k=0; k<N; k++)
    {
        int i, m, n;
        int ki, nm = 0, ki2;

        double max_speed, max_speed_total=0 , length, inv_area, zl, zr;
        timestep[k] = evolve_max_timestep;

#ifdef USING_ORIGINAL_CUDA
        double ql[3], qr[3], edgeflux[3];
#else
        double ql0, ql1, ql2,
               qr0, qr1, qr2,
               edgeflux0, edgeflux1, edgeflux2;
#endif

        // Loop through neighbours and compute edge flux for each
        for (i = 0; i < 3; i++) 
        {
            ki = k * 3 + i; // Linear index to edge i of triangle k

            n = neighbours[ki];

#ifdef USING_ORIGINAL_CUDA
            ql[0] = stage_edge_values[ki];
            ql[1] = xmom_edge_values[ki];
            ql[2] = ymom_edge_values[ki];
#else   
            ql0 = stage_edge_values[ki];
            ql1 = xmom_edge_values[ki];
            ql2 = ymom_edge_values[ki];
#endif
            zl = bed_edge_values[ki];

            if (n < 0) {
                m = -n - 1; // Convert negative flag to boundary index

#ifdef USING_ORIGINAL_CUDA
                qr[0] = stage_boundary_values[m];
                qr[1] = xmom_boundary_values[m];
                qr[2] = ymom_boundary_values[m];
#else   
                qr0 = stage_boundary_values[m];
                qr1 = xmom_boundary_values[m];
                qr2 = ymom_boundary_values[m];
#endif
                zr = zl; // Extend bed elevation to boundary
            } else {
                m = neighbour_edges[ki];
                nm = n * 3 + m; // Linear index (triangle n, edge m)

#ifdef USING_ORIGINAL_CUDA
                qr[0] = stage_edge_values[nm];
                qr[1] = xmom_edge_values[nm];
                qr[2] = ymom_edge_values[nm];
#else   
                qr0 = stage_edge_values[nm];
                qr1 = xmom_edge_values[nm];
                qr2 = ymom_edge_values[nm];
#endif
                zr = bed_edge_values[nm];
            }

            if (optimise_dry_cells){
#ifdef USING_ORIGINAL_CUDA
                if ( fabs(ql[0] - zl) < epsilon &&
                        fabs(qr[0] - zr) < epsilon) {
#else   
                if ( fabs(ql0 - zl) < epsilon &&
                        fabs(qr0 - zr) < epsilon) {
#endif
                    max_speed = 0.0;
                    continue;
                }
            }

            ki2 = 2 * ki; //k*6 + i*2
#ifdef USING_ORIGINAL_CUDA
            _flux_function_central(ql, qr, zl, zr,
                    normals[ki2], normals[ki2 + 1],
                    epsilon, h0, limiting_threshold, g,
                    edgeflux, &max_speed);

            length = edgelengths[ki];
            edgeflux[0] *= length;
            edgeflux[1] *= length;
            edgeflux[2] *= length;


            stage_explicit_update[k] -= edgeflux[0];
            xmom_explicit_update[k] -= edgeflux[1];
            ymom_explicit_update[k] -= edgeflux[2];
#else   
            _flux_function_central(
                    ql0, ql1, ql2, 
                    qr0, qr1, qr2,
                    zl, 
                    zr,
                    normals[ki2], normals[ki2 + 1],
                    epsilon, h0, limiting_threshold, g,
                    &edgeflux0,
                    &edgeflux1,
                    &edgeflux2,
                    &max_speed);

            length = edgelengths[ki];
            edgeflux0 *= length;
            edgeflux1 *= length;
            edgeflux2 *= length;


            stage_explicit_update[k] -= edgeflux0;
            xmom_explicit_update[k] -= edgeflux1;
            ymom_explicit_update[k] -= edgeflux2;
#endif
            if (tri_full_flag[k] == 1) {
                //if (max_speed > elements[Depsilon]) {
                if ( max_speed > epsilon) {
                    timestep[k] = fmin(timestep[k], radii[k] / max_speed);
                    if (n >= 0) {
                        timestep[k] = fmin(timestep[k], radii[n] / max_speed);
                    }
                }
            }

            if (n < 0 ||  n > k){
                max_speed_total = fmax(max_speed_total, max_speed);
            }
        } // End edge i (and neighbour n)

        inv_area = 1.0 / areas[k];
        stage_explicit_update[k] *= inv_area;
        xmom_explicit_update[k] *= inv_area;
        ymom_explicit_update[k] *= inv_area;

        max_speed_array[k] =  max_speed_total;
    }
}



//#pragma hmpp cf_central_single codelet, target=CUDA args[*].transfer=atcall
void compute_fluxes_central_structure_cuda_single(
        int N,
        int N3,
        int N6,
        int N2,

        double timestep[N],
        int neighbours[N3],
        int neighbour_edges[N3],
        double normals[N6],
        double edgelengths[N3],
        double radii[N],
        double areas[N],
        int tri_full_flag[N],
        double stage_edge_values[N3],
        double xmom_edge_values[N3],
        double ymom_edge_values[N3],
        double bed_edge_values[N3],
        double stage_boundary_values[N2],
        double xmom_boundary_values[N2],
        double ymom_boundary_values[N2],
        double stage_explicit_update[N],
        double xmom_explicit_update[N],
        double ymom_explicit_update[N],
        double max_speed_array[N],

        double evolve_max_timestep,
        double g,
        double epsilon,
        double h0,
        double limiting_threshold,
        int optimise_dry_cells)
{
    int k; 


    for (k=0; k < N; k++)
    {
        double w_left, h_left, uh_left, vh_left, u_left;
        double w_right, h_right, uh_right, vh_right, u_right;
        double s_min, s_max, soundspeed_left, soundspeed_right;
        double denom, inverse_denominator, z;
        double temp;


        double q_left_0, q_left_1, q_left_2;
        double q_right_0, q_right_1, q_right_2;
        double q_left_rotated_0, q_left_rotated_1, q_left_rotated_2;
        double q_right_rotated_0, q_right_rotated_1, q_right_rotated_2;
        double flux_left_0, flux_left_1, flux_left_2;
        double flux_right_0, flux_right_1, flux_right_2;
        double edgeflux_0, edgeflux_1, edgeflux_2;


        double max_speed, max_speed_total;
        double z_left, z_right;
        double n1, n2;
        double length, inv_area;


        int i, m, ni;
        int ki, nm;

        max_speed_total = 0;
        stage_explicit_update[k] = 0;
        xmom_explicit_update[k] = 0;
        ymom_explicit_update[k] = 0;
        
        for (i=0; i< 3; i++)
        {
            ki= k*3 + i;
            
            q_left_0 = stage_edge_values[ki];
            q_left_1 = xmom_edge_values[ki];
            q_left_2 = ymom_edge_values[ki];
            z_left = bed_edge_values[ki];
            ni = neighbours[ki];
            if (ni<0) {
                m= -ni -1;
                
                q_right_0 = stage_boundary_values[m];
                q_right_1 = xmom_boundary_values[m];
                q_right_2 = ymom_boundary_values[m];
                z_right = z_left;
            } else {
                m = neighbour_edges[ki];
                nm = ni *3 + m;
                
                q_right_0 = stage_edge_values[nm];
                q_right_1 = xmom_edge_values[nm];
                q_right_2 = ymom_edge_values[nm];
                z_right = bed_edge_values[nm];
            }
            
            if(optimise_dry_cells) {
                if(fabs(q_left_0 - z_left) < epsilon &&
                   fabs(q_right_0 - z_right) < epsilon ) {
                    max_speed = 0.0;
                    continue;
                }
            }
            n1 = normals[2*ki];
            n2 = normals[2*ki + 1];
            /////////////////////////////////////////////////////////
            
            
            // Copy conserved quantities to protect from modification
            q_left_rotated_0 = q_left_0;
            q_right_rotated_0 = q_right_0;
            q_left_rotated_1 = q_left_1;
            q_right_rotated_1 = q_right_1;
            q_left_rotated_2 = q_left_2;
            q_right_rotated_2 = q_right_2;
            
            // Align x- and y-momentum with x-axis
            //_rotate(q_left_rotated, n1, n2);
            q_left_rotated_1 = n1*q_left_1 + n2*q_left_2;
            q_left_rotated_2 = -n2*q_left_1 + n1*q_left_2;
            
            //_rotate(q_right_rotated, n1, n2);
            q_right_rotated_1 = n1*q_right_1 + n2*q_right_2;
            q_right_rotated_2 = -n2*q_right_1 + n1*q_right_2;
            
            
            if (fabs(z_left - z_right) > 1.0e-10) {
                //report_python_error(AT, "Discontinuous Elevation");
                //return 0.0;
            }
            z = 0.5 * (z_left + z_right); // Average elevation values.
            
            // Compute speeds in x-direction
            w_left = q_left_rotated_0;
            h_left = w_left - z;
            uh_left = q_left_rotated_1;
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
            
            w_right = q_right_rotated_0;
            h_right = w_right - z;
            uh_right = q_right_rotated_1;
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
            vh_left = q_left_rotated_2;
            vh_right = q_right_rotated_2;
            
            soundspeed_left = sqrt(g * h_left);
            soundspeed_right = sqrt(g * h_right);
            
            s_max = u_left + soundspeed_left;
            if(s_max < u_right + soundspeed_right)
                s_max = u_right + soundspeed_right;
                

            if (s_max < 0.0) {
                s_max = 0.0;
            }
            
            s_min = u_left - soundspeed_left;
            if (s_min > u_right - soundspeed_right)
                s_min = u_right - soundspeed_right;

            if (s_min > 0.0) {
                s_min = 0.0;
            }
            
            // Flux formulas
            flux_left_0 = u_left*h_left;
            flux_left_1 = u_left * uh_left + 0.5 * g * h_left*h_left;
            flux_left_2 = u_left*vh_left;
            
            flux_right_0 = u_right*h_right;
            flux_right_1 = u_right * uh_right + 0.5 * g * h_right*h_right;
            flux_right_2 = u_right*vh_right;
            
            // Flux computation
            denom = s_max - s_min;
            if (denom < epsilon) { // FIXME (Ole): Try using h0 here
                edgeflux_0 = 0;
                edgeflux_1 = 0;
                edgeflux_2 = 0;

                max_speed = 0.0;
            }
            else {
                inverse_denominator = 1.0 / denom;
                //for (j = 0; j < 3; j++) {
                //    edgeflux[j] = s_max *flux_left[j] -s_min *flux_right[j];
                //    edgeflux[j] += s_max * s_min * (
                //        q_right_rotated[j] - q_left_rotated[j]);
                //    edgeflux[j] *= inverse_denominator;
                //}
                edgeflux_0 = s_max *flux_left_0 -s_min *flux_right_0;
                edgeflux_0 += s_max * s_min * (
                        q_right_rotated_0 - q_left_rotated_0);
                edgeflux_0 *= inverse_denominator;
                
                edgeflux_1 = s_max *flux_left_1 -s_min *flux_right_1;
                edgeflux_1 += s_max * s_min * (
                        q_right_rotated_1 - q_left_rotated_1);
                edgeflux_1 *= inverse_denominator;
                
                edgeflux_2 = s_max *flux_left_2 -s_min *flux_right_2;
                edgeflux_2 += s_max * s_min * (
                        q_right_rotated_2 - q_left_rotated_2);
                edgeflux_2 *= inverse_denominator;
                

                
                // Maximal wavespeed
                //max_speed = max(fabs(s_max), fabs(s_min));
                max_speed = fabs(s_max);
                if (max_speed < fabs(s_min))
                    max_speed = fabs(s_min);
                
                // Rotate back
                //_rotate(edgeflux, n1, -n2);
                temp = edgeflux_1;
                edgeflux_1 = n1* temp -n2*edgeflux_2;
                edgeflux_2 = n2*temp + n1*edgeflux_2;
            }
            
            //////////////////////////////////////////////////////

            length = edgelengths[ki];
            edgeflux_0 *= length;
            edgeflux_1 *= length;
            edgeflux_2 *= length;
            
            stage_explicit_update[k]-= edgeflux_0;
            xmom_explicit_update[k] -= edgeflux_1;
            ymom_explicit_update[k] -= edgeflux_2;
            
            if (tri_full_flag[k] == 1) {
                if (max_speed > epsilon) {
                    //timestep[k] = min(timestep[k], radii[k]/ max_speed);
                    if (timestep[k] > radii[k]/max_speed)
                        timestep[k] = radii[k]/ max_speed;
                    
                    if (ni>=0){
                        //timestep[k] =min(timestep[k], radii[ni]/max_speed);
                        if (timestep[k] > radii[ni]/max_speed)
                            timestep[k] = radii[ni]/ max_speed;

                    }
                }
            }
            if (ni < 0 ||  ni > k){
                max_speed_total = fmax(max_speed_total, max_speed);
            }
        }

        inv_area = 1.0 / areas[k];
        stage_explicit_update[k] *= inv_area;
        xmom_explicit_update[k] *= inv_area;
        ymom_explicit_update[k] *= inv_area;
        
        max_speed_array[k] = max_speed_total;
    }
}


#ifdef USING_MAIN
int main(int argc, char *argv[])
{
    int N, N2, optimise_dry_cells;
    int cnt_s =0, cnt_x =0, cnt_y =0, cnt_m = 0;

    double evolve_max_timestep, 
            g,
            epsilon,
            h0,
            limiting_threshold;

    double * sv, // stage_vertex_values
            * sb, // stage_boundary_values
            * se, // stage_edge_values, 
            * sc, // stage_centroid_values,
            * su, // stage_explicit_update

            * be, // bed_edge_values, 
            * bc, // bed_centroid_values, 
            
            * xb, // xmom_boundary_values
            * xe, // xmom_edge_values
            * xu, // xmom_explicit_update

            * yb, // ymom_boundary_values
            * ye, // ymom_edge_values
            * yu, // ymom_explicit_update

            * normals, 
            * el, // edgelengths;
            * radii,
            * a, // areas, 
            * vc; // vertex_coordinates, 

    int * tri; // tri_full_flag
    int * n, // neighbours
         * ne; // neighbour_edges

    double * timestep, * max_speed_array;

    // Testing group
    double * test_stage_explicit_update,
            * test_xmom_explicit_update,
            * test_ymom_explicit_update,
            * test_max_speed_array,
            * test_timestep;
    int i;

    char file_name[]="merimbula_cf.dat";

    freopen(file_name, "r", stdin);
    
    scanf("%d\n", &N);
    scanf("%d\n", &N2);
    scanf("%d\n", &optimise_dry_cells);
    scanf("%lf\n",&evolve_max_timestep);
    scanf("%lf\n",&g);
    scanf("%lf\n",&epsilon);
    scanf("%lf\n",&h0);
    scanf("%lf\n",&limiting_threshold);


    printf("%d %lf\n", N, g);

    
    sv = (double *)malloc( N*3*sizeof(double));
    assert(sv != NULL);
    sb = (double *)malloc( N2*sizeof(double));
    assert(sb != NULL);
    se = (double *)malloc( N*3*sizeof(double));
    assert(se != NULL);
    sc = (double *)malloc( N*sizeof(double));
    assert(sc != NULL);
    su = (double *)malloc( N*sizeof(double));
    assert(su != NULL);
    
    
    be = (double *)malloc( N*3*sizeof(double));
    assert(be != NULL);
    bc = (double *)malloc( N*sizeof(double));
    assert(bc != NULL);
    
    
    xb = (double *)malloc( N2*sizeof(double));
    assert(xb != NULL);
    xe = (double *)malloc( N*3*sizeof(double));
    assert(xe != NULL);
    xu = (double *)malloc( N*sizeof(double));
    assert(xu != NULL);

    yb = (double *)malloc( N2*sizeof(double));
    assert(yb != NULL);
    ye = (double *)malloc( N*3*sizeof(double));
    assert(ye != NULL);
    yu = (double *)malloc( N*sizeof(double));
    assert(yu != NULL);

    
    n  = (int *)malloc( N*3*sizeof(int));
    assert(n != NULL);
    ne  = (int *)malloc( N*3*sizeof(int));
    assert(ne != NULL);


    normals  = (double *)malloc( N*6*sizeof(double));
    assert(normals != NULL);
    el = (double *)malloc( N*3*sizeof(double));
    assert(el != NULL);
    radii  = (double *)malloc( N*sizeof(double));
    assert(radii != NULL);
    a  = (double *)malloc( N*sizeof(double));
    assert(a != NULL);
    
    tri  = (int *)malloc( N*sizeof(int));
    assert(tri != NULL);
    
    vc = (double *)malloc( N*6*sizeof(double));
    assert(vc != NULL);
    
    timestep = (double *)malloc(N*sizeof(double));
    assert(timestep != NULL);
    max_speed_array = (double *)malloc(N*sizeof(double));
    assert(max_speed_array != NULL);

    test_stage_explicit_update = (double *)malloc(N*sizeof(double));
    assert(test_stage_explicit_update != NULL);
    test_xmom_explicit_update = (double *)malloc(N*sizeof(double));
    assert(test_xmom_explicit_update != NULL);
    test_ymom_explicit_update = (double *)malloc(N*sizeof(double));
    assert(test_ymom_explicit_update != NULL);
    test_max_speed_array = (double *)malloc(N*sizeof(double));
    assert(test_max_speed_array != NULL);
    test_timestep = (double *)malloc(N*sizeof(double));
    assert(test_timestep != NULL);


    for(i=0; i < N; i++)
    {
        scanf("%lf %lf %lf\n", sv+i*3, sv+i*3 +1, sv+i*3 +2);
        scanf("%lf %lf %lf\n", se+i*3, se+i*3 +1, se+i*3 +2);
        scanf("%lf\n", sc+i);
        scanf("%lf\n", su+i);

        scanf("%lf %lf %lf\n",be+i*3, be+i*3 +1, be+i*3 +2);
        scanf("%lf\n", bc+i);

        scanf("%lf %lf %lf\n", xe+i*3, xe+i*3 +1, xe+i*3 +2);
        scanf("%lf\n", xu+i);
        scanf("%lf %lf %lf\n", ye+i*3, ye+i*3 +1, ye+i*3 +2);
        scanf("%lf\n", yu+i);

        scanf("%d %d %d\n", n+i*3, n+i*3 +1, n+i*3 +2);

        scanf("%d %d %d\n", ne+i*3, ne+i*3 +1, ne+i*3 +2);


        scanf("%lf %lf %lf %lf %lf %lf\n", normals+i*6, normals+i*6+1, normals+i*6+2,normals+i*6+3, normals+i*6+4, normals+i*6+5);

        scanf("%lf %lf %lf\n", el+i*3, el+i*3 +1, el+i*3 +2);

        scanf("%lf\n", radii+i);
        scanf("%lf\n", a+i);
        scanf("%d\n", tri+i);

        scanf("%lf %lf\n", vc+i*6, vc+i*6 +1);
        scanf("%lf %lf\n", vc+i*6+2, vc+i*6 +3);
        scanf("%lf %lf\n", vc+i*6+4, vc+i*6 +5);
    }

    for(i=0; i < N2; i++)
    {
        scanf("%lf\n", sb+i);
        scanf("%lf\n", xb+i);
        scanf("%lf\n", yb+i);
    }

    cnt_m = 0;
    for(i=0; i < N; i++)
    {
        if (max_speed_array[i] != test_max_speed_array[i])
            cnt_m ++;
    }
    printf("\n --> Test initial input %d\n", cnt_m);
    
    printf(" --> Enter Kernel\n");
    #pragma hmpp cf_central_single callsite
    compute_fluxes_central_structure_cuda_single(
            N, 
            N*3,
            N*6,
            N2,

            timestep,
            n,
            ne,
            normals,
            el,
            radii,
            a,
            tri,
            se,
            xe,
            ye,
            be,
            sb,
            xb,
            yb,
            su,
            xu,
            yu,
            max_speed_array,

            evolve_max_timestep, 
            g, 
            epsilon,
            h0,
            limiting_threshold,
            optimise_dry_cells);

    /*
    printf(" --> Enter C\n");
    compute_fluxes_central_structure_cuda_single( N, N*3, N*6, N2, timestep, n, ne, normals, el, radii, a, tri, se, xe, ye, be, sb, xb, yb, su, xu, yu, max_speed_array, evolve_max_timestep, g, epsilon, h0, limiting_threshold, optimise_dry_cells);

    */


    printf(" --> Enter Kernel\n");
    #pragma hmpp cf_central callsite
    compute_fluxes_central_structure_CUDA(
            N, 
            N*3,
            N*6,
            N2,

            timestep,
            n,
            ne,
            normals,
            el,
            radii,
            a,
            tri,
            se,
            xe,
            ye,
            be,
            sb,
            xb,
            yb,
            su,
            xu,
            yu,
            max_speed_array,

            evolve_max_timestep, 
            g, 
            epsilon,
            h0,
            limiting_threshold,
            optimise_dry_cells);



    printf(" --> Enter original C\n");
    compute_fluxes_central_structure_CUDA( 
            N, N*3, N*6, N2, 
            test_timestep, n, ne, normals, el, radii, a, tri, 
            se, xe, ye, be, sb, xb, yb, 
            test_stage_explicit_update, 
            test_xmom_explicit_update, 
            test_ymom_explicit_update, 
            test_max_speed_array, 
            evolve_max_timestep, g, epsilon, h0, limiting_threshold, 
            optimise_dry_cells);



    for (cnt_s=cnt_x=cnt_y=cnt_m=i=0; i < N; i++)
    {
        if ( fabs(su[i] - test_stage_explicit_update[i]) >= TOLERANCE)
        {
            cnt_s++;
            if (cnt_s <= 10)
                printf(" sta %d %lf %lf\n", i, su[i], test_stage_explicit_update[i]);
        }
        if ( fabs(xu[i] - test_xmom_explicit_update[i]) >= TOLERANCE)
        {
            cnt_x++;
            if (cnt_x <= 10)
                printf(" xmom %d %lf %lf\n", i, xu[i], test_xmom_explicit_update[i]);
        }
        if ( fabs(yu[i] - test_ymom_explicit_update[i]) >= TOLERANCE)
        {
            cnt_y++;
            if (cnt_y <= 10)
                printf(" ymom %d %lf %lf\n", i, yu[i], test_ymom_explicit_update[i]);
        }
        if ( fabs(max_speed_array[i] -test_max_speed_array[i]) >=TOLERANCE)
        {    
            cnt_m++;
            if (cnt_m <= 10)
                printf(" max %d %lf %lf\n", i, max_speed_array[i], test_max_speed_array[i]);
        }
    }
    printf("se:%d  xe:%d  ye:%d  max:%d errors found\n", cnt_s, cnt_x, cnt_y, cnt_m);
}
#endif

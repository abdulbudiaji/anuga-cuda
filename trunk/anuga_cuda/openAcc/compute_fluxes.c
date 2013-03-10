#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<assert.h>

double max(double a, double b)
{ return (a >= b) ? a : b; }

double min(double a, double b)
{ return (a <= b) ? a : b; }


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


int _flux_function_central(double *q_left, double *q_right,
        double z_left, double z_right,
        double n1, double n2,
        double epsilon,
        double h0,
        double limiting_threshold,
        double g,
        double *edgeflux, 
        double *max_speed) 
{
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


void compute_fluxes_central_structure_CUDA(
        int N,
        int N2,
        double evolve_max_timestep,
        double  g,
        double epsilon,
        double h0,
        double limiting_threshold,
        int optimise_dry_cells,

        double * restrict timestep,
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
        double * restrict stage_explicit_update,
        double * restrict xmom_explicit_update,
        double * restrict ymom_explicit_update, 
        double * restrict max_speed_array)
{


    double max_speed, max_speed_total=0 , length, inv_area, zl, zr;

    //double h0 = elements[DH0] * elements[DH0]; // This ensures a good balance when h approaches H0.

    //double limiting_threshold = 10 * elements[DH0]; // Avoid applying limiter below this

    int i, m, n, k;
    int ki, nm = 0, ki2; // Index shorthands

    double ql[3], qr[3], edgeflux[3];


    timestep[k] = evolve_max_timestep;
    
    memset((char*) stage_explicit_update, 0, N * sizeof (double));
    memset((char*) xmom_explicit_update, 0, N * sizeof (double));
    memset((char*) ymom_explicit_update, 0, N * sizeof (double));


//#ifdef UNSORTED_DOMAIN
//    int b[3]={0,1,2}, j;
//    spe_bubble_sort( b, neighbours+k*3, k);
//    for (j = 0; j < 3; j++) {
//        i = b[j];
//#else
//    for ( i = 0; i < 3; i++) {
//#endif
    for(k=0; k<N; k++)
    {
        max_speed_total = 0.0;
        // Loop through neighbours and compute edge flux for each
        for (i = 0; i < 3; i++) 
        {
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

            if (optimise_dry_cells){
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

            if (n < 0 ||  n > k){
                max_speed_total = max(max_speed_total, max_speed);
            }
        } // End edge i (and neighbour n)

            inv_area = 1.0 / areas[k];
            stage_explicit_update[k] *= inv_area;
            xmom_explicit_update[k] *= inv_area;
            ymom_explicit_update[k] *= inv_area;

            max_speed_array[k] =  max_speed_total;
    }
    
}




void compute_fluxes_central_structure_cuda_single(
        int N,
        int N2,
        double evolve_max_timestep,
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
        double * max_speed_array){
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
    #pragma acc kernels loop copyin(neighbours[0:3*N], neighbour_edges[0:3*N], normals[0:6*N],edgelengths[0:3*N], radii[0:N], areas[0:N], tri_full_flag[0:N], stage_edge_values[0:3*N], xmom_edge_values[0:3*N], ymom_edge_values[0:3*N], bed_edge_values[0:3*N], stage_boundary_values[0:N2], xmom_boundary_values[0:N2], ymom_boundary_values[0:N2]) copy(timestep[0:N], stage_explicit_update[0:N], xmom_explicit_update[0:N], ymom_explicit_update[0:N], max_speed_array[0:N]) 
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
            flux_left[0] = u_left*h_left;
            flux_left[1] = u_left * uh_left + 0.5 * g * h_left*h_left;
            flux_left[2] = u_left*vh_left;
            
            flux_right[0] = u_right*h_right;
            flux_right[1] = u_right * uh_right + 0.5 * g * h_right*h_right;
            flux_right[2] = u_right*vh_right;
            
            // Flux computation
            denom = s_max - s_min;
            if (denom < epsilon) { // FIXME (Ole): Try using h0 here
                edgeflux[0] = 0;
                edgeflux[1] = 0;
                edgeflux[2] = 0;

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
                //max_speed = max(fabs(s_max), fabs(s_min));
                max_speed = fabs(s_max);
                if (max_speed < fabs(s_min))
                    max_speed = fabs(s_min);
                
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
                    //timestep[k] = min(timestep[k], radii[k]/ max_speed);
                    if (timestep[k] > radii[k]/max_speed)
                        timestep[k] = radii[k]/ max_speed;
                    
                    if (n>=0){
                        //timestep[k] =min(timestep[k], radii[n]/max_speed);
                        if (timestep[k] > radii[n]/max_speed)
                            timestep[k] = radii[n]/ max_speed;
                    }
                }
            }
            //max_speed_total = max(max_speed_total, max_speed);
            if (max_speed_total < max_speed)
            max_speed_total = max_speed;
        }
        
        
        
        inv_area = 1.0 / areas[k];
        stage_explicit_update[k] *= inv_area;
        xmom_explicit_update[k] *= inv_area;
        ymom_explicit_update[k] *= inv_area;
        
        max_speed_array[k] = max_speed_total;
    }
}


int main(int argc, char *argv[])
{
    int N, N2, optimise_dry_cells;

    long * tri; // tri_full_flag
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

    long * n, // neighbours
         * ne; // neighbour_edges

    double * timestep, * max_speed_array;
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
    xe = (double *)malloc( N*sizeof(double));
    assert(xe != NULL);
    xu = (double *)malloc( N*3*sizeof(double));
    assert(xu != NULL);

    yb = (double *)malloc( N2*sizeof(double));
    assert(yb != NULL);
    ye = (double *)malloc( N*sizeof(double));
    assert(ye != NULL);
    yu = (double *)malloc( N*3*sizeof(double));
    assert(yu != NULL);

    
    n  = (long *)malloc( N*3*sizeof(long));
    assert(n != NULL);
    ne  = (long *)malloc( N*3*sizeof(long));
    assert(ne != NULL);
    normals  = (double *)malloc( N*6*sizeof(double));
    assert(normals != NULL);
    el = (double *)malloc( N*3*sizeof(double));
    assert(el != NULL);
    radii  = (double *)malloc( N*sizeof(double));
    assert(radii != NULL);
    a  = (double *)malloc( N*sizeof(double));
    assert(a != NULL);
    
    tri  = (long *)malloc( N*sizeof(long));
    assert(tri != NULL);
    
    vc = (double *)malloc( N*6*sizeof(double));
    assert(vc != NULL);
    
    timestep = (double *)malloc(N*sizeof(double));
    max_speed_array = (double *)malloc(N*sizeof(double));

    for(i=0; i < N; i++)
    {
        scanf("%lf %lf %lf\n", sv+i*3, sv+i*3 +1, sv+i*3 +2);
        scanf("%lf %lf %lf\n", se+i*3, se+i*3 +1, se+i*3 +2);
        scanf("%lf\n", sc+i);
        scanf("%lf\n", su+i);

        scanf("%lf %lf %lf\n",be+i*3, be+i*3 +1, be+i*3 +2);
        scanf("%lf\n", bc+i);

        scanf("%lf %lf %lf\n", xe+i*3, xe+i*3 +1, xe+i*3 +2);
        scanf("%lf\n", xe+i);
        scanf("%lf %lf %lf\n", ye+i*3, ye+i*3 +1, ye+i*3 +2);
        scanf("%lf\n", ye+i);

        scanf("%ld %ld %ld\n", n+i*3, n+i*3 +1, n+i*3 +2);

        scanf("%ld %ld %ld\n", ne+i*3, ne+i*3 +1, ne+i*3 +2);


        scanf("%lf %lf %lf %lf %lf %lf\n", normals+i*6, normals+i*6+1, normals+i*6+2,normals+i*6+3, normals+i*6+4, normals+i*6+5);

        scanf("%lf %lf %lf\n", el+i*3, el+i*3 +1, el+i*3 +2);

        scanf("%lf\n", radii+i);
        scanf("%lf\n", a+i);
        scanf("%ld\n", tri+i);

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

    compute_fluxes_central_structure_cuda_single(
            N, 
            N2,
            evolve_max_timestep, 
            g, 
            epsilon,
            h0,
            limiting_threshold,
            optimise_dry_cells,

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
            max_speed_array);
}


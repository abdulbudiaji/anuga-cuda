#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define AT __FILE__ ":" TOSTRING(__LINE__)
#define P_ERROR_BUFFER_SIZE 65


/*
 * this function can be used to detect kernel executing error, with cuda_memcopy in/out the 
 * pre-reserved error var
*/

__device__ void report_python_error(const char *location, const char *msg)
{
    //char buf[P_ERROR_BUFFER_SIZE];
    
    //snprintf(buf, P_ERROR_BUFFER_SIZE, "Error at %s: %s\n", location, msg);

    //PyErr_SetString(PyExc_RuntimeError, buf);
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


__device__ int _flux_function_central_wb(double *q_left, double *q_right,
        double z_left, double h_left, double h1_left, double h2_left,
        double z_right, double h_right, double h1_right, double h2_right,
        double n1, double n2,
        double epsilon,
        double h0,
        double limiting_threshold,
        double g,
        double *edgeflux,
        double *max_speed) 
{

    /*Compute fluxes between volumes for the shallow water wave equation
     * cast in terms of the 'stage', w = h+z using
     * the 'central scheme' as described in
     *
     * Kurganov, Noelle, Petrova. 'Semidiscrete Central-Upwind Schemes For
     * Hyperbolic Conservation Laws and Hamilton-Jacobi Equations'.
     * Siam J. Sci. Comput. Vol. 23, No. 3, pp. 707-740.
     *
     * The implemented formula is given in equation (3.15) on page 714.
     *
     * The pressure term is calculated using simpson's rule so that it
     * will well balance with the standard ghz_x gravity term
     */


    int i;
    //double hl, hr;
    double uh_left, vh_left, u_left, v_left;
    double uh_right, vh_right, u_right, v_right;
    double s_min, s_max, soundspeed_left, soundspeed_right;
    double denom, inverse_denominator;
    double p_left, p_right;

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


    uh_left = q_left_rotated[1];
    u_left = _compute_speed(&uh_left, &h_left,
            epsilon, h0, limiting_threshold);

    //w_right = q_right_rotated[0];
    //h_right = w_right - z;
    uh_right = q_right_rotated[1];
    u_right = _compute_speed(&uh_right, &h_right,
            epsilon, h0, limiting_threshold);

    // Momentum in y-direction
    vh_left = q_left_rotated[2];
    v_left = _compute_speed(&vh_left, &h_left,
            epsilon, h0, limiting_threshold);

    vh_right = q_right_rotated[2];
    v_right = _compute_speed(&vh_right, &h_right,
            epsilon, h0, limiting_threshold);

    // Limit y-momentum if necessary
    // Leaving this out, improves speed significantly (Ole 27/5/2009)
    // All validation tests pass, so do we really need it anymore?
    //_compute_speed(&vh_left, &h_left,
    //       epsilon, h0, limiting_threshold);
    //_compute_speed(&vh_right, &h_right,
    //      epsilon, h0, limiting_threshold);





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


    p_left = 0.5 * g / 6.0 * (h1_left * h1_left + 4.0 * h_left * h_left + h2_left * h2_left);
    p_right = 0.5 * g / 6.0 * (h1_right * h1_right + 4.0 * h_right * h_right + h2_right * h2_right);



    // Flux formulas
    flux_left[0] = u_left*h_left;
    flux_left[1] = u_left * uh_left + p_left;
    flux_left[2] = u_left*vh_left;


    flux_right[0] = u_right*h_right;
    flux_right[1] = u_right * uh_right + p_right;
    flux_right[2] = u_right*vh_right;

    
    // Flux computation
    denom = s_max - s_min;
    if (denom < epsilon) { // FIXME (Ole): Try using h0 here
        memset(edgeflux, 0, 3 * sizeof (double));
        *max_speed = 0.0;
    } else {
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

__device__ int _flux_function_central_wb_3(
        //struct edge *E_left,
        //struct edge *E_right,
        int * E_left_cell_id,
        int * E_left_edge_id,
        double * E_left_w,
        double * E_left_h,
        double * E_left_z,
        double * E_left_uh,
        double * E_left_vh,
        double * E_left_u,
        double * E_left_v,
        double * E_left_w1,
        double * E_left_w2,
        double * E_left_h1,
        double * E_left_h2,
        double * E_left_z1,
        double * E_left_z2,
        double * E_left_uh1,
        double * E_left_uh2,
        double * E_left_vh1,
        double * E_left_vh2,
        double * E_left_u1,
        double * E_left_u2,
        double * E_left_v1,
        double * E_left_v2,

        int * E_right_cell_id,
        int * E_right_edge_id,
        double * E_right_w,
        double * E_right_h,
        double * E_right_z,
        double * E_right_uh,
        double * E_right_vh,
        double * E_right_u,
        double * E_right_v,
        double * E_right_w1,
        double * E_right_w2,
        double * E_right_h1,
        double * E_right_h2,
        double * E_right_z1,
        double * E_right_z2,
        double * E_right_uh1,
        double * E_right_uh2,
        double * E_right_vh1,
        double * E_right_vh2,
        double * E_right_u1,
        double * E_right_u2,
        double * E_right_v1,
        double * E_right_v2,
        double n1,
        double n2,
        double epsilon,
        double h0,
        double limiting_threshold,
        double g,
        double *edgeflux,
        double *max_speed) {

    /*Compute fluxes between volumes for the shallow water wave equation
      cast in terms of the 'stage', w = h+z using
      the 'central scheme' as described in

      Kurganov, Noelle, Petrova. 'Semidiscrete Central-Upwind Schemes For
      Hyperbolic Conservation Laws and Hamilton-Jacobi Equations'.
      Siam J. Sci. Comput. Vol. 23, No. 3, pp. 707-740.

      The implemented formula is given in equation (3.15) on page 714
     */

    int i;
    double h_left, h_right;
    double uh_left, vh_left, u_left;
    double uh_right, vh_right, u_right;
    double h1_left, h2_left, h1_right, h2_right;
    double uh1_left, uh2_left, uh1_right, uh2_right;
    double vh1_left, vh2_left, vh1_right, vh2_right;
    double u1_left, u2_left, u1_right, u2_right;
    double p_left, p_right;

    double s_min, s_max, soundspeed_left, soundspeed_right;
    double denom, inverse_denominator;

    double q_left[3], q_right[3], flux_right[3], flux_left[3];


    // Align x- and y-momentum with x-axis
    //_rotate_edge(E_left, n1, n2);
    (*E_left_uh) = n1 * (*E_left_uh) + n2 * (*E_left_vh);
    (*E_left_vh) = -n2 * (*E_left_uh) + n1 * (*E_left_vh);

    (*E_left_uh1) = n1 * (*E_left_uh1) + n2 * (*E_left_vh1);
    (*E_left_vh1) = -n2 * (*E_left_uh1) + n1 * (*E_left_vh1);

    (*E_left_uh2) = n1 * (*E_left_uh2) + n2 * (*E_left_vh2);
    (*E_left_vh2) = -n2 * (*E_left_uh2) + n1 * (*E_left_vh2);
    

    //_rotate_edge(E_right, n1, n2);
    (*E_right_uh) = n1 * (*E_right_uh) + n2 * (*E_right_vh);
    (*E_right_vh) = -n2 * (*E_right_uh) + n1 * (*E_right_vh);

    (*E_right_uh1) = n1 * (*E_right_uh1) + n2 * (*E_right_vh1);
    (*E_right_vh1) = -n2 * (*E_right_uh1) + n1 * (*E_right_vh1);

    (*E_right_uh2) = n1 * (*E_right_uh2) + n2 * (*E_right_vh2);
    (*E_right_vh2) = -n2 * (*E_right_uh2) + n1 * (*E_right_vh2);



    q_left[0] = *E_left_w;
    q_left[1] = *E_left_uh;
    q_left[2] = *E_left_vh;


    q_right[0] = *E_right_w;
    q_right[1] = *E_right_uh;
    q_right[2] = *E_right_vh;

    //printf("========== wb_3 ==============\n");
    //printf("E_left %i %i \n",E_left->cell_id, E_left->edge_id);
    //printf("E_right %i %i \n",E_right->cell_id, E_right->edge_id);
    
    //printf("q_left %f %f %f \n", q_left[0], q_left[1], q_left[2]);
    //printf("q_right %f %f %f \n", q_right[0], q_right[1], q_right[2]);

    //z = 0.5*(E_left->z + E_right->z); // Average elevation values.
    // Even though this will nominally allow
    // for discontinuities in the elevation data,
    // there is currently no numerical support for
    // this so results may be strange near
    // jumps in the bed.

    // Compute speeds in x-direction
    uh_left = *E_left_uh;
    h_left = *E_left_h;
    u_left = _compute_speed(&uh_left, &h_left,
            epsilon, h0, limiting_threshold);

    h_right = *E_right_h;
    uh_right = *E_right_uh;
    u_right = _compute_speed(&uh_right, &h_right,
            epsilon, h0, limiting_threshold);


    //printf("uh_left, u_left, h_left %g %g %g \n", uh_left, u_left, h_left);
    //printf("uh_right, u_right, h_right %g %g %g \n", uh_right, u_right, h_right);

    // Momentum in y-direction
    vh_left = *E_left_vh;
    vh_right = *E_right_vh;

    // Maximal and minimal wave speeds
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


    // Check if boundary edge or not
    if (*E_right_cell_id >= 0) {
        // Interior Edge
        // First vertex left edge data
        h1_left = *E_left_h1;
        uh1_left = *E_left_uh1;
        vh1_left = *E_left_vh1;
        u1_left = _compute_speed(&uh1_left, &h1_left, epsilon, h0, limiting_threshold);

        // Second vertex left edge data
        h2_left = *E_left_h2;
        uh2_left = *E_left_uh2;
        vh2_left = *E_left_vh2;
        u2_left = _compute_speed(&uh2_left, &h2_left, epsilon, h0, limiting_threshold);

        // First vertex right edge data (needs interchange)
        h1_right = *E_right_h2;
        uh1_right = *E_right_uh2;
        vh1_right = *E_right_vh2;
        u1_right = _compute_speed(&uh1_right, &h1_right, epsilon, h0, limiting_threshold);

        // Second vertex right edge data (needs interchange)
        h2_right = *E_right_h1;
        uh2_right = *E_right_uh1;
        vh2_right = *E_right_vh1;
        u2_right = _compute_speed(&uh2_right, &h2_right, epsilon, h0, limiting_threshold);
    } else {
        // Boundary Edge
        // First vertex left edge data
        h1_left  = h_left;
        uh1_left = uh_left;
        vh1_left = vh_left;
        u1_left  = u_left;

        // Second vertex left edge data
        h2_left  = h_left;
        uh2_left = uh_left;
        vh2_left = vh_left;
        u2_left  = u_left;

        // First vertex right edge data (needs interchange)
        h1_right  = h_right;
        uh1_right = uh_right;
        vh1_right = vh_right;
        u1_right  = u_right;

        // Second vertex right edge data (needs interchange)
        h2_right  = h_right;
        uh2_right = uh_right;
        vh2_right = vh_right;
        u2_right  = u_right;
    }

    //printf("h1_left, h2_left %g %g \n", h1_left, h2_left);
    //printf("h1_right, h2_right %g %g \n", h1_right, h2_right);

    p_left  = 0.5*g*h_left*h_left;
    p_right = 0.5*g*h_right*h_right;


    //p_left = 0.5 * g / 6.0 * (h1_left * h1_left + 4.0 * h_left * h_left + h2_left * h2_left);
    //p_right = 0.5 * g / 6.0 * (h1_right * h1_right + 4.0 * h_right * h_right + h2_right * h2_right);


    //printf("p_left, p_right %g %g \n", p_left, p_right);

    // Flux formulas
    //flux_left[0] = u_left*h_left;
    //flux_left[1] = u_left * uh_left + p_left;
    //flux_left[2] = u_left*vh_left;


    flux_left[0] = (u1_left*h1_left  + 4.0*u_left*h_left  + u2_left*h2_left)/6.0;
    flux_left[1] = (u1_left*uh1_left + 4.0*u_left*uh_left + u2_left*uh2_left)/6.0 + p_left;
    flux_left[2] = (u1_left*vh1_left + 4.0*u_left*vh_left + u2_left*vh2_left)/6.0;


    //printf("flux_left %g %g %g \n",flux_left[0],flux_left[1],flux_left[2]);

    //flux_right[0] = u_right*h_right;
    //flux_right[1] = u_right*uh_right + p_right;
    //flux_right[2] = u_right*vh_right;


    flux_right[0] = (u1_right*h1_right  + 4.0*u_right*h_right  + u2_right*h2_right)/6.0;
    flux_right[1] = (u1_right*uh1_right + 4.0*u_right*uh_right + u2_right*uh2_right)/6.0 + p_right;
    flux_right[2] = (u1_right*vh1_right + 4.0*u_right*vh_right + u2_right*vh2_right)/6.0;

    //printf("flux_right %g %g %g \n",flux_right[0],flux_right[1],flux_right[2]);

    // Flux computation
    denom = s_max - s_min;
    if (denom < epsilon) { // FIXME (Ole): Try using h0 here
        memset(edgeflux, 0, 3 * sizeof (double));
        *max_speed = 0.0;
    } else {
        inverse_denominator = 1.0 / denom;
        for (i = 0; i < 3; i++) {
            edgeflux[i] = s_max * flux_left[i] - s_min * flux_right[i];
            edgeflux[i] += s_max * s_min * (q_right[i] - q_left[i]);
            edgeflux[i] *= inverse_denominator;
        }

        // Maximal wavespeed
        *max_speed = max(fabs(s_max), fabs(s_min));

        // Rotate back
        _rotate(edgeflux, n1, -n2);
    }


    //printf("edgeflux %g %g %g \n",edgeflux[0],edgeflux[1],edgeflux[2]);
    
    return 0;
}


#define Dtimestep 0
#define Depsilon 1
#define DH0 2
#define Dg 3
#define Doptimise_dry_cells 4
#define Dcall 5

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
            report_python_error(AT, "Discontinuous Elevation");
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
}

#define Devolve_max_timestep 6
/*
 * elements[0] = _timestep
 * elements[1] = _epsilon
 * elements[2] = _H0
 * elements[3] = _g
 * elements[4] = _optimise_dry_cells
 * elements[5] = _compute_fluxes_central_structure_call
 * elements[6] = _evolve_max_timestep
 */

static long compute_fluxes_central_structure_call = 1;

__global__ void _compute_fluxes_central_structure(
        double * elements,
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
        double * max_speed_array)
{
    int k ;//= threadIdx.x + threadIdx.y + blockIdx.x * blockDim.x + blockIdx.y *blockDim.y;

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

        if (already_computed_flux[ki] == elements[Dcall]) {
            // We've already computed the flux across this edge
            //continue;
            return;
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
                //continue;
                return;
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
                //  timestep = min(timestep, (radii[k]+radii[n])/max_speed);
                //} else {
                //  timestep = min(timestep, radii[k]/max_speed);
                // }
            }
        }

    } // End edge i (and neighbour n)

}
    // Normalise triangle k by area and store for when all conserved
    // quantities get updated
    inv_area = 1.0 / areas[k];
    stage_explicit_update[k] *= inv_area;
    xmom_explicit_update[k] *= inv_area;
    ymom_explicit_update[k] *= inv_area;


    // Keep track of maximal speeds
    max_speed_array[k] = max_speed;

    //return timestep;
}

__global__ void compute_fluxes_central_structure_cuda(
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
        long * already_computed_flux,
        double * max_speed_array)
{
    int k = threadIdx.x + threadIdx.y + blockIdx.x * blockDim.x + blockIdx.y *blockDim.y;


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

    double temp=0;


    double w_left, h_left, uh_left, vh_left, u_left;
    double w_right, h_right, uh_right, vh_right, u_right;
    double s_min, s_max, soundspeed_left, soundspeed_right;
    double denom, inverse_denominator, z;
    
    double q_left_rotated[3], q_right_rotated[3], flux_right[3], flux_left[3];

     // Loop through neighbours and compute edge flux for each
    for (i = 0; i < 3; i++) {
        ki = k * 3 + i; // Linear index to edge i of triangle k

        n = neighbours[ki];

        // Get left hand side values from triangle k, edge i
        ql[0] = stage_edge_values[ki];
        ql[1] = xmom_edge_values[ki];
        ql[2] = ymom_edge_values[ki];
        zl = bed_edge_values[ki];

        // Get right hand side values either from neighbouring triangle
        // or from boundary array (Quantities at neighbour on nearest face).
        
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

                max_speed = 0.0;
                continue;
            }
        }



        // Outward pointing normal vector (domain.normals[k, 2*i:2*i+2])
        ki2 = 2 * ki; //k*6 + i*2

        /*
        // Edge flux computation (triangle k, edge i)
        _flux_function_central(ql, qr, zl, zr,
               normals[ki2], normals[ki2 + 1],
                elements[Depsilon], h0, limiting_threshold, elements[Dg],
                edgeflux, &max_speed);
        */


        q_left_rotated[0] = ql[0];
        q_right_rotated[0] = qr[0];
    q_left_rotated[1] = ql[1];
    q_right_rotated[1] = qr[1];
    q_left_rotated[2] = ql[2];
    q_right_rotated[2] = qr[2];
        //_rotate(double *q, double n1, double n2)

        
        temp = q_left_rotated[1];
        q_left_rotated[1] = normals[ki2] * q_left_rotated[1] + normals[ki2+1]*q_left_rotated[2];
        q_left_rotated[2] = -normals[ki2+1] * temp + normals[ki2]*q_left_rotated[2];

        temp = q_right_rotated[1];
        q_right_rotated[1] = normals[ki2] * q_right_rotated[1] + normals[ki2+1]*q_right_rotated[2];
        q_right_rotated[2] = -normals[ki2+1] * temp + normals[ki2]*q_right_rotated[2];
        
        z = 0.5 * (zl + zr);

        
    h_left =  q_left_rotated[0] - z;
    uh_left = q_left_rotated[1];
    //u_left = _compute_speed(&uh_left, &h_left,
    //        epsilon, h0, limiting_threshold);
    if( h_left < limiting_threshold)
    {
        if(h_left < elements[Depsilon] )
        {
            h_left = 0.0;
            u_left = 0.0;
        }
        else 
            u_left = uh_left/ (h_left + h0/h_left);
        uh_left = u_left * h_left;
    }
    else
        u_left = uh_left / h_left;



    
    h_right = q_right_rotated[0] - z;
    uh_right = q_right_rotated[1];
    //u_right = _compute_speed(&uh_right, &h_right,
    //        epsilon, h0, limiting_threshold);
    if( h_right < limiting_threshold)
    {
        if(h_right < elements[Depsilon] )
        {
            h_right = 0.0;
            u_right = 0.0;
        }
        else 
            u_right = uh_right/ (h_right + h0/h_right);
        uh_right = u_right * h_right;
    }
    else
        u_right = uh_right / h_right;


    
    vh_left = q_left_rotated[2];
    vh_right = q_right_rotated[2];
    
    //_compute_speed(&vh_left, &h_left,
    //        epsilon, h0, limiting_threshold);
    if( h_left < limiting_threshold)
    {
        if(h_left < elements[Depsilon] )
        {
            h_left = 0.0;
            temp = 0.0;
        }
        else 
            temp = vh_left/ (h_left + h0/h_left);
        vh_left = u_left * h_left;
    }

    //_compute_speed(&vh_right, &h_right,
    //        epsilon, h0, limiting_threshold);
    if( h_right < limiting_threshold)
    {
        if(h_right < elements[Depsilon] )
        {
            h_right = 0.0;
            temp = 0.0;
        }
        else 
            temp = vh_right/ (h_right + h0/h_right);
        vh_right = temp * h_right;
    }


    soundspeed_left = sqrt(elements[Dg] * h_left);
    soundspeed_right = sqrt(elements[Dg] * h_right);

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
    flux_left[1] = u_left * uh_left + 0.5 * elements[Dg] * h_left*h_left;
    flux_left[2] = u_left*vh_left;

    flux_right[0] = u_right*h_right;
    flux_right[1] = u_right * uh_right + 0.5 * elements[Dg] * h_right*h_right;
    flux_right[2] = u_right*vh_right;

    // Flux computation
    denom = s_max - s_min;
    if (denom < elements[Depsilon]) { // FIXME (Ole): Try using h0 here
        memset(edgeflux, 0, 3 * sizeof (double));
        max_speed = 0.0;
    }
    else {
        inverse_denominator = 1.0 / denom;
        for (i = 0; i < 3; i++) {
            edgeflux[i] = s_max * flux_left[i] - s_min * flux_right[i];
            edgeflux[i] += s_max * s_min * (q_right_rotated[i] - q_left_rotated[i]);
            edgeflux[i] *= inverse_denominator;
        }

        // Maximal wavespeed
        max_speed = max(fabs(s_max), fabs(s_min));

        // Rotate back
        //_rotate(edgeflux, n1, -n2);
        temp = edgeflux[1];
        edgeflux[1] = normals[ki2] * edgeflux[1] + normals[ki2+1]*edgeflux[2];
        edgeflux[2] = -normals[ki2+1] * temp + normals[ki2]*edgeflux[2];

    }



        // Multiply edgeflux by edgelength
        length = edgelengths[ki];
        edgeflux[0] *= length;
        edgeflux[1] *= length;
        edgeflux[2] *= length;


        // Update triangle k with flux from edge i
        stage_explicit_update[k] -= edgeflux[0];
        xmom_explicit_update[k] -= edgeflux[1];
        ymom_explicit_update[k] -= edgeflux[2];




        if (tri_full_flag[k] == 1) {
            if (max_speed > elements[Depsilon]) {
                timestep[k] = radii[k]/max_speed;
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

    max_speed_array[k] = max_speed;

}




    /*****************************************/
    /* Rearrange the domain variable order   */
    /* so as to achieve memory corelasing    */
    /*****************************************/

#define COMPUTE_FLUXES_CENTRAL_STRUCTURE_SHARED_MEMORY_LENGTH

__global__ void compute_fluxes_central_structure_MeCo(
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
        long * already_computed_flux,
        double * max_speed_array)
{
    const int k = threadIdx.x + (blockIdx.x )*blockDim.x;


    double max_speed, length, inv_area, zl, zr;

    double h0 = elements[DH0] * elements[DH0]; // This ensures a good balance when h approaches H0.

    double limiting_threshold = 10 * elements[DH0]; // Avoid applying limiter below this

    int i, m, n;
    int ki, nm , ki2; // Index shorthands

//    double ql[3], qr[3], edgeflux[3]; // Work array for summing up fluxes


    __shared__ double shared_var[COMPUTE_FLUXES_CENTRAL_STRUCTURE_SHARED_MEMORY_LENGTH];


        for (i = 0; i < 3; i++) {


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
    



/*
 * elements[0] = _timestep
 * elements[1] = _epsilon
 * elements[2] = _H0
 * elements[3] = _g
 * elements[4] = _optimise_dry_cells
 * elements[5] = _compute_fluxes_central_wb_call
 * elements[6] = _evolve_max_timestep
 */

static long compute_fluxes_central_wb_call = 1;


__global__ void _compute_fluxes_central_wb(
        double * elements,
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
        double * stage_vertex_values,
        double * bed_vertex_values,
        double * ymom_boundary_values,
        double * stage_explicit_update,
        double * xmom_explicit_update,
        double * ymom_explicit_update, 
        long * already_computed_flux,
        double * max_speed_array)

{
    int k = threadIdx.x + threadIdx.y + blockIdx.x * blockDim.x + blockIdx.y *blockDim.y;

    // Local variables
    double max_speed, length, inv_area;
    double timestep = 1.0e30;
    double h0 = elements[DH0] * elements[DH0]; // This ensures a good balance when h approaches H0.

    double limiting_threshold = 10 * elements[DH0]; // Avoid applying limiter below this
    // threshold for performance reasons.
    // See ANUGA manual under flux limiting
    int i, m, n, nn;
    int k3, k3i, k3i1, k3i2, k2i; // Index short hands

    int n3m = 0, n3m1, n3m2; // Index short hand for neightbours

    double hl, hl1, hl2;
    double hr, hr1, hr2;
    double zl, zl1, zl2;
    double zr;
    double wr;

    // Workspace (making them static actually made function slightly slower (Ole))
    double ql[3], qr[3], edgeflux[3]; // Work array for summing up fluxes

    //static long call = 1; // Static local variable flagging already computed flux

    // Start computation
    //call++; // Flag 'id' of flux calculation for this timestep


    timestep = elements[Devolve_max_timestep];

    // Set explicit_update to zero for all conserved_quantities.
    // This assumes compute_fluxes called before forcing terms
    //memset((char*) D->stage_explicit_update, 0, D->number_of_elements * sizeof (double));
    //memset((char*) D->xmom_explicit_update, 0, D->number_of_elements * sizeof (double));
    //memset((char*) D->ymom_explicit_update, 0, D->number_of_elements * sizeof (double));

    // For all triangles
    //for (k = 0; k < D->number_of_elements; k++) {
    // Loop through neighbours and compute edge flux for each
    for (i = 0; i < 3; i++) {
        k3 = 3 * k;
        k3i = k3 + i; // Linear index to edge i of triangle k
        k3i1 = k3 + (i + 1) % 3;
        k3i2 = k3 + (i + 2) % 3;

        if (already_computed_flux[k3i] == elements[Dcall]) {
            // We've already computed the flux across this edge
            continue;
        }


        // Get the inside values at the vertices from the triangle k, edge i
        ql[0] = stage_edge_values[k3i];
        ql[1] = xmom_edge_values[k3i];
        ql[2] = ymom_edge_values[k3i];

        zl = bed_edge_values[k3i];

        zl1 = bed_vertex_values[k3i1];
        zl2 = bed_vertex_values[k3i2];


        hl = stage_edge_values[k3i] - zl;

        hl1 = stage_vertex_values[k3i1] - zl1;
        hl2 = stage_vertex_values[k3i2] - zl2;


        // Get right hand side values either from neighbouring triangle
        // or from boundary array (Quantities at neighbour on nearest face).
        n = neighbours[k3i];
        if (n < 0) {
            // Neighbour is a boundary condition
            nn = -n - 1; // Convert negative flag to boundary index
            m = 0;

            qr[0] = stage_boundary_values[nn];
            qr[1] = xmom_boundary_values[nn];
            qr[2] = ymom_boundary_values[nn];

            zr = zl; // Extend bed elevation to boundary and assume flat stage

            wr = stage_boundary_values[nn];

            hr = wr - zr;
            hr1 = wr - zl2;
            hr2 = wr - zl1;


        } else {
            // Neighbour is a real triangle
            m = neighbour_edges[k3i];
            n3m = 3 * n + m; // Linear index (triangle n, edge m)
            n3m1 = 3 * n + (m + 1) % 3;
            n3m2 = 3 * n + (m + 2) % 3;

            qr[0] = stage_edge_values[n3m];
            qr[1] = xmom_edge_values[n3m];
            qr[2] = ymom_edge_values[n3m];

            zr = bed_edge_values[n3m];
            hr = stage_edge_values[n3m] - zr;

            hr1 = stage_vertex_values[n3m2] - bed_vertex_values[n3m2];
            hr2 = stage_vertex_values[n3m1] - bed_vertex_values[n3m1];

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

                already_computed_flux[k3i] = elements[Dcall]; // #k Done
                if (n >= 0) {
                    already_computed_flux[n3m] = elements[Dcall]; // #n Done
                }

                max_speed = 0.0;
                continue;
            }
        }


        if (fabs(zl - zr) > 1.0e-10) {
            report_python_error(AT, "Discontinuous Elevation");
            //return 0.0;
            return;
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
                normals[k2i], normals[k2i + 1],
                elements[Depsilon], h0, limiting_threshold, elements[Dg],
                edgeflux, &max_speed);



        // Multiply edgeflux by edgelength
        length = edgelengths[k3i];
        edgeflux[0] *= length;
        edgeflux[1] *= length;
        edgeflux[2] *= length;


        // Update triangle k with flux from edge i
        stage_explicit_update[k] -= edgeflux[0];
        xmom_explicit_update[k] -= edgeflux[1];
        ymom_explicit_update[k] -= edgeflux[2];

        already_computed_flux[k3i] = elements[Dcall]; // #k Done


        // Update neighbour n with same flux but reversed sign
        if (n >= 0) {
            stage_explicit_update[n] += edgeflux[0];
            xmom_explicit_update[n] += edgeflux[1];
            ymom_explicit_update[n] += edgeflux[2];

            already_computed_flux[n3m] = elements[Dcall]; // #n Done
        }

        // Update timestep based on edge i and possibly neighbour n
        if (tri_full_flag[k] == 1) {
            if (max_speed > elements[Depsilon]) {
                // Apply CFL condition for triangles joining this edge (triangle k and triangle n)

                // CFL for triangle k
                timestep = min(timestep, radii[k] / max_speed);

                if (n >= 0) {
                    // Apply CFL condition for neigbour n (which is on the ith edge of triangle k)
                    timestep = min(timestep, radii[n] / max_speed);
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
    inv_area = 1.0 / areas[k];
    stage_explicit_update[k] *= inv_area;
    xmom_explicit_update[k] *= inv_area;
    ymom_explicit_update[k] *= inv_area;


    // Keep track of maximal speeds
    max_speed_array[k] = max_speed;

    //} // End triangle k


    //return timestep;
}


/*
 * elements[0] = _timestep
 * elements[1] = _epsilon
 * elements[2] = _H0
 * elements[3] = _g
 * elements[4] = _optimise_dry_cells
 * elements[5] = _compute_fluxes_central_wb_3_call
 * elements[6] = _evolve_max_timestep
 */

static long compute_fluxes_central_wb_3_call = 1;


__global__ void compute_fluxes_central_wb_3(
        double* elements,
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
        double * stage_vertex_values,
        double * xmom_vertex_values,
        double * ymom_vertex_values,
        double * bed_vertex_values,
        double * ymom_boundary_values,
        double * stage_explicit_update,
        double * xmom_explicit_update,
        double * ymom_explicit_update, 
        long * already_computed_flux,
        double * max_speed_array) 
{
    int k = threadIdx.x + threadIdx.y + blockIdx.x * blockDim.x + blockIdx.y *blockDim.y;

    // Local variables
    double max_speed, length, inv_area;
    double timestep = 1.0e30;
    double h0 = elements[DH0] * elements[DH0]; // This ensures a good balance when h approaches H0.

    double limiting_threshold = 10 * elements[DH0]; // Avoid applying limiter below this
    // threshold for performance reasons.
    // See ANUGA manual under flux limiting
    int i, m, n;
    int k3, k3i, k2i; // Index short hands
    int n3m = 0; // Index short hand for neightbours

    //struct edge E_left, E_right;

    int E_left_cell_id, E_left_edge_id;
    int E_right_cell_id, E_right_edge_id;

    double E_left_w, E_left_h, E_left_z, E_left_uh, E_left_vh, E_left_u, E_left_v;
    double E_left_w1, E_left_h1, E_left_z1, E_left_uh1, E_left_vh1, E_left_u1, E_left_v1;
    double E_left_w2, E_left_h2, E_left_z2, E_left_uh2, E_left_vh2, E_left_u2, E_left_v2;

    double E_right_w, E_right_h, E_right_z, E_right_uh, E_right_vh, E_right_u, E_right_v;
    double E_right_w1, E_right_h1, E_right_z1, E_right_uh1, E_right_vh1, E_right_u1, E_right_v1;
    double E_right_w2, E_right_h2, E_right_z2, E_right_uh2, E_right_vh2, E_right_u2, E_right_v2;

    double edgeflux[3]; // Work array for summing up fluxes

    //static long call = 1; // Static local variable flagging already computed flux

    // Start computation
    // call++; // Flag 'id' of flux calculation for this timestep


    timestep = elements[Devolve_max_timestep];

    // Set explicit_update to zero for all conserved_quantities.
    // This assumes compute_fluxes called before forcing terms

    //    memset((char*) D->stage_explicit_update, 0, D->number_of_elements * sizeof (double));
    //    memset((char*) D->xmom_explicit_update,  0, D->number_of_elements * sizeof (double));
    //    memset((char*) D->ymom_explicit_update,  0, D->number_of_elements * sizeof (double));

    // For all triangles
    //    for (k = 0; k < D->number_of_elements; k++) {
    // Loop through neighbours and compute edge flux for each
    for (i = 0; i < 3; i++) {
        k3 = 3 * k;
        k3i = k3 + i; // Linear index to edge i of triangle k

        if (already_computed_flux[k3i] == elements[Dcall]) {
            // We've already computed the flux across this edge
            //continue;
            return;
        }

        // Get the "left" values at the edge vertices and midpoint
        // from the triangle k, edge i
        // get_edge_data(&E_left, D, k, i);
        E_left_cell_id = k;
        E_left_edge_id = i;

        E_left_w = stage_edge_values[k3 + i];
        E_left_z = bed_edge_values[k3 + i];
        E_left_h = E_left_w - E_left_z;
        E_left_uh = xmom_edge_values[k3 + i];
        E_left_vh = ymom_edge_values[k3 + i];

        E_left_w1 = stage_vertex_values[k3 + (i+1)%3];
        E_left_z1 = bed_vertex_values[k3 + (i+1)%3];
        E_left_h1 = E_left_w1 - E_left_z1;
        E_left_uh1 = xmom_vertex_values[k3 + (i+1)%3];
        E_left_vh1 = ymom_vertex_values[k3 + (i+1)%3];

        E_left_w2 = stage_vertex_values[k3 + (i+2)%3];
        E_left_z2 = bed_vertex_values[k3 + (i+2)%3];
        E_left_h2 = E_left_w2 - E_left_z2;
        E_left_uh2 = xmom_vertex_values[k3 + (i+2)%3];
        E_left_vh2 = ymom_vertex_values[k3 + (i+2)%3];


        // Get right hand side values either from neighbouring triangle
        // or from boundary array (Quantities at neighbour on nearest face).
        n = neighbours[k3i];
        if (n < 0) {
            // Neighbour is a boundary condition
            m = -n - 1; // Convert negative flag to boundary index

            E_right_cell_id = n;
            E_right_edge_id = 0;

            // Midpoint Values provided by the boundary conditions
            E_right_w = stage_boundary_values[m];
            E_right_uh = xmom_boundary_values[m];
            E_right_vh = ymom_boundary_values[m];

            // Set bed and height as equal to neighbour
            E_right_z = E_left_z;
            E_right_h = E_right_w - E_right_z;

            // vertex values same as midpoint values
            E_right_w1  = E_right_w;
            E_right_h1  = E_right_h;

            E_right_uh1 = E_left_uh2;
            E_right_vh1 = E_left_vh2;
            E_right_z1  = E_left_z;

            E_right_w2  = E_right_w;
            E_right_h2  = E_right_h;

            E_right_uh2 = E_right_uh;
            E_right_vh2 = E_right_vh;
            E_right_z2  = E_left_z;

        } else {
            // Neighbour is a real triangle
            m = neighbour_edges[k3i];
            n3m  = 3*n+m;

            //get_edge_data(&E_right_ D, n, m);

            E_right_cell_id = k;
            E_right_edge_id = i;

            E_right_w = stage_edge_values[n*3 + m];
            E_right_z = bed_edge_values[n*3 + m];
            E_right_h = E_right_w - E_right_z;
            E_right_uh = xmom_edge_values[n*3 + m];
            E_right_vh = ymom_edge_values[n*3 + m];

            E_right_w1 = stage_vertex_values[n*3 + (m+1)%3];
            E_right_z1 = bed_vertex_values[n*3 + (m+1)%3];
            E_right_h1 = E_right_w1 - E_right_z1;
            E_right_uh1 = xmom_vertex_values[n*3 + (m+1)%3];
            E_right_vh1 = ymom_vertex_values[n*3 + (m+1)%3];

            E_right_w2 = stage_vertex_values[n*3 + (m+2)%3];
            E_right_z2 = bed_vertex_values[n*3 + (m+2)%3];
            E_right_h2 = E_right_w2 - E_right_z2;
            E_right_uh2 = xmom_vertex_values[n*3 + (m+2)%3];
            E_right_vh2 = ymom_vertex_values[n*3 + (m+2)%3];


        }



        // Now we have values for this edge - both from left and right side.

        if (elements[Doptimise_dry_cells]) {
            // Check if flux calculation is necessary across this edge
            // This check will exclude dry cells.
            // This will also optimise cases where zl != zr as
            // long as both are dry

            if (fabs(E_left_h) < elements[Depsilon] &&
                    fabs(E_right_h) < elements[Depsilon]) {
                // Cell boundary is dry

                already_computed_flux[k3i] = elements[Dcall]; // #k Done
                if (n >= 0) {
                    already_computed_flux[n3m] = elements[Dcall]; // #n Done
                }

                max_speed = 0.0;
                continue;
            }
        }


        if (fabs(E_left_z - E_right_z) > 1.0e-10) {
            report_python_error(AT, "Discontinuous Elevation");
            //return 0.0;
            return;
        }

        // Outward pointing normal vector (domain.normals[k, 2*i:2*i+2])
        k2i = 2 * k3i; //k*6 + i*2



        // Edge flux computation (triangle k, edge i)
        _flux_function_central_wb_3(
                &E_left_cell_id,
                &E_left_edge_id,
                &E_left_w,
                &E_left_h,
                &E_left_z,
                &E_left_uh,
                &E_left_vh,
                &E_left_u,
                &E_left_v,
                &E_left_w1,
                &E_left_w2,
                &E_left_h1,
                &E_left_h2,
                &E_left_z1,
                &E_left_z2,
                &E_left_uh1,
                &E_left_uh2,
                &E_left_vh1,
                &E_left_vh2,
                &E_left_u1,
                &E_left_u2,
                &E_left_v1,
                &E_left_v2,

                &E_right_cell_id,
                &E_right_edge_id,
                &E_right_w,
                &E_right_h,
                &E_right_z,
                &E_right_uh,
                &E_right_vh,
                &E_right_u,
                &E_right_v,
                &E_right_w1,
                &E_right_w2,
                &E_right_h1,
                &E_right_h2,
                &E_right_z1,
                &E_right_z2,
                &E_right_uh1,
                &E_right_uh2,
                &E_right_vh1,
                &E_right_vh2,
                &E_right_u1,
                &E_right_u2,
                &E_right_v1,
                &E_right_v2,
                normals[k2i], normals[k2i + 1],
                elements[Depsilon], h0, limiting_threshold, elements[Dg],
                edgeflux, &max_speed);



        // Multiply edgeflux by edgelength
        length = edgelengths[k3i];
        edgeflux[0] *= length;
        edgeflux[1] *= length;
        edgeflux[2] *= length;


        // Update triangle k with flux from edge i
        stage_explicit_update[k] -= edgeflux[0];
        xmom_explicit_update[k]  -= edgeflux[1];
        ymom_explicit_update[k]  -= edgeflux[2];

        already_computed_flux[k3i] = elements[Dcall]; // #k Done

        //printf("k i n m %i %i %i %i \n",k,i,n,m);

        // Update neighbour n with same flux but reversed sign
        if (n >= 0) {

            stage_explicit_update[n] += edgeflux[0];
            xmom_explicit_update[n]  += edgeflux[1];
            ymom_explicit_update[n]  += edgeflux[2];


            already_computed_flux[n3m] = elements[Dcall]; // #n Done
        }

        // Update timestep based on edge i and possibly neighbour n
        if (tri_full_flag[k] == 1) {
            if (max_speed > elements[Depsilon]) {
                // Apply CFL condition for triangles joining this edge (triangle k and triangle n)

                // CFL for triangle k
                timestep = min(timestep, radii[k] / max_speed);

                if (n >= 0) {
                    // Apply CFL condition for neigbour n (which is on the ith edge of triangle k)
                    timestep = min(timestep, radii[n] / max_speed);
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
    inv_area = 1.0 / areas[k];
    stage_explicit_update[k] *= inv_area;
    xmom_explicit_update[k] *= inv_area;
    ymom_explicit_update[k] *= inv_area;


    // Keep track of maximal speeds
    max_speed_array[k] = max_speed;

    //} // End triangle k


    //return timestep;
}



#ifndef MAIN_COMPUTE_FLUEXS
#define MAIN_COMPUTE_FLUEXS
int main()
{
    compute_fluxes_central_call = compute_fluxes_central_call;
    compute_fluxes_central_structure_call = compute_fluxes_central_structure_call;
    compute_fluxes_central_wb_call = compute_fluxes_central_wb_call;
    compute_fluxes_central_wb_3_call = compute_fluxes_central_wb_3_call;
}
#endif

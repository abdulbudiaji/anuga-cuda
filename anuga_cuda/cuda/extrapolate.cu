//#include<Python.h>
//#include "extrapolate.h"

#define    W = 16
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define AT __FILE__ ":" TOSTRING(__LINE__)
#define P_ERROR_BUFFER_SIZE 65

__device__ void report_python_error(const char *location, const char *msg)
{

    //    char buf[P_ERROR_BUFFER_SIZE];

    //    snprintf(buf, P_ERROR_BUFFER_SIZE, "Error at %s: %s\n", location, msg);

    //    PyErr_SetString(PyExc_RuntimeError, buf);
}

__device__ int limit_gradient(double *dqv, double qmin, double qmax, double beta_w) {
    // Given provisional jumps dqv from the FV triangle centroid to its
    // vertices and jumps qmin (qmax) between the centroid of the FV
    // triangle and the minimum (maximum) of the values at the centroid of
    // the FV triangle and the auxiliary triangle vertices,
    // calculate a multiplicative factor phi by which the provisional
    // vertex jumps are to be limited

    int i;
    double r = 1000.0, r0 = 1.0, phi = 1.0;
    //static double TINY = 1.0e-100; // to avoid machine accuracy problems.
    double TINY = 1.0e-100; // to avoid machine accuracy problems.
    // FIXME: Perhaps use the epsilon used elsewhere.

    // Any provisional jump with magnitude < TINY does not contribute to
    // the limiting process.
    for (i = 0; i < 3; i++) {
        if (dqv[i]<-TINY)
            r0 = qmin / dqv[i];

        if (dqv[i] > TINY)
            r0 = qmax / dqv[i];

        r = min(r0, r);
    }

    phi = min(r*beta_w, 1.0);
    //for (i=0;i<3;i++)
    dqv[0] = dqv[0] * phi;
    dqv[1] = dqv[1] * phi;
    dqv[2] = dqv[2] * phi;

    return 0;
}

__device__ int find_qmin_and_qmax(double dq0, double dq1, double dq2,
        double *qmin, double *qmax) {
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

    if (dq0 >= 0.0) {
        if (dq1 >= dq2) {
            if (dq1 >= 0.0)
                *qmax = dq0 + dq1;
            else
                *qmax = dq0;

            *qmin = dq0 + dq2;
            if (*qmin >= 0.0) *qmin = 0.0;
        } else {// dq1<dq2
            if (dq2 > 0)
                *qmax = dq0 + dq2;
            else
                *qmax = dq0;

            *qmin = dq0 + dq1;
            if (*qmin >= 0.0) *qmin = 0.0;
        }
    } else {//dq0<0
        if (dq1 <= dq2) {
            if (dq1 < 0.0)
                *qmin = dq0 + dq1;
            else
                *qmin = dq0;

            *qmax = dq0 + dq2;
            if (*qmax <= 0.0) *qmax = 0.0;
        } else {// dq1>dq2
            if (dq2 < 0.0)
                *qmin = dq0 + dq2;
            else
                *qmin = dq0;

            *qmax = dq0 + dq1;
            if (*qmax <= 0.0) *qmax = 0.0;
        }
    }
    return 0;
}


#define Eepsilon                0
#define Eminimum_allowed_height 1
#define Ebeta_w                 2
#define Ebeta_w_dry             3
#define Ebeta_uh                4
#define Ebeta_uh_dry            5
#define Ebeta_vh                6
#define Ebeta_vh_dry            7
#define Eoptimise_dry_cells     8

/*
 * elements[0] = _epsilon
 * elements[1] = _minimum_allowed_height
 * elements[2] = _beta_w
 * elements[3] = _beta_w_dry
 * elements[4] = _beta_uh
 * elements[5] = _beta_uh_dry
 * elements[6] = _beta_vh
 * elements[7] = _beta_vh_dry
 * elements[8] = _optimise_dry_cells
 */


__global__ void _extrapolate_second_order_sw_TRUE(
        double * elements,
        long * surrogate_neighbours,
        long * number_of_boundaries,
        double * centroid_coordinates,
        double * bed_centroid_coordinates,
        double * bed_centroid_values,
        double * stage_centroid_values,
        double * xmom_centroid_values,
        double * ymom_centroid_values,
        double * elevation_centroid_values,
        double * vertex_coordinates,
        double * bed_vertex_values,
        double * stage_vertex_values,
        double * xmom_vertex_values,
        double * ymom_vertex_values,
        double * elevation_vertex_values,

        double * xmom_centroid_store, 
        double * ymom_centroid_store, 
        double * stage_centroid_store)
{
    int k = threadIdx.x + threadIdx.y + blockIdx.x * blockDim.x + blockIdx.y *blockDim.y;

    // Local variables
    double a, b; // Gradient vector used to calculate edge values from centroids
    int k0, k1, k2, k3, k6, coord_index, i;
    double x, y, x0, y0, x1, y1, x2, y2, xv0, yv0, xv1, yv1, xv2, yv2; // Vertices of the auxiliary triangle
    double dx1, dx2, dy1, dy2, dxv0, dxv1, dxv2, dyv0, dyv1, dyv2, dq0, dq1, dq2, area2, inv_area2;
    double dqv[3], qmin, qmax, hmin, hmax;
    double hc, h0, h1, h2, beta_tmp, hfactor;
    //double dk, dv0, dv1, dv2, de[3], demin, dcmax, r0scale;
    double dk, dv0, dv1, dv2;


    //weng_true if (extrapolate_velocity_second_order == 1) {
    // Replace momentum centroid with velocity centroid to allow velocity
    // extrapolation This will be changed back at the end of the routine
    //weng_true for (k = 0; k < number_of_elements; k++) {

    dk = max(stage_centroid_values[k] - bed_centroid_values[k], elements[Eminimum_allowed_height]);
    xmom_centroid_store[k] = xmom_centroid_values[k];
    xmom_centroid_values[k] = xmom_centroid_values[k] / dk;

    ymom_centroid_store[k] = ymom_centroid_values[k];
    ymom_centroid_values[k] = ymom_centroid_values[k] / dk;
    //weng_true }
    //weng_true }

    // Begin extrapolation routine
    //weng_true for (k = 0; k < number_of_elements; k++) {
    k3 = k * 3;
    k6 = k * 6;

    if (number_of_boundaries[k] == 3) {
        // No neighbours, set gradient on the triangle to zero

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
    } else {
        // Triangle k has one or more neighbours.
        // Get centroid and vertex coordinates of the triangle

        // Get the vertex coordinates
        xv0 = vertex_coordinates[k6];
        yv0 = vertex_coordinates[k6 + 1];
        xv1 = vertex_coordinates[k6 + 2];
        yv1 = vertex_coordinates[k6 + 3];
        xv2 = vertex_coordinates[k6 + 4];
        yv2 = vertex_coordinates[k6 + 5];

        // Get the centroid coordinates
        coord_index = 2 * k;
        x = centroid_coordinates[coord_index];
        y = centroid_coordinates[coord_index + 1];

        // Store x- and y- differentials for the vertices of
        // triangle k relative to the centroid
        dxv0 = xv0 - x;
        dxv1 = xv1 - x;
        dxv2 = xv2 - x;
        dyv0 = yv0 - y;
        dyv1 = yv1 - y;
        dyv2 = yv2 - y;
    }




    if (number_of_boundaries[k] <= 1) {
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

        // If the mesh is 'weird' near the boundary,
        // the triangle might be flat or clockwise
        // Default to zero gradient
        if (area2 <= 0) {
            //printf("Error negative triangle area \n");
            //report_python_error(AT, "Negative triangle area");
            //return -1;

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

        if (elements[Eoptimise_dry_cells]) {
            // Check if linear reconstruction is necessary for triangle k
            // This check will exclude dry cells.

            hmax = max(h0, max(h1, h2));
            if (hmax < elements[Eepsilon]) {
                //continue;
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
        //beta_tmp = beta_w_dry;
        //if (hmin>minimum_allowed_height)
        beta_tmp = elements[Ebeta_w_dry] + (elements[Ebeta_w] - elements[Ebeta_w_dry]) * hfactor;

        //printf("min_alled_height = %f\n",minimum_allowed_height);
        //printf("hmin = %f\n",hmin);
        //printf("beta_w = %f\n",beta_w);
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
        //beta_tmp = beta_uh;
        //if (hmin<minimum_allowed_height)
        //beta_tmp = beta_uh_dry;
        beta_tmp = elements[Ebeta_uh_dry] + (elements[Ebeta_uh] - elements[Ebeta_uh_dry]) * hfactor;

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

        //beta_tmp = beta_vh;
        //
        //if (hmin<minimum_allowed_height)
        //beta_tmp = beta_vh_dry;
        beta_tmp = elements[Ebeta_vh_dry] + (elements[Ebeta_vh] - elements[Ebeta_vh_dry]) * hfactor;

        // Limit the gradient
        limit_gradient(dqv, qmin, qmax, beta_tmp);

        for (i = 0; i < 3; i++) {
            ymom_vertex_values[k3 + i] = ymom_centroid_values[k] + dqv[i];
        }
    }// End number_of_boundaries <=1
    else {

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
            //report_python_error(AT, "Internal neighbour not found");
            //return -1;
            return;
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
        limit_gradient(dqv, qmin, qmax, elements[Ebeta_w]);

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
        limit_gradient(dqv, qmin, qmax, elements[Ebeta_w]);

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
        limit_gradient(dqv, qmin, qmax, elements[Ebeta_w]);

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




    //weng_true } // for k=0 to number_of_elements-1

    //weng_true if (extrapolate_velocity_second_order == 1) {
    // Convert back from velocity to momentum
    //weng_true for (k = 0; k < number_of_elements; k++) {
    k3 = 3 * k;
    //dv0 = max(stage_vertex_values[k3]-bed_vertex_values[k3],minimum_allowed_height);
    //dv1 = max(stage_vertex_values[k3+1]-bed_vertex_values[k3+1],minimum_allowed_height);
    //dv2 = max(stage_vertex_values[k3+2]-bed_vertex_values[k3+2],minimum_allowed_height);
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

    //weng_true }
    //weng_true }

}




/*
 * elements[0] = _epsilon
 * elements[1] = _minimum_allowed_height
 * elements[2] = _beta_w
 * elements[3] = _beta_w_dry
 * elements[4] = _beta_uh
 * elements[5] = _beta_uh_dry
 * elements[6] = _beta_vh
 * elements[7] = _beta_vh_dry
 * elements[8] = _optimise_dry_cells
 */


// Computational routine
//__global__ int _extrapolate_second_order_edge_sw(struct domain *D) 
__global__ void _extrapolate_second_order_edge_sw_TRUE(
        double * elements,
        long * surrogate_neighbours,
        long * number_of_boundaries,
        double * centroid_coordinates,
        double * stage_centroid_values,
        double * xmom_centroid_values,
        double * ymom_centroid_values,
        double * bed_centroid_values,
        double * vertex_coordinates,
        double * bed_vertex_values,
        double * stage_vertex_values,
        double * xmom_vertex_values,
        double * ymom_vertex_values,
        double * elevation_vertex_values,
        double * edge_coordinates,
        double * bed_edge_values,
        double * stage_edge_values,
        double * xmom_edge_values,
        double * ymom_edge_values,
        double * xmom_centroid_store, 
        double * ymom_centroid_store,
        double * stage_centroid_store)
{
    // Local variables
    int k = threadIdx.x + threadIdx.y + blockIdx.x * blockDim.x + blockIdx.y *blockDim.y;

    double a, b; // Gradient vector used to calculate edge values from centroids
    int k0, k1, k2, k3, k6, coord_index, i;
    double x, y, x0, y0, x1, y1, x2, y2, xv0, yv0, xv1, yv1, xv2, yv2; // Vertices of the auxiliary triangle
    double dx1, dx2, dy1, dy2, dxv0, dxv1, dxv2, dyv0, dyv1, dyv2, dq0, dq1, dq2, area2, inv_area2;
    double dqv[3], qmin, qmax, hmin, hmax;
    double hc, h0, h1, h2, beta_tmp, hfactor;
    //double dk, dv0, dv1, dv2, de[3], demin, dcmax, r0scale;
    double dk, de[3];



    //weng_true if(extrapolate_velocity_second_order==1){
    // Replace momentum centroid with velocity centroid to allow velocity
    // extrapolation This will be changed back at the end of the routine
    //weng_true for (k=0; k<number_of_elements; k++){

    dk = max(stage_centroid_values[k]-bed_centroid_values[k], elements[Eminimum_allowed_height]);
    xmom_centroid_store[k] = xmom_centroid_values[k];
    xmom_centroid_values[k] = xmom_centroid_values[k]/dk;

    ymom_centroid_store[k] = ymom_centroid_values[k];
    ymom_centroid_values[k] = ymom_centroid_values[k]/dk;

    //weng_true }
    //weng_true }

    // Begin extrapolation routine
    //weng_true for (k = 0; k < number_of_elements; k++){
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

        //continue;
        return;
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

        // Test to see whether we accept the surrogate neighbours
        // Note that if ki is replaced with k in more than 1 neighbour, then the
        // triangle area will be zero, and a first order extrapolation will be
        // used
        if(stage_centroid_values[k2]<=bed_centroid_values[k2]){
            k2 = k ;
        }
        if(stage_centroid_values[k0]<=bed_centroid_values[k0]){
            k0 = k ;
        }
        if(stage_centroid_values[k1]<=bed_centroid_values[k1]){
            k1 = k ;
        }
        // Alternative approach (BRUTAL) -- if the max neighbour bed elevation is greater
        // than the min neighbour stage, then we use first order extrapolation
        //bedmax = max(bed_centroid_values[k],
        //             max(bed_centroid_values[k0],
        //                 max(bed_centroid_values[k1], bed_centroid_values[k2])));
        //stagemin = min(stage_centroid_values[k],
        //             min(stage_centroid_values[k0],
        //                 min(stage_centroid_values[k1], stage_centroid_values[k2])));
        //
        //if(stagemin < bedmax){
        //   // This will cause first order extrapolation
        //   k2 = k;
        //   k0 = k;
        //   k1 = k;
        //}

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

        // compute the maximum distance from the centroid to a neighbouring
        // centroid
        //dcmax=max( (x0-x)*(x0-x) + (y0-y)*(y0-y),
        //           max((x1-x)*(x1-x) + (y1-y)*(y1-y),
        //               (x2-x)*(x2-x) + (y2-y)*(y2-y)));
        //dcmax=sqrt(dcmax);
        //// Ratio of centroid to edge distance -- useful in attempting to adapt limiter
        //if(dcmax>0.){
        //    r0scale=demin/dcmax;
        //    //printf("%f \n", r0scale);
        //}else{
        //    r0scale=0.5;
        //}

        // Store x- and y- differentials for the vertices
        // of the auxiliary triangle
        dx1 = x1 - x0;
        dx2 = x2 - x0;
        dy1 = y1 - y0;
        dy2 = y2 - y0;

        // Calculate 2*area of the auxiliary triangle
        // The triangle is guaranteed to be counter-clockwise
        area2 = dy2*dx1 - dy1*dx2;

        // If the mesh is 'weird' near the boundary,
        // the triangle might be flat or clockwise
        // Default to zero gradient
        if (area2 <= 0)
        {
            //printf("Error negative triangle area \n");
            //report_python_error(AT, "Negative triangle area");
            //return -1;

            stage_edge_values[k3]   = stage_centroid_values[k];
            stage_edge_values[k3+1] = stage_centroid_values[k];
            stage_edge_values[k3+2] = stage_centroid_values[k];
            xmom_edge_values[k3]    = xmom_centroid_values[k];
            xmom_edge_values[k3+1]  = xmom_centroid_values[k];
            xmom_edge_values[k3+2]  = xmom_centroid_values[k];
            ymom_edge_values[k3]    = ymom_centroid_values[k];
            ymom_edge_values[k3+1]  = ymom_centroid_values[k];
            ymom_edge_values[k3+2]  = ymom_centroid_values[k];

            //continue;
            return;
        }

        // Calculate heights of neighbouring cells
        hc = stage_centroid_values[k]  - bed_centroid_values[k];
        h0 = stage_centroid_values[k0] - bed_centroid_values[k0];
        h1 = stage_centroid_values[k1] - bed_centroid_values[k1];
        h2 = stage_centroid_values[k2] - bed_centroid_values[k2];
        hmin = min(min(h0, min(h1, h2)), hc);
        //hmin = min(h0, min(h1, h2));
        //hmin = max(hmin, 0.0);
        //hfactor = hc/(hc + 1.0);

        hfactor = 0.0;
        //if (hmin > 0.001)
        if (hmin > 0.)
            //if (hc>0.0)
        {
            hfactor = 1.0 ;//hmin/(hmin + 0.004);
            //hfactor=hmin/(hmin + 0.004);
        }

        if (elements[Eoptimise_dry_cells])
        {
            // Check if linear reconstruction is necessary for triangle k
            // This check will exclude dry cells.

            //hmax = max(h0, max(h1, max(h2, hc)));
            hmax = max(h0, max(h1, h2));
            if (hmax < elements[Eepsilon])
            {
                //continue;
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

        beta_tmp = elements[Ebeta_w_dry] + (elements[Ebeta_w] - elements[Ebeta_w_dry]) * hfactor;

        // Limit the gradient
        limit_gradient(dqv, qmin, qmax, beta_tmp);
        //limit_gradient2(dqv, qmin, qmax, beta_tmp,r0scale);

        //for (i=0;i<3;i++)
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

        beta_tmp = elements[Ebeta_uh_dry] + (elements[Ebeta_uh] - elements[Ebeta_uh_dry]) * hfactor;

        // Limit the gradient
        limit_gradient(dqv, qmin, qmax, beta_tmp);
        //limit_gradient2(dqv, qmin, qmax, beta_tmp,r0scale);


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

        beta_tmp = elements[Ebeta_vh_dry] + (elements[Ebeta_vh] - elements[Ebeta_vh_dry]) * hfactor;

        // Limit the gradient
        limit_gradient(dqv, qmin, qmax, beta_tmp);
        //limit_gradient2(dqv, qmin, qmax, beta_tmp,r0scale);

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
            // report_python_error(AT, "Internal neighbour not found");
            //return -1;
            return;
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
        limit_gradient(dqv, qmin, qmax, elements[Ebeta_w]);

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
        limit_gradient(dqv, qmin, qmax, elements[Ebeta_w]);

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
        limit_gradient(dqv, qmin, qmax, elements[Ebeta_w]);

        for (i=0;i<3;i++)
        {
            ymom_edge_values[k3 + i] = ymom_centroid_values[k] + dqv[i];
        }
    } // else [number_of_boundaries==2]
    //weng_true } // for k=0 to number_of_elements-1

    // Compute vertex values of quantities
    //weng_true for (k=0; k<number_of_elements; k++){
    k3 = 3*k;

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
    //weng_true if(extrapolate_velocity_second_order==1){
    //Convert velocity back to momenta at centroids
    xmom_centroid_values[k] = xmom_centroid_store[k];
    ymom_centroid_values[k] = ymom_centroid_store[k];

    // Re-compute momenta at edges
    for (i=0; i<3; i++){
        de[i] = max(stage_edge_values[k3+i]-bed_edge_values[k3+i],0.0);
        xmom_edge_values[k3+i]=xmom_edge_values[k3+i]*de[i];
        ymom_edge_values[k3+i]=ymom_edge_values[k3+i]*de[i];
    }

    // Re-compute momenta at vertices
    for (i=0; i<3; i++){
        de[i] = max(stage_vertex_values[k3+i]-bed_vertex_values[k3+i],0.0);
        xmom_vertex_values[k3+i]=xmom_vertex_values[k3+i]*de[i];
        ymom_vertex_values[k3+i]=ymom_vertex_values[k3+i]*de[i];
    }
    //weng_true }
    //weng_true }
}


#ifndef MAIN_EXTRAPOLATE
#define MAIN_EXTRAPOLATE
int main()
{}
#endif

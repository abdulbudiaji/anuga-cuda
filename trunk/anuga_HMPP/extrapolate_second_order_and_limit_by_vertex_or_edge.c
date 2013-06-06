#include "fun.h"

//#define NON_DIRECTIVES_EXTRA2_VERTEX
//#define NON_DIRECTIVES_EXTRA2_VERTEX_EXTRA_FROM_GRA

double fmax(double x, double y) {  
    //Return maximum of two doubles

    if (x > y) return x;
    else return y;
}


double fmin(double x, double y) {  
    //Return minimum of two doubles

    if (x < y) return x;
    else return y;
}



int gradient(double x0, double y0, 
        double x1, double y1, 
        double x2, double y2, 
        double q0, double q1, double q2, 
        double *a, double *b) {

    /*Compute gradient (a,b) based on three points (x0,y0), (x1,y1) and (x2,y2) 
      with values q0, q1 and q2.

      Extrapolation formula (q0 is selected as an arbitrary origin)
      q(x,y) = q0 + a*(x-x0) + b*(y-y0)                    (1)

      Substituting the known values for q1 and q2 into (1) yield the 
      equations for a and b 

      q1-q0 = a*(x1-x0) + b*(y1-y0)                      (2)
      q2-q0 = a*(x2-x0) + b*(y2-y0)                      (3)      

      or in matrix form

      /               \  /   \   /       \  
      |  x1-x0  y1-y0 |  | a |   | q1-q0 |
      |               |  |   | = |       | 
      |  x2-x0  y2-y0 |  | b |   | q2-q0 |
      \               /  \   /   \       /

      which is solved using the standard determinant technique    

     */


    double det;

    det = (y2-y0)*(x1-x0) - (y1-y0)*(x2-x0);

    *a = (y2-y0)*(q1-q0) - (y1-y0)*(q2-q0);
    *a /= det;

    *b = (x1-x0)*(q2-q0) - (x2-x0)*(q1-q0);
    *b /= det;

    return 0;
}




int gradient2(double x0, double y0, 
        double x1, double y1, 
        double q0, double q1, 
        double *a, double *b) {
    /*Compute gradient (a,b) between two points (x0,y0) and (x1,y1) 
      with values q0 and q1 such that the plane is constant in the direction 
      orthogonal to (x1-x0, y1-y0).

      Extrapolation formula
      q(x,y) = q0 + a*(x-x0) + b*(y-y0)                    (1)

      Substituting the known values for q1 into (1) yields an 
      under determined  equation for a and b 
      q1-q0 = a*(x1-x0) + b*(y1-y0)                      (2)


      Now add the additional requirement that the gradient in the direction 
      orthogonal to (x1-x0, y1-y0) should be zero. The orthogonal direction 
      is given by the vector (y0-y1, x1-x0).

      Define the point (x2, y2) = (x0 + y0-y1, y0 + x1-x0) on the orthognal line. 
      Then we know that the corresponding value q2 should be equal to q0 in order 
      to obtain the zero gradient, hence applying (1) again    
      q0 = q2 = q(x2, y2) = q0 + a*(x2-x0) + b*(y2-y0)
      = q0 + a*(x0 + y0-y1-x0) + b*(y0 + x1-x0 - y0)
      = q0 + a*(y0-y1) + b*(x1-x0)

      leads to the orthogonality constraint
      a*(y0-y1) + b*(x1-x0) = 0                           (3) 

      which closes the system and yields

      /               \  /   \   /       \  
      |  x1-x0  y1-y0 |  | a |   | q1-q0 |
      |               |  |   | = |       | 
      |  y0-y1  x1-x0 |  | b |   |   0   |
      \               /  \   /   \       /

      which is solved using the standard determinant technique    

     */

    double det, xx, yy, qq;

    xx = x1-x0;
    yy = y1-y0;
    qq = q1-q0;

    det = xx*xx + yy*yy;  //FIXME  catch det == 0
    *a = xx*qq/det;
    *b = yy*qq/det;

    return 0;
}


#ifndef NON_DIRECTIVES_EXTRA2_VERTEX_CPTGRA
#pragma hmpp cptGradients codelet, target=CUDA args[*].transfer=atcall
#endif
void _compute_gradients(
        int N,
        int N2,
        int N3,
        double centroids[N2],
        double centroid_values[N],
        long number_of_boundaries[N],
        long surrogate_neighbours[N3],
        double a[N],
        double b[N])
{
    int i, k, k0, k1, k2, k3;
    double x0, x1, x2, y0, y1, y2, q0, q1, q2; //, det;


    #pragma hmppcg gridify(k), &
    #pragma hmppcg & private( i, k0, k1, k2, k3, &
    #pragma hmppcg & x0, x1, x2, y0, y1, y2, q0, q1, q2), &
    #pragma hmppcg & global( centroids, centroid_values, &
    #pragma hmppcg & number_of_boundaries, surrogate_neighbours, a, b)
    for (k=0; k<N; k++) {
        k3 = 3*k;

        if (number_of_boundaries[k] < 2) {
            // Two or three true neighbours

            // Get indices of neighbours (or self when used as surrogate)
            // k0, k1, k2 = surrogate_neighbours[k,:]

            k0 = surrogate_neighbours[k3 + 0];
            k1 = surrogate_neighbours[k3 + 1];
            k2 = surrogate_neighbours[k3 + 2];


            // FIXME
            //if (k0 == k1 || k1 == k2) return -1;

            // Get data
            q0 = centroid_values[k0];
            q1 = centroid_values[k1];
            q2 = centroid_values[k2];

            x0 = centroids[k0*2]; y0 = centroids[k0*2+1];
            x1 = centroids[k1*2]; y1 = centroids[k1*2+1];
            x2 = centroids[k2*2]; y2 = centroids[k2*2+1];

            // Gradient
            gradient(x0, y0, x1, y1, x2, y2, q0, q1, q2, &a[k], &b[k]);

        } else if (number_of_boundaries[k] == 2) {
            // One true neighbour

            // FIXME
            // Get index of the one neighbour
            //i=0; k0 = k;
            //#pragma hmppcg noParallel
            //while (i<3 && k0==k) {
            //    k0 = surrogate_neighbours[k3 + i];
            //    i++;
            //}
            k0 = surrogate_neighbours[k3];
            if( k0 == k)
            {
                k0 = surrogate_neighbours[k3 + 1];
                if ( k0 == k)
                {
                    k0 = surrogate_neighbours[k3 + 2];
                }
            }
            

            // FIXME
            //if (k0 == k) return -1;

            k1 = k; //self

            // Get data
            q0 = centroid_values[k0];
            q1 = centroid_values[k1];

            x0 = centroids[k0*2]; y0 = centroids[k0*2+1];
            x1 = centroids[k1*2]; y1 = centroids[k1*2+1];

            // Two point gradient
            gradient2(x0, y0, x1, y1, q0, q1, &a[k], &b[k]);

        }
        //    else:
        //        #No true neighbours -
        //        #Fall back to first order scheme
    }
}


#ifndef NON_DIRECTIVES_EXTRA2_VERTEX_EXTRA_FROM_GRA
#pragma hmpp extraFromGradient codelet, target=CUDA args[*].transfer=atcall
#endif
void _extrapolate_from_gradient(
        int N,
        int N2,
        int N3,
        int N6,
        double centroids[N2],
        double centroid_values[N],
        double vertex_coordinates[N6],
        double vertex_values[N3],
        double edge_values[N3],
        double a[N],
        double b[N]) 
{
    int k, k2, k3, k6;
    double x, y, x0, y0, x1, y1, x2, y2;

    #pragma hmppcg gridify(k), &
    #pragma hmppcg & private( k2, k3, k6, &
    #pragma hmppcg & x, y, x0, x1, x2, y0, y1, y2), &
    #pragma hmppcg & global( centroids, centroid_values, &
    #pragma hmppcg & vertex_coordinates, vertex_values, edge_values, a, b)
    for (k=0; k<N; k++){
        k6 = 6*k;
        k3 = 3*k;
        k2 = 2*k;

        // Centroid coordinates
        x = centroids[k2]; y = centroids[k2+1];

        // vertex coordinates
        // x0, y0, x1, y1, x2, y2 = X[k,:]
        x0 = vertex_coordinates[k6 + 0];
        y0 = vertex_coordinates[k6 + 1];
        x1 = vertex_coordinates[k6 + 2];
        y1 = vertex_coordinates[k6 + 3];
        x2 = vertex_coordinates[k6 + 4];
        y2 = vertex_coordinates[k6 + 5];

        // Extrapolate to Vertices
        vertex_values[k3+0] = centroid_values[k] + a[k]*(x0-x) + b[k]*(y0-y);
        vertex_values[k3+1] = centroid_values[k] + a[k]*(x1-x) + b[k]*(y1-y);
        vertex_values[k3+2] = centroid_values[k] + a[k]*(x2-x) + b[k]*(y2-y);

        // Extrapolate to Edges (midpoints)
        edge_values[k3+0] = k3;//0.5*(vertex_values[k3 + 1]+vertex_values[k3 + 2]);
        edge_values[k3+1] = k3;//0.5*(vertex_values[k3 + 2]+vertex_values[k3 + 0]);
        edge_values[k3+2] = k3;//0.5*(vertex_values[k3 + 0]+vertex_values[k3 + 1]);

    }
}



#ifndef NON_DIRECTIVES_EXTRA2_VERTEX_LMT_V
#pragma hmpp lmtVByNeigh codelet, target=CUDA args[*].transfer=atcall
#endif
void _limit_vertices_by_all_neighbours(
        int N, 
        int N3,
        double beta,
        double centroid_values[N],
        double vertex_values[N3],
        double edge_values[N3],
        long   neighbours[N3],
        double x_gradient[N],
        double y_gradient[N]) 
{
    int i, k, k3, k6;
    long n;
    double qmin, qmax, qn, qc;
    double dq, dqa[3], phi, r;

    #pragma hmppcg gridify(k), &
    #pragma hmppcg & private( i, k3, k6, n,&
    #pragma hmppcg & qmin, qmax, qn, qc, dq, dqa, phi, r), &
    #pragma hmppcg & global( centroid_values, vertex_values, edge_values, &
    #pragma hmppcg & neighbours, x_gradient, y_gradient)
    for (k=0; k<N; k++){
        k6 = 6*k;
        k3 = 3*k;

        qc = centroid_values[k];
        qmin = qc;
        qmax = qc;

        #pragma hmppcg noparallel
        for (i=0; i<3; i++) {
            n = neighbours[k3+i];
            if (n >= 0) {
                qn = centroid_values[n]; //Neighbour's centroid value

                qmin = fmin(qmin, qn);
                qmax = fmax(qmax, qn);
            }
        }

        phi = 1.0;
        #pragma hmppcg noparallel
        for (i=0; i<3; i++) {    
            r = 1.0;

            dq = vertex_values[k3+i] - qc;    //Delta between vertex and centroid values
            dqa[i] = dq;                      //Save dq for use in updating vertex values

            if (dq > 0.0) r = (qmax - qc)/dq;
            if (dq < 0.0) r = (qmin - qc)/dq;      


            phi = fmin( fmin(r*beta, 1.0), phi);    
        }

        //Update gradient, vertex and edge values using phi limiter
        x_gradient[k] = x_gradient[k]*phi;
        y_gradient[k] = y_gradient[k]*phi;

        vertex_values[k3+0] = qc + phi*dqa[0];
        vertex_values[k3+1] = qc + phi*dqa[1];
        vertex_values[k3+2] = qc + phi*dqa[2];

        edge_values[k3+0] = 0.5f*(vertex_values[k3+1] + vertex_values[k3+2]);
        edge_values[k3+1] = 0.5f*(vertex_values[k3+2] + vertex_values[k3+0]);
        edge_values[k3+2] = 0.5f*(vertex_values[k3+0] + vertex_values[k3+1]);

    }
}




#ifndef NON_DIRECTIVES_EXTRA2_VERTEX_LMT_E
#pragma hmpp lmtEByNeigh codelet, target=CUDA args[*].transfer=atcall
#endif
void _limit_edges_by_all_neighbours(
        int N, 
        int N3,
        double beta,
        double centroid_values[N],
        double vertex_values[N3],
        double edge_values[N3],
        long   neighbours[N3],
        double x_gradient[N],
        double y_gradient[N]) 
{
    int i, k, k2, k3, k6;
    long n;
    double qmin, qmax, qn, qc, sign;
    double dq, dqa[3], phi, r;

    #pragma hmppcg gridify(k), &
    #pragma hmppcg & private( i, k2, k3, k6, n,&
    #pragma hmppcg & qmin, qmax, qn, qc, sign, dq, dqa, phi, r), &
    #pragma hmppcg & global( beta, &
    #pragma hmppcg & centroid_values, vertex_values, edge_values, &
    #pragma hmppcg & neighbours, x_gradient, y_gradient)
    for (k=0; k<N; k++){
        k6 = 6*k;
        k3 = 3*k;
        k2 = 2*k;

        qc = centroid_values[k];
        qmin = qc;
        qmax = qc;

        for (i=0; i<3; i++) {
            n = neighbours[k3+i];
            if (n >= 0) {
                qn = centroid_values[n]; //Neighbour's centroid value

                qmin = fmin(qmin, qn);
                qmax = fmax(qmax, qn);
            }
        }

        sign = 0.0;
        if (qmin > 0.0) {
            sign = 1.0;
        } else if (qmax < 0) {
            sign = -1.0;
        }

        phi = 1.0;
        for (i=0; i<3; i++) {  
            dq = edge_values[k3+i] - qc;      //Delta between edge and centroid values
            dqa[i] = dq;                      //Save dq for use in updating vertex values  


            // Just limit non boundary edges so that we can reconstruct a linear function
            // FIXME Problem with stability on edges 
            //if (neighbours[k3+i] >= 0) {
            r = 1.0;

            if (dq > 0.0) r = (qmax - qc)/dq;
            if (dq < 0.0) r = (qmin - qc)/dq;      

            phi = fmin( fmin(r*beta, 1.0), phi);
            //	}

            //
            /* if (neighbours[k3+i] < 0) { */
            /* 	r = 1.0; */

            /* 	if (dq > 0.0 && (sign == -1.0 || sign == 0.0 )) r = (0.0 - qc)/dq; */
            /* 	if (dq < 0.0 && (sign ==  1.0 || sign == 0.0 )) r = (0.0 - qc)/dq; */

            /* 	phi = min( min(r*beta, 1.0), phi); */
            /* 	} */

        }

        //Update gradient, vertex and edge values using phi limiter
        x_gradient[k] = x_gradient[k]*phi;
        y_gradient[k] = y_gradient[k]*phi;

        edge_values[k3+0] = qc + phi*dqa[0];
        edge_values[k3+1] = qc + phi*dqa[1];
        edge_values[k3+2] = qc + phi*dqa[2];

        vertex_values[k3+0] = edge_values[k3+1] + edge_values[k3+2] - edge_values[k3+0];
        vertex_values[k3+1] = edge_values[k3+2] + edge_values[k3+0] - edge_values[k3+1];
        vertex_values[k3+2] = edge_values[k3+0] + edge_values[k3+1] - edge_values[k3+2];

    }

}



void extrapolate_second_order_and_limit_by_vertex(
        int N,
        int N2,
        int N3,
        int N6,
        double beta,

        double domain_centroid_coordinates[N2],
        double domain_vertex_coordinates[N6],
        long domain_number_of_boundaries[N],
        long domain_surrogate_neighbours[N3],
        long domain_neighbours[N3],

        double quantity_centroid_values[N],
        double quantity_vertex_values[N3],
        double quantity_edge_values[N3],
        double quantity_x_gradient[N],
        double quantity_y_gradient[N]
        )
{
    #ifndef NON_DIRECTIVES_EXTRA2_VERTEX_CPTGRA
    #pragma hmpp cptGradients callsite
    #endif
    _compute_gradients(
            N,
            N2,
            N3,
            domain_centroid_coordinates,
            quantity_centroid_values,
            domain_number_of_boundaries,
            domain_surrogate_neighbours,
            quantity_x_gradient,
            quantity_y_gradient);

    #ifndef NON_DIRECTIVES_EXTRA2_VERTEX_CPTGRA
    #pragma hmpp cptGradients synchronize
    #endif


    #ifndef NON_DIRECTIVES_EXTRA2_VERTEX_EXTRA_FROM_GRA
    #pragma hmpp extraFromGradient callsite
    #endif
    _extrapolate_from_gradient(
            N,
            N2,
            N3,
            N6,
            domain_centroid_coordinates,
            quantity_centroid_values,
            domain_vertex_coordinates,
            quantity_vertex_values,
            quantity_edge_values,
            quantity_x_gradient,
            quantity_y_gradient);

    #ifndef NON_DIRECTIVES_EXTRA2_VERTEX_EXTRA_FROM_GRA
    #pragma hmpp extraFromGradient synchronize
    #endif

/*
    #ifndef NON_DIRECTIVES_EXTRA2_VERTEX_LMT_V
    #pragma hmpp lmtVByNeigh callsite
    #endif
    _limit_vertices_by_all_neighbours(
            N,
            N3,
            beta,
            quantity_centroid_values,
            quantity_vertex_values,
            quantity_edge_values,
            domain_neighbours,
            quantity_x_gradient,
            quantity_y_gradient);
    
    #ifndef NON_DIRECTIVES_EXTRA2_VERTEX_LMT_V
    #pragma hmpp lmtVByNeigh synchronize
    #endif
*/

}


void extrapolate_second_order_and_limit_by_vertex_normal(
        int N,
        int N2,
        int N3,
        int N6,
        double beta,

        double domain_centroid_coordinates[N2],
        double domain_vertex_coordinates[N6],
        long domain_number_of_boundaries[N],
        long domain_surrogate_neighbours[N3],
        long domain_neighbours[N3],

        double quantity_centroid_values[N],
        double quantity_vertex_values[N3],
        double quantity_edge_values[N3],
        double quantity_x_gradient[N],
        double quantity_y_gradient[N]
        )
{
    _compute_gradients(
            N,
            N2,
            N3,
            domain_centroid_coordinates,
            quantity_centroid_values,
            domain_number_of_boundaries,
            domain_surrogate_neighbours,
            quantity_x_gradient,
            quantity_y_gradient);

    _extrapolate_from_gradient(
            N,
            N2,
            N3,
            N6,
            domain_centroid_coordinates,
            quantity_centroid_values,
            domain_vertex_coordinates,
            quantity_vertex_values,
            quantity_edge_values,
            quantity_x_gradient,
            quantity_y_gradient);

/*
    _limit_vertices_by_all_neighbours(
            N,
            N3,
            beta,
            quantity_centroid_values,
            quantity_vertex_values,
            quantity_edge_values,
            domain_neighbours,
            quantity_x_gradient,
            quantity_y_gradient);

*/
}

void extrapolate_second_order_and_limit_by_edge(
        int N,
        int N2,
        int N3,
        int N6,
        double beta,

        double domain_centroid_coordinates[N2],
        double domain_vertex_coordinates[N6],
        long domain_number_of_boundaries[N],
        long domain_surrogate_neighbours[N3],
        long domain_neighbours[N3],

        double quantity_centroid_values[N],
        double quantity_vertex_values[N3],
        double quantity_edge_values[N3],
        double quantity_x_gradient[N],
        double quantity_y_gradient[N]
        )
{
    #ifndef NON_DIRECTIVES_EXTRA2_VERTEX_CPTGRA
    #pragma hmpp cptGradients callsite
    #endif
    _compute_gradients(
            N,
            N2,
            N3,
            domain_centroid_coordinates,
            quantity_centroid_values,
            domain_number_of_boundaries,
            domain_surrogate_neighbours,
            quantity_x_gradient,
            quantity_y_gradient);


    #ifndef NON_DIRECTIVES_EXTRA2_VERTEX_EXTRA_FROM_GRA
    #pragma hmpp extraFromGradient callsite
    #endif
    _extrapolate_from_gradient(
            N,
            N2,
            N3,
            N6,
            domain_centroid_coordinates,
            quantity_centroid_values,
            domain_vertex_coordinates,
            quantity_vertex_values,
            quantity_edge_values,
            quantity_x_gradient,
            quantity_y_gradient);

    #ifndef NON_DIRECTIVES_EXTRA2_VERTEX_LMT_E
    #pragma hmpp lmtEByNeigh callsite
    #endif
    _limit_edges_by_all_neighbours(
            N,
            N3,
            beta,
            quantity_centroid_values,
            quantity_vertex_values,
            quantity_edge_values,
            domain_neighbours,
            quantity_x_gradient,
            quantity_y_gradient);

}




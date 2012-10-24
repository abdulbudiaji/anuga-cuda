__device__ int _gradient(double x0, double y0, 
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


__device__ int _gradient2(double x0, double y0, 
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



__device__ int _limit_vertices_by_all_neighbours(
        //int N, 
        double beta,
        double* centroid_values,
        double* vertex_values,
        double* edge_values,
        long*   neighbours,
        double* x_gradient,
        double* y_gradient) {

    const int k = threadIdx.x + blockIdx.x * blockDim.x;

    int i, k3;
    long n;
    double qmin, qmax, qn, qc;
    double dq, dqa[3], phi, r;

    //for (k=0; k<N; k++){
        k3 = 3*k;

        qc = centroid_values[k];
        qmin = qc;
        qmax = qc;

        for (i=0; i<3; i++) {
            n = neighbours[k3+i];
            if (n >= 0) {
                qn = centroid_values[n]; //Neighbour's centroid value

                qmin = min(qmin, qn);
                qmax = max(qmax, qn);
            }
        }

        phi = 1.0;
        for (i=0; i<3; i++) {    
            r = 1.0;

            dq = vertex_values[k3+i] - qc;    //Delta between vertex and centroid values
            dqa[i] = dq;                      //Save dq for use in updating vertex values

            if (dq > 0.0) r = (qmax - qc)/dq;
            if (dq < 0.0) r = (qmin - qc)/dq;      


            phi = min( min(r*beta, 1.0), phi);    
        }

        //Update gradient, vertex and edge values using phi limiter
        x_gradient[k] = x_gradient[k]*phi;
        y_gradient[k] = y_gradient[k]*phi;

        vertex_values[k3+0] = qc + phi*dqa[0];
        vertex_values[k3+1] = qc + phi*dqa[1];
        vertex_values[k3+2] = qc + phi*dqa[2];

        edge_values[k3+0] = 0.5*(vertex_values[k3+1] + vertex_values[k3+2]);
        edge_values[k3+1] = 0.5*(vertex_values[k3+2] + vertex_values[k3+0]);
        edge_values[k3+2] = 0.5*(vertex_values[k3+0] + vertex_values[k3+1]);

    //}

    return 0;
}

__device__ int _extrapolate_from_gradient(
        //int N,
        double* centroids,
        double* centroid_values,
        double* vertex_coordinates,
        double* vertex_values,
        double* edge_values,
        double* a,
        double* b) {
    
    const int k = threadIdx.x + blockIdx.x * blockDim.x;

    int k2, k3, k6;
    double x, y, x0, y0, x1, y1, x2, y2;

    //for (k=0; k<N; k++){
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
        edge_values[k3+0] = 0.5*(vertex_values[k3 + 1]+vertex_values[k3 + 2]);
        edge_values[k3+1] = 0.5*(vertex_values[k3 + 2]+vertex_values[k3 + 0]);
        edge_values[k3+2] = 0.5*(vertex_values[k3 + 0]+vertex_values[k3 + 1]);

    //}
    return 0;
}


__device__ int _compute_gradients(
        //int N,
        double* centroids,
        double* centroid_values,
        long* number_of_boundaries,
        long* surrogate_neighbours,
        double* a,
        double* b)
{
    const int k = threadIdx.x + blockIdx.x * blockDim.x;
    int i, k0, k1, k2;
    double x0, x1, x2, y0, y1, y2, q0, q1, q2;
    double det;


//    for (k=0; k<N; k++) {

        if (number_of_boundaries[k] < 2) {
            // Two or three true neighbours

            // Get indices of neighbours (or self when used as surrogate)
            // k0, k1, k2 = surrogate_neighbours[k,:]

            k0 = surrogate_neighbours[3*k + 0];
            k1 = surrogate_neighbours[3*k + 1];
            k2 = surrogate_neighbours[3*k + 2];


            if (k0 == k1 || k1 == k2) return -1;

            // Get data
            q0 = centroid_values[k0];
            q1 = centroid_values[k1];
            q2 = centroid_values[k2];

            x0 = centroids[k0*2]; y0 = centroids[k0*2+1];
            x1 = centroids[k1*2]; y1 = centroids[k1*2+1];
            x2 = centroids[k2*2]; y2 = centroids[k2*2+1];

            // Gradient
            //_gradient(x0, y0, x1, y1, x2, y2, q0, q1, q2, &a[k], &b[k]);

            det = (y2-y0)*(x1-x0) - (y1-y0)*(x2-x0);

            a[k] = (y2-y0)*(q1-q0) - (y1-y0)*(q2-q0);
            a[k] /= det;

            b[k] = (x1-x0)*(q2-q0) - (x2-x0)*(q1-q0);
            b[k] /= det;

        } else if (number_of_boundaries[k] == 2) {
            // One true neighbour

            // Get index of the one neighbour
            i=0; k0 = k;
            while (i<3 && k0==k) {
                k0 = surrogate_neighbours[3*k + i];
                i++;
            }
            if (k0 == k) return -1;

            k1 = k; //self

            // Get data
            q0 = centroid_values[k0];
            q1 = centroid_values[k1];

            x0 = centroids[k0*2]; y0 = centroids[k0*2+1];
            x1 = centroids[k1*2]; y1 = centroids[k1*2+1];

            // Two point gradient
            //_gradient2(x0, y0, x1, y1, q0, q1, &a[k], &b[k]);

            det = (x1-x0) * (x1-x0) + (y1-y0)*(y1-y0);
            a[k] = (x1-x0)*(q1-q0) / det;
            b[k] = (y1-y0)*(q1-q0) / det;

        }
        //    else:
        //        #No true neighbours -
        //        #Fall back to first order scheme
//    }
    return 0;
}




__global__ void extrapolate_second_order_and_limit_by_vertex(
        double beta,
        double * domain_centroid_coordinates,
        double * domain_vertex_coordinates,
        long * domain_number_of_boundaries,
        long * domain_surrogate_neighbours,
        long * domain_neighbours,
        
        double * quantity_centroid_values,
        double * quantity_vertex_values,
        double * quantity_edge_values,
        double * quantity_x_gradient,
        double * quantity_y_gradient
        )
{
    _compute_gradients(
            domain_centroid_coordinates,
            quantity_centroid_values,
            domain_number_of_boundaries,
            domain_surrogate_neighbours,
            quantity_x_gradient,
            quantity_y_gradient);

    _extrapolate_from_gradient(
            domain_centroid_coordinates,
            quantity_centroid_values,
            domain_vertex_coordinates,
            quantity_vertex_values,
            quantity_edge_values,
            quantity_x_gradient,
            quantity_y_gradient);

    _limit_vertices_by_all_neighbours(
            beta,
            quantity_centroid_values,
            quantity_vertex_values,
            quantity_edge_values,
            domain_neighbours,
            quantity_x_gradient,
            quantity_y_gradient);
}


#ifndef MAIN_EXTRO
#define MAIN_EXTRO
int main()
{}
#endif

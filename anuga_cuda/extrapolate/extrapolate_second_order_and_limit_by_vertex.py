#!/usr/bin/env python

"""
    This function is from quantity_ext.c
    In the evolve procedure:

    * shallow_water/shallow_water_domain.py -> Domain.evolve
        
        * abstruct_2d_finite_volumes/generic_domain -> Domain.evolve
            
            * shallow_water/shallow_water_domain -> Domain.update_other_quantities
                
                * abstruct_2d_finite_volumes/quantity.py -> \
                            Quantity.extrapolate_second_order_and_limit_by_vertex

                    * quantity_ext.c -> extrapolate_second_order_and_limit_by_vertex

                        * _compute_gradients
                        * _extrapolate_from_gradient
                        * _limit_vertices_by_all_neighbours

                        
"""

def extrapolate_second_order_and_limit_by_vertex(domain, Q, step=3):
    
    mod = SourceModule("""
#ifndef USING_UNION
    __global__ void _limit_vertices_by_all_neighbours(
#else
    __device__ int _limit_vertices_by_all_neighbours(
#endif
            //int N, 
            double beta,
            double* centroid_values,
            double* vertex_values,
            double* edge_values,
            long*   neighbours,
            double* x_gradient,
            double* y_gradient) 
    {
        const int k = 
            threadIdx.x+threadIdx.y*blockDim.x+
            (blockIdx.x+blockIdx.y*gridDim.x)*blockDim.x*blockDim.y;
    
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
    
            dq = vertex_values[k3+i] - qc;    
                    //Delta between vertex and centroid values
            dqa[i] = dq;                      
                    //Save dq for use in updating vertex values
    
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
    }
    
#ifndef USING_UNION
    __global__ void _extrapolate_from_gradient(
#else
    __device__ int _extrapolate_from_gradient(
#endif
            //int N,
            double* centroids,
            double* centroid_values,
            double* vertex_coordinates,
            double* vertex_values,
            double* edge_values,
            double* a,
            double* b) 
    {
    
        const int k = 
            threadIdx.x+threadIdx.y*blockDim.x+
            (blockIdx.x+blockIdx.y*gridDim.x)*blockDim.x*blockDim.y;
    
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
    }
    
#ifndef USING_UNION 
    __global__ void _compute_gradients(
#else   
    __device__ int _compute_gradients(
#endif
            //int N,
            double* centroids,
            double* centroid_values,
            long* number_of_boundaries,
            long* surrogate_neighbours,
            double* a,
            double* b)
    {
        const int k = 
            threadIdx.x+threadIdx.y*blockDim.x+
            (blockIdx.x+blockIdx.y*gridDim.x)*blockDim.x*blockDim.y;
        
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
    
    
            if (k0 == k1 || k1 == k2) 
                //return -1;
                return;
    
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
    
            a[k] = ((y2-y0)*(q1-q0) - (y1-y0)*(q2-q0)) /det;
            //a[k] /= det;
    
            b[k] = ((x1-x0)*(q2-q0) - (x2-x0)*(q1-q0)) /det;
            //b[k] /= det;
    
        } else if (number_of_boundaries[k] == 2) {
            // One true neighbour
    
            // Get index of the one neighbour
            i=0; k0 = k;
            while (i<3 && k0==k) {
                k0 = surrogate_neighbours[3*k + i];
                i++;
            }
            if (k0 == k) 
                //return -1;
                return;
    
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
    }
    
    
    
#ifdef USING_UNION    
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
#endif
    """)
       
    N = Q.centroid_values.shape[0]
    
        
    #extra_func = mod.get_function("extrapolate_second_order_and_limit_by_vertex")
    #extra_func( 
    #        cuda.In(beta_gpu ),
    #        cuda.In( domain.centroid_coordinates ),
    #        cuda.In( domain.vertex_coordinates ),
    #        cuda.In( domain.number_of_boundaries ),
    #        cuda.In( domain.surrogate_neighbours ),
    #        cuda.In( domain.neighbours ),

    #        cuda.InOut( Q.centroid_values ),
    #        cuda.InOut( Q.vertex_values ),
    #        cuda.InOut( Q.edge_values ),
    #        cuda.InOut( Q.x_gradient ),
    #        cuda.InOut( Q.y_gradient ),
    #        block = ( W1, 1, 1),
    #        grid = ( (N + W1 -1 ) / W1, 1) )
    
    extra_func = mod.get_function("_compute_gradients")
    W1 = 32
    #W1 = extra_func.max_threads_per_block
    #print W1
    extra_func( 
            cuda.In( domain.centroid_coordinates ),
            cuda.InOut( Q.centroid_values ),
            cuda.In( domain.number_of_boundaries ),
            cuda.In( domain.surrogate_neighbours ),
            cuda.InOut( Q.x_gradient ),
            cuda.InOut( Q.y_gradient ),
            block = ( W1, 1, 1),
            grid = ( (N + W1 -1 ) / W1, 1) )
    
    if step > 1:
        extra_func = mod.get_function("_extrapolate_from_gradient")
        #W1 = extra_func.max_threads_per_block
        #print W1
        extra_func( 
                cuda.In( domain.centroid_coordinates ),
                cuda.InOut( Q.centroid_values ),
                cuda.In( domain.vertex_coordinates ),
                cuda.InOut( Q.vertex_values ),
                cuda.InOut( Q.edge_values ),
                cuda.InOut( Q.x_gradient ),
                cuda.InOut( Q.y_gradient ),
                block = ( W1, 1, 1),
                grid = ( (N + W1 -1 ) / W1, 1) )

    if step > 2:
        extra_func = mod.get_function("_limit_vertices_by_all_neighbours")
        #W1 = extra_func.max_threads_per_block
        #print W1
        extra_func( 
                numpy.float64(Q.beta),
                cuda.InOut( Q.centroid_values ),
                cuda.InOut( Q.vertex_values ),
                cuda.InOut( Q.edge_values ),
                cuda.In( domain.neighbours ),
                cuda.InOut( Q.x_gradient ),
                cuda.InOut( Q.y_gradient ),
                block = ( W1, 1, 1),
                grid = ( (N + W1 -1 ) / W1, 1) )



def extrapolate_second_order_and_limit_by_vertex_single(domain, Q):
    
    mod = SourceModule("""
    __global__ void extrapolate_second_order_and_limit_by_vertex_single(
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
    
        const int k = threadIdx.x + blockIdx.x * blockDim.x;
        int i, k0, k1, k2, k3, k6;
        double x0, x1, x2, y0, y1, y2, q0, q1, q2;
        double det;
        double x, y;
    
        long n;
        double qmin, qmax, qn, qc;
        double dq, dqa[3], phi, r;
    
    
        if (domain_number_of_boundaries[k] < 2) {
            // Two or three true neighbours
    
            // Get indices of neighbours (or self when used as surrogate)
            // k0, k1, k2 = domain_surrogate_neighbours[k,:]
    
            k0 = domain_surrogate_neighbours[3*k + 0];
            k1 = domain_surrogate_neighbours[3*k + 1];
            k2 = domain_surrogate_neighbours[3*k + 2];
    
    
            if (k0 == k1 || k1 == k2) 
                return;
    
            // Get data
            q0 = quantity_centroid_values[k0];
            q1 = quantity_centroid_values[k1];
            q2 = quantity_centroid_values[k2];
    
            x0 = domain_centroid_coordinates[k0*2]; y0 = domain_centroid_coordinates[k0*2+1];
            x1 = domain_centroid_coordinates[k1*2]; y1 = domain_centroid_coordinates[k1*2+1];
            x2 = domain_centroid_coordinates[k2*2]; y2 = domain_centroid_coordinates[k2*2+1];
    
            // Gradient
            //_gradient(x0, y0, x1, y1, x2, y2, q0, q1, q2, &a[k], &b[k]);
    
            det = (y2-y0)*(x1-x0) - (y1-y0)*(x2-x0);
    
            quantity_x_gradient[k] = (y2-y0)*(q1-q0) - (y1-y0)*(q2-q0);
            quantity_x_gradient[k] /= det;
    
            quantity_y_gradient[k] = (x1-x0)*(q2-q0) - (x2-x0)*(q1-q0);
            quantity_y_gradient[k] /= det;
    
        } else if (domain_number_of_boundaries[k] == 2) {
            // One true neighbour
    
            // Get index of the one neighbour
            i=0; k0 = k;
            while (i<3 && k0==k) {
                k0 = domain_surrogate_neighbours[3*k + i];
                i++;
            }
            if (k0 == k) return;
    
            k1 = k; //self
    
            // Get data
            q0 = quantity_centroid_values[k0];
            q1 = quantity_centroid_values[k1];
    
            x0 = domain_centroid_coordinates[k0*2]; y0 = domain_centroid_coordinates[k0*2+1];
            x1 = domain_centroid_coordinates[k1*2]; y1 = domain_centroid_coordinates[k1*2+1];
    
            // Two point gradient
            //_gradient2(x0, y0, x1, y1, q0, q1, &a[k], &b[k]);
    
            det = (x1-x0) * (x1-x0) + (y1-y0)*(y1-y0);
            quantity_x_gradient[k] = (x1-x0)*(q1-q0) / det;
            quantity_y_gradient[k] = (y1-y0)*(q1-q0) / det;
    
        }
    
        k6 = 6*k;
        k3 = 3*k;
        k2 = 2*k;
    
        // Centroid coordinates
        x = domain_centroid_coordinates[k2]; y = domain_centroid_coordinates[k2+1];
    
        // vertex coordinates
        // x0, y0, x1, y1, x2, y2 = X[k,:]
        x0 = domain_vertex_coordinates[k6 + 0];
        y0 = domain_vertex_coordinates[k6 + 1];
        x1 = domain_vertex_coordinates[k6 + 2];
        y1 = domain_vertex_coordinates[k6 + 3];
        x2 = domain_vertex_coordinates[k6 + 4];
        y2 = domain_vertex_coordinates[k6 + 5];
    
        // Extrapolate to Vertices
        quantity_vertex_values[k3+0] = quantity_centroid_values[k] + quantity_x_gradient[k]*(x0-x) + quantity_y_gradient[k]*(y0-y);
        quantity_vertex_values[k3+1] = quantity_centroid_values[k] + quantity_x_gradient[k]*(x1-x) + quantity_y_gradient[k]*(y1-y);
        quantity_vertex_values[k3+2] = quantity_centroid_values[k] + quantity_x_gradient[k]*(x2-x) + quantity_y_gradient[k]*(y2-y);
    
        // Extrapolate to Edges (midpoints)
        quantity_edge_values[k3+0] = 0.5*(quantity_vertex_values[k3 + 1]+quantity_vertex_values[k3 + 2]);
        quantity_edge_values[k3+1] = 0.5*(quantity_vertex_values[k3 + 2]+quantity_vertex_values[k3 + 0]);
        quantity_edge_values[k3+2] = 0.5*(quantity_vertex_values[k3 + 0]+quantity_vertex_values[k3 + 1]);
    
    
        k3 = 3*k;
    
        qc = quantity_centroid_values[k];
        qmin = qc;
        qmax = qc;
    
        for (i=0; i<3; i++) {
            n = domain_neighbours[k3+i];
            if (n >= 0) {
                qn = quantity_centroid_values[n]; //Neighbour's centroid value
    
                qmin = min(qmin, qn);
                qmax = max(qmax, qn);
            }
        }
    
        phi = 1.0;
        for (i=0; i<3; i++) {    
            r = 1.0;
    
            dq = quantity_vertex_values[k3+i] - qc;    //Delta between vertex and centroid values
            dqa[i] = dq;                      //Save dq for use in updating vertex values
    
            if (dq > 0.0) r = (qmax - qc)/dq;
            if (dq < 0.0) r = (qmin - qc)/dq;      
    
    
            phi = min( min(r*beta, 1.0), phi);    
        }
    
        //Update gradient, vertex and edge values using phi limiter
        quantity_x_gradient[k] = quantity_x_gradient[k]*phi;
        quantity_y_gradient[k] = quantity_y_gradient[k]*phi;
    
        quantity_vertex_values[k3+0] = qc + phi*dqa[0];
        quantity_vertex_values[k3+1] = qc + phi*dqa[1];
        quantity_vertex_values[k3+2] = qc + phi*dqa[2];
    
        quantity_edge_values[k3+0] = 0.5*(quantity_vertex_values[k3+1] + quantity_vertex_values[k3+2]);
        quantity_edge_values[k3+1] = 0.5*(quantity_vertex_values[k3+2] + quantity_vertex_values[k3+0]);
        quantity_edge_values[k3+2] = 0.5*(quantity_vertex_values[k3+0] + quantity_vertex_values[k3+1]);
    
    } 
    """)
       
    N = Q.centroid_values.shape[0]
    if (N % 32 == 0):
        W1 = 32
    else:
        raise Exception('N can not be splited')

    beta_gpu = numpy.zeros(1, dtype=numpy.float64)
    beta_gpu[0] = Q.beta

        
    extra_func = mod.get_function("extrapolate_second_order_and_limit_by_vertex_single")

    extra_func( 
            cuda.In(beta_gpu ),
            cuda.In( domain.centroid_coordinates ),
            cuda.In( domain.vertex_coordinates ),
            cuda.In( domain.number_of_boundaries ),
            cuda.In( domain.surrogate_neighbours ),
            cuda.In( domain.neighbours ),

            cuda.InOut( Q.centroid_values ),
            cuda.InOut( Q.vertex_values ),
            cuda.InOut( Q.edge_values ),
            cuda.InOut( Q.x_gradient ),
            cuda.InOut( Q.y_gradient ),
            block = ( W1, 1, 1),
            grid = ( (N + W1 -1 ) / W1, 1) )
  



def distribute_using_vertex_limiter(domain,kind = 'cuda', step=3):
    domain.protect_against_infinitesimal_and_negative_heights()
    for name in ['height', 'xvelocity', 'yvelocity']:
        Q = domain.quantities[name]
        if domain._order_ == 1:
            Q.extrapolate_first_order()
        elif domain._order_ == 2:
            if domain.use_edge_limiter:
                extrapolate_second_order_and_limit_by_edge()
            else:
                if  kind == 'single':
                    print name + "  in single function"
                    extrapolate_second_order_and_limit_by_vertex_single(domain,Q)
                elif kind == 'python':
                    print name + "  in python"
                    compute_gradients(domain, Q)
                    if step > 1:
                        extrapolate_from_gradient(domain, Q)
                    if step > 2:
                        limit_vertices_by_all_neighbours(domain, Q)
                else:
                    print name + "  in cuda"
                    extrapolate_second_order_and_limit_by_vertex(domain, Q, step)
        else:
            raise Exception('Unknown order')



def approx_cmp(a,b):
    if abs(a-b) > abs(a)*pow(10,-6):
        return True
    else:
        return False

def _gradient(x0, y0, x1, y1, x2, y2, q0, q1, q2, a, b):
    det = (y2-y0)*(x1-x0) - (y1-y0)*(x2-x0);

    a = (y2-y0)*(q1-q0) - (y1-y0)*(q2-q0);
    a /= det;

    b = (x1-x0)*(q2-q0) - (x2-x0)*(q1-q0);
    b /= det;

def _gradient2(x0, y0, x1, y1, q0, q1, a, b):
    xx = x1-x0
    yy = y1-y0
    qq = q1-q0
            
    det = xx*xx + yy*yy
    a = xx*qq/det
    b = yy*qq/det


def limit_vertices_by_all_neighbours(domain, Q):
    N = Q.centroid_values.shape[0]
    a = Q.x_gradient
    b = Q.y_gradient
    for k in range(N):
        qc = Q.centroid_values[k]
        qmin = qc
        qmax = qc

        for i in range(3):
            n = domain.neighbours[k][i]
            if n >= 0:
                qn = Q.centroid_values[n]
                qmin = min(qmin, qn)
                qmax = max(qmax, qn)
        phi = 1.0
        dqa = numpy.zeros(3, dtype=numpy.float64)
        for i in range(3):
            r = 1.0
            dq = Q.vertex_values[k][i] - qc
            dqa[i] = dq

            if dq > 0.0 :
                r = (qmax - qc) / dq
            if dq < 0.0:
                r = (qmin - qc) / dq

            phi = min ( min(r*Q.beta, 1.0), phi)

        Q.x_gradient[k] *= phi
        Q.y_gradient[k] *= phi

        Q.vertex_values[k][0] = qc + phi * dqa[0]
        Q.vertex_values[k][1] = qc + phi * dqa[1]
        Q.vertex_values[k][2] = qc + phi * dqa[2]

        Q.edge_values[k][0] = 0.5*(Q.vertex_values[k][1]+Q.vertex_values[k][2])
        Q.edge_values[k][1] = 0.5*(Q.vertex_values[k][2]+Q.vertex_values[k][0])
        Q.edge_values[k][2] = 0.5*(Q.vertex_values[k][0]+Q.vertex_values[k][1])

def extrapolate_from_gradient(domain, Q):
    N = Q.centroid_values.shape[0]
    a = Q.x_gradient
    b = Q.y_gradient
    for k in range(N):
        x = domain.centroid_coordinates[k][0]
        y = domain.centroid_coordinates[k][1]
        
        x0 = domain.vertex_coordinates[k*3][0]
        y0 = domain.vertex_coordinates[k*3][1]
        x1 = domain.vertex_coordinates[k*3+1][0]
        y1 = domain.vertex_coordinates[k*3+1][1]
        x2 = domain.vertex_coordinates[k*3+2][0]
        y2 = domain.vertex_coordinates[k*3+2][1]

        Q.vertex_values[k][0] = \
            Q.centroid_values[k] +a[k]*(x0-x)+b[k]*(y0-y)
        Q.vertex_values[k][1] = \
            Q.centroid_values[k] +a[k]*(x1-x)+b[k]*(y1-y)
        Q.vertex_values[k][2] = \
            Q.centroid_values[k] +a[k]*(x2-x)+b[k]*(y2-y)

        Q.edge_values[k][0] = \
                0.5*(Q.vertex_values[k][1] + Q.vertex_values[k][2])
        Q.edge_values[k][1] = \
                0.5*(Q.vertex_values[k][2] + Q.vertex_values[k][0])
        Q.edge_values[k][2] = \
                0.5*(Q.vertex_values[k][0] + Q.vertex_values[k][1])

def compute_gradients(domain, Q):
    N = Q.centroid_values.shape[0]
    a = Q.x_gradient
    b = Q.y_gradient
    for k in range(N):
        if domain.number_of_boundaries[k] < 2:
            k0 = domain.surrogate_neighbours[k][0];
            k1 = domain.surrogate_neighbours[k][1];
            k2 = domain.surrogate_neighbours[k][2];
            
            assert not (k0 == k1 and k1 == k2)

            q0 = Q.centroid_values[k0]
            q1 = Q.centroid_values[k1]
            q2 = Q.centroid_values[k2]
                
            x0 = domain.centroid_coordinates[k0][0]
            y0 = domain.centroid_coordinates[k0][1]
            x1 = domain.centroid_coordinates[k1][0]
            y1 = domain.centroid_coordinates[k1][1]
            x2 = domain.centroid_coordinates[k2][0]
            y2 = domain.centroid_coordinates[k2][1]

            _gradient(x0, y0, x1, y1, x2, y2, q0, q1, q2, a[k], b[k])

        elif domain.number_of_boundaries[k] == 2:
            for i in range(3):
                k0 = domain.surrogate_neighbours[k][i]
                if k0 != k:
                    break
            if k0 == k:
                return -1

            k1 = k
            q0 = Q.centroid_values[k0]
            q1 = Q.centroid_values[k1]

            x0 = domain.centroid_coordinates[k0][0]
            y0 = domain.centroid_coordinates[k0][1]
            x1 = domain.centroid_coordinates[k1][0]
            y1 = domain.centroid_coordinates[k1][1]

            _gradient2(x0, y0, x1, y1, q0, q1, a[k], b[k])

if __name__ == '__main__':
    from anuga_cuda import *

    testing_gpu_domain= True
    testing_2 = True
    testing_3 = False

    domain1 = generate_merimbula_domain()
    domain2 = generate_merimbula_domain(True)

    print "~~~~~~~ domain 1 ~~~~~~~"
    domain1.protect_against_infinitesimal_and_negative_heights()
    for name in domain1.conserved_quantities:
        Q = domain1.quantities[name]
        Q.extrapolate_second_order_and_limit_by_vertex()

    import pycuda.driver as cuda
    import pycuda.autoinit
    from pycuda.compiler import SourceModule
    import numpy

    step = 3

    print "~~~~~~~ domain 2 ~~~~~~~"
    if testing_gpu_domain:
        N = domain2.number_of_elements
        W1 = 32
        W2 = 1
        W3 = 1
        for name in ['height', 'xvelocity', 'yvelocity']:
            Q = domain2.quantities[name]
            domain2.extrapolate_second_order_and_limit_by_vertex_func(
                numpy.int32( N ),
                numpy.float64( Q.beta ),
                cuda.In( domain2.centroid_coordinates ),
                cuda.In( domain2.vertex_coordinates ),
                cuda.In( domain2.number_of_boundaries ),
                cuda.In( domain2.surrogate_neighbours ),
                cuda.In( domain2.neighbours ),

                cuda.InOut( Q.centroid_values ),
                cuda.InOut( Q.vertex_values ),
                cuda.InOut( Q.edge_values ),
                cuda.InOut( Q.x_gradient ),
                cuda.InOut( Q.y_gradient ),
                block = (W1, W2, W3),
                grid = (( N +W1-1)/W1, 1)
                )
            
    else:
        distribute_using_vertex_limiter(domain2, 'cuda', step)

    if testing_2:
        print "\n~~~~~~~ compare 1 2~~~~~~~"
        for name in ['height', 'xvelocity', 'yvelocity']:
            Q1 = domain1.quantities[name]
            Q2 = domain2.quantities[name]
            print name
            
            cnt_cv = 0
            cnt_vv = 0
            cnt_ev = 0
            cnt_xg = 0
            cnt_yg = 0
            for i in range(Q1.centroid_values.shape[0]):
                if Q1.centroid_values[i] != Q2.centroid_values[i]:
                    cnt_cv += 1
                    if cnt_cv < 10:
                        print "cv %d, %lf, %lf" % \
                                (i, Q1.centroid_values[i], Q2.centroid_values[i])
                #if (Q1.vertex_values[i] != Q2.vertex_values[i]).any():
                if approx_cmp(Q1.vertex_values[i][0],Q2.vertex_values[i][0]) or \
                         approx_cmp(Q1.vertex_values[i][1],Q2.vertex_values[i][1]) or \
                         approx_cmp(Q1.vertex_values[i][2],Q2.vertex_values[i][2]) :
                    cnt_vv += 1
                    if cnt_vv < 10:
                        print "vv %d, %s, %s" % \
                                (i, Q1.vertex_values[i], Q2.vertex_values[i])
                #if (Q1.edge_values[i] != Q2.edge_values[i]).any():
                if approx_cmp(Q1.edge_values[i][0] , Q2.edge_values[i][0]) or\
                        approx_cmp(Q1.edge_values[i][1] , Q2.edge_values[i][1]) or\
                        approx_cmp(Q1.edge_values[i][2] , Q2.edge_values[i][2]) :
                    cnt_ev += 1
                    if cnt_ev < 10:
                        print "ev %d, %s, %s" % \
                                (i, Q1.edge_values[i], Q2.edge_values[i])
                #if Q1.x_gradient[i] != Q2.x_gradient[i]:
                if approx_cmp(Q1.x_gradient[i], Q2.x_gradient[i]) :
                    cnt_xg += 1
                    if cnt_xg < 10:
                        print "xg %d, %lf, %lf" % \
                                (i, Q1.x_gradient[i], Q2.x_gradient[i])
                #if Q1.y_gradient[i] != Q2.y_gradient[i]:
                if approx_cmp(Q1.y_gradient[i], Q2.y_gradient[i]) :
                    cnt_yg += 1
                    if cnt_yg < 10:
                        print "yg %d, %lf, %lf" % \
                                (i, Q1.y_gradient[i], Q2.y_gradient[i])

            print "~~ # of diff %d, %d, %d, %d, %d" % \
                    (cnt_cv, cnt_vv, cnt_ev, cnt_xg, cnt_yg)


    if testing_3:
        domain3 = generate_merimbula_domain()
        print "~~~~~~~ domain 3 ~~~~~~~"
        #distribute_using_vertex_limiter(domain3, True, 0)
        distribute_using_vertex_limiter(domain3, 'python', step)

        print "\n~~~~~~~ compare 1 3~~~~~~~"
        for name in ['height', 'xvelocity', 'yvelocity']:
            Q1 = domain2.quantities[name]
            Q2 = domain3.quantities[name]
            print name
            
            cnt_cv = 0
            cnt_vv = 0
            cnt_ev = 0
            cnt_xg = 0
            cnt_yg = 0
            for i in range(Q1.centroid_values.shape[0]):
                if Q1.centroid_values[i] != Q2.centroid_values[i]:
                    cnt_cv += 1
                    if cnt_cv < 3:
                        print "cv %d, %lf, %lf" % \
                                (i, Q1.centroid_values[i], Q2.centroid_values[i])
                #if (Q1.vertex_values[i] != Q2.vertex_values[i]).any():
                if approx_cmp(Q1.vertex_values[i][0],Q2.vertex_values[i][0]) or \
                         approx_cmp(Q1.vertex_values[i][1],Q2.vertex_values[i][1]) or \
                         approx_cmp(Q1.vertex_values[i][2],Q2.vertex_values[i][2]) :
                    cnt_vv += 1
                    if cnt_vv < 3:
                        print "vv %d, %s, %s" % \
                                (i, Q1.vertex_values[i], Q2.vertex_values[i])
                #if (Q1.edge_values[i] != Q2.edge_values[i]).any():
                if approx_cmp(Q1.edge_values[i][0] , Q2.edge_values[i][0]) or\
                        approx_cmp(Q1.edge_values[i][1] , Q2.edge_values[i][1]) or\
                        approx_cmp(Q1.edge_values[i][2] , Q2.edge_values[i][2]) :
                    cnt_ev += 1
                    if cnt_ev < 3:
                        print "ev %d, %s, %s" % \
                                (i, Q1.edge_values[i], Q2.edge_values[i])
                #if Q1.x_gradient[i] != Q2.x_gradient[i]:
                if approx_cmp(Q1.x_gradient[i], Q2.x_gradient[i]) :
                    cnt_xg += 1
                    if cnt_xg < 3:
                        print "xg %d, %lf, %lf" % \
                                (i, Q1.x_gradient[i], Q2.x_gradient[i])
                #if Q1.y_gradient[i] != Q2.y_gradient[i]:
                if approx_cmp(Q1.y_gradient[i], Q2.y_gradient[i]) :
                    cnt_yg += 1
                    if cnt_yg < 10:
                        print "yg %d, %lf, %lf" % \
                                (i, Q1.y_gradient[i], Q2.y_gradient[i])

            print "~~ # of diff %d, %d, %d, %d, %d" % \
                    (cnt_cv, cnt_vv, cnt_ev, cnt_xg, cnt_yg)


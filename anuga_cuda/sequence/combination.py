#!/usr/bin/env python
    
import pycuda.driver as cuda
import pycuda.autoinit
from pycuda.compiler import SourceModule
import numpy

from anuga_cuda.merimbula_data.generate_domain import domain_create    
 
domain1 = domain_create()
domain2 = domain_create()




def mem_all_cpy(a):
    a_gpu = cuda.mem_alloc(a.nbytes)
    cuda.memcpy_htod(a_gpu, a)
    return a_gpu


def compute_fluxes_and_grivity(domain, parallelFlag = 1):

    
    print " in the compute_fluxes_central_structure cuda function"
    N = domain.number_of_elements

    

    mod = SourceModule("""
    #include "math.h"
    #define Dg 0
    #define Depsilon 1
    #define DH0 2
    #define Doptimise_dry_cells 3
    #define Dcall 5
    

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
            //report_python_error(AT, "Discontinuous Elevation");
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

    __device__ void spe_bubble_sort(int* _list , long* neighbours, int k)
    {
        int temp;
        if ( neighbours[ _list[2] ]>=0 and neighbours[ _list[2] ] < k and (neighbours[ _list[1] ]<0 or neighbours[ _list[2] ] <neighbours [_list[1] ]) )
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


    // The parallel CUDA compute fluex function
    __global__ void compute_fluxes_central_structure_CUDA(
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
            double * max_speed_array)
    {
        const int k = threadIdx.x + (blockIdx.x )*blockDim.x;
    
    
        double max_speed, length, inv_area, zl, zr;
    
        double h0 = elements[DH0] * elements[DH0]; // This ensures a good balance when h approaches H0.
    
        double limiting_threshold = 10 * elements[DH0]; // Avoid applying limiter below this
        
        int i, m, n, j;
        int ki, nm = 0, ki2; // Index shorthands
    
        double ql[3], qr[3], edgeflux[3]; // Work array for summing up fluxes
        

        int b[3]={0,1,2};

        spe_bubble_sort( b, neighbours+k*3, k);

        for (j = 0; j < 3; j++) {
            i = b[j];

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
    
    
        
        inv_area = 1.0 / areas[k];
        stage_explicit_update[k] *= inv_area;
        xmom_explicit_update[k] *= inv_area;
        ymom_explicit_update[k] *= inv_area;
    
        max_speed_array[k] =  max_speed;
    }
    
    
     __global__ void gravity_wb(
                double * stage_vertex_values, 
                double * stage_edge_values, 
                double * stage_centroid_values, 
                double * bed_edge_values, 
                double * bed_centroid_values, 
                double * vertex_coordinates, 
                double * xmom_explicit_update, 
                double * ymom_explicit_update, 
                double * normals, 
                double * areas, 
                double * edgelengths,
                double * g,
                double * temp_array)
        {
            const int k = threadIdx.x + (blockIdx.x )*blockDim.x;

            int i;
            
            double w0, w1, w2, 
                    x0, y0, x1, y1, x2, y2,
                    avg_h;
                    
            double wx, wy, det,
                hh[3];
            double sidex, sidey, area, n0, n1, fact;

    
            w0 = stage_vertex_values[3*k + 0];
            w1 = stage_vertex_values[3*k + 1];
            w2 = stage_vertex_values[3*k + 2];

            x0 = vertex_coordinates[k*6 + 0];
            y0 = vertex_coordinates[k*6 + 1];
            x1 = vertex_coordinates[k*6 + 2];
            y1 = vertex_coordinates[k*6 + 3];
            x2 = vertex_coordinates[k*6 + 4];
            y2 = vertex_coordinates[k*6 + 5];

        
            //_gradient(x0, y0, x1, y1, x2, y2, w0, w1, w2, &wx, &wy);

            det = (y2 - y0)*(x1 - x0) - (y1 - y0)*(x2 - x0);
    
            wx = (y2 -y0)*(w1 - w0) - (y1 - y0)*(w2 -w0);
            wx /= det;
    
            wy = (x1 - x0)*(w2 - w0) - (x2 - x0)*(w1 -w0);
            wy /= det;
    

            avg_h = stage_centroid_values[k] - bed_centroid_values[k];

            xmom_explicit_update[k] += -g[0] * wx * avg_h;
            ymom_explicit_update[k] += -g[0] * wy * avg_h;
    

            hh[0] = stage_edge_values[k*3] - bed_edge_values[k*3];
            hh[1] = stage_edge_values[k*3+1] - bed_edge_values[k*3+1];
            hh[2] = stage_edge_values[k*3+2] - bed_edge_values[k*3+2];

            sidex = 0.0;
            sidey = 0.0;
            area = areas[k];

            for ( i = 0 ; i < 3 ; i++ )
            {
                n0 = normals[k*6 + 2*i];
                n1 = normals[k*6 + 2*i + 1];
                
                fact =  -0.5 * g[0] * hh[i] * hh[i] * edgelengths[k*3 + i];
                
                sidex += fact*n0;
                sidey += fact*n1;
                    
                temp_array[k*6 + 2*i] = fact*n0;
                temp_array[k*6 + 2*i + 1] = fact*n1;
            }
            
            //xmom_explicit_update[k] += -sidex / area;
            //ymom_explicit_update[k] += -sidey / area;
                
        }


    
        """)
    
    elements = numpy.random.randn(5)
    elements = elements.astype(numpy.float64)
    elements[0] = domain.g
    elements[1] = domain.epsilon
    elements[2] = domain.H0
    elements[3] = domain.optimise_dry_cells

    
    temp_array = numpy.zeros((domain1.number_of_elements, 3,2), dtype=numpy.float64)

    timestep_array = numpy.zeros( (domain.number_of_elements) , dtype=numpy.float64)

    # Necessary data initialization
    for i in range(N):
        timestep_array[i] = domain.evolve_max_timestep
        domain.quantities['stage'].explicit_update[i] = 0.0
        domain.quantities['xmomentum'].explicit_update[i]=0.0
        domain.quantities['ymomentum'].explicit_update[i] =0.0
    

    
    if (N % 256 == 0):
        W1 = 16
        W2 = 16
    elif (N % 32 ==0):
        W1 = 32
        W2 = 1
    else:
        raise Exception('N can not be splited')

    if 1:
        compute_fluxes_central_function = mod.get_function('compute_fluxes_central_structure_CUDA')

         # facilities vertors
        timestep_array_gpu = mem_all_cpy(timestep_array )
        temp_array_gpu = mem_all_cpy(temp_array)
        elements_gpu = mem_all_cpy(elements)

        # neighbours vertors
        neighbours_gpu = mem_all_cpy(domain.neighbours)
        neighbour_edges_gpu = mem_all_cpy(domain.neighbour_edges)

        # 
        normals_gpu = mem_all_cpy(domain.normals)
        edgelengths_gpu = mem_all_cpy(domain.edgelengths)
        radii_gpu = mem_all_cpy(domain.radii)
        areas_gpu = mem_all_cpy(domain.areas)
        tri_full_flag_gpu = mem_all_cpy(domain.tri_full_flag)
        vertex_coordinates_gpu = mem_all_cpy(domain.vertex_coordinates)

        # edge_values vertors
        stage_edge_values_gpu = mem_all_cpy(domain.quantities['stage'].edge_values)
        xmom_edge_values_gpu = mem_all_cpy(domain.quantities['xmomentum'].edge_values)
        ymom_edge_values_gpu = mem_all_cpy(domain.quantities['ymomentum'].edge_values)
        bed_edge_values_gpu = mem_all_cpy(domain.quantities['elevation'].edge_values)

        # boundary_values vertors
        stage_boundary_values_gpu = mem_all_cpy(domain.quantities['stage'].boundary_values)
        xmom_boundary_values_gpu = mem_all_cpy(domain.quantities['xmomentum'].boundary_values)
        ymom_boundary_values_gpu = mem_all_cpy(domain.quantities['ymomentum'].boundary_values)

        # vertex_values
        stage_vertex_values_gpu = mem_all_cpy(domain.quantities['stage'].vertex_values)

        # centroid_values
        stage_centroid_values_gpu = mem_all_cpy(domain.quantities['stage'].centroid_values)
        bed_centroid_values_gpu = mem_all_cpy(domain.quantities['elevation'].centroid_values)

        # explicit_update vertors
        stage_explicit_update_gpu = mem_all_cpy(domain.quantities['stage'].explicit_update)
        xmom_explicit_update_gpu = mem_all_cpy(domain.quantities['xmomentum'].explicit_update)
        ymom_explicit_update_gpu = mem_all_cpy(domain.quantities['ymomentum'].explicit_update)

        # max_speed vertors
        max_speed_gpu = mem_all_cpy(domain.max_speed)


        compute_fluxes_central_function(
            elements_gpu, 
            timestep_array_gpu, 
            neighbours_gpu, 
            neighbour_edges_gpu,
            normals_gpu, 
            edgelengths_gpu, 
            radii_gpu, 
            areas_gpu, 
            tri_full_flag_gpu,
            stage_edge_values_gpu, 
            xmom_edge_values_gpu, 
            ymom_edge_values_gpu, 
            bed_edge_values_gpu, 
            stage_boundary_values_gpu, 
            xmom_boundary_values_gpu, 
            ymom_boundary_values_gpu, 
            stage_explicit_update_gpu, 
            xmom_explicit_update_gpu, 
            ymom_explicit_update_gpu, 
            max_speed_gpu, 
            block = ( W1, 1, 1),
            grid = ( (N + W1 - 1)/W1, 1) 
            )

#        cuda.memcpy_dtoh(domain.quantities['stage'].explicit_update, stage_explicit_update_gpu)
#        cuda.memcpy_dtoh(domain.quantities['xmomentum'].explicit_update, xmom_explicit_update_gpu)
#        cuda.memcpy_dtoh(domain.quantities['ymomentum'].explicit_update, ymom_explicit_update_gpu)

        gravity_wb_func = mod.get_function("gravity_wb")

        
        gravity_wb_func( 
            stage_vertex_values_gpu, 
            stage_edge_values_gpu, 
            stage_centroid_values_gpu, 
  			bed_edge_values_gpu, 
            bed_centroid_values_gpu, 
            vertex_coordinates_gpu, 
            xmom_explicit_update_gpu,
            ymom_explicit_update_gpu,
            normals_gpu,
            areas_gpu,
            edgelengths_gpu,
            elements_gpu,
            temp_array_gpu,
            block = ( W1, 1, 1),
            grid = ( (N + W1 -1 ) / W1, 1) )
        
        
        cuda.memcpy_dtoh(temp_array, temp_array_gpu)
        cuda.memcpy_dtoh(domain.quantities['stage'].explicit_update, stage_explicit_update_gpu)
        cuda.memcpy_dtoh(domain.quantities['xmomentum'].explicit_update, xmom_explicit_update_gpu)
        cuda.memcpy_dtoh(domain.quantities['ymomentum'].explicit_update, ymom_explicit_update_gpu)
        cuda.memcpy_dtoh(timestep_array, timestep_array_gpu)

        b = numpy.argsort(timestep_array)
        domain.flux_timestep = timestep_array[b[0]] 

        for i in range(domain.number_of_elements):
            area = domain.areas[i]
            temp = temp_array[i][0][0] + temp_array[i][1][0] + temp_array[i][2][0]
            domain.quantities['xmomentum'].explicit_update[i] += -temp / area 
            temp = temp_array[i][0][1] + temp_array[i][1][1] + temp_array[i][2][1]
            domain.quantities['ymomentum'].explicit_update[i] += -temp / area



def approx_cmp(a,b):
    if abs(a-b) > abs(a)*pow(10,-6):
        return True
    else:
        return False



if __name__ == '__main__':
    
    compute_fluxes_and_grivity(domain2)       
    
    from shallow_water_ext import compute_fluxes_ext_central_structure
    from shallow_water_ext import gravity_wb as gravity_wb_c

    domain1.flux_timestep = compute_fluxes_ext_central_structure(domain1)
    gravity_wb_c(domain1)


    print "\n~~~~~~~~~~~~~ domain 2 ~~~~~~~~~~~~"
    print "******* flux_timestep : %lf %lf %d" % \
            (domain1.flux_timestep, domain2.flux_timestep, \
                domain1.flux_timestep == domain2.flux_timestep)
       

    stage_h1= domain1.quantities['stage']
    xmom_h1 = domain1.quantities['xmomentum']
    ymom_h1 = domain1.quantities['ymomentum']

    stage_h2= domain2.quantities['stage']
    xmom_h2 = domain2.quantities['xmomentum']
    ymom_h2 = domain2.quantities['ymomentum']

    counter_s_vv = 0
    counter_x_cv = 0
    counter_y_cv = 0
    counter_x_vv = 0
    counter_y_vv = 0
    
#    print "------------- stage_vertex_values ---------"
#    print stage_h1.vertex_values
#    print stage_h2.vertex_values
#    print "------------- xmom_centroid_values ---------"
#    print xmom_h1.centroid_values
#    print xmom_h2.centroid_values
#    print "------------- xmom_vertex_values ---------"
#    print xmom_h1.vertex_values
#    print xmom_h2.vertex_values
#    print "------------- ymom_centroid_values ---------"
#    print ymom_h1.centroid_values
#    print ymom_h2.centroid_values
#    print "------------- ymom_vertex_values ---------"
#    print ymom_h1.vertex_values
#    print ymom_h2.vertex_values
    svv1 = stage_h1.vertex_values
    svv2 = stage_h2.vertex_values

    cnt_seu = 0
    cnt_xeu = 0
    cnt_yeu = 0
    for i in range(domain1.number_of_elements):
        if approx_cmp(svv1[i][0], svv2[i][0]) or \
                approx_cmp(svv1[i][1], svv2[i][1]) or \
                approx_cmp(svv1[i][2], svv2[i][2]):
            counter_s_vv += 1
            if counter_s_vv < 10:
                print i, stage_h1.vertex_values[i], stage_h2.vertex_values[i],\
                        (stage_h1.vertex_values[i] == stage_h2.vertex_values[i])

        if xmom_h1.centroid_values[i] != xmom_h2.centroid_values[i]:
            counter_x_cv += 1

        if (xmom_h1.vertex_values[i] != xmom_h2.vertex_values[i]).any():
            counter_x_vv += 1

        if ymom_h1.centroid_values[i] != ymom_h2.centroid_values[i]:
            counter_y_cv += 1

        if (ymom_h1.vertex_values[i] != ymom_h2.vertex_values[i]).any():
            counter_y_vv += 1

        if approx_cmp( stage_h1.explicit_update[i], stage_h2.explicit_update[i]):
            cnt_seu += 1

        if approx_cmp( xmom_h1.explicit_update[i], xmom_h2.explicit_update[i] ):
            cnt_xeu += 1

        if approx_cmp( ymom_h2.explicit_update[i], ymom_h2.explicit_update[i] ):
            cnt_yeu += 1

    print "*** # diff svv : %d" % counter_s_vv
    print "*** # diff xcv : %d" % counter_x_cv
    print "*** # diff xvv : %d" % counter_x_vv
    print "*** # diff ycv : %d" % counter_y_cv
    print "*** # diff yvv : %d" % counter_y_vv
    print "*** # diff seu : %d" % cnt_seu
    print "*** # diff xeu : %d" % cnt_xeu
    print "*** # diff yeu : %d" % cnt_yeu


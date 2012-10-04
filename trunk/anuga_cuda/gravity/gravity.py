#!/usr/bin/env python


def gravity(domain):
    
    W = 32

    print domain.quantities['xmomentum'].explicit_update
    print domain.quantities['ymomentum'].explicit_update
    
    mod = SourceModule("""
        __global__ void gravity(double *bed_vertex_values, double *stage_centroid_values, double *bed_centroid_values, double *vertex_coordinates, double * xmom_explicit_update, double *ymom_explicit_update, float g)
        {
            int i = threadIdx.x + threadIdx.y + blockIdx.x * 32*32;
            int k3 = 3*i,
                k6 = 6*i;
            double 
                z0 = bed_vertex_values[k3],
                z1 = bed_vertex_values[k3 + 1],
                z2 = bed_vertex_values[k3 + 2],
    
                avg_h = stage_centroid_values[i] - bed_centroid_values[i],
    
                x0 = vertex_coordinates[k6],
                y0 = vertex_coordinates[k6 + 1],
                x1 = vertex_coordinates[k6 + 2],
                y1 = vertex_coordinates[k6 + 3],
                x2 = vertex_coordinates[k6 + 4],
                y2 = vertex_coordinates[k6 + 5],
                
                zx, zy, det;
    
            det = (y2 - y0)*(x1 - x0) - (y1 - y0)*(x2 - x0);
    
            zx = (y2 -y0)*(z1 - z0) - (y1 - y0)*(z2 -z0);
            zx /= det;
    
            zy = (x1 - x0)*(z2 - z0) - (x2 - x0)*(z1 -z0);
            zy /= det;
    
            xmom_explicit_update[i] += -g * zx * avg_h;
            ymom_explicit_update[i] += -g * zy * avg_h;
        }
        """)


    
    gravity_func = mod.get_function("gravity")
    gravity_func( \
            cuda.In( domain.quantities['elevation'].vertex_values), \
            cuda.In( domain.quantities['stage'].centroid_values), \
            cuda.In( domain.quantities['elevation'].centroid_values), \
            cuda.In( domain.vertex_coordinates), \
            cuda.InOut( domain.quantities['xmomentum'].explicit_update),\
            cuda.InOut( domain.quantities['ymomentum'].explicit_update),\
            cuda.In( domain.g),\
            block = ( W, W, 1),\
            grid = ( (N + W*W -1 ) / (W*W), 1) )
    
    print domain.quantities['xmomentum'].explicit_update
    print domain.quantities['ymomentum'].explicit_update
 
def gravity_wb(domain):
    
    W = 32

    print domain.quantities['xmomentum'].explicit_update
    print domain.quantities['ymomentum'].explicit_update
    
    
    mod = SourceModule( """
        __global__ void gravity_wb(double *stage_vertex_values, double *stage_edge_values, double *stage_centroid_values, double *bed_centroid_values, double *vertex_coordinates, double * xmom_explicit_update, double *ymom_explicit_update, double * normals, double *areas, double * edgelength,float g)
        {
            int k = threadIdx.x + threadIdx.y + blockIdx.x * 32*32;
            int k3 = 3*k,
                k6 = 6*k;
                
            int j;
            
            double 
                w0 = stage_vertex_values[k3],
                w1 = stage_vertex_values[k3 + 1],
                w2 = stage_vertex_values[k3 + 2],
    
                avg_h = stage_centroid_values[k] - bed_centroid_values[k],
    
                x0 = vertex_coordinates[k6],
                y0 = vertex_coordinates[k6 + 1],
                x1 = vertex_coordinates[k6 + 2],
                y1 = vertex_coordinates[k6 + 3],
                x2 = vertex_coordinates[k6 + 4],
                y2 = vertex_coordinates[k6 + 5],
                    
                wx, wy, det,
                hh[3],
                sidex, sidey, area, n0, n1, fact;
    
            det = (y2 - y0)*(x1 - x0) - (y1 - y0)*(x2 - x0);
    
            wx = (y2 -y0)*(w1 - w0) - (y1 - y0)*(w2 -w0);
            wx /= det;
    
            wy = (x1 - x0)*(w2 - w0) - (x2 - x0)*(w1 -w0);
            wy /= det;
    
            xmom_explicit_update[k] += -g * wx * avg_h;
            ymom_explicit_update[k] += -g * wy * avg_h;
    
            hh[0] = stage_edge_values[k3] - bed_edge_values[k3];
            hh[1] = stage_edge_values[k3+1] - bed_edge_values[k3+1];
            hh[2] = stage_edge_values[k3+2] - bed_edge_values[k3+2];
            
            sidex = 0.0
            sidey = 0.0
            
            for ( i = 0 ; i < 3 ; i++ )
            {
                n0 = normal[k6 + 2*i];
                n1 = normal[k6 + 2*i + 1];
                
                fact =  -0.5 * g * hh[i] * hh[i] * edgelengths[k3 + i];
                sidex = sidex + fact*n0;
                sidey = sidey + fact*n1;
            }
            
            area = areas[k];
            xmom_explicit_update[k] += -sidex / area;
            ymom_explicit_update[k] += -sidey / area;
                
        }
        """)
 
    
    gravity_wb_func = mod.get_function("gravity_wb")
    gravity_wb_func( \
            cuda.In( domain.quantities['stage'].vertex_values), \
            cuda.In( domain.quantities['stage'].edge_values), \
            cuda.In( domain.quantities['stage'].centroid_values), \
            cuda.In( domain.quantities['elevation'].centroid_values), \
            cuda.In( domain.vertex_coordinates), \
            cuda.In( domain.quantities['xmomentum'].explicit_update),\
            cuda.In( domain.quantities['ymomentum'].explicit_update),\
            cuda.In( domain.normals),\
            cuda.In( domain.areas),\
            cuda.In( domain.edgelength),\
            cuda.In( domain.g),\
            block = ( W, W, 1),\
            grid = ( (N + W*W -1 ) / (W*W), 1) )
    
    print domain.quantities['xmomentum'].explicit_update
    print domain.quantities['ymomentum'].explicit_update



def gravity_old( domain ):
    N = domain.number_of_elements
    W = 32
    
    mod = SourceModule("""
        __global__ void gravity(double *bed_vertex_values, double *stage_centroid_values, double *bed_centroid_values, double *vertex_coordinates, double * xmom_explicit_update, double *ymom_explicit_update, float g)
        {
            int i = threadIdx.x + threadIdx.y + blockIdx.x * 32*32;
            int k3 = 3*i,
                k6 = 6*i;
            double 
                z0 = bed_vertex_values[k3],
                z1 = bed_vertex_values[k3 + 1],
                z2 = bed_vertex_values[k3 + 2],
    
                avg_h = stage_centroid_values[i] - bed_centroid_values[i],
    
                x0 = vertex_coordinates[k6],
                y0 = vertex_coordinates[k6 + 1],
                x1 = vertex_coordinates[k6 + 2],
                y1 = vertex_coordinates[k6 + 3],
                x2 = vertex_coordinates[k6 + 4],
                y2 = vertex_coordinates[k6 + 5],
                
                zx, zy, det;
    
            det = (y2 - y0)*(x1 - x0) - (y1 - y0)*(x2 - x0);
    
            zx = (y2 -y0)*(z1 - z0) - (y1 - y0)*(z2 -z0);
            zx /= det;
    
            zy = (x1 - x0)*(z2 - z0) - (x2 - x0)*(z1 -z0);
            zy /= det;
    
            xmom_explicit_update[i] += -g * zx * avg_h;
            ymom_explicit_update[i] += -g * zy * avg_h;
        }
        """)
    
    bed_vertex_values_gpu = cuda.mem_alloc(domain.quantities['elevation'].vertex_values.nbytes)
    cuda.memcpy_htod(bed_vertex_values_gpu, domain.quantities['elevation'].vertex_values)

    stage_centroid_values_gpu = cuda.mem_alloc(domain.quantities['stage'].centroid_values.nbytes)
    cuda.memcpy_htod(stage_centroid_values_gpu, domain.quantities['stage'].centroid_values)

    bed_centroid_values_gpu = cuda.mem_alloc(domain.quantities['elevation'].centroid_values.nbytes)
    cuda.memcpy_htod(bed_centroid_values_gpu, domain.quantities['elevation'].centroid_values)

    vertex_coordinates_gpu = cuda.mem_alloc(domain.vertex_coordinates.nbytes)
    cuda.memcpy_htod(vertex_coordinates_gpu, domain.vertex_coordinates)

    xmom_explicit_update_gpu = cuda.mem_alloc(domain.quantities['xmomentum'].explicit_update.nbytes)
    cuda.memcpy_htod(xmom_explicit_update_gpu, domain.quantities['xmomentum'].explicit_update)
    
    ymom_explicit_update_gpu = cuda.mem_alloc(domain.quantities['ymomentum'].explicit_update.nbytes)
    cuda.memcpy_htod(ymom_explicit_update_gpu, domain.quantities['ymomentum'].explicit_update)

    g_gpu = cuda.mem_alloc(4)
    cuda.memcpy_htod(g_gpu, domain.g)

    threadInBlock_gpu = cuda.mem_alloc(numpy.float32.nbytes)
    cuda.memcpy_htod(threadInBlock_gpu, W*W)
	

   	
	#print domain.quantities['xmomentum'].explicit_update
	#print domain.quantities['ymomentum'].explicit_update

    gravity_func = mod.get_function("gravity")
    gravity_func( bed_vertex_values_gpu, stage_centroid_values_gpu, bed_centroid_values_gpu, vertex_coordinates_gpu, xmom_explicit_update_gpu, ymom_explicit_update_gpu, g_gpu, block = ( W, W, 1), grid = ( (N + W*W -1 ) / (W*W), 1) )
    

    cuda.memcpy_dtoh(domain.quantities['xmomentum'].explicit_update, xmom_explicit_update_gpu)
    cuda.memcpy_dtoh(domain.quantities['ymomentum'].explicit_update, ymom_explicit_update_gpu)
	
	#print "--------- after cuda copy ---------------"

	#print domain.quantities['xmomentum'].explicit_update

	#print "-----------------------------------------"

	#print domain.quantities['ymomentum'].explicit_update



if __name__ == '__main__':
    from anuga_cuda.merimbula_data.generate_domain import *
    domain=domain_create()
    print domain.quantities['xmomentum'].explicit_update
    print domain.quantities['ymomentum'].explicit_update
    import pycuda.driver as cuda
    import pycuda.autoinit
    from pycuda.compiler import SourceModule
    import numpy
    gravity_old(domain)
    print domain.quantities['xmomentum'].explicit_update
    print domain.quantities['ymomentum'].explicit_update

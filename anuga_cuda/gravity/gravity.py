#!/usr/bin/env python

def gravity(domain):
    
    N = domain.number_of_elements
    W = 16

    #print domain.quantities['xmomentum'].explicit_update
    #print domain.quantities['ymomentum'].explicit_update
    
    mod = SourceModule("""
        __global__ void gravity(double *bed_vertex_values, double *stage_centroid_values, double *bed_centroid_values, double *vertex_coordinates, double * xmom_explicit_update, double *ymom_explicit_update, float *g)
        {
            int k = threadIdx.x + blockIdx.x * blockDim.x;
            int k3 = 3*k,
                k6 = 6*k;
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
    
            xmom_explicit_update[i] += -g[0] * zx * avg_h;
            ymom_explicit_update[i] += -g[0] * zy * avg_h;
        }
        """)

    g_gpu = numpy.random.rand(1)
    g_gpu = g_gpu.astype(numpy.float32)
    g_gpu[0] = domain.g

    if N % 32 == 0:
        W = 32
    
    gravity_func = mod.get_function("gravity")
    gravity_func( \
            cuda.In( domain.quantities['elevation'].vertex_values), \
            cuda.In( domain.quantities['stage'].centroid_values), \
            cuda.In( domain.quantities['elevation'].centroid_values), \
            cuda.In( domain.vertex_coordinates), \
            cuda.InOut( domain.quantities['xmomentum'].explicit_update),\
            cuda.InOut( domain.quantities['ymomentum'].explicit_update),\
            cuda.In( g_gpu),\
            block = ( W, 1, 1),\
            grid = ( (N + W -1 ) / W, 1) )
    
	
    print "--------------Gravity-----------------------"
    print domain.quantities['xmomentum'].explicit_update
    print domain.quantities['ymomentum'].explicit_update
 
def gravity_wb(domain):
    N = domain.number_of_elements
    
    import pycuda.driver as cuda                            
    import pycuda.autoinit
    from pycuda.compiler import SourceModule  
    import numpy
    
    mod = SourceModule( """
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
                double g,
                long N
                )
        {
            const int k = threadIdx.x + (blockIdx.x )*blockDim.x;

            int i;
            
            double w0, w1, w2, 
                    x0, y0, x1, y1, x2, y2,
                    avg_h;
                    
            double wx, wy, det,
                hh[3];
            double sidex, sidey, area, n0, n1, fact;

            __shared__ double sh_data[32*6];

            if (k >= N)
                return;
    
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

            xmom_explicit_update[k] += -g * wx * avg_h;
            ymom_explicit_update[k] += -g * wy * avg_h;
    

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
                
                fact =  -0.5 * g * hh[i] * hh[i] * edgelengths[k*3 + i];
                
                //sidex += fact*n0;
                //sidey += fact*n1;
                
                sh_data[threadIdx.x + i*blockDim.x] = fact*n0;
                sh_data[threadIdx.x + (i+3)*blockDim.x] = fact*n1;
            }
            
            //xmom_explicit_update[k] += -sidex / area;
            //ymom_explicit_update[k] += -sidey / area;

            xmom_explicit_update[k] += -(sh_data[threadIdx.x] + sh_data[threadIdx.x + blockDim.x] + sh_data[threadIdx.x+2*blockDim.x]) / area;

            ymom_explicit_update[k] += -(sh_data[threadIdx.x+3*blockDim.x] + sh_data[threadIdx.x + 4*blockDim.x] + sh_data[threadIdx.x+5*blockDim.x]) / area;
                
        }
        """)

    #if (N % 32 == 0):
    #    W1 = 32
    #elif (N % 256 ==0):
    #    W1 = 32
    #else:
    #    raise Exception('N can not be splited')

    gravity_wb_func = mod.get_function("gravity_wb")
    W2 = 1
    W3 = 1

    gravity_wb_func( 
            cuda.In( domain.quantities['stage'].vertex_values), 
            cuda.In( domain.quantities['stage'].edge_values), 
            cuda.In( domain.quantities['stage'].centroid_values), 
			cuda.In( domain.quantities['elevation'].edge_values), 
            cuda.In( domain.quantities['elevation'].centroid_values), 
            cuda.In( domain.vertex_coordinates), 
            cuda.InOut( domain.quantities['xmomentum'].explicit_update),
            cuda.InOut( domain.quantities['ymomentum'].explicit_update),
            cuda.In( domain.normals),
            cuda.In( domain.areas),
            cuda.In( domain.edgelengths),
            numpy.float64( domain.g),
            numpy.uint64( N),
            block = ( W1, W2, W3),
            grid = ( N/(W1*W2*W3), 1) )
    
    #print "----------Gravity_wb------------------------"
    #for i in range(domain.number_of_elements):
    #    area = domain.areas[i]
    #    temp = temp_array[i][0][0] + temp_array[i][1][0] + temp_array[i][2][0]
    #    domain.quantities['xmomentum'].explicit_update[i] += -temp / area 
    #    temp = temp_array[i][0][1] + temp_array[i][1][1] + temp_array[i][2][1]
    #    domain.quantities['ymomentum'].explicit_update[i] += -temp / area
    





def gravity_old( domain ):
    N = domain.number_of_elements
    W = 16
    
    mod = SourceModule("""
        __global__ void gravity(double *bed_vertex_values, double *stage_centroid_values, double *bed_centroid_values, double *vertex_coordinates, double * xmom_explicit_update, double *ymom_explicit_update, float *g)
        {
            int i = threadIdx.x + threadIdx.y + blockIdx.x * 16 * 16;
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
    
            xmom_explicit_update[i] += -g[0] * zx * avg_h;
            ymom_explicit_update[i] += -g[0] * zy * avg_h;
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

	# Since transimiting to gpu needs array
    a = numpy.random.randn(1)
    a = a.astype(numpy.float32)
    a[0] = domain.g
    g_gpu = cuda.mem_alloc(a.nbytes)
    cuda.memcpy_htod(g_gpu, a)

   	
	#print domain.quantities['xmomentum'].explicit_update
	#print domain.quantities['ymomentum'].explicit_update

    gravity_func = mod.get_function("gravity")
    gravity_func( 
            bed_vertex_values_gpu, 
            stage_centroid_values_gpu, 
            bed_centroid_values_gpu, 
            vertex_coordinates_gpu, 
            xmom_explicit_update_gpu, 
            ymom_explicit_update_gpu, 
            g_gpu, 
            block = ( W, W, 1), 
            grid = ( (N + W*W -1 ) / (W*W), 1) )
    

    cuda.memcpy_dtoh(domain.quantities['xmomentum'].explicit_update, xmom_explicit_update_gpu)
    cuda.memcpy_dtoh(domain.quantities['ymomentum'].explicit_update, ymom_explicit_update_gpu)
	
    print "---------gravity_old ---------------"
    print domain.quantities['xmomentum'].explicit_update
    print domain.quantities['ymomentum'].explicit_update






def gravity_single(domain, k=0, flag = 4):
    import numpy
    w0 = domain.quantities['stage'].vertex_values[k][0]
    w1 = domain.quantities['stage'].vertex_values[k][1]
    w2 = domain.quantities['stage'].vertex_values[k][2]

    x0 = domain.vertex_coordinates[k*3][0]
    y0 = domain.vertex_coordinates[k*3][1]
    x1 = domain.vertex_coordinates[k*3+1][0]
    y1 = domain.vertex_coordinates[k*3+1][1]
    x2 = domain.vertex_coordinates[k*3+2][0]
    y2 = domain.vertex_coordinates[k*3+2][1]

    #_gradient(x0, y0, x1, y1, x2, y2, w0, w1, w2, &wx, &wy)
    det = (y2 - y0)*(x1 - x0) - (y1 - y0)*(x2 - x0)
    
    wx = (y2 -y0)*(w1 - w0) - (y1 - y0)*(w2 -w0)
    wx /= det
    
    wy = (x1 - x0)*(w2 - w0) - (x2 - x0)*(w1 -w0)
    wy /= det

    avg_h = domain.quantities['stage'].centroid_values[k] -\
            domain.quantities['elevation'].centroid_values[k]

    
    domain.quantities['xmomentum'].explicit_update[k] += -domain.g * wx*avg_h
    domain.quantities['ymomentum'].explicit_update[k] += -domain.g * wy*avg_h
    
    if flag == 3:
        return

    hh = numpy.zeros(3, dtype=numpy.float64)

    hh[0] = domain.quantities['stage'].edge_values[k][0] -\
            domain.quantities['elevation'].edge_values[k][0]
    hh[1] = domain.quantities['stage'].edge_values[k][1] -\
            domain.quantities['elevation'].edge_values[k][1]
    hh[2] = domain.quantities['stage'].edge_values[k][2] -\
            domain.quantities['elevation'].edge_values[k][2]

    
    
    sidex = 0.0
    sidey = 0.0
    for i in range(3): 
        n0 = domain.normals[k][2 * i]
        n1 = domain.normals[k][2 * i + 1]

        fact = -0.5 * domain.g * hh[i] * hh[i] * domain.edgelengths[k][i]
        
        
        sidex = sidex + fact*n0
        sidey = sidey + fact*n1

        if i == flag or ( i == 2 and flag == 5):
            domain.quantities['xmomentum'].explicit_update[k] = fact*n0
            domain.quantities['ymomentum'].explicit_update[k] = fact*n1
            if flag == 5:
                area = domain.areas[k]
                domain.quantities['xmomentum'].explicit_update[k] = sidex / area
                domain.quantities['ymomentum'].explicit_update[k] = sidey / area
            return
            
    area = domain.areas[k]
    domain.quantities['xmomentum'].explicit_update[k] += -sidex / area
    domain.quantities['ymomentum'].explicit_update[k] += -sidey / area


if __name__ == '__main__':
    from anuga_cuda import *
    from time import time

    #domain1 = generate_merimbula_domain()
    #domain2 = generate_merimbula_domain(True)
    domain1 = generate_channel3_domain()
    domain2 = generate_channel3_domain(True)


    domain2.equip_kernel_functions()

    testing_domain3 = False


    if domain1.compute_fluxes_method == 'original':
        from shallow_water_ext import gravity as gravity_c
        print "original"
        gravity_c(domain1)

    elif domain1.compute_fluxes_method == 'wb_1':
        from shallow_water_ext import gravity as gravity_c
        print "wb_1"
        gravity_c(domain1)

    elif domain1.compute_fluxes_method == 'wb_2':
        from shallow_water_ext import gravity_wb as gravity_wb_c
        print "wb_2"
        t1=time()
        gravity_wb_c(domain1)
        print " gravity duration : %lf" % (time()-t1)

    elif domain1.compute_fluxes_method == 'wb_3':
        from shallow_water_ext import gravity_wb as gravity_wb_c
        print "wb_3"
        gravity_wb_c(domain1)

    else:
        raise Exception('unknown compute_fluxes_method')


    import pycuda.driver as cuda
    import pycuda.autoinit
    from pycuda.compiler import SourceModule
    import numpy

    


    if domain2.compute_fluxes_method == 'original':
        gravity(domain2)

    elif domain2.compute_fluxes_method == 'wb_1':
        gravity(domain2)

    elif domain2.compute_fluxes_method == 'wb_2':
        #gravity_wb(domain2)

        import sys
        W1 = 0
        W2 = 0
        for i in range( len(sys.argv)):
            if sys.argv[i] == "-b":
                W1 = int(sys.argv[i+1])
            elif sys.argv[i] == "-b1":
                W2 = int(sys.argv[i+1])

        if not W1:
            W1 = domain2.gravity_wb_func.max_threads_per_block
        if not W2:
            W2 = 1 
        print W1, W2
        W3 = 1
        get_kernel_function_info(domain2.gravity_wb_func, W1, W2, W3)
        
        N = domain2.number_of_elements
        domain2.gravity_wb_func(
            numpy.uint( N),
            numpy.float64( domain2.g),
            cuda.In( domain2.quantities['stage'].vertex_values), 
            cuda.In( domain2.quantities['stage'].edge_values), 
            cuda.In( domain2.quantities['stage'].centroid_values), 
			cuda.In( domain2.quantities['elevation'].edge_values), 
            cuda.In( domain2.quantities['elevation'].centroid_values), 
            cuda.In( domain2.vertex_coordinates), 
            cuda.InOut( domain2.quantities['xmomentum'].explicit_update),
            cuda.InOut( domain2.quantities['ymomentum'].explicit_update),
            cuda.In( domain2.normals),
            cuda.In( domain2.areas),
            cuda.In( domain2.edgelengths),
            block = ( W1, W2, W3),
            grid = ( N/(W1*W2*W3), 1) )



    elif domain2.compute_fluxes_method == 'wb_3':

        gravity_wb(domain2)

    else:
        raise Exception('unknown compute_fluxes_method')





    print "******* xmom_explicit_update"
    #print domain1.quantities['xmomentum'].explicit_update
    #print domain2.quantities['xmomentum'].explicit_update
    counter = 0
    for i in range(domain1.number_of_elements):
        if  abs(domain1.quantities['xmomentum'].explicit_update[i] -\
                domain2.quantities['xmomentum'].explicit_update[i]) >\
                    abs(domain1.quantities['xmomentum'].explicit_update[i])*pow(10,-6):
            counter += 1
            if counter < 10:
                print i, domain1.quantities['xmomentum'].explicit_update[i],\
                        domain2.quantities['xmomentum'].explicit_update[i]
    print "---------> # of differences: %d" % counter


    print "******* ymom_explicit_update"
    #print domain1.quantities['ymomentum'].explicit_update
    #print domain2.quantities['ymomentum'].explicit_update
    counter = 0
    for i in range(domain1.number_of_elements):
        if abs( domain1.quantities['ymomentum'].explicit_update[i] -\
                domain2.quantities['ymomentum'].explicit_update[i]) >\
                    abs(domain1.quantities['ymomentum'].explicit_update[i])*pow(10,-6):
            counter += 1
    print "---------> # of differences: %d" % counter

    
    
    
    if testing_domain3:
        domain3 = generate_merimbula_domain()
        print "~~~~~~~~~~~~~ domain3  mom_explicit_update ~~~~~~~"
        flag = 4

        counter=0 
        for i in range(domain2.number_of_elements):
            gravity_single(domain3,i,flag)
            #if temp_array[i][flag][0]!= domain3.quantities['xmomentum'].explicit_update[i]:
            if abs(domain2.quantities['xmomentum'].explicit_update[i] -\
                    domain3.quantities['xmomentum'].explicit_update[i]) >\
                    abs(domain1.quantities['xmomentum'].explicit_update[i])*pow(10,-6):
                counter +=1
                if counter < 10:
                    print i , domain2.quantities['xmomentum'].explicit_update[i] ,\
                            domain3.quantities['xmomentum'].explicit_update[i]
        print  "---------> # of differences: %d" % counter


        counter = 0
        for i in range(domain2.number_of_elements):
            #gravity_single(domain3,i,flag)
            #if temp_array[i][flag][1] != domain3.quantities['ymomentum'].explicit_update[i]:
            if abs(domain2.quantities['ymomentum'].explicit_update[i] -\
                    domain3.quantities['ymomentum'].explicit_update[i]) >\
                    abs(domain1.quantities['ymomentum'].explicit_update[i])*pow(10,-6):
                counter +=1
                if counter < 10:
                    print i , domain2.quantities['ymomentum'].explicit_update[i] ,\
                            domain3.quantities['ymomentum'].explicit_update[i]
        print  "---------> # of differences: %d" % counter

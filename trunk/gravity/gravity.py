def gravity( domain ):
    N = domain.number_of_elements
    W = 32
    
    
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

    g_gpu = cuda.mem_alloc(numpy.float32.nbytes)
    cuda.memcpy_htod(g_gpu, domain.g)

    threadInBlock_gpu = cuda.mem_alloc(numpy.float32.nbytes)
    cuda.memcpy_htod(threadInBlock_gpu, W*W)

    gravity_mod = SourceModule("""
        __global__ void gravity(double *bed_vertex_values_gpu, double *stage_centroid_values_gpu, double *bed_centroid_values_gpu, double *vertex_coordinates_gpu, double * xmom_explicit_update_gpu, double *ymom_explicit_update_gpu, g_gpu, threadInBlock_gpu)
        {
            int i = threadIdx.x + threadIdx.y + blockIdx.x *threadInBlock_gpu;
            int k3 = 3*i,
                k6 = 6*i;
            double 
                z0 = bed_vertex_values_gpu[k3],
                z1 = bed_vertex_values_gpu[k3 + 1],
                z2 = bed_vertex_values_gpu[k3 + 2],

                avg_h = stage_centroid_values_gpu[i] - bed_centroid_values_gpu[i],

                x0 = vertex_coordinates_gpu[k6],
                y0 = vertex_coordinates_gpu[k6 + 1],
                x1 = vertex_coordinates_gpu[k6 + 2],
                y1 = vertex_coordinates_gpu[k6 + 3],
                x2 = vertex_coordinates_gpu[k6 + 4],
                y2 = vertex_coordinates_gpu[k6 + 5],
                
                zx, zy, det;

            det = (y2 - y0)*(x1 - x0) - (y1 - y0)*(x2 - x0);

            zx = (y2 -y0)*(z1 - z0) - (y1 - y0)*(z2 -z0);
            zx /= det;

            zy = (x1 - x0)*(z2 - z0) - (x2 - x0)*(z1 -z0);
            zy /= det;

            xmom_explicit_update_gpu[i] += -g * zx * avg_h;
            ymom_explicit_update_gpu[i] += -g * zy * avg_h;
        }
    """)

    gravity_func = gravity_mod.get_function("gravity")
    gravity_func(
            bed_vertex_values_gpu, 
            stage_centroid_values_gpu,
            bed_centroid_values_gpu,
            vertex_coordinates_gpu,
            xmom_explicit_update_gpu,
            ymom_explicit_update_gpu,
            g_gpu,
            threadInBlock_gpu,
            block = ( W, W, 1),
            grid = ( (N + W*W -1 ) / (W*W), 1)
            )
    
    cuda.memcpy_dtoh(domain.quantities['xmomentum'].explicit_update, xmom_explicit_update_gpu)
    cuda.memcpy_dtoh(domain.quantities['ymomentum'].explicit_update, ymom_explicit_update_gpu)



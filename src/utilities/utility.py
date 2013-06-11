#!/usr/bin/env python

def get_kernel_function_info(a, W1=0, W2=1, W3=1):
    """Show kernel information
    
    Including 
        1. max #threads per block, 
        2. active warps per MP, 
        3. thread block per MP, 
        4. usage of shared memory, 
        5. const memory , 
        6. local memory 
        7. registers
        8. hardware occupancy
        9. limitation of the hardware occupancy
    """

    import pycuda.tools as tl
    import pycuda.driver as dri
    dev = dri.Device(0)
    td = tl.DeviceData()
    if not W1:
        W1 = a.max_threads_per_block
    to = tl.OccupancyRecord(td, W1*W2*W3, a.shared_size_bytes, a.num_regs)

    print "***************************************"
    print "  Function Info    "
    print "   -> max threads per block: %d / %d / %d" % \
                (a.max_threads_per_block, 
                        dev.max_threads_per_block,
                        dev.max_threads_per_multiprocessor)
    print "   -> shared mem : %d / %d" % (a.shared_size_bytes, 
            td.shared_memory)
    print "   -> const mem : %d" % a.const_size_bytes
    print "   -> local mem : %d" % a.local_size_bytes
    print "   -> register : %d / %d" % (a.num_regs, td.registers)
    print "   -> thread block per MP %d / %d" % \
            (to.tb_per_mp, td.thread_blocks_per_mp)
    print "   -> warps per MP %d / %d" % (to.warps_per_mp, td.warps_per_mp)
    print "   -> occupancy %f" % to.occupancy
    print "   -> limitation %s" % to.limited_by
    print "  Block size : %dx%dx%d" % (W1, W2, W3)
    print "***************************************"




def get_sourceModule(k_dir, k_name, rearranged_domain=False):
    """Compile kernel code and return the PyCUDA function object"""

    from pycuda.compiler import SourceModule
    from anuga_cuda import archM, codeM

    if rearranged_domain:
        defince_macro = "#define REARRANGED_DOMAIN\n"
    else:
        defince_macro = ""
    return SourceModule(
            defince_macro + open( k_dir + k_name, "r").read(),
            arch = 'compute_20',
            code = 'sm_20',
            options =['-use_fast_math', '--compiler-options', '-O3'],
            include_dirs=[ k_dir ]
            )



def get_page_locked_array(a):
    """Replace the pageable array to page-locked array"""
    
    import pycuda.driver as drv

    temp_page_lock_p = drv.pagelocked_zeros_like(a,
            mem_flags=drv.host_alloc_flags.DEVICEMAP)
    if len(a.shape) == 1:
        temp_page_lock_p[:] = a
    else:
        temp_page_lock_p[:, :] = a
    assert numpy.allclose(a, temp_page_lock_p)           
    return temp_page_lock_p


def get_device_array(a):
    """Allocate device memory"""

    import pycuda.driver as drv

    return drv.mem_alloc(a.nbytes)



def asy_cpy(a, a_gpu, auto_init_context= True):
    """Data transfer from host to device.
    
    Asynchronous will be enabled when auto_init_context is True, otherwise
    use normal transfer.
    """

    import pycuda.driver as drv

    if auto_init_context:
        strm = drv.Stream()
        drv.memcpy_htod_async(a_gpu, a, strm)
        
        # Test correctness
        #ctx.synchronize()
        #b= numpy.zeros_like(a, a.dtype)
        #drv.memcpy_dtoh(b, a_gpu)
        #print numpy.allclose(a, b)
        return strm
    else:
        drv.memcpy_htod(a_gpu, a)



def cpy_back(a, a_gpu, auto_init_context=True):
    """Data transfer from device to host.
    
    Asynchronous will be enabled when auto_init_context is True, otherwise
    use normal transfer.
    """

    import pycuda.driver as drv

    if auto_init_context:
        strm = drv.Stream()
        drv.memcpy_dtoh_async(a, a_gpu, strm)
        return strm
    else:
        drv.memcpy_dtoh(a, a_gpu)



def cpy_back_and_cmp(a, b, value_type, gpu = True, rg = False):
    """Download mesh information and check result."""
    
    import numpy 

    if gpu:
        if value_type is "centroid_values":
            cpy_back(a.centroid_values, a.centroid_values_gpu)
            return numpy.allclose(a.centroid_values, b.centroid_values)
        elif value_type is "vertex_values":
            cpy_back(a.vertex_values, a.vertex_values_gpu)
            if rg:
                return check_rearranged_array(
                        b.vertex_values, a.vertex_values, 3)
            return numpy.allclose(a.vertex_values, b.vertex_values)
        elif value_type is "boundary_values":
            cpy_back(a.boundary_values, a.boundary_values_gpu)
            return numpy.allclose(a.boundary_values, b.boundary_values)
        elif value_type is "edge_values":
            cpy_back(a.edge_values, a.edge_values_gpu)
            if rg:
                return check_rearranged_array(
                        b.edge_values, a.edge_values, 3)
            return numpy.allclose(a.edge_values, b.edge_values)
        elif value_type is "x_gradient_values":
            cpy_back(a.x_gradient, a.x_gradient_gpu)
            return numpy.allclose(a.x_gradient, b.x_gradient)
        elif value_type is "y_gradient_values":
            cpy_back(a.y_gradient, a.y_gradient_gpu)
            return numpy.allclose(a.y_gradient, b.y_gradient)
        elif value_type is "explicit_update":
            cpy_back(a.explicit_update, a.explicit_update_gpu)
            return numpy.allclose(a.explicit_update, b.explicit_update)
        elif value_type is "semi_implicit_update":
            cpy_back(a.semi_implicit_update, a.semi_implicit_update_gpu)
            return numpy.allclose(a.semi_implicit_update,
                b.semi_implicit_update)
        elif value_type is "areas":
            cpy_back(a.areas, a.areas_gpu)
            return numpy.allclose(a.areas, b.areas)
        elif value_type is "surrogate_neighbours":
            cpy_back(a.surrogate_neighbours, a.surrogate_neighbours_gpu)
            if rg:
                return check_rearranged_array(
                        b.surrogate_neighbours, a.surrogate_neighbours, 3)
            return numpy.allclose(a.surrogate_neighbours, b.surrogate_neighbours)
        elif value_type is "number_of_boundaries":
            cpy_back(a.number_of_boundaries, a.number_of_boundaries_gpu)
            return numpy.allclose(a.number_of_boundaries, 
                    b.number_of_boundaries)
        elif value_type is "centroid_coordinates":
            cpy_back(a.centroid_coordinates, a.centroid_coordinates_gpu)
            if rg:
                return check_rearranged_array(
                        b.centroid_coordinates, a.centroid_coordinates, 2)
            return numpy.allclose(a.centroid_coordinates, 
                    b.centroid_coordinates)
        elif value_type is "vertex_coordinates":
            cpy_back(a.vertex_coordinates, a.vertex_coordinates_gpu)
            if rg:
                return check_rearranged_array(
                        b.vertex_coordinates, a.vertex_coordinates, 32)
            return numpy.allclose(a.vertex_coordinates, 
                    b.vertex_coordinates)
        elif value_type is "edge_coordinates":
            cpy_back(a.edge_coordinates, a.edge_coordinates_gpu)
            if rg:
                return check_rearranged_array(
                        b.edge_coordinates, a.edge_coordinates, 32)
            return numpy.allclose(a.edge_coordinates, 
                    b.edge_coordinates)
        else:
            raise Exception('Unknown value_type %s' % value_type)
    else:
        if value_type is "centroid_values":
            return numpy.allclose(a.centroid_values, b.centroid_values)
        elif value_type is "vertex_values":
            return numpy.allclose(a.vertex_values, b.vertex_values)
        elif value_type is "boundary_values":
            return numpy.allclose(a.boundary_values, b.boundary_values)
        elif value_type is "edge_values":
            return numpy.allclose(a.edge_values, b.edge_values)
        elif value_type is "x_gradient_values":
            return numpy.allclose(a.x_gradient, b.x_gradient)
        elif value_type is "y_gradient_values":
            return numpy.allclose(a.y_gradient, b.y_gradient)
        elif value_type is "explicit_update":
            return numpy.allclose(a.explicit_update, b.explicit_update)
        elif value_type is "semi_implicit_update":
            return numpy.allclose(
                a.semi_implicit_update, b.semi_implicit_update)
        elif value_type is "vertex_coordinates":
            return numpy.allclose(
                a.vertex_coordinates, b.vertex_coordinates)
        elif value_type is "areas":
            return numpy.allclose(a.areas, b.areas)
        elif value_type is "surrogate_neighbours":
            return numpy.allclose(
                a.surrogate_neighbours, b.surrogate_neighbours)
        elif value_type is "number_of_boundaries":
            return numpy.allclose(
                a.number_of_boundaries, b.number_of_boundaries)
        elif value_type is "centroid_coordinates":
            return numpy.allclose(
                a.centroid_coordinates, b.centroid_coordinates)
        elif value_type is "vertex_coordinates":
            return numpy.allclose(
                a.vertex_coordinates, b.vertex_coordinates)
        else:
            raise Exception('Unknown value_type %s' % value_type)



def number_domain_method(domain):
    """Convert mesh information stored in text string to ingeger.
    
    This is used in hmpp_pyhton_glue.

    Return value:
        (compute_fluxes_method, flow_algorithm, timestepping_method)
    """

    if domain.timestepping_method == 'euler':
        timestepping_method = 1
    elif domain.timestepping_method == 'rk2':
        timestepping_method = 2
    elif domain.timestepping_method == 'rk3':
        timestepping_method = 3
    else:
        timestepping_method = 4
    print " The timestepping_method is '%s' %d" % (domain.timestepping_method, timestepping_method)
    
    
    if domain.flow_algorithm == 'tsunami':
        flow_algorithm = 1
    elif domain.flow_algorithm == 'yusuke':
        flow_algorithm = 2
    else:
        flow_algorithm = 3
    print " The flow_algorithm us '%s' %d" % (domain.flow_algorithm, flow_algorithm)
    
    
    if domain.compute_fluxes_method == 'original':
        compute_fluxes_method = 0
    elif domain.compute_fluxes_method == 'wb_1':
        compute_fluxes_method = 1
    elif domain.compute_fluxes_method == 'wb_2':
        compute_fluxes_method = 2
    elif domain.compute_fluxes_method == 'wb_3':
        compute_fluxes_method = 3
    elif domain.compute_fluxes_method == 'tsunami':
        compute_fluxes_method = 4
    else:
        compute_fluxes_method = 5
    print " The compute_fluxes_method is '%s' %d" % (domain.compute_fluxes_method, compute_fluxes_method)


    return (compute_fluxes_method, flow_algorithm, timestepping_method)





#!/usr/bin/env python

from time import time
import numpy

auto_init_context = False
show_func_info = True

using_page_locked_array = True
testing_async = True
using_page_locked_mapped_pointer = True
testing_pl_mapped = True

testing_aligned_memory = False

testing = False


"""
Pycuda modules
"""
import pycuda.driver as drv
if auto_init_context:
    import pycuda.autoinit
else:
    drv.init()
    dev = drv.Device(0)
    ctx = dev.make_context( drv.ctx_flags.MAP_HOST)

from pycuda.compiler import SourceModule
import pycuda.tools as tl
td = tl.DeviceData()
cuda = drv


"""
ANUGA modules
"""
from anuga_cuda.config import compute_fluxes_dir
from anuga_cuda.merimbula_data.utility import approx_cmp
from anuga_cuda.merimbula_data.generate_domain import domain_create    
from anuga_cuda.merimbula_data.sort_domain import \
            sort_domain, rearrange_domain


domain1 = domain_create()
domain2 = rearrange_domain(domain1)

N = domain2.number_of_elements
timestep_array = numpy.zeros( N, dtype=numpy.float64)

# Necessary data initialization
for i in range(N):
    timestep_array[i] = domain2.evolve_max_timestep
    
if not (N % 320 == 0):
    W1 = 32
    W2 = 10
    W3 = 1
elif (N % 256 == 0):
    W1 = 16
    W2 = 16
    W3 = 1
elif (N % 32 ==0):
    W1 = 32
    W2 = 1
    W3 = 1
else:
    raise Exception('N can not be splited')
    
macro ="#define THREAD_BLOCK_SIZE %d\n" % (W1*W2*W3) 
mod = SourceModule(
        macro + open("compute_fluxes.cu","r").read(),
        #options=["--ptxas-options=-v"],
        include_dirs=[compute_fluxes_dir]
        )

compute_fluxes_central_function = mod.get_function(
            "compute_fluxes_central_structure_cuda_single")


"""********** ********** ********** ********** **********
    Get function info
********** ********** ********** ********** **********"""
def get_func_info(a, W1, W2, W3):
    to = tl.OccupancyRecord(td, W1*W2*W3, a.shared_size_bytes, a.num_regs)
    print "***************************************"
    print "  Function Info    "
    print "   -> max threads per block: %d" % a.max_threads_per_block
    print "   -> shared mem : %d / %d" % (a.shared_size_bytes, td.shared_memory)
    print "   -> const mem : %d" % a.const_size_bytes
    print "   -> local mem : %d" % a.local_size_bytes
    print "   -> register : %d / %d" % (a.num_regs, td.registers)
    print "   -> thread block per MP %d / %d" % (to.tb_per_mp, td.thread_blocks_per_mp)
    print "   -> warps per MP %d / %d" % (to.warps_per_mp, td.warps_per_mp)
    print "   -> occupancy %f" % to.occupancy
    print "   -> limitation %s" % to.limited_by
    print "  Block size : %dx%dx%d" % (W1, W2, W3)
    print "***************************************"

if show_func_info:
    get_func_info( compute_fluxes_central_function, W1, W2, W3)


strm = drv.Stream()
start = drv.Event()
end = drv.Event()


"""********** ********** ********** ********** **********
    Starting performance testing
********** ********** ********** ********** **********"""

"""********** ********** ********** ********** **********
     Pageable Memory
********** ********** ********** ********** **********"""
start.record()
compute_fluxes_central_function( 
        numpy.uint64( N ),
        numpy.float64( domain2.g ),
        numpy.float64( domain2.epsilon ),
        numpy.float64( domain2.H0 * domain2.H0 ),
        numpy.float64( domain2.H0*10 ),
        numpy.uint32( domain2.optimise_dry_cells),
        cuda.InOut( timestep_array ), 
        cuda.In( domain2.neighbours ), 
        cuda.In( domain2.neighbour_edges ),
        cuda.In( domain2.normals ), 
        cuda.In( domain2.edgelengths ), 
        cuda.In( domain2.radii ), 
        cuda.In( domain2.areas ), 
        cuda.In( domain2.tri_full_flag ),
        cuda.In( domain2.quantities['stage'].edge_values ), 
        cuda.In( domain2.quantities['xmomentum'].edge_values ), 
        cuda.In( domain2.quantities['ymomentum'].edge_values ), 
        cuda.In( domain2.quantities['elevation'].edge_values ), 
        cuda.In( domain2.quantities['stage'].boundary_values ), 
        cuda.In( domain2.quantities['xmomentum'].boundary_values ), 
        cuda.In( domain2.quantities['ymomentum'].boundary_values ), 
        cuda.InOut( domain2.quantities['stage'].explicit_update ), 
        cuda.InOut( domain2.quantities['xmomentum'].explicit_update ), 
        cuda.InOut( domain2.quantities['ymomentum'].explicit_update ),  
        cuda.InOut( domain2.max_speed), 
        block = ( W1, W2, W3),
        grid = (N/(W1*W2*W3), 1) )
end.record()
end.synchronize()
secs = start.time_till(end)*1e-3
print "Pageable memory exe time: %3.7f" % secs




"""********** ********** ********** ********** **********
    Page-locked array
********** ********** ********** ********** **********"""
def page_locked_array(a):
    a_pl = drv.pagelocked_zeros_like(a, mem_flags=drv.host_alloc_flags.DEVICEMAP)
    if len(a.shape) == 1:
        a_pl[:] = a
    else:
        a_pl[:,:] = a
    return a_pl


if using_page_locked_array:
    timestep_pla = page_locked_array( timestep_array )
    neighbours_pla = page_locked_array( domain2.neighbours )
    nei_edges_pla = page_locked_array( domain2.neighbour_edges )
    normals_pla = page_locked_array( domain2.normals )
    edgelen_pla = page_locked_array( domain2.edgelengths )
    radii_pla = page_locked_array( domain2.radii )
    areas_pla = page_locked_array( domain2.areas )
    tri_pla = page_locked_array( domain2.tri_full_flag )
    stage_edge_val_pla = page_locked_array( domain2.quantities['stage'].edge_values )
    xmom_edge_val_pla = page_locked_array( domain2.quantities['xmomentum'].edge_values )
    ymom_edge_val_pla = page_locked_array( domain2.quantities['ymomentum'].edge_values )
    bed_edge_val_pla = page_locked_array( domain2.quantities['elevation'].edge_values )
    stage_bou_val_pla = page_locked_array( domain2.quantities['stage'].boundary_values )
    xmom_bou_val_pla = page_locked_array( domain2.quantities['xmomentum'].boundary_values )
    ymom_bou_val_pla = page_locked_array( domain2.quantities['ymomentum'].boundary_values )
    stage_up_pla = page_locked_array( domain2.quantities['stage'].explicit_update )
    xmom_up_pla = page_locked_array( domain2.quantities['xmomentum'].explicit_update )
    ymom_up_pla = page_locked_array( domain2.quantities['ymomentum'].explicit_update )
    maxspeed_pla = page_locked_array( domain2.max_speed )


"""********** ********** ********** ********** **********
    Asynchronized data transmition
********** ********** ********** ********** **********"""
def async(a):
    strm = drv.Stream()
    gpu_ptr = drv.mem_alloc(a.nbytes)
    drv.memcpy_htod_async(gpu_ptr, a, strm)
    return gpu_ptr

if testing_async:
    start.record()
    
    timestep_gpuptr = async( timestep_pla)
    neighbours_gpuptr = async( neighbours_pla )
    nei_edges_gpuptr = async( nei_edges_pla )
    normals_gpuptr = async( normals_pla )
    edgelen_gpuptr = async( edgelen_pla )
    radii_gpuptr = async( radii_pla )
    areas_gpuptr = async( areas_pla )
    tri_gpuptr = async( tri_pla )
    stage_edge_val_gpuptr = async( stage_edge_val_pla )
    xmom_edge_val_gpuptr = async( xmom_edge_val_pla )
    ymom_edge_val_gpuptr = async( ymom_edge_val_pla )
    bed_edge_val_gpuptr = async( bed_edge_val_pla )
    stage_bou_val_gpuptr = async( stage_bou_val_pla )
    xmom_bou_val_gpuptr = async( xmom_bou_val_pla )
    ymom_bou_val_gpuptr = async( ymom_bou_val_pla )
    stage_up_gpuptr = async( stage_up_pla )
    xmom_up_gpuptr = async( xmom_up_pla )
    ymom_up_gpuptr = async( ymom_up_pla )
    maxspeed_gpuptr = async( maxspeed_pla )
    
    compute_fluxes_central_function( 
            numpy.uint64( N ),
            numpy.float64( domain2.g ),
            numpy.float64( domain2.epsilon ),
            numpy.float64( domain2.H0 * domain2.H0 ),
            numpy.float64( domain2.H0*10 ),
            numpy.uint32( domain2.optimise_dry_cells),
            timestep_gpuptr, 
            neighbours_gpuptr, 
            nei_edges_gpuptr,
            normals_gpuptr, 
            edgelen_gpuptr,
            radii_gpuptr,
            areas_gpuptr,
            tri_gpuptr,
            stage_edge_val_gpuptr,
            xmom_edge_val_gpuptr,
            ymom_edge_val_gpuptr,
            bed_edge_val_gpuptr,
            stage_bou_val_gpuptr,
            xmom_bou_val_gpuptr,
            ymom_bou_val_gpuptr,
            stage_up_gpuptr,
            xmom_up_gpuptr,
            ymom_up_gpuptr,
            maxspeed_gpuptr,
            block = ( W1, W2, W3),
            grid = (N/(W1*W2*W3), 1) )
    end.record()
    end.synchronize()
    secs = start.time_till(end)*1e-3
    print "Pinned memory aysnc exe time: %3.7f" % secs


"""********** ********** ********** ********** **********
    Get device mapped pointer
********** ********** ********** ********** **********"""
def get_dev_p(a):
    return numpy.intp(a.base.get_device_pointer())

if using_page_locked_mapped_pointer:
    timestep_plptr = get_dev_p( timestep_pla )
    neighbours_plptr = get_dev_p( neighbours_pla )
    nei_edges_plptr = get_dev_p( nei_edges_pla )
    normals_plptr = get_dev_p( normals_pla )
    edgelen_plptr = get_dev_p( edgelen_pla )
    radii_plptr = get_dev_p( radii_pla )
    areas_plptr = get_dev_p( areas_pla )
    tri_plptr = get_dev_p( tri_pla )
    stage_edge_val_plptr = get_dev_p( stage_edge_val_pla )
    xmom_edge_val_plptr = get_dev_p( xmom_edge_val_pla )
    ymom_edge_val_plptr = get_dev_p( ymom_edge_val_pla )
    bed_edge_val_plptr = get_dev_p( bed_edge_val_pla )
    stage_bou_val_plptr = get_dev_p( stage_bou_val_pla )
    xmom_bou_val_plptr = get_dev_p( xmom_bou_val_pla )
    ymom_bou_val_plptr = get_dev_p( ymom_bou_val_pla )
    stage_up_plptr = get_dev_p( stage_up_pla )
    xmom_up_plptr = get_dev_p( xmom_up_pla )
    ymom_up_plptr = get_dev_p( ymom_up_pla )
    maxspeed_plptr = get_dev_p( maxspeed_pla )



"""********** ********** ********** ********** **********
    Memory mapped memory
********** ********** ********** ********** **********"""
if testing_pl_mapped :
    start.record()
    compute_fluxes_central_function( 
            numpy.uint64( N ),
            numpy.float64( domain2.g ),
            numpy.float64( domain2.epsilon ),
            numpy.float64( domain2.H0 * domain2.H0 ),
            numpy.float64( domain2.H0*10 ),
            numpy.uint32( domain2.optimise_dry_cells),
            timestep_plptr, 
            neighbours_plptr, 
            nei_edges_plptr,
            normals_plptr, 
            edgelen_plptr,
            radii_plptr,
            areas_plptr,
            tri_plptr,
            stage_edge_val_plptr,
            xmom_edge_val_plptr,
            ymom_edge_val_plptr,
            bed_edge_val_plptr,
            stage_bou_val_plptr,
            xmom_bou_val_plptr,
            ymom_bou_val_plptr,
            stage_up_plptr,
            xmom_up_plptr,
            ymom_up_plptr,
            maxspeed_plptr,
            block = ( W1, W2, W3),
            grid = (N/(W1*W2*W3), 1) )
    end.record()
    end.synchronize()
    secs = start.time_till(end)*1e-3
    print "Pinned mapped memory exe time: %3.7f" % secs


"""********** ********** ********** ********** **********
    Page-locked aligned memory
********** ********** ********** ********** **********"""
def aligned_mem(a):
    temp = drv.aligned_empty(a.shape, dtype=a.dtype, order='C')
    if len(a.shape) == 1:
        temp[:] = a
    else:
        temp[:,:] = a
    return temp

if testing_aligned_memory:    
    timestep_alm = aligned_mem( timestep_array )
    neighbours_alm = aligned_mem( domain2.neighbours )
    nei_edges_alm = aligned_mem( domain2.neighbour_edges )
    normals_alm = aligned_mem( domain2.normals )
    edgelen_alm = aligned_mem( domain2.edgelengths )
    radii_alm = aligned_mem( domain2.radii )
    areas_alm = aligned_mem( domain2.areas )
    tri_alm = aligned_mem( domain2.tri_full_flag )
    stage_edge_val_alm = aligned_mem( domain2.quantities['stage'].edge_values )
    xmom_edge_val_alm = aligned_mem( domain2.quantities['xmomentum'].edge_values )
    ymom_edge_val_alm = aligned_mem( domain2.quantities['ymomentum'].edge_values )
    bed_edge_val_alm = aligned_mem( domain2.quantities['elevation'].edge_values )
    stage_bou_val_alm = aligned_mem( domain2.quantities['stage'].boundary_values )
    xmom_bou_val_alm = aligned_mem( domain2.quantities['xmomentum'].boundary_values )
    ymom_bou_val_alm = aligned_mem( domain2.quantities['ymomentum'].boundary_values )
    stage_up_alm = aligned_mem( domain2.quantities['stage'].explicit_update )
    xmom_up_alm = aligned_mem( domain2.quantities['xmomentum'].explicit_update )
    ymom_up_alm = aligned_mem( domain2.quantities['ymomentum'].explicit_update )
    maxspeed_alm = aligned_mem( domain2.max_speed )
    
    
    """********** ********** ********** ********** **********
        Register page-locked aligned memory
    ********** ********** ********** ********** **********"""
    def reg_mem(a):
        temp = drv.register_host_memory( a, 
                    flags=drv.mem_host_register_flags.DEVICEMAP)
        return numpy.intp(temp.base.get_device_pointer())
    
    timestep_plalm = reg_mem( timestep_alm )
    neighbours_plalm = reg_mem( neighbours_alm )
    nei_edges_plalm = reg_mem( nei_edges_alm )
    normals_plalm = reg_mem( normals_alm )
    edgelen_plalm = reg_mem( edgelen_alm )
    radii_plalm = reg_mem( radii_alm )
    areas_plalm = reg_mem( areas_alm )
    tri_plalm = reg_mem( tri_alm )
    stage_edge_val_plalm = reg_mem( stage_edge_val_alm )
    xmom_edge_val_plalm = reg_mem( xmom_edge_val_alm )
    ymom_edge_val_plalm = reg_mem( ymom_edge_val_alm )
    bed_edge_val_plalm = reg_mem( bed_edge_val_alm )
    stage_bou_val_plalm = reg_mem( stage_bou_val_alm )
    xmom_bou_val_plalm = reg_mem( xmom_bou_val_alm )
    ymom_bou_val_plalm = reg_mem( ymom_bou_val_alm )
    stage_up_plalm = reg_mem( stage_up_alm )
    xmom_up_plalm = reg_mem( xmom_up_alm )
    ymom_up_plalm = reg_mem( ymom_up_alm )
    maxspeed_plalm = reg_mem( maxspeed_alm )
    
    start.record()
    compute_fluxes_central_function( 
            numpy.uint64( N ),
            numpy.float64( domain2.g ),
            numpy.float64( domain2.epsilon ),
            numpy.float64( domain2.H0 * domain2.H0 ),
            numpy.float64( domain2.H0*10 ),
            numpy.uint32( domain2.optimise_dry_cells),
            timestep_plalm, 
            neighbours_plalm, 
            nei_edges_plalm,
            normals_plalm, 
            edgelen_plalm,
            radii_plalm,
            areas_plalm,
            tri_plalm,
            stage_edge_val_plalm,
            xmom_edge_val_plalm,
            ymom_edge_val_plalm,
            bed_edge_val_plalm,
            stage_bou_val_plalm,
            xmom_bou_val_plalm,
            ymom_bou_val_plalm,
            stage_up_plalm,
            xmom_up_plalm,
            ymom_up_plalm,
            maxspeed_plalm,
            block = ( W1, W2, W3),
            grid = (N/(W1*W2*W3), 1) )
    end.record()
    end.synchronize()
    secs = start.time_till(end)*1e-3
    print "Aligned memory exe time: %3.7f" % secs



"""********** ********** ********** ********** **********
    Get the result
********** ********** ********** ********** **********"""
#b = numpy.argsort(timestep_array)
#domain2.flux_timestep = timestep_array[b[0]] 



"""********** ********** ********** ********** **********
    ANUGA C
********** ********** ********** ********** **********"""
from shallow_water_ext import compute_fluxes_ext_central_structure
start_time = time()
domain1.flux_timestep = compute_fluxes_ext_central_structure(domain1)
duration = time() - start_time
print "ANUGA C exe time: %3.7f" % duration



if not auto_init_context:
    ctx.pop()

"""********** ********** ********** ********** **********
     Start validation testing
********** ********** ********** ********** **********"""
if testing:
    """
    ANUGA C Function
    """
    if domain1.compute_fluxes_method == 'original':
        from shallow_water_ext import compute_fluxes_ext_central_structure
        from shallow_water_ext import gravity as gravity_c
        
        print "original"
        domain1.flux_timestep = compute_fluxes_ext_central_structure(domain1)
        #gravity_c(domain1)
        
    elif domain1.compute_fluxes_method == 'wb_1':
        # Calc pressure terms using Simpson rule in flux
        # computations. Then they match up exactly with
        # standard gravity term - g h grad(z)
        from shallow_water_ext import compute_fluxes_ext_wb
        from shallow_water_ext import gravity as gravity_c

        print "wb_1"
        domain1.flux_timestep = compute_fluxes_ext_wb(domain1)
        #gravity_c(domain1)

    elif domain1.compute_fluxes_method == 'wb_2':
        # Use standard flux calculation, but calc gravity
        # as -g h grad(w) - sum midpoint edge pressure terms
        from shallow_water_ext import compute_fluxes_ext_central_structure
        from shallow_water_ext import gravity_wb as gravity_wb_c

        #print "wb_2"
        #t1 = time()
        domain1.flux_timestep = compute_fluxes_ext_central_structure(domain1)
        #gravity_wb_c(domain1)
        #print "\n C compute flux duration is %lf\n" % (time()-t1)
    elif domain1.compute_fluxes_method == 'wb_3':
        # Calculate pure flux terms with simpsons rule, and
        # gravity flux and gravity forcing via
        # as -g h grad(w) - sum midpoint edge pressure terms
        from shallow_water_ext import compute_fluxes_ext_wb_3
        from shallow_water_ext import gravity_wb as gravity_wb_c

        print "wb_3"
        domain1.flux_timestep = compute_fluxes_ext_wb_3(domain1)
        #gravity_wb_c(domain1)

    else:
        raise Exception('unknown compute_fluxes_method')
    



    print "\n~~~~~~~~~~~~~ domain 2 ~~~~~~~~~~~~"
    print "******* flux_timestep : %lf %lf %d" % \
            (domain1.flux_timestep, domain2.flux_timestep, \
                domain1.flux_timestep == domain2.flux_timestep)

    print "******* max_speed"
    counter = 0
    for i in range(domain1.number_of_elements):
        if ( domain1.max_speed[i] != domain2.max_speed[i]):
            counter += 1
    print "---------> # of differences: %d" % counter


    print "******* stage_explicit_update"
    #print domain1.quantities['stage'].explicit_update
    #print domain2.quantities['stage'].explicit_update
    counter = 0
    for i in range(domain1.number_of_elements):
        if approx_cmp( domain1.quantities['stage'].explicit_update[i] ,
                domain2.quantities['stage'].explicit_update[i]):
            counter += 1
    print "---------> # of differences: %d" % counter


    print "******* xmom_explicit_update"
    #print domain1.quantities['xmomentum'].explicit_update
    #print domain2.quantities['xmomentum'].explicit_update
    counter = 0
    for i in range(domain1.number_of_elements):
        if approx_cmp( domain1.quantities['xmomentum'].explicit_update[i] ,
                domain2.quantities['xmomentum'].explicit_update[i]):
            counter += 1
            if counter < 10:
                print i, domain1.quantities['xmomentum'].explicit_update[i], \
                        domain2.quantities['xmomentum'].explicit_update[i]
    print "---------> # of differences: %d" % counter

    
    print "******* ymom_explicit_update"
    #print domain1.quantities['ymomentum'].explicit_update
    #print domain2.quantities['ymomentum'].explicit_update
    counter = 0
    for i in range(domain1.number_of_elements):
        if approx_cmp( domain1.quantities['ymomentum'].explicit_update[i],
                domain2.quantities['ymomentum'].explicit_update[i]):
            counter += 1
    print "---------> # of differences: %d" % counter


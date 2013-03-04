#!/usr/bin/env python

def approx_equal(a,b, approx=True):
    if approx:
        if abs(a-b) > abs(a)*pow(10,-6):
            return False
        else: 
            return True
    else:
        if a != b:
            return False
        else:
            return True


def mem_all_cpy(a):
    import pycuda.driver as cuda
    a_gpu = cuda.mem_alloc(a.nbytes)
    cuda_memcpy_htod(a_gpu, a)
    return a_gpu




"""********** ********** ********** ********** **********
    Get kernel function info
********** ********** ********** ********** **********"""
def get_kernel_function_info(a, W1=0, W2=1, W3=1):
    import pycuda.tools as tl
    import pycuda.driver as dri
    dev = dri.Device(0)
    td = tl.DeviceData()
    if not W1:
        W1 = a.max_threads_per_block
    to = tl.OccupancyRecord(td, W1*W2*W3, a.shared_size_bytes, a.num_regs)

    print "***************************************"
    print "  Function Info    "
    print "   -> max threads per block: %d / %d" % \
                (a.max_threads_per_block, dev.max_threads_per_block)
    print "   -> shared mem : %d / %d" % (a.shared_size_bytes, td.shared_memory)
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


def configure_kernel_function_block_dimension(a, tl=None, td=None, dev=None):
    if dev is None:
        import pycuda.autoinit
        dev = pycuda.autoinit.device
    if tl is None:
        import pycuda.tools as tl
        if td is None:
            td = tl.DeviceData()

    threadBlocks = min( td.shared_memory / a.shared_size_bytes,
                        td.registers / a.num_regs ,
                        dev.max_threads_per_block)

    to = tl.OccupancyRecord(td, a.max_threads_per_block, 
                        a.shared_size_bytes, a.num_regs)

    print "\n***************************************"
    #print "* %s :" % a.__name__
    print "   -> max threads per block: %d / %d" % \
                (a.max_threads_per_block, dev.max_threads_per_block)
    print "   -> shared mem : %d / %d" % (a.shared_size_bytes, td.shared_memory)
    print "   -> register : %d / %d" % (a.num_regs, td.registers)
    print "   -> const mem : %d" % a.const_size_bytes
    print "   -> local mem : %d" % a.local_size_bytes
    print "   -> occupancy %f" % to.occupancy
    print "   -> limitation %s" % to.limited_by
    print "  Block size : %d" % threadBlocks
    print "***************************************"


    return threadBlocks
    

if __name__ == '__main__':
    from anuga_cuda import generate_merimbula_domain
    d = generate_merimbula_domain(True)
    d.equip_kernel_functions()

    configure_kernel_function_block_dimension(d.compute_fluxes_func)

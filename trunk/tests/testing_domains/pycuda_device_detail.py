#!/usr/bin/env python

from pycuda import tools as tl
from pycuda import driver as dr
import pycuda.autoinit

td = tl.DeviceData()
to = tl.OccupancyRecord

def list_device_detail(td):
    print "  **  max_threads %d" % td.max_threads
    print "  **  warps per mp %d" % td.warps_per_mp
    print "  **  thread_blocks_per_mp %d" % td.thread_blocks_per_mp
    print "  **  registers %d" % td.registers
    print "  **  shared_memory %d" % td.shared_memory

def calculate_occupancy(threads, smem, registers):
    occ = to(td, threads, smem, registers)
    print "  **  Occupancy %lf" % occ.occupancy
    print "  **  thread block per SM %d" % occ.tb_per_mp
    print "  **  warps per SM %d" % occ.warps_per_mp
    print "  **  limited %s" % occ.limited_by

def list_function_detail(func):
    print "  **  max_threads _per_block %d " % func.max_threads_per_block
    print "  **  shared mem %d " % func.shared_size_bytes
    print "  **  const mem %d " % func.const_size_bytes
    print "  **  local mem %d " % func.local_size_bytes
    print "  **  register %d " % func.num_regs
    print "  **  ptx and binary version %d , %d" % (func.ptx_version, func.binary_version )



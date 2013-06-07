#!/usr/bin/env python

"""
Compare performance of using host-registered pinned and unpinned host memory,
with more than one block for larger arrays, and with unpinned tried first.
"""

import numpy as np

import pycuda.autoinit
import pycuda.driver as drv
from pycuda.compiler import SourceModule

from time import time

increment_mod = SourceModule("""
__global__ void increment(double *a, int N)
{
    int idx = threadIdx.x + blockIdx.x*blockDim.x;
    if (idx < N)
        a[idx] = a[idx]+1;
}
""")
increment = increment_mod.get_function("increment")

N = 20 # breaks. Works if <= 22
M = 3

#from anuga_cuda.merimbula_data.generate_domain import domain_create
#domain = domain_create()

 Time use of pageable host memory:
x = np.empty((N, N), np.float64)
#x = domain.radii
#N = domain.number_of_elements

times = np.empty(M)
for i in xrange(M):
    x[:, :] = np.random.rand(N, N)
    x_orig = x.copy()
    start = time()
    increment(
            drv.InOut(x), 
            np.uint32(N), 
            block=(320, 1, 1),
            grid=(N/320,1,1)
            )
    times[i] = time()-start
    np.allclose(x_orig + 1, x), "%r %r" % (x_orig, x)

print "Average kernel execution time with pageable memory: %3.7f" % np.mean(times)

# Time use of pinned host memory:
#x = drv.aligned_empty((N, N), dtype=np.float64, order='C')
x = drv.register_host_memory(x, flags=drv.mem_host_register_flags.DEVICEMAP)

x_gpu_ptr = np.intp(x.base.get_device_pointer())

times = np.empty(M)
for i in xrange(M):
    #x[:, :] = np.random.rand(N, N)
    x_orig = x.copy()
    start = time()
    increment(
            x_gpu_ptr, 
            np.uint32(N), 
            block=(320, 1, 1), 
            grid=(N/320,1,1)
            )
    times[i] = time()-start
print "Average kernel execution time with pinned memory:   %3.7f" % np.mean(times)

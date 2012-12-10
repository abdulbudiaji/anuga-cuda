#!/usr/bin/env python

"""
Compare performance of using pinned and unpinned host memory.
"""

import numpy as np

import pycuda.autoinit
import pycuda.driver as drv
from pycuda.compiler import SourceModule

from time import time

increment_mod = SourceModule("""
__global__ void increment(double *a, int N)
{
    int idx = threadIdx.x;
    if (idx < N)
        a[idx] = a[idx]+1;
}
""")
increment = increment_mod.get_function("increment")

N = 20
M = 3

# Time use of pinned host memory:
x = drv.pagelocked_empty((N, N), np.float64, mem_flags=drv.host_alloc_flags.DEVICEMAP)
x_gpu_ptr = np.intp(x.base.get_device_pointer())

times = np.empty(M)
for i in xrange(M):
    x[:, :] = np.random.rand(N, N)
    x_orig = x.copy()
    start = time()
    increment(x_gpu_ptr, np.uint32(x.size), block=(512, 1, 1))
    times[i] = time()-start
    np.allclose(x_orig + 1, x)

print "Average kernel execution time with pinned memory:   %3.7f" % np.mean(times)

# Time use of pageable host memory:
x = np.empty((N, N), np.float64)

times = np.empty(M)
for i in xrange(M):
    x[:, :] = np.random.rand(N, N)
    x_orig = x.copy()
    start = time()
    increment(drv.InOut(x), np.uint32(x.size), block=(512, 1, 1))
    times[i] = time()-start
    np.allclose(x_orig + 1, x)

print "Average kernel execution time with pageable memory: %3.7f" % np.mean(times)

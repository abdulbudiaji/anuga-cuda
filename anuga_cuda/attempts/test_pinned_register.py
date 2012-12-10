#!/usr/bin/env python

"""
Compare performance of using host-registered pinned and unpinned host memory.
"""

import numpy as np

import pycuda.autoinit
import pycuda.driver as drv
from pycuda.compiler import SourceModule

from time import time


#drv.init()
#dev = drv.Device(0)
#ctx = dev.make_context(drv.ctx_flags.SCHED_AUTO | drv.ctx_flags.MAP_HOST)

increment_mod = SourceModule("""
__global__ void increment(double *a, int N)
{
    int idx = threadIdx.x;
    if (idx < N)
        a[idx] = a[idx]+1;
}
""")
increment_mod2 = SourceModule("""
__global__ void increment(int N, double *a)
{
    int idx = threadIdx.x;
    if (idx < N)
        a[idx] = a[idx]+1;
}
""")
increment = increment_mod.get_function("increment")
increment2 = increment_mod2.get_function("increment")

N = 24 
M = 3

# Time use of pinned host memory:
#x = drv.aligned_empty((N, N), dtype=np.float64, order='C')
from anuga_cuda.merimbula_data.generate_domain import domain_create
domain = domain_create()

#x = drv.register_host_memory(x, flags=drv.mem_host_register_flags.DEVICEMAP)
x = drv.register_host_memory(domain.radii, flags=drv.mem_host_register_flags.DEVICEMAP)
#x_gpu_ptr = np.intp(x.base.get_device_pointer())

a = drv.pagelocked_zeros(x.size, dtype=np.float64)
a[:] = domain.radii
x_gpu_ptr = drv.mem_alloc(x.nbytes)
y_gpu_ptr = drv.mem_alloc(x.nbytes)
z_gpu_ptr = drv.mem_alloc(x.nbytes)



strm1 = drv.Stream()
strm2 = drv.Stream()
strm3 = drv.Stream()

start = drv.Event()
end = drv.Event()



start.record()
drv.memcpy_htod_async(x_gpu_ptr, domain.radii, strm1)
increment(x_gpu_ptr, np.uint32(x.size), block=(320, 1, 1), grid=(x.size/320,1), stream=strm1)
drv.memcpy_htod_async(y_gpu_ptr, a, strm2)
increment(y_gpu_ptr, np.uint32(x.size), block=(320, 1, 1), grid=(x.size/320,1), stream=strm2)
drv.memcpy_htod_async(z_gpu_ptr, a, strm3)
increment(z_gpu_ptr, np.uint32(x.size), block=(320, 1, 1), grid=(x.size/320,1), stream=strm3)
end.record()
end.synchronize()
secs = start.time_till(end)*1e-3
print "3 kernel execution time with async memcpy pinned:   %3.7f" % secs

start.record()
increment(drv.InOut(a), np.uint32(x.size), block=(320, 1, 1), grid=(x.size/320, 1))
increment(drv.InOut(a), np.uint32(x.size), block=(320, 1, 1), grid=(x.size/320, 1))
increment(drv.InOut(a), np.uint32(x.size), block=(320, 1, 1), grid=(x.size/320, 1))
end.record()
end.synchronize()
secs = start.time_till(end)*1e-3
print "3 kernel execution time with memcpy pinned:   %3.7f" % secs

start.record()
increment(drv.InOut(domain.areas), np.uint32(x.size), block=(320, 1, 1), grid=(x.size/320, 1))
increment(drv.InOut(domain.areas), np.uint32(x.size), block=(320, 1, 1), grid=(x.size/320, 1))
increment(drv.InOut(domain.areas), np.uint32(x.size), block=(320, 1, 1), grid=(x.size/320, 1))
end.record()
end.synchronize()
secs = start.time_till(end)*1e-3
print "3 kernel execution time with memcpy:   %3.7f" % secs


x_gpu_ptr = np.intp(x.base.get_device_pointer())
start.record()
increment(x_gpu_ptr, np.uint32(x.size), block=(320, 1, 1), grid=(x.size/320,1))
increment(x_gpu_ptr, np.uint32(x.size), block=(320, 1, 1), grid=(x.size/320,1))
increment(x_gpu_ptr, np.uint32(x.size), block=(320, 1, 1), grid=(x.size/320,1))
end.record()
end.synchronize()
secs = start.time_till(end)*1e-3
print "3 kernel execution time with pinned:   %3.7f measured in event" % secs


start = time()
increment(x_gpu_ptr, np.uint32(x.size), block=(320, 1, 1), grid=(x.size/320,1))
increment(x_gpu_ptr, np.uint32(x.size), block=(320, 1, 1), grid=(x.size/320,1))
increment(x_gpu_ptr, np.uint32(x.size), block=(320, 1, 1), grid=(x.size/320,1))
secs = time()-start
print "3 kernel execution time with pinned:   %3.7f measured in time" % secs


areas_handler = domain.areas
start = time()
for i in range(domain.number_of_elements):
    areas_handler[i] = areas_handler[i] +1
for i in range(domain.number_of_elements):
    areas_handler[i] = areas_handler[i] +1
for i in range(domain.number_of_elements):
    areas_handler[i] = areas_handler[i] +1
secs = time()-start
print "3 kernel execution time python:   %3.7f" % secs




times = np.empty(M)
for i in xrange(M):
    #x[:, :] = np.random.rand(N, N)
    #x_orig = x.copy()
    #print x, x.size
    start = time()
    increment(x_gpu_ptr, np.uint32(x.size), block=(320, 1, 1), grid=(x.size/320,1))
    times[i] = time()-start
    #np.allclose(x_orig + 1, x)
    
print "3 kernel execution time with pinned memory:   %3.7f" % np.mean(times)


# Time use of pageable host memory:
N = x.size
x = np.empty(N, np.float64)
times = np.empty(M)
for i in xrange(M):
    x[:] = np.random.rand(N)
#    x_orig = x.copy()
    start = time()
    increment(drv.InOut(a), np.uint32(x.size), block=(320, 1, 1), grid=(x.size/320, 1))
    times[i] = time()-start
#    np.allclose(x_orig + 1, x)
#
print "Average kernel execution time with pageable memory: %3.7f" % np.mean(times)


#ctx.pop()

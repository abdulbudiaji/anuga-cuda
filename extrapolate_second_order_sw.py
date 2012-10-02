def extrapolate_second_order_sw_cuda(domain=None):
    import pycuda.driver as cuda
    import pycuda.autoinit
    from pycuda.compiler import SourceModule
    import numpy

    a = numpy.random.randn(4,4)
    a = a.astype(numpy.float32)
    a_gpu = cuda.mem_alloc(a.nbytes)
    cuda.memcpy_htod(a_gpu, a)

    mod = SourceModule("""
        __global__ void doublify(float *a)
        {
            int idx = threadIdx.x + threadIdx.y*4;
            a[idx] *= 2;
        }
    """)

    func = mod.get_function("doublify")
    func(a_gpu, block=(4,4,1))

    a_doubled = numpy.empty_like(a)
    cuda.memcpy_dtoh(a_doubled, a_gpu)

    print a_doubled
    print a


def double_64(a):
    import pycuda.driver as cuda
    import pycuda.autoinit
    from pycuda.compiler import SourceModule
    import numpy

    a_gpu = cuda.mem_alloc(a.nbytes)
    cuda.memcpy_htod(a_gpu, a)

    mod = SourceModule("""
        __global__ void doublify(double *a)
        {
            int idx = threadIdx.x + threadIdx.y*4;
            a[idx] *= 2;
        }
    """)

    func = mod.get_function("doublify")
    func(a_gpu, block=(4,4,1))

    a_doubled = numpy.empty_like(a)
    cuda.memcpy_dtoh(a_doubled, a_gpu)
    print a_doubled
    


if __name__ == '__main__':
    import numpy
    a = numpy.random.randn(4,4)
    a = a.astype(numpy.float64)
    print a
    double_64(a)
    a = a.astype(numpy.float32)
    print a 
    double_64(a)

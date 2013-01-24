#!/usr/bin/env python

def approx_cmp(a,b, approx=True):
    if approx:
        if abs(a-b) > abs(a)*pow(10,-6):
            return True
        else: 
            return False
    else:
        if a != b:
            return True
        else:
            return False


def mem_all_cpy(a):
    import pycuda.driver as cuda
    a_gpu = cuda.mem_alloc(a.nbytes)
    cuda_memcpy_htod(a_gpu, a)
    return a_gpu

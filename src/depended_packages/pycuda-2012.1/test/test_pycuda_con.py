#!/usr/bin/env python

import pycuda.driver as cuda
import pycuda.autoinit
from pycuda.compiler import SourceModule

import numpy

def pycuda_Doubled():
    a = numpy.array([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16])

    a = a.astype(numpy.float64)
    mod=SourceModule("""

        __device__ double cuda_atomicAdd(double *address, double val)
        {
            double assumed,old=*address;
            //do {
                assumed=old;
                old= __longlong_as_double(atomicCAS((unsigned long long int*)address,
                    __double_as_longlong(assumed),
                    __double_as_longlong(val+assumed)));
            //}while (assumed!=old);

            return old;
        }
            
        __global__ void doublify(double *a)  
        {    
            //float b = max(1.0,fabs(-3.1));
            int idx = threadIdx.x + threadIdx.y*4;
            //if (a[idx]== idx+1)
			//    a[idx] = b;
            
            cuda_atomicAdd(a, a[idx]);
		}  
		""")
	
    print a
	#func = mod.get_function("doublify")
	#func(cuda.InOut(a), block=(4,4,1))
    a_gpu = cuda.mem_alloc(a.nbytes)
    cuda.memcpy_htod(a_gpu, a)
    print a
    cuda.memcpy_dtoh(a, a_gpu)
    print a
    cuda.memset_d32(a_gpu, 0, 16*2)
    cuda.memcpy_dtoh(a, a_gpu)
    print a
    b = numpy.float64(0)
    print b==a[0]

	
	#func(cuda.InOut(a), block=(4,4,1))
	#print a

if __name__ == "__main__":
	pycuda_Doubled()
	

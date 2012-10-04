#!/usr/bin/env python

import pycuda.driver as cuda
import pycuda.autoinit
from pycuda.compiler import SourceModule

import numpy

def pycuda_Doubled():
	a = numpy.random.randn(4,4)

	a = a.astype(numpy.float32)

	mod=SourceModule("""
		__global__ void doublify(float *a)  
		{    
			int idx = threadIdx.x + threadIdx.y*4;    
			a[idx] *= 2;  
		}  
		""")
	
	print a
	func = mod.get_function("doublify")
	func(cuda.InOut(a), block=(4,4,1))


	print a

if __name__ == "__main__":
	pycuda_Doubled()
	

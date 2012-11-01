
CC = gcc
INCLUDES = -I"/usr/lib64/python2.6/site-packages/numpy/core/include"  -I"/usr/include/python2.6" -I"../anuga/utilities/"
CXXFLAGS += -O3 -Wall -fPIC ${INCLUDES}

NVCC = nvcc
NVCCFLAGS = -O3 --gpu-architecture sm_13 ${IPATH} -I${CUDA_INSTALL_PATH}/include -Xcompiler -fpic --ptxas-options="-v"
# CUDA_INSTALL_PATH should evaluate to cuda home folder (where bin/ include/ lib/ etc. are located);
# ptxas -verbose adds some kernel information output;
# Be sure to adjust these environment variables:
#	export PATH=$PATH:${CUDA_INSTALL_PATH}/bin:${CUDA_INSTALL_PATH}/open64/bin
#	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${CUDA_INSTALL_PATH}/lib64/:${CUDA_INSTALL_PATH}/computeprof/bin/

LIBS = -L$(CUDA_INSTALL_PATH)/lib64 -lcudart -lm

.PHONY : clean

gpu_python_glue.so : gpu_python_glue.o cudafun.o
	${CC} -shared $< cudafun.o -o $@ ${LIBS}

%.o : %.c
	${CC} ${CXXFLAGS} -c $< -o $@

%.o : %.cu
	$(NVCC) $(NVCCFLAGS) -o $@ -c $<

clean : 
	rm -f *.o
	rm -f *.so 

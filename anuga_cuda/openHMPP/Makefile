gravity = gravity
compute_fluxes = compute_fluxes
hmpp_python_glue = hmpp_python_glue

HMPPRT_PATH = /home/devres/john/openHMPP/

INCLUDES = -I"/usr/lib/python2.7/dist-packages/numpy/core/include/" -I"/usr/include/python2.7" -I"/home/u5044856/anuga_core/source/anuga/utilities/"
CC = gcc
CXXFLAGS += -O0 -fno-inline -g -std=c99 -O3 -Wall -fPIC ${INCLUDES}
CXX = g++
HMPP_FLAGS13 = -e --nvcc-options -Xptxas=-v,-arch=sm_13 -c --force
HMPP_FLAGS20 = -e --nvcc-options -Xptxas=-v,-arch=sm_15 -c --force
#HMPP_FLAGS += -m

LIBS = -lm

CUDA_ENV_VAR = CUDA_PROFILE=1 
HMPP_ENV_VAR = HMPPRT_LOG_LEVEL=info HMPPRT_NO_FALLBACK=0


all: gravity compute_fluxes hmpp_python_glue.so

hmpp_python_glue.so : hmpp_python_glue.o gravity.o
	${CC} -shared $< gravity.o gravity.o -o $@ ${LIBS} 

%.o : %.c
	${CC} ${CXXFLAGS} -c $< -o $@


gravity: g.13
	#$(CUDA_ENV_VAR) $(HMPP_ENV_VAR) $(HMPPRT_PATH)g.13
	$(CUDA_ENV_VAR) $(HMPPRT_PATH)g.13

g.13: $(gravity).c
	capsmc $(HMPP_FLAGS13) $(CC) $(CXXFLAGS) $^ -o $@

compute_fluxes: cf.13
	#$(CUDA_ENV_VAR) $(HMPP_ENV_VAR) $(HMPPRT_PATH)cf.13
	$(CUDA_ENV_VAR) $(HMPPRT_PATH)cf.13

cf.13: $(compute_fluxes).c
	capsmc $(HMPP_FLAGS13) $(CC) $(CXXFLAGS) $^ -o $@



clean:
	rm -f example?/app_ported_cuda.hmf.o example?/app_ported_cuda.hmf.fatbin hmpprt-quickstart.html example?/app_ported.preproc.ii \
           hmppcg_functions.dpil example1/app_orig.exe example?/app_ported_cuda.hmf example?/app_ported_cuda.hmf.cu example?/app_ported.exe 

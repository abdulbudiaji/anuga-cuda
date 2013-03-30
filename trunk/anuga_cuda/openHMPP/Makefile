gravity = gravity
compute_fluxes = compute_fluxes

HMPPRT_PATH = /home/devres/john/openHMPP/

CFLAGS = -O0 -fno-inline -g -std=c99
CC = gcc
CXX = g++
HMPP_FLAGS13 = -e --nvcc-options -Xptxas=-v,-arch=sm_13 -c --force
HMPP_FLAGS20 = -e --nvcc-options -Xptxas=-v,-arch=sm_15 -c --force
#HMPP_FLAGS += -m

CUDA_ENV_VAR = CUDA_PROFILE=1 
HMPP_ENV_VAR = HMPPRT_LOG_LEVEL=info HMPPRT_NO_FALLBACK=0


all: gravity compute_fluxes

gravity: g.13
	#$(CUDA_ENV_VAR) $(HMPP_ENV_VAR) $(HMPPRT_PATH)g.13
	$(CUDA_ENV_VAR) $(HMPPRT_PATH)g.13

g.13: $(gravity).c
	capsmc $(HMPP_FLAGS13) $(CC) $(CFLAGS) $^ -o $@

compute_fluxes: cf.13
	#$(CUDA_ENV_VAR) $(HMPP_ENV_VAR) $(HMPPRT_PATH)cf.13
	$(CUDA_ENV_VAR) $(HMPPRT_PATH)cf.13

cf.13: $(compute_fluxes).c
	capsmc $(HMPP_FLAGS13) $(CC) $(CFLAGS) $^ -o $@



clean:
	rm -f example?/app_ported_cuda.hmf.o example?/app_ported_cuda.hmf.fatbin hmpprt-quickstart.html example?/app_ported.preproc.ii \
           hmppcg_functions.dpil example1/app_orig.exe example?/app_ported_cuda.hmf example?/app_ported_cuda.hmf.cu example?/app_ported.exe 

gravity = gravity

HMPPRT_PATH = /home/devres/john/openHMPP/

CFLAGS = -O0 -fno-inline -g -std=c99
CC = gcc
CXX = g++
HMPP_FLAGS13 = -e --nvcc-options -Xptxas=-v,-arch=sm_13 -c --force
HMPP_FLAGS20 = -e --nvcc-options -Xptxas=-v,-arch=sm_15 -c --force
#HMPP_FLAGS += -m

CUDA_ENV_VAR = CUDA_PROFILE=1 
HMPP_ENV_VAR = HMPPRT_LOG_LEVEL=info HMPPRT_NO_FALLBACK=0

all: gravity

gravity: $(gravity).exe
	$(CUDA_ENV_VAR) $(HMPPRT_PATH)$(gravity).exe

$(gravity).exe: $(gravity).c
	capsmc $(HMPP_FLAGS13) $(CC) $(CFLAGS) $^ -o $@


clean:
	rm -f example?/app_ported_cuda.hmf.o example?/app_ported_cuda.hmf.fatbin hmpprt-quickstart.html example?/app_ported.preproc.ii \
           hmppcg_functions.dpil example1/app_orig.exe example?/app_ported_cuda.hmf example?/app_ported_cuda.hmf.cu example?/app_ported.exe 

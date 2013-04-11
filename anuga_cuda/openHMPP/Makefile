
CC = gcc
INCLUDES = -I"/usr/lib/python2.7/dist-packages/numpy/core/include"  -I"/usr/include/python2.7" -I"/home/u5044856/anuga_core/source/anuga/utilities/"

CXXFLAGS += -O0 -Wall -fno-inline -g -std=c99 -fPIC ${INCLUDES}

HMPP_FLAGS13 = -e --nvcc-options -Xptxas=-v,-arch=sm_13 -c --force  

LIBS = -lm -lpython2.7 

.PHONY : clean

CUDA_ENV_VAR = CUDA_PROFILE=1
HMPP_ENV_VAR = HMPPRT_LOG_LEVEL=info HMPPRT_NO_FALLBACK=0



#hmpp_python_glue.so : hmpp_python_glue.o gravity.o
#	capsmc ${HMPP_FLAGS13} ${CC} -shared $^ -o $@ ${LIBS}
	#${CC} -shared $< -o $@ ${LIBS}


#hmpp_python_glue.o : hmpp_python_glue.c
	#${CC} ${CXXFLAGS} -c $< -o $@
	#capsmc ${HMPP_FLAGS13} ${CC} ${CXXFLAGS} -c $^ -o $@


hmpp_evolve.13 : hmpp_evolve.o gravity.o compute_fluxes.o
	capsmc ${HMPP_FLAGS13} ${CC} $^ -o $@ ${LIBS}


%.o : %.c   
	capsmc ${HMPP_FLAGS13} ${CC} ${CXXFLAGS} -c $^ -o $@
	
    
#hmpp_python_glue.so : hmpp_python_glue.o gravity.o
#	capsmc ${HMPP_FLAGS13} ${CC} -shared $^ -o $@ ${LIBS}
#
#gravity.o : gravity.c
#	capsmc ${HMPP_FLAGS13} ${CC} ${CXXFLAGS} -c $< -o $@ ${LIBS}
#
#%.o : %.c
#	${CC} ${CXXFLAGS} -c $< -o $@
#




clean : 
	rm -f *.o
	rm -f *.so 


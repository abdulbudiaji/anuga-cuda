
CC = gcc
INCLUDES = -I"/usr/lib/python2.7/dist-packages/numpy/core/include"  -I"/usr/include/python2.7" -I"/home/u5044856/anuga_core/source/anuga/utilities/"

CXXFLAGS += -O3 -Wall -fno-inline -g -std=c99 -fPIC ${INCLUDES}

HMPP_FLAGS13 = -e --nvcc-options -Xptxas=-v,-arch=sm_13 -c --force  

LIBS = -lm -lpython2.7 -L/usr/lib/python2.7

.PHONY : clean

CUDA_ENV_VAR = CUDA_PROFILE=1
HMPP_ENV_VAR = HMPPRT_LOG_LEVEL=info HMPPRT_NO_FALLBACK=0



hmpp_python_glue.so : hmpp_python_glue.o gravity.o
	#capsmc ${HMPP_FLAGS13} ${CC} $< -o $@ ${LIBS}
	${CC} -shared $< -o $@ ${LIBS}

#hmpp_python_glue.o : hmpp_python_glue.c
#	${CC} ${CXXFLAGS} -c $< -o $@

#main : main.o gravity.o
#	capsmc ${HMPP_FLAGS13} ${CC} $< -o $@ ${LIBS}


%.o : %.c   
	#capsmc ${HMPP_FLAGS13} ${CC} ${CXXFLAGS} -c $^ -o $@
	${CC} ${CXXFLAGS} -c $^ -o $@
	
    






clean : 
	rm -f *.o
	rm -f *.so 


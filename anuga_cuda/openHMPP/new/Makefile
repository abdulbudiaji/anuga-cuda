
CC = gcc
INCLUDES = -I"/usr/lib/python2.7/dist-packages/numpy/core/include"  -I"/usr/include/python2.7" -I"/home/u5044856/anuga_core/source/anuga/utilities/"

CXXFLAGS += -O0 -Wall -fno-inline -g -std=c99 -fPIC ${INCLUDES}

HMPP_FLAGS13 = -e --nvcc-options -Xptxas=-v,-arch=sm_13 -c --force  

LIBS = -lm -lpython2.7 

.PHONY : clean

CUDA_ENV_VAR = CUDA_PROFILE=1
HMPP_ENV_VAR = HMPPRT_LOG_LEVEL=info HMPPRT_NO_FALLBACK=0



hmpp_python_glue.so : hmpp_python_glue.o evolve.o gravity.o compute_fluxes.o extrapolate_first_order.o extrapolate_second_order_sw.o extrapolate_second_order_and_limit_by_vertex_or_edge.o swb2_domain_ext.o balance_deep_and_shallow.o protect.o set_boundary.o update_centroids_of_velocities_and_height.o update.o interpolate_from_vertices_to_edges.o manning_friction.o saxpy_centroid_values.o evaluate_segment.o
	capsmc ${HMPP_FLAGS13} ${CC} -shared $^ -o $@ ${LIBS}


%.o : %.c
	capsmc ${HMPP_FLAGS13} ${CC} ${CXXFLAGS} -c $^ -o $@


#gravity.o : gravity.c
#	capsmc ${HMPP_FLAGS13} ${CC} ${CXXFLAGS} -c $^ -o $@


hmpp_python_glue.o : hmpp_python_glue.c sw_domain_fun.o
	${CC} ${CXXFLAGS} -c $^ -o $@

sw_domain_fun.o : sw_domain_fun.h  
	${CC} ${CXXFLAGS} -c $^ -o $@

#%.o : %.c   
#	${CC} ${CXXFLAGS} -c $^ -o $@




clean : 
	rm -f *.o *.so *.hmc *.cu

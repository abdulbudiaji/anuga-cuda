// cudafun.cu -- 
// CUDA memory allocation & shallow water kernel routines
// To compile (Cuda 3.2):
// 	nvcc -c --gpu_architecture sm_13 -I${CUDA_INSTALL_PATH}/include 
//       -Xcompiler -fpic
//
// ! Incomplete kernel call !
//
// Matthias Griessinger, University of Erlangen, 2011.


#include <stdio.h>
//#include <cuda_runtime_api.h>
#ifndef _CUDA_MACROS_H_
#define _CUDA_MACROS_H_

#define safecall(call) do{\
  cudaError_t err = call ;\
  if( cudaSuccess != err ){\
    fprintf(stdout, "cuda error at %s:%d, %s\n",\
      __FILE__, __LINE__, cudaGetErrorString(err));\
    fflush(stdout);\
  }\
  } while(0)

#define BLOCKDIM 96

#ifdef _SHARED_WRITEBACK_3_
#define SHARED_MEM_MULT 3
#else
#define SHARED_MEM_MULT 1
#endif

#endif


/* *********** DEVICE SELECTION ************************* */

extern "C" void getDeviceInfo( int rank, int size, const char* hostname) {
  /* Print device information */
  int deviceCount, device;
  cudaDeviceProp deviceProp;

  cudaGetDeviceCount(&deviceCount);

  if ( 0 == rank ) {
    printf("## rank %i/%i on %s --\t Device Test: No. Cards: %d\n",
      rank, size-1, hostname, deviceCount);
    for( device = 0; device < deviceCount; ++device) {
      cudaGetDeviceProperties(&deviceProp, device);
      printf("## rank %i/%i on %s --\t Device %d: %s\n",
        rank, size-1, hostname, device, deviceProp.name);
    }
  }
}


extern "C" int selectDevice( int rank, int size, const char* hostname ) {
  /* Select GPU device (for multiple cards);
     call before any GPU memory/kernel calls,
     otherwise no effect (default card selected)
  */
  int deviceCount, takedevice, device;
  cudaDeviceProp deviceProp;

  cudaGetDeviceCount(&deviceCount);

  takedevice = rank%deviceCount;
  cudaSetDevice(takedevice);
  cudaGetDevice(&device);
  cudaGetDeviceProperties(&deviceProp, device);

  printf("rank %i/%i on %s --\t Selecting Device %d: %s\n",
    rank, size, hostname, device, deviceProp.name);

  return device;
}


/* *********** KERNEL LAUNCH PARAMETERS ***************** */

typedef struct {
  int gridDim;
  int blockDim;
} KERNEL_LAUNCHER;

KERNEL_LAUNCHER _launcher_;

extern "C" void setKernelDims( const int gridDim, const int blockDim ) {
  _launcher_.gridDim  = gridDim;
  _launcher_.blockDim = blockDim;
}

extern "C" void printKernelDims() {
	printf(" kernel dims: %i x %i\n", 
		_launcher_.gridDim, _launcher_.blockDim);
   fflush(stdout);
}


/* *********** CUDA MEMORY **************************** */

extern "C" void* allocDeviceMemory( size_t bytesize ) {
    char* mem = NULL;
    safecall(cudaMalloc( (void**)&mem, bytesize ));
	 fprintf(stdout,"allocDevice: allocating %lu bytes at %p\n", bytesize, mem);fflush(stdout);

    return (void*)mem;
}

extern "C" void* allocHostMemory( size_t bytesize ) {
    /* returns aligned CPU memory for faster transfer */
    char* mem = NULL;
    safecall(cudaHostAlloc( (void**)&mem, bytesize, 0 ));

    return (void*)mem;
}


extern "C" void copyDeviceToHost( void* hostmem, void* devicemem, size_t bytesize ) {
  /* copy bytesize bytes from GPU to host */
  safecall(cudaMemcpy( hostmem, devicemem, bytesize, cudaMemcpyDeviceToHost ));
}

extern "C" void copyHostToDevice( void* devmem, void* hostmem, size_t bytesize ) {
  /* copy bytesize bytes from host to GPU */
  safecall(cudaMemcpy( devmem, hostmem, bytesize, cudaMemcpyHostToDevice ));
}


extern "C" void freeDeviceMemory( void* mem ) {
	 fprintf(stdout,"freeDevice: freeing at %p\n", mem);fflush(stdout);
    safecall(cudaFree( mem ));
}

extern "C" void freeHostMemory( void* mem ) {
    safecall(cudaFreeHost( mem ));
}


extern "C" void dummy( ) {
    fprintf(stdout, "dummy 1 at %s:%d\n",__FILE__, __LINE__);
    fflush(stdout);
	double* mem;
	//double* mem = (double*)allocDeviceMemory( 128*sizeof(double) );
	cudaMalloc( (void**)&mem, 128*sizeof(double));

    fprintf(stdout, "dummy 2 at %s:%d\n",__FILE__, __LINE__);
    fflush(stdout);
	cudaFree(mem);
    fprintf(stdout, "dummy 3 at %s:%d\n",__FILE__, __LINE__);
    fflush(stdout);
}

/* *********** GPU KERNELS ************************** */


typedef struct {
	double edgeflux_s;
	double edgeflux_x;
	double edgeflux_y;
	double max_speed;
} Fluxes;

typedef struct {
	double u;
	double uh;
	double h;
} Velocity;

__global__ void __set_to_default__(double* edge, 
												size_t N, double def) {
	/* set input array edge of length N to value def */

	size_t k;
	for( k = blockIdx.x * blockDim.x + threadIdx.x; k < N; k += gridDim.x * blockDim.x ) {
		edge[k] = def;
 	}
}

__global__ void __set_arrays_to_default__(double* edge, 
												double* xmom, 
												double* ymom, 
												size_t N, double def) {
	/* set input arrays edge, xmom, ymom of length N to value def */

	size_t k;
	for( k = blockIdx.x * blockDim.x + threadIdx.x; k < N; k += gridDim.x * blockDim.x ) {
		edge[k] = def;
		xmom[k] = def;
		ymom[k] = def;
 	}
}


__global__ void __compute_time_step__( const long* tri_full_flag, 
							const double* max_speed_array, 
							const double* radii,
							const size_t number_of_elements,
							const double epsilon,
							const double time0,
							double* times_out ) {
    /* Computes minimal timestep in each triangle k, and finds 
       minimal timestep in each block of threads.

       Output is written to times_out, final reduction must
       be performed on CPU.*/
	
    // shared memory size defined at kernel launch,
    // set according to blockDim
	 extern __shared__ double sdata[];
	
	 unsigned int tid = threadIdx.x;

    // initialize thread with default (previous) timestep
    double mytime = time0;
   
    // For all triangles 
	 for( size_t k = blockIdx.x * blockDim.x + threadIdx.x; k < number_of_elements; k += gridDim.x*blockDim.x) {
        if( 1 == tri_full_flag[k] && max_speed_array[k] > epsilon ) {
				mytime = fmin( mytime, radii[k] / max_speed_array[k] );
		  }
	 }
    
    // each thread in block writes to shared memory
    sdata[tid] = mytime;
    __syncthreads();

    // Reduce to one value per thread block by successively 
    // comparing value pairs; in first sweep, first half of 
    // threads compares first half of values to second half and
    // writes min to first half; 2nd sweep, first fourth compares
    // first and second fourth a.s.o.
    for(unsigned int s=blockDim.x/2; s>0; s>>=1) {
		if(tid < s) {
			sdata[tid] = fmin( sdata[tid + s], sdata[tid]);
		}
		__syncthreads();
	}

    // Lead thread writes min for this block to global mem
    if (tid == 0) {
		times_out[blockIdx.x] = sdata[0];
	}
}


__device__ double2 __rotate__(double q1, double q2, double n1, double n2) {
  /*Rotate the momentum component q (q1, q2)
    from x,y coordinates to coordinates based on normal vector (n1, n2).

    Result is returned in array 3x1 r
    To rotate in opposite direction, call rotate with (q, n1, -n2)

    Contents of q are changed by this function */

  double2 q;

  // Rotate
  q.x =  n1*q1 + n2*q2;
  q.y = -n2*q1 + n1*q2;

  return q;
}

__device__ Velocity __compute_speed__(
              double uh, 
              double h,
              const double epsilon,
              const double h0,
              const double limiting_threshold) {
  
  Velocity result;

  if (h < limiting_threshold) {   
    // Apply limiting of speeds according to the ANUGA manual
    if (h < epsilon) {
      h = 0.0;  // Could have been negative
      result.u = 0.0;
    } else {
      result.u = uh/(h + h0/ h);    
    }
  

    // Adjust momentum to be consistent with speed
    uh = result.u * h;
  } else {
    // We are in deep water - no need for limiting
    result.u = uh/ h;
  }
  
  result.uh = uh; 
  result.h = h; 
  return result;
}

__device__ Fluxes __flux_function_central__(double2 stage, double2 xmom,
					 double2 ymom, double2 z_lr,
					 double2 n,
					 const double epsilon, 
					 const double h0,
					 const double limiting_threshold,
					 const double g) {
					 //double *edgeflux,
					 //double *max_speed) {

  /*Compute fluxes between volumes for the shallow water wave equation
    cast in terms of the 'stage', w = h+z using
    the 'central scheme' as described in

    Kurganov, Noelle, Petrova. 'Semidiscrete Central-Upwind Schemes For
    Hyperbolic Conservation Laws and Hamilton-Jacobi Equations'.
    Siam J. Sci. Comput. Vol. 23, No. 3, pp. 707-740.

    The implemented formula is given in equation (3.15) on page 714
  */
  double2 tmp;
  double h_left, uh_left, vh_left, u_left;
  double h_right, uh_right, vh_right, u_right;
  double s_min, s_max, soundspeed_left, soundspeed_right;
  double denom, inverse_denominator, z;

  // Cuda doesn't do static arrays
  double fs_l, fs_r;
  double2 fv_l, fv_r;
  Velocity velos;
  Fluxes fluxes;

  // Align x- and y-momentum with x-axis
  // Do not be confused: xmom.x and xmom.y are 
  // left and right momenti in x direction
  tmp = __rotate__( xmom.x, ymom.x, n.x, n.y);
  xmom.x = tmp.x; ymom.x = tmp.y;
  tmp = __rotate__( xmom.y, ymom.y, n.x, n.y);
  xmom.y = tmp.x; ymom.y = tmp.y;

  z = 0.5*(z_lr.x + z_lr.y); // Average elevation values. 
                              // Even though this will nominally allow 
                  // for discontinuities in the elevation data,
                  // there is currently no numerical support for
                  // this so results may be strange near
                  // jumps in the bed.

  // Compute speeds in x-direction
  h_left = stage.x - z;

  uh_left = xmom.x; // q_left_rotated[1];
  velos = __compute_speed__(uh_left, h_left, 
              epsilon, h0, limiting_threshold);
  u_left = velos.u;
  uh_left = velos.uh;
  h_left = velos.h;

  h_right = stage.y - z;
  uh_right = xmom.y; // q_right_rotated[1];
  velos = __compute_speed__(uh_right, h_right, 
               epsilon, h0, limiting_threshold);
  u_right = velos.u;
  uh_right = velos.uh;
  h_right = velos.h;

  // Momentum in y-direction
  vh_left  = ymom.x; //q_left_rotated[2];
  vh_right = ymom.y; //q_right_rotated[2];

  // Limit y-momentum if necessary 
  // Leaving this out, improves speed significantly (Ole 27/5/2009)
  // All validation tests pass, so do we really need it anymore? 
  velos = __compute_speed__(vh_left, h_left, 
         epsilon, h0, limiting_threshold);
  vh_left = velos.uh;
  h_left = velos.h;
  velos = __compute_speed__(vh_right, h_right, 
         epsilon, h0, limiting_threshold);
  vh_right = velos.uh;
  h_right = velos.h;

  // Maximal and minimal wave speeds
  soundspeed_left  = sqrt(g*h_left);
  soundspeed_right = sqrt(g*h_right);  

  s_max = fmax(u_left + soundspeed_left, u_right + soundspeed_right);
  if (s_max < 0.0) 
  {
    s_max = 0.0;
  }

  s_min = fmin(u_left - soundspeed_left, u_right - soundspeed_right);
  if (s_min > 0.0)
  {
    s_min = 0.0;
  }
  
  // Flux formulas
  fs_l 	= u_left*h_left;
  fv_l.x = u_left*uh_left + 0.5*g*h_left*h_left;
  fv_l.y = u_left*vh_left;

  fs_r   = u_right*h_right;
  fv_r.x = u_right*uh_right + 0.5*g*h_right*h_right;
  fv_r.y = u_right*vh_right;

  // Flux computation
  denom = s_max - s_min;
  if (denom < epsilon) 
  { // FIXME (Ole): Try using h0 here
	 fluxes.edgeflux_s = 0.0;
	 fluxes.edgeflux_x = 0.0;
	 fluxes.edgeflux_y = 0.0;
    fluxes.max_speed = 0.0;
  } 
  else 
  {
    inverse_denominator = 1.0/denom;

      fluxes.edgeflux_s = s_max*fs_l - s_min*fs_r;
      fluxes.edgeflux_s += s_max*s_min*(stage.y - stage.x);
      fluxes.edgeflux_s *= inverse_denominator;
      
      fluxes.edgeflux_x = s_max*fv_l.x - s_min*fv_r.x;
      fluxes.edgeflux_x += s_max*s_min*(xmom.y - xmom.x);
      fluxes.edgeflux_x *= inverse_denominator;
      
      fluxes.edgeflux_y = s_max*fv_l.y - s_min*fv_r.y;
      fluxes.edgeflux_y += s_max*s_min*(ymom.y - ymom.x);
      fluxes.edgeflux_y *= inverse_denominator;

    	// Maximal wavespeed
    	fluxes.max_speed = fmax(fabs(s_max), fabs(s_min));

    // Rotate back
    tmp = __rotate__( fluxes.edgeflux_x, fluxes.edgeflux_y, n.x, -n.y);
	 fluxes.edgeflux_x = tmp.x; fluxes.edgeflux_y = tmp.y;
  }
 
  return fluxes; 
}


__global__ void __compute_fluxes_central_kernel__(
		const int number_of_elements,
        double timestep,
        const double epsilon,
        const double H0,
        const double g,
        const long* neighbours,
        const long* neighbour_edges,
        const double* normals,
        double* edgelengths,
        const double* areas,
        const double* stage_edge_values,
        const double* xmom_edge_values,
        const double* ymom_edge_values,
        const double* bed_edge_values,
        const double* stage_boundary_values,
        const double* xmom_boundary_values,
        const double* ymom_boundary_values,
        double* stage_explicit_update,
        double* xmom_explicit_update,
        double* ymom_explicit_update,
        double* max_speed_array,
        const int optimise_dry_cells) {
    /* __global__: called by CPU and executed on GPU;
                   must return void, scalar variables, structs (and scalar members)
                   are copied automatically, arrays must be allocated on device.
    */

    // Local variables
    double length; //, zl, zr;
    double h0 = H0*H0; // This ensures a good balance when h approaches H0.

    double limiting_threshold = 10 * H0; // Avoid applying limiter below this
    // threshold for performance reasons.
    // See ANUGA manual under flux limiting
    int i, k, m, n;
    int ki, nm = 0, ki2; // Index shorthands

    Fluxes fluxes;
	
    double2 stage, xmom, ymom, z_lr, normvec;
  
    // Shared memory for explicit update quantity reduction
    extern __shared__ double update_shared[]; // empty [] array:(byte)size defined by kernel call
	 //__shared__ double update_shared[BLOCKDIM*SHARED_MEM_MULT]; // OR static size (not both)

    // For all edges; 
    for( ki = blockIdx.x * blockDim.x + threadIdx.x; ki < 3*number_of_elements; ki += gridDim.x * blockDim.x ) {

            // Get left hand side values from triangle k, edge i
            stage.x = stage_edge_values[ki];
            xmom.x = xmom_edge_values[ki];
            ymom.x = ymom_edge_values[ki];
            z_lr.x = bed_edge_values[ki];

            // Get right hand side values either from neighbouring triangle
            // or from boundary array (Quantities at neighbour on nearest face).
            n = neighbours[ki];
            if (n < 0) {
                // Neighbour is a boundary condition
                m = -n - 1; // Convert negative flag to boundary index

                // Bad access order consider binding boundary_values to Texture cache
                stage.y = stage_boundary_values[m];
                xmom.y = xmom_boundary_values[m];
                ymom.y = ymom_boundary_values[m];
                z_lr.y = z_lr.x; // Extend bed elevation to boundary
            }
            else {
                // Neighbour is a real triangle
                m = neighbour_edges[ki];
                nm = n * 3 + m; // Linear index (triangle n, edge m)

                // Again, bind to Texture cache
                stage.y = stage_edge_values[nm];
                xmom.y = xmom_edge_values[nm];
                ymom.y = ymom_edge_values[nm];
                z_lr.y = bed_edge_values[nm];
            }

            // Now we have values for this edge - both from left and right side.

            /*if (optimise_dry_cells) {
                // Check if flux calculation is necessary across this edge
                // This check will exclude dry cells.
                // This will also optimise cases where zl != zr as
                // long as both are dry

                if (fabs(ql[0] - zl) < epsilon &&
                        fabs(qr[0] - zr) < epsilon) {
                    // Cell boundary is dry

                    already_computed_flux[ki] = call; // #k Done
                    if (n >= 0) {
                        already_computed_flux[nm] = call; // #n Done
                    }

                    max_speed = 0.0;
                    continue;
                }
            }*/


            /*if (fabs(zl-zr)>1.0e-10) {TODO:
                report_python_error(AT,"Discontinuous Elevation");
                return 0.0;
            }*/
            
            // Outward pointing normal vector (domain.normals[k, 2*i:2*i+2])
            ki2 = 2 * ki; //k*6 + i*2
            // Bad access order, use Texture or shared mem
				normvec.x = normals[ki2];
            normvec.y = normals[ki2+1];

            // Edge flux computation (triangle k, edge i);
				// TODO: subroutine causes unspecified launch failure
            //fluxes = __flux_function_central__( stage, xmom, ymom, z_lr, normvec,
            //        epsilon, h0, limiting_threshold, g);
  

            // Multiply edgeflux by edgelength
            length = edgelengths[ki];
            fluxes.edgeflux_s *= length;
            fluxes.edgeflux_x *= length;
            fluxes.edgeflux_y *= length;

            
				// Use shared memory to accumulate flux for one triangle;
				// requires blockDim.x to be multiple of 3 (should also be multiple of 32)

				#ifdef _SHARED_WRITEBACK_3_

					// Accumulate all update arrays in one sweep;
					// if shared memory is large enough

            	update_shared[threadIdx.x] = fluxes.edgeflux_s;
            	update_shared[threadIdx.x+blockDim.x] = fluxes.edgeflux_x;
            	update_shared[threadIdx.x+2*blockDim.x] = fluxes.edgeflux_y;
					__syncthreads();
					

					// Each third of the threads in block write back an update array
					// with contiguous access, each thread sums fluxes in a triangle

					if( threadIdx.x < blockDim.x/3) {
                  // Calculate contiguous index
						i = threadIdx.x;
						k = blockIdx.x*(blockDim.x/3) + i;
						stage_explicit_update[k] = ( update_shared[3*i]
															+ update_shared[3*i+1]
															+ update_shared[3*i+2]) / areas[k];
					} else if( threadIdx.x < 2*(blockDim.x/3) ) {
						i = threadIdx.x - (blockDim.x/3);
						k = blockIdx.x*(blockDim.x/3) + i;
						xmom_explicit_update[k] = ( update_shared[blockDim.x+3*i]
															+ update_shared[blockDim.x+3*i+1]
															+ update_shared[blockDim.x+3*i+2]) / areas[k];
					} else {
						i = threadIdx.x - 2*(blockDim.x/3);
						k = blockIdx.x*(blockDim.x/3) + i;
						ymom_explicit_update[k] = ( update_shared[2*blockDim.x+3*i]
															+ update_shared[2*blockDim.x+3*i+1]
															+ update_shared[2*blockDim.x+3*i+2] ) / areas[k];

					}
					__syncthreads();
            
            #else   
            
					// Write each update array back by itself;
					// only first third of threads busy
	            update_shared[threadIdx.x]= fluxes.edgeflux_s;
					__syncthreads();
               if( threadIdx.x < blockDim.x/3 ) {
						i = threadIdx.x;
						k = blockIdx.x*(blockDim.x/3) + i;
               	stage_explicit_update[k] = ( update_shared[3*i] 
												+ update_shared[3*i+1]
												+ update_shared[3*i+2] ) / areas[k];
					}
					__syncthreads();
            
				
               update_shared[threadIdx.x] = fluxes.edgeflux_x;
					__syncthreads();
               if( threadIdx.x < blockDim.x/3 ) {
						i = threadIdx.x;
						k = blockIdx.x*(blockDim.x/3) + i;
               	xmom_explicit_update[k] = ( update_shared[3*i] 
												+ update_shared[3*i+1]
												+ update_shared[3*i+2] ) / areas[k];
					}
					__syncthreads();
               
					update_shared[threadIdx.x] = fluxes.edgeflux_y;
					__syncthreads();
               if( threadIdx.x < blockDim.x/3 ) {
						i = threadIdx.x;
						k = blockIdx.x*(blockDim.x/3) + i;
               	ymom_explicit_update[k] = ( update_shared[3*i] 
												+ update_shared[3*i+1]
												+ update_shared[3*i+2] ) / areas[k];
					}
					__syncthreads();
				
				#endif
			   
				// Likewise, get and write maximum speed within triangle
            // update_shared[threadIdx.x] = fluxes.max_speed;
				__syncthreads();
            if( threadIdx.x < blockDim.x/3 ) {
					i = threadIdx.x;
					k = blockIdx.x*(blockDim.x/3) + i;
					fluxes.max_speed = fmax( update_shared[3*i], update_shared[3*i+1] );
					max_speed_array[k] = fmax( fluxes.max_speed, update_shared[3*i+2] );
				}
           
    } // End edge ki

    // computation of timestep in seperate routine because of triangle-wise access
}


/* *********** KERNEL WRAPPER FUNCTIONS ************************** */

extern "C" void _set_to_default( double* edge, double* xmom, double* ymom, size_t N, double def) {
	
	 __set_arrays_to_default__ <<< _launcher_.gridDim, _launcher_.blockDim >>> ( edge, xmom, ymom, N, def);
    safecall(cudaThreadSynchronize());
}

/*extern "C" void _set_to_default( double* edge, size_t N, double def) {
	
	 __set_to_default__ <<< _launcher_.gridDim, _launcher_.blockDim >>> ( edge, N, def);
    safecall(cudaThreadSynchronize());
}*/

extern "C" double _compute_fluxes_central(
		int number_of_elements,
        double timestep,
        double epsilon,
        double H0,
        double g,
        long* neighbours,
        long* neighbour_edges,
        double* normals,
        double* edgelengths,
        double* radii,
        double* areas,
        long* tri_full_flag,
        double* stage_edge_values,
        double* xmom_edge_values,
        double* ymom_edge_values,
        double* bed_edge_values,
        double* stage_boundary_values,
        double* xmom_boundary_values,
        double* ymom_boundary_values,
        double* stage_explicit_update,
        double* xmom_explicit_update,
        double* ymom_explicit_update,
        double* max_speed_array,
        int optimise_dry_cells) {

    static long call = 1; // Static local variable flagging already computed flux
    int i;
    // Start computation
    call++; // Flag 'id' of flux calculation for this timestep

    // prepare memory for timestep reduction (TODO: (de)allocate only once)
    const size_t reduction_size = _launcher_.gridDim*sizeof(double);
    double* times_out = (double*)allocHostMemory( reduction_size  );
    double* times_out_gpu = (double*)allocDeviceMemory( reduction_size );

    printf("shared mum mult: %i\n", SHARED_MEM_MULT);

    if( 0 != _launcher_.blockDim%3 ) {
		fprintf(stderr,"error: blockDim required to be multiple of 3!\n");
    }	
	 
    __compute_fluxes_central_kernel__ <<< _launcher_.gridDim, _launcher_.blockDim, 
	  											 		_launcher_.blockDim*sizeof(double)*SHARED_MEM_MULT >>> (
    //__compute_fluxes_central_kernel__ <<< _launcher_.gridDim, BLOCKDIM >>> ( // for static shared memory size
		number_of_elements,
        timestep,
        epsilon,
        H0,
        g,
        neighbours,
        neighbour_edges,
        normals,
        edgelengths,
        areas,
        stage_edge_values,
        xmom_edge_values,
        ymom_edge_values,
        bed_edge_values,
        stage_boundary_values,
        xmom_boundary_values,
        ymom_boundary_values,
        stage_explicit_update,
        xmom_explicit_update,
        ymom_explicit_update,
        max_speed_array,
        optimise_dry_cells);
	 

    safecall(cudaThreadSynchronize()); // prevents overlap of kernels
    
	 // Some timestepping debug: (timestep 1.0)
	 //printKernelDims();
	 //_set_to_default( max_speed_array, radii, areas, number_of_elements, 3.3); // 
	
    __compute_time_step__ <<< _launcher_.gridDim, _launcher_.blockDim, _launcher_.blockDim*sizeof(double) >>> (
	                  tri_full_flag, 
							max_speed_array, 
							radii,
							number_of_elements,
							epsilon,
							timestep,
							times_out_gpu );

    safecall(cudaThreadSynchronize());
    
    copyDeviceToHost( times_out, times_out_gpu, reduction_size );
    for( i=0; i < _launcher_.gridDim; ++i) {
		timestep = min( timestep, times_out[i] );
    }
	 //printf("\ntimestep = %f\n",timestep);
	 //fflush(stdout);

    freeDeviceMemory( times_out_gpu );
    freeHostMemory( times_out );

    return timestep;
}

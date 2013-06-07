
#include <stdio.h>

#ifndef __CUDACC__
#include <stdlib.h>
#include <math.h>

#include <hmpprt/Grouplet.h>
#include <hmpprt/HostTypes.h>
#include <hmpprt/Context.h>
#include <hmpprt/CUDAGrid.h>
#include <hmpprt/CUDAModule.h>
#include <hmpprt/DeviceManager.h>
#include <hmpperr/hmpperr.h>

#ifdef _WIN32
#  define CDLT_API __declspec(dllexport)
#else /* ! _WIN32 */
#  define CDLT_API
#endif /* _WIN32 */



#else // ! __CUDACC__

#include <hmpprt/HostTypes.h>
#include <hmpprt/CUDAIntrinsics.h>

extern __shared__ int64_t hmpp_sharedmem[];
#endif // __CUDACC__



#ifndef __CUDACC__

#else


#endif

#define HMPPCG_SIMD_LENGTH 1

# 234 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"

#ifndef __CUDACC__
extern "C" CDLT_API  void _extrapolate_from_gradient(hmpprt::s32 N_1, hmpprt::s32 N2_1, hmpprt::s32 N3_2, hmpprt::s32 N6_1, double* centroids_1, double* centroid_values_1, double* vertex_coordinates_1, double* vertex_values, double* edge_values, double* a, double* b)
;
#endif // __CUDACC__



# 234 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"

#ifndef __CUDACC__
void _extrapolate_from_gradient_internal_1(hmpprt::s32 N_11, hmpprt::s32 N2_11, hmpprt::s32 N3_21, hmpprt::s32 N6_11, double* centroids_11, double* centroid_values_11, double* vertex_coordinates_11, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  vertex_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  edge_values, double* a, double* b)
;
#endif // __CUDACC__



# 234 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"

#ifndef __CUDACC__
static hmpprt::CUDAGrid * _extrapolate_from_gradient_loop1D_1 = 0;
#else

extern "C" __global__ void _extrapolate_from_gradient_loop1D_1(hmpprt::s32 N_2, double* vertex_values_2, double* edge_values_2);
#endif // __CUDACC__




# 234 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"

#ifdef __CUDACC__

extern "C" __global__ void _extrapolate_from_gradient_loop1D_1(hmpprt::s32 N_2, double* vertex_values_2, double* edge_values_2)
{
 # 247 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 hmpprt::s32 k3_1;
 # 247 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 hmpprt::s32 k_1;
 # 257 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 k_1 = (hmpprt::gr_atidf());
 # 257 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 if (k_1 > N_2 - 1)
 {
  # 257 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  goto __hmppcg_label_1;
 }
 # 257 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 k3_1 = 3 * k_1;
 # 273 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 *(vertex_values_2 + k3_1) = (double) 1.0;
 # 274 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 *(vertex_values_2 + (k3_1 + 1)) = (double) 1.0;
 # 275 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 *(vertex_values_2 + (k3_1 + 2)) = (double) 1.0;
 # 278 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 *(edge_values_2 + k3_1) = (double) 0.5 * (*(vertex_values_2 + (k3_1 + 1)) + *(vertex_values_2 + (k3_1 + 2)));
 # 279 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 *(edge_values_2 + (k3_1 + 1)) = (double) 0.5 * (*(vertex_values_2 + (k3_1 + 2)) + *(vertex_values_2 + k3_1));
 # 280 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 *(edge_values_2 + (k3_1 + 2)) = (double) 0.5 * (*(vertex_values_2 + k3_1) + *(vertex_values_2 + (k3_1 + 1)));
 # 234 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 __hmppcg_label_1:;
}
#endif // __CUDACC__



# 234 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"

#ifndef __CUDACC__
void _extrapolate_from_gradient_internal_1(hmpprt::s32 N_11, hmpprt::s32 N2_11, hmpprt::s32 N3_21, hmpprt::s32 N6_11, double* centroids_11, double* centroid_values_11, double* vertex_coordinates_11, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  vertex_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  edge_values, double* a, double* b)
{
 # 234 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 if (N_11 - 1 >= 0)
 {
  hmpprt::CUDAGridCall __hmppcg_call;
  __hmppcg_call.setSizeX((N_11 - 1) / 128 + 1);
  __hmppcg_call.setSizeY(1);
  __hmppcg_call.setBlockSizeX(32);
  __hmppcg_call.setBlockSizeY(4);
  __hmppcg_call.addLocalParameter((hmpprt::s32) (N_11), "N_2");
  __hmppcg_call.addLocalParameter(&vertex_values, 8, "vertex_values_2");
  __hmppcg_call.addLocalParameter(&edge_values, 8, "edge_values_2");
  __hmppcg_call.launch(_extrapolate_from_gradient_loop1D_1, hmpprt::Context::getInstance()->getCUDADevice());
 }
 ;
}
#endif // __CUDACC__



# 234 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"

#ifndef __CUDACC__
extern "C" CDLT_API  void _extrapolate_from_gradient(hmpprt::s32 N_1, hmpprt::s32 N2_1, hmpprt::s32 N3_2, hmpprt::s32 N6_1, double* centroids_1, double* centroid_values_1, double* vertex_coordinates_1, double* vertex_values, double* edge_values, double* a, double* b)
{
 # 1 "<preprocessor>"
 (_extrapolate_from_gradient_internal_1(N_1, N2_1, N3_2, N6_1, centroids_1, centroid_values_1, vertex_coordinates_1, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (vertex_values), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (edge_values), a, b));
}
#endif // __CUDACC__




#ifndef __CUDACC__
extern "C" const char * hmpprt_cuda_get_gpu_code();

static hmpprt::CUDAModule * hmpprt_module = 0;
static int hmpprt_uses = 0;

extern "C" CDLT_API void * hmpprt_init()
{
  try
  {
    if (hmpprt_uses++ == 0)
    {
      hmpprt_module = new hmpprt::CUDAModule(hmpprt_cuda_get_gpu_code());
      _extrapolate_from_gradient_loop1D_1 = new hmpprt::CUDAGrid(hmpprt_module, "_extrapolate_from_gradient_loop1D_1");

    }
    hmpprt::Context::getInstance()->getGrouplet()->setTarget(hmpprt::CUDA);
    hmpprt::Context::getInstance()->getGrouplet()->addSignature("_extrapolate_from_gradient", "prototype _extrapolate_from_gradient(N: s32, N2: s32, N3: s32, N6: s32, centroids: unused ^host double, centroid_values: unused ^host double, vertex_coordinates: unused ^host double, vertex_values: ^cudaglob double, edge_values: ^cudaglob double, a: unused ^host double, b: unused ^host double)");

  }
  catch (hmpperr::Error & e)
  {
    return e.clone();
  }
  catch(...)
  {
    fprintf(stderr,"Unexpected error in hmpprt_init()\n");
    abort();
  }
  return 0;
}
#endif // __CUDACC__

#ifndef __CUDACC__
extern "C" CDLT_API void * hmpprt_fini()
{
  try
  {
    if (--hmpprt_uses == 0)
    {
      delete _extrapolate_from_gradient_loop1D_1;

      delete hmpprt_module;
      hmpprt_module = 0;
    }
  }
  catch (hmpperr::Error & e)
  {
    return e.clone();
  }
  catch(...)
  {
    fprintf(stderr,"Unexpected error in hmpprt_fini()\n");
    abort();
  }
  return 0;
}
#endif // __CUDACC__

// footer

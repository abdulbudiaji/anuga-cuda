
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

# 8 "extrapolate_first_order.c"

#ifndef __CUDACC__
extern "C" CDLT_API  void extrapolate_first_order(hmpprt::s32 N_1, hmpprt::s32 N3_1, double* centroid_values_1, double* edge_values_1, double* vertex_values_2)
;
#endif // __CUDACC__



# 8 "extrapolate_first_order.c"

#ifndef __CUDACC__
void extrapolate_first_order_internal_1(hmpprt::s32 N_11, hmpprt::s32 N3_11, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  centroid_values_11, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  edge_values_11, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  vertex_values_21)
;
#endif // __CUDACC__



# 8 "extrapolate_first_order.c"

#ifndef __CUDACC__
static hmpprt::CUDAGrid * extrapolate_first_order_loop1D_1 = 0;
#else

extern "C" __global__ void extrapolate_first_order_loop1D_1(hmpprt::s32 N_2, double* centroid_values_2, double* edge_values_2, double* vertex_values_1);
#endif // __CUDACC__




# 8 "extrapolate_first_order.c"

#ifdef __CUDACC__

extern "C" __global__ void extrapolate_first_order_loop1D_1(hmpprt::s32 N_2, double* centroid_values_2, double* edge_values_2, double* vertex_values_1)
{
 # 15 "extrapolate_first_order.c"
 hmpprt::s32 k3_1;
 # 15 "extrapolate_first_order.c"
 hmpprt::s32 k_1;
 # 23 "extrapolate_first_order.c"
 k_1 = (hmpprt::gr_atidf());
 # 23 "extrapolate_first_order.c"
 if (k_1 > N_2 - 1)
 {
  # 23 "extrapolate_first_order.c"
  goto __hmppcg_label_1;
 }
 # 23 "extrapolate_first_order.c"
 k3_1 = k_1 * 3;
 # 24 "extrapolate_first_order.c"
 *(edge_values_2 + k3_1) = *(centroid_values_2 + k_1);
 # 25 "extrapolate_first_order.c"
 *(edge_values_2 + (k3_1 + 1)) = *(centroid_values_2 + k_1);
 # 26 "extrapolate_first_order.c"
 *(edge_values_2 + (k3_1 + 2)) = *(centroid_values_2 + k_1);
 # 28 "extrapolate_first_order.c"
 *(vertex_values_1 + k3_1) = *(centroid_values_2 + k_1);
 # 29 "extrapolate_first_order.c"
 *(vertex_values_1 + (k3_1 + 1)) = *(centroid_values_2 + k_1);
 # 30 "extrapolate_first_order.c"
 *(vertex_values_1 + (k3_1 + 2)) = *(centroid_values_2 + k_1);
 # 8 "extrapolate_first_order.c"
 __hmppcg_label_1:;
}
#endif // __CUDACC__



# 8 "extrapolate_first_order.c"

#ifndef __CUDACC__
void extrapolate_first_order_internal_1(hmpprt::s32 N_11, hmpprt::s32 N3_11, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  centroid_values_11, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  edge_values_11, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  vertex_values_21)
{
 # 8 "extrapolate_first_order.c"
 if (N_11 - 1 >= 0)
 {
  hmpprt::CUDAGridCall __hmppcg_call;
  __hmppcg_call.setSizeX((N_11 - 1) / 128 + 1);
  __hmppcg_call.setSizeY(1);
  __hmppcg_call.setBlockSizeX(32);
  __hmppcg_call.setBlockSizeY(4);
  __hmppcg_call.addLocalParameter((hmpprt::s32) (N_11), "N_2");
  __hmppcg_call.addLocalParameter(&centroid_values_11, 8, "centroid_values_2");
  __hmppcg_call.addLocalParameter(&edge_values_11, 8, "edge_values_2");
  __hmppcg_call.addLocalParameter(&vertex_values_21, 8, "vertex_values_1");
  __hmppcg_call.launch(extrapolate_first_order_loop1D_1, hmpprt::Context::getInstance()->getCUDADevice());
 }
 ;
}
#endif // __CUDACC__



# 8 "extrapolate_first_order.c"

#ifndef __CUDACC__
extern "C" CDLT_API  void extrapolate_first_order(hmpprt::s32 N_1, hmpprt::s32 N3_1, double* centroid_values_1, double* edge_values_1, double* vertex_values_2)
{
 # 1 "<preprocessor>"
 (extrapolate_first_order_internal_1(N_1, N3_1, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (centroid_values_1), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (edge_values_1), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (vertex_values_2)));
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
      extrapolate_first_order_loop1D_1 = new hmpprt::CUDAGrid(hmpprt_module, "extrapolate_first_order_loop1D_1");

    }
    hmpprt::Context::getInstance()->getGrouplet()->setTarget(hmpprt::CUDA);
    hmpprt::Context::getInstance()->getGrouplet()->addSignature("extrapolate_first_order", "prototype extrapolate_first_order(N: s32, N3: s32, centroid_values: ^cudaglob double, edge_values: ^cudaglob double, vertex_values: ^cudaglob double)");

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
      delete extrapolate_first_order_loop1D_1;

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

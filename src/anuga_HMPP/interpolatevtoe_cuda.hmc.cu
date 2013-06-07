
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

# 8 "interpolate_from_vertices_to_edges.c"

#ifndef __CUDACC__
extern "C" CDLT_API  void interpolate_from_vertices_to_edges(hmpprt::s32 N_2, hmpprt::s32 N3_2, double* vertex_values_2, double* edge_values_2)
;
#endif // __CUDACC__



# 8 "interpolate_from_vertices_to_edges.c"

#ifndef __CUDACC__
void interpolate_from_vertices_to_edges_internal_1(hmpprt::s32 N_21, hmpprt::s32 N3_21, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  vertex_values_21, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  edge_values_21)
;
#endif // __CUDACC__



# 8 "interpolate_from_vertices_to_edges.c"

#ifndef __CUDACC__
static hmpprt::CUDAGrid * interpolate_from_vertices_to_edges_loop1D_1 = 0;
#else

extern "C" __global__ void interpolate_from_vertices_to_edges_loop1D_1(hmpprt::s32 N_1, double* vertex_values_1, double* edge_values_1);
#endif // __CUDACC__




# 8 "interpolate_from_vertices_to_edges.c"

#ifdef __CUDACC__

extern "C" __global__ void interpolate_from_vertices_to_edges_loop1D_1(hmpprt::s32 N_1, double* vertex_values_1, double* edge_values_1)
{
 # 14 "interpolate_from_vertices_to_edges.c"
 hmpprt::s32 k3_1;
 # 14 "interpolate_from_vertices_to_edges.c"
 double q1_1;
 # 14 "interpolate_from_vertices_to_edges.c"
 double q2_1;
 # 14 "interpolate_from_vertices_to_edges.c"
 double q0_1;
 # 14 "interpolate_from_vertices_to_edges.c"
 hmpprt::s32 k_1;
 # 23 "interpolate_from_vertices_to_edges.c"
 k_1 = (hmpprt::gr_atidf());
 # 23 "interpolate_from_vertices_to_edges.c"
 if (k_1 > N_1 - 1)
 {
  # 23 "interpolate_from_vertices_to_edges.c"
  goto __hmppcg_label_1;
 }
 # 23 "interpolate_from_vertices_to_edges.c"
 k3_1 = k_1 * 3;
 # 24 "interpolate_from_vertices_to_edges.c"
 q0_1 = *(vertex_values_1 + k3_1);
 # 25 "interpolate_from_vertices_to_edges.c"
 q1_1 = *(vertex_values_1 + (k3_1 + 1));
 # 26 "interpolate_from_vertices_to_edges.c"
 q2_1 = *(vertex_values_1 + (k3_1 + 2));
 # 28 "interpolate_from_vertices_to_edges.c"
 *(edge_values_1 + k3_1) = (double) 0.5 * (q1_1 + q2_1);
 # 29 "interpolate_from_vertices_to_edges.c"
 *(edge_values_1 + (k3_1 + 1)) = (double) 0.5 * (q0_1 + q2_1);
 # 30 "interpolate_from_vertices_to_edges.c"
 *(edge_values_1 + (k3_1 + 2)) = (double) 0.5 * (q0_1 + q1_1);
 # 8 "interpolate_from_vertices_to_edges.c"
 __hmppcg_label_1:;
}
#endif // __CUDACC__



# 8 "interpolate_from_vertices_to_edges.c"

#ifndef __CUDACC__
void interpolate_from_vertices_to_edges_internal_1(hmpprt::s32 N_21, hmpprt::s32 N3_21, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  vertex_values_21, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  edge_values_21)
{
 # 8 "interpolate_from_vertices_to_edges.c"
 if (N_21 - 1 >= 0)
 {
  hmpprt::CUDAGridCall __hmppcg_call;
  __hmppcg_call.setSizeX((N_21 - 1) / 128 + 1);
  __hmppcg_call.setSizeY(1);
  __hmppcg_call.setBlockSizeX(32);
  __hmppcg_call.setBlockSizeY(4);
  __hmppcg_call.addLocalParameter((hmpprt::s32) (N_21), "N_1");
  __hmppcg_call.addLocalParameter(&vertex_values_21, 8, "vertex_values_1");
  __hmppcg_call.addLocalParameter(&edge_values_21, 8, "edge_values_1");
  __hmppcg_call.launch(interpolate_from_vertices_to_edges_loop1D_1, hmpprt::Context::getInstance()->getCUDADevice());
 }
 ;
}
#endif // __CUDACC__



# 8 "interpolate_from_vertices_to_edges.c"

#ifndef __CUDACC__
extern "C" CDLT_API  void interpolate_from_vertices_to_edges(hmpprt::s32 N_2, hmpprt::s32 N3_2, double* vertex_values_2, double* edge_values_2)
{
 # 1 "<preprocessor>"
 (interpolate_from_vertices_to_edges_internal_1(N_2, N3_2, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (vertex_values_2), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (edge_values_2)));
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
      interpolate_from_vertices_to_edges_loop1D_1 = new hmpprt::CUDAGrid(hmpprt_module, "interpolate_from_vertices_to_edges_loop1D_1");

    }
    hmpprt::Context::getInstance()->getGrouplet()->setTarget(hmpprt::CUDA);
    hmpprt::Context::getInstance()->getGrouplet()->addSignature("interpolate_from_vertices_to_edges", "prototype interpolate_from_vertices_to_edges(N: s32, N3: s32, vertex_values: ^cudaglob double, edge_values: ^cudaglob double)");

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
      delete interpolate_from_vertices_to_edges_loop1D_1;

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

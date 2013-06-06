
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

# 8 "set_boundary.c"

#ifndef __CUDACC__
extern "C" CDLT_API  void set_boundary_values_from_edges(hmpprt::s32 Nb, hmpprt::s32 N3, hmpprt::s64* vol_id, hmpprt::s64* edge_id, double* boundary_values_1, double* edge_values_1)
;
#endif // __CUDACC__



# 8 "set_boundary.c"

#ifndef __CUDACC__
void set_boundary_values_from_edges_internal_1(hmpprt::s32 Nb, hmpprt::s32 N3, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s64>  vol_id, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s64>  edge_id, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  boundary_values_12, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  edge_values_11)
;
#endif // __CUDACC__



# 8 "set_boundary.c"

#ifndef __CUDACC__
static hmpprt::CUDAGrid * set_boundary_values_from_edges_loop1D_1 = 0;
#else

extern "C" __global__ void set_boundary_values_from_edges_loop1D_1(hmpprt::s32 Nb, hmpprt::s64* vol_id, hmpprt::s64* edge_id, double* boundary_values_11, double* edge_values_12);
#endif // __CUDACC__




# 8 "set_boundary.c"

#ifdef __CUDACC__

extern "C" __global__ void set_boundary_values_from_edges_loop1D_1(hmpprt::s32 Nb, hmpprt::s64* vol_id, hmpprt::s64* edge_id, double* boundary_values_11, double* edge_values_12)
{
 # 16 "set_boundary.c"
 hmpprt::s32 id_1;
 # 16 "set_boundary.c"
 hmpprt::s32 k_1;
 # 38 "set_boundary.c"
 k_1 = (hmpprt::gr_atidf());
 # 38 "set_boundary.c"
 if (k_1 > Nb - 1)
 {
  # 38 "set_boundary.c"
  goto __hmppcg_label_1;
 }
 # 38 "set_boundary.c"
 id_1 = (hmpprt::s32 ) (3LL * *(vol_id + k_1) + *(edge_id + k_1));
 # 40 "set_boundary.c"
 *(boundary_values_11 + k_1) = *(edge_values_12 + id_1);
 # 8 "set_boundary.c"
 __hmppcg_label_1:;
}
#endif // __CUDACC__



# 8 "set_boundary.c"

#ifndef __CUDACC__
void set_boundary_values_from_edges_internal_1(hmpprt::s32 Nb, hmpprt::s32 N3, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s64>  vol_id, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s64>  edge_id, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  boundary_values_12, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  edge_values_11)
{
 # 8 "set_boundary.c"
 if (Nb - 1 >= 0)
 {
  hmpprt::CUDAGridCall __hmppcg_call;
  __hmppcg_call.setSizeX((Nb - 1) / 128 + 1);
  __hmppcg_call.setSizeY(1);
  __hmppcg_call.setBlockSizeX(32);
  __hmppcg_call.setBlockSizeY(4);
  __hmppcg_call.addLocalParameter((hmpprt::s32) (Nb), "Nb");
  __hmppcg_call.addLocalParameter(&vol_id, 8, "vol_id");
  __hmppcg_call.addLocalParameter(&edge_id, 8, "edge_id");
  __hmppcg_call.addLocalParameter(&boundary_values_12, 8, "boundary_values_11");
  __hmppcg_call.addLocalParameter(&edge_values_11, 8, "edge_values_12");
  __hmppcg_call.launch(set_boundary_values_from_edges_loop1D_1, hmpprt::Context::getInstance()->getCUDADevice());
 }
 ;
}
#endif // __CUDACC__



# 8 "set_boundary.c"

#ifndef __CUDACC__
extern "C" CDLT_API  void set_boundary_values_from_edges(hmpprt::s32 Nb, hmpprt::s32 N3, hmpprt::s64* vol_id, hmpprt::s64* edge_id, double* boundary_values_1, double* edge_values_1)
{
 # 1 "<preprocessor>"
 (set_boundary_values_from_edges_internal_1(Nb, N3, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s64> (vol_id), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s64> (edge_id), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (boundary_values_1), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (edge_values_1)));
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
      set_boundary_values_from_edges_loop1D_1 = new hmpprt::CUDAGrid(hmpprt_module, "set_boundary_values_from_edges_loop1D_1");

    }
    hmpprt::Context::getInstance()->getGrouplet()->setTarget(hmpprt::CUDA);
    hmpprt::Context::getInstance()->getGrouplet()->addSignature("set_boundary_values_from_edges", "prototype set_boundary_values_from_edges(Nb: s32, N3: s32, vol_id: ^cudaglob s64, edge_id: ^cudaglob s64, boundary_values: ^cudaglob double, edge_values: ^cudaglob double)");

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
      delete set_boundary_values_from_edges_loop1D_1;

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

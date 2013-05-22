
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

# 9 "update.c"

#ifndef __CUDACC__
extern "C" CDLT_API  void update(hmpprt::s32 N_1, double timestep_1, double* centroid_values_2, double* explicit_update_1, double* semi_implicit_update_1)
;
#endif // __CUDACC__



# 9 "update.c"

#ifndef __CUDACC__
void update_internal_1(hmpprt::s32 N_11, double timestep_11, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  centroid_values_21, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  explicit_update_11, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  semi_implicit_update_11)
;
#endif // __CUDACC__



# 9 "update.c"

#ifndef __CUDACC__
static hmpprt::CUDAGrid * update_loop1D_1 = 0;
#else

extern "C" __global__ void update_loop1D_1(hmpprt::s32 N_2, double timestep_2, double* centroid_values_1, double* explicit_update_2, double* semi_implicit_update_2);
#endif // __CUDACC__




# 9 "update.c"

#ifdef __CUDACC__

extern "C" __global__ void update_loop1D_1(hmpprt::s32 N_2, double timestep_2, double* centroid_values_1, double* explicit_update_2, double* semi_implicit_update_2)
{
 # 16 "update.c"
 double x_1;
 # 16 "update.c"
 double denominator_1;
 # 16 "update.c"
 hmpprt::s32 k_1;
 # 26 "update.c"
 k_1 = (hmpprt::gr_atidf());
 # 26 "update.c"
 if (k_1 > N_2 - 1)
 {
  # 26 "update.c"
  goto __hmppcg_label_1;
 }
 # 26 "update.c"
 x_1 = *(centroid_values_1 + k_1);
 # 27 "update.c"
 if (x_1 == (double) 0.0)
 {
  # 28 "update.c"
  *(semi_implicit_update_2 + k_1) = (double) 0.0;
 }
 else
 {
  # 30 "update.c"
  *(semi_implicit_update_2 + k_1) = *(semi_implicit_update_2 + k_1) / x_1;
 }
 # 37 "update.c"
 *(centroid_values_1 + k_1) = *(centroid_values_1 + k_1) + timestep_2 * *(explicit_update_2 + k_1);
 # 44 "update.c"
 denominator_1 = (double) 1.0 - timestep_2 * *(semi_implicit_update_2 + k_1);
 # 45 "update.c"
 if (denominator_1 <= (double) 0.0)
 {
  }
 else
 {
  # 50 "update.c"
  *(centroid_values_1 + k_1) = *(centroid_values_1 + k_1) / denominator_1;
 }
 # 9 "update.c"
 __hmppcg_label_1:;
}
#endif // __CUDACC__



# 9 "update.c"

#ifndef __CUDACC__
void update_internal_1(hmpprt::s32 N_11, double timestep_11, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  centroid_values_21, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  explicit_update_11, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  semi_implicit_update_11)
{
 # 9 "update.c"
 if (N_11 - 1 >= 0)
 {
  hmpprt::CUDAGridCall __hmppcg_call;
  __hmppcg_call.setSizeX((N_11 - 1) / 128 + 1);
  __hmppcg_call.setSizeY(1);
  __hmppcg_call.setBlockSizeX(32);
  __hmppcg_call.setBlockSizeY(4);
  __hmppcg_call.addLocalParameter((hmpprt::s32) (N_11), "N_2");
  __hmppcg_call.addLocalParameter(&timestep_11, 8, "timestep_2");
  __hmppcg_call.addLocalParameter(&centroid_values_21, 8, "centroid_values_1");
  __hmppcg_call.addLocalParameter(&explicit_update_11, 8, "explicit_update_2");
  __hmppcg_call.addLocalParameter(&semi_implicit_update_11, 8, "semi_implicit_update_2");
  __hmppcg_call.launch(update_loop1D_1, hmpprt::Context::getInstance()->getCUDADevice());
 }
 ;
}
#endif // __CUDACC__



# 9 "update.c"

#ifndef __CUDACC__
extern "C" CDLT_API  void update(hmpprt::s32 N_1, double timestep_1, double* centroid_values_2, double* explicit_update_1, double* semi_implicit_update_1)
{
 # 1 "<preprocessor>"
 (update_internal_1(N_1, timestep_1, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (centroid_values_2), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (explicit_update_1), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (semi_implicit_update_1)));
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
      update_loop1D_1 = new hmpprt::CUDAGrid(hmpprt_module, "update_loop1D_1");

    }
    hmpprt::Context::getInstance()->getGrouplet()->setTarget(hmpprt::CUDA);
    hmpprt::Context::getInstance()->getGrouplet()->addSignature("update", "prototype update(N: s32, timestep: double, centroid_values: ^cudaglob double, explicit_update: ^cudaglob double, semi_implicit_update: ^cudaglob double)");

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
      delete update_loop1D_1;

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

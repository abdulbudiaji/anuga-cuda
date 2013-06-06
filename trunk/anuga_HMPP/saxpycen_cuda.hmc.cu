
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

# 8 "saxpy_centroid_values.c"

#ifndef __CUDACC__
extern "C" CDLT_API  void saxpy_centroid_values(hmpprt::s32 N, double a, double b, double* centroid_values, double* centroid_backup_values)
;
#endif // __CUDACC__



# 8 "saxpy_centroid_values.c"

#ifndef __CUDACC__
void saxpy_centroid_values_internal_1(hmpprt::s32 N, double a, double b, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  centroid_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  centroid_backup_values)
;
#endif // __CUDACC__



# 8 "saxpy_centroid_values.c"

#ifndef __CUDACC__
static hmpprt::CUDAGrid * saxpy_centroid_values_loop1D_1 = 0;
#else

extern "C" __global__ void saxpy_centroid_values_loop1D_1(hmpprt::s32 N, double a, double b, double* centroid_values, double* centroid_backup_values);
#endif // __CUDACC__




# 8 "saxpy_centroid_values.c"

#ifdef __CUDACC__

extern "C" __global__ void saxpy_centroid_values_loop1D_1(hmpprt::s32 N, double a, double b, double* centroid_values, double* centroid_backup_values)
{
 # 16 "saxpy_centroid_values.c"
 hmpprt::s32 k_1;
 # 18 "saxpy_centroid_values.c"
 k_1 = (hmpprt::gr_atidf());
 # 18 "saxpy_centroid_values.c"
 if (k_1 > N - 1)
 {
  # 18 "saxpy_centroid_values.c"
  goto __hmppcg_label_1;
 }
 # 18 "saxpy_centroid_values.c"
 *(centroid_values + k_1) = a * *(centroid_values + k_1) + b * *(centroid_backup_values + k_1);
 # 8 "saxpy_centroid_values.c"
 __hmppcg_label_1:;
}
#endif // __CUDACC__



# 8 "saxpy_centroid_values.c"

#ifndef __CUDACC__
void saxpy_centroid_values_internal_1(hmpprt::s32 N, double a, double b, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  centroid_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  centroid_backup_values)
{
 # 8 "saxpy_centroid_values.c"
 if (N - 1 >= 0)
 {
  hmpprt::CUDAGridCall __hmppcg_call;
  __hmppcg_call.setSizeX((N - 1) / 128 + 1);
  __hmppcg_call.setSizeY(1);
  __hmppcg_call.setBlockSizeX(32);
  __hmppcg_call.setBlockSizeY(4);
  __hmppcg_call.addLocalParameter((hmpprt::s32) (N), "N");
  __hmppcg_call.addLocalParameter(&a, 8, "a");
  __hmppcg_call.addLocalParameter(&b, 8, "b");
  __hmppcg_call.addLocalParameter(&centroid_values, 8, "centroid_values");
  __hmppcg_call.addLocalParameter(&centroid_backup_values, 8, "centroid_backup_values");
  __hmppcg_call.launch(saxpy_centroid_values_loop1D_1, hmpprt::Context::getInstance()->getCUDADevice());
 }
 ;
}
#endif // __CUDACC__



# 8 "saxpy_centroid_values.c"

#ifndef __CUDACC__
extern "C" CDLT_API  void saxpy_centroid_values(hmpprt::s32 N, double a, double b, double* centroid_values, double* centroid_backup_values)
{
 # 1 "<preprocessor>"
 (saxpy_centroid_values_internal_1(N, a, b, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (centroid_values), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (centroid_backup_values)));
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
      saxpy_centroid_values_loop1D_1 = new hmpprt::CUDAGrid(hmpprt_module, "saxpy_centroid_values_loop1D_1");

    }
    hmpprt::Context::getInstance()->getGrouplet()->setTarget(hmpprt::CUDA);
    hmpprt::Context::getInstance()->getGrouplet()->addSignature("saxpy_centroid_values", "prototype saxpy_centroid_values(N: s32, a: double, b: double, centroid_values: ^cudaglob double, centroid_backup_values: ^cudaglob double)");

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
      delete saxpy_centroid_values_loop1D_1;

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

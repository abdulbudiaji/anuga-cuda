
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

# 7 "swb2_domain_ext.c"

#ifndef __CUDACC__
extern "C" CDLT_API  void protect_swb2(hmpprt::s64 N, hmpprt::s64 N3, double minimum_allowed_height, double maximum_allowed_speed, double epsilon, double* wc, double* wv, double* zc, double* zv, double* xmomc, double* ymomc, double* areas)
;
#endif // __CUDACC__



# 7 "swb2_domain_ext.c"

#ifndef __CUDACC__
void protect_swb2_internal_1(hmpprt::s64 N, hmpprt::s64 N3, double minimum_allowed_height, double maximum_allowed_speed, double epsilon, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  wc, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  wv, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  zc, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  zv, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  xmomc, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  ymomc, double* areas)
;
#endif // __CUDACC__



# 7 "swb2_domain_ext.c"

#ifndef __CUDACC__
static hmpprt::CUDAGrid * protect_swb2_loop1D_1 = 0;
#else

extern "C" __global__ void protect_swb2_loop1D_1(hmpprt::s64 N, double minimum_allowed_height, double* wc, double* wv, double* zc, double* zv, double* xmomc, double* ymomc);
#endif // __CUDACC__




# 7 "swb2_domain_ext.c"

#ifdef __CUDACC__

extern "C" __global__ void protect_swb2_loop1D_1(hmpprt::s64 N, double minimum_allowed_height, double* wc, double* wv, double* zc, double* zv, double* xmomc, double* ymomc)
{
 # 22 "swb2_domain_ext.c"
 double hc_1;
 # 22 "swb2_domain_ext.c"
 double bmax_1;
 # 22 "swb2_domain_ext.c"
 hmpprt::s32 k_1;
 # 34 "swb2_domain_ext.c"
 k_1 = (hmpprt::gr_atidf());
 # 34 "swb2_domain_ext.c"
 if (k_1 > (hmpprt::s32 ) (N) - 1)
 {
  # 34 "swb2_domain_ext.c"
  goto __hmppcg_label_1;
 }
 # 34 "swb2_domain_ext.c"
 hc_1 = *(wc + k_1) - *(zc + k_1);
 # 36 "swb2_domain_ext.c"
 bmax_1 = (double) 0.5 * fmax(*(zv + 3 * k_1) + *(zv + (3 * k_1 + 1)), fmax(*(zv + (3 * k_1 + 1)) + *(zv + (3 * k_1 + 2)), *(zv + (3 * k_1 + 2)) + *(zv + 3 * k_1)));
 # 38 "swb2_domain_ext.c"
 if (hc_1 < fmax((double) 0.1000000000000000056 * (bmax_1 - *(zc + k_1)), minimum_allowed_height))
 {
  # 43 "swb2_domain_ext.c"
  *(xmomc + k_1) = (double) 0.0;
  # 44 "swb2_domain_ext.c"
  *(ymomc + k_1) = (double) 0.0;
  # 47 "swb2_domain_ext.c"
  if (hc_1 <= (double) 0.0)
  {
   # 53 "swb2_domain_ext.c"
   double bmin_1;
   # 53 "swb2_domain_ext.c"
   bmin_1 = fmin(*(zv + 3 * k_1), fmin(*(zv + (3 * k_1 + 1)), *(zv + (3 * k_1 + 2)))) - minimum_allowed_height;
   # 57 "swb2_domain_ext.c"
   if (*(wc + k_1) < bmin_1)
   {
    # 62 "swb2_domain_ext.c"
    *(wc + k_1) = fmax(*(wc + k_1), bmin_1);
    # 69 "swb2_domain_ext.c"
    *(wv + 3 * k_1) = fmax(*(wv + 3 * k_1), bmin_1);
    # 70 "swb2_domain_ext.c"
    *(wv + (3 * k_1 + 1)) = fmax(*(wv + (3 * k_1 + 1)), bmin_1);
    # 71 "swb2_domain_ext.c"
    *(wv + (3 * k_1 + 2)) = fmax(*(wv + (3 * k_1 + 2)), bmin_1);
   }
  }
 }
 # 7 "swb2_domain_ext.c"
 __hmppcg_label_1:;
}
#endif // __CUDACC__



# 7 "swb2_domain_ext.c"

#ifndef __CUDACC__
void protect_swb2_internal_1(hmpprt::s64 N, hmpprt::s64 N3, double minimum_allowed_height, double maximum_allowed_speed, double epsilon, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  wc, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  wv, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  zc, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  zv, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  xmomc, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  ymomc, double* areas)
{
 # 7 "swb2_domain_ext.c"
 if ((hmpprt::s32 ) (N) - 1 >= 0)
 {
  hmpprt::CUDAGridCall __hmppcg_call;
  __hmppcg_call.setSizeX(((hmpprt::s32 ) (N) - 1) / 128 + 1);
  __hmppcg_call.setSizeY(1);
  __hmppcg_call.setBlockSizeX(32);
  __hmppcg_call.setBlockSizeY(4);
  __hmppcg_call.addLocalParameter(&N, 8, "N");
  __hmppcg_call.addLocalParameter(&minimum_allowed_height, 8, "minimum_allowed_height");
  __hmppcg_call.addLocalParameter(&wc, 8, "wc");
  __hmppcg_call.addLocalParameter(&wv, 8, "wv");
  __hmppcg_call.addLocalParameter(&zc, 8, "zc");
  __hmppcg_call.addLocalParameter(&zv, 8, "zv");
  __hmppcg_call.addLocalParameter(&xmomc, 8, "xmomc");
  __hmppcg_call.addLocalParameter(&ymomc, 8, "ymomc");
  __hmppcg_call.launch(protect_swb2_loop1D_1, hmpprt::Context::getInstance()->getCUDADevice());
 }
 ;
}
#endif // __CUDACC__



# 7 "swb2_domain_ext.c"

#ifndef __CUDACC__
extern "C" CDLT_API  void protect_swb2(hmpprt::s64 N, hmpprt::s64 N3, double minimum_allowed_height, double maximum_allowed_speed, double epsilon, double* wc, double* wv, double* zc, double* zv, double* xmomc, double* ymomc, double* areas)
{
 # 1 "<preprocessor>"
 (protect_swb2_internal_1(N, N3, minimum_allowed_height, maximum_allowed_speed, epsilon, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (wc), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (wv), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (zc), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (zv), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (xmomc), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (ymomc), areas));
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
      protect_swb2_loop1D_1 = new hmpprt::CUDAGrid(hmpprt_module, "protect_swb2_loop1D_1");

    }
    hmpprt::Context::getInstance()->getGrouplet()->setTarget(hmpprt::CUDA);
    hmpprt::Context::getInstance()->getGrouplet()->addSignature("protect_swb2", "prototype protect_swb2(N: s64, N3: s64, minimum_allowed_height: double, maximum_allowed_speed: double, epsilon: double, wc: ^cudaglob double, wv: ^cudaglob double, zc: ^cudaglob double, zv: ^cudaglob double, xmomc: ^cudaglob double, ymomc: ^cudaglob double, areas: unused ^host double)");

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
      delete protect_swb2_loop1D_1;

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

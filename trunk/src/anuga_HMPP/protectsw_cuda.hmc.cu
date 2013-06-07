
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

# 13 "protect.c"

#ifndef __CUDACC__
extern "C" CDLT_API  void protect_sw(hmpprt::s32 N_3, hmpprt::s32 N3_1, double minimum_allowed_height_1, double maximum_allowed_speed_1, double epsilon_1, double* wc_11, double* zc_3, double* xmomc_3, double* ymomc_3)
;
#endif // __CUDACC__



# 13 "protect.c"

#ifndef __CUDACC__
void protect_sw_internal_1(hmpprt::s32 N_31, hmpprt::s32 N3_11, double minimum_allowed_height_11, double maximum_allowed_speed_11, double epsilon_11, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  wc_1, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  zc_31, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  xmomc_31, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  ymomc_31)
;
#endif // __CUDACC__



# 13 "protect.c"

#ifndef __CUDACC__
static hmpprt::CUDAGrid * protect_sw_loop1D_2 = 0;
#else

extern "C" __global__ void protect_sw_loop1D_2(hmpprt::s32 N_2, double minimum_allowed_height_3, double maximum_allowed_speed_2, double* wc_2, double* zc_1, double* xmomc_1, double* ymomc_2);
#endif // __CUDACC__




# 13 "protect.c"

#ifndef __CUDACC__
static hmpprt::CUDAGrid * protect_sw_loop1D_1 = 0;
#else

extern "C" __global__ void protect_sw_loop1D_1(hmpprt::s32 N_1, double minimum_allowed_height_2, double* wc_3, double* zc_2, double* xmomc_2, double* ymomc_1);
#endif // __CUDACC__




# 13 "protect.c"

#ifdef __CUDACC__

extern "C" __global__ void protect_sw_loop1D_1(hmpprt::s32 N_1, double minimum_allowed_height_2, double* wc_3, double* zc_2, double* xmomc_2, double* ymomc_1)
{
 # 26 "protect.c"
 double hc_1;
 # 26 "protect.c"
 hmpprt::s32 k_1;
 # 38 "protect.c"
 k_1 = (hmpprt::gr_atidf());
 # 38 "protect.c"
 if (k_1 > N_1 - 1)
 {
  # 38 "protect.c"
  goto __hmppcg_label_1;
 }
 # 38 "protect.c"
 hc_1 = *(wc_3 + k_1) - *(zc_2 + k_1);
 # 39 "protect.c"
 if (hc_1 < minimum_allowed_height_2)
 {
  # 41 "protect.c"
  *(xmomc_2 + k_1) = (double) 0.0;
  # 42 "protect.c"
  *(ymomc_1 + k_1) = (double) 0.0;
  # 43 "protect.c"
  if (hc_1 <= (double) 0.0)
  {
   # 43 "protect.c"
   *(wc_3 + k_1) = *(zc_2 + k_1);
  }
 }
 # 13 "protect.c"
 __hmppcg_label_1:;
}
#endif // __CUDACC__



# 13 "protect.c"

#ifdef __CUDACC__

extern "C" __global__ void protect_sw_loop1D_2(hmpprt::s32 N_2, double minimum_allowed_height_3, double maximum_allowed_speed_2, double* wc_2, double* zc_1, double* xmomc_1, double* ymomc_2)
{
 # 26 "protect.c"
 double hc_2;
 # 26 "protect.c"
 hmpprt::s32 k_2;
 # 54 "protect.c"
 k_2 = (hmpprt::gr_atidf());
 # 54 "protect.c"
 if (k_2 > N_2 - 1)
 {
  # 54 "protect.c"
  goto __hmppcg_label_2;
 }
 # 54 "protect.c"
 hc_2 = *(wc_2 + k_2) - *(zc_1 + k_2);
 # 55 "protect.c"
 if (hc_2 < minimum_allowed_height_3)
 {
  # 58 "protect.c"
  if (hc_2 <= (double) 0.0)
  {
   # 59 "protect.c"
   *(wc_2 + k_2) = *(zc_1 + k_2);
   # 60 "protect.c"
   *(xmomc_1 + k_2) = (double) 0.0;
   # 61 "protect.c"
   *(ymomc_2 + k_2) = (double) 0.0;
  }
  else
  {
   # 67 "protect.c"
   double u_1;
   # 67 "protect.c"
   double v_1;
   # 67 "protect.c"
   u_1 = *(xmomc_1 + k_2) / hc_2;
   # 68 "protect.c"
   if (fabs(u_1) > maximum_allowed_speed_2)
   {
    # 69 "protect.c"
    double reduced_speed_1;
    # 69 "protect.c"
    reduced_speed_1 = maximum_allowed_speed_2 * u_1 / fabs(u_1);
    # 70 "protect.c"
    *(xmomc_1 + k_2) = reduced_speed_1 * hc_2;
   }
   # 73 "protect.c"
   v_1 = *(ymomc_2 + k_2) / hc_2;
   # 74 "protect.c"
   if (fabs(v_1) > maximum_allowed_speed_2)
   {
    # 75 "protect.c"
    double reduced_speed_2;
    # 75 "protect.c"
    reduced_speed_2 = maximum_allowed_speed_2 * v_1 / fabs(v_1);
    # 76 "protect.c"
    *(ymomc_2 + k_2) = reduced_speed_2 * hc_2;
   }
  }
 }
 # 13 "protect.c"
 __hmppcg_label_2:;
}
#endif // __CUDACC__



# 13 "protect.c"

#ifndef __CUDACC__
void protect_sw_internal_1(hmpprt::s32 N_31, hmpprt::s32 N3_11, double minimum_allowed_height_11, double maximum_allowed_speed_11, double epsilon_11, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  wc_1, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  zc_31, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  xmomc_31, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  ymomc_31)
{
 # 32 "protect.c"
 if (maximum_allowed_speed_11 < epsilon_11)
 {
  # 47 "protect.c"
  if (N_31 - 1 >= 0)
  {
   hmpprt::CUDAGridCall __hmppcg_call;
   __hmppcg_call.setSizeX((N_31 - 1) / 128 + 1);
   __hmppcg_call.setSizeY(1);
   __hmppcg_call.setBlockSizeX(32);
   __hmppcg_call.setBlockSizeY(4);
   __hmppcg_call.addLocalParameter((hmpprt::s32) (N_31), "N_1");
   __hmppcg_call.addLocalParameter(&minimum_allowed_height_11, 8, "minimum_allowed_height_2");
   __hmppcg_call.addLocalParameter(&wc_1, 8, "wc_3");
   __hmppcg_call.addLocalParameter(&zc_31, 8, "zc_2");
   __hmppcg_call.addLocalParameter(&xmomc_31, 8, "xmomc_2");
   __hmppcg_call.addLocalParameter(&ymomc_31, 8, "ymomc_1");
   __hmppcg_call.launch(protect_sw_loop1D_1, hmpprt::Context::getInstance()->getCUDADevice());
  }
  ;
 }
 else
 {
  # 13 "protect.c"
  if (N_31 - 1 >= 0)
  {
   hmpprt::CUDAGridCall __hmppcg_call;
   __hmppcg_call.setSizeX((N_31 - 1) / 128 + 1);
   __hmppcg_call.setSizeY(1);
   __hmppcg_call.setBlockSizeX(32);
   __hmppcg_call.setBlockSizeY(4);
   __hmppcg_call.addLocalParameter((hmpprt::s32) (N_31), "N_2");
   __hmppcg_call.addLocalParameter(&minimum_allowed_height_11, 8, "minimum_allowed_height_3");
   __hmppcg_call.addLocalParameter(&maximum_allowed_speed_11, 8, "maximum_allowed_speed_2");
   __hmppcg_call.addLocalParameter(&wc_1, 8, "wc_2");
   __hmppcg_call.addLocalParameter(&zc_31, 8, "zc_1");
   __hmppcg_call.addLocalParameter(&xmomc_31, 8, "xmomc_1");
   __hmppcg_call.addLocalParameter(&ymomc_31, 8, "ymomc_2");
   __hmppcg_call.launch(protect_sw_loop1D_2, hmpprt::Context::getInstance()->getCUDADevice());
  }
  ;
 }
}
#endif // __CUDACC__



# 13 "protect.c"

#ifndef __CUDACC__
extern "C" CDLT_API  void protect_sw(hmpprt::s32 N_3, hmpprt::s32 N3_1, double minimum_allowed_height_1, double maximum_allowed_speed_1, double epsilon_1, double* wc_11, double* zc_3, double* xmomc_3, double* ymomc_3)
{
 # 1 "<preprocessor>"
 (protect_sw_internal_1(N_3, N3_1, minimum_allowed_height_1, maximum_allowed_speed_1, epsilon_1, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (wc_11), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (zc_3), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (xmomc_3), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (ymomc_3)));
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
      protect_sw_loop1D_1 = new hmpprt::CUDAGrid(hmpprt_module, "protect_sw_loop1D_1");
      protect_sw_loop1D_2 = new hmpprt::CUDAGrid(hmpprt_module, "protect_sw_loop1D_2");

    }
    hmpprt::Context::getInstance()->getGrouplet()->setTarget(hmpprt::CUDA);
    hmpprt::Context::getInstance()->getGrouplet()->addSignature("protect_sw", "prototype protect_sw(N: s32, N3: s32, minimum_allowed_height: double, maximum_allowed_speed: double, epsilon: double, wc: ^cudaglob double, zc: ^cudaglob double, xmomc: ^cudaglob double, ymomc: ^cudaglob double)");

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
      delete protect_sw_loop1D_1;
      delete protect_sw_loop1D_2;

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

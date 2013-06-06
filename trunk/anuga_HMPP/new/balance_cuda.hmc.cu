
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

# 7 "balance_deep_and_shallow.c"

#ifndef __CUDACC__
extern "C" CDLT_API  void balance_deep_and_shallow(hmpprt::s32 N_21, hmpprt::s32 N3, double H0_21, double alpha_balance_2, hmpprt::s32 tight_slope_limiters_21, hmpprt::s32 use_centroid_velocities, double* wc_2, double* zc, double* wv, double* zv, double* xmomc_2, double* ymomc, double* xmomv, double* ymomv)
;
#endif // __CUDACC__



# 7 "balance_deep_and_shallow.c"

#ifndef __CUDACC__
void balance_deep_and_shallow_internal_1(hmpprt::s32 N_2, hmpprt::s32 N3, double H0_2, double alpha_balance_21, hmpprt::s32 tight_slope_limiters_2, hmpprt::s32 use_centroid_velocities, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  wc_21, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  zc, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  wv, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  zv, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  xmomc_21, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  ymomc, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  xmomv, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  ymomv)
;
#endif // __CUDACC__



# 7 "balance_deep_and_shallow.c"

#ifndef __CUDACC__
static hmpprt::CUDAGrid * balance_deep_and_shallow_loop1D_1 = 0;
#else

extern "C" __global__ void balance_deep_and_shallow_loop1D_1(hmpprt::s32 N, double H0, double alpha_balance, hmpprt::s32 tight_slope_limiters, hmpprt::s32 use_centroid_velocities_2, double* wc, double* zc_2, double* wv_2, double* zv_2, double* xmomc, double* ymomc_2, double* xmomv_2, double* ymomv_2);
#endif // __CUDACC__




# 7 "balance_deep_and_shallow.c"

#ifdef __CUDACC__

extern "C" __global__ void balance_deep_and_shallow_loop1D_1(hmpprt::s32 N, double H0, double alpha_balance, hmpprt::s32 tight_slope_limiters, hmpprt::s32 use_centroid_velocities_2, double* wc, double* zc_2, double* wv_2, double* zv_2, double* xmomc, double* ymomc_2, double* xmomv_2, double* ymomv_2)
{
 # 30 "balance_deep_and_shallow.c"
 double hv_1[3uLL];
 # 26 "balance_deep_and_shallow.c"
 double dz_1;
 # 26 "balance_deep_and_shallow.c"
 hmpprt::s32 k3_1;
 # 26 "balance_deep_and_shallow.c"
 double hmin_1;
 # 26 "balance_deep_and_shallow.c"
 double hc_k_1;
 # 26 "balance_deep_and_shallow.c"
 double alpha_1;
 # 26 "balance_deep_and_shallow.c"
 hmpprt::s32 k_1;
 # 50 "balance_deep_and_shallow.c"
 k_1 = (hmpprt::gr_atidf());
 # 50 "balance_deep_and_shallow.c"
 if (k_1 > N - 1)
 {
  # 50 "balance_deep_and_shallow.c"
  goto __hmppcg_label_1;
 }
 # 50 "balance_deep_and_shallow.c"
 k3_1 = 3 * k_1;
 # 51 "balance_deep_and_shallow.c"
 hc_k_1 = *(wc + k_1) - *(zc_2 + k_1);
 # 53 "balance_deep_and_shallow.c"
 dz_1 = (double) 0.0;
 # 54 "balance_deep_and_shallow.c"
 if (tight_slope_limiters == 0)
 {
  # 56 "balance_deep_and_shallow.c"
  hmpprt::s32 i_1;
  # 56 "balance_deep_and_shallow.c"
  # 56 "balance_deep_and_shallow.c"
  for (i_1 = 0 ; i_1 <= 2 ; i_1 = i_1 + 1)
  {
   # 57 "balance_deep_and_shallow.c"
   dz_1 = fmax(dz_1, fabs(*(zv_2 + (k3_1 + i_1)) - *(zc_2 + k_1)));
  }
  # 62 "balance_deep_and_shallow.c"
 }
 # 62 "balance_deep_and_shallow.c"
 hv_1[0] = *(wv_2 + k3_1) - *(zv_2 + k3_1);
 # 63 "balance_deep_and_shallow.c"
 hv_1[1] = *(wv_2 + (k3_1 + 1)) - *(zv_2 + (k3_1 + 1));
 # 64 "balance_deep_and_shallow.c"
 hv_1[2] = *(wv_2 + (k3_1 + 2)) - *(zv_2 + (k3_1 + 2));
 # 67 "balance_deep_and_shallow.c"
 hmin_1 = fmin(hv_1[0], fmin(hv_1[1], hv_1[2]));
 # 77 "balance_deep_and_shallow.c"
 if (tight_slope_limiters == 0)
 {
  # 84 "balance_deep_and_shallow.c"
  if (dz_1 > (double) 0.0)
  {
   # 85 "balance_deep_and_shallow.c"
   alpha_1 = fmax(fmin(alpha_balance * hmin_1 / dz_1, (double) 1.0), (double) 0.0);
  }
  else
  {
   # 87 "balance_deep_and_shallow.c"
   alpha_1 = (double) 1.0;
  }
 }
 else
 {
  # 99 "balance_deep_and_shallow.c"
  if (hmin_1 < H0)
  {
   # 100 "balance_deep_and_shallow.c"
   alpha_1 = (double) 1.0;
   # 101 "balance_deep_and_shallow.c"
   hmpprt::s32 i_2;
   # 101 "balance_deep_and_shallow.c"
   # 101 "balance_deep_and_shallow.c"
   for (i_2 = 0 ; i_2 <= 2 ; i_2 = i_2 + 1)
   {
    # 103 "balance_deep_and_shallow.c"
    double h_diff_1;
    # 103 "balance_deep_and_shallow.c"
    h_diff_1 = hc_k_1 - hv_1[i_2];
    # 104 "balance_deep_and_shallow.c"
    if (h_diff_1 <= (double) 0.0)
    {
     }
    else
    {
     # 112 "balance_deep_and_shallow.c"
     alpha_1 = fmin(alpha_1, (hc_k_1 - H0) / h_diff_1);
    }
   }
   # 117 "balance_deep_and_shallow.c"
   # 117 "balance_deep_and_shallow.c"
   if (alpha_1 > (double) 1.0)
   {
    # 117 "balance_deep_and_shallow.c"
    alpha_1 = (double) 1.0;
   }
   # 118 "balance_deep_and_shallow.c"
   if (alpha_1 < (double) 0.0)
   {
    # 118 "balance_deep_and_shallow.c"
    alpha_1 = (double) 0.0;
   }
  }
  else
  {
   # 122 "balance_deep_and_shallow.c"
   alpha_1 = (double) 1.0;
  }
 }
 # 145 "balance_deep_and_shallow.c"
 if (alpha_1 < (double) 1.0)
 {
  # 146 "balance_deep_and_shallow.c"
  hmpprt::s32 i_3;
  # 146 "balance_deep_and_shallow.c"
  # 146 "balance_deep_and_shallow.c"
  for (i_3 = 0 ; i_3 <= 2 ; i_3 = i_3 + 1)
  {
   # 148 "balance_deep_and_shallow.c"
   *(wv_2 + (k3_1 + i_3)) = *(zv_2 + (k3_1 + i_3)) + ((double) 1.0 - alpha_1) * hc_k_1 + alpha_1 * hv_1[i_3];
   # 151 "balance_deep_and_shallow.c"
   if (use_centroid_velocities_2 == 1)
   {
    # 157 "balance_deep_and_shallow.c"
    double uc_1;
    # 157 "balance_deep_and_shallow.c"
    double vc_1;
    # 157 "balance_deep_and_shallow.c"
    if (hc_k_1 > (double) 9.999999999999999547e-07)
    {
     # 158 "balance_deep_and_shallow.c"
     uc_1 = *(xmomc + k_1) / hc_k_1;
     # 159 "balance_deep_and_shallow.c"
     vc_1 = *(ymomc_2 + k_1) / hc_k_1;
    }
    else
    {
     # 161 "balance_deep_and_shallow.c"
     uc_1 = (double) 0.0;
     # 162 "balance_deep_and_shallow.c"
     vc_1 = (double) 0.0;
    }
    # 167 "balance_deep_and_shallow.c"
    hv_1[i_3] = *(wv_2 + (k3_1 + i_3)) - *(zv_2 + (k3_1 + i_3));
    # 168 "balance_deep_and_shallow.c"
    *(xmomv_2 + (k3_1 + i_3)) = uc_1 * hv_1[i_3];
    # 169 "balance_deep_and_shallow.c"
    *(ymomv_2 + (k3_1 + i_3)) = vc_1 * hv_1[i_3];
   }
   else
   {
    # 183 "balance_deep_and_shallow.c"
    *(xmomv_2 + (k3_1 + i_3)) = ((double) 1.0 - alpha_1) * *(xmomc + k_1) + alpha_1 * *(xmomv_2 + (k3_1 + i_3));
    # 184 "balance_deep_and_shallow.c"
    *(ymomv_2 + (k3_1 + i_3)) = ((double) 1.0 - alpha_1) * *(ymomc_2 + k_1) + alpha_1 * *(ymomv_2 + (k3_1 + i_3));
   }
  }
  # 7 "balance_deep_and_shallow.c"
 }
 # 7 "balance_deep_and_shallow.c"
 __hmppcg_label_1:;
}
#endif // __CUDACC__



# 7 "balance_deep_and_shallow.c"

#ifndef __CUDACC__
void balance_deep_and_shallow_internal_1(hmpprt::s32 N_2, hmpprt::s32 N3, double H0_2, double alpha_balance_21, hmpprt::s32 tight_slope_limiters_2, hmpprt::s32 use_centroid_velocities, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  wc_21, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  zc, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  wv, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  zv, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  xmomc_21, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  ymomc, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  xmomv, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  ymomv)
{
 # 7 "balance_deep_and_shallow.c"
 if (N_2 - 1 >= 0)
 {
  hmpprt::CUDAGridCall __hmppcg_call;
  __hmppcg_call.setSizeX((N_2 - 1) / 128 + 1);
  __hmppcg_call.setSizeY(1);
  __hmppcg_call.setBlockSizeX(32);
  __hmppcg_call.setBlockSizeY(4);
  __hmppcg_call.addLocalParameter((hmpprt::s32) (N_2), "N");
  __hmppcg_call.addLocalParameter(&H0_2, 8, "H0");
  __hmppcg_call.addLocalParameter(&alpha_balance_21, 8, "alpha_balance");
  __hmppcg_call.addLocalParameter((hmpprt::s32) (tight_slope_limiters_2), "tight_slope_limiters");
  __hmppcg_call.addLocalParameter((hmpprt::s32) (use_centroid_velocities), "use_centroid_velocities_2");
  __hmppcg_call.addLocalParameter(&wc_21, 8, "wc");
  __hmppcg_call.addLocalParameter(&zc, 8, "zc_2");
  __hmppcg_call.addLocalParameter(&wv, 8, "wv_2");
  __hmppcg_call.addLocalParameter(&zv, 8, "zv_2");
  __hmppcg_call.addLocalParameter(&xmomc_21, 8, "xmomc");
  __hmppcg_call.addLocalParameter(&ymomc, 8, "ymomc_2");
  __hmppcg_call.addLocalParameter(&xmomv, 8, "xmomv_2");
  __hmppcg_call.addLocalParameter(&ymomv, 8, "ymomv_2");
  __hmppcg_call.launch(balance_deep_and_shallow_loop1D_1, hmpprt::Context::getInstance()->getCUDADevice());
 }
 ;
}
#endif // __CUDACC__



# 7 "balance_deep_and_shallow.c"

#ifndef __CUDACC__
extern "C" CDLT_API  void balance_deep_and_shallow(hmpprt::s32 N_21, hmpprt::s32 N3, double H0_21, double alpha_balance_2, hmpprt::s32 tight_slope_limiters_21, hmpprt::s32 use_centroid_velocities, double* wc_2, double* zc, double* wv, double* zv, double* xmomc_2, double* ymomc, double* xmomv, double* ymomv)
{
 # 1 "<preprocessor>"
 (balance_deep_and_shallow_internal_1(N_21, N3, H0_21, alpha_balance_2, tight_slope_limiters_21, use_centroid_velocities, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (wc_2), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (zc), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (wv), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (zv), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (xmomc_2), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (ymomc), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (xmomv), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (ymomv)));
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
      balance_deep_and_shallow_loop1D_1 = new hmpprt::CUDAGrid(hmpprt_module, "balance_deep_and_shallow_loop1D_1");

    }
    hmpprt::Context::getInstance()->getGrouplet()->setTarget(hmpprt::CUDA);
    hmpprt::Context::getInstance()->getGrouplet()->addSignature("balance_deep_and_shallow", "prototype balance_deep_and_shallow(N: s32, N3: s32, H0: double, alpha_balance: double, tight_slope_limiters: s32, use_centroid_velocities: s32, wc: ^cudaglob double, zc: ^cudaglob double, wv: ^cudaglob double, zv: ^cudaglob double, xmomc: ^cudaglob double, ymomc: ^cudaglob double, xmomv: ^cudaglob double, ymomv: ^cudaglob double)");

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
      delete balance_deep_and_shallow_loop1D_1;

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

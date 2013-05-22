
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

# 67 "manning_friction.c"

#ifndef __CUDACC__
extern "C" CDLT_API  void manning_friction_sloped(hmpprt::s32 N_2, hmpprt::s32 N3, hmpprt::s32 N6, double g, double eps, double* x, double* w, double* zv, double* uh, double* vh, double* eta, double* xmom_update, double* ymom_update)
;
#endif // __CUDACC__



# 67 "manning_friction.c"

#ifndef __CUDACC__
void manning_friction_sloped_internal_1(hmpprt::s32 N_21, hmpprt::s32 N3, hmpprt::s32 N6, double g, double eps, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  x, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  w, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  zv, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  uh, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  vh, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  eta, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  xmom_update, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  ymom_update)
;
#endif // __CUDACC__



# 67 "manning_friction.c"

#ifndef __CUDACC__
static hmpprt::CUDAGrid * manning_friction_sloped_loop1D_1 = 0;
#else

extern "C" __global__ void manning_friction_sloped_loop1D_1(hmpprt::s32 N, double g_2, double eps_2, double* x_2, double* w_2, double* zv_2, double* uh_2, double* vh_2, double* eta_2, double* xmom_update_2, double* ymom_update_2);
#endif // __CUDACC__




# 67 "manning_friction.c"

#ifdef __CUDACC__

extern "C" __global__ void manning_friction_sloped_loop1D_1(hmpprt::s32 N, double g_2, double eps_2, double* x_2, double* w_2, double* zv_2, double* uh_2, double* vh_2, double* eta_2, double* xmom_update_2, double* ymom_update_2)
{
 # 83 "manning_friction.c"
 hmpprt::s32 k_1;
 # 97 "manning_friction.c"
 k_1 = (hmpprt::gr_atidf());
 # 97 "manning_friction.c"
 if (k_1 > N - 1)
 {
  # 97 "manning_friction.c"
  goto __hmppcg_label_1;
 }
 # 97 "manning_friction.c"
 if (*(eta_2 + k_1) > eps_2)
 {
  # 99 "manning_friction.c"
  hmpprt::s32 k3_1;
  # 99 "manning_friction.c"
  hmpprt::s32 k6_1;
  # 99 "manning_friction.c"
  double y2_1;
  # 99 "manning_friction.c"
  double y0_1;
  # 99 "manning_friction.c"
  double x1_1;
  # 99 "manning_friction.c"
  double x0_1;
  # 99 "manning_friction.c"
  double y1_1;
  # 99 "manning_friction.c"
  double x2_1;
  # 99 "manning_friction.c"
  double z1_1;
  # 99 "manning_friction.c"
  double z0_1;
  # 99 "manning_friction.c"
  double z2_1;
  # 99 "manning_friction.c"
  double zx_1;
  # 99 "manning_friction.c"
  double det_1;
  # 99 "manning_friction.c"
  double zy_1;
  # 99 "manning_friction.c"
  double zx_2;
  # 99 "manning_friction.c"
  double zy_2;
  # 99 "manning_friction.c"
  double z_1;
  # 99 "manning_friction.c"
  double h_1;
  # 99 "manning_friction.c"
  double zs_1;
  # 99 "manning_friction.c"
  k3_1 = 3 * k_1;
  # 100 "manning_friction.c"
  z0_1 = *(zv_2 + k3_1);
  # 101 "manning_friction.c"
  z1_1 = *(zv_2 + (k3_1 + 1));
  # 102 "manning_friction.c"
  z2_1 = *(zv_2 + (k3_1 + 2));
  # 105 "manning_friction.c"
  k6_1 = 6 * k_1;
  # 107 "manning_friction.c"
  x0_1 = *(x_2 + k6_1);
  # 108 "manning_friction.c"
  y0_1 = *(x_2 + (k6_1 + 1));
  # 109 "manning_friction.c"
  x1_1 = *(x_2 + (k6_1 + 2));
  # 110 "manning_friction.c"
  y1_1 = *(x_2 + (k6_1 + 3));
  # 111 "manning_friction.c"
  x2_1 = *(x_2 + (k6_1 + 4));
  # 112 "manning_friction.c"
  y2_1 = *(x_2 + (k6_1 + 5));
  # 128 "manning_friction.c"
  det_1 = (y2_1 - y0_1) * (x1_1 - x0_1) - (y1_1 - y0_1) * (x2_1 - x0_1);
  # 130 "manning_friction.c"
  zx_1 = (y2_1 - y0_1) * (z1_1 - z0_1) - (y1_1 - y0_1) * (z2_1 - z0_1);
  # 131 "manning_friction.c"
  zx_2 = zx_1 / det_1;
  # 133 "manning_friction.c"
  zy_1 = (x1_1 - x0_1) * (z2_1 - z0_1) - (x2_1 - x0_1) * (z1_1 - z0_1);
  # 134 "manning_friction.c"
  zy_2 = zy_1 / det_1;
  # 137 "manning_friction.c"
  zs_1 = sqrt((double) 1.0 + zx_2 * zx_2 + zy_2 * zy_2);
  # 138 "manning_friction.c"
  z_1 = (z0_1 + z1_1 + z2_1) / (double) 3.0;
  # 139 "manning_friction.c"
  h_1 = *(w_2 + k_1) - z_1;
  # 140 "manning_friction.c"
  if (h_1 >= eps_2)
  {
   # 141 "manning_friction.c"
   double S_1;
   # 141 "manning_friction.c"
   double S_2;
   # 141 "manning_friction.c"
   S_1 =  - g_2 * *(eta_2 + k_1) * *(eta_2 + k_1) * zs_1 * sqrt(*(uh_2 + k_1) * *(uh_2 + k_1) + *(vh_2 + k_1) * *(vh_2 + k_1));
   # 142 "manning_friction.c"
   S_2 = S_1 / pow(h_1, (double) 2.333333333333333481);
   # 149 "manning_friction.c"
   *(xmom_update_2 + k_1) = *(xmom_update_2 + k_1) + S_2 * *(uh_2 + k_1);
   # 150 "manning_friction.c"
   *(ymom_update_2 + k_1) = *(ymom_update_2 + k_1) + S_2 * *(vh_2 + k_1);
  }
 }
 # 67 "manning_friction.c"
 __hmppcg_label_1:;
}
#endif // __CUDACC__



# 67 "manning_friction.c"

#ifndef __CUDACC__
void manning_friction_sloped_internal_1(hmpprt::s32 N_21, hmpprt::s32 N3, hmpprt::s32 N6, double g, double eps, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  x, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  w, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  zv, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  uh, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  vh, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  eta, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  xmom_update, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  ymom_update)
{
 # 67 "manning_friction.c"
 if (N_21 - 1 >= 0)
 {
  hmpprt::CUDAGridCall __hmppcg_call;
  __hmppcg_call.setSizeX((N_21 - 1) / 128 + 1);
  __hmppcg_call.setSizeY(1);
  __hmppcg_call.setBlockSizeX(32);
  __hmppcg_call.setBlockSizeY(4);
  __hmppcg_call.addLocalParameter((hmpprt::s32) (N_21), "N");
  __hmppcg_call.addLocalParameter(&g, 8, "g_2");
  __hmppcg_call.addLocalParameter(&eps, 8, "eps_2");
  __hmppcg_call.addLocalParameter(&x, 8, "x_2");
  __hmppcg_call.addLocalParameter(&w, 8, "w_2");
  __hmppcg_call.addLocalParameter(&zv, 8, "zv_2");
  __hmppcg_call.addLocalParameter(&uh, 8, "uh_2");
  __hmppcg_call.addLocalParameter(&vh, 8, "vh_2");
  __hmppcg_call.addLocalParameter(&eta, 8, "eta_2");
  __hmppcg_call.addLocalParameter(&xmom_update, 8, "xmom_update_2");
  __hmppcg_call.addLocalParameter(&ymom_update, 8, "ymom_update_2");
  __hmppcg_call.launch(manning_friction_sloped_loop1D_1, hmpprt::Context::getInstance()->getCUDADevice());
 }
 ;
}
#endif // __CUDACC__



# 67 "manning_friction.c"

#ifndef __CUDACC__
extern "C" CDLT_API  void manning_friction_sloped(hmpprt::s32 N_2, hmpprt::s32 N3, hmpprt::s32 N6, double g, double eps, double* x, double* w, double* zv, double* uh, double* vh, double* eta, double* xmom_update, double* ymom_update)
{
 # 1 "<preprocessor>"
 (manning_friction_sloped_internal_1(N_2, N3, N6, g, eps, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (x), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (w), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (zv), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (uh), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (vh), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (eta), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (xmom_update), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (ymom_update)));
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
      manning_friction_sloped_loop1D_1 = new hmpprt::CUDAGrid(hmpprt_module, "manning_friction_sloped_loop1D_1");

    }
    hmpprt::Context::getInstance()->getGrouplet()->setTarget(hmpprt::CUDA);
    hmpprt::Context::getInstance()->getGrouplet()->addSignature("manning_friction_sloped", "prototype manning_friction_sloped(N: s32, N3: s32, N6: s32, g: double, eps: double, x: ^cudaglob double, w: ^cudaglob double, zv: ^cudaglob double, uh: ^cudaglob double, vh: ^cudaglob double, eta: ^cudaglob double, xmom_update: ^cudaglob double, ymom_update: ^cudaglob double)");

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
      delete manning_friction_sloped_loop1D_1;

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

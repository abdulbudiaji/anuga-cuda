
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

# 8 "manning_friction.c"

#ifndef __CUDACC__
extern "C" CDLT_API  void manning_friction_flat(hmpprt::s32 N_11, hmpprt::s32 N3_11, double g_1, double eps_1, double* w_1, double* zv_1, double* uh_1, double* vh_2, double* eta_11, double* xmom_1, double* ymom_2)
;
#endif // __CUDACC__



# 8 "manning_friction.c"

#ifndef __CUDACC__
void manning_friction_flat_internal_1(hmpprt::s32 N_1, hmpprt::s32 N3_1, double g_11, double eps_11, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  w_11, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  zv_11, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  uh_11, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  vh_21, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  eta_1, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  xmom_11, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  ymom_21)
;
#endif // __CUDACC__



# 8 "manning_friction.c"

#ifndef __CUDACC__
static hmpprt::CUDAGrid * manning_friction_flat_loop1D_1 = 0;
#else

extern "C" __global__ void manning_friction_flat_loop1D_1(hmpprt::s32 N_2, double g_2, double eps_2, double* w_2, double* zv_2, double* uh_2, double* vh_1, double* eta_2, double* xmom_2, double* ymom_1);
#endif // __CUDACC__




# 8 "manning_friction.c"

#ifdef __CUDACC__

extern "C" __global__ void manning_friction_flat_loop1D_1(hmpprt::s32 N_2, double g_2, double eps_2, double* w_2, double* zv_2, double* uh_2, double* vh_1, double* eta_2, double* xmom_2, double* ymom_1)
{
 # 22 "manning_friction.c"
 hmpprt::s32 k_1;
 # 32 "manning_friction.c"
 k_1 = (hmpprt::gr_atidf());
 # 32 "manning_friction.c"
 if (k_1 > N_2 - 1)
 {
  # 32 "manning_friction.c"
  goto __hmppcg_label_1;
 }
 # 32 "manning_friction.c"
 if (*(eta_2 + k_1) > eps_2)
 {
  # 34 "manning_friction.c"
  hmpprt::s32 k3_1;
  # 34 "manning_friction.c"
  double z0_1;
  # 34 "manning_friction.c"
  double z1_1;
  # 34 "manning_friction.c"
  double z2_1;
  # 34 "manning_friction.c"
  double z_1;
  # 34 "manning_friction.c"
  double h_1;
  # 34 "manning_friction.c"
  k3_1 = 3 * k_1;
  # 35 "manning_friction.c"
  z0_1 = *(zv_2 + k3_1);
  # 36 "manning_friction.c"
  z1_1 = *(zv_2 + (k3_1 + 1));
  # 37 "manning_friction.c"
  z2_1 = *(zv_2 + (k3_1 + 2));
  # 44 "manning_friction.c"
  z_1 = (z0_1 + z1_1 + z2_1) / (double) 3.0;
  # 45 "manning_friction.c"
  h_1 = *(w_2 + k_1) - z_1;
  # 46 "manning_friction.c"
  if (h_1 >= eps_2)
  {
   # 47 "manning_friction.c"
   double S_1;
   # 47 "manning_friction.c"
   double S_2;
   # 47 "manning_friction.c"
   S_1 =  - g_2 * *(eta_2 + k_1) * *(eta_2 + k_1) * sqrt(*(uh_2 + k_1) * *(uh_2 + k_1) + *(vh_1 + k_1) * *(vh_1 + k_1));
   # 48 "manning_friction.c"
   S_2 = S_1 / pow(h_1, (double) 2.333333333333333481);
   # 55 "manning_friction.c"
   *(xmom_2 + k_1) = *(xmom_2 + k_1) + S_2 * *(uh_2 + k_1);
   # 56 "manning_friction.c"
   *(ymom_1 + k_1) = *(ymom_1 + k_1) + S_2 * *(vh_1 + k_1);
  }
 }
 # 8 "manning_friction.c"
 __hmppcg_label_1:;
}
#endif // __CUDACC__



# 8 "manning_friction.c"

#ifndef __CUDACC__
void manning_friction_flat_internal_1(hmpprt::s32 N_1, hmpprt::s32 N3_1, double g_11, double eps_11, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  w_11, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  zv_11, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  uh_11, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  vh_21, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  eta_1, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  xmom_11, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  ymom_21)
{
 # 8 "manning_friction.c"
 if (N_1 - 1 >= 0)
 {
  hmpprt::CUDAGridCall __hmppcg_call;
  __hmppcg_call.setSizeX((N_1 - 1) / 128 + 1);
  __hmppcg_call.setSizeY(1);
  __hmppcg_call.setBlockSizeX(32);
  __hmppcg_call.setBlockSizeY(4);
  __hmppcg_call.addLocalParameter((hmpprt::s32) (N_1), "N_2");
  __hmppcg_call.addLocalParameter(&g_11, 8, "g_2");
  __hmppcg_call.addLocalParameter(&eps_11, 8, "eps_2");
  __hmppcg_call.addLocalParameter(&w_11, 8, "w_2");
  __hmppcg_call.addLocalParameter(&zv_11, 8, "zv_2");
  __hmppcg_call.addLocalParameter(&uh_11, 8, "uh_2");
  __hmppcg_call.addLocalParameter(&vh_21, 8, "vh_1");
  __hmppcg_call.addLocalParameter(&eta_1, 8, "eta_2");
  __hmppcg_call.addLocalParameter(&xmom_11, 8, "xmom_2");
  __hmppcg_call.addLocalParameter(&ymom_21, 8, "ymom_1");
  __hmppcg_call.launch(manning_friction_flat_loop1D_1, hmpprt::Context::getInstance()->getCUDADevice());
 }
 ;
}
#endif // __CUDACC__



# 8 "manning_friction.c"

#ifndef __CUDACC__
extern "C" CDLT_API  void manning_friction_flat(hmpprt::s32 N_11, hmpprt::s32 N3_11, double g_1, double eps_1, double* w_1, double* zv_1, double* uh_1, double* vh_2, double* eta_11, double* xmom_1, double* ymom_2)
{
 # 1 "<preprocessor>"
 (manning_friction_flat_internal_1(N_11, N3_11, g_1, eps_1, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (w_1), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (zv_1), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (uh_1), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (vh_2), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (eta_11), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (xmom_1), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (ymom_2)));
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
      manning_friction_flat_loop1D_1 = new hmpprt::CUDAGrid(hmpprt_module, "manning_friction_flat_loop1D_1");

    }
    hmpprt::Context::getInstance()->getGrouplet()->setTarget(hmpprt::CUDA);
    hmpprt::Context::getInstance()->getGrouplet()->addSignature("manning_friction_flat", "prototype manning_friction_flat(N: s32, N3: s32, g: double, eps: double, w: ^cudaglob double, zv: ^cudaglob double, uh: ^cudaglob double, vh: ^cudaglob double, eta: ^cudaglob double, xmom: ^cudaglob double, ymom: ^cudaglob double)");

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
      delete manning_friction_flat_loop1D_1;

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

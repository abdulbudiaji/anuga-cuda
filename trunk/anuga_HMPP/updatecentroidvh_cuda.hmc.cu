
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

# 8 "update_centroids_of_velocities_and_height.c"

#ifndef __CUDACC__
extern "C" CDLT_API  void _update_centroids_of_velocities_and_height(hmpprt::s32 N_c_11, hmpprt::s32 N_b_1, double* w_C_1, double* uh_C_11, double* vh_C_1, double* h_C_12, double* z_C_12, double* u_C_1, double* v_C_11, double* w_B_11, double* uh_B_11, double* vh_B_11, double* h_B_11, double* z_B_11, double* u_B_12, double* v_B_12)
;
#endif // __CUDACC__



# 8 "update_centroids_of_velocities_and_height.c"

#ifndef __CUDACC__
void _update_centroids_of_velocities_and_height_internal_1(hmpprt::s32 N_c_12, hmpprt::s32 N_b_11, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  w_C_11, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  uh_C_1, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  vh_C_11, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  h_C_1, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  z_C_1, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  u_C_11, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  v_C_1, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  w_B_1, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  uh_B_1, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  vh_B_1, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  h_B_1, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  z_B_1, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  u_B_1, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  v_B_11)
;
#endif // __CUDACC__



# 8 "update_centroids_of_velocities_and_height.c"

#ifndef __CUDACC__
static hmpprt::CUDAGrid * _update_centroids_of_velocities_and_height_loop1D_1 = 0;
#else

extern "C" __global__ void _update_centroids_of_velocities_and_height_loop1D_1(hmpprt::s32 N_c_1, hmpprt::s32 N_b_12, double* w_C_12, double* uh_C_12, double* vh_C_12, double* h_C_11, double* z_C_11, double* u_C_12, double* v_C_12, double* w_B_12, double* uh_B_12, double* vh_B_12, double* h_B_12, double* z_B_12, double* u_B_11, double* v_B_1);
#endif // __CUDACC__




# 8 "update_centroids_of_velocities_and_height.c"

#ifdef __CUDACC__

extern "C" __global__ void _update_centroids_of_velocities_and_height_loop1D_1(hmpprt::s32 N_c_1, hmpprt::s32 N_b_12, double* w_C_12, double* uh_C_12, double* vh_C_12, double* h_C_11, double* z_C_11, double* u_C_12, double* v_C_12, double* w_B_12, double* uh_B_12, double* vh_B_12, double* h_B_12, double* z_B_12, double* u_B_11, double* v_B_1)
{
 # 28 "update_centroids_of_velocities_and_height.c"
 hmpprt::s32 k_1;
 # 35 "update_centroids_of_velocities_and_height.c"
 k_1 = (hmpprt::gr_atidf());
 # 35 "update_centroids_of_velocities_and_height.c"
 if (k_1 > N_c_1 - 1)
 {
  # 35 "update_centroids_of_velocities_and_height.c"
  goto __hmppcg_label_1;
 }
 # 35 "update_centroids_of_velocities_and_height.c"
 if (k_1 < N_c_1)
 {
  # 37 "update_centroids_of_velocities_and_height.c"
  double factor_1;
  # 37 "update_centroids_of_velocities_and_height.c"
  *(h_C_11 + k_1) = *(w_C_12 + k_1) - *(z_C_11 + k_1);
  # 38 "update_centroids_of_velocities_and_height.c"
  if (*(h_C_11 + k_1) < (double) 0.0)
  {
   # 39 "update_centroids_of_velocities_and_height.c"
   *(h_C_11 + k_1) = (double) 0.0;
  }
  # 42 "update_centroids_of_velocities_and_height.c"
  factor_1 = *(h_C_11 + k_1) / (*(h_C_11 + k_1) * *(h_C_11 + k_1) + (double) 1.000000000000000021e-08);
  # 43 "update_centroids_of_velocities_and_height.c"
  *(u_C_12 + k_1) = *(uh_C_12 + k_1) * factor_1;
  # 44 "update_centroids_of_velocities_and_height.c"
  *(v_C_12 + k_1) = *(vh_C_12 + k_1) * factor_1;
 }
 # 48 "update_centroids_of_velocities_and_height.c"
 if (k_1 < N_b_12)
 {
  # 50 "update_centroids_of_velocities_and_height.c"
  double factor_2;
  # 50 "update_centroids_of_velocities_and_height.c"
  *(h_B_12 + k_1) = *(w_B_12 + k_1) - *(z_B_12 + k_1);
  # 51 "update_centroids_of_velocities_and_height.c"
  if (*(h_B_12 + k_1) < (double) 0.0)
  {
   # 52 "update_centroids_of_velocities_and_height.c"
   *(h_B_12 + k_1) = (double) 0.0;
  }
  # 55 "update_centroids_of_velocities_and_height.c"
  factor_2 = *(h_B_12 + k_1) / (*(h_B_12 + k_1) * *(h_B_12 + k_1) + (double) 1.000000000000000021e-08);
  # 56 "update_centroids_of_velocities_and_height.c"
  *(u_B_11 + k_1) = *(uh_B_12 + k_1) * factor_2;
  # 57 "update_centroids_of_velocities_and_height.c"
  *(v_B_1 + k_1) = *(vh_B_12 + k_1) * factor_2;
 }
 # 8 "update_centroids_of_velocities_and_height.c"
 __hmppcg_label_1:;
}
#endif // __CUDACC__



# 8 "update_centroids_of_velocities_and_height.c"

#ifndef __CUDACC__
void _update_centroids_of_velocities_and_height_internal_1(hmpprt::s32 N_c_12, hmpprt::s32 N_b_11, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  w_C_11, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  uh_C_1, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  vh_C_11, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  h_C_1, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  z_C_1, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  u_C_11, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  v_C_1, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  w_B_1, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  uh_B_1, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  vh_B_1, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  h_B_1, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  z_B_1, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  u_B_1, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  v_B_11)
{
 # 8 "update_centroids_of_velocities_and_height.c"
 if (N_c_12 - 1 >= 0)
 {
  hmpprt::CUDAGridCall __hmppcg_call;
  __hmppcg_call.setSizeX((N_c_12 - 1) / 128 + 1);
  __hmppcg_call.setSizeY(1);
  __hmppcg_call.setBlockSizeX(32);
  __hmppcg_call.setBlockSizeY(4);
  __hmppcg_call.addLocalParameter((hmpprt::s32) (N_c_12), "N_c_1");
  __hmppcg_call.addLocalParameter((hmpprt::s32) (N_b_11), "N_b_12");
  __hmppcg_call.addLocalParameter(&w_C_11, 8, "w_C_12");
  __hmppcg_call.addLocalParameter(&uh_C_1, 8, "uh_C_12");
  __hmppcg_call.addLocalParameter(&vh_C_11, 8, "vh_C_12");
  __hmppcg_call.addLocalParameter(&h_C_1, 8, "h_C_11");
  __hmppcg_call.addLocalParameter(&z_C_1, 8, "z_C_11");
  __hmppcg_call.addLocalParameter(&u_C_11, 8, "u_C_12");
  __hmppcg_call.addLocalParameter(&v_C_1, 8, "v_C_12");
  __hmppcg_call.addLocalParameter(&w_B_1, 8, "w_B_12");
  __hmppcg_call.addLocalParameter(&uh_B_1, 8, "uh_B_12");
  __hmppcg_call.addLocalParameter(&vh_B_1, 8, "vh_B_12");
  __hmppcg_call.addLocalParameter(&h_B_1, 8, "h_B_12");
  __hmppcg_call.addLocalParameter(&z_B_1, 8, "z_B_12");
  __hmppcg_call.addLocalParameter(&u_B_1, 8, "u_B_11");
  __hmppcg_call.addLocalParameter(&v_B_11, 8, "v_B_1");
  __hmppcg_call.launch(_update_centroids_of_velocities_and_height_loop1D_1, hmpprt::Context::getInstance()->getCUDADevice());
 }
 ;
}
#endif // __CUDACC__



# 8 "update_centroids_of_velocities_and_height.c"

#ifndef __CUDACC__
extern "C" CDLT_API  void _update_centroids_of_velocities_and_height(hmpprt::s32 N_c_11, hmpprt::s32 N_b_1, double* w_C_1, double* uh_C_11, double* vh_C_1, double* h_C_12, double* z_C_12, double* u_C_1, double* v_C_11, double* w_B_11, double* uh_B_11, double* vh_B_11, double* h_B_11, double* z_B_11, double* u_B_12, double* v_B_12)
{
 # 1 "<preprocessor>"
 (_update_centroids_of_velocities_and_height_internal_1(N_c_11, N_b_1, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (w_C_1), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (uh_C_11), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (vh_C_1), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (h_C_12), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (z_C_12), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (u_C_1), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (v_C_11), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (w_B_11), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (uh_B_11), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (vh_B_11), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (h_B_11), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (z_B_11), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (u_B_12), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (v_B_12)));
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
      _update_centroids_of_velocities_and_height_loop1D_1 = new hmpprt::CUDAGrid(hmpprt_module, "_update_centroids_of_velocities_and_height_loop1D_1");

    }
    hmpprt::Context::getInstance()->getGrouplet()->setTarget(hmpprt::CUDA);
    hmpprt::Context::getInstance()->getGrouplet()->addSignature("_update_centroids_of_velocities_and_height", "prototype _update_centroids_of_velocities_and_height(N_c: s32, N_b: s32, w_C: ^cudaglob double, uh_C: ^cudaglob double, vh_C: ^cudaglob double, h_C: ^cudaglob double, z_C: ^cudaglob double, u_C: ^cudaglob double, v_C: ^cudaglob double, w_B: ^cudaglob double, uh_B: ^cudaglob double, vh_B: ^cudaglob double, h_B: ^cudaglob double, z_B: ^cudaglob double, u_B: ^cudaglob double, v_B: ^cudaglob double)");

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
      delete _update_centroids_of_velocities_and_height_loop1D_1;

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

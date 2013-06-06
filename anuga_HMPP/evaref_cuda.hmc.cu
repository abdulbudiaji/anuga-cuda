
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

# 9 "evaluate_segment.c"

#ifndef __CUDACC__
extern "C" CDLT_API  void evaluate_segment_reflective(hmpprt::s32 N1_21, hmpprt::s32 N2_2, hmpprt::s32 N3, hmpprt::s32 N6, hmpprt::s64* ids, hmpprt::s64* vol_ids_21, hmpprt::s64* edge_ids, double* normals_2, double* stage_edge_values, double* bed_edge_values, double* height_edge_values, double* xmom_edge_values, double* ymom_edge_values, double* xvel_edge_values, double* yvel_edge_values, double* stage_boundary_values_21, double* bed_boundary_values, double* height_boundary_values, double* xmom_boundary_values_21, double* ymom_boundary_values, double* xvel_boundary_values, double* yvel_boundary_values)
;
#endif // __CUDACC__



# 9 "evaluate_segment.c"

#ifndef __CUDACC__
void evaluate_segment_reflective_internal_1(hmpprt::s32 N1_2, hmpprt::s32 N2_21, hmpprt::s32 N3, hmpprt::s32 N6, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s64>  ids, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s64>  vol_ids_2, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s64>  edge_ids, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  normals_21, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  stage_edge_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  bed_edge_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  height_edge_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  xmom_edge_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  ymom_edge_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  xvel_edge_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  yvel_edge_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  stage_boundary_values_2, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  bed_boundary_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  height_boundary_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  xmom_boundary_values_2, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  ymom_boundary_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  xvel_boundary_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  yvel_boundary_values)
;
#endif // __CUDACC__



# 9 "evaluate_segment.c"

#ifndef __CUDACC__
static hmpprt::CUDAGrid * evaluate_segment_reflective_loop1D_1 = 0;
#else

extern "C" __global__ void evaluate_segment_reflective_loop1D_1(hmpprt::s32 N1, hmpprt::s64* ids_2, hmpprt::s64* vol_ids, hmpprt::s64* edge_ids_2, double* normals, double* stage_edge_values_2, double* bed_edge_values_2, double* height_edge_values_2, double* xmom_edge_values_2, double* ymom_edge_values_2, double* xvel_edge_values_2, double* yvel_edge_values_2, double* stage_boundary_values, double* bed_boundary_values_2, double* height_boundary_values_2, double* xmom_boundary_values, double* ymom_boundary_values_2, double* xvel_boundary_values_2, double* yvel_boundary_values_2);
#endif // __CUDACC__




# 9 "evaluate_segment.c"

#ifdef __CUDACC__

extern "C" __global__ void evaluate_segment_reflective_loop1D_1(hmpprt::s32 N1, hmpprt::s64* ids_2, hmpprt::s64* vol_ids, hmpprt::s64* edge_ids_2, double* normals, double* stage_edge_values_2, double* bed_edge_values_2, double* height_edge_values_2, double* xmom_edge_values_2, double* ymom_edge_values_2, double* xvel_edge_values_2, double* yvel_edge_values_2, double* stage_boundary_values, double* bed_boundary_values_2, double* height_boundary_values_2, double* xmom_boundary_values, double* ymom_boundary_values_2, double* xvel_boundary_values_2, double* yvel_boundary_values_2)
{
 # 37 "evaluate_segment.c"
 hmpprt::s32 id_vol_1;
 # 37 "evaluate_segment.c"
 hmpprt::s32 id_edge_1;
 # 37 "evaluate_segment.c"
 double q1_1;
 # 37 "evaluate_segment.c"
 double n1_1;
 # 37 "evaluate_segment.c"
 double q2_1;
 # 37 "evaluate_segment.c"
 double n2_1;
 # 37 "evaluate_segment.c"
 double r1_1;
 # 37 "evaluate_segment.c"
 double r2_1;
 # 37 "evaluate_segment.c"
 double q1_2;
 # 37 "evaluate_segment.c"
 double q2_2;
 # 37 "evaluate_segment.c"
 double r1_2;
 # 37 "evaluate_segment.c"
 double r2_2;
 # 37 "evaluate_segment.c"
 hmpprt::s32 k_1;
 # 56 "evaluate_segment.c"
 k_1 = (hmpprt::gr_atidf());
 # 56 "evaluate_segment.c"
 if (k_1 > N1 - 1)
 {
  # 56 "evaluate_segment.c"
  goto __hmppcg_label_1;
 }
 # 56 "evaluate_segment.c"
 id_vol_1 = (hmpprt::s32 ) (*(vol_ids + k_1));
 # 57 "evaluate_segment.c"
 id_edge_1 = (hmpprt::s32 ) (*(edge_ids_2 + k_1));
 # 60 "evaluate_segment.c"
 n1_1 = *(normals + (id_vol_1 * 6 + id_edge_1 * 2));
 # 61 "evaluate_segment.c"
 n2_1 = *(normals + (id_vol_1 * 6 + id_edge_1 * 2 + 1));
 # 62 "evaluate_segment.c"
 q1_1 = *(xmom_edge_values_2 + (id_vol_1 * 3 + id_edge_1));
 # 63 "evaluate_segment.c"
 q2_1 = *(ymom_edge_values_2 + (id_vol_1 * 3 + id_edge_1));
 # 71 "evaluate_segment.c"
 r1_1 =  - q1_1 * n1_1 - q2_1 * n2_1;
 # 72 "evaluate_segment.c"
 r2_1 =  - q1_1 * n2_1 + q2_1 * n1_1;
 # 75 "evaluate_segment.c"
 *(stage_boundary_values + k_1) = *(stage_edge_values_2 + (id_vol_1 * 3 + id_edge_1));
 # 76 "evaluate_segment.c"
 *(bed_boundary_values_2 + k_1) = *(bed_edge_values_2 + (id_vol_1 * 3 + id_edge_1));
 # 77 "evaluate_segment.c"
 *(height_boundary_values_2 + k_1) = *(height_edge_values_2 + (id_vol_1 * 3 + id_edge_1));
 # 79 "evaluate_segment.c"
 q1_2 = *(xvel_edge_values_2 + (id_vol_1 * 3 + id_edge_1));
 # 80 "evaluate_segment.c"
 q2_2 = *(yvel_edge_values_2 + (id_vol_1 * 3 + id_edge_1));
 # 90 "evaluate_segment.c"
 *(xmom_boundary_values + k_1) = n1_1 * r1_1 - n2_1 * r2_1;
 # 91 "evaluate_segment.c"
 *(ymom_boundary_values_2 + k_1) = n2_1 * r1_1 + n1_1 * r2_1;
 # 94 "evaluate_segment.c"
 r1_2 = q1_2 * n1_1 + q2_2 * n2_1;
 # 95 "evaluate_segment.c"
 r2_2 = q1_2 * n2_1 - q2_2 * n1_1;
 # 97 "evaluate_segment.c"
 *(xvel_boundary_values_2 + k_1) = n1_1 * r1_2 - n2_1 * r2_2;
 # 98 "evaluate_segment.c"
 *(yvel_boundary_values_2 + k_1) = n2_1 * r1_2 + n1_1 * r2_2;
 # 100 "evaluate_segment.c"
 *(vol_ids + k_1) = *(ids_2 + k_1);
 # 9 "evaluate_segment.c"
 __hmppcg_label_1:;
}
#endif // __CUDACC__



# 9 "evaluate_segment.c"

#ifndef __CUDACC__
void evaluate_segment_reflective_internal_1(hmpprt::s32 N1_2, hmpprt::s32 N2_21, hmpprt::s32 N3, hmpprt::s32 N6, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s64>  ids, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s64>  vol_ids_2, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s64>  edge_ids, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  normals_21, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  stage_edge_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  bed_edge_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  height_edge_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  xmom_edge_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  ymom_edge_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  xvel_edge_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  yvel_edge_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  stage_boundary_values_2, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  bed_boundary_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  height_boundary_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  xmom_boundary_values_2, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  ymom_boundary_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  xvel_boundary_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  yvel_boundary_values)
{
 # 9 "evaluate_segment.c"
 if (N1_2 - 1 >= 0)
 {
  hmpprt::CUDAGridCall __hmppcg_call;
  __hmppcg_call.setSizeX((N1_2 - 1) / 128 + 1);
  __hmppcg_call.setSizeY(1);
  __hmppcg_call.setBlockSizeX(32);
  __hmppcg_call.setBlockSizeY(4);
  __hmppcg_call.addLocalParameter((hmpprt::s32) (N1_2), "N1");
  __hmppcg_call.addLocalParameter(&ids, 8, "ids_2");
  __hmppcg_call.addLocalParameter(&vol_ids_2, 8, "vol_ids");
  __hmppcg_call.addLocalParameter(&edge_ids, 8, "edge_ids_2");
  __hmppcg_call.addLocalParameter(&normals_21, 8, "normals");
  __hmppcg_call.addLocalParameter(&stage_edge_values, 8, "stage_edge_values_2");
  __hmppcg_call.addLocalParameter(&bed_edge_values, 8, "bed_edge_values_2");
  __hmppcg_call.addLocalParameter(&height_edge_values, 8, "height_edge_values_2");
  __hmppcg_call.addLocalParameter(&xmom_edge_values, 8, "xmom_edge_values_2");
  __hmppcg_call.addLocalParameter(&ymom_edge_values, 8, "ymom_edge_values_2");
  __hmppcg_call.addLocalParameter(&xvel_edge_values, 8, "xvel_edge_values_2");
  __hmppcg_call.addLocalParameter(&yvel_edge_values, 8, "yvel_edge_values_2");
  __hmppcg_call.addLocalParameter(&stage_boundary_values_2, 8, "stage_boundary_values");
  __hmppcg_call.addLocalParameter(&bed_boundary_values, 8, "bed_boundary_values_2");
  __hmppcg_call.addLocalParameter(&height_boundary_values, 8, "height_boundary_values_2");
  __hmppcg_call.addLocalParameter(&xmom_boundary_values_2, 8, "xmom_boundary_values");
  __hmppcg_call.addLocalParameter(&ymom_boundary_values, 8, "ymom_boundary_values_2");
  __hmppcg_call.addLocalParameter(&xvel_boundary_values, 8, "xvel_boundary_values_2");
  __hmppcg_call.addLocalParameter(&yvel_boundary_values, 8, "yvel_boundary_values_2");
  __hmppcg_call.launch(evaluate_segment_reflective_loop1D_1, hmpprt::Context::getInstance()->getCUDADevice());
 }
 ;
}
#endif // __CUDACC__



# 9 "evaluate_segment.c"

#ifndef __CUDACC__
extern "C" CDLT_API  void evaluate_segment_reflective(hmpprt::s32 N1_21, hmpprt::s32 N2_2, hmpprt::s32 N3, hmpprt::s32 N6, hmpprt::s64* ids, hmpprt::s64* vol_ids_21, hmpprt::s64* edge_ids, double* normals_2, double* stage_edge_values, double* bed_edge_values, double* height_edge_values, double* xmom_edge_values, double* ymom_edge_values, double* xvel_edge_values, double* yvel_edge_values, double* stage_boundary_values_21, double* bed_boundary_values, double* height_boundary_values, double* xmom_boundary_values_21, double* ymom_boundary_values, double* xvel_boundary_values, double* yvel_boundary_values)
{
 # 1 "<preprocessor>"
 (evaluate_segment_reflective_internal_1(N1_21, N2_2, N3, N6, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s64> (ids), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s64> (vol_ids_21), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s64> (edge_ids), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (normals_2), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (stage_edge_values), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (bed_edge_values), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (height_edge_values), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (xmom_edge_values), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (ymom_edge_values), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (xvel_edge_values), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (yvel_edge_values), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (stage_boundary_values_21), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (bed_boundary_values), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (height_boundary_values), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (xmom_boundary_values_21), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (ymom_boundary_values), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (xvel_boundary_values), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (yvel_boundary_values)));
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
      evaluate_segment_reflective_loop1D_1 = new hmpprt::CUDAGrid(hmpprt_module, "evaluate_segment_reflective_loop1D_1");

    }
    hmpprt::Context::getInstance()->getGrouplet()->setTarget(hmpprt::CUDA);
    hmpprt::Context::getInstance()->getGrouplet()->addSignature("evaluate_segment_reflective", "prototype evaluate_segment_reflective(N1: s32, N2: s32, N3: s32, N6: s32, ids: ^cudaglob s64, vol_ids: ^cudaglob s64, edge_ids: ^cudaglob s64, normals: ^cudaglob double, stage_edge_values: ^cudaglob double, bed_edge_values: ^cudaglob double, height_edge_values: ^cudaglob double, xmom_edge_values: ^cudaglob double, ymom_edge_values: ^cudaglob double, xvel_edge_values: ^cudaglob double, yvel_edge_values: ^cudaglob double, stage_boundary_values: ^cudaglob double, bed_boundary_values: ^cudaglob double, height_boundary_values: ^cudaglob double, xmom_boundary_values: ^cudaglob double, ymom_boundary_values: ^cudaglob double, xvel_boundary_values: ^cudaglob double, yvel_boundary_values: ^cudaglob double)");

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
      delete evaluate_segment_reflective_loop1D_1;

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

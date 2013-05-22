
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

# 49 "gravity.c"

#ifndef __CUDACC__
extern "C" CDLT_API  void gravity_wb(hmpprt::s32 n, hmpprt::s32 n3, hmpprt::s32 n6, double* xmom_explicit_update, double* ymom_explicit_update, double* stage_vertex_values, double* stage_edge_values, double* stage_centroid_values, double* bed_edge_values, double* bed_centroid_values, double* vertex_coordinates, double* normals, double* areas, double* edgelengths, double g)
;
#endif // __CUDACC__



# 49 "gravity.c"

#ifndef __CUDACC__
void gravity_wb_internal_1(hmpprt::s32 n, hmpprt::s32 n3, hmpprt::s32 n6, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  xmom_explicit_update, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  ymom_explicit_update, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  stage_vertex_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  stage_edge_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  stage_centroid_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  bed_edge_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  bed_centroid_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  vertex_coordinates, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  normals, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  areas, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  edgelengths, double g)
;
#endif // __CUDACC__



# 49 "gravity.c"

#ifndef __CUDACC__
static hmpprt::CUDAGrid * gravity_wb_loop1D_1 = 0;
#else

extern "C" __global__ void gravity_wb_loop1D_1(hmpprt::s32 n, double* xmom_explicit_update, double* ymom_explicit_update, double* stage_vertex_values, double* stage_edge_values, double* stage_centroid_values, double* bed_edge_values, double* bed_centroid_values, double* vertex_coordinates, double* normals, double* areas, double* edgelengths, double g);
#endif // __CUDACC__




# 49 "gravity.c"

#ifdef __CUDACC__

extern "C" __global__ void gravity_wb_loop1D_1(hmpprt::s32 n, double* xmom_explicit_update, double* ymom_explicit_update, double* stage_vertex_values, double* stage_edge_values, double* stage_centroid_values, double* bed_edge_values, double* bed_centroid_values, double* vertex_coordinates, double* normals, double* areas, double* edgelengths, double g)
{
 # 69 "gravity.c"
 hmpprt::s32 k3_1;
 # 69 "gravity.c"
 hmpprt::s32 k6_1;
 # 69 "gravity.c"
 double y2_1;
 # 69 "gravity.c"
 double y0_1;
 # 69 "gravity.c"
 double x1_1;
 # 69 "gravity.c"
 double x0_1;
 # 69 "gravity.c"
 double y1_1;
 # 69 "gravity.c"
 double x2_1;
 # 69 "gravity.c"
 double w1_1;
 # 69 "gravity.c"
 double w0_1;
 # 69 "gravity.c"
 double w2_1;
 # 69 "gravity.c"
 double wx_1;
 # 69 "gravity.c"
 double det_1;
 # 69 "gravity.c"
 double wy_1;
 # 69 "gravity.c"
 double wx_2;
 # 69 "gravity.c"
 double avg_h_1;
 # 69 "gravity.c"
 double wy_2;
 # 69 "gravity.c"
 double hh0_1;
 # 69 "gravity.c"
 double hh1_1;
 # 69 "gravity.c"
 double hh2_1;
 # 69 "gravity.c"
 double sidex_0_1;
 # 69 "gravity.c"
 double sidex_1_1;
 # 69 "gravity.c"
 double sidex_2_1;
 # 69 "gravity.c"
 double area_1;
 # 69 "gravity.c"
 double sidey_0_1;
 # 69 "gravity.c"
 double sidey_1_1;
 # 69 "gravity.c"
 double sidey_2_1;
 # 69 "gravity.c"
 hmpprt::s32 k_1;
 # 84 "gravity.c"
 k_1 = (hmpprt::gr_atidf());
 # 84 "gravity.c"
 if (k_1 > n - 1)
 {
  # 84 "gravity.c"
  goto __hmppcg_label_1;
 }
 # 84 "gravity.c"
 k3_1 = k_1 * 3;
 # 85 "gravity.c"
 w0_1 = *(stage_vertex_values + k3_1);
 # 86 "gravity.c"
 w1_1 = *(stage_vertex_values + (k3_1 + 1));
 # 87 "gravity.c"
 w2_1 = *(stage_vertex_values + (k3_1 + 2));
 # 89 "gravity.c"
 k6_1 = k_1 * 6;
 # 90 "gravity.c"
 x0_1 = *(vertex_coordinates + k6_1);
 # 91 "gravity.c"
 y0_1 = *(vertex_coordinates + (k6_1 + 1));
 # 92 "gravity.c"
 x1_1 = *(vertex_coordinates + (k6_1 + 2));
 # 93 "gravity.c"
 y1_1 = *(vertex_coordinates + (k6_1 + 3));
 # 94 "gravity.c"
 x2_1 = *(vertex_coordinates + (k6_1 + 4));
 # 95 "gravity.c"
 y2_1 = *(vertex_coordinates + (k6_1 + 5));
 # 98 "gravity.c"
 det_1 = (y2_1 - y0_1) * (x1_1 - x0_1) - (y1_1 - y0_1) * (x2_1 - x0_1);
 # 100 "gravity.c"
 wx_1 = (y2_1 - y0_1) * (w1_1 - w0_1) - (y1_1 - y0_1) * (w2_1 - w0_1);
 # 101 "gravity.c"
 wx_2 = wx_1 / det_1;
 # 103 "gravity.c"
 wy_1 = (x1_1 - x0_1) * (w2_1 - w0_1) - (x2_1 - x0_1) * (w1_1 - w0_1);
 # 104 "gravity.c"
 wy_2 = wy_1 / det_1;
 # 106 "gravity.c"
 avg_h_1 = *(stage_centroid_values + k_1) - *(bed_centroid_values + k_1);
 # 108 "gravity.c"
 *(xmom_explicit_update + k_1) = *(xmom_explicit_update + k_1) +  - g * wx_2 * avg_h_1;
 # 109 "gravity.c"
 *(ymom_explicit_update + k_1) = *(ymom_explicit_update + k_1) +  - g * wy_2 * avg_h_1;
 # 114 "gravity.c"
 hh0_1 = *(stage_edge_values + k3_1) - *(bed_edge_values + k3_1);
 # 115 "gravity.c"
 hh1_1 = *(stage_edge_values + (k3_1 + 1)) - *(bed_edge_values + (k3_1 + 1));
 # 116 "gravity.c"
 hh2_1 = *(stage_edge_values + (k3_1 + 2)) - *(bed_edge_values + (k3_1 + 2));
 # 118 "gravity.c"
 area_1 = *(areas + k_1);
 # 120 "gravity.c"
 sidex_0_1 = hh0_1 * hh0_1 * *(edgelengths + k3_1) * *(normals + k6_1);
 # 121 "gravity.c"
 sidex_1_1 = hh1_1 * hh1_1 * *(edgelengths + (k3_1 + 1)) * *(normals + (k6_1 + 2));
 # 122 "gravity.c"
 sidex_2_1 = hh2_1 * hh2_1 * *(edgelengths + (k3_1 + 2)) * *(normals + (k6_1 + 4));
 # 123 "gravity.c"
 *(xmom_explicit_update + k_1) = *(xmom_explicit_update + k_1) + (double) 0.5 * g * (sidex_0_1 + sidex_1_1 + sidex_2_1) / area_1;
 # 126 "gravity.c"
 sidey_0_1 = hh0_1 * hh0_1 * *(edgelengths + k3_1) * *(normals + (k6_1 + 1));
 # 127 "gravity.c"
 sidey_1_1 = hh1_1 * hh1_1 * *(edgelengths + (k3_1 + 1)) * *(normals + (k6_1 + 3));
 # 128 "gravity.c"
 sidey_2_1 = hh2_1 * hh2_1 * *(edgelengths + (k3_1 + 2)) * *(normals + (k6_1 + 5));
 # 129 "gravity.c"
 *(ymom_explicit_update + k_1) = *(ymom_explicit_update + k_1) + (double) 0.5 * g * (sidey_0_1 + sidey_1_1 + sidey_2_1) / area_1;
 # 49 "gravity.c"
 __hmppcg_label_1:;
}
#endif // __CUDACC__



# 49 "gravity.c"

#ifndef __CUDACC__
void gravity_wb_internal_1(hmpprt::s32 n, hmpprt::s32 n3, hmpprt::s32 n6, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  xmom_explicit_update, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  ymom_explicit_update, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  stage_vertex_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  stage_edge_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  stage_centroid_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  bed_edge_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  bed_centroid_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  vertex_coordinates, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  normals, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  areas, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  edgelengths, double g)
{
 # 49 "gravity.c"
 if (n - 1 >= 0)
 {
  hmpprt::CUDAGridCall __hmppcg_call;
  __hmppcg_call.setSizeX((n - 1) / 128 + 1);
  __hmppcg_call.setSizeY(1);
  __hmppcg_call.setBlockSizeX(32);
  __hmppcg_call.setBlockSizeY(4);
  __hmppcg_call.addLocalParameter((hmpprt::s32) (n), "n");
  __hmppcg_call.addLocalParameter(&xmom_explicit_update, 8, "xmom_explicit_update");
  __hmppcg_call.addLocalParameter(&ymom_explicit_update, 8, "ymom_explicit_update");
  __hmppcg_call.addLocalParameter(&stage_vertex_values, 8, "stage_vertex_values");
  __hmppcg_call.addLocalParameter(&stage_edge_values, 8, "stage_edge_values");
  __hmppcg_call.addLocalParameter(&stage_centroid_values, 8, "stage_centroid_values");
  __hmppcg_call.addLocalParameter(&bed_edge_values, 8, "bed_edge_values");
  __hmppcg_call.addLocalParameter(&bed_centroid_values, 8, "bed_centroid_values");
  __hmppcg_call.addLocalParameter(&vertex_coordinates, 8, "vertex_coordinates");
  __hmppcg_call.addLocalParameter(&normals, 8, "normals");
  __hmppcg_call.addLocalParameter(&areas, 8, "areas");
  __hmppcg_call.addLocalParameter(&edgelengths, 8, "edgelengths");
  __hmppcg_call.addLocalParameter(&g, 8, "g");
  __hmppcg_call.launch(gravity_wb_loop1D_1, hmpprt::Context::getInstance()->getCUDADevice());
 }
 ;
}
#endif // __CUDACC__



# 49 "gravity.c"

#ifndef __CUDACC__
extern "C" CDLT_API  void gravity_wb(hmpprt::s32 n, hmpprt::s32 n3, hmpprt::s32 n6, double* xmom_explicit_update, double* ymom_explicit_update, double* stage_vertex_values, double* stage_edge_values, double* stage_centroid_values, double* bed_edge_values, double* bed_centroid_values, double* vertex_coordinates, double* normals, double* areas, double* edgelengths, double g)
{
 # 1 "<preprocessor>"
 (gravity_wb_internal_1(n, n3, n6, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (xmom_explicit_update), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (ymom_explicit_update), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (stage_vertex_values), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (stage_edge_values), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (stage_centroid_values), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (bed_edge_values), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (bed_centroid_values), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (vertex_coordinates), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (normals), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (areas), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (edgelengths), g));
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
      gravity_wb_loop1D_1 = new hmpprt::CUDAGrid(hmpprt_module, "gravity_wb_loop1D_1");

    }
    hmpprt::Context::getInstance()->getGrouplet()->setTarget(hmpprt::CUDA);
    hmpprt::Context::getInstance()->getGrouplet()->addSignature("gravity_wb", "prototype gravity_wb(n: s32, n3: s32, n6: s32, xmom_explicit_update: ^cudaglob double, ymom_explicit_update: ^cudaglob double, stage_vertex_values: ^cudaglob double, stage_edge_values: ^cudaglob double, stage_centroid_values: ^cudaglob double, bed_edge_values: ^cudaglob double, bed_centroid_values: ^cudaglob double, vertex_coordinates: ^cudaglob double, normals: ^cudaglob double, areas: ^cudaglob double, edgelengths: ^cudaglob double, g: double)");

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
      delete gravity_wb_loop1D_1;

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

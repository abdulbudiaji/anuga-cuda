
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

# 290 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"

#ifndef __CUDACC__
extern "C" CDLT_API  void _limit_vertices_by_all_neighbours(hmpprt::s32 N_21, hmpprt::s32 N3, double beta, double* centroid_values, double* vertex_values, double* edge_values_2, hmpprt::s64* neighbours_21, double* x_gradient_21, double* y_gradient_21)
;
#endif // __CUDACC__



# 290 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"

#ifndef __CUDACC__
void _limit_vertices_by_all_neighbours_internal_1(hmpprt::s32 N_2, hmpprt::s32 N3, double beta, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  centroid_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  vertex_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  edge_values_21, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s64>  neighbours_2, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  x_gradient_2, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  y_gradient_2)
;
#endif // __CUDACC__



# 14 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"

#ifndef __CUDACC__
static hmpprt::CUDAGrid * _limit_vertices_by_all_neighbours_loop1D_1 = 0;
#else

extern "C" __global__ void _limit_vertices_by_all_neighbours_loop1D_1(hmpprt::s32 N, double beta_2, double* centroid_values_2, double* vertex_values_2, double* edge_values, hmpprt::s64* neighbours, double* x_gradient, double* y_gradient);
#endif // __CUDACC__




# 14 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"

#ifdef __CUDACC__
__device__ double fmax_cudalocal_1(double x, double y)
;
#endif // __CUDACC__



# 22 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"

#ifdef __CUDACC__
__device__ double fmin_cudalocal_1(double x_3, double y_3)
;
#endif // __CUDACC__



# 22 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"

#ifndef __CUDACC__
extern "C" CDLT_API  double fmin(double x_31, double y_32)
;
#endif // __CUDACC__



# 22 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"

#ifndef __CUDACC__
double fmin_internal_1(double x_32, double y_31)
;
#endif // __CUDACC__



# 14 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"

#ifndef __CUDACC__
extern "C" CDLT_API  double fmax(double x, double y)
;
#endif // __CUDACC__



# 14 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"

#ifndef __CUDACC__
double fmax_internal_1(double x, double y)
;
#endif // __CUDACC__



# 14 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"

#ifndef __CUDACC__
double fmax_internal_1(double x, double y)
{
 # 17 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 double rvalue_11;
 # 17 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 if (x > y)
 {
  # 17 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  rvalue_11 = x;
  # 17 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  goto endf_1;
 }
 else
 {
  # 14 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  rvalue_11 = y;
  # 14 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  goto endf_1;
 }
 # 14 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 endf_1:;
 # 14 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 return rvalue_11;
}
#endif // __CUDACC__



# 14 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"

#ifndef __CUDACC__
extern "C" CDLT_API  double fmax(double x, double y)
{
 # 22 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 (fmax_internal_1(x, y));
}
#endif // __CUDACC__



# 22 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"

#ifndef __CUDACC__
double fmin_internal_1(double x_32, double y_31)
{
 # 25 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 double rvalue_21;
 # 25 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 if (x_32 < y_31)
 {
  # 25 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  rvalue_21 = x_32;
  # 25 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  goto endf_2;
 }
 else
 {
  # 22 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  rvalue_21 = y_31;
  # 22 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  goto endf_2;
 }
 # 22 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 endf_2:;
 # 22 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 return rvalue_21;
}
#endif // __CUDACC__



# 22 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"

#ifndef __CUDACC__
extern "C" CDLT_API  double fmin(double x_31, double y_32)
{
 # 22 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 (fmin_internal_1(x_31, y_32));
}
#endif // __CUDACC__



# 22 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"

#ifdef __CUDACC__
__device__ double fmin_cudalocal_1(double x_3, double y_3)
{
 # 25 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 double rvalue_2;
 # 25 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 if (x_3 < y_3)
 {
  # 25 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  rvalue_2 = x_3;
  # 25 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  goto endf_21;
 }
 else
 {
  # 14 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  rvalue_2 = y_3;
  # 14 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  goto endf_21;
 }
 # 14 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 endf_21:;
 # 14 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 return rvalue_2;
}
#endif // __CUDACC__



# 14 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"

#ifdef __CUDACC__
__device__ double fmax_cudalocal_1(double x, double y)
{
 # 17 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 double rvalue_1;
 # 17 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 if (x > y)
 {
  # 17 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  rvalue_1 = x;
  # 17 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  goto endf_11;
 }
 else
 {
  # 290 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  rvalue_1 = y;
  # 290 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  goto endf_11;
 }
 # 290 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 endf_11:;
 # 290 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 return rvalue_1;
}
#endif // __CUDACC__



# 290 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"

#ifdef __CUDACC__

extern "C" __global__ void _limit_vertices_by_all_neighbours_loop1D_1(hmpprt::s32 N, double beta_2, double* centroid_values_2, double* vertex_values_2, double* edge_values, hmpprt::s64* neighbours, double* x_gradient, double* y_gradient)
{
 # 304 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 double dqa_1[3uLL];
 # 301 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 double qc_1;
 # 301 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 hmpprt::s32 k3_1;
 # 301 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 double qmin_1;
 # 301 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 double qmax_1;
 # 301 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 double phi_1;
 # 301 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 hmpprt::s32 k_1;
 # 313 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 k_1 = (hmpprt::gr_atidf());
 # 313 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 if (k_1 > N - 1)
 {
  # 313 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  goto __hmppcg_label_1;
 }
 # 313 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 k3_1 = 3 * k_1;
 # 315 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 qc_1 = *(centroid_values_2 + k_1);
 # 316 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 qmin_1 = qc_1;
 # 317 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 qmax_1 = qc_1;
 # 320 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 hmpprt::s32 i_1;
 # 320 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 # 320 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 for (i_1 = 0 ; i_1 <= 2 ; i_1 = i_1 + 1)
 {
  # 321 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  hmpprt::s64 n_1;
  # 321 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  n_1 = *(neighbours + (k3_1 + i_1));
  # 322 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  if (n_1 >= 0LL)
  {
   # 323 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
   double qn_1;
   # 323 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
   qn_1 = *(centroid_values_2 + n_1);
   # 325 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
   qmin_1 = (fmin_cudalocal_1(qmin_1, qn_1));
   # 326 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
   qmax_1 = (fmax_cudalocal_1(qmax_1, qn_1));
  }
 }
 # 330 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 # 330 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 phi_1 = (double) 1.0;
 # 332 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 hmpprt::s32 i_2;
 # 332 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 # 332 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 for (i_2 = 0 ; i_2 <= 2 ; i_2 = i_2 + 1)
 {
  # 333 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  double dq_1;
  # 333 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  double r_1;
  # 333 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  double returntype_2;
  # 333 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  r_1 = (double) 1.0;
  # 335 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  dq_1 = *(vertex_values_2 + (k3_1 + i_2)) - qc_1;
  # 336 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  dqa_1[i_2] = dq_1;
  # 338 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  if (dq_1 > (double) 0.0)
  {
   # 338 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
   r_1 = (qmax_1 - qc_1) / dq_1;
  }
  # 339 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  if (dq_1 < (double) 0.0)
  {
   # 339 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
   r_1 = (qmin_1 - qc_1) / dq_1;
  }
  # 342 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  returntype_2 = (fmin_cudalocal_1(r_1 * beta_2, (double) 1.0));
  # 342 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  phi_1 = (fmin_cudalocal_1(returntype_2, phi_1));
 }
 # 346 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 # 346 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 *(x_gradient + k_1) = *(x_gradient + k_1) * phi_1;
 # 347 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 *(y_gradient + k_1) = *(y_gradient + k_1) * phi_1;
 # 349 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 *(vertex_values_2 + k3_1) = qc_1 + phi_1 * dqa_1[0];
 # 350 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 *(vertex_values_2 + (k3_1 + 1)) = qc_1 + phi_1 * dqa_1[1];
 # 351 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 *(vertex_values_2 + (k3_1 + 2)) = qc_1 + phi_1 * dqa_1[2];
 # 353 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 *(edge_values + k3_1) = (double) 0.5 * (*(vertex_values_2 + (k3_1 + 1)) + *(vertex_values_2 + (k3_1 + 2)));
 # 354 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 *(edge_values + (k3_1 + 1)) = (double) 0.5 * (*(vertex_values_2 + (k3_1 + 2)) + *(vertex_values_2 + k3_1));
 # 355 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 *(edge_values + (k3_1 + 2)) = (double) 0.5 * (*(vertex_values_2 + k3_1) + *(vertex_values_2 + (k3_1 + 1)));
 # 290 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 __hmppcg_label_1:;
}
#endif // __CUDACC__



# 290 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"

#ifndef __CUDACC__
void _limit_vertices_by_all_neighbours_internal_1(hmpprt::s32 N_2, hmpprt::s32 N3, double beta, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  centroid_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  vertex_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  edge_values_21, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s64>  neighbours_2, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  x_gradient_2, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  y_gradient_2)
{
 # 290 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 if (N_2 - 1 >= 0)
 {
  hmpprt::CUDAGridCall __hmppcg_call;
  __hmppcg_call.setSizeX((N_2 - 1) / 128 + 1);
  __hmppcg_call.setSizeY(1);
  __hmppcg_call.setBlockSizeX(32);
  __hmppcg_call.setBlockSizeY(4);
  __hmppcg_call.addLocalParameter((hmpprt::s32) (N_2), "N");
  __hmppcg_call.addLocalParameter(&beta, 8, "beta_2");
  __hmppcg_call.addLocalParameter(&centroid_values, 8, "centroid_values_2");
  __hmppcg_call.addLocalParameter(&vertex_values, 8, "vertex_values_2");
  __hmppcg_call.addLocalParameter(&edge_values_21, 8, "edge_values");
  __hmppcg_call.addLocalParameter(&neighbours_2, 8, "neighbours");
  __hmppcg_call.addLocalParameter(&x_gradient_2, 8, "x_gradient");
  __hmppcg_call.addLocalParameter(&y_gradient_2, 8, "y_gradient");
  __hmppcg_call.launch(_limit_vertices_by_all_neighbours_loop1D_1, hmpprt::Context::getInstance()->getCUDADevice());
 }
 ;
}
#endif // __CUDACC__



# 290 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"

#ifndef __CUDACC__
extern "C" CDLT_API  void _limit_vertices_by_all_neighbours(hmpprt::s32 N_21, hmpprt::s32 N3, double beta, double* centroid_values, double* vertex_values, double* edge_values_2, hmpprt::s64* neighbours_21, double* x_gradient_21, double* y_gradient_21)
{
 # 1 "<preprocessor>"
 (_limit_vertices_by_all_neighbours_internal_1(N_21, N3, beta, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (centroid_values), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (vertex_values), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (edge_values_2), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s64> (neighbours_21), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (x_gradient_21), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (y_gradient_21)));
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
      _limit_vertices_by_all_neighbours_loop1D_1 = new hmpprt::CUDAGrid(hmpprt_module, "_limit_vertices_by_all_neighbours_loop1D_1");

    }
    hmpprt::Context::getInstance()->getGrouplet()->setTarget(hmpprt::CUDA);
    hmpprt::Context::getInstance()->getGrouplet()->addSignature("fmax", "prototype fmax(x: double, y: double) : double");
    hmpprt::Context::getInstance()->getGrouplet()->addSignature("fmin", "prototype fmin(x: double, y: double) : double");
    hmpprt::Context::getInstance()->getGrouplet()->addSignature("_limit_vertices_by_all_neighbours", "prototype _limit_vertices_by_all_neighbours(N: s32, N3: s32, beta: double, centroid_values: ^cudaglob double, vertex_values: ^cudaglob double, edge_values: ^cudaglob double, neighbours: ^cudaglob s64, x_gradient: ^cudaglob double, y_gradient: ^cudaglob double)");

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
      delete _limit_vertices_by_all_neighbours_loop1D_1;

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

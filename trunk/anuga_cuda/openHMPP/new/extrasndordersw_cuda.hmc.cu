
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

# 117 "extrapolate_second_order_sw.c"

#ifndef __CUDACC__
extern "C" CDLT_API  void extrapolate_second_order_sw(hmpprt::s32 N_41, hmpprt::s32 N2, hmpprt::s32 N3_31, hmpprt::s32 N6, double epsilon, double minimum_allowed_height, double beta_w_31, double beta_w_dry, double beta_uh, double beta_uh_dry_21, double beta_vh, double beta_vh_dry, hmpprt::s32 optimise_dry_cells, hmpprt::s32 extrapolate_velocity_second_order, hmpprt::s64* surrogate_neighbours, hmpprt::s64* number_of_boundaries, double* centroid_coordinates, double* stage_centroid_values, double* bed_centroid_values, double* xmom_centroid_values_41, double* ymom_centroid_values_4, double* vertex_coordinates, double* stage_vertex_values, double* bed_vertex_values, double* xmom_vertex_values, double* ymom_vertex_values, double* stage_centroid_store, double* xmom_centroid_store, double* ymom_centroid_store_31)
;
#endif // __CUDACC__



# 117 "extrapolate_second_order_sw.c"

#ifndef __CUDACC__
void extrapolate_second_order_sw_internal_1(hmpprt::s32 N_4, hmpprt::s32 N2, hmpprt::s32 N3_3, hmpprt::s32 N6, double epsilon, double minimum_allowed_height, double beta_w_3, double beta_w_dry, double beta_uh, double beta_uh_dry_2, double beta_vh, double beta_vh_dry, hmpprt::s32 optimise_dry_cells, hmpprt::s32 extrapolate_velocity_second_order, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s64>  surrogate_neighbours, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s64>  number_of_boundaries, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  centroid_coordinates, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  stage_centroid_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  bed_centroid_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  xmom_centroid_values_4, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  ymom_centroid_values_41, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  vertex_coordinates, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  stage_vertex_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  bed_vertex_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  xmom_vertex_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  ymom_vertex_values, double* stage_centroid_store, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  xmom_centroid_store, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  ymom_centroid_store_3)
;
#endif // __CUDACC__



# 71 "extrapolate_second_order_sw.c"

#ifndef __CUDACC__
static hmpprt::CUDAGrid * extrapolate_second_order_sw_loop1D_3 = 0;
#else

extern "C" __global__ void extrapolate_second_order_sw_loop1D_3(hmpprt::s32 N, double* xmom_centroid_values_2, double* ymom_centroid_values_2, double* stage_vertex_values_3, double* bed_vertex_values_2, double* xmom_vertex_values_3, double* ymom_vertex_values_3, double* xmom_centroid_store_2, double* ymom_centroid_store_2);
#endif // __CUDACC__




# 71 "extrapolate_second_order_sw.c"

#ifndef __CUDACC__
static hmpprt::CUDAGrid * extrapolate_second_order_sw_loop1D_2 = 0;
#else

extern "C" __global__ void extrapolate_second_order_sw_loop1D_2(hmpprt::s32 N_3, double epsilon_2, double beta_w_4, double beta_w_dry_2, double beta_uh_2, double beta_uh_dry, double beta_vh_2, double beta_vh_dry_2, hmpprt::s32 optimise_dry_cells_2, hmpprt::s64* surrogate_neighbours_2, hmpprt::s64* number_of_boundaries_2, double* centroid_coordinates_2, double* stage_centroid_values_2, double* bed_centroid_values_2, double* xmom_centroid_values_3, double* ymom_centroid_values, double* vertex_coordinates_2, double* stage_vertex_values_2, double* xmom_vertex_values_2, double* ymom_vertex_values_2);
#endif // __CUDACC__




# 71 "extrapolate_second_order_sw.c"

#ifdef __CUDACC__
__device__ hmpprt::s32 limit_gradient_old_cudalocal_1(double* dqv, double qmin, double qmax, double beta_w)
;
#endif // __CUDACC__



# 11 "extrapolate_second_order_sw.c"

#ifdef __CUDACC__
__device__ hmpprt::s32 find_qmin_and_qmax_cudalocal_1(double dq0, double dq1, double dq2, double* qmin_42, double* qmax_41)
;
#endif // __CUDACC__



# 71 "extrapolate_second_order_sw.c"

#ifndef __CUDACC__
static hmpprt::CUDAGrid * extrapolate_second_order_sw_loop1D_1 = 0;
#else

extern "C" __global__ void extrapolate_second_order_sw_loop1D_1(hmpprt::s32 N_2, double minimum_allowed_height_2, double* stage_centroid_values_3, double* bed_centroid_values_3, double* xmom_centroid_values, double* ymom_centroid_values_3, double* xmom_centroid_store_3, double* ymom_centroid_store);
#endif // __CUDACC__




# 71 "extrapolate_second_order_sw.c"

#ifndef __CUDACC__
extern "C" CDLT_API  hmpprt::s32 limit_gradient_old(double* dqv, double qmin, double qmax, double beta_w)
;
#endif // __CUDACC__



# 71 "extrapolate_second_order_sw.c"

#ifndef __CUDACC__
hmpprt::s32 limit_gradient_old_internal_1(double* dqv, double qmin, double qmax, double beta_w)
;
#endif // __CUDACC__



# 11 "extrapolate_second_order_sw.c"

#ifndef __CUDACC__
extern "C" CDLT_API  hmpprt::s32 find_qmin_and_qmax(double dq0, double dq1, double dq2, double* qmin_41, double* qmax_42)
;
#endif // __CUDACC__



# 11 "extrapolate_second_order_sw.c"

#ifndef __CUDACC__
hmpprt::s32 find_qmin_and_qmax_internal_1(double dq0, double dq1, double dq2, double* qmin_4, double* qmax_4)
;
#endif // __CUDACC__



# 11 "extrapolate_second_order_sw.c"

#ifndef __CUDACC__
hmpprt::s32 find_qmin_and_qmax_internal_1(double dq0, double dq1, double dq2, double* qmin_4, double* qmax_4)
{
 # 29 "extrapolate_second_order_sw.c"
 if (dq0 >= (double) 0.0)
 {
  # 30 "extrapolate_second_order_sw.c"
  if (dq1 >= dq2)
  {
   # 31 "extrapolate_second_order_sw.c"
   if (dq1 >= (double) 0.0)
   {
    # 32 "extrapolate_second_order_sw.c"
    *qmax_4 = dq0 + dq1;
   }
   else
   {
    # 34 "extrapolate_second_order_sw.c"
    *qmax_4 = dq0;
   }
   # 36 "extrapolate_second_order_sw.c"
   *qmin_4 = dq0 + dq2;
   # 37 "extrapolate_second_order_sw.c"
   if (*qmin_4 >= (double) 0.0)
   {
    # 37 "extrapolate_second_order_sw.c"
    *qmin_4 = (double) 0.0;
   }
  }
  else
  {
   # 39 "extrapolate_second_order_sw.c"
   if (dq2 > (double) 0.0)
   {
    # 40 "extrapolate_second_order_sw.c"
    *qmax_4 = dq0 + dq2;
   }
   else
   {
    # 42 "extrapolate_second_order_sw.c"
    *qmax_4 = dq0;
   }
   # 44 "extrapolate_second_order_sw.c"
   *qmin_4 = dq0 + dq1;
   # 45 "extrapolate_second_order_sw.c"
   if (*qmin_4 >= (double) 0.0)
   {
    # 45 "extrapolate_second_order_sw.c"
    *qmin_4 = (double) 0.0;
   }
  }
 }
 else
 {
  # 48 "extrapolate_second_order_sw.c"
  if (dq1 <= dq2)
  {
   # 49 "extrapolate_second_order_sw.c"
   if (dq1 < (double) 0.0)
   {
    # 50 "extrapolate_second_order_sw.c"
    *qmin_4 = dq0 + dq1;
   }
   else
   {
    # 52 "extrapolate_second_order_sw.c"
    *qmin_4 = dq0;
   }
   # 54 "extrapolate_second_order_sw.c"
   *qmax_4 = dq0 + dq2;
   # 55 "extrapolate_second_order_sw.c"
   if (*qmax_4 <= (double) 0.0)
   {
    # 55 "extrapolate_second_order_sw.c"
    *qmax_4 = (double) 0.0;
   }
  }
  else
  {
   # 57 "extrapolate_second_order_sw.c"
   if (dq2 < (double) 0.0)
   {
    # 58 "extrapolate_second_order_sw.c"
    *qmin_4 = dq0 + dq2;
   }
   else
   {
    # 60 "extrapolate_second_order_sw.c"
    *qmin_4 = dq0;
   }
   # 62 "extrapolate_second_order_sw.c"
   *qmax_4 = dq0 + dq1;
   # 63 "extrapolate_second_order_sw.c"
   if (*qmax_4 <= (double) 0.0)
   {
    # 63 "extrapolate_second_order_sw.c"
    *qmax_4 = (double) 0.0;
   }
  }
 }
 # 11 "extrapolate_second_order_sw.c"
 return 0;
}
#endif // __CUDACC__



# 11 "extrapolate_second_order_sw.c"

#ifndef __CUDACC__
extern "C" CDLT_API  hmpprt::s32 find_qmin_and_qmax(double dq0, double dq1, double dq2, double* qmin_41, double* qmax_42)
{
 # 71 "extrapolate_second_order_sw.c"
 (find_qmin_and_qmax_internal_1(dq0, dq1, dq2, qmin_41, qmax_42));
}
#endif // __CUDACC__



# 71 "extrapolate_second_order_sw.c"

#ifndef __CUDACC__
hmpprt::s32 limit_gradient_old_internal_1(double* dqv, double qmin, double qmax, double beta_w)
{
 # 85 "extrapolate_second_order_sw.c"
 double r;
 # 85 "extrapolate_second_order_sw.c"
 double r0;
 # 85 "extrapolate_second_order_sw.c"
 double phi;
 # 85 "extrapolate_second_order_sw.c"
 r = (double) 1000.0;
 # 85 "extrapolate_second_order_sw.c"
 r0 = (double) 1.0;
 # 93 "extrapolate_second_order_sw.c"
 hmpprt::s32 i_2;
 # 93 "extrapolate_second_order_sw.c"
 # 93 "extrapolate_second_order_sw.c"
 for (i_2 = 0 ; i_2 <= 2 ; i_2 = i_2 + 1)
 {
  # 94 "extrapolate_second_order_sw.c"
  if (*(dqv + i_2) <  - (double) 1.00000000000000002e-100)
  {
   # 95 "extrapolate_second_order_sw.c"
   r0 = qmin / *(dqv + i_2);
  }
  # 97 "extrapolate_second_order_sw.c"
  if (*(dqv + i_2) > (double) 1.00000000000000002e-100)
  {
   # 98 "extrapolate_second_order_sw.c"
   r0 = qmax / *(dqv + i_2);
  }
  # 100 "extrapolate_second_order_sw.c"
  r = fmin(r0, r);
 }
 # 103 "extrapolate_second_order_sw.c"
 # 103 "extrapolate_second_order_sw.c"
 phi = fmin(r * beta_w, (double) 1.0);
 # 105 "extrapolate_second_order_sw.c"
 *dqv = *dqv * phi;
 # 106 "extrapolate_second_order_sw.c"
 *(dqv + 1) = *(dqv + 1) * phi;
 # 107 "extrapolate_second_order_sw.c"
 *(dqv + 2) = *(dqv + 2) * phi;
 # 71 "extrapolate_second_order_sw.c"
 return 0;
}
#endif // __CUDACC__



# 71 "extrapolate_second_order_sw.c"

#ifndef __CUDACC__
extern "C" CDLT_API  hmpprt::s32 limit_gradient_old(double* dqv, double qmin, double qmax, double beta_w)
{
 # 11 "extrapolate_second_order_sw.c"
 (limit_gradient_old_internal_1(dqv, qmin, qmax, beta_w));
}
#endif // __CUDACC__



# 11 "extrapolate_second_order_sw.c"

#ifdef __CUDACC__

extern "C" __global__ void extrapolate_second_order_sw_loop1D_1(hmpprt::s32 N_2, double minimum_allowed_height_2, double* stage_centroid_values_3, double* bed_centroid_values_3, double* xmom_centroid_values, double* ymom_centroid_values_3, double* xmom_centroid_store_3, double* ymom_centroid_store)
{
 # 155 "extrapolate_second_order_sw.c"
 double dk_1;
 # 155 "extrapolate_second_order_sw.c"
 hmpprt::s32 k_1;
 # 171 "extrapolate_second_order_sw.c"
 k_1 = (hmpprt::gr_atidf());
 # 171 "extrapolate_second_order_sw.c"
 if (k_1 > N_2 - 1)
 {
  # 171 "extrapolate_second_order_sw.c"
  goto __hmppcg_label_1;
 }
 # 171 "extrapolate_second_order_sw.c"
 dk_1 = fmax(*(stage_centroid_values_3 + k_1) - *(bed_centroid_values_3 + k_1), minimum_allowed_height_2);
 # 172 "extrapolate_second_order_sw.c"
 *(xmom_centroid_store_3 + k_1) = *(xmom_centroid_values + k_1);
 # 173 "extrapolate_second_order_sw.c"
 *(xmom_centroid_values + k_1) = *(xmom_centroid_values + k_1) / dk_1;
 # 175 "extrapolate_second_order_sw.c"
 *(ymom_centroid_store + k_1) = *(ymom_centroid_values_3 + k_1);
 # 176 "extrapolate_second_order_sw.c"
 *(ymom_centroid_values_3 + k_1) = *(ymom_centroid_values_3 + k_1) / dk_1;
 # 11 "extrapolate_second_order_sw.c"
 __hmppcg_label_1:;
}
#endif // __CUDACC__



# 11 "extrapolate_second_order_sw.c"

#ifdef __CUDACC__
__device__ hmpprt::s32 find_qmin_and_qmax_cudalocal_1(double dq0, double dq1, double dq2, double* qmin_42, double* qmax_41)
{
 # 29 "extrapolate_second_order_sw.c"
 if (dq0 >= (double) 0.0)
 {
  # 30 "extrapolate_second_order_sw.c"
  if (dq1 >= dq2)
  {
   # 31 "extrapolate_second_order_sw.c"
   if (dq1 >= (double) 0.0)
   {
    # 32 "extrapolate_second_order_sw.c"
    *qmax_41 = dq0 + dq1;
   }
   else
   {
    # 34 "extrapolate_second_order_sw.c"
    *qmax_41 = dq0;
   }
   # 36 "extrapolate_second_order_sw.c"
   *qmin_42 = dq0 + dq2;
   # 37 "extrapolate_second_order_sw.c"
   if (*qmin_42 >= (double) 0.0)
   {
    # 37 "extrapolate_second_order_sw.c"
    *qmin_42 = (double) 0.0;
   }
  }
  else
  {
   # 39 "extrapolate_second_order_sw.c"
   if (dq2 > (double) 0.0)
   {
    # 40 "extrapolate_second_order_sw.c"
    *qmax_41 = dq0 + dq2;
   }
   else
   {
    # 42 "extrapolate_second_order_sw.c"
    *qmax_41 = dq0;
   }
   # 44 "extrapolate_second_order_sw.c"
   *qmin_42 = dq0 + dq1;
   # 45 "extrapolate_second_order_sw.c"
   if (*qmin_42 >= (double) 0.0)
   {
    # 45 "extrapolate_second_order_sw.c"
    *qmin_42 = (double) 0.0;
   }
  }
 }
 else
 {
  # 48 "extrapolate_second_order_sw.c"
  if (dq1 <= dq2)
  {
   # 49 "extrapolate_second_order_sw.c"
   if (dq1 < (double) 0.0)
   {
    # 50 "extrapolate_second_order_sw.c"
    *qmin_42 = dq0 + dq1;
   }
   else
   {
    # 52 "extrapolate_second_order_sw.c"
    *qmin_42 = dq0;
   }
   # 54 "extrapolate_second_order_sw.c"
   *qmax_41 = dq0 + dq2;
   # 55 "extrapolate_second_order_sw.c"
   if (*qmax_41 <= (double) 0.0)
   {
    # 55 "extrapolate_second_order_sw.c"
    *qmax_41 = (double) 0.0;
   }
  }
  else
  {
   # 57 "extrapolate_second_order_sw.c"
   if (dq2 < (double) 0.0)
   {
    # 58 "extrapolate_second_order_sw.c"
    *qmin_42 = dq0 + dq2;
   }
   else
   {
    # 60 "extrapolate_second_order_sw.c"
    *qmin_42 = dq0;
   }
   # 62 "extrapolate_second_order_sw.c"
   *qmax_41 = dq0 + dq1;
   # 63 "extrapolate_second_order_sw.c"
   if (*qmax_41 <= (double) 0.0)
   {
    # 63 "extrapolate_second_order_sw.c"
    *qmax_41 = (double) 0.0;
   }
  }
 }
 # 71 "extrapolate_second_order_sw.c"
 return 0;
}
#endif // __CUDACC__



# 71 "extrapolate_second_order_sw.c"

#ifdef __CUDACC__
__device__ hmpprt::s32 limit_gradient_old_cudalocal_1(double* dqv, double qmin, double qmax, double beta_w)
{
 # 85 "extrapolate_second_order_sw.c"
 double r;
 # 85 "extrapolate_second_order_sw.c"
 double r0;
 # 85 "extrapolate_second_order_sw.c"
 double phi;
 # 85 "extrapolate_second_order_sw.c"
 r = (double) 1000.0;
 # 85 "extrapolate_second_order_sw.c"
 r0 = (double) 1.0;
 # 93 "extrapolate_second_order_sw.c"
 hmpprt::s32 i_21;
 # 93 "extrapolate_second_order_sw.c"
 # 93 "extrapolate_second_order_sw.c"
 for (i_21 = 0 ; i_21 <= 2 ; i_21 = i_21 + 1)
 {
  # 94 "extrapolate_second_order_sw.c"
  if (*(dqv + i_21) <  - (double) 1.00000000000000002e-100)
  {
   # 95 "extrapolate_second_order_sw.c"
   r0 = qmin / *(dqv + i_21);
  }
  # 97 "extrapolate_second_order_sw.c"
  if (*(dqv + i_21) > (double) 1.00000000000000002e-100)
  {
   # 98 "extrapolate_second_order_sw.c"
   r0 = qmax / *(dqv + i_21);
  }
  # 100 "extrapolate_second_order_sw.c"
  r = fmin(r0, r);
 }
 # 103 "extrapolate_second_order_sw.c"
 # 103 "extrapolate_second_order_sw.c"
 phi = fmin(r * beta_w, (double) 1.0);
 # 105 "extrapolate_second_order_sw.c"
 *dqv = *dqv * phi;
 # 106 "extrapolate_second_order_sw.c"
 *(dqv + 1) = *(dqv + 1) * phi;
 # 107 "extrapolate_second_order_sw.c"
 *(dqv + 2) = *(dqv + 2) * phi;
 # 117 "extrapolate_second_order_sw.c"
 return 0;
}
#endif // __CUDACC__



# 117 "extrapolate_second_order_sw.c"

#ifdef __CUDACC__

extern "C" __global__ void extrapolate_second_order_sw_loop1D_2(hmpprt::s32 N_3, double epsilon_2, double beta_w_4, double beta_w_dry_2, double beta_uh_2, double beta_uh_dry, double beta_vh_2, double beta_vh_dry_2, hmpprt::s32 optimise_dry_cells_2, hmpprt::s64* surrogate_neighbours_2, hmpprt::s64* number_of_boundaries_2, double* centroid_coordinates_2, double* stage_centroid_values_2, double* bed_centroid_values_2, double* xmom_centroid_values_3, double* ymom_centroid_values, double* vertex_coordinates_2, double* stage_vertex_values_2, double* xmom_vertex_values_2, double* ymom_vertex_values_2)
{
 # 158 "extrapolate_second_order_sw.c"
 double dqv_3[3uLL];
 # 158 "extrapolate_second_order_sw.c"
 double qmin_5;
 # 158 "extrapolate_second_order_sw.c"
 double qmax_5;
 # 155 "extrapolate_second_order_sw.c"
 hmpprt::s32 k3_1;
 # 155 "extrapolate_second_order_sw.c"
 hmpprt::s32 k6_1;
 # 155 "extrapolate_second_order_sw.c"
 double x_1;
 # 155 "extrapolate_second_order_sw.c"
 double y_1;
 # 155 "extrapolate_second_order_sw.c"
 double dxv0_1;
 # 155 "extrapolate_second_order_sw.c"
 double dyv0_1;
 # 155 "extrapolate_second_order_sw.c"
 double dxv1_1;
 # 155 "extrapolate_second_order_sw.c"
 double dyv1_1;
 # 155 "extrapolate_second_order_sw.c"
 double dxv2_1;
 # 155 "extrapolate_second_order_sw.c"
 double dyv2_1;
 # 155 "extrapolate_second_order_sw.c"
 hmpprt::s32 k_2;
 # 194 "extrapolate_second_order_sw.c"
 k_2 = (hmpprt::gr_atidf());
 # 194 "extrapolate_second_order_sw.c"
 if (k_2 > N_3 - 1)
 {
  # 194 "extrapolate_second_order_sw.c"
  goto __hmppcg_label_2;
 }
 # 194 "extrapolate_second_order_sw.c"
 k3_1 = k_2 * 3;
 # 195 "extrapolate_second_order_sw.c"
 k6_1 = k_2 * 6;
 # 197 "extrapolate_second_order_sw.c"
 if (*(number_of_boundaries_2 + k_2) == 3LL)
 {
  # 200 "extrapolate_second_order_sw.c"
  *(stage_vertex_values_2 + k3_1) = *(stage_centroid_values_2 + k_2);
  # 201 "extrapolate_second_order_sw.c"
  *(stage_vertex_values_2 + (k3_1 + 1)) = *(stage_centroid_values_2 + k_2);
  # 202 "extrapolate_second_order_sw.c"
  *(stage_vertex_values_2 + (k3_1 + 2)) = *(stage_centroid_values_2 + k_2);
  # 203 "extrapolate_second_order_sw.c"
  *(xmom_vertex_values_2 + k3_1) = *(xmom_centroid_values_3 + k_2);
  # 204 "extrapolate_second_order_sw.c"
  *(xmom_vertex_values_2 + (k3_1 + 1)) = *(xmom_centroid_values_3 + k_2);
  # 205 "extrapolate_second_order_sw.c"
  *(xmom_vertex_values_2 + (k3_1 + 2)) = *(xmom_centroid_values_3 + k_2);
  # 206 "extrapolate_second_order_sw.c"
  *(ymom_vertex_values_2 + k3_1) = *(ymom_centroid_values + k_2);
  # 207 "extrapolate_second_order_sw.c"
  *(ymom_vertex_values_2 + (k3_1 + 1)) = *(ymom_centroid_values + k_2);
  # 208 "extrapolate_second_order_sw.c"
  *(ymom_vertex_values_2 + (k3_1 + 2)) = *(ymom_centroid_values + k_2);
  # 211 "extrapolate_second_order_sw.c"
  goto __hmppcg_label_2;
 }
 else
 {
  # 216 "extrapolate_second_order_sw.c"
  hmpprt::s32 coord_index_1;
  # 216 "extrapolate_second_order_sw.c"
  double xv0_1;
  # 216 "extrapolate_second_order_sw.c"
  double xv1_1;
  # 216 "extrapolate_second_order_sw.c"
  double xv2_1;
  # 216 "extrapolate_second_order_sw.c"
  double yv0_1;
  # 216 "extrapolate_second_order_sw.c"
  double yv1_1;
  # 216 "extrapolate_second_order_sw.c"
  double yv2_1;
  # 216 "extrapolate_second_order_sw.c"
  xv0_1 = *(vertex_coordinates_2 + k6_1);
  # 217 "extrapolate_second_order_sw.c"
  yv0_1 = *(vertex_coordinates_2 + (k6_1 + 1));
  # 218 "extrapolate_second_order_sw.c"
  xv1_1 = *(vertex_coordinates_2 + (k6_1 + 2));
  # 219 "extrapolate_second_order_sw.c"
  yv1_1 = *(vertex_coordinates_2 + (k6_1 + 3));
  # 220 "extrapolate_second_order_sw.c"
  xv2_1 = *(vertex_coordinates_2 + (k6_1 + 4));
  # 221 "extrapolate_second_order_sw.c"
  yv2_1 = *(vertex_coordinates_2 + (k6_1 + 5));
  # 224 "extrapolate_second_order_sw.c"
  coord_index_1 = 2 * k_2;
  # 225 "extrapolate_second_order_sw.c"
  x_1 = *(centroid_coordinates_2 + coord_index_1);
  # 226 "extrapolate_second_order_sw.c"
  y_1 = *(centroid_coordinates_2 + (coord_index_1 + 1));
  # 230 "extrapolate_second_order_sw.c"
  dxv0_1 = xv0_1 - x_1;
  # 231 "extrapolate_second_order_sw.c"
  dxv1_1 = xv1_1 - x_1;
  # 232 "extrapolate_second_order_sw.c"
  dxv2_1 = xv2_1 - x_1;
  # 233 "extrapolate_second_order_sw.c"
  dyv0_1 = yv0_1 - y_1;
  # 234 "extrapolate_second_order_sw.c"
  dyv1_1 = yv1_1 - y_1;
  # 235 "extrapolate_second_order_sw.c"
  dyv2_1 = yv2_1 - y_1;
 }
 # 241 "extrapolate_second_order_sw.c"
 if (*(number_of_boundaries_2 + k_2) <= 1LL)
 {
  # 252 "extrapolate_second_order_sw.c"
  hmpprt::s32 k0_1;
  # 252 "extrapolate_second_order_sw.c"
  hmpprt::s32 coord_index_2;
  # 252 "extrapolate_second_order_sw.c"
  hmpprt::s32 k1_1;
  # 252 "extrapolate_second_order_sw.c"
  hmpprt::s32 coord_index_3;
  # 252 "extrapolate_second_order_sw.c"
  hmpprt::s32 k2_1;
  # 252 "extrapolate_second_order_sw.c"
  hmpprt::s32 coord_index_4;
  # 252 "extrapolate_second_order_sw.c"
  double x1_1;
  # 252 "extrapolate_second_order_sw.c"
  double x0_1;
  # 252 "extrapolate_second_order_sw.c"
  double x2_1;
  # 252 "extrapolate_second_order_sw.c"
  double y1_1;
  # 252 "extrapolate_second_order_sw.c"
  double y0_1;
  # 252 "extrapolate_second_order_sw.c"
  double y2_1;
  # 252 "extrapolate_second_order_sw.c"
  double dy2_1;
  # 252 "extrapolate_second_order_sw.c"
  double dx1_1;
  # 252 "extrapolate_second_order_sw.c"
  double dy1_1;
  # 252 "extrapolate_second_order_sw.c"
  double dx2_1;
  # 252 "extrapolate_second_order_sw.c"
  double area2_1;
  # 252 "extrapolate_second_order_sw.c"
  double h0_1;
  # 252 "extrapolate_second_order_sw.c"
  double h1_1;
  # 252 "extrapolate_second_order_sw.c"
  double h2_1;
  # 252 "extrapolate_second_order_sw.c"
  double hc_1;
  # 252 "extrapolate_second_order_sw.c"
  double hmin_1;
  # 252 "extrapolate_second_order_sw.c"
  double dq1_3;
  # 252 "extrapolate_second_order_sw.c"
  double dq2_3;
  # 252 "extrapolate_second_order_sw.c"
  double a_1;
  # 252 "extrapolate_second_order_sw.c"
  double inv_area2_1;
  # 252 "extrapolate_second_order_sw.c"
  double b_1;
  # 252 "extrapolate_second_order_sw.c"
  double a_2;
  # 252 "extrapolate_second_order_sw.c"
  double b_2;
  # 252 "extrapolate_second_order_sw.c"
  double dq0_3;
  # 252 "extrapolate_second_order_sw.c"
  double hfactor_1;
  # 252 "extrapolate_second_order_sw.c"
  double beta_tmp_1;
  # 252 "extrapolate_second_order_sw.c"
  double dq1_4;
  # 252 "extrapolate_second_order_sw.c"
  double dq2_4;
  # 252 "extrapolate_second_order_sw.c"
  double a_3;
  # 252 "extrapolate_second_order_sw.c"
  double b_3;
  # 252 "extrapolate_second_order_sw.c"
  double a_4;
  # 252 "extrapolate_second_order_sw.c"
  double b_4;
  # 252 "extrapolate_second_order_sw.c"
  double dq0_4;
  # 252 "extrapolate_second_order_sw.c"
  double beta_tmp_2;
  # 252 "extrapolate_second_order_sw.c"
  double dq1_5;
  # 252 "extrapolate_second_order_sw.c"
  double dq2_5;
  # 252 "extrapolate_second_order_sw.c"
  double a_5;
  # 252 "extrapolate_second_order_sw.c"
  double b_5;
  # 252 "extrapolate_second_order_sw.c"
  double a_6;
  # 252 "extrapolate_second_order_sw.c"
  double b_6;
  # 252 "extrapolate_second_order_sw.c"
  double dq0_5;
  # 252 "extrapolate_second_order_sw.c"
  double beta_tmp_3;
  # 252 "extrapolate_second_order_sw.c"
  k0_1 = (hmpprt::s32 ) (*(surrogate_neighbours_2 + k3_1));
  # 253 "extrapolate_second_order_sw.c"
  k1_1 = (hmpprt::s32 ) (*(surrogate_neighbours_2 + (k3_1 + 1)));
  # 254 "extrapolate_second_order_sw.c"
  k2_1 = (hmpprt::s32 ) (*(surrogate_neighbours_2 + (k3_1 + 2)));
  # 258 "extrapolate_second_order_sw.c"
  coord_index_2 = 2 * k0_1;
  # 259 "extrapolate_second_order_sw.c"
  x0_1 = *(centroid_coordinates_2 + coord_index_2);
  # 260 "extrapolate_second_order_sw.c"
  y0_1 = *(centroid_coordinates_2 + (coord_index_2 + 1));
  # 262 "extrapolate_second_order_sw.c"
  coord_index_3 = 2 * k1_1;
  # 263 "extrapolate_second_order_sw.c"
  x1_1 = *(centroid_coordinates_2 + coord_index_3);
  # 264 "extrapolate_second_order_sw.c"
  y1_1 = *(centroid_coordinates_2 + (coord_index_3 + 1));
  # 266 "extrapolate_second_order_sw.c"
  coord_index_4 = 2 * k2_1;
  # 267 "extrapolate_second_order_sw.c"
  x2_1 = *(centroid_coordinates_2 + coord_index_4);
  # 268 "extrapolate_second_order_sw.c"
  y2_1 = *(centroid_coordinates_2 + (coord_index_4 + 1));
  # 272 "extrapolate_second_order_sw.c"
  dx1_1 = x1_1 - x0_1;
  # 273 "extrapolate_second_order_sw.c"
  dx2_1 = x2_1 - x0_1;
  # 274 "extrapolate_second_order_sw.c"
  dy1_1 = y1_1 - y0_1;
  # 275 "extrapolate_second_order_sw.c"
  dy2_1 = y2_1 - y0_1;
  # 279 "extrapolate_second_order_sw.c"
  area2_1 = dy2_1 * dx1_1 - dy1_1 * dx2_1;
  # 284 "extrapolate_second_order_sw.c"
  if (area2_1 <= (double) 0.0)
  {
   # 289 "extrapolate_second_order_sw.c"
   *(stage_vertex_values_2 + k3_1) = *(stage_centroid_values_2 + k_2);
   # 290 "extrapolate_second_order_sw.c"
   *(stage_vertex_values_2 + (k3_1 + 1)) = *(stage_centroid_values_2 + k_2);
   # 291 "extrapolate_second_order_sw.c"
   *(stage_vertex_values_2 + (k3_1 + 2)) = *(stage_centroid_values_2 + k_2);
   # 292 "extrapolate_second_order_sw.c"
   *(xmom_vertex_values_2 + k3_1) = *(xmom_centroid_values_3 + k_2);
   # 293 "extrapolate_second_order_sw.c"
   *(xmom_vertex_values_2 + (k3_1 + 1)) = *(xmom_centroid_values_3 + k_2);
   # 294 "extrapolate_second_order_sw.c"
   *(xmom_vertex_values_2 + (k3_1 + 2)) = *(xmom_centroid_values_3 + k_2);
   # 295 "extrapolate_second_order_sw.c"
   *(ymom_vertex_values_2 + k3_1) = *(ymom_centroid_values + k_2);
   # 296 "extrapolate_second_order_sw.c"
   *(ymom_vertex_values_2 + (k3_1 + 1)) = *(ymom_centroid_values + k_2);
   # 297 "extrapolate_second_order_sw.c"
   *(ymom_vertex_values_2 + (k3_1 + 2)) = *(ymom_centroid_values + k_2);
   # 303 "extrapolate_second_order_sw.c"
   goto __hmppcg_label_2;
  }
  # 303 "extrapolate_second_order_sw.c"
  hc_1 = *(stage_centroid_values_2 + k_2) - *(bed_centroid_values_2 + k_2);
  # 304 "extrapolate_second_order_sw.c"
  h0_1 = *(stage_centroid_values_2 + k0_1) - *(bed_centroid_values_2 + k0_1);
  # 305 "extrapolate_second_order_sw.c"
  h1_1 = *(stage_centroid_values_2 + k1_1) - *(bed_centroid_values_2 + k1_1);
  # 306 "extrapolate_second_order_sw.c"
  h2_1 = *(stage_centroid_values_2 + k2_1) - *(bed_centroid_values_2 + k2_1);
  # 307 "extrapolate_second_order_sw.c"
  hmin_1 = fmin(fmin(h0_1, fmin(h1_1, h2_1)), hc_1);
  # 310 "extrapolate_second_order_sw.c"
  hfactor_1 = (double) 0.0;
  # 311 "extrapolate_second_order_sw.c"
  if (hmin_1 > (double) 0.001000000000000000021)
  {
   # 312 "extrapolate_second_order_sw.c"
   hfactor_1 = (hmin_1 - (double) 0.001000000000000000021) / (hmin_1 + (double) 0.004000000000000000083);
  }
  # 315 "extrapolate_second_order_sw.c"
  if (optimise_dry_cells_2)
  {
   # 319 "extrapolate_second_order_sw.c"
   double hmax_1;
   # 319 "extrapolate_second_order_sw.c"
   hmax_1 = fmax(h0_1, fmax(h1_1, h2_1));
   # 320 "extrapolate_second_order_sw.c"
   if (hmax_1 < epsilon_2)
   {
    # 331 "extrapolate_second_order_sw.c"
    goto __hmppcg_label_2;
   }
  }
  # 331 "extrapolate_second_order_sw.c"
  dq0_3 = *(stage_centroid_values_2 + k0_1) - *(stage_centroid_values_2 + k_2);
  # 335 "extrapolate_second_order_sw.c"
  dq1_3 = *(stage_centroid_values_2 + k1_1) - *(stage_centroid_values_2 + k0_1);
  # 336 "extrapolate_second_order_sw.c"
  dq2_3 = *(stage_centroid_values_2 + k2_1) - *(stage_centroid_values_2 + k0_1);
  # 338 "extrapolate_second_order_sw.c"
  inv_area2_1 = (double) 1.0 / area2_1;
  # 340 "extrapolate_second_order_sw.c"
  a_1 = dy2_1 * dq1_3 - dy1_1 * dq2_3;
  # 341 "extrapolate_second_order_sw.c"
  a_2 = a_1 * inv_area2_1;
  # 342 "extrapolate_second_order_sw.c"
  b_1 = dx1_1 * dq2_3 - dx2_1 * dq1_3;
  # 343 "extrapolate_second_order_sw.c"
  b_2 = b_1 * inv_area2_1;
  # 347 "extrapolate_second_order_sw.c"
  dqv_3[0] = a_2 * dxv0_1 + b_2 * dyv0_1;
  # 348 "extrapolate_second_order_sw.c"
  dqv_3[1] = a_2 * dxv1_1 + b_2 * dyv1_1;
  # 349 "extrapolate_second_order_sw.c"
  dqv_3[2] = a_2 * dxv2_1 + b_2 * dyv2_1;
  # 354 "extrapolate_second_order_sw.c"
  (find_qmin_and_qmax_cudalocal_1(dq0_3, dq1_3, dq2_3, &qmin_5, &qmax_5));
  # 360 "extrapolate_second_order_sw.c"
  beta_tmp_1 = beta_w_dry_2 + (beta_w_4 - beta_w_dry_2) * hfactor_1;
  # 367 "extrapolate_second_order_sw.c"
  (limit_gradient_old_cudalocal_1(&dqv_3[0], qmin_5, qmax_5, beta_tmp_1));
  # 370 "extrapolate_second_order_sw.c"
  *(stage_vertex_values_2 + k3_1) = *(stage_centroid_values_2 + k_2) + dqv_3[0];
  # 371 "extrapolate_second_order_sw.c"
  *(stage_vertex_values_2 + (k3_1 + 1)) = *(stage_centroid_values_2 + k_2) + dqv_3[1];
  # 372 "extrapolate_second_order_sw.c"
  *(stage_vertex_values_2 + (k3_1 + 2)) = *(stage_centroid_values_2 + k_2) + dqv_3[2];
  # 381 "extrapolate_second_order_sw.c"
  dq0_4 = *(xmom_centroid_values_3 + k0_1) - *(xmom_centroid_values_3 + k_2);
  # 385 "extrapolate_second_order_sw.c"
  dq1_4 = *(xmom_centroid_values_3 + k1_1) - *(xmom_centroid_values_3 + k0_1);
  # 386 "extrapolate_second_order_sw.c"
  dq2_4 = *(xmom_centroid_values_3 + k2_1) - *(xmom_centroid_values_3 + k0_1);
  # 389 "extrapolate_second_order_sw.c"
  a_3 = dy2_1 * dq1_4 - dy1_1 * dq2_4;
  # 390 "extrapolate_second_order_sw.c"
  a_4 = a_3 * inv_area2_1;
  # 391 "extrapolate_second_order_sw.c"
  b_3 = dx1_1 * dq2_4 - dx2_1 * dq1_4;
  # 392 "extrapolate_second_order_sw.c"
  b_4 = b_3 * inv_area2_1;
  # 396 "extrapolate_second_order_sw.c"
  dqv_3[0] = a_4 * dxv0_1 + b_4 * dyv0_1;
  # 397 "extrapolate_second_order_sw.c"
  dqv_3[1] = a_4 * dxv1_1 + b_4 * dyv1_1;
  # 398 "extrapolate_second_order_sw.c"
  dqv_3[2] = a_4 * dxv2_1 + b_4 * dyv2_1;
  # 403 "extrapolate_second_order_sw.c"
  (find_qmin_and_qmax_cudalocal_1(dq0_4, dq1_4, dq2_4, &qmin_5, &qmax_5));
  # 407 "extrapolate_second_order_sw.c"
  beta_tmp_2 = beta_uh_dry + (beta_uh_2 - beta_uh_dry) * hfactor_1;
  # 410 "extrapolate_second_order_sw.c"
  (limit_gradient_old_cudalocal_1(&dqv_3[0], qmin_5, qmax_5, beta_tmp_2));
  # 412 "extrapolate_second_order_sw.c"
  hmpprt::s32 i_3;
  # 412 "extrapolate_second_order_sw.c"
  # 412 "extrapolate_second_order_sw.c"
  for (i_3 = 0 ; i_3 <= 2 ; i_3 = i_3 + 1)
  {
   # 413 "extrapolate_second_order_sw.c"
   *(xmom_vertex_values_2 + (k3_1 + i_3)) = *(xmom_centroid_values_3 + k_2) + dqv_3[i_3];
  }
  # 422 "extrapolate_second_order_sw.c"
  # 422 "extrapolate_second_order_sw.c"
  dq0_5 = *(ymom_centroid_values + k0_1) - *(ymom_centroid_values + k_2);
  # 426 "extrapolate_second_order_sw.c"
  dq1_5 = *(ymom_centroid_values + k1_1) - *(ymom_centroid_values + k0_1);
  # 427 "extrapolate_second_order_sw.c"
  dq2_5 = *(ymom_centroid_values + k2_1) - *(ymom_centroid_values + k0_1);
  # 430 "extrapolate_second_order_sw.c"
  a_5 = dy2_1 * dq1_5 - dy1_1 * dq2_5;
  # 431 "extrapolate_second_order_sw.c"
  a_6 = a_5 * inv_area2_1;
  # 432 "extrapolate_second_order_sw.c"
  b_5 = dx1_1 * dq2_5 - dx2_1 * dq1_5;
  # 433 "extrapolate_second_order_sw.c"
  b_6 = b_5 * inv_area2_1;
  # 437 "extrapolate_second_order_sw.c"
  dqv_3[0] = a_6 * dxv0_1 + b_6 * dyv0_1;
  # 438 "extrapolate_second_order_sw.c"
  dqv_3[1] = a_6 * dxv1_1 + b_6 * dyv1_1;
  # 439 "extrapolate_second_order_sw.c"
  dqv_3[2] = a_6 * dxv2_1 + b_6 * dyv2_1;
  # 444 "extrapolate_second_order_sw.c"
  (find_qmin_and_qmax_cudalocal_1(dq0_5, dq1_5, dq2_5, &qmin_5, &qmax_5));
  # 450 "extrapolate_second_order_sw.c"
  beta_tmp_3 = beta_vh_dry_2 + (beta_vh_2 - beta_vh_dry_2) * hfactor_1;
  # 453 "extrapolate_second_order_sw.c"
  (limit_gradient_old_cudalocal_1(&dqv_3[0], qmin_5, qmax_5, beta_tmp_3));
  # 455 "extrapolate_second_order_sw.c"
  hmpprt::s32 i_4;
  # 455 "extrapolate_second_order_sw.c"
  # 455 "extrapolate_second_order_sw.c"
  for (i_4 = 0 ; i_4 <= 2 ; i_4 = i_4 + 1)
  {
   # 456 "extrapolate_second_order_sw.c"
   *(ymom_vertex_values_2 + (k3_1 + i_4)) = *(ymom_centroid_values + k_2) + dqv_3[i_4];
  }
  # 459 "extrapolate_second_order_sw.c"
 }
 else
 {
  # 469 "extrapolate_second_order_sw.c"
  hmpprt::s32 k2_2;
  # 469 "extrapolate_second_order_sw.c"
  hmpprt::s32 k1_2;
  # 469 "extrapolate_second_order_sw.c"
  hmpprt::s32 coord_index_5;
  # 469 "extrapolate_second_order_sw.c"
  double x1_2;
  # 469 "extrapolate_second_order_sw.c"
  double y1_2;
  # 469 "extrapolate_second_order_sw.c"
  double dx1_2;
  # 469 "extrapolate_second_order_sw.c"
  double dy1_2;
  # 469 "extrapolate_second_order_sw.c"
  double area2_2;
  # 469 "extrapolate_second_order_sw.c"
  double dx2_2;
  # 469 "extrapolate_second_order_sw.c"
  double dq1_6;
  # 469 "extrapolate_second_order_sw.c"
  double dx2_3;
  # 469 "extrapolate_second_order_sw.c"
  double dy2_2;
  # 469 "extrapolate_second_order_sw.c"
  double a_7;
  # 469 "extrapolate_second_order_sw.c"
  double b_7;
  # 469 "extrapolate_second_order_sw.c"
  double dq1_7;
  # 469 "extrapolate_second_order_sw.c"
  double a_8;
  # 469 "extrapolate_second_order_sw.c"
  double b_8;
  # 469 "extrapolate_second_order_sw.c"
  double dq1_8;
  # 469 "extrapolate_second_order_sw.c"
  double a_9;
  # 469 "extrapolate_second_order_sw.c"
  double b_9;
  # 469 "extrapolate_second_order_sw.c"
  # 469 "extrapolate_second_order_sw.c"
  for (k2_2 = k3_1 ; k2_2 < k3_1 + 3 ; k2_2 = k2_2 + 1)
  {
   # 473 "extrapolate_second_order_sw.c"
   if (*(surrogate_neighbours_2 + k2_2) != (hmpprt::s64 ) (k_2))
   {
    # 474 "extrapolate_second_order_sw.c"
    break;
   }
  }
  # 474 "extrapolate_second_order_sw.c"
  # 484 "extrapolate_second_order_sw.c"
  k1_2 = (hmpprt::s32 ) (*(surrogate_neighbours_2 + k2_2));
  # 488 "extrapolate_second_order_sw.c"
  coord_index_5 = 2 * k1_2;
  # 489 "extrapolate_second_order_sw.c"
  x1_2 = *(centroid_coordinates_2 + coord_index_5);
  # 490 "extrapolate_second_order_sw.c"
  y1_2 = *(centroid_coordinates_2 + (coord_index_5 + 1));
  # 494 "extrapolate_second_order_sw.c"
  dx1_2 = x1_2 - x_1;
  # 495 "extrapolate_second_order_sw.c"
  dy1_2 = y1_2 - y_1;
  # 498 "extrapolate_second_order_sw.c"
  area2_2 = dx1_2 * dx1_2 + dy1_2 * dy1_2;
  # 504 "extrapolate_second_order_sw.c"
  dx2_2 = (double) 1.0 / area2_2;
  # 505 "extrapolate_second_order_sw.c"
  dy2_2 = dx2_2 * dy1_2;
  # 506 "extrapolate_second_order_sw.c"
  dx2_3 = dx2_2 * dx1_2;
  # 514 "extrapolate_second_order_sw.c"
  dq1_6 = *(stage_centroid_values_2 + k1_2) - *(stage_centroid_values_2 + k_2);
  # 518 "extrapolate_second_order_sw.c"
  a_7 = dq1_6 * dx2_3;
  # 519 "extrapolate_second_order_sw.c"
  b_7 = dq1_6 * dy2_2;
  # 522 "extrapolate_second_order_sw.c"
  dqv_3[0] = a_7 * dxv0_1 + b_7 * dyv0_1;
  # 523 "extrapolate_second_order_sw.c"
  dqv_3[1] = a_7 * dxv1_1 + b_7 * dyv1_1;
  # 524 "extrapolate_second_order_sw.c"
  dqv_3[2] = a_7 * dxv2_1 + b_7 * dyv2_1;
  # 527 "extrapolate_second_order_sw.c"
  if (dq1_6 >= (double) 0.0)
  {
   # 528 "extrapolate_second_order_sw.c"
   qmin_5 = (double) 0.0;
   # 529 "extrapolate_second_order_sw.c"
   qmax_5 = dq1_6;
  }
  else
  {
   # 531 "extrapolate_second_order_sw.c"
   qmin_5 = dq1_6;
   # 532 "extrapolate_second_order_sw.c"
   qmax_5 = (double) 0.0;
  }
  # 536 "extrapolate_second_order_sw.c"
  (limit_gradient_old_cudalocal_1(&dqv_3[0], qmin_5, qmax_5, beta_w_4));
  # 540 "extrapolate_second_order_sw.c"
  *(stage_vertex_values_2 + k3_1) = *(stage_centroid_values_2 + k_2) + dqv_3[0];
  # 541 "extrapolate_second_order_sw.c"
  *(stage_vertex_values_2 + (k3_1 + 1)) = *(stage_centroid_values_2 + k_2) + dqv_3[1];
  # 542 "extrapolate_second_order_sw.c"
  *(stage_vertex_values_2 + (k3_1 + 2)) = *(stage_centroid_values_2 + k_2) + dqv_3[2];
  # 550 "extrapolate_second_order_sw.c"
  dq1_7 = *(xmom_centroid_values_3 + k1_2) - *(xmom_centroid_values_3 + k_2);
  # 554 "extrapolate_second_order_sw.c"
  a_8 = dq1_7 * dx2_3;
  # 555 "extrapolate_second_order_sw.c"
  b_8 = dq1_7 * dy2_2;
  # 558 "extrapolate_second_order_sw.c"
  dqv_3[0] = a_8 * dxv0_1 + b_8 * dyv0_1;
  # 559 "extrapolate_second_order_sw.c"
  dqv_3[1] = a_8 * dxv1_1 + b_8 * dyv1_1;
  # 560 "extrapolate_second_order_sw.c"
  dqv_3[2] = a_8 * dxv2_1 + b_8 * dyv2_1;
  # 563 "extrapolate_second_order_sw.c"
  if (dq1_7 >= (double) 0.0)
  {
   # 564 "extrapolate_second_order_sw.c"
   qmin_5 = (double) 0.0;
   # 565 "extrapolate_second_order_sw.c"
   qmax_5 = dq1_7;
  }
  else
  {
   # 567 "extrapolate_second_order_sw.c"
   qmin_5 = dq1_7;
   # 568 "extrapolate_second_order_sw.c"
   qmax_5 = (double) 0.0;
  }
  # 572 "extrapolate_second_order_sw.c"
  (limit_gradient_old_cudalocal_1(&dqv_3[0], qmin_5, qmax_5, beta_w_4));
  # 579 "extrapolate_second_order_sw.c"
  hmpprt::s32 i_5;
  # 579 "extrapolate_second_order_sw.c"
  # 579 "extrapolate_second_order_sw.c"
  for (i_5 = 0 ; i_5 <= 2 ; i_5 = i_5 + 1)
  {
   # 580 "extrapolate_second_order_sw.c"
   *(xmom_vertex_values_2 + (k3_1 + i_5)) = *(xmom_centroid_values_3 + k_2) + dqv_3[i_5];
  }
  # 588 "extrapolate_second_order_sw.c"
  # 588 "extrapolate_second_order_sw.c"
  dq1_8 = *(ymom_centroid_values + k1_2) - *(ymom_centroid_values + k_2);
  # 592 "extrapolate_second_order_sw.c"
  a_9 = dq1_8 * dx2_3;
  # 593 "extrapolate_second_order_sw.c"
  b_9 = dq1_8 * dy2_2;
  # 596 "extrapolate_second_order_sw.c"
  dqv_3[0] = a_9 * dxv0_1 + b_9 * dyv0_1;
  # 597 "extrapolate_second_order_sw.c"
  dqv_3[1] = a_9 * dxv1_1 + b_9 * dyv1_1;
  # 598 "extrapolate_second_order_sw.c"
  dqv_3[2] = a_9 * dxv2_1 + b_9 * dyv2_1;
  # 601 "extrapolate_second_order_sw.c"
  if (dq1_8 >= (double) 0.0)
  {
   # 602 "extrapolate_second_order_sw.c"
   qmin_5 = (double) 0.0;
   # 603 "extrapolate_second_order_sw.c"
   qmax_5 = dq1_8;
  }
  else
  {
   # 606 "extrapolate_second_order_sw.c"
   qmin_5 = dq1_8;
   # 607 "extrapolate_second_order_sw.c"
   qmax_5 = (double) 0.0;
  }
  # 611 "extrapolate_second_order_sw.c"
  (limit_gradient_old_cudalocal_1(&dqv_3[0], qmin_5, qmax_5, beta_w_4));
  # 618 "extrapolate_second_order_sw.c"
  hmpprt::s32 i_6;
  # 618 "extrapolate_second_order_sw.c"
  # 618 "extrapolate_second_order_sw.c"
  for (i_6 = 0 ; i_6 <= 2 ; i_6 = i_6 + 1)
  {
   # 619 "extrapolate_second_order_sw.c"
   *(ymom_vertex_values_2 + (k3_1 + i_6)) = *(ymom_centroid_values + k_2) + dqv_3[i_6];
  }
  # 117 "extrapolate_second_order_sw.c"
 }
 # 117 "extrapolate_second_order_sw.c"
 __hmppcg_label_2:;
}
#endif // __CUDACC__



# 117 "extrapolate_second_order_sw.c"

#ifdef __CUDACC__

extern "C" __global__ void extrapolate_second_order_sw_loop1D_3(hmpprt::s32 N, double* xmom_centroid_values_2, double* ymom_centroid_values_2, double* stage_vertex_values_3, double* bed_vertex_values_2, double* xmom_vertex_values_3, double* ymom_vertex_values_3, double* xmom_centroid_store_2, double* ymom_centroid_store_2)
{
 # 155 "extrapolate_second_order_sw.c"
 hmpprt::s32 k3_2;
 # 155 "extrapolate_second_order_sw.c"
 double dv0_1;
 # 155 "extrapolate_second_order_sw.c"
 double dv1_1;
 # 155 "extrapolate_second_order_sw.c"
 double dv2_1;
 # 155 "extrapolate_second_order_sw.c"
 hmpprt::s32 k_3;
 # 640 "extrapolate_second_order_sw.c"
 k_3 = (hmpprt::gr_atidf());
 # 640 "extrapolate_second_order_sw.c"
 if (k_3 > N - 1)
 {
  # 640 "extrapolate_second_order_sw.c"
  goto __hmppcg_label_3;
 }
 # 640 "extrapolate_second_order_sw.c"
 k3_2 = 3 * k_3;
 # 644 "extrapolate_second_order_sw.c"
 dv0_1 = fmax(*(stage_vertex_values_3 + k3_2) - *(bed_vertex_values_2 + k3_2), (double) 0.0);
 # 645 "extrapolate_second_order_sw.c"
 dv1_1 = fmax(*(stage_vertex_values_3 + (k3_2 + 1)) - *(bed_vertex_values_2 + (k3_2 + 1)), (double) 0.0);
 # 646 "extrapolate_second_order_sw.c"
 dv2_1 = fmax(*(stage_vertex_values_3 + (k3_2 + 2)) - *(bed_vertex_values_2 + (k3_2 + 2)), (double) 0.0);
 # 649 "extrapolate_second_order_sw.c"
 *(xmom_centroid_values_2 + k_3) = *(xmom_centroid_store_2 + k_3);
 # 650 "extrapolate_second_order_sw.c"
 *(xmom_vertex_values_3 + k3_2) = *(xmom_vertex_values_3 + k3_2) * dv0_1;
 # 651 "extrapolate_second_order_sw.c"
 *(xmom_vertex_values_3 + (k3_2 + 1)) = *(xmom_vertex_values_3 + (k3_2 + 1)) * dv1_1;
 # 652 "extrapolate_second_order_sw.c"
 *(xmom_vertex_values_3 + (k3_2 + 2)) = *(xmom_vertex_values_3 + (k3_2 + 2)) * dv2_1;
 # 654 "extrapolate_second_order_sw.c"
 *(ymom_centroid_values_2 + k_3) = *(ymom_centroid_store_2 + k_3);
 # 655 "extrapolate_second_order_sw.c"
 *(ymom_vertex_values_3 + k3_2) = *(ymom_vertex_values_3 + k3_2) * dv0_1;
 # 656 "extrapolate_second_order_sw.c"
 *(ymom_vertex_values_3 + (k3_2 + 1)) = *(ymom_vertex_values_3 + (k3_2 + 1)) * dv1_1;
 # 657 "extrapolate_second_order_sw.c"
 *(ymom_vertex_values_3 + (k3_2 + 2)) = *(ymom_vertex_values_3 + (k3_2 + 2)) * dv2_1;
 # 117 "extrapolate_second_order_sw.c"
 __hmppcg_label_3:;
}
#endif // __CUDACC__



# 117 "extrapolate_second_order_sw.c"

#ifndef __CUDACC__
void extrapolate_second_order_sw_internal_1(hmpprt::s32 N_4, hmpprt::s32 N2, hmpprt::s32 N3_3, hmpprt::s32 N6, double epsilon, double minimum_allowed_height, double beta_w_3, double beta_w_dry, double beta_uh, double beta_uh_dry_2, double beta_vh, double beta_vh_dry, hmpprt::s32 optimise_dry_cells, hmpprt::s32 extrapolate_velocity_second_order, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s64>  surrogate_neighbours, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s64>  number_of_boundaries, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  centroid_coordinates, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  stage_centroid_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  bed_centroid_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  xmom_centroid_values_4, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  ymom_centroid_values_41, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  vertex_coordinates, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  stage_vertex_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  bed_vertex_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  xmom_vertex_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  ymom_vertex_values, double* stage_centroid_store, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  xmom_centroid_store, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  ymom_centroid_store_3)
{
 # 164 "extrapolate_second_order_sw.c"
 if (extrapolate_velocity_second_order == 1)
 {
  # 631 "extrapolate_second_order_sw.c"
  if (N_4 - 1 >= 0)
  {
   hmpprt::CUDAGridCall __hmppcg_call;
   __hmppcg_call.setSizeX((N_4 - 1) / 128 + 1);
   __hmppcg_call.setSizeY(1);
   __hmppcg_call.setBlockSizeX(32);
   __hmppcg_call.setBlockSizeY(4);
   __hmppcg_call.addLocalParameter((hmpprt::s32) (N_4), "N_2");
   __hmppcg_call.addLocalParameter(&minimum_allowed_height, 8, "minimum_allowed_height_2");
   __hmppcg_call.addLocalParameter(&stage_centroid_values, 8, "stage_centroid_values_3");
   __hmppcg_call.addLocalParameter(&bed_centroid_values, 8, "bed_centroid_values_3");
   __hmppcg_call.addLocalParameter(&xmom_centroid_values_4, 8, "xmom_centroid_values");
   __hmppcg_call.addLocalParameter(&ymom_centroid_values_41, 8, "ymom_centroid_values_3");
   __hmppcg_call.addLocalParameter(&xmom_centroid_store, 8, "xmom_centroid_store_3");
   __hmppcg_call.addLocalParameter(&ymom_centroid_store_3, 8, "ymom_centroid_store");
   __hmppcg_call.launch(extrapolate_second_order_sw_loop1D_1, hmpprt::Context::getInstance()->getCUDADevice());
  }
  ;
 }
 # 631 "extrapolate_second_order_sw.c"
 if (N_4 - 1 >= 0)
 {
  hmpprt::CUDAGridCall __hmppcg_call;
  __hmppcg_call.setSizeX((N_4 - 1) / 128 + 1);
  __hmppcg_call.setSizeY(1);
  __hmppcg_call.setBlockSizeX(32);
  __hmppcg_call.setBlockSizeY(4);
  __hmppcg_call.addLocalParameter((hmpprt::s32) (N_4), "N_3");
  __hmppcg_call.addLocalParameter(&epsilon, 8, "epsilon_2");
  __hmppcg_call.addLocalParameter(&beta_w_3, 8, "beta_w_4");
  __hmppcg_call.addLocalParameter(&beta_w_dry, 8, "beta_w_dry_2");
  __hmppcg_call.addLocalParameter(&beta_uh, 8, "beta_uh_2");
  __hmppcg_call.addLocalParameter(&beta_uh_dry_2, 8, "beta_uh_dry");
  __hmppcg_call.addLocalParameter(&beta_vh, 8, "beta_vh_2");
  __hmppcg_call.addLocalParameter(&beta_vh_dry, 8, "beta_vh_dry_2");
  __hmppcg_call.addLocalParameter((hmpprt::s32) (optimise_dry_cells), "optimise_dry_cells_2");
  __hmppcg_call.addLocalParameter(&surrogate_neighbours, 8, "surrogate_neighbours_2");
  __hmppcg_call.addLocalParameter(&number_of_boundaries, 8, "number_of_boundaries_2");
  __hmppcg_call.addLocalParameter(&centroid_coordinates, 8, "centroid_coordinates_2");
  __hmppcg_call.addLocalParameter(&stage_centroid_values, 8, "stage_centroid_values_2");
  __hmppcg_call.addLocalParameter(&bed_centroid_values, 8, "bed_centroid_values_2");
  __hmppcg_call.addLocalParameter(&xmom_centroid_values_4, 8, "xmom_centroid_values_3");
  __hmppcg_call.addLocalParameter(&ymom_centroid_values_41, 8, "ymom_centroid_values");
  __hmppcg_call.addLocalParameter(&vertex_coordinates, 8, "vertex_coordinates_2");
  __hmppcg_call.addLocalParameter(&stage_vertex_values, 8, "stage_vertex_values_2");
  __hmppcg_call.addLocalParameter(&xmom_vertex_values, 8, "xmom_vertex_values_2");
  __hmppcg_call.addLocalParameter(&ymom_vertex_values, 8, "ymom_vertex_values_2");
  __hmppcg_call.launch(extrapolate_second_order_sw_loop1D_2, hmpprt::Context::getInstance()->getCUDADevice());
 }
 ;
 # 631 "extrapolate_second_order_sw.c"
 if (extrapolate_velocity_second_order == 1)
 {
  # 117 "extrapolate_second_order_sw.c"
  if (N_4 - 1 >= 0)
  {
   hmpprt::CUDAGridCall __hmppcg_call;
   __hmppcg_call.setSizeX((N_4 - 1) / 128 + 1);
   __hmppcg_call.setSizeY(1);
   __hmppcg_call.setBlockSizeX(32);
   __hmppcg_call.setBlockSizeY(4);
   __hmppcg_call.addLocalParameter((hmpprt::s32) (N_4), "N");
   __hmppcg_call.addLocalParameter(&xmom_centroid_values_4, 8, "xmom_centroid_values_2");
   __hmppcg_call.addLocalParameter(&ymom_centroid_values_41, 8, "ymom_centroid_values_2");
   __hmppcg_call.addLocalParameter(&stage_vertex_values, 8, "stage_vertex_values_3");
   __hmppcg_call.addLocalParameter(&bed_vertex_values, 8, "bed_vertex_values_2");
   __hmppcg_call.addLocalParameter(&xmom_vertex_values, 8, "xmom_vertex_values_3");
   __hmppcg_call.addLocalParameter(&ymom_vertex_values, 8, "ymom_vertex_values_3");
   __hmppcg_call.addLocalParameter(&xmom_centroid_store, 8, "xmom_centroid_store_2");
   __hmppcg_call.addLocalParameter(&ymom_centroid_store_3, 8, "ymom_centroid_store_2");
   __hmppcg_call.launch(extrapolate_second_order_sw_loop1D_3, hmpprt::Context::getInstance()->getCUDADevice());
  }
  ;
 }
}
#endif // __CUDACC__



# 117 "extrapolate_second_order_sw.c"

#ifndef __CUDACC__
extern "C" CDLT_API  void extrapolate_second_order_sw(hmpprt::s32 N_41, hmpprt::s32 N2, hmpprt::s32 N3_31, hmpprt::s32 N6, double epsilon, double minimum_allowed_height, double beta_w_31, double beta_w_dry, double beta_uh, double beta_uh_dry_21, double beta_vh, double beta_vh_dry, hmpprt::s32 optimise_dry_cells, hmpprt::s32 extrapolate_velocity_second_order, hmpprt::s64* surrogate_neighbours, hmpprt::s64* number_of_boundaries, double* centroid_coordinates, double* stage_centroid_values, double* bed_centroid_values, double* xmom_centroid_values_41, double* ymom_centroid_values_4, double* vertex_coordinates, double* stage_vertex_values, double* bed_vertex_values, double* xmom_vertex_values, double* ymom_vertex_values, double* stage_centroid_store, double* xmom_centroid_store, double* ymom_centroid_store_31)
{
 # 1 "<preprocessor>"
 (extrapolate_second_order_sw_internal_1(N_41, N2, N3_31, N6, epsilon, minimum_allowed_height, beta_w_31, beta_w_dry, beta_uh, beta_uh_dry_21, beta_vh, beta_vh_dry, optimise_dry_cells, extrapolate_velocity_second_order, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s64> (surrogate_neighbours), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s64> (number_of_boundaries), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (centroid_coordinates), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (stage_centroid_values), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (bed_centroid_values), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (xmom_centroid_values_41), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (ymom_centroid_values_4), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (vertex_coordinates), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (stage_vertex_values), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (bed_vertex_values), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (xmom_vertex_values), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (ymom_vertex_values), stage_centroid_store, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (xmom_centroid_store), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (ymom_centroid_store_31)));
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
      extrapolate_second_order_sw_loop1D_1 = new hmpprt::CUDAGrid(hmpprt_module, "extrapolate_second_order_sw_loop1D_1");
      extrapolate_second_order_sw_loop1D_2 = new hmpprt::CUDAGrid(hmpprt_module, "extrapolate_second_order_sw_loop1D_2");
      extrapolate_second_order_sw_loop1D_3 = new hmpprt::CUDAGrid(hmpprt_module, "extrapolate_second_order_sw_loop1D_3");

    }
    hmpprt::Context::getInstance()->getGrouplet()->setTarget(hmpprt::CUDA);
    hmpprt::Context::getInstance()->getGrouplet()->addSignature("find_qmin_and_qmax", "prototype find_qmin_and_qmax(dq0: double, dq1: double, dq2: double, qmin: ^host double, qmax: ^host double) : s32");
    hmpprt::Context::getInstance()->getGrouplet()->addSignature("limit_gradient_old", "prototype limit_gradient_old(dqv: ^host double, qmin: double, qmax: double, beta_w: double) : s32");
    hmpprt::Context::getInstance()->getGrouplet()->addSignature("extrapolate_second_order_sw", "prototype extrapolate_second_order_sw(N: s32, N2: s32, N3: s32, N6: s32, epsilon: double, minimum_allowed_height: double, beta_w: double, beta_w_dry: double, beta_uh: double, beta_uh_dry: double, beta_vh: double, beta_vh_dry: double, optimise_dry_cells: s32, extrapolate_velocity_second_order: s32, surrogate_neighbours: ^cudaglob s64, number_of_boundaries: ^cudaglob s64, centroid_coordinates: ^cudaglob double, stage_centroid_values: ^cudaglob double, bed_centroid_values: ^cudaglob double, xmom_centroid_values: ^cudaglob double, ymom_centroid_values: ^cudaglob double, vertex_coordinates: ^cudaglob double, stage_vertex_values: ^cudaglob double, bed_vertex_values: ^cudaglob double, xmom_vertex_values: ^cudaglob double, ymom_vertex_values: ^cudaglob double, stage_centroid_store: unused ^host double, xmom_centroid_store: ^cudaglob double, ymom_centroid_store: ^cudaglob double)");

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
      delete extrapolate_second_order_sw_loop1D_1;
      delete extrapolate_second_order_sw_loop1D_2;
      delete extrapolate_second_order_sw_loop1D_3;

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


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

# 176 "swb2_domain_ext.c"

#ifndef __CUDACC__
extern "C" CDLT_API  void extrapolate_second_order_edge_sw(hmpprt::s32 number_of_elements, hmpprt::s32 optimise_dry_cells, hmpprt::s32 extrapolate_velocity_second_order, double epsilon, double minimum_allowed_height, double beta_w, double beta_w_dry, double beta_uh, double beta_uh_dry, double beta_vh, double beta_vh_dry, hmpprt::s64* surrogate_neighbours, hmpprt::s64* number_of_boundaries, double* centroid_coordinates, double* stage_centroid_values, double* elevation_centroid_values, double* xmom_centroid_values, double* ymom_centroid_values, double* edge_coordinates, double* stage_edge_values, double* elevation_edge_values, double* xmom_edge_values, double* ymom_edge_values, double* stage_vertex_values, double* xmom_vertex_values, double* ymom_vertex_values, double* elevation_vertex_values, double* stage_centroid_store, double* xmom_centroid_store, double* ymom_centroid_store, double* min_elevation_edgevalue, double* max_elevation_edgevalue, hmpprt::s32* count_wet_neighbours)
;
#endif // __CUDACC__



# 176 "swb2_domain_ext.c"

#ifndef __CUDACC__
void extrapolate_second_order_edge_sw_internal_1(hmpprt::s32 number_of_elements, hmpprt::s32 optimise_dry_cells, hmpprt::s32 extrapolate_velocity_second_order, double epsilon, double minimum_allowed_height, double beta_w, double beta_w_dry, double beta_uh, double beta_uh_dry, double beta_vh, double beta_vh_dry, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s64>  surrogate_neighbours, hmpprt::s64* number_of_boundaries, double* centroid_coordinates, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  stage_centroid_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  elevation_centroid_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  xmom_centroid_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  ymom_centroid_values, double* edge_coordinates, double* stage_edge_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  elevation_edge_values, double* xmom_edge_values, double* ymom_edge_values, double* stage_vertex_values, double* xmom_vertex_values, double* ymom_vertex_values, double* elevation_vertex_values, double* stage_centroid_store, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  xmom_centroid_store, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  ymom_centroid_store, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  min_elevation_edgevalue, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  max_elevation_edgevalue, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s32>  count_wet_neighbours)
;
#endif // __CUDACC__



# 108 "swb2_domain_ext.c"

#ifndef __CUDACC__
static hmpprt::CUDAGrid * extrapolate_second_order_edge_sw_loop1D_2 = 0;
#else

extern "C" __global__ void extrapolate_second_order_edge_sw_loop1D_2(hmpprt::s32 number_of_elements, hmpprt::s64* surrogate_neighbours, double* stage_centroid_values, double* max_elevation_edgevalue, hmpprt::s32* count_wet_neighbours);
#endif // __CUDACC__




# 108 "swb2_domain_ext.c"

#ifndef __CUDACC__
static hmpprt::CUDAGrid * extrapolate_second_order_edge_sw_loop1D_1 = 0;
#else

extern "C" __global__ void extrapolate_second_order_edge_sw_loop1D_1(hmpprt::s32 number_of_elements, double minimum_allowed_height, double* stage_centroid_values, double* elevation_centroid_values, double* xmom_centroid_values, double* ymom_centroid_values, double* elevation_edge_values, double* xmom_centroid_store, double* ymom_centroid_store, double* min_elevation_edgevalue, double* max_elevation_edgevalue);
#endif // __CUDACC__




# 108 "swb2_domain_ext.c"

#ifndef __CUDACC__
extern "C" CDLT_API  hmpprt::s32 _limit_gradient(double* dqv, double qmin_41, double qmax_41, double beta_w_31)
;
#endif // __CUDACC__



# 108 "swb2_domain_ext.c"

#ifndef __CUDACC__
hmpprt::s32 _limit_gradient_internal_1(double* dqv, double qmin_4, double qmax_4, double beta_w_3)
;
#endif // __CUDACC__



# 86 "swb2_domain_ext.c"

#ifndef __CUDACC__
extern "C" CDLT_API  hmpprt::s32 _find_qmin_and_qmax(double dq0, double dq1, double dq2, double* qmin, double* qmax)
;
#endif // __CUDACC__



# 86 "swb2_domain_ext.c"

#ifndef __CUDACC__
hmpprt::s32 _find_qmin_and_qmax_internal_1(double dq0, double dq1, double dq2, double* qmin, double* qmax)
;
#endif // __CUDACC__



# 86 "swb2_domain_ext.c"

#ifndef __CUDACC__
hmpprt::s32 _find_qmin_and_qmax_internal_1(double dq0, double dq1, double dq2, double* qmin, double* qmax)
{
 # 100 "swb2_domain_ext.c"
 *qmax = fmax(fmax(dq0, fmax(dq0 + dq1, dq0 + dq2)), (double) 0.0);
 # 101 "swb2_domain_ext.c"
 *qmin = fmin(fmin(dq0, fmin(dq0 + dq1, dq0 + dq2)), (double) 0.0);
 # 86 "swb2_domain_ext.c"
 return 0;
}
#endif // __CUDACC__



# 86 "swb2_domain_ext.c"

#ifndef __CUDACC__
extern "C" CDLT_API  hmpprt::s32 _find_qmin_and_qmax(double dq0, double dq1, double dq2, double* qmin, double* qmax)
{
 # 108 "swb2_domain_ext.c"
 (_find_qmin_and_qmax_internal_1(dq0, dq1, dq2, qmin, qmax));
}
#endif // __CUDACC__



# 108 "swb2_domain_ext.c"

#ifndef __CUDACC__
hmpprt::s32 _limit_gradient_internal_1(double* dqv, double qmin_4, double qmax_4, double beta_w_3)
{
 # 117 "swb2_domain_ext.c"
 double r;
 # 117 "swb2_domain_ext.c"
 double r0;
 # 117 "swb2_domain_ext.c"
 double phi;
 # 117 "swb2_domain_ext.c"
 r0 = (double) 1.0;
 # 118 "swb2_domain_ext.c"
 double TINY;
 # 136 "swb2_domain_ext.c"
 if (*dqv <  - TINY)
 {
  # 137 "swb2_domain_ext.c"
  r0 = qmin_4 / *dqv;
 }
 # 139 "swb2_domain_ext.c"
 if (*dqv > TINY)
 {
  # 140 "swb2_domain_ext.c"
  r0 = qmax_4 / *dqv;
 }
 # 142 "swb2_domain_ext.c"
 r = fmin(r0, (double) 1000.0);
 # 144 "swb2_domain_ext.c"
 if (*(dqv + 1) <  - TINY)
 {
  # 145 "swb2_domain_ext.c"
  r0 = qmin_4 / *(dqv + 1);
 }
 # 147 "swb2_domain_ext.c"
 if (*(dqv + 1) > TINY)
 {
  # 148 "swb2_domain_ext.c"
  r0 = qmax_4 / *(dqv + 1);
 }
 # 150 "swb2_domain_ext.c"
 r = fmin(r0, r);
 # 152 "swb2_domain_ext.c"
 if (*(dqv + 2) <  - TINY)
 {
  # 153 "swb2_domain_ext.c"
  r0 = qmin_4 / *(dqv + 2);
 }
 # 155 "swb2_domain_ext.c"
 if (*(dqv + 2) > TINY)
 {
  # 156 "swb2_domain_ext.c"
  r0 = qmax_4 / *(dqv + 2);
 }
 # 158 "swb2_domain_ext.c"
 r = fmin(r0, r);
 # 161 "swb2_domain_ext.c"
 phi = fmin(r * beta_w_3, (double) 1.0);
 # 163 "swb2_domain_ext.c"
 *dqv = *dqv * phi;
 # 164 "swb2_domain_ext.c"
 *(dqv + 1) = *(dqv + 1) * phi;
 # 165 "swb2_domain_ext.c"
 *(dqv + 2) = *(dqv + 2) * phi;
 # 108 "swb2_domain_ext.c"
 return 0;
}
#endif // __CUDACC__



# 108 "swb2_domain_ext.c"

#ifndef __CUDACC__
extern "C" CDLT_API  hmpprt::s32 _limit_gradient(double* dqv, double qmin_41, double qmax_41, double beta_w_31)
{
 # 176 "swb2_domain_ext.c"
 (_limit_gradient_internal_1(dqv, qmin_41, qmax_41, beta_w_31));
}
#endif // __CUDACC__



# 176 "swb2_domain_ext.c"

#ifdef __CUDACC__

extern "C" __global__ void extrapolate_second_order_edge_sw_loop1D_1(hmpprt::s32 number_of_elements, double minimum_allowed_height, double* stage_centroid_values, double* elevation_centroid_values, double* xmom_centroid_values, double* ymom_centroid_values, double* elevation_edge_values, double* xmom_centroid_store, double* ymom_centroid_store, double* min_elevation_edgevalue, double* max_elevation_edgevalue)
{
 # 222 "swb2_domain_ext.c"
 double dk_1;
 # 222 "swb2_domain_ext.c"
 hmpprt::s32 k_1;
 # 236 "swb2_domain_ext.c"
 k_1 = (hmpprt::gr_atidf());
 # 236 "swb2_domain_ext.c"
 if (k_1 > number_of_elements - 1)
 {
  # 236 "swb2_domain_ext.c"
  goto __hmppcg_label_1;
 }
 # 236 "swb2_domain_ext.c"
 dk_1 = fmax(*(stage_centroid_values + k_1) - *(elevation_centroid_values + k_1), minimum_allowed_height);
 # 237 "swb2_domain_ext.c"
 *(xmom_centroid_store + k_1) = *(xmom_centroid_values + k_1);
 # 238 "swb2_domain_ext.c"
 *(xmom_centroid_values + k_1) = *(xmom_centroid_values + k_1) / dk_1;
 # 240 "swb2_domain_ext.c"
 *(ymom_centroid_store + k_1) = *(ymom_centroid_values + k_1);
 # 241 "swb2_domain_ext.c"
 *(ymom_centroid_values + k_1) = *(ymom_centroid_values + k_1) / dk_1;
 # 245 "swb2_domain_ext.c"
 *(min_elevation_edgevalue + k_1) = fmin(*(elevation_edge_values + 3 * k_1), fmin(*(elevation_edge_values + (3 * k_1 + 1)), *(elevation_edge_values + (3 * k_1 + 2))));
 # 248 "swb2_domain_ext.c"
 *(max_elevation_edgevalue + k_1) = fmax(*(elevation_edge_values + 3 * k_1), fmax(*(elevation_edge_values + (3 * k_1 + 1)), *(elevation_edge_values + (3 * k_1 + 2))));
 # 176 "swb2_domain_ext.c"
 __hmppcg_label_1:;
}
#endif // __CUDACC__



# 176 "swb2_domain_ext.c"

#ifdef __CUDACC__

extern "C" __global__ void extrapolate_second_order_edge_sw_loop1D_2(hmpprt::s32 number_of_elements, hmpprt::s64* surrogate_neighbours, double* stage_centroid_values, double* max_elevation_edgevalue, hmpprt::s32* count_wet_neighbours)
{
 # 222 "swb2_domain_ext.c"
 hmpprt::s32 k_2;
 # 255 "swb2_domain_ext.c"
 k_2 = (hmpprt::gr_atidf());
 # 255 "swb2_domain_ext.c"
 if (k_2 > number_of_elements - 1)
 {
  # 255 "swb2_domain_ext.c"
  goto __hmppcg_label_2;
 }
 # 255 "swb2_domain_ext.c"
 *(count_wet_neighbours + k_2) = 0;
 # 256 "swb2_domain_ext.c"
 hmpprt::s32 i_1;
 # 256 "swb2_domain_ext.c"
 # 256 "swb2_domain_ext.c"
 for (i_1 = 0 ; i_1 <= 2 ; i_1 = i_1 + 1)
 {
  # 257 "swb2_domain_ext.c"
  hmpprt::s32 ktmp_1;
  # 257 "swb2_domain_ext.c"
  ktmp_1 = (hmpprt::s32 ) (*(surrogate_neighbours + (3 * k_2 + i_1)));
  # 258 "swb2_domain_ext.c"
  if (*(stage_centroid_values + ktmp_1) > *(max_elevation_edgevalue + ktmp_1))
  {
   # 259 "swb2_domain_ext.c"
   *(count_wet_neighbours + k_2) = *(count_wet_neighbours + k_2) + 1;
  }
 }
 # 176 "swb2_domain_ext.c"
 # 176 "swb2_domain_ext.c"
 __hmppcg_label_2:;
}
#endif // __CUDACC__



# 176 "swb2_domain_ext.c"

#ifndef __CUDACC__
void extrapolate_second_order_edge_sw_internal_1(hmpprt::s32 number_of_elements, hmpprt::s32 optimise_dry_cells, hmpprt::s32 extrapolate_velocity_second_order, double epsilon, double minimum_allowed_height, double beta_w, double beta_w_dry, double beta_uh, double beta_uh_dry, double beta_vh, double beta_vh_dry, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s64>  surrogate_neighbours, hmpprt::s64* number_of_boundaries, double* centroid_coordinates, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  stage_centroid_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  elevation_centroid_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  xmom_centroid_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  ymom_centroid_values, double* edge_coordinates, double* stage_edge_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  elevation_edge_values, double* xmom_edge_values, double* ymom_edge_values, double* stage_vertex_values, double* xmom_vertex_values, double* ymom_vertex_values, double* elevation_vertex_values, double* stage_centroid_store, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  xmom_centroid_store, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  ymom_centroid_store, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  min_elevation_edgevalue, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  max_elevation_edgevalue, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s32>  count_wet_neighbours)
{
 # 225 "swb2_domain_ext.c"
 double dqv_2[3uLL];
 # 225 "swb2_domain_ext.c"
 double qmin_3;
 # 225 "swb2_domain_ext.c"
 double qmax_3;
 # 227 "swb2_domain_ext.c"
 double de[3uLL];
 # 231 "swb2_domain_ext.c"
 if (extrapolate_velocity_second_order == 1)
 {
  # 265 "swb2_domain_ext.c"
  if (number_of_elements - 1 >= 0)
  {
   hmpprt::CUDAGridCall __hmppcg_call;
   __hmppcg_call.setSizeX((number_of_elements - 1) / 128 + 1);
   __hmppcg_call.setSizeY(1);
   __hmppcg_call.setBlockSizeX(32);
   __hmppcg_call.setBlockSizeY(4);
   __hmppcg_call.addLocalParameter((hmpprt::s32) (number_of_elements), "number_of_elements");
   __hmppcg_call.addLocalParameter(&minimum_allowed_height, 8, "minimum_allowed_height");
   __hmppcg_call.addLocalParameter(&stage_centroid_values, 8, "stage_centroid_values");
   __hmppcg_call.addLocalParameter(&elevation_centroid_values, 8, "elevation_centroid_values");
   __hmppcg_call.addLocalParameter(&xmom_centroid_values, 8, "xmom_centroid_values");
   __hmppcg_call.addLocalParameter(&ymom_centroid_values, 8, "ymom_centroid_values");
   __hmppcg_call.addLocalParameter(&elevation_edge_values, 8, "elevation_edge_values");
   __hmppcg_call.addLocalParameter(&xmom_centroid_store, 8, "xmom_centroid_store");
   __hmppcg_call.addLocalParameter(&ymom_centroid_store, 8, "ymom_centroid_store");
   __hmppcg_call.addLocalParameter(&min_elevation_edgevalue, 8, "min_elevation_edgevalue");
   __hmppcg_call.addLocalParameter(&max_elevation_edgevalue, 8, "max_elevation_edgevalue");
   __hmppcg_call.launch(extrapolate_second_order_edge_sw_loop1D_1, hmpprt::Context::getInstance()->getCUDADevice());
  }
  ;
 }
 # 265 "swb2_domain_ext.c"
 if (number_of_elements - 1 >= 0)
 {
  hmpprt::CUDAGridCall __hmppcg_call;
  __hmppcg_call.setSizeX((number_of_elements - 1) / 128 + 1);
  __hmppcg_call.setSizeY(1);
  __hmppcg_call.setBlockSizeX(32);
  __hmppcg_call.setBlockSizeY(4);
  __hmppcg_call.addLocalParameter((hmpprt::s32) (number_of_elements), "number_of_elements");
  __hmppcg_call.addLocalParameter(&surrogate_neighbours, 8, "surrogate_neighbours");
  __hmppcg_call.addLocalParameter(&stage_centroid_values, 8, "stage_centroid_values");
  __hmppcg_call.addLocalParameter(&max_elevation_edgevalue, 8, "max_elevation_edgevalue");
  __hmppcg_call.addLocalParameter(&count_wet_neighbours, 8, "count_wet_neighbours");
  __hmppcg_call.launch(extrapolate_second_order_edge_sw_loop1D_2, hmpprt::Context::getInstance()->getCUDADevice());
 }
 ;
 # 265 "swb2_domain_ext.c"
 hmpprt::s32 k_3;
 # 265 "swb2_domain_ext.c"
 hmpprt::s32 end_6;
 # 265 "swb2_domain_ext.c"
 # 266 "swb2_domain_ext.c"
 # 266 "swb2_domain_ext.c"
 for (k_3 = 0, end_6 = number_of_elements - 1 ; k_3 <= end_6 ; k_3 = k_3 + 1)
 {
  # 267 "swb2_domain_ext.c"
  hmpprt::s32 k3_1;
  # 267 "swb2_domain_ext.c"
  hmpprt::s32 k6_1;
  # 267 "swb2_domain_ext.c"
  double x_1;
  # 267 "swb2_domain_ext.c"
  double y_1;
  # 267 "swb2_domain_ext.c"
  double dxv0_1;
  # 267 "swb2_domain_ext.c"
  double dyv0_1;
  # 267 "swb2_domain_ext.c"
  double dxv1_1;
  # 267 "swb2_domain_ext.c"
  double dyv1_1;
  # 267 "swb2_domain_ext.c"
  double dxv2_1;
  # 267 "swb2_domain_ext.c"
  double dyv2_1;
  # 267 "swb2_domain_ext.c"
  k3_1 = k_3 * 3;
  # 268 "swb2_domain_ext.c"
  k6_1 = k_3 * 6;
  # 270 "swb2_domain_ext.c"
  if (*(number_of_boundaries + k_3) == 3LL)
  {
   # 275 "swb2_domain_ext.c"
   *(stage_edge_values + k3_1) = *(stage_centroid_values + k_3);
   # 276 "swb2_domain_ext.c"
   *(stage_edge_values + (k3_1 + 1)) = *(stage_centroid_values + k_3);
   # 277 "swb2_domain_ext.c"
   *(stage_edge_values + (k3_1 + 2)) = *(stage_centroid_values + k_3);
   # 278 "swb2_domain_ext.c"
   *(xmom_edge_values + k3_1) = *(xmom_centroid_values + k_3);
   # 279 "swb2_domain_ext.c"
   *(xmom_edge_values + (k3_1 + 1)) = *(xmom_centroid_values + k_3);
   # 280 "swb2_domain_ext.c"
   *(xmom_edge_values + (k3_1 + 2)) = *(xmom_centroid_values + k_3);
   # 281 "swb2_domain_ext.c"
   *(ymom_edge_values + k3_1) = *(ymom_centroid_values + k_3);
   # 282 "swb2_domain_ext.c"
   *(ymom_edge_values + (k3_1 + 1)) = *(ymom_centroid_values + k_3);
   # 283 "swb2_domain_ext.c"
   *(ymom_edge_values + (k3_1 + 2)) = *(ymom_centroid_values + k_3);
   # 285 "swb2_domain_ext.c"
   continue;
  }
  else
  {
   # 293 "swb2_domain_ext.c"
   hmpprt::s32 coord_index_1;
   # 293 "swb2_domain_ext.c"
   double xv0_1;
   # 293 "swb2_domain_ext.c"
   double xv1_1;
   # 293 "swb2_domain_ext.c"
   double xv2_1;
   # 293 "swb2_domain_ext.c"
   double yv0_1;
   # 293 "swb2_domain_ext.c"
   double yv1_1;
   # 293 "swb2_domain_ext.c"
   double yv2_1;
   # 293 "swb2_domain_ext.c"
   xv0_1 = *(edge_coordinates + k6_1);
   # 294 "swb2_domain_ext.c"
   yv0_1 = *(edge_coordinates + (k6_1 + 1));
   # 295 "swb2_domain_ext.c"
   xv1_1 = *(edge_coordinates + (k6_1 + 2));
   # 296 "swb2_domain_ext.c"
   yv1_1 = *(edge_coordinates + (k6_1 + 3));
   # 297 "swb2_domain_ext.c"
   xv2_1 = *(edge_coordinates + (k6_1 + 4));
   # 298 "swb2_domain_ext.c"
   yv2_1 = *(edge_coordinates + (k6_1 + 5));
   # 301 "swb2_domain_ext.c"
   coord_index_1 = 2 * k_3;
   # 302 "swb2_domain_ext.c"
   x_1 = *(centroid_coordinates + coord_index_1);
   # 303 "swb2_domain_ext.c"
   y_1 = *(centroid_coordinates + (coord_index_1 + 1));
   # 307 "swb2_domain_ext.c"
   dxv0_1 = xv0_1 - x_1;
   # 308 "swb2_domain_ext.c"
   dxv1_1 = xv1_1 - x_1;
   # 309 "swb2_domain_ext.c"
   dxv2_1 = xv2_1 - x_1;
   # 310 "swb2_domain_ext.c"
   dyv0_1 = yv0_1 - y_1;
   # 311 "swb2_domain_ext.c"
   dyv1_1 = yv1_1 - y_1;
   # 312 "swb2_domain_ext.c"
   dyv2_1 = yv2_1 - y_1;
  }
  # 320 "swb2_domain_ext.c"
  if (*(number_of_boundaries + k_3) <= 1LL)
  {
   # 332 "swb2_domain_ext.c"
   hmpprt::s32 k0_1;
   # 332 "swb2_domain_ext.c"
   hmpprt::s32 k1_1;
   # 332 "swb2_domain_ext.c"
   hmpprt::s32 k2_1;
   # 332 "swb2_domain_ext.c"
   double stagemin_1;
   # 332 "swb2_domain_ext.c"
   double bedmax_1;
   # 332 "swb2_domain_ext.c"
   hmpprt::s32 coord_index_2;
   # 332 "swb2_domain_ext.c"
   hmpprt::s32 coord_index_3;
   # 332 "swb2_domain_ext.c"
   hmpprt::s32 coord_index_4;
   # 332 "swb2_domain_ext.c"
   double x1_1;
   # 332 "swb2_domain_ext.c"
   double x0_1;
   # 332 "swb2_domain_ext.c"
   double x2_1;
   # 332 "swb2_domain_ext.c"
   double y1_1;
   # 332 "swb2_domain_ext.c"
   double y0_1;
   # 332 "swb2_domain_ext.c"
   double y2_1;
   # 332 "swb2_domain_ext.c"
   double dy2_1;
   # 332 "swb2_domain_ext.c"
   double dx1_1;
   # 332 "swb2_domain_ext.c"
   double dy1_1;
   # 332 "swb2_domain_ext.c"
   double dx2_1;
   # 332 "swb2_domain_ext.c"
   double area2_1;
   # 332 "swb2_domain_ext.c"
   double h0_1;
   # 332 "swb2_domain_ext.c"
   double h1_1;
   # 332 "swb2_domain_ext.c"
   double h2_1;
   # 332 "swb2_domain_ext.c"
   double hc_1;
   # 332 "swb2_domain_ext.c"
   double hmin_1;
   # 332 "swb2_domain_ext.c"
   double dq1_3;
   # 332 "swb2_domain_ext.c"
   double dq2_3;
   # 332 "swb2_domain_ext.c"
   double a_1;
   # 332 "swb2_domain_ext.c"
   double inv_area2_1;
   # 332 "swb2_domain_ext.c"
   double b_1;
   # 332 "swb2_domain_ext.c"
   double a_2;
   # 332 "swb2_domain_ext.c"
   double b_2;
   # 332 "swb2_domain_ext.c"
   double dq0_3;
   # 332 "swb2_domain_ext.c"
   double hfactor_1;
   # 332 "swb2_domain_ext.c"
   double beta_tmp_1;
   # 332 "swb2_domain_ext.c"
   double dq1_4;
   # 332 "swb2_domain_ext.c"
   double dq2_4;
   # 332 "swb2_domain_ext.c"
   double a_3;
   # 332 "swb2_domain_ext.c"
   double b_3;
   # 332 "swb2_domain_ext.c"
   double a_4;
   # 332 "swb2_domain_ext.c"
   double b_4;
   # 332 "swb2_domain_ext.c"
   double dq0_4;
   # 332 "swb2_domain_ext.c"
   double beta_tmp_2;
   # 332 "swb2_domain_ext.c"
   double dq1_5;
   # 332 "swb2_domain_ext.c"
   double dq2_5;
   # 332 "swb2_domain_ext.c"
   double a_5;
   # 332 "swb2_domain_ext.c"
   double b_5;
   # 332 "swb2_domain_ext.c"
   double a_6;
   # 332 "swb2_domain_ext.c"
   double b_6;
   # 332 "swb2_domain_ext.c"
   double dq0_5;
   # 332 "swb2_domain_ext.c"
   double beta_tmp_3;
   # 332 "swb2_domain_ext.c"
   k0_1 = (hmpprt::s32 ) (*(surrogate_neighbours + k3_1));
   # 333 "swb2_domain_ext.c"
   k1_1 = (hmpprt::s32 ) (*(surrogate_neighbours + (k3_1 + 1)));
   # 334 "swb2_domain_ext.c"
   k2_1 = (hmpprt::s32 ) (*(surrogate_neighbours + (k3_1 + 2)));
   # 341 "swb2_domain_ext.c"
   bedmax_1 = fmax(*(elevation_centroid_values + k_3), fmax(*(elevation_centroid_values + k0_1), fmax(*(elevation_centroid_values + k1_1), *(elevation_centroid_values + k2_1))));
   # 350 "swb2_domain_ext.c"
   stagemin_1 = fmin(fmax(*(stage_centroid_values + k_3), *(elevation_centroid_values + k_3)), fmin(fmax(*(stage_centroid_values + k0_1), *(elevation_centroid_values + k0_1)), fmin(fmax(*(stage_centroid_values + k1_1), *(elevation_centroid_values + k1_1)), fmax(*(stage_centroid_values + k2_1), *(elevation_centroid_values + k2_1)))));
   # 351 "swb2_domain_ext.c"
   if (stagemin_1 < bedmax_1)
   {
    # 353 "swb2_domain_ext.c"
    k2_1 = k_3;
    # 354 "swb2_domain_ext.c"
    k0_1 = k_3;
    # 355 "swb2_domain_ext.c"
    k1_1 = k_3;
   }
   # 360 "swb2_domain_ext.c"
   coord_index_2 = 2 * k0_1;
   # 361 "swb2_domain_ext.c"
   x0_1 = *(centroid_coordinates + coord_index_2);
   # 362 "swb2_domain_ext.c"
   y0_1 = *(centroid_coordinates + (coord_index_2 + 1));
   # 364 "swb2_domain_ext.c"
   coord_index_3 = 2 * k1_1;
   # 365 "swb2_domain_ext.c"
   x1_1 = *(centroid_coordinates + coord_index_3);
   # 366 "swb2_domain_ext.c"
   y1_1 = *(centroid_coordinates + (coord_index_3 + 1));
   # 368 "swb2_domain_ext.c"
   coord_index_4 = 2 * k2_1;
   # 369 "swb2_domain_ext.c"
   x2_1 = *(centroid_coordinates + coord_index_4);
   # 370 "swb2_domain_ext.c"
   y2_1 = *(centroid_coordinates + (coord_index_4 + 1));
   # 374 "swb2_domain_ext.c"
   dx1_1 = x1_1 - x0_1;
   # 375 "swb2_domain_ext.c"
   dx2_1 = x2_1 - x0_1;
   # 376 "swb2_domain_ext.c"
   dy1_1 = y1_1 - y0_1;
   # 377 "swb2_domain_ext.c"
   dy2_1 = y2_1 - y0_1;
   # 381 "swb2_domain_ext.c"
   area2_1 = dy2_1 * dx1_1 - dy1_1 * dx2_1;
   # 384 "swb2_domain_ext.c"
   if (area2_1 <= (double) 0.0)
   {
    # 389 "swb2_domain_ext.c"
    *(stage_edge_values + k3_1) = *(stage_centroid_values + k_3);
    # 390 "swb2_domain_ext.c"
    *(stage_edge_values + (k3_1 + 1)) = *(stage_centroid_values + k_3);
    # 391 "swb2_domain_ext.c"
    *(stage_edge_values + (k3_1 + 2)) = *(stage_centroid_values + k_3);
    # 393 "swb2_domain_ext.c"
    *(xmom_edge_values + k3_1) = *(xmom_centroid_values + k_3);
    # 394 "swb2_domain_ext.c"
    *(xmom_edge_values + (k3_1 + 1)) = *(xmom_centroid_values + k_3);
    # 395 "swb2_domain_ext.c"
    *(xmom_edge_values + (k3_1 + 2)) = *(xmom_centroid_values + k_3);
    # 396 "swb2_domain_ext.c"
    *(ymom_edge_values + k3_1) = *(ymom_centroid_values + k_3);
    # 397 "swb2_domain_ext.c"
    *(ymom_edge_values + (k3_1 + 1)) = *(ymom_centroid_values + k_3);
    # 398 "swb2_domain_ext.c"
    *(ymom_edge_values + (k3_1 + 2)) = *(ymom_centroid_values + k_3);
    # 400 "swb2_domain_ext.c"
    continue;
   }
   # 404 "swb2_domain_ext.c"
   hc_1 = *(stage_centroid_values + k_3) - *(elevation_centroid_values + k_3);
   # 405 "swb2_domain_ext.c"
   h0_1 = *(stage_centroid_values + k0_1) - *(elevation_centroid_values + k0_1);
   # 406 "swb2_domain_ext.c"
   h1_1 = *(stage_centroid_values + k1_1) - *(elevation_centroid_values + k1_1);
   # 407 "swb2_domain_ext.c"
   h2_1 = *(stage_centroid_values + k2_1) - *(elevation_centroid_values + k2_1);
   # 408 "swb2_domain_ext.c"
   hmin_1 = fmin(fmin(h0_1, fmin(h1_1, h2_1)), hc_1);
   # 410 "swb2_domain_ext.c"
   hfactor_1 = (double) 0.0;
   # 412 "swb2_domain_ext.c"
   if (hmin_1 > (double) 0.0)
   {
    # 415 "swb2_domain_ext.c"
    hfactor_1 = (double) 1.0;
   }
   # 426 "swb2_domain_ext.c"
   dq0_3 = *(stage_centroid_values + k0_1) - *(stage_centroid_values + k_3);
   # 430 "swb2_domain_ext.c"
   dq1_3 = *(stage_centroid_values + k1_1) - *(stage_centroid_values + k0_1);
   # 431 "swb2_domain_ext.c"
   dq2_3 = *(stage_centroid_values + k2_1) - *(stage_centroid_values + k0_1);
   # 433 "swb2_domain_ext.c"
   inv_area2_1 = (double) 1.0 / area2_1;
   # 435 "swb2_domain_ext.c"
   a_1 = dy2_1 * dq1_3 - dy1_1 * dq2_3;
   # 436 "swb2_domain_ext.c"
   a_2 = a_1 * inv_area2_1;
   # 437 "swb2_domain_ext.c"
   b_1 = dx1_1 * dq2_3 - dx2_1 * dq1_3;
   # 438 "swb2_domain_ext.c"
   b_2 = b_1 * inv_area2_1;
   # 442 "swb2_domain_ext.c"
   dqv_2[0] = a_2 * dxv0_1 + b_2 * dyv0_1;
   # 443 "swb2_domain_ext.c"
   dqv_2[1] = a_2 * dxv1_1 + b_2 * dyv1_1;
   # 444 "swb2_domain_ext.c"
   dqv_2[2] = a_2 * dxv2_1 + b_2 * dyv2_1;
   # 449 "swb2_domain_ext.c"
   (_find_qmin_and_qmax_internal_1(dq0_3, dq1_3, dq2_3, &qmin_3, &qmax_3));
   # 451 "swb2_domain_ext.c"
   beta_tmp_1 = beta_w_dry + (beta_w - beta_w_dry) * hfactor_1;
   # 455 "swb2_domain_ext.c"
   (_limit_gradient_internal_1(&dqv_2[0], qmin_3, qmax_3, beta_tmp_1));
   # 456 "swb2_domain_ext.c"
   *(stage_edge_values + k3_1) = *(stage_centroid_values + k_3) + dqv_2[0];
   # 457 "swb2_domain_ext.c"
   *(stage_edge_values + (k3_1 + 1)) = *(stage_centroid_values + k_3) + dqv_2[1];
   # 458 "swb2_domain_ext.c"
   *(stage_edge_values + (k3_1 + 2)) = *(stage_centroid_values + k_3) + dqv_2[2];
   # 467 "swb2_domain_ext.c"
   dq0_4 = *(xmom_centroid_values + k0_1) - *(xmom_centroid_values + k_3);
   # 471 "swb2_domain_ext.c"
   dq1_4 = *(xmom_centroid_values + k1_1) - *(xmom_centroid_values + k0_1);
   # 472 "swb2_domain_ext.c"
   dq2_4 = *(xmom_centroid_values + k2_1) - *(xmom_centroid_values + k0_1);
   # 475 "swb2_domain_ext.c"
   a_3 = dy2_1 * dq1_4 - dy1_1 * dq2_4;
   # 476 "swb2_domain_ext.c"
   a_4 = a_3 * inv_area2_1;
   # 477 "swb2_domain_ext.c"
   b_3 = dx1_1 * dq2_4 - dx2_1 * dq1_4;
   # 478 "swb2_domain_ext.c"
   b_4 = b_3 * inv_area2_1;
   # 482 "swb2_domain_ext.c"
   dqv_2[0] = a_4 * dxv0_1 + b_4 * dyv0_1;
   # 483 "swb2_domain_ext.c"
   dqv_2[1] = a_4 * dxv1_1 + b_4 * dyv1_1;
   # 484 "swb2_domain_ext.c"
   dqv_2[2] = a_4 * dxv2_1 + b_4 * dyv2_1;
   # 490 "swb2_domain_ext.c"
   (_find_qmin_and_qmax_internal_1(dq0_4, dq1_4, dq2_4, &qmin_3, &qmax_3));
   # 492 "swb2_domain_ext.c"
   beta_tmp_2 = beta_uh_dry + (beta_uh - beta_uh_dry) * hfactor_1;
   # 494 "swb2_domain_ext.c"
   (_limit_gradient_internal_1(&dqv_2[0], qmin_3, qmax_3, beta_tmp_2));
   # 497 "swb2_domain_ext.c"
   hmpprt::s32 i_2;
   # 497 "swb2_domain_ext.c"
   # 498 "swb2_domain_ext.c"
   for (i_2 = 0 ; i_2 <= 2 ; i_2 = i_2 + 1)
   {
    # 499 "swb2_domain_ext.c"
    *(xmom_edge_values + (k3_1 + i_2)) = *(xmom_centroid_values + k_3) + dqv_2[i_2];
   }
   # 508 "swb2_domain_ext.c"
   # 508 "swb2_domain_ext.c"
   dq0_5 = *(ymom_centroid_values + k0_1) - *(ymom_centroid_values + k_3);
   # 512 "swb2_domain_ext.c"
   dq1_5 = *(ymom_centroid_values + k1_1) - *(ymom_centroid_values + k0_1);
   # 513 "swb2_domain_ext.c"
   dq2_5 = *(ymom_centroid_values + k2_1) - *(ymom_centroid_values + k0_1);
   # 516 "swb2_domain_ext.c"
   a_5 = dy2_1 * dq1_5 - dy1_1 * dq2_5;
   # 517 "swb2_domain_ext.c"
   a_6 = a_5 * inv_area2_1;
   # 518 "swb2_domain_ext.c"
   b_5 = dx1_1 * dq2_5 - dx2_1 * dq1_5;
   # 519 "swb2_domain_ext.c"
   b_6 = b_5 * inv_area2_1;
   # 523 "swb2_domain_ext.c"
   dqv_2[0] = a_6 * dxv0_1 + b_6 * dyv0_1;
   # 524 "swb2_domain_ext.c"
   dqv_2[1] = a_6 * dxv1_1 + b_6 * dyv1_1;
   # 525 "swb2_domain_ext.c"
   dqv_2[2] = a_6 * dxv2_1 + b_6 * dyv2_1;
   # 531 "swb2_domain_ext.c"
   (_find_qmin_and_qmax_internal_1(dq0_5, dq1_5, dq2_5, &qmin_3, &qmax_3));
   # 533 "swb2_domain_ext.c"
   beta_tmp_3 = beta_vh_dry + (beta_vh - beta_vh_dry) * hfactor_1;
   # 535 "swb2_domain_ext.c"
   (_limit_gradient_internal_1(&dqv_2[0], qmin_3, qmax_3, beta_tmp_3));
   # 537 "swb2_domain_ext.c"
   hmpprt::s32 i_3;
   # 537 "swb2_domain_ext.c"
   # 538 "swb2_domain_ext.c"
   for (i_3 = 0 ; i_3 <= 2 ; i_3 = i_3 + 1)
   {
    # 539 "swb2_domain_ext.c"
    *(ymom_edge_values + (k3_1 + i_3)) = *(ymom_centroid_values + k_3) + dqv_2[i_3];
   }
   # 544 "swb2_domain_ext.c"
  }
  else
  {
   # 553 "swb2_domain_ext.c"
   hmpprt::s32 k2_2;
   # 553 "swb2_domain_ext.c"
   hmpprt::s32 k1_2;
   # 553 "swb2_domain_ext.c"
   hmpprt::s32 coord_index_5;
   # 553 "swb2_domain_ext.c"
   double x1_2;
   # 553 "swb2_domain_ext.c"
   double y1_2;
   # 553 "swb2_domain_ext.c"
   double dx1_2;
   # 553 "swb2_domain_ext.c"
   double dy1_2;
   # 553 "swb2_domain_ext.c"
   double area2_2;
   # 553 "swb2_domain_ext.c"
   double dx2_2;
   # 553 "swb2_domain_ext.c"
   double dq1_6;
   # 553 "swb2_domain_ext.c"
   double dx2_3;
   # 553 "swb2_domain_ext.c"
   double dy2_2;
   # 553 "swb2_domain_ext.c"
   double a_7;
   # 553 "swb2_domain_ext.c"
   double b_7;
   # 553 "swb2_domain_ext.c"
   double dq1_7;
   # 553 "swb2_domain_ext.c"
   double a_8;
   # 553 "swb2_domain_ext.c"
   double b_8;
   # 553 "swb2_domain_ext.c"
   double dq1_8;
   # 553 "swb2_domain_ext.c"
   double a_9;
   # 553 "swb2_domain_ext.c"
   double b_9;
   # 553 "swb2_domain_ext.c"
   # 553 "swb2_domain_ext.c"
   for (k2_2 = k3_1 ; k2_2 < k3_1 + 3 ; k2_2 = k2_2 + 1)
   {
    # 558 "swb2_domain_ext.c"
    if (*(surrogate_neighbours + k2_2) != (hmpprt::s64 ) (k_3))
    {
     # 560 "swb2_domain_ext.c"
     break;
    }
   }
   # 560 "swb2_domain_ext.c"
   # 571 "swb2_domain_ext.c"
   k1_2 = (hmpprt::s32 ) (*(surrogate_neighbours + k2_2));
   # 575 "swb2_domain_ext.c"
   coord_index_5 = 2 * k1_2;
   # 576 "swb2_domain_ext.c"
   x1_2 = *(centroid_coordinates + coord_index_5);
   # 577 "swb2_domain_ext.c"
   y1_2 = *(centroid_coordinates + (coord_index_5 + 1));
   # 581 "swb2_domain_ext.c"
   dx1_2 = x1_2 - x_1;
   # 582 "swb2_domain_ext.c"
   dy1_2 = y1_2 - y_1;
   # 585 "swb2_domain_ext.c"
   area2_2 = dx1_2 * dx1_2 + dy1_2 * dy1_2;
   # 591 "swb2_domain_ext.c"
   dx2_2 = (double) 1.0 / area2_2;
   # 592 "swb2_domain_ext.c"
   dy2_2 = dx2_2 * dy1_2;
   # 593 "swb2_domain_ext.c"
   dx2_3 = dx2_2 * dx1_2;
   # 601 "swb2_domain_ext.c"
   dq1_6 = *(stage_centroid_values + k1_2) - *(stage_centroid_values + k_3);
   # 605 "swb2_domain_ext.c"
   a_7 = dq1_6 * dx2_3;
   # 606 "swb2_domain_ext.c"
   b_7 = dq1_6 * dy2_2;
   # 609 "swb2_domain_ext.c"
   dqv_2[0] = a_7 * dxv0_1 + b_7 * dyv0_1;
   # 610 "swb2_domain_ext.c"
   dqv_2[1] = a_7 * dxv1_1 + b_7 * dyv1_1;
   # 611 "swb2_domain_ext.c"
   dqv_2[2] = a_7 * dxv2_1 + b_7 * dyv2_1;
   # 614 "swb2_domain_ext.c"
   if (dq1_6 >= (double) 0.0)
   {
    # 616 "swb2_domain_ext.c"
    qmin_3 = (double) 0.0;
    # 617 "swb2_domain_ext.c"
    qmax_3 = dq1_6;
   }
   else
   {
    # 621 "swb2_domain_ext.c"
    qmin_3 = dq1_6;
    # 622 "swb2_domain_ext.c"
    qmax_3 = (double) 0.0;
   }
   # 626 "swb2_domain_ext.c"
   (_limit_gradient_internal_1(&dqv_2[0], qmin_3, qmax_3, beta_w));
   # 630 "swb2_domain_ext.c"
   *(stage_edge_values + k3_1) = *(stage_centroid_values + k_3) + dqv_2[0];
   # 631 "swb2_domain_ext.c"
   *(stage_edge_values + (k3_1 + 1)) = *(stage_centroid_values + k_3) + dqv_2[1];
   # 632 "swb2_domain_ext.c"
   *(stage_edge_values + (k3_1 + 2)) = *(stage_centroid_values + k_3) + dqv_2[2];
   # 640 "swb2_domain_ext.c"
   dq1_7 = *(xmom_centroid_values + k1_2) - *(xmom_centroid_values + k_3);
   # 644 "swb2_domain_ext.c"
   a_8 = dq1_7 * dx2_3;
   # 645 "swb2_domain_ext.c"
   b_8 = dq1_7 * dy2_2;
   # 648 "swb2_domain_ext.c"
   dqv_2[0] = a_8 * dxv0_1 + b_8 * dyv0_1;
   # 649 "swb2_domain_ext.c"
   dqv_2[1] = a_8 * dxv1_1 + b_8 * dyv1_1;
   # 650 "swb2_domain_ext.c"
   dqv_2[2] = a_8 * dxv2_1 + b_8 * dyv2_1;
   # 653 "swb2_domain_ext.c"
   if (dq1_7 >= (double) 0.0)
   {
    # 655 "swb2_domain_ext.c"
    qmin_3 = (double) 0.0;
    # 656 "swb2_domain_ext.c"
    qmax_3 = dq1_7;
   }
   else
   {
    # 660 "swb2_domain_ext.c"
    qmin_3 = dq1_7;
    # 661 "swb2_domain_ext.c"
    qmax_3 = (double) 0.0;
   }
   # 665 "swb2_domain_ext.c"
   (_limit_gradient_internal_1(&dqv_2[0], qmin_3, qmax_3, beta_w));
   # 672 "swb2_domain_ext.c"
   hmpprt::s32 i_4;
   # 672 "swb2_domain_ext.c"
   # 673 "swb2_domain_ext.c"
   for (i_4 = 0 ; i_4 <= 2 ; i_4 = i_4 + 1)
   {
    # 674 "swb2_domain_ext.c"
    *(xmom_edge_values + (k3_1 + i_4)) = *(xmom_centroid_values + k_3) + dqv_2[i_4];
   }
   # 682 "swb2_domain_ext.c"
   # 682 "swb2_domain_ext.c"
   dq1_8 = *(ymom_centroid_values + k1_2) - *(ymom_centroid_values + k_3);
   # 686 "swb2_domain_ext.c"
   a_9 = dq1_8 * dx2_3;
   # 687 "swb2_domain_ext.c"
   b_9 = dq1_8 * dy2_2;
   # 690 "swb2_domain_ext.c"
   dqv_2[0] = a_9 * dxv0_1 + b_9 * dyv0_1;
   # 691 "swb2_domain_ext.c"
   dqv_2[1] = a_9 * dxv1_1 + b_9 * dyv1_1;
   # 692 "swb2_domain_ext.c"
   dqv_2[2] = a_9 * dxv2_1 + b_9 * dyv2_1;
   # 695 "swb2_domain_ext.c"
   if (dq1_8 >= (double) 0.0)
   {
    # 697 "swb2_domain_ext.c"
    qmin_3 = (double) 0.0;
    # 698 "swb2_domain_ext.c"
    qmax_3 = dq1_8;
   }
   else
   {
    # 702 "swb2_domain_ext.c"
    qmin_3 = dq1_8;
    # 703 "swb2_domain_ext.c"
    qmax_3 = (double) 0.0;
   }
   # 707 "swb2_domain_ext.c"
   (_limit_gradient_internal_1(&dqv_2[0], qmin_3, qmax_3, beta_w));
   # 709 "swb2_domain_ext.c"
   hmpprt::s32 i_5;
   # 709 "swb2_domain_ext.c"
   # 710 "swb2_domain_ext.c"
   for (i_5 = 0 ; i_5 <= 2 ; i_5 = i_5 + 1)
   {
    # 711 "swb2_domain_ext.c"
    *(ymom_edge_values + (k3_1 + i_5)) = *(ymom_centroid_values + k_3) + dqv_2[i_5];
   }
   # 717 "swb2_domain_ext.c"
  }
 }
 # 717 "swb2_domain_ext.c"
 # 717 "swb2_domain_ext.c"
 hmpprt::s32 k_4;
 # 717 "swb2_domain_ext.c"
 hmpprt::s32 end_9;
 # 717 "swb2_domain_ext.c"
 # 717 "swb2_domain_ext.c"
 # 717 "swb2_domain_ext.c"
 for (k_4 = 0, end_9 = number_of_elements - 1 ; k_4 <= end_9 ; k_4 = k_4 + 1)
 {
  # 718 "swb2_domain_ext.c"
  hmpprt::s32 k3_2;
  # 718 "swb2_domain_ext.c"
  k3_2 = 3 * k_4;
  # 721 "swb2_domain_ext.c"
  *(stage_vertex_values + k3_2) = *(stage_edge_values + (k3_2 + 1)) + *(stage_edge_values + (k3_2 + 2)) - *(stage_edge_values + k3_2);
  # 722 "swb2_domain_ext.c"
  *(stage_vertex_values + (k3_2 + 1)) = *(stage_edge_values + k3_2) + *(stage_edge_values + (k3_2 + 2)) - *(stage_edge_values + (k3_2 + 1));
  # 723 "swb2_domain_ext.c"
  *(stage_vertex_values + (k3_2 + 2)) = *(stage_edge_values + k3_2) + *(stage_edge_values + (k3_2 + 1)) - *(stage_edge_values + (k3_2 + 2));
  # 726 "swb2_domain_ext.c"
  *(xmom_vertex_values + k3_2) = *(xmom_edge_values + (k3_2 + 1)) + *(xmom_edge_values + (k3_2 + 2)) - *(xmom_edge_values + k3_2);
  # 727 "swb2_domain_ext.c"
  *(xmom_vertex_values + (k3_2 + 1)) = *(xmom_edge_values + k3_2) + *(xmom_edge_values + (k3_2 + 2)) - *(xmom_edge_values + (k3_2 + 1));
  # 728 "swb2_domain_ext.c"
  *(xmom_vertex_values + (k3_2 + 2)) = *(xmom_edge_values + k3_2) + *(xmom_edge_values + (k3_2 + 1)) - *(xmom_edge_values + (k3_2 + 2));
  # 731 "swb2_domain_ext.c"
  *(ymom_vertex_values + k3_2) = *(ymom_edge_values + (k3_2 + 1)) + *(ymom_edge_values + (k3_2 + 2)) - *(ymom_edge_values + k3_2);
  # 732 "swb2_domain_ext.c"
  *(ymom_vertex_values + (k3_2 + 1)) = *(ymom_edge_values + k3_2) + *(ymom_edge_values + (k3_2 + 2)) - *(ymom_edge_values + (k3_2 + 1));
  # 733 "swb2_domain_ext.c"
  *(ymom_vertex_values + (k3_2 + 2)) = *(ymom_edge_values + k3_2) + *(ymom_edge_values + (k3_2 + 1)) - *(ymom_edge_values + (k3_2 + 2));
  # 736 "swb2_domain_ext.c"
  if (extrapolate_velocity_second_order == 1)
  {
   # 738 "swb2_domain_ext.c"
   *(xmom_centroid_values + k_4) = *(xmom_centroid_store + k_4);
   # 739 "swb2_domain_ext.c"
   *(ymom_centroid_values + k_4) = *(ymom_centroid_store + k_4);
   # 742 "swb2_domain_ext.c"
   hmpprt::s32 i_6;
   # 742 "swb2_domain_ext.c"
   # 742 "swb2_domain_ext.c"
   for (i_6 = 0 ; i_6 <= 2 ; i_6 = i_6 + 1)
   {
    # 743 "swb2_domain_ext.c"
    de[i_6] = fmax(*(stage_edge_values + (k3_2 + i_6)) - *(elevation_edge_values + (k3_2 + i_6)), (double) 0.0);
    # 744 "swb2_domain_ext.c"
    *(xmom_edge_values + (k3_2 + i_6)) = *(xmom_edge_values + (k3_2 + i_6)) * de[i_6];
    # 745 "swb2_domain_ext.c"
    *(ymom_edge_values + (k3_2 + i_6)) = *(ymom_edge_values + (k3_2 + i_6)) * de[i_6];
   }
   # 749 "swb2_domain_ext.c"
   # 749 "swb2_domain_ext.c"
   hmpprt::s32 i_7;
   # 749 "swb2_domain_ext.c"
   # 749 "swb2_domain_ext.c"
   for (i_7 = 0 ; i_7 <= 2 ; i_7 = i_7 + 1)
   {
    # 750 "swb2_domain_ext.c"
    de[i_7] = fmax(*(stage_vertex_values + (k3_2 + i_7)) - *(elevation_vertex_values + (k3_2 + i_7)), (double) 0.0);
    # 751 "swb2_domain_ext.c"
    *(xmom_vertex_values + (k3_2 + i_7)) = *(xmom_vertex_values + (k3_2 + i_7)) * de[i_7];
    # 752 "swb2_domain_ext.c"
    *(ymom_vertex_values + (k3_2 + i_7)) = *(ymom_vertex_values + (k3_2 + i_7)) * de[i_7];
   }
   # 176 "swb2_domain_ext.c"
  }
 }
 # 176 "swb2_domain_ext.c"
}
#endif // __CUDACC__



# 176 "swb2_domain_ext.c"

#ifndef __CUDACC__
extern "C" CDLT_API  void extrapolate_second_order_edge_sw(hmpprt::s32 number_of_elements, hmpprt::s32 optimise_dry_cells, hmpprt::s32 extrapolate_velocity_second_order, double epsilon, double minimum_allowed_height, double beta_w, double beta_w_dry, double beta_uh, double beta_uh_dry, double beta_vh, double beta_vh_dry, hmpprt::s64* surrogate_neighbours, hmpprt::s64* number_of_boundaries, double* centroid_coordinates, double* stage_centroid_values, double* elevation_centroid_values, double* xmom_centroid_values, double* ymom_centroid_values, double* edge_coordinates, double* stage_edge_values, double* elevation_edge_values, double* xmom_edge_values, double* ymom_edge_values, double* stage_vertex_values, double* xmom_vertex_values, double* ymom_vertex_values, double* elevation_vertex_values, double* stage_centroid_store, double* xmom_centroid_store, double* ymom_centroid_store, double* min_elevation_edgevalue, double* max_elevation_edgevalue, hmpprt::s32* count_wet_neighbours)
{
 # 1 "<preprocessor>"
 (extrapolate_second_order_edge_sw_internal_1(number_of_elements, optimise_dry_cells, extrapolate_velocity_second_order, epsilon, minimum_allowed_height, beta_w, beta_w_dry, beta_uh, beta_uh_dry, beta_vh, beta_vh_dry, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s64> (surrogate_neighbours), number_of_boundaries, centroid_coordinates, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (stage_centroid_values), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (elevation_centroid_values), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (xmom_centroid_values), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (ymom_centroid_values), edge_coordinates, stage_edge_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (elevation_edge_values), xmom_edge_values, ymom_edge_values, stage_vertex_values, xmom_vertex_values, ymom_vertex_values, elevation_vertex_values, stage_centroid_store, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (xmom_centroid_store), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (ymom_centroid_store), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (min_elevation_edgevalue), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (max_elevation_edgevalue), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s32> (count_wet_neighbours)));
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
      extrapolate_second_order_edge_sw_loop1D_1 = new hmpprt::CUDAGrid(hmpprt_module, "extrapolate_second_order_edge_sw_loop1D_1");
      extrapolate_second_order_edge_sw_loop1D_2 = new hmpprt::CUDAGrid(hmpprt_module, "extrapolate_second_order_edge_sw_loop1D_2");

    }
    hmpprt::Context::getInstance()->getGrouplet()->setTarget(hmpprt::CUDA);
    hmpprt::Context::getInstance()->getGrouplet()->addSignature("_find_qmin_and_qmax", "prototype _find_qmin_and_qmax(dq0: double, dq1: double, dq2: double, qmin: ^host double, qmax: ^host double) : s32");
    hmpprt::Context::getInstance()->getGrouplet()->addSignature("_limit_gradient", "prototype _limit_gradient(dqv: ^host double, qmin: double, qmax: double, beta_w: double) : s32");
    hmpprt::Context::getInstance()->getGrouplet()->addSignature("extrapolate_second_order_edge_sw", "prototype extrapolate_second_order_edge_sw(number_of_elements: s32, optimise_dry_cells: s32, extrapolate_velocity_second_order: s32, epsilon: double, minimum_allowed_height: double, beta_w: double, beta_w_dry: double, beta_uh: double, beta_uh_dry: double, beta_vh: double, beta_vh_dry: double, surrogate_neighbours: ^cudaglob s64, number_of_boundaries: ^host s64, centroid_coordinates: ^host double, stage_centroid_values: ^cudaglob double, elevation_centroid_values: ^cudaglob double, xmom_centroid_values: ^cudaglob double, ymom_centroid_values: ^cudaglob double, edge_coordinates: ^host double, stage_edge_values: ^host double, elevation_edge_values: ^cudaglob double, xmom_edge_values: ^host double, ymom_edge_values: ^host double, stage_vertex_values: ^host double, xmom_vertex_values: ^host double, ymom_vertex_values: ^host double, elevation_vertex_values: ^host double, stage_centroid_store: unused ^host double, xmom_centroid_store: ^cudaglob double, ymom_centroid_store: ^cudaglob double, min_elevation_edgevalue: ^cudaglob double, max_elevation_edgevalue: ^cudaglob double, count_wet_neighbours: ^cudaglob s32)");

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
      delete extrapolate_second_order_edge_sw_loop1D_1;
      delete extrapolate_second_order_edge_sw_loop1D_2;

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

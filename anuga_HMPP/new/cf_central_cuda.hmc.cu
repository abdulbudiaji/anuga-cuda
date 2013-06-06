
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

# 349 "compute_fluxes.c"

#ifndef __CUDACC__
extern "C" CDLT_API  void compute_fluxes_central_structure_CUDA(hmpprt::s32 N, hmpprt::s32 N3, hmpprt::s32 N6_21, hmpprt::s32 N2_2, double* timestep, hmpprt::s64* neighbours_2, hmpprt::s64* neighbour_edges, double* normals_21, double* edgelengths_2, double* radii, double* areas, hmpprt::s64* tri_full_flag_2, double* stage_edge_values, double* xmom_edge_values_21, double* ymom_edge_values_21, double* bed_edge_values_21, double* stage_boundary_values_2, double* xmom_boundary_values, double* ymom_boundary_values_2, double* stage_explicit_update, double* xmom_explicit_update, double* ymom_explicit_update_21, double* max_speed_array, double evolve_max_timestep, double g_21, double epsilon_2, double h0_2, double limiting_threshold_41, hmpprt::s32 optimise_dry_cells_2)
;
#endif // __CUDACC__



# 349 "compute_fluxes.c"

#ifndef __CUDACC__
void compute_fluxes_central_structure_CUDA_internal_1(hmpprt::s32 N, hmpprt::s32 N3, hmpprt::s32 N6_2, hmpprt::s32 N2_21, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  timestep, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s64>  neighbours_21, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s64>  neighbour_edges, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  normals_2, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  edgelengths_21, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  radii, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  areas, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s64>  tri_full_flag_21, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  stage_edge_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  xmom_edge_values_2, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  ymom_edge_values_2, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  bed_edge_values_2, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  stage_boundary_values_21, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  xmom_boundary_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  ymom_boundary_values_21, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  stage_explicit_update, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  xmom_explicit_update, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  ymom_explicit_update_2, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  max_speed_array, double evolve_max_timestep, double g_2, double epsilon_21, double h0_21, double limiting_threshold_4, hmpprt::s32 optimise_dry_cells_21)
;
#endif // __CUDACC__



# 124 "compute_fluxes.c"

#ifndef __CUDACC__
static hmpprt::CUDAGrid * compute_fluxes_central_structure_CUDA_loop1D_1 = 0;
#else

extern "C" __global__ void compute_fluxes_central_structure_CUDA_loop1D_1(hmpprt::s32 N_2, double* timestep_2, hmpprt::s64* neighbours, hmpprt::s64* neighbour_edges_2, double* normals, double* edgelengths, double* radii_2, double* areas_2, hmpprt::s64* tri_full_flag, double* stage_edge_values_2, double* xmom_edge_values, double* ymom_edge_values, double* bed_edge_values, double* stage_boundary_values, double* xmom_boundary_values_2, double* ymom_boundary_values, double* stage_explicit_update_2, double* xmom_explicit_update_2, double* ymom_explicit_update, double* max_speed_array_2, double evolve_max_timestep_2, double g_3, double epsilon_4, double h0_4, double limiting_threshold_2, hmpprt::s32 optimise_dry_cells);
#endif // __CUDACC__




# 124 "compute_fluxes.c"

#ifdef __CUDACC__
__device__ hmpprt::s32 _flux_function_central_cudalocal_1(double* q_left, double* q_right, double z_left, double z_right, double n1, double n2, double epsilon, double h0, double limiting_threshold, double g, double* edgeflux, double* max_speed)
;
#endif // __CUDACC__



# 95 "compute_fluxes.c"

#ifdef __CUDACC__
__device__ double _compute_speed_cudalocal_1(double* uh, double* h, double epsilon_32, double h0_31, double limiting_threshold_3)
;
#endif // __CUDACC__



# 124 "compute_fluxes.c"

#ifndef __CUDACC__
extern "C" CDLT_API  hmpprt::s32 _flux_function_central(double* q_left, double* q_right, double z_left, double z_right, double n1, double n2, double epsilon, double h0, double limiting_threshold, double g, double* edgeflux, double* max_speed)
;
#endif // __CUDACC__



# 124 "compute_fluxes.c"

#ifndef __CUDACC__
hmpprt::s32 _flux_function_central_internal_1(double* q_left, double* q_right, double z_left, double z_right, double n1, double n2, double epsilon, double h0, double limiting_threshold, double g, double* edgeflux, double* max_speed)
;
#endif // __CUDACC__



# 95 "compute_fluxes.c"

#ifndef __CUDACC__
extern "C" CDLT_API  double _compute_speed(double* uh, double* h, double epsilon_31, double h0_32, double limiting_threshold_32)
;
#endif // __CUDACC__



# 95 "compute_fluxes.c"

#ifndef __CUDACC__
double _compute_speed_internal_1(double* uh, double* h, double epsilon_3, double h0_3, double limiting_threshold_31)
;
#endif // __CUDACC__



# 95 "compute_fluxes.c"

#ifndef __CUDACC__
double _compute_speed_internal_1(double* uh, double* h, double epsilon_3, double h0_3, double limiting_threshold_31)
{
 # 101 "compute_fluxes.c"
 double u;
 # 103 "compute_fluxes.c"
 if (*h < limiting_threshold_31)
 {
  # 105 "compute_fluxes.c"
  if (*h < epsilon_3)
  {
   # 106 "compute_fluxes.c"
   *h = (double) 0.0;
   # 107 "compute_fluxes.c"
   u = (double) 0.0;
  }
  else
  {
   # 109 "compute_fluxes.c"
   u = *uh / (*h + h0_3 / *h);
  }
  # 114 "compute_fluxes.c"
  *uh = u * *h;
 }
 else
 {
  # 117 "compute_fluxes.c"
  u = *uh / *h;
 }
 # 95 "compute_fluxes.c"
 return u;
}
#endif // __CUDACC__



# 95 "compute_fluxes.c"

#ifndef __CUDACC__
extern "C" CDLT_API  double _compute_speed(double* uh, double* h, double epsilon_31, double h0_32, double limiting_threshold_32)
{
 # 124 "compute_fluxes.c"
 (_compute_speed_internal_1(uh, h, epsilon_31, h0_32, limiting_threshold_32));
}
#endif // __CUDACC__



# 124 "compute_fluxes.c"

#ifndef __CUDACC__
hmpprt::s32 _flux_function_central_internal_1(double* q_left, double* q_right, double z_left, double z_right, double n1, double n2, double epsilon, double h0, double limiting_threshold, double g, double* edgeflux, double* max_speed)
{
 # 160 "compute_fluxes.c"
 double w_left;
 # 160 "compute_fluxes.c"
 double h_left;
 # 160 "compute_fluxes.c"
 double uh_left;
 # 160 "compute_fluxes.c"
 double vh_left;
 # 160 "compute_fluxes.c"
 double u_left;
 # 161 "compute_fluxes.c"
 double w_right;
 # 161 "compute_fluxes.c"
 double h_right;
 # 161 "compute_fluxes.c"
 double uh_right;
 # 161 "compute_fluxes.c"
 double vh_right;
 # 161 "compute_fluxes.c"
 double u_right;
 # 162 "compute_fluxes.c"
 double s_min;
 # 162 "compute_fluxes.c"
 double s_max;
 # 162 "compute_fluxes.c"
 double soundspeed_left;
 # 162 "compute_fluxes.c"
 double soundspeed_right;
 # 163 "compute_fluxes.c"
 double denom;
 # 163 "compute_fluxes.c"
 double z;
 # 167 "compute_fluxes.c"
 double q_left_rotated[3uLL];
 # 167 "compute_fluxes.c"
 double q_right_rotated[3uLL];
 # 167 "compute_fluxes.c"
 double flux_right[3uLL];
 # 167 "compute_fluxes.c"
 double flux_left[3uLL];
 # 171 "compute_fluxes.c"
 hmpprt::s32 rvalue_1;
 # 171 "compute_fluxes.c"
 q_left_rotated[0] = *q_left;
 # 172 "compute_fluxes.c"
 q_right_rotated[0] = *q_right;
 # 173 "compute_fluxes.c"
 q_left_rotated[1] = *(q_left + 1);
 # 174 "compute_fluxes.c"
 q_right_rotated[1] = *(q_right + 1);
 # 175 "compute_fluxes.c"
 q_left_rotated[2] = *(q_left + 2);
 # 176 "compute_fluxes.c"
 q_right_rotated[2] = *(q_right + 2);
 # 180 "compute_fluxes.c"
 q_left_rotated[1] = n1 * *(q_left + 1) + n2 * *(q_left + 2);
 # 181 "compute_fluxes.c"
 q_left_rotated[2] =  - n2 * *(q_left + 1) + n1 * *(q_left + 2);
 # 184 "compute_fluxes.c"
 q_right_rotated[1] = n1 * *(q_right + 1) + n2 * *(q_right + 2);
 # 185 "compute_fluxes.c"
 q_right_rotated[2] =  - n2 * *(q_right + 1) + n1 * *(q_right + 2);
 # 205 "compute_fluxes.c"
 if (fabs(z_left - z_right) > (double) 1.000000000000000036e-10)
 {
  # 209 "compute_fluxes.c"
  rvalue_1 = 0;
  # 209 "compute_fluxes.c"
  goto endf_11;
 }
 # 209 "compute_fluxes.c"
 z = (double) 0.5 * (z_left + z_right);
 # 212 "compute_fluxes.c"
 w_left = q_left_rotated[0];
 # 213 "compute_fluxes.c"
 h_left = w_left - z;
 # 214 "compute_fluxes.c"
 uh_left = q_left_rotated[1];
 # 216 "compute_fluxes.c"
 u_left = (_compute_speed_internal_1(&uh_left, &h_left, epsilon, h0, limiting_threshold));
 # 217 "compute_fluxes.c"
 w_right = q_right_rotated[0];
 # 218 "compute_fluxes.c"
 h_right = w_right - z;
 # 219 "compute_fluxes.c"
 uh_right = q_right_rotated[1];
 # 221 "compute_fluxes.c"
 u_right = (_compute_speed_internal_1(&uh_right, &h_right, epsilon, h0, limiting_threshold));
 # 223 "compute_fluxes.c"
 vh_left = q_left_rotated[2];
 # 224 "compute_fluxes.c"
 vh_right = q_right_rotated[2];
 # 230 "compute_fluxes.c"
 (_compute_speed_internal_1(&vh_left, &h_left, epsilon, h0, limiting_threshold));
 # 232 "compute_fluxes.c"
 (_compute_speed_internal_1(&vh_right, &h_right, epsilon, h0, limiting_threshold));
 # 235 "compute_fluxes.c"
 soundspeed_left = sqrt(g * h_left);
 # 236 "compute_fluxes.c"
 soundspeed_right = sqrt(g * h_right);
 # 254 "compute_fluxes.c"
 s_max = fmax(u_left + soundspeed_left, u_right + soundspeed_right);
 # 255 "compute_fluxes.c"
 if (s_max < (double) 0.0)
 {
  # 256 "compute_fluxes.c"
  s_max = (double) 0.0;
 }
 # 259 "compute_fluxes.c"
 s_min = fmin(u_left - soundspeed_left, u_right - soundspeed_right);
 # 260 "compute_fluxes.c"
 if (s_min > (double) 0.0)
 {
  # 261 "compute_fluxes.c"
  s_min = (double) 0.0;
 }
 # 265 "compute_fluxes.c"
 flux_left[0] = u_left * h_left;
 # 266 "compute_fluxes.c"
 flux_left[1] = u_left * uh_left + (double) 0.5 * g * h_left * h_left;
 # 267 "compute_fluxes.c"
 flux_left[2] = u_left * vh_left;
 # 269 "compute_fluxes.c"
 flux_right[0] = u_right * h_right;
 # 270 "compute_fluxes.c"
 flux_right[1] = u_right * uh_right + (double) 0.5 * g * h_right * h_right;
 # 271 "compute_fluxes.c"
 flux_right[2] = u_right * vh_right;
 # 274 "compute_fluxes.c"
 denom = s_max - s_min;
 # 275 "compute_fluxes.c"
 if (denom < epsilon)
 {
  # 277 "compute_fluxes.c"
  *edgeflux = (double) 0.0;
  # 278 "compute_fluxes.c"
  *(edgeflux + 1) = (double) 0.0;
  # 279 "compute_fluxes.c"
  *(edgeflux + 2) = (double) 0.0;
  # 280 "compute_fluxes.c"
  *max_speed = (double) 0.0;
 }
 else
 {
  # 283 "compute_fluxes.c"
  double inverse_denominator_11;
  # 283 "compute_fluxes.c"
  double temp_11;
  # 283 "compute_fluxes.c"
  inverse_denominator_11 = (double) 1.0 / denom;
  # 291 "compute_fluxes.c"
  *edgeflux = s_max * flux_left[0] - s_min * flux_right[0];
  # 292 "compute_fluxes.c"
  *edgeflux = *edgeflux + s_max * s_min * (q_right_rotated[0] - q_left_rotated[0]);
  # 293 "compute_fluxes.c"
  *edgeflux = *edgeflux * inverse_denominator_11;
  # 295 "compute_fluxes.c"
  *(edgeflux + 1) = s_max * flux_left[1] - s_min * flux_right[1];
  # 296 "compute_fluxes.c"
  *(edgeflux + 1) = *(edgeflux + 1) + s_max * s_min * (q_right_rotated[1] - q_left_rotated[1]);
  # 297 "compute_fluxes.c"
  *(edgeflux + 1) = *(edgeflux + 1) * inverse_denominator_11;
  # 299 "compute_fluxes.c"
  *(edgeflux + 2) = s_max * flux_left[2] - s_min * flux_right[2];
  # 300 "compute_fluxes.c"
  *(edgeflux + 2) = *(edgeflux + 2) + s_max * s_min * (q_right_rotated[2] - q_left_rotated[2]);
  # 301 "compute_fluxes.c"
  *(edgeflux + 2) = *(edgeflux + 2) * inverse_denominator_11;
  # 317 "compute_fluxes.c"
  *max_speed = fmax(fabs(s_max), fabs(s_min));
  # 322 "compute_fluxes.c"
  temp_11 = *(edgeflux + 1);
  # 323 "compute_fluxes.c"
  *(edgeflux + 1) = n1 * temp_11 - n2 * *(edgeflux + 2);
  # 324 "compute_fluxes.c"
  *(edgeflux + 2) = n2 * temp_11 + n1 * *(edgeflux + 2);
 }
 # 124 "compute_fluxes.c"
 rvalue_1 = 0;
 # 124 "compute_fluxes.c"
 endf_11:;
 # 124 "compute_fluxes.c"
 return rvalue_1;
}
#endif // __CUDACC__



# 124 "compute_fluxes.c"

#ifndef __CUDACC__
extern "C" CDLT_API  hmpprt::s32 _flux_function_central(double* q_left, double* q_right, double z_left, double z_right, double n1, double n2, double epsilon, double h0, double limiting_threshold, double g, double* edgeflux, double* max_speed)
{
 # 95 "compute_fluxes.c"
 (_flux_function_central_internal_1(q_left, q_right, z_left, z_right, n1, n2, epsilon, h0, limiting_threshold, g, edgeflux, max_speed));
}
#endif // __CUDACC__



# 95 "compute_fluxes.c"

#ifdef __CUDACC__
__device__ double _compute_speed_cudalocal_1(double* uh, double* h, double epsilon_32, double h0_31, double limiting_threshold_3)
{
 # 101 "compute_fluxes.c"
 double u;
 # 103 "compute_fluxes.c"
 if (*h < limiting_threshold_3)
 {
  # 105 "compute_fluxes.c"
  if (*h < epsilon_32)
  {
   # 106 "compute_fluxes.c"
   *h = (double) 0.0;
   # 107 "compute_fluxes.c"
   u = (double) 0.0;
  }
  else
  {
   # 109 "compute_fluxes.c"
   u = *uh / (*h + h0_31 / *h);
  }
  # 114 "compute_fluxes.c"
  *uh = u * *h;
 }
 else
 {
  # 117 "compute_fluxes.c"
  u = *uh / *h;
 }
 # 124 "compute_fluxes.c"
 return u;
}
#endif // __CUDACC__



# 124 "compute_fluxes.c"

#ifdef __CUDACC__
__device__ hmpprt::s32 _flux_function_central_cudalocal_1(double* q_left, double* q_right, double z_left, double z_right, double n1, double n2, double epsilon, double h0, double limiting_threshold, double g, double* edgeflux, double* max_speed)
{
 # 160 "compute_fluxes.c"
 double w_left;
 # 160 "compute_fluxes.c"
 double h_left;
 # 160 "compute_fluxes.c"
 double uh_left;
 # 160 "compute_fluxes.c"
 double vh_left;
 # 160 "compute_fluxes.c"
 double u_left;
 # 161 "compute_fluxes.c"
 double w_right;
 # 161 "compute_fluxes.c"
 double h_right;
 # 161 "compute_fluxes.c"
 double uh_right;
 # 161 "compute_fluxes.c"
 double vh_right;
 # 161 "compute_fluxes.c"
 double u_right;
 # 162 "compute_fluxes.c"
 double s_min;
 # 162 "compute_fluxes.c"
 double s_max;
 # 162 "compute_fluxes.c"
 double soundspeed_left;
 # 162 "compute_fluxes.c"
 double soundspeed_right;
 # 163 "compute_fluxes.c"
 double denom;
 # 163 "compute_fluxes.c"
 double z;
 # 167 "compute_fluxes.c"
 double q_left_rotated[3uLL];
 # 167 "compute_fluxes.c"
 double q_right_rotated[3uLL];
 # 167 "compute_fluxes.c"
 double flux_right[3uLL];
 # 167 "compute_fluxes.c"
 double flux_left[3uLL];
 # 171 "compute_fluxes.c"
 hmpprt::s32 rvalue_11;
 # 171 "compute_fluxes.c"
 q_left_rotated[0] = *q_left;
 # 172 "compute_fluxes.c"
 q_right_rotated[0] = *q_right;
 # 173 "compute_fluxes.c"
 q_left_rotated[1] = *(q_left + 1);
 # 174 "compute_fluxes.c"
 q_right_rotated[1] = *(q_right + 1);
 # 175 "compute_fluxes.c"
 q_left_rotated[2] = *(q_left + 2);
 # 176 "compute_fluxes.c"
 q_right_rotated[2] = *(q_right + 2);
 # 180 "compute_fluxes.c"
 q_left_rotated[1] = n1 * *(q_left + 1) + n2 * *(q_left + 2);
 # 181 "compute_fluxes.c"
 q_left_rotated[2] =  - n2 * *(q_left + 1) + n1 * *(q_left + 2);
 # 184 "compute_fluxes.c"
 q_right_rotated[1] = n1 * *(q_right + 1) + n2 * *(q_right + 2);
 # 185 "compute_fluxes.c"
 q_right_rotated[2] =  - n2 * *(q_right + 1) + n1 * *(q_right + 2);
 # 205 "compute_fluxes.c"
 if (fabs(z_left - z_right) > (double) 1.000000000000000036e-10)
 {
  # 209 "compute_fluxes.c"
  rvalue_11 = 0;
  # 209 "compute_fluxes.c"
  goto endf_1;
 }
 # 209 "compute_fluxes.c"
 z = (double) 0.5 * (z_left + z_right);
 # 212 "compute_fluxes.c"
 w_left = q_left_rotated[0];
 # 213 "compute_fluxes.c"
 h_left = w_left - z;
 # 214 "compute_fluxes.c"
 uh_left = q_left_rotated[1];
 # 216 "compute_fluxes.c"
 u_left = (_compute_speed_cudalocal_1(&uh_left, &h_left, epsilon, h0, limiting_threshold));
 # 217 "compute_fluxes.c"
 w_right = q_right_rotated[0];
 # 218 "compute_fluxes.c"
 h_right = w_right - z;
 # 219 "compute_fluxes.c"
 uh_right = q_right_rotated[1];
 # 221 "compute_fluxes.c"
 u_right = (_compute_speed_cudalocal_1(&uh_right, &h_right, epsilon, h0, limiting_threshold));
 # 223 "compute_fluxes.c"
 vh_left = q_left_rotated[2];
 # 224 "compute_fluxes.c"
 vh_right = q_right_rotated[2];
 # 230 "compute_fluxes.c"
 (_compute_speed_cudalocal_1(&vh_left, &h_left, epsilon, h0, limiting_threshold));
 # 232 "compute_fluxes.c"
 (_compute_speed_cudalocal_1(&vh_right, &h_right, epsilon, h0, limiting_threshold));
 # 235 "compute_fluxes.c"
 soundspeed_left = sqrt(g * h_left);
 # 236 "compute_fluxes.c"
 soundspeed_right = sqrt(g * h_right);
 # 254 "compute_fluxes.c"
 s_max = fmax(u_left + soundspeed_left, u_right + soundspeed_right);
 # 255 "compute_fluxes.c"
 if (s_max < (double) 0.0)
 {
  # 256 "compute_fluxes.c"
  s_max = (double) 0.0;
 }
 # 259 "compute_fluxes.c"
 s_min = fmin(u_left - soundspeed_left, u_right - soundspeed_right);
 # 260 "compute_fluxes.c"
 if (s_min > (double) 0.0)
 {
  # 261 "compute_fluxes.c"
  s_min = (double) 0.0;
 }
 # 265 "compute_fluxes.c"
 flux_left[0] = u_left * h_left;
 # 266 "compute_fluxes.c"
 flux_left[1] = u_left * uh_left + (double) 0.5 * g * h_left * h_left;
 # 267 "compute_fluxes.c"
 flux_left[2] = u_left * vh_left;
 # 269 "compute_fluxes.c"
 flux_right[0] = u_right * h_right;
 # 270 "compute_fluxes.c"
 flux_right[1] = u_right * uh_right + (double) 0.5 * g * h_right * h_right;
 # 271 "compute_fluxes.c"
 flux_right[2] = u_right * vh_right;
 # 274 "compute_fluxes.c"
 denom = s_max - s_min;
 # 275 "compute_fluxes.c"
 if (denom < epsilon)
 {
  # 277 "compute_fluxes.c"
  *edgeflux = (double) 0.0;
  # 278 "compute_fluxes.c"
  *(edgeflux + 1) = (double) 0.0;
  # 279 "compute_fluxes.c"
  *(edgeflux + 2) = (double) 0.0;
  # 280 "compute_fluxes.c"
  *max_speed = (double) 0.0;
 }
 else
 {
  # 283 "compute_fluxes.c"
  double inverse_denominator_1;
  # 283 "compute_fluxes.c"
  double temp_1;
  # 283 "compute_fluxes.c"
  inverse_denominator_1 = (double) 1.0 / denom;
  # 291 "compute_fluxes.c"
  *edgeflux = s_max * flux_left[0] - s_min * flux_right[0];
  # 292 "compute_fluxes.c"
  *edgeflux = *edgeflux + s_max * s_min * (q_right_rotated[0] - q_left_rotated[0]);
  # 293 "compute_fluxes.c"
  *edgeflux = *edgeflux * inverse_denominator_1;
  # 295 "compute_fluxes.c"
  *(edgeflux + 1) = s_max * flux_left[1] - s_min * flux_right[1];
  # 296 "compute_fluxes.c"
  *(edgeflux + 1) = *(edgeflux + 1) + s_max * s_min * (q_right_rotated[1] - q_left_rotated[1]);
  # 297 "compute_fluxes.c"
  *(edgeflux + 1) = *(edgeflux + 1) * inverse_denominator_1;
  # 299 "compute_fluxes.c"
  *(edgeflux + 2) = s_max * flux_left[2] - s_min * flux_right[2];
  # 300 "compute_fluxes.c"
  *(edgeflux + 2) = *(edgeflux + 2) + s_max * s_min * (q_right_rotated[2] - q_left_rotated[2]);
  # 301 "compute_fluxes.c"
  *(edgeflux + 2) = *(edgeflux + 2) * inverse_denominator_1;
  # 317 "compute_fluxes.c"
  *max_speed = fmax(fabs(s_max), fabs(s_min));
  # 322 "compute_fluxes.c"
  temp_1 = *(edgeflux + 1);
  # 323 "compute_fluxes.c"
  *(edgeflux + 1) = n1 * temp_1 - n2 * *(edgeflux + 2);
  # 324 "compute_fluxes.c"
  *(edgeflux + 2) = n2 * temp_1 + n1 * *(edgeflux + 2);
 }
 # 349 "compute_fluxes.c"
 rvalue_11 = 0;
 # 349 "compute_fluxes.c"
 endf_1:;
 # 349 "compute_fluxes.c"
 return rvalue_11;
}
#endif // __CUDACC__



# 349 "compute_fluxes.c"

#ifdef __CUDACC__

extern "C" __global__ void compute_fluxes_central_structure_CUDA_loop1D_1(hmpprt::s32 N_2, double* timestep_2, hmpprt::s64* neighbours, hmpprt::s64* neighbour_edges_2, double* normals, double* edgelengths, double* radii_2, double* areas_2, hmpprt::s64* tri_full_flag, double* stage_edge_values_2, double* xmom_edge_values, double* ymom_edge_values, double* bed_edge_values, double* stage_boundary_values, double* xmom_boundary_values_2, double* ymom_boundary_values, double* stage_explicit_update_2, double* xmom_explicit_update_2, double* ymom_explicit_update, double* max_speed_array_2, double evolve_max_timestep_2, double g_3, double epsilon_4, double h0_4, double limiting_threshold_2, hmpprt::s32 optimise_dry_cells)
{
 # 386 "compute_fluxes.c"
 double max_speed_1;
 # 386 "compute_fluxes.c"
 double max_speed_total_1;
 # 389 "compute_fluxes.c"
 double ql_1[3uLL];
 # 389 "compute_fluxes.c"
 double qr[3uLL];
 # 389 "compute_fluxes.c"
 double edgeflux_1[3uLL];
 # 382 "compute_fluxes.c"
 double inv_area_1;
 # 382 "compute_fluxes.c"
 hmpprt::s32 k_1;
 # 408 "compute_fluxes.c"
 k_1 = (hmpprt::gr_atidf());
 # 408 "compute_fluxes.c"
 if (k_1 > N_2 - 1)
 {
  # 408 "compute_fluxes.c"
  goto __hmppcg_label_1;
 }
 # 408 "compute_fluxes.c"
 *(timestep_2 + k_1) = evolve_max_timestep_2;
 # 411 "compute_fluxes.c"
 hmpprt::s32 i_2;
 # 411 "compute_fluxes.c"
 # 412 "compute_fluxes.c"
 for (i_2 = 0 ; i_2 <= 2 ; i_2 = i_2 + 1)
 {
  # 413 "compute_fluxes.c"
  hmpprt::s32 ki_1;
  # 413 "compute_fluxes.c"
  hmpprt::s32 n_1;
  # 413 "compute_fluxes.c"
  double zl_1;
  # 413 "compute_fluxes.c"
  double zr_1;
  # 413 "compute_fluxes.c"
  hmpprt::s32 ki2_1;
  # 413 "compute_fluxes.c"
  double length_1;
  # 413 "compute_fluxes.c"
  ki_1 = k_1 * 3 + i_2;
  # 415 "compute_fluxes.c"
  n_1 = (hmpprt::s32 ) (*(neighbours + ki_1));
  # 418 "compute_fluxes.c"
  ql_1[0] = *(stage_edge_values_2 + ki_1);
  # 419 "compute_fluxes.c"
  ql_1[1] = *(xmom_edge_values + ki_1);
  # 420 "compute_fluxes.c"
  ql_1[2] = *(ymom_edge_values + ki_1);
  # 426 "compute_fluxes.c"
  zl_1 = *(bed_edge_values + ki_1);
  # 428 "compute_fluxes.c"
  if (n_1 < 0)
  {
   # 429 "compute_fluxes.c"
   hmpprt::s32 m_1;
   # 429 "compute_fluxes.c"
   m_1 =  - (n_1 + 1);
   # 432 "compute_fluxes.c"
   qr[0] = *(stage_boundary_values + m_1);
   # 433 "compute_fluxes.c"
   qr[1] = *(xmom_boundary_values_2 + m_1);
   # 434 "compute_fluxes.c"
   qr[2] = *(ymom_boundary_values + m_1);
   # 440 "compute_fluxes.c"
   zr_1 = zl_1;
  }
  else
  {
   # 442 "compute_fluxes.c"
   hmpprt::s32 m_2;
   # 442 "compute_fluxes.c"
   hmpprt::s32 nm_1;
   # 442 "compute_fluxes.c"
   m_2 = (hmpprt::s32 ) (*(neighbour_edges_2 + ki_1));
   # 443 "compute_fluxes.c"
   nm_1 = n_1 * 3 + m_2;
   # 446 "compute_fluxes.c"
   qr[0] = *(stage_edge_values_2 + nm_1);
   # 447 "compute_fluxes.c"
   qr[1] = *(xmom_edge_values + nm_1);
   # 448 "compute_fluxes.c"
   qr[2] = *(ymom_edge_values + nm_1);
   # 454 "compute_fluxes.c"
   zr_1 = *(bed_edge_values + nm_1);
  }
  # 457 "compute_fluxes.c"
  if (optimise_dry_cells)
  {
   # 459 "compute_fluxes.c"
   if (fabs(ql_1[0] - zl_1) < epsilon_4 && fabs(qr[0] - zr_1) < epsilon_4)
   {
    # 465 "compute_fluxes.c"
    max_speed_1 = (double) 0.0;
    # 466 "compute_fluxes.c"
    continue;
   }
  }
  # 470 "compute_fluxes.c"
  ki2_1 = 2 * ki_1;
  # 475 "compute_fluxes.c"
  (_flux_function_central_cudalocal_1(&ql_1[0], &qr[0], zl_1, zr_1, *(normals + ki2_1), *(normals + (ki2_1 + 1)), epsilon_4, h0_4, limiting_threshold_2, g_3, &edgeflux_1[0], &max_speed_1));
  # 477 "compute_fluxes.c"
  length_1 = *(edgelengths + ki_1);
  # 478 "compute_fluxes.c"
  edgeflux_1[0] = edgeflux_1[0] * length_1;
  # 479 "compute_fluxes.c"
  edgeflux_1[1] = edgeflux_1[1] * length_1;
  # 480 "compute_fluxes.c"
  edgeflux_1[2] = edgeflux_1[2] * length_1;
  # 483 "compute_fluxes.c"
  *(stage_explicit_update_2 + k_1) = *(stage_explicit_update_2 + k_1) - edgeflux_1[0];
  # 484 "compute_fluxes.c"
  *(xmom_explicit_update_2 + k_1) = *(xmom_explicit_update_2 + k_1) - edgeflux_1[1];
  # 485 "compute_fluxes.c"
  *(ymom_explicit_update + k_1) = *(ymom_explicit_update + k_1) - edgeflux_1[2];
  # 509 "compute_fluxes.c"
  if (*(tri_full_flag + k_1) == 1LL)
  {
   # 511 "compute_fluxes.c"
   if (max_speed_1 > epsilon_4)
   {
    # 512 "compute_fluxes.c"
    *(timestep_2 + k_1) = fmin(*(timestep_2 + k_1), *(radii_2 + k_1) / max_speed_1);
    # 513 "compute_fluxes.c"
    if (n_1 >= 0)
    {
     # 514 "compute_fluxes.c"
     *(timestep_2 + k_1) = fmin(*(timestep_2 + k_1), *(radii_2 + n_1) / max_speed_1);
    }
   }
  }
  # 519 "compute_fluxes.c"
  if (n_1 < 0 || n_1 > k_1)
  {
   # 520 "compute_fluxes.c"
   max_speed_total_1 = fmax(max_speed_total_1, max_speed_1);
  }
 }
 # 524 "compute_fluxes.c"
 # 524 "compute_fluxes.c"
 inv_area_1 = (double) 1.0 / *(areas_2 + k_1);
 # 525 "compute_fluxes.c"
 *(stage_explicit_update_2 + k_1) = *(stage_explicit_update_2 + k_1) * inv_area_1;
 # 526 "compute_fluxes.c"
 *(xmom_explicit_update_2 + k_1) = *(xmom_explicit_update_2 + k_1) * inv_area_1;
 # 527 "compute_fluxes.c"
 *(ymom_explicit_update + k_1) = *(ymom_explicit_update + k_1) * inv_area_1;
 # 529 "compute_fluxes.c"
 *(max_speed_array_2 + k_1) = max_speed_total_1;
 # 349 "compute_fluxes.c"
 __hmppcg_label_1:;
}
#endif // __CUDACC__



# 349 "compute_fluxes.c"

#ifndef __CUDACC__
void compute_fluxes_central_structure_CUDA_internal_1(hmpprt::s32 N, hmpprt::s32 N3, hmpprt::s32 N6_2, hmpprt::s32 N2_21, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  timestep, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s64>  neighbours_21, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s64>  neighbour_edges, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  normals_2, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  edgelengths_21, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  radii, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  areas, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s64>  tri_full_flag_21, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  stage_edge_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  xmom_edge_values_2, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  ymom_edge_values_2, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  bed_edge_values_2, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  stage_boundary_values_21, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  xmom_boundary_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  ymom_boundary_values_21, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  stage_explicit_update, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  xmom_explicit_update, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  ymom_explicit_update_2, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  max_speed_array, double evolve_max_timestep, double g_2, double epsilon_21, double h0_21, double limiting_threshold_4, hmpprt::s32 optimise_dry_cells_21)
{
 # 349 "compute_fluxes.c"
 if (N - 1 >= 0)
 {
  hmpprt::CUDAGridCall __hmppcg_call;
  __hmppcg_call.setSizeX((N - 1) / 128 + 1);
  __hmppcg_call.setSizeY(1);
  __hmppcg_call.setBlockSizeX(32);
  __hmppcg_call.setBlockSizeY(4);
  __hmppcg_call.addLocalParameter((hmpprt::s32) (N), "N_2");
  __hmppcg_call.addLocalParameter(&timestep, 8, "timestep_2");
  __hmppcg_call.addLocalParameter(&neighbours_21, 8, "neighbours");
  __hmppcg_call.addLocalParameter(&neighbour_edges, 8, "neighbour_edges_2");
  __hmppcg_call.addLocalParameter(&normals_2, 8, "normals");
  __hmppcg_call.addLocalParameter(&edgelengths_21, 8, "edgelengths");
  __hmppcg_call.addLocalParameter(&radii, 8, "radii_2");
  __hmppcg_call.addLocalParameter(&areas, 8, "areas_2");
  __hmppcg_call.addLocalParameter(&tri_full_flag_21, 8, "tri_full_flag");
  __hmppcg_call.addLocalParameter(&stage_edge_values, 8, "stage_edge_values_2");
  __hmppcg_call.addLocalParameter(&xmom_edge_values_2, 8, "xmom_edge_values");
  __hmppcg_call.addLocalParameter(&ymom_edge_values_2, 8, "ymom_edge_values");
  __hmppcg_call.addLocalParameter(&bed_edge_values_2, 8, "bed_edge_values");
  __hmppcg_call.addLocalParameter(&stage_boundary_values_21, 8, "stage_boundary_values");
  __hmppcg_call.addLocalParameter(&xmom_boundary_values, 8, "xmom_boundary_values_2");
  __hmppcg_call.addLocalParameter(&ymom_boundary_values_21, 8, "ymom_boundary_values");
  __hmppcg_call.addLocalParameter(&stage_explicit_update, 8, "stage_explicit_update_2");
  __hmppcg_call.addLocalParameter(&xmom_explicit_update, 8, "xmom_explicit_update_2");
  __hmppcg_call.addLocalParameter(&ymom_explicit_update_2, 8, "ymom_explicit_update");
  __hmppcg_call.addLocalParameter(&max_speed_array, 8, "max_speed_array_2");
  __hmppcg_call.addLocalParameter(&evolve_max_timestep, 8, "evolve_max_timestep_2");
  __hmppcg_call.addLocalParameter(&g_2, 8, "g_3");
  __hmppcg_call.addLocalParameter(&epsilon_21, 8, "epsilon_4");
  __hmppcg_call.addLocalParameter(&h0_21, 8, "h0_4");
  __hmppcg_call.addLocalParameter(&limiting_threshold_4, 8, "limiting_threshold_2");
  __hmppcg_call.addLocalParameter((hmpprt::s32) (optimise_dry_cells_21), "optimise_dry_cells");
  __hmppcg_call.launch(compute_fluxes_central_structure_CUDA_loop1D_1, hmpprt::Context::getInstance()->getCUDADevice());
 }
 ;
}
#endif // __CUDACC__



# 349 "compute_fluxes.c"

#ifndef __CUDACC__
extern "C" CDLT_API  void compute_fluxes_central_structure_CUDA(hmpprt::s32 N, hmpprt::s32 N3, hmpprt::s32 N6_21, hmpprt::s32 N2_2, double* timestep, hmpprt::s64* neighbours_2, hmpprt::s64* neighbour_edges, double* normals_21, double* edgelengths_2, double* radii, double* areas, hmpprt::s64* tri_full_flag_2, double* stage_edge_values, double* xmom_edge_values_21, double* ymom_edge_values_21, double* bed_edge_values_21, double* stage_boundary_values_2, double* xmom_boundary_values, double* ymom_boundary_values_2, double* stage_explicit_update, double* xmom_explicit_update, double* ymom_explicit_update_21, double* max_speed_array, double evolve_max_timestep, double g_21, double epsilon_2, double h0_2, double limiting_threshold_41, hmpprt::s32 optimise_dry_cells_2)
{
 # 1 "<preprocessor>"
 (compute_fluxes_central_structure_CUDA_internal_1(N, N3, N6_21, N2_2, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (timestep), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s64> (neighbours_2), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s64> (neighbour_edges), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (normals_21), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (edgelengths_2), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (radii), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (areas), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s64> (tri_full_flag_2), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (stage_edge_values), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (xmom_edge_values_21), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (ymom_edge_values_21), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (bed_edge_values_21), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (stage_boundary_values_2), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (xmom_boundary_values), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (ymom_boundary_values_2), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (stage_explicit_update), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (xmom_explicit_update), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (ymom_explicit_update_21), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (max_speed_array), evolve_max_timestep, g_21, epsilon_2, h0_2, limiting_threshold_41, optimise_dry_cells_2));
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
      compute_fluxes_central_structure_CUDA_loop1D_1 = new hmpprt::CUDAGrid(hmpprt_module, "compute_fluxes_central_structure_CUDA_loop1D_1");

    }
    hmpprt::Context::getInstance()->getGrouplet()->setTarget(hmpprt::CUDA);
    hmpprt::Context::getInstance()->getGrouplet()->addSignature("_compute_speed", "prototype _compute_speed(uh: ^host double, h: ^host double, epsilon: double, h0: double, limiting_threshold: double) : double");
    hmpprt::Context::getInstance()->getGrouplet()->addSignature("_flux_function_central", "prototype _flux_function_central(q_left: ^host double, q_right: ^host double, z_left: double, z_right: double, n1: double, n2: double, epsilon: double, h0: double, limiting_threshold: double, g: double, edgeflux: ^host double, max_speed: ^host double) : s32");
    hmpprt::Context::getInstance()->getGrouplet()->addSignature("compute_fluxes_central_structure_CUDA", "prototype compute_fluxes_central_structure_CUDA(N: s32, N3: s32, N6: s32, N2: s32, timestep: ^cudaglob double, neighbours: ^cudaglob s64, neighbour_edges: ^cudaglob s64, normals: ^cudaglob double, edgelengths: ^cudaglob double, radii: ^cudaglob double, areas: ^cudaglob double, tri_full_flag: ^cudaglob s64, stage_edge_values: ^cudaglob double, xmom_edge_values: ^cudaglob double, ymom_edge_values: ^cudaglob double, bed_edge_values: ^cudaglob double, stage_boundary_values: ^cudaglob double, xmom_boundary_values: ^cudaglob double, ymom_boundary_values: ^cudaglob double, stage_explicit_update: ^cudaglob double, xmom_explicit_update: ^cudaglob double, ymom_explicit_update: ^cudaglob double, max_speed_array: ^cudaglob double, evolve_max_timestep: double, g: double, epsilon: double, h0: double, limiting_threshold: double, optimise_dry_cells: s32)");

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
      delete compute_fluxes_central_structure_CUDA_loop1D_1;

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

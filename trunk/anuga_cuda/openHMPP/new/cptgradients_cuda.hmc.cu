
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

# 137 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"

#ifndef __CUDACC__
extern "C" CDLT_API  void _compute_gradients(hmpprt::s32 N, hmpprt::s32 N2, hmpprt::s32 N3, double* centroids, double* centroid_values, hmpprt::s64* number_of_boundaries, hmpprt::s64* surrogate_neighbours, double* a, double* b)
;
#endif // __CUDACC__



# 137 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"

#ifndef __CUDACC__
void _compute_gradients_internal_1(hmpprt::s32 N, hmpprt::s32 N2, hmpprt::s32 N3, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  centroids, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  centroid_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s64>  number_of_boundaries, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s64>  surrogate_neighbours, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  a, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  b)
;
#endif // __CUDACC__



# 78 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"

#ifndef __CUDACC__
static hmpprt::CUDAGrid * _compute_gradients_loop1D_1 = 0;
#else

extern "C" __global__ void _compute_gradients_loop1D_1(hmpprt::s32 N_2, double* centroids_2, double* centroid_values_2, hmpprt::s64* number_of_boundaries_2, hmpprt::s64* surrogate_neighbours_2, double* a_4, double* b_4);
#endif // __CUDACC__




# 78 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"

#ifdef __CUDACC__
__device__ hmpprt::s32 gradient2_cudalocal_1(double x0_22, double y0_21, double x1_21, double y1_2, double q0_2, double q1_2, double* a_32, double* b_3)
;
#endif // __CUDACC__



# 31 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"

#ifdef __CUDACC__
__device__ hmpprt::s32 gradient_cudalocal_1(double x0_12, double y0_12, double x1_12, double y1_12, double x2_12, double y2_12, double q0_1, double q1_12, double q2_12, double* a_22, double* b_21)
;
#endif // __CUDACC__



# 78 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"

#ifndef __CUDACC__
extern "C" CDLT_API  hmpprt::s32 gradient2(double x0_21, double y0_22, double x1_22, double y1_22, double q0_22, double q1_22, double* a_31, double* b_32)
;
#endif // __CUDACC__



# 78 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"

#ifndef __CUDACC__
hmpprt::s32 gradient2_internal_1(double x0_2, double y0_2, double x1_2, double y1_21, double q0_21, double q1_21, double* a_3, double* b_31)
;
#endif // __CUDACC__



# 31 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"

#ifndef __CUDACC__
extern "C" CDLT_API  hmpprt::s32 gradient(double x0_1, double y0_1, double x1_1, double y1_11, double x2_11, double y2_11, double q0_12, double q1_11, double q2_11, double* a_21, double* b_22)
;
#endif // __CUDACC__



# 31 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"

#ifndef __CUDACC__
hmpprt::s32 gradient_internal_1(double x0_11, double y0_11, double x1_11, double y1_1, double x2_1, double y2_1, double q0_11, double q1_1, double q2_1, double* a_2, double* b_2)
;
#endif // __CUDACC__



# 31 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"

#ifndef __CUDACC__
hmpprt::s32 gradient_internal_1(double x0_11, double y0_11, double x1_11, double y1_1, double x2_1, double y2_1, double q0_11, double q1_1, double q2_1, double* a_2, double* b_2)
{
 # 62 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 double det_1;
 # 64 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 det_1 = (y2_1 - y0_11) * (x1_11 - x0_11) - (y1_1 - y0_11) * (x2_1 - x0_11);
 # 66 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 *a_2 = (y2_1 - y0_11) * (q1_1 - q0_11) - (y1_1 - y0_11) * (q2_1 - q0_11);
 # 67 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 *a_2 = *a_2 / det_1;
 # 69 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 *b_2 = (x1_11 - x0_11) * (q2_1 - q0_11) - (x2_1 - x0_11) * (q1_1 - q0_11);
 # 70 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 *b_2 = *b_2 / det_1;
 # 31 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 return 0;
}
#endif // __CUDACC__



# 31 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"

#ifndef __CUDACC__
extern "C" CDLT_API  hmpprt::s32 gradient(double x0_1, double y0_1, double x1_1, double y1_11, double x2_11, double y2_11, double q0_12, double q1_11, double q2_11, double* a_21, double* b_22)
{
 # 78 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 (gradient_internal_1(x0_1, y0_1, x1_1, y1_11, x2_11, y2_11, q0_12, q1_11, q2_11, a_21, b_22));
}
#endif // __CUDACC__



# 78 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"

#ifndef __CUDACC__
hmpprt::s32 gradient2_internal_1(double x0_2, double y0_2, double x1_2, double y1_21, double q0_21, double q1_21, double* a_3, double* b_31)
{
 # 120 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 double det;
 # 120 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 double xx;
 # 120 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 double yy;
 # 120 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 double qq;
 # 122 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 xx = x1_2 - x0_2;
 # 123 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 yy = y1_21 - y0_2;
 # 124 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 qq = q1_21 - q0_21;
 # 126 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 det = xx * xx + yy * yy;
 # 127 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 *a_3 = xx * qq / det;
 # 128 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 *b_31 = yy * qq / det;
 # 78 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 return 0;
}
#endif // __CUDACC__



# 78 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"

#ifndef __CUDACC__
extern "C" CDLT_API  hmpprt::s32 gradient2(double x0_21, double y0_22, double x1_22, double y1_22, double q0_22, double q1_22, double* a_31, double* b_32)
{
 # 31 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 (gradient2_internal_1(x0_21, y0_22, x1_22, y1_22, q0_22, q1_22, a_31, b_32));
}
#endif // __CUDACC__



# 31 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"

#ifdef __CUDACC__
__device__ hmpprt::s32 gradient_cudalocal_1(double x0_12, double y0_12, double x1_12, double y1_12, double x2_12, double y2_12, double q0_1, double q1_12, double q2_12, double* a_22, double* b_21)
{
 # 62 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 double det_11;
 # 64 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 det_11 = (y2_12 - y0_12) * (x1_12 - x0_12) - (y1_12 - y0_12) * (x2_12 - x0_12);
 # 66 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 *a_22 = (y2_12 - y0_12) * (q1_12 - q0_1) - (y1_12 - y0_12) * (q2_12 - q0_1);
 # 67 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 *a_22 = *a_22 / det_11;
 # 69 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 *b_21 = (x1_12 - x0_12) * (q2_12 - q0_1) - (x2_12 - x0_12) * (q1_12 - q0_1);
 # 70 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 *b_21 = *b_21 / det_11;
 # 78 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 return 0;
}
#endif // __CUDACC__



# 78 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"

#ifdef __CUDACC__
__device__ hmpprt::s32 gradient2_cudalocal_1(double x0_22, double y0_21, double x1_21, double y1_2, double q0_2, double q1_2, double* a_32, double* b_3)
{
 # 120 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 double det;
 # 120 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 double xx;
 # 120 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 double yy;
 # 120 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 double qq;
 # 122 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 xx = x1_21 - x0_22;
 # 123 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 yy = y1_2 - y0_21;
 # 124 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 qq = q1_2 - q0_2;
 # 126 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 det = xx * xx + yy * yy;
 # 127 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 *a_32 = xx * qq / det;
 # 128 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 *b_3 = yy * qq / det;
 # 137 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 return 0;
}
#endif // __CUDACC__



# 137 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"

#ifdef __CUDACC__

extern "C" __global__ void _compute_gradients_loop1D_1(hmpprt::s32 N_2, double* centroids_2, double* centroid_values_2, hmpprt::s64* number_of_boundaries_2, hmpprt::s64* surrogate_neighbours_2, double* a_4, double* b_4)
{
 # 148 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 hmpprt::s32 k3_1;
 # 148 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 hmpprt::s32 k_1;
 # 158 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 k_1 = (hmpprt::gr_atidf());
 # 158 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 if (k_1 > N_2 - 1)
 {
  # 158 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  goto __hmppcg_label_1;
 }
 # 158 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 k3_1 = 3 * k_1;
 # 160 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 if (*(number_of_boundaries_2 + k_1) < 2LL)
 {
  # 166 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  hmpprt::s32 k0_1;
  # 166 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  hmpprt::s32 k1_1;
  # 166 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  hmpprt::s32 k2_1;
  # 166 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  double x0_3;
  # 166 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  double y0_3;
  # 166 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  double x1_3;
  # 166 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  double y1_3;
  # 166 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  double x2_2;
  # 166 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  double y2_2;
  # 166 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  double q0_3;
  # 166 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  double q1_3;
  # 166 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  double q2_2;
  # 166 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  k0_1 = (hmpprt::s32 ) (*(surrogate_neighbours_2 + k3_1));
  # 167 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  k1_1 = (hmpprt::s32 ) (*(surrogate_neighbours_2 + (k3_1 + 1)));
  # 168 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  k2_1 = (hmpprt::s32 ) (*(surrogate_neighbours_2 + (k3_1 + 2)));
  # 175 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  q0_3 = *(centroid_values_2 + k0_1);
  # 176 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  q1_3 = *(centroid_values_2 + k1_1);
  # 177 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  q2_2 = *(centroid_values_2 + k2_1);
  # 179 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  x0_3 = *(centroids_2 + k0_1 * 2);
  # 179 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  y0_3 = *(centroids_2 + (k0_1 * 2 + 1));
  # 180 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  x1_3 = *(centroids_2 + k1_1 * 2);
  # 180 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  y1_3 = *(centroids_2 + (k1_1 * 2 + 1));
  # 181 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  x2_2 = *(centroids_2 + k2_1 * 2);
  # 181 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  y2_2 = *(centroids_2 + (k2_1 * 2 + 1));
  # 184 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  (gradient_cudalocal_1(x0_3, y0_3, x1_3, y1_3, x2_2, y2_2, q0_3, q1_3, q2_2, a_4 + k_1, b_4 + k_1));
 }
 else
 {
  # 186 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
  if (*(number_of_boundaries_2 + k_1) == 2LL)
  {
   # 197 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
   hmpprt::s32 k0_2;
   # 197 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
   double x0_4;
   # 197 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
   double y0_4;
   # 197 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
   double x1_4;
   # 197 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
   double y1_4;
   # 197 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
   double q0_4;
   # 197 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
   double q1_4;
   # 197 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
   k0_2 = (hmpprt::s32 ) (*(surrogate_neighbours_2 + k3_1));
   # 198 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
   if (k0_2 == k_1)
   {
    # 200 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
    k0_2 = (hmpprt::s32 ) (*(surrogate_neighbours_2 + (k3_1 + 1)));
    # 201 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
    if (k0_2 == k_1)
    {
     # 203 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
     k0_2 = (hmpprt::s32 ) (*(surrogate_neighbours_2 + (k3_1 + 2)));
    }
   }
   # 214 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
   q0_4 = *(centroid_values_2 + k0_2);
   # 215 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
   q1_4 = *(centroid_values_2 + k_1);
   # 217 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
   x0_4 = *(centroids_2 + k0_2 * 2);
   # 217 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
   y0_4 = *(centroids_2 + (k0_2 * 2 + 1));
   # 218 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
   x1_4 = *(centroids_2 + k_1 * 2);
   # 218 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
   y1_4 = *(centroids_2 + (k_1 * 2 + 1));
   # 221 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
   (gradient2_cudalocal_1(x0_4, y0_4, x1_4, y1_4, q0_4, q1_4, a_4 + k_1, b_4 + k_1));
  }
 }
 # 137 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 __hmppcg_label_1:;
}
#endif // __CUDACC__



# 137 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"

#ifndef __CUDACC__
void _compute_gradients_internal_1(hmpprt::s32 N, hmpprt::s32 N2, hmpprt::s32 N3, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  centroids, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  centroid_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s64>  number_of_boundaries, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s64>  surrogate_neighbours, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  a, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  b)
{
 # 137 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"
 if (N - 1 >= 0)
 {
  hmpprt::CUDAGridCall __hmppcg_call;
  __hmppcg_call.setSizeX((N - 1) / 128 + 1);
  __hmppcg_call.setSizeY(1);
  __hmppcg_call.setBlockSizeX(32);
  __hmppcg_call.setBlockSizeY(4);
  __hmppcg_call.addLocalParameter((hmpprt::s32) (N), "N_2");
  __hmppcg_call.addLocalParameter(&centroids, 8, "centroids_2");
  __hmppcg_call.addLocalParameter(&centroid_values, 8, "centroid_values_2");
  __hmppcg_call.addLocalParameter(&number_of_boundaries, 8, "number_of_boundaries_2");
  __hmppcg_call.addLocalParameter(&surrogate_neighbours, 8, "surrogate_neighbours_2");
  __hmppcg_call.addLocalParameter(&a, 8, "a_4");
  __hmppcg_call.addLocalParameter(&b, 8, "b_4");
  __hmppcg_call.launch(_compute_gradients_loop1D_1, hmpprt::Context::getInstance()->getCUDADevice());
 }
 ;
}
#endif // __CUDACC__



# 137 "extrapolate_second_order_and_limit_by_vertex_or_edge.c"

#ifndef __CUDACC__
extern "C" CDLT_API  void _compute_gradients(hmpprt::s32 N, hmpprt::s32 N2, hmpprt::s32 N3, double* centroids, double* centroid_values, hmpprt::s64* number_of_boundaries, hmpprt::s64* surrogate_neighbours, double* a, double* b)
{
 # 1 "<preprocessor>"
 (_compute_gradients_internal_1(N, N2, N3, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (centroids), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (centroid_values), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s64> (number_of_boundaries), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s64> (surrogate_neighbours), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (a), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (b)));
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
      _compute_gradients_loop1D_1 = new hmpprt::CUDAGrid(hmpprt_module, "_compute_gradients_loop1D_1");

    }
    hmpprt::Context::getInstance()->getGrouplet()->setTarget(hmpprt::CUDA);
    hmpprt::Context::getInstance()->getGrouplet()->addSignature("gradient", "prototype gradient(x0: double, y0: double, x1: double, y1: double, x2: double, y2: double, q0: double, q1: double, q2: double, a: ^host double, b: ^host double) : s32");
    hmpprt::Context::getInstance()->getGrouplet()->addSignature("gradient2", "prototype gradient2(x0: double, y0: double, x1: double, y1: double, q0: double, q1: double, a: ^host double, b: ^host double) : s32");
    hmpprt::Context::getInstance()->getGrouplet()->addSignature("_compute_gradients", "prototype _compute_gradients(N: s32, N2: s32, N3: s32, centroids: ^cudaglob double, centroid_values: ^cudaglob double, number_of_boundaries: ^cudaglob s64, surrogate_neighbours: ^cudaglob s64, a: ^cudaglob double, b: ^cudaglob double)");

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
      delete _compute_gradients_loop1D_1;

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

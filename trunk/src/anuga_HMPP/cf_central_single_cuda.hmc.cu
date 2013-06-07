
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

# 538 "compute_fluxes.c"

#ifndef __CUDACC__
extern "C" CDLT_API  void compute_fluxes_central_structure_cuda_single(hmpprt::s32 N, hmpprt::s32 N3, hmpprt::s32 N6, hmpprt::s32 N2, double* timestep, hmpprt::s32* neighbours, hmpprt::s32* neighbour_edges, double* normals, double* edgelengths, double* radii, double* areas, hmpprt::s32* tri_full_flag, double* stage_edge_values, double* xmom_edge_values, double* ymom_edge_values, double* bed_edge_values, double* stage_boundary_values, double* xmom_boundary_values, double* ymom_boundary_values, double* stage_explicit_update, double* xmom_explicit_update, double* ymom_explicit_update, double* max_speed_array, double evolve_max_timestep, double g, double epsilon, double h0, double limiting_threshold, hmpprt::s32 optimise_dry_cells)
;
#endif // __CUDACC__



# 538 "compute_fluxes.c"

#ifndef __CUDACC__
void compute_fluxes_central_structure_cuda_single_internal_1(hmpprt::s32 N, hmpprt::s32 N3, hmpprt::s32 N6, hmpprt::s32 N2, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  timestep, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s32>  neighbours, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s32>  neighbour_edges, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  normals, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  edgelengths, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  radii, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  areas, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s32>  tri_full_flag, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  stage_edge_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  xmom_edge_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  ymom_edge_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  bed_edge_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  stage_boundary_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  xmom_boundary_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  ymom_boundary_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  stage_explicit_update, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  xmom_explicit_update, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  ymom_explicit_update, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  max_speed_array, double evolve_max_timestep, double g, double epsilon, double h0, double limiting_threshold, hmpprt::s32 optimise_dry_cells)
;
#endif // __CUDACC__



# 538 "compute_fluxes.c"

#ifndef __CUDACC__
static hmpprt::CUDAGrid * compute_fluxes_central_structure_cuda_single_loop1D_1 = 0;
#else

extern "C" __global__ void compute_fluxes_central_structure_cuda_single_loop1D_1(double g, double epsilon, double h0, double limiting_threshold, hmpprt::s32 optimise_dry_cells, hmpprt::s32 N, double* timestep, hmpprt::s32* neighbours, hmpprt::s32* neighbour_edges, double* normals, double* edgelengths, double* radii, double* areas, hmpprt::s32* tri_full_flag, double* stage_edge_values, double* xmom_edge_values, double* ymom_edge_values, double* bed_edge_values, double* stage_boundary_values, double* xmom_boundary_values, double* ymom_boundary_values, double* stage_explicit_update, double* xmom_explicit_update, double* ymom_explicit_update, double* max_speed_array);
#endif // __CUDACC__




# 538 "compute_fluxes.c"

#ifdef __CUDACC__

extern "C" __global__ void compute_fluxes_central_structure_cuda_single_loop1D_1(double g, double epsilon, double h0, double limiting_threshold, hmpprt::s32 optimise_dry_cells, hmpprt::s32 N, double* timestep, hmpprt::s32* neighbours, hmpprt::s32* neighbour_edges, double* normals, double* edgelengths, double* radii, double* areas, hmpprt::s32* tri_full_flag, double* stage_edge_values, double* xmom_edge_values, double* ymom_edge_values, double* bed_edge_values, double* stage_boundary_values, double* xmom_boundary_values, double* ymom_boundary_values, double* stage_explicit_update, double* xmom_explicit_update, double* ymom_explicit_update, double* max_speed_array)
{
 # 592 "compute_fluxes.c"
 double max_speed_total;
 # 595 "compute_fluxes.c"
 double inv_area;
 # 571 "compute_fluxes.c"
 hmpprt::s32 k_1;
 # 601 "compute_fluxes.c"
 k_1 = (hmpprt::gr_atidf());
 # 601 "compute_fluxes.c"
 if (k_1 > N - 1)
 {
  # 601 "compute_fluxes.c"
  goto __hmppcg_label_1;
 }
 # 601 "compute_fluxes.c"
 max_speed_total = (double) 0.0;
 # 602 "compute_fluxes.c"
 *(stage_explicit_update + k_1) = (double) 0.0;
 # 603 "compute_fluxes.c"
 *(xmom_explicit_update + k_1) = (double) 0.0;
 # 604 "compute_fluxes.c"
 *(ymom_explicit_update + k_1) = (double) 0.0;
 # 606 "compute_fluxes.c"
 hmpprt::s32 i_1;
 # 606 "compute_fluxes.c"
 # 607 "compute_fluxes.c"
 for (i_1 = 0 ; i_1 <= 2 ; i_1 = i_1 + 1)
 {
  # 608 "compute_fluxes.c"
  hmpprt::s32 ki_1;
  # 608 "compute_fluxes.c"
  hmpprt::s32 ni_1;
  # 608 "compute_fluxes.c"
  double z_left_1;
  # 608 "compute_fluxes.c"
  double q_left_0_1;
  # 608 "compute_fluxes.c"
  double q_right_0_1;
  # 608 "compute_fluxes.c"
  double z_right_1;
  # 608 "compute_fluxes.c"
  double n1_1;
  # 608 "compute_fluxes.c"
  double q_left_1_1;
  # 608 "compute_fluxes.c"
  double n2_1;
  # 608 "compute_fluxes.c"
  double q_left_2_1;
  # 608 "compute_fluxes.c"
  double q_right_1_1;
  # 608 "compute_fluxes.c"
  double q_right_2_1;
  # 608 "compute_fluxes.c"
  double z_1;
  # 608 "compute_fluxes.c"
  double q_left_rotated_1_1;
  # 608 "compute_fluxes.c"
  double h_left_1;
  # 608 "compute_fluxes.c"
  double u_left_1;
  # 608 "compute_fluxes.c"
  double q_right_rotated_1_1;
  # 608 "compute_fluxes.c"
  double h_right_1;
  # 608 "compute_fluxes.c"
  double u_right_1;
  # 608 "compute_fluxes.c"
  double soundspeed_left_1;
  # 608 "compute_fluxes.c"
  double s_max_1;
  # 608 "compute_fluxes.c"
  double soundspeed_right_1;
  # 608 "compute_fluxes.c"
  double s_min_1;
  # 608 "compute_fluxes.c"
  double uh_left_1;
  # 608 "compute_fluxes.c"
  double q_left_rotated_2_1;
  # 608 "compute_fluxes.c"
  double uh_right_1;
  # 608 "compute_fluxes.c"
  double q_right_rotated_2_1;
  # 608 "compute_fluxes.c"
  double denom_1;
  # 608 "compute_fluxes.c"
  double flux_left_0_1;
  # 608 "compute_fluxes.c"
  double flux_right_0_1;
  # 608 "compute_fluxes.c"
  double edgeflux_0_1;
  # 608 "compute_fluxes.c"
  double edgeflux_0_2;
  # 608 "compute_fluxes.c"
  double flux_left_1_1;
  # 608 "compute_fluxes.c"
  double flux_right_1_1;
  # 608 "compute_fluxes.c"
  double edgeflux_1_1;
  # 608 "compute_fluxes.c"
  double edgeflux_1_2;
  # 608 "compute_fluxes.c"
  double flux_left_2_1;
  # 608 "compute_fluxes.c"
  double flux_right_2_1;
  # 608 "compute_fluxes.c"
  double edgeflux_2_1;
  # 608 "compute_fluxes.c"
  double edgeflux_2_2;
  # 608 "compute_fluxes.c"
  double max_speed_1;
  # 608 "compute_fluxes.c"
  double edgeflux_1_3;
  # 608 "compute_fluxes.c"
  double edgeflux_2_3;
  # 608 "compute_fluxes.c"
  double edgeflux_0_3;
  # 608 "compute_fluxes.c"
  double length_1;
  # 608 "compute_fluxes.c"
  double edgeflux_1_4;
  # 608 "compute_fluxes.c"
  double edgeflux_2_4;
  # 608 "compute_fluxes.c"
  double edgeflux_0_4;
  # 608 "compute_fluxes.c"
  double edgeflux_1_5;
  # 608 "compute_fluxes.c"
  double edgeflux_2_5;
  # 608 "compute_fluxes.c"
  ki_1 = k_1 * 3 + i_1;
  # 610 "compute_fluxes.c"
  q_left_0_1 = *(stage_edge_values + ki_1);
  # 611 "compute_fluxes.c"
  q_left_1_1 = *(xmom_edge_values + ki_1);
  # 612 "compute_fluxes.c"
  q_left_2_1 = *(ymom_edge_values + ki_1);
  # 613 "compute_fluxes.c"
  z_left_1 = *(bed_edge_values + ki_1);
  # 614 "compute_fluxes.c"
  ni_1 = *(neighbours + ki_1);
  # 615 "compute_fluxes.c"
  if (ni_1 < 0)
  {
   # 616 "compute_fluxes.c"
   hmpprt::s32 m_1;
   # 616 "compute_fluxes.c"
   m_1 =  - (ni_1 + 1);
   # 618 "compute_fluxes.c"
   q_right_0_1 = *(stage_boundary_values + m_1);
   # 619 "compute_fluxes.c"
   q_right_1_1 = *(xmom_boundary_values + m_1);
   # 620 "compute_fluxes.c"
   q_right_2_1 = *(ymom_boundary_values + m_1);
   # 621 "compute_fluxes.c"
   z_right_1 = z_left_1;
  }
  else
  {
   # 623 "compute_fluxes.c"
   hmpprt::s32 m_2;
   # 623 "compute_fluxes.c"
   hmpprt::s32 nm_1;
   # 623 "compute_fluxes.c"
   m_2 = *(neighbour_edges + ki_1);
   # 624 "compute_fluxes.c"
   nm_1 = ni_1 * 3 + m_2;
   # 626 "compute_fluxes.c"
   q_right_0_1 = *(stage_edge_values + nm_1);
   # 627 "compute_fluxes.c"
   q_right_1_1 = *(xmom_edge_values + nm_1);
   # 628 "compute_fluxes.c"
   q_right_2_1 = *(ymom_edge_values + nm_1);
   # 629 "compute_fluxes.c"
   z_right_1 = *(bed_edge_values + nm_1);
  }
  # 632 "compute_fluxes.c"
  if (optimise_dry_cells)
  {
   # 633 "compute_fluxes.c"
   if (fabs(q_left_0_1 - z_left_1) < epsilon && fabs(q_right_0_1 - z_right_1) < epsilon)
   {
    # 636 "compute_fluxes.c"
    continue;
   }
  }
  # 639 "compute_fluxes.c"
  n1_1 = *(normals + 2 * ki_1);
  # 640 "compute_fluxes.c"
  n2_1 = *(normals + (2 * ki_1 + 1));
  # 654 "compute_fluxes.c"
  q_left_rotated_1_1 = n1_1 * q_left_1_1 + n2_1 * q_left_2_1;
  # 655 "compute_fluxes.c"
  q_left_rotated_2_1 =  - n2_1 * q_left_1_1 + n1_1 * q_left_2_1;
  # 658 "compute_fluxes.c"
  q_right_rotated_1_1 = n1_1 * q_right_1_1 + n2_1 * q_right_2_1;
  # 659 "compute_fluxes.c"
  q_right_rotated_2_1 =  - n2_1 * q_right_1_1 + n1_1 * q_right_2_1;
  # 666 "compute_fluxes.c"
  z_1 = (double) 0.5 * (z_left_1 + z_right_1);
  # 670 "compute_fluxes.c"
  h_left_1 = q_left_0_1 - z_1;
  # 671 "compute_fluxes.c"
  uh_left_1 = q_left_rotated_1_1;
  # 675 "compute_fluxes.c"
  if (h_left_1 < limiting_threshold)
  {
   # 677 "compute_fluxes.c"
   if (h_left_1 < epsilon)
   {
    # 678 "compute_fluxes.c"
    h_left_1 = (double) 0.0;
    # 679 "compute_fluxes.c"
    u_left_1 = (double) 0.0;
   }
   else
   {
    # 681 "compute_fluxes.c"
    u_left_1 = q_left_rotated_1_1 / (h_left_1 + h0 / h_left_1);
   }
   # 684 "compute_fluxes.c"
   uh_left_1 = u_left_1 * h_left_1;
  }
  else
  {
   # 686 "compute_fluxes.c"
   u_left_1 = q_left_rotated_1_1 / h_left_1;
  }
  # 690 "compute_fluxes.c"
  h_right_1 = q_right_0_1 - z_1;
  # 691 "compute_fluxes.c"
  uh_right_1 = q_right_rotated_1_1;
  # 695 "compute_fluxes.c"
  if (h_right_1 < limiting_threshold)
  {
   # 697 "compute_fluxes.c"
   if (h_right_1 < epsilon)
   {
    # 698 "compute_fluxes.c"
    h_right_1 = (double) 0.0;
    # 699 "compute_fluxes.c"
    u_right_1 = (double) 0.0;
   }
   else
   {
    # 701 "compute_fluxes.c"
    u_right_1 = q_right_rotated_1_1 / (h_right_1 + h0 / h_right_1);
   }
   # 704 "compute_fluxes.c"
   uh_right_1 = u_right_1 * h_right_1;
  }
  else
  {
   # 706 "compute_fluxes.c"
   u_right_1 = q_right_rotated_1_1 / h_right_1;
  }
  # 716 "compute_fluxes.c"
  soundspeed_left_1 = sqrt(g * h_left_1);
  # 717 "compute_fluxes.c"
  soundspeed_right_1 = sqrt(g * h_right_1);
  # 719 "compute_fluxes.c"
  s_max_1 = u_left_1 + soundspeed_left_1;
  # 720 "compute_fluxes.c"
  if (s_max_1 < u_right_1 + soundspeed_right_1)
  {
   # 721 "compute_fluxes.c"
   s_max_1 = u_right_1 + soundspeed_right_1;
  }
  # 724 "compute_fluxes.c"
  if (s_max_1 < (double) 0.0)
  {
   # 725 "compute_fluxes.c"
   s_max_1 = (double) 0.0;
  }
  # 728 "compute_fluxes.c"
  s_min_1 = u_left_1 - soundspeed_left_1;
  # 729 "compute_fluxes.c"
  if (s_min_1 > u_right_1 - soundspeed_right_1)
  {
   # 730 "compute_fluxes.c"
   s_min_1 = u_right_1 - soundspeed_right_1;
  }
  # 732 "compute_fluxes.c"
  if (s_min_1 > (double) 0.0)
  {
   # 733 "compute_fluxes.c"
   s_min_1 = (double) 0.0;
  }
  # 737 "compute_fluxes.c"
  flux_left_0_1 = u_left_1 * h_left_1;
  # 738 "compute_fluxes.c"
  flux_left_1_1 = u_left_1 * uh_left_1 + (double) 0.5 * g * h_left_1 * h_left_1;
  # 739 "compute_fluxes.c"
  flux_left_2_1 = u_left_1 * q_left_rotated_2_1;
  # 741 "compute_fluxes.c"
  flux_right_0_1 = u_right_1 * h_right_1;
  # 742 "compute_fluxes.c"
  flux_right_1_1 = u_right_1 * uh_right_1 + (double) 0.5 * g * h_right_1 * h_right_1;
  # 743 "compute_fluxes.c"
  flux_right_2_1 = u_right_1 * q_right_rotated_2_1;
  # 746 "compute_fluxes.c"
  denom_1 = s_max_1 - s_min_1;
  # 747 "compute_fluxes.c"
  if (denom_1 < epsilon)
  {
   # 748 "compute_fluxes.c"
   edgeflux_0_3 = (double) 0.0;
   # 749 "compute_fluxes.c"
   edgeflux_1_4 = (double) 0.0;
   # 750 "compute_fluxes.c"
   edgeflux_2_4 = (double) 0.0;
   # 752 "compute_fluxes.c"
   max_speed_1 = (double) 0.0;
  }
  else
  {
   # 755 "compute_fluxes.c"
   double inverse_denominator_1;
   # 755 "compute_fluxes.c"
   inverse_denominator_1 = (double) 1.0 / denom_1;
   # 762 "compute_fluxes.c"
   edgeflux_0_1 = s_max_1 * flux_left_0_1 - s_min_1 * flux_right_0_1;
   # 764 "compute_fluxes.c"
   edgeflux_0_2 = edgeflux_0_1 + s_max_1 * s_min_1 * (q_right_0_1 - q_left_0_1);
   # 765 "compute_fluxes.c"
   edgeflux_0_3 = edgeflux_0_2 * inverse_denominator_1;
   # 767 "compute_fluxes.c"
   edgeflux_1_1 = s_max_1 * flux_left_1_1 - s_min_1 * flux_right_1_1;
   # 769 "compute_fluxes.c"
   edgeflux_1_2 = edgeflux_1_1 + s_max_1 * s_min_1 * (q_right_rotated_1_1 - q_left_rotated_1_1);
   # 770 "compute_fluxes.c"
   edgeflux_1_3 = edgeflux_1_2 * inverse_denominator_1;
   # 772 "compute_fluxes.c"
   edgeflux_2_1 = s_max_1 * flux_left_2_1 - s_min_1 * flux_right_2_1;
   # 774 "compute_fluxes.c"
   edgeflux_2_2 = edgeflux_2_1 + s_max_1 * s_min_1 * (q_right_rotated_2_1 - q_left_rotated_2_1);
   # 775 "compute_fluxes.c"
   edgeflux_2_3 = edgeflux_2_2 * inverse_denominator_1;
   # 781 "compute_fluxes.c"
   max_speed_1 = fabs(s_max_1);
   # 782 "compute_fluxes.c"
   if (max_speed_1 < fabs(s_min_1))
   {
    # 783 "compute_fluxes.c"
    max_speed_1 = fabs(s_min_1);
   }
   # 788 "compute_fluxes.c"
   edgeflux_1_4 = n1_1 * edgeflux_1_3 - n2_1 * edgeflux_2_3;
   # 789 "compute_fluxes.c"
   edgeflux_2_4 = n2_1 * edgeflux_1_3 + n1_1 * edgeflux_2_3;
  }
  # 794 "compute_fluxes.c"
  length_1 = *(edgelengths + ki_1);
  # 795 "compute_fluxes.c"
  edgeflux_0_4 = edgeflux_0_3 * length_1;
  # 796 "compute_fluxes.c"
  edgeflux_1_5 = edgeflux_1_4 * length_1;
  # 797 "compute_fluxes.c"
  edgeflux_2_5 = edgeflux_2_4 * length_1;
  # 799 "compute_fluxes.c"
  *(stage_explicit_update + k_1) = *(stage_explicit_update + k_1) - edgeflux_0_4;
  # 800 "compute_fluxes.c"
  *(xmom_explicit_update + k_1) = *(xmom_explicit_update + k_1) - edgeflux_1_5;
  # 801 "compute_fluxes.c"
  *(ymom_explicit_update + k_1) = *(ymom_explicit_update + k_1) - edgeflux_2_5;
  # 803 "compute_fluxes.c"
  if (*(tri_full_flag + k_1) == 1)
  {
   # 804 "compute_fluxes.c"
   if (max_speed_1 > epsilon)
   {
    # 806 "compute_fluxes.c"
    if (*(timestep + k_1) > *(radii + k_1) / max_speed_1)
    {
     # 807 "compute_fluxes.c"
     *(timestep + k_1) = *(radii + k_1) / max_speed_1;
    }
    # 809 "compute_fluxes.c"
    if (ni_1 >= 0)
    {
     # 811 "compute_fluxes.c"
     if (*(timestep + k_1) > *(radii + ni_1) / max_speed_1)
     {
      # 812 "compute_fluxes.c"
      *(timestep + k_1) = *(radii + ni_1) / max_speed_1;
     }
    }
   }
  }
  # 817 "compute_fluxes.c"
  if (ni_1 < 0 || ni_1 > k_1)
  {
   # 818 "compute_fluxes.c"
   max_speed_total = fmax(max_speed_total, max_speed_1);
  }
 }
 # 822 "compute_fluxes.c"
 # 822 "compute_fluxes.c"
 inv_area = (double) 1.0 / *(areas + k_1);
 # 823 "compute_fluxes.c"
 *(stage_explicit_update + k_1) = *(stage_explicit_update + k_1) * inv_area;
 # 824 "compute_fluxes.c"
 *(xmom_explicit_update + k_1) = *(xmom_explicit_update + k_1) * inv_area;
 # 825 "compute_fluxes.c"
 *(ymom_explicit_update + k_1) = *(ymom_explicit_update + k_1) * inv_area;
 # 827 "compute_fluxes.c"
 *(max_speed_array + k_1) = max_speed_total;
 # 538 "compute_fluxes.c"
 __hmppcg_label_1:;
}
#endif // __CUDACC__



# 538 "compute_fluxes.c"

#ifndef __CUDACC__
void compute_fluxes_central_structure_cuda_single_internal_1(hmpprt::s32 N, hmpprt::s32 N3, hmpprt::s32 N6, hmpprt::s32 N2, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  timestep, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s32>  neighbours, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s32>  neighbour_edges, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  normals, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  edgelengths, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  radii, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  areas, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s32>  tri_full_flag, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  stage_edge_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  xmom_edge_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  ymom_edge_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  bed_edge_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  stage_boundary_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  xmom_boundary_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  ymom_boundary_values, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  stage_explicit_update, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  xmom_explicit_update, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  ymom_explicit_update, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double>  max_speed_array, double evolve_max_timestep, double g, double epsilon, double h0, double limiting_threshold, hmpprt::s32 optimise_dry_cells)
{
 # 538 "compute_fluxes.c"
 if (N - 1 >= 0)
 {
  hmpprt::CUDAGridCall __hmppcg_call;
  __hmppcg_call.setSizeX((N - 1) / 128 + 1);
  __hmppcg_call.setSizeY(1);
  __hmppcg_call.setBlockSizeX(32);
  __hmppcg_call.setBlockSizeY(4);
  __hmppcg_call.addLocalParameter(&g, 8, "g");
  __hmppcg_call.addLocalParameter(&epsilon, 8, "epsilon");
  __hmppcg_call.addLocalParameter(&h0, 8, "h0");
  __hmppcg_call.addLocalParameter(&limiting_threshold, 8, "limiting_threshold");
  __hmppcg_call.addLocalParameter((hmpprt::s32) (optimise_dry_cells), "optimise_dry_cells");
  __hmppcg_call.addLocalParameter((hmpprt::s32) (N), "N");
  __hmppcg_call.addLocalParameter(&timestep, 8, "timestep");
  __hmppcg_call.addLocalParameter(&neighbours, 8, "neighbours");
  __hmppcg_call.addLocalParameter(&neighbour_edges, 8, "neighbour_edges");
  __hmppcg_call.addLocalParameter(&normals, 8, "normals");
  __hmppcg_call.addLocalParameter(&edgelengths, 8, "edgelengths");
  __hmppcg_call.addLocalParameter(&radii, 8, "radii");
  __hmppcg_call.addLocalParameter(&areas, 8, "areas");
  __hmppcg_call.addLocalParameter(&tri_full_flag, 8, "tri_full_flag");
  __hmppcg_call.addLocalParameter(&stage_edge_values, 8, "stage_edge_values");
  __hmppcg_call.addLocalParameter(&xmom_edge_values, 8, "xmom_edge_values");
  __hmppcg_call.addLocalParameter(&ymom_edge_values, 8, "ymom_edge_values");
  __hmppcg_call.addLocalParameter(&bed_edge_values, 8, "bed_edge_values");
  __hmppcg_call.addLocalParameter(&stage_boundary_values, 8, "stage_boundary_values");
  __hmppcg_call.addLocalParameter(&xmom_boundary_values, 8, "xmom_boundary_values");
  __hmppcg_call.addLocalParameter(&ymom_boundary_values, 8, "ymom_boundary_values");
  __hmppcg_call.addLocalParameter(&stage_explicit_update, 8, "stage_explicit_update");
  __hmppcg_call.addLocalParameter(&xmom_explicit_update, 8, "xmom_explicit_update");
  __hmppcg_call.addLocalParameter(&ymom_explicit_update, 8, "ymom_explicit_update");
  __hmppcg_call.addLocalParameter(&max_speed_array, 8, "max_speed_array");
  __hmppcg_call.launch(compute_fluxes_central_structure_cuda_single_loop1D_1, hmpprt::Context::getInstance()->getCUDADevice());
 }
 ;
}
#endif // __CUDACC__



# 538 "compute_fluxes.c"

#ifndef __CUDACC__
extern "C" CDLT_API  void compute_fluxes_central_structure_cuda_single(hmpprt::s32 N, hmpprt::s32 N3, hmpprt::s32 N6, hmpprt::s32 N2, double* timestep, hmpprt::s32* neighbours, hmpprt::s32* neighbour_edges, double* normals, double* edgelengths, double* radii, double* areas, hmpprt::s32* tri_full_flag, double* stage_edge_values, double* xmom_edge_values, double* ymom_edge_values, double* bed_edge_values, double* stage_boundary_values, double* xmom_boundary_values, double* ymom_boundary_values, double* stage_explicit_update, double* xmom_explicit_update, double* ymom_explicit_update, double* max_speed_array, double evolve_max_timestep, double g, double epsilon, double h0, double limiting_threshold, hmpprt::s32 optimise_dry_cells)
{
 # 1 "<preprocessor>"
 (compute_fluxes_central_structure_cuda_single_internal_1(N, N3, N6, N2, hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (timestep), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s32> (neighbours), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s32> (neighbour_edges), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (normals), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (edgelengths), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (radii), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (areas), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,hmpprt::s32> (tri_full_flag), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (stage_edge_values), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (xmom_edge_values), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (ymom_edge_values), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (bed_edge_values), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (stage_boundary_values), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (xmom_boundary_values), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (ymom_boundary_values), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (stage_explicit_update), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (xmom_explicit_update), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (ymom_explicit_update), hmpprt::DevicePtr<hmpprt::MS_CUDA_GLOB,double> (max_speed_array), evolve_max_timestep, g, epsilon, h0, limiting_threshold, optimise_dry_cells));
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
      compute_fluxes_central_structure_cuda_single_loop1D_1 = new hmpprt::CUDAGrid(hmpprt_module, "compute_fluxes_central_structure_cuda_single_loop1D_1");

    }
    hmpprt::Context::getInstance()->getGrouplet()->setTarget(hmpprt::CUDA);
    hmpprt::Context::getInstance()->getGrouplet()->addSignature("compute_fluxes_central_structure_cuda_single", "prototype compute_fluxes_central_structure_cuda_single(N: s32, N3: s32, N6: s32, N2: s32, timestep: ^cudaglob double, neighbours: ^cudaglob s32, neighbour_edges: ^cudaglob s32, normals: ^cudaglob double, edgelengths: ^cudaglob double, radii: ^cudaglob double, areas: ^cudaglob double, tri_full_flag: ^cudaglob s32, stage_edge_values: ^cudaglob double, xmom_edge_values: ^cudaglob double, ymom_edge_values: ^cudaglob double, bed_edge_values: ^cudaglob double, stage_boundary_values: ^cudaglob double, xmom_boundary_values: ^cudaglob double, ymom_boundary_values: ^cudaglob double, stage_explicit_update: ^cudaglob double, xmom_explicit_update: ^cudaglob double, ymom_explicit_update: ^cudaglob double, max_speed_array: ^cudaglob double, evolve_max_timestep: double, g: double, epsilon: double, h0: double, limiting_threshold: double, optimise_dry_cells: s32)");

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
      delete compute_fluxes_central_structure_cuda_single_loop1D_1;

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

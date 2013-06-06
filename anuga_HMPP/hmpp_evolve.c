
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<assert.h>

#include "domain.h"
#include "gravity.c"
#include "compute_fluxes.c"

//#define ON_XE

#define TOLERANCE_RANGE 1e-14


extern int errs_x, errs_y;


int main(int argc, char *argv[])
{
    struct domain sw_domain;
    int cnt_s =0, cnt_x =0, cnt_y =0, cnt_m = 0;
    int i, j, N;
    
    char file_name[]="merimbula_cf.dat";
    
    get_domain_from_file( sw_domain, file_name);

    printf("%d %lf\n", N, g);

    N = sw_domain.number_of_elements;
    
    printf(" --> Enter Kernel with Single Function\n");
    #pragma hmpp cf_central_single callsite
    compute_fluxes_central_structure_cuda_single(
            N, 
            N*3,
            N*6,
            N2,

            timestep,
            neighbours,
            neighbour_edges,
            normals,
            edgelengths,
            radii,
            areas,
            tri,
            se,
            xe,
            ye,
            be,
            sb,
            xb,
            yb,
            su,
            xu,
            yu,
            max_speed_array,

            evolve_max_timestep, 
            g, 
            epsilon,
            h0,
            limiting_threshold,
            optimise_dry_cells);

    /*
    printf(" --> Enter C\n");
    compute_fluxes_central_structure_cuda_single( N, N*3, N*6, N2, timestep, n, ne, normals, el, radii, a, tri, se, xe, ye, be, sb, xb, yb, su, xu, yu, max_speed_array, evolve_max_timestep, g, epsilon, h0, limiting_threshold, optimise_dry_cells);



    printf(" --> Enter Kernel with Multi-Function\n");
    #pragma hmpp cf_central callsite
    compute_fluxes_central_structure_CUDA(
            N, 
            N*3,
            N*6,
            N2,

            timestep,
            neighbours,
            neighbour_edges,
            normals,
            edgelengths,
            radii,
            areas,
            tri,
            se,
            xe,
            ye,
            be,
            sb,
            xb,
            yb,
            su,
            xu,
            yu,
            max_speed_array,

            evolve_max_timestep, 
            g, 
            epsilon,
            h0,
            limiting_threshold,
            optimise_dry_cells);
    */



    printf(" --> Enter original C\n");
    compute_fluxes_central_structure_CUDA( 
            N, N*3, N*6, N2, 
            test_timestep, 
            neighbours, neighbour_edges, normals, 
            edgelengths, radii, areas, tri, 
            se, xe, ye, be, sb, xb, yb, 
            test_stage_explicit_update, 
            test_xmom_explicit_update, 
            test_ymom_explicit_update, 
            test_max_speed_array, 
            evolve_max_timestep, g, epsilon, h0, limiting_threshold, 
            optimise_dry_cells);



    for (cnt_s=cnt_x=cnt_y=cnt_m=i=0; i < N; i++)
    {
        if ( fabs(su[i] - test_stage_explicit_update[i]) >= TOLERANCE_RANGE)
        {
            cnt_s++;
            if (cnt_s <= 10)
                printf(" sta %d %lf %lf\n", i, su[i], test_stage_explicit_update[i]);
        }
        if ( fabs(xu[i] - test_xmom_explicit_update[i]) >= TOLERANCE_RANGE)
        {
            cnt_x++;
            if (cnt_x <= 10)
                printf(" xmom %d %lf %lf\n", i, xu[i], test_xmom_explicit_update[i]);
        }
        if ( fabs(yu[i] - test_ymom_explicit_update[i]) >= TOLERANCE_RANGE)
        {
            cnt_y++;
            if (cnt_y <= 10)
                printf(" ymom %d %lf %lf\n", i, yu[i], test_ymom_explicit_update[i]);
        }
        if ( fabs(max_speed_array[i] -test_max_speed_array[i]) >=TOLERANCE_RANGE)
        {    
            cnt_m++;
            if (cnt_m <= 10)
                printf(" max %d %lf %lf\n", i, max_speed_array[i], test_max_speed_array[i]);
        }
    }
    printf("se:%d  xe:%d  ye:%d  max:%d errors found\n", 
        cnt_s, cnt_x, cnt_y, cnt_m);


    #pragma hmpp gravity callsite
    gravity_wb(N, N*3, N*6, xu, yu, sv, se, sc, be, bc, vc, 
        normals, areas, edgelengths, g);

    gravity_wb_orig( xu, yu, sv, se, sc, be, bc, vc, normals, areas,
        edgelengths, test_xmom_explicit_update, test_ymom_explicit_update,
        N, g);

    printf("\nxe:%d ye:%d\n", errs_x, errs_y);
}


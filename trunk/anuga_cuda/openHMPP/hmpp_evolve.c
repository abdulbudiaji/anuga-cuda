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
    int N, N2, optimise_dry_cells;
    int cnt_s =0, cnt_x =0, cnt_y =0, cnt_m = 0;

    double evolve_max_timestep, 
            g,
            epsilon,
            h0,
            limiting_threshold;

    double * sv, // stage_vertex_values
            * sb, // stage_boundary_values
            * se, // stage_edge_values, 
            * sc, // stage_centroid_values,
            * su, // stage_explicit_update

            * be, // bed_edge_values, 
            * bc, // bed_centroid_values, 
            
            * xb, // xmom_boundary_values
            * xe, // xmom_edge_values
            * xu, // xmom_explicit_update

            * yb, // ymom_boundary_values
            * ye, // ymom_edge_values
            * yu, // ymom_explicit_update

            * normals, 
            * edgelengths, // edgelengths;
            * radii,
            * areas, // areas, 
            * vc; // vertex_coordinates, 

    int * tri; // tri_full_flag
    int * neighbours, // neighbours
         * neighbour_edges; // neighbour_edges

    double * timestep, * max_speed_array;

    // Testing group
    double * test_stage_explicit_update,
            * test_xmom_explicit_update,
            * test_ymom_explicit_update,
            * test_max_speed_array,
            * test_timestep;
    int i;

    char file_name[]="merimbula_cf.dat";

    freopen(file_name, "r", stdin);
    
    scanf("%d\n", &N);
    scanf("%d\n", &N2);
    scanf("%d\n", &optimise_dry_cells);
    scanf("%lf\n",&evolve_max_timestep);
    scanf("%lf\n",&g);
    scanf("%lf\n",&epsilon);
    scanf("%lf\n",&h0);
    scanf("%lf\n",&limiting_threshold);


    printf("%d %lf\n", N, g);

    
    sv = (double *)malloc( N*3*sizeof(double));
    assert(sv != NULL);
    sb = (double *)malloc( N2*sizeof(double));
    assert(sb != NULL);
    se = (double *)malloc( N*3*sizeof(double));
    assert(se != NULL);
    sc = (double *)malloc( N*sizeof(double));
    assert(sc != NULL);
    su = (double *)malloc( N*sizeof(double));
    assert(su != NULL);
    
    
    be = (double *)malloc( N*3*sizeof(double));
    assert(be != NULL);
    bc = (double *)malloc( N*sizeof(double));
    assert(bc != NULL);
    
    
    xb = (double *)malloc( N2*sizeof(double));
    assert(xb != NULL);
    xe = (double *)malloc( N*3*sizeof(double));
    assert(xe != NULL);
    xu = (double *)malloc( N*sizeof(double));
    assert(xu != NULL);

    yb = (double *)malloc( N2*sizeof(double));
    assert(yb != NULL);
    ye = (double *)malloc( N*3*sizeof(double));
    assert(ye != NULL);
    yu = (double *)malloc( N*sizeof(double));
    assert(yu != NULL);

    
    neighbours  = (int *)malloc( N*3*sizeof(int));
    assert(neighbours != NULL);
    neighbour_edges  = (int *)malloc( N*3*sizeof(int));
    assert(neighbour_edges != NULL);


    normals  = (double *)malloc( N*6*sizeof(double));
    assert(normals != NULL);
    edgelengths = (double *)malloc( N*3*sizeof(double));
    assert(edgelengths != NULL);
    radii  = (double *)malloc( N*sizeof(double));
    assert(radii != NULL);
    areas  = (double *)malloc( N*sizeof(double));
    assert(areas != NULL);
    
    tri  = (int *)malloc( N*sizeof(int));
    assert(tri != NULL);
    
    vc = (double *)malloc( N*6*sizeof(double));
    assert(vc != NULL);
    
    timestep = (double *)malloc(N*sizeof(double));
    assert(timestep != NULL);
    max_speed_array = (double *)malloc(N*sizeof(double));
    assert(max_speed_array != NULL);

    test_stage_explicit_update = (double *)malloc(N*sizeof(double));
    assert(test_stage_explicit_update != NULL);
    test_xmom_explicit_update = (double *)malloc(N*sizeof(double));
    assert(test_xmom_explicit_update != NULL);
    test_ymom_explicit_update = (double *)malloc(N*sizeof(double));
    assert(test_ymom_explicit_update != NULL);
    test_max_speed_array = (double *)malloc(N*sizeof(double));
    assert(test_max_speed_array != NULL);
    test_timestep = (double *)malloc(N*sizeof(double));
    assert(test_timestep != NULL);


    for(i=0; i < N; i++)
    {
        scanf("%lf %lf %lf\n", sv+i*3, sv+i*3 +1, sv+i*3 +2);
        scanf("%lf %lf %lf\n", se+i*3, se+i*3 +1, se+i*3 +2);
        scanf("%lf\n", sc+i);
        scanf("%lf\n", su+i);

        scanf("%lf %lf %lf\n",be+i*3, be+i*3 +1, be+i*3 +2);
        scanf("%lf\n", bc+i);

        scanf("%lf %lf %lf\n", xe+i*3, xe+i*3 +1, xe+i*3 +2);
        scanf("%lf\n", xu+i);
        scanf("%lf %lf %lf\n", ye+i*3, ye+i*3 +1, ye+i*3 +2);
        scanf("%lf\n", yu+i);

        scanf("%d %d %d\n", neighbours+i*3, neighbours+i*3 +1, neighbours+i*3 +2);

        scanf("%d %d %d\n", neighbour_edges+i*3, neighbour_edges+i*3 +1, neighbour_edges+i*3 +2);


        scanf("%lf %lf %lf %lf %lf %lf\n", normals+i*6, normals+i*6+1, normals+i*6+2,normals+i*6+3, normals+i*6+4, normals+i*6+5);

        scanf("%lf %lf %lf\n", edgelengths+i*3, edgelengths+i*3 +1, edgelengths+i*3 +2);

        scanf("%lf\n", radii+i);
        scanf("%lf\n", areas+i);
        scanf("%d\n", tri+i);

        scanf("%lf %lf\n", vc+i*6, vc+i*6 +1);
        scanf("%lf %lf\n", vc+i*6+2, vc+i*6 +3);
        scanf("%lf %lf\n", vc+i*6+4, vc+i*6 +5);
    }

    for(i=0; i < N2; i++)
    {
        scanf("%lf\n", sb+i);
        scanf("%lf\n", xb+i);
        scanf("%lf\n", yb+i);
    }

    cnt_m = 0;
    for(i=0; i < N; i++)
    {
        if (max_speed_array[i] != test_max_speed_array[i])
            cnt_m ++;
    }

    printf("\n --> Test initial input %d\n", cnt_m);
    
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


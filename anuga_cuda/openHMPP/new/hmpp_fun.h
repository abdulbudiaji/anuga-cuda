#define USING_DOUBLE
#define ON_XE

#ifdef USING_CPP
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>
using namespace std;
#else
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#endif

#ifdef USING_DOUBLE
#define TOLERANCE 0.000000000000001
#define DATA_TYPE double 

#else
#define TOLERANCE 0.0000001
#define DATA_TYPE float
#endif


// Shallow_water domain structure
#include "sw_domain.h"


int check_tolerance(DATA_TYPE a,DATA_TYPE b);


int evolve(struct domain D, double yieldstep, 
            double finaltime, double duration, int skip_initial_step);


#pragma hmpp gravity codelet, target=CUDA args[*].transfer=atcall
void gravity_wb( 
        int n, int n3, int n6, 
        DATA_TYPE xmom_explicit_update[n], 
        DATA_TYPE ymom_explicit_update[n], 

        DATA_TYPE stage_vertex_values[n3],
        DATA_TYPE stage_edge_values[n3],
        DATA_TYPE stage_centroid_values[n],

        DATA_TYPE bed_edge_values[n3],
        DATA_TYPE bed_centroid_values[n],

        DATA_TYPE vertex_coordinates[n6],

        DATA_TYPE normals[n6],
        DATA_TYPE areas[n],
        DATA_TYPE edgelengths[n3],

        DATA_TYPE g );



#pragma hmpp cf_central codelet, target=CUDA args[*].transfer=atcall
void compute_fluxes_central_structure_CUDA(
        int N,
        int N3,
        int N6,
        int N2,

        double timestep[N],
        int neighbours[N3],
        int neighbour_edges[N3],
        double normals[N6],
        double edgelengths[N3],
        double radii[N],
        double areas[N],
        int tri_full_flag[N],
        double stage_edge_values[N3],
        double xmom_edge_values[N3],
        double ymom_edge_values[N3],
        double bed_edge_values[N3],
        double stage_boundary_values[N2],
        double xmom_boundary_values[N2],
        double ymom_boundary_values[N2],
        double stage_explicit_update[N],
        double xmom_explicit_update[N],
        double ymom_explicit_update[N],
        double max_speed_array[N],

        double evolve_max_timestep,
        double g,
        double epsilon,
        double h0,
        double limiting_threshold,
        int optimise_dry_cells);



#pragma hmpp cf_central_single codelet, target=CUDA args[*].transfer=atcall
void compute_fluxes_central_structure_cuda_single(
        int N,
        int N3,
        int N6,
        int N2,

        double timestep[N],
        int neighbours[N3],
        int neighbour_edges[N3],
        double normals[N6],
        double edgelengths[N3],
        double radii[N],
        double areas[N],
        int tri_full_flag[N],
        double stage_edge_values[N3],
        double xmom_edge_values[N3],
        double ymom_edge_values[N3],
        double bed_edge_values[N3],
        double stage_boundary_values[N2],
        double xmom_boundary_values[N2],
        double ymom_boundary_values[N2],
        double stage_explicit_update[N],
        double xmom_explicit_update[N],
        double ymom_explicit_update[N],
        double max_speed_array[N],

        double evolve_max_timestep,
        double g,
        double epsilon,
        double h0,
        double limiting_threshold,
        int optimise_dry_cells);



void gravity_wb_orig(
        DATA_TYPE * xmom_explicit_update, 
        DATA_TYPE * ymom_explicit_update, 
        DATA_TYPE * stage_vertex_values, 
        DATA_TYPE * stage_edge_values, 
        DATA_TYPE * stage_centroid_values, 
        DATA_TYPE * bed_edge_values, 
        DATA_TYPE * bed_centroid_values, 
        DATA_TYPE * vertex_coordinates, 
        DATA_TYPE * normals, 
        DATA_TYPE * areas, 
        DATA_TYPE * edgelengths,
        DATA_TYPE * test_xe,
        DATA_TYPE * test_ye,
        int N,
        DATA_TYPE g);



void gravity_call(
        int n, int n3, int n6, 
        DATA_TYPE xmom_explicit_update[n], 
        DATA_TYPE ymom_explicit_update[n], 

        DATA_TYPE stage_vertex_values[n3],
        DATA_TYPE stage_edge_values[n3],
        DATA_TYPE stage_centroid_values[n],

        DATA_TYPE bed_edge_values[n3],
        DATA_TYPE bed_centroid_values[n],

        DATA_TYPE vertex_coordinates[n6],

        DATA_TYPE normals[n6],
        DATA_TYPE areas[n],
        DATA_TYPE edgelengths[n3],

        DATA_TYPE g );



void test_call();



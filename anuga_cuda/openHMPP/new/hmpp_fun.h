// Using double floating variable
#define USING_DOUBLE
// When porting to Xe
//#define ON_XE
// Putting directives along with the implementation
//#define USING_LOCAL_DIRECTIVES
#define USING_GLOBAL_DIRECTIVES
// Do not back and force between Python and C
#define EVOLVE_ALL_IN_C

#define USING_MIRROR_DATA

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
#include <string.h>
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

double evolve(struct domain * D, double yieldstep, 
            double finaltime, double duration,
            double epsilon, int skip_initial_step,
            int step);

void allocate_upload_data(struct domain *D);

void download_data(struct domain *D);

int _distribute_to_vertices_and_edges(struct domain * D);

int _extrapolate_second_order_sw(struct domain * D);

void test_extrapolate_second_order_and_limit_by_vertex( struct domain *D);

void test_extrapolate_second_order_and_limit_by_vertex_normal( struct domain *D);

double compute_fluxes(struct domain * D);

int update_boundary(struct domain * D);

int update_ghosts(struct domain * D);

int update_extrema(struct domain * D);




#ifdef USING_GLOBAL_DIRECTIVES
#ifdef USING_MIRROR_DATA
#pragma hmpp gravity codelet, target=CUDA, transfer=atcall, & 
#pragma hmpp & args[xmom_explicit_update, ymom_explicit_update, &
#pragma hmpp & stage_vertex_values, stage_edge_values, stage_centroid_values, &
#pragma hmpp & bed_edge_values, bed_centroid_values, vertex_coordinates, &
#pragma hmpp & normals, areas, edgelengths].mirror, &
#pragma hmpp & args[xmom_explicit_update, ymom_explicit_update, &
#pragma hmpp & stage_vertex_values, stage_edge_values, stage_centroid_values, &
#pragma hmpp & bed_edge_values, bed_centroid_values, vertex_coordinates, &
#pragma hmpp & normals, areas, edgelengths].transfer=manual
#else
#pragma hmpp gravity codelet, target=CUDA args[*].transfer=atcall
#endif
#endif
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



#ifdef USING_GLOBAL_DIRECTIVES
#ifdef USING_MIRROR_DATA
#pragma hmpp cf_central codelet, target=CUDA, transfer=atcall, &
#pragma hmpp & args[timestep, neighbours, neighbour_edges, normals, &
#pragma hmpp & edgelengths, radii, areas, tri_full_flag, stage_edge_values, &
#pragma hmpp & xmom_edge_values, ymom_edge_values, bed_edge_values, &
#pragma hmpp & stage_boundary_values, xmom_boundary_values, &
#pragma hmpp & ymom_boundary_values, stage_explicit_update, &
#pragma hmpp & xmom_explicit_update, ymom_explicit_update, &
#pragma hmpp & max_speed_array].mirror, &
#pragma hmpp & args[timestep, neighbours, neighbour_edges, normals, &
#pragma hmpp & edgelengths, radii, areas, tri_full_flag, stage_edge_values, &
#pragma hmpp & xmom_edge_values, ymom_edge_values, bed_edge_values, &
#pragma hmpp & stage_boundary_values, xmom_boundary_values, &
#pragma hmpp & ymom_boundary_values, stage_explicit_update, &
#pragma hmpp & xmom_explicit_update, ymom_explicit_update, &
#pragma hmpp & max_speed_array].transfer=manual
#else
#pragma hmpp cf_central codelet, target=CUDA args[*].transfer=atcall
#endif
#endif
void compute_fluxes_central_structure_CUDA(
        int N,
        int N3,
        int N6,
        int N2,

        double timestep[N],
        long neighbours[N3],
        long neighbour_edges[N3],
        double normals[N6],
        double edgelengths[N3],
        double radii[N],
        double areas[N],
        long tri_full_flag[N],
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



#ifdef USING_GLOBAL_DIRECTIVES
#ifdef USING_MIRROR_DATA
#pragma hmpp cf_central_single codelet, target=CUDA, transfer=atcall &
#pragma hmpp & args[timestep, neighbours, neighbour_edges, normals, &
#pragma hmpp & edgelengths, radii, areas, tri_full_flag, stage_edge_values, &
#pragma hmpp & xmom_edge_values, ymom_edge_values, bed_edge_values, &
#pragma hmpp & stage_boundary_values, xmom_boundary_values, &
#pragma hmpp & ymom_boundary_values, stage_explicit_update, &
#pragma hmpp & xmom_explicit_update, ymom_explicit_update, &
#pragma hmpp & max_speed_array].mirror &
#pragma hmpp & args[timestep, neighbours, neighbour_edges, normals, &
#pragma hmpp & edgelengths, radii, areas, tri_full_flag, stage_edge_values, &
#pragma hmpp & xmom_edge_values, ymom_edge_values, bed_edge_values, &
#pragma hmpp & stage_boundary_values, xmom_boundary_values, &
#pragma hmpp & ymom_boundary_values, stage_explicit_update, &
#pragma hmpp & xmom_explicit_update, ymom_explicit_update, &
#pragma hmpp & max_speed_array].transfer=manual
#else
#pragma hmpp cf_central_single codelet, target=CUDA args[*].transfer=atcall
#endif
#endif
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



#ifdef USING_GLOBAL_DIRECTIVES
#ifdef USING_MIRROR_DATA
#pragma hmpp extraFstOrder codelet, target=CUDA, transfer=atcall, &
#pragma hmpp & args[centroid_values, edge_values, vertex_values].mirror, &
#pragma hmpp & args[centroid_values, edge_values, vertex_values].transfer=manual
#else
#pragma hmpp extraFstOrder codelet, target=CUDA args[*].transfer=atcall
#endif
#endif
void extrapolate_first_order(
        int N,
        int N3,
        double centroid_values[N],
        double edge_values[N3],
        double vertex_values[N3]
        );



// swb2_domain.c
// Not used so far
#ifdef USING_GLOBAL_DIRECTIVES
#pragma hmpp extraSndOrderEdge codelet, target=CUDA args[*].transfer=atcall
#endif
void extrapolate_second_order_edge_sw(
        int number_of_elements,
        int optimise_dry_cells, 
        int extrapolate_velocity_second_order, 

        double epsilon,
        double minimum_allowed_height,
        double beta_w,
        double beta_w_dry,
        double beta_uh,
        double beta_uh_dry,
        double beta_vh,
        double beta_vh_dry,

        long* surrogate_neighbours,
        long* number_of_boundaries,

        double* centroid_coordinates,

        double* stage_centroid_values,
        double* elevation_centroid_values,
        double* xmom_centroid_values,
        double* ymom_centroid_values,

        double* edge_coordinates,

        double* stage_edge_values,
        double* elevation_edge_values,
        double* xmom_edge_values,
        double* ymom_edge_values,

        double* stage_vertex_values,
        double* xmom_vertex_values,
        double* ymom_vertex_values,
        double* elevation_vertex_values,

        double* stage_centroid_store,
        double* xmom_centroid_store,
        double* ymom_centroid_store,
        double* min_elevation_edgevalue,
        double* max_elevation_edgevalue,
        int* count_wet_neighbours
        );



//#pragma hmpp extraSndOrderLmtV codelet, target=CUDA args[*].transfer=atcall
void extrapolate_second_order_and_limit_by_vertex(
        int N,
        int N2,
        int N3,
        int N6,
        double beta,
        double * domain_centroid_coordinates,
        double * domain_vertex_coordinates,
        long * domain_number_of_boundaries,
        long * domain_surrogate_neighbours,
        long * domain_neighbours,

        double * quantity_centroid_values,
        double * quantity_vertex_values,
        double * quantity_edge_values,
        double * quantity_x_gradient,
        double * quantity_y_gradient
        );



void extrapolate_second_order_and_limit_by_vertex_normal(
        int N,
        int N2,
        int N3,
        int N6,
        double beta,

        double domain_centroid_coordinates[N2],
        double domain_vertex_coordinates[N6],
        long domain_number_of_boundaries[N],
        long domain_surrogate_neighbours[N3],
        long domain_neighbours[N3],

        double quantity_centroid_values[N],
        double quantity_vertex_values[N3],
        double quantity_edge_values[N3],
        double quantity_x_gradient[N],
        double quantity_y_gradient[N]
        );


#ifdef TESTING_EXTRA_2_LV_GROUP
#pragma hmpp <extra2LV> group, target=CUDA

#pragma hmpp <extra2LV> map, &
#pragma hmpp & args[cptGradients::N; extraFromGradient::N;lmtVByNeigh::N]
#pragma hmpp <extra2LV> map, &
#pragma hmpp & args[cptGradients::N2; extraFromGradient::N2]
#pragma hmpp <extra2LV> map, &
#pragma hmpp & args[cptGradients::N3; extraFromGradient::N3; lmtVByNeigh::N3]
#pragma hmpp <extra2LV> map, &
#pragma hmpp & args[cptGradients::centroids; extraFromGradient::centroids]
#pragma hmpp <extra2LV> map, &
#pragma hmpp & args[cptGradients::centroid_values; extraFromGradient::centroid_values; lmtVByNeigh::centroid_values]
#pragma hmpp <extra2LV> map, &
#pragma hmpp & args[cptGradients::a; extraFromGradient::a; lmtVByNeigh::x_gradient]
#pragma hmpp <extra2LV> map, &
#pragma hmpp & args[cptGradients::b; extraFromGradient::b; lmtVByNeigh::y_gradient]
#pragma hmpp <extra2LV> map, &
#pragma hmpp & args[extraFromGradient::vertex_values; lmtVByNeigh::vertex_values]
#pragma hmpp <extra2LV> map, &
#pragma hmpp & args[extraFromGradient::edge_values; lmtVByNeigh::edge_values]
#endif


#ifndef NON_DIRECTIVES_EXTRA2_VERTEX_CPTGRA
#pragma hmpp cptGradients codelet, target=CUDA, transfer=atcall, &
#pragma hmpp & args[ centroids, centroid_values, number_of_boundaries, &
#pragma hmpp & surrogate_neighbours, a, b].mirror, &
#pragma hmpp & args[ centroids, centroid_values, number_of_boundaries, &
#pragma hmpp & surrogate_neighbours, a, b].transfer=manual
#endif
void _compute_gradients(
        int N,
        int N2,
        int N3,
        double centroids[N2],
        double centroid_values[N],
        long number_of_boundaries[N],
        long surrogate_neighbours[N3],
        double a[N],
        double b[N]);



#ifndef NON_DIRECTIVES_EXTRA2_VERTEX_EXTRA_FROM_GRA
#pragma hmpp extraFromGradient codelet, target=CUDA, transfer=atcall, &
#pragma hmpp & args[ centroids, centroid_values, vertex_coordinates, vertex_values, a, &
#pragma hmpp & b].mirror, &
#pragma hmpp & args[ centroids, centroid_values, vertex_coordinates, vertex_values, a, &
#pragma hmpp & b].transfer=manual
#endif
void _extrapolate_from_gradient(
        int N,
        int N2,
        int N3,
        int N6,
        double centroids[N2],
        double centroid_values[N],
        double vertex_coordinates[N6],
        double vertex_values[N3],
        double edge_values[N3],
        double a[N],
        double b[N]);



#ifndef NON_DIRECTIVES_EXTRA2_VERTEX
#pragma hmpp lmtVByNeigh codelet, target=CUDA, transfer=atcall &
#pragma hmpp & args[centroid_values, vertex_values, edge_values, &
#pragma hmpp & neighbours, x_gradient, y_gradient].mirror &
#pragma hmpp & args[centroid_values, vertex_values, edge_values, &
#pragma hmpp & neighbours, x_gradient, y_gradient].transfer=manual
#endif
void _limit_vertices_by_all_neighbours(
        int N, 
        int N3,
        double beta,
        double centroid_values[N],
        double vertex_values[N3],
        double edge_values[N3],
        long   neighbours[N3],
        double x_gradient[N],
        double y_gradient[N]);



//#pragma hmpp extraSndOrderLmtE codelet, target=CUDA args[*].transfer=atcall
void extrapolate_second_order_and_limit_by_edge(
        int N,
        int N2,
        int N3,
        int N6,
        double beta,
        double * domain_centroid_coordinates,
        double * domain_vertex_coordinates,
        long * domain_number_of_boundaries,
        long * domain_surrogate_neighbours,
        long * domain_neighbours,

        double * quantity_centroid_values,
        double * quantity_vertex_values,
        double * quantity_edge_values,
        double * quantity_x_gradient,
        double * quantity_y_gradient
        );



#ifdef USING_GLOBAL_DIRECTIVES
#ifdef USING_MIRROR_DATA
#pragma hmpp balance codelet, target=CUDA, transfer=atcall &
#pragma hmpp & args[wc, zc, wv, zv, xmomc, ymomc, xmomv, ymomv].mirror &
#pragma hmpp & args[wc, zc, wv, zv, xmomc, ymomc, xmomv, ymomv].transfer=manual
#else
#pragma hmpp balance codelet, target=CUDA args[*].transfer=atcall
#endif
#endif
void balance_deep_and_shallow(
        int N,
        int N3,
        double H0,
        double alpha_balance,
        int tight_slope_limiters,
        int use_centroid_velocities,

        double wc[N],   // stage_centroid_values
        double zc[N],   // elevation_centroid_values
        double wv[N3],  // stage_vertex_values
        double zv[N3],  // elevation_vertex_values
        //double* hvbar,// Retire this
        double xmomc[N],  // xmom_centroid_values
        double ymomc[N],  // ymom_centroid_values
        double xmomv[N3],  // xmom_vertex_values
        double ymomv[N3]   // ymom_vertex_values
        ); 



#ifdef USING_GLOBAL_DIRECTIVES
#ifdef USING_MIRROR_DATA
#pragma hmpp setBoundaryE codelet, target=CUDA, transfer=atcall &
#pragma hmpp & args[vol_id, edge_id, boundary_values, edge_values].mirror &
#pragma hmpp & args[vol_id, edge_id, boundary_values, edge_values].transfer=manual
#else
#pragma hmpp setBoundaryE codelet, target=CUDA args[*].transfer=atcall
#endif
#endif
void set_boundary_values_from_edges(
        int Nb,
        int N3,
        long vol_id[Nb],
        long edge_id[Nb],
        double boundary_values[Nb],
        double edge_values[N3]
        );


// swb2_domain.c
// Not used so far
#ifdef USING_GLOBAL_DIRECTIVES
#ifdef USING_MIRROR_DATA
#pragma hmpp protectSWB2 codelet, target=CUDA, transfer=atcall &
#pragma hmpp & args[wc, wv, zc, zv, xmomc, ymomc, areas].mirror &
#pragma hmpp & args[wc, wv, zc, zv, xmomc, ymomc, areas].transfer=manual
#else
#pragma hmpp protectSWB2 codelet, target=CUDA args[*].transfer=atcall
#endif
#endif
void protect_swb2(
        long N,
        long N3,

        double minimum_allowed_height,
        double maximum_allowed_speed,
        double epsilon,
        
        double* wc,
        double* wv,
        double* zc,
        double* zv,
        double* xmomc,
        double* ymomc,
        double* areas);



#ifdef USING_GLOBAL_DIRECTIVES
#ifdef USING_MIRROR_DATA
#pragma hmpp protectSW codelet, target=CUDA, transfer=atcall &
#pragma hmpp & args[wc, zc, xmomc, ymomc].mirror &
#pragma hmpp & args[wc, zc, xmomc, ymomc].transfer=manual
#else
#pragma hmpp protectSW codelet, target=CUDA args[*].transfer=atcall
#endif
#endif
void protect_sw(
        int N,
        int N3,
        double minimum_allowed_height,
        double maximum_allowed_speed,
        double epsilon,

        double wc[N],   // stage_centroid_values
        double zc[N],   // bed_centroid_values
        double xmomc[N],// xmom_centroid_values
        double ymomc[N] //ymom_centroid_values 
        );


#ifdef USING_GLOBAL_DIRECTIVES
#ifdef USING_MIRROR_DATA
#pragma hmpp interpolateVtoE codelet, target=CUDA, transfer=atcall &
#pragma hmpp & args[vertex_values, edge_values].mirror &
#pragma hmpp & args[vertex_values, edge_values].transfer=manual
#else
#pragma hmpp interpolateVtoE codelet, target=CUDA args[*].transfer=atcall
#endif
#endif
void interpolate_from_vertices_to_edges(
        int N,
        int N3,
        double vertex_values[N3],
        double edge_values[N3]
        ); 


        
#ifdef USING_GLOBAL_DIRECTIVES
#ifdef USING_MIRROR_DATA
#pragma hmpp updateCentroidVH codelet, target=CUDA,  transfer=atcall &
#pragma hmpp & args[w_C, uh_C, vh_C, h_C, z_C, u_C, v_C, &
#pragma hmpp &      w_B, uh_B, vh_B, h_B, z_B, u_B, v_B].mirror &
#pragma hmpp & args[w_C, uh_C, vh_C, h_C, z_C, u_C, v_C, &
#pragma hmpp &      w_B, uh_B, vh_B, h_B, z_B, u_B, v_B].transfer=manual
#else
#pragma hmpp updateCentroidVH codelet, target=CUDA args[*].transfer=atcall
#endif
#endif
void _update_centroids_of_velocities_and_height(
        int N_c,
        int N_b,
        double w_C[N_c], // stage_centroid_values
        double uh_C[N_c],// xmomentum_centroid_values
        double vh_C[N_c],// ymomentum_centroid_values
        double h_C[N_c], // height_centroid_values
        double z_C[N_c], // elevation_centroid_values
        double u_C[N_c], // xvelocity_centroid_values
        double v_C[N_c], // yvelocity_centroid_values

        double w_B[N_b], // stage_boundary_values
        double uh_B[N_b],// xmomentum_boundary_values
        double vh_B[N_b],// ymomentum_boundary_values
        double h_B[N_b], // height_boundary_values
        double z_B[N_b], // elevation_boundary_values
        double u_B[N_b], // xvelocity_boundary_values
        double v_B[N_b] // yvelocity_boundary_values
        );



#ifdef USING_GLOBAL_DIRECTIVES
#ifdef USING_MIRROR_DATA
#pragma hmpp manFrictionFlat codelet, target=CUDA, transfer=atcall &
#pragma hmpp & args[w, zv, uh, vh, eta, xmom, ymom].mirror &
#pragma hmpp & args[w, zv, uh, vh, eta, xmom, ymom].transfer=manual
#else
#pragma hmpp manFrictionFlat codelet, target=CUDA args[*].transfer=atcall
#endif
#endif
void manning_friction_flat(
        int N,
        int N3,
        double g, 
        double eps, // minimum_allowed_height 

        double w[N],  // stage_centroid_values
        double zv[N3], // elevation_vertex_values
        double uh[N], // xmom_centroid_values
        double vh[N], // ymom_centroid_values
        double eta[N],// friction_centroid_values
        double xmom[N],//xmom_semi_implicit_update 
        double ymom[N]//ymom_semi_implicit_update 
        );



#ifdef USING_GLOBAL_DIRECTIVES
#ifdef USING_MIRROR_DATA
#pragma hmpp manFrictionSloped codelet, target=CUDA, transfer=atcall &
#pragma hmpp & args[x, w, zv, uh, vh, eta, xmom_update, ymom_update].mirror &
#pragma hmpp & args[x, w, zv, uh, vh, eta, xmom_update, &
#pragma hmpp & ymom_update].transfer=manual
#else
#pragma hmpp manFrictionSloped codelet, target=CUDA args[*].transfer=atcall
#endif
#endif
void manning_friction_sloped(
        int N,
        int N3,
        int N6,
        double g, 
        double eps, // minimum_allowed_height

        double x[N6],  // vertex_coordinates
        double w[N],  // stage_centroid_values
        double zv[N3], // elevation_vertex_values
        double uh[N], // xmom_centroid_values
        double vh[N], // ymom_centroid_values
        double eta[N],// friction_centroid_values
        double xmom_update[N],    // xmom_semi_implicit_update
        double ymom_update[N]    // ymom_semi_implicit_update
        );



#ifdef USING_GLOBAL_DIRECTIVES
#ifdef USING_MIRROR_DATA
#pragma hmpp extraSndVelocity codelet, target=CUDA, transfer=atcall &
#pragma hmpp & args[stage_centroid_values, bed_centroid_values, &
#pragma hmpp & xmom_centroid_values, xmom_centroid_store, &
#pragma hmpp & ymom_centroid_values, ymom_centroid_store].mirror &
#pragma hmpp & args[stage_centroid_values, bed_centroid_values, &
#pragma hmpp & xmom_centroid_values, xmom_centroid_store, &
#pragma hmpp & ymom_centroid_values, ymom_centroid_store].transfer=manual
#else
#pragma hmpp extraSndVelocity codelet, target=CUDA args[*].transfer=atcall
#endif
#endif
void extrapolate_second_order_velocity_true(
            int N,
            double minimum_allowed_height,
            double stage_centroid_values[N],
            double bed_centroid_values[N],
            double xmom_centroid_values[N],
            double xmom_centroid_store[N],
            double ymom_centroid_values[N],
            double ymom_centroid_store[N]
            );


            
#ifdef USING_GLOBAL_DIRECTIVES
#ifdef USING_MIRROR_DATA
#pragma hmpp extraSndOrderSWT codelet, target=CUDA, transfer=atcall &
#pragma hmpp & args[surrogate_neighbours, number_of_boundaries, &
#pragma hmpp & centroid_coordinates, stage_centroid_values, &
#pragma hmpp & bed_centroid_values, xmom_centroid_values, &
#pragma hmpp & ymom_centroid_values, vertex_coordinates, &
#pragma hmpp & stage_vertex_values, bed_vertex_values, &
#pragma hmpp & xmom_vertex_values, ymom_vertex_values].mirror &
#pragma hmpp & args[surrogate_neighbours, number_of_boundaries, &
#pragma hmpp & centroid_coordinates, stage_centroid_values, &
#pragma hmpp & bed_centroid_values, xmom_centroid_values, &
#pragma hmpp & ymom_centroid_values, vertex_coordinates, &
#pragma hmpp & stage_vertex_values, bed_vertex_values, &
#pragma hmpp & xmom_vertex_values, ymom_vertex_values].transfer=manual
#else
#pragma hmpp extraSndOrderSWT codelet, target=CUDA args[*].transfer=atcall
#endif
#endif
void extrapolate_second_order_sw_true (
        int N,
        int N2,
        int N3,
        int N6,
        double epsilon,
        double minimum_allowed_height,
        double beta_w,
        double beta_w_dry,
        double beta_uh,
        double beta_uh_dry,
        double beta_vh,
        double beta_vh_dry,
        int optimise_dry_cells,

        long surrogate_neighbours[N3],
        long number_of_boundaries[N],
        double centroid_coordinates[N2],

        double stage_centroid_values[N],
        double bed_centroid_values[N],
        double xmom_centroid_values[N],
        double ymom_centroid_values[N],

        double vertex_coordinates[N6],
        
        double stage_vertex_values[N3],
        double bed_vertex_values[N3],
        double xmom_vertex_values[N3],
        double ymom_vertex_values[N3]
        );



#ifdef USING_GLOBAL_DIRECTIVES
#ifdef USING_MIRROR_DATA
#pragma hmpp extraSndOrderSW codelet, target=CUDA, transfer=atcall &
#pragma hmpp & args[surrogate_neighbours, number_of_boundaries, &
#pragma hmpp & centroid_coordinates, stage_centroid_values, &
#pragma hmpp & bed_centroid_values, xmom_centroid_values, &
#pragma hmpp & ymom_centroid_values, vertex_coordinates, &
#pragma hmpp & stage_vertex_values, bed_vertex_values, &
#pragma hmpp & xmom_vertex_values, ymom_vertex_values, &
#pragma hmpp & stage_centroid_store, xmom_centroid_store, &
#pragma hmpp & ymom_centroid_store].mirror &
#pragma hmpp & args[surrogate_neighbours, number_of_boundaries, &
#pragma hmpp & centroid_coordinates, stage_centroid_values, &
#pragma hmpp & bed_centroid_values, xmom_centroid_values, &
#pragma hmpp & ymom_centroid_values, vertex_coordinates, &
#pragma hmpp & stage_vertex_values, bed_vertex_values, &
#pragma hmpp & xmom_vertex_values, ymom_vertex_values, &
#pragma hmpp & stage_centroid_store, xmom_centroid_store, &
#pragma hmpp & ymom_centroid_store].transfer=manual
#else
#pragma hmpp extraSndOrderSW codelet, target=CUDA args[*].transfer=atcall
#endif
#endif
void extrapolate_second_order_sw( 
        int N,
        int N2,
        int N3,
        int N6,
        double epsilon,
        double minimum_allowed_height,
        double beta_w,
        double beta_w_dry,
        double beta_uh,
        double beta_uh_dry,
        double beta_vh,
        double beta_vh_dry,
        int optimise_dry_cells,
        int extrapolate_velocity_second_order,

        long surrogate_neighbours[N3],
        long number_of_boundaries[N],
        double centroid_coordinates[N2],

        double stage_centroid_values[N],
        double bed_centroid_values[N],
        double xmom_centroid_values[N],
        double ymom_centroid_values[N],

        double vertex_coordinates[N6],
        
        double stage_vertex_values[N3],
        double bed_vertex_values[N3],
        double xmom_vertex_values[N3],
        double ymom_vertex_values[N3],
        double stage_centroid_store[N],
        double xmom_centroid_store[N],
        double ymom_centroid_store[N]
        );



#ifdef USING_GLOBAL_DIRECTIVES
#ifdef USING_MIRROR_DATA
#pragma hmpp extraSndOrderSWF codelet, target=CUDA, transfer=atcall &
#pragma hmpp & args[surrogate_neighbours, number_of_boundaries, &
#pragma hmpp & centroid_coordinates, stage_centroid_values, &
#pragma hmpp & bed_centroid_values, xmom_centroid_values, &
#pragma hmpp & ymom_centroid_values, vertex_coordinates, &
#pragma hmpp & stage_vertex_values, bed_vertex_values, &
#pragma hmpp & xmom_vertex_values, ymom_vertex_values].mirror &
#pragma hmpp & args[surrogate_neighbours, number_of_boundaries, &
#pragma hmpp & centroid_coordinates, stage_centroid_values, &
#pragma hmpp & bed_centroid_values, xmom_centroid_values, &
#pragma hmpp & ymom_centroid_values, vertex_coordinates, &
#pragma hmpp & stage_vertex_values, bed_vertex_values, &
#pragma hmpp & xmom_vertex_values, ymom_vertex_values].transfer=manual
#else
#pragma hmpp extraSndOrderSWF codelet, target=CUDA args[*].transfer=atcall
#endif
#endif
void extrapolate_second_order_sw_false (
        int N,
        int N3,
        int N6,
        double epsilon,
        double minimum_allowed_height,
        double beta_w,
        double beta_w_dry,
        double beta_uh,
        double beta_uh_dry,
        double beta_vh,
        double beta_vh_dry,
        int optimise_dry_cells,

        long surrogate_neighbours[N3],
        long number_of_boundaries[N],
        double centroid_coordinates[N3],

        double stage_centroid_values[N],
        double bed_centroid_values[N],
        double xmom_centroid_values[N],
        double ymom_centroid_values[N],

        double vertex_coordinates[N6],
        
        double stage_vertex_values[N3],
        double bed_vertex_values[N3],
        double xmom_vertex_values[N3],
        double ymom_vertex_values[N3]
        );



#ifdef USING_GLOBAL_DIRECTIVES
#ifdef USING_MIRROR_DATA
#pragma hmpp update codelet, target=CUDA, transfer=atcall &
#pragma hmpp & args[centroid_values, explicit_update, &
#pragma hmpp & semi_implicit_update].mirror &
#pragma hmpp & args[centroid_values, explicit_update, &
#pragma hmpp & semi_implicit_update].transfer=manual
#else
#pragma hmpp update codelet, target=CUDA args[*].transfer=atcall
#endif
#endif
void update(
        int N,
        double timestep,
        double centroid_values[N],
        double explicit_update[N],
        double semi_implicit_update[N]
        );
    


#ifdef USING_GLOBAL_DIRECTIVES
#ifdef USING_MIRROR_DATA
#pragma hmpp saxpyCen codelet, target=CUDA, transfer=atcall &
#pragma hmpp & args[centroid_values, centroid_backup_values].mirror &
#pragma hmpp & args[centroid_values, centroid_backup_values].transfer=manual
#else
#pragma hmpp saxpyCen codelet, target=CUDA args[*].transfer=atcall
#endif
#endif
void saxpy_centroid_values(
        int N,
        double a,
        double b,
        double centroid_values[N],
        double centroid_backup_values[N]
        );




#ifdef USING_GLOBAL_DIRECTIVES
#ifdef USING_MIRROR_DATA
#pragma hmpp evaRef codelet, target=CUDA, transfer=atcall, &
#pragma hmpp & args[stage_edge_values, bed_edge_values, height_edge_values, &
#pragma hmpp & xmom_edge_values, xvel_edge_values, yvel_edge_values, &
#pragma hmpp & stage_boundary_values, bed_boundary_values, &
#pragma hmpp & height_boundary_values, xmom_boundary_values, &
#pragma hmpp & ymom_boundary_values, xvel_boundary_values, &
#pragma hmpp & yvel_boundary_values].mirror, &
#pragma hmpp & args[stage_edge_values, bed_edge_values, height_edge_values, &
#pragma hmpp & xmom_edge_values, xvel_edge_values, yvel_edge_values, &
#pragma hmpp & stage_boundary_values, bed_boundary_values, &
#pragma hmpp & height_boundary_values, xmom_boundary_values, &
#pragma hmpp & ymom_boundary_values, xvel_boundary_values, &
#pragma hmpp & yvel_boundary_values].transfer=manual
#else
#pragma hmpp evaRef codelet, target=CUDA args[*].transfer=atcall
#endif
#endif
void evaluate_segment_reflective(
    int N1,   // Nids
    int N2,     // Nb
    int N3,
    int N6,

    long ids[N1],
    long vol_ids[N2],   // domain.boundary_cells
    long edge_ids[N2],  // domain.boundary_edges
    double normals[N6], 
    
    double stage_edge_values[N3],
    double bed_edge_values[N3],
    double height_edge_values[N3],
    double xmom_edge_values[N3],
    double ymom_edge_values[N3],
    double xvel_edge_values[N3],
    double yvel_edge_values[N3],

    double stage_boundary_values[N2],
    double bed_boundary_values[N2],
    double height_boundary_values[N2],
    double xmom_boundary_values[N2],
    double ymom_boundary_values[N2],
    double xvel_boundary_values[N2],
    double yvel_boundary_values[N2]
    );



// swb2
int _find_qmin_and_qmax(double dq0, double dq1, double dq2, 
               double *qmin, double *qmax);


int _limit_gradient(double *dqv, double qmin, double qmax, double beta_w);



// extrapolate_second_order_sw
int limit_gradient(
        double *dqv0, 
        double *dqv1, 
        double *dqv2, 
        double qmin, 
        double qmax, 
        double beta_w); 


int limit_gradient_old(
        double *dqv, 
        double qmin, 
        double qmax, 
        double beta_w);

int find_qmin_and_qmax(
        double dq0, 
        double dq1, 
        double dq2,
        double *qmin, 
        double *qmax); 



void test_call();
// protect_sw.c
void test_protect_sw();



void test_single( struct domain *D);

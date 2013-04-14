
// OpenHMPP ANUGA evolve function
//
// Zhe Weng 2013
#include <hmpp_fun.h>

int evolve( struct domain D, 
            double yieldstep, 
            double finaltime,
            double duration,
            int skip_initial_step
            )
{
    assert( D.beta_w >= 0 && D.beta_w <= 2.0 );

    #pragma hmpp cf_central callsite
    compute_fluxes_central_structure_CUDA(
            D.number_of_elements, 
            D.number_of_elements*3,
            D.number_of_elements*6,
            D.number_of_boundary_elements,

            D.timestep,
            D.neighbours,
            D.neighbour_edges,
            D.normals,
            D.edgelengths,
            D.radii,
            D.areas,
            D.tri_full_flag,
            
            D.stage_edge_values,
            D.xmom_edge_values,
            D.ymom_edge_values,
            D.bed_edge_values,
            
            D.stage_boundary_values,
            D.xmom_boundary_values,
            D.ymom_boundary_values,
            
            D.stage_explicit_update,
            D.xmom_explicit_update,
            D.ymom_explicit_update,

            D.max_speed,

            D.evolve_max_timestep, 
            D.g, 
            D.epsilon,
            D.h0,
            D.limiting_threshold,
            D.optimise_dry_cells);



    #pragma hmpp gravity callsite
    gravity_wb(
        D.number_of_elements,
        D.number_of_elements * 3,
        D.number_of_elements * 6,
        D.xmom_explicit_update,
        D.ymom_explicit_update,

        D.stage_vertex_values,
        D.stage_edge_values,
        D.stage_centroid_values,

        D.bed_edge_values,
        D.bed_centroid_values,

        D.vertex_coordinates,

        D.normals,
        D.areas,
        D.edgelengths,

        D.g
        );
}

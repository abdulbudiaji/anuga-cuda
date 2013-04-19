// C struct for domain and quantities
//
// Stephen Roberts 2012
// John Weng 2013

#include <stdlib.h>
#include <assert.h>


// Shared code snippets
#include "util_ext.h"

#ifndef SW_DOMAIN 
#define SW_DOMAIN
#include "sw_domain.h"
#endif


#define SAFE_MALLOC(a, n, type) \
    (assert((a= (type *) malloc ((n) * sizeof (type))) != NULL))

void get_edge_data(struct edge *E, struct domain *D, int k, int i) {
    // fill edge data (conserved and bed) for ith edge of kth triangle

    int k3i, k3i1, k3i2;

    k3i = 3 * k + i;
    k3i1 = 3 * k + (i + 1) % 3;
    k3i2 = 3 * k + (i + 2) % 3;

    E->cell_id = k;
    E->edge_id = i;

    E->w = D->stage_edge_values[k3i];
    E->z = D->bed_edge_values[k3i];
    E->h = E->w - E->z;
    E->uh = D->xmom_edge_values[k3i];
    E->vh = D->ymom_edge_values[k3i];

    E->w1 = D->stage_vertex_values[k3i1];
    E->z1 = D->bed_vertex_values[k3i1];
    E->h1 = E->w1 - E->z1;
    E->uh1 = D->xmom_vertex_values[k3i1];
    E->vh1 = D->ymom_vertex_values[k3i1];


    E->w2 = D->stage_vertex_values[k3i2];
    E->z2 = D->bed_vertex_values[k3i2];
    E->h2 = E->w2 - E->z2;
    E->uh2 = D->xmom_vertex_values[k3i2];
    E->vh2 = D->ymom_vertex_values[k3i2];

}


struct domain* get_python_domain(struct domain *D, PyObject *domain) {
    PyArrayObject
            *neighbours,
            *neighbour_edges,
            *normals,
            *edgelengths,
            *radii,
            *areas,
            *tri_full_flag,
            *already_computed_flux,
            *vertex_coordinates,
            *edge_coordinates,
            *centroid_coordinates,
            *number_of_boundaries,
            *surrogate_neighbours,
            *max_speed,
            *boundary;
            //*min_bed_edge_values,
            //*max_bed_edge_values,
            //*count_wet_neighbours;

    PyObject *quantities, *stage;

    D->number_of_elements   = get_python_integer(domain, "number_of_elements");
    D->epsilon              = get_python_double(domain, "epsilon");
    D->H0                   = get_python_double(domain, "H0");
    D->h0 = D->H0 * D->H0;
    D->limiting_threshold = 10 * D->H0;
    D->g                    = get_python_double(domain, "g");
    D->optimise_dry_cells   = get_python_integer(domain, "optimise_dry_cells");
    D->evolve_max_timestep  = get_python_double(domain, "evolve_max_timestep");
    D->extrapolate_velocity_second_order  = get_python_integer(domain, "extrapolate_velocity_second_order");
    D->minimum_allowed_height = get_python_double(domain, "minimum_allowed_height");

    
    D->CFL                  = get_python_double(domain, "CFL");
    D->flux_timestep        = get_python_double(domain, "flux_timestep");

    D->beta_w      = get_python_double(domain, "beta_w");;
    D->beta_w_dry  = get_python_double(domain, "beta_w_dry");
    D->beta_uh     = get_python_double(domain, "beta_uh");
    D->beta_uh_dry = get_python_double(domain, "beta_uh_dry");
    D->beta_vh     = get_python_double(domain, "beta_vh");
    D->beta_vh_dry = get_python_double(domain, "beta_vh_dry");


    
    neighbours = get_consecutive_array(domain, "neighbours");
    D->neighbours = (long *) neighbours->data;

    surrogate_neighbours = get_consecutive_array(domain, "surrogate_neighbours");
    D->surrogate_neighbours = (long *) surrogate_neighbours->data;

    neighbour_edges = get_consecutive_array(domain, "neighbour_edges");
    D->neighbour_edges = (long *) neighbour_edges->data;

    normals = get_consecutive_array(domain, "normals");
    D->normals = (double *) normals->data;

    edgelengths = get_consecutive_array(domain, "edgelengths");
    D->edgelengths = (double *) edgelengths->data;

    radii = get_consecutive_array(domain, "radii");
    D->radii = (double *) radii->data;

    areas = get_consecutive_array(domain, "areas");
    D->areas = (double *) areas->data;



    tri_full_flag = get_consecutive_array(domain, "tri_full_flag");
    D->tri_full_flag = (long *) tri_full_flag->data;

    already_computed_flux = get_consecutive_array(domain, "already_computed_flux");
    D->already_computed_flux = (long *) already_computed_flux->data;

    max_speed = get_consecutive_array(domain, "max_speed");
    D->max_speed = (double *) max_speed->data;

    SAFE_MALLOC( D->timestep_array, D->number_of_elements, double);



    vertex_coordinates = get_consecutive_array(domain, "vertex_coordinates");
    D->vertex_coordinates = (double *) vertex_coordinates->data;

    edge_coordinates = get_consecutive_array(domain, "edge_coordinates");
    D->edge_coordinates = (double *) edge_coordinates->data;

    centroid_coordinates = get_consecutive_array(domain, "centroid_coordinates");
    D->centroid_coordinates = (double *) centroid_coordinates->data;


    
    number_of_boundaries = get_consecutive_array(domain, "number_of_boundaries");
    D->number_of_boundaries = (long *) number_of_boundaries->data;



    // Some others
    SAFE_MALLOC( D->min_bed_edge_values, D->number_of_elements, double);
    SAFE_MALLOC( D->max_bed_edge_values, D->number_of_elements, double);
    SAFE_MALLOC( D->count_wet_neighbours, D->number_of_elements, int);

    

    // Quantities
    quantities = get_python_object(domain, "quantities");
    stage = PyDict_GetItemString(quantities,   "stage");
    boundary = get_consecutive_array( stage, "boundary_values");
    // Number of boundary elements
    D->number_of_boundary_elements = boundary->dimensions[0];


    // Edge values    
    D->stage_edge_values    = get_python_array_data_from_dict(quantities, "stage",      "edge_values");
    D->xmom_edge_values     = get_python_array_data_from_dict(quantities, "xmomentum",  "edge_values");
    D->ymom_edge_values     = get_python_array_data_from_dict(quantities, "ymomentum",  "edge_values");
    D->bed_edge_values      = get_python_array_data_from_dict(quantities, "elevation",  "edge_values");
    D->height_edge_values   = get_python_array_data_from_dict(quantities, "height",     "edge_values");
    D->xvelocity_edge_values= get_python_array_data_from_dict(quantities, "xvelocity",  "edge_values");
    D->yvelocity_edge_values= get_python_array_data_from_dict(quantities, "yvelocity",  "edge_values");


    // Centroid values
    D->stage_centroid_values    = get_python_array_data_from_dict(quantities, "stage",      "centroid_values");
    D->xmom_centroid_values     = get_python_array_data_from_dict(quantities, "xmomentum",  "centroid_values");
    D->ymom_centroid_values     = get_python_array_data_from_dict(quantities, "ymomentum",  "centroid_values");
    D->bed_centroid_values      = get_python_array_data_from_dict(quantities, "elevation",  "centroid_values");
    D->height_centroid_values   = get_python_array_data_from_dict(quantities, "height",     "centroid_values");
    D->xvelocity_centroid_values= get_python_array_data_from_dict(quantities, "xvelocity",  "centroid_values");
    D->yvelocity_centroid_values= get_python_array_data_from_dict(quantities, "yvelocity",  "centroid_values");
    // Store values
    //D->stage_centroid_store    = get_python_array_data_from_dict(quantities, "stage",      "centroid_store");
    //D->xmom_centroid_store     = get_python_array_data_from_dict(quantities, "xmomentum",  "centroid_store");
    //D->ymom_centroid_store     = get_python_array_data_from_dict(quantities, "ymomentum",  "centroid_store");
    SAFE_MALLOC( D->stage_centroid_store, D->number_of_elements, double);
    SAFE_MALLOC( D->xmom_centroid_store, D->number_of_elements, double);
    SAFE_MALLOC( D->ymom_centroid_store, D->number_of_elements, double);
    
    // Backup values
    //D->stage_centroid_backup    = get_python_array_data_from_dict(quantities, "stage",      "centroid_backup");
    //D->xmom_centroid_backup     = get_python_array_data_from_dict(quantities, "xmomentum",  "centroid_backup");
    //D->ymom_centroid_backup     = get_python_array_data_from_dict(quantities, "ymomentum",  "centroid_backup");
    SAFE_MALLOC( D->stage_centroid_backup, D->number_of_elements, double);
    SAFE_MALLOC( D->xmom_centroid_backup, D->number_of_elements, double);
    SAFE_MALLOC( D->ymom_centroid_backup, D->number_of_elements, double);



    // Vertex values
    D->stage_vertex_values      = get_python_array_data_from_dict(quantities, "stage",      "vertex_values");
    D->xmom_vertex_values       = get_python_array_data_from_dict(quantities, "xmomentum",  "vertex_values");
    D->ymom_vertex_values       = get_python_array_data_from_dict(quantities, "ymomentum",  "vertex_values");
    D->bed_vertex_values        = get_python_array_data_from_dict(quantities, "elevation",  "vertex_values");
    D->height_vertex_values     = get_python_array_data_from_dict(quantities, "height",     "vertex_values");
    D->xvelocity_vertex_values  = get_python_array_data_from_dict(quantities, "xvelocity",  "vertex_values");
    D->yvelocity_vertex_values  = get_python_array_data_from_dict(quantities, "yvelocity",  "vertex_values");



    // Boundary values    
    D->stage_boundary_values = get_python_array_data_from_dict(quantities, "stage",     "boundary_values");
    D->xmom_boundary_values  = get_python_array_data_from_dict(quantities, "xmomentum", "boundary_values");
    D->ymom_boundary_values  = get_python_array_data_from_dict(quantities, "ymomentum", "boundary_values");
    D->bed_boundary_values   = get_python_array_data_from_dict(quantities, "elevation", "boundary_values");
    D->height_boundary_values   = get_python_array_data_from_dict(quantities, "height",     "boundary_values");
    D->xvelocity_boundary_values= get_python_array_data_from_dict(quantities, "xvelocity",  "boundary_values");
    D->yvelocity_boundary_values= get_python_array_data_from_dict(quantities, "yvelocity",  "boundary_values");



    // Explicit update values
    D->stage_explicit_update = get_python_array_data_from_dict(quantities, "stage",     "explicit_update");
    D->xmom_explicit_update  = get_python_array_data_from_dict(quantities, "xmomentum", "explicit_update");
    D->ymom_explicit_update  = get_python_array_data_from_dict(quantities, "ymomentum", "explicit_update");



    // Semi_implicit update values 
    D->stage_semi_implicit_update = get_python_array_data_from_dict(quantities, "stage",     "semi_implicit_update");
    D->xmom_semi_implicit_update  = get_python_array_data_from_dict(quantities, "xmomentum", "semi_implicit_update");
    D->ymom_semi_implicit_update  = get_python_array_data_from_dict(quantities, "ymomentum", "semi_implicit_update");



    // Gradient values
    D->stage_x_gradient = get_python_array_data_from_dict(quantities, "stage",     "x_gradient");
    D->xmom_x_gradient  = get_python_array_data_from_dict(quantities, "xmomentum", "x_gradient");
    D->ymom_x_gradient  = get_python_array_data_from_dict(quantities, "ymomentum", "x_gradient");
    D->height_x_gradient   = get_python_array_data_from_dict(quantities, "height",     "x_gradient");
    D->xvelocity_x_gradient= get_python_array_data_from_dict(quantities, "xvelocity",  "x_gradient");
    D->yvelocity_x_gradient= get_python_array_data_from_dict(quantities, "yvelocity",  "x_gradient");

    D->stage_y_gradient = get_python_array_data_from_dict(quantities, "stage",     "y_gradient");
    D->xmom_y_gradient  = get_python_array_data_from_dict(quantities, "xmomentum", "y_gradient");
    D->ymom_y_gradient  = get_python_array_data_from_dict(quantities, "ymomentum", "y_gradient");
    D->height_y_gradient   = get_python_array_data_from_dict(quantities, "height",     "y_gradient");
    D->xvelocity_y_gradient= get_python_array_data_from_dict(quantities, "xvelocity",  "y_gradient");
    D->yvelocity_y_gradient= get_python_array_data_from_dict(quantities, "yvelocity",  "y_gradient");
    
    
    return D;
}



int print_domain_struct(struct domain *D) {
    printf("D->number_of_elements     %ld  \n", D->number_of_elements);
    printf("D->number_of_boundary_elements  %ld  \n", D->number_of_boundary_elements);
    printf("D->epsilon                %g \n", D->epsilon);
    printf("D->H0                     %g \n", D->H0);
    printf("D->h0                     %g \n", D->h0);
    printf("D->limiting_threshold     %g \n", D->limiting_threshold);
    printf("D->g                      %g \n", D->g);
    printf("D->optimise_dry_cells     %ld \n", D->optimise_dry_cells);
    printf("D->evolve_max_timestep    %g \n", D->evolve_max_timestep);
    printf("D->minimum_allowed_height %g \n", D->minimum_allowed_height);
    printf("D->extrapolate_velocity_second_order %ld \n", D->extrapolate_velocity_second_order);
    printf("D->beta_w                 %g \n", D->beta_w);
    printf("D->beta_w_dry             %g \n", D->beta_w_dry);
    printf("D->beta_uh                %g \n", D->beta_uh);
    printf("D->beta_uh_dry            %g \n", D->beta_uh_dry);
    printf("D->beta_vh                %g \n", D->beta_vh);
    printf("D->beta_vh_dry            %g \n", D->beta_vh_dry);



    printf("D->neighbours             %p \n", D->neighbours);
    printf("D->surrogate_neighbours   %p \n", D->surrogate_neighbours);
    printf("D->neighbour_edges        %p \n", D->neighbour_edges);
    printf("D->normals                %p \n", D->normals);
    printf("D->edgelengths            %p \n", D->edgelengths);
    printf("D->radii                  %p \n", D->radii);
    printf("D->areas                  %p \n", D->areas);
    printf("D->tri_full_flag          %p \n", D->tri_full_flag);
    printf("D->already_computed_flux  %p \n", D->already_computed_flux);
    printf("D->vertex_coordinates     %p \n", D->vertex_coordinates);
    printf("D->edge_coordinates       %p \n", D->edge_coordinates);
    printf("D->centroid_coordinates   %p \n", D->centroid_coordinates);
    printf("D->max_speed              %p \n", D->max_speed);
    printf("D->number_of_boundaries   %p \n", D->number_of_boundaries);
    printf("D->stage_edge_values      %p \n", D->stage_edge_values);
    printf("D->xmom_edge_values       %p \n", D->xmom_edge_values);
    printf("D->ymom_edge_values       %p \n", D->ymom_edge_values);
    printf("D->bed_edge_values        %p \n", D->bed_edge_values);
    printf("D->stage_centroid_values  %p \n", D->stage_centroid_values);
    printf("D->xmom_centroid_values   %p \n", D->xmom_centroid_values);
    printf("D->ymom_centroid_values   %p \n", D->ymom_centroid_values);
    printf("D->bed_centroid_values    %p \n", D->bed_centroid_values);
    printf("D->stage_vertex_values    %p \n", D->stage_vertex_values);
    printf("D->xmom_vertex_values     %p \n", D->xmom_vertex_values);
    printf("D->ymom_vertex_values     %p \n", D->ymom_vertex_values);
    printf("D->bed_vertex_values      %p \n", D->bed_vertex_values);
    printf("D->stage_boundary_values  %p \n", D->stage_boundary_values);
    printf("D->xmom_boundary_values   %p \n", D->xmom_boundary_values);
    printf("D->ymom_boundary_values   %p \n", D->ymom_boundary_values);
    printf("D->bed_boundary_values    %p \n", D->bed_boundary_values);
    printf("D->stage_explicit_update  %p \n", D->stage_explicit_update);
    printf("D->xmom_explicit_update   %p \n", D->xmom_explicit_update);
    printf("D->ymom_explicit_update   %p \n", D->ymom_explicit_update);


    return 0;
}


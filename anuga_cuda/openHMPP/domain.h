// C struct for domain and quantities
//
// Stephen Roberts 2012
#include <stdio.h>
#include <stdlib.h>


#define SAFE_MALLOC(a, n, type) \
    (assert((a= (type *) malloc ((n) * sizeof (type))) != NULL))

/*
#define MALLOC(n, type) \
    ((type *) malloc ((n) * sizeof (type)))

#define SAFECALL(a, b) \
    (assert ( (a) != (b) ))
*/


// structures
struct domain {
    // Changing these don't change the data in python object
    long    number_of_elements;
    long    number_of_boundary_elements;
    long    number_of_tags;
    double  epsilon;
    double  H0;
    double  g;
    long    optimise_dry_cells;
    double  evolve_max_timestep;
    long    extrapolate_velocity_second_order;
    double  minimum_allowed_height;

    double beta_w;
    double beta_w_dry;
    double beta_uh;
    double beta_uh_dry;
    double beta_vh;
    double beta_vh_dry;

    // Compute fluxes method
    //char * compute_fluxes_method;
    int compute_fluxes_method;

    // Changing values in these arrays will change the values in the python object
    long*   neighbours;
    long*   neighbour_edges;
    long*   surrogate_neighbours;
    double* normals;
    double* edgelengths;
    double* radii;
    double* areas;



    long*   tri_full_flag;
    long*   already_computed_flux; // useless to GPU
    double* max_speed;

    double* vertex_coordinates;
    double* edge_coordinates;
    double* centroid_coordinates;

    // Edge values
    long*   number_of_boundaries;
    double* stage_edge_values;
    double* xmom_edge_values;
    double* ymom_edge_values;
    double* bed_edge_values;
    double* height_edge_values;
    double* xvelocity_edge_values;
    double* yvelocity_edge_values;


    // Centroid values
    double* stage_centroid_values;
    double* xmom_centroid_values;
    double* ymom_centroid_values;
    double* bed_centroid_values;
    double* friction_centroid_values;
    double* height_centroid_values;
    double* xvelocity_centroid_values;
    double* yvelocity_centroid_values;
    // Store values
    double* stage_centroid_store;
    double* xmom_centroid_store;
    double* ymom_centroid_store;
    // Backup values
    double* stage_centroid_backup;
    double* xmom_centroid_backup;
    double* ymom_centroid_backup;


    // Vertex values
    double* stage_vertex_values;
    double* xmom_vertex_values;
    double* ymom_vertex_values;
    double* bed_vertex_values;
    double* height_vertex_values;
    double* xvelocity_vertex_values;
    double* yvelocity_vertex_values;


    // Explicit update values
    double* stage_explicit_update;
    double* xmom_explicit_update;
    double* ymom_explicit_update;


    // Semi_implicit_update
    double* stage_semi_implicit;
    double* xmom_semi_implicit;
    double* ymom_semi_implicit;

    
    // Gradient values
    double* stage_x_gradient;
    double* xmom_x_gradient;
    double* ymom_x_gradient;
    double* height_x_gradient;
    double* xvelocity_x_gradient;
    double* yvelocity_x_gradient;

    double* stage_y_gradient;
    double* xmom_y_gradient;
    double* ymom_y_gradient;
    double* height_y_gradient;
    double* xvelocity_y_gradient;
    double* yvelocity_y_gradient;

    
    // Some others
    double* min_bed_edge_values;
    double* max_bed_edge_values;
    int * count_wet_neighbours;


    // Boundary values
    double* stage_boundary_values;
    double* xmom_boundary_values;
    double* ymom_boundary_values;
    double* bed_boundary_values;
    double* height_boundary_values;
    double* xvelocity_boundary_values;
    double* yvelocity_boundary_values;


    
    // Tags
    long* tags;
    long** boundary_ids;
    long** boundary_cells;
    long** boundary_edges;
};

    
struct edge {
    int cell_id;
    int edge_id;

    // mid point values
    double w;
    double h;
    double z;
    double uh;
    double vh;
    double u;
    double v;

    // vertex values
    double w1;
    double h1;
    double z1;
    double uh1;
    double vh1;
    double u1;
    double v1;

    double w2;
    double h2;
    double z2;
    double uh2;
    double vh2;
    double u2;
    double v2;
};



struct domain * get_domain_from_file( struct domain * sw_domain, char * file_name)
{
    int i, j;
    freopen(file_name, "r", stdin);
    
    scanf("%ld %ld %ld %lf %lf %lf %ld %lf %ld %lf", 
        &sw_domain->number_of_elements, 
        &sw_domain->number_of_boundary_elements,
        &sw_domain->number_of_tags,
        &sw_domain->epsilon, &sw_domain->H0, &sw_domain->g,
        &sw_domain->optimise_dry_cells,
        &sw_domain->evolve_max_timestep,
        &sw_domain->extrapolate_velocity_second_order,
        &sw_domain->minimum_allowed_height);


    scanf("%lf %lf %lf %lf %lf %lf",
        &sw_domain->beta_w,
        &sw_domain->beta_w_dry,
        &sw_domain->beta_uh,
        &sw_domain->beta_uh_dry,
        &sw_domain->beta_vh,
        &sw_domain->beta_vh_dry);
    
    // FIX: maybe better to use number 
    SAFE_MALLOC( sw_domain->compute_fluxes_method, 20, char);
    scanf("%s", sw_domain->compute_fluxes_method);


    i = sw_domain->number_of_elements;
    SAFE_MALLOC( sw_domain->neighbours, i*3, long);
    SAFE_MALLOC( sw_domain->neighbour_edges, i*3, long);
    SAFE_MALLOC( sw_domain->surrogate_neighbours, i*3, long);

    SAFE_MALLOC( sw_domain->normals, i*6, double);
    SAFE_MALLOC( sw_domain->edgelengths, i*3, double);
    SAFE_MALLOC( sw_domain->radii, i, double);
    SAFE_MALLOC( sw_domain->areas, i, double);

    SAFE_MALLOC( sw_domain->tri_full_flag, i, long);
    SAFE_MALLOC( sw_domain->max_speed, i, double);
    

    // Coordinates values
    SAFE_MALLOC( sw_domain->vertex_coordinates, i*6, double);
    SAFE_MALLOC( sw_domain->edge_coordinates, i*6, double);
    SAFE_MALLOC( sw_domain->centroid_coordinates, i*2, double);

    
    SAFE_MALLOC( sw_domain->number_of_boundaries, i, long);


    // Edge values
    SAFE_MALLOC( sw_domain->stage_edge_values, i*3, double);
    SAFE_MALLOC( sw_domain->xmom_edge_values, i*3, double);
    SAFE_MALLOC( sw_domain->ymom_edge_values, i*3, double);
    SAFE_MALLOC( sw_domain->bed_edge_values, i*3, double);
    SAFE_MALLOC( sw_domain->height_edge_values, i*3, double);
    SAFE_MALLOC( sw_domain->xvelocity_edge_values, i*3, double);
    SAFE_MALLOC( sw_domain->yvelocity_edge_values, i*3, double);
    

    // Centroid values
    SAFE_MALLOC( sw_domain->stage_centroid_values, i, double);
    SAFE_MALLOC( sw_domain->xmom_centroid_values, i, double);
    SAFE_MALLOC( sw_domain->ymom_centroid_values, i, double);
    SAFE_MALLOC( sw_domain->bed_centroid_values, i, double);
    SAFE_MALLOC( sw_domain->friction_centroid_values, i, double);
    SAFE_MALLOC( sw_domain->height_centroid_values, i, double);
    SAFE_MALLOC( sw_domain->xvelocity_centroid_values, i, double);
    SAFE_MALLOC( sw_domain->yvelocity_centroid_values, i, double);


    // Store values ?? backup values as well?
    SAFE_MALLOC( sw_domain->stage_centroid_store, i, double);
    SAFE_MALLOC( sw_domain->xmom_centroid_store, i, double);
    SAFE_MALLOC( sw_domain->ymom_centroid_store, i, double);

    SAFE_MALLOC( sw_domain->stage_centroid_backup, i, double);
    SAFE_MALLOC( sw_domain->xmom_centroid_backup, i, double);
    SAFE_MALLOC( sw_domain->ymom_centroid_backup, i, double);


    // Vertex values
    SAFE_MALLOC( sw_domain->stage_vertex_values, i*3, double);
    SAFE_MALLOC( sw_domain->xmom_vertex_values, i*3, double);
    SAFE_MALLOC( sw_domain->ymom_vertex_values, i*3, double);
    SAFE_MALLOC( sw_domain->bed_vertex_values, i*3, double);
    SAFE_MALLOC( sw_domain->height_vertex_values, i*3, double);
    SAFE_MALLOC( sw_domain->xvelocity_vertex_values, i*3, double);
    SAFE_MALLOC( sw_domain->yvelocity_vertex_values, i*3, double);


    // Explicit update values
    SAFE_MALLOC( sw_domain->stage_explicit_update, i, double);
    SAFE_MALLOC( sw_domain->xmom_explicit_update, i, double);
    SAFE_MALLOC( sw_domain->ymom_explicit_update, i, double);


    // Semi_implicit_update
    SAFE_MALLOC( sw_domain->stage_semi_implicit, i, double);
    SAFE_MALLOC( sw_domain->xmom_semi_implicit, i, double);
    SAFE_MALLOC( sw_domain->ymom_semi_implicit, i, double);

    // Gradient values
    SAFE_MALLOC( sw_domain->stage_x_gradient, i, double);
    SAFE_MALLOC( sw_domain->xmom_x_gradient, i, double);
    SAFE_MALLOC( sw_domain->ymom_x_gradient, i, double);
    SAFE_MALLOC( sw_domain->height_x_gradient, i, double);
    SAFE_MALLOC( sw_domain->xvelocity_x_gradient, i, double);
    SAFE_MALLOC( sw_domain->yvelocity_x_gradient, i, double);

    SAFE_MALLOC( sw_domain->stage_y_gradient, i, double);
    SAFE_MALLOC( sw_domain->xmom_y_gradient, i, double);
    SAFE_MALLOC( sw_domain->ymom_y_gradient, i, double);
    SAFE_MALLOC( sw_domain->height_y_gradient, i, double);
    SAFE_MALLOC( sw_domain->xvelocity_y_gradient, i, double);
    SAFE_MALLOC( sw_domain->yvelocity_y_gradient, i, double);


    // Some others
    SAFE_MALLOC( sw_domain->min_bed_edge_values, i, double);
    SAFE_MALLOC( sw_domain->max_bed_edge_values, i, double);
    SAFE_MALLOC( sw_domain->count_wet_neighbouts, i, int);

    
    // Boundary values
    i = sw_domain->number_of_boundary_elements;
    SAFE_MALLOC( sw_domain->stage_boundary_values, i, double);
    SAFE_MALLOC( sw_domain->xmom_boundary_values, i, double);
    SAFE_MALLOC( sw_domain->ymom_boundary_values, i, double);
    SAFE_MALLOC( sw_domain->bed_boundary_values, i, double);
    SAFE_MALLOC( sw_domain->height_boundary_values, i, double);
    SAFE_MALLOC( sw_domain->xvelocity_boundary_values, i, double);
    SAFE_MALLOC( sw_domain->yvelocity_boundary_values, i, double);


    
    // Storing the number of elements for each tags
    SAFE_MALLOC( sw_domain->tags, sw_domain->number_of_tags, long);
    SAFE_MALLOC( sw_domain->boundary_ids, sw_domain->number_of_tags, long *);
    SAFE_MALLOC( sw_domain->boundary_cells, sw_domain->number_of_tags, long *);
    SAFE_MALLOC( sw_domain->boundary_edges, sw_domain->number_of_tags, long *);

    for (i=0; i< sw_domain->number_of_tags; i++)
    {
        scanf("%ld", sw_domain->tags+i);
        SAFE_MALLOC( sw_domain->boundary_ids[i], sw_domain->tags[i], long);
        SAFE_MALLOC( sw_domain->boundary_cells[i], sw_domain->tags[i], long);
        SAFE_MALLOC( sw_domain->boundary_edges[i], sw_domain->tags[i], long);
        
        for (j=0; j < sw_domain->tags[i]; i++)
        {
            scanf("%ld %ld %ld", 
                sw_domain->boundary_ids[i]+j,
                sw_domain->boundary_cells[i]+j,
                sw_domain->boundary_edges[i]+j);
        }
    }


    for (i=0; i< sw_domain->number_of_elements; i++)
    {
        
    }
}

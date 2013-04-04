// C struct for domain and quantities
//
// Stephen Roberts 2012




// structures
struct domain {
    // Changing these don't change the data in python object
    long    number_of_elements;
    long    number_of_boundarie_elements;
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

    // 
    char * compute_fluxes_method;

    // Changing values in these arrays will change the values in the python object
    long*   neighbours;
    long*   neighbour_edges;
    long*   surrogate_neighbours;
    double* normals;
    double* edgelengths;
    double* radii;
    double* areas;



    long*   tri_full_flag;
    long*   already_computed_flux;
    double* max_speed;

    double* vertex_coordinates;
    double* edge_coordinates;
    double* centroid_coordinates;

    long*   number_of_boundaries;
    double* stage_edge_values;
    double* xmom_edge_values;
    double* ymom_edge_values;
    double* bed_edge_values;

    double* stage_centroid_values;
    double* xmom_centroid_values;
    double* ymom_centroid_values;
    double* bed_centroid_values;

    double* stage_vertex_values;
    double* xmom_vertex_values;
    double* ymom_vertex_values;
    double* bed_vertex_values;


    double* stage_boundary_values;
    double* xmom_boundary_values;
    double* ymom_boundary_values;
    double* bed_boundary_values;

    double* stage_explicit_update;
    double* xmom_explicit_update;
    double* ymom_explicit_update;
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



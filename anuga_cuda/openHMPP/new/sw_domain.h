// C struct for domain and quantities
//
// Stephen Roberts 2012
// John Weng 2013

#define SW_DOMAIN

// structures
struct domain {
    // Changing these don't change the data in python object
    long    number_of_elements;
    long    number_of_boundary_elements;
    double  epsilon;
    double  H0;
    double  h0;
    double  limiting_threshold;
    double  g;
    long    optimise_dry_cells;
    double  evolve_max_timestep;
    double  evolve_min_timestep;
    long    extrapolate_velocity_second_order;
    double  minimum_allowed_height;
    double  time;
    double  starttime;
    double  finaltime;
    double  yieldtime;
    double  timestep;

    int     smallsteps;
    int     max_smallsteps;
    int     number_of_steps;
    int     number_of_first_order_steps;
    int     timestepping_method;
    int     _order_;
    int     default_order;
    int     use_sloped_mannings;
    int     use_centroid_velocities;
    int     use_edge_limiter;

    int     flow_algorithm;
    int     compute_fluxes_method;

    long*   boundary_cells;
    long*   boundary_edges;


    double CFL;
    double flux_timestep;
    double recorded_max_timestep;
    double recorded_min_timestep;
    double maximum_allowed_speed;
    double optimised_gradient_limiter;
    double alpha_balance;
    double tight_slope_limiters;

    double beta_w;
    double beta_w_dry;
    double beta_uh;
    double beta_uh_dry;
    double beta_vh;
    double beta_vh_dry;

    double stage_beta;
    double xmom_beta;
    double ymom_beta;

    
    // long or int ???
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
    double* timestep_array;

    double* vertex_coordinates;
    double* edge_coordinates;
    double* centroid_coordinates;

    long*   number_of_boundaries;

    // Edge Values
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


    // Semi_implicit update values
    double* stage_semi_implicit_update;
    double* xmom_semi_implicit_update;
    double* ymom_semi_implicit_update;


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
    int*    count_wet_neighbours;

    // Boundary values
    double* stage_boundary_values;
    double* xmom_boundary_values;
    double* ymom_boundary_values;
    double* bed_boundary_values;
    double* height_boundary_values;
    double* xvelocity_boundary_values;
    double* yvelocity_boundary_values;

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

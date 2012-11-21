__global__ void _extrapolate_second_order_sw_old(int number_of_elements,
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
        double* xmom_centroid_values,
        double* ymom_centroid_values,
        double* elevation_centroid_values,
        double* vertex_coordinates,
        double* stage_vertex_values,
        double* xmom_vertex_values,
        double* ymom_vertex_values,
        double* elevation_vertex_values,
        int optimise_dry_cells,
        int extrapolate_velocity_second_order);
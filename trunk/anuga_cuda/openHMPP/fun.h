
void extrapolate_second_order_and_limit_by_vertex(
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

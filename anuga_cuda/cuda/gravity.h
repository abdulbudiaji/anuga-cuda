__global__ void gravity(double *bed_vertex_values, double *stage_centroid_values, double *bed_centroid_values, double *vertex_coordinates, double * xmom_explicit_update, double *ymom_explicit_update, float g);

__global__ void gravity_wb(double *stage_vertex_values, double *stage_edge_values, double *stage_centroid_values, double *bed_centroid_values, double *bed_edge_values,  double *vertex_coordinates, double * xmom_explicit_update, double *ymom_explicit_update, double * normals, double *areas, double * edgelengths,float g);

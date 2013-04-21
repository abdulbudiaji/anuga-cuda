//#define REARRANGED_DOMAIN

void interpolate_from_vertices_to_edges(
        int N,
        double* vertex_values,
        double* edge_values) 
{
    int k;

#ifndef REARRANGED_DOMAIN
    int k3;
#endif
    double q0, q1, q2;

    for (k=0; k<N; k++) {

#ifndef REARRANGED_DOMAIN
        k3 = k*3;
        q0 = vertex_values[k3 + 0];
        q1 = vertex_values[k3 + 1];
        q2 = vertex_values[k3 + 2];

        edge_values[k3 + 0] = 0.5*(q1+q2);
        edge_values[k3 + 1] = 0.5*(q0+q2);
        edge_values[k3 + 2] = 0.5*(q0+q1);
#else
        q0 = vertex_values[k];
        q1 = vertex_values[k + N];
        q2 = vertex_values[k + 2*N];

        edge_values[k] = 0.5*(q1+q2);
        edge_values[k + N] = 0.5*(q0+q2);
        edge_values[k + 2*N] = 0.5*(q0+q1);
#endif
    }
}
#include "hmpp_fun.h"


void interpolate_from_vertices_to_edges(
        int N,
        int N3,
        double vertex_values[N3],
        double edge_values[N3]) 
{
    int k;

    double q0, q1, q2;

    for (k=0; k<N; k++) {

#ifndef REARRANGED_DOMAIN
        q0 = vertex_values[k*3 + 0];
        q1 = vertex_values[k*3 + 1];
        q2 = vertex_values[k*3 + 2];

        edge_values[k*3 + 0] = 0.5*(q1+q2);
        edge_values[k*3 + 1] = 0.5*(q0+q2);
        edge_values[k*3 + 2] = 0.5*(q0+q1);
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

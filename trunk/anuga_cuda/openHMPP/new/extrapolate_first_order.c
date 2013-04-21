
#include "hmpp_fun.h"


void extrapolate_first_order(
        int N,
        int N3,
        double centroid_values[N],
        double edge_values[N3],
        double vertex_values[N3])
{
    int k;

    for (k=0; k < N; k++)
    {
#ifndef REARRANGED_DOMAIN
        edge_values[k*3] = centroid_values[k];
        edge_values[k*3+ 1] = centroid_values[k];
        edge_values[k*3+ 2] = centroid_values[k];

        vertex_values[k*3] = centroid_values[k];
        vertex_values[k*3+ 1] = centroid_values[k];
        vertex_values[k*3+ 2] = centroid_values[k];
#else
        edge_values[k] = centroid_values[k];
        edge_values[k + N] = centroid_values[k];
        edge_values[k + 2*N] = centroid_values[k];

        vertex_values[k] = centroid_values[k];
        vertex_values[k + N] = centroid_values[k];
        vertex_values[k + 2*N] = centroid_values[k];
#endif
    }
}

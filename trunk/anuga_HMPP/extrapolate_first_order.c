
#include "hmpp_fun.h"


#ifdef USING_LOCAL_DIRECTIVES
#pragma hmpp extraFstOrder codelet, target=CUDA args[*].transfer=atcall
#endif
void extrapolate_first_order(
        int N,
        int N3,
        double centroid_values[N],
        double edge_values[N3],
        double vertex_values[N3])
{
    int k, k3;

    #pragma hmppcg gridify(k), &
    #pragma hmppcg & private(k3), &
    #pragma hmppcg & global( centroid_values, edge_values, vertex_values)
    for (k=0; k < N; k++)
    {
#ifndef REARRANGED_DOMAIN
        k3 = k*3;
        edge_values[k3] = centroid_values[k];
        edge_values[k3+ 1] = centroid_values[k];
        edge_values[k3+ 2] = centroid_values[k];

        vertex_values[k3] = centroid_values[k];
        vertex_values[k3+ 1] = centroid_values[k];
        vertex_values[k3+ 2] = centroid_values[k];
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

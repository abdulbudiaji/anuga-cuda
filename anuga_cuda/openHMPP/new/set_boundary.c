#include "hmpp_fun.h"


void set_boundary_values_from_edges(
        int N,
        int N3,
        long vol_id[N],
        long edge_id[N],
        double boundary_values[N],
        double edge_values[N3])
{
    int k, id;
    
    for(k=0; k<N; k++)
    {
        id= 3*vol_id[k] + edge_id[k];

        boundary_values[k] = edge_values[id];
    }
}

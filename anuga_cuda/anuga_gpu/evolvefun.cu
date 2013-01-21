#include "evolvefun.h"

__host__ void distribute_to_vertices_and_edges(struct D *D)
{
    //FIXME:can not find swb2_D_ext module
    if (D->compute_fluxes_method[0] == 't') //tsunami
    {}
    else if (D->use_edge_limiter)
        distribute_using_edge_limiter(D);
    else
        distribute_using_vertex_limiter(D);
}

__host__ void distribute_using_edge_limiter(struct D *D)
{
    
}

__host__ void distribute_using_vertex_limiter(struct D *D)
{
    
}

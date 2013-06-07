//#define REARRANGED_DOMAIN
#include "hmpp_fun.h"



#ifdef USING_LOCAL_DIRECTIVES
#pragma hmpp evaRef codelet, target=CUDA args[*].transfer=atcall
#endif
void evaluate_segment_reflective(
    int N1, // Nids
    int N2,
    int N3,
    int N6,

    long ids[N1],
    long vol_ids[N2],   // domain.boundary_cells
    long edge_ids[N2],  // domain.boundary_edges
    double normals[N6], 
    
    double stage_edge_values[N3],
    double bed_edge_values[N3],
    double height_edge_values[N3],
    double xmom_edge_values[N3],
    double ymom_edge_values[N3],
    double xvel_edge_values[N3],
    double yvel_edge_values[N3],

    double stage_boundary_values[N2],
    double bed_boundary_values[N2],
    double height_boundary_values[N2],
    double xmom_boundary_values[N2],
    double ymom_boundary_values[N2],
    double xvel_boundary_values[N2],
    double yvel_boundary_values[N2]
    )
{
    int k, id, id_vol, id_edge;

    double n1, n2, q1, q2, r1, r2;


    #pragma hmppcg gridify(k), &
    #pragma hmppcg & private(id, id_vol, id_edge, n1, n2, q1, q2, r1, r2), &
    #pragma hmppcg & global( edge_ids, normals, &
    #pragma hmppcg & stage_edge_values, bed_edge_values, &
    #pragma hmppcg & height_edge_values, xmom_edge_values, &
    #pragma hmppcg & ymom_edge_values, xvel_edge_values, &
    #pragma hmppcg & yvel_edge_values, stage_boundary_values, &
    #pragma hmppcg & bed_boundary_values, height_boundary_values, &
    #pragma hmppcg & xmom_boundary_values, ymom_boundary_values, &
    #pragma hmppcg & xvel_boundary_values, yvel_boundary_values)
    for(k=0; k<N1; k++)
    {
        //id = ids[k];
        id = k;
        id_vol = vol_ids[ id ];
        id_edge= edge_ids[ id ];

#ifndef REARRANGED_DOMAIN
        n1 = normals[id_vol*6 + id_edge*2];
        n2 = normals[id_vol*6 + id_edge*2 + 1];
        q1 = xmom_edge_values[id_vol*3 + id_edge];
        q2 = ymom_edge_values[id_vol*3 + id_edge];
#else
        n1 = normals[id_vol + id_edge*2*N];
        n2 = normals[id_vol + (id_edge*2+1)*N];
        q1 = xmom_edge_values[id_vol + id_edge*N];
        q2 = ymom_edge_values[id_vol + id_edge*N];
#endif

        r1 = -q1*n1 - q2*n2;
        r2 = -q1*n2 + q2*n1;

#ifndef REARRANGED_DOMAIN
        stage_boundary_values[id] = stage_edge_values[id_vol*3 + id_edge];
        bed_boundary_values[id] = bed_edge_values[id_vol*3 + id_edge];
        height_boundary_values[id] = height_edge_values[id_vol*3 + id_edge];

        q1 = xvel_edge_values[id_vol*3 + id_edge];
        q2 = yvel_edge_values[id_vol*3 + id_edge];
#else
        stage_boundary_values[id] = stage_edge_values[id_vol + id_edge*N];
        bed_boundary_values[id] = bed_edge_values[id_vol + id_edge*N];
        height_boundary_values[id] = height_edge_values[id_vol + id_edge*N];

        q1 = xvel_edge_values[id_vol + id_edge*N];
        q2 = yvel_edge_values[id_vol + id_edge*N];
#endif

        xmom_boundary_values[id] = n1*r1 - n2*r2;
        ymom_boundary_values[id] = n2*r1 + n1*r2;


        r1 = q1*n1 + q2*n2;
        r2 = q1*n2 - q2*n1;

        xvel_boundary_values[id] = n1*r1 - n2*r2;
        yvel_boundary_values[id] = n2*r1 + n1*r2;
    
        vol_ids[k] = ids[k];
    }
}


/*
void evaluate_segment_dirichlet_1(
    int N,
    long * ids,
    long * vol_ids,
    long * edge_ids,
    double * boundary_values,
    double *edge_values)
{
    const int k = 
        threadIdx.x+threadIdx.y*blockDim.x+
        (blockIdx.x+blockIdx.y*gridDim.x)*blockDim.x*blockDim.y;
    
    if (k >= N)
        return;

    int id = ids[k],
        id_vol = vol_ids[k],
        id_edge = edge_ids[k];

    boundary_values[id] = edge_values[id_vol*3 + id_edge];
}


void evaluate_segment_dirichlet_2(
    int N,
    double  q_bdry,
    long * ids,
    double * boundary_values)
{
    const int k = 
        threadIdx.x+threadIdx.y*blockDim.x+
        (blockIdx.x+blockIdx.y*gridDim.x)*blockDim.x*blockDim.y;
    
    int id = ids[k];
    
    if (k >= N)
        return;

    boundary_values[id] = q_bdry;
}

*/

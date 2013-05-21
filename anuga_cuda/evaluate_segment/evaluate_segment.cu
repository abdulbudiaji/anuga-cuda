#define REARRANGED_DOMAIN

__global__ void evaluate_segment_reflective(
    int N,
    int Nids,
    long * ids,
    long * vol_ids,
    long * edge_ids,
    double * normals,
    
    double * stage_edge_values,
    double * bed_edge_values,
    double * height_edge_values,
    double * xmom_edge_values,
    double * ymom_edge_values,
    double * xvel_edge_values,
    double * yvel_edge_values,

    double * stage_boundary_values,
    double * bed_boundary_values,
    double * height_boundary_values,
    double * xmom_boundary_values,
    double * ymom_boundary_values,
    double * xvel_boundary_values,
    double * yvel_boundary_values)
{
    const int k = 
        threadIdx.x+threadIdx.y*blockDim.x+
        (blockIdx.x+blockIdx.y*gridDim.x)*blockDim.x*blockDim.y;
    
    if (k >= Nids)
        return;

    long id = ids[k],
         id_vol = vol_ids[k],   // Note here should be k, since we pass
         id_edge= edge_ids[k];  // the vol_ids and edge_ids from CPU

#ifndef REARRANGED_DOMAIN
    double n1 = normals[id_vol*6 + id_edge*2],
           n2 = normals[id_vol*6 + id_edge*2 + 1];
    double q1 = xmom_edge_values[id_vol*3 + id_edge],
           q2 = ymom_edge_values[id_vol*3 + id_edge];
#else
    double n1 = normals[id_vol + id_edge*2*N],
           n2 = normals[id_vol + (id_edge*2+1)*N];
    double q1 = xmom_edge_values[id_vol + id_edge*N],
           q2 = ymom_edge_values[id_vol + id_edge*N];
#endif
    
    double r1 = -q1*n1 - q2*n2,
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
}

__global__ void evaluate_segment_dirichlet_1(
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


__global__ void evaluate_segment_dirichlet_2(
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

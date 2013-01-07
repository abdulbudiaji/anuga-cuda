__global__ void evaluate_segment(
    int N,
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
    
    if (k >= N)
        return;

    long id = ids[k],
         id_vol = vol_ids[k],
         id_edge= edge_ids[k];

    double n1 = normals[id_vol*6 + id_edge*2],
            n2 = normals[id_vol*6 + id_edge*2 + 1];
    
    double q1 = xmom_edge_values[id_vol*3 + id_edge],
            q2 = ymom_edge_values[id_vol*3 + id_edge];
    
    double r1 = -q1*n1 - q2*n2,
            r2 = -q1*n2 + q2*n1;

    stage_boundary_values[id] = stage_edge_values[id_vol*3 + id_edge];
    bed_boundary_values[id] = bed_edge_values[id_vol*3 + id_edge];
    height_boundary_values[id] = height_edge_values[id_vol*3 + id_edge];

    xmom_boundary_values[id] = n1*r1 - n2*r2;
    ymom_boundary_values[id] = n2*r1 + n1*r2;

    q1 = xvel_edge_values[id_vol*3 + id_edge];
    q2 = yvel_edge_values[id_vol*3 + id_edge];

    r1 = q1*n1 + q2*n2;
    r2 = q1*n2 - q2*n1;

    xvel_boundary_values[id] = n1*r1 - n2*r2;
    yvel_boundary_values[id] = n2*r1 + n1*r2;
}



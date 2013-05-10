//#define REARRANGED_DOMAIN 

#define BLOCK_SIZE 960

__global__ void gravity_wb(
        int N,
        double g,
        double * stage_vertex_values, 
        double * stage_edge_values, 
        double * stage_centroid_values, 
        double * bed_edge_values, 
        double * bed_centroid_values, 
        double * vertex_coordinates, 
        double * xmom_explicit_update, 
        double * ymom_explicit_update, 
        double * normals, 
        double * areas, 
        double * edgelengths
        )
{
    const int k = 
            threadIdx.x+threadIdx.y*blockDim.x+
            (blockIdx.x+blockIdx.y*gridDim.x)*blockDim.x*blockDim.y;

    int i;
#ifndef REARRANGED_DOMAIN
    int k3=k*3;
#endif

    double w0, w1, w2, 
           x0, y0, x1, y1, x2, y2,
           avg_h;

    double wx, wy, det,
           hh[3];
    //double sidex, sidey;
    double area, n0, n1, fact;

    __shared__ double sh_data[ BLOCK_SIZE *6];

    if (k >= N)
        return;
#ifndef REARRANGED_DOMAIN
    w0 = stage_vertex_values[k3];
    w1 = stage_vertex_values[k3 + 1];
    w2 = stage_vertex_values[k3 + 2];

    x0 = vertex_coordinates[k*6];
    y0 = vertex_coordinates[k*6 + 1];
    x1 = vertex_coordinates[k*6 + 2];
    y1 = vertex_coordinates[k*6 + 3];
    x2 = vertex_coordinates[k*6 + 4];
    y2 = vertex_coordinates[k*6 + 5];
#else
    w0 = stage_vertex_values[k];
    w1 = stage_vertex_values[k + N];
    w2 = stage_vertex_values[k + 2*N];
    
    x0 = vertex_coordinates[k];
    y0 = vertex_coordinates[k + N];
    x1 = vertex_coordinates[k + 2*N];
    y1 = vertex_coordinates[k + 3*N];
    x2 = vertex_coordinates[k + 4*N];
    y2 = vertex_coordinates[k + 5*N];
#endif

    //_gradient(x0, y0, x1, y1, x2, y2, w0, w1, w2, &wx, &wy);

    det = (y2 - y0)*(x1 - x0) - (y1 - y0)*(x2 - x0);

    wx = (y2 -y0)*(w1 - w0) - (y1 - y0)*(w2 -w0);
    wx /= det;

    wy = (x1 - x0)*(w2 - w0) - (x2 - x0)*(w1 -w0);
    wy /= det;


    avg_h = stage_centroid_values[k] - bed_centroid_values[k];

    xmom_explicit_update[k] += -g * wx * avg_h;
    ymom_explicit_update[k] += -g * wy * avg_h;

#ifndef REARRANGED_DOMAIN
    hh[0] = stage_edge_values[k3] - bed_edge_values[k3];
    hh[1] = stage_edge_values[k3+1] - bed_edge_values[k3+1];
    hh[2] = stage_edge_values[k3+2] - bed_edge_values[k3+2];
#else
    hh[0] = stage_edge_values[k] - bed_edge_values[k];
    hh[1] = stage_edge_values[k+N] - bed_edge_values[k+N];
    hh[2] = stage_edge_values[k+2*N] - bed_edge_values[k+2*N];
#endif

    //sidex = 0.0;
    //sidey = 0.0;
    area = areas[k];

    for ( i = 0 ; i < 3 ; i++ )
    {
#ifndef REARRANGED_DOMAIN
        n0 = normals[k*6 + 2*i];
        n1 = normals[k*6 + 2*i + 1];

        fact =  -0.5 * g * hh[i] * hh[i] * edgelengths[k3 + i];
#else
        n0 = normals[k + 2*i*N];
        n1 = normals[k + (2*i + 1)*N];

        fact =  -0.5 * g * hh[i] * hh[i] * edgelengths[k + i*N];
#endif

        //sidex += fact*n0;
        //sidey += fact*n1;

        sh_data[threadIdx.x + i*blockDim.x] = fact*n0;
        sh_data[threadIdx.x + (i+3)*blockDim.x] = fact*n1;
    }

    //xmom_explicit_update[k] += -sidex / area;
    //ymom_explicit_update[k] += -sidey / area;

    xmom_explicit_update[k] += -(sh_data[threadIdx.x] + sh_data[threadIdx.x + blockDim.x] + sh_data[threadIdx.x+2*blockDim.x]) / area;

    ymom_explicit_update[k] += -(sh_data[threadIdx.x+3*blockDim.x] + sh_data[threadIdx.x + 4*blockDim.x] + sh_data[threadIdx.x+5*blockDim.x]) / area;
}

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

    double w0, w1, w2, 
           x0, y0, x1, y1, x2, y2,
           avg_h;

    double wx, wy, det,
           hh[3];
    double sidex, sidey, area, n0, n1, fact;

    __shared__ double sh_data[32*6];

    if (k >= N)
        return;

    w0 = stage_vertex_values[3*k + 0];
    w1 = stage_vertex_values[3*k + 1];
    w2 = stage_vertex_values[3*k + 2];

    x0 = vertex_coordinates[k*6 + 0];
    y0 = vertex_coordinates[k*6 + 1];
    x1 = vertex_coordinates[k*6 + 2];
    y1 = vertex_coordinates[k*6 + 3];
    x2 = vertex_coordinates[k*6 + 4];
    y2 = vertex_coordinates[k*6 + 5];


    //_gradient(x0, y0, x1, y1, x2, y2, w0, w1, w2, &wx, &wy);

    det = (y2 - y0)*(x1 - x0) - (y1 - y0)*(x2 - x0);

    wx = (y2 -y0)*(w1 - w0) - (y1 - y0)*(w2 -w0);
    wx /= det;

    wy = (x1 - x0)*(w2 - w0) - (x2 - x0)*(w1 -w0);
    wy /= det;


    avg_h = stage_centroid_values[k] - bed_centroid_values[k];

    xmom_explicit_update[k] += -g * wx * avg_h;
    ymom_explicit_update[k] += -g * wy * avg_h;


    hh[0] = stage_edge_values[k*3] - bed_edge_values[k*3];
    hh[1] = stage_edge_values[k*3+1] - bed_edge_values[k*3+1];
    hh[2] = stage_edge_values[k*3+2] - bed_edge_values[k*3+2];

    sidex = 0.0;
    sidey = 0.0;
    area = areas[k];

    for ( i = 0 ; i < 3 ; i++ )
    {
        n0 = normals[k*6 + 2*i];
        n1 = normals[k*6 + 2*i + 1];

        fact =  -0.5 * g * hh[i] * hh[i] * edgelengths[k*3 + i];

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

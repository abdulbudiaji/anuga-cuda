#include "gravity.h"
#define    W = 32

 __global__ void gravity(double *bed_vertex_values, double *stage_centroid_values, double *bed_centroid_values, double *vertex_coordinates, double * xmom_explicit_update, double *ymom_explicit_update, float g)
{
    int i = threadIdx.x + threadIdx.y + blockIdx.x * 32*32;
    int k3 = 3*i,
        k6 = 6*i;
	double 
        z0 = bed_vertex_values[k3],
        z1 = bed_vertex_values[k3 + 1],
        z2 = bed_vertex_values[k3 + 2],
    
		avg_h = stage_centroid_values[i] - bed_centroid_values[i],
    
		x0 = vertex_coordinates[k6],
        y0 = vertex_coordinates[k6 + 1],
        x1 = vertex_coordinates[k6 + 2],
        y1 = vertex_coordinates[k6 + 3],
        x2 = vertex_coordinates[k6 + 4],
        y2 = vertex_coordinates[k6 + 5],
                
        zx, zy, det;
    
    det = (y2 - y0)*(x1 - x0) - (y1 - y0)*(x2 - x0);
    
    zx = (y2 -y0)*(z1 - z0) - (y1 - y0)*(z2 -z0);
    zx /= det;
    
    zy = (x1 - x0)*(z2 - z0) - (x2 - x0)*(z1 -z0);
    zy /= det;
    
    xmom_explicit_update[i] += -g * zx * avg_h;
    ymom_explicit_update[i] += -g * zy * avg_h;
}
        
__global__ void gravity_wb(double *stage_vertex_values, double *stage_edge_values, double *stage_centroid_values, double *bed_centroid_values, double *bed_edge_values,  double *vertex_coordinates, double * xmom_explicit_update, double *ymom_explicit_update, double * normals, double *areas, double * edgelengths,float g)
{            
	int k = threadIdx.x + threadIdx.y + blockIdx.x * 32*32;
    int k3 = 3*k,
        k6 = 6*k,
        i=0;
            
    double 
        w0 = stage_vertex_values[k3],
        w1 = stage_vertex_values[k3 + 1],
        w2 = stage_vertex_values[k3 + 2],
    
        avg_h = stage_centroid_values[k] - bed_centroid_values[k],
    
        x0 = vertex_coordinates[k6],
        y0 = vertex_coordinates[k6 + 1],
        x1 = vertex_coordinates[k6 + 2],
        y1 = vertex_coordinates[k6 + 3],
        x2 = vertex_coordinates[k6 + 4],
        y2 = vertex_coordinates[k6 + 5],
                    
        wx, wy, det,
        hh[3],
        sidex, sidey, area, n0, n1, fact;
               
    det = (y2 - y0)*(x1 - x0) - (y1 - y0)*(x2 - x0);
    
    wx = (y2 -y0)*(w1 - w0) - (y1 - y0)*(w2 -w0);
    wx /= det;
    
    wy = (x1 - x0)*(w2 - w0) - (x2 - x0)*(w1 -w0);
    wy /= det;
    
    xmom_explicit_update[k] += -g * wx * avg_h;
    ymom_explicit_update[k] += -g * wy * avg_h;
    
    hh[0] = stage_edge_values[k3] - bed_edge_values[k3];
    hh[1] = stage_edge_values[k3+1] - bed_edge_values[k3+1];
    hh[2] = stage_edge_values[k3+2] - bed_edge_values[k3+2];
            
    sidex = 0.0;
    sidey = 0.0;
            
    for ( i = 0 ; i < 3 ; i++ )
    {
        n0 = normals[k6 + 2*i];
        n1 = normals[k6 + 2*i + 1];
                
        fact =  -0.5 * g * hh[i] * hh[i] * edgelengths[k3 + i];
        sidex = sidex + fact*n0;
        sidey = sidey + fact*n1;
    }
            
    area = areas[k];
    xmom_explicit_update[k] += -sidex / area;
    ymom_explicit_update[k] += -sidey / area;
                
}

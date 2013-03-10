#include<stdio.h>
#include<stdlib.h>
#include<assert.h>

void gravity_wb(
        double *restrict xmom_explicit_update, 
        double *restrict ymom_explicit_update, 
        int N,
        double g,
        double * stage_vertex_values, 
        double * stage_edge_values, 
        double * stage_centroid_values, 
        double * bed_edge_values, 
        double * bed_centroid_values, 
        double * vertex_coordinates, 
        double * normals, 
        double * areas, 
        double * edgelengths
        )
{
    int i, k, k3, k6;

    double w0, w1, w2, 
           x0, y0, x1, y1, x2, y2,
           avg_h;

    double wx, wy, det,
           hh[3];
    double sidex, sidey, area, n0, n1, fact;

    #pragma acc kernels loop copyin(stage_vertex_values[0:3*N], stage_edge_values[0:3*N], &
    #pragma acc & stage_centroid_values[0:N], bed_edge_values[0:3*N], bed_centroid_values[0:N], &
    #pragma acc & vertex_coordinates[0:6*N], normals[0:6*N], area[0:N], edgelengths[0:3*N]) &
    #pragma acc copy(xmom_explicit_update[0:N], ymom_explicit_update[0:N])
    for (k = 0; k < N; k++)
    {
        k3 = 3*k;
        w0 = stage_vertex_values[3*k + 0];
        w1 = stage_vertex_values[3*k + 1];
        w2 = stage_vertex_values[3*k + 2];


        k6 = 6*k;
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

        for ( i = 0 ; i < 3 ; i++ )
        {
            n0 = normals[k*6 + 2*i];
            n1 = normals[k*6 + 2*i + 1];

            fact =  -0.5 * g * hh[i] * hh[i] * edgelengths[k*3 + i];

            sidex += fact*n0;
            sidey += fact*n1;
        }

        area = areas[k];
        xmom_explicit_update[k] += -sidex / area;
        ymom_explicit_update[k] += -sidey / area;

    }
}


void gravity_wb_orig(
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
    int i, k, k3, k6;

    double w0, w1, w2, 
           x0, y0, x1, y1, x2, y2,
           avg_h;

    double wx, wy, det,
           hh[3];
    double sidex, sidey, area, n0, n1, fact;

    for (k = 0; k < N; k++)
    {
        k3 = 3*k;
        w0 = stage_vertex_values[3*k + 0];
        w1 = stage_vertex_values[3*k + 1];
        w2 = stage_vertex_values[3*k + 2];


        k6 = 6*k;
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

        for ( i = 0 ; i < 3 ; i++ )
        {
            n0 = normals[k*6 + 2*i];
            n1 = normals[k*6 + 2*i + 1];

            fact =  -0.5 * g * hh[i] * hh[i] * edgelengths[k*3 + i];

            sidex += fact*n0;
            sidey += fact*n1;
        }

        area = areas[k];
        xmom_explicit_update[k] += -sidex / area;
        ymom_explicit_update[k] += -sidey / area;

    }
}


int main(int argc, char *argv[])
{
    int N, cnt1, cnt2;
    double g;
    double * sv, // stage_vertex_values
            * se, // stage_edge_values, 
            * sc, // stage_centroid_values, 
            * be, // bed_edge_values, 
            * bc, // bed_centroid_values, 
            * vc, // vertex_coordinates, 
            * xe, // xmom_explicit_update, 
            * ye, // ymom_explicit_update, 
            * n, // normals, 
            * a, // areas, 
            * el; // edgelengths;

    double *test_xe, *test_ye;
    int i;

    char file_name[]="merimbula.dat";

    freopen(file_name, "r", stdin);
    scanf("%d\n", &N);
    scanf("%lf\n",&g);
    printf("%d %lf\n", N, g);

    
    sv = (double *)malloc( N*3*sizeof(double));
    assert(sv != NULL);
    se = (double *)malloc( N*3*sizeof(double));
    assert(se != NULL);
    sc = (double *)malloc( N*sizeof(double));
    assert(sc != NULL);
    
    
    be = (double *)malloc( N*3*sizeof(double));
    assert(be != NULL);
    bc = (double *)malloc( N*sizeof(double));
    assert(bc != NULL);
    
    
    xe = (double *)malloc( N*sizeof(double));
    assert(xe != NULL);
    ye = (double *)malloc( N*sizeof(double));
    assert(ye != NULL);

    
    vc = (double *)malloc( N*6*sizeof(double));
    assert(vc != NULL);
    n  = (double *)malloc( N*6*sizeof(double));
    assert(n != NULL);
    a  = (double *)malloc( N*sizeof(double));
    assert(a != NULL);
    el = (double *)malloc( N*3*sizeof(double));
    assert(el != NULL);
    

    test_xe = (double *)malloc( N*sizeof(double));
    assert(test_xe != NULL);
    test_ye = (double *)malloc( N*sizeof(double));
    assert(test_ye != NULL);
    

    for(i=0; i < N; i++)
    {
        scanf("%lf %lf %lf\n", sv+i*3, sv+i*3 +1, sv+i*3 +2);
        scanf("%lf %lf %lf\n", se+i*3, se+i*3 +1, se+i*3 +2);
        scanf("%lf\n", sc+i);

        scanf("%lf %lf %lf\n",be+i*3, be+i*3 +1, be+i*3 +2);
        scanf("%lf\n", bc+i);

        scanf("%lf %lf\n", vc+i*6, vc+i*6 +1);
        scanf("%lf %lf\n", vc+i*6+2, vc+i*6 +3);
        scanf("%lf %lf\n", vc+i*6+4, vc+i*6 +5);

        scanf("%lf\n", xe+i);
        scanf("%lf\n", ye+i);

        scanf("%lf %lf %lf %lf %lf %lf\n", n+i*6, n+i*6+1, n+i*6+2,n+i*6+3, n+i*6+4, n+i*6+5);

        scanf("%lf\n", a+i);

        scanf("%lf %lf %lf\n", el+i*3, el+i*3 +1, el+i*3 +2);
    }
    gravity_wb( xe, ye, N, g, sv, se, sc, be, bc, vc, n, a, el);
    gravity_wb_orig( N, g, sv, se, sc, be, bc, vc, test_xe, test_ye, n, a, el);


    cnt1 = 0;
    cnt2 = 0;
    for (i=0; i < N; i++)
    {
        if (xe[i] != test_xe[i])
            cnt1 ++;
        if (ye[i] != test_ye[i])
            cnt2 ++;
    }
    printf("errors: %d %d\n", cnt1, cnt2);
}

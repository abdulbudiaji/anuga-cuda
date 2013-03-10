#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>


#define TOLERANCE 0.0000001

int errs_x=0, errs_y=0;


int check_tolerance(float a, float b)
{
    if (fabs(a) <= 0 && fabs(a-b) > TOLERANCE)
        return 1;
    else if(fabs(a) < 10 && fabs(a-b) > TOLERANCE*10)
        return 1;
    else if(fabs(a) < 100 && fabs(a-b) > TOLERANCE*100)
        return 1;
    else if(fabs(a) < 1000 && fabs(a-b) > TOLERANCE*1000)
        return 1;
    else if(fabs(a) < 10000 && fabs(a-b) > TOLERANCE*10000)
        return 1;
    else return 0;
}
        
void vecaddgpu( float *xmom_explicit_update, 
                float *ymom_explicit_update, 
                
                float * stage_vertex_values,
                float * stage_edge_values,
                float * stage_centroid_values,

                float * bed_edge_values,
                float * bed_centroid_values,

                float * vertex_coordinates,

                float * normals,
                float * areas,
                float * edgelengths,
                
                int n,
                int n3,
                int n6,
                float g )
{
    int k;

    #pragma acc kernels loop copyin(stage_vertex_values[0:n3], &
    #pragma acc & stage_edge_values[0:n3], &
    #pragma acc & stage_centroid_values[0:n], &
    #pragma acc & bed_edge_values[0:n3], &
    #pragma acc & bed_centroid_values[0:n], &
    #pragma acc & vertex_coordinates[0:n6], &
    #pragma acc & normals[0:n6], &
    #pragma acc & areas[0:n], &
    #pragma acc & edgelengths[0:n3]) &
    #pragma acc & copy(xmom_explicit_update[0:n], &
    #pragma acc & ymom_explicit_update[0:n])
    for( int k = 0; k < n; ++k ){
        int i, k3, k6;
        float w0, w1, w2,
              x0, y0, x1, y1, x2, y2,
              avg_h;

        float wx, wy, det, hh[3];//hh0, hh1,hh2;
        float sidex_0, sidex_1, sidex_2;
        float sidey_0, sidey_1, sidey_2;
        float area, n0, n1, fact;

        k3 = k*3;
        w0 = stage_vertex_values[k3];
        w1 = stage_vertex_values[k3 + 1];
        w2 = stage_vertex_values[k3 + 2];

        k6 = k*6;
        x0 = vertex_coordinates[k6 + 0];
        y0 = vertex_coordinates[k6 + 1];
        x1 = vertex_coordinates[k6 + 2];
        y1 = vertex_coordinates[k6 + 3];
        x2 = vertex_coordinates[k6 + 4];
        y2 = vertex_coordinates[k6 + 5];

        det = (y2 - y0)*(x1 - x0) - (y1 - y0)*(x2 - x0);

        wx = (y2 -y0)*(w1 - w0) - (y1 - y0)*(w2 -w0);
        wx /= det;

        wy = (x1 - x0)*(w2 - w0) - (x2 - x0)*(w1 -w0);
        wy /= det;

        avg_h = stage_centroid_values[k] - bed_centroid_values[k];

        xmom_explicit_update[k] += -g * wx * avg_h;
        ymom_explicit_update[k] += -g * wy * avg_h;


        hh[0] = stage_edge_values[k3] - bed_edge_values[k3];
        hh[1] = stage_edge_values[k3+1] - bed_edge_values[k3+1];
        hh[2] = stage_edge_values[k3+2] - bed_edge_values[k3+2];

        //sidex = 0.0;
        //sidey = 0.0;
        area = areas[k];
        
        sidex_0 = hh[0]*hh[0]*edgelengths[k3] * normals[k6];
        sidex_1 = hh[1]*hh[1]*edgelengths[k3+1] * normals[k6+2];
        sidex_2 = hh[2]*hh[2]*edgelengths[k3+2] * normals[k6+4];
        xmom_explicit_update[k] += 0.5*g*(sidex_0 + sidex_1 + sidex_2)/area;


        sidey_0 = hh[0]*hh[0]*edgelengths[k3] * normals[k6+1];
        sidey_1 = hh[1]*hh[1]*edgelengths[k3+1] * normals[k6+3];
        sidey_2 = hh[2]*hh[2]*edgelengths[k3+2] * normals[k6+5];
        ymom_explicit_update[k] += 0.5*g*(sidey_0 + sidey_1 + sidey_2)/area;
        
    }
}



void gravity_wb_orig(
        float * xmom_explicit_update, 
        float * ymom_explicit_update, 
        float * stage_vertex_values, 
        float * stage_edge_values, 
        float * stage_centroid_values, 
        float * bed_edge_values, 
        float * bed_centroid_values, 
        float * vertex_coordinates, 
        float * normals, 
        float * areas, 
        float * edgelengths,
        float * test_xe,
        float * test_ye,
        int N,
        float g
        )
{
    int i, k, k3, k6;

    float w0, w1, w2, 
           x0, y0, x1, y1, x2, y2,
           avg_h;

    float wx, wy, det,
           hh[3];
    float sidex, sidey, area, n0, n1, fact;

    for (k = 0; k < N; k++)
    {
        k3 = 3*k;
        w0 = stage_vertex_values[k3 + 0];
        w1 = stage_vertex_values[k3 + 1];
        w2 = stage_vertex_values[k3 + 2];


        k6 = 6*k;
        x0 = vertex_coordinates[k6 + 0];
        y0 = vertex_coordinates[k6 + 1];
        x1 = vertex_coordinates[k6 + 2];
        y1 = vertex_coordinates[k6 + 3];
        x2 = vertex_coordinates[k6 + 4];
        y2 = vertex_coordinates[k6 + 5];


        //_gradient(x0, y0, x1, y1, x2, y2, w0, w1, w2, &wx, &wy);

        det = (y2 - y0)*(x1 - x0) - (y1 - y0)*(x2 - x0);

        wx = (y2 -y0)*(w1 - w0) - (y1 - y0)*(w2 -w0);
        wx /= det;

        wy = (x1 - x0)*(w2 - w0) - (x2 - x0)*(w1 -w0);
        wy /= det;


        avg_h = stage_centroid_values[k] - bed_centroid_values[k];

        test_xe[k] += -g * wx * avg_h;
        test_ye[k] += -g * wy * avg_h;


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

            sidex += fact*n0;
            sidey += fact*n1;
        }

        area = areas[k];
        test_xe[k] += -sidex / area;
        test_ye[k] += -sidey / area;
        // For testing purpose
        if ( check_tolerance(xmom_explicit_update[k], test_xe[k]) )
        {
            errs_x += 1;
            if (errs_x <= 10)
                printf("   Errors on xe: %d %.9f %.9f\n",
                        k, xmom_explicit_update[k], test_xe[k]);
        }
        if ( check_tolerance(ymom_explicit_update[k], test_ye[k]) )
        {
            errs_y += 1;
            if (errs_y <= 10)
                printf("   Errors on ye: %d %f %f\n", 
                        k, ymom_explicit_update[k], test_ye[k]);
        }

    }
}



int main( int argc, char* argv[] ){
    int n;   /* vector length */
    float g = 9.8;
    int i;
    
    float * xe, //xmom_explicit_update, 
          * ye; //ymom_explicit_update;
    float * test_xe, //xmom_explicit_update, 
          * test_ye; //ymom_explicit_update;


    float * sv, //stage_vertex_values, 
          * se, //stage_edge_values, 
          * sc; //stage_centroid_values;

    float * be, //bed_edge_values, 
          * bc; //bed_centroid_values;

    float * vc, //vertex_coordinates,
          * normals, 
          * a, //areas, 
          * el; //edgelengths;


    char file_name[] = "merimbula.dat";
    freopen(file_name, "r", stdin);


    scanf("%d\n %f\n", &n, &g);
    printf("\n\nThe number of elements is %d\n", n);
    //if( argc > 1 ) n = atoi( argv[1] );
    //else  n = 100000;  /* default vector length */
    //if( n <= 0 ) n = 100000;
    
    
    xe =    (float*)malloc( n*sizeof(float));
    assert( xe != NULL);
    ye =    (float*)malloc( n*sizeof(float));
    assert( ye != NULL);
    test_xe =    (float*)malloc( n*sizeof(float));
    assert( test_xe != NULL);
    test_ye =    (float*)malloc( n*sizeof(float));
    assert( test_ye != NULL);


    sv =    (float*)malloc( n*sizeof(float) *3);
    assert( sv != NULL);
    se =    (float*)malloc( n*sizeof(float) *3);
    assert( se != NULL);
    sc =    (float*)malloc( n*sizeof(float) );
    assert( sc != NULL);
    
    be =    (float*)malloc( n*sizeof(float) *3);
    assert( be != NULL);
    bc =    (float*)malloc( n*sizeof(float) );
    assert( bc != NULL);

    vc =    (float*)malloc( n*sizeof(float) *6);
    assert( vc != NULL);
    normals=(float*)malloc( n*sizeof(float) *6);
    assert( normals != NULL);
    a =     (float*)malloc( n*sizeof(float) );
    assert( a != NULL);
    el =    (float*)malloc( n*sizeof(float) *3);
    assert( el != NULL);



    for( i = 0; i < n; ++i ){
        scanf("%f %f %f\n", sv+i*3, sv+i*3 +1, sv+i*3 +2);
        scanf("%f %f %f\n", se+i*3, se+i*3 +1, se+i*3 +2);
        scanf("%f\n", sc+i);

        scanf("%f %f %f\n",be+i*3, be+i*3 +1, be+i*3 +2);
        scanf("%f\n", bc+i);

        scanf("%f %f\n", vc+i*6, vc+i*6 +1);
        scanf("%f %f\n", vc+i*6+2, vc+i*6 +3);
        scanf("%f %f\n", vc+i*6+4, vc+i*6 +5);

        scanf("%f\n", xe+i);
        scanf("%f\n", ye+i);

        scanf("%f %f %f %f %f %f\n", normals+i*6, normals+i*6+1, normals+i*6+2,normals+i*6+3, normals+i*6+4, normals+i*6+5);

        scanf("%f\n", a+i);

        scanf("%f %f %f\n", el+i*3, el+i*3 +1, el+i*3 +2);
    }

    memcpy(test_xe, xe, n*sizeof(float));
    memcpy(test_ye, ye, n*sizeof(float));
    errs_x = 0;
    errs_y = 0;
    for( i = 0; i < n; ++i ){
        if (test_xe[i] != xe[i])
            errs_x +=1;
        if (test_ye[i] != ye[i])
            errs_y +=1;
    }
    printf("Verifying initial input: %d %d \n", errs_x, errs_y);

    /* compute on the GPU */
    printf(" --> Entering GPU .... \n");
    vecaddgpu( xe, ye, sv, se, sc, be, bc, vc, normals, a, el, n,n*3,n*6, g);
    
    /* compare results */
    printf(" --> Entering CPU .... \n");
    errs_x = 0;
    errs_y = 0;
    gravity_wb_orig( xe, ye, sv, se, sc, be, bc, vc, normals, a, el, 
            test_xe, test_ye, n, g);

    
    
    printf( "%d  %derrors found\n", errs_x, errs_y );
    return errs_x + errs_y;
}

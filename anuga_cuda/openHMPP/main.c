#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "fun.h"


#define DATA_TYPE double 

#define CHECK_RESULT(a, b) ( fabs(a-b) > 1e-11)



int main( int argc, char* argv[] ){
    int n;   /* vector length */
    int i, temp;

    DATA_TYPE * sv, //stage_vertex_values, 
          * se, //stage_edge_values, 
          * sc, //stage_centroid_values,
          * s_x, // stage_x_gradient
          * s_y; // stage_y_gradient
    DATA_TYPE beta;

    DATA_TYPE * centroid_coordinates,
            * vertex_coordinates;
    long * number_of_boundaries,
            * surrogate_neighbours,
            * neighbours;

    DATA_TYPE *test_c,
            *test_v,
            *test_e,
            *test_x, 
            *test_y;

    int cnt_c, cnt_v, cnt_e, cnt_x, cnt_y;

    char file_name[] = "merimbula_all.dat";
    freopen(file_name, "r", stdin);

    scanf("%d\n %lf\n", &n, &beta);
    
    centroid_coordinates =    (DATA_TYPE*)malloc( n*sizeof(DATA_TYPE)*2);
    //assert( xe != NULL);
    vertex_coordinates =    (DATA_TYPE*)malloc( n*sizeof(DATA_TYPE) * 6);
    //assert( ye != NULL);
    number_of_boundaries =    (long*)malloc( n*sizeof(long));
    //assert( test_xe != NULL);
    surrogate_neighbours =    (long*)malloc( n*sizeof(long)*3);
    //assert( test_ye != NULL);
    neighbours =    (long*)malloc( n*sizeof(long)*3);
    //assert( test_ye != NULL);


    sv =    (DATA_TYPE*)malloc( n*sizeof(DATA_TYPE) *3);
    //assert( sv != NULL);
    se =    (DATA_TYPE*)malloc( n*sizeof(DATA_TYPE) *3);
    //assert( se != NULL);
    sc =    (DATA_TYPE*)malloc( n*sizeof(DATA_TYPE) );
    //assert( sc != NULL);
    s_x =    (DATA_TYPE*)malloc( n*sizeof(DATA_TYPE) );
    //assert( be != NULL);
    s_y =    (DATA_TYPE*)malloc( n*sizeof(DATA_TYPE) );
    //assert( bc != NULL);
    

    test_v =    (DATA_TYPE*)malloc( n*sizeof(DATA_TYPE) *3);
    //assert( sv != NULL);
    test_e =    (DATA_TYPE*)malloc( n*sizeof(DATA_TYPE) *3);
    //assert( se != NULL);
    test_c =    (DATA_TYPE*)malloc( n*sizeof(DATA_TYPE) );
    //assert( sc != NULL);
    test_x =     (DATA_TYPE*)malloc( n*sizeof(DATA_TYPE) );
    //assert( a != NULL);
    test_y =    (DATA_TYPE*)malloc( n*sizeof(DATA_TYPE) );
    //assert( el != NULL);



    for( i = 0; i < n; ++i ){
        scanf("%lf %lf\n", 
            centroid_coordinates+i*2, centroid_coordinates+i*2 +1);

        scanf("%lf %lf %lf %lf %lf %lf\n", 
            vertex_coordinates+i*6, 
            vertex_coordinates+i*6 +1, 
            vertex_coordinates+i*6 +2,
            vertex_coordinates+i*6 +3, 
            vertex_coordinates+i*6 +4, 
            vertex_coordinates+i*6 +5);

        scanf("%d\n", &temp);
        number_of_boundaries[i] = temp;

        scanf("%ld %ld %ld\n",
            surrogate_neighbours+i*3, 
            surrogate_neighbours+i*3 +1, 
            surrogate_neighbours+i*3 +2);

        scanf("%ld %ld %ld\n",
            neighbours+i*3, 
            neighbours+i*3 +1, 
            neighbours+i*3 +2);

        scanf("%lf\n", sc+i);
        scanf("%lf %lf %lf\n", sv+i*3, sv+i*3 +1, sv+i*3 +2);
        scanf("%lf %lf %lf\n", se+i*3, se+i*3 +1, se+i*3 +2);

    }


    for( i = 0; i < n; ++i ){
        test_c[i] = sc[i];

        test_v[i*3] = sv[i*3];
        test_v[i*3+1] = sv[i*3+1];
        test_v[i*3+2] = sv[i*3+2];

        test_e[i*3] = se[i*3];
        test_e[i*3+1] = se[i*3+1];
        test_e[i*3+2] = se[i*3+2];
    }


    /* compute on the GPU */
    printf(" --> Entering GPU .... \n");
    
    extrapolate_second_order_and_limit_by_vertex(
            n,
            n*2,
            n*3,
            n*6,
            beta,

            centroid_coordinates,
            vertex_coordinates,
            number_of_boundaries,
            surrogate_neighbours,
            neighbours,

            sc,
            sv,
            se,
            s_x,
            s_y
            );

    /* compare results */
    printf(" --> Entering CPU .... \n");
    extrapolate_second_order_and_limit_by_vertex_normal(
            n,
            n*2,
            n*3,
            n*6,
            beta,

            centroid_coordinates,
            vertex_coordinates,
            number_of_boundaries,
            surrogate_neighbours,
            neighbours,

            test_c,
            test_v,
            test_e,
            test_x,
            test_y
            );
    
    cnt_c = cnt_v = cnt_e = cnt_x = cnt_y = 0;
    for(i=0; i<n; i++)
    {
        if(sc[i] != test_c[i])
        {
            cnt_c += 1;
            if (cnt_c < 10)
                printf("sc: %d  %lf  %lf\n", i, sc[i], test_c[i]);
        }

        if( CHECK_RESULT( sv[i*3+0], test_v[i*3+0]) || 
            CHECK_RESULT( sv[i*3+1], test_v[i*3+1]) || 
            CHECK_RESULT( sv[i*3+2], test_v[i*3+2]) )
        {
            cnt_v += 1;
            if (cnt_v < 10)
                printf("sv: %d [%lf %lf %lf]  [%lf %lf %lf]\n", 
                    i, sv[i*3+0],sv[i*3+1],sv[i*3+2], test_v[i*3+0], 
                    test_v[i*3+1], test_v[i*3+2]);
        }

        if(se[i*3+0] != test_e[i*3+0] || se[i*3+1] !=test_e[i*3+1] || se[i*3+2]!=test_e[i*3+2])
            cnt_e += 1;

        if( CHECK_RESULT( s_x[i], test_x[i] ))
        {
            cnt_x += 1;
            if (cnt_x < 10)
                printf("sx: %d  %lf  %lf\n", i, s_x[i], test_x[i]);
        }

        if( CHECK_RESULT( s_y[i], test_y[i] ))
            cnt_y += 1;
    }
    
    if ( cnt_c || cnt_v || cnt_e || cnt_x || cnt_y )
        printf( "%d  %d  %d  %d  %derrors found\n", 
            cnt_c, cnt_v, cnt_e, cnt_x, cnt_y);
}

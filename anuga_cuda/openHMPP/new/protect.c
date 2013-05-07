//#define USING_MAIN

#ifndef USING_MAIN
#include "hmpp_fun.h"
#else
#include <math.h>
#endif

//#ifdef USING_MAIN
#ifdef USING_LOCAL_DIRECTIVES
#pragma hmpp protectSW codelet, target=CUDA args[*].transfer=atcall
#endif
void protect_sw(
        int N,
        int N3,
        double minimum_allowed_height,
        double maximum_allowed_speed,
        double epsilon,

        double wc[N],
        double zc[N],
        double xmomc[N],
        double ymomc[N]) 
{
    int k;
    double hc;
    double u, v, reduced_speed;


    // Protect against initesimal and negative heights
    if (maximum_allowed_speed < epsilon) {
        #pragma hmppcg gridify(k), &
        #pragma hmppcg & private( hc, u, v, reduced_speed), &
        #pragma hmppcg & global( minimum_allowed_height, maximum_allowed_speed, &
        #pragma hmppcg & epsilon, wc, zc, xmomc, ymomc)
        for (k = 0; k < N; k++) {
            hc = wc[k] - zc[k];
            if (hc < minimum_allowed_height) {
                // Set momentum to zero and ensure h is non negative
                xmomc[k] = 0.0;
                ymomc[k] = 0.0;
                if (hc <= 0.0) wc[k] = zc[k];
            }
        }

    } else {
        // Protect against initesimal and negative heights
        #pragma hmppcg gridify(k), &
        #pragma hmppcg & private( hc, u, v, reduced_speed), &
        #pragma hmppcg & global( minimum_allowed_height, maximum_allowed_speed, &
        #pragma hmppcg & epsilon, wc, zc, xmomc, ymomc)
        for (k = 0; k < N; k++) {
            hc = wc[k] - zc[k];
            if (hc < minimum_allowed_height) {
                //New code: Adjust momentum to guarantee speeds are physical
                //          ensure h is non negative
                if (hc <= 0.0) {
                    wc[k] = zc[k];
                    xmomc[k] = 0.0;
                    ymomc[k] = 0.0;
                } else {
                    //Reduce excessive speeds derived from division by small hc
                    //FIXME (Ole): This may be unnecessary with new slope limiters
                    //in effect.

                    u = xmomc[k] / hc;
                    if (fabs(u) > maximum_allowed_speed) {
                        reduced_speed = maximum_allowed_speed * u / fabs(u);
                        xmomc[k] = reduced_speed * hc;
                    }

                    v = ymomc[k] / hc;
                    if (fabs(v) > maximum_allowed_speed) {
                        reduced_speed = maximum_allowed_speed * v / fabs(v);
                        ymomc[k] = reduced_speed * hc;
                    }
                }
            }
        }
    }
}



#ifdef USING_MAIN
int main()
#else
void test_protect_sw()
#endif
{
    int N=4;
    double minimum_allowed_height=0.001, 
            maximum_allowed_speed=0, 
            epsilon=1e-12;

    double  wc[]={1,2,3,4}, 
            zc[]={1,2,3,4}, 
            xmomc[]={1,2,3,4}, 
            ymomc[]={1,2,3,4};

    #pragma hmpp protectSW callsite
    protect_sw(N, N*3, minimum_allowed_height,
            maximum_allowed_speed, epsilon, wc, zc, xmomc, ymomc);
}

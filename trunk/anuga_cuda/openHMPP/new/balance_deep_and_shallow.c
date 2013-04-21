#include "hmpp_fun.h"

void balance_deep_and_shallow(
        int N,
        int N3,
        double H0,
        double alpha_balance,
        int tight_slope_limiters,
        int use_centroid_velocities,
        double wc[N],   // stage_centroid_values
        double zc[N],   // elevation_centroid_values
        double wv[N3],  // stage_vertex_values
        double zv[N3],  // elevation_vertex_values
        //double* hvbar,// Retire this
        double xmomc[N],  // xmom_centroid_values
        double ymomc[N],  // ymom_centroid_values
        double xmomv[N3],  // xmom_vertex_values
        double ymomv[N3]   // ymom_vertex_values
        ) 
{
    int k, i;

    double dz, hmin, alpha, h_diff, hc_k;
    double epsilon = 1.0e-6; // FIXME: Temporary measure
#ifdef USING_ORIGINAL_FUNCTION
    double hv[3]; // Depths at vertices
#else 
    double hv0, hv1, hv2;
#endif
    double uc, vc; // Centroid speeds

    // Compute linear combination between w-limited stages and
    // h-limited stages close to the bed elevation.

    for (k = 0; k < N; k++) {
        // Compute maximal variation in bed elevation
        // This quantitiy is
        //     dz = max_i abs(z_i - z_c)
        // and it is independent of dimension
        // In the 1d case zc = (z0+z1)/2
        // In the 2d case zc = (z0+z1+z2)/3
#ifndef REARRANGED_DOMAIN
        int k3 = 3 * k;
#endif

        hc_k = wc[k] - zc[k]; // Centroid value at triangle k

        dz = 0.0;
        if (tight_slope_limiters == 0) {
            // FIXME: Try with this one precomputed
            //for (i = 0; i < 3; i++) {
            //    dz = fmax(dz, fabs(zv[k3 + i] - zc[k]));
            //}
#ifndef REARRANGED_DOMAIN
                dz = fmax(dz, fabs(zv[k3] - zc[k]));
                dz = fmax(dz, fabs(zv[k3 + 1] - zc[k]));
                dz = fmax(dz, fabs(zv[k3 + 2] - zc[k]));
#else   
                dz = fmax(dz, fabs(zv[k] -zc[k]));
                dz = fmax(dz, fabs(zv[k + N] -zc[k]));
                dz = fmax(dz, fabs(zv[k + 2*N] -zc[k]));
#endif
        }

        // Calculate depth at vertices (possibly negative here!)
#ifndef REARRANGED_DOMAIN
#ifdef USING_ORIGINAL_FUNCTION
        hv[0] = wv[k3] - zv[k3];
        hv[1] = wv[k3 + 1] - zv[k3 + 1];
        hv[2] = wv[k3 + 2] - zv[k3 + 2];
#else
        hv0 = wv[k3] - zv[k3];
        hv1 = wv[k3 + 1] - zv[k3 + 1];
        hv2 = wv[k3 + 2] - zv[k3 + 2];
#endif
#else   
#ifdef USING_ORIGINAL_FUNCTION
        hv[0] = wv[k] - zv[k];
        hv[1] = wv[k + N] - zv[k + N];
        hv[2] = wv[k + 2*N] - zv[k + 2*N];
#else
        hv0 = wv[k] - zv[k];
        hv1 = wv[k + N] - zv[k + N];
        hv2 = wv[k + 2*N] - zv[k + 2*N];
#endif
#endif

        // Calculate minimal depth across all three vertices
#ifdef USING_ORIGINAL_FUNCTION
        hmin = fmin(hv[0], fmin(hv[1], hv[2]));
#else
        hmin = fmin(hv0, fmin(hv1, hv2));
#endif

        //if (hmin < 0.0 ) {
        //  printf("hmin = %f\n",hmin);
        //}


        // Create alpha in [0,1], where alpha==0 means using the h-limited
        // stage and alpha==1 means using the w-limited stage as
        // computed by the gradient limiter (both 1st or 2nd order)
        if (tight_slope_limiters == 0) {
            // If hmin > dz/alpha_balance then alpha = 1 and the bed will have no
            // effect
            // If hmin < 0 then alpha = 0 reverting to constant height above bed.
            // The parameter alpha_balance==2 by default


            if (dz > 0.0) {
                alpha = fmax(fmin(alpha_balance * hmin / dz, 1.0), 0.0);
            } else {
                alpha = 1.0; // Flat bed
            }
            //printf("Using old style limiter\n");

        } else {

            // Tight Slope Limiter (2007)

            // Make alpha as large as possible but still ensure that
            // final depth is positive and that velocities at vertices
            // are controlled

            if (hmin < H0) {
                alpha = 1.0;
#ifdef USING_ORIGINAL_FUNCTION
                for (i = 0; i < 3; i++) {
                    h_diff = hc_k - hv[i];
                    if (h_diff <= 0) {
                        // Deep water triangle is further away from bed than
                        // shallow water (hbar < h). Any alpha will do

                    } else {
                        // Denominator is positive which means that we need some of the
                        // h-limited stage.

                        alpha = fmin(alpha, (hc_k - H0) / h_diff);
                    }
                }
#else
                h_diff = hc_k - hv0;
                if (h_diff > 0) 
                    alpha = fmin(alpha, (hc_k - H0) / h_diff);
                h_diff = hc_k - hv1;
                if (h_diff > 0) 
                    alpha = fmin(alpha, (hc_k - H0) / h_diff);
                h_diff = hc_k - hv2;
                if (h_diff > 0) 
                    alpha = fmin(alpha, (hc_k - H0) / h_diff);
#endif

                // Ensure alpha in [0,1]
                if (alpha > 1.0) alpha = 1.0;
                if (alpha < 0.0) alpha = 0.0;

            } else {
                // Use w-limited stage exclusively in deeper water.
                alpha = 1.0;
            }
        }


        //  Let
        //
        //    wvi be the w-limited stage (wvi = zvi + hvi)
        //    wvi- be the h-limited state (wvi- = zvi + hvi-)
        //
        //
        //  where i=0,1,2 denotes the vertex ids
        //
        //  Weighted balance between w-limited and h-limited stage is
        //
        //    wvi := (1-alpha)*(zvi+hvi-) + alpha*(zvi+hvi)
        //
        //  It follows that the updated wvi is
        //    wvi := zvi + (1-alpha)*hvi- + alpha*hvi
        //
        //   Momentum is balanced between constant and limited


        if (alpha < 1) {
#ifdef USING_ORIGINAL_FUNCTION
            for (i = 0; i < 3; i++) {
#ifndef REARRANGED_DOMAIN
                wv[k3 + i] = zv[k3 + i] + (1 - alpha) * hc_k + alpha * hv[i];
#else   
                wv[k + i*N] = zv[k + i*N] + (1 - alpha) * hc_k + alpha * hv[i];
#endif

                // Update momentum at vertices
                if (use_centroid_velocities == 1) {
                    // This is a simple, efficient and robust option
                    // It uses first order approximation of velocities, but retains
                    // the order used by stage.

                    // Speeds at centroids
                    if (hc_k > epsilon) {
                        uc = xmomc[k] / hc_k;
                        vc = ymomc[k] / hc_k;
                    } else {
                        uc = 0.0;
                        vc = 0.0;
                    }

                    // controlled speed
                    // Recompute (balanced) vertex depth
#ifndef REARRANGED_DOMAIN
                    hv[i] = wv[k3 + i] - zv[k3 + i]; 
                    xmomv[k3 + i] = uc * hv[i];
                    ymomv[k3 + i] = vc * hv[i];
#else   
                    hv[i] = wv[k + i*N] - zv[k + i*N]; 
                    xmomv[k + i*N] = uc * hv[i];
                    ymomv[k + i*N] = vc * hv[i];
#endif

                } else {
                    // Update momentum as a linear combination of
                    // xmomc and ymomc (shallow) and momentum
                    // from extrapolator xmomv and ymomv (deep).
                    // This assumes that values from xmomv and ymomv have
                    // been established e.g. by the gradient limiter.

                    // FIXME (Ole): I think this should be used with vertex momenta
                    // computed above using centroid_velocities instead of xmomc
                    // and ymomc as they'll be more representative first order
                    // values.

#ifndef REARRANGED_DOMAIN
                    xmomv[k3 + i] = (1 - alpha) * xmomc[k] + alpha *xmomv[k3 + i];
                    ymomv[k3 + i] = (1 - alpha) * ymomc[k] + alpha *ymomv[k3 + i];
#else   
                    xmomv[k + i*N] = (1-alpha) * xmomc[k] + alpha *xmomv[k + i*N];
                    ymomv[k + i*N] = (1-alpha) * ymomc[k] + alpha *ymomv[k + i*N];
#endif

                }
            }


#else // USING_ORIGINAL_FUNCTION
            //for (i = 0; i < 3; i++) {
            // i = 0
#ifndef REARRANGED_DOMAIN
//            wv[k3] = zv[k3] + (1 - alpha) * hc_k + alpha * hv0;
#else   
            wv[k] = zv[k] + (1 - alpha) * hc_k + alpha * hv0;
#endif
            if (use_centroid_velocities == 1) {
                if (hc_k > epsilon) {
                    uc = xmomc[k] / hc_k;
                    vc = ymomc[k] / hc_k;
                } else {
                    uc = 0.0;
                    vc = 0.0;
                }
#ifndef REARRANGED_DOMAIN
                hv0 = wv[k3] - zv[k3]; 
                xmomv[k3] = uc * hv0;
                ymomv[k3] = vc * hv0;
#else   
                hv0 = wv[k] - zv[k]; 
                xmomv[k] = uc * hv0;
                ymomv[k] = vc * hv0;
#endif
            } else {
#ifndef REARRANGED_DOMAIN
//                xmomv[k3] = (1 - alpha) * xmomc[k] + alpha *xmomv[k3];
//                ymomv[k3] = (1 - alpha) * ymomc[k] + alpha *ymomv[k3];
#else   
                xmomv[k] = (1-alpha) * xmomc[k] + alpha *xmomv[k];
                ymomv[k] = (1-alpha) * ymomc[k] + alpha *ymomv[k];
#endif
            }
/*
            // i = 1
#ifndef REARRANGED_DOMAIN
            wv[k3 + 1] = zv[k3 + 1] + (1 - alpha) * hc_k + alpha * hv1;
#else   
            wv[k + N] = zv[k + N] + (1 - alpha) * hc_k + alpha * hv1;
#endif
            if (use_centroid_velocities == 1) {
                if (hc_k > epsilon) {
                    uc = xmomc[k] / hc_k;
                    vc = ymomc[k] / hc_k;
                } else {
                    uc = 0.0;
                    vc = 0.0;
                }
#ifndef REARRANGED_DOMAIN
                hv1 = wv[k3 + 1] - zv[k3 + 1]; 
                xmomv[k3 + 1] = uc * hv1;
                ymomv[k3 + 1] = vc * hv1;
#else   
                hv1 = wv[k + N] - zv[k + N]; 
                xmomv[k + N] = uc * hv1;
                ymomv[k + N] = vc * hv1;
#endif

            } else {
#ifndef REARRANGED_DOMAIN
                xmomv[k3 + 1] = (1 - alpha) * xmomc[k] + alpha *xmomv[k3 + 1];
                ymomv[k3 + 1] = (1 - alpha) * ymomc[k] + alpha *ymomv[k3 + 1];
#else   
                xmomv[k + N] = (1-alpha) * xmomc[k] + alpha *xmomv[k + N];
                ymomv[k + N] = (1-alpha) * ymomc[k] + alpha *ymomv[k + N];
#endif
            }
            // i = 2
#ifndef REARRANGED_DOMAIN
            wv[k3 + 2] = zv[k3 + 2] + (1 - alpha) * hc_k + alpha * hv2;
#else   
            wv[k + 2*N] = zv[k + 2*N] + (1 - alpha) * hc_k + alpha * hv2;
#endif

            // Update momentum at vertices
            if (use_centroid_velocities == 1) {
                // This is a simple, efficient and robust option
                // It uses first order approximation of velocities, but retains
                // the order used by stage.

                // Speeds at centroids
                if (hc_k > epsilon) {
                    uc = xmomc[k] / hc_k;
                    vc = ymomc[k] / hc_k;
                } else {
                    uc = 0.0;
                    vc = 0.0;
                }

                // controlled speed
                // Recompute (balanced) vertex depth
#ifndef REARRANGED_DOMAIN
                hv2 = wv[k3 + 2] - zv[k3 + 2]; 
                xmomv[k3 + 2] = uc * hv2;
                ymomv[k3 + 2] = vc * hv2;
#else   
                hv2 = wv[k + 2*N] - zv[k + 2*N]; 
                xmomv[k + 2*N] = uc * hv2;
                ymomv[k + 2*N] = vc * hv2;
#endif

            } else {
                // Update momentum as a linear combination of
                // xmomc and ymomc (shallow) and momentum
                // from extrapolator xmomv and ymomv (deep).
                // This assumes that values from xmomv and ymomv have
                // been established e.g. by the gradient limiter.

                // FIXME (Ole): I think this should be used with vertex momenta
                // computed above using centroid_velocities instead of xmomc
                // and ymomc as they'll be more representative first order
                // values.

#ifndef REARRANGED_DOMAIN
                xmomv[k3 + 2] = (1 - alpha) * xmomc[k] + alpha *xmomv[k3 + 2];
                ymomv[k3 + 2] = (1 - alpha) * ymomc[k] + alpha *ymomv[k3 + 2];
#else   
                xmomv[k + 2*N] = (1-alpha) * xmomc[k] + alpha *xmomv[k + 2*N];
                ymomv[k + 2*N] = (1-alpha) * ymomc[k] + alpha *ymomv[k + 2*N];
#endif
            }
*/
#endif // USING_ORIGINAL_FUNCTION
            //} //for (i = 0; i < 3; i++) 
        }// if (alpha < 1) {
    }// for (k = 0; k < N; k++) {
}


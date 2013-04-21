#include "hmpp_fun.h"


void saxpy_centroid_values(
        int N,
        double a,
        double b,
        double centroid_values[N],
        double centroid_backup_values[N]
        )
{
    int k;
    for(k=0; k<N; k++)
        centroid_values[k]= a*centroid_values[k] + b*centroid_backup_values[k];
}

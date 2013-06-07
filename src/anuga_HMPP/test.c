#include <stdio.h>

#pragma hmpp myCall codelet, target=CUDA, transfer=atcall, &
#pragma hmpp & args[A, B].mirror, args[A, B].transfer=manual
void myFunc(int n, int *A, int *B)
{
    int i, i3;

    int hh[2];

    #pragma hmppcg gridify(i), private(i3, hh)
    for (i=0; i<n ; ++i)
    {
        i3 = i*1;
        hh[0] = 3*A[i];
        hh[1] = 4*A[i]; 
        if (i < 3)
            B[i3] =  hh[0];
        else
            B[i] = A[i] + hh[1];
    }
}



#pragma hmpp myCall2 codelet, target=CUDA, transfer=atcall, &
#pragma hmpp & args[C, D].mirror, args[C, D].transfer=manual
void myFunc2(int n, int C[n], int D[n])
{
    int i, i3;

    int hh[2];

    #pragma hmppcg gridify(i), private(i3, hh)
    for (i=0; i<n ; ++i)
    {
        i3 = i*1;
        hh[0] = 3*C[i];
        hh[1] = 4*C[i]; 
        if (i < 3)
            D[i3] =  hh[0];
        else
            D[i] = C[i] + hh[1];
    }
}

void temp(int n, int *E, int *F)
{
    printf(" In GPU\n");
#pragma hmpp myCall2 callsite
    myFunc2(n, E, F);
}


void main(void) {
    int n = 10, i;
    int X[10], Y[10], Z[10];
    for(i = 0; i< n; i++)
    {
        X[i] = i;
        Y[i] = 0;
        Z[i] = 0;
    }

#pragma hmpp myCall allocate, data[X, Y, Z], size={n}

#pragma hmpp advancedload data[X]

    printf(" In GPU\n");
#pragma hmpp myCall callsite
    myFunc(n, X, Y);
#pragma hmpp delegatedstore data[Y]
    for(i = 0; i< n; i++)
        printf("%d ", Y[i]);
    printf("\n");
    
    for(i = 0; i< n; i++)
        printf("%d ", Z[i]);
    printf("\n");

    temp(n, X, Z);
    for(i = 0; i< n; i++)
        printf("%d ", Z[i]);
    printf("\n");
#pragma hmpp delegatedstore data[Z]
    for(i = 0; i< n; i++)
        printf("%d ", Z[i]);
    printf("\n");



#pragma hmpp myCall release
#pragma hmpp myCall2 release


    for(i = 0; i< n; i++)
        if (Y[i] != Z[i])
            printf(" %d %d %d A=%d\n", i, Y[i], Z[i], X[i]);
    printf(" In CPU\n");
    myFunc(n, X, Z);

    for(i = 0; i< n; i++)
        if (Y[i] != Z[i])
            printf(" %d %d %d A=%d\n", i, Y[i], Z[i], X[i]);
}

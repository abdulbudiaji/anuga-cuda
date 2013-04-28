

#pragma hmpp myCall codelet, args[*].transfer=atcall, &
#pragma hmpp & target=CUDA
void myFunc(int n, int A[n], int B[n])
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



void main(void) {
    int n = 10, i;
    int X[10], Y[10], Z[10];
    for(i = 0; i< n; i++)
    {
        X[i] = i;
        Y[i] = 0;
        Z[i] = 0;
    }

    printf(" In GPU\n");
#pragma hmpp myCall callsite
    myFunc(n, X, Y);

    for(i = 0; i< n; i++)
        if (Y[i] != Z[i])
            printf(" %d %d %d A=%d\n", i, Y[i], Z[i], X[i]);
    printf(" In CPU\n");
    myFunc(n, X, Z);

    for(i = 0; i< n; i++)
        if (Y[i] != Z[i])
            printf(" %d %d %d A=%d\n", i, Y[i], Z[i], X[i]);
}

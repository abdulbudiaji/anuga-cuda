#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

int main()
{
	int i =0 ;
	assert(i == 0);
    if (i != 1)
        exit(EXIT_FAILURE);
    else
        printf("!");
}

#ifndef _CUDAFUN_H_
#define _CUDAFUN_H_
#include <stdlib.h>

void getDeviceInfo( int, int, const char* );
int  selectDevice( int, int, const char* );

void setKernelDims( const int, const int );
void printKernelDims(); 

#ifdef TEXCACHE
void bindTexRefToPtr();
void bindMemoryToTexCache( double*, int );
#endif

void* allocDeviceMemory( size_t );
void* allocHostMemory( size_t );
void copyDeviceToHost( void*, void*, size_t );
void copyHostToDevice( void*, void*, size_t );
void freeDeviceMemory( void* );
void freeHostMemory( void* );

void dummy();

void _set_to_default( double*, double*, double*, size_t, double);
double _compute_fluxes_central( int, double, double, double, double,
                      long*, long*, double*, double*, double*, double*,
                      long*, double*, double*, double*, double*, double*,
                      double*, double*, double*, double *, double *, double*,
                      int );

#endif


// this file contains some useful inline functions, mostly for managing memory
#include<stdio.h>
#include<iostream>
#include<stdlib.h>
typedef int Vec[3];

static const int X = 0;
static const int  Y = 1;
static const int  Z = 2;
static const int Inf = 1024*1024*1024;

using namespace std;

static inline int*** alloc_arr(size_t x, size_t y) {
    int ***arr = (int***) malloc( x * sizeof(int**) );
    if (arr == NULL) return NULL;
    for (size_t i=0 ; i<x ; i++) {
        arr[i] = (int**) malloc( y * sizeof(int*) );
        for(size_t j=0; j<y; j++)
            arr[i][j]= (int*) malloc(3*sizeof(int));
    }
    return arr;
}

static inline void free_arr(int ***arr, size_t sizeX, size_t sizeY) {
    if (arr != NULL) {
        for(size_t i=0; i<sizeX; i++)
            for(size_t j=0; j<sizeY; j++)
                free(arr[i][j]);
        for (size_t i=0 ; i<sizeX ; i++)
            free(arr[i]);
        free(arr);
    }
}

static inline void fill_arr(int value, int ***arr, size_t pos, size_t imin, size_t imax, size_t jmin, size_t jmax) {
    for (size_t i=imin ; i<=imax ; i++) {
        for (size_t j=jmin ; j<=jmax ; j++) {
            arr[i][j][pos] = value;
            arr[i][j][pos] = value;
            arr[i][j][pos] = value;

            //cout<<"i j value: "<<i<<" "<<j<<" "<<value<<endl;
        }
    }
}

static inline int maximum2(int arg0, int arg1, int arg2, int *idx) {

    if (arg0 >= arg1 && arg0 >= arg2) {
        *idx = X;
        return arg0;
    }
    if (arg1 >= arg0 && arg1 >= arg2) {
        *idx = Y;
        return arg1;
    }

    *idx = Z;
    return arg2;
}


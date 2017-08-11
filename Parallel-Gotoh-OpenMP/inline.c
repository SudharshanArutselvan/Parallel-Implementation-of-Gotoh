// this file contains some useful inline functions, mostly for managing memory

typedef int Vec[3];

static const int X = 0;
static const int Y = 1;
static const int Z = 2;
static const int Inf = 1024*1024*1024;

static inline Vec** alloc_arr(int x, int y) {
    Vec **arr = (Vec**) malloc( x * sizeof(Vec*) );
    if (arr == NULL) return NULL;
	// #pragma omp parallel for 
    for (int i=0 ; i<x ; i++) {
        arr[i] = (Vec*) malloc( y * sizeof(Vec) );
    }
    return arr;
}

static inline void free_arr(Vec **arr, int size) {
    if (arr != NULL) {
        for (int i=0 ; i<size ; i++) free(arr[i]);
        free(arr);
    }
}

static inline void fill_arr(int value, Vec **arr, int pos, int imin, int imax, int jmin, int jmax) {
    // #pragma omp parallel for collapse(2)
	for (int i=imin ; i<=imax ; i++) {
        for (int j=jmin ; j<=jmax ; j++) {
            arr[i][j][pos] = value;
            arr[i][j][pos] = value;
            arr[i][j][pos] = value;
        }
    }
}

static inline int maximum(int arg0, int arg1, int arg2, int *idx) {

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


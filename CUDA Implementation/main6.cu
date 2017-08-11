#include <stdio.h>
#include "gotoh.h"
#include "inline.cc"
#include<iostream>
#include "timer.h"
using namespace std;

__device__ int maximum(int arg0, int arg1, int arg2, int *idx)
{

    if (arg0 >= arg1 && arg0 >= arg2) {
        *idx = 0;
        return arg0;
    }
    if (arg1 >= arg0 && arg1 >= arg2) {
       *idx = 1;
        return arg1;
    }

    *idx = 2;
    return arg2;
}

__global__ void block_max(int *a, int no,int N, int *b, int *gapX, int *gapY, const int n)

{
    __shared__ int mat[3][3][3];
    int index_i=0,index_j=0;

    if(no<N/n)
    {
        index_i=n*blockIdx.x;
        index_j=n*no-index_i;
    }
    else
    {
        index_i=(blockIdx.x+no-N/n+1)*n;
        index_j=n*no-index_i;
    }

    int index = index_j*(N+1)*3+index_i*3;

    //column by column copying; per column, 3 rows (tid.x), each tid.x has 3 tid.y 
    for(int i=0;i<n+1;i++)
    {    
        mat[threadIdx.x][i][threadIdx.y]=a[index+i*(N+1)*3+3*threadIdx.x+threadIdx.y];
    }

    for(int i=0;i<2*n-1;i++)
    {
        if(i<n)
        {
            if(threadIdx.x<=i)
            {
                if(threadIdx.y==0)
                    mat[threadIdx.x+1][i-threadIdx.x+1][0]+=maximum(mat[threadIdx.x][i-threadIdx.x+1][0], mat[threadIdx.x][i-threadIdx.x+1][1] + gapY[index_j+1], mat[threadIdx.x][i-threadIdx.x+1][2] + gapY[index_j+1], &b[index+(i-threadIdx.x+1)*(N+1)*3+(threadIdx.x+1)*3]);

                if(threadIdx.y==1) 
                    mat[threadIdx.x+1][i-threadIdx.x+1][1]+=maximum(mat[threadIdx.x+1][i-threadIdx.x][0] + gapX[index_i+1], mat[threadIdx.x+1][i-threadIdx.x][1], mat[threadIdx.x+1][i-threadIdx.x][2] + gapX[index_i+1], &b[index+(i-threadIdx.x+1)*(N+1)*3+(threadIdx.x+1)*3+1]);

                if(threadIdx.x==2)
                    mat[threadIdx.x+1][i-threadIdx.x+1][2]+=maximum(mat[threadIdx.x][i-threadIdx.x][0], mat[threadIdx.x][i-threadIdx.x][1], mat[threadIdx.x][i-threadIdx.x][2], &b[index+(i-threadIdx.x+1)*(N+1)*3+(threadIdx.x+1)*3+2]);
            }
        }
        else
        {
            if(threadIdx.x<i-n+1)
            {
                if(threadIdx.y==0)
                    mat[n-threadIdx.x][i+2-n+threadIdx.x][0]+=maximum(mat[n-threadIdx.x-1][i+2-n+threadIdx.x][0], mat[n-threadIdx.x-1][i+2-n+threadIdx.x][1] + gapY[index_j+1], mat[n-threadIdx.x-1][i+2-n+threadIdx.x][2] + gapY[index_j+1], &b[index+(i+2-n+threadIdx.x)*(N+1)*3+(n-threadIdx.x)*3]);

                if(threadIdx.y==1)
                     mat[n-threadIdx.x][i+2-n+threadIdx.x][1]+=maximum(mat[n-threadIdx.x][i+2-n+threadIdx.x-1][0] + gapX[index_i+1], mat[n-threadIdx.x][i+2-n+threadIdx.x-1][1], mat[n-threadIdx.x][i+2-n+threadIdx.x-1][2] + gapX[index_i+1], &b[index+(i+2-n+threadIdx.x)*(N+1)*3+(n-threadIdx.x)*3+1]);

                if(threadIdx.y==2)
                     mat[n-threadIdx.x][i+2-n+threadIdx.x][2]+=maximum(mat[n-threadIdx.x-1][i+2-n+threadIdx.x-1][0] , mat[n-threadIdx.x-1][i+2-n+threadIdx.x-1][1], mat[n-threadIdx.x-1][i+2-n+threadIdx.x-1][2], &b[index+(i+2-n+threadIdx.x)*(N+1)*3+(n-threadIdx.x)*3+2]);     
            }
        }
    }
 
__syncthreads();

    for(int j=0;j<n+1;j++)
    {
        a[index+3*(N+1)*j+3*threadIdx.x+threadIdx.y]=mat[threadIdx.x][j][threadIdx.y];
    }
}

void recur(gth_Arr arr)
{
//assuming two strings of equal length
int N=arr.lenX;
int size=(N+1)*(N+1)*3*sizeof(int); 

//Creating reduntant copies of data, copying to 1D array explicitly
int *data_copy=(int*)malloc(size);
int *path_copy=(int*)malloc(size);
for(int i=0;i<=arr.lenX;i++)
{
    for(int j=0;j<=arr.lenY;j++){
        for(int k=0;k<3;k++)
        {    data_copy[j*(arr.lenX+1)*3+i*3+k]=arr.data[i][j][k];
             path_copy[j*(arr.lenX+1)*3+i*3+k]=arr.path[i][j][k];
        }
       }  
}

int *gapX=(int*)malloc((N+1)*sizeof(int));
int *gapY=(int*)malloc((N+1)*sizeof(int));

for(int i=0;i<N+1;i++)
 {
    gapX[i]=arr.gapX[i];
    gapY[i]=arr.gapY[i];
}

//Creating variables for device memory
 int *dev_a, *dev_b, *dev_gapX, *dev_gapY;
 int n=2;

    cudaMalloc(&dev_a, size);
    cudaMalloc(&dev_b, size);
    cudaMalloc(&dev_gapX, (N+1)*sizeof(int));
    cudaMalloc(&dev_gapY, (N+1)*sizeof(int));

    cudaMemcpy(dev_a,data_copy,size,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_b, path_copy,size,  cudaMemcpyHostToDevice);
    cudaMemcpy(dev_gapX, gapX, (N+1)*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_gapY, gapY, (N+1)*sizeof(int), cudaMemcpyHostToDevice);

    dim3 thread_count(n+1,3);

    struct stopwatch_t* timer=NULL;
    long double t_gpu;
    stopwatch_init();

    timer=stopwatch_create();
    stopwatch_start(timer);
    for(int i=0;i<(2*N-1)/n;i++)
    {
        if(i<N/2)
        {
            block_max<<<i+1,thread_count>>>(dev_a,i,N,dev_b,dev_gapX,dev_gapY,n);
        }
       else
        {
            block_max<<<(2*N-1)/n-i-N%n,thread_count>>>(dev_a,i,N,dev_b,dev_gapX,dev_gapY,n);
        }
    }
    t_gpu=stopwatch_stop(timer);
    cout<<"GPU: "<<t_gpu;
    cudaMemcpy(data_copy,dev_a,size,cudaMemcpyDeviceToHost);
    cudaMemcpy(path_copy,dev_b,size,cudaMemcpyDeviceToHost);
    cudaFree(dev_a);
    cudaFree(dev_b);
    cudaFree(dev_gapX);
    cudaFree(dev_gapY);

for(int i=0;i<=N;i++)
{
    for(int j=0;j<=N;j++){
        for(int k=0;k<3;k++)
        {   arr.data[i][j][k]= data_copy[j*(arr.lenX+1)*3+i*3+k];
            arr.path[i][j][k]=path_copy[j*(arr.lenX+1)*3+i*3+k];
        }
       }
}

if(N%2!=0)
{
	for(int j=1;j<N;j++)
	{
		arr.data[N][j][X] += maximum2(
                arr.data[N-1][j  ][X],
                arr.data[N-1][j  ][Y] + arr.gapY[j],
                arr.data[N-1][j  ][Z] + arr.gapY[j],
            &arr.path[N][j][X]);

            arr.data[N][j][Y] += maximum2(
                arr.data[N  ][j-1][X] + arr.gapX[N],
                arr.data[N  ][j-1][Y],
                arr.data[N  ][j-1][Z] + arr.gapX[N],
            &arr.path[N][j][Y]);

            arr.data[N][j][Z] += maximum2(
                arr.data[N-1][j-1][X],
                arr.data[N-1][j-1][Y],
                arr.data[N-1][j-1][Z],
            &arr.path[N][j][Z]);

	    arr.data[j][N][X] += maximum2(
                arr.data[j-1][N  ][X],
                arr.data[j-1][N  ][Y] + arr.gapY[N],
                arr.data[j-1][N  ][Z] + arr.gapY[N],
            &arr.path[j][N][X]);

            arr.data[j][N][Y] += maximum2(
                arr.data[j  ][N-1][X] + arr.gapX[j],
                arr.data[j  ][N-1][Y],
                arr.data[j  ][N-1][Z] + arr.gapX[j],
            &arr.path[j][N][Y]);

            arr.data[j][N][Z] += maximum2(
                arr.data[j-1][N-1][X],
                arr.data[j-1][N-1][Y],
                arr.data[j-1][N-1][Z],
            &arr.path[j][N][Z]);

	}
	int i=N,j=N;
	arr.data[i][j][X] += maximum2(
                arr.data[i-1][j  ][X],
                arr.data[i-1][j  ][Y] + arr.gapY[j],
                arr.data[i-1][j  ][Z] + arr.gapY[j],
            &arr.path[i][j][X]);

            arr.data[i][j][Y] += maximum2(
                arr.data[i  ][j-1][X] + arr.gapX[i],
                arr.data[i  ][j-1][Y],
                arr.data[i  ][j-1][Z] + arr.gapX[i],
            &arr.path[i][j][Y]);

            arr.data[i][j][Z] += maximum2(
                arr.data[i-1][j-1][X],
                arr.data[i-1][j-1][Y],
                arr.data[i-1][j-1][Z],
            &arr.path[i][j][Z]);
}

}

int gth_align(gth_Arr arr) {

    for (size_t i=1 ; i<=arr.lenX ; i++) {
        arr.gapX[i] = arr.data[i][0][Y];
        arr.data[i][0][X] += arr.data[i-1][0][X];
    }

    for (size_t j=1 ; j<=arr.lenY ; j++) {
        arr.gapY[j] = arr.data[0][j][X];
        arr.data[0][j][Y] += arr.data[0][j-1][Y];
    }

    fill_arr(X, arr.path, X, 1, arr.lenX, 0,        0);
    fill_arr(Y, arr.path, Y, 0,        0, 1, arr.lenY);

    fill_arr(-Inf, arr.data, Y, 1, arr.lenX, 0,        0);
    fill_arr(-Inf, arr.data, Z, 1, arr.lenX, 0,        0);
    fill_arr(-Inf, arr.data, X, 0,        0, 1, arr.lenY);
    fill_arr(-Inf, arr.data, Z, 0,        0, 1, arr.lenY);

    recur(arr);

    for (size_t i=0 ; i<=arr.lenX ; i++) arr.gapX[i] = 0;
    for (size_t j=0 ; j<=arr.lenY ; j++) arr.gapY[j] = 0;

    // backtracking
    int K;
    int max = maximum2(
        arr.data[arr.lenX][arr.lenY][X],
        arr.data[arr.lenX][arr.lenY][Y],
        arr.data[arr.lenX][arr.lenY][Z],
    &K);

    size_t i = arr.lenX, j = arr.lenY;
    while (i > 0 || j > 0) {
        if (K == X) {
            arr.gapY[j]++;
            K = arr.path[i--][j  ][X];
        }
        else if (K == Y) {
            arr.gapX[i]++;
            K = arr.path[i  ][j--][Y];
        }
        else if (K == Z) {
            K = arr.path[i--][j--][Z];
        }

    }

    return max;
}

// default BLOSUM62 matrix from the EMBOSS package. See end of the file for its initialization.
//const gth_Sub BLOSUM62;
const gth_Sub BLOSUM62 = {
    .alpha = "ARNDCQEGHILKMFPSTWYVBZX", .score = 
    {{  4, -2,  0, -2, -1, -2,  0, -2, -1, -4, -1, -1, -1, -2, -4, -1, -1, -1,  1,  0, -4,  0, -3,  0, -2, -1},
     { -2,  4, -3,  4,  1, -3, -1,  0, -3, -4,  0, -4, -3,  3, -4, -2,  0, -1,  0, -1, -4, -3, -4, -1, -3,  1},
     {  0, -3,  9, -3, -4, -2, -3, -3, -1, -4, -3, -1, -1, -3, -4, -3, -3, -3, -1, -1, -4, -1, -2, -2, -2, -3},
     { -2,  4, -3,  6,  2, -3, -1, -1, -3, -4, -1, -4, -3,  1, -4, -1,  0, -2,  0, -1, -4, -3, -4, -1, -3,  1},
     { -1,  1, -4,  2,  5, -3, -2,  0, -3, -4,  1, -3, -2,  0, -4, -1,  2,  0,  0, -1, -4, -2, -3, -1, -2,  4},
     { -2, -3, -2, -3, -3,  6, -3, -1,  0, -4, -3,  0,  0, -3, -4, -4, -3, -3, -2, -2, -4, -1,  1, -1,  3, -3},
     {  0, -1, -3, -1, -2, -3,  6, -2, -4, -4, -2, -4, -3,  0, -4, -2, -2, -2,  0, -2, -4, -3, -2, -1, -3, -2},
     { -2,  0, -3, -1,  0, -1, -2,  8, -3, -4, -1, -3, -2,  1, -4, -2,  0,  0, -1, -2, -4, -3, -2, -1,  2,  0},
     { -1, -3, -1, -3, -3,  0, -4, -3,  4, -4, -3,  2,  1, -3, -4, -3, -3, -3, -2, -1, -4,  3, -3, -1, -1, -3},
     { -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4},
     { -1,  0, -3, -1,  1, -3, -2, -1, -3, -4,  5, -2, -1,  0, -4, -1,  1,  2,  0, -1, -4, -2, -3, -1, -2,  1},
     { -1, -4, -1, -4, -3,  0, -4, -3,  2, -4, -2,  4,  2, -3, -4, -3, -2, -2, -2, -1, -4,  1, -2, -1, -1, -3},
     { -1, -3, -1, -3, -2,  0, -3, -2,  1, -4, -1,  2,  5, -2, -4, -2,  0, -1, -1, -1, -4,  1, -1, -1, -1, -1},
     { -2,  3, -3,  1,  0, -3,  0,  1, -3, -4,  0, -3, -2,  6, -4, -2,  0,  0,  1,  0, -4, -3, -4, -1, -2,  0},
     { -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4},
     { -1, -2, -3, -1, -1, -4, -2, -2, -3, -4, -1, -3, -2, -2, -4,  7, -1, -2, -1, -1, -4, -2, -4, -2, -3, -1},
     { -1,  0, -3,  0,  2, -3, -2,  0, -3, -4,  1, -2,  0,  0, -4, -1,  5,  1,  0, -1, -4, -2, -2, -1, -1,  3},
     { -1, -1, -3, -2,  0, -3, -2,  0, -3, -4,  2, -2, -1,  0, -4, -2,  1,  5, -1, -1, -4, -3, -3, -1, -2,  0},
     {  1,  0, -1,  0,  0, -2,  0, -1, -2, -4,  0, -2, -1,  1, -4, -1,  0, -1,  4,  1, -4, -2, -3,  0, -2,  0},
     {  0, -1, -1, -1, -1, -2, -2, -2, -1, -4, -1, -1, -1,  0, -4, -1, -1, -1,  1,  5, -4,  0, -2,  0, -2, -1},
     { -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4},
     {  0, -3, -1, -3, -2, -1, -3, -3,  3, -4, -2,  1,  1, -3, -4, -2, -2, -3, -2,  0, -4,  4, -3, -1, -1, -2},
     { -3, -4, -2, -4, -3,  1, -2, -2, -3, -4, -3, -2, -1, -4, -4, -4, -2, -3, -3, -2, -4, -3, 11, -2,  2, -3},
     {  0, -1, -2, -1, -1, -1, -1, -1, -1, -4, -1, -1, -1, -1, -4, -2, -1, -1,  0,  0, -4, -1, -2, -1, -1, -1},
     { -2, -3, -2, -3, -2,  3, -3,  2, -1, -4, -2, -1, -1, -2, -4, -3, -1, -2, -2, -2, -4, -1,  2, -1,  7, -2},
     { -1,  1, -3,  1,  4, -3, -2,  0, -3, -4,  1, -3, -1,  0, -4, -1,  3,  0,  0, -1, -4, -2, -3, -1, -2,  4}}
};


int main(int argc, char **argv) {

    // arguments parsing
    
    int error = 0;
    int dump = 0;
    int quiet = 0;
    char *seqX_path = NULL;
    char *seqY_path = NULL;
    char *matrix_path = NULL;
    char *arr_paths[] = {NULL, NULL, NULL};
    double gapopen = 9.5;
    double gapextend = 0.5;
    double endopen = 0.0;
    double endextend = 0.0;

    for (size_t i=1 ; i<argc && !error ; i++) {
        if (!strcmp(argv[i], "-dump")) dump = 1;
        else if (!strcmp(argv[i], "-quiet")) quiet = 1;
        else if (argv[i][0] == '-') {
            if (i+1 >= argc) {
                fprintf(stderr, "ERROR: missing argument for option '%s'\n\n", argv[i]);
                error = 1;
            }
            else if (!strcmp(argv[i], "-matrix"))    matrix_path  = argv[i+1];
            else if (!strcmp(argv[i], "-gapopen"))   gapopen      = atof(argv[i+1]);
            else if (!strcmp(argv[i], "-gapextend")) gapextend    = atof(argv[i+1]);
            else if (!strcmp(argv[i], "-endopen"))   endopen      = atof(argv[i+1]);
            else if (!strcmp(argv[i], "-endextend")) endextend    = atof(argv[i+1]);
            else if (!strcmp(argv[i], "-arrxfile"))  arr_paths[0] = argv[i+1];
            else if (!strcmp(argv[i], "-arryfile"))  arr_paths[1] = argv[i+1];
            else if (!strcmp(argv[i], "-arrzfile"))  arr_paths[2] = argv[i+1];
            else {
                fprintf(stderr, "ERROR: unknown option '%s'\n\n", argv[i]);
                error = 1;
            }
            i++;
        }
        else if (seqX_path == NULL) seqX_path = argv[i];
        else if (seqY_path == NULL) seqY_path = argv[i];
        else {
            fprintf(stderr, "ERROR: too many arguments\n\n");
            error = 1;
        }
    }

    if (argc <= 1) error = 1;
    else if ((seqX_path == NULL || seqY_path == NULL) && !error) {
        fprintf(stderr, "ERROR: you must provide two sequence files\n\n");
        error = 1;
    }


    // help message

    if (error) {
        printf("Needleman-Wunsch global alignment of two sequences\n\n");
        printf("Usage: %s [OPTIONS] SEQUENCE1.fasta SEQUENCE2.fasta\n\n", argv[0]);
        printf("Options and their defaults:\n");
        printf("    -matrix     BLOSUM62    matrix file in NCBI/EMBOSS format\n");
        printf("    -gapopen    9.5         opening gap penalty\n");
        printf("    -gapextend  0.5         extending gap penalty\n");
        printf("    -endopen    0.0         opening end gap penalty\n");
        printf("    -endextend  0.0         extending end gap penalty\n");
        printf("    -quiet                  decrease verbosity\n\n");
        printf("Advanced:\n");
        printf("    -arrxfile   arrx.txt    load initial array X from this file\n");
        printf("    -arryfile   arry.txt    load initial array Y from this file\n");
        printf("    -arrzfile   arrz.txt    load initial array Z from this file\n");
        printf("    -dump                   output final arrays containing partial scores\n\n");
        printf("Default parameters are the same as those of EMBOSS Needle.\n");
        printf("Feeding corrupted file formats to the program may result in undefined behaviour\n");
        return 0;
    }


    // load sequences

    gth_Seq seqX = gth_read_fasta(seqX_path);
    if (seqX.len == 0) {
        fprintf(stderr, "ERROR: problem reading file '%s'\n", seqX_path);
        return -1;
    }

    gth_Seq seqY = gth_read_fasta(seqY_path);
    if (seqY.len == 0) {
        fprintf(stderr, "ERROR: problem reading file '%s'\n", seqY_path);
        return -1;
    }


    // load matrix

    gth_Sub matrix = BLOSUM62;
    if (matrix_path != NULL) {
        matrix = gth_read_matrix(matrix_path);
        if (matrix.alpha[0] == '\0') {
            fprintf(stderr, "ERROR: problem reading file '%s'\n", matrix_path);
            return -1;
        }
    }
    for (size_t i=0 ; i<26 ; i++) {
        for (size_t j=0 ; j<26 ; j++)
        {    matrix.score[i][j] *= 10;
        }
    }

    // create, fill, and backtrack the arrays
    gth_Arr array = gth_init(seqX.len, seqY.len);
    gth_set_sub(array, seqX.res, seqY.res, matrix.score);

    gth_set_gap(array, (int)(gapopen*10), (int)(gapextend*10), (int)(endopen*10), (int)(endextend*10));
    for (int k=0 ; k<3 ; k++) {
        if (arr_paths[k] == NULL) continue;
        FILE *file = fopen(arr_paths[k], "r");
        if (!file) {
            fprintf(stderr, "ERROR: problem reading file '%s'\n", arr_paths[k]);
            gth_free(array);
            free(seqX.res);
            free(seqY.res);
            return -1;
        }
        for (size_t i=0 ; i<=array.lenX ; i++) {
            for (size_t j=0 ; j<=array.lenY ; j++) {
                fscanf(file, "%d", &array.data[i][j][k]);
                array.data[i][j][k] *= 10;
            }
            fscanf(file, "%*[^\n]");
        }
        fclose(file);
    }
    double score = (double)(gth_align(array)) / 10;


    // output

    if (!quiet) {
        printf("# Needleman-Wunsch global alignment of two sequences\n");
        printf("#\n");
        if (arr_paths[2] == NULL) {
            printf("# matrix file: %s\n", (matrix_path == NULL) ? "none specified, using BLOSUM62" : matrix_path);
        }
        else {
            printf("# substitution scores array provided by user\n");
        }
        if (arr_paths[0] == NULL && arr_paths[1] == NULL) {
            printf("# gap opening penalty: %.1f\n", gapopen);
            printf("# gap extending penalty: %.1f\n", gapextend);
            printf("# end gap opening penalty: %.1f\n", endopen);
            printf("# end gap extending penalty: %.1f\n", endextend);
        }
        else {
            printf("# gap scores array(s) provided by user\n");
        }
        printf("#\n");
        printf("# score: %.1f\n\n", score);
    }
    if (dump) {
        for (int k=0 ; k<3 ; k++) {
            printf("# scores in array %c:\n", 'X'+k);
            for (size_t i=0 ; i<=array.lenX ; i++) {
                for (size_t j=0 ; j<=array.lenY ; j++) {
                    printf("% 15d ", array.data[i][j][k]);
                }
                printf("\n");
            }
            printf("\n");
        }
    }
    printf(">%s\n", seqX.name);
    gth_putseq(stdout, seqX.res, array.gapX);
    printf("\n\n>%s\n", seqY.name);
    gth_putseq(stdout, seqY.res, array.gapY);
    printf("\n");


    // cleanup

    gth_free(array);
    free(seqX.res);
    free(seqY.res);
    return 0;
}


// default BLOSUM62 matrix from the EMBOSS package.


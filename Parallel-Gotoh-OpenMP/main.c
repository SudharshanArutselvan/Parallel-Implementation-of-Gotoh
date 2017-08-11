#include <stdio.h>
#include "gotoh.h"
#include <time.h>
#include <omp.h>
#include "timer.c"
#include "timer.h"

// default BLOSUM62 matrix from the EMBOSS package. See end of the file for its initialization.
const gth_Sub BLOSUM62;


int main(int argc, char **argv) {
	
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
	gth_Seq seqY;
	gth_Seq seqX;
	double score;
	gth_Arr array;
	gth_Sub matrix;
	
	
    for (int i=1 ; i<argc && !error ; i++) {
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

    
	seqX = gth_read_fasta(seqX_path);
	printf("Size of X:%d\n",seqX.len);
	if (seqX.len == 0) {
		fprintf(stderr, "ERROR: problem reading file '%s'\n", seqX_path);
		exit(0);
	}
		
	seqY = gth_read_fasta(seqY_path);
	printf("Size of Y:%d\n",seqY.len);
	if (seqY.len == 0) {
		fprintf(stderr, "ERROR: problem reading file '%s'\n", seqY_path);
		exit(0);
	}
	
	matrix = BLOSUM62;
	array = gth_init(seqX.len, seqY.len);
		
	gth_set_sub(array, seqX.res, seqY.res, matrix.score);
		
	gth_set_gap(array, (int)(gapopen*10), (int)(gapextend*10), (int)(endopen*10), (int)(endextend*10));
		
	score = (double)(gth_align(array)) / 10;
		
	printf("# score: %.1f\n\n", score);
	printf("\n");
			
    // cleanup

    gth_free(array);
    free(seqX.res);
    free(seqY.res);
    return 0;
}


// default BLOSUM62 matrix from the EMBOSS package.
const gth_Sub BLOSUM62 = {
    .alpha = "ARNDCQEGHILKMFPSTWYVBZX", .score = 
    {{  40, -20,  00, -20, -10, -20,  00, -20, -10, -40, -10, -10, -10, -20, -40, -10, -10, -10,  10,  00, -40,  00, -30,  00, -20, -10},
     { -20,  40, -30,  40,  10, -30, -10,  00, -30, -40,  00, -40, -30,  30, -40, -20,  00, -10,  00, -10, -40, -30, -40, -10, -30,  10},
     {  00, -30,  90, -30, -40, -20, -30, -30, -10, -40, -30, -10, -10, -30, -40, -30, -30, -30, -10, -10, -40, -10, -20, -20, -20, -30},
     { -20,  40, -30,  60,  20, -30, -10, -10, -30, -40, -10, -40, -30,  10, -40, -10,  00, -20,  00, -10, -40, -30, -40, -10, -30,  10},
     { -10,  10, -40,  20,  50, -30, -20,  00, -30, -40,  10, -30, -20,  00, -40, -10,  20,  00,  00, -10, -40, -20, -30, -10, -20,  40},
     { -20, -30, -20, -30, -30,  60, -30, -10,  00, -40, -30,  00,  00, -30, -40, -40, -30, -30, -20, -20, -40, -10,  10, -10,  30, -30},
     {  00, -10, -30, -10, -20, -30,  60, -20, -40, -40, -20, -40, -30,  00, -40, -20, -20, -20,  00, -20, -40, -30, -20, -10, -30, -20},
     { -20,  00, -30, -10,  00, -10, -20,  80, -30, -40, -10, -30, -20,  10, -40, -20,  00,  00, -10, -20, -40, -30, -20, -10,  20,  00},
     { -10, -30, -10, -30, -30,  00, -40, -30,  40, -40, -30,  20,  10, -30, -40, -30, -30, -30, -20, -10, -40,  30, -30, -10, -10, -30},
     { -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40},
     { -10,  00, -30, -10,  10, -30, -20, -10, -30, -40,  50, -20, -10,  00, -40, -10,  10,  20,  00, -10, -40, -20, -30, -10, -20,  10},
     { -10, -40, -10, -40, -30,  00, -40, -30,  20, -40, -20,  40,  20, -30, -40, -30, -20, -20, -20, -10, -40,  10, -20, -10, -10, -30},
     { -10, -30, -10, -30, -20,  00, -30, -20,  10, -40, -10,  20,  50, -20, -40, -20,  00, -10, -10, -10, -40,  10, -10, -10, -10, -10},
     { -20,  30, -30,  10,  00, -30,  00,  10, -30, -40,  00, -30, -20,  60, -40, -20,  00,  00,  10,  00, -40, -30, -40, -10, -20,  00},
     { -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40},
     { -10, -20, -30, -10, -10, -40, -20, -20, -30, -40, -10, -30, -20, -20, -40,  70, -10, -20, -10, -10, -40, -20, -40, -20, -30, -10},
     { -10,  00, -30,  00,  20, -30, -20,  00, -30, -40,  10, -20,  00,  00, -40, -10,  50,  10,  00, -10, -40, -20, -20, -10, -10,  30},
     { -10, -10, -30, -20,  00, -30, -20,  00, -30, -40,  20, -20, -10,  00, -40, -20,  10,  50, -10, -10, -40, -30, -30, -10, -20,  00},
     {  10,  00, -10,  00,  00, -20,  00, -10, -20, -40,  00, -20, -10,  10, -40, -10,  00, -10,  40,  10, -40, -20, -30,  00, -20,  00},
     {  00, -10, -10, -10, -10, -20, -20, -20, -10, -40, -10, -10, -10,  00, -40, -10, -10, -10,  10,  50, -40,  00, -20,  00, -20, -10},
     { -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40, -40},
     {  00, -30, -10, -30, -20, -10, -30, -30,  30, -40, -20,  10,  10, -30, -40, -20, -20, -30, -20,  00, -40,  40, -30, -10, -10, -20},
     { -30, -40, -20, -40, -30,  10, -20, -20, -30, -40, -30, -20, -10, -40, -40, -40, -20, -30, -30, -20, -40, -30, 110, -20,  20, -30},
     {  00, -10, -20, -10, -10, -10, -10, -10, -10, -40, -10, -10, -10, -10, -40, -20, -10, -10,  00,  00, -40, -10, -20, -10, -10, -10},
     { -20, -30, -20, -30, -20,  30, -30,  20, -10, -40, -20, -10, -10, -20, -40, -30, -10, -20, -20, -20, -40, -10,  20, -10,  70, -20},
     { -10,  10, -30,  10,  40, -30, -20,  00, -30, -40,  10, -30, -10,  00, -40, -10,  30,  00,  00, -10, -40, -20, -30, -10, -20,  40}}
};


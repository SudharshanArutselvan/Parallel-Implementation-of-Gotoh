# Parallel-Gotoh

Parallel Implementation of the Gotoh Algorithm using OpenMP

To compile the program:

gcc -std=c99 -O3 -fopenmp gotoh.c main.c -o align

To run the OpenMP Gotoh Algorithm:

If running in the HPC cluster run "qsub gotoh.sh"
	This command runs the program 5 times and the number of cores can be changed in it.

If running in another machine, run
	./align sequence1.fasta sequence2.fasta